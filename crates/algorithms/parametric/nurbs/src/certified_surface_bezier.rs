//! Interval-aware tensor-product NURBS refinement for certified surface queries.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Scalar;
use axiolid_surface::BSplineSurface;

use crate::{
    certified_bezier::{Cell, HomogeneousPoint, Interval},
    certified_refinement::{insert_homogeneous_once, RefinementBudget},
    periodic_surface::{exact_scalar_add, exact_scalar_subtract},
};

#[derive(Debug, Clone)]
pub(crate) struct Patch {
    pub(crate) controls: Vec<Vec<HomogeneousPoint>>,
    pub(crate) u_start: Scalar,
    pub(crate) u_end: Scalar,
    pub(crate) v_start: Scalar,
    pub(crate) v_end: Scalar,
}

#[derive(Debug, Clone, Copy)]
struct AxisPlan {
    degree: usize,
    domain: (Scalar, Scalar),
    periodic: bool,
    initial_controls: usize,
    unique_controls: usize,
    segments: usize,
    refined_controls: usize,
    work_per_line: u128,
}

impl AxisPlan {
    fn required_insertions(self, knot: Scalar, multiplicity: u32) -> GeomResult<usize> {
        if knot < self.domain.0 || knot > self.domain.1 {
            return Ok(0);
        }
        let multiplicity = usize::try_from(multiplicity).map_err(|_| allocation_error())?;
        if multiplicity > self.degree {
            if !self.periodic
                && (knot == self.domain.0 || knot == self.domain.1)
                && multiplicity == self.degree + 1
            {
                return Ok(0);
            }
            return Err(GeomError::InvalidInput(
                "certified surface refinement active knot multiplicity exceeds its degree"
                    .to_owned(),
            ));
        }
        Ok(self.degree - multiplicity)
    }
}

fn allocation_error() -> GeomError {
    GeomError::BudgetExceeded {
        resource: "certified surface refinement allocation",
    }
}

fn work_overflow() -> GeomError {
    GeomError::BudgetExceeded {
        resource: "certified surface refinement work arithmetic",
    }
}

fn reserve<T>(target: &mut Vec<T>, additional: usize) -> GeomResult<()> {
    target
        .try_reserve_exact(additional)
        .map_err(|_| allocation_error())
}

fn clone_controls(input: &[HomogeneousPoint]) -> GeomResult<Vec<HomogeneousPoint>> {
    let mut output = Vec::new();
    reserve(&mut output, input.len())?;
    output.extend(input.iter().cloned());
    Ok(output)
}

fn expanded_knot_at(knots: &[Scalar], multiplicities: &[u32], target: usize) -> Option<Scalar> {
    let mut offset = 0usize;
    for (&knot, &multiplicity) in knots.iter().zip(multiplicities) {
        let multiplicity = usize::try_from(multiplicity).ok()?;
        let end = offset.checked_add(multiplicity)?;
        if target < end {
            return Some(knot);
        }
        offset = end;
    }
    None
}

fn axis_plan(
    degree: u16,
    control_count: usize,
    knots: &[Scalar],
    multiplicities: &[u32],
    periodic: bool,
    name: &str,
) -> GeomResult<AxisPlan> {
    if knots.len() < 2 || knots.len() != multiplicities.len() {
        return Err(GeomError::InvalidInput(format!(
            "B-spline surface {name} knots and multiplicities must be nonempty and aligned"
        )));
    }
    if knots.iter().any(|knot| !knot.is_finite()) || knots.windows(2).any(|pair| pair[0] >= pair[1])
    {
        return Err(GeomError::InvalidInput(format!(
            "B-spline surface {name} knots must be finite and strictly increasing"
        )));
    }
    let degree = usize::from(degree);
    let order = degree.checked_add(1).ok_or_else(allocation_error)?;
    let order_u32 = u32::try_from(order).map_err(|_| allocation_error())?;
    if control_count < order {
        return Err(GeomError::InvalidInput(format!(
            "B-spline surface {name} degree exceeds its control count"
        )));
    }
    if !periodic
        && (multiplicities.first().copied() != Some(order_u32)
            || multiplicities.last().copied() != Some(order_u32))
    {
        return Err(GeomError::InvalidInput(format!(
            "certified surface refinement requires a clamped {name} axis"
        )));
    }
    if multiplicities.iter().any(|&m| m == 0 || m > order_u32) {
        return Err(GeomError::InvalidInput(format!(
            "B-spline surface {name} multiplicity is outside 1..={order}"
        )));
    }
    if !periodic
        && multiplicities[1..multiplicities.len() - 1]
            .iter()
            .any(|&m| usize::try_from(m).map_or(true, |m| m > degree))
    {
        return Err(GeomError::InvalidInput(format!(
            "certified surface refinement does not support a full-multiplicity internal {name} knot; internal multiplicity must be between one and the degree"
        )));
    }
    let expanded = multiplicities
        .iter()
        .try_fold(0_usize, |sum, &value| sum.checked_add(value as usize));
    let expected = control_count
        .checked_add(degree)
        .and_then(|value| value.checked_add(1));
    if expanded != expected {
        return Err(GeomError::InvalidInput(format!(
            "B-spline surface {name} compact knot count is inconsistent"
        )));
    }
    let domain = (
        expanded_knot_at(knots, multiplicities, degree).ok_or_else(allocation_error)?,
        expanded_knot_at(knots, multiplicities, control_count).ok_or_else(allocation_error)?,
    );
    if !domain.0.is_finite() || !domain.1.is_finite() || domain.1 <= domain.0 {
        return Err(GeomError::InvalidInput(format!(
            "B-spline surface {name} active domain must be finite and positive"
        )));
    }
    let segments = knots
        .windows(2)
        .filter(|pair| pair[0] >= domain.0 && pair[1] <= domain.1)
        .count();
    if segments == 0 {
        return Err(GeomError::InvalidInput(format!(
            "B-spline surface {name} active domain has no spans"
        )));
    }
    let mut insertions = 0usize;
    for (&knot, &multiplicity) in knots.iter().zip(multiplicities) {
        if knot < domain.0 || knot > domain.1 {
            continue;
        }
        let multiplicity = usize::try_from(multiplicity).map_err(|_| {
            GeomError::InvalidInput(format!(
                "B-spline surface {name} multiplicity does not fit usize"
            ))
        })?;
        let missing = if multiplicity <= degree {
            degree - multiplicity
        } else if !periodic && (knot == domain.0 || knot == domain.1) && multiplicity == order {
            0
        } else {
            return Err(GeomError::InvalidInput(format!(
                "B-spline surface {name} active multiplicity exceeds its degree"
            )));
        };
        insertions = insertions
            .checked_add(missing)
            .ok_or_else(allocation_error)?;
    }
    let unique_controls = if periodic {
        control_count
            .checked_sub(degree)
            .filter(|&count| count > degree)
            .ok_or_else(allocation_error)?
    } else {
        control_count
    };
    let extension = if periodic {
        degree.checked_mul(2).ok_or_else(allocation_error)?
    } else {
        0
    };
    let extended_controls = control_count
        .checked_add(extension)
        .ok_or_else(allocation_error)?;
    let refined_controls = extended_controls
        .checked_add(insertions)
        .ok_or_else(allocation_error)?;
    let extended_knots = expanded
        .ok_or_else(allocation_error)?
        .checked_add(extension)
        .ok_or_else(allocation_error)?;
    let mut work = u128::try_from(extended_knots).map_err(|_| allocation_error())?;
    let extension_work = u128::try_from(extension).map_err(|_| allocation_error())?;
    work = work
        .checked_add(extension_work)
        .ok_or_else(allocation_error)?;
    for insertion in 0..insertions {
        let controls_before = extended_controls
            .checked_add(insertion)
            .ok_or_else(allocation_error)?;
        work = work
            .checked_add(controls_before as u128)
            .and_then(|value| value.checked_add((degree + 2) as u128))
            .ok_or_else(allocation_error)?;
    }
    Ok(AxisPlan {
        degree,
        domain,
        periodic,
        initial_controls: control_count,
        unique_controls,
        segments,
        refined_controls,
        work_per_line: work,
    })
}

fn validate_surface(
    surface: &BSplineSurface,
    u_periodic: bool,
    v_periodic: bool,
) -> GeomResult<(usize, usize, AxisPlan, AxisPlan)> {
    let rows = surface.control_points.len();
    let Some(first_row) = surface.control_points.first() else {
        return Err(GeomError::InvalidInput(
            "B-spline surface has no control points".to_owned(),
        ));
    };
    let cols = first_row.len();
    if cols == 0 || surface.control_points.iter().any(|row| row.len() != cols) {
        return Err(GeomError::InvalidInput(
            "B-spline surface control net must be nonempty and rectangular".to_owned(),
        ));
    }
    if surface
        .control_points
        .iter()
        .flatten()
        .any(|point| !point.is_finite())
    {
        return Err(GeomError::InvalidInput(
            "B-spline surface control points must be finite".to_owned(),
        ));
    }
    if let Some(weights) = &surface.weights {
        if weights.len() != rows || weights.iter().any(|row| row.len() != cols) {
            return Err(GeomError::InvalidInput(
                "B-spline surface weight net does not match the control net".to_owned(),
            ));
        }
        if weights
            .iter()
            .flatten()
            .any(|weight| !weight.is_finite() || *weight <= 0.0)
        {
            return Err(GeomError::InvalidInput(
                "B-spline surface weights must be finite and strictly positive".to_owned(),
            ));
        }
    }
    let u = axis_plan(
        surface.u_degree,
        rows,
        &surface.u_knots,
        &surface.u_multiplicities,
        u_periodic,
        "u",
    )?;
    let v = axis_plan(
        surface.v_degree,
        cols,
        &surface.v_knots,
        &surface.v_multiplicities,
        v_periodic,
        "v",
    )?;
    Ok((rows, cols, u, v))
}

fn homogeneous_control(
    surface: &BSplineSurface,
    u: usize,
    v: usize,
) -> GeomResult<HomogeneousPoint> {
    let weight = surface.weights.as_ref().map_or(1.0, |rows| rows[u][v]);
    let point = surface.control_points[u][v];
    Ok(HomogeneousPoint {
        numerator: [
            Interval::exact(point.x)?.multiply(Interval::exact(weight)?)?,
            Interval::exact(point.y)?.multiply(Interval::exact(weight)?)?,
            Interval::exact(point.z)?.multiply(Interval::exact(weight)?)?,
        ],
        weight: Interval::exact(weight)?,
    })
}

fn expanded_knots(knots: &[Scalar], multiplicities: &[u32]) -> GeomResult<Vec<Scalar>> {
    let count = multiplicities
        .iter()
        .try_fold(0_usize, |sum, &value| sum.checked_add(value as usize));
    let mut expanded = Vec::new();
    reserve(&mut expanded, count.ok_or_else(allocation_error)?)?;
    for (&knot, &multiplicity) in knots.iter().zip(multiplicities) {
        expanded.extend(std::iter::repeat_n(knot, multiplicity as usize));
    }
    Ok(expanded)
}

fn cyclically_extend_controls(
    plan: AxisPlan,
    controls: Vec<HomogeneousPoint>,
) -> GeomResult<Vec<HomogeneousPoint>> {
    if !plan.periodic {
        return Ok(controls);
    }
    if controls.len() != plan.initial_controls {
        return Err(GeomError::Degenerate(
            "periodic control line length changed before cyclic extension".to_owned(),
        ));
    }
    let extension = plan.degree.checked_mul(2).ok_or_else(allocation_error)?;
    let target = controls
        .len()
        .checked_add(extension)
        .ok_or_else(allocation_error)?;
    let mut extended = Vec::new();
    reserve(&mut extended, target)?;
    let prefix = plan
        .unique_controls
        .checked_sub(plan.degree)
        .ok_or_else(allocation_error)?;
    extended.extend(controls[prefix..plan.unique_controls].iter().cloned());
    extended.extend(controls.iter().cloned());
    for offset in 0..plan.degree {
        let source = plan
            .degree
            .checked_add(offset)
            .ok_or_else(allocation_error)?
            % plan.unique_controls;
        extended.push(controls[source].clone());
    }
    Ok(extended)
}

fn cyclically_extend_knots(plan: AxisPlan, expanded: &[Scalar]) -> GeomResult<Vec<Scalar>> {
    if !plan.periodic {
        let mut cloned = Vec::new();
        reserve(&mut cloned, expanded.len())?;
        cloned.extend(expanded.iter().copied());
        return Ok(cloned);
    }
    let extension = plan.degree.checked_mul(2).ok_or_else(allocation_error)?;
    let target = expanded
        .len()
        .checked_add(extension)
        .ok_or_else(allocation_error)?;
    let period = exact_scalar_subtract(plan.domain.1, plan.domain.0).ok_or_else(|| {
        GeomError::Degenerate("periodic knot period is not exact in binary64".to_owned())
    })?;
    if period <= 0.0 {
        return Err(GeomError::Degenerate(
            "periodic knot extension lost its positive period".to_owned(),
        ));
    }
    let mut extended = Vec::new();
    reserve(&mut extended, target)?;
    let prefix_source = plan
        .unique_controls
        .checked_sub(plan.degree)
        .ok_or_else(allocation_error)?;
    for offset in 0..plan.degree {
        extended.push(
            exact_scalar_subtract(expanded[prefix_source + offset], period).ok_or_else(|| {
                GeomError::Degenerate(
                    "periodic prefix knot translation is not exact in binary64".to_owned(),
                )
            })?,
        );
    }
    extended.extend(expanded.iter().copied());
    let suffix_source = plan
        .degree
        .checked_mul(2)
        .and_then(|value| value.checked_add(1))
        .ok_or_else(allocation_error)?;
    for offset in 0..plan.degree {
        extended.push(
            exact_scalar_add(expanded[suffix_source + offset], period).ok_or_else(|| {
                GeomError::Degenerate(
                    "periodic suffix knot translation is not exact in binary64".to_owned(),
                )
            })?,
        );
    }
    if extended.iter().any(|value| !value.is_finite())
        || extended.windows(2).any(|pair| pair[1] < pair[0])
    {
        return Err(GeomError::Degenerate(
            "periodic knot extension is non-finite or unordered".to_owned(),
        ));
    }
    Ok(extended)
}

fn refined_axis_knots(
    plan: AxisPlan,
    knots: &[Scalar],
    multiplicities: &[u32],
) -> GeomResult<Vec<Scalar>> {
    let base = expanded_knots(knots, multiplicities)?;
    let mut expanded = cyclically_extend_knots(plan, &base)?;
    let mut insertions = 0usize;
    for (&knot, &multiplicity) in knots.iter().zip(multiplicities) {
        insertions = insertions
            .checked_add(plan.required_insertions(knot, multiplicity)?)
            .ok_or_else(allocation_error)?;
    }
    reserve(&mut expanded, insertions)?;
    for (&knot, &multiplicity) in knots.iter().zip(multiplicities) {
        for _ in 0..plan.required_insertions(knot, multiplicity)? {
            let position = expanded.partition_point(|&existing| existing <= knot);
            expanded.insert(position, knot);
        }
    }
    Ok(expanded)
}

fn span_control_first(expanded: &[Scalar], degree: usize, start: Scalar) -> GeomResult<usize> {
    let after_start = expanded.partition_point(|&knot| knot <= start);
    let span = after_start.checked_sub(1).ok_or_else(allocation_error)?;
    span.checked_sub(degree).ok_or_else(|| {
        GeomError::Degenerate("refined periodic span escaped its control support".to_owned())
    })
}

fn refine_line(
    controls: Vec<HomogeneousPoint>,
    plan: AxisPlan,
    knots: &[Scalar],
    multiplicities: &[u32],
) -> GeomResult<Vec<HomogeneousPoint>> {
    let base_knots = expanded_knots(knots, multiplicities)?;
    let mut expanded = cyclically_extend_knots(plan, &base_knots)?;
    let mut controls = cyclically_extend_controls(plan, controls)?;
    for (&knot, &multiplicity) in knots.iter().zip(multiplicities) {
        let needed = plan.required_insertions(knot, multiplicity)?;
        for _ in 0..needed {
            (controls, expanded) = insert_homogeneous_once(controls, expanded, plan.degree, knot)?;
        }
    }
    Ok(controls)
}

fn collect_original_net(
    surface: &BSplineSurface,
    rows: usize,
    cols: usize,
) -> GeomResult<Vec<Vec<HomogeneousPoint>>> {
    let mut net = Vec::new();
    reserve(&mut net, rows)?;
    for u in 0..rows {
        let mut row = Vec::new();
        reserve(&mut row, cols)?;
        for v in 0..cols {
            row.push(homogeneous_control(surface, u, v)?);
        }
        net.push(row);
    }
    Ok(net)
}

impl Patch {
    pub(crate) fn midpoint_point(&self) -> GeomResult<HomogeneousPoint> {
        let mut collapsed_v = Vec::new();
        collapsed_v
            .try_reserve(self.controls.len())
            .map_err(|_| allocation_error())?;
        for row in &self.controls {
            collapsed_v.push(
                Cell {
                    controls: clone_controls(row)?,
                    start: self.v_start,
                    end: self.v_end,
                    depth: 0,
                }
                .midpoint_point()?,
            );
        }
        Cell {
            controls: collapsed_v,
            start: self.u_start,
            end: self.u_end,
            depth: 0,
        }
        .midpoint_point()
    }

    pub(crate) fn point_at(&self, u: Scalar, v: Scalar) -> GeomResult<HomogeneousPoint> {
        if !u.is_finite()
            || !v.is_finite()
            || u < self.u_start
            || u > self.u_end
            || v < self.v_start
            || v > self.v_end
        {
            return Err(GeomError::InvalidInput(
                "Bezier patch evaluation parameter is outside the native box".to_owned(),
            ));
        }
        let mut collapsed_v = Vec::new();
        reserve(&mut collapsed_v, self.controls.len())?;
        for row in &self.controls {
            collapsed_v.push(
                Cell {
                    controls: clone_controls(row)?,
                    start: self.v_start,
                    end: self.v_end,
                    depth: 0,
                }
                .point_at(v)?,
            );
        }
        Cell {
            controls: collapsed_v,
            start: self.u_start,
            end: self.u_end,
            depth: 0,
        }
        .point_at(u)
    }

    /// Conservative work units for restricting this complete tensor patch to an
    /// arbitrary nonempty subbox and then scanning its coordinate hull.
    ///
    /// Restricting one Bézier line may require two de Casteljau splits. A split
    /// allocates one copied level, two outputs, and every shorter intermediate
    /// level. The checked formula counts all of those homogeneous-control slots,
    /// the row/column copies used by tensor restriction, and one full hull scan.
    pub(crate) fn restriction_bound_work(&self) -> GeomResult<u128> {
        let u_count = u128::try_from(self.controls.len()).map_err(|_| work_overflow())?;
        let v_count_usize = self.controls.first().map_or(0, Vec::len);
        let v_count = u128::try_from(v_count_usize).map_err(|_| work_overflow())?;
        if u_count == 0
            || v_count == 0
            || self.controls.iter().any(|row| row.len() != v_count_usize)
        {
            return Err(GeomError::InvalidInput(
                "Bezier patch control net must be nonempty and rectangular".to_owned(),
            ));
        }

        let line_restriction_work = |count: u128| {
            count.checked_mul(count).and_then(|square| {
                count
                    .checked_mul(5)
                    .and_then(|linear| square.checked_add(linear))
            })
        };
        let u_line = u_count
            .checked_add(line_restriction_work(u_count).ok_or_else(work_overflow)?)
            .ok_or_else(work_overflow)?;
        let v_line = v_count
            .checked_add(line_restriction_work(v_count).ok_or_else(work_overflow)?)
            .ok_or_else(work_overflow)?;
        let restriction = v_count
            .checked_mul(u_line)
            .and_then(|columns| {
                u_count
                    .checked_mul(v_line)
                    .and_then(|rows| columns.checked_add(rows))
            })
            .ok_or_else(work_overflow)?;
        let hull_scan = u_count.checked_mul(v_count).ok_or_else(work_overflow)?;
        restriction.checked_add(hull_scan).ok_or_else(work_overflow)
    }

    /// Conservative work units for one interval representative and one hull scan.
    pub(crate) fn representative_bound_work(&self) -> GeomResult<u128> {
        let u_count = u128::try_from(self.controls.len()).map_err(|_| work_overflow())?;
        let v_count_usize = self.controls.first().map_or(0, Vec::len);
        let v_count = u128::try_from(v_count_usize).map_err(|_| work_overflow())?;
        if u_count == 0
            || v_count == 0
            || self.controls.iter().any(|row| row.len() != v_count_usize)
        {
            return Err(GeomError::InvalidInput(
                "Bezier patch control net must be nonempty and rectangular".to_owned(),
            ));
        }

        let split_work = |count: u128| {
            count.checked_mul(3).and_then(|linear| {
                count
                    .checked_mul(count.checked_sub(1)?)
                    .and_then(|product| product.checked_div(2))
                    .and_then(|triangle| linear.checked_add(triangle))
            })
        };
        let point_work = u_count
            .checked_mul(v_count)
            .and_then(|row_copies| {
                u_count
                    .checked_mul(split_work(v_count)?)
                    .and_then(|row_splits| row_copies.checked_add(row_splits))
            })
            .and_then(|work| work.checked_add(u_count))
            .and_then(|work| work.checked_add(split_work(u_count)?))
            .ok_or_else(work_overflow)?;
        point_work
            .checked_add(u_count.checked_mul(v_count).ok_or_else(work_overflow)?)
            .ok_or_else(work_overflow)
    }

    pub(crate) fn coordinate_intervals(&self) -> GeomResult<[Interval; 3]> {
        let mut lo = [Scalar::INFINITY; 3];
        let mut hi = [Scalar::NEG_INFINITY; 3];
        for control in self.controls.iter().flatten() {
            let point = control.euclidean()?;
            for axis in 0..3 {
                lo[axis] = lo[axis].min(point[axis].lower());
                hi[axis] = hi[axis].max(point[axis].upper());
            }
        }
        Ok([
            Interval::bounds(lo[0], hi[0])?,
            Interval::bounds(lo[1], hi[1])?,
            Interval::bounds(lo[2], hi[2])?,
        ])
    }

    pub(crate) fn partial_u_intervals(&self) -> GeomResult<[Interval; 3]> {
        self.partial_intervals(true)
    }

    pub(crate) fn partial_v_intervals(&self) -> GeomResult<[Interval; 3]> {
        self.partial_intervals(false)
    }

    fn partial_intervals(&self, along_u: bool) -> GeomResult<[Interval; 3]> {
        let u_count = self.controls.len();
        let v_count = self.controls.first().map_or(0, Vec::len);
        let degree = if along_u {
            u_count.checked_sub(1)
        } else {
            v_count.checked_sub(1)
        }
        .ok_or_else(|| GeomError::InvalidInput("empty Bezier patch control net".to_owned()))?;
        if degree == 0 {
            return Ok([Interval::exact(0.0)?; 3]);
        }
        let (start, end) = if along_u {
            (self.u_start, self.u_end)
        } else {
            (self.v_start, self.v_end)
        };
        let span = Interval::exact(end)?.subtract(Interval::exact(start)?)?;
        let scale = Interval::exact(degree as Scalar)?.divide(span)?;
        let mut derivatives = Vec::new();
        let derivative_count = if along_u {
            degree.checked_mul(v_count)
        } else {
            u_count.checked_mul(degree)
        }
        .ok_or_else(allocation_error)?;
        reserve(&mut derivatives, derivative_count)?;
        let mut push_difference = |before: &HomogeneousPoint, after: &HomogeneousPoint| {
            derivatives.push(HomogeneousPoint {
                numerator: [
                    after.numerator[0]
                        .subtract(before.numerator[0])?
                        .multiply(scale)?,
                    after.numerator[1]
                        .subtract(before.numerator[1])?
                        .multiply(scale)?,
                    after.numerator[2]
                        .subtract(before.numerator[2])?
                        .multiply(scale)?,
                ],
                weight: after.weight.subtract(before.weight)?.multiply(scale)?,
            });
            GeomResult::Ok(())
        };
        if along_u {
            for rows in self.controls.windows(2) {
                for (before, after) in rows[0].iter().zip(&rows[1]) {
                    push_difference(before, after)?;
                }
            }
        } else {
            for row in &self.controls {
                for pair in row.windows(2) {
                    push_difference(&pair[0], &pair[1])?;
                }
            }
        }
        let weights = Interval::hull(self.controls.iter().flatten().map(|point| point.weight))?;
        let weight_partial = Interval::hull(derivatives.iter().map(|point| point.weight))?;
        let denominator = weights.multiply(weights)?;
        let mut result = [Interval::exact(0.0)?; 3];
        for (axis, value) in result.iter_mut().enumerate() {
            let numerator = Interval::hull(
                self.controls
                    .iter()
                    .flatten()
                    .map(|point| point.numerator[axis]),
            )?;
            let numerator_partial =
                Interval::hull(derivatives.iter().map(|point| point.numerator[axis]))?;
            *value = numerator_partial
                .multiply(weights)?
                .subtract(numerator.multiply(weight_partial)?)?
                .divide(denominator)?;
        }
        Ok(result)
    }

    pub(crate) fn restrict(
        &self,
        u_start: Scalar,
        u_end: Scalar,
        v_start: Scalar,
        v_end: Scalar,
    ) -> GeomResult<Self> {
        use crate::certified_bezier::Cell;

        if u_start < self.u_start
            || u_end > self.u_end
            || v_start < self.v_start
            || v_end > self.v_end
            || u_start >= u_end
            || v_start >= v_end
        {
            return Err(GeomError::InvalidInput(
                "Bezier patch restriction must be a finite nonempty parameter box".to_owned(),
            ));
        }
        let v_count = self.controls[0].len();
        let mut columns = Vec::new();
        reserve(&mut columns, v_count)?;
        for v in 0..v_count {
            let mut column_controls = Vec::new();
            reserve(&mut column_controls, self.controls.len())?;
            column_controls.extend(self.controls.iter().map(|row| row[v].clone()));
            let column = Cell {
                controls: column_controls,
                start: self.u_start,
                end: self.u_end,
                depth: 0,
            }
            .restrict(u_start, u_end)?
            .controls;
            columns.push(column);
        }
        let mut controls = Vec::new();
        reserve(&mut controls, self.controls.len())?;
        for u in 0..self.controls.len() {
            let mut row_controls = Vec::new();
            reserve(&mut row_controls, columns.len())?;
            row_controls.extend(columns.iter().map(|column| column[u].clone()));
            let row = Cell {
                controls: row_controls,
                start: self.v_start,
                end: self.v_end,
                depth: 0,
            }
            .restrict(v_start, v_end)?
            .controls;
            controls.push(row);
        }
        Ok(Self {
            controls,
            u_start,
            u_end,
            v_start,
            v_end,
        })
    }
}

pub(crate) fn piecewise_bezier_patches(
    surface: &BSplineSurface,
    budget: &mut RefinementBudget,
) -> GeomResult<Vec<Patch>> {
    piecewise_bezier_patches_with_modes(surface, false, false, budget)
}

pub(crate) fn piecewise_periodic_bezier_patches(
    surface: &BSplineSurface,
    u_periodic: bool,
    v_periodic: bool,
    budget: &mut RefinementBudget,
) -> GeomResult<Vec<Patch>> {
    piecewise_bezier_patches_with_modes(surface, u_periodic, v_periodic, budget)
}

fn piecewise_bezier_patches_with_modes(
    surface: &BSplineSurface,
    u_periodic: bool,
    v_periodic: bool,
    budget: &mut RefinementBudget,
) -> GeomResult<Vec<Patch>> {
    let (rows, cols, u_plan, v_plan) = validate_surface(surface, u_periodic, v_periodic)?;
    let patch_count = u_plan
        .segments
        .checked_mul(v_plan.segments)
        .ok_or_else(allocation_error)?;
    let controls_per_patch = u_plan
        .degree
        .checked_add(1)
        .and_then(|u| v_plan.degree.checked_add(1).and_then(|v| u.checked_mul(v)))
        .ok_or_else(allocation_error)?;
    let original = rows.checked_mul(cols).ok_or_else(allocation_error)? as u128;
    let u_work = u_plan
        .work_per_line
        .checked_mul(cols as u128)
        .ok_or_else(allocation_error)?;
    let transpose = u_plan
        .refined_controls
        .checked_mul(cols)
        .ok_or_else(allocation_error)? as u128;
    let v_work = v_plan
        .work_per_line
        .checked_mul(u_plan.refined_controls as u128)
        .ok_or_else(allocation_error)?;
    let emitted = patch_count
        .checked_mul(controls_per_patch)
        .ok_or_else(allocation_error)? as u128;
    let total_work = original
        .checked_add(u_work)
        .and_then(|value| value.checked_add(transpose))
        .and_then(|value| value.checked_add(v_work))
        .and_then(|value| value.checked_add(emitted));
    budget.charge(total_work)?;

    let refined_u_knots = refined_axis_knots(u_plan, &surface.u_knots, &surface.u_multiplicities)?;
    let refined_v_knots = refined_axis_knots(v_plan, &surface.v_knots, &surface.v_multiplicities)?;
    let original_net = collect_original_net(surface, rows, cols)?;
    let mut u_columns = Vec::new();
    reserve(&mut u_columns, cols)?;
    for v in 0..cols {
        let mut column = Vec::new();
        reserve(&mut column, rows)?;
        column.extend(original_net.iter().map(|row| row[v].clone()));
        u_columns.push(refine_line(
            column,
            u_plan,
            &surface.u_knots,
            &surface.u_multiplicities,
        )?);
    }

    let mut u_rows = Vec::new();
    reserve(&mut u_rows, u_plan.refined_controls)?;
    for u in 0..u_plan.refined_controls {
        let mut row = Vec::new();
        reserve(&mut row, cols)?;
        row.extend(u_columns.iter().map(|column| column[u].clone()));
        u_rows.push(row);
    }

    let mut refined = Vec::new();
    reserve(&mut refined, u_plan.refined_controls)?;
    for row in u_rows {
        refined.push(refine_line(
            row,
            v_plan,
            &surface.v_knots,
            &surface.v_multiplicities,
        )?);
    }
    if refined
        .iter()
        .any(|row| row.len() != v_plan.refined_controls)
    {
        return Err(GeomError::Degenerate(
            "certified surface refinement produced an inconsistent control net".to_owned(),
        ));
    }

    let mut patches = Vec::new();
    reserve(&mut patches, patch_count)?;
    for u_span in surface
        .u_knots
        .windows(2)
        .filter(|pair| pair[0] >= u_plan.domain.0 && pair[1] <= u_plan.domain.1)
    {
        let u_first = span_control_first(&refined_u_knots, u_plan.degree, u_span[0])?;
        for v_span in surface
            .v_knots
            .windows(2)
            .filter(|pair| pair[0] >= v_plan.domain.0 && pair[1] <= v_plan.domain.1)
        {
            let v_first = span_control_first(&refined_v_knots, v_plan.degree, v_span[0])?;
            let u_last = u_first
                .checked_add(u_plan.degree)
                .ok_or_else(allocation_error)?;
            let v_last = v_first
                .checked_add(v_plan.degree)
                .ok_or_else(allocation_error)?;
            let source_rows = refined.get(u_first..=u_last).ok_or_else(|| {
                GeomError::Degenerate("refined surface U span escaped its control net".to_owned())
            })?;
            let mut controls = Vec::new();
            reserve(&mut controls, u_plan.degree + 1)?;
            for source_row in source_rows {
                let source_controls = source_row.get(v_first..=v_last).ok_or_else(|| {
                    GeomError::Degenerate(
                        "refined surface V span escaped its control net".to_owned(),
                    )
                })?;
                let mut row = Vec::new();
                reserve(&mut row, v_plan.degree + 1)?;
                row.extend(source_controls.iter().cloned());
                controls.push(row);
            }
            patches.push(Patch {
                controls,
                u_start: u_span[0],
                u_end: u_span[1],
                v_start: v_span[0],
                v_end: v_span[1],
            });
        }
    }
    if patches.len() != patch_count {
        return Err(GeomError::Degenerate(
            "certified surface refinement emitted an unexpected patch count".to_owned(),
        ));
    }
    Ok(patches)
}

#[cfg(test)]
mod tests {
    use axiolid_core::Point3;
    use axiolid_curve::KnotSpec;
    use axiolid_evaluate::surface::bspline_jet;

    use super::*;

    fn plane() -> BSplineSurface {
        BSplineSurface {
            u_degree: 1,
            v_degree: 1,
            control_points: vec![
                vec![Point3::new(-1.0, -1.0, 0.0), Point3::new(-1.0, 1.0, 0.0)],
                vec![Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
            ],
            u_knots: vec![2.0, 5.0],
            u_multiplicities: vec![2, 2],
            v_knots: vec![-1.0, 2.0],
            v_multiplicities: vec![2, 2],
            weights: None,
            u_closed: false,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: Some(false),
        }
    }

    #[test]
    fn bilinear_surface_refines_to_one_native_domain_patch() {
        let mut budget = RefinementBudget::new(1_000, "surface refinement test");
        let patches =
            piecewise_bezier_patches(&plane(), &mut budget).expect("bilinear surface refines");
        assert_eq!(patches.len(), 1);
        let patch = &patches[0];
        assert_eq!((patch.u_start, patch.u_end), (2.0, 5.0));
        assert_eq!((patch.v_start, patch.v_end), (-1.0, 2.0));
        assert_eq!(patch.controls.len(), 2);
        assert!(patch.controls.iter().all(|row| row.len() == 2));
        assert!(patch.controls[0][0].numerator[0].lower() <= -1.0);
        assert!(patch.controls[0][0].numerator[0].upper() >= -1.0);
        assert!(patch.controls[1][1].numerator[1].lower() <= 1.0);
        assert!(patch.controls[1][1].numerator[1].upper() >= 1.0);
    }

    #[test]
    fn bilinear_restriction_bound_accounts_for_temporaries_and_hull_scan() {
        let mut budget = RefinementBudget::new(1_000, "surface refinement test");
        let patch = piecewise_bezier_patches(&plane(), &mut budget)
            .expect("bilinear surface refines")
            .pop()
            .expect("one patch");
        assert_eq!(patch.restriction_bound_work().expect("work is bounded"), 68);
        assert_eq!(
            patch
                .representative_bound_work()
                .expect("representative work is bounded"),
            31
        );
    }

    fn periodic_u_strip() -> BSplineSurface {
        let ring = [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        BSplineSurface {
            u_degree: 1,
            v_degree: 1,
            control_points: ring
                .into_iter()
                .map(|point| vec![point, point + Point3::new(0.0, 0.0, 1.0)])
                .collect(),
            u_knots: vec![-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            u_multiplicities: vec![1; 7],
            v_knots: vec![0.0, 1.0],
            v_multiplicities: vec![2, 2],
            weights: None,
            u_closed: true,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: Some(false),
        }
    }

    #[test]
    fn periodic_active_spans_refine_to_scalar_oracle_enclosures() {
        let surface = periodic_u_strip();
        let mut budget = RefinementBudget::new(100_000, "periodic surface refinement test");
        let patches = piecewise_periodic_bezier_patches(&surface, true, false, &mut budget)
            .expect("periodic strip refines");
        assert_eq!(patches.len(), 4);
        for patch in patches {
            let u = (patch.u_start + patch.u_end) * 0.5;
            let v = (patch.v_start + patch.v_end) * 0.5;
            let enclosure = patch
                .point_at(u, v)
                .expect("patch evaluates")
                .euclidean()
                .expect("positive weight");
            let oracle = bspline_jet(&surface, u, v)
                .expect("surface evaluates")
                .point;
            for (interval, coordinate) in enclosure.into_iter().zip(oracle.to_array()) {
                assert!(interval.contains(coordinate));
            }
        }
    }

    fn periodic_u_quadratic_strip() -> BSplineSurface {
        let ring = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 3.0, 0.0),
            Point3::new(7.0, -2.0, 0.0),
            Point3::new(-4.0, 11.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 3.0, 0.0),
        ];
        BSplineSurface {
            u_degree: 2,
            v_degree: 1,
            control_points: ring
                .into_iter()
                .map(|point| vec![point, point + Point3::new(0.0, 0.0, 1.0)])
                .collect(),
            u_knots: vec![-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            u_multiplicities: vec![1; 9],
            v_knots: vec![0.0, 1.0],
            v_multiplicities: vec![2, 2],
            weights: None,
            u_closed: true,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: Some(false),
        }
    }

    #[test]
    fn quadratic_periodic_boundary_refinement_matches_scalar_oracle() {
        let surface = periodic_u_quadratic_strip();
        let mut budget = RefinementBudget::new(100_000, "periodic surface refinement test");
        let patches = piecewise_periodic_bezier_patches(&surface, true, false, &mut budget)
            .expect("quadratic periodic strip refines");
        assert_eq!(patches.len(), 4);
        for patch in patches {
            for fraction in [0.125, 0.5, 0.875] {
                let u = patch.u_start + (patch.u_end - patch.u_start) * fraction;
                let v = 0.375;
                let enclosure = patch
                    .point_at(u, v)
                    .expect("patch evaluates")
                    .euclidean()
                    .expect("positive weight");
                let oracle = bspline_jet(&surface, u, v)
                    .expect("surface evaluates")
                    .point;
                for (interval, coordinate) in enclosure.into_iter().zip(oracle.to_array()) {
                    assert!(
                        interval.contains(coordinate),
                        "span [{}, {}], fraction {fraction}, interval {interval:?}, coordinate {coordinate}",
                        patch.u_start,
                        patch.u_end,
                    );
                }
            }
        }
    }

    #[test]
    fn periodic_cyclic_extension_work_is_precharged() {
        let surface = periodic_u_strip();
        let mut insufficient = RefinementBudget::new(89, "periodic surface conversion test");
        assert!(
            piecewise_periodic_bezier_patches(&surface, true, false, &mut insufficient).is_err()
        );

        let mut exact = RefinementBudget::new(90, "periodic surface conversion test");
        let patches = piecewise_periodic_bezier_patches(&surface, true, false, &mut exact)
            .expect("exact conversion work");
        assert_eq!(patches.len(), 4);
    }

    #[test]
    fn native_span_partials_are_conservatively_scaled() {
        let mut budget = RefinementBudget::new(1_000, "surface refinement test");
        let patch = piecewise_bezier_patches(&plane(), &mut budget)
            .expect("bilinear surface refines")
            .pop()
            .expect("one patch");
        let u = patch.partial_u_intervals().expect("u partials");
        let v = patch.partial_v_intervals().expect("v partials");
        assert!(u[0].contains(2.0 / 3.0));
        assert!(u[1].contains(0.0));
        assert!(u[2].contains(0.0));
        assert!(v[0].contains(0.0));
        assert!(v[1].contains(2.0 / 3.0));
        assert!(v[2].contains(0.0));
    }

    #[test]
    fn restricted_patch_encloses_only_requested_native_subbox() {
        let mut budget = RefinementBudget::new(1_000, "surface refinement test");
        let patch = piecewise_bezier_patches(&plane(), &mut budget)
            .expect("bilinear surface refines")
            .pop()
            .expect("one patch")
            .restrict(3.0, 4.0, 0.0, 1.0)
            .expect("patch restricts");
        let bounds = patch.coordinate_intervals().expect("coordinate hull");
        assert!(bounds[0].contains(-1.0 / 3.0));
        assert!(bounds[0].contains(1.0 / 3.0));
        assert!(bounds[1].contains(-1.0 / 3.0));
        assert!(bounds[1].contains(1.0 / 3.0));
        assert!(bounds[2].contains(0.0));
    }

    #[test]
    fn rational_internal_knot_refinement_encloses_original_surface() {
        let surface = BSplineSurface {
            u_degree: 2,
            v_degree: 1,
            control_points: vec![
                vec![Point3::new(-1.0, -1.0, 0.0), Point3::new(-1.0, 1.0, 0.0)],
                vec![Point3::new(-0.3, -1.0, 0.2), Point3::new(-0.3, 1.0, 0.2)],
                vec![Point3::new(0.3, -1.0, -0.1), Point3::new(0.3, 1.0, -0.1)],
                vec![Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
            ],
            u_knots: vec![0.0, 0.5, 1.0],
            u_multiplicities: vec![3, 1, 3],
            v_knots: vec![-2.0, 2.0],
            v_multiplicities: vec![2, 2],
            weights: Some(vec![
                vec![1.0, 1.2],
                vec![0.8, 1.5],
                vec![1.4, 0.9],
                vec![1.1, 1.3],
            ]),
            u_closed: false,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: Some(false),
        };
        let mut budget = RefinementBudget::new(10_000, "surface refinement test");
        let patches = piecewise_bezier_patches(&surface, &mut budget).expect("surface refines");
        assert_eq!(patches.len(), 2);
        for patch in patches {
            let u = patch.u_start * 0.5 + patch.u_end * 0.5;
            let v = patch.v_start * 0.5 + patch.v_end * 0.5;
            let expected = bspline_jet(&surface, u, v)
                .expect("surface evaluates")
                .point;
            let enclosed = patch
                .midpoint_point()
                .expect("midpoint encloses")
                .euclidean()
                .expect("positive weight projects");
            assert!(enclosed[0].contains(expected.x));
            assert!(enclosed[1].contains(expected.y));
            assert!(enclosed[2].contains(expected.z));
        }
    }

    #[test]
    fn multispan_axis_emits_corresponding_tensor_patches() {
        let mut surface = plane();
        surface.control_points = vec![
            vec![Point3::new(-1.0, -1.0, 0.0), Point3::new(-1.0, 1.0, 0.0)],
            vec![Point3::new(0.0, -1.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        surface.u_knots = vec![2.0, 3.0, 4.0];
        surface.u_multiplicities = vec![2, 1, 2];
        let mut budget = RefinementBudget::new(2_000, "surface refinement test");
        let patches =
            piecewise_bezier_patches(&surface, &mut budget).expect("multispan surface refines");
        assert_eq!(patches.len(), 2);
        assert_eq!((patches[0].u_start, patches[0].u_end), (2.0, 3.0));
        assert_eq!((patches[1].u_start, patches[1].u_end), (3.0, 4.0));
        assert!(patches
            .iter()
            .all(|patch| (patch.v_start, patch.v_end) == (-1.0, 2.0)));
    }
}
