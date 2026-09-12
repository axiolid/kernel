//! Certified transverse intersections between internally continuous clamped 3D NURBS
//! curves and surfaces.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar};
use axiolid_curve::BSplineCurve3;
use axiolid_evaluate::{curve::bspline_jet3, surface::bspline_jet};
use axiolid_surface::BSplineSurface;

use crate::{
    certified_bezier::{Cell, Interval},
    certified_projection::ParameterInterval,
    certified_refinement::{piecewise_bezier_cells, RefinementBudget},
    certified_surface_bezier::{piecewise_bezier_patches, Patch},
};

const MAX_NODES: u32 = 100_000;
const MAX_DEPTH: u16 = 64;

/// Bounded policy for certified curve/surface root isolation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertifiedCurveSurfaceIntersectionOptions {
    parameter_tolerance: Scalar,
    max_nodes: u32,
    max_depth: u16,
}

impl CertifiedCurveSurfaceIntersectionOptions {
    /// Construct a policy with a positive finite native-parameter resolution.
    pub fn new(parameter_tolerance: Scalar, max_nodes: u32, max_depth: u16) -> GeomResult<Self> {
        if !parameter_tolerance.is_finite() || parameter_tolerance <= 0.0 {
            return Err(GeomError::InvalidInput(
                "curve/surface parameter tolerance must be finite and positive".to_owned(),
            ));
        }
        if max_nodes == 0 || max_nodes > MAX_NODES {
            return Err(GeomError::InvalidInput(format!(
                "curve/surface max_nodes must be in 1..={MAX_NODES}"
            )));
        }
        if max_depth == 0 || max_depth > MAX_DEPTH {
            return Err(GeomError::InvalidInput(format!(
                "curve/surface max_depth must be in 1..={MAX_DEPTH}"
            )));
        }
        Ok(Self {
            parameter_tolerance,
            max_nodes,
            max_depth,
        })
    }

    fn parameter_tolerance(self) -> Scalar {
        self.parameter_tolerance
    }
    fn max_nodes(self) -> u32 {
        self.max_nodes
    }
    fn max_depth(self) -> u16 {
        self.max_depth
    }
}

impl Default for CertifiedCurveSurfaceIntersectionOptions {
    fn default() -> Self {
        Self {
            parameter_tolerance: 1.0e-8,
            max_nodes: MAX_NODES,
            max_depth: MAX_DEPTH,
        }
    }
}

/// Native parameter box for a curve and a tensor-product surface patch.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct CurveSurfaceParameterBox {
    /// Curve parameter enclosure.
    pub curve: ParameterInterval,
    /// Surface `u` parameter enclosure.
    pub surface_u: ParameterInterval,
    /// Surface `v` parameter enclosure.
    pub surface_v: ParameterInterval,
}

/// Existence-and-uniqueness certificate for one transverse intersection.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct TransverseCurveSurfaceIntersection3 {
    /// Certified curve parameter enclosure.
    pub curve_parameter: ParameterInterval,
    /// Certified surface `u` parameter enclosure.
    pub surface_u_parameter: ParameterInterval,
    /// Certified surface `v` parameter enclosure.
    pub surface_v_parameter: ParameterInterval,
    /// Representative midpoint of the two evaluated images.
    pub point: Point3,
    /// Conservative residual norm upper bound over the certified parameter box.
    pub residual_upper_bound: Scalar,
    /// Positive lower bound for the absolute 3x3 Jacobian determinant.
    pub jacobian_determinant_lower_bound: Scalar,
}

/// Certified curve/surface query outcome.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum CertifiedCurveSurfaceIntersection3 {
    /// Every candidate box was excluded or proved to contain one transverse root.
    Complete {
        /// Pairwise-disjoint transverse-root certificates.
        intersections: Vec<TransverseCurveSurfaceIntersection3>,
        /// Number of parameter boxes processed or generated.
        visited_nodes: u32,
    },
    /// One or more boxes could not be proved or excluded within policy.
    Unresolved {
        /// Transverse roots that were proved before another box remained unresolved.
        intersections: Vec<TransverseCurveSurfaceIntersection3>,
        /// Conservative boxes that may contain singular, tangential, boundary, or unresolved roots.
        candidate_boxes: Vec<CurveSurfaceParameterBox>,
        /// Number of parameter boxes processed or generated.
        visited_nodes: u32,
    },
}

#[derive(Debug, Clone, Copy)]
struct Pending {
    curve_index: usize,
    patch_index: usize,
    parameters: CurveSurfaceParameterBox,
    depth: u16,
}

/// Certify all isolated transverse roots of `curve(t) = surface(u, v)`.
///
/// The implementation accepts finite polynomial inputs and rational inputs with
/// finite, strictly positive weights. All axes must be clamped and have continuous
/// internal span joins (internal multiplicity `1..=degree`). Full-multiplicity internal knots are valid NURBS representation but are not
/// supported by this certified query. A successful certificate uses outward rational Bézier bounds and a strict
/// interior 3D Krawczyk image. Tangential, singular, patch-boundary, and
/// proof-insufficient cases are returned conservatively as [`Unresolved`](CertifiedCurveSurfaceIntersection3::Unresolved).
pub fn intersect_curve_surface_certified(
    curve: &BSplineCurve3,
    surface: &BSplineSurface,
    options: CertifiedCurveSurfaceIntersectionOptions,
) -> GeomResult<CertifiedCurveSurfaceIntersection3> {
    let options = CertifiedCurveSurfaceIntersectionOptions::new(
        options.parameter_tolerance,
        options.max_nodes,
        options.max_depth,
    )?;
    let mut budget = RefinementBudget::new(
        options.max_nodes(),
        "certified curve/surface intersection budget",
    );
    let curve_cells =
        piecewise_bezier_cells(curve, |point| [point.x, point.y, point.z], &mut budget)?;
    let patches = piecewise_bezier_patches(surface, &mut budget)?;
    let initial =
        curve_cells
            .len()
            .checked_mul(patches.len())
            .ok_or(GeomError::BudgetExceeded {
                resource: "certified curve/surface intersection budget",
            })?;
    budget.charge(u128::try_from(initial).ok())?;
    let mut visited_nodes = u32::try_from(initial).map_err(|_| GeomError::BudgetExceeded {
        resource: "certified curve/surface intersection budget",
    })?;

    let mut intersections = Vec::new();
    let mut unresolved = Vec::new();
    // The surface's own domain, captured before subdivision so a true
    // edge stays distinguishable from a face introduced by splitting.
    let domain = SurfaceDomain {
        u_start: *surface
            .u_knots
            .first()
            .ok_or_else(|| GeomError::InvalidInput("surface has no u knots".to_owned()))?,
        u_end: *surface
            .u_knots
            .last()
            .ok_or_else(|| GeomError::InvalidInput("surface has no u knots".to_owned()))?,
        v_start: *surface
            .v_knots
            .first()
            .ok_or_else(|| GeomError::InvalidInput("surface has no v knots".to_owned()))?,
        v_end: *surface
            .v_knots
            .last()
            .ok_or_else(|| GeomError::InvalidInput("surface has no v knots".to_owned()))?,
    };
    let mut pending = Vec::new();
    pending
        .try_reserve(usize::from(options.max_depth()) + 1)
        .map_err(|_| GeomError::BudgetExceeded {
            resource: "certified curve/surface pending allocation",
        })?;

    for curve_index in 0..curve_cells.len() {
        for patch_index in 0..patches.len() {
            pending.push(Pending {
                curve_index,
                patch_index,
                parameters: base_box(&curve_cells[curve_index], &patches[patch_index]),
                depth: 0,
            });
            while let Some(current) = pending.pop() {
                let curve_cell = curve_cells[current.curve_index]
                    .restrict(current.parameters.curve.start, current.parameters.curve.end)?;
                let patch = patches[current.patch_index].restrict(
                    current.parameters.surface_u.start,
                    current.parameters.surface_u.end,
                    current.parameters.surface_v.start,
                    current.parameters.surface_v.end,
                )?;
                if residual_excludes_zero(&curve_cell, &patch)? {
                    continue;
                }
                if let Some(root) = krawczyk_root(curve, surface, &curve_cell, &patch)? {
                    if certificate_meets_resolution(&root, options.parameter_tolerance()) {
                        push_result(&mut intersections, root)?;
                        continue;
                    }
                    if current.depth >= options.max_depth() {
                        push_result(&mut unresolved, certificate_box(&root))?;
                        continue;
                    }
                    let contracted = CurveSurfaceParameterBox {
                        curve: contract_interval(
                            current.parameters.curve,
                            root.curve_parameter,
                            options.parameter_tolerance(),
                        ),
                        surface_u: contract_interval(
                            current.parameters.surface_u,
                            root.surface_u_parameter,
                            options.parameter_tolerance(),
                        ),
                        surface_v: contract_interval(
                            current.parameters.surface_v,
                            root.surface_v_parameter,
                            options.parameter_tolerance(),
                        ),
                    };
                    if contracted == current.parameters {
                        push_result(&mut unresolved, contracted)?;
                        continue;
                    }
                    budget.charge(Some(1))?;
                    visited_nodes = checked_nodes(visited_nodes, 1, options.max_nodes())?;
                    pending.push(Pending {
                        parameters: contracted,
                        depth: current.depth.checked_add(1).ok_or_else(|| {
                            GeomError::Degenerate(
                                "curve/surface intersection depth overflow".to_owned(),
                            )
                        })?,
                        ..current
                    });
                    continue;
                }
                // The 3x3 operator could not certify this box. Before
                // subdividing, try the reduced solve on any DOMAIN edge the
                // box touches: a root sitting exactly on a patch edge is
                // unreachable for strict containment and would otherwise
                // halve forever until depth runs out.
                //
                // Only true domain edges qualify. A face created by
                // subdivision has interior geometry on both sides, so a root
                // there is genuinely interior and belongs to the 3x3 path;
                // pinning it would report a root the caller could also find
                // by subdividing, i.e. a duplicate.
                if let Some(root) = edge_root(
                    curve,
                    surface,
                    &curve_cell,
                    &patch,
                    &domain,
                    options.parameter_tolerance(),
                )? {
                    push_result(&mut intersections, root)?;
                    continue;
                }
                if current.depth >= options.max_depth() {
                    push_result(&mut unresolved, current.parameters)?;
                    continue;
                }
                budget.charge(Some(2))?;
                visited_nodes = checked_nodes(visited_nodes, 2, options.max_nodes())?;
                split_pending(current, &mut pending)?;
            }
        }
    }

    if unresolved.is_empty() {
        Ok(CertifiedCurveSurfaceIntersection3::Complete {
            intersections,
            visited_nodes,
        })
    } else {
        Ok(CertifiedCurveSurfaceIntersection3::Unresolved {
            intersections,
            candidate_boxes: unresolved,
            visited_nodes,
        })
    }
}

fn checked_nodes(current: u32, additional: u32, maximum: u32) -> GeomResult<u32> {
    current
        .checked_add(additional)
        .filter(|&value| value <= maximum)
        .ok_or(GeomError::BudgetExceeded {
            resource: "certified curve/surface intersection budget",
        })
}

fn push_result<T>(target: &mut Vec<T>, value: T) -> GeomResult<()> {
    target
        .try_reserve(1)
        .map_err(|_| GeomError::BudgetExceeded {
            resource: "certified curve/surface result allocation",
        })?;
    target.push(value);
    Ok(())
}

fn base_box(curve: &Cell, patch: &Patch) -> CurveSurfaceParameterBox {
    CurveSurfaceParameterBox {
        curve: ParameterInterval {
            start: curve.start,
            end: curve.end,
        },
        surface_u: ParameterInterval {
            start: patch.u_start,
            end: patch.u_end,
        },
        surface_v: ParameterInterval {
            start: patch.v_start,
            end: patch.v_end,
        },
    }
}

fn certificate_box(root: &TransverseCurveSurfaceIntersection3) -> CurveSurfaceParameterBox {
    CurveSurfaceParameterBox {
        curve: root.curve_parameter,
        surface_u: root.surface_u_parameter,
        surface_v: root.surface_v_parameter,
    }
}

fn certificate_meets_resolution(
    root: &TransverseCurveSurfaceIntersection3,
    tolerance: Scalar,
) -> bool {
    root.curve_parameter.end - root.curve_parameter.start <= tolerance
        && root.surface_u_parameter.end - root.surface_u_parameter.start <= tolerance
        && root.surface_v_parameter.end - root.surface_v_parameter.start <= tolerance
}

fn residual_excludes_zero(curve: &Cell, patch: &Patch) -> GeomResult<bool> {
    let curve = curve.coordinate_intervals()?;
    let surface = patch.coordinate_intervals()?;
    Ok((0..3).any(|axis| {
        curve[axis].upper() < surface[axis].lower() || surface[axis].upper() < curve[axis].lower()
    }))
}

fn split_interval(
    interval: ParameterInterval,
) -> GeomResult<(ParameterInterval, ParameterInterval)> {
    let midpoint = interval.start * 0.5 + interval.end * 0.5;
    if midpoint <= interval.start || midpoint >= interval.end {
        return Err(GeomError::Degenerate(
            "certified curve/surface parameter split did not advance".to_owned(),
        ));
    }
    Ok((
        ParameterInterval {
            start: interval.start,
            end: midpoint,
        },
        ParameterInterval {
            start: midpoint,
            end: interval.end,
        },
    ))
}

fn split_pending(current: Pending, pending: &mut Vec<Pending>) -> GeomResult<()> {
    let widths = [
        current.parameters.curve.end - current.parameters.curve.start,
        current.parameters.surface_u.end - current.parameters.surface_u.start,
        current.parameters.surface_v.end - current.parameters.surface_v.start,
    ];
    let axis = if widths[0] >= widths[1] && widths[0] >= widths[2] {
        0
    } else if widths[1] >= widths[2] {
        1
    } else {
        2
    };
    let source = match axis {
        0 => current.parameters.curve,
        1 => current.parameters.surface_u,
        _ => current.parameters.surface_v,
    };
    let (left, right) = split_interval(source)?;
    let depth = current.depth.checked_add(1).ok_or_else(|| {
        GeomError::Degenerate("curve/surface intersection depth overflow".to_owned())
    })?;
    let with_interval = |interval| {
        let mut parameters = current.parameters;
        match axis {
            0 => parameters.curve = interval,
            1 => parameters.surface_u = interval,
            _ => parameters.surface_v = interval,
        }
        Pending {
            parameters,
            depth,
            ..current
        }
    };
    pending.push(with_interval(left));
    pending.push(with_interval(right));
    Ok(())
}

fn stable_start(cell: Scalar, root: Scalar, desired: Scalar) -> Scalar {
    if desired > cell {
        desired.min(root)
    } else {
        root
    }
}

fn stable_end(cell: Scalar, root: Scalar, desired: Scalar) -> Scalar {
    if desired < cell {
        desired.max(root)
    } else {
        root
    }
}

fn contract_interval(
    cell: ParameterInterval,
    root: ParameterInterval,
    tolerance: Scalar,
) -> ParameterInterval {
    if cell.end - cell.start <= tolerance {
        return cell;
    }
    let center = root.start * 0.5 + root.end * 0.5;
    let half = tolerance * 0.5;
    ParameterInterval {
        start: stable_start(cell.start, root.start, center - half),
        end: stable_end(cell.end, root.end, center + half),
    }
}

fn krawczyk_root(
    curve: &BSplineCurve3,
    surface: &BSplineSurface,
    curve_cell: &Cell,
    patch: &Patch,
) -> GeomResult<Option<TransverseCurveSurfaceIntersection3>> {
    let center = [
        curve_cell.start * 0.5 + curve_cell.end * 0.5,
        patch.u_start * 0.5 + patch.u_end * 0.5,
        patch.v_start * 0.5 + patch.v_end * 0.5,
    ];
    let curve_jet = bspline_jet3(curve, center[0])?;
    let surface_jet = bspline_jet(surface, center[1], center[2])?;
    let point_jacobian = [
        [curve_jet.first.x, -surface_jet.du.x, -surface_jet.dv.x],
        [curve_jet.first.y, -surface_jet.du.y, -surface_jet.dv.y],
        [curve_jet.first.z, -surface_jet.du.z, -surface_jet.dv.z],
    ];
    let Some(inverse) = inverse3(point_jacobian) else {
        return Ok(None);
    };

    let curve_midpoint = curve_cell.midpoint_point()?.euclidean()?;
    let surface_midpoint = patch.midpoint_point()?.euclidean()?;
    let residual = [
        curve_midpoint[0].subtract(surface_midpoint[0])?,
        curve_midpoint[1].subtract(surface_midpoint[1])?,
        curve_midpoint[2].subtract(surface_midpoint[2])?,
    ];
    let curve_derivative = curve_cell.derivative_intervals()?;
    let surface_u = patch.partial_u_intervals()?;
    let surface_v = patch.partial_v_intervals()?;
    let minus_one = Interval::exact(-1.0)?;
    let jacobian = [
        [
            curve_derivative[0],
            surface_u[0].multiply(minus_one)?,
            surface_v[0].multiply(minus_one)?,
        ],
        [
            curve_derivative[1],
            surface_u[1].multiply(minus_one)?,
            surface_v[1].multiply(minus_one)?,
        ],
        [
            curve_derivative[2],
            surface_u[2].multiply(minus_one)?,
            surface_v[2].multiply(minus_one)?,
        ],
    ];

    let zero = Interval::exact(0.0)?;
    let one = Interval::exact(1.0)?;
    let mut corrected = [zero; 3];
    for row in 0..3 {
        corrected[row] = Interval::exact(center[row])?.subtract(dot3(inverse[row], residual)?)?;
    }
    let mut matrix = [[zero; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            let jacobian_column = [
                jacobian[0][column],
                jacobian[1][column],
                jacobian[2][column],
            ];
            let identity = if row == column { one } else { zero };
            matrix[row][column] = identity.subtract(dot3(inverse[row], jacobian_column)?)?;
        }
    }
    let bounds = [
        ParameterInterval {
            start: curve_cell.start,
            end: curve_cell.end,
        },
        ParameterInterval {
            start: patch.u_start,
            end: patch.u_end,
        },
        ParameterInterval {
            start: patch.v_start,
            end: patch.v_end,
        },
    ];
    let mut delta = [zero; 3];
    for axis in 0..3 {
        delta[axis] = Interval::hull([
            Interval::exact(bounds[axis].start)?.subtract(Interval::exact(center[axis])?)?,
            Interval::exact(bounds[axis].end)?.subtract(Interval::exact(center[axis])?)?,
        ])?;
    }
    let mut image = [zero; 3];
    for row in 0..3 {
        image[row] = corrected[row].add(
            matrix[row][0]
                .multiply(delta[0])?
                .add(matrix[row][1].multiply(delta[1])?)?
                .add(matrix[row][2].multiply(delta[2])?)?,
        )?;
    }
    if !(image[0].lower() > bounds[0].start
        && image[0].upper() < bounds[0].end
        && image[1].lower() > bounds[1].start
        && image[1].upper() < bounds[1].end
        && image[2].lower() > bounds[2].start
        && image[2].upper() < bounds[2].end)
    {
        return Ok(None);
    }
    let determinant_lower = determinant3_interval(jacobian)?.absolute_lower_bound();
    if determinant_lower == 0.0 {
        return Ok(None);
    }
    let curve_parameter = ParameterInterval {
        start: image[0].lower(),
        end: image[0].upper(),
    };
    let surface_u_parameter = ParameterInterval {
        start: image[1].lower(),
        end: image[1].upper(),
    };
    let surface_v_parameter = ParameterInterval {
        start: image[2].lower(),
        end: image[2].upper(),
    };
    let root_curve = curve_cell.restrict(curve_parameter.start, curve_parameter.end)?;
    let root_patch = patch.restrict(
        surface_u_parameter.start,
        surface_u_parameter.end,
        surface_v_parameter.start,
        surface_v_parameter.end,
    )?;
    let residual_upper_bound = residual_norm_upper(&root_curve, &root_patch)?;
    let curve_value = bspline_jet3(curve, interval_midpoint(curve_parameter))?.point;
    let surface_value = bspline_jet(
        surface,
        interval_midpoint(surface_u_parameter),
        interval_midpoint(surface_v_parameter),
    )?
    .point;
    let point = Point3::new(
        curve_value.x * 0.5 + surface_value.x * 0.5,
        curve_value.y * 0.5 + surface_value.y * 0.5,
        curve_value.z * 0.5 + surface_value.z * 0.5,
    );
    if !point.is_finite() || !residual_upper_bound.is_finite() {
        return Err(GeomError::Degenerate(
            "certified curve/surface representative overflowed".to_owned(),
        ));
    }
    Ok(Some(TransverseCurveSurfaceIntersection3 {
        curve_parameter,
        surface_u_parameter,
        surface_v_parameter,
        point,
        residual_upper_bound,
        jacobian_determinant_lower_bound: determinant_lower,
    }))
}

fn dot3(coefficients: [Scalar; 3], values: [Interval; 3]) -> GeomResult<Interval> {
    Interval::exact(coefficients[0])?
        .multiply(values[0])?
        .add(Interval::exact(coefficients[1])?.multiply(values[1])?)?
        .add(Interval::exact(coefficients[2])?.multiply(values[2])?)
}

fn inverse3(matrix: [[Scalar; 3]; 3]) -> Option<[[Scalar; 3]; 3]> {
    if matrix.iter().flatten().any(|value| !value.is_finite()) {
        return None;
    }
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    let determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if determinant == 0.0 || !determinant.is_finite() {
        return None;
    }
    let inverse = [
        [
            (e * i - f * h) / determinant,
            (c * h - b * i) / determinant,
            (b * f - c * e) / determinant,
        ],
        [
            (f * g - d * i) / determinant,
            (a * i - c * g) / determinant,
            (c * d - a * f) / determinant,
        ],
        [
            (d * h - e * g) / determinant,
            (b * g - a * h) / determinant,
            (a * e - b * d) / determinant,
        ],
    ];
    inverse
        .iter()
        .flatten()
        .all(|value| value.is_finite())
        .then_some(inverse)
}

fn determinant3_interval(matrix: [[Interval; 3]; 3]) -> GeomResult<Interval> {
    let first = matrix[0][0].multiply(
        matrix[1][1]
            .multiply(matrix[2][2])?
            .subtract(matrix[1][2].multiply(matrix[2][1])?)?,
    )?;
    let second = matrix[0][1].multiply(
        matrix[1][0]
            .multiply(matrix[2][2])?
            .subtract(matrix[1][2].multiply(matrix[2][0])?)?,
    )?;
    let third = matrix[0][2].multiply(
        matrix[1][0]
            .multiply(matrix[2][1])?
            .subtract(matrix[1][1].multiply(matrix[2][0])?)?,
    )?;
    first.subtract(second)?.add(third)
}

fn residual_norm_upper(curve: &Cell, patch: &Patch) -> GeomResult<Scalar> {
    let curve = curve.coordinate_intervals()?;
    let surface = patch.coordinate_intervals()?;
    let mut maximum = [0.0; 3];
    for axis in 0..3 {
        let difference = curve[axis].subtract(surface[axis])?;
        maximum[axis] = difference.lower().abs().max(difference.upper().abs());
    }
    let norm = maximum[0].hypot(maximum[1]).hypot(maximum[2]);
    if !norm.is_finite() {
        return Err(GeomError::Degenerate(
            "certified curve/surface residual bound overflowed".to_owned(),
        ));
    }
    Ok(next_up(norm))
}

fn interval_midpoint(interval: ParameterInterval) -> Scalar {
    interval.start * 0.5 + interval.end * 0.5
}

fn next_up(value: Scalar) -> Scalar {
    if value == Scalar::INFINITY {
        return value;
    }
    if value == 0.0 {
        return Scalar::from_bits(1);
    }
    let bits = value.to_bits();
    if value > 0.0 {
        Scalar::from_bits(bits + 1)
    } else {
        Scalar::from_bits(bits - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_unbounded_or_zero_policy() {
        assert!(CertifiedCurveSurfaceIntersectionOptions::new(0.0, 1, 1).is_err());
        assert!(CertifiedCurveSurfaceIntersectionOptions::new(1.0e-6, 0, 1).is_err());
        assert!(CertifiedCurveSurfaceIntersectionOptions::new(1.0e-6, MAX_NODES + 1, 1).is_err());
        assert!(CertifiedCurveSurfaceIntersectionOptions::new(1.0e-6, 1, MAX_DEPTH + 1).is_err());
    }

    #[test]
    fn interval_determinant_encloses_identity() {
        let zero = Interval::exact(0.0).expect("zero");
        let one = Interval::exact(1.0).expect("one");
        let determinant =
            determinant3_interval([[one, zero, zero], [zero, one, zero], [zero, zero, one]])
                .expect("determinant");
        assert!(determinant.lower() <= 1.0 && determinant.upper() >= 1.0);
        assert!(determinant.absolute_lower_bound() > 0.0);
    }
}

/// Which surface parameter is pinned to a domain edge, and where.
///
/// A root lying exactly ON a patch edge cannot be certified by the 3x3
/// operator: strict containment `image.lower() > bounds.start` is
/// unsatisfiable when the root equals `bounds.start`. Pinning that one
/// parameter turns the 3-unknown system into a 2-unknown one whose root
/// IS interior, so the same Krawczyk argument applies unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PinnedEdge {
    /// Surface `u` is fixed at the domain edge.
    U,
    /// Surface `v` is fixed at the domain edge.
    V,
}
/// Certify a root pinned to `edge` at parameter `pinned_value`.
///
/// Solves the reduced 2x2 system in the two free unknowns with the same
/// interval-Krawczyk argument the 3x3 path uses: build the operator at the
/// box centre, require its image to land strictly inside the free-parameter
/// box, and require a nonzero Jacobian determinant lower bound. The pinned
/// parameter is reported as a degenerate interval, which is the honest
/// enclosure: it is known exactly, not bracketed.
fn krawczyk_root_on_edge(
    curve: &BSplineCurve3,
    surface: &BSplineSurface,
    curve_cell: &Cell,
    patch: &Patch,
    edge: PinnedEdge,
    pinned_value: Scalar,
) -> GeomResult<Option<TransverseCurveSurfaceIntersection3>> {
    // Free axes: curve t is always free; the other surface parameter is
    // whichever one is not pinned.
    let (free_start, free_end) = match edge {
        PinnedEdge::U => (patch.v_start, patch.v_end),
        PinnedEdge::V => (patch.u_start, patch.u_end),
    };
    let center_t = curve_cell.start * 0.5 + curve_cell.end * 0.5;
    let center_free = free_start * 0.5 + free_end * 0.5;
    let (center_u, center_v) = match edge {
        PinnedEdge::U => (pinned_value, center_free),
        PinnedEdge::V => (center_free, pinned_value),
    };

    let curve_jet = bspline_jet3(curve, center_t)?;
    let surface_jet = bspline_jet(surface, center_u, center_v)?;
    // Column 2 is the derivative along the FREE surface parameter only.
    let free_partial = match edge {
        PinnedEdge::U => surface_jet.dv,
        PinnedEdge::V => surface_jet.du,
    };
    let point_jacobian = [
        [curve_jet.first.x, -free_partial.x],
        [curve_jet.first.y, -free_partial.y],
        [curve_jet.first.z, -free_partial.z],
    ];

    // Three equations, two unknowns. Pick the two rows whose 2x2 minor has
    // the largest magnitude: that is the best-conditioned choice available,
    // and any nonzero minor is a valid certificate basis. The discarded row
    // is not ignored -- the residual bound below is taken over all three.
    let rows = [(0usize, 1usize), (0, 2), (1, 2)];
    let mut best: Option<RowSelection> = None;
    let mut best_magnitude = 0.0;
    for (first, second) in rows {
        let candidate = [
            [point_jacobian[first][0], point_jacobian[first][1]],
            [point_jacobian[second][0], point_jacobian[second][1]],
        ];
        let Some(inverse) = inverse2(candidate) else {
            continue;
        };
        let magnitude =
            (candidate[0][0] * candidate[1][1] - candidate[0][1] * candidate[1][0]).abs();
        if magnitude > best_magnitude {
            best_magnitude = magnitude;
            best = Some((candidate, inverse, (first, second)));
        }
    }
    let Some((_, inverse, (row_a, row_b))) = best else {
        return Ok(None);
    };

    let curve_midpoint = curve_cell.midpoint_point()?.euclidean()?;
    let surface_midpoint = patch.midpoint_point()?.euclidean()?;
    let residual_all = [
        curve_midpoint[0].subtract(surface_midpoint[0])?,
        curve_midpoint[1].subtract(surface_midpoint[1])?,
        curve_midpoint[2].subtract(surface_midpoint[2])?,
    ];
    let residual = [residual_all[row_a], residual_all[row_b]];

    let curve_derivative = curve_cell.derivative_intervals()?;
    let free_intervals = match edge {
        PinnedEdge::U => patch.partial_v_intervals()?,
        PinnedEdge::V => patch.partial_u_intervals()?,
    };
    let minus_one = Interval::exact(-1.0)?;
    let jacobian = [
        [
            curve_derivative[row_a],
            free_intervals[row_a].multiply(minus_one)?,
        ],
        [
            curve_derivative[row_b],
            free_intervals[row_b].multiply(minus_one)?,
        ],
    ];

    let zero = Interval::exact(0.0)?;
    let one = Interval::exact(1.0)?;
    let center = [center_t, center_free];
    let mut corrected = [zero; 2];
    for row in 0..2 {
        corrected[row] = Interval::exact(center[row])?.subtract(dot2(inverse[row], residual)?)?;
    }
    let mut matrix = [[zero; 2]; 2];
    for row in 0..2 {
        for column in 0..2 {
            let jacobian_column = [jacobian[0][column], jacobian[1][column]];
            let identity = if row == column { one } else { zero };
            matrix[row][column] = identity.subtract(dot2(inverse[row], jacobian_column)?)?;
        }
    }
    let bounds = [
        ParameterInterval {
            start: curve_cell.start,
            end: curve_cell.end,
        },
        ParameterInterval {
            start: free_start,
            end: free_end,
        },
    ];
    let mut delta = [zero; 2];
    for axis in 0..2 {
        delta[axis] = Interval::hull([
            Interval::exact(bounds[axis].start)?.subtract(Interval::exact(center[axis])?)?,
            Interval::exact(bounds[axis].end)?.subtract(Interval::exact(center[axis])?)?,
        ])?;
    }
    let mut image = [zero; 2];
    for row in 0..2 {
        image[row] = corrected[row].add(
            matrix[row][0]
                .multiply(delta[0])?
                .add(matrix[row][1].multiply(delta[1])?)?,
        )?;
    }
    // Same strict-containment demand as the 3x3 path, now over the two free
    // axes only. The pinned axis needs no containment: it is fixed exactly.
    if !(image[0].lower() > bounds[0].start
        && image[0].upper() < bounds[0].end
        && image[1].lower() > bounds[1].start
        && image[1].upper() < bounds[1].end)
    {
        return Ok(None);
    }
    let determinant_lower = determinant2_interval(jacobian)?.absolute_lower_bound();
    if determinant_lower == 0.0 {
        return Ok(None);
    }

    let curve_parameter = ParameterInterval {
        start: image[0].lower(),
        end: image[0].upper(),
    };
    let free_parameter = ParameterInterval {
        start: image[1].lower(),
        end: image[1].upper(),
    };
    let pinned = ParameterInterval {
        start: pinned_value,
        end: pinned_value,
    };
    let (surface_u_parameter, surface_v_parameter) = match edge {
        PinnedEdge::U => (pinned, free_parameter),
        PinnedEdge::V => (free_parameter, pinned),
    };

    // `restrict` rejects a zero-width box, and the pinned axis is exactly
    // zero-width by construction. Widen it to the smallest sliver that still
    // lies inside the patch. This is sound in both directions a caller
    // depends on: a SUPERSET of the true degenerate box can only inflate
    // `residual_norm_upper` (conservative) and can only make the
    // disjointness test below harder to pass (conservative). It never
    // certifies a root that the exact degenerate box would reject.
    let root_curve = {
        let (start, end) = sliver(
            curve_parameter.start,
            curve_parameter.end,
            curve_cell.start,
            curve_cell.end,
        );
        curve_cell.restrict(start, end)?
    };
    let root_patch = {
        let (u_start, u_end) = sliver(
            surface_u_parameter.start,
            surface_u_parameter.end,
            patch.u_start,
            patch.u_end,
        );
        let (v_start, v_end) = sliver(
            surface_v_parameter.start,
            surface_v_parameter.end,
            patch.v_start,
            patch.v_end,
        );
        patch.restrict(u_start, u_end, v_start, v_end)?
    };
    // The 2x2 solve only enforced two of the three coordinate equations.
    // The third must be CHECKED, not assumed: a curve can meet the edge line
    // in the chosen two coordinates while missing it in the third. The
    // residual bound is conservative over the whole certified box, so
    // requiring it to enclose zero is a sound test.
    //
    // NOT COVERED BY TESTS: every off-patch case reachable today is rejected
    // earlier by `residual_excludes_zero` during subdivision, so this branch
    // never fires in the current corpus -- deleting it does not fail any
    // test. It is retained as defence for inputs that reach here with a
    // satisfied 2x2 system and a violated third row, which the planar-only
    // certified path cannot yet construct. Do not treat it as verified.
    let residual_upper_bound = residual_norm_upper(&root_curve, &root_patch)?;
    let discarded = 3 - row_a - row_b;
    let curve_span = root_curve.coordinate_intervals()?[discarded];
    let patch_span = root_patch.coordinate_intervals()?[discarded];
    if curve_span.upper() < patch_span.lower() || patch_span.upper() < curve_span.lower() {
        return Ok(None);
    }

    let curve_value = bspline_jet3(curve, interval_midpoint(curve_parameter))?.point;
    let surface_value = bspline_jet(
        surface,
        interval_midpoint(surface_u_parameter),
        interval_midpoint(surface_v_parameter),
    )?
    .point;
    let point = Point3::new(
        curve_value.x * 0.5 + surface_value.x * 0.5,
        curve_value.y * 0.5 + surface_value.y * 0.5,
        curve_value.z * 0.5 + surface_value.z * 0.5,
    );
    if !point.is_finite() || !residual_upper_bound.is_finite() {
        return Err(GeomError::Degenerate(
            "certified curve/surface edge representative overflowed".to_owned(),
        ));
    }
    Ok(Some(TransverseCurveSurfaceIntersection3 {
        curve_parameter,
        surface_u_parameter,
        surface_v_parameter,
        point,
        residual_upper_bound,
        jacobian_determinant_lower_bound: determinant_lower,
    }))
}

/// Widen a possibly-degenerate parameter range into a nonempty one.
///
/// The pinned axis of a boundary-restricted certificate has `start == end`,
/// which `restrict` rejects. Returns the smallest nonempty range containing
/// `[start, end]` and clamped inside `[low, high]`, so the result is always
/// a superset of the requested range that the Bezier restriction accepts.
fn sliver(start: Scalar, end: Scalar, low: Scalar, high: Scalar) -> (Scalar, Scalar) {
    if end > start {
        return (start, end);
    }
    if high <= low {
        return (low, high);
    }
    let widened_end = next_up(end);
    if widened_end < high {
        return (start, widened_end);
    }
    // At the upper limit there is no room above, so grow downwards instead.
    let widened_start = -next_up(-start);
    if widened_start > low {
        (widened_start, end)
    } else {
        (low, high)
    }
}

/// A certified 2x2 row selection for a boundary-restricted solve: the chosen
/// submatrix, its numeric inverse, and which two coordinate rows were kept.
type RowSelection = ([[Scalar; 2]; 2], [[Scalar; 2]; 2], (usize, usize));

fn dot2(coefficients: [Scalar; 2], values: [Interval; 2]) -> GeomResult<Interval> {
    Interval::exact(coefficients[0])?
        .multiply(values[0])?
        .add(Interval::exact(coefficients[1])?.multiply(values[1])?)
}

fn inverse2(matrix: [[Scalar; 2]; 2]) -> Option<[[Scalar; 2]; 2]> {
    if matrix.iter().flatten().any(|value| !value.is_finite()) {
        return None;
    }
    let [[a, b], [c, d]] = matrix;
    let determinant = a * d - b * c;
    if determinant == 0.0 || !determinant.is_finite() {
        return None;
    }
    let inverse = [
        [d / determinant, -b / determinant],
        [-c / determinant, a / determinant],
    ];
    inverse
        .iter()
        .flatten()
        .all(|value| value.is_finite())
        .then_some(inverse)
}

fn determinant2_interval(matrix: [[Interval; 2]; 2]) -> GeomResult<Interval> {
    matrix[0][0]
        .multiply(matrix[1][1])?
        .subtract(matrix[0][1].multiply(matrix[1][0])?)
}

/// The surface's own parameter domain, used to tell a true edge from a
/// face created by subdivision.
#[derive(Debug, Clone, Copy)]
struct SurfaceDomain {
    u_start: Scalar,
    u_end: Scalar,
    v_start: Scalar,
    v_end: Scalar,
}

/// Try a boundary-restricted certificate on every domain edge this box
/// touches.
///
/// Returns the first certified root. At most four edges are examined and
/// each is a bounded 2x2 solve, so this adds constant work per node.
fn edge_root(
    curve: &BSplineCurve3,
    surface: &BSplineSurface,
    curve_cell: &Cell,
    patch: &Patch,
    domain: &SurfaceDomain,
    tolerance: Scalar,
) -> GeomResult<Option<TransverseCurveSurfaceIntersection3>> {
    let candidates = [
        (
            patch.u_start == domain.u_start,
            PinnedEdge::U,
            patch.u_start,
        ),
        (patch.u_end == domain.u_end, PinnedEdge::U, patch.u_end),
        (
            patch.v_start == domain.v_start,
            PinnedEdge::V,
            patch.v_start,
        ),
        (patch.v_end == domain.v_end, PinnedEdge::V, patch.v_end),
    ];
    for (on_domain_edge, edge, value) in candidates {
        if !on_domain_edge {
            continue;
        }
        if let Some(root) = krawczyk_root_on_edge(curve, surface, curve_cell, patch, edge, value)? {
            if certificate_meets_resolution(&root, tolerance) {
                return Ok(Some(root));
            }
        }
    }
    Ok(None)
}
