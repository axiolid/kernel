//! Curve interpolation and lofting (#33).
//!
//! # Interpolation, not approximation
//!
//! `interpolate_curve3` produces a curve that passes THROUGH its input points,
//! not near them. That distinction is the contract: a caller handing over
//! survey points or a section outline needs those points on the curve, and an
//! approximating fit that misses them by a tolerance is a different operation
//! with different uses.
//!
//! The test that matters therefore evaluates the result at each computed
//! parameter and requires the input point back to near machine precision.
//!
//! # Continuity is documented, not assumed
//!
//! A cubic interpolant through n points is C2 across interior joins: the
//! natural consequence of a single global system with continuity built into
//! the basis. A piecewise fit stitched segment by segment would only be C0,
//! which is why this solves globally rather than locally.
//!
//! # Lofting
//!
//! `loft_surface` runs a surface through an ordered set of section curves. The
//! sections become rows of the control net, so the surface interpolates each
//! section exactly. Sections must agree in degree and control-point count:
//! reconciling mismatched sections means knot merging and degree elevation,
//! and doing it implicitly would hide a shape change inside a construction.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar};
use axiolid_curve::{BSplineCurve3, KnotSpec};
use axiolid_surface::BSplineSurface;

/// Interpolate a cubic B-spline through `points` in order.
///
/// The curve passes through every input point. Parameters are assigned by
/// chord length, which keeps the parameterisation proportional to distance --
/// uniform spacing on unevenly spaced points produces visible overshoot.
///
/// Fewer than two points cannot define a curve, and repeated coincident points
/// make chord length degenerate, so both are refused.
pub fn interpolate_curve3(points: &[Point3]) -> GeomResult<BSplineCurve3> {
    if points.len() < 2 {
        return Err(GeomError::InvalidInput(
            "interpolation needs at least two points".to_owned(),
        ));
    }
    if !points.iter().all(|p| p.is_finite()) {
        return Err(GeomError::InvalidInput(
            "interpolation points must be finite".to_owned(),
        ));
    }
    let parameters = chord_parameters(points)?;
    interpolate_with(points, &parameters)
}

/// Chord-length parameters normalised to `[0, 1]`.
pub(crate) fn chord_parameters(points: &[Point3]) -> GeomResult<Vec<Scalar>> {
    let mut distances = Vec::with_capacity(points.len());
    distances.push(0.0);
    let mut total = 0.0;
    for pair in points.windows(2) {
        let step = (pair[1] - pair[0]).length();
        if step <= 0.0 {
            return Err(GeomError::Degenerate(
                "consecutive interpolation points coincide, so chord length is undefined"
                    .to_owned(),
            ));
        }
        total += step;
        distances.push(total);
    }
    Ok(distances.into_iter().map(|d| d / total).collect())
}

/// Cubic degree, or lower when there are too few points to support it.
///
/// Three points cannot define a cubic, so the degree drops rather than the
/// call failing: an interpolant through two or three points is still a useful
/// answer, and padding with invented points would fabricate shape.
fn degree_for(count: usize) -> u16 {
    match count {
        0 | 1 => 0,
        2 => 1,
        3 => 2,
        _ => 3,
    }
}

/// Expanded (repeated) knot vector for interpolation, averaging interior knots.
///
/// Averaging is what keeps the system well-conditioned: it guarantees every
/// basis function has a parameter inside its support, so the matrix is
/// non-singular.
fn averaged_knots(parameters: &[Scalar], degree: usize) -> Vec<Scalar> {
    let n = parameters.len() - 1;
    let mut knots = vec![0.0; degree + 1];
    for j in 1..=n.saturating_sub(degree) {
        let sum: Scalar = parameters[j..j + degree].iter().sum();
        knots.push(sum / degree as Scalar);
    }
    knots.extend(core::iter::repeat_n(1.0, degree + 1));
    knots
}

/// Knot span containing `t`, clamped to the last non-empty span.
pub(crate) fn span_of(knots: &[Scalar], n: usize, degree: usize, t: Scalar) -> usize {
    if t >= knots[n + 1] {
        return n;
    }
    let (mut lo, mut hi) = (degree, n + 1);
    let mut mid = lo.midpoint(hi);
    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = lo.midpoint(hi);
    }
    mid
}

/// Non-zero basis functions at `t`, by the Cox-de Boor recurrence.
pub(crate) fn basis_at(span: usize, t: Scalar, degree: usize, knots: &[Scalar]) -> Vec<Scalar> {
    let mut basis = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    basis[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = basis[r] / (right[r + 1] + left[j - r]);
            basis[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        basis[j] = saved;
    }
    basis
}

/// Solve `matrix * x = rhs` by Gaussian elimination with partial pivoting.
///
/// Pivoting is not optional here: the interpolation matrix is banded but its
/// diagonal is not guaranteed dominant for arbitrary point spacing, and
/// without pivoting a small pivot amplifies rounding into visible wobble.
/// A singular matrix is reported rather than producing infinities.
fn solve(mut matrix: Vec<Vec<Scalar>>, mut rhs: Vec<[Scalar; 3]>) -> GeomResult<Vec<[Scalar; 3]>> {
    let n = matrix.len();
    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&a, &b| matrix[a][column].abs().total_cmp(&matrix[b][column].abs()))
            .expect("range is non-empty");
        if matrix[pivot][column].abs() < 1e-12 {
            return Err(GeomError::Degenerate(
                "interpolation system is singular for these points".to_owned(),
            ));
        }
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);

        for row in (column + 1)..n {
            let factor = matrix[row][column] / matrix[column][column];
            if factor == 0.0 {
                continue;
            }
            for k in column..n {
                matrix[row][k] -= factor * matrix[column][k];
            }
            for axis in 0..3 {
                rhs[row][axis] -= factor * rhs[column][axis];
            }
        }
    }

    let mut solution = vec![[0.0; 3]; n];
    for row in (0..n).rev() {
        let mut accumulated = rhs[row];
        for (k, solved) in solution.iter().enumerate().take(n).skip(row + 1) {
            for (axis, value) in accumulated.iter_mut().enumerate() {
                *value -= matrix[row][k] * solved[axis];
            }
        }
        let pivot = matrix[row][row];
        for (axis, value) in solution[row].iter_mut().enumerate() {
            *value = accumulated[axis] / pivot;
        }
    }
    Ok(solution)
}

/// Assemble and solve the global interpolation system.
fn interpolate_with(points: &[Point3], parameters: &[Scalar]) -> GeomResult<BSplineCurve3> {
    let count = points.len();
    let degree = usize::from(degree_for(count));
    let expanded = averaged_knots(parameters, degree);
    let n = count - 1;

    // One row per point: the basis functions at its parameter must combine
    // the unknown control points into exactly that point.
    let mut matrix = vec![vec![0.0; count]; count];
    let mut rhs = vec![[0.0; 3]; count];
    for (row, (&t, point)) in parameters.iter().zip(points).enumerate() {
        let span = span_of(&expanded, n, degree, t);
        let basis = basis_at(span, t, degree, &expanded);
        for (offset, value) in basis.iter().enumerate() {
            matrix[row][span - degree + offset] = *value;
        }
        rhs[row] = [point.x, point.y, point.z];
    }

    let solved = solve(matrix, rhs)?;
    let control_points: Vec<Point3> = solved
        .into_iter()
        .map(|c| Point3::new(c[0], c[1], c[2]))
        .collect();

    // Collapse the expanded vector into distinct knots plus multiplicities,
    // which is the representation BSplineCurve carries.
    let (knots, multiplicities) = collapse(&expanded);

    Ok(BSplineCurve3 {
        degree: u16::try_from(degree)
            .map_err(|_| GeomError::InvalidInput("degree overflows".to_owned()))?,
        control_points,
        knots,
        multiplicities,
        weights: None,
        knot_spec: KnotSpec::Unspecified,
        closed: false,
        self_intersect: None,
    })
}

/// Group a repeated knot vector into distinct values and multiplicities.
pub(crate) fn collapse(expanded: &[Scalar]) -> (Vec<Scalar>, Vec<u32>) {
    let mut knots: Vec<Scalar> = Vec::new();
    let mut multiplicities: Vec<u32> = Vec::new();
    for &knot in expanded {
        if knots.last().is_some_and(|&last| last == knot) {
            *multiplicities
                .last_mut()
                .expect("knots and counts stay in step") += 1;
        } else {
            knots.push(knot);
            multiplicities.push(1);
        }
    }
    (knots, multiplicities)
}

/// Loft a surface through ordered section curves.
///
/// Each section becomes a row of the control net, so the surface passes
/// exactly through every section: at the section's own `u` parameter the
/// surface reduces to that curve. Section spacing along `u` is chord length
/// between corresponding control points, matching the curve case.
///
/// Sections must share a degree and control-point count. Reconciling
/// mismatched sections requires knot merging and degree elevation, and doing
/// that implicitly would change the caller's curves inside what looks like a
/// pure construction -- so it is refused, and the caller elevates explicitly
/// with `elevate_degree3`.
pub fn loft_surface(sections: &[BSplineCurve3]) -> GeomResult<BSplineSurface> {
    if sections.len() < 2 {
        return Err(GeomError::InvalidInput(
            "lofting needs at least two sections".to_owned(),
        ));
    }

    let first = &sections[0];
    let width = first.control_points.len();
    for (index, section) in sections.iter().enumerate() {
        if section.degree != first.degree {
            return Err(GeomError::InvalidInput(format!(
                "section {index} has degree {} but section 0 has degree {}; \
                 elevate explicitly rather than having the loft change your curves",
                section.degree, first.degree
            )));
        }
        if section.control_points.len() != width {
            return Err(GeomError::InvalidInput(format!(
                "section {index} has {} control points but section 0 has {width}; \
                 sections must share a control net width",
                section.control_points.len()
            )));
        }
        if section.weights.is_some() {
            return Err(GeomError::Unsupported {
                backend: axiolid_contracts::BackendId::new("nurbs"),
                operation: axiolid_contracts::Operation::SurfaceEvaluation,
            });
        }
    }

    // Space sections along u by the average chord between corresponding
    // control points: a section that sits far from its neighbour gets a
    // proportionally longer parameter interval, as in the curve case.
    let mut spans = vec![0.0];
    let mut total = 0.0;
    for pair in sections.windows(2) {
        let mean: Scalar = pair[0]
            .control_points
            .iter()
            .zip(&pair[1].control_points)
            .map(|(a, b)| (*b - *a).length())
            .sum::<Scalar>()
            / width as Scalar;
        if mean <= 0.0 {
            return Err(GeomError::Degenerate(
                "consecutive sections coincide, so loft spacing is undefined".to_owned(),
            ));
        }
        total += mean;
        spans.push(total);
    }
    let u_parameters: Vec<Scalar> = spans.into_iter().map(|d| d / total).collect();

    // Interpolate down each column of control points, so the surface passes
    // through every section rather than merely near it. Using the sections as
    // raw control rows would only approximate the interior ones.
    let u_degree = usize::from(degree_for(sections.len()));
    let u_expanded = averaged_knots(&u_parameters, u_degree);
    let rows = sections.len();
    let n = rows - 1;

    let mut matrix = vec![vec![0.0; rows]; rows];
    for (row, &t) in u_parameters.iter().enumerate() {
        let span = span_of(&u_expanded, n, u_degree, t);
        let basis = basis_at(span, t, u_degree, &u_expanded);
        for (offset, value) in basis.iter().enumerate() {
            matrix[row][span - u_degree + offset] = *value;
        }
    }

    let mut net: Vec<Vec<Point3>> = vec![Vec::with_capacity(width); rows];
    for column in 0..width {
        let rhs: Vec<[Scalar; 3]> = sections
            .iter()
            .map(|s| {
                let p = s.control_points[column];
                [p.x, p.y, p.z]
            })
            .collect();
        let solved = solve(matrix.clone(), rhs)?;
        for (row, coordinate) in solved.into_iter().enumerate() {
            net[row].push(Point3::new(coordinate[0], coordinate[1], coordinate[2]));
        }
    }

    let (u_knots, u_multiplicities) = collapse(&u_expanded);

    Ok(BSplineSurface {
        u_degree: u16::try_from(u_degree)
            .map_err(|_| GeomError::InvalidInput("u degree overflows".to_owned()))?,
        v_degree: first.degree,
        control_points: net,
        u_knots,
        u_multiplicities,
        // The v direction is the sections' own parameterisation, carried
        // through unchanged so the surface reproduces each section exactly.
        v_knots: first.knots.clone(),
        v_multiplicities: first.multiplicities.clone(),
        weights: None,
        u_closed: false,
        v_closed: first.closed,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    })
}
