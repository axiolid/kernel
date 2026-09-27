//! Knot removal, degree change and iso-curves of tensor-product surfaces
//! (#141).
//!
//! Each operation along `u` applies the curve operation to every column of
//! the control net (every row, along `v`), in homogeneous coordinates, so
//! rational surfaces are handled like polynomial ones. The `v` operations
//! transpose the surface, work along `u`, and transpose back.
//!
//! Degree elevation and iso-curves are exact. Knot removal and degree
//! reduction are not in general, so each bounds the deviation it introduced
//! and keeps it within the caller's tolerance. The bound is not sampled: the
//! result is refined back onto the original's knot vector (by knot insertion
//! or degree elevation, both exact), and the two control nets bound the
//! distance between the surfaces everywhere, by the convex-hull property of
//! the B-spline basis.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Point3, Scalar};
use axiolid_curve::{BSplineCurve, BSplineCurve3};
use axiolid_evaluate::surface::bspline_jet;
use axiolid_surface::BSplineSurface;

use crate::surface_transform::u_curve;

/// A homogeneous control point: `[w x, w y, w z, w]`.
type Homogeneous = [Scalar; 4];

/// A surface changed within a stated bound.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedSurface {
    /// The transformed surface.
    pub surface: BSplineSurface,
    /// Upper bound on the distance between the original and transformed
    /// surfaces at every parameter, in model units (computed in `f64`).
    pub deviation_upper_bound: Scalar,
}

/// The outcome of removing one knot value from a surface up to some number
/// of times.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceKnotRemoval {
    /// The surface with `removed` copies of the knot taken out; the input
    /// surface unchanged when `removed` is zero.
    pub surface: BSplineSurface,
    /// How many copies were removed.
    pub removed: u32,
    /// Upper bound on the distance between the input and `surface` at every
    /// parameter, in model units (computed in `f64`): the sum of the bounds
    /// of the single removals.
    pub deviation_upper_bound: Scalar,
}

/// Remove the `u` knot `parameter` up to `times` times, keeping the surface
/// within `tolerance` of the input.
///
/// A copy is removed only when it can be removed from every column of the
/// control net with the accumulated deviation bound still within
/// `tolerance`; removal stops at the first copy that cannot. A knot that
/// carries shape is therefore left in place (`removed == 0`), never
/// approximated away.
///
/// # Errors
///
/// An invalid surface, a tolerance that is negative or not finite, or a
/// `parameter` that is not an interior `u` knot.
pub fn remove_surface_knot_u(
    surface: &BSplineSurface,
    parameter: Scalar,
    times: u32,
    tolerance: Scalar,
) -> GeomResult<SurfaceKnotRemoval> {
    bspline_jet(surface, 0.0, 0.0)?;
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(GeomError::InvalidInput(
            "tolerance must be finite and non-negative".to_owned(),
        ));
    }
    let index = surface
        .u_knots
        .iter()
        .position(|&k| k == parameter)
        .ok_or_else(|| GeomError::InvalidInput(format!("{parameter} is not a u knot")))?;
    if index == 0 || index + 1 == surface.u_knots.len() {
        return Err(GeomError::InvalidInput(
            "end knots bound the domain and cannot be removed".to_owned(),
        ));
    }
    let rational = surface.weights.is_some();
    let mut columns = homogeneous_columns(surface);
    let mut removed = 0;
    let mut bound = 0.0;
    while removed < times && columns[0].knots.contains(&parameter) {
        let Some((candidate, step)) = remove_once(&columns, parameter, rational) else {
            break;
        };
        if bound + step > tolerance {
            break;
        }
        columns = candidate;
        bound += step;
        removed += 1;
    }
    let surface = if removed == 0 {
        surface.clone()
    } else {
        from_homogeneous_columns(surface, &columns)?
    };
    Ok(SurfaceKnotRemoval {
        surface,
        removed,
        deviation_upper_bound: bound,
    })
}

/// Remove the `v` knot `parameter` up to `times` times; as
/// [`remove_surface_knot_u`].
///
/// # Errors
///
/// As [`remove_surface_knot_u`].
pub fn remove_surface_knot_v(
    surface: &BSplineSurface,
    parameter: Scalar,
    times: u32,
    tolerance: Scalar,
) -> GeomResult<SurfaceKnotRemoval> {
    let mut result = remove_surface_knot_u(&transpose(surface), parameter, times, tolerance)?;
    result.surface = transpose(&result.surface);
    Ok(result)
}

/// Raise the surface's `u` degree by one without changing it.
///
/// Exact. Interior `u` knots come out with multiplicity equal to the new
/// degree (Bezier form), as for curves: collapsing them would be lossy, and
/// is [`remove_surface_knot_u`]'s job.
///
/// # Errors
///
/// An invalid surface.
pub fn elevate_surface_degree_u(surface: &BSplineSurface) -> GeomResult<BSplineSurface> {
    bspline_jet(surface, 0.0, 0.0)?;
    let columns = (0..surface.control_points[0].len())
        .map(|v| crate::degree::elevate_degree3(&u_curve(surface, v)))
        .collect::<GeomResult<Vec<_>>>()?;
    assemble_columns(surface, &columns)
}

/// Raise the surface's `v` degree by one; as [`elevate_surface_degree_u`].
///
/// # Errors
///
/// An invalid surface.
pub fn elevate_surface_degree_v(surface: &BSplineSurface) -> GeomResult<BSplineSurface> {
    Ok(transpose(&elevate_surface_degree_u(&transpose(surface))?))
}

/// Lower the surface's `u` degree by one, within `tolerance`, or refuse.
///
/// Only a surface already representable at the lower degree reduces
/// cleanly; one that needs its degree is refused. The result is in Bezier
/// form along `u`. Rational surfaces are refused, as for curves.
///
/// # Errors
///
/// An invalid or rational surface, a `u` degree below 2, a tolerance that is
/// negative or not finite, or a deviation bound above `tolerance`
/// (`GeomError::Degenerate`).
pub fn reduce_surface_degree_u(
    surface: &BSplineSurface,
    tolerance: Scalar,
) -> GeomResult<BoundedSurface> {
    bspline_jet(surface, 0.0, 0.0)?;
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(GeomError::InvalidInput(
            "tolerance must be finite and non-negative".to_owned(),
        ));
    }
    if surface.weights.is_some() {
        return Err(GeomError::Unsupported {
            backend: BackendId::new("nurbs"),
            operation: Operation::SurfaceEvaluation,
        });
    }
    if surface.u_degree < 2 {
        return Err(GeomError::InvalidInput(
            "a u degree of 1 cannot be reduced".to_owned(),
        ));
    }
    let mut reduced = Vec::with_capacity(surface.control_points[0].len());
    let mut bound: Scalar = 0.0;
    for v in 0..surface.control_points[0].len() {
        let column = u_curve(surface, v);
        // The per-column curve check samples; the bound used here is the
        // control-net one below, so the curve step accepts any deviation.
        let candidate = crate::degree::reduce_degree3(&column, Scalar::MAX)?.curve;
        // Elevating back is exact and gives Bezier form at the original
        // degree; refining the original to Bezier form is exact too. The
        // nets then share a knot vector, and their largest difference
        // bounds the curves' distance everywhere.
        let back = crate::degree::elevate_degree3(&candidate)?;
        let original = bezier_form(&column)?;
        if back.knots != original.knots
            || back.multiplicities != original.multiplicities
            || back.control_points.len() != original.control_points.len()
        {
            return Err(GeomError::Degenerate(
                "reduced surface does not refine onto the original knots".to_owned(),
            ));
        }
        for (a, b) in back.control_points.iter().zip(&original.control_points) {
            bound = bound.max(a.distance(*b));
        }
        reduced.push(candidate);
    }
    if bound.is_nan() || bound > tolerance {
        return Err(GeomError::Degenerate(format!(
            "u degree is not reducible within tolerance: deviation {bound:.3e} exceeds {tolerance:.3e}"
        )));
    }
    Ok(BoundedSurface {
        surface: assemble_columns(surface, &reduced)?,
        deviation_upper_bound: bound,
    })
}

/// Lower the surface's `v` degree by one; as [`reduce_surface_degree_u`].
///
/// # Errors
///
/// As [`reduce_surface_degree_u`].
pub fn reduce_surface_degree_v(
    surface: &BSplineSurface,
    tolerance: Scalar,
) -> GeomResult<BoundedSurface> {
    let mut result = reduce_surface_degree_u(&transpose(surface), tolerance)?;
    result.surface = transpose(&result.surface);
    Ok(result)
}

/// The iso-curve `u = parameter`, parametrised by `v`.
///
/// Exact: its control points are the columns' homogeneous points blended by
/// the `u` basis at `parameter`, so it has the surface's `v` degree and
/// knots, and is rational exactly when the surface is.
///
/// # Errors
///
/// An invalid surface, or a `parameter` outside the `u` domain.
pub fn iso_curve_at_u(surface: &BSplineSurface, parameter: Scalar) -> GeomResult<BSplineCurve3> {
    bspline_jet(surface, 0.0, 0.0)?;
    let degree = usize::from(surface.u_degree);
    let knots = expand(&surface.u_knots, &surface.u_multiplicities);
    let count = surface.control_points.len();
    let (lo, hi) = (knots[degree], knots[count]);
    if !(parameter >= lo && parameter <= hi) {
        return Err(GeomError::InvalidInput(format!(
            "u = {parameter} is outside the domain [{lo}, {hi}]"
        )));
    }
    let span = find_span(&knots, count, degree, parameter);
    let basis = basis_functions(&knots, span, degree, parameter);
    let rows = surface.control_points[0].len();
    let mut points = Vec::with_capacity(rows);
    let mut weights = Vec::with_capacity(rows);
    for v in 0..rows {
        let mut h = [0.0; 4];
        for (k, &n) in basis.iter().enumerate() {
            let i = span - degree + k;
            let p = surface.control_points[i][v];
            let w = surface.weights.as_ref().map_or(1.0, |rows| rows[i][v]);
            h[0] += n * w * p.x;
            h[1] += n * w * p.y;
            h[2] += n * w * p.z;
            h[3] += n * w;
        }
        if surface.weights.is_some() {
            if !(h[3].is_finite() && h[3] > 0.0) {
                return Err(GeomError::Degenerate(
                    "iso-curve weight is not positive and finite".to_owned(),
                ));
            }
            points.push(Point3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]));
            weights.push(h[3]);
        } else {
            points.push(Point3::new(h[0], h[1], h[2]));
        }
    }
    Ok(BSplineCurve {
        degree: surface.v_degree,
        control_points: points,
        knots: surface.v_knots.clone(),
        multiplicities: surface.v_multiplicities.clone(),
        weights: surface.weights.as_ref().map(|_| weights),
        knot_spec: surface.knot_spec,
        closed: surface.v_closed,
        self_intersect: surface.self_intersect,
    })
}

/// The iso-curve `v = parameter`, parametrised by `u`; as
/// [`iso_curve_at_u`].
///
/// # Errors
///
/// As [`iso_curve_at_u`].
pub fn iso_curve_at_v(surface: &BSplineSurface, parameter: Scalar) -> GeomResult<BSplineCurve3> {
    iso_curve_at_u(&transpose(surface), parameter)
}

/// Remove one copy of `parameter` from every column, and bound the
/// deviation; `None` when some column cannot take the removal.
fn remove_once(
    columns: &[BSplineCurve<Homogeneous>],
    parameter: Scalar,
    rational: bool,
) -> Option<(Vec<BSplineCurve<Homogeneous>>, Scalar)> {
    let mut candidate = Vec::with_capacity(columns.len());
    let mut refined = Vec::with_capacity(columns.len());
    for column in columns {
        let removed = crate::degree::remove(column, parameter, |h| *h, |h| h).ok()?;
        // A removal that drives a weight to zero or below leaves no valid
        // rational surface: the knot is not removable.
        if rational
            && removed
                .control_points
                .iter()
                .any(|h| !(h[3].is_finite() && h[3] > 0.0))
        {
            return None;
        }
        // Inserting the knot back is exact and restores the original knot
        // vector, so the two nets can be compared point for point.
        let back = crate::transform::insert(&removed, parameter, |h| *h, |h| h).ok()?;
        if back.control_points.len() != column.control_points.len() {
            return None;
        }
        candidate.push(removed);
        refined.push(back);
    }
    let bound = net_deviation(columns, &refined, rational)?;
    Some((candidate, bound))
}

/// A bound on the distance between two surfaces on one knot vector, from
/// their homogeneous nets.
///
/// With `H`, `h` the homogeneous points and `W`, `w` the weight sums of the
/// first and second surface, and any centre `c`:
/// `S1 - S2 = sum N (H - c W1 - (h - c w2) - (S2 - c)(W1 - w2)) / W1`, so
/// `|S1 - S2| <= max (|dH_c| + R |dw|) / min W1`, with `R` the radius about
/// `c` of the second net, which contains `S2`. Polynomial nets reduce to the
/// largest control-point difference.
fn net_deviation(
    first: &[BSplineCurve<Homogeneous>],
    second: &[BSplineCurve<Homogeneous>],
    rational: bool,
) -> Option<Scalar> {
    let pairs = || {
        first
            .iter()
            .zip(second)
            .flat_map(|(a, b)| a.control_points.iter().zip(&b.control_points))
    };
    if !rational {
        return Some(
            pairs()
                .map(|(a, b)| {
                    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
                    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
                })
                .fold(0.0, Scalar::max),
        );
    }
    let euclid = |h: &Homogeneous| [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
    let mut min = [Scalar::INFINITY; 3];
    let mut max = [Scalar::NEG_INFINITY; 3];
    for (_, b) in pairs() {
        let p = euclid(b);
        for k in 0..3 {
            min[k] = min[k].min(p[k]);
            max[k] = max[k].max(p[k]);
        }
    }
    let c = [0, 1, 2].map(|k| 0.5 * (min[k] + max[k]));
    let radius = pairs()
        .map(|(_, b)| {
            let p = euclid(b);
            ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt()
        })
        .fold(0.0, Scalar::max);
    let mut worst: Scalar = 0.0;
    let mut least_weight = Scalar::INFINITY;
    for (a, b) in pairs() {
        least_weight = least_weight.min(a[3]);
        let d = [0, 1, 2].map(|k| (a[k] - c[k] * a[3]) - (b[k] - c[k] * b[3]));
        let term = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() + radius * (a[3] - b[3]).abs();
        worst = worst.max(term);
    }
    if least_weight.is_nan() || least_weight <= 0.0 {
        return None;
    }
    let bound = worst / least_weight;
    bound.is_finite().then_some(bound)
}

/// Refine a curve to Bezier form: every interior knot at multiplicity equal
/// to the degree.
fn bezier_form(curve: &BSplineCurve3) -> GeomResult<BSplineCurve3> {
    let degree = u32::from(curve.degree);
    let mut refined = curve.clone();
    let interior = &curve.knots[1..curve.knots.len() - 1];
    for &knot in interior {
        loop {
            let index = refined
                .knots
                .iter()
                .position(|&k| k == knot)
                .expect("insertion keeps the knot");
            if refined.multiplicities[index] >= degree {
                break;
            }
            refined = crate::transform::insert_knot3(&refined, knot)?;
        }
    }
    Ok(refined)
}

/// Every column of the net, along `u`, in homogeneous coordinates.
fn homogeneous_columns(surface: &BSplineSurface) -> Vec<BSplineCurve<Homogeneous>> {
    (0..surface.control_points[0].len())
        .map(|v| BSplineCurve {
            degree: surface.u_degree,
            control_points: surface
                .control_points
                .iter()
                .enumerate()
                .map(|(u, row)| {
                    let p = row[v];
                    let w = surface.weights.as_ref().map_or(1.0, |rows| rows[u][v]);
                    [w * p.x, w * p.y, w * p.z, w]
                })
                .collect(),
            knots: surface.u_knots.clone(),
            multiplicities: surface.u_multiplicities.clone(),
            weights: None,
            knot_spec: surface.knot_spec,
            closed: surface.u_closed,
            self_intersect: surface.self_intersect,
        })
        .collect()
}

/// The surface from homogeneous columns sharing one knot vector.
fn from_homogeneous_columns(
    template: &BSplineSurface,
    columns: &[BSplineCurve<Homogeneous>],
) -> GeomResult<BSplineSurface> {
    let rational = template.weights.is_some();
    let u_count = columns[0].control_points.len();
    let mut control_points = vec![Vec::with_capacity(columns.len()); u_count];
    let mut weights = vec![Vec::with_capacity(columns.len()); u_count];
    for column in columns {
        for (u, h) in column.control_points.iter().enumerate() {
            if rational {
                if !(h[3].is_finite() && h[3] > 0.0) {
                    return Err(GeomError::Degenerate(
                        "a weight is not positive and finite".to_owned(),
                    ));
                }
                control_points[u].push(Point3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]));
                weights[u].push(h[3]);
            } else {
                // Polynomial: the weight coordinate is one, up to rounding.
                control_points[u].push(Point3::new(h[0], h[1], h[2]));
            }
        }
    }
    Ok(BSplineSurface {
        u_degree: columns[0].degree,
        v_degree: template.v_degree,
        control_points,
        u_knots: columns[0].knots.clone(),
        u_multiplicities: columns[0].multiplicities.clone(),
        v_knots: template.v_knots.clone(),
        v_multiplicities: template.v_multiplicities.clone(),
        weights: rational.then_some(weights),
        knot_spec: template.knot_spec,
        u_closed: template.u_closed,
        v_closed: template.v_closed,
        self_intersect: template.self_intersect,
    })
}

/// The surface from Euclidean (possibly rational) columns sharing one knot
/// vector and degree.
fn assemble_columns(
    template: &BSplineSurface,
    columns: &[BSplineCurve3],
) -> GeomResult<BSplineSurface> {
    let first = &columns[0];
    if columns.iter().any(|c| {
        c.degree != first.degree
            || c.knots != first.knots
            || c.multiplicities != first.multiplicities
            || c.control_points.len() != first.control_points.len()
    }) {
        return Err(GeomError::Degenerate(
            "columns disagree on their knot vector".to_owned(),
        ));
    }
    let u_count = first.control_points.len();
    let control_points = (0..u_count)
        .map(|u| columns.iter().map(|c| c.control_points[u]).collect())
        .collect();
    let weights = template.weights.as_ref().map(|_| {
        (0..u_count)
            .map(|u| {
                columns
                    .iter()
                    .map(|c| c.weights.as_ref().map_or(1.0, |w| w[u]))
                    .collect()
            })
            .collect()
    });
    Ok(BSplineSurface {
        u_degree: first.degree,
        v_degree: template.v_degree,
        control_points,
        u_knots: first.knots.clone(),
        u_multiplicities: first.multiplicities.clone(),
        v_knots: template.v_knots.clone(),
        v_multiplicities: template.v_multiplicities.clone(),
        weights,
        knot_spec: template.knot_spec,
        u_closed: template.u_closed,
        v_closed: template.v_closed,
        self_intersect: template.self_intersect,
    })
}

/// The same surface with `u` and `v` exchanged.
fn transpose(surface: &BSplineSurface) -> BSplineSurface {
    BSplineSurface {
        u_degree: surface.v_degree,
        v_degree: surface.u_degree,
        control_points: flip(&surface.control_points),
        u_knots: surface.v_knots.clone(),
        u_multiplicities: surface.v_multiplicities.clone(),
        v_knots: surface.u_knots.clone(),
        v_multiplicities: surface.u_multiplicities.clone(),
        weights: surface.weights.as_deref().map(flip),
        knot_spec: surface.knot_spec,
        u_closed: surface.v_closed,
        v_closed: surface.u_closed,
        self_intersect: surface.self_intersect,
    }
}

fn flip<T: Copy>(net: &[Vec<T>]) -> Vec<Vec<T>> {
    (0..net[0].len())
        .map(|v| net.iter().map(|row| row[v]).collect())
        .collect()
}

fn expand(knots: &[Scalar], multiplicities: &[u32]) -> Vec<Scalar> {
    let mut out = Vec::new();
    for (&k, &m) in knots.iter().zip(multiplicities) {
        out.extend(core::iter::repeat_n(k, m as usize));
    }
    out
}

/// The span `i` with `knots[i] <= t < knots[i + 1]`, the last non-empty
/// one at the domain's end.
fn find_span(knots: &[Scalar], count: usize, degree: usize, t: Scalar) -> usize {
    if t >= knots[count] {
        let mut span = count - 1;
        while span > degree && knots[span] == knots[count] {
            span -= 1;
        }
        return span;
    }
    let mut span = degree;
    while span + 1 < count && knots[span + 1] <= t {
        span += 1;
    }
    span
}

/// The `degree + 1` non-zero basis functions on `span` at `t` (Piegl and
/// Tiller A2.2).
fn basis_functions(knots: &[Scalar], span: usize, degree: usize, t: Scalar) -> Vec<Scalar> {
    let mut n = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}
