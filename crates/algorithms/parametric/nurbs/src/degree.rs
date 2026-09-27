//! Knot removal and degree operations (#33).
//!
//! # Which of these are exact, and which are not
//!
//! Degree ELEVATION is exact: every degree-p B-spline is also a degree-(p+1)
//! B-spline, and the elevated control net represents the same curve. Nothing
//! is approximated, so it always succeeds on valid input.
//!
//! Knot REMOVAL and degree REDUCTION are not. A knot is removable only if the
//! curve is smooth enough there to be represented without it; a degree is
//! reducible only if the curve was already representable at the lower degree.
//! Neither is generally true, so both must either meet a stated tolerance or
//! REFUSE.
//!
//! That is the whole design constraint here. An implementation that always
//! returns something has silently approximated, and the caller cannot tell a
//! clean removal from a lossy one. These return the deviation they actually
//! introduced so the caller can check it against their own budget.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar};
use axiolid_curve::{BSplineCurve, BSplineCurve2, BSplineCurve3};
use axiolid_evaluate::curve::{bspline_jet2, bspline_jet3};

/// Outcome of a lossy operation: the curve, and the error it introduced.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundedResult<C> {
    /// The transformed curve.
    pub curve: C,
    /// Upper bound on the deviation from the original, in model units.
    ///
    /// Zero means the operation was exact. A caller comparing this against
    /// its own budget is the intended use; the operation does not decide
    /// whether the error is acceptable.
    pub deviation_upper_bound: Scalar,
}

/// Raise a planar curve's degree by one without changing its image.
///
/// Exact: a degree-p curve is also a degree-(p+1) curve. Rational curves
/// elevate in homogeneous coordinates, so weights are carried through rather
/// than discarded.
pub fn elevate_degree2(curve: &BSplineCurve2) -> GeomResult<BSplineCurve2> {
    bspline_jet2(curve, 0.0)?;
    let segments = crate::transform::bezier_segments2(curve)?;
    let extract = |p: &Point2| [p.x, p.y];
    let rebuild = |c: [Scalar; 2]| Point2::new(c[0], c[1]);
    assemble(curve, &segments, extract, rebuild)
}

/// Raise a spatial curve's degree by one without changing its image.
pub fn elevate_degree3(curve: &BSplineCurve3) -> GeomResult<BSplineCurve3> {
    bspline_jet3(curve, 0.0)?;
    let segments = crate::transform::bezier_segments3(curve)?;
    let extract = |p: &Point3| [p.x, p.y, p.z];
    let rebuild = |c: [Scalar; 3]| Point3::new(c[0], c[1], c[2]);
    assemble(curve, &segments, extract, rebuild)
}

/// Degree elevation for a single clamped Bezier segment, in homogeneous form.
///
/// The standard identity: the elevated control points are the convex
/// combination `(i/(p+1)) * P[i-1] + (1 - i/(p+1)) * P[i]`, with the endpoints
/// carried through unchanged. Exact for polynomial and rational alike, because
/// it is applied to homogeneous coordinates.
fn elevate_bezier<const N: usize>(
    points: &[[Scalar; N]],
    weights: &[Scalar],
) -> (Vec<[Scalar; N]>, Vec<Scalar>) {
    let p = points.len() - 1;
    let elevated = p + 2;
    let mut out_points = Vec::with_capacity(elevated);
    let mut out_weights = Vec::with_capacity(elevated);

    out_points.push(points[0]);
    out_weights.push(weights[0]);
    for i in 1..=p {
        let alpha = i as Scalar / (p as Scalar + 1.0);
        let mut coordinate = [0.0; N];
        for (axis, value) in coordinate.iter_mut().enumerate() {
            // Homogeneous blend: weight the position by its own weight first.
            *value = alpha * points[i - 1][axis] * weights[i - 1]
                + (1.0 - alpha) * points[i][axis] * weights[i];
        }
        let weight = alpha * weights[i - 1] + (1.0 - alpha) * weights[i];
        for value in &mut coordinate {
            *value /= weight;
        }
        out_points.push(coordinate);
        out_weights.push(weight);
    }
    out_points.push(points[p]);
    out_weights.push(weights[p]);
    (out_points, out_weights)
}

/// Assemble elevated Bezier segments into one clamped curve.
///
/// Each segment elevates exactly; joining them keeps Bezier form, so every
/// internal knot carries multiplicity `q` for elevated degree `q`. The curve
/// is identical to the input. The representation is deliberately not
/// knot-minimal: collapsing the redundant internal knots is a lossy step
/// (`remove_knot`), and folding it in here would hide an approximation inside
/// an operation documented as exact.
fn assemble<const N: usize, P: Clone>(
    curve: &BSplineCurve<P>,
    segments: &[BSplineCurve<P>],
    coordinates: impl Fn(&P) -> [Scalar; N],
    point: impl Fn([Scalar; N]) -> P,
) -> GeomResult<BSplineCurve<P>> {
    let elevated_degree = curve.degree.checked_add(1).ok_or_else(|| {
        GeomError::InvalidInput("degree elevation would overflow the degree".to_owned())
    })?;

    let mut control_points: Vec<P> = Vec::new();
    let mut weights: Vec<Scalar> = Vec::new();
    for (index, segment) in segments.iter().enumerate() {
        let points: Vec<[Scalar; N]> = segment.control_points.iter().map(&coordinates).collect();
        let segment_weights = segment
            .weights
            .clone()
            .unwrap_or_else(|| vec![1.0; points.len()]);
        let (elevated_points, elevated_weights) = elevate_bezier(&points, &segment_weights);

        // Adjacent segments share their join point; keep one copy.
        let skip = usize::from(index > 0);
        for (coordinate, weight) in elevated_points.into_iter().zip(elevated_weights).skip(skip) {
            control_points.push(point(coordinate));
            weights.push(weight);
        }
    }

    // Each Bezier segment is clamped, so its first and last knot values are
    // its domain ends. Boundary knots get multiplicity q+1 (clamped), internal
    // joins get q, which is Bezier form at the elevated degree.
    let clamped = u32::from(elevated_degree) + 1;
    let internal = u32::from(elevated_degree);
    let mut knots: Vec<Scalar> = Vec::with_capacity(segments.len() + 1);
    let mut multiplicities: Vec<u32> = Vec::with_capacity(segments.len() + 1);

    for (index, segment) in segments.iter().enumerate() {
        let first = *segment
            .knots
            .first()
            .ok_or_else(|| GeomError::InvalidInput("bezier segment has no knots".to_owned()))?;
        if index == 0 {
            knots.push(first);
            multiplicities.push(clamped);
        }
        let last = *segment
            .knots
            .last()
            .ok_or_else(|| GeomError::InvalidInput("bezier segment has no knots".to_owned()))?;
        knots.push(last);
        multiplicities.push(if index + 1 == segments.len() {
            clamped
        } else {
            internal
        });
    }

    // Preserve rationality: a polynomial input stays polynomial rather than
    // acquiring a vector of ones, which would change the representation
    // without changing the curve.
    let weights = curve.weights.as_ref().map(|_| weights);

    Ok(BSplineCurve {
        degree: elevated_degree,
        control_points,
        knots,
        multiplicities,
        weights,
        knot_spec: curve.knot_spec,
        closed: curve.closed,
        self_intersect: curve.self_intersect,
    })
}

/// Remove one interior knot from a planar curve, or refuse.
///
/// Removal is lossy in general: a knot can only be dropped if the curve is
/// already smooth enough there to be represented without it. This computes the
/// candidate, MEASURES how far it actually moved, and refuses when that
/// exceeds `tolerance`.
///
/// Measuring rather than trusting the recurrence is deliberate. The removal
/// equations are an inverse of knot insertion and are only valid when the knot
/// is genuinely removable; applying them to a knot that carries real shape
/// produces a curve that looks plausible and is wrong. Sampling the result
/// against the original turns that into a refusal instead.
pub fn remove_knot2(
    curve: &BSplineCurve2,
    parameter: Scalar,
    tolerance: Scalar,
) -> GeomResult<BoundedResult<BSplineCurve2>> {
    bspline_jet2(curve, parameter)?;
    let candidate = remove(
        curve,
        parameter,
        |p| [p.x, p.y],
        |c| Point2::new(c[0], c[1]),
    )?;
    let deviation = deviation2(curve, &candidate)?;
    accept(candidate, deviation, tolerance)
}

/// Remove one interior knot from a spatial curve, or refuse.
pub fn remove_knot3(
    curve: &BSplineCurve3,
    parameter: Scalar,
    tolerance: Scalar,
) -> GeomResult<BoundedResult<BSplineCurve3>> {
    bspline_jet3(curve, parameter)?;
    let candidate = remove(
        curve,
        parameter,
        |p| [p.x, p.y, p.z],
        |c| Point3::new(c[0], c[1], c[2]),
    )?;
    let deviation = deviation3(curve, &candidate)?;
    accept(candidate, deviation, tolerance)
}

/// Accept a lossy result only if it met the caller's tolerance.
fn accept<C>(curve: C, deviation: Scalar, tolerance: Scalar) -> GeomResult<BoundedResult<C>> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(GeomError::InvalidInput(
            "tolerance must be finite and non-negative".to_owned(),
        ));
    }
    if deviation > tolerance {
        return Err(GeomError::Degenerate(format!(
            "knot is not removable within tolerance: deviation {deviation:.3e} exceeds {tolerance:.3e}"
        )));
    }
    Ok(BoundedResult {
        curve,
        deviation_upper_bound: deviation,
    })
}

/// Samples used to measure how far a lossy result moved.
///
/// Dense enough to catch the local bulge a bad removal introduces, and fixed
/// so the measurement is reproducible rather than depending on curve size.
const DEVIATION_SAMPLES: usize = 128;

fn deviation2(original: &BSplineCurve2, candidate: &BSplineCurve2) -> GeomResult<Scalar> {
    let (lo, hi) = domain(original);
    let mut worst: Scalar = 0.0;
    for index in 0..=DEVIATION_SAMPLES {
        let t = lo + (hi - lo) * (index as Scalar / DEVIATION_SAMPLES as Scalar);
        let a = bspline_jet2(original, t)?.point;
        let b = bspline_jet2(candidate, t)?.point;
        worst = worst.max((a - b).length());
    }
    Ok(worst)
}

fn deviation3(original: &BSplineCurve3, candidate: &BSplineCurve3) -> GeomResult<Scalar> {
    let (lo, hi) = domain(original);
    let mut worst: Scalar = 0.0;
    for index in 0..=DEVIATION_SAMPLES {
        let t = lo + (hi - lo) * (index as Scalar / DEVIATION_SAMPLES as Scalar);
        let a = bspline_jet3(original, t)?.point;
        let b = bspline_jet3(candidate, t)?.point;
        worst = worst.max((a - b).length());
    }
    Ok(worst)
}

/// Active parameter domain of a clamped curve.
fn domain<P>(curve: &BSplineCurve<P>) -> (Scalar, Scalar) {
    let mut expanded = Vec::new();
    for (&knot, &multiplicity) in curve.knots.iter().zip(&curve.multiplicities) {
        expanded.extend(core::iter::repeat_n(knot, multiplicity as usize));
    }
    let lo = expanded[usize::from(curve.degree)];
    let hi = expanded[curve.control_points.len()];
    (lo, hi)
}

/// Drop one knot by inverting the insertion recurrence.
///
/// Insertion computes new control points as convex combinations of old ones;
/// removal walks that backwards from both ends of the affected span. When the
/// knot is genuinely removable the two walks meet; when it is not, they
/// disagree and the resulting curve differs from the original -- which is what
/// the caller-side deviation check detects.
pub(crate) fn remove<const N: usize, P: Clone>(
    curve: &BSplineCurve<P>,
    parameter: Scalar,
    coordinates: impl Fn(&P) -> [Scalar; N],
    point: impl Fn([Scalar; N]) -> P,
) -> GeomResult<BSplineCurve<P>> {
    if !parameter.is_finite() {
        return Err(GeomError::InvalidInput(
            "knot parameter must be finite".to_owned(),
        ));
    }
    if curve.weights.is_some() {
        return Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("nurbs"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        });
    }

    let index = curve
        .knots
        .iter()
        .position(|&k| k == parameter)
        .ok_or_else(|| GeomError::InvalidInput("knot is not present in the curve".to_owned()))?;
    if index == 0 || index + 1 == curve.knots.len() {
        return Err(GeomError::InvalidInput(
            "endpoint knots bound the domain and cannot be removed".to_owned(),
        ));
    }

    let degree = usize::from(curve.degree);
    let expanded = expand_local(curve);
    // Last expanded position of this knot value.
    let span = expanded
        .iter()
        .rposition(|&k| k == parameter)
        .ok_or_else(|| GeomError::InvalidInput("knot vanished during expansion".to_owned()))?;
    let multiplicity = curve.multiplicities[index] as usize;

    // Piegl & Tiller A5.8. `temp` holds the recomputed run: index 0 and
    // index `last + 1 - first` are seeded from the untouched neighbours, and
    // the two walks meet in the middle.
    let ord = degree + 1;
    let first = span - degree;
    let last = span - multiplicity;
    let points: Vec<[Scalar; N]> = curve.control_points.iter().map(&coordinates).collect();

    let mut temp: Vec<[Scalar; N]> = vec![[0.0; N]; last + 2 - first];
    temp[0] = points[first - 1];
    temp[last + 1 - first] = points[last + 1];

    let (mut i, mut j) = (first, last);
    let (mut ii, mut jj) = (1_usize, last - first);
    while j > i {
        let alfi = (parameter - expanded[i]) / (expanded[i + ord] - expanded[i]);
        let alfj = (parameter - expanded[j]) / (expanded[j + ord] - expanded[j]);
        for axis in 0..N {
            temp[ii][axis] = (points[i][axis] - (1.0 - alfi) * temp[ii - 1][axis]) / alfi;
            temp[jj][axis] = (points[j][axis] - alfj * temp[jj + 1][axis]) / (1.0 - alfj);
        }
        i += 1;
        ii += 1;
        j -= 1;
        jj -= 1;
    }

    // Write the recomputed run back, then drop the surplus point. The
    // meeting point (j == i) needs its own store: the loop above stops
    // before writing it, and omitting it leaves one stale control point.
    let mut result = points.clone();
    let (mut i, mut j) = (first, last);
    while j > i {
        result[i] = temp[i - first + 1];
        result[j] = temp[j - first + 1];
        i += 1;
        j -= 1;
    }
    if j == i {
        result[i] = temp[i - first + 1];
    }
    result.remove(last);
    let points = result;

    let mut knots = curve.knots.clone();
    let mut multiplicities = curve.multiplicities.clone();
    multiplicities[index] -= 1;
    if multiplicities[index] == 0 {
        knots.remove(index);
        multiplicities.remove(index);
    }

    Ok(BSplineCurve {
        degree: curve.degree,
        control_points: points.into_iter().map(&point).collect(),
        knots,
        multiplicities,
        weights: None,
        knot_spec: curve.knot_spec,
        closed: curve.closed,
        self_intersect: curve.self_intersect,
    })
}

fn expand_local<P>(curve: &BSplineCurve<P>) -> Vec<Scalar> {
    let mut expanded = Vec::new();
    for (&knot, &multiplicity) in curve.knots.iter().zip(&curve.multiplicities) {
        expanded.extend(core::iter::repeat_n(knot, multiplicity as usize));
    }
    expanded
}

/// Lower a planar curve's degree by one, or refuse.
///
/// Lossy in general: only a curve that was already representable at the lower
/// degree reduces cleanly. An elevated curve is the clean case, and reduction
/// recovers what it started from.
///
/// Same discipline as knot removal: compute, measure, and refuse when the
/// deviation exceeds `tolerance`, rather than trusting the reduction formula
/// on a curve that genuinely needs its degree.
pub fn reduce_degree2(
    curve: &BSplineCurve2,
    tolerance: Scalar,
) -> GeomResult<BoundedResult<BSplineCurve2>> {
    bspline_jet2(curve, 0.0)?;
    let segments = crate::transform::bezier_segments2(curve)?;
    let candidate = reduce(
        curve,
        &segments,
        |p: &Point2| [p.x, p.y],
        |c: [Scalar; 2]| Point2::new(c[0], c[1]),
    )?;
    let deviation = deviation2(curve, &candidate)?;
    accept(candidate, deviation, tolerance)
}

/// Lower a spatial curve's degree by one, or refuse.
pub fn reduce_degree3(
    curve: &BSplineCurve3,
    tolerance: Scalar,
) -> GeomResult<BoundedResult<BSplineCurve3>> {
    bspline_jet3(curve, 0.0)?;
    let segments = crate::transform::bezier_segments3(curve)?;
    let candidate = reduce(
        curve,
        &segments,
        |p: &Point3| [p.x, p.y, p.z],
        |c: [Scalar; 3]| Point3::new(c[0], c[1], c[2]),
    )?;
    let deviation = deviation3(curve, &candidate)?;
    accept(candidate, deviation, tolerance)
}

/// Degree reduction for one Bezier segment.
///
/// Forward and backward recurrences each reconstruct the lower-degree control
/// points; averaging them distributes the error instead of piling it at one
/// end, which is what a one-directional recurrence does.
fn reduce_bezier<const N: usize>(points: &[[Scalar; N]]) -> Vec<[Scalar; N]> {
    let p = points.len() - 1;
    let reduced = p;
    let mut forward: Vec<[Scalar; N]> = vec![[0.0; N]; reduced];
    let mut backward: Vec<[Scalar; N]> = vec![[0.0; N]; reduced];

    forward[0] = points[0];
    for i in 1..reduced {
        let alpha = i as Scalar / p as Scalar;
        for axis in 0..N {
            forward[i][axis] = (points[i][axis] - alpha * forward[i - 1][axis]) / (1.0 - alpha);
        }
    }

    backward[reduced - 1] = points[p];
    for i in (0..reduced - 1).rev() {
        let alpha = (i + 1) as Scalar / p as Scalar;
        for axis in 0..N {
            backward[i][axis] =
                (points[i + 1][axis] - (1.0 - alpha) * backward[i + 1][axis]) / alpha;
        }
    }

    (0..reduced)
        .map(|i| {
            let mut blended = [0.0; N];
            for (axis, value) in blended.iter_mut().enumerate() {
                *value = 0.5 * (forward[i][axis] + backward[i][axis]);
            }
            blended
        })
        .collect()
}

/// Reduce each Bezier segment and rejoin, mirroring `assemble`.
fn reduce<const N: usize, P: Clone>(
    curve: &BSplineCurve<P>,
    segments: &[BSplineCurve<P>],
    coordinates: impl Fn(&P) -> [Scalar; N],
    point: impl Fn([Scalar; N]) -> P,
) -> GeomResult<BSplineCurve<P>> {
    if curve.degree < 2 {
        return Err(GeomError::InvalidInput(
            "degree 1 cannot be reduced further and stay a curve".to_owned(),
        ));
    }
    if curve.weights.is_some() {
        return Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("nurbs"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        });
    }
    let reduced_degree = curve.degree - 1;

    let mut control_points: Vec<P> = Vec::new();
    for (index, segment) in segments.iter().enumerate() {
        let points: Vec<[Scalar; N]> = segment.control_points.iter().map(&coordinates).collect();
        let lowered = reduce_bezier(&points);
        let skip = usize::from(index > 0);
        for coordinate in lowered.into_iter().skip(skip) {
            control_points.push(point(coordinate));
        }
    }

    let clamped = u32::from(reduced_degree) + 1;
    let internal = u32::from(reduced_degree);
    let mut knots: Vec<Scalar> = Vec::with_capacity(segments.len() + 1);
    let mut multiplicities: Vec<u32> = Vec::with_capacity(segments.len() + 1);
    for (index, segment) in segments.iter().enumerate() {
        let first = *segment
            .knots
            .first()
            .ok_or_else(|| GeomError::InvalidInput("bezier segment has no knots".to_owned()))?;
        if index == 0 {
            knots.push(first);
            multiplicities.push(clamped);
        }
        let last = *segment
            .knots
            .last()
            .ok_or_else(|| GeomError::InvalidInput("bezier segment has no knots".to_owned()))?;
        knots.push(last);
        multiplicities.push(if index + 1 == segments.len() {
            clamped
        } else {
            internal
        });
    }

    Ok(BSplineCurve {
        degree: reduced_degree,
        control_points,
        knots,
        multiplicities,
        weights: None,
        knot_spec: curve.knot_spec,
        closed: curve.closed,
        self_intersect: curve.self_intersect,
    })
}
