//! Exact shape-preserving NURBS transformations.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar};
use axiolid_curve::{BSplineCurve, BSplineCurve2, BSplineCurve3};
use axiolid_evaluate::curve::{bspline_jet2, bspline_jet3};

/// Reverse a planar B-spline curve without changing its image.
pub fn reverse2(curve: &BSplineCurve2) -> GeomResult<BSplineCurve2> {
    bspline_jet2(curve, 0.0)?;
    reverse(curve)
}

/// Reverse a spatial B-spline curve without changing its image.
pub fn reverse3(curve: &BSplineCurve3) -> GeomResult<BSplineCurve3> {
    bspline_jet3(curve, 0.0)?;
    reverse(curve)
}

/// Insert one interior knot into a planar curve in homogeneous coordinates.
///
/// The represented rational or polynomial curve is unchanged. Endpoint knots
/// and knots whose multiplicity already reaches the degree are rejected.
pub fn insert_knot2(curve: &BSplineCurve2, parameter: Scalar) -> GeomResult<BSplineCurve2> {
    bspline_jet2(curve, parameter)?;
    insert(
        curve,
        parameter,
        |p| [p.x, p.y],
        |p| Point2::new(p[0], p[1]),
    )
}

/// Insert one interior knot into a spatial curve in homogeneous coordinates.
///
/// The represented rational or polynomial curve is unchanged. Endpoint knots
/// and knots whose multiplicity already reaches the degree are rejected.
pub fn insert_knot3(curve: &BSplineCurve3, parameter: Scalar) -> GeomResult<BSplineCurve3> {
    bspline_jet3(curve, parameter)?;
    insert(
        curve,
        parameter,
        |p| [p.x, p.y, p.z],
        |p| Point3::new(p[0], p[1], p[2]),
    )
}

/// Split a planar curve exactly at an interior parameter.
///
/// Both output curves include the shared cut point and are marked open.
pub fn split2(
    curve: &BSplineCurve2,
    parameter: Scalar,
) -> GeomResult<(BSplineCurve2, BSplineCurve2)> {
    bspline_jet2(curve, parameter)?;
    let mut refined = curve.clone();
    check_interior(&refined, parameter)?;
    while multiplicity(&refined, parameter) < usize::from(refined.degree) {
        refined = insert_knot2(&refined, parameter)?;
    }
    split_ready(&refined, parameter)
}

/// Split a spatial curve exactly at an interior parameter.
///
/// Both output curves include the shared cut point and are marked open.
pub fn split3(
    curve: &BSplineCurve3,
    parameter: Scalar,
) -> GeomResult<(BSplineCurve3, BSplineCurve3)> {
    bspline_jet3(curve, parameter)?;
    let mut refined = curve.clone();
    check_interior(&refined, parameter)?;
    while multiplicity(&refined, parameter) < usize::from(refined.degree) {
        refined = insert_knot3(&refined, parameter)?;
    }
    split_ready(&refined, parameter)
}

/// Decompose a planar B-spline into exact rational/polynomial Bézier segments.
pub fn bezier_segments2(curve: &BSplineCurve2) -> GeomResult<Vec<BSplineCurve2>> {
    bspline_jet2(curve, 0.0)?;
    decompose(curve, split2)
}

/// Decompose a spatial B-spline into exact rational/polynomial Bézier segments.
pub fn bezier_segments3(curve: &BSplineCurve3) -> GeomResult<Vec<BSplineCurve3>> {
    bspline_jet3(curve, 0.0)?;
    decompose(curve, split3)
}

fn reverse<P: Clone>(curve: &BSplineCurve<P>) -> GeomResult<BSplineCurve<P>> {
    let (knots, multiplicities) = crate::axis::reverse_axis(&curve.knots, &curve.multiplicities)?;
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let weights = curve.weights.as_ref().map(|weights| {
        let mut reversed = weights.clone();
        reversed.reverse();
        reversed
    });
    Ok(BSplineCurve {
        degree: curve.degree,
        control_points,
        knots,
        multiplicities,
        weights,
        knot_spec: curve.knot_spec,
        closed: curve.closed,
        self_intersect: curve.self_intersect,
    })
}

fn check_interior<P>(curve: &BSplineCurve<P>, parameter: Scalar) -> GeomResult<()> {
    if !parameter.is_finite() {
        return Err(GeomError::InvalidInput(
            "split parameter must be finite".to_owned(),
        ));
    }
    let expanded = expand(curve);
    let lo = expanded[usize::from(curve.degree)];
    let hi = expanded[curve.control_points.len()];
    if parameter <= lo || parameter >= hi {
        return Err(GeomError::InvalidInput(
            "split parameter must lie strictly inside the active domain".to_owned(),
        ));
    }
    Ok(())
}

fn expand<P>(curve: &BSplineCurve<P>) -> Vec<Scalar> {
    let mut expanded =
        Vec::with_capacity(curve.control_points.len() + usize::from(curve.degree) + 1);
    for (&knot, &multiplicity) in curve.knots.iter().zip(&curve.multiplicities) {
        expanded.extend(core::iter::repeat_n(knot, multiplicity as usize));
    }
    expanded
}

fn multiplicity<P>(curve: &BSplineCurve<P>, parameter: Scalar) -> usize {
    curve
        .knots
        .iter()
        .position(|&k| k == parameter)
        .map_or(0, |i| curve.multiplicities[i] as usize)
}
fn split_ready<P: Clone>(
    curve: &BSplineCurve<P>,
    parameter: Scalar,
) -> GeomResult<(BSplineCurve<P>, BSplineCurve<P>)> {
    let expanded = expand(curve);
    let p = usize::from(curve.degree);
    let n = curve.control_points.len() - 1;
    let k = find_span(&expanded, n, p, parameter);
    let shared = k
        .checked_sub(p)
        .ok_or_else(|| GeomError::InvalidInput("split control index underflows".to_owned()))?;
    if shared == 0 || shared >= curve.control_points.len() - 1 {
        return Err(GeomError::InvalidInput(
            "split would create an empty segment".to_owned(),
        ));
    }
    let ki = curve
        .knots
        .iter()
        .position(|&value| value == parameter)
        .ok_or_else(|| {
            GeomError::InvalidInput("split knot is absent after refinement".to_owned())
        })?;
    let mut lm = curve.multiplicities[..=ki].to_vec();
    let mut rm = curve.multiplicities[ki..].to_vec();
    lm[ki] = u32::from(curve.degree) + 1;
    rm[0] = u32::from(curve.degree) + 1;
    let si = if curve.self_intersect == Some(false) {
        Some(false)
    } else {
        None
    };
    let make = |control_points: Vec<P>,
                knots: Vec<Scalar>,
                multiplicities: Vec<u32>,
                weights: Option<Vec<Scalar>>| BSplineCurve {
        degree: curve.degree,
        control_points,
        knots,
        multiplicities,
        weights,
        knot_spec: curve.knot_spec,
        closed: false,
        self_intersect: si,
    };
    let lw = curve.weights.as_ref().map(|w| w[..=shared].to_vec());
    let rw = curve.weights.as_ref().map(|w| w[shared..].to_vec());
    Ok((
        make(
            curve.control_points[..=shared].to_vec(),
            curve.knots[..=ki].to_vec(),
            lm,
            lw,
        ),
        make(
            curve.control_points[shared..].to_vec(),
            curve.knots[ki..].to_vec(),
            rm,
            rw,
        ),
    ))
}
type SplitFn<P> = fn(&BSplineCurve<P>, Scalar) -> GeomResult<(BSplineCurve<P>, BSplineCurve<P>)>;

fn decompose<P: Clone>(
    curve: &BSplineCurve<P>,
    split: SplitFn<P>,
) -> GeomResult<Vec<BSplineCurve<P>>> {
    let expanded = expand(curve);
    let lo = expanded[usize::from(curve.degree)];
    let hi = expanded[curve.control_points.len()];
    let internal: Vec<_> = curve
        .knots
        .iter()
        .copied()
        .filter(|&k| k > lo && k < hi)
        .collect();
    let mut result = Vec::with_capacity(internal.len() + 1);
    let mut remainder = curve.clone();
    for parameter in internal {
        let (left, right) = split(&remainder, parameter)?;
        result.push(left);
        remainder = right;
    }
    result.push(remainder);
    Ok(result)
}

pub(crate) fn insert<const N: usize, P: Clone>(
    curve: &BSplineCurve<P>,
    parameter: Scalar,
    coordinates: impl Fn(&P) -> [Scalar; N],
    point: impl Fn([Scalar; N]) -> P,
) -> GeomResult<BSplineCurve<P>> {
    if !parameter.is_finite() {
        return Err(GeomError::InvalidInput(
            "inserted knot must be finite".to_owned(),
        ));
    }
    let expanded = expand_knots(curve);
    let p = usize::from(curve.degree);
    let n = curve.control_points.len() - 1;
    let lo = expanded[p];
    let hi = expanded[n + 1];
    if parameter <= lo || parameter >= hi {
        return Err(GeomError::InvalidInput(format!(
            "inserted knot {parameter} must be strictly inside ({lo}, {hi})"
        )));
    }
    let k = find_span(&expanded, n, p, parameter);
    let s = expanded.iter().filter(|&&knot| knot == parameter).count();
    if s >= p {
        return Err(GeomError::InvalidInput(format!(
            "knot multiplicity {s} already reaches degree {p}"
        )));
    }
    let weights = curve.weights.clone().unwrap_or_else(|| vec![1.0; n + 1]);
    let homogeneous: Vec<_> = curve
        .control_points
        .iter()
        .zip(&weights)
        .map(|(control, &weight)| {
            let mut h = coordinates(control);
            for value in &mut h {
                *value *= weight;
            }
            (h, weight)
        })
        .collect();
    let mut output = vec![([0.0; N], 0.0); n + 2];
    output[..=k - p].clone_from_slice(&homogeneous[..=k - p]);
    output[k - s + 1..n + 2].copy_from_slice(&homogeneous[k - s..n + 1]);
    for i in k - p + 1..=k - s {
        let denominator = expanded[i + p] - expanded[i];
        if denominator == 0.0 {
            return Err(GeomError::Degenerate(
                "knot insertion denominator is zero".to_owned(),
            ));
        }
        let alpha = (parameter - expanded[i]) / denominator;
        let mut h = [0.0; N];
        for (d, value) in h.iter_mut().enumerate() {
            *value = alpha * homogeneous[i].0[d] + (1.0 - alpha) * homogeneous[i - 1].0[d];
        }
        output[i] = (
            h,
            alpha * homogeneous[i].1 + (1.0 - alpha) * homogeneous[i - 1].1,
        );
    }
    let mut new_expanded = expanded;
    new_expanded.insert(k + 1, parameter);
    let (knots, multiplicities) = compact(&new_expanded)?;
    let mut controls = Vec::with_capacity(output.len());
    let mut new_weights = Vec::with_capacity(output.len());
    for (mut h, weight) in output {
        if !weight.is_finite() || weight <= 0.0 {
            return Err(GeomError::Degenerate(
                "inserted homogeneous weight is not positive and finite".to_owned(),
            ));
        }
        for value in &mut h {
            *value /= weight;
        }
        controls.push(point(h));
        new_weights.push(weight);
    }
    Ok(BSplineCurve {
        degree: curve.degree,
        control_points: controls,
        knots,
        multiplicities,
        weights: curve.weights.as_ref().map(|_| new_weights),
        knot_spec: curve.knot_spec,
        closed: curve.closed,
        self_intersect: curve.self_intersect,
    })
}

fn expand_knots<P>(curve: &BSplineCurve<P>) -> Vec<Scalar> {
    let expected = curve.control_points.len() + usize::from(curve.degree) + 1;
    let mut expanded = Vec::with_capacity(expected);
    for (&knot, &multiplicity) in curve.knots.iter().zip(&curve.multiplicities) {
        expanded.extend(core::iter::repeat_n(knot, multiplicity as usize));
    }
    expanded
}

fn find_span(knots: &[Scalar], n: usize, degree: usize, parameter: Scalar) -> usize {
    if parameter >= knots[n + 1] {
        return n;
    }
    let mut low = degree;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;
    while parameter < knots[mid] || parameter >= knots[mid + 1] {
        if parameter < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

fn compact(expanded: &[Scalar]) -> GeomResult<(Vec<Scalar>, Vec<u32>)> {
    let mut knots = Vec::new();
    let mut multiplicities = Vec::new();
    for &knot in expanded {
        if knots.last().copied() == Some(knot) {
            *multiplicities
                .last_mut()
                .expect("knot and multiplicity stay parallel") += 1;
        } else {
            knots.push(knot);
            multiplicities.push(1);
        }
    }
    Ok((knots, multiplicities))
}
