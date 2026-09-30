//! Hausdorff distance between polylines, and a certified Fréchet decision
//! (#147, ledger row H4).
//!
//! The Hausdorff distance treats both curves as continuous piecewise-linear
//! point sets (not just their vertices). The one-sided distance from `a` to
//! `b` is the supremum, over every point of `a`, of that point's distance to
//! the nearest point of `b`; the two-sided distance is the larger of the two
//! one-sided distances. Unlike the Fréchet distance, it has no notion of a
//! monotone walk, so it can be much smaller than the Fréchet distance for
//! curves that pass close to each other out of step, and it is not
//! symmetric in general -- see [`one_sided_polyline_hausdorff_distance`].
//!
//! # Algorithm
//!
//! For one segment of `a`, parametrised `a(t) = p0 + t * d` for `t` in
//! `[0, 1]`, the distance to `b` is `f(t) = min` over every feature of `b`
//! (a vertex, or a segment's infinite line restricted to where the foot of
//! `a(t)` falls inside the segment) of `dist(a(t), feature)`. Distance
//! squared to a fixed point, or to an infinite line, is a quadratic function
//! of `t` (an affine image of a convex distance), so on any sub-interval
//! where a single feature is uniquely the minimum, `f` is convex there and
//! has no interior local maximum: its maximum on such a sub-interval sits at
//! an endpoint. The candidates for the segment's supremum are therefore the
//! segment's own endpoints, `t = 0` and `t = 1`, plus every `t` where two
//! features of `b` are simultaneously the minimum -- a tie found by solving
//! a (at most) quadratic equation per pair of features, restricted to any
//! line feature's activity window. At every candidate, the true minimum
//! distance over all of `b`'s features is evaluated (ordinary clamped
//! point-to-segment distance, minimised over every segment of `b`), and the
//! largest of these is the segment's contribution. This is `O(m^2)` per
//! segment of `a` against `m` segments of `b`, so `O(n * m^2)` for one
//! direction -- worth revisiting if a crate consumer needs it faster, but
//! correctness came first here.
//!
//! # Rounding
//!
//! Everything is `f64`. Tie candidates are roots of quadratics solved in
//! floating point, so a genuinely exact tie can be missed, or found a few
//! ULPs from its true parameter; the evaluation step then reads the true
//! minimum at that (slightly off) point, so the reported supremum can be a
//! few ULPs short of the exact answer. Nothing here is certified.

use axiolid_core::{Point2, Point3};

use crate::frechet::{check, decide, distance, foot, lift, FrechetError};

/// The supremum, over every point of polyline `a`, of that point's distance
/// to the nearest point of polyline `b`. Both curves are treated as
/// continuous piecewise-linear point sets, not just their vertices.
///
/// This is not symmetric: `a` can lie entirely close to `b` while `b` has
/// stretches far from `a` (for example, a short segment sitting directly
/// above the middle of a much longer one). Use [`polyline_hausdorff_distance`] for
/// the symmetric, two-sided distance.
///
/// # Errors
///
/// An empty polyline or a non-finite coordinate.
pub fn one_sided_polyline_hausdorff_distance(
    a: &[Point3],
    b: &[Point3],
) -> Result<f64, FrechetError> {
    check(a)?;
    check(b)?;
    if a.len() == 1 {
        return Ok(point_to_polyline_distance(a[0], b));
    }
    let sup = a
        .windows(2)
        .map(|s| segment_sup(s[0], s[1], b))
        .fold(f64::NEG_INFINITY, f64::max);
    Ok(sup)
}

/// [`one_sided_polyline_hausdorff_distance`] in the plane.
///
/// # Errors
///
/// As [`one_sided_polyline_hausdorff_distance`].
pub fn one_sided_polyline_hausdorff_distance_2d(
    a: &[Point2],
    b: &[Point2],
) -> Result<f64, FrechetError> {
    one_sided_polyline_hausdorff_distance(&lift(a), &lift(b))
}

/// The symmetric Hausdorff distance between polylines `a` and `b`: the
/// larger of the two one-sided distances (see [`one_sided_polyline_hausdorff_distance`]).
///
/// # Errors
///
/// An empty polyline or a non-finite coordinate.
pub fn polyline_hausdorff_distance(a: &[Point3], b: &[Point3]) -> Result<f64, FrechetError> {
    let a_to_b = one_sided_polyline_hausdorff_distance(a, b)?;
    let b_to_a = one_sided_polyline_hausdorff_distance(b, a)?;
    Ok(a_to_b.max(b_to_a))
}

/// [`polyline_hausdorff_distance`] in the plane.
///
/// # Errors
///
/// As [`polyline_hausdorff_distance`].
pub fn polyline_hausdorff_distance_2d(a: &[Point2], b: &[Point2]) -> Result<f64, FrechetError> {
    polyline_hausdorff_distance(&lift(a), &lift(b))
}

/// Distance from `p` to the nearest point of polyline `b` (a single point
/// when `b` has one vertex).
fn point_to_polyline_distance(p: Point3, b: &[Point3]) -> f64 {
    if b.len() == 1 {
        return distance(p, b[0]);
    }
    b.windows(2)
        .map(|s| point_segment_distance(p, s[0], s[1]))
        .fold(f64::INFINITY, f64::min)
}

/// Ordinary clamped point-to-segment distance.
fn point_segment_distance(p: Point3, a: Point3, b: Point3) -> f64 {
    match foot(p, a, b) {
        None => distance(p, a),
        Some((t, h)) => {
            if t <= 0.0 {
                distance(p, a)
            } else if t >= 1.0 {
                distance(p, b)
            } else {
                h
            }
        }
    }
}

/// `(c1, c0)` such that `dist(a(t), v)^2 = c2 * t^2 + c1 * t + c0`, where
/// `a(t) = p0 + t * d` and `c2 = d.dot(d)` is shared by every point feature
/// of the same segment.
fn point_linear_coeffs(p0: Point3, d: Point3, v: Point3) -> (f64, f64) {
    let diff = p0 - v;
    (2.0 * d.dot(diff), diff.dot(diff))
}

/// `(c2, c1, c0)` such that the squared distance from `a(t) = p0 + t * d` to
/// the infinite line through segment `(u, w)` is `c2 * t^2 + c1 * t + c0`.
/// `None` for a degenerate (zero-length) segment, which has no line feature.
fn line_quadratic_coeffs(p0: Point3, d: Point3, u: Point3, w: Point3) -> Option<(f64, f64, f64)> {
    let e = w - u;
    let e2 = e.dot(e);
    if e2 == 0.0 {
        return None;
    }
    let q0 = p0 - u;
    // |q(t)|^2 coefficients, q(t) = q0 + t * d.
    let big_a = q0.dot(q0);
    let big_b = 2.0 * d.dot(q0);
    let big_c = d.dot(d);
    // dot(q(t), e) coefficients (linear).
    let big_e = q0.dot(e);
    let big_f = d.dot(e);
    Some((
        big_c - big_f * big_f / e2,
        big_b - 2.0 * big_e * big_f / e2,
        big_a - big_e * big_e / e2,
    ))
}

/// Whether the foot of `a(t) = p0 + t * d` onto segment `(u, w)` lands
/// inside the segment (parameter in `[0, 1]`). `false` for a degenerate
/// segment, which has no line feature to be active.
fn foot_active(p0: Point3, d: Point3, t: f64, u: Point3, w: Point3) -> bool {
    match foot(p0 + d * t, u, w) {
        Some((s, _)) => (0.0..=1.0).contains(&s),
        None => false,
    }
}

/// Roots of `c2 * t^2 + c1 * t + c0 = 0`, handling the near-linear and
/// no-real-root cases. Deliberately permissive: any spurious root is just an
/// extra evaluation point, which cannot corrupt the supremum.
fn quadratic_roots(c2: f64, c1: f64, c0: f64) -> Vec<f64> {
    let scale = c2.abs().max(c1.abs()).max(c0.abs());
    if scale == 0.0 {
        return Vec::new();
    }
    let threshold = 1e-12 * scale;
    if c2.abs() <= threshold {
        return if c1.abs() <= threshold {
            Vec::new()
        } else {
            vec![-c0 / c1]
        };
    }
    let discriminant = c1 * c1 - 4.0 * c2 * c0;
    if discriminant < 0.0 {
        return Vec::new();
    }
    let root = discriminant.sqrt();
    vec![(-c1 - root) / (2.0 * c2), (-c1 + root) / (2.0 * c2)]
}

/// The supremum of the distance to polyline `b`, over the segment from `p0`
/// to `p1` of `a`. See the module's `# Algorithm` section.
fn segment_sup(p0: Point3, p1: Point3, b: &[Point3]) -> f64 {
    let d = p1 - p0;
    let mut candidates: Vec<f64> = vec![0.0, 1.0];

    if b.len() >= 2 {
        let c2 = d.dot(d);
        // Vertex-vertex ties: same c2, so the quadratic collapses to linear.
        for (i, &vi) in b.iter().enumerate() {
            let (c1i, c0i) = point_linear_coeffs(p0, d, vi);
            for &vj in b.iter().skip(i + 1) {
                let (c1j, c0j) = point_linear_coeffs(p0, d, vj);
                let slope = c1i - c1j;
                if slope != 0.0 {
                    let t = (c0j - c0i) / slope;
                    if (0.0..=1.0).contains(&t) {
                        candidates.push(t);
                    }
                }
            }
        }
        let segments: Vec<(Point3, Point3)> = b.windows(2).map(|s| (s[0], s[1])).collect();
        let lines: Vec<Option<(f64, f64, f64)>> = segments
            .iter()
            .map(|&(u, w)| line_quadratic_coeffs(p0, d, u, w))
            .collect();
        for (j, &(u, w)) in segments.iter().enumerate() {
            let Some((c2l, c1l, c0l)) = lines[j] else {
                continue;
            };
            // Vertex-vs-line ties.
            for &v in b {
                let (c1v, c0v) = point_linear_coeffs(p0, d, v);
                for t in quadratic_roots(c2 - c2l, c1v - c1l, c0v - c0l) {
                    if (0.0..=1.0).contains(&t) && foot_active(p0, d, t, u, w) {
                        candidates.push(t);
                    }
                }
            }
            // Line-vs-line ties.
            for (k, &(u2, w2)) in segments.iter().enumerate().skip(j + 1) {
                let Some((c2l2, c1l2, c0l2)) = lines[k] else {
                    continue;
                };
                for t in quadratic_roots(c2l - c2l2, c1l - c1l2, c0l - c0l2) {
                    if (0.0..=1.0).contains(&t)
                        && foot_active(p0, d, t, u, w)
                        && foot_active(p0, d, t, u2, w2)
                    {
                        candidates.push(t);
                    }
                }
            }
        }
    }

    candidates
        .into_iter()
        .map(|t| point_to_polyline_distance(p0 + d * t.clamp(0.0, 1.0), b))
        .fold(f64::NEG_INFINITY, f64::max)
}

/// A certified answer to the Fréchet decision problem: whether the
/// continuous Fréchet distance is at most `eps`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrechetDecision {
    /// The Fréchet distance is certainly at most `eps`.
    AtMost,
    /// The Fréchet distance is certainly more than `eps`.
    MoreThan,
    /// `eps` is within the computed error margin of the true distance: the
    /// sign cannot be certified from `f64` arithmetic alone.
    Undecided,
}

/// A conservative bound on the absolute error in the decision's distance and
/// foot computations, for coordinates up to magnitude `r` and a query up to
/// `eps`. `k = 64.0` is a deliberately generous constant covering the
/// handful of subtractions, dot products and one square root in `distance`
/// and `foot`'s short, fixed-length chains, not a tight one.
fn error_margin(a: &[Point3], b: &[Point3], eps: f64) -> f64 {
    const K: f64 = 64.0;
    let r = a
        .iter()
        .chain(b.iter())
        .flat_map(|p| [p.x, p.y, p.z])
        .fold(0.0_f64, |m, c| m.max(c.abs()));
    let r = if r == 0.0 { 1.0 } else { r };
    K * (r + eps.abs() + 1.0) * f64::EPSILON
}

/// A certified decision for whether the continuous Fréchet distance between
/// `a` and `b` is at most `eps`, unlike [`frechet_at_most`](crate::frechet_at_most)
/// (cheaper, but its module docs say plainly it is not certified).
///
/// # Design
///
/// The decision is monotone non-decreasing in `eps`, and every distance or
/// foot value it reads is within a conservative bound `margin` of the true
/// value. If the *computed* free space already permits a monotone
/// path at the smaller leash `eps - margin`, the *true* free space at `eps`
/// -- a superset, since a true distance is at most `margin` above its
/// computed value -- permits it too, so the answer is certainly
/// [`FrechetDecision::AtMost`]. Symmetrically, if even the *computed* free
/// space at the larger leash `eps + margin` -- a superset of the true free
/// space at `eps` -- fails to permit a path, the true one certainly does
/// not, so the answer is certainly [`FrechetDecision::MoreThan`]. Otherwise
/// the true distance is provably within `margin` of `eps` and the decision
/// is [`FrechetDecision::Undecided`].
///
/// # Errors
///
/// An empty polyline, a non-finite coordinate, or a negative or non-finite
/// `eps`.
pub fn frechet_decide_certified(
    a: &[Point3],
    b: &[Point3],
    eps: f64,
) -> Result<FrechetDecision, FrechetError> {
    check(a)?;
    check(b)?;
    if !eps.is_finite() || eps < 0.0 {
        return Err(FrechetError::InvalidBound);
    }
    let margin = error_margin(a, b, eps);
    if eps >= margin && decide(a, b, eps - margin) {
        return Ok(FrechetDecision::AtMost);
    }
    if !decide(a, b, eps + margin) {
        return Ok(FrechetDecision::MoreThan);
    }
    Ok(FrechetDecision::Undecided)
}

/// [`frechet_decide_certified`] in the plane.
///
/// # Errors
///
/// As [`frechet_decide_certified`].
pub fn frechet_decide_certified_2d(
    a: &[Point2],
    b: &[Point2],
    eps: f64,
) -> Result<FrechetDecision, FrechetError> {
    frechet_decide_certified(&lift(a), &lift(b), eps)
}
