//! Chord and derivative bounds of intrinsic curves from their curvature
//! laws (#317).
//!
//! An intrinsic curve is parameterised by its arc length `s`. In its start
//! frame's coordinates its tangent `T` is a unit vector and `T' = k N`
//! (Frenet: also for a space curve, whose torsion turns `N` but not the
//! size of `T'`), so `|c''| = |k|`; in 2D `c''' = k' N - k^2 T` gives
//! `|c'''| <= |k'| + k^2`.
//! The start frame maps those coordinates to the curve's, stretching by at
//! most [`frame_stretch2`] (`frame_stretch3`), so the suprema of `|k|` and
//! `|k'|` over a span (`crate::chain::curvature_bound`, the law's own
//! closed-form bound, term by term) give certified derivative suprema, and
//! `h^2 / 8 sup |c''|` the chord bound.
//!
//! A seam of a piecewise law is where `k` may jump. The heading, its
//! integral, stays continuous, so `c` is `C^1` with `c''` bounded across a
//! seam and the chord bound crosses it. `c'''` does not exist there: the
//! third-derivative bound refuses a span with a seam strictly inside.
//!
//! These bound the curve the law defines. Its points are read by
//! quadrature (`crate::arc_length`, `crate::frenet`), whose error is the
//! evaluator's stated estimate, as an arc-length chain's intrinsic pieces'
//! are (`crate::chain`): a caller that flattens one says so.

use axiolid_core::Scalar;
use axiolid_curve::{CurvatureLaw, Intrinsic2, Intrinsic3};

use super::{frame_stretch2, frame_stretch3, inflate, ordered, taylor, CurveDerivativeBounds};
use crate::chain::curvature_bound;

/// The span `[lo, hi]` of an intrinsic curve of `length`, or `None` outside
/// `[0, length]`.
fn span(a: Scalar, b: Scalar, length: Scalar) -> Option<(Scalar, Scalar)> {
    if !(a.is_finite() && b.is_finite() && length.is_finite() && length > 0.0) {
        return None;
    }
    let (lo, hi) = ordered(a, b);
    (lo >= 0.0 && hi <= length).then_some((lo, hi))
}

/// Certified derivative suprema of an intrinsic law placed by a frame of
/// `stretch` over `[lo, hi]`: `|c'| <= stretch`, `|c''| <= stretch K`,
/// `|c'''| <= stretch (K' + K^2)`, the last infinite across a seam.
fn law_bounds(
    law: &CurvatureLaw,
    length: Scalar,
    stretch: Scalar,
    lo: Scalar,
    hi: Scalar,
) -> Option<CurveDerivativeBounds> {
    let k = curvature_bound(law, lo, hi)?;
    let seamed = law
        .seams_within(length)
        .iter()
        .any(|&seam| seam > lo && seam < hi);
    let third = if seamed {
        Scalar::INFINITY
    } else {
        curvature_bound(&law.derivative(), lo, hi)? + k * k
    };
    Some(CurveDerivativeBounds {
        first: inflate(stretch),
        second: inflate(stretch * k),
        third: inflate(stretch * third),
    })
}

/// Certified derivative suprema of a 2D intrinsic curve over `[a, b]` of
/// its arc length; `None` outside `[0, length]` or for a law with no
/// closed-form bound. `third` is infinite over a span with a seam of the
/// law strictly inside.
#[must_use]
pub fn intrinsic_derivative_bounds2(
    curve: &Intrinsic2,
    a: Scalar,
    b: Scalar,
) -> Option<CurveDerivativeBounds> {
    let (lo, hi) = span(a, b, curve.length)?;
    law_bounds(
        &curve.curvature,
        curve.length,
        frame_stretch2(&curve.start),
        lo,
        hi,
    )
}

/// Certified chord bound of a 2D intrinsic curve over `[a, b]`:
/// `h^2 / 8 sigma sup |k|`, across seams of its law too. `None` outside
/// `[0, length]` or for a law with no closed-form bound.
#[must_use]
pub fn intrinsic_chord_bound2(curve: &Intrinsic2, a: Scalar, b: Scalar) -> Option<Scalar> {
    let bounds = intrinsic_derivative_bounds2(curve, a, b)?;
    let out = inflate(taylor(b - a, bounds.second));
    out.is_finite().then_some(out)
}

/// Certified chord bound of a 3D intrinsic curve over `[a, b]`: `h^2 / 8
/// sigma sup |k|`, whatever its torsion (`|c''| = |k|`). `None` outside
/// `[0, length]` or for a curvature law with no closed-form bound.
#[must_use]
pub fn intrinsic_chord_bound3(curve: &Intrinsic3, a: Scalar, b: Scalar) -> Option<Scalar> {
    let (lo, hi) = span(a, b, curve.length)?;
    let k = curvature_bound(&curve.curvature, lo, hi)?;
    let out = inflate(taylor(hi - lo, frame_stretch3(&curve.start) * k));
    out.is_finite().then_some(out)
}
