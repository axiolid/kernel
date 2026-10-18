//! Points and chord bounds of constant offsets at their basis's native
//! parameter (#317), from the derivative bounds #294 and #307 derive.
//!
//! # The curve
//!
//! A constant offset of a basis `B` is `P(t) = B(t) + d u(t)`, `u = w /
//! |w|` a unit vector built from the basis's velocity:
//!
//! - beside a 2D curve, `w = Z x B'`, the left normal (`d` positive to the
//!   left of the curve's direction, as an offset curve 2D is);
//! - beside a 3D curve, `w = V x B'`, `V` the unit reference direction (an
//!   offset curve 3D), as `crate::station`'s directed offset law reads it.
//!
//! `P` takes the basis's parameter, so a trim or a composite of the basis
//! offsets span by span. A tangent that vanishes, or lies along `V`, leaves
//! `u` undefined: such a point is refused.
//!
//! # The chord bound
//!
//! `P'' = B'' + d u''`. On `[a, b]`, `h = b - a`, `|w|` falls at most `h /
//! 2 sup |B''|` from its midpoint value, which bounds it below by `r`, and
//! `|w'| <= |B''|`, `|w''| <= |B'''|` (a cross product with a unit vector
//! never lengthens). The unit-vector rule the certified station lengths use
//! (`crate::station::certify`'s `unit_jet`: `|u''| <= |w''| / r + 3 |w'|^2
//! / r^2`) then gives
//!
//! `sup |P''| <= sup |B''| + |d| (sup |B'''| / r + 3 sup |B''|^2 / r^2)`,
//!
//! and the chord bound is `h^2 / 8` times that, as for any curve
//! ([`super::chord_bound2`]). The basis's suprema come from
//! [`super::curve_derivative_bounds2`] (`3`), and from
//! [`super::intrinsic_derivative_bounds2`] beside a 2D intrinsic curve.
//! `P` is `C^1` only where `B` is `C^2`, so a span with such a break of the
//! basis strictly inside (a B-spline knot, a polyline vertex, a seam of a
//! curvature law) is refused; split there. A span of one polyline segment
//! offsets to a segment: its bound is `0`.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_curve::{Curve2, Curve3};

use super::{
    breaks_inside, continuity_breaks2, continuity_breaks3, curve_derivative_bounds2,
    curve_derivative_bounds3, inflate, intrinsic_derivative_bounds2, taylor, CurveDerivativeBounds,
    ROUNDING,
};
use crate::curve::{derivative2, derivative3, evaluate2, evaluate3};

/// The unit vector `w / |w|`, refused where `w` vanishes.
fn unit(w: Vec3) -> GeomResult<Vec3> {
    let size = w.length();
    if !(size.is_finite() && size > 0.0) {
        return Err(GeomError::Degenerate(
            "an offset curve whose basis tangent vanishes, or lies along its reference \
             direction: the offset direction is undefined there"
                .into(),
        ));
    }
    Ok(w / size)
}

/// `V` normalised, refused unless finite and non-zero.
fn reference(direction: Vec3) -> Option<Vec3> {
    let size = direction.length();
    (size.is_finite() && size > 0.0).then(|| direction / size)
}

/// The point of the offset by `distance` to the left of a 2D curve at its
/// parameter `t`: `B(t) + d normalise(Z x B'(t))`.
///
/// # Errors
///
/// The basis's evaluation errors, and [`GeomError::Degenerate`] where its
/// tangent vanishes.
pub fn offset_point2(curve: &Curve2, distance: Scalar, t: Scalar) -> GeomResult<Point2> {
    let p = evaluate2(curve, t)?;
    let d = derivative2(curve, t)?;
    let left = unit(Vec3::new(-d.y, d.x, 0.0))?;
    Ok(Point2::new(
        p.x + distance * left.x,
        p.y + distance * left.y,
    ))
}

/// The point of the offset by `distance` along `normalise(V x B'(t))` of a
/// 3D curve at its parameter `t`, `V` the reference direction.
///
/// # Errors
///
/// The basis's evaluation errors, and [`GeomError::Degenerate`] where its
/// tangent vanishes or lies along `V`, or `V` is zero or not finite.
pub fn offset_point3(
    curve: &Curve3,
    distance: Scalar,
    reference_direction: Vec3,
    t: Scalar,
) -> GeomResult<Point3> {
    let v = reference(reference_direction).ok_or_else(|| {
        GeomError::Degenerate("an offset curve's reference direction is zero".into())
    })?;
    let p = evaluate3(curve, t)?;
    let d = derivative3(curve, t)?;
    Ok(p + distance * unit(v.cross(d))?)
}

/// `h^2 / 8 sup |P''|` from the basis's suprema and `|w|` at the midpoint.
fn offset_chord(
    h: Scalar,
    bounds: CurveDerivativeBounds,
    middle: Scalar,
    distance: Scalar,
) -> Option<Scalar> {
    let h = h.abs();
    let r = middle * (1.0 - ROUNDING) - inflate(0.5 * h * bounds.second);
    if r.is_nan() || r <= 0.0 {
        return None;
    }
    let (x, y) = (bounds.second / r, bounds.third / r);
    let second = bounds.second + distance.abs() * (y + 3.0 * x * x);
    let out = inflate(taylor(h, inflate(second)));
    out.is_finite().then_some(out)
}

/// Certified chord bound of the offset of a 2D curve by `distance` (see
/// [`offset_point2`]) over `[a, b]` of the basis's parameter. `None` for a
/// basis with no derivative bound, a span with a break of the basis's
/// `C^2` continuity strictly inside, or one whose tangent may vanish.
#[must_use]
pub fn offset_chord_bound2(
    curve: &Curve2,
    distance: Scalar,
    a: Scalar,
    b: Scalar,
) -> Option<Scalar> {
    if !(a.is_finite() && b.is_finite() && distance.is_finite()) {
        return None;
    }
    if a == b {
        return Some(0.0);
    }
    if breaks_inside(&continuity_breaks2(curve, 2), a, b) {
        return None;
    }
    let bounds = match curve {
        Curve2::Polyline(_) => return Some(0.0),
        Curve2::Intrinsic(intrinsic) => intrinsic_derivative_bounds2(intrinsic, a, b)?,
        _ => curve_derivative_bounds2(curve, a, b)?,
    };
    let middle = derivative2(curve, 0.5 * (a + b)).ok()?.length();
    offset_chord(b - a, bounds, middle, distance)
}

/// Certified chord bound of the offset of a 3D curve by `distance` along
/// `V x B'` (see [`offset_point3`]) over `[a, b]` of the basis's
/// parameter; `None` as for [`offset_chord_bound2`], or where the tangent
/// may lie along `V`.
#[must_use]
pub fn offset_chord_bound3(
    curve: &Curve3,
    distance: Scalar,
    reference_direction: Vec3,
    a: Scalar,
    b: Scalar,
) -> Option<Scalar> {
    if !(a.is_finite() && b.is_finite() && distance.is_finite()) {
        return None;
    }
    let v = reference(reference_direction)?;
    if a == b {
        return Some(0.0);
    }
    if breaks_inside(&continuity_breaks3(curve, 2), a, b) {
        return None;
    }
    let bounds = match curve {
        Curve3::Polyline(_) => return Some(0.0),
        _ => curve_derivative_bounds3(curve, a, b)?,
    };
    let middle = v.cross(derivative3(curve, 0.5 * (a + b)).ok()?).length();
    offset_chord(b - a, bounds, middle, distance)
}
