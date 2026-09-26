//! Sections of analytic surfaces with no closed form, as implicit curves
//! (ADR 0077).
//!
//! The section, read in one surface's parameters, is the zero set of the
//! other surface's equation composed with the first one's point: an exact
//! [`Field2`]. [`implicit_surface_intersection`] traces every component of
//! it in a parameter window with certified topology
//! (`implicit_trace`), each as an [`ImplicitSection3`] on the first
//! surface. A torus against a cylinder, cone or torus off its axis takes
//! this path from `exact_surface_intersection`.

use axiolid_core::{Interval, Point2, Scalar, Vec3};
use axiolid_curve::{Carrier, Curve3, Field2, ImplicitSection3};
use axiolid_surface::Surface;
use core::f64::consts::{FRAC_PI_2, PI};

use crate::exact_surface_intersection::{
    Derivation, ExactIntersectionCurve, ExactIntersectionRefusal,
};
use crate::field::{carrier_of, section_field};
use crate::implicit_trace::{trace, Periodic, TraceRefusal};
use axiolid_curve::implicit::Cell;

/// `other`'s implicit equation read in `carrier`'s parameters, where both
/// are analytic: zero exactly where `carrier`'s point lies on `other` (for
/// a cone, on either nappe).
#[must_use]
pub fn section_field_of(carrier: &Surface, other: &Surface) -> Option<Field2> {
    section_field(&carrier_of(carrier)?, other)
}

/// The refusal a trace's own refusal amounts to: touching only at isolated
/// points is `NotRegularCurve`, as for the closed forms; anything else is
/// `Undecided`.
pub(crate) fn refusal_of(refusal: TraceRefusal) -> ExactIntersectionRefusal {
    match refusal {
        TraceRefusal::Touching(_) => ExactIntersectionRefusal::NotRegularCurve,
        TraceRefusal::Singular(_) | TraceRefusal::Many(_) | TraceRefusal::Budget => {
            ExactIntersectionRefusal::Undecided
        }
    }
}

/// Every component of the section of `carrier` by `other` in a window of
/// `carrier`'s parameters, as implicit curves on `carrier`.
///
/// Without a `window`, the whole of a compact carrier (a sphere or torus)
/// is searched, or the stretch of an unbounded one that a compact `other`
/// can reach.
///
/// # Errors
///
/// - `UnsupportedPair`: a B-spline operand, or an unbounded pair with no
///   window.
/// - `Undecided`: the trace ran out of budget, or met a singular point it
///   could not match to the field's sign changes about it.
/// - `NotRegularCurve`: the surfaces only touch in the window (a
///   singular point of the section), or the trace exhausted its budget.
/// - `Disjoint`: no section in the window.
pub fn implicit_surface_intersection(
    carrier: &Surface,
    other: &Surface,
    window: Option<(Point2, Point2)>,
) -> Result<Vec<ImplicitSection3>, ExactIntersectionRefusal> {
    let form = carrier_of(carrier).ok_or(ExactIntersectionRefusal::UnsupportedPair)?;
    let field = section_field(&form, other).ok_or(ExactIntersectionRefusal::UnsupportedPair)?;
    let (periodic_u, periodic_v) = form.periodic();
    let (lo, hi) = match window {
        Some(w) => w,
        None => default_window(&form, other).ok_or(ExactIntersectionRefusal::UnsupportedPair)?,
    };
    // A window a whole turn wide in a periodic parameter wraps.
    let wraps = |periodic: bool, a: Scalar, b: Scalar| {
        periodic && ((b - a) - 2.0 * PI).abs() <= 1e-12 * (1.0 + a.abs())
    };
    let periodic = Periodic {
        u: wraps(periodic_u, lo.x, hi.x),
        v: wraps(periodic_v, lo.y, hi.y),
    };
    let curves = trace(&field, Cell { lo, hi }, periodic).map_err(refusal_of)?;
    // A B-spline ends at its domain: its sections do too.
    let curves: Vec<_> = match &form {
        Carrier::Spline(b) => match b.domain() {
            Some(((u0, u1), (v0, v1))) => curves
                .iter()
                .flat_map(|c| c.clipped(Point2::new(u0, v0), Point2::new(u1, v1)))
                .collect(),
            None => curves,
        },
        _ => curves,
    };
    if curves.is_empty() {
        return Err(ExactIntersectionRefusal::Disjoint);
    }
    Ok(curves
        .into_iter()
        .map(|curve| ImplicitSection3 {
            carrier: form.clone(),
            curve,
        })
        .collect())
}

/// Where a whole-turn window starts, past `-pi`: an irrational fraction of
/// a radian, so the window's edge is never where a symmetric section
/// crosses or turns.
pub(crate) const TURN_OFFSET: Scalar = 0.123_456_789_012_345_67;

/// A centre and radius holding a compact surface.
fn ball(surface: &Surface) -> Option<(Vec3, Scalar)> {
    match surface {
        Surface::Sphere(s) => Some((s.frame.origin, s.radius)),
        Surface::Torus(t) => Some((t.frame.origin, t.major_radius + t.minor_radius)),
        _ => None,
    }
}

/// The window a trace needs on `carrier` to find every section with
/// `other`.
fn default_window(carrier: &Carrier, other: &Surface) -> Option<(Point2, Point2)> {
    // A whole turn, starting off the round angles a symmetric input puts
    // its special points on.
    let o = TURN_OFFSET;
    match carrier {
        Carrier::Torus(_) => Some((Point2::new(-PI + o, -PI + o), Point2::new(PI + o, PI + o))),
        Carrier::Sphere { .. } => Some((
            Point2::new(-PI + o, -FRAC_PI_2),
            Point2::new(PI + o, FRAC_PI_2),
        )),
        Carrier::Ruled(k) => {
            let (centre, radius) = ball(other)?;
            let z = k.frame.z.normalize();
            let along = (centre - k.frame.origin).dot(z);
            let margin = 1e-6 * (1.0 + radius + along.abs());
            let (mut v0, mut v1) = (along - radius - margin, along + radius + margin);
            // One nappe of a cone: stop short of the apex.
            if k.slope != 0.0 {
                let apex = -k.x_radius / k.slope;
                let clear = 1e-9 * (1.0 + apex.abs());
                if k.slope > 0.0 {
                    v0 = v0.max(apex + clear);
                } else {
                    v1 = v1.min(apex - clear);
                }
                if v0 >= v1 {
                    return None;
                }
            }
            Some((Point2::new(-PI + o, v0), Point2::new(PI + o, v1)))
        }
        Carrier::Spline(b) => {
            // A little past the domain, so a section reaching its edge is
            // found there and not on the window's own boundary.
            let ((u0, u1), (v0, v1)) = b.domain()?;
            let (pu, pv) = (1e-6 * (u1 - u0), 1e-6 * (v1 - v0));
            Some((Point2::new(u0 - pu, v0 - pv), Point2::new(u1 + pu, v1 + pv)))
        }
        Carrier::Plane(f) => {
            let (centre, radius) = ball(other)?;
            let d = centre - f.origin;
            let (x, y) = (d.dot(f.x.normalize()), d.dot(f.y.normalize()));
            let r = radius * (1.0 + 1e-6);
            Some((Point2::new(x - r, y - r), Point2::new(x + r, y + r)))
        }
    }
}

/// The traced fallback of `exact_surface_intersection`: analytic pairs with
/// no closed form, carried on the compact surface (a torus first). `None`
/// when the pair has no compact surface to carry the trace.
pub(crate) fn traced_section(
    first: &Surface,
    second: &Surface,
) -> Result<Option<ExactIntersectionCurve>, ExactIntersectionRefusal> {
    if carrier_of(first).is_none() || carrier_of(second).is_none() {
        return Ok(None);
    }
    let rank = |s: &Surface| match s {
        Surface::Torus(_) => 2,
        Surface::Sphere(_) => 1,
        _ => 0,
    };
    let (carrier, other) = if rank(first) >= rank(second) {
        (first, second)
    } else {
        (second, first)
    };
    if rank(carrier) == 0 {
        return Ok(None);
    }
    let sections = implicit_surface_intersection(carrier, other, None)?;
    let spans = sections
        .iter()
        .map(|s| Some(Interval::new(0.0, s.curve.end())))
        .collect();
    Ok(Some(ExactIntersectionCurve {
        branches: sections.into_iter().map(Curve3::ImplicitSection).collect(),
        derivation: Derivation::ImplicitTrace,
        spans,
    }))
}
