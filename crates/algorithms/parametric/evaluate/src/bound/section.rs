//! Derivative bounds of the space curves a lifted pcurve reads (#414).
//!
//! The families with a derivative bound in [`curve_derivative_bounds3`]
//! are covered by it. The sections of an analytic carrier are a plane
//! curve `q(t)` read through the surface: `c = S(q)`, so
//! `|c'| <= S_u |q_u'| + S_v |q_v'|` and
//! `|c''| <= S_uu q_u'^2 + 2 S_uv |q_u' q_v'| + S_vv q_v'^2 + S_u |q_u''| +
//! S_v |q_v''|`, each with the suprema of the surface's partials over the
//! box of `q` and those of `q`'s own derivatives over the span:
//!
//! - a [`RuledSection3`] is `q = (t, w)` and a [`TorusSection3`] `q = (w,
//!   t)`, `w` the graph's root (`graph`);
//! - an [`ImplicitSection3`] is one cell of its implicit curve, whose rate
//!   `r` along its free parameter is constant (`implicit`).
//!
//! A second derivative needs the carrier to be smooth over the box; on a
//! B-spline carrier crossing a crease it is unknown (`second: None`), and
//! the first-order bound stands alone. The box of a spline carrier is
//! clamped to the knot domain, which holds the section by construction.
//! Every other family has no bound here.

use axiolid_core::Scalar;
use axiolid_curve::{Carrier, Curve3};

use super::graph::{angle_bounds, quadratic_bounds, PlaneBounds};
use super::interval::Iv;
use super::{curve_derivative_bounds3, implicit, SurfaceBoundOracle};

/// Suprema of the first two derivatives of a space curve over a span.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct SpaceBounds {
    /// `sup |c'|`.
    pub first: Scalar,
    /// `sup |c''|`, unknown where the carrier is not smooth.
    pub second: Option<Scalar>,
}

/// The bounds of `curve` over `[lo, hi]`, or `None` for a family without
/// one or where a bound cannot be certified.
pub(super) fn space_bounds(curve: &Curve3, lo: Scalar, hi: Scalar) -> Option<SpaceBounds> {
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return None;
    }
    match curve {
        Curve3::RuledSection(s) => {
            let carrier = Carrier::Ruled(s.carrier);
            through(&carrier, &quadratic_bounds(&s.graph, lo, hi)?)
        }
        Curve3::TorusSection(s) => {
            let carrier = Carrier::Torus(s.torus);
            through(&carrier, &angle_bounds(&s.graph, lo, hi)?)
        }
        Curve3::ImplicitSection(s) => {
            through(&s.carrier, &implicit::derivative_bounds(&s.curve, lo, hi)?)
        }
        _ => curve_derivative_bounds3(curve, lo, hi).map(|d| SpaceBounds {
            first: d.first,
            second: Some(d.second),
        }),
    }
}

/// `sum a_i b_i`, rounded up.
fn up(terms: &[(Scalar, Scalar)]) -> Scalar {
    terms
        .iter()
        .map(|&(a, b)| Iv::point(a).mul(Iv::point(b)))
        .fold(Iv::point(0.0), Iv::add)
        .hi
}

/// The bounds of `S(q(t))` from those of the plane curve `q` over the box
/// of `bounds`.
fn through(carrier: &Carrier, q: &PlaneBounds) -> Option<SpaceBounds> {
    let oracle = SurfaceBoundOracle::for_carrier(carrier)?;
    let clamp = |lo: Scalar, hi: Scalar, range: Option<(Scalar, Scalar)>| match range {
        Some((a, b)) => (lo.clamp(a, b), hi.clamp(a, b)),
        None => (lo, hi),
    };
    let domain = oracle.domain();
    let u = clamp(q.lo.x, q.hi.x, domain.map(|d| d.0));
    let v = clamp(q.lo.y, q.hi.y, domain.map(|d| d.1));
    let b = oracle.bounds(u, v)?;
    let [fu, fv] = q.first;
    let [su, sv] = q.second;
    let first = up(&[(b.du, fu), (b.dv, fv)]);
    let second = b.smooth.then(|| {
        up(&[
            (b.duu, up(&[(fu, fu)])),
            (2.0 * b.duv, up(&[(fu, fv)])),
            (b.dvv, up(&[(fv, fv)])),
            (b.du, su),
            (b.dv, sv),
        ])
    });
    (first.is_finite() && second.is_none_or(|s| s.is_finite()))
        .then_some(SpaceBounds { first, second })
}

#[cfg(test)]
mod tests;
