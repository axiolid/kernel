//! Chord bounds of lifted pcurves on analytic carriers (#414).
//!
//! A [`LiftedCurve2`] reads a space curve `c(t)` in the parameters `p` of an
//! analytic surface `S`: `S(p(t)) = c(t)`. Per piece `[a, b]` of the span
//! (`h = b - a`) take the evaluated ends `p_a`, `p_b` and the box `P` of
//! their hull widened by `rho`, and the separation `sigma` of `S` over `P`
//! in parameters scaled by `D = diag(|S_u|, |S_v|)` at the box's centre
//! (`SurfaceBoundOracle::separation`): `|J w| >= sigma |D w|` at every point
//! of `P`, and `|S(x) - S(y)| >= sigma |D (x - y)|` on it.
//!
//! - **Speed.** `J p' = c'` gives `|D p'| <= C1 / sigma` for `C1 >= sup
//!   |c'|` (`section`), so `|p'| <= C1 / (sigma m)`, `m = min(|S_u|,
//!   |S_v|)`.
//! - **Reach.** The lift of an exact end moves at most `C1 h / (2 sigma m)`
//!   up to the middle of the piece; it starts within `e` of the evaluated
//!   end. If that is less than `rho` the lift stays in `P` up to the middle
//!   from either end, and injectivity on `P` makes the two lifts one: the
//!   evaluated ends lie on one continuous branch. `rho` grows from zero
//!   until it does, or the piece is refused. A larger box only has a
//!   smaller separation, so a box without one ends the search.
//! - **Ends.** The evaluated end `p_i` is off its exact lift by at most
//!   `|S(p_i) - c_i| / (sigma m)`, `c_i` the space curve's point as
//!   evaluated: the residual is enclosed in outward-rounded arithmetic
//!   (`carrier::point_enclosure`), and injectivity turns it into a
//!   parameter distance. The space curve is read as the evaluator computes
//!   it, and lies on the carrier.
//! - **Bend.** Differentiating twice, `J p'' = c'' - H[p', p']`, so `|D p''|
//!   <= (C2 + |D p'|^2 (S_uu / |S_u|^2 + 2 S_uv / (|S_u| |S_v|) + S_vv /
//!   |S_v|^2)) / sigma`, `C2 >= sup |c''|`. Where `C2` is unknown only the
//!   speed bounds the chord.
//!
//! An angle is evaluated at the turn nearest the curve's guide. The
//! evaluated pcurve is the continuous lift only if that turn does not
//! change inside the piece: the guide over the piece must stay within half
//! a turn of the box. Otherwise the piece is refused.

use axiolid_core::{Point2, Scalar};
use axiolid_curve::{Carrier, LiftedCurve2};

use super::carrier::point_enclosure;
use super::continuation::{split_chord, PieceBound};
use super::interval::Iv;
use super::section::space_bounds;
use super::SurfaceBoundOracle;
use crate::curve::{evaluate3, pair_side};

/// Growth steps of `rho` before the bound is refused.
const RHO_STEPS: usize = 64;

/// How far inside half a turn the guide must stay of the box, so that the
/// nearest turn is never a tie.
const TURN_MARGIN: Scalar = 1e-9;

/// A certified bound on the distance from the lifted pcurve over `[lo, hi]`
/// to the chord joining its evaluated ends, or `None` where a piece cannot
/// be certified (a B-spline carrier outside a pair section, a pole, a
/// guide that changes the turn, a space curve without derivative bounds).
pub(super) fn chord_bound(l: &LiftedCurve2, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    // A section of two surfaces read on one of them is evaluated by the
    // section's own solve, not by inverting the carrier.
    if !l.is_finite() || matches!(l.carrier, Carrier::Spline(_)) || pair_side(l).is_some() {
        return None;
    }
    let oracle = SurfaceBoundOracle::for_carrier(&l.carrier)?;
    split_chord(lo, hi, |a, b| piece(l, &oracle, a, b))
}

/// An evaluated end, and the distance of its image from the space curve's
/// point.
struct End {
    at: Point2,
    residual: Scalar,
}

fn end(l: &LiftedCurve2, t: Scalar) -> Option<End> {
    let space = evaluate3(&l.curve, t).ok()?;
    let at = l.unwrap_at(t, space)?;
    let image = point_enclosure(&l.carrier, at.x, at.y)?;
    let c = [space.x, space.y, space.z];
    let squares = image
        .iter()
        .zip(c)
        .map(|(s, c)| s.add(Iv::point(-c)).sq())
        .fold(Iv::point(0.0), Iv::add);
    let residual = squares.sqrt().hi;
    (at.is_finite() && residual.is_finite()).then_some(End { at, residual })
}

/// The range of the guide over `[a, b]`: its values at the ends and at the
/// nodes between.
fn guide_range(l: &LiftedCurve2, a: Scalar, b: Scalar) -> Option<(Point2, Point2)> {
    let n = l.guide.len();
    let mut lo = l.guide_at(a)?.min(l.guide_at(b)?);
    let mut hi = l.guide_at(a)?.max(l.guide_at(b)?);
    if n > 1 && l.end != l.start {
        for (i, g) in l.guide.iter().enumerate() {
            let t = l.start + (l.end - l.start) * (i as Scalar / (n - 1) as Scalar);
            if t > a && t < b {
                lo = lo.min(*g);
                hi = hi.max(*g);
            }
        }
    }
    Some((lo, hi))
}

/// Whether the guide stays within half a turn of the box on every periodic
/// axis, so the turn read at each `t` is the box's.
fn guide_holds(l: &LiftedCurve2, a: Scalar, b: Scalar, boxed: [(Scalar, Scalar); 2]) -> bool {
    let (periodic_u, periodic_v) = l.carrier.periodic();
    let Some((g_lo, g_hi)) = guide_range(l, a, b) else {
        return false;
    };
    let half = Iv::point(core::f64::consts::PI)
        .add(Iv::point(-TURN_MARGIN))
        .lo;
    [
        (periodic_u, g_lo.x, g_hi.x, boxed[0]),
        (periodic_v, g_lo.y, g_hi.y, boxed[1]),
    ]
    .into_iter()
    .all(|(periodic, g_lo, g_hi, (x_lo, x_hi))| {
        !periodic
            || (Iv::point(g_hi).add(Iv::point(-x_lo)).hi < half
                && Iv::point(x_hi).add(Iv::point(-g_lo)).hi < half)
    })
}

/// One piece `[a, b]`; see the module docs.
fn piece(
    l: &LiftedCurve2,
    oracle: &SurfaceBoundOracle<'_>,
    a: Scalar,
    b: Scalar,
) -> Option<PieceBound> {
    let (start, finish) = (end(l, a)?, end(l, b)?);
    let space = space_bounds(&l.curve, a, b)?;
    let h = Iv::point(b).add(Iv::point(-a)).hi;
    let hull = |x: Scalar, y: Scalar| (x.min(y), x.max(y));
    let (u, v) = (hull(start.at.x, finish.at.x), hull(start.at.y, finish.at.y));
    let residual = start.residual.max(finish.residual);
    let mut rho: Scalar = 0.0;
    for _ in 0..RHO_STEPS {
        let grow = |(lo, hi): (Scalar, Scalar)| {
            (
                Iv::point(lo).add(Iv::point(-rho)).lo,
                Iv::point(hi).add(Iv::point(rho)).hi,
            )
        };
        let boxed = [grow(u), grow(v)];
        if !guide_holds(l, a, b, boxed) {
            return None;
        }
        let sep = oracle.separation(boxed[0], boxed[1])?;
        let m = sep.scale[0].min(sep.scale[1]);
        let sigma = Iv::point(sep.sigma);
        let floor = sigma.mul(Iv::point(m));
        // `|D p'| <= C1 / sigma`, `|p'| <= C1 / (sigma m)`.
        let speed = Iv::point(space.first).div(sigma).hi;
        let slope = Iv::point(space.first).div(floor).hi;
        let error = Iv::point(residual).div(floor).hi;
        let reach = Iv::point(slope)
            .mul(Iv::point(h))
            .scale(0.5)
            .add(Iv::point(error))
            .hi;
        if reach < rho {
            let bend = match space.second {
                Some(second) => Some(bend(oracle, boxed, sep.scale, speed, second, sigma, m)?),
                None => None,
            };
            return Some(PieceBound {
                slope,
                bend,
                start: Iv::point(start.residual).div(floor).hi,
                end: Iv::point(finish.residual).div(floor).hi,
            });
        }
        rho = (1.1 * reach).max(1.1 * rho).max(Scalar::MIN_POSITIVE);
    }
    None
}

/// `sup |p''|` from `|D p''| <= (C2 + |D p'|^2 coefficient) / sigma`.
fn bend(
    oracle: &SurfaceBoundOracle<'_>,
    boxed: [(Scalar, Scalar); 2],
    scale: [Scalar; 2],
    speed: Scalar,
    second: Scalar,
    sigma: Iv,
    m: Scalar,
) -> Option<Scalar> {
    let b = oracle.bounds(boxed[0], boxed[1])?;
    let (al, be) = (Iv::point(scale[0]), Iv::point(scale[1]));
    let coefficient = Iv::point(b.duu)
        .div(al.mul(al))
        .add(Iv::point(2.0 * b.duv).div(al.mul(be)))
        .add(Iv::point(b.dvv).div(be.mul(be)));
    let bend = Iv::point(speed)
        .mul(Iv::point(speed))
        .mul(coefficient)
        .add(Iv::point(second))
        .div(sigma)
        .div(Iv::point(m))
        .hi;
    bend.is_finite().then_some(bend)
}
