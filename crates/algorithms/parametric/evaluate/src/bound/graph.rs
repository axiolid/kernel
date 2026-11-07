//! Chord bounds of quadratic-graph and angle-graph pcurves (#414).
//!
//! Both families are a closed-form root `w(t)` of an implicit `F(t, w) = 0`:
//! `F = a w^2 + b w + c` for a [`QuadraticGraph2`] (the point is `(t, w)`)
//! and `F = a cos w + b sin w - c` for an [`AngleGraph2`] (the point is
//! `(w, t)`), the coefficients being [`Trig2`]. Either way the curve's
//! second derivative is `(0, w'')` or `(w'', 0)`, so the chord of `[lo, hi]`
//! is within `(hi - lo)^2 / 8 sup |w''|` (the module's chord bound), and the
//! implicit function theorem gives `w'' = -(F_tt + 2 F_tw w' + F_ww w'^2) /
//! F_w`, `|w'| <= L = sup |F_t| / inf |F_w|`.
//!
//! The suprema are taken over the box of the span by an interval `W` of
//! heights, which must hold the graph. The ends give `[w(lo), w(hi)]`,
//! widened by `rho` to `W`. The computed ends are not exact, so a root is
//! certified at each: `F` must change sign, in outward-rounded intervals,
//! across `[w_end - d, w_end + d]`, which holds the true root (`F_w` has no
//! zero on the box, so `F` is monotone in `w` there). If `F_w` has no zero
//! on the box and `rho > d + L (hi - lo) / 2`, with `L` taken on that same
//! box, the graph cannot leave `W`: from either end it moves at most `L` per
//! unit of `t`, so it stays within `rho` of the nearer end up to the middle
//! of the span. The largest `d` of the span's two ends is added to the bound
//! (the chord joins the computed ends).
//!
//! An [`AngleGraph2`] is evaluated in `(-pi, pi]`: where its closed form
//! crosses that cut the evaluated curve jumps by a turn, though the
//! implicit solution does not. A piece whose box `W` reaches the cut (within
//! a few ulps of it) is refused, so the closed form is continuous across a
//! certified piece.
//! Over the box the root is then unique (two roots would put a zero of
//! `F_w` between them), so the implicit solution, which exists across the
//! span for the same reason, is the branch the closed form names. `rho`
//! starts at `d` and is grown until that holds, or the bound is refused
//! (`None`).
//!
//! Every enclosure is outward-rounded: each interval operation widens its
//! result by one ulp each way, and the trigonometric ranges by a few.

use axiolid_core::Scalar;
use axiolid_curve::{AngleGraph2, QuadraticGraph2, Trig2};

use super::interval::{cos_iv, sin_iv, Iv};

/// Growth steps of `rho` before the bound is refused.
const RHO_STEPS: usize = 64;

/// A trigonometric polynomial and its first two derivatives over a span.
#[derive(Debug, Clone, Copy)]
struct TrigEnclosure {
    value: Iv,
    first: Iv,
    second: Iv,
}

/// Enclosures of `t`, `t'` and `t''` of a [`Trig2`] over `[lo, hi]`.
///
/// Each of the four waves is enclosed on its own (its extremes at the ends
/// or at the lattice of its turning points), and the terms are summed.
fn trig_enclose(t: &Trig2, lo: Scalar, hi: Scalar) -> TrigEnclosure {
    let (c1, s1) = (cos_iv(lo, hi), sin_iv(lo, hi));
    let (c2, s2) = (cos_iv(2.0 * lo, 2.0 * hi), sin_iv(2.0 * lo, 2.0 * hi));
    let sum = |terms: [Iv; 4]| terms.into_iter().fold(Iv::point(0.0), Iv::add);
    TrigEnclosure {
        value: sum([
            c1.scale(t.cos),
            s1.scale(t.sin),
            c2.scale(t.cos2),
            s2.scale(t.sin2),
        ])
        .add(Iv::point(t.constant)),
        first: sum([
            s1.scale(-t.cos),
            c1.scale(t.sin),
            s2.scale(-2.0 * t.cos2),
            c2.scale(2.0 * t.sin2),
        ]),
        second: sum([
            c1.scale(-t.cos),
            s1.scale(-t.sin),
            c2.scale(-4.0 * t.cos2),
            s2.scale(-4.0 * t.sin2),
        ]),
    }
}

/// The enclosures of the three coefficients over the span.
struct Coefficients {
    a: TrigEnclosure,
    b: TrigEnclosure,
    c: TrigEnclosure,
}

/// The partials of `F` over a box.
struct Partials {
    t: Iv,
    w: Iv,
    tt: Iv,
    tw: Iv,
    ww: Iv,
}

/// A family `F(t, w) = 0` with a closed form `w(t)`.
trait Field {
    /// What the family encloses over a span, the input of `value` and
    /// `partials`.
    type Piece;

    /// Whether every number of the family is finite.
    fn is_finite(&self) -> bool;

    /// The enclosure of the family over `[lo, hi]`, or `None` where `F(t,
    /// .)` certainly has no root for some `t` of the span.
    fn piece(&self, lo: Scalar, hi: Scalar) -> Option<Self::Piece>;

    /// The closed-form `w(t)`.
    fn height(&self, t: Scalar) -> Option<Scalar>;

    /// Whether `heights` reaches a cut of the closed form, across which it
    /// jumps though the implicit solution does not.
    fn cut(&self, heights: Iv) -> bool;

    /// `F` at `t` in the span and `w` in `heights`.
    fn value(&self, k: &Self::Piece, heights: Iv) -> Iv;

    /// The partials of `F` over `t` in the span and `w` in `heights`.
    fn partials(&self, k: &Self::Piece, heights: Iv) -> Partials;
}

/// The enclosures of the three coefficients of a [`Trig2`] family over a
/// span, or `None` where `undefined` holds of them.
fn enclose_family(
    [a, b, c]: [&Trig2; 3],
    lo: Scalar,
    hi: Scalar,
    undefined: impl Fn(&Coefficients) -> bool,
) -> Option<Coefficients> {
    let [a, b, c] = [a, b, c].map(|t| trig_enclose(t, lo, hi));
    let k = Coefficients { a, b, c };
    (!undefined(&k)).then_some(k)
}

impl Field for QuadraticGraph2 {
    type Piece = Coefficients;

    fn is_finite(&self) -> bool {
        [&self.a, &self.b, &self.c].iter().all(|t| t.is_finite())
    }

    fn piece(&self, lo: Scalar, hi: Scalar) -> Option<Coefficients> {
        enclose_family([&self.a, &self.b, &self.c], lo, hi, |k| {
            // `b^2 - 4ac < 0`.
            let d = k.b.value.sq().add(k.a.value.mul(k.c.value).scale(-4.0));
            d.hi < 0.0
        })
    }

    fn height(&self, t: Scalar) -> Option<Scalar> {
        QuadraticGraph2::height(self, t)
    }

    fn cut(&self, _heights: Iv) -> bool {
        false
    }

    fn value(&self, k: &Coefficients, w: Iv) -> Iv {
        k.a.value.mul(w.sq()).add(k.b.value.mul(w)).add(k.c.value)
    }

    fn partials(&self, k: &Coefficients, w: Iv) -> Partials {
        // `F = a w^2 + b w + c`.
        let w2 = w.sq();
        let poly = |a: Iv, b: Iv, c: Iv| a.mul(w2).add(b.mul(w)).add(c);
        Partials {
            t: poly(k.a.first, k.b.first, k.c.first),
            w: k.a.value.mul(w).scale(2.0).add(k.b.value),
            tt: poly(k.a.second, k.b.second, k.c.second),
            tw: k.a.first.mul(w).scale(2.0).add(k.b.first),
            ww: k.a.value.scale(2.0),
        }
    }
}

impl Field for AngleGraph2 {
    type Piece = Coefficients;

    fn is_finite(&self) -> bool {
        [&self.a, &self.b, &self.c].iter().all(|t| t.is_finite())
    }

    fn piece(&self, lo: Scalar, hi: Scalar) -> Option<Coefficients> {
        enclose_family([&self.a, &self.b, &self.c], lo, hi, |k| {
            // `a^2 + b^2 - c^2 < 0`.
            let e = k.a.value.sq().add(k.b.value.sq()).add(k.c.value.sq().neg());
            e.hi < 0.0
        })
    }

    fn height(&self, t: Scalar) -> Option<Scalar> {
        self.angle(t)
    }

    fn cut(&self, heights: Iv) -> bool {
        // The angle is returned in `(-pi, pi]`, within a few ulps of it.
        let margin = 64.0 * Scalar::EPSILON * core::f64::consts::PI;
        heights.hi >= core::f64::consts::PI - margin
            || heights.lo <= -core::f64::consts::PI + margin
    }

    fn value(&self, k: &Coefficients, w: Iv) -> Iv {
        let (cw, sw) = (cos_iv(w.lo, w.hi), sin_iv(w.lo, w.hi));
        k.a.value
            .mul(cw)
            .add(k.b.value.mul(sw))
            .add(k.c.value.neg())
    }

    fn partials(&self, k: &Coefficients, w: Iv) -> Partials {
        // `F = a cos w + b sin w - c`.
        let (cw, sw) = (cos_iv(w.lo, w.hi), sin_iv(w.lo, w.hi));
        let mix = |a: Iv, b: Iv| a.mul(cw).add(b.mul(sw));
        let turn = |a: Iv, b: Iv| b.mul(cw).add(a.mul(sw).neg());
        Partials {
            t: mix(k.a.first, k.b.first).add(k.c.first.neg()),
            w: turn(k.a.value, k.b.value),
            tt: mix(k.a.second, k.b.second).add(k.c.second.neg()),
            tw: turn(k.a.first, k.b.first),
            ww: mix(k.a.value, k.b.value).neg(),
        }
    }
}

/// A certified bound on the distance from the quadratic graph over
/// `[lo, hi]` to the chord joining its ends, or `None` where the root does
/// not certainly exist across the span or `F_w` is not bounded away from
/// zero (a branch end nearby).
pub(crate) fn quadratic_chord_bound(
    graph: &QuadraticGraph2,
    lo: Scalar,
    hi: Scalar,
) -> Option<Scalar> {
    certify(graph, lo, hi)
}

/// As [`quadratic_chord_bound`], for an angle graph.
pub(crate) fn angle_chord_bound(graph: &AngleGraph2, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    certify(graph, lo, hi)
}

/// The finest split of a span into equal pieces that is tried.
const MAX_PIECES: usize = 64;

/// The chord bound over `[lo, hi]`, the tightest of the splits into `1, 2,
/// 4, ...` equal pieces that all certify.
///
/// `sup |w''|` over the span is the largest of the pieces', and each piece
/// is certified from its own ends (a piece's box is smaller, so the
/// interval enclosures, which ignore the correlation of `t` and `w`, are
/// tighter). A piece's end is the closed-form height there, which is the
/// next piece's start; the pieces chain into one continuous solution, whose
/// `F_w` has the branch's sign at every end, so it is the branch.
fn certify<F: Field>(field: &F, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return None;
    }
    if !field.is_finite() {
        return None;
    }
    let span = Iv::point(hi).add(Iv::point(-lo)).hi;
    let mut best: Option<Scalar> = None;
    let mut pieces = 1;
    while pieces <= MAX_PIECES {
        if let Some((bend, err)) = split_bend(field, lo, hi, pieces) {
            // The chord joins the computed ends, whose certified distance
            // from the roots the bound carries.
            let bound = Iv::point(span)
                .mul(Iv::point(span))
                .scale(0.125)
                .mul(Iv::point(bend))
                .add(Iv::point(err))
                .hi;
            if bound.is_finite() && best.is_none_or(|b| bound < b) {
                best = Some(bound);
            }
        }
        pieces *= 2;
    }
    best
}

/// The rounding of two computed heights.
fn rounding(a: Scalar, b: Scalar) -> Scalar {
    64.0 * Scalar::EPSILON * a.abs().max(b.abs()).max(1.0)
}

/// A certified piece: `sup |w''|`, and how far each computed end is from its
/// root.
#[derive(Debug, Clone, Copy)]
struct Enclosure {
    bend: Scalar,
    start: Scalar,
    end: Scalar,
}

/// A distance `d` such that the root of `F(t, .)` lies within `d` of `w`,
/// certified by a sign change of `F` across `[w - d, w + d]`, or `None`.
fn end_error<F: Field>(field: &F, t: Scalar, w: Scalar) -> Option<Scalar> {
    let k = field.piece(t, t)?;
    let mut d = rounding(w, w);
    for _ in 0..RHO_STEPS {
        let at = |x: Scalar| field.value(&k, Iv::point(x).widen());
        let (below, above) = (at(w - d), at(w + d));
        if (below.hi < 0.0 && above.lo > 0.0) || (below.lo > 0.0 && above.hi < 0.0) {
            return Some(d);
        }
        d *= 2.0;
    }
    None
}

/// `sup |w''|` over `[lo, hi]` split into `pieces`, with the larger distance
/// of the span's two computed ends from their roots, or `None` if a piece
/// does not certify.
fn split_bend<F: Field>(
    field: &F,
    lo: Scalar,
    hi: Scalar,
    pieces: usize,
) -> Option<(Scalar, Scalar)> {
    let at = |i: usize| {
        if i == pieces {
            hi
        } else {
            lo + (hi - lo) * (i as Scalar / pieces as Scalar)
        }
    };
    let mut bend: Scalar = 0.0;
    let (mut start, mut end) = (0.0, 0.0);
    for i in 0..pieces {
        let (a, b) = (at(i), at(i + 1));
        if a >= b {
            return None;
        }
        let piece = piece_enclosure(field, a, b)?;
        bend = bend.max(piece.bend);
        if i == 0 {
            start = piece.start;
        }
        end = piece.end;
    }
    Some((bend, start.max(end)))
}

/// `sup |w''|` over one piece, from the box of the piece by the heights
/// `rho` about the ends' (see the module docs).
fn piece_enclosure<F: Field>(field: &F, lo: Scalar, hi: Scalar) -> Option<Enclosure> {
    let (w_lo, w_hi) = (field.height(lo)?, field.height(hi)?);
    let k = field.piece(lo, hi)?;
    let span = Iv::point(hi).add(Iv::point(-lo)).hi;
    let (w_min, w_max) = (w_lo.min(w_hi), w_lo.max(w_hi));
    let (d_lo, d_hi) = (end_error(field, lo, w_lo)?, end_error(field, hi, w_hi)?);
    let d = d_lo.max(d_hi);
    let mut rho = d;
    for _ in 0..RHO_STEPS {
        let heights = Iv {
            lo: w_min - rho,
            hi: w_max + rho,
        }
        .widen();
        if field.cut(heights) {
            return None;
        }
        let p = field.partials(&k, heights);
        if ![p.t, p.w, p.tt, p.tw, p.ww].iter().all(|i| i.finite()) || p.w.straddles_zero() {
            return None;
        }
        let floor = p.w.mig();
        let lipschitz = (p.t.mag() / floor).next_up();
        let reach = Iv::point(lipschitz).mul(Iv::point(span)).scale(0.5).hi;
        if Iv::point(reach).add(Iv::point(d)).hi < rho {
            let l = Iv::point(lipschitz);
            let numerator = Iv::point(p.tt.mag())
                .add(Iv::point(p.tw.mag()).mul(l).scale(2.0))
                .add(Iv::point(p.ww.mag()).mul(l).mul(l));
            let bend = (numerator.hi / floor).next_up();
            return bend.is_finite().then_some(Enclosure {
                bend,
                start: d_lo,
                end: d_hi,
            });
        }
        rho = (1.1 * reach).max(1.1 * rho);
    }
    None
}

#[cfg(test)]
mod tests;
