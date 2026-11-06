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
//! widened by `rho` to `W`. If `F_w` has no zero on the box and `rho >
//! L (hi - lo) / 2`, with `L` taken on that same box, the graph cannot
//! leave `W`: from either end it moves at most `L` per unit of `t`, so
//! it stays within `rho` of the nearer end up to the middle of the span.
//! Over the box the root is then unique (two roots would put a zero of
//! `F_w` between them), so the implicit solution, which exists across the
//! span for the same reason, is the branch the closed form names. `rho`
//! starts at the ends' rounding and is grown until that holds, or the
//! bound is refused (`None`).
//!
//! Every enclosure is outward-rounded: each interval operation widens its
//! result by one ulp each way, and the trigonometric ranges by a few.

use axiolid_core::Scalar;
use axiolid_curve::{AngleGraph2, QuadraticGraph2, Trig2};

use super::{contains_lattice, cos_range};

/// Absolute padding of a sine or cosine range, covering the library's
/// rounding of the endpoint values.
const TRIG_PAD: Scalar = 8.0 * Scalar::EPSILON;

/// Growth steps of `rho` before the bound is refused.
const RHO_STEPS: usize = 64;

/// A closed interval whose operations round outward.
#[derive(Debug, Clone, Copy)]
struct Iv {
    lo: Scalar,
    hi: Scalar,
}

impl Iv {
    fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    fn widen(self) -> Self {
        Self {
            lo: self.lo.next_down(),
            hi: self.hi.next_up(),
        }
    }

    fn add(self, o: Self) -> Self {
        Self {
            lo: self.lo + o.lo,
            hi: self.hi + o.hi,
        }
        .widen()
    }

    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self {
            lo: p.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            hi: p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        }
        .widen()
    }

    fn scale(self, k: Scalar) -> Self {
        self.mul(Self::point(k))
    }

    fn sq(self) -> Self {
        let (a, b) = (self.lo * self.lo, self.hi * self.hi);
        let lo = if self.lo <= 0.0 && self.hi >= 0.0 {
            0.0
        } else {
            a.min(b).next_down().max(0.0)
        };
        Self {
            lo,
            hi: a.max(b).next_up(),
        }
    }

    fn straddles_zero(self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }

    /// `sup |x|`.
    fn mag(self) -> Scalar {
        self.lo.abs().max(self.hi.abs())
    }

    /// `inf |x|`, for an interval holding no zero.
    fn mig(self) -> Scalar {
        self.lo.abs().min(self.hi.abs())
    }

    fn finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }
}

/// The range of `cos` over `[lo, hi]`.
fn cos_iv(lo: Scalar, hi: Scalar) -> Iv {
    let (min, max) = cos_range(lo, hi);
    Iv {
        lo: min - TRIG_PAD,
        hi: max + TRIG_PAD,
    }
}

/// The range of `sin` over `[lo, hi]`.
fn sin_iv(lo: Scalar, hi: Scalar) -> Iv {
    use core::f64::consts::{FRAC_PI_2, TAU};
    let max = if contains_lattice(lo, hi, FRAC_PI_2, TAU) {
        1.0
    } else {
        lo.sin().max(hi.sin())
    };
    let min = if contains_lattice(lo, hi, -FRAC_PI_2, TAU) {
        -1.0
    } else {
        lo.sin().min(hi.sin())
    };
    Iv {
        lo: min - TRIG_PAD,
        hi: max + TRIG_PAD,
    }
}

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
    fn coefficients(&self) -> [&Trig2; 3];

    /// The closed-form `w(t)`.
    fn height(&self, t: Scalar) -> Option<Scalar>;

    /// Whether `F(t, .)` certainly has no root for some `t` of the span.
    fn undefined(&self, k: &Coefficients) -> bool;

    /// The partials of `F` over `t` in the span and `w` in `heights`.
    fn partials(&self, k: &Coefficients, heights: Iv) -> Partials;
}

impl Field for QuadraticGraph2 {
    fn coefficients(&self) -> [&Trig2; 3] {
        [&self.a, &self.b, &self.c]
    }

    fn height(&self, t: Scalar) -> Option<Scalar> {
        QuadraticGraph2::height(self, t)
    }

    fn undefined(&self, k: &Coefficients) -> bool {
        // `b^2 - 4ac < 0`.
        let d = k.b.value.sq().add(k.a.value.mul(k.c.value).scale(-4.0));
        d.hi < 0.0
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
    fn coefficients(&self) -> [&Trig2; 3] {
        [&self.a, &self.b, &self.c]
    }

    fn height(&self, t: Scalar) -> Option<Scalar> {
        self.angle(t)
    }

    fn undefined(&self, k: &Coefficients) -> bool {
        // `a^2 + b^2 - c^2 < 0`.
        let e = k.a.value.sq().add(k.b.value.sq()).add(k.c.value.sq().neg());
        e.hi < 0.0
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
    if field.coefficients().iter().any(|t| !t.is_finite()) {
        return None;
    }
    let span = Iv::point(hi).add(Iv::point(-lo)).hi;
    // The chord joins the computed ends, whose rounding the bound carries.
    let err = rounding(field.height(lo)?, field.height(hi)?);
    let mut best: Option<Scalar> = None;
    let mut pieces = 1;
    while pieces <= MAX_PIECES {
        if let Some(bend) = split_bend(field, lo, hi, pieces) {
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

/// `sup |w''|` over `[lo, hi]` split into `pieces`, or `None` if a piece
/// does not certify.
fn split_bend<F: Field>(field: &F, lo: Scalar, hi: Scalar, pieces: usize) -> Option<Scalar> {
    let at = |i: usize| {
        if i == pieces {
            hi
        } else {
            lo + (hi - lo) * (i as Scalar / pieces as Scalar)
        }
    };
    let mut bend: Scalar = 0.0;
    for i in 0..pieces {
        let (a, b) = (at(i), at(i + 1));
        if a >= b {
            return None;
        }
        bend = bend.max(piece_bend(field, a, b)?);
    }
    Some(bend)
}

/// `sup |w''|` over one piece, from the box of the piece by the heights
/// `rho` about the ends' (see the module docs).
fn piece_bend<F: Field>(field: &F, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    let (w_lo, w_hi) = (field.height(lo)?, field.height(hi)?);
    let [a, b, c] = field.coefficients().map(|t| trig_enclose(t, lo, hi));
    let k = Coefficients { a, b, c };
    if field.undefined(&k) {
        return None;
    }
    let span = Iv::point(hi).add(Iv::point(-lo)).hi;
    let (w_min, w_max) = (w_lo.min(w_hi), w_lo.max(w_hi));
    let mut rho = rounding(w_lo, w_hi);
    for _ in 0..RHO_STEPS {
        let heights = Iv {
            lo: w_min - rho,
            hi: w_max + rho,
        }
        .widen();
        let p = field.partials(&k, heights);
        if ![p.t, p.w, p.tt, p.tw, p.ww].iter().all(|i| i.finite()) || p.w.straddles_zero() {
            return None;
        }
        let floor = p.w.mig();
        let lipschitz = (p.t.mag() / floor).next_up();
        let reach = Iv::point(lipschitz).mul(Iv::point(span)).scale(0.5).hi;
        if reach < rho {
            let l = Iv::point(lipschitz);
            let numerator = Iv::point(p.tt.mag())
                .add(Iv::point(p.tw.mag()).mul(l).scale(2.0))
                .add(Iv::point(p.ww.mag()).mul(l).mul(l));
            return Some((numerator.hi / floor).next_up()).filter(|b| b.is_finite());
        }
        rho = (1.1 * reach).max(1.1 * rho);
    }
    None
}

#[cfg(test)]
mod tests;
