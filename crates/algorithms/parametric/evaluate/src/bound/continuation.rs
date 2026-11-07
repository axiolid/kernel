//! Chord bounds from piecewise certified derivative bounds (#414).
//!
//! A lifted pcurve, or one side of a pair section, is the solution `p(t)`
//! of an implicit equation. Over a piece of the span the implicit function
//! theorem, with a certified separation of the surface's parameters over a
//! box holding the piece (`lifted`, `pair`), bounds `sup |p'|` and `sup
//! |p''|`; the pieces' bounds are combined here, as `graph` does for its
//! graphs.
//!
//! For `p` with `|p'| <= L` the chord of its ends on `[a, b]`, `h = b - a`,
//! is within `L h / 2` of it: `p(t)` is within `L (t - a)` of `p(a)` and
//! `L (b - t)` of `p(b)`, and the chord's point is their convex mix, giving
//! `2 L h s (1 - s) <= L h / 2` at `t = a + s h`. With `|p''| <= B` it is
//! within `B h^2 / 8`; the smaller of the two stands, and where `B` is
//! unknown the first alone.
//!
//! The chord joins the *evaluated* ends, not the exact ones: each end is
//! off by a certified distance, and the bound carries the larger of the
//! two (the chord of the exact ends and that of the evaluated ends differ
//! by at most that, everywhere).
//!
//! The supremum over the span is the largest over the pieces of an equal
//! split, each certified from its own box, so a finer split is tighter
//! where the box shrinks faster than the second derivative grows. The
//! tightest of the splits into 1, 2, 4, ... pieces that all certify stands.

use axiolid_core::Scalar;

use super::interval::Iv;

/// The finest split of a span into equal pieces that is tried.
const MAX_PIECES: usize = 64;

/// A certified piece of a pcurve, in the parameters' own units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct PieceBound {
    /// `sup |p'|` over the piece.
    pub slope: Scalar,
    /// `sup |p''|`, where the curve's own second derivative is bounded.
    pub bend: Option<Scalar>,
    /// How far the evaluated start is from the exact one.
    pub start: Scalar,
    /// How far the evaluated end is from the exact one.
    pub end: Scalar,
}

/// The bound on the distance of the curve over `[lo, hi]` from the chord of
/// its evaluated ends, from the pieces `piece` certifies, or `None` if no
/// split certifies.
pub(super) fn split_chord(
    lo: Scalar,
    hi: Scalar,
    mut piece: impl FnMut(Scalar, Scalar) -> Option<PieceBound>,
) -> Option<Scalar> {
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return None;
    }
    let span = Iv::point(hi).add(Iv::point(-lo)).hi;
    let mut best: Option<Scalar> = None;
    let mut pieces = 1;
    while pieces <= MAX_PIECES {
        if let Some(bound) = split(lo, hi, span, pieces, &mut piece) {
            if bound.is_finite() && best.is_none_or(|b| bound < b) {
                best = Some(bound);
            }
        }
        pieces *= 2;
    }
    best
}

/// The bound of one split into `pieces` equal pieces.
fn split(
    lo: Scalar,
    hi: Scalar,
    span: Scalar,
    pieces: usize,
    piece: &mut impl FnMut(Scalar, Scalar) -> Option<PieceBound>,
) -> Option<Scalar> {
    let at = |i: usize| {
        if i == pieces {
            hi
        } else {
            lo + (hi - lo) * (i as Scalar / pieces as Scalar)
        }
    };
    let mut slope: Scalar = 0.0;
    let mut bend: Option<Scalar> = Some(0.0);
    let (mut start, mut end) = (0.0, 0.0);
    for i in 0..pieces {
        let (a, b) = (at(i), at(i + 1));
        if a >= b {
            return None;
        }
        let p = piece(a, b)?;
        slope = slope.max(p.slope);
        bend = bend.zip(p.bend).map(|(x, y)| x.max(y));
        if i == 0 {
            start = p.start;
        }
        end = p.end;
    }
    let first = Iv::point(span).mul(Iv::point(slope)).scale(0.5).hi;
    let second = bend.map(|bend| {
        Iv::point(span)
            .mul(Iv::point(span))
            .scale(0.125)
            .mul(Iv::point(bend))
            .hi
    });
    let within = second.map_or(first, |s| s.min(first));
    Some(Iv::point(within).add(Iv::point(start.max(end))).hi)
}
