//! Chord bounds of implicit curves (ADR 0077, #249).
//!
//! An [`ImplicitCurve2`] is a chain of cells. In a regular cell the free
//! parameter `f` runs linearly and the solved one `w(f)` is the field's
//! unique zero in the cell's bracket, where the field is strictly monotone
//! in it. By the implicit function theorem `|w'| <= L`, `L` the largest
//! `|F_free|` over the smallest `|F_solved|` on the box of the stretch's
//! free range by the bracket: interval bounds of the field's partials
//! (`axiolid_curve::implicit::bound_simple`, which cover their own
//! rounding). Between its ends `(f0, w0)` and `(f1, w1)` such a graph
//! leaves its chord, of slope `m`, by at most `(L + |m|) |f1 - f0| / 2`
//! along `w`, and its distance to the chord is no more. Halving the
//! stretch halves the bound.
//!
//! A bridge cell (the last stretch into a point where two branches cross)
//! is a cubic in its local parameter, so it lies in the convex hull of its
//! Bezier control points: the largest distance of those points from the
//! chord bounds it, and shrinks quadratically under halving.
//!
//! Cells meet with matching value and slope, but each is bounded on its
//! own, so a stretch over more than one cell is refused: the cell joins
//! are named as continuity breaks (`continuity_breaks2`).

use axiolid_core::{Point2, Scalar};
use axiolid_curve::implicit::{bound_simple, partial, Cell, Range};
use axiolid_curve::{Axis, ImplicitCurve2};

/// The cell joins of an implicit curve: `1, 2, ...` up to its last cell.
pub(crate) fn joins(curve: &ImplicitCurve2) -> Vec<Scalar> {
    (1..curve.cells.len()).map(|k| k as Scalar).collect()
}

/// A certified bound on the distance from the curve over `[lo, hi]` to the
/// chord joining its ends, or `None` over more than one cell, where the
/// field's partial along the solved parameter may vanish, or where a point
/// cannot be evaluated.
pub(crate) fn chord_bound(curve: &ImplicitCurve2, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    let count = curve.cells.len();
    if count == 0 || !(lo.is_finite() && hi.is_finite()) || lo < 0.0 || hi > count as Scalar {
        return None;
    }
    let index = (lo.floor() as usize).min(count - 1);
    if hi > (index + 1) as Scalar {
        return None;
    }
    let cell = curve.cells[index];
    let local = |t: Scalar| (t - index as Scalar).clamp(0.0, 1.0);
    let part = cell.part(local(lo), local(hi));
    let place = |free: Scalar, solved: Scalar| match cell.axis {
        Axis::U => Point2::new(free, solved),
        Axis::V => Point2::new(solved, free),
    };
    let span = part.to - part.from;
    if let Some((m0, m1)) = part.bridge {
        // The cubic's Bezier control points, in `(free, solved)`.
        let controls = [
            place(part.from, part.low),
            place(part.from + span / 3.0, part.low + m0 * span / 3.0),
            place(part.to - span / 3.0, part.high - m1 * span / 3.0),
            place(part.to, part.high),
        ];
        return Some(
            controls
                .iter()
                .map(|&p| distance_to_line(p, controls[0], controls[3]))
                .fold(0.0, Scalar::max),
        );
    }
    let (start, end) = (curve.point(lo)?, curve.point(hi)?);
    let solved = |p: Point2| match cell.axis {
        Axis::U => p.y,
        Axis::V => p.x,
    };
    if span == 0.0 {
        return Some((solved(end) - solved(start)).abs());
    }
    let along_u = cell.axis == Axis::U;
    let (a, b) = (
        place(part.from.min(part.to), cell.low.min(cell.high)),
        place(part.from.max(part.to), cell.low.max(cell.high)),
    );
    let corner = Cell {
        lo: a.min(b),
        hi: a.max(b),
    };
    let d_free = partial(&curve.field, along_u);
    let d_solved = partial(&curve.field, !along_u);
    let solved_rate = bound_simple(&d_solved, &corner);
    if solved_rate.straddles_zero() {
        return None;
    }
    let most = |r: Range| r.lo.abs().max(r.hi.abs());
    let floor = solved_rate.lo.abs().min(solved_rate.hi.abs());
    let lipschitz = most(bound_simple(&d_free, &corner)) / floor;
    let slope = ((solved(end) - solved(start)) / span).abs();
    let first = 0.5 * (lipschitz + slope) * span.abs();
    // `w'' = -(F_ff + 2 F_fw w' + F_ww w'^2) / F_w`, so the interpolant of
    // the ends is within `span^2 / 8 sup |w''|` of the graph.
    let ff = most(bound_simple(&partial(&d_free, along_u), &corner));
    let fw = most(bound_simple(&partial(&d_free, !along_u), &corner));
    let ww = most(bound_simple(&partial(&d_solved, !along_u), &corner));
    let bend = (ff + 2.0 * fw * lipschitz + ww * lipschitz * lipschitz) / floor;
    let second = span * span / 8.0 * bend;
    let bound = first.min(second);
    bound.is_finite().then_some(bound)
}

/// The distance from `p` to the line through `a` and `b` (to `a` when they
/// coincide).
fn distance_to_line(p: Point2, a: Point2, b: Point2) -> Scalar {
    let d = b - a;
    let length = d.length();
    if length == 0.0 {
        return (p - a).length();
    }
    (d.perp_dot(p - a) / length).abs()
}
