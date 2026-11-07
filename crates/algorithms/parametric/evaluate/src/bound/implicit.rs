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

use super::graph::PlaneBounds;
use axiolid_core::{Point2, Scalar};
use axiolid_curve::implicit::{bound_simple, partial, Cell, Range};
use axiolid_curve::{Axis, Field2, ImplicitCell, ImplicitCurve2};

/// The cell joins of an implicit curve: `1, 2, ...` up to its last cell.
pub(crate) fn joins(curve: &ImplicitCurve2) -> Vec<Scalar> {
    (1..curve.cells.len()).map(|k| k as Scalar).collect()
}

/// A certified bound on the distance from the curve over `[lo, hi]` to the
/// chord joining its ends, or `None` over more than one cell, where the
/// field's partial along the solved parameter may vanish, or where a point
/// cannot be evaluated.
pub(crate) fn chord_bound(curve: &ImplicitCurve2, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    let (cell, part, _) = stretch(curve, lo, hi)?;
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
    let (_, lipschitz, bend) = regular_rates(&curve.field, &cell, &corner(&cell, &part))?;
    let slope = ((solved(end) - solved(start)) / span).abs();
    let first = 0.5 * (lipschitz + slope) * span.abs();
    // The interpolant of the ends is within `span^2 / 8 sup |w''|` of the
    // graph.
    let second = span * span / 8.0 * bend;
    let bound = first.min(second);
    bound.is_finite().then_some(bound)
}

/// The box of an implicit curve over `[lo, hi]` within one cell, and the
/// suprema of the derivatives of its two coordinates in `t`, or `None`
/// over more than one cell or where the field's partial along the solved
/// parameter may vanish.
///
/// A regular cell moves its free parameter at the constant rate `r`, so
/// `q' = r (1, w')` and `q'' = r^2 (0, w'')` in `(free, solved)`, with the
/// bounds of [`regular_rates`]. A bridge's solved parameter is a cubic
/// whose Bezier control values `c` bound `|w_s| <= 3 max |c_{i+1} - c_i|`
/// and `|w_ss| <= 6 max |c_{i+2} - 2 c_{i+1} + c_i|` in its local
/// parameter, which spans `dt` of the curve's.
pub(super) fn derivative_bounds(
    curve: &ImplicitCurve2,
    lo: Scalar,
    hi: Scalar,
) -> Option<PlaneBounds> {
    let (cell, part, dt) = stretch(curve, lo, hi)?;
    if dt.is_nan() || dt <= 0.0 {
        return None;
    }
    let rate = (part.to - part.from).abs() / dt;
    let place = |free: Scalar, solved: Scalar| match cell.axis {
        Axis::U => [free, solved],
        Axis::V => [solved, free],
    };
    let (solved_first, solved_second, solved_box) = if let Some((m0, m1)) = part.bridge {
        let third = (part.to - part.from) / 3.0;
        let c = [
            part.low,
            part.low + m0 * third,
            part.high - m1 * third,
            part.high,
        ];
        let d1 = (0..3)
            .map(|i| (c[i + 1] - c[i]).abs())
            .fold(0.0, Scalar::max);
        let d2 = (0..2)
            .map(|i| (c[i + 2] - 2.0 * c[i + 1] + c[i]).abs())
            .fold(0.0, Scalar::max);
        (3.0 * d1 / dt, 6.0 * d2 / (dt * dt), part.solved_range())
    } else {
        let b = corner(&cell, &part);
        let (_, lipschitz, bend) = regular_rates(&curve.field, &cell, &b)?;
        let range = match cell.axis {
            Axis::U => (b.lo.y, b.hi.y),
            Axis::V => (b.lo.x, b.hi.x),
        };
        (rate * lipschitz, rate * rate * bend, range)
    };
    let first = place(rate, solved_first);
    let second = place(0.0, solved_second);
    let free = (part.from.min(part.to), part.from.max(part.to));
    let (lo_box, hi_box) = (place(free.0, solved_box.0), place(free.1, solved_box.1));
    let finite = first.iter().chain(&second).all(|x| x.is_finite());
    finite.then_some(PlaneBounds {
        lo: Point2::new(lo_box[0], lo_box[1]),
        hi: Point2::new(hi_box[0], hi_box[1]),
        first: first.map(super::inflate),
        second: second.map(super::inflate),
    })
}

/// The cell holding `[lo, hi]`, the part of it over that range, and the
/// range's length in the cell's local parameter, or `None` over more than
/// one cell.
fn stretch(
    curve: &ImplicitCurve2,
    lo: Scalar,
    hi: Scalar,
) -> Option<(ImplicitCell, ImplicitCell, Scalar)> {
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
    Some((cell, part, local(hi) - local(lo)))
}

/// The box of a regular stretch: its free range by the cell's bracket.
fn corner(cell: &ImplicitCell, part: &ImplicitCell) -> Cell {
    let place = |free: Scalar, solved: Scalar| match cell.axis {
        Axis::U => Point2::new(free, solved),
        Axis::V => Point2::new(solved, free),
    };
    let (a, b) = (
        place(part.from.min(part.to), cell.low.min(cell.high)),
        place(part.from.max(part.to), cell.low.max(cell.high)),
    );
    Cell {
        lo: a.min(b),
        hi: a.max(b),
    }
}

/// `sup |w'|` and `sup |w''|` of the solved parameter `w` of a regular
/// stretch against the free one over `corner`, by the implicit function
/// theorem: `|w'| <= L = sup |F_free| / inf |F_w|` and `w'' = -(F_ff + 2
/// F_fw w' + F_ww w'^2) / F_w`; `None` where `F_w` may vanish. The first
/// number is the floor of `|F_w|`.
fn regular_rates(
    field: &Field2,
    cell: &ImplicitCell,
    corner: &Cell,
) -> Option<(Scalar, Scalar, Scalar)> {
    let along_u = cell.axis == Axis::U;
    let d_free = partial(field, along_u);
    let d_solved = partial(field, !along_u);
    let solved_rate = bound_simple(&d_solved, corner);
    if solved_rate.straddles_zero() {
        return None;
    }
    let most = |r: Range| r.lo.abs().max(r.hi.abs());
    let floor = solved_rate.lo.abs().min(solved_rate.hi.abs());
    let lipschitz = most(bound_simple(&d_free, corner)) / floor;
    let ff = most(bound_simple(&partial(&d_free, along_u), corner));
    let fw = most(bound_simple(&partial(&d_free, !along_u), corner));
    let ww = most(bound_simple(&partial(&d_solved, !along_u), corner));
    let bend = (ff + 2.0 * fw * lipschitz + ww * lipschitz * lipschitz) / floor;
    Some((floor, lipschitz, bend))
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
