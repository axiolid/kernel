//! Every component of a field's zero set in a parameter box, as
//! [`ImplicitCurve2`]s (ADR 0077).
//!
//! The box is subdivided until each piece either certainly misses the
//! curve (the field's bound excludes zero) or is *regular*: one partial
//! derivative's bound excludes zero, so along that parameter the field is
//! strictly monotone and the curve crosses each line across the piece at
//! most once. A regular piece holds the curve as a graph over the other
//! parameter, existing exactly where the field has opposite signs on the
//! piece's two sides -- decided by isolating the field's roots along those
//! sides, each root certified simple. Those stretches are the curve's
//! cells. Cells of neighbouring pieces meet where the curve crosses their
//! common side, so chaining them by their end points gives the components,
//! closed loops or open chains leaving the box. Nothing is marched and no
//! step size is guessed: a component is found because some regular piece
//! must contain part of it, and pieces cover the box.
//!
//! A piece where neither partial can be bounded away from zero while the
//! field may vanish, down to the smallest size allowed, may hold a singular
//! point of the curve (where the surfaces touch): the trace is refused
//! there by name rather than guessed through.

use axiolid_core::{Point2, Scalar};
use axiolid_curve::{Axis, Field2, ImplicitCell, ImplicitCurve2};
use core::f64::consts::TAU;

use axiolid_curve::implicit::{bound, bound_simple, partial, Cell, Range};

/// Why a trace was refused.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum TraceRefusal {
    /// The field and both partials may vanish together near this point: a
    /// singular point of the curve (the surfaces touch), or two branches
    /// closer than the trace resolves.
    Singular(Point2),
    /// The subdivision exceeded its work budget.
    Budget,
}

/// Periodicity of the box's parameters: a component leaving through one
/// side of a periodic parameter comes back through the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Periodic {
    pub(crate) u: bool,
    pub(crate) v: bool,
}

/// Every component of `field = 0` in `domain`.
pub(crate) fn trace(
    field: &Field2,
    domain: Cell,
    periodic: Periodic,
) -> Result<Vec<ImplicitCurve2>, TraceRefusal> {
    let du = partial(field, true);
    let dv = partial(field, false);
    let extent = (domain.hi - domain.lo).abs();
    let smallest = extent * 1e-9;
    // Start from a grid so that the bounds are tight from the outset; its
    // lines sit off the round fractions a symmetric input would put a
    // crossing on.
    let mut queue = split_grid(domain, 8);
    let mut all: Vec<ImplicitCell> = Vec::new();
    let mut work = 0usize;
    while let Some((cell, depth)) = queue.pop() {
        work += 1;
        if work > 400_000 {
            return Err(TraceRefusal::Budget);
        }
        if !bound(field, &du, &dv, &cell).straddles_zero() {
            continue;
        }
        let fu = bound_simple(&du, &cell);
        let fv = bound_simple(&dv, &cell);
        // Monotone in v: a graph over u. Where both hold, the steeper.
        let axis = match (!fu.straddles_zero(), !fv.straddles_zero()) {
            (true, true) => {
                if steeper(fv, &cell, false) >= steeper(fu, &cell, true) {
                    Some(Axis::U)
                } else {
                    Some(Axis::V)
                }
            }
            (false, true) => Some(Axis::U),
            (true, false) => Some(Axis::V),
            (false, false) => None,
        };
        match axis {
            Some(axis) => match cells_in(field, &du, &dv, &cell, axis) {
                Ok(found) => all.extend(found),
                Err(()) => {
                    // A root on a side could not be certified simple: move
                    // the side by splitting off-centre.
                    if depth > 60 {
                        return Err(TraceRefusal::Singular(cell.centre()));
                    }
                    for piece in split(cell, 0.4129) {
                        queue.push((piece, depth + 1));
                    }
                }
            },
            None => {
                if (cell.hi.x - cell.lo.x) <= smallest.x && (cell.hi.y - cell.lo.y) <= smallest.y {
                    return Err(TraceRefusal::Singular(cell.centre()));
                }
                for piece in split(cell, 0.5) {
                    queue.push((piece, depth + 1));
                }
            }
        }
    }
    Ok(chain(field, all, domain, periodic))
}

/// How strongly a partial's bound keeps the field monotone, scaled by the
/// cell's size along that parameter.
fn steeper(r: Range, cell: &Cell, along_u: bool) -> Scalar {
    let width = if along_u {
        cell.hi.x - cell.lo.x
    } else {
        cell.hi.y - cell.lo.y
    };
    r.lo.abs().min(r.hi.abs()) * width
}

fn split_grid(domain: Cell, n: usize) -> Vec<(Cell, u32)> {
    let line = |k: usize, lo: Scalar, hi: Scalar| -> Scalar {
        if k == 0 {
            lo
        } else if k == n {
            hi
        } else {
            lo + (hi - lo) * ((k as Scalar + 0.0371) / n as Scalar)
        }
    };
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let lo = Point2::new(
                line(i, domain.lo.x, domain.hi.x),
                line(j, domain.lo.y, domain.hi.y),
            );
            let hi = Point2::new(
                line(i + 1, domain.lo.x, domain.hi.x),
                line(j + 1, domain.lo.y, domain.hi.y),
            );
            out.push((Cell { lo, hi }, 0));
        }
    }
    out
}

/// Four pieces, split at `fraction` of each side.
fn split(cell: Cell, fraction: Scalar) -> [Cell; 4] {
    let m = cell.lo + (cell.hi - cell.lo) * fraction;
    [
        Cell { lo: cell.lo, hi: m },
        Cell {
            lo: Point2::new(m.x, cell.lo.y),
            hi: Point2::new(cell.hi.x, m.y),
        },
        Cell {
            lo: Point2::new(cell.lo.x, m.y),
            hi: Point2::new(m.x, cell.hi.y),
        },
        Cell { lo: m, hi: cell.hi },
    ]
}

/// The field along one side of a regular piece, as a function of the free
/// parameter: its value at `x` and a bound over `[a, b]`.
struct Side<'a> {
    field: &'a Field2,
    d_free: &'a Field2,
    axis: Axis,
    fixed: Scalar,
}

impl Side<'_> {
    fn at(&self, x: Scalar) -> Point2 {
        match self.axis {
            Axis::U => Point2::new(x, self.fixed),
            Axis::V => Point2::new(self.fixed, x),
        }
    }

    fn value(&self, x: Scalar) -> Scalar {
        self.field.value(self.at(x))
    }

    fn cell(&self, a: Scalar, b: Scalar) -> Cell {
        let (p, q) = (self.at(a), self.at(b));
        Cell {
            lo: p.min(q),
            hi: p.max(q),
        }
    }

    fn range(&self, a: Scalar, b: Scalar) -> Range {
        bound_simple(self.field, &self.cell(a, b))
    }

    fn slope(&self, a: Scalar, b: Scalar) -> Range {
        bound_simple(self.d_free, &self.cell(a, b))
    }

    /// The roots in `(a, b)`, each certified simple; `Err` where a root
    /// cannot be separated or certified.
    fn roots(&self, a: Scalar, b: Scalar, out: &mut Vec<Scalar>, depth: u32) -> Result<(), ()> {
        if !self.range(a, b).straddles_zero() {
            return Ok(());
        }
        let (fa, fb) = (self.value(a), self.value(b));
        if !self.slope(a, b).straddles_zero() {
            // Monotone: one root exactly where the signs differ.
            if fa == 0.0 || fb == 0.0 {
                // A root at an end: it belongs to the neighbouring piece's
                // side as much as to this one, and is found as an end there.
                return Err(());
            }
            if (fa < 0.0) != (fb < 0.0) {
                out.push(self.refine(a, b, fa));
            }
            return Ok(());
        }
        if depth > 60 || (b - a).abs() <= 1e-13 * (1.0 + a.abs().max(b.abs())) {
            return Err(());
        }
        let m = 0.5 * (a + b);
        self.roots(a, m, out, depth + 1)?;
        self.roots(m, b, out, depth + 1)
    }

    /// The root in `[a, b]`, where the field is monotone and changes sign.
    fn refine(&self, mut a: Scalar, mut b: Scalar, fa: Scalar) -> Scalar {
        let negative_at_a = fa < 0.0;
        for _ in 0..200 {
            let m = 0.5 * (a + b);
            if m <= a.min(b) || m >= a.max(b) {
                break;
            }
            let fm = self.value(m);
            if fm == 0.0 {
                return m;
            }
            if (fm < 0.0) == negative_at_a {
                a = m;
            } else {
                b = m;
            }
        }
        0.5 * (a + b)
    }
}

/// The cells of the curve in a regular piece, where the field is monotone
/// along the solved parameter (the one that is not `axis`).
fn cells_in(
    field: &Field2,
    du: &Field2,
    dv: &Field2,
    cell: &Cell,
    axis: Axis,
) -> Result<Vec<ImplicitCell>, ()> {
    let (d_free, lo_free, hi_free, lo_solved, hi_solved) = match axis {
        Axis::U => (du, cell.lo.x, cell.hi.x, cell.lo.y, cell.hi.y),
        Axis::V => (dv, cell.lo.y, cell.hi.y, cell.lo.x, cell.hi.x),
    };
    let low = Side {
        field,
        d_free,
        axis,
        fixed: lo_solved,
    };
    let high = Side {
        field,
        d_free,
        axis,
        fixed: hi_solved,
    };
    let mut breaks = vec![lo_free, hi_free];
    low.roots(lo_free, hi_free, &mut breaks, 0)?;
    high.roots(lo_free, hi_free, &mut breaks, 0)?;
    breaks.sort_by(Scalar::total_cmp);
    let mut out = Vec::new();
    for pair in breaks.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b - a <= 1e-14 * (1.0 + a.abs().max(b.abs())) {
            continue;
        }
        let m = 0.5 * (a + b);
        let (f0, f1) = (low.value(m), high.value(m));
        if (f0 < 0.0) != (f1 < 0.0) && f0 != 0.0 && f1 != 0.0 {
            out.push(ImplicitCell {
                axis,
                from: a,
                to: b,
                low: lo_solved,
                high: hi_solved,
            });
        }
    }
    Ok(out)
}

/// A cell's end points in `(u, v)`, at its `from` and `to`.
fn ends(field: &Field2, cell: &ImplicitCell) -> Option<(Point2, Point2)> {
    let curve = ImplicitCurve2 {
        field: field.clone(),
        cells: vec![*cell],
    };
    Some((curve.point(0.0)?, curve.point(1.0)?))
}

/// Shift a cell by whole periods.
fn shifted(cell: &ImplicitCell, du: Scalar, dv: Scalar) -> ImplicitCell {
    let (df, ds) = match cell.axis {
        Axis::U => (du, dv),
        Axis::V => (dv, du),
    };
    ImplicitCell {
        axis: cell.axis,
        from: cell.from + df,
        to: cell.to + df,
        low: cell.low + ds,
        high: cell.high + ds,
    }
}

fn reversed(cell: &ImplicitCell) -> ImplicitCell {
    ImplicitCell {
        from: cell.to,
        to: cell.from,
        ..*cell
    }
}

/// Chain cells end to end into components.
fn chain(
    field: &Field2,
    cells: Vec<ImplicitCell>,
    domain: Cell,
    periodic: Periodic,
) -> Vec<ImplicitCurve2> {
    let points: Vec<(Point2, Point2)> = cells
        .iter()
        .map(|c| ends(field, c).unwrap_or((Point2::splat(Scalar::NAN), Point2::splat(Scalar::NAN))))
        .collect();
    let extent = (domain.hi - domain.lo).abs();
    let eps = 1e-8 * (1.0 + extent.x.max(extent.y));
    // Wrap a point into the domain for matching across periods.
    let wrap = |p: Point2| -> Point2 {
        let mut q = p;
        if periodic.u {
            q.x = domain.lo.x + (q.x - domain.lo.x).rem_euclid(TAU);
            if (q.x - domain.lo.x - TAU).abs() <= eps {
                q.x = domain.lo.x;
            }
        }
        if periodic.v {
            q.y = domain.lo.y + (q.y - domain.lo.y).rem_euclid(TAU);
            if (q.y - domain.lo.y - TAU).abs() <= eps {
                q.y = domain.lo.y;
            }
        }
        q
    };
    let near = |p: Point2, q: Point2| {
        let (a, b) = (wrap(p), wrap(q));
        let mut d = a - b;
        if periodic.u && d.x.abs() > TAU - eps {
            d.x = d.x.abs() - TAU;
        }
        if periodic.v && d.y.abs() > TAU - eps {
            d.y = d.y.abs() - TAU;
        }
        d.length() <= eps
    };
    // Which other cell end meets this one.
    let find = |me: usize, p: Point2, used: &[bool]| -> Option<(usize, bool)> {
        let mut best: Option<(usize, bool, Scalar)> = None;
        for (k, (a, b)) in points.iter().enumerate() {
            if k == me || used[k] {
                continue;
            }
            for (at_start, q) in [(true, *a), (false, *b)] {
                if near(p, q) {
                    let d = (wrap(p) - wrap(q)).length();
                    if best.is_none_or(|(_, _, e)| d < e) {
                        best = Some((k, at_start, d));
                    }
                }
            }
        }
        best.map(|(k, s, _)| (k, s))
    };
    let mut used = vec![false; cells.len()];
    let mut out = Vec::new();
    for first in 0..cells.len() {
        if used[first] || points[first].0.x.is_nan() {
            continue;
        }
        used[first] = true;
        let mut chain = vec![cells[first]];
        let mut tail = points[first].1;
        let mut closed = false;
        // Forward.
        loop {
            if near(tail, points[first].0) && chain.len() > 1 {
                closed = true;
                break;
            }
            let Some((k, at_start)) = find(usize::MAX, tail, &used) else {
                break;
            };
            used[k] = true;
            let (next, (start, end)) = if at_start {
                (cells[k], points[k])
            } else {
                (reversed(&cells[k]), (points[k].1, points[k].0))
            };
            // Unwrap: shift the cell so it starts where the chain is.
            let du = if periodic.u {
                ((tail.x - start.x) / TAU).round() * TAU
            } else {
                0.0
            };
            let dv = if periodic.v {
                ((tail.y - start.y) / TAU).round() * TAU
            } else {
                0.0
            };
            chain.push(shifted(&next, du, dv));
            tail = Point2::new(end.x + du, end.y + dv);
        }
        if !closed {
            // Backward from the first cell's start.
            let mut head = points[first].0;
            loop {
                let Some((k, at_start)) = find(usize::MAX, head, &used) else {
                    break;
                };
                used[k] = true;
                // The cell must end at `head`.
                let (prev, (start, end)) = if at_start {
                    (reversed(&cells[k]), (points[k].1, points[k].0))
                } else {
                    (cells[k], points[k])
                };
                let du = if periodic.u {
                    ((head.x - end.x) / TAU).round() * TAU
                } else {
                    0.0
                };
                let dv = if periodic.v {
                    ((head.y - end.y) / TAU).round() * TAU
                } else {
                    0.0
                };
                chain.insert(0, shifted(&prev, du, dv));
                head = Point2::new(start.x + du, start.y + dv);
            }
        }
        out.push(ImplicitCurve2 {
            field: field.clone(),
            cells: chain,
        });
    }
    out
}
