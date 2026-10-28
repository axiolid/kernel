//! Parameters along an edge and the interval algebra on them.
//!
//! The edge runs from vertex `start` (parameter 0) to `end` (parameter 1).
//! A constraint is affine in the parameter, `L(t) = (1-t)·a + t·b` with
//! `a`, `b` its values at the two ends, so it changes sign at
//! `t = a / (a - b)` and nowhere else. A root is kept as its constraint and
//! never evaluated: two roots are ordered by the sign of a polynomial in
//! the constraints' values, decided exactly. Zero-measure intersections are
//! ignored, so open and closed ends never matter; only an identically zero
//! constraint needs its strictness.

use std::cmp::Ordering;

use axiolid_contracts::GeomResult;

use crate::decide::decide;
use crate::recipe::{Constraint, Delta, Kind, RootOrder, ValueSign};
use crate::view::Vert;

/// The root of a constraint whose end values have opposite signs, so it lies
/// strictly inside the edge.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Root {
    pub constraint: Constraint,
    /// Sign of `a - b`: positive when the constraint holds before the root.
    pub diff: Ordering,
}

/// A parameter in `[0, 1]`.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Param {
    Zero,
    One,
    Root(Root),
}

/// A parameter interval of positive length.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Span {
    pub lo: Param,
    pub hi: Param,
}

impl Span {
    pub(crate) const FULL: Self = Self {
        lo: Param::Zero,
        hi: Param::One,
    };
}

/// Where along the edge a constraint holds.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Holds {
    Nowhere,
    Everywhere,
    /// From the start to the root.
    Before(Root),
    /// From the root to the end.
    After(Root),
}

/// Where a constraint holds, from its signs at the two ends.
fn classify(constraint: &Constraint, a: Ordering, b: Ordering) -> Holds {
    let root = |diff| Root {
        constraint: *constraint,
        diff,
    };
    match (a, b) {
        (Ordering::Equal, Ordering::Equal) if constraint.strict => Holds::Nowhere,
        (Ordering::Equal, Ordering::Equal) => Holds::Everywhere,
        (Ordering::Greater, Ordering::Less) => Holds::Before(root(Ordering::Greater)),
        (Ordering::Less, Ordering::Greater) => Holds::After(root(Ordering::Less)),
        (Ordering::Less | Ordering::Equal, Ordering::Less | Ordering::Equal) => Holds::Nowhere,
        (Ordering::Greater | Ordering::Equal, Ordering::Greater | Ordering::Equal) => {
            Holds::Everywhere
        }
    }
}

/// Decides signs and orders parameters along one edge.
pub(crate) struct Solver<'a> {
    pub verts: &'a [Vert],
    pub start: u32,
    pub end: u32,
    /// Decisions the interval filter could not settle.
    pub escalations: usize,
}

impl<'a> Solver<'a> {
    pub(crate) fn new(verts: &'a [Vert], start: u32, end: u32) -> Self {
        Self {
            verts,
            start,
            end,
            escalations: 0,
        }
    }

    /// The signs of `constraint` at the edge's two ends.
    fn sides(&mut self, constraint: &Constraint) -> GeomResult<(Ordering, Ordering)> {
        let mut at = |vertex| {
            decide(
                &ValueSign {
                    verts: self.verts,
                    constraint,
                    at: vertex,
                },
                &mut self.escalations,
            )
        };
        Ok((at(self.start)?, at(self.end)?))
    }

    /// Where `constraint` holds along the edge.
    pub(crate) fn holds(&mut self, constraint: &Constraint) -> GeomResult<Holds> {
        let (a, b) = self.sides(constraint)?;
        Ok(classify(constraint, a, b))
    }

    /// Whether `constraint` is zero along the whole edge.
    pub(crate) fn identically_zero(&mut self, constraint: &Constraint) -> GeomResult<bool> {
        Ok(self.sides(constraint)? == (Ordering::Equal, Ordering::Equal))
    }

    /// `span ∩ {constraint holds}` for an image-side constraint, with the
    /// image point perturbed by `δ = (ε, ε²)` where the constraint is zero
    /// along the whole edge. The perturbed point is strictly on one side of
    /// the line `p'q'`: `orient2d(p, q, P + δ) = orient2d(p, q, P) - dy·ε +
    /// dx·ε²` with `(dx, dy) = q' - p'`, so the side is `-sign(dy)` and, when
    /// `dy` is zero, `sign(dx)`.
    pub(crate) fn clip_perturbed(
        &mut self,
        span: Span,
        constraint: &Constraint,
    ) -> GeomResult<Option<Span>> {
        let (a, b) = self.sides(constraint)?;
        let Kind::ImageSide { p, q, neg } = constraint.kind else {
            unreachable!("only image-side constraints are perturbed");
        };
        if (a, b) != (Ordering::Equal, Ordering::Equal) {
            return self.apply(span, classify(constraint, a, b));
        }
        let mut delta = |axis| {
            decide(
                &Delta {
                    verts: self.verts,
                    from: p,
                    to: q,
                    axis,
                },
                &mut self.escalations,
            )
        };
        let (dx, dy) = (delta(0)?, delta(1)?);
        let side = if dy == Ordering::Equal {
            dx
        } else {
            dy.reverse()
        };
        let side = if neg { side.reverse() } else { side };
        Ok((side == Ordering::Greater).then_some(span))
    }

    /// The exact order of two parameters.
    pub(crate) fn cmp(&mut self, x: Param, y: Param) -> GeomResult<Ordering> {
        fn rank(p: Param) -> u8 {
            match p {
                Param::Zero => 0,
                Param::Root(_) => 1,
                Param::One => 2,
            }
        }
        let (Param::Root(rx), Param::Root(ry)) = (x, y) else {
            return Ok(rank(x).cmp(&rank(y)));
        };
        if rx.constraint.kind == ry.constraint.kind {
            return Ok(Ordering::Equal);
        }
        // r_x - r_y = (a_y·b_x - a_x·b_y) / ((a_x - b_x)·(a_y - b_y)).
        let mut order = decide(
            &RootOrder {
                verts: self.verts,
                x: &rx.constraint,
                y: &ry.constraint,
                start: self.start,
                end: self.end,
            },
            &mut self.escalations,
        )?;
        for diff in [rx.diff, ry.diff] {
            if diff == Ordering::Less {
                order = order.reverse();
            }
        }
        Ok(order)
    }

    fn later(&mut self, x: Param, y: Param) -> GeomResult<Param> {
        Ok(if self.cmp(x, y)? == Ordering::Less {
            y
        } else {
            x
        })
    }

    fn earlier(&mut self, x: Param, y: Param) -> GeomResult<Param> {
        Ok(if self.cmp(x, y)? == Ordering::Greater {
            y
        } else {
            x
        })
    }

    /// `lo < hi` exactly.
    fn positive(&mut self, lo: Param, hi: Param) -> GeomResult<bool> {
        Ok(self.cmp(lo, hi)? == Ordering::Less)
    }

    /// `span ∩ {constraint holds}`, if it has positive length.
    pub(crate) fn clip(&mut self, span: Span, constraint: &Constraint) -> GeomResult<Option<Span>> {
        let holds = self.holds(constraint)?;
        self.apply(span, holds)
    }

    fn apply(&mut self, span: Span, holds: Holds) -> GeomResult<Option<Span>> {
        let (lo, hi) = match holds {
            Holds::Nowhere => return Ok(None),
            Holds::Everywhere => return Ok(Some(span)),
            Holds::Before(root) => (span.lo, self.earlier(span.hi, Param::Root(root))?),
            Holds::After(root) => (self.later(span.lo, Param::Root(root))?, span.hi),
        };
        Ok(self.positive(lo, hi)?.then_some(Span { lo, hi }))
    }

    /// Each span of `list` clipped by `constraint`.
    pub(crate) fn clip_all(
        &mut self,
        list: Vec<Span>,
        constraint: &Constraint,
    ) -> GeomResult<Vec<Span>> {
        let mut out = Vec::with_capacity(list.len());
        for span in list {
            out.extend(self.clip(span, constraint)?);
        }
        Ok(out)
    }

    /// Add `span` to the sorted, disjoint, non-touching `list`, merging
    /// every span it overlaps or touches.
    pub(crate) fn insert(&mut self, list: &mut Vec<Span>, span: Span) -> GeomResult<()> {
        // First span ending at or after the new start.
        let (mut low, mut high) = (0, list.len());
        while low < high {
            let mid = low + (high - low) / 2;
            if self.cmp(list[mid].hi, span.lo)? == Ordering::Less {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        let first = low;
        let mut merged = span;
        let mut next = first;
        while next < list.len() && self.cmp(list[next].lo, merged.hi)? != Ordering::Greater {
            merged.lo = self.earlier(merged.lo, list[next].lo)?;
            merged.hi = self.later(merged.hi, list[next].hi)?;
            next += 1;
        }
        list.splice(first..next, [merged]);
        Ok(())
    }

    /// The positive-length overlaps of two sorted, disjoint lists.
    pub(crate) fn intersect(&mut self, a: &[Span], b: &[Span]) -> GeomResult<Vec<Span>> {
        let mut out = Vec::new();
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            let lo = self.later(a[i].lo, b[j].lo)?;
            let hi = self.earlier(a[i].hi, b[j].hi)?;
            if self.positive(lo, hi)? {
                out.push(Span { lo, hi });
            }
            if self.cmp(a[i].hi, b[j].hi)? == Ordering::Less {
                i += 1;
            } else {
                j += 1;
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use axiolid_exact::{Dyadic, Interval};

    use super::*;
    use crate::recipe::Kind;

    fn vert(u: f64, v: f64, w: f64) -> Vert {
        let exact = [u, v, w].map(|x| Dyadic::try_from_f64(x).expect("finite"));
        Vert {
            exact,
            enc: [u, v, w].map(Interval::point),
        }
    }

    /// An edge along U from `from` to `to` at V = 0, W = `w`.
    fn edge(from: f64, to: f64, w: f64) -> [Vert; 2] {
        [vert(from, 0.0, w), vert(to, 0.0, w)]
    }

    fn lower(bound: f64) -> Constraint {
        Constraint::new(
            Kind::Region {
                axis: 0,
                upper: false,
                bound,
            },
            false,
        )
    }

    fn upper(bound: f64) -> Constraint {
        Constraint::new(
            Kind::Region {
                axis: 0,
                upper: true,
                bound,
            },
            false,
        )
    }

    fn root(solver: &mut Solver<'_>, constraint: &Constraint) -> Param {
        match solver.holds(constraint).expect("decided") {
            Holds::Before(root) | Holds::After(root) => Param::Root(root),
            other => panic!("expected a root, got {other:?}"),
        }
    }

    #[test]
    fn roots_order_by_position_whatever_the_denominator_signs() {
        // Along U = 0..4. `lower(b)` has a negative denominator (a < b);
        // `upper(b)` a positive one. Roots at t = bound / 4.
        let verts = edge(0.0, 4.0, 0.0);
        let mut solver = Solver::new(&verts, 0, 1);
        let cases = [
            (lower(1.0), upper(3.0), Ordering::Less),
            (upper(3.0), lower(1.0), Ordering::Greater),
            (lower(3.0), upper(1.0), Ordering::Greater),
            (upper(1.0), lower(3.0), Ordering::Less),
            (lower(1.0), lower(3.0), Ordering::Less),
            (upper(3.0), upper(1.0), Ordering::Greater),
            (lower(2.0), upper(2.0), Ordering::Equal),
        ];
        for (x, y, expected) in cases {
            let (px, py) = (root(&mut solver, &x), root(&mut solver, &y));
            assert_eq!(solver.cmp(px, py).unwrap(), expected, "{x:?} vs {y:?}");
            assert_eq!(solver.cmp(py, px).unwrap(), expected.reverse());
        }
    }

    #[test]
    fn roots_order_when_the_edge_runs_backwards() {
        let verts = edge(4.0, 0.0, 0.0);
        let mut solver = Solver::new(&verts, 0, 1);
        // Along t, U = 4 - 4t: lower(1) holds for t <= 3/4, upper(3) for t >= 1/4.
        let (px, py) = (
            root(&mut solver, &lower(1.0)),
            root(&mut solver, &upper(3.0)),
        );
        assert_eq!(solver.cmp(px, py).unwrap(), Ordering::Greater);
        assert!(matches!(
            solver.holds(&lower(1.0)).unwrap(),
            Holds::Before(Root {
                diff: Ordering::Greater,
                ..
            })
        ));
    }

    #[test]
    fn endpoints_bracket_every_root() {
        let verts = edge(0.0, 4.0, 0.0);
        let mut solver = Solver::new(&verts, 0, 1);
        let p = root(&mut solver, &lower(1.0));
        assert_eq!(solver.cmp(Param::Zero, p).unwrap(), Ordering::Less);
        assert_eq!(solver.cmp(p, Param::One).unwrap(), Ordering::Less);
        assert_eq!(solver.cmp(Param::One, Param::One).unwrap(), Ordering::Equal);
    }

    #[test]
    fn an_identically_zero_constraint_follows_its_strictness() {
        // W is 1 along the whole edge, so `near - W` with near = 1 is zero.
        let verts = edge(0.0, 4.0, 1.0);
        let mut solver = Solver::new(&verts, 0, 1);
        let closed = Constraint::new(Kind::BelowNear(1.0), false);
        let open = Constraint::new(Kind::BelowNear(1.0), true);
        assert!(matches!(solver.holds(&closed).unwrap(), Holds::Everywhere));
        assert!(matches!(solver.holds(&open).unwrap(), Holds::Nowhere));
        assert!(solver.clip(Span::FULL, &closed).unwrap().is_some());
        assert!(solver.clip(Span::FULL, &open).unwrap().is_none());
    }

    #[test]
    fn a_constraint_touching_zero_at_one_end_has_zero_measure_or_everything() {
        // lower(0) on U = 0..4: zero at the start, positive after.
        let verts = edge(0.0, 4.0, 0.0);
        let mut solver = Solver::new(&verts, 0, 1);
        assert!(matches!(
            solver.holds(&lower(0.0)).unwrap(),
            Holds::Everywhere
        ));
        // upper(0): zero at the start, negative after.
        assert!(matches!(solver.holds(&upper(0.0)).unwrap(), Holds::Nowhere));
    }

    #[test]
    fn insertion_merges_touching_and_overlapping_spans() {
        let verts = edge(0.0, 8.0, 0.0);
        let mut solver = Solver::new(&verts, 0, 1);
        let at = |solver: &mut Solver<'_>, bound| root(solver, &lower(bound));
        let (p1, p2, p3, p4) = (
            at(&mut solver, 1.0),
            at(&mut solver, 2.0),
            at(&mut solver, 3.0),
            at(&mut solver, 5.0),
        );
        let mut list = Vec::new();
        solver.insert(&mut list, Span { lo: p3, hi: p4 }).unwrap();
        solver.insert(&mut list, Span { lo: p1, hi: p2 }).unwrap();
        assert_eq!(list.len(), 2);
        // [2,3] touches both neighbours.
        solver.insert(&mut list, Span { lo: p2, hi: p3 }).unwrap();
        assert_eq!(list.len(), 1);
        assert_eq!(solver.cmp(list[0].lo, p1).unwrap(), Ordering::Equal);
        assert_eq!(solver.cmp(list[0].hi, p4).unwrap(), Ordering::Equal);
    }

    #[test]
    fn intersection_keeps_only_positive_length_overlaps() {
        let verts = edge(0.0, 8.0, 0.0);
        let mut solver = Solver::new(&verts, 0, 1);
        let (p1, p2, p3) = (
            root(&mut solver, &lower(1.0)),
            root(&mut solver, &lower(2.0)),
            root(&mut solver, &lower(3.0)),
        );
        let a = [Span { lo: p1, hi: p2 }];
        let touching = [Span { lo: p2, hi: p3 }];
        assert!(solver.intersect(&a, &touching).unwrap().is_empty());
        let overlapping = [Span { lo: p1, hi: p3 }];
        assert_eq!(solver.intersect(&a, &overlapping).unwrap().len(), 1);
    }
}
