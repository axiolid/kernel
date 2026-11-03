//! A certified, random-access position table of intrinsic curves (#361).
//!
//! [`IntrinsicPositions`](super::IntrinsicPositions) (#317) marches the
//! midpoint rule forward from the last point it read: its answers depend on
//! the order of the queries, a backward read starts again from zero, and
//! its `O(h^2)` rule needs about a million panels for a 10 km alignment at
//! a thousandth of a millimetre. [`IntrinsicPositionTable`] reads any point
//! of the same curve in any order, in at most `CHECKPOINT_PANELS` (64)
//! panels of a fourth-order rule, and bounds the error of every point it
//! returns.
//!
//! # The rule
//!
//! A point is the integral of the unit tangent `f = e^(i theta)` in the
//! curve's start frame, `theta' = k`, `theta'' = k'`, `theta''' = k''`,
//! `theta'''' = k'''` the closed-form derivatives of the curvature law. By
//! the chain rule
//!
//! `f'''' = (i theta'''' - 4 theta' theta''' - 3 theta''^2 - 6 i theta'^2 theta'' + theta'^4) f`,
//!
//! so on a seam-free piece of the law
//!
//! `|f''''| <= B4 = sup|k'''| + 4 sup|k| sup|k''| + 3 sup k'^2 + 6 sup k^2 sup|k'| + sup k^4`,
//!
//! each supremum from the law's closed-form bound term by term
//! (`crate::chain::curvature_bound`). The two-point Gauss-Legendre rule on
//! a panel of width `w` has the error `w^5 / 4320 f''''(xi)` for a real
//! integrand. Its Peano kernel keeps one sign, so the vector form
//! `|E| <= w^5 / 4320 B4` holds without a factor `sqrt 2`.
//!
//! # The grid
//!
//! Each seam-free piece of `[0, end]` (the law may jump at a seam, so no
//! panel straddles one) is cut into equal panels no wider than `step`,
//! chosen so the remainders over the whole curve sum to half the requested
//! error. Every 64 panels the table stores the running sums and the
//! remainder so far. A point at `s` starts from the last checkpoint at or
//! before the panel holding `s`, sums the full panels up to it, and adds
//! one partial panel from the grid line to `s` by the same rule, with its
//! own remainder `(s - grid)^5 / 4320 B4`. A query is a pure function of
//! `s`.
//!
//! The rounding terms are those of the midpoint positions: the heading's,
//! at most `64 eps (1 + K L) (1 + seams)` radians (`K` bounding `|k|`, `L`
//! the end), which turns each tangent by as much so a point at `s` moves by
//! at most `s` times it; the compensated sums' `16 eps s`; and the start
//! frame's stretch, applied to all of these, with `8 eps (|origin| +
//! stretch s)` for placing the point.

use axiolid_core::{Point2, Point3, Scalar};
use axiolid_curve::{Intrinsic2, Intrinsic3};

use super::positions::{add, unit_frame, Placement};
use super::{frame_stretch2, frame_stretch3, inflate};
use crate::chain::curvature_bound;

/// Most panels one table may use.
pub const MAX_TABLE_PANELS: usize = 1 << 22;

/// Panels between two stored checkpoints, and so most full panels a query
/// sums.
const CHECKPOINT_PANELS: usize = 64;

/// How far past the curve's length an end may lie, relative to `max(1,
/// length)`, and still be read as the length.
const END_SLACK: Scalar = 1e-9;

/// The running sums at a panel boundary.
#[derive(Debug, Clone, Copy)]
struct Checkpoint {
    sum: [(Scalar, Scalar); 2],
    remainder: Scalar,
}

/// A seam-free piece of `[0, end]`, cut into equal panels.
#[derive(Debug, Clone, Copy)]
struct Span {
    lo: Scalar,
    hi: Scalar,
    width: Scalar,
    panels: usize,
    /// The index of its first panel in the whole grid.
    base: usize,
    /// `sup |f''''|` over it.
    b4: Scalar,
}

/// Certified positions of an intrinsic curve over `[0, end]`, read in any
/// order: the unit tangent of its closed-form heading integrated by the
/// composite two-point Gauss-Legendre rule, each point with a certified
/// bound on its error (see the module documentation).
#[derive(Debug, Clone)]
pub struct IntrinsicPositionTable {
    /// The curve in its own start frame.
    local: Intrinsic2,
    placement: Placement,
    stretch: Scalar,
    origin: Scalar,
    end: Scalar,
    spans: Vec<Span>,
    checkpoints: Vec<Checkpoint>,
    heading_rounding: Scalar,
    panels: usize,
}

/// The remainder of one panel of width `w`.
fn remainder(width: Scalar, b4: Scalar) -> Scalar {
    width.powi(5) / 4320.0 * b4
}

impl IntrinsicPositionTable {
    /// The table of a 2D intrinsic curve over `[0, end]`, its points
    /// within about `error`, half of it spent on the remainders; `None`
    /// for an error or an end that is not finite and positive, an end past
    /// the curve's length, a law with no closed-form bound of `k` to `k'''`
    /// on a seam-free piece, or more than [`MAX_TABLE_PANELS`] panels.
    #[must_use]
    pub fn new2(curve: &Intrinsic2, error: Scalar, end: Scalar) -> Option<Self> {
        let local = Intrinsic2::new(unit_frame(), curve.curvature.clone(), curve.length);
        Self::build(
            local,
            Placement::Plane(curve.start),
            frame_stretch2(&curve.start),
            curve.start.origin.abs().max_element(),
            error,
            end,
        )
    }

    /// The table of a planar space curve (its torsion law identically zero,
    /// [`Intrinsic3::is_planar`]); `None` as for [`Self::new2`], or for a
    /// curve with torsion.
    #[must_use]
    pub fn new3(curve: &Intrinsic3, error: Scalar, end: Scalar) -> Option<Self> {
        if !curve.is_planar() {
            return None;
        }
        let local = Intrinsic2::new(unit_frame(), curve.curvature.clone(), curve.length);
        Self::build(
            local,
            Placement::Space(curve.start),
            frame_stretch3(&curve.start),
            curve.start.origin.abs().max_element(),
            error,
            end,
        )
    }

    fn build(
        local: Intrinsic2,
        placement: Placement,
        stretch: Scalar,
        origin: Scalar,
        error: Scalar,
        end: Scalar,
    ) -> Option<Self> {
        let length = local.length;
        if !(error.is_finite() && error > 0.0 && length.is_finite() && length > 0.0) {
            return None;
        }
        if !(end.is_finite() && end > 0.0 && end <= length + END_SLACK * length.max(1.0)) {
            return None;
        }
        let end = end.min(length);
        let law = &local.curvature;
        let rate = law.derivative();
        let second = rate.derivative();
        let third = second.derivative();
        let seams: Vec<Scalar> = law
            .seams_within(end)
            .into_iter()
            .filter(|&s| s > 0.0 && s < end)
            .collect();
        let mut ends = vec![0.0];
        ends.extend(seams.iter().copied());
        ends.push(end);
        let mut pieces = Vec::with_capacity(ends.len() - 1);
        let mut worst: Scalar = 0.0;
        for pair in ends.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            let k = curvature_bound(law, lo, hi)?;
            let k1 = curvature_bound(&rate, lo, hi)?;
            let k2 = curvature_bound(&second, lo, hi)?;
            let k3 = curvature_bound(&third, lo, hi)?;
            let b4 = inflate(k3 + 4.0 * k * k2 + 3.0 * k1 * k1 + 6.0 * k * k * k1 + k.powi(4));
            if !b4.is_finite() {
                return None;
            }
            worst = worst.max(b4);
            pieces.push((lo, hi, b4));
        }
        let budget = 0.5 * error / stretch;
        let step = if worst > 0.0 {
            (4320.0 * budget / (end * worst)).powf(0.25).min(end)
        } else {
            end
        };
        if !(step.is_finite() && step > 0.0) {
            return None;
        }
        let mut spans = Vec::with_capacity(pieces.len());
        let mut total = 0usize;
        for (lo, hi, b4) in pieces {
            let count = ((hi - lo) / step).ceil().max(1.0);
            if count > MAX_TABLE_PANELS as Scalar || total + count as usize > MAX_TABLE_PANELS {
                return None;
            }
            let panels = count as usize;
            spans.push(Span {
                lo,
                hi,
                width: (hi - lo) / count,
                panels,
                base: total,
                b4,
            });
            total += panels;
        }
        let turning = curvature_bound(law, 0.0, end)? * end;
        let heading_rounding =
            64.0 * Scalar::EPSILON * (1.0 + turning) * (1.0 + seams.len() as Scalar);
        let mut table = Self {
            local,
            placement,
            stretch,
            origin,
            end,
            spans,
            checkpoints: Vec::with_capacity(total / CHECKPOINT_PANELS + 1),
            heading_rounding,
            panels: total,
        };
        let mut state = Checkpoint {
            sum: [(0.0, 0.0); 2],
            remainder: 0.0,
        };
        for index in 0..total {
            if index % CHECKPOINT_PANELS == 0 {
                table.checkpoints.push(state);
            }
            table.add_panel(index, &mut state)?;
        }
        Some(table)
    }

    /// The number of panels the table is cut into.
    #[must_use]
    pub fn panels(&self) -> usize {
        self.panels
    }

    /// The span holding grid panel `index`.
    fn span_of(&self, index: usize) -> &Span {
        let at = self.spans.partition_point(|s| s.base + s.panels <= index);
        &self.spans[at.min(self.spans.len() - 1)]
    }

    /// The start of panel `i` of `span`.
    fn grid(span: &Span, i: usize) -> Scalar {
        if i == 0 {
            span.lo
        } else {
            span.lo + i as Scalar * span.width
        }
    }

    /// Add the integral of `f` over `[a, a + w]` by two Gauss points to
    /// `sum`.
    fn rule(&self, a: Scalar, w: Scalar, sum: &mut [(Scalar, Scalar); 2]) -> Option<()> {
        let half = 0.5 * w;
        let offset = half / (3.0 as Scalar).sqrt();
        for x in [a + half - offset, a + half + offset] {
            let theta = self.local.heading_at(x)?;
            add(&mut sum[0], half * theta.cos());
            add(&mut sum[1], half * theta.sin());
        }
        Some(())
    }

    /// Add grid panel `index` to `state`.
    fn add_panel(&self, index: usize, state: &mut Checkpoint) -> Option<()> {
        let span = self.span_of(index);
        let a = Self::grid(span, index - span.base);
        self.rule(a, span.width, &mut state.sum)?;
        state.remainder += remainder(span.width, span.b4);
        Some(())
    }

    /// The point at `s` in the start frame's coordinates, and its error
    /// there.
    fn local_point(&self, s: Scalar) -> Option<(Point2, Scalar)> {
        if !(s.is_finite() && (0.0..=self.end).contains(&s)) {
            return None;
        }
        let at = self.spans.partition_point(|span| span.hi < s);
        let span = self.spans.get(at)?;
        let i = (((s - span.lo) / span.width).floor().max(0.0) as usize).min(span.panels - 1);
        let index = span.base + i;
        let first = index - index % CHECKPOINT_PANELS;
        let mut state = self.checkpoints[index / CHECKPOINT_PANELS];
        for panel in first..index {
            self.add_panel(panel, &mut state)?;
        }
        let from = Self::grid(span, i);
        let rest = (s - from).max(0.0);
        if rest > 0.0 {
            self.rule(from, rest, &mut state.sum)?;
            state.remainder += remainder(rest, span.b4);
        }
        let p = Point2::new(
            state.sum[0].0 + state.sum[0].1,
            state.sum[1].0 + state.sum[1].1,
        );
        let error = state.remainder + s * (self.heading_rounding + 16.0 * Scalar::EPSILON);
        Some((p, error))
    }

    /// A local error through the start frame, with the placement's
    /// rounding.
    fn placed_error(&self, local: Scalar, s: Scalar) -> Scalar {
        inflate(self.stretch * local + 8.0 * Scalar::EPSILON * (self.origin + self.stretch * s))
    }

    /// The point at arc length `s` of a 2D curve and a certified bound on
    /// its distance from the exact point: the same bits whatever was read
    /// before. `None` outside `[0, end]`, for a law that does not integrate
    /// there, or for a space curve.
    #[must_use]
    pub fn point2(&self, s: Scalar) -> Option<(Point2, Scalar)> {
        let Placement::Plane(frame) = self.placement else {
            return None;
        };
        let (p, local) = self.local_point(s)?;
        let out = Point2::new(
            frame.origin.x + frame.x.x * p.x + frame.y.x * p.y,
            frame.origin.y + frame.x.y * p.x + frame.y.y * p.y,
        );
        Some((out, self.placed_error(local, s)))
    }

    /// [`Self::point2`] for a planar space curve.
    #[must_use]
    pub fn point3(&self, s: Scalar) -> Option<(Point3, Scalar)> {
        let Placement::Space(frame) = self.placement else {
            return None;
        };
        let (p, local) = self.local_point(s)?;
        Some((
            frame.origin + frame.x * p.x + frame.y * p.y,
            self.placed_error(local, s),
        ))
    }

    /// A certified bound on the rounding of the heading, in radians: what a
    /// tangent read from it may be turned by.
    #[must_use]
    pub fn heading_rounding(&self) -> Scalar {
        self.heading_rounding
    }
}
