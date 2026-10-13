//! Opt-in dropping of degenerate input rings, every drop reported (#256).
//!
//! By default a degenerate input ring refuses the whole operation with a
//! named [`OverlayError`]. [`OverlayOptions::dropping_degenerate_rings`]
//! asks [`crate::overlay_with`] and [`crate::union_soup_with`] instead to
//! leave such a ring out, to report it as a [`DroppedRing`], and to run
//! the exact boolean over the rings that remain. Nothing is dropped
//! silently and nothing is repaired: a ring is either used as given or
//! left out whole, and the report says which, why and how much area went
//! with it.
//!
//! # Which faults are droppable
//!
//! A fault is droppable when leaving the ring out has a stated, bounded
//! effect, and the ring is otherwise well formed:
//!
//! - [`DropReason::TooShort`]: fewer than three points, all finite. Such
//!   a ring encloses nothing, so its area is exactly zero.
//! - [`DropReason::ZeroArea`]: the crate's exact `ZeroArea` rule, an area
//!   of at most `tolerance.linear()^2`, decided exactly from the fan of
//!   triangles at the first vertex (#274). It is decided after the
//!   self-intersection test, so the ring is simple and its area is the
//!   area it encloses. The threshold is that rule, not a second one.
//! - [`DropReason::RepeatedVertex`]: an edge no longer than
//!   `tolerance.linear()`. It is droppable only when the ring with each
//!   run of such vertices merged into its first vertex does not intersect
//!   itself, so its signed area is the area it encloses and the reported
//!   area is what leaving it out removes. That area can be large: a
//!   square with a doubled corner is dropped whole, and the report says
//!   how much went.
//! - [`DropReason::OuterDropped`]: a hole whose outer ring was dropped
//!   goes with it, since a hole has no meaning without its outer. The
//!   hole must still be droppable or valid on its own; its own fault is
//!   reported in preference to this reason. Whether it lies inside the
//!   dropped outer is not tested.
//!
//! Never droppable, refused as without the option: a non-finite point
//! (no area can be stated), a ring that intersects itself (its interior
//! is ambiguous and may be large, so no bounded removal can be named), a
//! hole outside its kept outer ring, and an invalid frame. A
//! [`DropReason::RepeatedVertex`] ring that still intersects itself once
//! its short edges are merged is refused as
//! [`OverlayError::SelfIntersection`].
//!
//! A degenerate hole of a kept outer ring is dropped alone: the operand
//! then covers that hole's area too, which the report states as the
//! ring's area.
//!
//! # The reported area
//!
//! [`DroppedRing::area`] is the magnitude of the ring's signed area, the
//! fan of triangles from its first vertex summed exactly and rounded once
//! for output. [`DroppedRing::area_bound`] is a proven upper bound on that
//! exact magnitude. A [`DropReason::ZeroArea`] ring's exact area is at
//! most `tolerance.linear()^2` by the rule that dropped it.
use std::borrow::Cow;

use axiolid_core::{Point2, Tolerance};
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;

use crate::{hole_outside, self_intersects, validate_frame, validate_ring};
use crate::{OverlayError, OverlayInput, OverlayResult, Polygon, Ring};

/// What an operation does with a degenerate input ring (#256).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum DegenerateRings {
    /// Refuse the whole operation with the ring's named error. The default.
    #[default]
    Refuse,
    /// Leave each droppable degenerate ring out and report it; refuse
    /// every other fault as [`DegenerateRings::Refuse`] does.
    DropReported,
}

/// Options for [`crate::overlay_with`] and [`crate::union_soup_with`].
///
/// The default is the behaviour of [`crate::overlay`] and
/// [`crate::union_soup`]: every degenerate ring is refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub struct OverlayOptions {
    /// What to do with a degenerate input ring.
    pub degenerate_rings: DegenerateRings,
}

impl OverlayOptions {
    /// The default options: degenerate rings are refused.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            degenerate_rings: DegenerateRings::Refuse,
        }
    }

    /// These options with degenerate rings dropped and reported
    /// ([`DegenerateRings::DropReported`]).
    #[must_use]
    pub const fn dropping_degenerate_rings(mut self) -> Self {
        self.degenerate_rings = DegenerateRings::DropReported;
        self
    }
}

/// Which operand a dropped ring came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum OverlayOperand {
    /// The subject, or the rings of a [`crate::union_soup_with`].
    Subject,
    /// The clip.
    Clip,
}

/// Which ring of its polygon a dropped ring was.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RingRole {
    /// The outer ring.
    Outer,
    /// The hole at this index in [`Polygon::holes`].
    Hole(usize),
}

/// Why a ring was dropped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DropReason {
    /// Fewer than three points ([`OverlayError::RingTooShort`]).
    TooShort,
    /// An edge no longer than the tolerance
    /// ([`OverlayError::RepeatedVertex`]).
    RepeatedVertex,
    /// An area of at most `tolerance.linear()^2`, decided exactly
    /// ([`OverlayError::ZeroArea`]).
    ZeroArea,
    /// A hole, valid on its own, whose outer ring was dropped.
    OuterDropped,
}

/// One input ring an operation left out, and why.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct DroppedRing {
    /// The operand the ring belonged to.
    pub operand: OverlayOperand,
    /// The index of its polygon in the operand; for
    /// [`crate::union_soup_with`], the index of the ring.
    pub polygon: usize,
    /// Which ring of that polygon it was.
    pub ring: RingRole,
    /// Why it was dropped.
    pub reason: DropReason,
    /// The magnitude of its signed area, computed exactly and rounded once
    /// for output. For a dropped outer ring, the area the operand lost; for
    /// a dropped hole, the area it gained.
    pub area: f64,
    /// A proven upper bound on the exact magnitude of its signed area.
    pub area_bound: f64,
}

/// An overlay's result and the rings it dropped (#256).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ReportedOverlay {
    /// The exact overlay of the rings that were kept. Its evidence counts
    /// those rings, not the dropped ones.
    pub result: OverlayResult,
    /// Every dropped ring, subject before clip, then by polygon, outer
    /// ring before holes. Empty when nothing was dropped.
    pub dropped: Vec<DroppedRing>,
}

/// A soup union's polygons and the rings it dropped (#256).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ReportedUnion {
    /// The union of the rings that were kept.
    pub polygons: Vec<Polygon>,
    /// Every dropped ring, in input order, as [`OverlayOperand::Subject`]
    /// with [`RingRole::Outer`].
    pub dropped: Vec<DroppedRing>,
}

/// Where a ring sits, for its report.
#[derive(Clone, Copy)]
pub(crate) struct Place {
    pub(crate) operand: OverlayOperand,
    pub(crate) polygon: usize,
    pub(crate) ring: RingRole,
}

/// The operand's polygons as the boolean should see them: refused as
/// [`crate::overlay`] refuses them, or with each droppable degenerate ring
/// left out and pushed onto `dropped`.
pub(crate) fn admit<'a>(
    input: &'a OverlayInput,
    operand: OverlayOperand,
    tolerance: Tolerance,
    options: OverlayOptions,
    dropped: &mut Vec<DroppedRing>,
) -> Result<Cow<'a, [Polygon]>, OverlayError> {
    if options.degenerate_rings == DegenerateRings::Refuse {
        crate::validate(input, tolerance)?;
        return Ok(Cow::Borrowed(&input.polygons));
    }
    validate_frame(input, tolerance)?;
    let mut kept = Vec::with_capacity(input.polygons.len());
    for (index, polygon) in input.polygons.iter().enumerate() {
        let place = |ring| Place {
            operand,
            polygon: index,
            ring,
        };
        let outer = checked(&polygon.outer, tolerance)?;
        if let Some(reason) = outer {
            dropped.push(report(&polygon.outer, place(RingRole::Outer), reason));
        }
        let mut holes = Vec::with_capacity(polygon.holes.len());
        for (h, hole) in polygon.holes.iter().enumerate() {
            let own = checked(hole, tolerance)?;
            let reason = match (own, outer) {
                (Some(reason), _) => reason,
                (None, Some(_)) => DropReason::OuterDropped,
                (None, None) => {
                    if hole_outside(&polygon.outer, hole, tolerance) {
                        return Err(OverlayError::HoleOutsideOuter);
                    }
                    holes.push(hole.clone());
                    continue;
                }
            };
            dropped.push(report(hole, place(RingRole::Hole(h)), reason));
        }
        if outer.is_none() {
            kept.push(Polygon {
                outer: polygon.outer.clone(),
                holes,
            });
        }
    }
    Ok(Cow::Owned(kept))
}

/// `None` for a valid ring, the reason for a droppable degenerate one,
/// and the error for any other fault.
pub(crate) fn checked(
    ring: &Ring,
    tolerance: Tolerance,
) -> Result<Option<DropReason>, OverlayError> {
    let error = match validate_ring(ring, tolerance) {
        Ok(()) => return Ok(None),
        Err(error) => error,
    };
    // No area can be stated for a ring with a non-finite point, and
    // `validate_ring` tests the point count first.
    if !ring.points.iter().all(|p| p.is_finite()) {
        return Err(OverlayError::NonFinitePoint);
    }
    match error {
        OverlayError::RingTooShort => Ok(Some(DropReason::TooShort)),
        OverlayError::ZeroArea => Ok(Some(DropReason::ZeroArea)),
        OverlayError::RepeatedVertex => {
            let merged = Ring {
                points: merge_short_edges(&ring.points, tolerance.linear()),
            };
            if merged.points.len() >= 3 && self_intersects(&merged, tolerance) {
                Err(OverlayError::SelfIntersection)
            } else {
                Ok(Some(DropReason::RepeatedVertex))
            }
        }
        other => Err(other),
    }
}

/// The ring with every vertex within `linear` of the last one kept left
/// out, and trailing vertices within `linear` of the first: no edge of
/// the result is that short, so the self-intersection test sees the ring's
/// shape rather than the short edges, whose neighbours always touch.
fn merge_short_edges(points: &[Point2], linear: f64) -> Vec<Point2> {
    let mut merged: Vec<Point2> = Vec::with_capacity(points.len());
    for &p in points {
        if merged.last().is_none_or(|&q| (p - q).length() > linear) {
            merged.push(p);
        }
    }
    while merged.len() > 1 && (merged[merged.len() - 1] - merged[0]).length() <= linear {
        merged.pop();
    }
    merged
}

pub(crate) fn report(ring: &Ring, place: Place, reason: DropReason) -> DroppedRing {
    let (area, area_bound) = exact_area(&ring.points);
    DroppedRing {
        operand: place.operand,
        polygon: place.polygon,
        ring: place.ring,
        reason,
        area,
        area_bound,
    }
}

/// The magnitude of a ring's signed area rounded once, and a proven upper
/// bound on it. The points must be finite.
fn exact_area(points: &[Point2]) -> (f64, f64) {
    let Some(o) = points.first() else {
        return (0.0, 0.0);
    };
    let f = Dyadic::from_f64;
    let local = |p: &Point2| (f(p.x).sub(&f(o.x)), f(p.y).sub(&f(o.y)));
    let mut twice = Dyadic::zero();
    for w in points.windows(2).skip(1) {
        let ((ax, ay), (bx, by)) = (local(&w[0]), local(&w[1]));
        twice = twice.add(&ax.mul(&by).sub(&ay.mul(&bx)));
    }
    if twice.sign() == Some(Sign::Negative) {
        twice = twice.neg();
    }
    let enclosure = twice.enclosure();
    let half_up = |v: f64| {
        let h = v * 0.5;
        // Halving is exact unless it leaves the normal range.
        if h * 2.0 == v {
            h
        } else {
            h.next_up()
        }
    };
    let bound = half_up(enclosure.hi());
    (twice.to_f64() * 0.5, bound)
}
