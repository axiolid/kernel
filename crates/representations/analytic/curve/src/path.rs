//! Pieces of curves laid end to end: the neutral form of a curve relation a
//! distance runs along (#285, #290, ADR 0082).
//!
//! A composite curve, a trim of one, and a run of segments each placed at a
//! station of another curve are graph relations in `axiolid-model`. A
//! reader that measures a distance along them -- a station, a point or a
//! frame at a distance -- does not need the graph: it needs the spans of
//! atomic curves the relation runs along, in order, each read forwards or
//! backwards and carried by a rigid placement. [`CurvePath`] is that list,
//! and [`PathPiece`] one entry of it. It owns its curves: this crate is part
//! of the geometry data plane, whose values carry no borrowed references so
//! that a native backend can copy them across an FFI boundary.
//!
//! This module holds data and its structural composition only (reversing a
//! path, placing it, joining paths). It measures nothing: where a span lies
//! on its curve, whether consecutive pieces meet and which distance runs
//! through them are read by an evaluator (`axiolid-evaluate`'s
//! `CompositeBasis`, through the curve-evaluation contract's `path_*`
//! queries), which refuses a path it cannot measure by name. A graph
//! compiler flattens a relation into a path (`axiolid-mesh-compile`'s
//! `station::curve_path`).
//!
//! # Reading a piece
//!
//! A piece is the span `[start, end]` of its curve in that curve's station
//! measure from its start: plan distance on an elevated or banked curve,
//! arc length on every other curve (a 2D curve lying in `z = 0`). A
//! reversed piece is that span read from `end` to `start`; a placed piece
//! is its curve carried by the rigid motion `placement`. The path's
//! distance runs end to end through its pieces, each in its own measure.
//! The rules an evaluator applies (joints as seams read by their
//! [`SeamSide`](crate::SeamSide), the common convention, exactness) are
//! ADR 0082's.
//!
//! # Exactness
//!
//! A frame read on a piece is exact, rounding aside, only on a line placed,
//! if at all, exactly ([`PathPiece::frame_is_exact`]); on a path, only
//! where every piece up to and including the one read is.
//!
//! # Offsets (#289)
//!
//! A piece may be an offset of another piece ([`PathCurve::Offset`]): the
//! curve beside a base piece at the displacement its [`OffsetLaw`] gives
//! in the base's frame, measured in its OWN length from the base piece's
//! start (its own plan length where the base is plan-measured). The base
//! is one span of an atomic curve, not itself an offset, with no seam of
//! its curve inside it: an offset of a curve with seams is several offset
//! pieces, one per span between seams, which an evaluator joins as it
//! joins any pieces (and refuses where they do not meet). An offset of a
//! line by a constant law is a line ([`PathCurve::is_line`]); the other
//! rules (closed forms, collapse and cusp refusals) are ADR 0082's.
//!
//! # Growth
//!
//! [`PathCurve`], [`PathPiece`], [`PathOffset`] and [`OffsetLaw`] are
//! `#[non_exhaustive]`: a later piece kind is a new variant an evaluator
//! that does not know it refuses by name.

use axiolid_core::{Scalar, Transform3, Vec3};

use crate::{Curve2, Curve3};

/// The curve under a [`PathPiece`].
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum PathCurve {
    /// An atomic 2D curve, lying in `z = 0`.
    Two(Curve2),
    /// An atomic 3D curve.
    Three(Curve3),
    /// An offset of a base piece (#289), measured in its own length; see
    /// the [module documentation](self#offsets-289).
    Offset(Box<PathOffset>),
}

impl PathCurve {
    /// Whether the curve is a line, the one family whose frame is exact
    /// at a distance, rounding aside: an atomic line, or an offset of a
    /// line by a constant law (a line beside it).
    #[must_use]
    pub fn is_line(&self) -> bool {
        match self {
            Self::Two(Curve2::Line(_)) | Self::Three(Curve3::Line(_)) => true,
            Self::Offset(offset) => offset.law.is_constant() && offset.base.curve.is_line(),
            _ => false,
        }
    }

    /// Whether every placement inside the curve (an offset's base
    /// piece's) is exact; `true` for an atomic curve.
    fn inner_placement_exact(&self) -> bool {
        match self {
            Self::Offset(offset) => {
                offset.base.placement_exact && offset.base.curve.inner_placement_exact()
            }
            _ => true,
        }
    }
}

/// The displacement of an offset at one place, in its base's frame:
/// `lateral` to the left, `vertical` along up, `longitudinal` along the
/// tangent (ADR 0082's station offsets).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PathOffsets {
    /// Along the lateral axis, positive to the left.
    pub lateral: Scalar,
    /// Along the section's up.
    pub vertical: Scalar,
    /// Along the tangent.
    pub longitudinal: Scalar,
}

impl PathOffsets {
    /// The displacement `(lateral, vertical, longitudinal)`.
    #[must_use]
    pub const fn new(lateral: Scalar, vertical: Scalar, longitudinal: Scalar) -> Self {
        Self {
            lateral,
            vertical,
            longitudinal,
        }
    }

    /// The displacement a fraction `u` of the way from `self` to `other`.
    #[must_use]
    pub fn lerp(self, other: Self, u: Scalar) -> Self {
        let mix = |a: Scalar, b: Scalar| a + (b - a) * u;
        Self::new(
            mix(self.lateral, other.lateral),
            mix(self.vertical, other.vertical),
            mix(self.longitudinal, other.longitudinal),
        )
    }
}

/// Which frame of the base an [`OffsetLaw::Linear`] displacement is read
/// in: ADR 0082's station frames.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum OffsetFrame {
    /// The base's section frame (planar, reference-up or banked).
    #[default]
    Section,
    /// The upright frame: horizontal tangent, horizontal left normal, `+Z`.
    Plan,
}

/// How far, and which way, an offset lies from its base piece.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OffsetLaw {
    /// A constant offset in the plane of a 2D base: `distance` along its
    /// left normal (an offset curve 2D, positive to the left).
    Planar {
        /// The offset distance, positive to the left.
        distance: Scalar,
    },
    /// A constant 3D offset: `distance` along `normalise(V x T)`, `V` the
    /// fixed `reference_direction` and `T` the base's unit tangent (an
    /// offset curve 3D). `T` must nowhere be parallel to `V`.
    Directed {
        /// The offset distance.
        distance: Scalar,
        /// The fixed reference direction `V`.
        reference_direction: Vec3,
    },
    /// A displacement in the base's `frame`, linear in the base's measure
    /// from `start` at the base piece's start to `end` at its end: one
    /// interval of an offset by distances at stations.
    Linear {
        /// The displacement at the base piece's start.
        start: PathOffsets,
        /// The displacement at the base piece's end.
        end: PathOffsets,
        /// The frame it is read in.
        frame: OffsetFrame,
    },
}

impl OffsetLaw {
    /// Whether the displacement is the same all along the base.
    #[must_use]
    pub fn is_constant(&self) -> bool {
        match self {
            Self::Planar { .. } | Self::Directed { .. } => true,
            Self::Linear { start, end, .. } => start == end,
        }
    }
}

/// An offset of a base piece (#289); see the
/// [module documentation](self#offsets-289).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct PathOffset {
    /// The base piece: a span of an atomic curve, reversed or placed, not
    /// itself an offset, with no seam of its curve inside it.
    pub base: PathPiece,
    /// The displacement from it.
    pub law: OffsetLaw,
}

impl PathOffset {
    /// The offset of `base` by `law`. Nothing is checked here: an
    /// evaluator refuses a base it cannot offset when it reads the path.
    #[must_use]
    pub fn new(base: PathPiece, law: OffsetLaw) -> Self {
        Self { base, law }
    }
}

/// One piece of a [`CurvePath`]: the span `[start, end]` of a curve in its
/// station measure, traversed forwards or backwards, optionally carried by
/// a rigid placement; see the [module documentation](self).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct PathPiece {
    /// The curve.
    pub curve: PathCurve,
    /// Where the span starts, in the curve's station measure from its
    /// start (negative only on a line).
    pub start: Scalar,
    /// Where the span ends; a measurable piece has `end > start`.
    pub end: Scalar,
    /// Whether the path runs from `end` to `start`.
    pub reversed: bool,
    /// The rigid motion carrying the curve, if it is placed.
    pub placement: Option<Transform3>,
    /// Whether that placement is exact, rounding aside (`true` with no
    /// placement).
    pub placement_exact: bool,
}

impl PathPiece {
    /// The span `[start, end]` of `curve`, forwards and unplaced.
    ///
    /// Nothing is checked here: an evaluator refuses a span that is empty,
    /// not finite or outside its curve when it reads the path.
    #[must_use]
    pub fn new(curve: PathCurve, start: Scalar, end: Scalar) -> Self {
        Self {
            curve,
            start,
            end,
            reversed: false,
            placement: None,
            placement_exact: true,
        }
    }

    /// The same span traversed the other way.
    #[must_use]
    pub fn reversed(mut self) -> Self {
        self.reversed = !self.reversed;
        self
    }

    /// The same piece carried further by the rigid motion `rigid`, applied
    /// after any placement it already has; `exact` when that motion is.
    #[must_use]
    pub fn placed(mut self, rigid: Transform3, exact: bool) -> Self {
        self.placement = Some(match self.placement {
            Some(inner) => rigid * inner,
            None => rigid,
        });
        self.placement_exact &= exact;
        self
    }

    /// The span's length in its curve's station measure, `end - start`.
    #[must_use]
    pub fn length(&self) -> Scalar {
        self.end - self.start
    }

    /// Whether a frame read on this piece is exact, rounding aside: a line
    /// (an offset of a line by a constant law included), placed, if at
    /// all, exactly.
    #[must_use]
    pub fn frame_is_exact(&self) -> bool {
        self.curve.is_line() && self.placement_exact && self.curve.inner_placement_exact()
    }
}

/// Pieces of curves laid end to end, in the order the path runs; see the
/// [module documentation](self).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CurvePath {
    pieces: Vec<PathPiece>,
}

impl CurvePath {
    /// The path through `pieces`, in order. Nothing is checked here (see
    /// the [module documentation](self)).
    #[must_use]
    pub fn new(pieces: Vec<PathPiece>) -> Self {
        Self { pieces }
    }

    /// The pieces, in order.
    #[must_use]
    pub fn pieces(&self) -> &[PathPiece] {
        &self.pieces
    }

    /// The pieces, in order, by value.
    #[must_use]
    pub fn into_pieces(self) -> Vec<PathPiece> {
        self.pieces
    }

    /// Whether the path has no pieces (an evaluator refuses it).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pieces.is_empty()
    }

    /// The same path traversed the other way: its pieces in reverse order,
    /// each reversed.
    #[must_use]
    pub fn reversed(self) -> Self {
        self.pieces
            .into_iter()
            .rev()
            .map(PathPiece::reversed)
            .collect()
    }

    /// The same path carried further by the rigid motion `rigid` (every
    /// piece, after any placement it has), `exact` when that motion is.
    #[must_use]
    pub fn placed(self, rigid: Transform3, exact: bool) -> Self {
        self.pieces
            .into_iter()
            .map(|piece| piece.placed(rigid, exact))
            .collect()
    }

    /// The sum of the pieces' lengths, in their station measure: the
    /// path's length once an evaluator has accepted it.
    #[must_use]
    pub fn length(&self) -> Scalar {
        self.pieces.iter().map(PathPiece::length).sum()
    }
}

impl From<PathPiece> for CurvePath {
    fn from(piece: PathPiece) -> Self {
        Self::new(vec![piece])
    }
}

impl FromIterator<PathPiece> for CurvePath {
    fn from_iter<I: IntoIterator<Item = PathPiece>>(pieces: I) -> Self {
        Self::new(pieces.into_iter().collect())
    }
}

impl Extend<PathPiece> for CurvePath {
    fn extend<I: IntoIterator<Item = PathPiece>>(&mut self, pieces: I) {
        self.pieces.extend(pieces);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Line3, Polyline3};
    use axiolid_core::{Point3, Vec3};

    #[test]
    fn reversing_a_path_reverses_its_order_and_each_piece() {
        let line = Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        });
        let corner = Curve3::Polyline(Polyline3 {
            points: vec![Point3::ZERO, Point3::new(1.0, 0.0, 0.0)],
            closed: false,
        });
        let a = PathPiece::new(PathCurve::Three(line), 0.0, 2.0);
        let b = PathPiece::new(PathCurve::Three(corner), 0.0, 1.0).reversed();
        let path = CurvePath::new(vec![a.clone(), b.clone()]).reversed();
        assert_eq!(path.pieces(), &[b.reversed(), a.reversed()]);
        assert!(!path.pieces()[0].reversed && path.pieces()[1].reversed);
        assert_eq!(path.length(), 3.0);
    }

    #[test]
    fn placements_compose_outermost_last_and_exactness_is_sticky() {
        let line = Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        });
        let inner = Transform3::from_translation(Vec3::X);
        let outer = Transform3::from_rotation_z(0.5);
        let piece = PathPiece::new(PathCurve::Three(line.clone()), 0.0, 1.0)
            .placed(inner, true)
            .placed(outer, false);
        assert_eq!(piece.placement, Some(outer * inner));
        assert!(!piece.placement_exact && !piece.frame_is_exact());
        let exact = PathPiece::new(PathCurve::Three(line), 0.0, 1.0).placed(inner, true);
        assert!(exact.frame_is_exact());
    }

    #[test]
    fn an_offset_of_a_line_by_a_constant_law_is_a_line() {
        let line = Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        });
        let directed = OffsetLaw::Directed {
            distance: 1.0,
            reference_direction: Vec3::Z,
        };
        let base = PathPiece::new(PathCurve::Three(line), 0.0, 4.0);
        let constant = PathCurve::Offset(Box::new(PathOffset::new(base.clone(), directed)));
        assert!(constant.is_line());
        let widening = PathCurve::Offset(Box::new(PathOffset::new(
            base.clone(),
            OffsetLaw::Linear {
                start: PathOffsets::new(1.0, 0.0, 0.0),
                end: PathOffsets::new(2.0, 0.0, 0.0),
                frame: OffsetFrame::Section,
            },
        )));
        assert!(!widening.is_line());
        // An inexactly placed base makes the offset's frame inexact.
        let placed = base.placed(Transform3::from_rotation_z(0.5), false);
        let offset = PathPiece::new(
            PathCurve::Offset(Box::new(PathOffset::new(placed, directed))),
            0.0,
            4.0,
        );
        assert!(offset.curve.is_line() && !offset.frame_is_exact());
        let half = PathOffsets::new(0.0, 0.0, 2.0).lerp(PathOffsets::new(2.0, 4.0, 0.0), 0.5);
        assert_eq!(half, PathOffsets::new(1.0, 2.0, 1.0));
    }
}
