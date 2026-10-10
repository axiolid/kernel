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
//! # Growth
//!
//! [`PathCurve`] and [`PathPiece`] are `#[non_exhaustive]`: a later piece
//! kind (an offset of a curve, #289) is a new variant an evaluator that
//! does not know it refuses by name.

use axiolid_core::{Scalar, Transform3};

use crate::{Curve2, Curve3};

/// The curve under a [`PathPiece`].
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum PathCurve {
    /// An atomic 2D curve, lying in `z = 0`.
    Two(Curve2),
    /// An atomic 3D curve.
    Three(Curve3),
}

impl PathCurve {
    /// Whether the curve is a line, the one family whose frame is exact
    /// at a distance, rounding aside.
    #[must_use]
    pub fn is_line(&self) -> bool {
        matches!(
            self,
            Self::Two(Curve2::Line(_)) | Self::Three(Curve3::Line(_))
        )
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

    /// Whether a frame read on this piece is exact, rounding aside: a line,
    /// placed, if at all, exactly.
    #[must_use]
    pub fn frame_is_exact(&self) -> bool {
        self.curve.is_line() && self.placement_exact
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
}
