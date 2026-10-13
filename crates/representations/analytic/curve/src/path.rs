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
//! # Open ends (#298)
//!
//! A caller cannot always state where a span ends: the length of an offset
//! beside a gradient curve, a clothoid or a B-spline, and of a B-spline
//! itself, has no closed form. A piece may therefore run to its curve's
//! end without stating it ([`PathPiece::to_end`], [`PathPiece::whole`]):
//! its [`open_end`](PathPiece::open_end) is set and its `end` holds
//! `Scalar::INFINITY`, which nothing reads. The evaluator resolves the end
//! with its own length of the curve, the number it would measure the
//! whole curve by: exact where the length is read from stored data, its
//! numerical estimate otherwise, so the joints after it fall where its
//! certified joint distances say. On an unbounded line the end is its
//! parameter domain's, `origin + direction`, as inside a composite. A
//! span from the curve's start needs no open start, since the station
//! measure starts there at `0` ([`PathPiece::from_start`]). Until an
//! evaluator has resolved it, an open piece has no length here:
//! [`PathPiece::length`] and [`CurvePath::length`] are infinite.
//!
//! # Bounds at a parameter or a point (#304)
//!
//! A trim of an ellipse or a B-spline at an interior native parameter, or
//! at a point, has a station-measure distance that is itself a
//! quadrature, so a caller without a provider cannot state it either. A
//! piece's start and end may therefore be given as a [`PathBound`]
//! ([`PathPiece::between_parameters`], [`PathPiece::bounded`],
//! [`PathPiece::with_start`], [`PathPiece::with_end`]): a native
//! parameter of the curve, or a point on it in the curve's own
//! coordinates (before the piece's placement; a 2D curve's points lie in
//! `z = 0`). The bound is held in [`start_bound`](PathPiece::start_bound)
//! or [`end_bound`](PathPiece::end_bound), and the stated `start` or `end`
//! beside it holds an infinity, which nothing reads. The evaluator
//! resolves a point to the one parameter naming it, refusing a point off
//! the curve or one several parameters name, and a parameter to its
//! distance in the station measure: exact where the curve's arc length to
//! it is read from stored data, the same numerical estimate a piece
//! between parameters is measured by otherwise. As both ends lie in the
//! curve's direction, a reversed piece is bounded as a forward one is.
//! Until an evaluator has resolved them, a bounded piece has no length
//! here, as an open one has none.
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
//! [`PathCurve`], [`PathPiece`], [`PathBound`], [`PathOffset`] and
//! [`OffsetLaw`] are `#[non_exhaustive]`: a later piece kind is a new variant an evaluator
//! that does not know it refuses by name.

use axiolid_core::{Point3, Scalar, Transform3, Vec3};

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

/// Where a [`PathPiece`] starts or ends when its distance is not stated
/// (#304): the evaluator resolves it to a distance in the curve's station
/// measure; see the
/// [module documentation](self#bounds-at-a-parameter-or-a-point-304).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PathBound {
    /// At the curve's native parameter (a closed conic's angle, a
    /// B-spline's knot-domain parameter, a polyline's segment index plus
    /// fraction, a line's `t` along its direction).
    Parameter(Scalar),
    /// At the point of the curve, in the curve's own coordinates (before
    /// the piece's placement; `z = 0` on a 2D curve). It must lie on the
    /// curve, and only one parameter may name it.
    Point(Point3),
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
    /// Where the span ends; a measurable piece has `end > start`. Not read
    /// where the end is open ([`Self::open_end`]).
    pub end: Scalar,
    /// Whether the path runs from `end` to `start`.
    pub reversed: bool,
    /// The rigid motion carrying the curve, if it is placed.
    pub placement: Option<Transform3>,
    /// Whether that placement is exact, rounding aside (`true` with no
    /// placement).
    pub placement_exact: bool,
    /// Whether the span runs to the curve's end, which the evaluator
    /// resolves with its own length (#298; see the
    /// [module documentation](self#open-ends-298)); `end` is then not read.
    pub open_end: bool,
    /// Where the span starts when that is given as a parameter or a point
    /// rather than a distance (#304; see the
    /// [module documentation](self#bounds-at-a-parameter-or-a-point-304));
    /// `start` is then not read.
    pub start_bound: Option<PathBound>,
    /// Where the span ends when that is given as a parameter or a point
    /// (#304); `end` is then not read. An end is bounded or open, not
    /// both: an evaluator refuses a piece whose end is both.
    pub end_bound: Option<PathBound>,
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
            open_end: false,
            start_bound: None,
            end_bound: None,
        }
    }

    /// The span of `curve` from `start` to its end, forwards and unplaced,
    /// the end left to the evaluator (#298; see the
    /// [module documentation](self#open-ends-298)).
    ///
    /// Nothing is checked here: an evaluator refuses a start that is not
    /// finite, or not before the curve's end, when it reads the path.
    #[must_use]
    pub fn to_end(curve: PathCurve, start: Scalar) -> Self {
        Self {
            open_end: true,
            ..Self::new(curve, start, Scalar::INFINITY)
        }
    }

    /// The span of `curve` from its start to `end`: `[0, end]`, the station
    /// measure starting at `0` on every curve.
    #[must_use]
    pub fn from_start(curve: PathCurve, end: Scalar) -> Self {
        Self::new(curve, 0.0, end)
    }

    /// The whole of `curve`, from its start to its end, no length stated
    /// (#298): [`Self::to_end`] from `0`. An unbounded line is its parameter
    /// domain, from its origin to `origin + direction`.
    #[must_use]
    pub fn whole(curve: PathCurve) -> Self {
        Self::to_end(curve, 0.0)
    }

    /// The span of `curve` from its native parameter `t0` to `t1`, `t0 <
    /// t1` (except on a closed conic, below), forwards and unplaced, no distance stated (#304; see the
    /// [module documentation](self#bounds-at-a-parameter-or-a-point-304)).
    /// A closed conic's angle runs round: its span runs from `t0` forwards
    /// to `t1` within one turn, across its parameter seam where `t1` comes
    /// before `t0`, as a trim of it in its sense does.
    ///
    /// Nothing is checked here: an evaluator refuses parameters that are
    /// not finite, not increasing or outside the curve when it reads the
    /// path.
    #[must_use]
    pub fn between_parameters(curve: PathCurve, t0: Scalar, t1: Scalar) -> Self {
        Self::bounded(curve, PathBound::Parameter(t0), PathBound::Parameter(t1))
    }

    /// The span of `curve` from `start` to `end`, each a parameter or a
    /// point, forwards and unplaced, no distance stated (#304). Both lie
    /// in the curve's direction; on a closed conic the span runs from the
    /// start forwards to the end within one turn.
    ///
    /// Nothing is checked here: an evaluator refuses a point off the
    /// curve or named by several parameters, and an empty span, when it
    /// reads the path.
    #[must_use]
    pub fn bounded(curve: PathCurve, start: PathBound, end: PathBound) -> Self {
        Self::new(curve, 0.0, 0.0).with_start(start).with_end(end)
    }

    /// The same piece starting at `bound` rather than at a stated distance
    /// (#304): `start` then holds `-Scalar::INFINITY`, which nothing reads.
    #[must_use]
    pub fn with_start(mut self, bound: PathBound) -> Self {
        self.start = -Scalar::INFINITY;
        self.start_bound = Some(bound);
        self
    }

    /// The same piece ending at `bound` rather than at a stated distance
    /// or at its curve's end (#304): `end` then holds `Scalar::INFINITY`,
    /// which nothing reads, and the end is no longer open.
    #[must_use]
    pub fn with_end(mut self, bound: PathBound) -> Self {
        self.end = Scalar::INFINITY;
        self.open_end = false;
        self.end_bound = Some(bound);
        self
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

    /// The span's length in its curve's station measure, `end - start`;
    /// infinite where the end is open or an end is bounded by a parameter
    /// or a point, which only an evaluator resolves.
    #[must_use]
    pub fn length(&self) -> Scalar {
        if !self.ends_stated() {
            return Scalar::INFINITY;
        }
        self.end - self.start
    }

    /// Whether both of its own ends are stated distances: neither open
    /// (#298) nor bounded by a parameter or a point (#304).
    fn ends_stated(&self) -> bool {
        !self.open_end && self.start_bound.is_none() && self.end_bound.is_none()
    }

    /// Whether every length the piece spans is stated: its own ends (not
    /// open, #298, nor bounded by a parameter or a point, #304), and its
    /// offset's base piece's.
    #[must_use]
    pub fn is_stated(&self) -> bool {
        self.ends_stated()
            && match &self.curve {
                PathCurve::Offset(offset) => offset.base.is_stated(),
                _ => true,
            }
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
    /// path's length once an evaluator has accepted it; infinite where a
    /// piece's end is open (#298), which only an evaluator resolves.
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

    #[test]
    fn an_open_end_is_left_to_the_evaluator_and_survives_composition() {
        let line = PathCurve::Three(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }));
        let open = PathPiece::to_end(line.clone(), 2.0);
        assert!(open.open_end && !open.is_stated());
        assert_eq!((open.start, open.end), (2.0, Scalar::INFINITY));
        assert_eq!(open.length(), Scalar::INFINITY);
        // An open piece equals itself (no NaN sentinel).
        assert_eq!(open, open.clone());
        let whole = PathPiece::whole(line.clone());
        assert!(whole.open_end && whole.start == 0.0);
        let from_start = PathPiece::from_start(line.clone(), 3.0);
        assert!(!from_start.open_end && from_start.is_stated());
        assert_eq!((from_start.start, from_start.length()), (0.0, 3.0));
        // Reversing and placing keep the open end.
        let path = CurvePath::new(vec![from_start, open])
            .reversed()
            .placed(Transform3::from_rotation_z(0.5), true);
        assert!(path.pieces()[0].open_end && path.pieces()[0].reversed);
        assert!(!path.pieces()[1].open_end);
        assert_eq!(path.length(), Scalar::INFINITY);
        // A stated offset over an open base is not stated.
        let offset = PathCurve::Offset(Box::new(PathOffset::new(
            whole,
            OffsetLaw::Planar { distance: 1.0 },
        )));
        let piece = PathPiece::new(offset, 0.0, 1.0);
        assert!(!piece.open_end && !piece.is_stated());
    }

    #[test]
    fn a_bound_is_left_to_the_evaluator_and_survives_composition() {
        let line = PathCurve::Three(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }));
        let between = PathPiece::between_parameters(line.clone(), 0.5, 2.0);
        assert_eq!(between.start_bound, Some(PathBound::Parameter(0.5)));
        assert_eq!(between.end_bound, Some(PathBound::Parameter(2.0)));
        assert!(!between.open_end && !between.is_stated());
        assert_eq!(
            (between.start, between.end),
            (-Scalar::INFINITY, Scalar::INFINITY)
        );
        assert_eq!(between.length(), Scalar::INFINITY);
        assert_eq!(between, between.clone());
        // One bound is enough to leave the length unstated.
        let point = PathBound::Point(Point3::new(3.0, 0.0, 0.0));
        let start = PathPiece::new(line.clone(), 0.0, 4.0).with_start(point);
        assert!(!start.is_stated() && start.end_bound.is_none());
        assert_eq!((start.end, start.length()), (4.0, Scalar::INFINITY));
        let end = PathPiece::new(line.clone(), 0.0, 4.0).with_end(point);
        assert!(!end.is_stated() && end.start_bound.is_none());
        assert_eq!((end.start, end.length()), (0.0, Scalar::INFINITY));
        // A bound end closes an open one; an open end keeps a bound start.
        let closed = PathPiece::whole(line.clone()).with_end(point);
        assert!(!closed.open_end && closed.end_bound == Some(point));
        let open = PathPiece::whole(line.clone()).with_start(point);
        assert!(open.open_end && open.start_bound == Some(point));
        // Reversing and placing keep the bounds, in the curve's direction.
        let path = CurvePath::new(vec![between.clone(), open])
            .reversed()
            .placed(Transform3::from_rotation_z(0.5), true);
        assert_eq!(path.pieces()[1].start_bound, between.start_bound);
        assert_eq!(path.pieces()[1].end_bound, between.end_bound);
        assert!(path.pieces()[1].reversed && path.pieces()[0].open_end);
        assert_eq!(path.length(), Scalar::INFINITY);
        // A stated offset over a bounded base is not stated.
        let offset = PathCurve::Offset(Box::new(PathOffset::new(
            between,
            OffsetLaw::Planar { distance: 1.0 },
        )));
        assert!(!PathPiece::new(offset, 0.0, 1.0).is_stated());
        let bounded = PathPiece::bounded(line, PathBound::Parameter(1.0), point);
        assert_eq!(bounded.start_bound, Some(PathBound::Parameter(1.0)));
        assert_eq!(bounded.end_bound, Some(point));
    }
}
