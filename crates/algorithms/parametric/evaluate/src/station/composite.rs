//! A composite station basis: stations along pieces of curves laid end to
//! end (#285, ADR 0082 amendment).
//!
//! A [`CompositeBasis`] is a run of [`StationPiece`]s. Each piece is a span
//! of an atomic 2D or 3D curve ([`StationCurve`]) between two distances in
//! that curve's own station measure, traversed forwards or backwards, and
//! optionally carried by a rigid placement; or of an offset of such a span
//! ([`StationCurve::Offset`], #289), in the offset's own length (see
//! [the offset module](super::offset)). This is the form of a composite
//! curve relation this module measures: a format adapter or a graph
//! compiler flattens its trims, sense flags, nested composites and curves
//! placed at stations into pieces, and this module measures stations along
//! them.
//!
//! The same pieces as an owned, neutral value are `axiolid-curve`'s
//! [`CurvePath`] (#290), which the curve-evaluation contract can see:
//! [`CompositeBasis::from_path`] measures a path with its pieces borrowed
//! from it, [`CompositeBasis::path`] hands the pieces out as one, and a
//! [`StationPiece`] converts to a [`PathPiece`] and back
//! ([`StationPiece::from_path_piece`]). The contract's `path_*` queries
//! read a path through the same [`CompositeBasis::section_on`] in the
//! reference provider. A [`StationPiece`] borrows its curve so that a
//! station on a graph relation, resolved again for every station of a run,
//! copies none.
//!
//! # Distance convention
//!
//! A composite measures stations as its pieces do, end to end: the
//! distance of a point is the sum of the lengths of the pieces before it
//! plus its distance into its own piece, each in its piece's station
//! measure. A piece measures in its atomic curve's convention ([the parent
//! module](super)): PLAN distance on an elevated or banked curve, arc
//! length on every other curve (a 2D curve's arc length is its plan
//! length). A rigid placement keeps arc length; it keeps plan distance only
//! when it keeps `+Z`, and a plan-measured piece placed by one that tilts
//! `+Z` is refused by name ([`SectionFrame::carried`]). So:
//!
//! - every piece plan-measured (elevated or banked curves, trimmed,
//!   reversed or placed upright): the composite is measured in plan
//!   distance, [`DistanceConvention::PlanDistance`];
//! - every piece arc-length-measured (lines, conics, polylines,
//!   B-splines, intrinsic curves and chains, 2D or 3D, placed or not): the
//!   composite is measured in arc length,
//!   [`DistanceConvention::ArcLength3d`]. A composite of segments placed at
//!   stations of an elevated curve is measured in ITS OWN pieces' arc
//!   length, whatever the base the segments were placed on;
//! - pieces of both kinds: refused by name, since no one distance runs
//!   through it.
//!
//! A piece with no length is refused by name, as is an unbounded one: an
//! untrimmed line inside a composite is its parameter domain `[0, 1]`, as
//! a composite directrix reads it ([`StationPiece::whole`]).
//!
//! # Joints
//!
//! Consecutive pieces must meet: the end of one within [`JOINT_TOLERANCE`]
//! `* max(1, |p|)` of the start of the next, `|p|` the largest coordinate
//! magnitude of the joint. A joint whose pieces do not meet is refused by
//! name: as a reversed piece where the next piece's END meets the previous
//! one (or the first piece's start meets the second), which is a piece
//! whose sense was not declared; as a gap otherwise.
//!
//! Every interior joint is a seam under the #263 rule
//! ([the seam module](super::seam)): it is never smooth, since a
//! declared continuity is not a guarantee in the data. A station within
//! [`ARC_LENGTH_TOLERANCE`]` * max(1, s)` of a joint is ON it and reads the
//! piece its [`SeamSide`] names at the joint's own distance:
//! [`SeamSide::Outgoing`] the piece that starts there, at its start;
//! [`SeamSide::Incoming`] the piece that ends there, at its end. Each side
//! reads its own piece's point: two pieces that meet within the joint
//! tolerance are not snapped onto one another. At the composite's start
//! and end only one piece exists, and both sides read it; likewise a piece
//! is always read from the inside at its own ends, so a trim that starts
//! on a seam of its curve never reads the part of the curve it trims away.
//! The seams of a piece's curve strictly inside the piece are the
//! composite's seams too, read by the curve's own rule.
//!
//! # Reversed and placed pieces
//!
//! A reversed piece is its span read from its end: the point at distance
//! `u` into the piece is the curve's at `end - u`, the side of a seam is
//! swapped (the incoming piece of the composite is the curve's outgoing
//! one), and the frame's tangent and lateral are negated while up is kept,
//! which is the section frame of the curve traversed the other way (the
//! reference-up frame and the planar frame are exactly that; a banked
//! curve's rolled section keeps its up and keeps the left rail on the
//! left).
//!
//! A placed piece's frame is its curve's carried by the placement
//! ([`SectionFrame::carried`]): moved rigidly where the placement keeps
//! `+Z`; where it tilts `+Z`, the reference-up frame of the placed curve's
//! own point and tangent, since the moved frame would not be the placed 3D
//! curve's own section frame. A plan-measured curve placed tilted is
//! refused by name.
//!
//! # Exactness
//!
//! A station frame on a composite is exact, rounding aside, only where
//! every piece up to and including the one it is read on is exact
//! ([`CompositeBasis::frame_is_exact_at`]): those pieces' lengths make the
//! distance, and the read piece makes the frame. A piece is exact when its
//! curve is a line and its placement, if any, was exact
//! ([`StationPiece::frame_is_exact`]). Elsewhere the frame carries the
//! tolerances of the parent module.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Point3, Scalar, Transform3};
use axiolid_curve::{Curve2, Curve3, CurvePath, PathCurve, PathOffset, PathPiece, SeamSide};
use axiolid_curve_evaluate_contract::DistanceConvention;

use super::offset::{OffsetTable, StationOffset};
use super::seam::tidy;
use super::{
    admitted, invalid, start2, start3, station_length2, station_length3, station_seams2,
    station_seams3, station_section2_on, station_section3_on, SectionFrame, StationSeam,
};
use crate::arc_parameter::{
    arc_length2, arc_length3, parameter_at_arc_length2, parameter_at_arc_length3,
    ARC_LENGTH_TOLERANCE,
};
use crate::curve::{evaluate2, evaluate3};

/// Relative tolerance within which two consecutive pieces of a
/// [`CompositeBasis`] meet: their joint points may be `JOINT_TOLERANCE *
/// max(1, |p|)` apart, `|p|` the joint's largest coordinate magnitude.
pub const JOINT_TOLERANCE: Scalar = 1e-9;

/// The arc-length tolerance around a distance: what counts as ON a seam or
/// a joint.
fn slack(distance: Scalar) -> Scalar {
    ARC_LENGTH_TOLERANCE * distance.abs().max(1.0)
}

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BackendId::new("axiolid-evaluate"),
        operation: Operation::CurveEvaluation,
        input,
    }
}

/// The atomic curve under a [`StationPiece`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StationCurve<'c> {
    /// A 2D curve, lying in `z = 0`.
    Two(&'c Curve2),
    /// A 3D curve.
    Three(&'c Curve3),
    /// An offset of a base piece, measured in its own length (#289; see
    /// [the offset module](super::offset)).
    Offset(StationOffset<'c>),
}

impl From<StationCurve<'_>> for PathCurve {
    /// The curve, copied.
    fn from(curve: StationCurve<'_>) -> Self {
        match curve {
            StationCurve::Two(curve) => PathCurve::Two(curve.clone()),
            StationCurve::Three(curve) => PathCurve::Three(curve.clone()),
            StationCurve::Offset(offset) => PathCurve::Offset(Box::new(PathOffset::new(
                PathPiece::from(offset.base()),
                offset.law(),
            ))),
        }
    }
}

impl<'c> StationCurve<'c> {
    /// The station curve of a path piece's curve, borrowed from it.
    ///
    /// # Errors
    ///
    /// A piece kind this evaluator does not know, by name.
    pub fn from_path_curve(curve: &'c PathCurve) -> GeomResult<Self> {
        match curve {
            PathCurve::Two(curve) => Ok(Self::Two(curve)),
            PathCurve::Three(curve) => Ok(Self::Three(curve)),
            PathCurve::Offset(offset) => Ok(Self::Offset(StationOffset::new(
                StationPiece::from_path_piece(&offset.base)?,
                offset.law,
            )?)),
            _ => Err(unsupported(
                "a curve path piece of a kind this evaluator does not know",
            )),
        }
    }
}

impl StationCurve<'_> {
    /// The distance convention a station on this curve is measured in:
    /// plan distance on an elevated or banked curve, arc length on every
    /// other (a 2D curve's arc length is its plan length).
    /// An offset is measured as its base is (#289).
    #[must_use]
    pub fn convention(self) -> DistanceConvention {
        match self {
            Self::Three(Curve3::Elevated(_) | Curve3::Banked(_)) => {
                DistanceConvention::PlanDistance
            }
            Self::Offset(offset) => offset.convention(),
            _ => DistanceConvention::ArcLength3d,
        }
    }

    /// The station length (the parent module's `station_length`).
    fn length(self) -> GeomResult<Option<Scalar>> {
        match self {
            Self::Two(curve) => station_length2(curve),
            Self::Three(curve) => station_length3(curve),
            Self::Offset(offset) => offset.length().map(Some),
        }
    }

    /// Whether a frame read on it is exact: a line, or an offset of a line
    /// by a constant law.
    fn frame_is_exact(self) -> bool {
        match self {
            Self::Offset(offset) => offset.frame_is_exact(),
            _ => self.is_line(),
        }
    }

    fn is_line(self) -> bool {
        matches!(
            self,
            Self::Two(Curve2::Line(_)) | Self::Three(Curve3::Line(_))
        )
    }

    /// A circle or an ellipse: its measure runs round, so a span may cross
    /// its parameter seam.
    fn is_closed_conic(self) -> bool {
        matches!(
            self,
            Self::Two(Curve2::Circle(_) | Curve2::Ellipse(_))
                | Self::Three(Curve3::Circle(_) | Curve3::Ellipse(_))
        )
    }

    /// Whether a distance on it is read from stored data alone: false for
    /// an ellipse and a B-spline, whose arc length is a quadrature, and an
    /// offset other than a line or a circle.
    fn measure_is_exact(self) -> bool {
        match self {
            Self::Offset(offset) => offset.measure_is_exact(),
            _ => !matches!(
                self,
                Self::Two(Curve2::Ellipse(_) | Curve2::BSpline(_))
                    | Self::Three(Curve3::Ellipse(_) | Curve3::BSpline(_))
            ),
        }
    }

    /// A line's direction length, `None` for any other family.
    fn line_speed(self) -> Option<Scalar> {
        match self {
            Self::Two(Curve2::Line(line)) => Some(line.direction.length()),
            Self::Three(Curve3::Line(line)) => Some(line.direction.length()),
            _ => None,
        }
    }

    /// A closed conic's measure `m` folded into one turn `[0, L]`.
    fn folded(self, m: Scalar) -> GeomResult<Scalar> {
        if !self.is_closed_conic() {
            return Ok(m);
        }
        let turn = self.length()?.unwrap_or(Scalar::INFINITY);
        Ok(if m > turn { m - turn } else { m })
    }

    /// The station measure at native parameter `t`, from the curve's
    /// start: `t |d|` on a line (any sign), `r t` on a circle, the arc
    /// length from the domain's start on a polyline, a B-spline or an
    /// ellipse, `t` itself on a curve parameterised by its measure.
    fn measure_at(self, t: Scalar) -> GeomResult<Scalar> {
        Ok(match self {
            Self::Two(Curve2::Line(line)) => t * line.direction.length(),
            Self::Three(Curve3::Line(line)) => t * line.direction.length(),
            Self::Two(Curve2::Circle(circle)) => t * circle.radius,
            Self::Three(Curve3::Circle(circle)) => t * circle.radius,
            Self::Two(curve @ (Curve2::Polyline(_) | Curve2::BSpline(_) | Curve2::Ellipse(_))) => {
                arc_length2(curve, start2(curve), t)?
            }
            Self::Three(
                curve @ (Curve3::Polyline(_) | Curve3::BSpline(_) | Curve3::Ellipse(_)),
            ) => arc_length3(curve, start3(curve), t)?,
            Self::Offset(_) => {
                return Err(unsupported(
                    "the native parameter of an offset piece: an offset takes its basis's \
                     parameter, not its own length",
                ))
            }
            _ => t,
        })
    }

    /// The point at measure `m` from the curve's start; an offset's
    /// through `table`.
    fn point_at(self, m: Scalar, table: Option<&OffsetTable>) -> GeomResult<Point3> {
        let m = self.folded(m)?;
        match self {
            Self::Offset(offset) => offset.point_at(m, table),
            Self::Two(curve) => {
                let t = match curve {
                    Curve2::Line(line) => m / line.direction.length(),
                    _ => parameter_at_arc_length2(curve, start2(curve), m)?,
                };
                let p = evaluate2(curve, t)?;
                Ok(Point3::new(p.x, p.y, 0.0))
            }
            Self::Three(curve) => match curve {
                Curve3::Elevated(elevated) => crate::arc_length::elevated_point(elevated, m),
                Curve3::Banked(banked) => crate::banked::banked_point(banked, m),
                Curve3::Line(line) => evaluate3(curve, m / line.direction.length()),
                _ => evaluate3(curve, parameter_at_arc_length3(curve, start3(curve), m)?),
            },
        }
    }

    /// The section frame at measure `m`, read from `side` on a seam; an
    /// offset's through `table`.
    fn section_at(
        self,
        m: Scalar,
        side: SeamSide,
        table: Option<&OffsetTable>,
    ) -> GeomResult<SectionFrame> {
        let m = self.folded(m)?;
        match self {
            Self::Offset(offset) => offset.section_at(m, table),
            // A line inside a composite may be trimmed before its origin, so
            // its measure may be negative: read it directly.
            Self::Two(Curve2::Line(line)) => {
                let point = line.origin + line.direction * (m / line.direction.length());
                super::planar_frame(point, line.direction)
            }
            Self::Three(curve @ Curve3::Line(line)) => super::reference_up_frame(
                evaluate3(curve, m / line.direction.length())?,
                line.direction,
            ),
            Self::Two(curve) => station_section2_on(curve, m, side),
            Self::Three(curve) => station_section3_on(curve, m, side),
        }
    }

    fn seams(self) -> GeomResult<Vec<StationSeam>> {
        if self.is_closed_conic() {
            return Ok(Vec::new());
        }
        match self {
            Self::Two(curve) => station_seams2(curve),
            Self::Three(curve) => station_seams3(curve),
            // Its base spans no seam of its curve.
            Self::Offset(_) => Ok(Vec::new()),
        }
    }
}

/// One piece of a [`CompositeBasis`]: the span `[start, end]` of an atomic
/// curve, in that curve's station measure, traversed forwards or
/// backwards, optionally placed by a rigid motion; see the
/// [module documentation](self).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StationPiece<'c> {
    /// The atomic curve.
    pub curve: StationCurve<'c>,
    /// Where the span starts, in the curve's station measure from its
    /// start (negative only on a line).
    pub start: Scalar,
    /// Where the span ends, greater than `start`.
    pub end: Scalar,
    /// Whether the composite runs from `end` to `start`.
    pub reversed: bool,
    /// The rigid motion carrying the curve, if it is placed.
    pub placement: Option<Transform3>,
    /// Whether that placement is exact, rounding aside (`true` with no
    /// placement).
    pub placement_exact: bool,
}

impl<'c> StationPiece<'c> {
    /// The whole curve: `[0, L]`, `L` its station length; an unbounded
    /// line is its parameter domain `[0, 1]`, from its origin to `origin +
    /// direction`.
    ///
    /// # Errors
    ///
    /// A curve whose length a station cannot be measured against, by name.
    pub fn whole(curve: StationCurve<'c>) -> GeomResult<Self> {
        let end = match curve.line_speed() {
            Some(speed) => speed,
            None => curve.length()?.ok_or_else(|| {
                invalid("station: a composite piece has no bounded length".into())
            })?,
        };
        Self::between(curve, 0.0, end)
    }

    /// The span `[start, end]` of `curve`, in its station measure.
    ///
    /// # Errors
    ///
    /// A non-finite or empty span, and one outside the curve's measure
    /// (before its start or past its length, beyond the arc-length
    /// tolerance; a closed conic's span may cross its parameter seam but
    /// not exceed one turn), by name.
    pub fn between(curve: StationCurve<'c>, start: Scalar, end: Scalar) -> GeomResult<Self> {
        if !(start.is_finite() && end.is_finite()) || end <= start {
            return Err(GeomError::Degenerate(format!(
                "station: a composite piece spans [{start}, {end}], which is empty or not finite"
            )));
        }
        if curve
            .line_speed()
            .is_some_and(|speed| !(speed.is_finite() && speed > 0.0))
        {
            return Err(GeomError::Degenerate(
                "station: a composite piece is a line with no direction".into(),
            ));
        }
        if !curve.is_line() {
            let length = curve.length()?.unwrap_or(Scalar::INFINITY);
            let reach = if curve.is_closed_conic() {
                start + length
            } else {
                length
            };
            if start < -slack(0.0) || start > length + slack(length) || end > reach + slack(reach) {
                return Err(invalid(format!(
                    "station: a composite piece spans [{start}, {end}], outside its curve's \
                     measure [0, {length}]"
                )));
            }
        }
        Ok(Self {
            curve,
            start,
            end,
            reversed: false,
            placement: None,
            placement_exact: true,
        })
    }

    /// The span between native parameters `t0 < t1` of `curve`, converted
    /// to its station measure. A closed conic's span may start anywhere
    /// and cross its parameter seam (at most one turn).
    ///
    /// # Errors
    ///
    /// As [`Self::between`], and the arc-length quadrature's refusals.
    pub fn between_parameters(curve: StationCurve<'c>, t0: Scalar, t1: Scalar) -> GeomResult<Self> {
        if !(t0.is_finite() && t1.is_finite()) || t1 <= t0 {
            return Err(GeomError::Degenerate(format!(
                "station: a composite piece spans parameters [{t0}, {t1}], which is empty or \
                 not finite"
            )));
        }
        if curve.is_closed_conic() {
            let turn = core::f64::consts::TAU;
            let from = t0.rem_euclid(turn);
            let to = from + (t1 - t0);
            let start = curve.measure_at(from)?;
            let end = if to <= turn {
                curve.measure_at(to)?
            } else {
                let length = curve.length()?.unwrap_or(Scalar::INFINITY);
                length + curve.measure_at(to - turn)?
            };
            return Self::between(curve, start, end);
        }
        Self::between(curve, curve.measure_at(t0)?, curve.measure_at(t1)?)
    }

    /// The same span traversed the other way.
    #[must_use]
    pub const fn reversed(mut self) -> Self {
        self.reversed = !self.reversed;
        self
    }

    /// The same piece carried further by the rigid motion `rigid` (applied
    /// after any placement it already has), `exact` when that motion is
    /// (as [`PathPiece::placed`]).
    #[must_use]
    pub fn placed(mut self, rigid: Transform3, exact: bool) -> Self {
        self.placement = Some(match self.placement {
            Some(inner) => rigid * inner,
            None => rigid,
        });
        self.placement_exact &= exact;
        self
    }

    /// The whole of an offset (#289): `[0, L]`, `L` its own length.
    ///
    /// # Errors
    ///
    /// The offset's refusals (a collapse, a cusp, a self-crossing), by
    /// name.
    pub fn offset(offset: StationOffset<'c>) -> GeomResult<Self> {
        let end = offset.length()?;
        if !(end.is_finite() && end > 0.0) {
            return Err(GeomError::Degenerate(
                "station: an offset piece has no length".into(),
            ));
        }
        Ok(Self {
            curve: StationCurve::Offset(offset),
            start: 0.0,
            end,
            reversed: false,
            placement: None,
            placement_exact: true,
        })
    }

    /// The offset's measure table, for an offset piece.
    fn table(&self) -> GeomResult<Option<OffsetTable>> {
        match self.curve {
            StationCurve::Offset(offset) => offset.table().map(Some),
            _ => Ok(None),
        }
    }

    /// The station piece of a path piece: its span checked against its
    /// curve as [`Self::between`] checks it, then reversed and placed as
    /// the path piece is.
    ///
    /// # Errors
    ///
    /// As [`Self::between`], and a piece kind this evaluator does not
    /// know, by name.
    pub fn from_path_piece(piece: &'c PathPiece) -> GeomResult<Self> {
        let curve = StationCurve::from_path_curve(&piece.curve)?;
        Ok(Self {
            reversed: piece.reversed,
            placement: piece.placement,
            placement_exact: piece.placement_exact,
            ..Self::between(curve, piece.start, piece.end)?
        })
    }

    /// The piece's length in its station measure.
    #[must_use]
    pub fn length(&self) -> Scalar {
        self.end - self.start
    }

    /// The distance convention the piece is measured in.
    #[must_use]
    pub fn convention(&self) -> DistanceConvention {
        self.curve.convention()
    }

    /// Whether a frame read on this piece is exact, rounding aside: a line
    /// (an offset of a line by a constant law included), placed, if at
    /// all, exactly (as [`PathPiece::frame_is_exact`]).
    #[must_use]
    pub fn frame_is_exact(&self) -> bool {
        self.curve.frame_is_exact() && self.placement_exact
    }

    /// The curve's measure at distance `u` into the piece, and the side of
    /// a seam of the curve that `side` (in the composite's direction) is.
    fn measure(&self, u: Scalar, side: SeamSide) -> (Scalar, SeamSide) {
        if self.reversed {
            let swapped = match side {
                SeamSide::Incoming => SeamSide::Outgoing,
                _ => SeamSide::Incoming,
            };
            (self.end - u, swapped)
        } else {
            (self.start + u, side)
        }
    }

    /// The point at distance `u` into the piece; an offset's through
    /// `table`.
    fn point_at(&self, u: Scalar, table: Option<&OffsetTable>) -> GeomResult<Point3> {
        let (m, _) = self.measure(u, SeamSide::Outgoing);
        let point = self.curve.point_at(m, table)?;
        Ok(match self.placement {
            Some(rigid) => rigid.transform_point3(point),
            None => point,
        })
    }

    /// The section frame at distance `u` into the piece, read from `side`
    /// on a seam of its curve; at the piece's own ends it is always read
    /// from inside the piece.
    ///
    /// # Errors
    ///
    /// The curve evaluators' refusals, and a plan-measured curve placed by
    /// a motion that tilts `+Z`, by name.
    pub fn section_on(&self, u: Scalar, side: SeamSide) -> GeomResult<SectionFrame> {
        self.section_with(u, side, None)
    }

    /// [`Self::section_on`], an offset read through `table`.
    fn section_with(
        &self,
        u: Scalar,
        side: SeamSide,
        table: Option<&OffsetTable>,
    ) -> GeomResult<SectionFrame> {
        let length = self.length();
        let u = u.clamp(0.0, length);
        let side = if u <= slack(u) {
            SeamSide::Outgoing
        } else if u >= length - slack(length) {
            SeamSide::Incoming
        } else {
            side
        };
        let (m, curve_side) = self.measure(u, side);
        let mut section = self.curve.section_at(m, curve_side, table)?;
        if self.reversed {
            section.tangent = -section.tangent;
            section.lateral = -section.lateral;
        }
        match self.placement {
            Some(rigid) => section.carried(rigid, self.convention()),
            None => Ok(section),
        }
    }

    /// The seams of the curve strictly inside the piece, at their distance
    /// into the piece in the composite's direction.
    pub(super) fn seams(&self) -> GeomResult<Vec<StationSeam>> {
        let mut out: Vec<StationSeam> = self
            .curve
            .seams()?
            .into_iter()
            .filter(|seam| {
                seam.distance > self.start + slack(self.start)
                    && seam.distance < self.end - slack(self.end)
            })
            .map(|seam| {
                let u = if self.reversed {
                    self.end - seam.distance
                } else {
                    seam.distance - self.start
                };
                StationSeam::new(u, u, seam.smooth, seam.exact)
            })
            .collect();
        out.sort_by(|a, b| a.distance.total_cmp(&b.distance));
        Ok(out)
    }
}

impl From<StationPiece<'_>> for PathPiece {
    /// The piece, its curve copied.
    fn from(piece: StationPiece<'_>) -> Self {
        let mut out = PathPiece::new(piece.curve.into(), piece.start, piece.end);
        out.reversed = piece.reversed;
        out.placement = piece.placement;
        out.placement_exact = piece.placement_exact;
        out
    }
}

/// Pieces of curves laid end to end as one station basis; see the
/// [module documentation](self).
#[derive(Debug, Clone, PartialEq)]
pub struct CompositeBasis<'c> {
    pieces: Vec<StationPiece<'c>>,
    /// Each offset piece's measure table, built once (#289).
    tables: Vec<Option<OffsetTable>>,
    /// Each piece's start distance, and the total length last.
    starts: Vec<Scalar>,
    convention: DistanceConvention,
}

impl<'c> CompositeBasis<'c> {
    /// The composite of `pieces`, in order.
    ///
    /// # Errors
    ///
    /// No pieces, a piece with no length, pieces measured in different
    /// conventions, and a joint whose pieces do not meet (a reversed piece
    /// or a gap), each by name; and the evaluators' refusals reading the
    /// joints' points.
    pub fn new(pieces: Vec<StationPiece<'c>>) -> GeomResult<Self> {
        let Some(first) = pieces.first() else {
            return Err(invalid("station: a composite basis has no pieces".into()));
        };
        let convention = first.convention();
        for (index, piece) in pieces.iter().enumerate() {
            let length = piece.length();
            if !length.is_finite() || length <= 0.0 {
                return Err(GeomError::Degenerate(format!(
                    "station: piece {index} of a composite basis has no length"
                )));
            }
            if piece.convention() != convention {
                return Err(unsupported(
                    "a composite station basis whose pieces measure stations differently (plan \
                     distance on an elevated or banked piece, arc length on another): no one \
                     distance runs through it",
                ));
            }
        }
        let tables = pieces
            .iter()
            .map(StationPiece::table)
            .collect::<GeomResult<Vec<_>>>()?;
        for index in 1..pieces.len() {
            check_joint(
                (&pieces[index - 1], tables[index - 1].as_ref()),
                (&pieces[index], tables[index].as_ref()),
                index,
            )?;
        }
        let mut starts = Vec::with_capacity(pieces.len() + 1);
        let mut run = 0.0;
        for piece in &pieces {
            starts.push(run);
            run += piece.length();
        }
        starts.push(run);
        Ok(Self {
            pieces,
            tables,
            starts,
            convention,
        })
    }

    /// The composite along `path` (#290), its curves borrowed from it: each
    /// piece checked against its curve ([`StationPiece::from_path_piece`]),
    /// then measured as [`Self::new`] measures them.
    ///
    /// # Errors
    ///
    /// As [`StationPiece::from_path_piece`] for each piece, naming its
    /// index, and as [`Self::new`].
    pub fn from_path(path: &'c CurvePath) -> GeomResult<Self> {
        let pieces = path
            .pieces()
            .iter()
            .map(StationPiece::from_path_piece)
            .collect::<GeomResult<Vec<_>>>()?;
        Self::new(pieces)
    }

    /// The composite's pieces as a neutral [`CurvePath`] (#290), their
    /// curves copied.
    #[must_use]
    pub fn path(&self) -> CurvePath {
        self.pieces
            .iter()
            .map(|piece| PathPiece::from(*piece))
            .collect()
    }

    /// The pieces, in order.
    #[must_use]
    pub fn pieces(&self) -> &[StationPiece<'c>] {
        &self.pieces
    }

    /// The total length in the composite's station measure.
    #[must_use]
    pub fn length(&self) -> Scalar {
        self.starts[self.pieces.len()]
    }

    /// The convention every piece, and so the composite, is measured in.
    #[must_use]
    pub fn convention(&self) -> DistanceConvention {
        self.convention
    }

    /// The piece a station at the admitted `distance` reads on `side`, and
    /// the distance into it.
    fn locate(&self, distance: Scalar, side: SeamSide) -> (usize, Scalar) {
        let count = self.pieces.len();
        for joint in 1..count {
            let at = self.starts[joint];
            if (at - distance).abs() <= slack(distance) {
                return match side {
                    SeamSide::Incoming => (joint - 1, self.pieces[joint - 1].length()),
                    _ => (joint, 0.0),
                };
            }
        }
        let k = self.starts[1..count].partition_point(|start| *start <= distance);
        (k, distance - self.starts[k])
    }

    /// The piece a station at `distance` reads on `side` (the one
    /// [`Self::section_on`] reads).
    ///
    /// # Errors
    ///
    /// As [`Self::section_on`] for the distance.
    pub(crate) fn piece_read(
        &self,
        distance: Scalar,
        side: SeamSide,
    ) -> GeomResult<&StationPiece<'c>> {
        let distance = admitted(distance, Some(self.length()))?;
        Ok(&self.pieces[self.locate(distance, side).0])
    }

    /// The section frame at `distance` along the composite, read from
    /// `side` on a joint or a seam inside a piece; see the
    /// [module documentation](self#joints).
    ///
    /// # Errors
    ///
    /// A distance that is not finite, negative or beyond the composite's
    /// length (by more than the arc-length tolerance), and the piece's
    /// refusals, by name.
    pub fn section_on(&self, distance: Scalar, side: SeamSide) -> GeomResult<SectionFrame> {
        let distance = admitted(distance, Some(self.length()))?;
        let (piece, into) = self.locate(distance, side);
        self.pieces[piece].section_with(into, side, self.tables[piece].as_ref())
    }

    /// Whether the frame at `distance`, read from `side`, is exact: every
    /// piece up to and including the one it is read on is
    /// ([`StationPiece::frame_is_exact`]). `false` for a distance outside
    /// the composite.
    #[must_use]
    pub fn frame_is_exact_at(&self, distance: Scalar, side: SeamSide) -> bool {
        let Ok(distance) = admitted(distance, Some(self.length())) else {
            return false;
        };
        let (piece, _) = self.locate(distance, side);
        self.pieces[..=piece]
            .iter()
            .all(StationPiece::frame_is_exact)
    }

    /// The composite's seams: every interior joint (never smooth) and the
    /// seams of each piece's curve strictly inside the piece, at their
    /// distance along the composite, which is also the parameter reported.
    /// A seam's distance is exact when every piece before it, and its own
    /// seam, is measured from stored data alone (not an ellipse or a
    /// B-spline).
    ///
    /// # Errors
    ///
    /// The seam readers' refusals.
    pub fn seams(&self) -> GeomResult<Vec<StationSeam>> {
        let mut out = Vec::new();
        let mut exact_before = true;
        for (k, piece) in self.pieces.iter().enumerate() {
            let at = self.starts[k];
            if k > 0 {
                out.push(StationSeam::new(at, at, false, exact_before));
            }
            for seam in piece.seams()? {
                let distance = at + seam.distance;
                out.push(StationSeam::new(
                    distance,
                    distance,
                    seam.smooth,
                    seam.exact && exact_before,
                ));
            }
            exact_before &= piece.curve.measure_is_exact();
        }
        Ok(tidy(out, Some(self.length())))
    }

    /// [`Self::seams`], every one exact or refused.
    ///
    /// # Errors
    ///
    /// As [`Self::seams`], and [`GeomError::UnsupportedInput`] for a seam
    /// whose distance needs a quadrature: one after an ellipse or a
    /// B-spline piece, or a B-spline's corner knot.
    pub fn exact_seams(&self) -> GeomResult<Vec<StationSeam>> {
        let seams = self.seams()?;
        if seams.iter().any(|seam| !seam.exact) {
            return Err(unsupported(
                "exact seam distances of a composite with an ellipse or a B-spline piece before \
                 a seam, or a B-spline's corner knot: their arc length is a quadrature",
            ));
        }
        Ok(seams)
    }

    /// The pieces between distances `start < end` along the composite,
    /// clipped: what a trim of the composite by its station measure keeps.
    ///
    /// # Errors
    ///
    /// A non-finite or empty interval, or one outside `[0, L]` beyond the
    /// arc-length tolerance, by name.
    pub fn pieces_between(&self, start: Scalar, end: Scalar) -> GeomResult<Vec<StationPiece<'c>>> {
        let length = self.length();
        if !(start.is_finite() && end.is_finite()) || end <= start {
            return Err(GeomError::Degenerate(format!(
                "station: a trim of a composite basis spans [{start}, {end}], which is empty or \
                 not finite"
            )));
        }
        if start < -slack(0.0) || end > length + slack(length) {
            return Err(invalid(format!(
                "station: a trim of a composite basis spans [{start}, {end}], outside its \
                 length {length}"
            )));
        }
        let (start, end) = (start.max(0.0), end.min(length));
        let mut out = Vec::new();
        for (k, piece) in self.pieces.iter().enumerate() {
            let (lo, hi) = (self.starts[k], self.starts[k + 1]);
            let (a, b) = (start.max(lo) - lo, end.min(hi) - lo);
            if b - a <= slack(b) {
                continue;
            }
            let mut clipped = *piece;
            if piece.reversed {
                clipped.start = piece.end - b;
                clipped.end = piece.end - a;
            } else {
                clipped.start = piece.start + a;
                clipped.end = piece.start + b;
            }
            out.push(clipped);
        }
        Ok(out)
    }
}

/// A piece with its offset table, if any.
type Read<'a, 'c> = (&'a StationPiece<'c>, Option<&'a OffsetTable>);

/// Refuse a joint whose pieces do not meet, naming an offset across a
/// corner, a reversed piece or the gap.
fn check_joint(
    (before, before_table): Read<'_, '_>,
    (after, after_table): Read<'_, '_>,
    index: usize,
) -> GeomResult<()> {
    let end = before.point_at(before.length(), before_table)?;
    let start = after.point_at(0.0, after_table)?;
    let tolerance = |p: Point3| JOINT_TOLERANCE * p.abs().max_element().max(1.0);
    let meets = |a: Point3, b: Point3| (a - b).length() <= tolerance(a);
    if meets(end, start) {
        return Ok(());
    }
    if [before, after]
        .iter()
        .any(|piece| matches!(piece.curve, StationCurve::Offset(_)))
    {
        return Err(invalid(format!(
            "station: an offset across a corner of its basis at joint {index}: the offsets of \
             the two sides do not meet (a gap of {}), so no distance runs across it (the \
             basis's tangent or roll jumps there)",
            (end - start).length()
        )));
    }
    let after_end = after.point_at(after.length(), after_table)?;
    let reversed = meets(end, after_end)
        || (index == 1 && {
            let before_start = before.point_at(0.0, before_table)?;
            meets(before_start, start) || meets(before_start, after_end)
        });
    if reversed {
        return Err(invalid(format!(
            "station: a reversed piece at joint {index} of a composite basis: its pieces meet \
             end to end or start to start, so one runs against the composite (its sense is not \
             declared)"
        )));
    }
    Err(invalid(format!(
        "station: a gap of {} at joint {index} of a composite basis: its pieces do not meet, so \
         no distance runs across it",
        (end - start).length()
    )))
}
