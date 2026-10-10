//! Stations: a position at a distance along a curve, framed by the curve's
//! section (#241, ADR 0082).
//!
//! A station names a place by a measure along a basis curve instead of by
//! coordinates, the way an alignment is staked out. Reading it back is
//! curve evaluation, so it lives here, beside the arc-length machinery it
//! reuses ([`crate::arc_parameter`]).
//!
//! # Distance convention
//!
//! The distance is measured in the basis curve's own convention, from its
//! start:
//!
//! - an [`Elevated3`](axiolid_curve::Elevated3) or a
//!   [`Banked3`](axiolid_curve::Banked3) is authored against PLAN distance,
//!   its own parameter, and a station on it is a plan distance (the
//!   curve-evaluation contract's `DistanceConvention::PlanDistance`);
//! - every other curve, 2D or 3D, is measured by its arc length (for a 2D
//!   curve, the plan length).
//!
//! The start is parameter `0` for a line, a circle, an ellipse, an
//! intrinsic curve and a chain, and the first parameter of the domain for
//! a polyline or a B-spline. A distance must lie in `[0, L]`, `L` the
//! curve's length ([`station_length2`], [`station_length3`]): one full turn
//! for a closed conic, no limit for a line. A negative or non-finite
//! distance, and one past `L` by more than
//! [`ARC_LENGTH_TOLERANCE`]` * max(1, L)`, is refused by name; one past it
//! by less is read at `L`.
//!
//! # Section frame
//!
//! [`SectionFrame`] carries the point and three unit axes: `tangent` along
//! the curve, `lateral` to the LEFT of it, and `up = tangent x lateral`.
//!
//! - A 2D curve lies in the plane `z = 0`: `lateral` is the left normal,
//!   `up` is `+Z`.
//! - A [`Banked3`](axiolid_curve::Banked3) is framed by its rolled section
//!   ([`crate::banked::banked_section`]): `lateral` towards the left rail
//!   head.
//! - Every other 3D curve, an [`Elevated3`](axiolid_curve::Elevated3)
//!   included, gets the reference-up frame against `+Z`, the one
//!   `ReferenceCurveEvaluator::frame_at` returns: `lateral` is horizontal,
//!   `up` leans back with the grade. A vertical tangent has no such frame
//!   and is refused by name.
//!
//! [`SectionFrame::frame`] returns the axes in the curve-evaluation
//! provider's layout: `x` the tangent, `y` up, `z = -lateral` (to the
//! right), the layout the contract's `frame_at` documents (#242).
//!
//! The curve-evaluation provider reads a seam side through this module:
//! `ReferenceCurveEvaluator`'s `point_at_on`, `tangent_at_on` and
//! `frame_at_on` answer a measure on a seam (below) with
//! [`station_section3_on`] (#286).
//!
//! [`SectionFrame::plan`] is the vertical alternative: the tangent's
//! horizontal projection, the horizontal left normal and `+Z`, with grade
//! and bank dropped.
//!
//! # Explicit orientation (#246)
//!
//! [`SectionFrame::oriented`] turns a section frame by an explicit axis and
//! reference direction given as components in that frame, `(tangent,
//! lateral, up)`, not in world coordinates: a linear placement's axes "are
//! relative to the curve used for linear referencing ..., maintaining the
//! relationship to the tangent of the curve" (buildingSMART IFC 4.3,
//! `IfcAxis2PlacementLinear`), whose local `X`, `Y`, `Z` are the tangent,
//! the left lateral and up. The axis is the oriented up, exact; the
//! reference direction is made perpendicular to it (Gram-Schmidt, axis
//! primary) and becomes the oriented tangent. The point does not move.
//!
//! # Seams (#263)
//!
//! A curve made of pieces has two frames where a polyline turns at a
//! vertex, a B-spline at a corner knot, an elevated curve's grade breaks
//! or a banked curve's cant jumps. A station within
//! [`ARC_LENGTH_TOLERANCE`]` * max(1, s)` of such a seam is ON it: it is
//! read at the seam's own distance, from the piece that starts there
//! (`SeamSide::Outgoing`, [`station_section2`] and [`station_section3`])
//! or from the side [`station_section2_on`] and [`station_section3_on`]
//! are given. [`station_seams2`] and [`station_seams3`] say where the
//! seams are without evaluating the curve, and [`Mitre`] is the plane a
//! run of sections crossing one stands in; see [`seam`] for both. At the
//! curve's start and end there is one piece, and both sides read it.
//!
//! # Placing at a station (#264)
//!
//! [`SectionFrame::placement`] is the rigid motion a node placed at a
//! station is moved by: local `x`, `y`, `z` onto the tangent, the left
//! lateral and up, the origin onto the point (a linear placement's
//! reading). [`station_frame_is_exact2`] and [`station_frame_is_exact3`]
//! say whether that frame is exact (a line) or carries the tolerances
//! below. [`SectionFrame::carried`] is the frame a station on the placed
//! curve has (#285): the source's moved where the placement keeps `+Z`,
//! the placed curve's own reference-up frame where it tilts it.
//!
//! # Composite bases (#285)
//!
//! [`CompositeBasis`] measures stations along pieces of atomic curves laid
//! end to end ([`StationPiece`]: a span of a curve, reversed or placed),
//! the neutral form of a composite curve relation. The distance runs
//! through the pieces in their common convention, every joint is a seam
//! read by the rule above, and a frame is exact only where every piece up
//! to the one read is; see [`composite`].
//!
//! # Offset bases (#289)
//!
//! A piece may be an offset of a base piece ([`StationOffset`]), measured
//! in its OWN length from its start: arc length, or its own plan length
//! beside an elevated or banked curve. [`offset_pieces`] offsets pieces
//! laid end to end, one offset piece per span between seams; see
//! [`offset`] for the closed forms, the numerical reading and the
//! refusals.
//!
//! # Accuracy contract
//!
//! The point and frame are the exact curve's at a measure within
//! [`ARC_LENGTH_TOLERANCE`]` * max(1, s)` of `s` where the arc length is
//! numerical (a conic, a polyline is exact, a B-spline), at `s` itself on a
//! line, an intrinsic curve, a chain and an elevated or banked curve (their
//! parameter is the measure), each evaluated to its own stated accuracy
//! (Gauss-Legendre quadrature for intrinsic positions, ADR 0060). A station
//! on a seam is read at the seam's distance, within the same tolerance of
//! its own.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame3, Point3, Scalar, Transform3, Vec3};
use axiolid_curve::{Curve2, Curve3, SeamSide};
/// The convention a station basis is measured in (#285): the
/// curve-evaluation contract's, re-exported so a caller can name it.
pub use axiolid_curve_evaluate_contract::DistanceConvention;

pub mod composite;
pub mod offset;
pub mod seam;

pub use composite::{CompositeBasis, StationCurve, StationPiece, JOINT_TOLERANCE};
pub use offset::{offset_pieces, StationOffset, CUSP_TOLERANCE, OFFSET_TOLERANCE};
pub use seam::{
    exact_station_seams2, exact_station_seams3, station_seams2, station_seams3, Mitre, StationSeam,
    MITRE_TOLERANCE, SEAM_TANGENT_TOLERANCE,
};

use crate::arc_parameter::{
    arc_length2, arc_length3, parameter_at_arc_length2, parameter_at_arc_length3,
    ARC_LENGTH_TOLERANCE,
};
use crate::curve::{derivative2, derivative3, domain2, domain3, evaluate2, evaluate3};

/// Largest sine of the angle between an orientation's axis and reference
/// direction that still counts as parallel ([`SectionFrame::oriented`]);
/// the graph refuses the same pairs when a station is pushed.
pub const ORIENTATION_TOLERANCE: Scalar = 1e-9;

/// How far a rigid placement's image of `+Z` may lean and still count as
/// keeping `+Z` ([`SectionFrame::carried`]): rounding only.
pub const KEEPS_UP_TOLERANCE: Scalar = 1e-12;

fn invalid(detail: String) -> GeomError {
    GeomError::InvalidInput(detail)
}

/// A point on a curve and the section axes there: tangent, left lateral
/// and up, unit and right-handed (`up = tangent x lateral`).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SectionFrame {
    /// The curve point.
    pub point: Point3,
    /// Unit tangent, in the direction of increasing distance.
    pub tangent: Vec3,
    /// Unit lateral axis, to the left of the tangent.
    pub lateral: Vec3,
    /// Unit section up, `tangent x lateral`.
    pub up: Vec3,
}

impl SectionFrame {
    /// The point offset from this frame's point: `lateral` along the
    /// lateral axis (positive to the left), `vertical` along up and
    /// `longitudinal` along the tangent.
    #[must_use]
    pub fn place(&self, lateral: Scalar, vertical: Scalar, longitudinal: Scalar) -> Point3 {
        self.point + lateral * self.lateral + vertical * self.up + longitudinal * self.tangent
    }

    /// The same point framed vertically: the tangent's horizontal
    /// projection, the horizontal left normal and `+Z`. Grade and bank are
    /// dropped.
    ///
    /// # Errors
    ///
    /// A vertical tangent, which has no horizontal projection.
    pub fn plan(&self) -> GeomResult<Self> {
        let horizontal = Vec3::new(self.tangent.x, self.tangent.y, 0.0);
        let length = horizontal.length();
        if !length.is_finite() || length <= 1e-12 {
            return Err(invalid(
                "station: the tangent is vertical, so it has no plan direction".into(),
            ));
        }
        let tangent = horizontal / length;
        Ok(Self {
            point: self.point,
            tangent,
            lateral: Vec3::new(-tangent.y, tangent.x, 0.0),
            up: Vec3::Z,
        })
    }

    /// The same point with its axes turned by an explicit orientation
    /// (#246): `axis` and `ref_direction` are components in THIS frame,
    /// `(tangent, lateral, up)`, defaulting to `(0, 0, 1)` and `(1, 0, 0)`.
    ///
    /// Gram-Schmidt with the axis primary: `up' = axis / |axis|`,
    /// `tangent' = normalise(r - (r . up') up')` for the unit reference
    /// direction `r`, `lateral' = up' x tangent'`. This is the reading of a
    /// linear placement's `Axis` and `RefDirection` relative to the curve
    /// (see the [module documentation](self#explicit-orientation-246)).
    ///
    /// # Errors
    ///
    /// A zero or non-finite vector, or an axis and reference direction
    /// parallel or anti-parallel within [`ORIENTATION_TOLERANCE`], by name.
    pub fn oriented(&self, axis: Option<Vec3>, ref_direction: Option<Vec3>) -> GeomResult<Self> {
        let unit = |vector: Vec3, what: &str| {
            let length = vector.length();
            if vector.is_finite() && length.is_finite() && length > 1e-12 {
                Ok(vector / length)
            } else {
                Err(invalid(format!(
                    "station: the orientation's {what} {vector:?} is zero or not finite"
                )))
            }
        };
        let axis = unit(axis.unwrap_or(Vec3::Z), "axis")?;
        let reference = unit(ref_direction.unwrap_or(Vec3::X), "reference direction")?;
        if axis.cross(reference).length() <= ORIENTATION_TOLERANCE {
            return Err(invalid(format!(
                "station: the orientation's axis {axis:?} and reference direction \
                 {reference:?} are parallel"
            )));
        }
        let x = (reference - reference.dot(axis) * axis).normalize();
        let y = axis.cross(x);
        let world = |v: Vec3| v.x * self.tangent + v.y * self.lateral + v.z * self.up;
        Ok(Self {
            point: self.point,
            tangent: world(x),
            lateral: world(y),
            up: world(axis),
        })
    }

    /// The rigid placement this frame stands for (#264): local `x` onto the
    /// tangent, `y` onto the lateral (to the left), `z` onto up, the local
    /// origin onto the point.
    ///
    /// This is how a linear placement maps its local axes (buildingSMART
    /// IFC 4.3, `IfcAxis2PlacementLinear`: `X` the tangent, `Y` the left
    /// lateral, `Z` up), and how a node placed at a station is moved. It is
    /// not [`Self::frame`]'s provider layout (`y` up, `z` right).
    #[must_use]
    pub fn placement(&self) -> Transform3 {
        Transform3::from_cols(self.tangent, self.lateral, self.up, self.point)
    }

    /// This frame carried by a rigid motion: the point by `rigid`, the
    /// axes by its linear part.
    ///
    /// A section frame on a curve moved rigidly is the moved curve's own
    /// only where its convention does not depend on the world: always for
    /// its tangent, and for the reference-up and plan frames only when
    /// `rigid` keeps `+Z`; checking that is the caller's.
    #[must_use]
    pub fn moved(&self, rigid: Transform3) -> Self {
        Self {
            point: rigid.transform_point3(self.point),
            tangent: rigid.transform_vector3(self.tangent),
            lateral: rigid.transform_vector3(self.lateral),
            up: rigid.transform_vector3(self.up),
        }
    }

    /// This frame, read on a curve, as the section frame of that curve
    /// carried by the rigid motion `rigid` (#264, #285): the frame a
    /// station on the placed curve has, `measure` the convention the
    /// source curve is measured in.
    ///
    /// Where `rigid` keeps `+Z` (within [`KEEPS_UP_TOLERANCE`]) this is
    /// [`Self::moved`]: the moved reference-up, planar or banked frame is
    /// the placed curve's own, and plan distance is kept. Where it tilts
    /// `+Z`, an arc-length-measured curve (arc length is kept by any rigid
    /// motion) gets the reference-up frame against `+Z` of its moved point
    /// and tangent, which is the placed 3D curve's own section frame; the
    /// moved frame would not be.
    ///
    /// # Errors
    ///
    /// A plan-measured curve (elevated or banked) under a motion that
    /// tilts `+Z`, whose plan distance and section frame are not the placed
    /// curve's, by a typed [`GeomError::UnsupportedInput`]; and a moved
    /// tangent with no reference-up frame (vertical), by name.
    pub fn carried(&self, rigid: Transform3, measure: DistanceConvention) -> GeomResult<Self> {
        if (rigid.transform_vector3(Vec3::Z) - Vec3::Z).length() <= KEEPS_UP_TOLERANCE {
            return Ok(self.moved(rigid));
        }
        if measure == DistanceConvention::PlanDistance {
            return Err(GeomError::UnsupportedInput {
                backend: axiolid_contracts::BackendId::new("axiolid-evaluate"),
                operation: axiolid_contracts::Operation::CurveEvaluation,
                input: PLAN_MEASURED_TILTED,
            });
        }
        reference_up_frame(
            rigid.transform_point3(self.point),
            rigid.transform_vector3(self.tangent),
        )
    }

    /// The frame in the curve-evaluation provider's layout: `x` the
    /// tangent, `y` up, `z = x x y`, which is `-lateral` (to the right).
    #[must_use]
    pub fn frame(&self) -> Frame3 {
        Frame3 {
            origin: self.point,
            x: self.tangent,
            y: self.up,
            z: -self.lateral,
        }
    }
}

/// The refusal of a plan-measured curve placed in a tilted frame.
const PLAN_MEASURED_TILTED: &str =
    "a station along an elevated or banked curve placed at a station whose frame tilts +Z: its \
     plan distance and its section frame carried by the placement are not the placed curve's own \
     (#264)";

/// A 2D curve's point and unit tangent lifted to `z = 0`, framed with the
/// left normal and `+Z`.
fn planar_frame(
    point: axiolid_core::Point2,
    tangent: axiolid_core::Vec2,
) -> GeomResult<SectionFrame> {
    let length = tangent.length();
    if !length.is_finite() || length <= 0.0 {
        return Err(invalid(
            "station: the curve has no tangent direction there".into(),
        ));
    }
    let t = tangent / length;
    Ok(SectionFrame {
        point: Point3::new(point.x, point.y, 0.0),
        tangent: Vec3::new(t.x, t.y, 0.0),
        lateral: Vec3::new(-t.y, t.x, 0.0),
        up: Vec3::Z,
    })
}

/// The reference-up frame against `+Z`, as `ReferenceCurveEvaluator`
/// builds it: `right = t x Z`, `up = right x t`, `lateral = -right`.
fn reference_up_frame(point: Point3, derivative: Vec3) -> GeomResult<SectionFrame> {
    let length = derivative.length();
    if !length.is_finite() || length <= 0.0 {
        return Err(invalid(
            "station: the curve has no tangent direction there".into(),
        ));
    }
    let tangent = derivative / length;
    let right = tangent.cross(Vec3::Z);
    let magnitude = right.length();
    if !magnitude.is_finite() || magnitude <= 1e-12 {
        return Err(invalid(
            "station: the tangent is vertical, so the section has no lateral axis".into(),
        ));
    }
    let right = right / magnitude;
    Ok(SectionFrame {
        point,
        tangent,
        lateral: -right,
        up: right.cross(tangent),
    })
}

/// A banked curve's rolled section at plan distance `d`, as a station
/// frame.
fn banked_frame(curve: &axiolid_curve::Banked3, d: Scalar) -> GeomResult<SectionFrame> {
    let section = crate::banked::banked_section(curve, d)?;
    Ok(SectionFrame {
        point: section.point,
        tangent: section.tangent,
        lateral: section.lateral,
        up: section.up,
    })
}

/// Check `distance` against the curve length `length` (`None`:
/// unbounded) and return the distance to read, clamped onto `L` within
/// tolerance.
fn admitted(distance: Scalar, length: Option<Scalar>) -> GeomResult<Scalar> {
    if !distance.is_finite() {
        return Err(invalid(format!(
            "station: distance {distance} along the curve is not finite"
        )));
    }
    if distance < 0.0 {
        return Err(invalid(format!(
            "station: distance {distance} is before the curve's start"
        )));
    }
    match length {
        None => Ok(distance),
        Some(length) => {
            let slack = ARC_LENGTH_TOLERANCE * length.abs().max(1.0);
            if distance > length + slack {
                Err(invalid(format!(
                    "station: distance {distance} is beyond the curve's length {length}"
                )))
            } else {
                Ok(distance.min(length))
            }
        }
    }
}

fn unmeasured(dimension: &str) -> GeomError {
    invalid(format!(
        "station: this {dimension} curve family has no distance a station can be measured by"
    ))
}

/// Where a 2D curve's measure starts: parameter `0`, or the domain's
/// first parameter for a polyline or a B-spline.
fn start2(curve: &Curve2) -> Scalar {
    match curve {
        Curve2::Polyline(_) | Curve2::BSpline(_) => {
            let domain = domain2(curve);
            domain.start.min(domain.end)
        }
        _ => 0.0,
    }
}

fn start3(curve: &Curve3) -> Scalar {
    match curve {
        Curve3::Polyline(_) | Curve3::BSpline(_) => {
            let domain = domain3(curve);
            domain.start.min(domain.end)
        }
        _ => 0.0,
    }
}

/// The length a station on a 2D curve is measured against: `None` for a
/// line (unbounded), one turn for a circle or an ellipse, the arc length
/// of its domain for a polyline or a B-spline, the stated length of an
/// intrinsic curve or a chain.
///
/// # Errors
///
/// Any other family, a chain without a valid length, and the quadrature's
/// refusals, by name.
pub fn station_length2(curve: &Curve2) -> GeomResult<Option<Scalar>> {
    match curve {
        Curve2::Line(_) => Ok(None),
        Curve2::Circle(_) | Curve2::Ellipse(_) => {
            arc_length2(curve, 0.0, core::f64::consts::TAU).map(Some)
        }
        Curve2::Polyline(_) | Curve2::BSpline(_) => {
            let domain = domain2(curve);
            let (lo, hi) = (domain.start.min(domain.end), domain.start.max(domain.end));
            if hi <= lo {
                return Err(invalid("station: the curve has an empty domain".into()));
            }
            arc_length2(curve, lo, hi).map(Some)
        }
        Curve2::Intrinsic(intrinsic) => {
            if intrinsic.length.is_finite() && intrinsic.length > 0.0 {
                Ok(Some(intrinsic.length))
            } else {
                Err(invalid(
                    "station: the intrinsic curve has no positive finite length".into(),
                ))
            }
        }
        Curve2::Chain(chain) => chain
            .length()
            .map(Some)
            .ok_or_else(|| invalid("station: the chain has no valid length".into())),
        _ => Err(unmeasured("2D")),
    }
}

/// The length a station on a 3D curve is measured against, as
/// [`station_length2`]; an elevated curve's is its plan's, a banked
/// curve's the span of its cant law.
///
/// # Errors
///
/// As [`station_length2`].
pub fn station_length3(curve: &Curve3) -> GeomResult<Option<Scalar>> {
    match curve {
        Curve3::Line(_) => Ok(None),
        Curve3::Circle(_) | Curve3::Ellipse(_) => {
            arc_length3(curve, 0.0, core::f64::consts::TAU).map(Some)
        }
        Curve3::Polyline(_) | Curve3::BSpline(_) => {
            let domain = domain3(curve);
            let (lo, hi) = (domain.start.min(domain.end), domain.start.max(domain.end));
            if hi <= lo {
                return Err(invalid("station: the curve has an empty domain".into()));
            }
            arc_length3(curve, lo, hi).map(Some)
        }
        Curve3::Intrinsic(intrinsic) => {
            if intrinsic.length.is_finite() && intrinsic.length > 0.0 {
                Ok(Some(intrinsic.length))
            } else {
                Err(invalid(
                    "station: the intrinsic curve has no positive finite length".into(),
                ))
            }
        }
        Curve3::Elevated(elevated) => station_length2(&elevated.plan),
        Curve3::Banked(banked) => {
            let domain = domain3(curve);
            if domain.end > domain.start {
                Ok(Some(domain.end))
            } else {
                Err(invalid(format!(
                    "station: the banked curve's cant law has no valid span ({})",
                    banked.span()
                )))
            }
        }
        _ => Err(unmeasured("3D")),
    }
}

/// Whether a station frame on this 2D curve is exact (#264): its point
/// and axes at the station's distance itself, in closed form, rounding
/// aside, with no quadrature or root find behind them.
///
/// Only a line: on every other family the distance is read by the
/// arc-length inverse ([`ARC_LENGTH_TOLERANCE`], an estimate) or the point
/// by quadrature (an intrinsic curve, a chain), so the frame is not exact
/// and a placement in it must not be reported as exact.
#[must_use]
pub fn station_frame_is_exact2(curve: &Curve2) -> bool {
    matches!(curve, Curve2::Line(_))
}

/// [`station_frame_is_exact2`] for a 3D curve: only a line.
#[must_use]
pub fn station_frame_is_exact3(curve: &Curve3) -> bool {
    matches!(curve, Curve3::Line(_))
}

/// The section of a 2D curve at arc length `distance` from its start; see
/// the [module documentation](self) for the conventions.
///
/// # Errors
///
/// A distance that is not finite, negative or beyond the curve's length,
/// an unmeasurable family, and the evaluators' refusals, by name.
pub fn station_section2(curve: &Curve2, distance: Scalar) -> GeomResult<SectionFrame> {
    station_section2_on(curve, distance, SeamSide::Outgoing)
}

/// [`station_section2`] reading `side` of a seam the station lies on (see
/// the [module documentation](self#seams-263)); off a seam the side does
/// not matter.
///
/// # Errors
///
/// As [`station_section2`].
pub fn station_section2_on(
    curve: &Curve2,
    distance: Scalar,
    side: SeamSide,
) -> GeomResult<SectionFrame> {
    let distance = admitted(distance, station_length2(curve)?)?;
    if let Some(section) = seam::seam_section2(curve, distance, side)? {
        return Ok(section);
    }
    let t = parameter_at_arc_length2(curve, start2(curve), distance)?;
    planar_frame(evaluate2(curve, t)?, derivative2(curve, t)?)
}

/// Whether a station at `distance` along a 3D curve lies ON a seam it
/// reads two ways, by the rule [`station_section3_on`] applies (#286):
/// the curve-evaluation provider's sided queries read a station there
/// through it and answer every other measure side-lessly.
///
/// A distance outside `[0, L]` lies on no interior seam: `false`, so the
/// side-less query keeps its own refusal or answer there.
///
/// # Errors
///
/// A curve whose length or seams cannot be read, by name.
pub(crate) fn station_on_seam3(curve: &Curve3, distance: Scalar) -> GeomResult<bool> {
    if !seam::sided3(curve) {
        return Ok(false);
    }
    let Ok(distance) = admitted(distance, station_length3(curve)?) else {
        return Ok(false);
    };
    Ok(seam::seam_at3(curve, distance)?.is_some())
}

/// The section of a 3D curve at `distance` from its start: plan distance
/// on an elevated or banked curve, arc length on any other; see the
/// [module documentation](self).
///
/// # Errors
///
/// As [`station_section2`], and a vertical tangent where the reference-up
/// frame is built.
pub fn station_section3(curve: &Curve3, distance: Scalar) -> GeomResult<SectionFrame> {
    station_section3_on(curve, distance, SeamSide::Outgoing)
}

/// [`station_section3`] reading `side` of a seam the station lies on (see
/// the [module documentation](self#seams-263)); off a seam the side does
/// not matter.
///
/// # Errors
///
/// As [`station_section3`].
pub fn station_section3_on(
    curve: &Curve3,
    distance: Scalar,
    side: SeamSide,
) -> GeomResult<SectionFrame> {
    let distance = admitted(distance, station_length3(curve)?)?;
    if let Some(section) = seam::seam_section3(curve, distance, side)? {
        return Ok(section);
    }
    match curve {
        Curve3::Banked(banked) => banked_frame(banked, distance),
        Curve3::Elevated(elevated) => reference_up_frame(
            crate::arc_length::elevated_point(elevated, distance)?,
            crate::arc_length::elevated_tangent(elevated, distance)?,
        ),
        _ => {
            let t = parameter_at_arc_length3(curve, start3(curve), distance)?;
            reference_up_frame(evaluate3(curve, t)?, derivative3(curve, t)?)
        }
    }
}
