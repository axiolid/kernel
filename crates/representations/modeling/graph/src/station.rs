//! Stations: positions named by a distance along a basis curve (#241,
//! ADR 0082).
//!
//! A station is exact stored data -- a distance and three offsets -- that a
//! kernel resolves by evaluating the basis curve; the graph never does.
//!
//! # Distance
//!
//! [`Station::distance`] is measured from the basis curve's start in that
//! curve's own convention: PLAN distance on an elevated or a banked curve
//! (their parameter), arc length on every other curve, 2D or 3D. The start
//! is parameter `0`, or the first parameter of the domain of a polyline or
//! a B-spline. A resolver refuses by name a distance before the start or
//! past the curve's length; the graph refuses a non-finite or negative one
//! when the node is pushed.
//!
//! # Frame and offsets
//!
//! The basis curve's section frame at the distance has three unit axes:
//! the tangent, the lateral axis to the LEFT of it, and `up = tangent x
//! lateral`. For a 2D curve (lying in `z = 0`) they are the tangent, the
//! left normal and `+Z`; for a banked curve, its rolled section (lateral
//! towards the left rail head); for every other 3D curve, an elevated one
//! included, the reference-up frame against `+Z` (lateral horizontal, up
//! leaning back with the grade). [`StationFrame::Plan`] instead takes the
//! tangent's horizontal projection, the horizontal left normal and `+Z`.
//!
//! The station's point is `point + lateral * l + up * v + tangent * g` for
//! the offsets `(l, v, g)` of [`StationOffsets`]. A resolver presents the
//! frame as `x` the tangent, `y` up, `z` to the right (`-lateral`), the
//! curve-evaluation provider's layout.
//!
//! # Sections between stations
//!
//! A profile placed at a station maps its `x` onto the lateral axis and its
//! `y` onto up, so its normal `x x y` is the tangent. Between two stations
//! everything is interpolated LINEARLY in distance: the offsets, and each
//! section point with its counterpart in the next section (a closed
//! profile's by ring and vertex index, an open section's by tag).

use axiolid_core::Scalar;

use crate::NodeId;

/// Offsets of a station from its basis curve, in the curve's section frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StationOffsets {
    /// Along the lateral axis, positive to the left of the tangent.
    pub lateral: Scalar,
    /// Along the section's up axis.
    pub vertical: Scalar,
    /// Along the tangent, positive in the direction of increasing distance.
    pub longitudinal: Scalar,
}

impl StationOffsets {
    /// No offset: the station is on the curve.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    /// Offsets along the lateral, up and tangent axes.
    #[must_use]
    pub const fn new(lateral: Scalar, vertical: Scalar, longitudinal: Scalar) -> Self {
        Self {
            lateral,
            vertical,
            longitudinal,
        }
    }

    /// Whether every offset is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.lateral.is_finite() && self.vertical.is_finite() && self.longitudinal.is_finite()
    }
}

/// A distance along a basis curve with offsets in its section frame.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Station {
    /// Distance from the basis curve's start, in its convention (see the
    /// [module documentation](self)).
    pub distance: Scalar,
    /// Offsets in the section frame at that distance.
    pub offsets: StationOffsets,
}

impl Station {
    /// A station at `distance` with `offsets`.
    #[must_use]
    pub const fn new(distance: Scalar, offsets: StationOffsets) -> Self {
        Self { distance, offsets }
    }

    /// A station on the curve at `distance`.
    #[must_use]
    pub const fn at(distance: Scalar) -> Self {
        Self::new(distance, StationOffsets::ZERO)
    }

    /// Whether the distance is finite and not negative and every offset is
    /// finite.
    #[must_use]
    pub fn is_well_formed(&self) -> bool {
        self.distance.is_finite() && self.distance >= 0.0 && self.offsets.is_finite()
    }
}

/// Which frame a station's offsets and orientation are read in.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum StationFrame {
    /// The basis curve's section frame: rolled on a banked curve, leaning
    /// with the grade on any 3D curve.
    #[default]
    Section,
    /// The vertical frame: the tangent's horizontal projection, the
    /// horizontal left normal and `+Z`. Grade and bank are dropped, so a
    /// section stands upright.
    Plan,
}

/// A point and frame at a station along a basis curve.
///
/// Resolving it is curve evaluation, which a kernel does; see the
/// [module documentation](self) for the conventions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveStation {
    /// Basis curve, 2D or 3D.
    pub basis: NodeId,
    /// Distance and offsets.
    pub station: Station,
    /// Frame the offsets are read in, and that the station presents.
    pub frame: StationFrame,
}

impl CurveStation {
    /// A station along `basis` in the basis curve's section frame.
    #[must_use]
    pub const fn new(basis: NodeId, station: Station) -> Self {
        Self {
            basis,
            station,
            frame: StationFrame::Section,
        }
    }
}

/// A closed profile standing at a station along a sectioned spine's
/// directrix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StationedSection {
    /// Area profile node.
    pub profile: NodeId,
    /// Distance along the directrix and offsets of the profile's origin.
    pub station: Station,
}

/// An open section standing at a station along a sectioned surface's
/// directrix, its points named by tags.
///
/// `tags` names the vertices of the section's polyline in order; two
/// consecutive sections are joined point to point by equal tags. Every
/// section of one surface carries the same tag sequence, or none (then the
/// vertices are joined by index).
#[derive(Debug, Clone, PartialEq)]
pub struct StationedOpenSection {
    /// Open profile node.
    pub profile: NodeId,
    /// One tag per vertex, or empty.
    pub tags: Vec<String>,
    /// Distance along the directrix and offsets of the profile's origin.
    pub station: Station,
}
