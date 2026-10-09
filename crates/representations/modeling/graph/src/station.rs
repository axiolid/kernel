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
//! # Explicit orientation (#246)
//!
//! A station may carry a [`StationOrientation`]: an optional `axis` and an
//! optional `ref_direction`, both given as COMPONENTS IN THE STATION'S BASE
//! FRAME at its distance (the [`StationFrame`] it names), never in world
//! coordinates: a vector `(a, b, c)` means `a * tangent + b * lateral + c *
//! up`. This is how a linear placement reads its own axes: they "are
//! relative to the curve used for linear referencing ..., maintaining the
//! relationship to the tangent of the curve" (buildingSMART IFC 4.3,
//! `IfcAxis2PlacementLinear`), whose local `X` is the tangent, `Y` the left
//! lateral and `Z` up, this module's `(tangent, lateral, up)`.
//!
//! - `axis` is the exact direction of the oriented frame's up (local `Z`);
//!   absent, it is `(0, 0, 1)`, the base frame's up.
//! - `ref_direction` fixes the oriented frame's tangent role (local `X`);
//!   absent, it is `(1, 0, 0)`, the curve tangent.
//! - The frame is orthonormalised by Gram-Schmidt with the axis primary:
//!   `up' = axis / |axis|`, `tangent' = normalise(r - (r . up') up')` for
//!   the unit `r = ref_direction / |ref_direction|`, and `lateral' = up' x
//!   tangent'`, so `(tangent', lateral', up')` stays right-handed.
//! - A zero or non-finite vector, and an axis and reference direction
//!   (given or defaulted) parallel or anti-parallel within
//!   [`ORIENTATION_TOLERANCE`] (the sine of the angle between them), are
//!   refused by name when the node is pushed.
//!
//! The oriented triad replaces the base one wherever a frame ORIENTS
//! something: the frame a resolved station presents, and the plane a
//! section's profile is placed in (profile `x` along `lateral'`, `y`
//! along `up'`). The OFFSETS stay in the base frame: they locate the
//! station's origin, which the orientation does not move.
//!
//! # Seams (#263)
//!
//! A basis curve made of pieces -- a polyline, a B-spline with a corner
//! knot, an elevated curve whose profile's grade breaks, a banked curve
//! whose cant or pivot law jumps -- has two frames where two pieces meet.
//! A station within the arc-length tolerance of such a seam is ON it and
//! reads the frame of one piece, chosen by [`SeamSide`]: the piece that
//! starts there ([`SeamSide::Outgoing`], the default and what every curve
//! evaluator reads) or the one that ends there ([`SeamSide::Incoming`]).
//! An [`OrientedCurveStation`] carries the side ([`OrientedCurveStation::seam`]);
//! a plain [`CurveStation`] reads the outgoing piece, and
//! [`CurveStation::with_seam_side`] turns it into an oriented station in
//! its base frame that reads the other one.
//!
//! A run of sections or offsets ([`Station`]s along one directrix) needs
//! no side: a section standing exactly on an interior seam where the
//! tangent turns is cut in the mitre plane, the bisector of the incoming
//! and outgoing tangents, and so is the run where it crosses such a seam
//! between two of its stations. The run's first station reads the
//! outgoing piece and its last the incoming one, the pieces the run lies
//! on. A seam where the tangent does not turn keeps the outgoing frame.
//!
//! # Placing a node at a station (#264)
//!
//! [`InstanceAtStation`] reuses a node -- a curve, a solid, a surface --
//! in the frame of an [`OrientedCurveStation`], the way [`Instance`]
//! reuses one under a fixed transform. The station stays symbolic: a
//! kernel resolves its frame when it evaluates the node, with the
//! station's [`StationFrame`], orientation and [`SeamSide`]. The source's
//! local axes map onto the station's oriented frame as a linear placement
//! reads them (buildingSMART IFC 4.3, `IfcAxis2PlacementLinear`):
//!
//! - local `x` onto the oriented tangent (`tangent'`),
//! - local `y` onto the oriented lateral, to the LEFT (`lateral'`),
//! - local `z` onto the oriented up (`up' = tangent' x lateral'`),
//! - the local origin onto the station's point, its offsets read in the
//!   base frame.
//!
//! This is a rigid motion (right-handed, no scale). It is not the profile
//! mapping above (a section's `x` along lateral and `y` along up), and not
//! the provider layout a resolved station presents (`x` tangent, `y` up,
//! `z` right): a placed curve is laid ALONG the curve, a section ACROSS
//! it. A 2D source curve lies in its local `z = 0`, so it is placed in the
//! plane of the oriented tangent and lateral, and the placed curve is a 3D
//! curve whatever the source's dimension.
//!
//! [`Instance`]: crate::Instance
//!
//! # Sections between stations
//!
//! A profile placed at a station maps its `x` onto the lateral axis and its
//! `y` onto up, so its normal `x x y` is the tangent. Between two stations
//! everything is interpolated LINEARLY in distance: the offsets, each
//! section point with its counterpart in the next section, and an explicit
//! orientation (the unit axis and the unit reference direction,
//! component-wise in the base frame, an absent one counting as `(0, 0, 1)`
//! and `(1, 0, 0)`, then orthonormalised as above; a resolver refuses by
//! name an interpolated pair that degenerates).
//!
//! # Matching points between sections
//!
//! Untagged sections are matched by position: a closed profile's points by
//! ring and vertex index, an open section's by vertex index. Tagged
//! sections ([`SectionAtStation::tags`], [`StationedOpenSection::tags`])
//! are matched by tag. One run of sections is tagged throughout or not at
//! all, and every section carries the same SET of tags, none repeated, so
//! the tags are a bijection between any two sections:
//!
//! - an open section's tags run along its polyline in the first section's
//!   order or in reverse (a reversed section is joined reversed); any other
//!   order would cross the sheet and is refused;
//! - a closed section's tags name its contour vertices, the outer ring's
//!   then each hole's, each ring from its first segment's start in its
//!   authored sense. Its profile must be polygonal (straight segments only,
//!   optionally under a `Profile::Derived` transform). The tags must map
//!   each ring of the first section onto one ring of every other, the
//!   outer onto the outer, in the same cyclic order once both are wound
//!   alike (outer counter-clockwise, holes clockwise): a section may start
//!   a ring at another vertex, list its holes in another order or wind a
//!   ring the other way, nothing else. A resolver checks this, since it
//!   needs the profile's vertices.

use axiolid_core::{Scalar, Vec3};
pub use axiolid_curve::SeamSide;

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

    /// This station, unturned, reading `side` of a seam it lies on (#263).
    ///
    /// A plain curve station always reads the outgoing piece; the side is
    /// carried by an [`OrientedCurveStation`] whose orientation is the base
    /// frame, so the node to push is
    /// [`GeometryNode::OrientedCurveStation`](crate::GeometryNode::OrientedCurveStation).
    #[must_use]
    pub const fn with_seam_side(self, side: SeamSide) -> OrientedCurveStation {
        OrientedCurveStation::new(self, StationOrientation::new(None, None)).with_seam_side(side)
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
/// section of one surface carries the same tags, in the first section's
/// order or in reverse, or none (then the vertices are joined by index);
/// see the [module documentation](self).
#[derive(Debug, Clone, PartialEq)]
pub struct StationedOpenSection {
    /// Open profile node.
    pub profile: NodeId,
    /// One tag per vertex, or empty.
    pub tags: Vec<String>,
    /// Distance along the directrix and offsets of the profile's origin.
    pub station: Station,
}

/// Largest sine of the angle between a [`StationOrientation`]'s axis and
/// reference direction that still counts as parallel.
pub const ORIENTATION_TOLERANCE: Scalar = 1e-9;

/// An explicit in-section orientation at a station (#246).
///
/// Both vectors are components in the station's base frame, `(tangent,
/// lateral, up)`; see the [module documentation](self) for their meaning,
/// defaults and orthonormalisation. Both absent is the base frame itself,
/// the [`Default`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct StationOrientation {
    /// The oriented frame's up (local `Z`), exact; `None` for `(0, 0, 1)`.
    pub axis: Option<Vec3>,
    /// The direction the oriented frame's tangent role (local `X`) is
    /// taken from, made perpendicular to the axis; `None` for `(1, 0, 0)`.
    pub ref_direction: Option<Vec3>,
}

impl StationOrientation {
    /// An orientation from an optional axis and reference direction, both
    /// in the base frame.
    #[must_use]
    pub const fn new(axis: Option<Vec3>, ref_direction: Option<Vec3>) -> Self {
        Self {
            axis,
            ref_direction,
        }
    }

    /// Whether neither vector is given, so the base frame is used as is.
    #[must_use]
    pub const fn is_base(&self) -> bool {
        self.axis.is_none() && self.ref_direction.is_none()
    }

    /// The unit axis and unit reference direction, defaults filled in.
    ///
    /// # Errors
    ///
    /// Why there are none, worded for
    /// [`GraphError::InvalidStation`](crate::GraphError::InvalidStation):
    /// a zero or non-finite vector, or the two parallel or anti-parallel
    /// within [`ORIENTATION_TOLERANCE`].
    pub fn unit_axes(&self) -> Result<(Vec3, Vec3), &'static str> {
        let unit = |vector: Vec3, what: &'static str| {
            let length = vector.length();
            if vector.is_finite() && length.is_finite() && length > 1e-12 {
                Ok(vector / length)
            } else {
                Err(what)
            }
        };
        let axis = unit(
            self.axis.unwrap_or(Vec3::Z),
            "the orientation's axis is zero or not finite",
        )?;
        let reference = unit(
            self.ref_direction.unwrap_or(Vec3::X),
            "the orientation's reference direction is zero or not finite",
        )?;
        if axis.cross(reference).length() <= ORIENTATION_TOLERANCE {
            return Err("the orientation's axis and reference direction are parallel");
        }
        Ok((axis, reference))
    }
}

/// A [`CurveStation`] with an explicit orientation (#246) and the side of
/// a seam it reads (#263).
///
/// Its point is the station's, offsets read in the base frame; its frame
/// is the base frame turned by `orientation` (see the
/// [module documentation](self)). On a seam of the basis curve the base
/// frame is the one of the piece `seam` names.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientedCurveStation {
    /// The basis curve, distance, offsets and base frame.
    pub station: CurveStation,
    /// The axis and reference direction, in the base frame.
    pub orientation: StationOrientation,
    /// The piece read when the station lies on a seam of its basis curve;
    /// [`SeamSide::Outgoing`] unless set by [`Self::with_seam_side`].
    pub seam: SeamSide,
}

impl OrientedCurveStation {
    /// `station` turned by `orientation`, reading the outgoing piece on a
    /// seam.
    #[must_use]
    pub const fn new(station: CurveStation, orientation: StationOrientation) -> Self {
        Self {
            station,
            orientation,
            seam: SeamSide::Outgoing,
        }
    }

    /// The same station reading `side` of a seam it lies on.
    #[must_use]
    pub const fn with_seam_side(mut self, side: SeamSide) -> Self {
        self.seam = side;
        self
    }
}

impl From<CurveStation> for OrientedCurveStation {
    /// The station in its base frame, reading the outgoing piece on a seam:
    /// what the plain station means.
    fn from(station: CurveStation) -> Self {
        Self::new(station, StationOrientation::default())
    }
}

/// A section standing at a station, with optional tags and an optional
/// explicit orientation (#246): the general form of [`StationedSection`]
/// (closed) and [`StationedOpenSection`] (open).
///
/// Built with [`SectionAtStation::new`] and the `with_` methods, so that
/// fields can be added without breaking callers.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct SectionAtStation {
    /// The profile node: an area profile in a sectioned spine, an open
    /// profile in a sectioned surface.
    pub profile: NodeId,
    /// Distance along the directrix and offsets of the profile's origin,
    /// in the base frame.
    pub station: Station,
    /// One tag per vertex, or empty to match by index; see the
    /// [module documentation](self).
    pub tags: Vec<String>,
    /// The section's orientation in the base frame; the default is the
    /// base frame.
    pub orientation: StationOrientation,
}

impl SectionAtStation {
    /// `profile` at `station`, untagged, in the base frame.
    #[must_use]
    pub const fn new(profile: NodeId, station: Station) -> Self {
        Self {
            profile,
            station,
            tags: Vec::new(),
            orientation: StationOrientation::new(None, None),
        }
    }

    /// The same section with `tags`, one per vertex.
    #[must_use]
    pub fn with_tags<I, T>(mut self, tags: I) -> Self
    where
        I: IntoIterator<Item = T>,
        T: Into<String>,
    {
        self.tags = tags.into_iter().map(Into::into).collect();
        self
    }

    /// The same section turned by `orientation`.
    #[must_use]
    pub const fn with_orientation(mut self, orientation: StationOrientation) -> Self {
        self.orientation = orientation;
        self
    }
}

impl From<StationedSection> for SectionAtStation {
    fn from(section: StationedSection) -> Self {
        Self::new(section.profile, section.station)
    }
}

impl From<StationedOpenSection> for SectionAtStation {
    fn from(section: StationedOpenSection) -> Self {
        Self::new(section.profile, section.station).with_tags(section.tags)
    }
}

/// A node reused in the frame of an [`OrientedCurveStation`] (#264): the
/// station-framed form of [`Instance`](crate::Instance).
///
/// The source's local `x`, `y`, `z` and origin are placed on the
/// station's oriented tangent, left lateral, up and point (see the
/// [module documentation](self#placing-a-node-at-a-station-264)); the
/// frame is resolved at evaluation, reading the side of a seam the station
/// names. A placed curve is a 3D curve, a placed solid or surface keeps
/// its family. Built with [`InstanceAtStation::new`], so fields can be
/// added without breaking callers.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InstanceAtStation {
    /// The reused node, in its own local coordinates.
    pub source: NodeId,
    /// The station whose oriented frame places it.
    pub station: OrientedCurveStation,
}

impl InstanceAtStation {
    /// `source` placed in the frame of `station`.
    #[must_use]
    pub const fn new(source: NodeId, station: OrientedCurveStation) -> Self {
        Self { source, station }
    }
}
