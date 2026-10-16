//! Bounding an unbounded half-space against a finite boundary.
//!
//! A half-space is one side of a plane, infinite in every direction, so it
//! has no mesh of its own. To take part in a finite operation it must first
//! be bounded, and the boundary profile is what bounds it.
//!
//! The region is the kept side of the plane within the prism of the
//! boundary: the lines through the boundary's region parallel to the prism
//! axis, infinite both ways ([`BoundaryPrism`]). The axis is the plane
//! normal ([`bounded_half_space`], [`bounded_half_space_in_frame`]) or any
//! direction that crosses the plane ([`bounded_half_space_prism`], #321),
//! so the prism may stand oblique to the plane. The slab is that region cut
//! off by a plane parallel to the clip plane at a depth on the kept side
//! ([`SlabDepth`]): relative to the boundary's own extent, which keeps the
//! result independent of model units (a boundary in millimetres and the
//! same boundary in metres produce the same shape), or the reach of a
//! subject the slab must cover.
//!
//! Its near cap is the boundary carried along the axis onto the clip plane,
//! its far cap the near one moved along the axis, and its walls the lines
//! between them: the boundary's own prism, exactly, for a polygon.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Plane3, PlaneFrame, Point2, Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_primitive::{ClipMargin, HalfSpace};

use crate::loft::{self, Frame, Station};
use crate::profile::Rings;

/// Where a bounded half-space's boundary stands and which way its prism
/// runs (#321).
///
/// The boundary point `(u, v)` lies at `origin + u x + v y`, and the prism
/// is the set of lines through the boundary's region parallel to `axis`,
/// infinite both ways. `x` and `y` are used as given: their lengths scale
/// the boundary, the angle between them shears it, and their handedness
/// mirrors it. Only the direction of `axis` counts, not its length or sign.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundaryPrism {
    /// Where the boundary's `(0, 0)` lies.
    pub origin: Point3,
    /// The image of the boundary's `+u`.
    pub x: Vec3,
    /// The image of the boundary's `+v`.
    pub y: Vec3,
    /// The direction of the prism's lines.
    pub axis: Vec3,
}

impl BoundaryPrism {
    /// The prism along `plane`'s normal of a boundary authored in `frame`
    /// ([`bounded_half_space_in_frame`]): `frame`'s `x` projected into the
    /// plane fixes the boundary's rotation about the normal, its `y` is the
    /// normal cross that `x`, and its origin is projected onto the plane.
    ///
    /// # Errors
    ///
    /// `InvalidInput`, by name, for a plane without a normal and for a
    /// frame whose `x` is parallel to it.
    pub fn along_normal(frame: PlaneFrame, plane: Plane3) -> GeomResult<Self> {
        let normal = plane_normal(plane)?;
        // Project the authored axes into the clip plane and re-orthonormalise
        // against the sweep normal. The authored frame fixes the rotation
        // about the normal; the sweep direction is still the clip plane's.
        let projected = frame.x_axis() - normal * frame.x_axis().dot(normal);
        let reference = projected.normalize_or_zero();
        if reference == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "authored boundary frame x axis is parallel to the clip plane normal, \
                 so it fixes no in-plane direction"
                    .to_owned(),
            ));
        }
        // The authored origin places the footprint in the plane (#164).
        // Anchoring at `plane.origin` instead silently moves the prism
        // whenever the boundary frame is offset from the clip plane's own
        // origin, which the frame's contract allows. Only the offset ALONG
        // the normal is dropped, exactly as for the axes: the sweep still
        // starts on the clip plane, so polarity and depth are unchanged.
        let offset = frame.origin() - plane.origin;
        let anchor = frame.origin() - normal * offset.dot(normal);
        let base = Frame::from_reference(anchor, normal, reference)?;
        Ok(Self {
            origin: base.origin,
            x: base.x,
            y: base.y,
            axis: normal,
        })
    }

    /// The prism of a boundary authored in `placement`'s `xy` plane and
    /// running along its `z` (#321): the point `(u, v)` at
    /// `placement * (u, v, 0)`, the axis `placement`'s `z` column. This is
    /// how IFC bounds an `IfcPolygonalBoundedHalfSpace`, by its
    /// `Position`.
    ///
    /// # Errors
    ///
    /// `InvalidInput`, by name, for a non-finite placement and for a
    /// degenerate one, whose `x`, `y` and `z` columns do not span space
    /// within the angular tolerance (parallel or zero axes): its boundary
    /// has no area, or its prism no direction off the boundary's plane.
    pub fn along_placement_z(placement: Transform3, tolerance: Tolerance) -> GeomResult<Self> {
        let [x, y, z] = [
            placement.matrix3.x_axis,
            placement.matrix3.y_axis,
            placement.matrix3.z_axis,
        ];
        if !(placement.translation.is_finite() && x.is_finite() && y.is_finite() && z.is_finite()) {
            return Err(GeomError::InvalidInput(
                "bounded half-space placement is not finite".to_owned(),
            ));
        }
        let volume = x.cross(y).dot(z).abs();
        if volume <= unit_floor(tolerance) * x.length() * y.length() * z.length() {
            return Err(GeomError::InvalidInput(
                "bounded half-space placement is degenerate: its axes are parallel or zero"
                    .to_owned(),
            ));
        }
        Ok(Self {
            origin: placement.translation,
            x,
            y,
            axis: z,
        })
    }

    /// The cosine of the prism's axis against `plane`'s unit normal, signed:
    /// how steeply, and from which side, the prism crosses the plane.
    ///
    /// # Errors
    ///
    /// `InvalidInput`, by name, for a plane without a normal, an axis that
    /// is zero or not finite, and an axis within the angular tolerance of
    /// the plane: its prism never crosses it.
    pub fn crossing(&self, plane: Plane3, tolerance: Tolerance) -> GeomResult<Scalar> {
        let normal = plane_normal(plane)?;
        let axis = self.axis.normalize_or_zero();
        if !axis.is_finite() || axis == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "bounded half-space prism axis must be a non-zero direction".to_owned(),
            ));
        }
        let crossing = axis.dot(normal);
        if crossing.abs() <= unit_floor(tolerance) {
            return Err(GeomError::InvalidInput(
                "bounded half-space prism axis is parallel to the clip plane, so its prism \
                 never crosses it"
                    .to_owned(),
            ));
        }
        Ok(crossing)
    }
}

/// How far a bounded half-space's slab reaches past its plane on the kept
/// side, measured perpendicular to the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SlabDepth {
    /// The boundary's own extent on the plane, its largest coordinate
    /// scaled by its longest axis there, times the margin: independent of
    /// the model's units, and of everything the slab meets.
    Relative(ClipMargin),
    /// This distance: the reach of the subject the slab must cover, so that
    /// a boolean with the slab is the boolean with the half-space.
    Reach(Scalar),
}

/// Build the finite solid that represents a half-space bounded by a profile.
///
/// `agreement` selects which side of the plane survives: `true` keeps the
/// normal side, `false` the opposite one. The boundary profile is assumed
/// to lie in the plane's own frame, which is how the graph stores it.
///
/// The margin scales the extrusion depth relative to the boundary's own
/// size. It is relative rather than absolute so the result does not depend
/// on the model's units.
pub fn bounded_half_space(
    boundary: &Rings,
    plane: Plane3,
    agreement: bool,
    margin: ClipMargin,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    bounded_half_space_framed(boundary, plane, None, agreement, margin, tolerance)
}

/// Build a bounded half-space with an explicitly authored in-plane frame.
///
/// [`bounded_half_space`] derives the boundary's in-plane axes from the clip
/// plane normal alone, which fixes only 2 of the 3 orientation degrees of
/// freedom: the rotation ABOUT the normal is picked by an internal heuristic.
/// That is fine when the profile carries no authored orientation of its own.
///
/// It is not fine when the source format authors one. IFC's
/// `IfcPolygonalBoundedHalfSpace.Position` is an placement explicitly
/// independent of `BaseSurface`, and its `PolygonalBoundary` is expressed in
/// THAT frame. Deriving the axes instead of honouring the authored ones
/// silently rotates the clipping profile about the plane normal, so the wrong
/// region of the subject is cut.
///
/// The frame's own axes are used as written. Its normal is not required to
/// match `plane.normal`: only the in-plane axes are taken from it, and the
/// sweep still runs along the clip plane's normal, so an authored frame that
/// is tilted relative to the clip plane still contributes only its rotation
/// about that normal.
///
/// The frame's origin anchors the boundary, projected onto the clip plane: an
/// in-plane offset from `plane.origin` moves the footprint by exactly that
/// offset, and the component along the normal is dropped so the sweep still
/// starts on the clip plane.
pub fn bounded_half_space_in_frame(
    boundary: &Rings,
    plane: Plane3,
    frame: PlaneFrame,
    agreement: bool,
    margin: ClipMargin,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    bounded_half_space_framed(boundary, plane, Some(frame), agreement, margin, tolerance)
}

fn bounded_half_space_framed(
    boundary: &Rings,
    plane: Plane3,
    authored: Option<PlaneFrame>,
    agreement: bool,
    margin: ClipMargin,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    // Frame the profile in the plane. The plane normal is the sweep
    // direction, so the profile's own x and y span the plane.
    //
    // An authored frame is used as written: it is the caller's statement of
    // where the profile's x and y point, and re-deriving them would discard
    // the one piece of information the caller supplied. Without one, the
    // rotation about the normal is unconstrained, so a reference axis is
    // picked -- the least parallel of X and Y, so the cross product stays
    // well conditioned.
    let prism = match authored {
        Some(frame) => BoundaryPrism::along_normal(frame, plane)?,
        None => {
            let normal = plane_normal(plane)?;
            let reference = if normal.x.abs() < 0.9 {
                Vec3::X
            } else {
                Vec3::Y
            };
            let base = Frame::from_reference(plane.origin, normal, reference)?;
            BoundaryPrism {
                origin: base.origin,
                x: base.x,
                y: base.y,
                axis: normal,
            }
        }
    };
    bounded_half_space_prism(
        boundary,
        plane,
        prism,
        agreement,
        SlabDepth::Relative(margin),
        tolerance,
    )
}

/// Build the slab of a half-space bounded by the prism of `boundary` along
/// any axis that crosses the plane (#321).
///
/// The region is the side of `plane` that `agreement` selects (`true` the
/// normal side) within `prism`, cut off at `depth` from the plane. The
/// prism may stand at any angle to the plane: the near cap is the boundary
/// carried along the axis onto the plane, the far cap the same moved along
/// the axis to the depth. The boundary's winding is its own affair: a ring
/// wound either way, under a frame of either handedness, gives the same
/// outward-facing solid, the loft's stations ordered to suit.
///
/// # Errors
///
/// `InvalidInput`, by name, for a plane without a normal, a boundary of
/// fewer than 3 vertices, a prism axis that is zero, not finite or within
/// the angular tolerance of the plane (its prism never crosses it), a frame
/// whose `x` and `y` carried onto the plane are parallel or zero, and a
/// depth that is not positive and finite. A boundary that bounds no region
/// is refused by the cap triangulation, naming the ring.
pub fn bounded_half_space_prism(
    boundary: &Rings,
    plane: Plane3,
    prism: BoundaryPrism,
    agreement: bool,
    depth: SlabDepth,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    let normal = plane_normal(plane)?;
    if boundary.outer.len() < 3 {
        return Err(GeomError::InvalidInput(format!(
            "half-space boundary needs at least 3 vertices, got {}",
            boundary.outer.len()
        )));
    }
    // The cosine of the axis against the kept side: how steeply the prism
    // crosses the plane. Its sign turns the axis to the kept side, which
    // is all a prism infinite both ways asks of it.
    let crossing = prism.crossing(plane, tolerance)?;
    let lean = if agreement { crossing } else { -crossing };
    let axis = prism.axis.normalize();
    let toward = axis * lean.signum();
    // Carry the frame along the axis onto the plane: a point `p` lands at
    // `p - axis ((p - origin) . n) / (axis . n)`, an affine map, so the
    // near cap is the boundary in the carried frame.
    let carry = |v: Vec3| v - axis * (v.dot(normal) / crossing);
    let origin = prism.origin - axis * ((prism.origin - plane.origin).dot(normal) / crossing);
    let (x, y) = (carry(prism.x), carry(prism.y));
    let face = x.cross(y);
    if !(origin.is_finite() && face.is_finite())
        || face.length() <= unit_floor(tolerance) * x.length() * y.length()
    {
        return Err(GeomError::InvalidInput(
            "bounded half-space boundary frame is degenerate on the clip plane: its axes \
             are parallel or zero"
                .to_owned(),
        ));
    }
    // Size the slab. Relative to the boundary, a profile spanning 10 units
    // gets a slab proportional to 10, so the construction carries no hidden
    // absolute length.
    let reach = match depth {
        SlabDepth::Relative(margin) => {
            let mut extent: Scalar = 0.0;
            for p in boundary.outer.iter().chain(boundary.holes.iter().flatten()) {
                extent = extent.max(p.x.abs()).max(p.y.abs());
            }
            // No zero-extent guard here: a boundary with no area cannot be
            // triangulated, so `loft` rejects it downstream with a message
            // naming the real problem.
            extent * x.length().max(y.length()) * margin.factor()
        }
        SlabDepth::Reach(reach) => reach,
    };
    if !(reach.is_finite() && reach > 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "bounded half-space slab depth must be positive and finite, got {reach}"
        )));
    }
    // Along the axis, the far cap stands `reach` from the plane.
    let step = toward * (reach / lean.abs());

    // `loft` winds its walls and caps outward when the stations run along
    // the face normal of a counter-clockwise ring, so a clockwise ring is
    // read reversed (the same region), and the stations run near to far
    // when that face normal points to the kept side, far to near when the
    // frame (or the side) turns it over. Swapping stations rather than
    // mirroring an axis keeps the boundary's footprint untouched.
    let rings = counter_clockwise(boundary);
    let near = Frame { origin, x, y };
    let far = Frame {
        origin: origin + step,
        x,
        y,
    };
    let ordered = if face.dot(toward) > 0.0 {
        [near, far]
    } else {
        [far, near]
    };
    let stations: Vec<Station> = ordered
        .iter()
        .map(|f| loft::place(&rings, |p| loft::at(f, p)))
        .collect();
    loft::loft(&rings, &stations, false)
}

/// The plane's unit normal, refused when it has none.
fn plane_normal(plane: Plane3) -> GeomResult<Vec3> {
    let normal = plane.normal.normalize_or_zero();
    if normal == Vec3::ZERO || !normal.is_finite() {
        return Err(GeomError::InvalidInput(
            "half-space boundary plane needs a non-zero normal".to_owned(),
        ));
    }
    Ok(normal)
}

/// The smallest sine (or cosine) a direction test accepts: the angular
/// tolerance, floored at a few ulps so `Tolerance::ZERO` still refuses an
/// exactly degenerate frame.
fn unit_floor(tolerance: Tolerance) -> Scalar {
    tolerance.angular().max(Scalar::EPSILON * 8.0)
}

/// `rings` with its outer ring counter-clockwise, every ring reversed when
/// it was not: the same region, wound as `loft` expects.
fn counter_clockwise(rings: &Rings) -> std::borrow::Cow<'_, Rings> {
    let twice_area: Scalar = (0..rings.outer.len())
        .map(|i| {
            let (a, b) = (rings.outer[i], rings.outer[(i + 1) % rings.outer.len()]);
            a.perp_dot(b)
        })
        .sum();
    if twice_area >= 0.0 {
        return std::borrow::Cow::Borrowed(rings);
    }
    let reversed = |ring: &Vec<Point2>| ring.iter().rev().copied().collect();
    std::borrow::Cow::Owned(Rings {
        outer: reversed(&rings.outer),
        holes: rings.holes.iter().map(reversed).collect(),
    })
}

pub fn for_subject(
    subject: &TriMesh,
    half_space: HalfSpace,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    let normal = half_space.boundary.normal.normalize_or_zero();
    if normal == Vec3::ZERO || subject.positions.is_empty() {
        return Err(GeomError::InvalidInput(
            "half-space boolean needs finite subject bounds".into(),
        ));
    }
    let reference = if normal.x.abs() < 0.9 {
        Vec3::X
    } else {
        Vec3::Y
    };
    let frame = Frame::from_reference(half_space.boundary.origin, normal, reference)?;
    let mut min = Point2::splat(Scalar::INFINITY);
    let mut max = Point2::splat(Scalar::NEG_INFINITY);
    let side = if half_space.agreement { 1.0 } else { -1.0 };
    let mut selected_depth: Scalar = 0.0;
    for &point in &subject.positions {
        if !point.is_finite() {
            return Err(GeomError::InvalidInput(
                "half-space subject contains non-finite points".into(),
            ));
        }
        let delta = point - half_space.boundary.origin;
        let uv = Point2::new(delta.dot(frame.x), delta.dot(frame.y));
        min = min.min(uv);
        max = max.max(uv);
        selected_depth = selected_depth.max(delta.dot(normal) * side);
    }
    let span = (max - min).max_element().max(tolerance.linear());
    let pad = span * 0.1 + tolerance.linear();
    min -= Point2::splat(pad);
    max += Point2::splat(pad);
    let boundary = Rings {
        outer: vec![
            min,
            Point2::new(max.x, min.y),
            max,
            Point2::new(min.x, max.y),
        ],
        holes: Vec::new(),
    };
    let extent = min.abs().max(max.abs()).max_element();
    let depth = selected_depth + pad;
    let factor = depth / extent;
    let margin = ClipMargin::new(factor).ok_or_else(|| {
        GeomError::InvalidInput("half-space subject has no finite clipping extent".into())
    })?;
    bounded_half_space(
        &boundary,
        half_space.boundary,
        half_space.agreement,
        margin,
        tolerance,
    )
}
