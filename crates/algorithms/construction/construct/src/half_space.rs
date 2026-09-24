//! Bounding an unbounded half-space against a finite boundary.
//!
//! A half-space is one side of a plane, infinite in every direction, so it
//! has no mesh of its own. To take part in a finite operation it must first
//! be bounded, and the boundary profile is what bounds it.
//!
//! The slab is built by extruding the boundary profile along the plane
//! normal, far enough in each direction to cover the boundary's own extent,
//! then keeping the selected side. Sizing from the boundary rather than
//! from a fixed constant is what keeps the result independent of model
//! units: a boundary in millimetres and the same boundary in metres
//! produce the same shape.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Plane3, PlaneFrame, Point2, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_primitive::{ClipMargin, HalfSpace};

use crate::loft::{self, Frame, Station};
use crate::profile::Rings;

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
    let normal = plane.normal.normalize_or_zero();
    if normal == Vec3::ZERO {
        return Err(GeomError::InvalidInput(
            "half-space boundary plane needs a non-zero normal".to_owned(),
        ));
    }
    if boundary.outer.len() < 3 {
        return Err(GeomError::InvalidInput(format!(
            "half-space boundary needs at least 3 vertices, got {}",
            boundary.outer.len()
        )));
    }

    // Size the slab from the boundary's own extent. A profile spanning 10
    // units gets a slab proportional to 10, so the construction carries no
    // hidden absolute length.
    let mut extent: Scalar = 0.0;
    for p in boundary.outer.iter().chain(boundary.holes.iter().flatten()) {
        extent = extent.max(p.x.abs()).max(p.y.abs());
    }
    // No zero-extent guard here: a boundary with no area cannot be
    // triangulated, so `loft` rejects it downstream with a message naming
    // the real problem. A guard here was unreachable, and unreachable
    // validation reads as a promise the code cannot keep.
    let depth = extent * margin.factor();

    // Frame the profile in the plane. The plane normal is the sweep
    // direction, so the profile's own x and y span the plane.
    //
    // An authored frame is used as written: it is the caller's statement of
    // where the profile's x and y point, and re-deriving them would discard
    // the one piece of information the caller supplied. Without one, the
    // rotation about the normal is unconstrained, so a reference axis is
    // picked -- the least parallel of X and Y, so the cross product stays
    // well conditioned.
    let base = match authored {
        Some(frame) => {
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
            Frame::from_reference(anchor, normal, reference)?
        }
        None => {
            let reference = if normal.x.abs() < 0.9 {
                Vec3::X
            } else {
                Vec3::Y
            };
            Frame::from_reference(plane.origin, normal, reference)?
        }
    };

    // Keeping the normal side sweeps from the plane along +normal; the
    // opposite side sweeps the other way.
    let step = if agreement { normal } else { -normal };
    let near = Frame {
        origin: base.origin,
        x: base.x,
        y: base.y,
    };
    let far = Frame {
        origin: base.origin + step * depth,
        x: base.x,
        y: base.y,
    };
    // Order the stations so traversal always runs along +normal. `loft`
    // derives its winding from station order, and the frame's own handedness
    // is fixed (x cross y == +normal), so the two must agree or the solid
    // comes out inside-out. Swapping the order rather than mirroring an
    // in-plane axis keeps the boundary's footprint untouched: `agreement`
    // selects a side, it never reshapes the profile.
    let ordered = if agreement { [near, far] } else { [far, near] };
    let stations: Vec<Station> = ordered
        .iter()
        .map(|f| loft::place(boundary, |p| loft::at(f, p)))
        .collect();
    let _ = tolerance;
    loft::loft(boundary, &stations, false)
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
