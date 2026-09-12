//! Chamfer and fillet on a straight vertical edge of an exact prism (#68).
//!
//! # The contract
//!
//! An edge is selected by the CORNER it belongs to, not by an opaque index:
//! a caller naming "the corner nearest this point" survives a topology change
//! that renumbers edges, while an index does not.
//!
//! # What is exact here
//!
//! Chamfering a vertical edge of a prism replaces one corner of its profile
//! with a straight cut, so the result is a prism over a polygon with one more
//! vertex. That is exactly the shape `extrude_polygon_rings` already builds --
//! no new surface families, no approximation.
//!
//! A constant-radius FILLET replaces the corner with a circular arc, so the
//! wall becomes a cylindrical face. That is representable, but stitching a
//! cylindrical wall into the prism assembly is not the same construction, so
//! it is refused here rather than approximated by a polyline. Refusing is the
//! honest answer: a caller asking for a fillet and receiving a many-segment
//! chamfer would have no way to tell.
//!
//! # Refusals
//!
//! Variable radius, edge loops, curved edges, and fillets all return typed
//! refusals naming the missing capability, never a silent no-op. A no-op is
//! the worst outcome: the caller believes the feature was applied.

use axiolid_brep::{EdgeName, ExactBRep, FaceName, SweptFace};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Point2, Scalar, Tolerance, Vec2, Vec3};
use axiolid_profile::{Profile, RectangleProfile};

use crate::extrude_exact::{extrude_polygon_rings, extrude_with_cylindrical_blend};
use crate::BACKEND_ID;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// Which edge a feature applies to.
///
/// Selecting by position rather than index keeps a caller's request stable
/// across a topology change that renumbers edges.
#[derive(Debug, Clone, PartialEq)]
pub enum EdgeSelector {
    /// The vertical edge at the profile corner nearest this point.
    NearestCorner(Point2),
    /// The edge carrying this persistent structural name.
    ///
    /// Unlike [`Self::NearestCorner`], this does not depend on the caller
    /// knowing where the edge currently is. A name obtained from a previous
    /// result still selects the same edge after the profile has been resized
    /// or the solid rebuilt, which is what re-applying a feature requires.
    Named(EdgeName),
}

impl EdgeSelector {
    /// Resolve this selector to a profile corner index.
    ///
    /// Both variants answer the same question -- which corner of the profile
    /// ring -- so the construction below stays one code path.
    fn corner_index(&self, corners: &[Point2]) -> GeomResult<usize> {
        match self {
            Self::NearestCorner(target) => corners
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    (**a - *target)
                        .length_squared()
                        .total_cmp(&(**b - *target).length_squared())
                })
                .map(|(index, _)| index)
                .ok_or_else(|| GeomError::Degenerate("profile has no corners".to_owned())),
            // A corner edge is where two consecutive walls meet. The name
            // states that pair, so resolving it is a search for the corner
            // whose adjacent side ordinals match -- no geometry involved, which
            // is exactly why it survives an edit that moves the corner.
            Self::Named(name) => {
                let count = corners.len();
                (0..count)
                    .find(|index| {
                        let previous = (index + count - 1) % count;
                        let (Ok(ordinal), Ok(previous_ordinal)) =
                            (u32::try_from(*index), u32::try_from(previous))
                        else {
                            return false;
                        };
                        let candidate = EdgeName::between(
                            FaceName::swept(SweptFace::Side(previous_ordinal)),
                            FaceName::swept(SweptFace::Side(ordinal)),
                        );
                        candidate == *name
                    })
                    .ok_or_else(|| {
                        unsupported("edge name does not resolve to a corner of this profile")
                    })
            }
        }
    }
}

/// How much material a feature removes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FeatureSize {
    /// Constant setback measured along both faces meeting at the edge.
    ConstantDistance(Scalar),
    /// Constant blend radius. Representable, not yet constructible here.
    ConstantRadius(Scalar),
}

/// Chamfer one vertical edge of an exact extruded rectangle.
///
/// Supported: a sharp filled rectangle extruded along +z, chamfered at one
/// corner by a constant distance. Everything else is refused by name.
pub fn chamfer_extruded_profile(
    profile: &Profile,
    direction: Vec3,
    depth: Scalar,
    edge: EdgeSelector,
    size: FeatureSize,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let distance = match size {
        FeatureSize::ConstantDistance(distance) => distance,
        // A fillet needs a cylindrical wall stitched into the prism assembly,
        // which is a different construction. Approximating it with a
        // multi-segment chamfer would be indistinguishable to the caller.
        FeatureSize::ConstantRadius(_) => {
            return Err(unsupported("constant-radius fillet on an extruded solid"))
        }
    };
    if !distance.is_finite() || distance <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "chamfer distance must be positive and finite, got {distance}"
        )));
    }

    let Profile::Rectangle(rectangle) = profile else {
        return Err(unsupported("chamfer on a non-rectangle profile"));
    };
    let RectangleProfile {
        x,
        y,
        thickness,
        outer_radius,
        inner_radius,
    } = *rectangle;
    if thickness.is_some() {
        return Err(unsupported("chamfer on a hollow profile"));
    }
    if outer_radius.is_some() || inner_radius.is_some() {
        return Err(unsupported("chamfer on an already-rounded profile"));
    }
    if !x.is_finite() || !y.is_finite() || x <= 0.0 || y <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "chamfer profile must have positive finite extents, got {x} x {y}"
        )));
    }
    if !depth.is_finite() || depth <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "chamfer extrusion depth must be positive and finite, got {depth}"
        )));
    }
    // Only a +z extrusion keeps the vertical edges vertical, which is what
    // makes the chamfered result another prism.
    if direction.normalize_or_zero().dot(Vec3::Z) < 1.0 - tolerance.linear() {
        return Err(unsupported("chamfer on an oblique extrusion"));
    }

    let (half_x, half_y) = (x / 2.0, y / 2.0);
    // Counter-clockwise, matching the outer ring extrude_rectangle builds.
    let corners = [
        Point2::new(-half_x, -half_y),
        Point2::new(half_x, -half_y),
        Point2::new(half_x, half_y),
        Point2::new(-half_x, half_y),
    ];

    // The setback is measured along each adjacent edge, so it cannot reach
    // the neighbouring corner: at exactly half an edge the two chamfers meet
    // and the edge vanishes, which is a different topology.
    if distance * 2.0 >= x || distance * 2.0 >= y {
        return Err(GeomError::Degenerate(format!(
            "chamfer distance {distance} consumes an entire edge of the {x} x {y} profile"
        )));
    }

    let index = edge.corner_index(&corners)?;

    // Replace the corner with two points, each set back along one adjacent
    // edge. The corner's own vertex disappears: that is the chamfer.
    let previous = corners[(index + corners.len() - 1) % corners.len()];
    let corner = corners[index];
    let next = corners[(index + 1) % corners.len()];
    let into_previous = (previous - corner).normalize();
    let into_next = (next - corner).normalize();

    let mut ring = Vec::with_capacity(corners.len() + 1);
    for (position, point) in corners.iter().enumerate() {
        if position == index {
            // Order matters: entering along the previous edge, leaving along
            // the next one, so the ring keeps its counter-clockwise winding.
            ring.push(corner + into_previous * distance);
            ring.push(corner + into_next * distance);
        } else {
            ring.push(*point);
        }
    }

    extrude_polygon_rings(&[ring], Vec3::Z * depth)
}

/// Fillet one vertical edge of an exact extruded rectangle.
///
/// The blend is a real cylindrical face, not a many-segment approximation.
/// A vertical edge of a +z prism sweeps its cross-section unchanged, so the
/// blend surface is `Surface::Cylinder` with its axis parallel to z, stitched
/// between the two planar walls it is tangent to.
///
/// # Tangency is constructed, not asserted
///
/// The blend axis sits on the internal angle bisector at distance
/// `r / sin(theta/2)` from the corner, where `theta` is the interior angle.
/// At that distance the perpendicular from the axis to each adjacent wall is
/// exactly `r`, so the cylinder meets both walls tangentially by
/// construction. For the right angle of a rectangle this is `r * sqrt(2)`.
/// Nothing is checked afterwards, so nothing can drift.
///
/// # Errors
///
/// Refuses a radius that reaches past either neighbouring corner, and the
/// same input families the chamfer refuses.
pub fn fillet_extruded_profile(
    profile: &Profile,
    direction: Vec3,
    depth: Scalar,
    edge: EdgeSelector,
    size: FeatureSize,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let radius = match size {
        FeatureSize::ConstantRadius(radius) => radius,
        FeatureSize::ConstantDistance(_) => {
            return Err(unsupported(
                "chamfer requested through the fillet entry point",
            ))
        }
    };
    if !radius.is_finite() || radius <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "fillet radius must be positive and finite, got {radius}"
        )));
    }

    let geometry = rectangle_prism(profile, direction, depth, tolerance)?;
    let (x, y) = geometry;

    // The tangent points sit `radius` back along each adjacent edge, so a
    // radius of half an edge reaches the neighbouring corner and the wall
    // between them vanishes -- a different topology, not a fillet.
    if radius * 2.0 >= x || radius * 2.0 >= y {
        return Err(GeomError::Degenerate(format!(
            "fillet radius {radius} reaches past a neighbouring corner of the {x} x {y} profile"
        )));
    }

    build_filleted_prism(x, y, depth, radius, edge)
}

/// Validate the profile family the fillet supports, returning its extents.
///
/// Same refusals as the chamfer, for the same reasons: a hollow or
/// already-rounded profile is a different corner, and an oblique extrusion
/// does not keep the vertical edges vertical, which is what makes the blend
/// a cylinder rather than a general swept surface.
fn rectangle_prism(
    profile: &Profile,
    direction: Vec3,
    depth: Scalar,
    tolerance: Tolerance,
) -> GeomResult<(Scalar, Scalar)> {
    let Profile::Rectangle(rectangle) = profile else {
        return Err(unsupported("fillet on a non-rectangle profile"));
    };
    let RectangleProfile {
        x,
        y,
        thickness,
        outer_radius,
        inner_radius,
    } = *rectangle;
    if thickness.is_some() {
        return Err(unsupported("fillet on a hollow profile"));
    }
    if outer_radius.is_some() || inner_radius.is_some() {
        return Err(unsupported("fillet on an already-rounded profile"));
    }
    if !x.is_finite() || !y.is_finite() || x <= 0.0 || y <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "fillet profile must have positive finite extents, got {x} x {y}"
        )));
    }
    if !depth.is_finite() || depth <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "fillet extrusion depth must be positive and finite, got {depth}"
        )));
    }
    if direction.normalize_or_zero().dot(Vec3::Z) < 1.0 - tolerance.linear() {
        return Err(unsupported("fillet on an oblique extrusion"));
    }
    Ok((x, y))
}

/// Corner geometry of the fillet: where the blend starts, ends, and centres.
///
/// For a right-angled rectangle corner the interior angle is pi/2, so the
/// centre lies on the bisector at `r / sin(pi/4)` = `r * sqrt(2)`. The
/// tangent points are the feet of the perpendiculars from that centre, which
/// land exactly `r` back along each adjacent edge.
pub(crate) struct BlendCorner {
    pub(crate) centre: Point2,
    pub(crate) start: Point2,
    pub(crate) end: Point2,
    pub(crate) sweep: Scalar,
}

/// Solve the blend for one corner of a closed polygon ring.
///
/// # Why the general angle, not the right-angle shortcut
///
/// The rectangle-only version summed the two tangent offsets, which lands on
/// the centre only when the edges meet at a right angle. For an arbitrary
/// interior angle `theta` the setback along each edge is `r / tan(theta/2)`
/// and the centre sits `r / sin(theta/2)` along the bisector. Both reduce to
/// the old expressions at `theta = pi/2`, so rectangles are unchanged.
///
/// Returns `None` when the corner is degenerate -- collinear or reversed
/// edges have no well-defined bisector, and a blend there is not geometry
/// this can name.
fn blend_corner(corners: &[Point2], index: usize, radius: Scalar) -> Option<BlendCorner> {
    let count = corners.len();
    let previous = corners[(index + count - 1) % count];
    let corner = corners[index];
    let next = corners[(index + 1) % count];

    let into_previous = (previous - corner).normalize_or_zero();
    let into_next = (next - corner).normalize_or_zero();
    if into_previous == Vec2::ZERO || into_next == Vec2::ZERO {
        return None;
    }

    // Interior angle at the corner.
    let cosine = into_previous.dot(into_next).clamp(-1.0, 1.0);
    let theta = cosine.acos();
    let half = theta / 2.0;
    let (sin_half, tan_half) = (half.sin(), half.tan());
    // Collinear (theta = pi) or folded back (theta = 0) has no blend.
    if sin_half.abs() <= f64::EPSILON || tan_half.abs() <= f64::EPSILON {
        return None;
    }

    // Tangent points sit `r / tan(theta/2)` back along each edge, which is
    // where the perpendicular from the centre meets the wall.
    let setback = radius / tan_half;
    let start = corner + into_previous * setback;
    let end = corner + into_next * setback;

    // The centre lies on the interior bisector at `r / sin(theta/2)`.
    let bisector = (into_previous + into_next).normalize_or_zero();
    if bisector == Vec2::ZERO {
        return None;
    }
    let centre = corner + bisector * (radius / sin_half);

    // Sweep is measured between the two tangent directions; the frame's own
    // x-axis is built from the start point, so no absolute start angle is
    // needed.
    let start_angle = (start.y - centre.y).atan2(start.x - centre.x);
    let end_angle = (end.y - centre.y).atan2(end.x - centre.x);
    // Shortest signed sweep, which for a convex corner is the minor arc.
    let mut sweep = end_angle - start_angle;
    while sweep > core::f64::consts::PI {
        sweep -= core::f64::consts::TAU;
    }
    while sweep < -core::f64::consts::PI {
        sweep += core::f64::consts::TAU;
    }

    Some(BlendCorner {
        centre,
        start,
        end,
        sweep,
    })
}

/// Build the filleted prism: planar walls plus one cylindrical blend face.
///
/// The ring is the rectangle with the filleted corner replaced by its two
/// tangent points, so the planar walls already stop exactly where the blend
/// begins. `extrude_polygon_rings` builds those walls and both caps; the
/// blend is then the one face spanning the gap, and it is a genuine
/// `Surface::Cylinder` rather than a fan of narrow planes.
fn build_filleted_prism(
    x: Scalar,
    y: Scalar,
    depth: Scalar,
    radius: Scalar,
    edge: EdgeSelector,
) -> GeomResult<ExactBRep> {
    let (half_x, half_y) = (x / 2.0, y / 2.0);
    let corners = [
        Point2::new(-half_x, -half_y),
        Point2::new(half_x, -half_y),
        Point2::new(half_x, half_y),
        Point2::new(-half_x, half_y),
    ];

    let index = edge.corner_index(&corners)?;

    let blend = blend_corner(&corners, index, radius)
        .ok_or_else(|| unsupported("fillet at a degenerate corner"))?;

    // The ring carries the tangent points in place of the sharp corner, in
    // winding order: enter along the previous edge, leave along the next.
    let mut ring = Vec::with_capacity(corners.len() + 1);
    for (position, point) in corners.iter().enumerate() {
        if position == index {
            ring.push(blend.start);
            ring.push(blend.end);
        } else {
            ring.push(*point);
        }
    }

    extrude_with_cylindrical_blend(&ring, Vec3::Z * depth, index, &blend, radius)
}

/// Fillet one corner of an arbitrary closed polygon profile.
///
/// # Why this exists next to [`fillet_extruded_profile`]
///
/// The profile-based entry point takes a `Profile::Rectangle` and so can
/// only ever round a four-corner box. This one takes the ring directly, so
/// an L-shape, a hexagon, or any polygon a boolean produced can be
/// filleted. Corners are addressed by index into `ring`, which is stable
/// under a rebuild in a way an arena id is not.
///
/// # What is still refused
///
/// - A radius whose setback does not fit on either adjacent edge. Trimming
///   past a neighbouring corner would silently delete that corner.
/// - A degenerate corner: collinear or folded edges have no bisector.
/// - A reflex corner. The blend there is a fillet on the outside of the
///   material, which is a different surface and is not derived here.
pub fn fillet_polygon_corner(
    ring: &[Point2],
    corner: usize,
    radius: Scalar,
    depth: Scalar,
) -> GeomResult<ExactBRep> {
    if ring.len() < 3 {
        return Err(GeomError::InvalidInput(
            "fillet needs a ring of at least three corners".to_owned(),
        ));
    }
    if corner >= ring.len() {
        return Err(GeomError::InvalidInput(format!(
            "fillet corner {corner} is outside a ring of {} corners",
            ring.len()
        )));
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "fillet radius must be positive and finite, got {radius}"
        )));
    }
    if !depth.is_finite() || depth <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "fillet extrusion depth must be positive and finite, got {depth}"
        )));
    }
    if !ring.iter().all(|p| p.x.is_finite() && p.y.is_finite()) {
        return Err(GeomError::InvalidInput(
            "fillet ring has a non-finite corner".to_owned(),
        ));
    }

    let count = ring.len();
    let previous = ring[(corner + count - 1) % count];
    let here = ring[corner];
    let next = ring[(corner + 1) % count];

    // A blend only makes sense on a convex corner of a counter-clockwise
    // ring. On a reflex corner the arc bulges into the material and is a
    // different surface, so it is refused rather than silently built wrong.
    let ring_is_ccw = signed_area(ring) > 0.0;
    let turn = (here - previous).perp_dot(next - here);
    let convex = if ring_is_ccw { turn > 0.0 } else { turn < 0.0 };
    if !convex {
        return Err(unsupported("fillet on a reflex corner"));
    }

    let blend = blend_corner(ring, corner, radius)
        .ok_or_else(|| unsupported("fillet at a degenerate corner"))?;

    // The setback must fit on both adjacent edges. If it ran past a
    // neighbour the blend would swallow that corner, changing the profile
    // rather than rounding it.
    let setback = (blend.start - here).length();
    let to_previous = (previous - here).length();
    let to_next = (next - here).length();
    if setback >= to_previous || setback >= to_next {
        return Err(unsupported("fillet radius larger than an adjacent edge"));
    }

    // Replace the corner with its two tangent points, so the planar walls
    // already stop exactly where the blend begins.
    let mut blended = Vec::with_capacity(count + 1);
    for (index, point) in ring.iter().enumerate() {
        if index == corner {
            blended.push(blend.start);
            blended.push(blend.end);
        } else {
            blended.push(*point);
        }
    }
    let blend_index = corner;
    extrude_with_cylindrical_blend(&blended, Vec3::Z * depth, blend_index, &blend, radius)
}

/// Twice the signed area of a closed ring; positive is counter-clockwise.
fn signed_area(ring: &[Point2]) -> Scalar {
    let count = ring.len();
    (0..count)
        .map(|index| {
            let a = ring[index];
            let b = ring[(index + 1) % count];
            a.x * b.y - b.x * a.y
        })
        .sum()
}
