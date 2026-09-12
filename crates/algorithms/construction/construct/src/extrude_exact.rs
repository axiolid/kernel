//! Exact analytic linear-extrusion constructors.

use std::f64::consts::TAU;

use axiolid_brep::{EdgeName, ExactBRep, ExactBRepBuilder, FaceName, SweptFace};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{
    audit_brep, Edge, EdgeId, EdgeUse, Face, FaceBound, FaceId, Loop, LoopId, Orientation, Shell,
    Solid, Vertex, VertexId,
};

use crate::BACKEND_ID;

#[derive(Debug)]
struct RingTopology {
    bottom_edges: Vec<EdgeId>,
    top_edges: Vec<EdgeId>,
    vertical_edges: Vec<EdgeId>,
    bottom_points: Vec<Point3>,
}

/// Extrude a supported profile into an exact, closed analytic B-rep.
///
/// Initial exact families are deliberately narrow: sharp filled/hollow
/// rectangles under a forward non-coplanar linear extrusion, and filled circles
/// along the positive profile normal. Every other family returns a typed refusal
/// instead of taking the existing tessellation path.
pub fn extrude_profile_exact(
    profile: &Profile,
    direction: Vec3,
    depth: Scalar,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let offset = extrusion_offset(direction, depth, tolerance)?;
    match profile {
        Profile::Rectangle(rectangle) => extrude_rectangle(rectangle, offset, tolerance),
        Profile::Circle(circle) => extrude_circle(circle, offset),
        Profile::Ellipse(_) => Err(unsupported("ellipse extrusion")),
        Profile::Section(_) => Err(unsupported("section-profile extrusion")),
        Profile::Contour(_) => Err(unsupported("contour extrusion")),
        Profile::Derived { .. } => Err(unsupported("derived-profile extrusion")),
        Profile::Composite(_) => Err(unsupported("composite-profile extrusion")),
        Profile::CenterLine(_) => Err(unsupported("center-line extrusion")),
        _ => Err(unsupported("unknown profile extrusion")),
    }
}

fn normalize_direction(direction: Vec3) -> Option<Vec3> {
    if !direction.is_finite() {
        return None;
    }

    // Scaling first avoids underflow/overflow in the squared norm while
    // preserving every finite nonzero direction, including subnormals.
    let scale = direction.abs().max_element();
    if scale == 0.0 {
        return None;
    }
    let scaled = direction / scale;
    Some(scaled / scaled.length())
}

// Compare IEEE-754 magnitudes so underflow cannot be optimized as real arithmetic.
const MAGNITUDE_BITS: u64 = 0x7fff_ffff_ffff_ffff;

fn loses_nonzero_component(input: Vec3, output: Vec3) -> bool {
    (input.x.to_bits() & MAGNITUDE_BITS != 0 && output.x.to_bits() & MAGNITUDE_BITS == 0)
        || (input.y.to_bits() & MAGNITUDE_BITS != 0 && output.y.to_bits() & MAGNITUDE_BITS == 0)
        || (input.z.to_bits() & MAGNITUDE_BITS != 0 && output.z.to_bits() & MAGNITUDE_BITS == 0)
}

fn extrusion_offset(direction: Vec3, depth: Scalar, tolerance: Tolerance) -> GeomResult<Vec3> {
    if !depth.is_finite() || depth <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "exact extrusion depth must be positive and finite, got {depth}"
        )));
    }
    let normalized = normalize_direction(direction).ok_or_else(|| {
        GeomError::InvalidInput("exact extrusion direction must be finite and non-zero".to_owned())
    })?;
    if loses_nonzero_component(direction, normalized) {
        return Err(GeomError::Degenerate(
            "exact extrusion direction cannot be normalized without losing a nonzero component"
                .to_owned(),
        ));
    }
    let offset = normalized * depth;
    if !offset.is_finite() || loses_nonzero_component(normalized, offset) {
        return Err(GeomError::Degenerate(
            "exact extrusion direction could not be scaled without range loss".to_owned(),
        ));
    }
    if offset.z <= tolerance.linear() {
        return Err(unsupported("non-forward planar extrusion"));
    }
    Ok(offset)
}

fn extrude_rectangle(
    rectangle: &RectangleProfile,
    offset: Vec3,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    if rectangle.outer_radius.is_some() || rectangle.inner_radius.is_some() {
        return Err(unsupported("rounded rectangle extrusion"));
    }
    if !rectangle.x.is_finite()
        || !rectangle.y.is_finite()
        || rectangle.x <= 0.0
        || rectangle.y <= 0.0
    {
        return Err(GeomError::InvalidInput(format!(
            "exact rectangle profile must have positive finite extents, got {} x {}",
            rectangle.x, rectangle.y
        )));
    }

    let (hx, hy) = (rectangle.x / 2.0, rectangle.y / 2.0);
    let mut rings = vec![vec![
        Point2::new(-hx, -hy),
        Point2::new(hx, -hy),
        Point2::new(hx, hy),
        Point2::new(-hx, hy),
    ]];
    if let Some(thickness) = rectangle.thickness {
        if !thickness.is_finite()
            || thickness <= 0.0
            || 2.0 * thickness >= rectangle.x
            || 2.0 * thickness >= rectangle.y
        {
            return Err(GeomError::InvalidInput(format!(
                "exact hollow rectangle wall thickness {thickness} does not fit inside {} x {}",
                rectangle.x, rectangle.y
            )));
        }
        let (ix, iy) = (hx - thickness, hy - thickness);
        if tolerance.eq(ix, 0.0) || tolerance.eq(iy, 0.0) {
            return Err(GeomError::Degenerate(
                "exact hollow rectangle has a collapsed inner ring".to_owned(),
            ));
        }
        // Clockwise inner ring: its side normals point into the through-passage.
        rings.push(vec![
            Point2::new(-ix, -iy),
            Point2::new(-ix, iy),
            Point2::new(ix, iy),
            Point2::new(ix, -iy),
        ]);
    }

    extrude_polygon_rings(&rings, offset)
}

/// Extrude closed polygon rings into an exact prism.
///
/// The ring-to-B-rep assembly is general over polygon rings, so it is shared
/// rather than duplicated: a chamfered box is a prism over a pentagon, and a
/// filleted one differs only in that an edge becomes an arc. Ring 0 is the
/// outer boundary; any further rings are through-holes.
pub(crate) fn extrude_polygon_rings(rings: &[Vec<Point2>], offset: Vec3) -> GeomResult<ExactBRep> {
    extrude_polygon_rings_named(rings, offset, &mut |_| None)
}

/// Extrude polygon rings, naming each side wall from its own geometry.
///
/// `name_side` receives the wall's two base corners and returns the name
/// that wall should carry, or `None` to leave it unnamed. Passing the
/// geometry rather than an index is deliberate: the caller recovering
/// provenance from a boolean has no index to match against, only position.
pub(crate) fn extrude_polygon_rings_named(
    rings: &[Vec<Point2>],
    offset: Vec3,
    name_side: &mut dyn FnMut((Point2, Point2)) -> Option<FaceName>,
) -> GeomResult<ExactBRep> {
    let ring_count = rings.len();
    let edge_count = ring_count * 12;
    let pcurve_count = ring_count * 24;
    let face_count = 2 + ring_count * 4;
    let mut builder = ExactBRepBuilder::default();
    reserve(
        &mut builder,
        ring_count * 8,
        edge_count,
        2 * ring_count + 4 * ring_count,
        face_count,
        edge_count,
        pcurve_count,
        face_count,
    )?;

    let ring_topology = rings
        .iter()
        .map(|ring| add_polygon_ring(&mut builder, ring.as_slice(), offset))
        .collect::<GeomResult<Vec<_>>>()?;

    let bottom_surface = builder.add_surface(Surface::Plane(Plane {
        frame: identity_frame3(Vec3::ZERO),
    }));
    let top_surface = builder.add_surface(Surface::Plane(Plane {
        frame: identity_frame3(offset),
    }));

    let mut bottom_bounds = Vec::with_capacity(ring_count);
    let mut top_bounds = Vec::with_capacity(ring_count);
    for (ring_index, ring) in ring_topology.iter().enumerate() {
        let bottom_loop = add_cap_loop(&mut builder, &ring.bottom_edges, &ring.bottom_points);
        let top_loop = add_cap_loop(&mut builder, &ring.top_edges, &ring.bottom_points);
        bottom_bounds.push(FaceBound {
            loop_id: bottom_loop,
            orientation: Orientation::Forward,
            outer: ring_index == 0,
        });
        top_bounds.push(FaceBound {
            loop_id: top_loop,
            orientation: Orientation::Forward,
            outer: ring_index == 0,
        });
    }

    let bottom_face = builder.topology_mut().add_face(Face {
        surface: Some(bottom_surface),
        bounds: bottom_bounds,
        orientation: Orientation::Reversed,
    });
    let top_face = builder.topology_mut().add_face(Face {
        surface: Some(top_surface),
        bounds: top_bounds,
        orientation: Orientation::Forward,
    });

    let mut faces = vec![bottom_face, top_face];
    for (ring_index, ring) in ring_topology.iter().enumerate() {
        let source = &rings[ring_index];
        for edge_index in 0..ring.bottom_edges.len() {
            let face = add_planar_side(&mut builder, ring, edge_index, offset)?;
            // The wall's own base corners are what the caller matches
            // against its inputs, so they are what gets handed over.
            let start = source[edge_index];
            let end = source[(edge_index + 1) % source.len()];
            if let Some(name) = name_side((start, end)) {
                builder.set_face_name(face, name);
            }
            faces.push(face);
        }
    }
    finish_closed(builder, faces)
}

fn add_polygon_ring(
    builder: &mut ExactBRepBuilder,
    ring: &[Point2],
    offset: Vec3,
) -> GeomResult<RingTopology> {
    let bottom_points: Vec<_> = ring
        .iter()
        .map(|point| Point3::new(point.x, point.y, 0.0))
        .collect();
    let top_points: Vec<_> = bottom_points.iter().map(|point| *point + offset).collect();
    let bottom_vertices: Vec<_> = bottom_points
        .iter()
        .map(|&position| builder.topology_mut().add_vertex(Vertex { position }))
        .collect();
    let top_vertices: Vec<_> = top_points
        .iter()
        .map(|&position| builder.topology_mut().add_vertex(Vertex { position }))
        .collect();

    let mut bottom_edges = Vec::with_capacity(ring.len());
    let mut top_edges = Vec::with_capacity(ring.len());
    let mut vertical_edges = Vec::with_capacity(ring.len());
    for index in 0..ring.len() {
        let next = (index + 1) % ring.len();
        bottom_edges.push(add_line_edge(
            builder,
            bottom_vertices[index],
            bottom_vertices[next],
            bottom_points[index],
            bottom_points[next] - bottom_points[index],
        ));
        top_edges.push(add_line_edge(
            builder,
            top_vertices[index],
            top_vertices[next],
            top_points[index],
            top_points[next] - top_points[index],
        ));
        vertical_edges.push(add_line_edge(
            builder,
            bottom_vertices[index],
            top_vertices[index],
            bottom_points[index],
            offset,
        ));
    }
    Ok(RingTopology {
        bottom_edges,
        top_edges,
        vertical_edges,
        bottom_points,
    })
}

pub(crate) fn add_line_edge(
    builder: &mut ExactBRepBuilder,
    start: VertexId,
    end: VertexId,
    origin: Point3,
    direction: Vec3,
) -> EdgeId {
    let curve = builder.add_curve3(Curve3::Line(Line3 { origin, direction }));
    let edge = builder.topology_mut().add_edge(Edge {
        start,
        end,
        curve: Some(curve),
    });
    builder.set_edge_interval(edge, Interval::UNIT);
    edge
}

fn add_cap_loop(builder: &mut ExactBRepBuilder, edges: &[EdgeId], points: &[Point3]) -> LoopId {
    let mut uses = Vec::with_capacity(edges.len());
    let mut intervals = Vec::with_capacity(edges.len());
    for (index, &edge) in edges.iter().enumerate() {
        let next = (index + 1) % edges.len();
        let origin = points[index].truncate();
        let direction = (points[next] - points[index]).truncate();
        let pcurve = builder.add_curve2(Curve2::Line(Line2 { origin, direction }));
        uses.push(EdgeUse {
            edge,
            orientation: Orientation::Forward,
            pcurve: Some(pcurve),
        });
        intervals.push(Interval::UNIT);
    }
    add_loop(builder, uses, intervals)
}

fn add_planar_side(
    builder: &mut ExactBRepBuilder,
    ring: &RingTopology,
    index: usize,
    offset: Vec3,
) -> GeomResult<FaceId> {
    let next = (index + 1) % ring.bottom_edges.len();
    let origin = ring.bottom_points[index];
    let edge_vector = ring.bottom_points[next] - origin;
    let x = edge_vector.normalize();
    let z = edge_vector.cross(offset).normalize();
    let y = z.cross(x);
    if !(x.is_finite() && y.is_finite() && z.is_finite()) {
        return Err(GeomError::Degenerate(
            "exact rectangle extrusion produced a degenerate side plane".to_owned(),
        ));
    }
    let surface = builder.add_surface(Surface::Plane(Plane {
        frame: Frame3 { origin, x, y, z },
    }));
    let edge_length = edge_vector.length();
    let offset_uv = Vec2::new(offset.dot(x), offset.dot(y));
    let pcurves = [
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(edge_length, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(edge_length, 0.0),
            direction: offset_uv,
        }),
        Curve2::Line(Line2 {
            origin: offset_uv,
            direction: Vec2::new(edge_length, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: offset_uv,
        }),
    ];
    let edge_uses = [
        (ring.bottom_edges[index], Orientation::Forward),
        (ring.vertical_edges[next], Orientation::Forward),
        (ring.top_edges[index], Orientation::Reversed),
        (ring.vertical_edges[index], Orientation::Reversed),
    ];
    let mut uses = Vec::with_capacity(4);
    let mut intervals = Vec::with_capacity(4);
    for ((edge, orientation), pcurve) in edge_uses.into_iter().zip(pcurves) {
        let pcurve = builder.add_curve2(pcurve);
        uses.push(EdgeUse {
            edge,
            orientation,
            pcurve: Some(pcurve),
        });
        intervals.push(match orientation {
            Orientation::Forward => Interval::UNIT,
            Orientation::Reversed => Interval::new(1.0, 0.0),
        });
    }
    let loop_id = add_loop(builder, uses, intervals);
    Ok(builder.topology_mut().add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    }))
}

fn extrude_circle(circle: &CircleProfile, offset: Vec3) -> GeomResult<ExactBRep> {
    if circle.thickness.is_some() {
        return Err(unsupported("annular circle extrusion"));
    }
    if !circle.radius.is_finite() || circle.radius <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "exact circle profile radius must be positive and finite, got {}",
            circle.radius
        )));
    }
    if offset.x != 0.0 || offset.y != 0.0 {
        return Err(unsupported("oblique circle extrusion"));
    }

    let depth = offset.z;
    let frame_bottom = identity_frame3(Vec3::ZERO);
    let frame_top = identity_frame3(offset);
    let frame2 = Frame2 {
        origin: Vec2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    };
    let bottom_point = Vec3::new(circle.radius, 0.0, 0.0);
    let top_point = bottom_point + offset;

    let mut builder = ExactBRepBuilder::default();
    reserve(&mut builder, 2, 3, 3, 3, 3, 8, 3)?;
    let bottom_vertex = builder.topology_mut().add_vertex(Vertex {
        position: bottom_point,
    });
    let top_vertex = builder.topology_mut().add_vertex(Vertex {
        position: top_point,
    });

    let bottom_curve = builder.add_curve3(Curve3::Circle(Circle3 {
        frame: frame_bottom,
        radius: circle.radius,
    }));
    let bottom_edge = builder.topology_mut().add_edge(Edge {
        start: bottom_vertex,
        end: bottom_vertex,
        curve: Some(bottom_curve),
    });
    builder.set_edge_interval(bottom_edge, Interval::new(0.0, TAU));

    let top_curve = builder.add_curve3(Curve3::Circle(Circle3 {
        frame: frame_top,
        radius: circle.radius,
    }));
    let top_edge = builder.topology_mut().add_edge(Edge {
        start: top_vertex,
        end: top_vertex,
        curve: Some(top_curve),
    });
    builder.set_edge_interval(top_edge, Interval::new(0.0, TAU));

    let seam_curve = builder.add_curve3(Curve3::Line(Line3 {
        origin: bottom_point,
        direction: offset,
    }));
    let seam_edge = builder.topology_mut().add_edge(Edge {
        start: bottom_vertex,
        end: top_vertex,
        curve: Some(seam_curve),
    });
    builder.set_edge_interval(seam_edge, Interval::UNIT);

    let bottom_loop = add_circle_cap_loop(&mut builder, bottom_edge, frame2, circle.radius);
    let top_loop = add_circle_cap_loop(&mut builder, top_edge, frame2, circle.radius);

    let side_curves = [
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::X,
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(TAU, 0.0),
            direction: Vec2::new(0.0, depth),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(0.0, depth),
            direction: Vec2::X,
        }),
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(0.0, depth),
        }),
    ];
    let side_edges = [
        (bottom_edge, Orientation::Forward, Interval::new(0.0, TAU)),
        (seam_edge, Orientation::Forward, Interval::UNIT),
        (top_edge, Orientation::Reversed, Interval::new(TAU, 0.0)),
        (seam_edge, Orientation::Reversed, Interval::new(1.0, 0.0)),
    ];
    let mut side_uses = Vec::with_capacity(4);
    let mut side_intervals = Vec::with_capacity(4);
    for ((edge, orientation, interval), curve) in side_edges.into_iter().zip(side_curves) {
        let pcurve = builder.add_curve2(curve);
        side_uses.push(EdgeUse {
            edge,
            orientation,
            pcurve: Some(pcurve),
        });
        side_intervals.push(interval);
    }
    let side_loop = add_loop(&mut builder, side_uses, side_intervals);

    let bottom_surface = builder.add_surface(Surface::Plane(Plane {
        frame: frame_bottom,
    }));
    let top_surface = builder.add_surface(Surface::Plane(Plane { frame: frame_top }));
    let side_surface = builder.add_surface(Surface::Cylinder(Cylinder {
        frame: frame_bottom,
        radius: circle.radius,
    }));
    let bottom_face = add_single_bound_face(
        &mut builder,
        bottom_surface,
        bottom_loop,
        Orientation::Reversed,
    );
    let top_face = add_single_bound_face(&mut builder, top_surface, top_loop, Orientation::Forward);
    let side_face =
        add_single_bound_face(&mut builder, side_surface, side_loop, Orientation::Forward);
    finish_closed(builder, vec![bottom_face, top_face, side_face])
}

fn add_circle_cap_loop(
    builder: &mut ExactBRepBuilder,
    edge: EdgeId,
    frame: Frame2,
    radius: Scalar,
) -> LoopId {
    let pcurve = builder.add_curve2(Curve2::Circle(Circle2 { frame, radius }));
    add_loop(
        builder,
        vec![EdgeUse {
            edge,
            orientation: Orientation::Forward,
            pcurve: Some(pcurve),
        }],
        vec![Interval::new(0.0, TAU)],
    )
}

pub(crate) fn add_single_bound_face(
    builder: &mut ExactBRepBuilder,
    surface: axiolid_brep::SurfaceId,
    loop_id: LoopId,
    orientation: Orientation,
) -> FaceId {
    builder.topology_mut().add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation,
    })
}

pub(crate) fn add_loop(
    builder: &mut ExactBRepBuilder,
    uses: Vec<EdgeUse<axiolid_brep::Curve2Id>>,
    intervals: Vec<Interval>,
) -> LoopId {
    debug_assert_eq!(uses.len(), intervals.len());
    let loop_id = builder.topology_mut().add_loop(Loop { edges: uses });
    for (use_index, interval) in intervals.into_iter().enumerate() {
        builder.set_pcurve_interval(loop_id, use_index, interval);
    }
    loop_id
}

pub(crate) fn finish_closed(
    mut builder: ExactBRepBuilder,
    faces: Vec<FaceId>,
) -> GeomResult<ExactBRep> {
    let shell = builder.topology_mut().add_shell(Shell {
        faces: faces
            .into_iter()
            .map(|face| (face, Orientation::Forward))
            .collect(),
        closed: true,
    });
    builder.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    let exact = builder
        .finish()
        .map_err(|error| GeomError::BackendContractViolation {
            backend: BACKEND_ID,
            detail: format!("exact extrusion assembly failed: {error}"),
        })?;
    let health = audit_brep(exact.topology());
    if !health.is_closed_manifold() {
        return Err(GeomError::BackendContractViolation {
            backend: BACKEND_ID,
            detail: format!("exact extrusion is not a closed manifold: {health:?}"),
        });
    }
    Ok(exact)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn reserve(
    builder: &mut ExactBRepBuilder,
    vertices: usize,
    edges: usize,
    loops: usize,
    faces: usize,
    curves3: usize,
    curves2: usize,
    surfaces: usize,
) -> GeomResult<()> {
    builder
        .topology_mut()
        .try_reserve(vertices, edges, loops, faces, 1, 1)
        .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
    builder
        .try_reserve(curves3, curves2, surfaces, edges, curves2)
        .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })
}

pub(crate) fn identity_frame3(origin: Point3) -> Frame3 {
    Frame3 {
        origin,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// Extrude a ring where one side is a cylindrical blend instead of a plane.
///
/// Mirrors `extrude_polygon_rings` for the planar walls and caps, then adds
/// the blend as a single face whose surface is `Surface::Cylinder`. The two
/// vertical edges bounding it are shared with the neighbouring planar walls,
/// which is what makes the shell closed and the blend tangent: same edges,
/// same vertices, no seam.
pub(crate) use crate::feature::BlendCorner;

pub(crate) fn extrude_with_cylindrical_blend(
    ring: &[Point2],
    offset: Vec3,
    blend_index: usize,
    blend: &BlendCorner,
    radius: Scalar,
) -> GeomResult<ExactBRep> {
    let mut builder = ExactBRepBuilder::default();
    let n = ring.len();
    reserve(
        &mut builder,
        n * 2,
        n * 3,
        2 + n,
        2 + n,
        n * 3,
        n * 6,
        2 + n,
    )?;

    let topology = add_polygon_ring(&mut builder, ring, offset)?;

    let bottom_surface = builder.add_surface(Surface::Plane(Plane {
        frame: identity_frame3(Vec3::ZERO),
    }));
    let top_surface = builder.add_surface(Surface::Plane(Plane {
        frame: identity_frame3(offset),
    }));
    let bottom_loop = add_cap_loop(
        &mut builder,
        &topology.bottom_edges,
        &topology.bottom_points,
    );
    let top_loop = add_cap_loop(&mut builder, &topology.top_edges, &topology.bottom_points);

    let bottom_face = builder.topology_mut().add_face(Face {
        surface: Some(bottom_surface),
        bounds: vec![FaceBound {
            loop_id: bottom_loop,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Reversed,
    });
    let top_face = builder.topology_mut().add_face(Face {
        surface: Some(top_surface),
        bounds: vec![FaceBound {
            loop_id: top_loop,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    builder.set_face_name(bottom_face, FaceName::swept(SweptFace::StartCap));
    builder.set_face_name(top_face, FaceName::swept(SweptFace::EndCap));

    let mut faces = vec![bottom_face, top_face];
    let side_count = topology.bottom_edges.len();
    for edge_index in 0..side_count {
        // `edge_index` is the profile edge ordinal this wall was swept from.
        // It survives a rebuild with a different face count, which is exactly
        // what an arena index does not.
        let ordinal = u32::try_from(edge_index).map_err(|_| {
            GeomError::Degenerate("profile edge count exceeds u32 capacity".to_owned())
        })?;
        if edge_index == blend_index {
            let face =
                add_cylindrical_blend(&mut builder, &topology, edge_index, offset, blend, radius)?;
            // The blend replaced the edge between the two walls it is tangent
            // to. Naming it after that edge -- rather than after its own
            // position -- is what lets the same fillet be re-applied after an
            // upstream edit: the edge name is still derivable from the profile.
            let previous = (edge_index + side_count - 1) % side_count;
            let previous_ordinal = u32::try_from(previous).map_err(|_| {
                GeomError::Degenerate("profile edge count exceeds u32 capacity".to_owned())
            })?;
            builder.set_face_name(
                face,
                FaceName::blend(EdgeName::between(
                    FaceName::swept(SweptFace::Side(previous_ordinal)),
                    FaceName::swept(SweptFace::Side(ordinal)),
                )),
            );
            faces.push(face);
        } else {
            let face = add_planar_side(&mut builder, &topology, edge_index, offset)?;
            builder.set_face_name(face, FaceName::swept(SweptFace::Side(ordinal)));
            faces.push(face);
        }
    }
    finish_closed(builder, faces)
}

/// The blend face: a cylinder tangent to both neighbouring walls.
///
/// The cylinder's axis is the blend centre extruded along the offset, and its
/// radius is the fillet radius. Tangency holds because the centre was placed
/// at the perpendicular distance `radius` from each adjacent wall -- it is a
/// property of the construction, not something verified here.
///
/// The face is bounded by the same four edges a planar wall would use, so the
/// blend shares vertices and edges with its neighbours and the shell closes.
fn add_cylindrical_blend(
    builder: &mut ExactBRepBuilder,
    ring: &RingTopology,
    index: usize,
    offset: Vec3,
    blend: &BlendCorner,
    radius: Scalar,
) -> GeomResult<FaceId> {
    let next = (index + 1) % ring.bottom_edges.len();
    let centre = Point3::new(blend.centre.x, blend.centre.y, 0.0);
    // Frame x-axis points at the arc start, so the surface parameter u is the
    // angle measured from there and the pcurve intervals are the sweep.
    let start = ring.bottom_points[index];
    let x = (start - centre).normalize();
    let z = offset.normalize();
    let y = z.cross(x);
    if !(x.is_finite() && y.is_finite() && z.is_finite()) {
        return Err(GeomError::Degenerate(
            "fillet blend produced a degenerate cylinder frame".to_owned(),
        ));
    }
    let surface = builder.add_surface(Surface::Cylinder(Cylinder {
        frame: Frame3 {
            origin: centre,
            x,
            y,
            z,
        },
        radius,
    }));

    let height = offset.length();
    let pcurves = [
        // Bottom arc: u sweeps the blend angle at v = 0.
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(blend.sweep, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(blend.sweep, 0.0),
            direction: Vec2::new(0.0, height),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(0.0, height),
            direction: Vec2::new(blend.sweep, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(0.0, height),
        }),
    ];
    let edge_uses = [
        (ring.bottom_edges[index], Orientation::Forward),
        (ring.vertical_edges[next], Orientation::Forward),
        (ring.top_edges[index], Orientation::Reversed),
        (ring.vertical_edges[index], Orientation::Reversed),
    ];

    let mut uses = Vec::with_capacity(4);
    let mut intervals = Vec::with_capacity(4);
    for ((edge, orientation), pcurve) in edge_uses.into_iter().zip(pcurves) {
        let pcurve = builder.add_curve2(pcurve);
        uses.push(EdgeUse {
            edge,
            orientation,
            pcurve: Some(pcurve),
        });
        intervals.push(match orientation {
            Orientation::Forward => Interval::UNIT,
            Orientation::Reversed => Interval::new(1.0, 0.0),
        });
    }
    let loop_id = add_loop(builder, uses, intervals);
    Ok(builder.topology_mut().add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    }))
}
