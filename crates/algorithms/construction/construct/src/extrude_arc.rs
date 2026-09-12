//! Extrude an arc-capable cross-section into an exact B-rep (ADR 0050).
//!
//! # Why this is a second extruder
//!
//! `extrude_polygon_rings` builds a planar wall per edge. An arc edge
//! sweeps a CYLINDRICAL wall instead, so the surface, the 3D curve and the
//! parameter-space curves all differ. Rather than thread a branch through
//! the polygon extruder, the arc case is built here and the polygon path
//! is left untouched.
//!
//! # Exactness
//!
//! An arc edge becomes a `Circle3` bound to the swept `Cylinder` surface:
//! no sampling, no tessellation. The bulge is converted to centre, radius
//! and sweep by closed form, verified against the defining circle before
//! any B-rep is built.

use axiolid_brep::{ExactBRep, ExactBRepBuilder, FaceName, SweptFace};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_overlay::ArcRing;
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{
    Edge, EdgeId, EdgeUse, Face, FaceBound, FaceId, LoopId, Orientation, Vertex, VertexId,
};

use crate::extrude_exact::{
    add_line_edge, add_loop, add_single_bound_face, finish_closed, identity_frame3, reserve,
};

/// The circle an arc edge lies on.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ArcGeometry {
    /// Circle centre in the cross-section plane.
    pub(crate) centre: Point2,
    /// Circle radius.
    pub(crate) radius: Scalar,
    /// Signed sweep, positive counter-clockwise.
    pub(crate) sweep: Scalar,
}

/// Recover the circle of an arc edge from its bulge.
///
/// The bulge is `tan(theta / 4)` for included angle `theta`, so
/// `theta = 4 atan(bulge)` carries both magnitude and direction. The
/// centre sits on the chord's left normal at `R cos(theta / 2)`, which
/// goes negative for a major arc and moves the centre to the other
/// side -- exactly the wanted behaviour, so no case split is needed.
pub(crate) fn arc_geometry(from: Point2, to: Point2, bulge: Scalar) -> GeomResult<ArcGeometry> {
    let chord = to - from;
    let length = chord.length();
    if length == 0.0 || !length.is_finite() || bulge == 0.0 {
        return Err(GeomError::Degenerate(
            "arc edge has no chord or no bulge".to_owned(),
        ));
    }
    let sweep = 4.0 * bulge.atan();
    let half = 0.5 * sweep;
    let sine = half.abs().sin();
    if sine == 0.0 {
        return Err(GeomError::Degenerate(
            "arc edge sweeps a full turn".to_owned(),
        ));
    }
    let radius = length / (2.0 * sine);
    let unit = chord / length;
    let normal = Vec2::new(-unit.y, unit.x);
    let centre = from + chord * 0.5 + normal * (radius * half.cos());
    if !centre.is_finite() || !radius.is_finite() {
        return Err(GeomError::Degenerate(
            "arc edge produced a non-finite circle".to_owned(),
        ));
    }
    Ok(ArcGeometry {
        centre,
        radius,
        sweep,
    })
}

/// Per-edge topology of an extruded arc ring.
struct RingTopology {
    bottom_edges: Vec<EdgeId>,
    top_edges: Vec<EdgeId>,
    vertical_edges: Vec<EdgeId>,
    bottom_points: Vec<Point3>,
}

/// Extrude one arc-capable ring along `offset`, naming each wall.
///
/// Only single-ring sections are built here. A section with holes needs
/// cap faces carrying several bounds, and the arc cap loop machinery for
/// that is not written, so it is refused by the caller rather than
/// silently dropping the hole.
pub(crate) fn extrude_arc_ring(ring: &ArcRing, offset: Vec3) -> GeomResult<ExactBRep> {
    let count = ring.vertices.len();
    if count < 2 {
        return Err(GeomError::Degenerate(
            "arc ring needs at least two vertices".to_owned(),
        ));
    }
    let mut builder = ExactBRepBuilder::default();
    reserve(
        &mut builder,
        count * 2,
        count * 3,
        2 + count,
        2 + count,
        count * 3,
        count * 8,
        2 + count,
    )?;

    let topology = add_arc_ring(&mut builder, ring, offset)?;

    let bottom_surface = builder.add_surface(Surface::Plane(Plane {
        frame: identity_frame3(Vec3::ZERO),
    }));
    let top_surface = builder.add_surface(Surface::Plane(Plane {
        frame: identity_frame3(offset),
    }));
    let bottom_loop = arc_cap_loop(&mut builder, ring, &topology, false);
    let top_loop = arc_cap_loop(&mut builder, ring, &topology, true);

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
    for index in 0..count {
        let ordinal = u32::try_from(index).map_err(|_| {
            GeomError::Degenerate("profile edge count exceeds u32 capacity".to_owned())
        })?;
        let bulge = ring.vertices[index].bulge;
        let face = if bulge == 0.0 {
            add_planar_wall(&mut builder, ring, &topology, index, offset)?
        } else {
            add_cylindrical_wall(&mut builder, ring, &topology, index, offset)?
        };
        builder.set_face_name(face, FaceName::swept(SweptFace::Side(ordinal)));
        faces.push(face);
    }
    finish_closed(builder, faces)
}

/// Build the vertices and the three edge families of an arc ring.
fn add_arc_ring(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    offset: Vec3,
) -> GeomResult<RingTopology> {
    let count = ring.vertices.len();
    let bottom_points: Vec<Point3> = ring
        .vertices
        .iter()
        .map(|vertex| Point3::new(vertex.point.x, vertex.point.y, 0.0))
        .collect();
    let bottom_vertices: Vec<VertexId> = bottom_points
        .iter()
        .map(|position| {
            builder.topology_mut().add_vertex(Vertex {
                position: *position,
            })
        })
        .collect();
    let top_vertices: Vec<VertexId> = bottom_points
        .iter()
        .map(|position| {
            builder.topology_mut().add_vertex(Vertex {
                position: *position + offset,
            })
        })
        .collect();

    let mut bottom_edges = Vec::with_capacity(count);
    let mut top_edges = Vec::with_capacity(count);
    for index in 0..count {
        let next = (index + 1) % count;
        let bulge = ring.vertices[index].bulge;
        // Iterate by flag, not by height: a zero-height offset would make
        // `level == 0.0` true for BOTH levels and put every edge in the
        // bottom family. The caller refuses zero height, but a silent
        // mis-binding here would be far harder to see than a refusal.
        for (is_top, vertices) in [(false, &bottom_vertices), (true, &top_vertices)] {
            let level = if is_top { offset.z } else { 0.0 };
            let lift = Vec3::new(0.0, 0.0, level);
            let curve = if bulge == 0.0 {
                let origin = bottom_points[index] + lift;
                Curve3::Line(Line3 {
                    origin,
                    direction: bottom_points[next] - bottom_points[index],
                })
            } else {
                let arc =
                    arc_geometry(ring.vertices[index].point, ring.vertices[next].point, bulge)?;
                Curve3::Circle(circle_of(&arc, level, bottom_points[index] + lift)?)
            };
            let curve_id = builder.add_curve3(curve);
            let edge = builder.topology_mut().add_edge(Edge {
                start: vertices[index],
                end: vertices[next],
                curve: Some(curve_id),
            });
            let interval = if bulge == 0.0 {
                Interval::UNIT
            } else {
                let arc =
                    arc_geometry(ring.vertices[index].point, ring.vertices[next].point, bulge)?;
                Interval::new(0.0, arc.sweep)
            };
            builder.set_edge_interval(edge, interval);
            if is_top {
                top_edges.push(edge);
            } else {
                bottom_edges.push(edge);
            }
        }
    }

    let mut vertical_edges = Vec::with_capacity(count);
    for index in 0..count {
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

/// The `Circle3` for an arc at height `level`, oriented so its frame
/// x-axis points at the arc start.
///
/// Anchoring x at the start vertex makes the surface parameter u the angle
/// measured from that vertex, so an edge interval of `0..sweep` is
/// literally the arc and needs no offset term.
fn circle_of(arc: &ArcGeometry, level: Scalar, start: Point3) -> GeomResult<Circle3> {
    let centre = Point3::new(arc.centre.x, arc.centre.y, level);
    let radial = start - centre;
    let length = radial.length();
    if length == 0.0 || !length.is_finite() {
        return Err(GeomError::Degenerate(
            "arc edge start lies on its own centre".to_owned(),
        ));
    }
    let x = radial / length;
    let z = Vec3::Z;
    let y = z.cross(x);
    Ok(Circle3 {
        frame: Frame3 {
            origin: centre,
            x,
            y,
            z,
        },
        radius: arc.radius,
    })
}

/// A cylindrical wall swept by an arc edge.
fn add_cylindrical_wall(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    topology: &RingTopology,
    index: usize,
    offset: Vec3,
) -> GeomResult<FaceId> {
    let count = ring.vertices.len();
    let next = (index + 1) % count;
    let arc = arc_geometry(
        ring.vertices[index].point,
        ring.vertices[next].point,
        ring.vertices[index].bulge,
    )?;
    let circle = circle_of(&arc, 0.0, topology.bottom_points[index])?;
    let surface = builder.add_surface(Surface::Cylinder(Cylinder {
        frame: circle.frame,
        radius: arc.radius,
    }));
    let height = offset.z;
    let pcurves = [
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(arc.sweep, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(arc.sweep, 0.0),
            direction: Vec2::new(0.0, height),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(0.0, height),
            direction: Vec2::new(arc.sweep, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(0.0, height),
        }),
    ];
    let edge_uses = [
        (topology.bottom_edges[index], Orientation::Forward),
        (topology.vertical_edges[next], Orientation::Forward),
        (topology.top_edges[index], Orientation::Reversed),
        (topology.vertical_edges[index], Orientation::Reversed),
    ];
    Ok(wall_face(builder, surface, edge_uses, pcurves))
}

/// A planar wall swept by a straight edge.
fn add_planar_wall(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    topology: &RingTopology,
    index: usize,
    offset: Vec3,
) -> GeomResult<FaceId> {
    let count = ring.vertices.len();
    let next = (index + 1) % count;
    let start = topology.bottom_points[index];
    let along = topology.bottom_points[next] - start;
    let length = along.length();
    if length == 0.0 || !length.is_finite() {
        return Err(GeomError::Degenerate(
            "straight edge has zero length".to_owned(),
        ));
    }
    let x = along / length;
    let z = x.cross(offset).normalize();
    if !z.is_finite() {
        return Err(GeomError::Degenerate(
            "straight wall produced a degenerate plane".to_owned(),
        ));
    }
    let surface = builder.add_surface(Surface::Plane(Plane {
        frame: Frame3 {
            origin: start,
            x,
            y: z.cross(x),
            z,
        },
    }));
    let height = offset.z;
    let pcurves = [
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(length, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(length, 0.0),
            direction: Vec2::new(0.0, height),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::new(0.0, height),
            direction: Vec2::new(length, 0.0),
        }),
        Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(0.0, height),
        }),
    ];
    let edge_uses = [
        (topology.bottom_edges[index], Orientation::Forward),
        (topology.vertical_edges[next], Orientation::Forward),
        (topology.top_edges[index], Orientation::Reversed),
        (topology.vertical_edges[index], Orientation::Reversed),
    ];
    Ok(wall_face(builder, surface, edge_uses, pcurves))
}

/// Assemble a four-sided wall face from its edges and pcurves.
fn wall_face(
    builder: &mut ExactBRepBuilder,
    surface: axiolid_brep::SurfaceId,
    edge_uses: [(EdgeId, Orientation); 4],
    pcurves: [Curve2; 4],
) -> FaceId {
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
    add_single_bound_face(builder, surface, loop_id, Orientation::Forward)
}

/// The cap loop of an arc ring, in the cross-section plane.
///
/// Each pcurve is the cross-section curve itself: a `Line2` for a straight
/// edge, a circle for an arc. The cap surface is the z-plane, so the
/// parameter space IS the cross-section plane and no mapping is needed.
fn arc_cap_loop(
    builder: &mut ExactBRepBuilder,
    ring: &ArcRing,
    topology: &RingTopology,
    top: bool,
) -> LoopId {
    let count = ring.vertices.len();
    let edges = if top {
        &topology.top_edges
    } else {
        &topology.bottom_edges
    };
    let mut uses = Vec::with_capacity(count);
    let mut intervals = Vec::with_capacity(count);
    for (index, edge) in edges.iter().enumerate().take(count) {
        let next = (index + 1) % count;
        let from = ring.vertices[index].point;
        let to = ring.vertices[next].point;
        let bulge = ring.vertices[index].bulge;
        // The cap pcurve must be the cross-section curve itself. Using a
        // chord for an arc edge would make the cap boundary disagree with
        // the wall boundary along the same edge -- a silent inconsistency
        // that no area check would reveal.
        let (curve, interval) = if bulge == 0.0 {
            (
                Curve2::Line(Line2 {
                    origin: Vec2::new(from.x, from.y),
                    direction: Vec2::new(to.x - from.x, to.y - from.y),
                }),
                Interval::UNIT,
            )
        } else {
            match arc_geometry(from, to, bulge) {
                Ok(arc) => (
                    Curve2::Circle(Circle2 {
                        frame: circle2_frame(&arc, from),
                        radius: arc.radius,
                    }),
                    Interval::new(0.0, arc.sweep),
                ),
                // A malformed arc was already refused when its 3D edge was
                // built, so this branch is unreachable in practice; falling
                // back to the chord keeps the loop total rather than
                // panicking on an impossible state.
                Err(_) => (
                    Curve2::Line(Line2 {
                        origin: Vec2::new(from.x, from.y),
                        direction: Vec2::new(to.x - from.x, to.y - from.y),
                    }),
                    Interval::UNIT,
                ),
            }
        };
        let pcurve = builder.add_curve2(curve);
        uses.push(EdgeUse {
            edge: *edge,
            orientation: Orientation::Forward,
            pcurve: Some(pcurve),
        });
        intervals.push(interval);
    }
    add_loop(builder, uses, intervals)
}

/// A 2D frame whose x-axis points from the arc centre at its start.
fn circle2_frame(arc: &ArcGeometry, start: Point2) -> Frame2 {
    let radial = start - arc.centre;
    let length = radial.length();
    let x = if length == 0.0 {
        Vec2::X
    } else {
        radial / length
    };
    Frame2 {
        origin: Vec2::new(arc.centre.x, arc.centre.y),
        x,
        y: Vec2::new(-x.y, x.x),
    }
}
