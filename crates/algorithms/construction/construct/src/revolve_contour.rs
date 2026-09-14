//! General contour revolution (ADR 0059).
//!
//! `revolve_rectangle` builds one shape from four known corners. This module
//! revolves an ARBITRARY closed section, which is what lets `Contour`,
//! `Section`, `CenterLine`, `Derived` and `Composite` revolve: they all lower
//! to a contour already, so the refusals they carried were about this module
//! not existing rather than about the geometry being impossible.
//!
//! # Which surface each segment sweeps
//!
//! Measured, not assumed (`revolve.py`):
//!
//! | profile segment | swept surface |
//! |---|---|
//! | parallel to the axis | cylinder |
//! | perpendicular to the axis | plane (annulus) |
//! | oblique | cone |
//! | circular arc | torus |
//!
//! The cone check fitted `r` as an affine function of `z` to a residual of
//! `2e-13`; the torus check evaluated the implicit form to `1.1e-15`.
//!
//! # Parameterisation
//!
//! Cylinder and cone both parameterise as `(angle, height)`, so a wall
//! between two coaxial circles uses the same four-edge seam loop in both
//! cases. That is why they share one code path here: giving the cone its own
//! loop builder would let the two drift.

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_surface::{Cone, Cylinder, Plane, Surface, Torus};
use axiolid_topology::{
    Edge, EdgeId, EdgeUse, Face, FaceBound, FaceId, Loop, LoopId, Orientation, Shell, Solid,
    Vertex, VertexId,
};

use axiolid_overlay::ArcRing;

use crate::contour_lower::arc_ring_signed_area;
use crate::BACKEND_ID;

const TAU: Scalar = core::f64::consts::TAU;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// One revolved section vertex: its distance from the axis and its height.
#[derive(Debug, Clone, Copy)]
struct Station {
    /// Signed distance from the axis, in the profile plane.
    radius: Scalar,
    /// Position along the axis.
    height: Scalar,
    /// Bulge of the segment LEAVING this station.
    bulge: Scalar,
}
/// Revolve a lowered section a full turn about the profile's local y axis.
///
/// The section must lie entirely on one side of the axis: a section touching
/// or crossing it produces a solid whose walls collapse onto the axis, which
/// is a different topology rather than a degenerate case of this one.
pub fn revolve_arc_ring(
    ring: &ArcRing,
    axis_origin: Point3,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let epsilon = tolerance.linear();
    let count = ring.vertices.len();
    if count < 3 {
        return Err(GeomError::Degenerate(format!(
            "a revolved section needs at least three vertices, got {count}"
        )));
    }

    // Counter-clockwise in the profile plane gives outward wall normals once
    // revolved, matching the extrusion convention.
    let area = arc_ring_signed_area(ring);
    if area == 0.0 {
        return Err(GeomError::Degenerate(
            "revolved section encloses no area".to_owned(),
        ));
    }
    let forward = area > 0.0;

    let stations: Vec<Station> = (0..count)
        .map(|index| {
            let source = if forward { index } else { count - 1 - index };
            let vertex = ring.vertices[source];
            let bulge = if forward {
                vertex.bulge
            } else {
                // Reversing the walk moves each bulge onto the segment that
                // now leaves this station, and flips its side.
                -ring.vertices[(source + count - 1) % count].bulge
            };
            Station {
                radius: vertex.point.x - axis_origin.x,
                height: vertex.point.y,
                bulge,
            }
        })
        .collect();

    let min = stations
        .iter()
        .map(|station| station.radius)
        .fold(Scalar::INFINITY, Scalar::min);
    let max = stations
        .iter()
        .map(|station| station.radius)
        .fold(Scalar::NEG_INFINITY, Scalar::max);
    if min * max < 0.0 || min.abs() <= epsilon || max.abs() <= epsilon {
        return Err(unsupported(
            "exact revolution of a section touching or crossing the axis",
        ));
    }
    // Mirror a section on the negative side so every radius is positive; a
    // full turn makes the two cases the same solid.
    let flip = max < 0.0;
    let stations: Vec<Station> = stations
        .iter()
        .map(|station| Station {
            radius: if flip {
                -station.radius
            } else {
                station.radius
            },
            height: station.height,
            bulge: if flip { -station.bulge } else { station.bulge },
        })
        .collect();

    build(&stations, axis_origin, epsilon)
}
/// Assemble the revolved solid: one circular edge per station, one wall face
/// per segment.
fn build(stations: &[Station], axis_origin: Point3, epsilon: Scalar) -> GeomResult<ExactBRep> {
    let count = stations.len();
    let mut builder = ExactBRepBuilder::default();
    builder
        .topology_mut()
        .try_reserve(count * 2, count * 2, count * 2, count, 1, 1)
        .map_err(|_| GeomError::BudgetExceeded {
            resource: "exact revolution topology",
        })?;

    let axis_frame = |height: Scalar| Frame3 {
        origin: Point3::new(axis_origin.x, height, 0.0),
        x: Vec3::X,
        y: Vec3::Z,
        // The profile plane is z = 0 and the axis is local y, so the surface
        // frames put their own z along the revolution axis.
        z: Vec3::Y,
    };

    // One circular edge per station, plus a seam so each wall loop closes.
    let mut circles = Vec::with_capacity(count);
    let mut vertices = Vec::with_capacity(count);
    let mut points = Vec::with_capacity(count);
    for station in stations {
        let position = Point3::new(axis_origin.x + station.radius, station.height, 0.0);
        let vertex = builder.topology_mut().add_vertex(Vertex { position });
        let curve = builder.add_curve3(Curve3::Circle(Circle3 {
            frame: axis_frame(station.height),
            radius: station.radius,
        }));
        let edge = builder.topology_mut().add_edge(Edge {
            start: vertex,
            end: vertex,
            curve: Some(curve),
        });
        builder.set_edge_interval(edge, Interval::new(0.0, TAU));
        circles.push(edge);
        vertices.push(vertex);
        points.push(position);
    }

    let mut faces = Vec::with_capacity(count);
    for index in 0..count {
        let next = (index + 1) % count;
        let here = stations[index];
        let there = stations[next];
        let face = if here.bulge != 0.0 {
            torus_wall(
                &mut builder,
                axis_origin.x,
                &circles,
                &vertices,
                stations,
                index,
            )?
        } else if (here.radius - there.radius).abs() <= epsilon {
            // Constant radius: a cylinder. A zero-height run would be a
            // degenerate face, so it is skipped rather than emitted.
            if (here.height - there.height).abs() <= epsilon {
                continue;
            }
            straight_wall(
                &mut builder,
                axis_origin.x,
                &circles,
                &vertices,
                &points,
                stations,
                index,
                false,
            )?
        } else if (here.height - there.height).abs() <= epsilon {
            // Constant height: a planar annulus.
            annulus_wall(&mut builder, axis_origin.x, &circles, stations, index)?
        } else {
            // Neither constant: a cone.
            straight_wall(
                &mut builder,
                axis_origin.x,
                &circles,
                &vertices,
                &points,
                stations,
                index,
                true,
            )?
        };
        faces.push(face);
    }

    if faces.len() < 3 {
        return Err(GeomError::Degenerate(format!(
            "a revolved solid needs at least three faces, got {}",
            faces.len()
        )));
    }
    finish(builder, faces)
}
/// A cylinder or cone wall between two coaxial circles.
///
/// Both parameterise as `(angle, height)`, so one seam loop serves both: only
/// the SURFACE differs. The seam is traversed forward and reversed so the
/// loop is vertex-connected without enclosing extra area.
#[allow(clippy::too_many_arguments)]
fn straight_wall(
    builder: &mut ExactBRepBuilder,
    axis_x: Scalar,
    circles: &[EdgeId],
    vertices: &[VertexId],
    points: &[Point3],
    stations: &[Station],
    index: usize,
    cone: bool,
) -> GeomResult<FaceId> {
    let count = stations.len();
    let next = (index + 1) % count;
    let here = stations[index];
    let there = stations[next];
    let rise = there.height - here.height;

    let seam_curve = builder.add_curve3(Curve3::Line(Line3 {
        origin: points[index],
        direction: points[next] - points[index],
    }));
    let seam = builder.topology_mut().add_edge(Edge {
        start: vertices[index],
        end: vertices[next],
        curve: Some(seam_curve),
    });
    builder.set_edge_interval(seam, Interval::UNIT);

    let surface = if cone {
        // `semi_angle` is measured so radius grows as `r + v * tan(angle)`,
        // which is exactly the slope of this segment.
        let semi_angle = ((there.radius - here.radius) / rise).atan();
        builder.add_surface(Surface::Cone(Cone {
            frame: frame_at(axis_x, here.height),
            radius: here.radius,
            semi_angle,
        }))
    } else {
        builder.add_surface(Surface::Cylinder(Cylinder {
            frame: frame_at(axis_x, here.height),
            radius: here.radius,
        }))
    };

    let loop_id = seam_loop(builder, circles[index], circles[next], seam, rise);
    Ok(add_face(builder, surface, loop_id, Orientation::Forward))
}

/// A planar annulus between two circles at the same height.
fn annulus_wall(
    builder: &mut ExactBRepBuilder,
    axis_x: Scalar,
    circles: &[EdgeId],
    stations: &[Station],
    index: usize,
) -> GeomResult<FaceId> {
    let count = stations.len();
    let next = (index + 1) % count;
    let here = stations[index];
    let there = stations[next];
    let surface = builder.add_surface(Surface::Plane(Plane {
        frame: frame_at(axis_x, here.height),
    }));

    // The larger circle bounds the face and the smaller one is its hole, so
    // the two loops carry opposite orientations.
    let (outer, inner) = if here.radius > there.radius {
        (index, next)
    } else {
        (next, index)
    };
    let outer_loop = circle_loop(
        builder,
        circles[outer],
        stations[outer].radius,
        Orientation::Forward,
    );
    let inner_loop = circle_loop(
        builder,
        circles[inner],
        stations[inner].radius,
        Orientation::Reversed,
    );
    let orientation = if here.radius > there.radius {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    Ok(builder.topology_mut().add_face(Face {
        surface: Some(surface),
        bounds: vec![
            FaceBound {
                loop_id: outer_loop,
                orientation: Orientation::Forward,
                outer: true,
            },
            FaceBound {
                loop_id: inner_loop,
                orientation: Orientation::Forward,
                outer: false,
            },
        ],
        orientation,
    }))
}
/// A torus wall: an arc segment revolved about the axis.
///
/// A torus parameterises as `(angle about the axis, angle about the tube)`,
/// so the trim region is the SAME rectangle the cylinder and cone use -- only
/// the `v` range is the arc's own sweep instead of a height.
fn torus_wall(
    builder: &mut ExactBRepBuilder,
    axis_x: Scalar,
    circles: &[EdgeId],
    vertices: &[VertexId],
    stations: &[Station],
    index: usize,
) -> GeomResult<FaceId> {
    let count = stations.len();
    let next = (index + 1) % count;
    let here = stations[index];
    let there = stations[next];

    let arc = crate::extrude_arc::arc_geometry(
        Point2::new(here.radius, here.height),
        Point2::new(there.radius, there.height),
        here.bulge,
    )?;
    let (centre_radius, centre_height) = (arc.centre.x, arc.centre.y);

    // A tube wider than its own ring distance is a SPINDLE torus: the surface
    // passes through the axis and self-intersects, so the solid it bounds is
    // not the one the section describes.
    if arc.radius >= centre_radius.abs() - Scalar::EPSILON {
        return Err(unsupported(
            "revolved arc whose tube reaches the axis (spindle torus)",
        ));
    }
    if centre_radius <= 0.0 {
        return Err(unsupported(
            "revolved arc centred on the far side of the axis",
        ));
    }

    let surface = builder.add_surface(Surface::Torus(Torus {
        frame: frame_at(axis_x, centre_height),
        major_radius: centre_radius,
        minor_radius: arc.radius,
    }));

    // `v` is measured from the outward radial direction in the tube plane.
    let angle_of =
        |station: Station| (station.height - centre_height).atan2(station.radius - centre_radius);
    let start = angle_of(here);
    // Walk the arc in its own direction rather than taking the shorter way
    // round: the sweep says which half of the tube is material.
    let end = start + arc.sweep;

    // The seam runs ACROSS the tube, and on a torus that path is an arc, not
    // a ruling: a straight chord sags below the surface by the sagitta
    // `r*(1 - cos(sweep/2))`, which the audit reports exactly. Cylinders and
    // cones get a straight seam because their v-direction genuinely is a
    // straight ruling; a torus does not.
    let seam_frame = Frame3 {
        origin: Point3::new(axis_x + centre_radius, centre_height, 0.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    };
    let seam_curve = builder.add_curve3(Curve3::Circle(Circle3 {
        frame: seam_frame,
        radius: arc.radius,
    }));
    let seam = builder.topology_mut().add_edge(Edge {
        start: vertices[index],
        end: vertices[next],
        curve: Some(seam_curve),
    });
    builder.set_edge_interval(seam, Interval::new(start, end));

    let loop_id = seam_loop_between(builder, circles[index], circles[next], seam, start, end);
    Ok(add_face(builder, surface, loop_id, Orientation::Forward))
}
// --- shared helpers ---------------------------------------------------------

/// A surface frame on the axis at `height`, with local z along the axis.
///
/// `axis_x` is not optional: every circle in the solid is centred on the axis,
/// so a frame at x = 0 disagrees with them by exactly the axis offset. The
/// audit reports that as a pcurve error equal to `|axis_x|`.
fn frame_at(axis_x: Scalar, height: Scalar) -> Frame3 {
    Frame3 {
        origin: Point3::new(axis_x, height, 0.0),
        x: Vec3::X,
        y: Vec3::Z,
        z: Vec3::Y,
    }
}

/// A closed loop around one circular edge, in the cap plane.
fn circle_loop(
    builder: &mut ExactBRepBuilder,
    edge: EdgeId,
    radius: Scalar,
    orientation: Orientation,
) -> LoopId {
    let pcurve = builder.add_curve2(Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Vec2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    }));
    let loop_id = builder.topology_mut().add_loop(Loop {
        edges: vec![EdgeUse {
            edge,
            orientation,
            pcurve: Some(pcurve),
        }],
    });
    // A reversed edge use walks the curve backwards, so its pcurve interval
    // must run backwards too. Leaving it forward puts the pcurve start
    // diametrically opposite the edge start: the audit reports the mismatch
    // as exactly twice the radius.
    let interval = match orientation {
        Orientation::Forward => Interval::new(0.0, TAU),
        Orientation::Reversed => Interval::new(TAU, 0.0),
    };
    builder.set_pcurve_interval(loop_id, 0, interval);
    loop_id
}

/// The four-edge seam loop for a wall spanning `v` from 0 to `rise`.
fn seam_loop(
    builder: &mut ExactBRepBuilder,
    lower: EdgeId,
    upper: EdgeId,
    seam: EdgeId,
    rise: Scalar,
) -> LoopId {
    seam_loop_between(builder, lower, upper, seam, 0.0, rise)
}

/// The four-edge seam loop for a wall spanning `v` from `from` to `to`.
///
/// A wall between two disjoint circles cannot close with two edges: a loop
/// must be vertex-connected, so the seam is traversed once each way. It adds
/// no area because the two traversals cancel.
fn seam_loop_between(
    builder: &mut ExactBRepBuilder,
    lower: EdgeId,
    upper: EdgeId,
    seam: EdgeId,
    from: Scalar,
    to: Scalar,
) -> LoopId {
    let lower_pcurve = builder.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(0.0, from),
        direction: Vec2::X,
    }));
    let upper_pcurve = builder.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(0.0, to),
        direction: Vec2::X,
    }));
    let up_pcurve = builder.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(TAU, from),
        direction: Vec2::new(0.0, to - from),
    }));
    let down_pcurve = builder.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(0.0, from),
        direction: Vec2::new(0.0, to - from),
    }));
    let loop_id = builder.topology_mut().add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge: lower,
                orientation: Orientation::Forward,
                pcurve: Some(lower_pcurve),
            },
            EdgeUse {
                edge: seam,
                orientation: Orientation::Forward,
                pcurve: Some(up_pcurve),
            },
            EdgeUse {
                edge: upper,
                orientation: Orientation::Reversed,
                pcurve: Some(upper_pcurve),
            },
            EdgeUse {
                edge: seam,
                orientation: Orientation::Reversed,
                pcurve: Some(down_pcurve),
            },
        ],
    });
    builder.set_pcurve_interval(loop_id, 0, Interval::new(0.0, TAU));
    builder.set_pcurve_interval(loop_id, 1, Interval::UNIT);
    builder.set_pcurve_interval(loop_id, 2, Interval::new(TAU, 0.0));
    builder.set_pcurve_interval(loop_id, 3, Interval::new(1.0, 0.0));
    loop_id
}

fn add_face(
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

fn finish(mut builder: ExactBRepBuilder, faces: Vec<FaceId>) -> GeomResult<ExactBRep> {
    let shell = builder.topology_mut().add_shell(Shell {
        faces: faces
            .into_iter()
            .map(|f| (f, Orientation::Forward))
            .collect(),
        closed: true,
    });
    builder.topology_mut().add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    builder.finish().map_err(|error| {
        GeomError::InvalidInput(format!("revolved solid failed validation: {error:?}"))
    })
}
