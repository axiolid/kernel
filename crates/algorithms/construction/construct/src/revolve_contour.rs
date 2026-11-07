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

use axiolid_guarantees::Sign;
use axiolid_overlay::ArcRing;

use crate::contour_lower::arc_ring_signed_area;
use crate::revolve::side_of_axis;
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
    /// Exactly on the axis (#379): the exact `side_of_axis` sign is zero.
    on_axis: bool,
}

/// What a section edge sweeps when one of its ends may lie on the axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Contact {
    /// Both ends on the axis, straight: sweeps nothing, so no face.
    Axis,
    /// One end on the axis, straight and exactly perpendicular to it: a disc.
    Cap,
    /// An arc with an end on the axis: a pole, or a horn torus's tangent
    /// point.
    Arc,
    /// A straight edge with an end on the axis that is not perpendicular to
    /// it: a cone's slant up to its apex.
    Cone,
}

/// Revolve a lowered section a full turn about the profile's local y axis.
///
/// The section may touch the axis (#379), classified with the exact
/// `side_of_axis` sign, no tolerance: an edge lying on the axis sweeps no
/// face, and a straight edge ending on it exactly perpendicularly (equal
/// heights) sweeps a disc cap, so a solid cylinder and a stepped shaft
/// revolve. A section with vertices on both sides of the axis is refused as
/// crossing it. A vertex on the axis is admitted when one of its two edges
/// lies on the axis and the other is another axis edge, a disc cap or a
/// straight edge to a cone's apex (an arc is refused until `TODO(#379)`
/// step 3). The apex is a vertex of the B-rep, not a degenerate edge: the
/// cone's wall loop runs `[rim circle], seam, [rim circle], seam reversed`
/// and leaves a gap along the pole line in its parameters, which consumers
/// close (ADR 0059). A vertex on the axis with no axis edge (a pinch point)
/// is refused by name, and so is a vertex off the axis by no more than the
/// tolerance, as neither on it nor clear of it, or a straight edge ending on
/// the axis within the tolerance of perpendicular to it but not exactly.
/// Several separate contacts (a cavity along the axis) give one solid per
/// outermost shell with the nested shells as its voids.
pub fn revolve_arc_ring(
    ring: &ArcRing,
    axis_origin: Point3,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    revolve_ring_about(ring, axis_origin, true, tolerance)
}

/// [`revolve_arc_ring`], told whether the axis direction is exactly the
/// profile's `y` (`x` and `z` zero), which a vertex on the axis needs.
pub(crate) fn revolve_ring_about(
    ring: &ArcRing,
    axis_origin: Point3,
    exact_direction: bool,
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

    let axis_point = Point2::new(axis_origin.x, axis_origin.y);
    let (mut positive, mut negative) = (false, false);
    let stations = (0..count)
        .map(|index| {
            let source = if forward { index } else { count - 1 - index };
            let vertex = ring.vertices[source];
            if !(vertex.point.x.is_finite() && vertex.point.y.is_finite()) {
                return Err(GeomError::InvalidInput(
                    "revolved section must have finite vertices".to_owned(),
                ));
            }
            let bulge = if forward {
                vertex.bulge
            } else {
                // Reversing the walk moves each bulge onto the segment that
                // now leaves this station, and flips its side.
                -ring.vertices[(source + count - 1) % count].bulge
            };
            let side = side_of_axis(axis_point, Vec3::Y, vertex.point);
            match side {
                Sign::Positive => positive = true,
                Sign::Negative => negative = true,
                _ => {}
            }
            Ok(Station {
                radius: vertex.point.x - axis_origin.x,
                height: vertex.point.y,
                bulge,
                on_axis: side == Sign::Zero,
            })
        })
        .collect::<GeomResult<Vec<Station>>>()?;

    if positive && negative {
        return Err(unsupported(
            "exact revolution of a section crossing the axis",
        ));
    }
    if stations
        .iter()
        .any(|station| !station.on_axis && station.radius.abs() <= epsilon)
    {
        return Err(unsupported(
            "exact revolution of a section within tolerance of the axis but not on it to the bit (an arc end computed onto the axis is needed, TODO(#379) step 3)",
        ));
    }
    check_axis_contact(&stations, epsilon)?;
    // The exact sign above is against the line through the origin along y.
    // The axis the caller named is that line only when it is exactly so: a
    // tilt or a z offset within the tolerance moves the real axis off a
    // vertex the sign called on it (#367 refuses the mesh the same way).
    if stations.iter().any(|station| station.on_axis) && !(exact_direction && axis_origin.z == 0.0)
    {
        return Err(unsupported(
            "exact revolution of a section touching an axis that is not exactly the profile's y axis (tilted or off the plane)",
        ));
    }
    // Mirror a section on the negative side so every radius is positive; a
    // full turn makes the two cases the same solid. `side_of_axis` is
    // positive to the left of the upward axis, which is the negative radius.
    // A mirror turns the walk clockwise, so it is walked backwards again:
    // each bulge moves onto the segment that now leaves its station, and
    // mirroring flips it back.
    let stations: Vec<Station> = if positive {
        (0..count)
            .map(|index| {
                let source = count - 1 - index;
                Station {
                    radius: -stations[source].radius,
                    height: stations[source].height,
                    bulge: stations[(source + count - 1) % count].bulge,
                    on_axis: stations[source].on_axis,
                }
            })
            .collect()
    } else {
        stations
    };

    build(&stations, axis_origin, epsilon)
}

/// What the edge leaving station `index` sweeps, judged exactly.
fn contact(stations: &[Station], index: usize) -> Option<Contact> {
    let here = stations[index];
    let there = stations[(index + 1) % stations.len()];
    match (here.on_axis, there.on_axis) {
        (false, false) => None,
        _ if here.bulge != 0.0 => Some(Contact::Arc),
        (true, true) => Some(Contact::Axis),
        _ if here.height == there.height => Some(Contact::Cap),
        _ => Some(Contact::Cone),
    }
}

/// Refuse every vertex on the axis that is not the end of an edge along the
/// axis (a pinch point, a non-manifold vertex), an arc ending on it (until
/// step 3), and a cone's edge so near perpendicular to the axis that it is
/// neither a disc nor a cone to the tolerance.
fn check_axis_contact(stations: &[Station], epsilon: Scalar) -> GeomResult<()> {
    let count = stations.len();
    for (index, station) in stations.iter().enumerate() {
        if !station.on_axis {
            continue;
        }
        let incoming = contact(stations, (index + count - 1) % count);
        let outgoing = contact(stations, index);
        match (incoming, outgoing) {
            (Some(Contact::Arc), _) | (_, Some(Contact::Arc)) => {
                return Err(unsupported(
                    "exact revolution of a section with an arc ending on the axis (pole or horn torus, TODO(#379) step 3)",
                ));
            }
            (Some(Contact::Axis), Some(_)) | (Some(_), Some(Contact::Axis)) => {}
            _ => {
                return Err(unsupported(
                    "exact revolution of a section touching the axis at an isolated vertex (a pinch point)",
                ));
            }
        }
        // A cone's rise is rounding, not geometry, when it is within the
        // tolerance: the edge is neither a disc (equal heights, exactly) nor
        // a cone the audit can hold to the tolerance.
        for (edge, other) in [
            (incoming, (index + count - 1) % count),
            (outgoing, (index + 1) % count),
        ] {
            if edge == Some(Contact::Cone)
                && (station.height - stations[other].height).abs() <= epsilon
            {
                return Err(unsupported(
                    "exact revolution of an edge ending on the axis within the tolerance of perpendicular to it but not exactly",
                ));
            }
        }
    }
    Ok(())
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

    let axis_frame = |height: Scalar| frame_at(axis_origin.x, height);

    // One circular edge per station off the axis, plus a seam so each wall
    // loop closes. A station on the axis traces no circle and has none.
    let mut circles = Vec::with_capacity(count);
    let mut vertices = Vec::with_capacity(count);
    let mut points = Vec::with_capacity(count);
    for (index, station) in stations.iter().enumerate() {
        let position = Point3::new(axis_origin.x + station.radius, station.height, 0.0);
        points.push(position);
        if station.on_axis {
            circles.push(None);
            // A cone's apex is a vertex its seam ends at; a station that only
            // ends axis edges and discs has none.
            let apex = [(index + count - 1) % count, index]
                .into_iter()
                .any(|edge| contact(stations, edge) == Some(Contact::Cone));
            vertices.push(apex.then(|| builder.topology_mut().add_vertex(Vertex { position })));
            continue;
        }
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
        circles.push(Some(edge));
        vertices.push(Some(vertex));
    }

    // Each face with the two stations of the edge it sweeps.
    let mut faces: Vec<(FaceId, usize, usize)> = Vec::with_capacity(count);
    for index in 0..count {
        let next = (index + 1) % count;
        let here = stations[index];
        let there = stations[next];
        let face = if here.on_axis || there.on_axis {
            // `check_axis_contact` admitted only an axis edge, a disc and a
            // cone's slant.
            if here.on_axis && there.on_axis {
                continue;
            }
            if here.height == there.height {
                cap_wall(&mut builder, axis_origin.x, &circles, stations, index)?
            } else {
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
            }
        } else if here.bulge != 0.0 {
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
        faces.push((face, index, next));
    }

    if faces.len() < 2 {
        return Err(GeomError::Degenerate(format!(
            "a revolved solid needs at least two faces, got {}",
            faces.len()
        )));
    }
    finish(builder, stations, faces)
}
/// A cylinder or cone wall between two coaxial circles.
///
/// Both parameterise as `(angle, height)`, so one seam loop serves both: only
/// the SURFACE differs. The seam is traversed forward and reversed so the
/// loop is vertex-connected without enclosing extra area.
///
/// A cone with an end on the axis is framed at its apex, where its radius is
/// exactly zero, so the apex is at `v = 0` and the loop omits the circle
/// there (#379).
#[allow(clippy::too_many_arguments)]
fn straight_wall(
    builder: &mut ExactBRepBuilder,
    axis_x: Scalar,
    circles: &[Option<EdgeId>],
    vertices: &[Option<VertexId>],
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
        start: slot(vertices, index)?,
        end: slot(vertices, next)?,
        curve: Some(seam_curve),
    });
    builder.set_edge_interval(seam, Interval::UNIT);

    // The station the surface is framed at: the apex when an end is on the
    // axis, else the start, as before.
    let basis = if there.on_axis { there } else { here };
    let surface = if cone {
        // `semi_angle` is measured so radius grows as `r + v * tan(angle)`,
        // which is exactly the slope of this segment.
        let semi_angle = ((there.radius - here.radius) / rise).atan();
        builder.add_surface(Surface::Cone(Cone {
            frame: frame_at(axis_x, basis.height),
            radius: basis.radius,
            semi_angle,
        }))
    } else {
        builder.add_surface(Surface::Cylinder(Cylinder {
            frame: frame_at(axis_x, here.height),
            radius: here.radius,
        }))
    };

    let loop_id = seam_loop_between(
        builder,
        circles[index],
        circles[next],
        seam,
        here.height - basis.height,
        there.height - basis.height,
    );
    Ok(add_face(builder, surface, loop_id, Orientation::Forward))
}

/// The circle or vertex of a station, which a station on the axis has none of.
fn slot<T: Copy>(slots: &[Option<T>], index: usize) -> GeomResult<T> {
    slots[index]
        .ok_or_else(|| GeomError::Degenerate("a wall references a station on the axis".to_owned()))
}

/// A planar disc: the straight edge between a station off the axis and one on
/// it, exactly perpendicular to the axis, swept a full turn.
///
/// The one circle of the off-axis station bounds it; it faces the way the
/// annulus does, outward on an edge that runs inward to the axis and the
/// other way on one that leaves it.
fn cap_wall(
    builder: &mut ExactBRepBuilder,
    axis_x: Scalar,
    circles: &[Option<EdgeId>],
    stations: &[Station],
    index: usize,
) -> GeomResult<FaceId> {
    let next = (index + 1) % stations.len();
    let (rim, orientation) = if stations[index].on_axis {
        (next, Orientation::Reversed)
    } else {
        (index, Orientation::Forward)
    };
    let surface = builder.add_surface(Surface::Plane(Plane {
        frame: frame_at(axis_x, stations[rim].height),
    }));
    let loop_id = circle_loop(
        builder,
        slot(circles, rim)?,
        stations[rim].radius,
        Orientation::Forward,
    );
    Ok(add_face(builder, surface, loop_id, orientation))
}

/// A planar annulus between two circles at the same height.
fn annulus_wall(
    builder: &mut ExactBRepBuilder,
    axis_x: Scalar,
    circles: &[Option<EdgeId>],
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
        slot(circles, outer)?,
        stations[outer].radius,
        Orientation::Forward,
    );
    let inner_loop = circle_loop(
        builder,
        slot(circles, inner)?,
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
    circles: &[Option<EdgeId>],
    vertices: &[Option<VertexId>],
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
        start: slot(vertices, index)?,
        end: slot(vertices, next)?,
        curve: Some(seam_curve),
    });
    builder.set_edge_interval(seam, Interval::new(start, end));

    let loop_id = seam_loop_between(
        builder,
        Some(slot(circles, index)?),
        Some(slot(circles, next)?),
        seam,
        start,
        end,
    );
    Ok(add_face(builder, surface, loop_id, Orientation::Forward))
}
// --- shared helpers ---------------------------------------------------------

/// A surface frame on the axis at `height`, with local z along the axis.
///
/// `axis_x` is not optional: every circle in the solid is centred on the axis,
/// so a frame at x = 0 disagrees with them by exactly the axis offset. The
/// audit reports that as a pcurve error equal to `|axis_x|`.
///
/// The profile plane is `z = 0` and the axis is world `y`, so each frame puts
/// its own `z` along the axis. Its `y` is world `-z`, not `+z`: `(X, Z, Y)` is
/// LEFT-handed, and on a left-handed frame a loop running anticlockwise in
/// `(u, v)` runs clockwise about `S_u x S_v`. Every loop here is built
/// anticlockwise in its parameters, so a left-handed frame turned the whole
/// solid inside out -- consistently enough that the topological and
/// geometric audits, which only compare faces with each other, passed it.
/// `exact_properties` measured `-2 pi R A` for every revolution.
pub(crate) fn frame_at(axis_x: Scalar, height: Scalar) -> Frame3 {
    Frame3 {
        origin: Point3::new(axis_x, height, 0.0),
        x: Vec3::X,
        y: Vec3::NEG_Z,
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

/// The seam loop for a wall spanning `v` from `from` to `to`: the lower
/// circle, the seam up, the upper circle back, the seam down.
///
/// A wall between two disjoint circles cannot close with two edges: a loop
/// must be vertex-connected, so the seam is traversed once each way. It adds
/// no area because the two traversals cancel. A circle that is `None` is a
/// pole or apex, a vertex: the loop omits its use and leaves a gap along the
/// pole line in the wall's parameters.
fn seam_loop_between(
    builder: &mut ExactBRepBuilder,
    lower: Option<EdgeId>,
    upper: Option<EdgeId>,
    seam: EdgeId,
    from: Scalar,
    to: Scalar,
) -> LoopId {
    let rise = to - from;
    let mut uses = Vec::with_capacity(4);
    if let Some(lower) = lower {
        uses.push((
            lower,
            Orientation::Forward,
            Line2 {
                origin: Vec2::new(0.0, from),
                direction: Vec2::X,
            },
            Interval::new(0.0, TAU),
        ));
    }
    uses.push((
        seam,
        Orientation::Forward,
        Line2 {
            origin: Vec2::new(TAU, from),
            direction: Vec2::new(0.0, rise),
        },
        Interval::UNIT,
    ));
    if let Some(upper) = upper {
        uses.push((
            upper,
            Orientation::Reversed,
            Line2 {
                origin: Vec2::new(0.0, to),
                direction: Vec2::X,
            },
            Interval::new(TAU, 0.0),
        ));
    }
    uses.push((
        seam,
        Orientation::Reversed,
        Line2 {
            origin: Vec2::new(0.0, from),
            direction: Vec2::new(0.0, rise),
        },
        Interval::new(1.0, 0.0),
    ));
    let (edges, intervals): (Vec<_>, Vec<_>) = uses
        .into_iter()
        .map(|(edge, orientation, line, interval)| {
            let pcurve = builder.add_curve2(Curve2::Line(line));
            (
                EdgeUse {
                    edge,
                    orientation,
                    pcurve: Some(pcurve),
                },
                interval,
            )
        })
        .unzip();
    let loop_id = builder.topology_mut().add_loop(Loop { edges });
    for (use_index, interval) in intervals.into_iter().enumerate() {
        builder.set_pcurve_interval(loop_id, use_index, interval);
    }
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

/// Close the faces into shells, one per connected piece of the revolved
/// section (two stations are connected by an edge between them off the
/// axis), and the shells into solids by nesting.
///
/// A section touching the axis in several places gives several pieces: a
/// cavity along the axis is a piece of its own. The pieces' ranges along the
/// axis (the heights of every station of each of its faces, an apex's
/// included) decide the nesting exactly, from heights taken from the input:
/// one strictly inside another's range is its void. A piece nested twice (a
/// solid inside a cavity) is refused. The faces already face out of the
/// material, so a void's faces face into its cavity as they are.
fn finish(
    mut builder: ExactBRepBuilder,
    stations: &[Station],
    faces: Vec<(FaceId, usize, usize)>,
) -> GeomResult<ExactBRep> {
    struct Piece {
        label: usize,
        faces: Vec<FaceId>,
        low: Scalar,
        high: Scalar,
    }
    let count = stations.len();
    let mut label: Vec<usize> = (0..count).collect();
    for index in 0..count {
        let next = (index + 1) % count;
        if stations[index].on_axis || stations[next].on_axis {
            continue;
        }
        let (a, b) = (label[index], label[next]);
        for l in label.iter_mut() {
            if *l == b {
                *l = a;
            }
        }
    }
    let mut pieces: Vec<Piece> = Vec::new();
    for (index, station) in stations.iter().enumerate() {
        if station.on_axis {
            continue;
        }
        match pieces.iter_mut().find(|piece| piece.label == label[index]) {
            Some(piece) => {
                piece.low = piece.low.min(station.height);
                piece.high = piece.high.max(station.height);
            }
            None => pieces.push(Piece {
                label: label[index],
                faces: Vec::new(),
                low: station.height,
                high: station.height,
            }),
        }
    }
    for (face, index, next) in faces {
        // The edge's station off the axis names its piece.
        let representative = if stations[index].on_axis { next } else { index };
        if let Some(piece) = pieces
            .iter_mut()
            .find(|piece| piece.label == label[representative])
        {
            piece.faces.push(face);
            // An apex or a disc's axis end lies off the piece's off-axis
            // stations; without it a cavity inside a spinning top would not
            // nest.
            for station in [index, next] {
                piece.low = piece.low.min(stations[station].height);
                piece.high = piece.high.max(stations[station].height);
            }
        }
    }
    let contains = |outer: &Piece, inner: &Piece| outer.low < inner.low && inner.high < outer.high;
    let depth: Vec<usize> = pieces
        .iter()
        .map(|piece| pieces.iter().filter(|other| contains(other, piece)).count())
        .collect();
    if depth.iter().any(|&d| d > 1) {
        return Err(unsupported(
            "exact revolution of a section with a solid nested inside a cavity",
        ));
    }
    let shells: Vec<_> = pieces
        .iter()
        .map(|piece| {
            builder.topology_mut().add_shell(Shell {
                faces: piece
                    .faces
                    .iter()
                    .map(|&face| (face, Orientation::Forward))
                    .collect(),
                closed: true,
            })
        })
        .collect();
    for outer in (0..pieces.len()).filter(|&piece| depth[piece] == 0) {
        let voids = (0..pieces.len())
            .filter(|&inner| depth[inner] == 1 && contains(&pieces[outer], &pieces[inner]))
            .map(|inner| shells[inner])
            .collect();
        builder.topology_mut().add_solid(Solid {
            outer: shells[outer],
            voids,
        });
    }
    builder.finish().map_err(|error| {
        GeomError::InvalidInput(format!("revolved solid failed validation: {error:?}"))
    })
}
