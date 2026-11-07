//! Exact revolution for the profile families exact extrusion already covers.
//!
//! # What is exactly representable
//!
//! Revolving a rectangle a full turn about an axis it does not cross produces
//! an annular tube: two cylinders (inner and outer radius) capped by two
//! annular planes. Every one of those is an elementary surface the kernel
//! already carries, so the result is exact -- no tessellation, no sampled
//! approximation.
//!
//! # Partial turns
//!
//! A PARTIAL turn is not this shape. It has two extra planar walls at the
//! start and end angles, its walls are bounded by circular arcs rather than
//! closed by seams, and nothing closes round the axis. That different
//! topology lives in [`crate::revolve_partial`] (#172); every profile that
//! revolves a full turn through the contour path also revolves partially,
//! and a partial turn may additionally touch the axis along straight
//! segments or at vertices.
//!
//! # What is refused, and why that is the honest answer
//!
//! A turn beyond a full turn sweeps through itself and is refused by name;
//! a zero or non-finite angle is invalid input.
//!
//! A FULL turn of a profile crossing the axis sweeps through itself and is
//! refused, as is a partial turn crossing it.
//!
//! A FULL turn of a profile touching the axis (#379) is exact when each
//! vertex on the axis ends an edge lying on the axis, which sweeps no face,
//! and its other edge is another such edge, a straight edge ending on the
//! axis exactly perpendicularly, which sweeps a disc cap, or any other
//! straight edge, which sweeps a cone up to its apex, or an arc ending on
//! the axis with its circle centred on it, which sweeps a sphere up to its
//! pole: a solid cylinder, a stepped shaft, a cone, a bolt, a spinning top, a
//! ball, a domed column. The contact is classified with the exact
//! `side_of_axis` sign, no tolerance; whether an arc's circle is centred on
//! the axis is decided within the rounding of its coordinates (the arc is a
//! bulge, which cannot encode a sphere's meridian exactly). The apex or pole
//! is a vertex of the B-rep, not a degenerate edge (ADR 0059). A vertex on
//! the axis with no edge along it (a pinch point), an arc whose circle is not
//! centred on the axis (a spindle or horn torus, a lemon or an apple), an arc
//! with both ends on the axis and an arc running across it are refused by
//! name, as is a straight edge ending on the axis within the tolerance of
//! perpendicular to it but not exactly.
//!
//! # The axis (#450)
//!
//! An axis exactly along +y or -y builds in place. Any other axis lying
//! exactly in the profile plane (`z` zero in its direction and origin) is
//! built exactly as given: the section is mapped rigidly into the axis's own
//! frame, revolved about +y there and placed back with
//! [`ExactBRep::transformed`], so a tilt of a milliradian is a tilt of a
//! milliradian, never snapped onto +y. A direction with a `z` component, or
//! an origin off the plane, sweeps hyperboloids and is refused by name: there
//! is no tolerance band to offer. Rejected: snapping the axis within the
//! tolerance and reporting it (a tilt of 1.4 mrad over 10 m moves geometry by
//! 14 mm, and construct has no reading channel).
//!
//! The mesh path in `revolve.rs` handles all of these. Refusing here means a
//! caller asking for exactness gets a typed refusal naming the gap instead of
//! a silently tessellated substitute.

use std::f64::consts::TAU;

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Interval, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_guarantees::Sign;
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{
    Edge, EdgeUse, Face, FaceBound, FaceId, Loop, Orientation, Shell, Solid, Vertex,
};

use crate::extrude_exact::extrude_profile_exact;
use crate::revolve::{along_axis, side_of_axis};
use crate::BACKEND_ID;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// Revolve a supported profile into an exact, closed analytic B-rep.
///
/// The axis must lie in the profile plane, exactly: `z` is zero in its
/// direction and its origin. Along y (either sense) it builds in place; any
/// other in-plane direction builds in the axis's own frame and is placed back
/// (#450). `angle` is in radians, by the right-hand rule about
/// `axis_direction`.
///
/// A full turn (`|angle|` within tolerance of `2 pi`) revolves every profile
/// that lowers to a contour, provided it does not cross the axis and touches
/// it only along edges or at perpendicular edge ends (#379): a sharp
/// rectangle becomes an annular tube of two cylinders and two annular caps,
/// a hole becomes a void shell. A partial turn (`0 < |angle| < 2 pi`)
/// revolves the same profiles into a solid capped by the profile at both
/// ends, each hole a tunnel between the caps, and may touch the axis along
/// straight segments (see [`crate::revolve_partial`]).
///
/// Refused by name: an ellipse, a turn beyond a full turn, a profile crossing
/// the axis, one touching it at an isolated vertex (a pinch point), an arc
/// whose circle reaches the axis without being centred on it (spindle or horn
/// torus), an axis out of the profile plane. A zero or non-finite angle is
/// invalid input.
pub fn revolve_profile_exact(
    profile: &Profile,
    axis_origin: Point3,
    axis_direction: Vec3,
    angle: Scalar,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let frame = axis_frame(axis_origin, axis_direction)?;
    if !angle.is_finite() || angle == 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "revolution angle must be finite and non-zero, got {angle}"
        )));
    }
    // More than a full turn sweeps the solid through itself.
    if angle.abs() > TAU + tolerance.linear() {
        return Err(unsupported("exact revolution beyond a full turn"));
    }
    // A partial turn has two extra planar walls and arc-bounded caps: a
    // different topology, not a different parameter (#172).
    if (angle.abs() - TAU).abs() > tolerance.linear() {
        return revolve_partial(profile, &frame, angle, tolerance);
    }

    match profile {
        // The dedicated rectangle builder knows only the axis along y.
        Profile::Rectangle(rectangle) if frame.placement.is_none() => {
            revolve_rectangle(rectangle, &frame, tolerance)
        }
        // A circle lowers to four quarter arcs, each sweeping a torus
        // quarter; a hollow circle's bore becomes a toroidal void (#111).
        Profile::Circle(_) => revolve_via_contour(profile, &frame, None, tolerance),
        Profile::Rectangle(_) => revolve_via_contour(profile, &frame, None, tolerance),
        Profile::Ellipse(_) => Err(unsupported("ellipse exact revolution")),
        // Every one of these lowers to a contour, and a contour revolves.
        // The refusals they carried described a missing module, not missing
        // geometry, so they route through the general path rather than
        // restating a limitation that no longer holds.
        Profile::Section(_)
        | Profile::Contour(_)
        | Profile::Derived { .. }
        | Profile::Composite(_)
        | Profile::CenterLine(_) => revolve_via_contour(profile, &frame, None, tolerance),
        _ => Err(unsupported("unknown profile exact revolution")),
    }
}

/// The axis of a revolution, in the profile plane (#450).
#[derive(Debug, Clone)]
struct AxisFrame {
    /// A point of the axis, as given.
    origin: Point2,
    /// The axis direction as given, in 3D (`z` zero): the exact sign tests
    /// run against it.
    given: Vec3,
    /// The unit direction in the plane, the frame's y.
    unit: Vec2,
    /// Where the frame sits in the profile plane: `None` for an axis exactly
    /// along y, built in place.
    placement: Option<Transform3>,
}

/// The axis `axis_origin` + `axis_direction`, which must lie exactly in the
/// profile plane.
///
/// A direction along +y or -y (`x` exactly zero) is the identity. Any other
/// in-plane direction `d` has the frame `(n, d, z)` with `n = (d.y, -d.x)`,
/// right-handed, whose origin is the axis origin: a section mapped into it is
/// revolved about +y there, and the placement maps the result back.
fn axis_frame(axis_origin: Point3, axis_direction: Vec3) -> GeomResult<AxisFrame> {
    if !(axis_origin.is_finite() && axis_direction.is_finite()) || axis_direction == Vec3::ZERO {
        return Err(GeomError::InvalidInput(
            "revolution axis must be a finite non-zero direction through a finite origin"
                .to_owned(),
        ));
    }
    // A component out of the plane sweeps hyperboloids: no exact family, and
    // no tolerance band to offer.
    if axis_direction.z != 0.0 || axis_origin.z != 0.0 {
        return Err(unsupported(
            "exact revolution about an axis out of the profile plane",
        ));
    }
    let origin = Point2::new(axis_origin.x, axis_origin.y);
    let length = axis_direction.x.hypot(axis_direction.y);
    let unit = Vec2::new(axis_direction.x / length, axis_direction.y / length);
    let placement = (axis_direction.x != 0.0).then(|| {
        let normal = Vec2::new(unit.y, -unit.x);
        Transform3::from_cols(
            Vec3::new(normal.x, normal.y, 0.0),
            Vec3::new(unit.x, unit.y, 0.0),
            Vec3::Z,
            Vec3::new(origin.x, origin.y, 0.0),
        )
    });
    Ok(AxisFrame {
        origin,
        given: axis_direction,
        unit,
        placement,
    })
}

/// `ring` in the axis's frame: `p' = ((p - o) . n, (p - o) . d)`, the bulges
/// kept (a rigid map keeps an arc's bulge).
///
/// The revolution classifies a vertex on the axis with an exact sign against
/// +y, so the map keeps the classification of the original: a vertex exactly
/// on the axis (the exact sign against the axis as given) maps to `x' = 0`,
/// and any other has `x'` of the sign of its side. A vertex whose `x'` rounds
/// to zero or the wrong sign is closer to the axis than the frame can tell
/// and is refused. An edge exactly perpendicular to the axis with one end on
/// it keeps equal heights, which makes a disc; an edge that is not, whose
/// heights come out equal, is refused.
fn into_axis_frame(ring: &ArcRing, frame: &AxisFrame) -> GeomResult<ArcRing> {
    let normal = Vec2::new(frame.unit.y, -frame.unit.x);
    let count = ring.vertices.len();
    let mut points = Vec::with_capacity(count);
    let mut on_axis = Vec::with_capacity(count);
    for vertex in &ring.vertices {
        let p = vertex.point;
        if !(p.x.is_finite() && p.y.is_finite()) {
            return Err(GeomError::InvalidInput(
                "revolved section must have finite vertices".to_owned(),
            ));
        }
        let offset = p - frame.origin;
        let side = side_of_axis(frame.origin, frame.given, p);
        let across = offset.x * normal.x + offset.y * normal.y;
        let along = offset.x * frame.unit.x + offset.y * frame.unit.y;
        // Positive is the left of the axis: negative `x'`.
        let consistent = match side {
            Sign::Zero => true,
            Sign::Positive => across < 0.0,
            Sign::Negative => across > 0.0,
            _ => false,
        };
        if !consistent {
            return Err(unsupported(
                "exact revolution of a section with a vertex closer to the axis than the rounding of the axis frame",
            ));
        }
        on_axis.push(side == Sign::Zero);
        points.push(Point2::new(
            if side == Sign::Zero { 0.0 } else { across },
            along,
        ));
    }
    for index in 0..count {
        let next = (index + 1) % count;
        if on_axis[index] == on_axis[next] || ring.vertices[index].bulge != 0.0 {
            continue;
        }
        let (on, off) = if on_axis[index] {
            (index, next)
        } else {
            (next, index)
        };
        let from = ring.vertices[index].point;
        let to = ring.vertices[next].point;
        if along_axis(frame.given, from, to) == Sign::Zero {
            points[on].y = points[off].y;
        } else if points[on].y == points[off].y {
            return Err(unsupported(
                "exact revolution of an edge ending on the axis within the tolerance of perpendicular to it but not exactly",
            ));
        }
    }
    Ok(ArcRing::new(
        ring.vertices
            .iter()
            .zip(points)
            .map(|(vertex, point)| ArcVertex {
                point,
                bulge: vertex.bulge,
            })
            .collect(),
    ))
}

/// `solid`, built in the axis's frame, placed where the axis is.
fn place_in_axis_frame(solid: ExactBRep, frame: &AxisFrame) -> GeomResult<ExactBRep> {
    match &frame.placement {
        None => Ok(solid),
        Some(placement) => {
            solid
                .transformed(placement)
                .map_err(|error| GeomError::BackendContractViolation {
                    backend: BACKEND_ID,
                    detail: format!("placing a revolution in its axis frame failed: {error}"),
                })
        }
    }
}

/// The dedicated rectangle builder, for an axis exactly along y.
fn revolve_rectangle(
    rectangle: &RectangleProfile,
    frame: &AxisFrame,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let axis_origin = Point3::new(frame.origin.x, frame.origin.y, 0.0);
    if rectangle.thickness.is_some()
        || rectangle.outer_radius.is_some()
        || rectangle.inner_radius.is_some()
    {
        // Rounded corners revolve into tori and a hollow section's opening
        // into a void, both of which the general contour revolver builds;
        // the dedicated path below only knows one straight-edged ring.
        return revolve_via_contour(&Profile::Rectangle(*rectangle), frame, None, tolerance);
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

    let (half_x, half_y) = (rectangle.x / 2.0, rectangle.y / 2.0);
    // Signed offsets of the rectangle's two vertical edges from the axis. The
    // rectangle is centred on the origin, so these are its x extents measured
    // in the axis's frame.
    let left = -half_x - axis_origin.x;
    let right = half_x - axis_origin.x;

    // Straddling the axis means the swept solid is not an annulus: the inner
    // wall collapses onto the axis and the caps stop being annuli. That is a
    // different topology, so refuse rather than emit a degenerate cylinder.
    // An edge exactly on the axis is a solid cylinder, which the contour path
    // builds (#379); crossing it, or lying within the tolerance of it without
    // touching it, is refused there and here.
    if left == 0.0 || right == 0.0 {
        return revolve_via_contour(&Profile::Rectangle(*rectangle), frame, None, tolerance);
    }
    if left * right < 0.0 {
        return Err(unsupported(
            "exact revolution of a section crossing the axis",
        ));
    }
    if left.abs() <= tolerance.linear() || right.abs() <= tolerance.linear() {
        return Err(unsupported(
            "exact revolution of a section within tolerance of the axis but not on it",
        ));
    }
    let inner = left.abs().min(right.abs());
    let outer = left.abs().max(right.abs());

    // The rectangle spans y in [-half_y, half_y] around the axis origin.
    let bottom = axis_origin.y - half_y;
    let top = axis_origin.y + half_y;
    // Right-handed, with `z` along the axis: see `revolve_contour::frame_at`.
    let frame_at = |y: Scalar| crate::revolve_contour::frame_at(axis_origin.x, y);
    let frame2 = Frame2 {
        origin: Vec2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    };

    let mut builder = ExactBRepBuilder::default();
    // 4 vertices, 4 circular edges, 4 loops, 4 faces.
    builder
        .topology_mut()
        .try_reserve(4, 6, 6, 4, 1, 1)
        .map_err(|_| GeomError::BudgetExceeded {
            resource: "exact revolution topology",
        })?;

    // One circular edge per (height, radius) corner of the revolved section.
    let circle_edge = |builder: &mut ExactBRepBuilder, y: Scalar, radius: Scalar| {
        let position = Point3::new(axis_origin.x + radius, y, 0.0);
        let vertex = builder.topology_mut().add_vertex(Vertex { position });
        let curve = builder.add_curve3(Curve3::Circle(Circle3 {
            frame: frame_at(y),
            radius,
        }));
        let edge = builder.topology_mut().add_edge(Edge {
            start: vertex,
            end: vertex,
            curve: Some(curve),
        });
        builder.set_edge_interval(edge, Interval::new(0.0, TAU));
        (edge, vertex, position)
    };

    let (bottom_outer, bottom_outer_vertex, bottom_outer_point) =
        circle_edge(&mut builder, bottom, outer);
    let (bottom_inner, bottom_inner_vertex, bottom_inner_point) =
        circle_edge(&mut builder, bottom, inner);
    let (top_outer, top_outer_vertex, _) = circle_edge(&mut builder, top, outer);
    let (top_inner, top_inner_vertex, _) = circle_edge(&mut builder, top, inner);

    // Seam edges close the cylinder walls; see `wall_loop`.
    let seam_edge = |builder: &mut ExactBRepBuilder, start, end, from: Point3| {
        let curve = builder.add_curve3(Curve3::Line(Line3 {
            origin: from,
            direction: Vec3::new(0.0, top - bottom, 0.0),
        }));
        let edge = builder.topology_mut().add_edge(Edge {
            start,
            end,
            curve: Some(curve),
        });
        builder.set_edge_interval(edge, Interval::UNIT);
        edge
    };
    let outer_seam = seam_edge(
        &mut builder,
        bottom_outer_vertex,
        top_outer_vertex,
        bottom_outer_point,
    );
    let inner_seam = seam_edge(
        &mut builder,
        bottom_inner_vertex,
        top_inner_vertex,
        bottom_inner_point,
    );

    // Cap loops: an annulus has an outer boundary and an inner hole, so each
    // cap face carries two loops rather than one.
    let circle_loop = |builder: &mut ExactBRepBuilder, edge, radius, orientation| {
        let pcurve = builder.add_curve2(Curve2::Circle(Circle2 {
            frame: frame2,
            radius,
        }));
        let loop_id = builder.topology_mut().add_loop(Loop {
            edges: vec![EdgeUse {
                edge,
                orientation,
                pcurve: Some(pcurve),
            }],
        });
        // A pcurve interval runs the way its use traverses the edge
        // (ADR 0024): a reversed use walks the circle backwards.
        let interval = match orientation {
            Orientation::Forward => Interval::new(0.0, TAU),
            Orientation::Reversed => Interval::new(TAU, 0.0),
        };
        builder.set_pcurve_interval(loop_id, 0, interval);
        loop_id
    };

    let bottom_outer_loop = circle_loop(&mut builder, bottom_outer, outer, Orientation::Forward);
    let bottom_inner_loop = circle_loop(&mut builder, bottom_inner, inner, Orientation::Reversed);
    let top_outer_loop = circle_loop(&mut builder, top_outer, outer, Orientation::Forward);
    let top_inner_loop = circle_loop(&mut builder, top_inner, inner, Orientation::Reversed);

    // A cylinder wall between two coaxial circles needs a SEAM edge. A loop
    // must be vertex-connected -- consecutive edges share an endpoint -- and
    // the two circles are disjoint, so a two-edge loop is an open loop and the
    // topology audit rejects it. The seam runs up the wall and is traversed
    // once forward and once reversed, closing the loop without adding area.
    let wall_loop = |builder: &mut ExactBRepBuilder, lower, upper, seam, height: Scalar| {
        let lower_pcurve = builder.add_curve2(Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::X,
        }));
        let upper_pcurve = builder.add_curve2(Curve2::Line(Line2 {
            origin: Vec2::new(0.0, height),
            direction: Vec2::X,
        }));
        let up_pcurve = builder.add_curve2(Curve2::Line(Line2 {
            origin: Vec2::new(TAU, 0.0),
            direction: Vec2::new(0.0, height),
        }));
        let down_pcurve = builder.add_curve2(Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(0.0, height),
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
    };

    let height = top - bottom;
    let outer_loop = wall_loop(&mut builder, bottom_outer, top_outer, outer_seam, height);
    let inner_loop = wall_loop(&mut builder, bottom_inner, top_inner, inner_seam, height);

    let bottom_surface = builder.add_surface(Surface::Plane(Plane {
        frame: frame_at(bottom),
    }));
    let top_surface = builder.add_surface(Surface::Plane(Plane {
        frame: frame_at(top),
    }));
    let outer_surface = builder.add_surface(Surface::Cylinder(Cylinder {
        frame: frame_at(bottom),
        radius: outer,
    }));
    let inner_surface = builder.add_surface(Surface::Cylinder(Cylinder {
        frame: frame_at(bottom),
        radius: inner,
    }));

    // Each cap is an annulus: outer boundary plus an inner hole.
    let annulus = |builder: &mut ExactBRepBuilder, surface, outer_loop, inner_loop, orientation| {
        builder.topology_mut().add_face(Face {
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
        })
    };

    let bottom_face = annulus(
        &mut builder,
        bottom_surface,
        bottom_outer_loop,
        bottom_inner_loop,
        Orientation::Reversed,
    );
    let top_face = annulus(
        &mut builder,
        top_surface,
        top_outer_loop,
        top_inner_loop,
        Orientation::Forward,
    );

    let wall_face = |builder: &mut ExactBRepBuilder, surface, loop_id, orientation| {
        builder.topology_mut().add_face(Face {
            surface: Some(surface),
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation,
        })
    };
    // The inner wall faces inward, so its orientation is reversed relative to
    // the outer wall: both must point out of the material.
    let outer_face = wall_face(
        &mut builder,
        outer_surface,
        outer_loop,
        Orientation::Forward,
    );
    let inner_face = wall_face(
        &mut builder,
        inner_surface,
        inner_loop,
        Orientation::Reversed,
    );

    finish_solid(builder, vec![bottom_face, top_face, outer_face, inner_face])
}

/// Assemble the faces into a closed solid and audit it.
///
/// The manifold audit is the real gate: it catches a wall wound the wrong way
/// or a cap loop that does not close, which are exactly the mistakes that
/// produce a B-rep that looks plausible and bounds nothing.
fn finish_solid(mut builder: ExactBRepBuilder, faces: Vec<FaceId>) -> GeomResult<ExactBRep> {
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
            detail: format!("exact revolution assembly failed: {error}"),
        })?;
    let health = axiolid_topology::audit_brep(exact.topology());
    if !health.is_closed_manifold() {
        return Err(GeomError::BackendContractViolation {
            backend: BACKEND_ID,
            detail: format!("exact revolution is not a closed manifold: {health:?}"),
        });
    }
    Ok(exact)
}

/// Sweep a supported profile along a straight path with a fixed reference.
///
/// A fixed-reference sweep along a STRAIGHT directrix keeps the profile's
/// orientation constant, so the swept solid is exactly a linear extrusion
/// along that segment. Rather than build a second implementation that could
/// drift from it, this delegates to `extrude_profile_exact`, which is the
/// same solid by construction.
///
/// A curved directrix is refused: the profile then rotates along the path and
/// the walls become general swept surfaces, not planes and cylinders. That is
/// the curved-surface work this issue explicitly puts out of scope.
pub fn fixed_reference_sweep_exact(
    profile: &Profile,
    path: &[Point3],
    reference_direction: Vec3,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    if path.len() < 2 {
        return Err(GeomError::InvalidInput(
            "a sweep path needs at least two points".to_owned(),
        ));
    }
    if !path.iter().all(|p| p.is_finite()) {
        return Err(GeomError::InvalidInput(
            "sweep path points must be finite".to_owned(),
        ));
    }
    if !reference_direction.is_finite() || reference_direction.length() <= 0.0 {
        return Err(GeomError::InvalidInput(
            "sweep reference direction must be finite and non-zero".to_owned(),
        ));
    }

    // Straightness is the precondition, so it is checked rather than assumed:
    // every interior point must lie on the segment from first to last.
    let start = path[0];
    let end = path[path.len() - 1];
    let span = end - start;
    let length = span.length();
    if length <= tolerance.linear() {
        return Err(GeomError::Degenerate(
            "sweep path start and end coincide, so the path has no direction".to_owned(),
        ));
    }
    let direction = span / length;
    for point in &path[1..path.len() - 1] {
        let offset = *point - start;
        let perpendicular = offset - direction * direction.dot(offset);
        if perpendicular.length() > tolerance.linear() {
            return Err(unsupported("exact sweep along a curved directrix"));
        }
    }

    extrude_profile_exact(profile, direction, length, tolerance)
}

/// Revolve a profile through a partial turn: every family the contour path
/// revolves a full turn, routed through [`crate::revolve_partial`].
fn revolve_partial(
    profile: &Profile,
    frame: &AxisFrame,
    angle: Scalar,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    match profile {
        Profile::Ellipse(_) => return Err(unsupported("ellipse exact revolution")),
        Profile::Rectangle(_)
        | Profile::Circle(_)
        | Profile::Section(_)
        | Profile::Contour(_)
        | Profile::Derived { .. }
        | Profile::Composite(_)
        | Profile::CenterLine(_) => {}
        _ => return Err(unsupported("unknown profile exact revolution")),
    }
    revolve_via_contour(profile, frame, Some(angle), tolerance)
}

/// Revolve any profile that lowers to a contour, a full turn (`angle` is
/// `None`) or partially.
///
/// About an axis exactly along y the sections are built in place, a turn
/// about -y running the other way about +y. About any other in-plane axis
/// each section is mapped into the axis's frame ([`into_axis_frame`]),
/// revolved about +y there with the angle as given, and the solid placed back
/// (#450).
fn revolve_via_contour(
    profile: &Profile,
    frame: &AxisFrame,
    angle: Option<Scalar>,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let placed = frame.placement.is_some();
    let (axis_origin, angle) = if placed {
        (Point3::ZERO, angle)
    } else {
        // The sections are built about world +y; an axis pointing down turns
        // the right-hand sweep the other way.
        (
            Point3::new(frame.origin.x, frame.origin.y, 0.0),
            angle.map(|angle| if frame.given.y < 0.0 { -angle } else { angle }),
        )
    };
    let section = |ring: &ArcRing| -> GeomResult<ArcRing> {
        if placed {
            into_axis_frame(ring, frame)
        } else {
            Ok(ring.clone())
        }
    };
    // The mapped coordinates are differences of the given ones and the
    // axis's: the rounding gates of the revolution read their magnitude.
    let magnitude = |rings: &[&ArcRing]| -> Scalar {
        if !placed {
            return 0.0;
        }
        rings
            .iter()
            .flat_map(|ring| ring.vertices.iter())
            .fold(frame.origin.x.abs().max(frame.origin.y.abs()), |m, v| {
                m.max(v.point.x.abs()).max(v.point.y.abs())
            })
    };

    // A composite is unioned exactly first; members that do not touch
    // revolve into separate solids of one `ExactBRep` (#111).
    if let Profile::Composite(members) = profile {
        let regions = crate::profile_lower::composite_regions(members, tolerance)?;
        let mut pieces = Vec::with_capacity(regions.len());
        for region in &regions {
            let holes = region
                .holes
                .iter()
                .map(&section)
                .collect::<GeomResult<Vec<_>>>()?;
            let rings: Vec<&ArcRing> = core::iter::once(&region.outer)
                .chain(region.holes.iter())
                .collect();
            pieces.push(revolve_region(
                &section(&region.outer)?,
                &holes,
                axis_origin,
                magnitude(&rings),
                angle,
                tolerance,
            )?);
        }
        return place_in_axis_frame(crate::assemble::merge_solids(pieces)?, frame);
    }

    let contour = crate::extrude_exact::profile_to_contour(profile, tolerance)?;
    let outer = crate::contour_lower::contour_to_arc_ring(&contour.outer, tolerance)?;
    let holes = contour
        .holes
        .iter()
        .map(|hole| crate::contour_lower::contour_to_arc_ring(hole, tolerance))
        .collect::<GeomResult<Vec<_>>>()?;
    let rings: Vec<&ArcRing> = core::iter::once(&outer).chain(holes.iter()).collect();
    let scale = magnitude(&rings);
    let holes = holes.iter().map(&section).collect::<GeomResult<Vec<_>>>()?;
    let solid = revolve_region(
        &section(&outer)?,
        &holes,
        axis_origin,
        scale,
        angle,
        tolerance,
    )?;
    place_in_axis_frame(solid, frame)
}

/// Revolve one section a full turn -- its outer ring as the solid, each hole
/// as a void -- or partially, each hole a tunnel between the caps.
fn revolve_region(
    outer: &ArcRing,
    holes: &[ArcRing],
    axis_origin: Point3,
    magnitude: Scalar,
    angle: Option<Scalar>,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    if let Some(angle) = angle {
        return crate::revolve_partial::revolve_section_partial(
            outer,
            holes,
            axis_origin,
            angle,
            tolerance,
        );
    }
    let solid =
        crate::revolve_contour::revolve_ring_scaled(outer, axis_origin, magnitude, tolerance)?;
    if holes.is_empty() {
        return Ok(solid);
    }
    if solid.topology().solids().len() != 1 || !solid.topology().solids()[0].voids.is_empty() {
        return Err(unsupported(
            "exact revolution of a section with holes and several contacts with the axis",
        ));
    }
    // A hole swept a full turn encloses a ring-shaped cavity: it touches
    // neither cap, so it is a second shell, not a second loop. Each hole
    // revolves on its own and joins the solid as a void, its faces used
    // reversed so they face into the cavity (#111).
    let mut builder = ExactBRepBuilder::default();
    let outer = builder.append(&solid, false);
    let mut voids = Vec::with_capacity(holes.len());
    for hole in holes {
        let cavity =
            crate::revolve_contour::revolve_ring_scaled(hole, axis_origin, magnitude, tolerance)?;
        voids.extend(builder.append(&cavity, true));
    }
    let [outer] = outer[..] else {
        return Err(GeomError::Degenerate(
            "a revolved section must have one outer shell".to_owned(),
        ));
    };
    builder.topology_mut().add_solid(Solid { outer, voids });
    builder.finish().map_err(|error| {
        GeomError::Degenerate(format!(
            "revolved section with holes did not assemble: {error}"
        ))
    })
}
