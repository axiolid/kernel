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
//! # What is refused, and why that is the honest answer
//!
//! A PARTIAL turn is not this shape. It has two extra planar walls at the
//! start and end angles, and its cap loops are not closed circles but
//! circular arcs joined by radial segments. That is a different topology, not
//! a parameter change, so it is refused rather than approximated.
//!
//! A profile crossing the axis degenerates: the inner cylinder collapses to
//! the axis line and the caps stop being annuli. Refused for the same reason.
//!
//! The mesh path in `revolve.rs` handles all of these. Refusing here means a
//! caller asking for exactness gets a typed refusal naming the gap instead of
//! a silently tessellated substitute.

use std::f64::consts::TAU;

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Interval, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{
    Edge, EdgeUse, Face, FaceBound, FaceId, Loop, Orientation, Shell, Solid, Vertex,
};

use crate::extrude_exact::extrude_profile_exact;
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
/// Supported: a sharp filled rectangle, revolved a full turn about an axis
/// parallel to the profile's local y and offset from it, so the profile does
/// not cross the axis. The result is an annular tube of two cylinders and two
/// annular planar caps.
///
/// Every other case -- partial turns, profiles crossing the axis, rounded or
/// hollow rectangles, other profile families -- returns a typed refusal.
pub fn revolve_profile_exact(
    profile: &Profile,
    axis_origin: Point3,
    axis_direction: Vec3,
    angle: Scalar,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    if !angle.is_finite() {
        return Err(GeomError::InvalidInput(
            "revolution angle must be finite".to_owned(),
        ));
    }
    // A partial turn has two extra planar walls and arc-bounded caps: a
    // different topology, not a different parameter.
    if (angle.abs() - TAU).abs() > tolerance.linear() {
        return Err(unsupported("partial-turn exact revolution"));
    }

    match profile {
        Profile::Rectangle(rectangle) => {
            revolve_rectangle(rectangle, axis_origin, axis_direction, tolerance)
        }
        Profile::Circle(_) => Err(unsupported("circle-profile exact revolution")),
        Profile::Ellipse(_) => Err(unsupported("ellipse exact revolution")),
        // Every one of these lowers to a contour, and a contour revolves.
        // The refusals they carried described a missing module, not missing
        // geometry, so they route through the general path rather than
        // restating a limitation that no longer holds.
        Profile::Section(_)
        | Profile::Contour(_)
        | Profile::Derived { .. }
        | Profile::Composite(_)
        | Profile::CenterLine(_) => {
            revolve_via_contour(profile, axis_origin, axis_direction, tolerance)
        }
        _ => Err(unsupported("unknown profile exact revolution")),
    }
}

fn revolve_rectangle(
    rectangle: &RectangleProfile,
    axis_origin: Point3,
    axis_direction: Vec3,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    if rectangle.thickness.is_some() {
        return Err(unsupported("hollow rectangle exact revolution"));
    }
    if rectangle.outer_radius.is_some() || rectangle.inner_radius.is_some() {
        // Rounded corners revolve into tori, which the general contour
        // revolver already builds; the dedicated path below only knows
        // straight edges.
        return revolve_via_contour(
            &Profile::Rectangle(*rectangle),
            axis_origin,
            axis_direction,
            tolerance,
        );
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

    let length = axis_direction.length();
    if !length.is_finite() || length <= 0.0 {
        return Err(GeomError::InvalidInput(
            "revolution axis must be a finite non-zero direction".to_owned(),
        ));
    }
    let axis = axis_direction / length;

    // The profile lives in the z = 0 plane. Only an axis parallel to local y
    // sweeps its edges into cylinders and its horizontal edges into annuli;
    // any other axis produces cones or general surfaces of revolution, which
    // is a different construction.
    if (axis.dot(Vec3::Y).abs() - 1.0).abs() > tolerance.linear() {
        return Err(unsupported(
            "exact revolution about an axis that is not the profile's local y",
        ));
    }
    if axis_origin.z.abs() > tolerance.linear() {
        return Err(unsupported(
            "exact revolution about an axis off the profile plane",
        ));
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
    if left * right < 0.0 || left.abs() <= tolerance.linear() || right.abs() <= tolerance.linear() {
        return Err(unsupported(
            "exact revolution of a profile touching or crossing the axis",
        ));
    }
    let inner = left.abs().min(right.abs());
    let outer = left.abs().max(right.abs());

    // The rectangle spans y in [-half_y, half_y] around the axis origin.
    let bottom = axis_origin.y - half_y;
    let top = axis_origin.y + half_y;
    let frame_at = |y: Scalar| Frame3 {
        origin: Point3::new(axis_origin.x, y, 0.0),
        x: Vec3::X,
        y: Vec3::Z,
        z: Vec3::Y,
    };
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
        builder.set_pcurve_interval(loop_id, 0, Interval::new(0.0, TAU));
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

/// Revolve any profile that lowers to a contour.
///
/// The axis must be the profile's local y and must lie in the profile plane,
/// the same restriction `revolve_rectangle` carries: any other axis sweeps a
/// general surface of revolution rather than the named quadrics this builds.
fn revolve_via_contour(
    profile: &Profile,
    axis_origin: Point3,
    axis_direction: Vec3,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let axis = axis_direction.normalize_or_zero();
    if (axis.dot(Vec3::Y).abs() - 1.0).abs() > tolerance.linear() {
        return Err(unsupported(
            "exact revolution about an axis that is not the profile's local y",
        ));
    }
    if axis_origin.z.abs() > tolerance.linear() {
        return Err(unsupported(
            "exact revolution about an axis off the profile plane",
        ));
    }

    let contour = crate::extrude_exact::profile_to_contour(profile, tolerance)?;
    if !contour.holes.is_empty() {
        // A hole in a revolved section makes an internal void, which is a
        // second shell rather than a second loop on a cap.
        return Err(unsupported("exact revolution of a section with holes"));
    }
    let ring = crate::contour_lower::contour_to_arc_ring(&contour.outer, tolerance)?;
    crate::revolve_contour::revolve_arc_ring(&ring, axis_origin, tolerance)
}
