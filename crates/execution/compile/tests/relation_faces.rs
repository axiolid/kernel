//! B-rep faces on swept and ruled surface relations (#411, ADR 0103).
//!
//! A face over a `Revolution`, `LinearExtrusion` or `Ruled` relation whose
//! edge uses carry pcurves is meshed with a bound certified from the
//! relation's closed-form derivative bounds (`relation_edge_inversion.rs`
//! covers uses without one). Closed solids are a wall on the relation with
//! two planar caps, each checked against the same solid on the atomic
//! surface it coincides with, for the fixture's sake. Open sheets are
//! checked against a closed-form point of the surface. Every case reads the
//! mesh back: closure, watertightness, volume and a dense sample of the
//! exact surface against the reported bound, which is certified within the
//! budget or the case fails.

mod common;

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, PI, TAU};

use axiolid_contracts::{GeomError, Operation};
use axiolid_core::{BooleanOperator, Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{Curve3, Polyline3};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, ReferenceMeshCompiler};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, SolidOperation};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_surface::{Cone, Cylinder, Sphere, Surface, Torus};

use common::*;

#[test]
fn atomic_controls_compile_certified() {
    let cylinder = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Cylinder(Cylinder {
                frame: frame_at(0.0),
                radius: RADIUS,
            }),
        )
    };
    let cone = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Cone(Cone {
                frame: frame_at(0.0),
                radius: 0.2,
                semi_angle: 0.2_f64.atan(),
            }),
        )
    };
    let sphere = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Sphere(Sphere {
                frame: frame_at(0.0),
                radius: 0.3,
            }),
        )
    };
    let family = |name| Solid3 {
        family: name,
        ..cylinder_solid()
    };
    assert_closed_wall("cylinder", &cylinder_wall(&cylinder), &family("cylinder"));
    assert_closed_wall(
        "cone",
        &frustum_wall(&cone),
        &Solid3 {
            family: "cone",
            ..frustum_solid()
        },
    );
    assert_closed_wall(
        "sphere",
        &zone_wall(&sphere),
        &Solid3 {
            family: "sphere",
            ..zone_solid()
        },
    );
}

// --- revolution: closed solids --------------------------------------------------

#[test]
fn a_line_revolved_about_a_parallel_axis_meshes_as_its_cylinder() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.3, 0.0, 0.0), Vec3::Z),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    assert_closed_wall("revolved line", &cylinder_wall(&support), &cylinder_solid());
}

#[test]
fn a_slanted_line_revolved_meshes_as_its_frustum() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.2, 0.0, 0.0), Vec3::new(0.2, 0.0, 1.0)),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    assert_closed_wall(
        "revolved slanted line",
        &frustum_wall(&support),
        &frustum_solid(),
    );
}

#[test]
fn an_arc_revolved_meshes_as_a_spherical_zone() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            circle_in(Point3::ZERO, Vec3::X, Vec3::Z, 0.3),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    assert_closed_wall("revolved arc", &zone_wall(&support), &zone_solid());
}

// --- revolution: open sheets ----------------------------------------------------

#[test]
fn a_revolved_arc_patch_like_revit_s_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(b, torus_profile(), Point3::new(1.0, 2.0, 0.0), Vec3::Z)
    };
    let sheet = torus_patch(&torus_point);
    assert_sheet(
        "revolved arc patch",
        &support,
        &sheet,
        "surface of revolution",
    );
    let control = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Torus(Torus {
                frame: Frame3 {
                    origin: Point3::new(1.0, 2.0, 0.1),
                    x: Vec3::X,
                    y: Vec3::Y,
                    z: Vec3::Z,
                },
                major_radius: 0.4,
                minor_radius: 0.1,
            }),
        )
    };
    assert_sheet("torus control", &control, &sheet, "torus");
}

#[test]
fn a_bezier_profile_revolved_like_a_revit_hermite_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(b, bezier(bezier_points()), Point3::ZERO, Vec3::Z)
    };
    let turned = bezier_points().map(|p| Point3::new(-p.x, -p.y, p.z));
    let sheet = Sheet {
        u: (0.0, PI),
        v: (0.0, 1.0),
        point: &bezier_point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, 0.3)),
            Some(bezier(turned)),
            Some(horizontal_circle(0.0, 0.0, 0.3, 0.3)),
            Some(bezier(bezier_points())),
        ],
        pcurves: [true; 4],
        window: None,
        hole: None,
    };
    assert_sheet("revolved Bezier", &support, &sheet, "surface of revolution");
}

/// A hole on a curved face is bounded above the budget today, on an atomic
/// surface as much as on a relation: the relation's bound is the atomic
/// one's, which this checks instead of the budget.
#[test]
fn a_revolved_sheet_with_a_hole_meshes_round_it() {
    let hole = || Hole {
        centre: (FRAC_PI_2, 0.25),
        radius: 0.1,
    };
    let (mesh, _, report) = check_sheet(
        "revolved sheet with a hole",
        &cylinder_revolution,
        &cylinder_sheet(true, Some(hole())),
        "surface of revolution",
        false,
    );
    let control = |b: &mut GeometryGraphBuilder| {
        atomic(
            b,
            Surface::Cylinder(Cylinder {
                frame: frame_at(0.0),
                radius: RADIUS,
            }),
        )
    };
    let (_, _, atomic_report) = check_sheet(
        "cylinder sheet with a hole",
        &control,
        &cylinder_sheet(true, Some(hole())),
        "cylinder",
        false,
    );
    let (bound, atomic_bound) = (report.bound.unwrap(), atomic_report.bound.unwrap());
    assert!(
        (bound - atomic_bound).abs() <= 1e-9 * atomic_bound,
        "{bound:e} against the atomic cylinder's {atomic_bound:e}"
    );
    // The hole is 0.1 radians wide, which is 0.03 along the wall: inside it,
    // shrunk by the budget, no triangle has its centre.
    for t in mesh.indices.as_chunks::<3>().0 {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        let m = (a + b + c) / 3.0;
        let (u, z) = (m.y.atan2(m.x), m.z);
        let inside = ((u - FRAC_PI_2) * RADIUS / (0.1 * RADIUS - BUDGET)).powi(2)
            + ((z - 0.25) / (0.1 - BUDGET)).powi(2);
        assert!(inside > 1.0, "a triangle centred in the hole at ({u}, {z})");
    }
}

// --- refusals ---------------------------------------------------------------------

#[test]
fn a_revolution_over_a_polyline_is_refused_by_name() {
    let support = |b: &mut GeometryGraphBuilder| {
        let polyline = Curve3::Polyline(Polyline3 {
            points: vec![Point3::new(0.3, 0.0, 0.0), Point3::new(0.3, 0.0, 1.0)],
            closed: false,
        });
        revolution(b, polyline, Point3::ZERO, Vec3::Z)
    };
    let error = refused(open_sheet(&support, &cylinder_sheet(true, None)));
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::SurfaceEvaluation,
                input: "a swept curve other than a line, circle, ellipse or B-spline",
                ..
            }
        ),
        "{error:?}"
    );
}

#[test]
fn a_revolution_about_a_zero_axis_is_refused() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.3, 0.0, 0.0), Vec3::Z),
            Point3::ZERO,
            Vec3::ZERO,
        )
    };
    let error = refused(open_sheet(&support, &cylinder_sheet(true, None)));
    assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
}

// --- linear extrusion ---------------------------------------------------------------

#[test]
fn a_circle_extruded_along_its_normal_meshes_as_its_cylinder() {
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, Vec3::Z);
    let expect = Solid3 {
        family: EXTRUSION,
        ..cylinder_solid()
    };
    assert_closed_wall("extruded circle", &cylinder_wall(&support), &expect);
}

#[test]
fn an_obliquely_extruded_circle_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, Vec3::new(0.2, 0.0, 1.0));
    let sheet = Sheet {
        u: (0.0, PI),
        v: (0.0, HEIGHT),
        point: &oblique_point,
        curves: [
            Some(horizontal_circle(0.0, 0.0, 0.0, RADIUS)),
            None,
            Some(horizontal_circle(0.1, 0.0, HEIGHT, RADIUS)),
            None,
        ],
        pcurves: [true; 4],
        window: None,
        hole: None,
    };
    let (mesh, _, _) = assert_sheet("oblique extrusion", &support, &sheet, EXTRUSION);
    for p in &mesh.positions {
        let off = ((p.x - 0.2 * p.z).hypot(p.y) - RADIUS).abs();
        assert!(off <= 1e-12, "a vertex {p:?} is {off:e} off the surface");
    }
}

#[test]
fn a_linear_extrusion_with_a_zero_direction_is_refused() {
    let support = |b: &mut GeometryGraphBuilder| extruded_circle(b, Vec3::ZERO);
    let error = refused(open_sheet(&support, &cylinder_sheet(true, None)));
    assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
}

// --- ruled surfaces -------------------------------------------------------------------

#[test]
fn two_circles_ruled_mesh_as_their_cylinder() {
    let support = |b: &mut GeometryGraphBuilder| {
        ruled(
            b,
            [
                (horizontal_circle(0.0, 0.0, 0.0, RADIUS), 0.0, TAU),
                (horizontal_circle(0.0, 0.0, HEIGHT, RADIUS), 0.0, TAU),
            ],
        )
    };
    let wall = Wall {
        u: (0.0, 1.0),
        v: (0.0, 1.0),
        ..cylinder_wall(&support)
    };
    let expect = Solid3 {
        family: RULED,
        ..cylinder_solid()
    };
    assert_closed_wall("ruled circles", &wall, &expect);
}

#[test]
fn an_arc_ruled_against_a_line_like_revit_s_meshes_within_its_bound() {
    let support = |b: &mut GeometryGraphBuilder| arc_line_support(b, false);
    let sheet = unit_sheet(
        &arc_line_point,
        [
            Some(horizontal_circle(0.0, 0.0, 0.0, 1.0)),
            None,
            Some(line3(Point3::new(1.0, 0.0, 1.0), Vec3::new(-1.0, 1.0, 0.0))),
            None,
        ],
    );
    assert_sheet("arc against a line", &support, &sheet, RULED);
}

fn arcs_point(u: Scalar, v: Scalar) -> Point3 {
    let (a, b) = (FRAC_PI_2 * u, FRAC_PI_3 * u);
    Point3::new(a.cos(), a.sin(), 0.0) * (1.0 - v)
        + Point3::new(0.5 * b.cos(), 0.5 * b.sin(), 1.0) * v
}

#[test]
fn two_arcs_of_different_sweeps_ruled_mesh_within_their_bound() {
    let support = |b: &mut GeometryGraphBuilder| {
        ruled(
            b,
            [
                (horizontal_circle(0.0, 0.0, 0.0, 1.0), 0.0, FRAC_PI_2),
                (horizontal_circle(0.0, 0.0, 1.0, 0.5), 0.0, FRAC_PI_3),
            ],
        )
    };
    let sheet = unit_sheet(
        &arcs_point,
        [
            Some(horizontal_circle(0.0, 0.0, 0.0, 1.0)),
            None,
            Some(horizontal_circle(0.0, 0.0, 1.0, 0.5)),
            None,
        ],
    );
    assert_sheet("arcs of different sweeps", &support, &sheet, RULED);
}

#[test]
fn a_reversed_side_is_read_backwards() {
    let support = |b: &mut GeometryGraphBuilder| arc_line_support(b, true);
    let mirrored = |u: Scalar, v: Scalar| arc_line_point(1.0 - u, v);
    // The arc now runs from the quarter turn back to the start.
    let sheet = unit_sheet(
        &mirrored,
        [
            Some(circle_in(Point3::ZERO, Vec3::Y, Vec3::X, 1.0)),
            None,
            Some(line3(Point3::new(0.0, 1.0, 1.0), Vec3::new(1.0, -1.0, 0.0))),
            None,
        ],
    );
    assert_sheet("reversed sides", &support, &sheet, RULED);
}

#[test]
fn a_boolean_over_a_relation_face_operand_meshes_with_a_named_unbounded_report() {
    let support = |b: &mut GeometryGraphBuilder| {
        revolution(
            b,
            line3(Point3::new(0.3, 0.0, 0.0), Vec3::Z),
            Point3::ZERO,
            Vec3::Z,
        )
    };
    let (mut builder, brep) = closed_wall_parts(&cylinder_wall(&support));
    let cylinder = builder.push(GeometryNode::BRep(brep)).unwrap();
    let profile = builder
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 0.2,
                y: 0.2,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let block = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 0.3,
        }))
        .unwrap();
    let root = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: cylinder,
            right: block,
            operator: BooleanOperator::Difference,
        }))
        .unwrap();
    let graph = builder.finish(vec![root]).unwrap();
    let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options())
        .expect("the boolean meshes");
    assert!(report.bound.is_none(), "{report:?}");
    assert!(
        report.contributions.iter().any(|c| matches!(
            c.bound,
            DeviationBound::Unbounded(reason)
                if reason.starts_with("exact boolean over a B-rep operand with a face that is neither a plane")
        )),
        "{report:?}"
    );
}
