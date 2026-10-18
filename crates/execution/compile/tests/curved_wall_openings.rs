//! A curved wall with a row of abutting openings (#332), certified, and a
//! box tangent to a curved wall's inner face (#327), exact.
//!
//! The wall is an annular sector, radii `6` to `6.25` about the origin,
//! `0` to `60` degrees, `3` high, less five wedge openings with radial
//! sides (radius `5` to `7`, `z` from `0.5` to `2.5`), cut one after
//! another as `ifc-geometry` lowers a wall's openings. Three abut, each
//! side reaching `15 um` past its neighbour's; one reaches `0.6 mm` past
//! another's; one stands on its own. The side planes of abutting openings
//! are turned by more than the angular tolerance and lie within the linear
//! one in the wall, and the general boolean refused the cut, so the
//! deviation report left the wall unbounded. The box of #327 touches a
//! quarter ring's inner cylinder along a ruling, its back face on the
//! tangent plane or a unit in the last place either side of it.
//!
//! The wall is now certified within the millimetre requested, and the
//! exact results have the closed-form volume. The box's exact result does
//! too; its mesh, touching itself along the contact, stays refused.

use std::f64::consts::PI;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Frame2, Interval, Point2, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{BooleanReport, ReferenceExactCompiler, ReferenceMeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

fn line(a: Point2, b: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn at(radius: f64, angle: f64) -> Point2 {
    Point2::new(radius * angle.cos(), radius * angle.sin())
}

/// An annular sector between radii `inner` and `outer` from angle `from`
/// to `to`, anticlockwise: outer arc, line in, inner arc back, line out,
/// as `ifc-geometry` lowers an `IfcIndexedPolyCurve` of two arcs and two
/// lines.
fn sector(outer: f64, inner: f64, from: f64, to: f64) -> Profile {
    let circle = |radius| {
        Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        })
    };
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve: circle(outer),
                domain: Interval::new(from, to),
                same_sense: true,
            },
            line(at(outer, to), at(inner, to)),
            ProfileSegment {
                curve: circle(inner),
                domain: Interval::new(from, to),
                same_sense: false,
            },
            line(at(inner, from), at(outer, from)),
        ]),
        holes: Vec::new(),
    })
}

/// A quadrilateral through `corners`, anticlockwise.
fn quad(corners: [Point2; 4]) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(
            (0..4)
                .map(|i| line(corners[i], corners[(i + 1) % 4]))
                .collect(),
        ),
        holes: Vec::new(),
    })
}

/// `profile` extruded up by `depth` and placed by `transform`.
fn placed(
    builder: &mut GeometryGraphBuilder,
    profile: Profile,
    depth: f64,
    transform: Transform3,
) -> NodeId {
    let profile = builder.push(GeometryNode::Profile(profile)).unwrap();
    let solid = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
        .unwrap();
    builder
        .push(GeometryNode::Instance(Instance {
            source: solid,
            transform,
        }))
        .unwrap()
}

fn minus(builder: &mut GeometryGraphBuilder, left: NodeId, right: NodeId) -> NodeId {
    builder
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator: BooleanOperator::Difference,
        }))
        .unwrap()
}

fn up(z: f64) -> Transform3 {
    Transform3::from_translation(Vec3::new(0.0, 0.0, z))
}

/// The volume of an annular sector extruded by `height`.
fn sector_volume(outer: f64, inner: f64, sweep: f64, height: f64) -> f64 {
    0.5 * sweep * (outer * outer - inner * inner) * height
}

/// How far an abutting opening's side reaches past its neighbour's, as an
/// angle: `15 um` and `0.6 mm` at the wall.
const DRIFT: f64 = 2.5e-6;
const NEAR: f64 = 1e-4;

fn deg(d: f64) -> f64 {
    d.to_radians()
}

/// The wall less its openings, and the closed-form volume.
fn wall_with_openings() -> (GeometryGraph, NodeId, f64) {
    let openings = [
        (deg(20.0), deg(26.0)),
        (deg(14.0), deg(20.0) + DRIFT),
        (deg(8.0), deg(14.0) + DRIFT),
        (deg(26.0) - NEAR, deg(32.0)),
        (deg(40.0), deg(48.0)),
    ];
    let mut builder = GeometryGraphBuilder::new();
    let mut root = placed(&mut builder, sector(6.25, 6.0, 0.0, PI / 3.0), 3.0, up(0.0));
    for (from, to) in openings {
        let wedge = quad([at(5.0, from), at(7.0, from), at(7.0, to), at(5.0, to)]);
        let opening = placed(&mut builder, wedge, 2.0, up(0.5));
        root = minus(&mut builder, root, opening);
    }
    // The openings cover `8` to `32` degrees together, and `40` to `48`.
    let removed = deg(24.0) + deg(8.0);
    let volume = sector_volume(6.25, 6.0, PI / 3.0, 3.0) - sector_volume(6.25, 6.0, removed, 2.0);
    (builder.finish(vec![root]).unwrap(), root, volume)
}

/// The quarter ring of #327 less a box in the tangent frame at `45`
/// degrees, from `back` to `back + 0.1` along the radius, `0.3` either
/// side, `z` from `1` to `2`; and the closed-form volume.
fn ring_less_tangent_box(back: f64) -> (GeometryGraph, NodeId, f64) {
    let mut builder = GeometryGraphBuilder::new();
    let ring = placed(&mut builder, sector(5.0, 4.8, 0.0, PI / 2.0), 3.0, up(0.0));
    let profile = quad([
        Point2::new(back, -0.3),
        Point2::new(back + 0.1, -0.3),
        Point2::new(back + 0.1, 0.3),
        Point2::new(back, 0.3),
    ]);
    let tool = placed(
        &mut builder,
        profile,
        1.0,
        Transform3::from_rotation_z(deg(45.0)) * up(1.0),
    );
    let root = minus(&mut builder, ring, tool);
    // The whole box lies in the wall but for the ruling its back face
    // touches the inner face along.
    let volume = sector_volume(5.0, 4.8, PI / 2.0, 3.0) - 0.1 * 0.6 * 1.0;
    (builder.finish(vec![root]).unwrap(), root, volume)
}

/// The exact result's report after checking its volume against `want`.
fn exact(what: &str, graph: &GeometryGraph, root: NodeId, want: f64) -> BooleanReport {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let (exact, report) = ReferenceExactCompiler::new()
        .compile_exact_with_report(graph, root, &options)
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    let volume = exact_properties(&exact, options.tolerance())
        .expect("measurable")
        .signed_volume;
    assert!(
        (volume - want).abs() <= 1e-9 * want,
        "{what}: expected {want}, got {volume}"
    );
    report
}

#[test]
fn a_curved_wall_with_abutting_openings_is_certified() {
    let (graph, root, want) = wall_with_openings();
    let report = exact("openings", &graph, root, want);
    // Cut at the finer tolerance where the openings abut.
    assert!(report.finer(), "{report:?}");
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let (_, deviation) = ReferenceMeshCompiler::new(BoolmeshBoolean)
        .compile_mesh_with_deviation(&graph, root, &options)
        .unwrap_or_else(|e| panic!("{e}"));
    let bound = deviation.bound.unwrap_or_else(|| panic!("{deviation:?}"));
    assert!(bound <= 1e-3, "{deviation:?}");
}

#[test]
fn a_box_tangent_to_a_curved_walls_inner_face_is_exact_and_its_mesh_refused() {
    // The exact boolean decides it at the caller's tolerance, at the
    // closed-form volume, wherever the back face lies. On the tangent
    // plane, or a unit in the last place past it (cutting the inner face
    // in two rulings `6e-8` apart), the result's material is as thin as
    // nothing along the contact, and the mesh compiler refuses a result
    // touching itself by name (#194): a known refusal of the mesh path,
    // pinned here. A unit in the last place short of it, the box touches
    // nothing but material, and its mesh is certified.
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    for (back, touching) in [
        (4.8, true),
        (4.8_f64.next_down(), true),
        (4.8_f64.next_up(), false),
    ] {
        let what = format!("back at {back}");
        let (graph, root, want) = ring_less_tangent_box(back);
        let report = exact(&what, &graph, root, want);
        assert!(!report.finer(), "{what}: {report:?}");
        let outcome = ReferenceMeshCompiler::new(BoolmeshBoolean)
            .compile_mesh_with_deviation(&graph, root, &options);
        if touching {
            let error = outcome.map(|(_, d)| d).expect_err(&what);
            assert!(
                error.to_string().contains("touches itself"),
                "{what}: {error}"
            );
        } else {
            let (_, deviation) = outcome.unwrap_or_else(|e| panic!("{what}: {e}"));
            let bound = deviation
                .bound
                .unwrap_or_else(|| panic!("{what}: {deviation:?}"));
            assert!(bound <= 1e-3, "{what}: {deviation:?}");
        }
    }
}
