//! A flat wall less openings whose faces coincide with the wall's, or with
//! each other's, within about `1e-8` (#334), certified.
//!
//! The wall is `5 x 0.25 x 3`, less a door on its bottom face whose caps
//! lie `e` outside the wall's faces, a window whose caps lie `e` inside
//! them, and two abutting windows whose shared side is `e` apart. The exact
//! boolean reads each such pair as one within tolerance and reports the
//! readings; its result's pcurve loops then close only within them, and
//! the deviation measurement refused the faces ("an exact face's loops do
//! not close in its parameters"), leaving the wall unbounded. It now
//! admits a gap up to the boolean's largest reading: the wall is
//! certified, its exact result has the closed-form volume up to what the
//! readings move, and at `e = 0` the bound stays at rounding (#314).

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Interval, Point2, Scalar, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, ReferenceExactCompiler, ReferenceMeshCompiler};
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

/// The rectangle from `(x0, y0)` to `(x1, y1)`, anticlockwise.
fn rect(x0: Scalar, y0: Scalar, x1: Scalar, y1: Scalar) -> Profile {
    let c = [
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ];
    Profile::Contour(ContourProfile {
        outer: Contour::new((0..4).map(|i| line(c[i], c[(i + 1) % 4])).collect()),
        holes: Vec::new(),
    })
}

/// `profile` extruded up by `depth` and placed by `transform`.
fn placed(
    builder: &mut GeometryGraphBuilder,
    profile: Profile,
    depth: Scalar,
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

/// A wall placement as a building model has one: turned about `z` and moved.
fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

/// An opening in the wall's elevation: `x` from `x0` to `x1`, `z` from
/// `z0` to `z1`, cut along `-y` from `back` by `depth`.
struct Opening {
    x: (Scalar, Scalar),
    z: (Scalar, Scalar),
    back: Scalar,
    depth: Scalar,
}

/// The wall less its openings under `world`, the closed-form volume with
/// the near-coincident faces read as one, and the area of every operand's
/// faces (how far a reading of `r` can change the volume is `r` times it).
fn wall(e: Scalar, world: Transform3) -> (GeometryGraph, NodeId, Scalar, Scalar) {
    let (length, thickness, height) = (5.0, 0.25, 3.0);
    let openings = [
        // A door on the bottom face, its caps `e` outside the wall's.
        Opening {
            x: (0.5, 1.5),
            z: (0.0, 2.2),
            back: thickness + e,
            depth: thickness + 2.0 * e,
        },
        // A window, its caps `e` inside them.
        Opening {
            x: (2.5, 3.5),
            z: (1.0, 2.0),
            back: thickness - e,
            depth: thickness - 2.0 * e,
        },
        // Two abutting windows through the wall, the shared side `e` apart.
        Opening {
            x: (3.8, 4.3),
            z: (0.5, 1.0),
            back: 0.3,
            depth: 0.35,
        },
        Opening {
            x: (4.3 - e, 4.8),
            z: (0.5, 1.0),
            back: 0.3,
            depth: 0.35,
        },
    ];
    let mut builder = GeometryGraphBuilder::new();
    let mut root = placed(
        &mut builder,
        rect(0.0, 0.0, length, thickness),
        height,
        world,
    );
    let mut area = 2.0 * (length * thickness + length * height + thickness * height);
    let across =
        |y: Scalar| Transform3::from_cols(Vec3::X, Vec3::Z, Vec3::NEG_Y, Vec3::new(0.0, y, 0.0));
    for o in &openings {
        let profile = rect(o.x.0, o.z.0, o.x.1, o.z.1);
        let opening = placed(&mut builder, profile, o.depth, world * across(o.back));
        root = minus(&mut builder, root, opening);
        let (w, h) = (o.x.1 - o.x.0, o.z.1 - o.z.0);
        area += 2.0 * (w * h + (w + h) * o.depth);
    }
    // Read as flush, the door and the first window go through; the
    // abutting windows remove `3.8..4.8` together.
    let removed = 1.0 * 2.2 + 1.0 * 1.0 + 1.0 * 0.5;
    let volume = (length * height - removed) * thickness;
    (builder.finish(vec![root]).unwrap(), root, volume, area)
}

#[test]
fn a_flat_wall_with_openings_coinciding_within_readings_is_certified() {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    for (name, world) in [("local", Transform3::IDENTITY), ("placed", building())] {
        for e in [0.0, 1e-8, 3e-8] {
            let what = format!("{name}, e = {e:e}");
            let (graph, root, want, area) = wall(e, world);
            let (exact, report) = ReferenceExactCompiler::new()
                .compile_exact_with_report(&graph, root, &options)
                .unwrap_or_else(|error| panic!("{what}: {error}"));
            let reading = report.linear();
            if e > 0.0 {
                // The readings the boolean took are of the offsets.
                assert!(reading >= e && reading <= 4.0 * e, "{what}: {report:?}");
            }
            // The closed-form volume, up to what the readings (and, placed,
            // the rounding) move the operands' faces by.
            let volume = exact_properties(&exact, options.tolerance())
                .expect("measurable")
                .signed_volume;
            let moved = (reading + report.rounding_floor()) * area + 1e-12 * want;
            assert!(
                (volume - want).abs() <= moved,
                "{what}: expected {want} within {moved:e}, got {volume}"
            );
            let (_, deviation) = ReferenceMeshCompiler::new(BoolmeshBoolean)
                .compile_mesh_with_deviation(&graph, root, &options)
                .unwrap_or_else(|error| panic!("{what}: {error}"));
            for c in &deviation.contributions {
                assert!(
                    matches!(c.bound, DeviationBound::Certified(_)),
                    "{what}: {deviation:?}"
                );
            }
            let bound = deviation
                .bound
                .unwrap_or_else(|| panic!("{what}: {deviation:?}"));
            eprintln!("{what}: reading {reading:e}, bound {bound:e}, volume {volume}");
            assert!(bound <= 1e-3, "{what}: {deviation:?}");
            if e == 0.0 {
                // The loops close exactly: the bound stays at rounding.
                assert!(bound <= 1e-11, "{what}: {deviation:?}");
            } else {
                // Within a few readings of the exact result: the slack holds
                // the gaps, no more.
                assert!(bound <= 10.0 * reading, "{what}: {deviation:?}");
            }
        }
    }
}
