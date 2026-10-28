//! A revolved contour whose arc ends on the axis meshes closed (#365).
//!
//! The arc's end is the joint with the next segment, stated exactly, so the
//! half- and quarter-disc of radius 0.5 weld to one vertex per pole instead
//! of landing beside the axis and pinching.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::volume_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, SolidOperation};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

const RADIUS: Scalar = 0.5;
const BUDGET: Scalar = 1e-3;

fn arc(from: Scalar, to: Scalar) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: RADIUS,
        }),
        domain: Interval::new(from, to),
        same_sense: true,
    }
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn contour(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

fn half_disc() -> Profile {
    contour(vec![
        arc(-FRAC_PI_2, FRAC_PI_2),
        line(Point2::new(0.0, RADIUS), Point2::new(0.0, -RADIUS)),
    ])
}

fn quarter_disc() -> Profile {
    contour(vec![
        arc(0.0, FRAC_PI_2),
        line(Point2::new(0.0, RADIUS), Point2::ZERO),
        line(Point2::ZERO, Point2::new(RADIUS, 0.0)),
    ])
}

/// Mesh `profile` revolved by `angle` about the Y axis: the volume, the
/// reported bound and the number of vertices exactly on the axis.
fn compile(profile: Profile, angle: Scalar) -> (Scalar, Scalar, usize) {
    let mut builder = GeometryGraphBuilder::new();
    let profile = builder.push(GeometryNode::Profile(profile)).unwrap();
    let root = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Revolution {
            profile,
            axis_origin: Point3::ZERO,
            axis_direction: Vec3::Y,
            angle,
        }))
        .unwrap();
    let graph = builder.finish(vec![root]).unwrap();
    let options = ExecutionOptions::new(Tolerance::new(BUDGET, 1e-9).unwrap());
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options)
        .expect("an arc profile touching the axis is meshed");
    // The profile's bound counts the distance an arc's end moved onto its
    // joint (3e-17 here), so it can top the budget by that rounding: the
    // report is "not certified", not "missed".
    let bound = report.bound.expect("a certified bound");
    assert!(bound <= report.requested * (1.0 + 1e-12), "{report:?}");
    let mesh = outcome.mesh;
    let exact = Tolerance::new(1e-9, 1e-9).unwrap();
    let health = axiolid_mesh::audit_mesh(&mesh, exact);
    assert!(health.is_closed_two_manifold(), "{health:?}");
    let volume = volume_properties(&mesh, exact)
        .expect("a closed solid")
        .signed_volume;
    let on_axis = mesh
        .positions
        .iter()
        .filter(|p| p.x == 0.0 && p.z == 0.0)
        .count();
    (volume, bound, on_axis)
}

#[test]
fn a_half_disc_arc_revolves_to_a_sphere_and_its_turns() {
    let sphere = 4.0 / 3.0 * PI * RADIUS.powi(3);
    for (angle, share) in [(TAU, 1.0), (PI, 0.5)] {
        let (volume, bound, on_axis) = compile(half_disc(), angle);
        assert_eq!(on_axis, 2, "angle {angle}");
        // Inscribed in the exact solid, and within the reported bound of
        // its surface: the shell between them is at most `bound` thick.
        let want = sphere * share;
        let shell = 4.0 * PI * RADIUS * RADIUS * bound * share;
        assert!(
            volume <= want && want - volume <= shell,
            "angle {angle}: {volume} vs {want}"
        );
    }
}

#[test]
fn a_quarter_disc_arc_revolves_to_a_hemisphere_and_its_turns() {
    let hemisphere = 2.0 / 3.0 * PI * RADIUS.powi(3);
    for (angle, share) in [(TAU, 1.0), (PI, 0.5)] {
        let (volume, bound, on_axis) = compile(quarter_disc(), angle);
        // The pole and the centre of the flat face.
        assert_eq!(on_axis, 2, "angle {angle}");
        let want = hemisphere * share;
        let shell = 2.0 * PI * RADIUS * RADIUS * bound * share;
        assert!(
            volume <= want && want - volume <= shell,
            "angle {angle}: {volume} vs {want}"
        );
    }
}
