//! A wall less flush openings on a national grid measures its closed-form
//! volume however the model's length unit rounds (#343).
//!
//! The wall is `4.7 x 0.25 x 2.85` less a door and a window, each
//! `1.01 x 2.26` and cut flush through it: their caps lie in the wall's
//! faces, and the door's sill in its bottom face. The placement is the
//! one #343 found the jumps at: about `(6.0e5, 5.6e6, 204.5)`, turned
//! `0.0396` about `z`. Every length, the placement's included, is scaled
//! by the length unit `s`, as a model read in millimetres is; `s` one ulp
//! off `1` is the same model with its unit factor rounded the other way.
//!
//! The mesh and the exact boolean were right at every `s`. The exact
//! measure was not: it integrated about the world origin, `5.6e6` m away,
//! and read `2.314 m^3` for this `2.2074 m^3` wall at `s = 1`; the
//! IFC fixture #343 reduces read `2.14` to `2.27 m^3`, a different number
//! in each rounding of its unit.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Interval, Point2, Scalar, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_measure::{exact_properties, volume_properties};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{ReferenceExactCompiler, ReferenceMeshCompiler};
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

/// The wall less its openings in length unit `s`, its closed-form volume,
/// and the area of every operand's faces.
fn wall(s: Scalar) -> (GeometryGraph, NodeId, Scalar, Scalar) {
    let (length, thickness, height) = (4.7 * s, 0.25 * s, 2.85 * s);
    let (width, tall) = (1.01 * s, 2.26 * s);
    let world = Transform3::from_translation(
        Vec3::new(600_013.917_809_576_6, 5_599_985.940_231_986, 204.5) * s,
    ) * Transform3::from_rotation_z(0.039_608_610_093_419_115);
    // An opening's profile lies in the wall's elevation and is cut along
    // `-y` from the wall's back face through to its front face.
    let across = Transform3::from_cols(
        Vec3::X,
        Vec3::Z,
        Vec3::NEG_Y,
        Vec3::new(0.0, thickness, 0.0),
    );
    let mut builder = GeometryGraphBuilder::new();
    let mut root = placed(
        &mut builder,
        rect(0.0, 0.0, length, thickness),
        height,
        world,
    );
    let mut area = 2.0 * (length * thickness + length * height + thickness * height);
    for (x, z) in [(0.5 * s, 0.0), (2.9 * s, 0.3 * s)] {
        let opening = placed(
            &mut builder,
            rect(x, z, x + width, z + tall),
            thickness,
            world * across,
        );
        root = builder
            .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: root,
                right: opening,
                operator: BooleanOperator::Difference,
            }))
            .unwrap();
        area += 2.0 * (width * tall + (width + tall) * thickness);
    }
    let volume = (length * height - 2.0 * width * tall) * thickness;
    (builder.finish(vec![root]).unwrap(), root, volume, area)
}

#[test]
fn a_wall_with_flush_openings_on_a_national_grid_measures_its_closed_form() {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let mut per_unit = Vec::new();
    for s in [1.0, 1.0 + Scalar::EPSILON, 1.0 - Scalar::EPSILON / 2.0] {
        let (graph, root, want, area) = wall(s);
        let (exact, report) = ReferenceExactCompiler::new()
            .compile_exact_with_report(&graph, root, &options)
            .unwrap_or_else(|error| panic!("s = {s}: {error}"));
        let volume = exact_properties(&exact, options.tolerance())
            .expect("measurable")
            .signed_volume;
        // The input's grid: at `5.6e6` adjacent f64 values are `9.3e-10`
        // apart, and placing each operand moves its faces by a few of them.
        let grid = area * 4.0 * 5.6e6 * s * Scalar::EPSILON;
        // Against the closed form, the boolean's readings and rounding
        // floor move the faces too (#244).
        let moved = grid + area * (report.linear() + report.rounding_floor());
        assert!(
            (volume - want).abs() <= moved,
            "s = {s}: expected {want} within {moved:e}, got {volume}"
        );
        // The mesh lies within its certified bound of the exact surface.
        let (outcome, deviation) = ReferenceMeshCompiler::new(BoolmeshBoolean)
            .compile_mesh_with_deviation(&graph, root, &options)
            .unwrap_or_else(|error| panic!("s = {s}: {error}"));
        let bound = deviation.bound.expect("certified");
        let meshed = volume_properties(&outcome.mesh, options.tolerance())
            .expect("closed")
            .signed_volume;
        let apart = grid + area * bound;
        assert!(
            (volume - meshed).abs() <= apart,
            "s = {s}: exact {volume} vs mesh {meshed}, within {apart:e}"
        );
        per_unit.push((s, volume / (s * s * s), grid / want));
    }
    // One model in three roundings of its unit: one volume per unit^3, up
    // to the grid. Before #343 the three differed by up to 5 %.
    let (_, first, _) = per_unit[0];
    for &(s, volume, relative) in &per_unit[1..] {
        assert!(
            (volume - first).abs() <= 2.0 * relative * first,
            "s = {s}: {volume} per unit^3 against {first} at s = 1"
        );
    }
}
