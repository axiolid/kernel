//! An extruded polyline profile with exactly collinear boundary vertices
//! meshes closed whatever its placement (#415).
//!
//! The profile is the issue's 23-point outer ring: 9 collinear vertices on
//! the `x = 125.549...` side, 5 on the `y = -3.9625` side. Its caps are
//! triangulated in the profile's own plane, so a run that is exactly
//! collinear there stays collinear and the placement never moves a cap
//! vertex off a line it was on.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Interval, Point2, Scalar, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_mesh::audit_mesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

/// The issue's ring, `IFCCARTESIANPOINT` #100..#122, counter-clockwise.
const RING: [[f64; 2]; 23] = [
    [132.8, -8.0],
    [132.4, -8.0],
    [132.4, -7.60000000000037],
    [132.8, -7.60000000000037],
    [132.8, -3.96250000000036],
    [127.4, -3.96250000000034],
    [127.3375, -3.96250000000034],
    [127.2622304, -3.96250000000034],
    [125.549671244606, -3.96250000000034],
    [125.549671244606, -5.36699040000169],
    [125.549671244606, -5.52000000000169],
    [125.549671244606, -9.98250000000169],
    [125.549671244606, -10.1355096000017],
    [125.549671244606, -12.2419904000017],
    [125.549671244606, -12.3950000000017],
    [125.549671244606, -14.2200000000018],
    [125.549671244606, -14.9450000000018],
    [125.549671244606, -15.6000000000008],
    [127.4, -15.6000000000008],
    [127.4, -14.4530096000007],
    [127.4, -14.2375000000007],
    [127.5377696, -14.2375000000007],
    [132.8, -14.2375000000008],
];

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

fn polygon(ring: &[[f64; 2]]) -> Profile {
    let p: Vec<Point2> = ring.iter().map(|q| Point2::new(q[0], q[1])).collect();
    Profile::Contour(ContourProfile {
        outer: Contour::new(
            (0..p.len())
                .map(|i| line(p[i], p[(i + 1) % p.len()]))
                .collect(),
        ),
        holes: Vec::new(),
    })
}

/// The turn about `z` whose x axis is `(cos, sin, 0)`, then `translation`.
fn placement(x_axis: [f64; 2], translation: [f64; 3]) -> Transform3 {
    Transform3::from_cols(
        Vec3::new(x_axis[0], x_axis[1], 0.0),
        Vec3::new(-x_axis[1], x_axis[0], 0.0),
        Vec3::Z,
        Vec3::new(translation[0], translation[1], translation[2]),
    )
}

/// Extrudes `ring` by `depth` along `+z`, places it, compiles it, and
/// asserts the mesh is a closed two-manifold with no degenerate triangle.
fn assert_closed(ring: &[[f64; 2]], depth: Scalar, transform: Transform3, what: &str) {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(polygon(ring))).unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
        .unwrap();
    let root = b
        .push(GeometryNode::Instance(Instance {
            source: solid,
            transform,
        }))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let tolerance = Tolerance::MILLIMETRE;
    let outcome = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_reported(&graph, root, &ExecutionOptions::new(tolerance))
        .unwrap_or_else(|e| panic!("{what}: refused: {e}"));
    let health = audit_mesh(&outcome.mesh, tolerance);
    assert!(
        health.is_closed_two_manifold() && health.degenerate_triangles == 0,
        "{what}: {health:?}"
    );
}

const X_AXIS: [f64; 2] = [0.999215271103513, 0.0396086100934191];

#[test]
fn the_issue_placement_is_closed() {
    assert_closed(
        &RING,
        3.6,
        placement(
            X_AXIS,
            [
                643692.128712773,
                5649156.75663248,
                197.600000000001 + 16.0599999999987,
            ],
        ),
        "georeferenced",
    );
}

#[test]
fn a_local_frame_placement_is_closed() {
    assert_closed(
        &RING,
        3.6,
        placement(X_AXIS, [3.2, -1.7, 0.4]),
        "local frame",
    );
}

#[test]
fn the_unturned_profile_is_closed() {
    assert_closed(
        &RING,
        3.6,
        placement([1.0, 0.0], [643692.128712773, 5649156.75663248, 213.66]),
        "unturned",
    );
}

#[test]
fn a_turned_rectangle_is_closed() {
    let rect = [[0.0, 0.0], [5.0, 0.0], [5.0, 3.0], [0.0, 3.0]];
    assert_closed(
        &rect,
        3.6,
        placement(X_AXIS, [643692.128712773, 5649156.75663248, 213.66]),
        "rectangle",
    );
}

/// A small deterministic generator (SplitMix64), no dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// An L-shaped ring whose every side carries a random run of vertices on
/// it, in order, counter-clockwise. Some run vertices sit a few ulps off the
/// side, as the issue's `y = -3.9625` side does (`2e-13`).
fn ring_with_runs(rng: &mut Rng) -> Vec<[f64; 2]> {
    let corners = [
        [0.0, 0.0],
        [6.0, 0.0],
        [6.0, 2.0],
        [2.5, 2.0],
        [2.5, 5.0],
        [0.0, 5.0],
    ];
    let mut ring = Vec::new();
    for i in 0..corners.len() {
        let (a, b) = (corners[i], corners[(i + 1) % corners.len()]);
        ring.push(a);
        let mut ts: Vec<f64> = (0..(rng.next() * 8.0) as usize)
            .map(|_| 0.05 + 0.9 * rng.next())
            .collect();
        ts.sort_by(f64::total_cmp);
        ts.dedup();
        let jitter = if rng.next() < 0.5 { 2e-13 } else { 0.0 };
        for t in ts {
            let off = if rng.next() < 0.5 { jitter } else { 0.0 };
            ring.push([
                a[0] + t * (b[0] - a[0]) + off,
                a[1] + t * (b[1] - a[1]) + off,
            ]);
        }
    }
    ring
}

#[test]
fn a_randomised_sweep_of_turns_and_translations_is_closed() {
    let mut rng = Rng(415);
    for case in 0..200 {
        let theta = 2.0 * std::f64::consts::PI * (0.001 + 0.998 * rng.next());
        let reach = 10f64.powf(7.0 * rng.next());
        let t = [
            reach * (rng.next() - 0.5),
            reach * (rng.next() - 0.5),
            reach * (rng.next() - 0.5),
        ];
        let ring = ring_with_runs(&mut rng);
        assert_closed(
            &ring,
            3.6,
            placement([theta.cos(), theta.sin()], t),
            &format!("case {case}, theta {theta}, t {t:?}"),
        );
    }
}

#[test]
fn a_cap_that_must_be_a_sliver_is_refused_by_name() {
    // One triangle, doubled area 1e-8, flat to rounding: no triangulation
    // avoids it, and the extrusion refuses it by name. A sliver that is not
    // flat (5e-8 high) is still meshed, as a fine chord's caps are thin.
    let sliver = [[0.0, 0.0], [10.0, 0.0], [5.0, 1e-9]];
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(polygon(&sliver))).unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let graph = b.finish(vec![solid]).unwrap();
    let error = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_reported(&graph, solid, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .expect_err("a sliver cap is refused");
    assert!(error.to_string().contains("cap triangle"), "{error}");
}
