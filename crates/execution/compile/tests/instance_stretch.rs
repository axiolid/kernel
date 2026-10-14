//! An instance's deviation bound is its source's, grown by the transform's
//! stretch bound, and its source is meshed to the budget shrunk by that
//! stretch so the product stays within the request (#258).
//!
//! The stretch bound sits a hair above one even for an identity or a
//! rotation, so a rigid instance of a body meshed to exactly its budget
//! must still report at most that budget, and `meets_requested` must hold.
//! Scaled instances must stay sound: every case samples the exact surface
//! densely (closed form, or the reference evaluator, then moved by the
//! instance's transform) and asserts no sample is further from the mesh
//! than the reported bound.

use std::collections::HashMap;
use std::f64::consts::TAU;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{
    Frame2, Frame3, Mat3, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3,
};
use axiolid_curve::{BSplineCurve3, Circle2, Curve2, Curve3, KnotSpec, Polyline3};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationReport, ReferenceMeshCompiler};
use axiolid_model::{
    GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation, SurfaceRelation,
};
use axiolid_primitive::Primitive;
use axiolid_surface::{Plane, Surface};

const BUDGET: Scalar = 1e-3;

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(budget, 1e-9).unwrap())
        .with_chord_error(budget)
        .unwrap()
}

/// Triangles binned on a uniform grid, each in every cell its box touches
/// once grown by `reach`.
struct TriangleGrid {
    triangles: Vec<[Point3; 3]>,
    cells: HashMap<[i64; 3], Vec<usize>>,
    cell: Scalar,
    reach: Scalar,
}

impl TriangleGrid {
    fn new(mesh: &TriMesh, reach: Scalar) -> Self {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .chunks_exact(3)
            .map(|t| [0, 1, 2].map(|i| mesh.positions[t[i] as usize]))
            .collect();
        let (lo, hi) = mesh.positions.iter().fold(
            (
                Vec3::splat(Scalar::INFINITY),
                Vec3::splat(Scalar::NEG_INFINITY),
            ),
            |(lo, hi), p| (lo.min(*p), hi.max(*p)),
        );
        let cell = ((hi - lo).max_element() / 64.0).max(4.0 * reach);
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (i, t) in triangles.iter().enumerate() {
            let a = t[0].min(t[1]).min(t[2]) - Vec3::splat(reach);
            let b = t[0].max(t[1]).max(t[2]) + Vec3::splat(reach);
            let (ka, kb) = (key(a, cell), key(b, cell));
            for x in ka[0]..=kb[0] {
                for y in ka[1]..=kb[1] {
                    for z in ka[2]..=kb[2] {
                        cells.entry([x, y, z]).or_default().push(i);
                    }
                }
            }
        }
        Self {
            triangles,
            cells,
            cell,
            reach,
        }
    }

    fn distance(&self, p: Point3) -> Scalar {
        let mut best = self.reach;
        for &i in self.cells.get(&key(p, self.cell)).into_iter().flatten() {
            if let Ok(q) = closest_point_on_triangle(p, self.triangles[i]) {
                best = best.min(p.distance(q));
            }
        }
        best
    }
}

fn key(p: Point3, cell: Scalar) -> [i64; 3] {
    [p.x, p.y, p.z].map(|v| (v / cell).floor() as i64)
}

fn span(a: Scalar, b: Scalar, n: usize) -> impl Iterator<Item = Scalar> + Clone {
    (0..n).map(move |i| a + (b - a) * ((i as Scalar) + 0.5) / n as Scalar)
}

/// Compile the node `build` pushes at [`BUDGET`], with the report.
fn compile(build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId) -> (TriMesh, DeviationReport) {
    compile_at(BUDGET, build)
}

/// Compile the node `build` pushes at `budget`, with the report.
fn compile_at(
    budget: Scalar,
    build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId,
) -> (TriMesh, DeviationReport) {
    let mut builder = GeometryGraphBuilder::new();
    let root = build(&mut builder);
    let graph = builder.finish(vec![root]).expect("a valid graph");
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options(budget))
        .expect("the body compiles");
    (outcome.mesh, report)
}

fn instance(b: &mut GeometryGraphBuilder, source: NodeId, transform: Transform3) -> NodeId {
    b.push(GeometryNode::Instance(Instance { source, transform }))
        .unwrap()
}

/// The report's bound, which must exist; asserts no sample is further from
/// the mesh than it, and returns it.
fn assert_sound(
    name: &str,
    mesh: &TriMesh,
    report: &DeviationReport,
    samples: &[Point3],
) -> Scalar {
    let bound = report
        .bound
        .unwrap_or_else(|| panic!("{name}: expected a certified bound, got {report:?}"));
    let grid = TriangleGrid::new(mesh, 2.0 * bound);
    let worst = samples
        .iter()
        .map(|p| grid.distance(*p))
        .fold(0.0, Scalar::max);
    eprintln!(
        "{name}: reported {bound:e}, measured {worst:e}, {} samples",
        samples.len()
    );
    assert!(
        worst <= bound,
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
    assert!(worst > 0.0, "{name}: nothing measured");
    bound
}

/// Sound, within the requested budget, and certified to meet it.
fn assert_meets(name: &str, mesh: &TriMesh, report: &DeviationReport, samples: &[Point3]) {
    let bound = assert_sound(name, mesh, report, samples);
    let budget = report.requested;
    assert!(bound <= budget, "{name}: bound {bound:e} over {budget:e}");
    assert!(report.meets_requested(), "{name}: {report:?}");
}

fn moved(transform: Transform3, samples: &[Point3]) -> Vec<Point3> {
    samples
        .iter()
        .map(|&p| transform.transform_point3(p))
        .collect()
}

/// The rigid placements every body is instanced with: identity, a rotation
/// about a skew axis, a translation, a mirror, and a rotation placed far
/// from the origin.
fn rigid() -> Vec<(&'static str, Transform3)> {
    let axis = Vec3::new(1.0, 2.0, 3.0).normalize();
    vec![
        ("identity", Transform3::IDENTITY),
        ("rotation", Transform3::from_axis_angle(axis, 0.7)),
        (
            "translation",
            Transform3::from_translation(Vec3::new(3.25, -1.5, 0.75)),
        ),
        ("mirror", Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0))),
        (
            "far placement",
            Transform3::from_translation(Vec3::new(1520.3, -2741.9, 388.1))
                * Transform3::from_axis_angle(Vec3::new(-0.3, 0.2, 0.9).normalize(), 2.3),
        ),
    ]
}

/// Scaled placements: uniform, non-uniform along rotated axes, and a shear
/// (columns not orthogonal, so the stretch is the Frobenius norm).
fn scaled() -> Vec<(&'static str, Transform3)> {
    vec![
        ("uniform scale", Transform3::from_scale(Vec3::splat(2.5))),
        ("uniform shrink", Transform3::from_scale(Vec3::splat(0.4))),
        (
            "non-uniform scale",
            Transform3::from_translation(Vec3::new(10.0, 4.0, -2.0))
                * Transform3::from_rotation_z(0.4)
                * Transform3::from_scale(Vec3::new(1.5, 0.75, 2.0)),
        ),
        (
            "shear",
            Transform3::from_mat3(Mat3::from_cols(
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.3, 1.0, 0.0),
                Vec3::new(0.0, -0.2, 1.1),
            )),
        ),
    ]
}

// --- tori -------------------------------------------------------------------

const MAJOR: Scalar = 0.6;
const MINOR: Scalar = 0.2;

fn torus(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push(GeometryNode::Primitive(Primitive::Torus {
        major_radius: MAJOR,
        minor_radius: MINOR,
    }))
    .unwrap()
}

fn torus_samples() -> Vec<Point3> {
    let mut out = Vec::new();
    for u in span(0.0, TAU, 1500) {
        for v in span(0.0, TAU, 400) {
            let rho = MAJOR + MINOR * v.cos();
            out.push(Point3::new(rho * u.cos(), rho * u.sin(), MINOR * v.sin()));
        }
    }
    out
}

#[test]
fn the_bare_torus_reports_its_budget() {
    let (mesh, report) = compile(torus);
    assert_eq!(report.bound, Some(BUDGET));
    assert_meets("torus", &mesh, &report, &torus_samples());
}

#[test]
fn a_rigid_instance_of_a_torus_meets_its_budget() {
    let samples = torus_samples();
    for (name, transform) in rigid() {
        let (mesh, report) = compile(|b| {
            let body = torus(b);
            instance(b, body, transform)
        });
        assert!(
            matches!(report.contributions[0].bound, DeviationBound::Proven(_)),
            "{report:?}"
        );
        assert_meets(
            &format!("torus, {name}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
}

#[test]
fn a_scaled_instance_of_a_torus_stays_sound() {
    let samples = torus_samples();
    for (name, transform) in scaled() {
        let (mesh, report) = compile(|b| {
            let body = torus(b);
            instance(b, body, transform)
        });
        // The local budget shrinks by the stretch and the bound grows back
        // by it: within the request, but never below what the stretch
        // makes of the source's budget.
        assert_meets(
            &format!("torus, {name}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
}

/// Scale factors at which the budget divided by the stretch and multiplied
/// back rounds one ulp above it: the local budget must be shrunk until the
/// product fits, not merely divided.
#[test]
fn a_torus_scaled_where_the_round_trip_rounds_up_meets_its_budget() {
    let budget = 1e-4;
    let samples = torus_samples();
    for factor in [0.33, 1.19, 2.38] {
        let transform = Transform3::from_scale(Vec3::splat(factor));
        let (mesh, report) = compile_at(budget, |b| {
            let body = torus(b);
            instance(b, body, transform)
        });
        assert_meets(
            &format!("torus, scale {factor}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
}

#[test]
fn nested_instances_of_a_torus_meet_their_budget() {
    let samples = torus_samples();
    let inner = Transform3::from_axis_angle(Vec3::new(0.0, 1.0, 1.0).normalize(), 1.1);
    let middle = Transform3::from_scale(Vec3::splat(1.75));
    let outer = Transform3::from_translation(Vec3::new(-812.5, 97.25, 4410.0))
        * Transform3::from_rotation_x(0.9);
    let (mesh, report) = compile(|b| {
        let body = torus(b);
        let one = instance(b, body, inner);
        let two = instance(b, one, middle);
        instance(b, two, outer)
    });
    assert_meets(
        "torus, nested",
        &mesh,
        &report,
        &moved(outer * middle * inner, &samples),
    );
}

// --- curve-bounded planes -----------------------------------------------------

const DISK: Scalar = 0.5;

/// A plane over the world XY frame bounded by a circle of radius [`DISK`].
fn disk(b: &mut GeometryGraphBuilder) -> NodeId {
    let basis = b
        .push(GeometryNode::Surface(Surface::Plane(Plane {
            frame: Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        })))
        .unwrap();
    let circle = b
        .push(GeometryNode::Curve2(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: DISK,
        })))
        .unwrap();
    b.push(GeometryNode::SurfaceRelation(
        SurfaceRelation::CurveBounded {
            basis,
            boundaries: vec![circle],
            implicit_outer: false,
        },
    ))
    .unwrap()
}

fn disk_samples() -> Vec<Point3> {
    let mut out = Vec::new();
    for t in span(0.0, TAU, 20_000) {
        out.push(Point3::new(DISK * t.cos(), DISK * t.sin(), 0.0));
    }
    for r in span(0.0, DISK, 60) {
        for t in span(0.0, TAU, 600) {
            out.push(Point3::new(r * t.cos(), r * t.sin(), 0.0));
        }
    }
    out
}

/// Uninstanced, the plane's own frame stretch (a hair above one) is all
/// that grows its boundary's bound: an instance's source leaves a share of
/// its budget for the placement's rounding (#295), which would hide it.
#[test]
fn the_bare_curve_bounded_plane_meets_its_budget() {
    let (mesh, report) = compile(disk);
    assert_meets("disk", &mesh, &report, &disk_samples());
}

#[test]
fn rigid_and_nested_instances_of_a_curve_bounded_plane_meet_their_budget() {
    let samples = disk_samples();
    for (name, transform) in rigid() {
        let (mesh, report) = compile(|b| {
            let body = disk(b);
            instance(b, body, transform)
        });
        assert_meets(
            &format!("disk, {name}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
    let inner = Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7);
    let outer = Transform3::from_translation(Vec3::new(1520.3, -2741.9, 388.1))
        * Transform3::from_rotation_y(2.3);
    let (mesh, report) = compile(|b| {
        let body = disk(b);
        let one = instance(b, body, inner);
        instance(b, one, outer)
    });
    assert_meets(
        "disk, nested",
        &mesh,
        &report,
        &moved(outer * inner, &samples),
    );
}

#[test]
fn a_scaled_instance_of_a_curve_bounded_plane_stays_sound() {
    let samples = disk_samples();
    for (name, transform) in scaled() {
        let (mesh, report) = compile(|b| {
            let body = disk(b);
            instance(b, body, transform)
        });
        assert_meets(
            &format!("disk, {name}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
}

// --- pipes ----------------------------------------------------------------------

const PIPE: Scalar = 0.05;

/// A disk swept along a straight polyline: its exact tube is a capped
/// cylinder along x, proven at the budget (`axiolid_construct::pipe`).
fn straight_pipe(b: &mut GeometryGraphBuilder) -> NodeId {
    let directrix = b
        .push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
            points: vec![Point3::ZERO, Point3::new(0.6, 0.0, 0.0)],
            closed: false,
        })))
        .unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix,
        radius: PIPE,
        inner_radius: None,
        parameter_range: None,
        fillet_radius: None,
    }))
    .unwrap()
}

fn straight_pipe_samples() -> Vec<Point3> {
    let mut out = Vec::new();
    for x in span(0.0, 0.6, 600) {
        for t in span(0.0, TAU, 2000) {
            out.push(Point3::new(x, PIPE * t.cos(), PIPE * t.sin()));
        }
    }
    out
}

/// A disk swept along a cubic B-spline: certified by branch and bound.
fn spline() -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.3, 0.0, 0.1),
            Point3::new(0.5, 0.3, 0.0),
            Point3::new(0.8, 0.4, 0.2),
            Point3::new(1.0, 0.2, 0.3),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

fn spline_pipe(b: &mut GeometryGraphBuilder) -> NodeId {
    let directrix = b.push(GeometryNode::Curve3(spline())).unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix,
        radius: PIPE,
        inner_radius: None,
        parameter_range: None,
        fillet_radius: None,
    }))
    .unwrap()
}

/// Samples of the spline pipe's walls and end disks.
fn spline_pipe_samples() -> Vec<Point3> {
    let curve = spline();
    let frame = |t: Scalar| {
        let c = axiolid_reference::curve::evaluate3(&curve, t).unwrap();
        let d = axiolid_reference::curve::derivative3(&curve, t)
            .unwrap()
            .normalize();
        let a = if d.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
        let n = (a - d * d.dot(a)).normalize();
        (c, n, d.cross(n))
    };
    let mut out = Vec::new();
    for t in span(0.0, 1.0, 2000) {
        let (c, n, b) = frame(t);
        for theta in span(0.0, TAU, 360) {
            out.push(c + (n * theta.cos() + b * theta.sin()) * PIPE);
        }
    }
    for t in [0.0, 1.0] {
        let (c, n, b) = frame(t);
        for r in span(0.0, PIPE, 40) {
            for theta in span(0.0, TAU, 360) {
                out.push(c + (n * theta.cos() + b * theta.sin()) * r);
            }
        }
    }
    out
}

#[test]
fn rigid_instances_of_a_proven_pipe_meet_their_budget() {
    let samples = straight_pipe_samples();
    for (name, transform) in rigid() {
        let (mesh, report) = compile(|b| {
            let body = straight_pipe(b);
            instance(b, body, transform)
        });
        assert!(
            matches!(report.contributions[0].bound, DeviationBound::Proven(_)),
            "{report:?}"
        );
        assert_meets(
            &format!("straight pipe, {name}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
    for (name, transform) in scaled() {
        let (mesh, report) = compile(|b| {
            let body = straight_pipe(b);
            instance(b, body, transform)
        });
        assert_meets(
            &format!("straight pipe, {name}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
}

#[test]
fn rigid_and_nested_instances_of_a_certified_pipe_meet_their_budget() {
    let samples = spline_pipe_samples();
    let rotation = Transform3::from_axis_angle(Vec3::new(1.0, -1.0, 0.5).normalize(), 0.8);
    let far = Transform3::from_translation(Vec3::new(-3017.0, 1288.5, 64.25));
    for (name, transform, nested) in [
        ("identity", Transform3::IDENTITY, false),
        ("far placement", far * rotation, false),
        ("nested", far * rotation, true),
    ] {
        let (mesh, report) = compile(|b| {
            let body = spline_pipe(b);
            if nested {
                let one = instance(b, body, rotation);
                instance(b, one, far)
            } else {
                instance(b, body, transform)
            }
        });
        assert!(
            matches!(report.contributions[0].bound, DeviationBound::Certified(_)),
            "{report:?}"
        );
        assert_meets(
            &format!("spline pipe, {name}"),
            &mesh,
            &report,
            &moved(transform, &samples),
        );
    }
}
