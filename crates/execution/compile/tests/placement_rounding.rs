//! An instance's deviation bound includes the rounding of applying its
//! `f64` transform far from the origin (#295).
//!
//! A torus instanced at `(6e6, 6e6, 100)` has its vertices moved by about
//! `5e-10` per coordinate by the placement itself: an absolute error no
//! relative stretch covers. The report grows every bound by `4 eps M`
//! (`deviation::PLACEMENT_EPSILONS`, `M` the largest coordinate magnitude
//! the placement sums), names it under `DeviationPath::Placement`, meshes
//! the source with room for it so the request is still met, and refuses a
//! budget at or below that floor by name. Soundness is checked in `f64`
//! against the mesh as placed: its positions are moved back by the
//! translation, which is exact there (Sterbenz), and compared with the
//! exact torus sampled densely at the origin.

use std::f64::consts::TAU;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::{
    BUDGET_BELOW_PLACEMENT_ROUNDING, PLACEMENT_EPSILONS, PLACEMENT_ROUNDING,
};
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, InstanceAtStation,
    NodeId, Station, StationFrame, StationOffsets,
};
use axiolid_primitive::Primitive;

/// A small torus, so a budget a few times the placement's rounding still
/// meshes in about a hundred thousand triangles.
const MAJOR: Scalar = 6e-4;
const MINOR: Scalar = 2e-4;

/// The issue's placement.
const SITE: Vec3 = Vec3::new(6e6, 6e6, 100.0);

/// The rounding floor at [`SITE`]: `4 eps` of its largest component.
fn floor() -> Scalar {
    PLACEMENT_EPSILONS * Scalar::EPSILON * 6e6
}

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(budget, 1e-9).unwrap())
        .with_chord_error(budget)
        .unwrap()
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn torus(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push(GeometryNode::Primitive(Primitive::Torus {
        major_radius: MAJOR,
        minor_radius: MINOR,
    }))
    .unwrap()
}

fn instance(b: &mut GeometryGraphBuilder, source: NodeId, transform: Transform3) -> NodeId {
    b.push(GeometryNode::Instance(Instance { source, transform }))
        .unwrap()
}

/// A torus instanced by `transforms`, innermost first.
fn placed_torus(transforms: &[Transform3]) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let mut root = torus(&mut b);
    for &transform in transforms {
        root = instance(&mut b, root, transform);
    }
    (b.finish(vec![root]).unwrap(), root)
}

fn compile(
    graph: &GeometryGraph,
    root: NodeId,
    budget: Scalar,
) -> GeomResult<(TriMesh, DeviationReport)> {
    compiler()
        .compile_mesh_with_deviation(graph, root, &options(budget))
        .map(|(outcome, report)| (outcome.mesh, report))
}

fn assert_refused<T: core::fmt::Debug>(result: GeomResult<T>) {
    match result {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert_eq!(input, BUDGET_BELOW_PLACEMENT_ROUNDING);
        }
        other => panic!("expected the placement rounding refusal, got {other:?}"),
    }
}

fn contribution(report: &DeviationReport, path: DeviationPath) -> Option<DeviationBound> {
    report
        .contributions
        .iter()
        .find(|c| c.path == path)
        .map(|c| c.bound)
}

/// The placement's rounding, as the report names it.
fn rounding(report: &DeviationReport) -> Scalar {
    let named = report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::Placement && c.detail == PLACEMENT_ROUNDING);
    match named.map(|c| c.bound) {
        Some(DeviationBound::Certified(value)) => value,
        other => panic!("no placement rounding named: {other:?} in {report:?}"),
    }
}

fn span(n: usize) -> impl Iterator<Item = Scalar> + Clone {
    (0..n).map(move |i| TAU * ((i as Scalar) + 0.5) / n as Scalar)
}

/// The exact torus at the origin, densely.
fn samples() -> Vec<Point3> {
    let mut out = Vec::new();
    for u in span(800) {
        for v in span(400) {
            let rho = MAJOR + MINOR * v.cos();
            out.push(Point3::new(rho * u.cos(), rho * u.sin(), MINOR * v.sin()));
        }
    }
    out
}

/// The farthest any sample lies from the mesh. Triangles are binned on a
/// grid, each in every cell its box grown by `REACH` touches, so a
/// sample's own cell holds every triangle within `REACH` of it.
fn worst_distance(mesh: &TriMesh, samples: &[Point3]) -> Scalar {
    use std::collections::HashMap;
    const REACH: Scalar = 1e-7;
    let cell = 2e-5;
    let key = |p: Point3| [p.x, p.y, p.z].map(|v| (v / cell).floor() as i64);
    let triangles: Vec<[Point3; 3]> = mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| [0, 1, 2].map(|i| mesh.positions[t[i] as usize]))
        .collect();
    let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    let grow = Vec3::splat(REACH);
    for (i, t) in triangles.iter().enumerate() {
        let a = key(t[0].min(t[1]).min(t[2]) - grow);
        let b = key(t[0].max(t[1]).max(t[2]) + grow);
        for x in a[0]..=b[0] {
            for y in a[1]..=b[1] {
                for z in a[2]..=b[2] {
                    cells.entry([x, y, z]).or_default().push(i);
                }
            }
        }
    }
    let mut worst: Scalar = 0.0;
    for &p in samples {
        let mut best = REACH;
        for &i in cells.get(&key(p)).into_iter().flatten() {
            if let Ok(q) = closest_point_on_triangle(p, triangles[i]) {
                best = best.min(p.distance(q));
            }
        }
        worst = worst.max(best);
    }
    worst
}

/// The mesh as placed at [`SITE`], moved back to the origin: each
/// coordinate within a factor two of the translation's, so the subtraction
/// is exact and the result is the placed mesh's own position, rounding and
/// all.
fn moved_back(mesh: &TriMesh) -> TriMesh {
    let positions = mesh
        .positions
        .iter()
        .map(|&p| {
            let back = p - SITE;
            assert_eq!(back + SITE, p, "the subtraction is exact");
            back
        })
        .collect();
    TriMesh::new(positions, mesh.indices.clone())
}

#[test]
fn a_torus_placed_at_six_million_reports_a_bound_that_holds_in_f64() {
    let budget = 2e-8;
    let (graph, root) = placed_torus(&[Transform3::from_translation(SITE)]);
    let (mesh, report) = compile(&graph, root, budget).expect("the budget is above the floor");

    // The placement's rounding is named, at `4 eps` of the placed
    // coordinates: about a quarter of the budget here.
    let r = rounding(&report);
    assert!(r >= floor() && r <= floor() * (1.0 + 1e-9), "{r:e}");
    // The torus's own bound carries it on top of its chord budget, and the
    // source was meshed with room for it: within the request, and above
    // what the chord budget alone would have reached.
    let primitive = contribution(&report, DeviationPath::Primitive)
        .and_then(DeviationBound::value)
        .expect("a proven torus bound");
    assert!(primitive <= budget, "{primitive:e}");
    assert!(
        primitive > budget - 0.5 * r,
        "the bound adds the rounding: {primitive:e}"
    );
    assert_eq!(report.bound, Some(primitive));
    assert!(report.meets_requested(), "{report:?}");

    // Dense samples of the exact torus, against the placed mesh moved back
    // exactly: the bound holds for the f64 mesh.
    let worst = worst_distance(&moved_back(&mesh), &samples());
    assert!(worst <= primitive, "{worst:e} over {primitive:e}");
    assert!(worst > 0.25 * budget, "nothing measured: {worst:e}");
}

#[test]
fn a_budget_below_the_placement_rounding_is_refused_by_name() {
    let (graph, root) = placed_torus(&[Transform3::from_translation(SITE)]);
    // The issue's budget: below `4 eps 6e6`, about 5.3e-9.
    assert!(1e-9 < floor());
    assert_refused(compile(&graph, root, 1e-9));
    // At the floor exactly, and through the plain compile path too.
    assert_refused(compile(&graph, root, floor()));
    assert_refused(compiler().compile_mesh(&graph, root, &options(1e-9)));
    // Placed near the origin, a budget below the far floor compiles and
    // is met.
    let near = Transform3::from_translation(Vec3::new(0.0, 0.0, 100.0));
    let (graph, root) = placed_torus(&[near]);
    let (_, report) = compile(&graph, root, 0.5 * floor()).expect("near the origin");
    assert!(report.meets_requested(), "{report:?}");
}

#[test]
fn a_source_far_from_the_origin_turned_there_is_meshed_with_room_for_its_rounding() {
    // The inner instance places the torus at the site; the outer one turns
    // it about the origin, with no translation, so the room left before
    // the source is meshed is only a sliver, and the source's own
    // coordinates round by about 7e-9 once turned: the source is meshed
    // once more with room for that.
    let budget = 3e-8;
    let turn = Transform3::from_rotation_z(0.3);
    let (graph, root) = placed_torus(&[Transform3::from_translation(SITE), turn]);
    let (mesh, report) = compile(&graph, root, budget).expect("compiles");
    let bound = report.bound.expect("bounded");
    assert!(bound <= budget, "{report:?}");
    assert!(report.meets_requested(), "{report:?}");
    // Both placements' rounding: the outer's at least `4 eps` of the turned
    // coordinates' summed magnitudes, the inner's grown by the turn.
    let turned = PLACEMENT_EPSILONS * Scalar::EPSILON * 6e6 * (0.3_f64.cos() + 0.3_f64.sin());
    assert!(rounding(&report) >= turned + floor(), "{report:?}");

    // Sampled in f64: the samples are placed by the same transforms, so
    // each is off by up to the same rounding; the bound must hold up to it.
    let world = turn * Transform3::from_translation(SITE);
    let sampled: Vec<Point3> = samples()
        .into_iter()
        .step_by(7)
        .map(|p| world.transform_point3(p))
        .collect();
    let back = TriMesh::new(mesh.positions.clone(), mesh.indices.clone());
    let worst = worst_distance(&back, &sampled);
    assert!(
        worst <= bound + 2.0 * rounding(&report),
        "{worst:e} over {bound:e}"
    );
}

#[test]
fn an_identity_or_a_mirror_adds_no_rounding() {
    let mirror = Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0));
    let (graph, root) = placed_torus(&[Transform3::IDENTITY, mirror]);
    let (_, report) = compile(&graph, root, 1e-5).expect("compiles");
    assert_eq!(
        contribution(&report, DeviationPath::Placement),
        None,
        "{report:?}"
    );
    assert!(report.meets_requested(), "{report:?}");
}

/// A plan line along `+y` from `(6e6, 6e6)`, its parameter the distance.
fn far_line(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::new(6e6, 6e6),
        direction: Vec2::Y,
    }))
    .unwrap()
}

/// A torus placed at a station `100` up a far line: an exact frame, so
/// only the placement's rounding is reported, and refused by name below it.
#[test]
fn a_torus_placed_at_a_far_station_carries_the_rounding_too() {
    let mut b = GeometryGraphBuilder::new();
    let basis = far_line(&mut b);
    let body = torus(&mut b);
    let station = CurveStation {
        basis,
        station: Station::new(2.0, StationOffsets::new(0.0, 100.0, 0.0)),
        frame: StationFrame::Section,
    };
    let placed = b
        .push_value(InstanceAtStation::new(body, station.into()))
        .unwrap();
    let graph = b.finish(vec![placed]).unwrap();

    let budget = 2e-8;
    let (_, report) = compile(&graph, placed, budget).expect("above the floor");
    assert_eq!(contribution(&report, DeviationPath::StationPlacement), None);
    let r = rounding(&report);
    assert!(r >= floor(), "{r:e}");
    let primitive = contribution(&report, DeviationPath::Primitive)
        .and_then(DeviationBound::value)
        .expect("a proven torus bound");
    assert!(
        primitive > budget - 0.5 * r && primitive <= budget,
        "{report:?}"
    );
    assert!(report.meets_requested(), "{report:?}");

    assert_refused(compile(&graph, placed, 1e-9));
}
