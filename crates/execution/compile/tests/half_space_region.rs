//! Intersections of half-spaces as clip tools (#326).
//!
//! A region is a half-space, or an intersection of regions, read through
//! instances. It is the tool of a difference or an intersection, and the
//! left operand of an intersection with a solid. The wall is
//! `6 x 0.3 x 3`, `x` in `[-3, 3]`, `y` in `[-0.15, 0.15]`, `z` in `[0, 3]`
//! in its own frame, under a rigid placement; `H1` is the side above
//! `z = 2 + 0.2 x`, `H2` the side `x <= 1`. Then `H1 & H2` takes
//! `0.3 * integral_{-3}^{1} (1 - 0.2 x) dx = 1.44` from the wall, whose
//! volume is `5.4`.
//!
//! Each case is checked against its closed form in the exact compiler and
//! the mesh compiler, and the mesh result against the exact one by the
//! deviation report.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Plane3, Point2, Point3, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Polyline2};
use axiolid_measure::exact_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::MESH_TO_EXACT;
use axiolid_mesh_compile::{DeviationBound, ReferenceExactCompiler, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Profile, RectangleProfile};

const L: f64 = 6.0;
const T: f64 = 0.3;
const H: f64 = 3.0;
const WALL: f64 = L * T * H;
/// What `H1 & H2` takes from the wall.
const CORNER: f64 = 1.44;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

fn site() -> Transform3 {
    Transform3::from_translation(Vec3::new(6.0e5, 5.6e6, 200.0)) * Transform3::from_rotation_z(0.6)
}

struct Graph {
    builder: GeometryGraphBuilder,
}

impl Graph {
    fn new() -> Self {
        Self {
            builder: GeometryGraphBuilder::new(),
        }
    }

    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.builder.push(node).expect("a valid node")
    }

    fn instance(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }

    fn minus(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.boolean(left, right, BooleanOperator::Difference)
    }

    fn meet(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.boolean(left, right, BooleanOperator::Intersection)
    }

    fn local_wall(&mut self) -> NodeId {
        let profile = self.push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: L,
                y: T,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: H,
        }))
    }

    fn wall(&mut self, placement: Transform3) -> NodeId {
        let wall = self.local_wall();
        self.instance(wall, placement)
    }

    /// The side `normal` points to of the plane through `origin`, both in
    /// the wall's frame, in world coordinates under `placement`.
    fn side(&mut self, origin: Point3, normal: Vec3, placement: Transform3) -> NodeId {
        self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: placement.transform_point3(origin),
                normal: placement.transform_vector3(normal),
            },
            agreement: true,
        }))
    }

    /// The same, in the wall's frame.
    fn local_side(&mut self, origin: Point3, normal: Vec3) -> NodeId {
        self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 { origin, normal },
            agreement: true,
        }))
    }

    /// `z >= 2 + 0.2 x`.
    fn h1(&mut self, placement: Transform3) -> NodeId {
        self.side(Point3::new(0.0, 0.0, 2.0), h1_normal(), placement)
    }

    /// `x <= 1`.
    fn h2(&mut self, placement: Transform3) -> NodeId {
        self.side(Point3::new(1.0, 0.0, 0.0), -Vec3::X, placement)
    }

    /// The side `z >= 0` within the square `[-1, 1]^2`.
    fn square(&mut self, placement: Transform3) -> NodeId {
        let half_space = self.local_side(Point3::ZERO, Vec3::Z);
        let boundary = self.push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::new(-1.0, -1.0),
                Point2::new(1.0, -1.0),
                Point2::new(1.0, 1.0),
                Point2::new(-1.0, 1.0),
            ],
            closed: true,
        })));
        let bounded = self.push(GeometryNode::SolidOperation(
            SolidOperation::BoundedHalfSpace {
                half_space,
                boundary,
                placement: Transform3::IDENTITY,
            },
        ));
        self.instance(bounded, placement)
    }

    fn finish(self, roots: Vec<NodeId>) -> GeometryGraph {
        self.builder.finish(roots).expect("a valid graph")
    }
}

fn h1_normal() -> Vec3 {
    Vec3::new(-0.2, 0.0, 1.0)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> Result<(ExactBRep, bool), GeomError> {
    ReferenceExactCompiler::new()
        .compile_exact_with_report(graph, root, &options())
        .map(|(body, report)| (body, report.is_exact()))
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

/// The exact result of `root`, audited, within `tolerance` (relative) of
/// `expected`; and whether its report is exact.
fn check(g: Graph, root: NodeId, expected: f64, tolerance: f64) -> (ExactBRep, bool) {
    let graph = g.finish(vec![root]);
    let (body, exact) = compile(&graph, root).expect("compiles exactly");
    let health = geometric_audit(&body, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let measured = volume(&body);
    assert!(
        (measured - expected).abs() <= tolerance * expected,
        "volume {measured}, expected {expected}"
    );
    (body, exact)
}

fn near(g: Graph, root: NodeId, expected: f64) -> (ExactBRep, bool) {
    check(g, root, expected, 1e-9)
}

fn refusal(g: Graph, root: NodeId) -> GeomError {
    let graph = g.finish(vec![root]);
    compile(&graph, root).expect_err("refused")
}

fn named(error: &GeomError, part: &str) {
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains(part)),
        "{error:?} must name `{part}`"
    );
}

/// A case: a graph, its root and the volume it must have.
struct Case {
    graph: Graph,
    root: NodeId,
    volume: f64,
}

/// `wall - (H1 & H2)`.
fn corner_removed(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let both = g.meet(h1, h2);
    let root = g.minus(wall, both);
    Case {
        graph: g,
        root,
        volume: WALL - CORNER,
    }
}

/// `wall & (H1 & H2)`.
fn corner_kept(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let both = g.meet(h1, h2);
    let root = g.meet(wall, both);
    Case {
        graph: g,
        root,
        volume: CORNER,
    }
}

/// `(H1 & H2) & wall`: the region as the left operand.
fn corner_kept_swapped(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let both = g.meet(h1, h2);
    let root = g.meet(both, wall);
    Case {
        graph: g,
        root,
        volume: CORNER,
    }
}

/// `H1 & wall`: a single half-space as the left operand.
fn roof_kept_swapped(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let root = g.meet(h1, wall);
    Case {
        graph: g,
        root,
        volume: 1.8,
    }
}

/// `wall - Instance(H1 & H2, placement)`, the leaves in the wall's frame.
fn region_under_an_instance(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.local_side(Point3::new(0.0, 0.0, 2.0), h1_normal());
    let h2 = g.local_side(Point3::new(1.0, 0.0, 0.0), -Vec3::X);
    let both = g.meet(h1, h2);
    let placed = g.instance(both, placement);
    let root = g.minus(wall, placed);
    Case {
        graph: g,
        root,
        volume: WALL - CORNER,
    }
}

/// `wall - ((H1 & H2) & {y <= 0})`: a region of three, nested.
fn region_of_three(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let corner = g.meet(h1, h2);
    let back = g.side(Point3::ZERO, -Vec3::Y, placement);
    let region = g.meet(corner, back);
    let root = g.minus(wall, region);
    Case {
        graph: g,
        root,
        volume: WALL - CORNER / 2.0,
    }
}

/// `wall - ({z >= 1} & {z <= 2})`: a band across the wall.
fn band_removed(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let above = g.side(Point3::new(0.0, 0.0, 1.0), Vec3::Z, placement);
    let below = g.side(Point3::new(0.0, 0.0, 2.0), -Vec3::Z, placement);
    let band = g.meet(above, below);
    let root = g.minus(wall, band);
    Case {
        graph: g,
        root,
        volume: WALL - L * T,
    }
}

/// `wall - (H1 & {x >= 10})`: the second leaf misses the wall.
fn leaf_beyond_the_wall(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let far = g.side(Point3::new(10.0, 0.0, 0.0), Vec3::X, placement);
    let region = g.meet(h1, far);
    let root = g.minus(wall, region);
    Case {
        graph: g,
        root,
        volume: WALL,
    }
}

/// `wall - (H1 & B)` or `wall - (B & H1)`, `B` a half-space bounded to
/// the square `[-1, 1]^2`.
fn bounded_leaf(placement: Transform3, bounded_first: bool) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let square = g.square(placement);
    let region = if bounded_first {
        g.meet(square, h1)
    } else {
        g.meet(h1, square)
    };
    let root = g.minus(wall, region);
    Case {
        graph: g,
        root,
        volume: WALL - 2.0 * T,
    }
}

/// `wall - ({z >= -5} & H1)`: the first leaf holds the whole wall.
fn leaf_holding_the_wall(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let floor = g.side(Point3::new(0.0, 0.0, -5.0), Vec3::Z, placement);
    let h1 = g.h1(placement);
    let region = g.meet(floor, h1);
    let root = g.minus(wall, region);
    Case {
        graph: g,
        root,
        volume: WALL - 1.8,
    }
}

/// `wall - ({z >= 2 + x/2} & {z <= 1 + x/2})`: parallel planes a unit
/// apart, each crossing the wall, with nothing in common.
fn stand_ins_with_nothing_in_common(placement: Transform3) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let above = g.side(
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(-0.5, 0.0, 1.0),
        placement,
    );
    let below = g.side(
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.5, 0.0, -1.0),
        placement,
    );
    let region = g.meet(above, below);
    let root = g.minus(wall, region);
    Case {
        graph: g,
        root,
        volume: WALL,
    }
}

fn every_case() -> Vec<(&'static str, Case)> {
    let p = building();
    vec![
        ("corner removed", corner_removed(p)),
        ("corner kept", corner_kept(p)),
        ("corner kept, swapped", corner_kept_swapped(p)),
        ("roof kept, swapped", roof_kept_swapped(p)),
        ("region under an instance", region_under_an_instance(p)),
        ("region of three", region_of_three(p)),
        ("band removed", band_removed(p)),
        ("leaf beyond the wall", leaf_beyond_the_wall(p)),
        ("bounded leaf second", bounded_leaf(p, false)),
        ("bounded leaf first", bounded_leaf(p, true)),
        ("leaf holding the wall", leaf_holding_the_wall(p)),
        (
            "stand-ins with nothing in common",
            stand_ins_with_nothing_in_common(p),
        ),
    ]
}

#[test]
fn the_exact_compiler_clips_by_every_region() {
    for (name, case) in every_case() {
        let Case {
            graph,
            root,
            volume,
        } = case;
        near_named(name, graph, root, volume);
    }
}

fn near_named(name: &str, g: Graph, root: NodeId, expected: f64) -> (ExactBRep, bool) {
    let graph = g.finish(vec![root]);
    let (body, exact) = compile(&graph, root).unwrap_or_else(|e| panic!("{name}: {e:?}"));
    let health = geometric_audit(&body, Tolerance::METRE);
    assert!(health.is_consistent(), "{name}: {:?}", health.defects());
    let measured = volume(&body);
    assert!(
        (measured - expected).abs() <= 1e-9 * expected,
        "{name}: volume {measured}, expected {expected}"
    );
    (body, exact)
}

#[test]
fn a_wall_less_two_half_spaces_is_exact_and_whole() {
    let Case { graph, root, .. } = corner_removed(building());
    let (solid, exact) = near(graph, root, WALL - CORNER);
    assert!(exact, "the report reads exact");
    assert_eq!(solid.topology().shells().len(), 1);
}

#[test]
fn a_wall_at_a_site_less_two_half_spaces_is_exact() {
    // The volume is read about the far origin: a loose bound.
    let Case { graph, root, .. } = corner_removed(site());
    let graph = graph.finish(vec![root]);
    let (solid, exact) = compile(&graph, root).expect("compiles exactly");
    assert!(
        (volume(&solid) - (WALL - CORNER)).abs() <= 1e-7,
        "{}",
        volume(&solid)
    );
    assert!(exact);
}

#[test]
fn a_band_cuts_the_wall_in_two() {
    let Case { graph, root, .. } = band_removed(building());
    let (solid, exact) = near(graph, root, WALL - L * T);
    assert!(exact);
    assert_eq!(solid.topology().shells().len(), 2);
}

#[test]
fn a_region_with_a_leaf_missing_the_wall_keeps_it_or_empties_the_intersection() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let far = g.side(Point3::new(10.0, 0.0, 0.0), Vec3::X, placement);
    let region = g.meet(h1, far);
    let root = g.meet(wall, region);
    let error = refusal(g, root);
    assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
}

#[test]
fn an_intersection_commutes_with_a_region() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let swapped = g.meet(h1, wall);
    let plain = g.meet(wall, h1);
    let graph = g.finish(vec![swapped, plain]);
    let (swapped_body, _) = compile(&graph, swapped).expect("compiles");
    let (plain_body, _) = compile(&graph, plain).expect("compiles");
    assert!((volume(&swapped_body) - 1.8).abs() <= 1e-9 * 1.8);
    assert_eq!(swapped_body, plain_body);
}

#[test]
fn a_region_with_no_solid_to_clip_is_refused_by_name() {
    let placement = building();

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let flipped = g.minus(h1, wall);
    named(&refusal(g, flipped), "subject is a half-space");

    let mut g = Graph::new();
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let alone = g.meet(h1, h2);
    named(
        &refusal(g, alone),
        "intersection of half-spaces without a subject",
    );

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let region = g.meet(h1, h2);
    let united = g.boolean(wall, region, BooleanOperator::Union);
    named(&refusal(g, united), "union with a half-space");

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let united = g.boolean(h1, wall, BooleanOperator::Union);
    named(&refusal(g, united), "union with a half-space");

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let region = g.meet(h1, h2);
    let xor = g.boolean(wall, region, BooleanOperator::SymmetricDifference);
    named(&refusal(g, xor), "symmetric difference with a half-space");
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

/// The mesh's volume about its first vertex, so far from the origin it is
/// not lost to cancellation.
fn mesh_volume(mesh: &TriMesh) -> f64 {
    let base = mesh.positions[0];
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize] - base);
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

#[test]
fn the_mesh_compiler_clips_by_every_region_and_is_certified() {
    for (name, case) in every_case() {
        let Case {
            graph,
            root,
            volume,
        } = case;
        let graph = graph.finish(vec![root]);
        let (outcome, report) = compiler()
            .compile_mesh_with_deviation(&graph, root, &options())
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        let measured = mesh_volume(&outcome.mesh);
        assert!(
            (measured - volume).abs() <= 1e-7 * volume,
            "{name}: volume {measured}, expected {volume}"
        );
        assert!(report.meets_requested(), "{name}: {report:?}");
    }
}

#[test]
fn a_mesh_clip_by_two_half_spaces_is_measured_both_ways() {
    let Case { graph, root, .. } = corner_removed(building());
    let graph = graph.finish(vec![root]);
    let (_, report) = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .expect("compiles");
    assert!(
        report
            .contributions
            .iter()
            .any(|c| c.detail == MESH_TO_EXACT && matches!(c.bound, DeviationBound::Certified(_))),
        "{report:?}"
    );
}

#[test]
fn the_mesh_compiler_clips_at_a_site() {
    let Case { graph, root, .. } = corner_removed(site());
    let graph = graph.finish(vec![root]);
    let outcome = compiler()
        .compile_mesh_reported(&graph, root, &options())
        .expect("compiles");
    // Volume about a local origin: it is meaningless far out.
    let measured = mesh_volume(&outcome.mesh);
    assert!((measured - (WALL - CORNER)).abs() <= 1e-6, "{measured}");
}

#[test]
fn the_mesh_compiler_refuses_a_region_with_no_solid_to_clip_by_name() {
    let placement = building();
    let refused = |g: Graph, root: NodeId, part: &str| {
        let graph = g.finish(vec![root]);
        let error = compiler()
            .compile_mesh(&graph, root, &options())
            .expect_err("refused");
        named(&error, part);
    };
    let alone_text = "boolean of a half-space that is not a clip";

    let mut g = Graph::new();
    let h1 = g.h1(placement);
    refused(g, h1, "half-space without a subject");

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let flipped = g.minus(h1, wall);
    refused(g, flipped, alone_text);

    let mut g = Graph::new();
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let alone = g.meet(h1, h2);
    refused(g, alone, alone_text);

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let region = g.meet(h1, h2);
    let united = g.boolean(wall, region, BooleanOperator::Union);
    refused(g, united, alone_text);
}

/// `Instance(wall, q) - (Instance(H1, q) & Instance({y <= 0.15}, q))` and
/// the same cut in the wall's own frame, `Instance(wall - (H1 & {y <= 0.15}), q)`.
/// The second plane lies on the wall's face, so the stand-in's face and the
/// wall's are a rounding residue apart once turned (#400, #326).
fn turned_pair(angle: f64) -> (Graph, NodeId, NodeId) {
    let q = Transform3::from_translation(Vec3::new(13.475, -15.95, 4.3))
        * Transform3::from_rotation_z(angle);
    let mut g = Graph::new();
    let wall = g.local_wall();
    let h1 = g.local_side(Point3::new(0.0, 0.0, 2.0), h1_normal());
    let face = g.local_side(Point3::new(0.0, 0.15, 0.0), -Vec3::Y);
    let local_region = g.meet(h1, face);
    let local_cut = g.minus(wall, local_region);
    let whole = g.instance(local_cut, q);

    let placed_wall = g.instance(wall, q);
    let placed_h1 = g.instance(h1, q);
    let placed_face = g.instance(face, q);
    let region = g.meet(placed_h1, placed_face);
    let world = g.minus(placed_wall, region);
    (g, world, whole)
}

#[test]
fn a_region_is_cut_in_the_subjects_frame() {
    for angle in [0.3, 0.6, 1.1, std::f64::consts::FRAC_PI_4, 2.5] {
        let (g, world, whole) = turned_pair(angle);
        let graph = g.finish(vec![world, whole]);
        let (placed, report) = compiler()
            .compile_mesh_with_deviation(&graph, world, &options())
            .unwrap_or_else(|error| panic!("{angle}: {error:?}"));
        let detail = report.contributions[0].detail;
        assert!(detail.starts_with("measured against"), "{angle}: {detail}");
        let measured = mesh_volume(&placed.mesh);
        assert!((measured - 3.6).abs() <= 1e-7 * 3.6, "{angle}: {measured}");
        let in_frame = compiler()
            .compile_mesh(&graph, whole, &options())
            .unwrap_or_else(|error| panic!("{angle}: {error:?}"));
        assert_eq!(
            placed.mesh, in_frame,
            "{angle}: cut in the wall's own frame"
        );
    }
}
