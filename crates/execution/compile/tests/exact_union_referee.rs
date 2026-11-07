//! A union the mesh cut cannot settle is decided by the exact compiler
//! (#391, ADR 0098).
//!
//! The #356 column and a corbel on one of its faces, the corbel placed by
//! the same turn composed another way: their linear parts differ in the
//! last bits, so no shared frame holds both exactly and the mesh boolean,
//! which has no tolerance, returns two bodies, a pinch, or one body with
//! sub-rounding wedge cracks. The exact compiler reads the faces as one
//! surface within the tolerance and fuses them; its union is tessellated.
//!
//! The referee is never worse than the mesh path: where it refuses, the
//! mesh result stands, or the mesh refusal with the referee's reason.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Mat3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_heal::self_intersections;
use axiolid_inspect::enclosed_volume;
use axiolid_mesh::{component_count, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::MeshBoolean;
use axiolid_mesh_compile::deviation::{
    BOOLEAN_ROUNDING_FLOOR, COINCIDENT_OPERANDS, REFEREED_UNION, REFEREED_UNION_EXACT,
    SNAPPED_OPERANDS, TOLERANCE_READINGS,
};
use axiolid_mesh_compile::{
    DeviationBound, DeviationReport, ReferenceMeshCompiler, UNMERGED_UNION_CONTACT,
};
use axiolid_mesh_compile_contract::{MeshClosure, MeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

const EPS: Scalar = 1e-6;

/// The column and one corbel.
const VOLUME: Scalar = 1.2 + 0.025;

/// The four faces: the corbel's centre in the column's frame (x, y) and its
/// profile size along x and y.
const FACES: [([Scalar; 2], [Scalar; 2]); 4] = [
    ([-0.325, 0.0], [0.25, 0.2]),
    ([0.325, 0.0], [0.25, 0.2]),
    ([0.0, -0.375], [0.2, 0.25]),
    ([0.0, 0.375], [0.2, 0.25]),
];

/// A site far from the origin, where the placement's rounding is a
/// micrometre's thousandth.
const FAR: Vec3 = Vec3::new(6.0e5, 5.6e6, 200.0);

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
        .with_chord_error(1e-3)
        .expect("a positive chord error")
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

    fn prism(&mut self, size: [Scalar; 2], depth: Scalar) -> NodeId {
        let profile = self.push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: size[0],
                y: size[1],
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn dowel(&mut self, radius: Scalar, depth: Scalar) -> NodeId {
        let profile = self.push(GeometryNode::Profile(Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        })));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn placed(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn union(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }

    fn finish(self, root: NodeId) -> (GeometryGraph, NodeId) {
        (
            self.builder.finish(vec![root]).expect("a valid graph"),
            root,
        )
    }
}

/// The turn of `theta` about Z at `site`, and the same turn computed as two
/// half turns, composed with the corbel's offset on `face`. The two linear
/// parts differ in bits for most `theta`.
fn placements(theta: Scalar, face: usize, site: Vec3) -> (Transform3, Transform3) {
    let ([cx, cy], _) = FACES[face];
    let column = Transform3::from_translation(site) * Transform3::from_rotation_z(theta);
    let half = Transform3::from_rotation_z(theta / 2.0);
    let corbel = Transform3::from_translation(site)
        * half
        * half
        * Transform3::from_translation(Vec3::new(cx, cy, 3.5));
    (column, corbel)
}

/// The column and the corbel on `face`, united, the turn composed two ways.
fn column_and_corbel(theta: Scalar, face: usize, site: Vec3) -> (GeometryGraph, NodeId) {
    let (column_at, corbel_at) = placements(theta, face, site);
    let mut g = Graph::new();
    let column = g.prism([0.4, 0.5], 6.0);
    let column = g.placed(column, column_at);
    let corbel = g.prism(FACES[face].1, 0.5);
    let corbel = g.placed(corbel, corbel_at);
    let root = g.union(column, corbel, BooleanOperator::Union);
    g.finish(root)
}

fn theta(k: usize) -> Scalar {
    0.001 + k as Scalar * 0.031_415_9
}

fn compile(graph: &GeometryGraph, root: NodeId) -> Result<(TriMesh, MeshClosure), GeomError> {
    let outcome = compiler().compile_mesh_reported(graph, root, &options())?;
    Ok((outcome.mesh, outcome.closure))
}

fn report(graph: &GeometryGraph, root: NodeId) -> DeviationReport {
    compiler()
        .compile_mesh_with_deviation(graph, root, &options())
        .expect("compiled")
        .1
}

/// The divergence volume: read through a mesh the certified measure refuses.
fn volume(mesh: &TriMesh) -> Scalar {
    let interval = enclosed_volume(mesh).expect("a measurable volume");
    0.5 * (interval.lower + interval.upper)
}

fn assert_fused(mesh: &TriMesh, closure: MeshClosure, expected: Scalar, what: &str) {
    assert_eq!(component_count(mesh), 1, "{what}: two bodies");
    assert_eq!(closure, MeshClosure::Solid, "{what}");
    assert!(
        self_intersections(mesh).is_empty(),
        "{what}: the body crosses itself"
    );
    let measured = volume(mesh);
    assert!(
        (measured - expected).abs() <= 1e-9 * expected.max(1.0),
        "{what}: volume {measured}"
    );
}

fn detail_bound(report: &DeviationReport, detail: &str) -> Option<DeviationBound> {
    report
        .contributions
        .iter()
        .find(|c| c.detail == detail)
        .map(|c| c.bound)
}

fn is_unmerged(error: &GeomError) -> bool {
    matches!(error, GeomError::Degenerate(text) if text.starts_with(UNMERGED_UNION_CONTACT))
}

fn is_pinch(error: &GeomError) -> bool {
    matches!(error, GeomError::Degenerate(text) if text.contains("touches itself"))
}

fn names_the_referee(error: &GeomError) -> bool {
    error.to_string().contains("; the exact referee refused: ")
}

/// Every differently composed turn at the origin, on every face, is one
/// solid that crosses nowhere (before #391 a third of them were refused and
/// some came back with wedge cracks). Far from the origin the exact
/// compiler refuses most turns by itself (`exact_placed_boolean.rs`); those
/// keep the mesh path's answer: a body, or its refusal naming the referee.
#[test]
fn a_union_under_differently_composed_turns_is_one_solid() {
    for face in 0..FACES.len() {
        for k in 0..50 {
            let (graph, root) = column_and_corbel(theta(k), face, Vec3::ZERO);
            let (mesh, closure) = compile(&graph, root)
                .unwrap_or_else(|e| panic!("face {face}, theta {}: {e}", theta(k)));
            assert_fused(
                &mesh,
                closure,
                VOLUME,
                &format!("face {face}, theta {}", theta(k)),
            );
        }
    }
    for face in 0..FACES.len() {
        for k in 0..50 {
            let (graph, root) = column_and_corbel(theta(k), face, FAR);
            match compile(&graph, root) {
                Ok((mesh, _)) => assert_eq!(
                    component_count(&mesh),
                    1,
                    "far, face {face}, theta {}: two bodies in contact",
                    theta(k)
                ),
                Err(error) => assert!(
                    (is_unmerged(&error) || is_pinch(&error)) && names_the_referee(&error),
                    "far, face {face}, theta {}: {error}",
                    theta(k)
                ),
            }
        }
    }
}

/// The #356 case the mesh cut answers wrongly with `Ok`: the column and the
/// corbel on face 3 at `theta(36)` come out of the mesh boolean as one body
/// that crosses itself where the faces meet. The referee replaces it.
#[test]
fn a_union_the_mesh_cut_cracks_is_decided_by_the_referee() {
    let (face, k) = (3, 36);
    let (column_at, corbel_at) = placements(theta(k), face, Vec3::ZERO);
    // The mesh boolean alone: one component, and it crosses itself.
    let mesh_of = |size: [Scalar; 2], depth: Scalar, at: Transform3| {
        let mut g = Graph::new();
        let solid = g.prism(size, depth);
        let root = g.placed(solid, at);
        let (graph, root) = g.finish(root);
        compile(&graph, root).expect("an operand").0
    };
    let column = mesh_of([0.4, 0.5], 6.0, column_at);
    let corbel = mesh_of(FACES[face].1, 0.5, corbel_at);
    let raw = BoolmeshBoolean::new()
        .boolean(&column, &corbel, BooleanOperator::Union, &options())
        .expect("the mesh boolean answers")
        .mesh;
    assert_eq!(component_count(&raw), 1);
    assert!(
        !self_intersections(&raw).is_empty(),
        "the mesh cut must be cracked for this case to mean anything"
    );
    // The compiler's answer is the exact union's.
    let (graph, root) = column_and_corbel(theta(k), face, Vec3::ZERO);
    let (mesh, closure) = compile(&graph, root).expect("fused");
    assert_fused(&mesh, closure, VOLUME, "the cracked case");
    let report = report(&graph, root);
    assert!(
        detail_bound(&report, REFEREED_UNION).is_some()
            || detail_bound(&report, REFEREED_UNION_EXACT).is_some(),
        "{:?}",
        report.contributions
    );
}

/// Operands that share a frame are cut exactly there and never refereed.
#[test]
fn a_union_in_a_shared_frame_is_not_refereed() {
    let mut g = Graph::new();
    let (column_at, _) = placements(0.7, 0, Vec3::ZERO);
    let column = g.prism([0.4, 0.5], 6.0);
    let column = g.placed(column, column_at);
    let corbel = g.prism(FACES[0].1, 0.5);
    let corbel = g.placed(
        corbel,
        column_at * Transform3::from_translation(Vec3::new(FACES[0].0[0], FACES[0].0[1], 3.5)),
    );
    let root = g.union(column, corbel, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let (mesh, closure) = compile(&graph, root).expect("fused in the shared frame");
    assert_fused(&mesh, closure, VOLUME, "shared frame");
    let report = report(&graph, root);
    assert!(detail_bound(&report, REFEREED_UNION).is_none());
    assert!(detail_bound(&report, REFEREED_UNION_EXACT).is_none());
}

/// Two boxes a gap apart, the second placed by the same turn composed
/// another way.
fn boxes(gap: Scalar, site: Vec3) -> (GeometryGraph, NodeId) {
    let a = Transform3::from_translation(site) * Transform3::from_rotation_z(0.3);
    let half = Transform3::from_rotation_z(0.15);
    let b = Transform3::from_translation(site) * half * half;
    assert_ne!(a.matrix3, b.matrix3, "the placements must differ in bits");
    let mut g = Graph::new();
    let first = g.prism([1.0, 1.0], 1.0);
    let first = g.placed(first, a);
    let second = g.prism([1.0, 1.0], 1.0);
    let second = g.placed(
        second,
        b * Transform3::from_translation(Vec3::new(1.0 + gap, 0.0, 0.0)),
    );
    let root = g.union(first, second, BooleanOperator::Union);
    g.finish(root)
}

#[test]
fn a_gap_of_twice_the_tolerance_stays_two_bodies() {
    let (graph, root) = boxes(2.0 * EPS, FAR);
    let (mesh, closure) = compile(&graph, root).expect("two bodies");
    assert_eq!(component_count(&mesh), 2);
    assert_eq!(closure, MeshClosure::Solid);
    let report = report(&graph, root);
    assert!(detail_bound(&report, REFEREED_UNION).is_none());
    assert!(detail_bound(&report, REFEREED_UNION_EXACT).is_none());
}

#[test]
fn a_referee_bound_is_widened_by_the_readings_and_the_floor_not_combined_by_worst() {
    // Two unit prisms overlapping by half a tolerance, at the origin: the
    // referee's own tessellation bound is proven, and the reading it takes
    // is added to it (#421), not the larger of the two.
    let mut g = Graph::new();
    let first = g.prism([1.0, 1.0], 1.0);
    let second = g.prism([1.0, 1.0], 1.0);
    let second = g.placed(
        second,
        Transform3::from_translation(Vec3::new(1.0 - 0.5 * EPS, 0.0, 0.0)),
    );
    let root = g.union(first, second, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let report = report(&graph, root);
    let number = |detail| detail_bound(&report, detail).and_then(|b| b.value());
    let reading = number(TOLERANCE_READINGS).expect("the reading is named");
    let floor = number(BOOLEAN_ROUNDING_FLOOR).unwrap_or(0.0);
    assert!(reading > 0.0, "{report:?}");
    let referee =
        detail_bound(&report, REFEREED_UNION).unwrap_or_else(|| panic!("referee: {report:?}"));
    assert!(
        matches!(referee, DeviationBound::Certified(v) if v >= reading + floor),
        "{referee:?} {reading:e} {floor:e}"
    );
}

#[test]
fn a_gap_within_the_tolerance_is_read_as_contact_and_says_so() {
    let (graph, root) = boxes(0.5 * EPS, FAR);
    let (mesh, closure) = compile(&graph, root).expect("fused");
    assert_eq!(component_count(&mesh), 1);
    assert_eq!(closure, MeshClosure::Solid);
    assert!((volume(&mesh) - 2.0).abs() < 1e-5, "{}", volume(&mesh));
    let report = report(&graph, root);
    // The exact union read the faces as one surface: the move is named, and
    // the union is the exact union of operands within the tolerance.
    let moved = detail_bound(&report, COINCIDENT_OPERANDS).expect("the reading is named");
    assert!(
        matches!(moved, DeviationBound::Certified(m) if m > 0.0 && m <= EPS),
        "{moved:?}"
    );
    // The largest reading of every kind and the rounding floor are named
    // too (#421), the reading covering the coincident supports' move.
    let readings = detail_bound(&report, TOLERANCE_READINGS).expect("every reading is named");
    assert!(
        matches!((moved, readings),
            (DeviationBound::Certified(m), DeviationBound::Certified(r)) if r >= m && r <= EPS),
        "{moved:?} {readings:?}"
    );
    assert!(detail_bound(&report, BOOLEAN_ROUNDING_FLOOR).is_some());
    assert!(detail_bound(&report, REFEREED_UNION).is_some());
    assert!(detail_bound(&report, REFEREED_UNION_EXACT).is_none());
    assert!(detail_bound(&report, SNAPPED_OPERANDS).is_none());
}

/// A sheared placement is no similarity, and the exact compiler refuses it:
/// the mesh cut's refusal stands and names the referee's.
#[test]
fn a_sheared_union_keeps_its_refusal_and_names_the_referee() {
    let shear = axiolid_core::Mat3::from_cols(Vec3::X, Vec3::new(0.3, 1.0, 0.0), Vec3::Z);
    let site = Vec3::new(1.0, 2.0, 0.0);
    let a = Transform3::from_mat3_translation(shear, site);
    let b = Transform3::from_mat3_translation(shear, site);
    let mut g = Graph::new();
    let first = g.prism([1.0, 1.0], 1.0);
    let first = g.placed(first, a);
    let second = g.prism([1.0, 1.0], 1.0);
    // One box-edge along the sheared x axis plus half the tolerance.
    let second = g.placed(
        second,
        b * Transform3::from_translation(Vec3::new(1.0 + 0.5 * EPS, 0.0, 0.0)),
    );
    let root = g.union(first, second, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let error = compile(&graph, root).expect_err("the contact is not merged");
    assert!(is_unmerged(&error), "{error}");
    assert!(names_the_referee(&error), "{error}");
}

/// Boxes that meet along one edge only: the pinch refusal stays, and the
/// referee's reason is appended.
#[test]
fn an_edge_contact_stays_a_refusal() {
    let mut g = Graph::new();
    let first = g.prism([1.0, 1.0], 1.0);
    let first = g.placed(first, Transform3::IDENTITY);
    let second = g.prism([1.0, 1.0], 1.0);
    let second = g.placed(
        second,
        Transform3::from_translation(Vec3::new(1.0, 1.0, 0.0)),
    );
    let root = g.union(first, second, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let error = compile(&graph, root).expect_err("an edge contact is refused");
    assert!(is_pinch(&error), "{error}");
    assert!(
        error
            .to_string()
            .ends_with("the exact referee refused: the tessellated exact union touches itself"),
        "{error}"
    );
}

/// A dowel standing on a slab, the turn composed two ways: the union is
/// certified within the chord budget of the closed-form volume.
#[test]
fn a_dowel_union_is_certified_against_its_closed_form() {
    let (radius, depth) = (0.3, 0.8);
    // An oblique turn, so no face of either is axis-aligned and the dowel's
    // foot meets the slab's top in a plane no frame holds exactly.
    let axis = Vec3::new(1.0, 2.0, 3.0).normalize();
    let site = Transform3::from_translation(Vec3::new(2.0, 3.0, 0.0));
    let at = site * Transform3::from_axis_angle(axis, 0.9);
    let half = Transform3::from_axis_angle(axis, 0.45);
    let dowel_at = site * half * half * Transform3::from_translation(Vec3::new(0.0, 0.0, 1.0));
    let mut g = Graph::new();
    let slab = g.prism([4.0, 4.0], 1.0);
    let slab = g.placed(slab, at);
    let dowel = g.dowel(radius, depth);
    let dowel = g.placed(dowel, dowel_at);
    let root = g.union(slab, dowel, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let (mesh, closure) = compile(&graph, root).expect("fused");
    assert_eq!(component_count(&mesh), 1);
    assert_eq!(closure, MeshClosure::Solid);
    let report = report(&graph, root);
    let bound = report.bound.expect("a certified bound");
    assert!(bound <= 1e-3, "{bound}");
    // The union's own bound is the tessellation's: positive for a round
    // dowel, and what the report's overall bound is made of.
    let own = detail_bound(&report, REFEREED_UNION)
        .or_else(|| detail_bound(&report, REFEREED_UNION_EXACT))
        .expect("a refereed union");
    assert!(
        matches!(own, DeviationBound::Certified(b) if b > 0.0 && b <= bound),
        "{own:?}"
    );
    assert!(report.meets_requested(), "{:?}", report.contributions);
    let expected = 16.0 + std::f64::consts::PI * radius * radius * depth;
    let wall = 2.0 * std::f64::consts::PI * radius * depth;
    assert!(
        (volume(&mesh) - expected).abs() <= bound * wall + 1e-9,
        "{} against {expected}",
        volume(&mesh)
    );
}

/// An authored mesh operand cannot be refereed (the exact compiler refuses
/// it, ADR 0080): an unsettled but valid mesh result is returned as before.
#[test]
fn an_authored_mesh_operand_returns_the_mesh_result() {
    // A unit cube as an authored mesh, standing on a slab turned another way.
    let cube = {
        let corners = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let triangles: [u32; 36] = [
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ];
        TriMesh::new(
            corners
                .iter()
                .map(|c| axiolid_core::Point3::new(c[0], c[1], c[2]))
                .collect(),
            triangles.to_vec(),
        )
    };
    let mut found = 0;
    for k in 0..40 {
        let theta = theta(k);
        let turn = Transform3::from_rotation_z(theta);
        let half = Transform3::from_rotation_z(theta / 2.0);
        let mut g = Graph::new();
        let slab = g.prism([4.0, 4.0], 1.0);
        let slab = g.placed(slab, turn);
        let authored = g.push(GeometryNode::TriMesh(cube.clone()));
        let authored = g.placed(
            authored,
            half * half * Transform3::from_translation(Vec3::new(-0.5, -0.5, 1.0)),
        );
        let root = g.union(slab, authored, BooleanOperator::Union);
        let (graph, root) = g.finish(root);
        match compile(&graph, root) {
            Ok((mesh, closure)) => {
                if component_count(&mesh) == 1 && closure == MeshClosure::Solid {
                    found += 1;
                    let report = report(&graph, root);
                    assert!(
                        detail_bound(&report, REFEREED_UNION).is_none()
                            && detail_bound(&report, REFEREED_UNION_EXACT).is_none(),
                        "theta {theta}: a mesh operand was refereed"
                    );
                }
            }
            Err(error) => assert!(
                (is_unmerged(&error) || is_pinch(&error)) && names_the_referee(&error),
                "theta {theta}: {error}"
            ),
        }
    }
    assert!(found > 0, "no turn fused as one body from the mesh cut");
}

/// The arm of a refused referee that keeps a valid mesh result: an authored
/// cube sunk `1e-7` into a slab tilted by `1e-8` (within the tolerance, and
/// no axis-aligned face for the snap to land on) leaves faces that nearly
/// share a plane, so the cut is unsettled; the exact compiler refuses the
/// mesh operand, and the mesh union, one valid body, stands. Pinned by
/// search: after #418 the exact compiler decides the far turns that once
/// reached this arm.
#[test]
fn a_refused_referee_keeps_a_valid_mesh_result() {
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let triangles: [u32; 36] = [
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    let cube = TriMesh::new(
        corners
            .iter()
            .map(|c| axiolid_core::Point3::new(c[0], c[1], c[2]))
            .collect(),
        triangles.to_vec(),
    );
    let mut g = Graph::new();
    let slab = g.prism([4.0, 4.0], 1.0);
    let slab = g.placed(slab, Transform3::from_rotation_x(1e-8));
    let authored = g.push(GeometryNode::TriMesh(cube));
    let authored = g.placed(
        authored,
        Transform3::from_translation(Vec3::new(-0.5, -0.5, 1.0 - 1e-7)),
    );
    let root = g.union(slab, authored, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let (mesh, closure) = compile(&graph, root).expect("the mesh union stands");
    assert_eq!(closure, MeshClosure::Solid);
    assert_eq!(component_count(&mesh), 1);
    let report = report(&graph, root);
    assert!(
        detail_bound(&report, REFEREED_UNION).is_none()
            && detail_bound(&report, REFEREED_UNION_EXACT).is_none(),
        "a mesh operand was refereed"
    );
}

/// Operands the snap lands on each other are settled, whatever the given
/// pair looks like: the pair actually cut decides. The two linear parts
/// differ (an exact quarter turn and none) and are exactly axis-aligned, so
/// the world-frame snap closes the gap.
#[test]
fn a_union_the_snap_settles_is_not_refereed() {
    let quarter = Mat3::from_cols(Vec3::Y, Vec3::NEG_X, Vec3::Z);
    let gap = 0.5 * EPS;
    let mut g = Graph::new();
    let first = g.prism([1.0, 1.0], 1.0);
    let first = g.placed(first, Transform3::IDENTITY);
    let second = g.prism([1.0, 1.0], 1.0);
    let second = g.placed(
        second,
        Transform3::from_mat3_translation(quarter, Vec3::new(0.0, 1.0 + gap, 0.0)),
    );
    let root = g.union(first, second, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let (mesh, closure) = compile(&graph, root).expect("fused");
    assert_eq!(component_count(&mesh), 1);
    assert_eq!(closure, MeshClosure::Solid);
    // The snap's: the tool's near face lands on the host's, gaining the gap.
    assert!(
        (volume(&mesh) - (2.0 + gap)).abs() < 1e-12,
        "{}",
        volume(&mesh)
    );
    let report = report(&graph, root);
    assert!(detail_bound(&report, SNAPPED_OPERANDS).is_some());
    assert!(detail_bound(&report, REFEREED_UNION).is_none());
    assert!(detail_bound(&report, REFEREED_UNION_EXACT).is_none());
}

/// The framed sources carry the differently composed turns themselves
/// (each wrapped in one shared translation), so the shared frame holds the
/// same near-contact: it yields, and the referee replaces the cracked body.
#[test]
fn a_shared_frame_that_holds_unsettled_contact_yields_to_the_referee() {
    let (face, k) = (3, 47);
    let wrap = Transform3::from_translation(Vec3::new(10.0, 0.0, 0.0));
    let mut g = Graph::new();
    let column = g.prism([0.4, 0.5], 6.0);
    let column = g.placed(column, Transform3::from_rotation_z(theta(k)));
    let half = Transform3::from_rotation_z(theta(k) / 2.0);
    let corbel = g.prism(FACES[face].1, 0.5);
    let corbel = g.placed(
        corbel,
        half * half
            * Transform3::from_translation(Vec3::new(FACES[face].0[0], FACES[face].0[1], 3.5)),
    );
    let column = g.placed(column, wrap);
    let corbel = g.placed(corbel, wrap);
    let root = g.union(column, corbel, BooleanOperator::Union);
    let (graph, root) = g.finish(root);
    let (mesh, closure) = compile(&graph, root).expect("fused");
    assert_fused(&mesh, closure, VOLUME, "wrapped sources");
    let report = report(&graph, root);
    assert!(
        detail_bound(&report, REFEREED_UNION).is_some()
            || detail_bound(&report, REFEREED_UNION_EXACT).is_some(),
        "{:?}",
        report.contributions
    );
}
