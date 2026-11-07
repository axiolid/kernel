//! Empty operands of a boolean, in both compilers (#326).
//!
//! A half-space can leave nothing of a wall: a plane above it keeps nothing
//! of an intersection, a plane below it removes all of a difference. Such an
//! empty intermediate is then an operand of another boolean, decided by the
//! identities of `crate::empty` (`S - 0 = S`, `0 - T = 0`, `S & 0 = 0`,
//! `S | 0 = S`), not refused. The exact compiler cannot hold an empty
//! `ExactBRep`, so an empty root stays a named `Degenerate` refusal; the
//! mesh compiler's root is an empty mesh.
//!
//! Each case is built once and run through both compilers. The wall is
//! `6 x 0.3 x 3` under a rigid placement; `Hmiss` is the side of
//! `z = 3.5` above it, `Hover` the side of `z = -0.5` above, which holds the
//! whole wall.

use axiolid_brep::ExactBRep;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Plane3, Point3, Tolerance, Transform3, Vec3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::EMPTY_RESULT;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_mesh_compile_contract::{MeshClosure, MeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

const L: f64 = 6.0;
const T: f64 = 0.3;
const H: f64 = 3.0;
const WALL: f64 = L * T * H;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
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

    fn extrusion(&mut self, profile: Profile, depth: f64) -> NodeId {
        let profile = self.push(GeometryNode::Profile(profile));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
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

    fn join(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.boolean(left, right, BooleanOperator::Union)
    }

    /// The wall in its own frame.
    fn local_wall(&mut self) -> NodeId {
        self.extrusion(rect(L, T), H)
    }

    fn wall(&mut self, placement: Transform3) -> NodeId {
        let wall = self.local_wall();
        self.instance(wall, placement)
    }

    /// The half-space through `origin` with `normal`, both given in the
    /// wall's frame, in world coordinates under `placement`.
    fn world_half_space(&mut self, origin: Point3, normal: Vec3, placement: Transform3) -> NodeId {
        self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: placement.transform_point3(origin),
                normal: placement.transform_vector3(normal),
            },
            agreement: true,
        }))
    }

    /// The same, in the wall's frame.
    fn local_half_space(&mut self, origin: Point3, normal: Vec3) -> NodeId {
        self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 { origin, normal },
            agreement: true,
        }))
    }

    /// The side of `z = 3.5` above the wall: it misses it.
    fn miss(&mut self, placement: Transform3) -> NodeId {
        self.world_half_space(Point3::new(0.0, 0.0, H + 0.5), Vec3::Z, placement)
    }

    /// The side of `z = -0.5` above: it holds the whole wall.
    fn over(&mut self, placement: Transform3) -> NodeId {
        self.world_half_space(Point3::new(0.0, 0.0, -0.5), Vec3::Z, placement)
    }

    fn finish(self, roots: Vec<NodeId>) -> GeometryGraph {
        self.builder.finish(roots).expect("a valid graph")
    }
}

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// What a case compiles to.
#[derive(Clone, Copy)]
enum Want {
    /// The wall under this placement, as the operand itself.
    Wall(Transform3),
    /// Nothing; the exact compiler's refusal of it as a root says this.
    Nothing(&'static str),
}

struct Case {
    graph: Graph,
    root: NodeId,
    want: Want,
}

const CLIP_EMPTIED: &str = "exact half-space clip leaves nothing of the subject";

fn apart() -> Transform3 {
    building() * Transform3::from_translation(Vec3::new(0.0, 2.0, 0.0))
}

/// `wall - (wall & Hmiss)`: the tool is the clip's empty intersection.
fn difference_by_empty_clip() -> Case {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let miss = g.miss(placement);
    let nothing = g.meet(wall, miss);
    let root = g.minus(wall, nothing);
    Case {
        graph: g,
        root,
        want: Want::Wall(placement),
    }
}

/// `(wall - Hover) | wall2`: an empty subject united with a wall.
fn union_with_empty_subject() -> Case {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let over = g.over(placement);
    let nothing = g.minus(wall, over);
    let wall2 = g.wall(apart());
    let root = g.join(nothing, wall2);
    Case {
        graph: g,
        root,
        want: Want::Wall(apart()),
    }
}

/// `wall - (wall - Hover)`.
fn difference_by_emptied_wall() -> Case {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let over = g.over(placement);
    let nothing = g.minus(wall, over);
    let root = g.minus(wall, nothing);
    Case {
        graph: g,
        root,
        want: Want::Wall(placement),
    }
}

/// `(wall - Hover) - wall2`: nothing, less a wall.
fn empty_subject_less_a_wall() -> Case {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let over = g.over(placement);
    let nothing = g.minus(wall, over);
    let wall2 = g.wall(apart());
    let root = g.minus(nothing, wall2);
    Case {
        graph: g,
        root,
        want: Want::Nothing(CLIP_EMPTIED),
    }
}

/// `wall - (wall & wall_far)`: the general boolean's empty result.
fn difference_by_disjoint_intersection() -> Case {
    let placement = building();
    let far = Transform3::from_translation(Vec3::new(100.0, 0.0, 0.0)) * placement;
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let away = g.wall(far);
    let nothing = g.meet(wall, away);
    let root = g.minus(wall, nothing);
    Case {
        graph: g,
        root,
        want: Want::Wall(placement),
    }
}

/// `wall & wall_far` on its own.
fn disjoint_intersection() -> Case {
    let placement = building();
    let far = Transform3::from_translation(Vec3::new(100.0, 0.0, 0.0)) * placement;
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let away = g.wall(far);
    let root = g.meet(wall, away);
    Case {
        graph: g,
        root,
        want: Want::Nothing("exact boolean intersection leaves nothing"),
    }
}

/// A column of radius 0.5 and a plane 0.6 from its axis, facing out: the
/// column's box reaches across the plane, the column does not, so the
/// clip's boolean is what finds nothing.
fn column_by_a_plane_it_misses(g: &mut Graph) -> NodeId {
    let placement = general();
    let profile = Profile::Circle(CircleProfile {
        radius: 0.5,
        thickness: None,
    });
    let column = g.extrusion(profile, H);
    let column = g.instance(column, placement);
    let d = 0.6 / 2.0_f64.sqrt();
    let plane = g.world_half_space(
        Point3::new(d, d, 0.0),
        Vec3::new(1.0, 1.0, 0.0).normalize(),
        placement,
    );
    g.meet(column, plane)
}

fn the_column_clip_alone() -> Case {
    let mut g = Graph::new();
    let root = column_by_a_plane_it_misses(&mut g);
    Case {
        graph: g,
        root,
        want: Want::Nothing(CLIP_EMPTIED),
    }
}

fn difference_by_the_column_clip() -> Case {
    let placement = building();
    let mut g = Graph::new();
    let nothing = column_by_a_plane_it_misses(&mut g);
    let wall = g.wall(placement);
    let root = g.minus(wall, nothing);
    Case {
        graph: g,
        root,
        want: Want::Wall(placement),
    }
}

/// `wall - Instance(wall_local & Hmiss_local, building)`.
fn difference_by_an_instance_of_nothing() -> Case {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let local = g.local_wall();
    let miss = g.local_half_space(Point3::new(0.0, 0.0, H + 0.5), Vec3::Z);
    let nothing = g.meet(local, miss);
    let placed = g.instance(nothing, placement);
    let root = g.minus(wall, placed);
    Case {
        graph: g,
        root,
        want: Want::Wall(placement),
    }
}

/// `wall & Hmiss` at the root.
fn empty_root() -> Case {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let miss = g.miss(placement);
    let root = g.meet(wall, miss);
    Case {
        graph: g,
        root,
        want: Want::Nothing(CLIP_EMPTIED),
    }
}

fn every_case() -> Vec<(&'static str, Case)> {
    vec![
        ("difference by an empty clip", difference_by_empty_clip()),
        ("union with an empty subject", union_with_empty_subject()),
        (
            "difference by an emptied wall",
            difference_by_emptied_wall(),
        ),
        ("an empty subject less a wall", empty_subject_less_a_wall()),
        (
            "difference by a disjoint intersection",
            difference_by_disjoint_intersection(),
        ),
        ("a disjoint intersection", disjoint_intersection()),
        ("the column clip alone", the_column_clip_alone()),
        (
            "difference by the column clip",
            difference_by_the_column_clip(),
        ),
        (
            "difference by an instance of nothing",
            difference_by_an_instance_of_nothing(),
        ),
        ("an empty root", empty_root()),
    ]
}

fn wall_only(placement: Transform3) -> (GeometryGraph, NodeId) {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    (g.finish(vec![wall]), wall)
}

fn exact(graph: &GeometryGraph, root: NodeId) -> Result<ExactBRep, GeomError> {
    ReferenceExactCompiler::new().compile_exact(graph, root, &options())
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

#[test]
fn the_exact_compiler_absorbs_an_empty_operand_and_refuses_an_empty_root() {
    for (name, case) in every_case() {
        let Case { graph, root, want } = case;
        let graph = graph.finish(vec![root]);
        match want {
            Want::Wall(placement) => {
                let (other, wall) = wall_only(placement);
                let wanted = exact(&other, wall).expect("the operand compiles");
                let solid = exact(&graph, root)
                    .unwrap_or_else(|error| panic!("{name}: the identity absorbs it: {error:?}"));
                assert!(
                    (volume(&solid) - volume(&wanted)).abs() <= 1e-9 * volume(&wanted),
                    "{name}: volume {}, expected {}",
                    volume(&solid),
                    volume(&wanted)
                );
                assert_eq!(solid, wanted, "{name}: the operand itself, not a copy");
            }
            Want::Nothing(why) => match exact(&graph, root) {
                Err(GeomError::Degenerate(text)) => {
                    assert!(text.contains(why), "{name}: {text:?} must say `{why}`");
                }
                other => panic!("{name}: expected Degenerate({why}), got {other:?}"),
            },
        }
    }
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn mesh_volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

#[test]
fn the_mesh_compiler_absorbs_an_empty_operand_and_returns_an_empty_root() {
    for (name, case) in every_case() {
        let Case { graph, root, want } = case;
        let graph = graph.finish(vec![root]);
        let outcome = compiler()
            .compile_mesh_reported(&graph, root, &options())
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        assert_eq!(outcome.closure, MeshClosure::Solid, "{name}");
        match want {
            Want::Wall(placement) => {
                let (other, wall) = wall_only(placement);
                let wanted = compiler()
                    .compile_mesh(&other, wall, &options())
                    .expect("the operand compiles");
                assert_eq!(outcome.mesh, wanted, "{name}: the operand's mesh itself");
                assert!(
                    (mesh_volume(&outcome.mesh) - WALL).abs() <= 1e-9 * WALL,
                    "{name}"
                );
            }
            Want::Nothing(_) => {
                assert!(outcome.mesh.indices.is_empty(), "{name}: no triangles");
            }
        }
    }
}

#[test]
fn an_authored_empty_mesh_is_still_a_refused_tool() {
    let mut g = Graph::new();
    let wall = g.wall(building());
    let empty = g.push(GeometryNode::TriMesh(TriMesh::default()));
    let root = g.minus(wall, empty);
    let graph = g.finish(vec![root]);
    let error = compiler()
        .compile_mesh(&graph, root, &options())
        .expect_err("a surface model encloses no volume");
    assert!(
        matches!(&error, GeomError::InvalidInput(text) if text.starts_with("boolean tool is a surface model")),
        "{error:?}"
    );
}

/// The deviation of `case` compiled to a mesh, within the request.
fn deviation_of(case: Case) -> axiolid_mesh_compile::DeviationReport {
    let Case { graph, root, .. } = case;
    let graph = graph.finish(vec![root]);
    compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .expect("the mesh compiles")
        .1
}

#[test]
fn an_empty_boolean_with_an_empty_mesh_is_certified_zero() {
    for (name, case) in [
        ("an empty root", empty_root()),
        ("the column clip alone", the_column_clip_alone()),
        ("an empty subject less a wall", empty_subject_less_a_wall()),
    ] {
        let report = deviation_of(case);
        assert_eq!(report.contributions.len(), 1, "{name}");
        let only = &report.contributions[0];
        assert_eq!(only.path, DeviationPath::Boolean, "{name}");
        assert_eq!(only.detail, EMPTY_RESULT, "{name}");
        assert_eq!(only.bound, DeviationBound::Certified(0.0), "{name}");
        assert_eq!(report.bound, Some(0.0), "{name}");
        assert!(report.meets_requested(), "{name}");
    }
}

#[test]
fn a_difference_by_an_empty_tool_is_certified_within_the_request() {
    let report = deviation_of(difference_by_empty_clip());
    assert!(report.bound.is_some(), "{report:?}");
    assert!(report.meets_requested(), "{report:?}");
}

/// `wall - (A & B)`, `A` and `B` blocks inside the wall that overlap by
/// `gap`: their intersection is a cavity `gap` wide.
fn wall_less_overlapping_blocks(gap: f64) -> (Graph, NodeId) {
    let mut g = Graph::new();
    let wall = g.wall(Transform3::IDENTITY);
    let a = g.extrusion(rect(1.0, T), 1.0);
    let a = g.instance(a, Transform3::from_translation(Vec3::new(0.0, 0.0, 1.0)));
    let b = g.extrusion(rect(1.0, T), 1.0);
    let b = g.instance(
        b,
        Transform3::from_translation(Vec3::new(1.0 - gap, 0.0, 1.0)),
    );
    let cavity = g.meet(a, b);
    let root = g.minus(wall, cavity);
    (g, root)
}

/// `wall - ({z >= 2 + 0.2 x} & {z <= 2 + 0.2 x + gap})`: a slab `gap` thick.
fn wall_less_a_slab(gap: f64) -> (Graph, NodeId) {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let above = g.world_half_space(
        Point3::new(0.0, 0.0, 2.0),
        Vec3::new(-0.2, 0.0, 1.0),
        placement,
    );
    let below = g.world_half_space(
        Point3::new(0.0, 0.0, 2.0 + gap),
        Vec3::new(0.2, 0.0, -1.0),
        placement,
    );
    let slab = g.meet(above, below);
    let root = g.minus(wall, slab);
    (g, root)
}

/// `wall - (wall & {z <= gap})`: the clip leaves a sliver `gap` thick.
fn wall_less_a_sliver_clip(gap: f64) -> (Graph, NodeId) {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let low = g.world_half_space(Point3::new(0.0, 0.0, gap), -Vec3::Z, placement);
    let sliver = g.meet(wall, low);
    let root = g.minus(wall, sliver);
    (g, root)
}

/// An emptiness the general boolean reads within the tolerance (operands
/// that overlap by less than it) is no identity: passing the wall through
/// would drop the cavity and call the result exact. It is refused as
/// `Degenerate`, as it was before empty bodies; above the tolerance the
/// cavity is cut.
#[test]
fn an_emptiness_read_within_the_tolerance_is_refused_not_absorbed() {
    // The tolerance is 1e-6.
    for (name, build) in [
        (
            "overlapping blocks",
            wall_less_overlapping_blocks as fn(f64) -> (Graph, NodeId),
        ),
        ("a slab between planes", wall_less_a_slab),
        ("a sliver clip", wall_less_a_sliver_clip),
    ] {
        let (g, root) = build(1e-7);
        let graph = g.finish(vec![root]);
        match exact(&graph, root) {
            Err(GeomError::Degenerate(_)) => {}
            other => panic!("{name}: a feature read away must be refused, got {other:?}"),
        }

        let (g, root) = build(1e-3);
        let graph = g.finish(vec![root]);
        let solid = exact(&graph, root).unwrap_or_else(|error| panic!("{name}: {error:?}"));
        let cavity = if name == "overlapping blocks" {
            T * 1e-3
        } else {
            L * T * 1e-3
        };
        assert!(
            (volume(&solid) - (WALL - cavity)).abs() <= 1e-7,
            "{name}: volume {}, expected {}",
            volume(&solid),
            WALL - cavity
        );
    }
}
