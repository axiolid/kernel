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
//! The wall is `6 x 0.3 x 3` under a rigid placement; `Hmiss` is the side of
//! `z = 3.5` above it, `Hover` the side of `z = -0.5` above, which holds the
//! whole wall.

use axiolid_brep::ExactBRep;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Plane3, Point3, Tolerance, Transform3, Vec3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh_compile::ReferenceExactCompiler;
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

fn exact(graph: Graph, root: NodeId) -> Result<ExactBRep, GeomError> {
    let graph = graph.finish(vec![root]);
    ReferenceExactCompiler::new().compile_exact(&graph, root, &options())
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

/// The exact result of `body`, which must be the wall's volume, and equal
/// to the exact result of `same`, built in a graph of its own.
fn is_the_wall(graph: Graph, body: NodeId, same: impl FnOnce(&mut Graph) -> NodeId) -> ExactBRep {
    let mut other = Graph::new();
    let wanted = same(&mut other);
    let wanted = exact(other, wanted).expect("the operand compiles");
    let solid = exact(graph, body).expect("the identity absorbs the empty operand");
    assert!(
        (volume(&solid) - volume(&wanted)).abs() <= 1e-9 * volume(&wanted),
        "volume {}, expected {}",
        volume(&solid),
        volume(&wanted)
    );
    assert_eq!(solid, wanted, "the operand itself, not a rebuilt copy");
    solid
}

fn degenerate(graph: Graph, root: NodeId, part: &str) {
    match exact(graph, root) {
        Err(GeomError::Degenerate(why)) => {
            assert!(why.contains(part), "{why:?} must say `{part}`");
        }
        other => panic!("expected Degenerate({part}), got {other:?}"),
    }
}

const CLIP_EMPTIED: &str = "exact half-space clip leaves nothing of the subject";

#[test]
fn a_difference_by_an_empty_intersection_is_the_wall() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let miss = g.miss(placement);
    let nothing = g.meet(wall, miss);
    let body = g.minus(wall, nothing);
    let solid = is_the_wall(g, body, |g| g.wall(placement));
    assert!((volume(&solid) - WALL).abs() <= 1e-9 * WALL);
}

#[test]
fn an_empty_difference_is_absorbed_by_a_union_and_by_a_difference() {
    let placement = building();
    let apart = placement * Transform3::from_translation(Vec3::new(0.0, 2.0, 0.0));

    // (wall - Hover) | wall2 is wall2.
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let over = g.over(placement);
    let nothing = g.minus(wall, over);
    let wall2 = g.wall(apart);
    let body = g.join(nothing, wall2);
    is_the_wall(g, body, |g| g.wall(apart));

    // wall - (wall - Hover) is the wall.
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let over = g.over(placement);
    let nothing = g.minus(wall, over);
    let body = g.minus(wall, nothing);
    is_the_wall(g, body, |g| g.wall(placement));

    // (wall - Hover) - wall2 is empty, and so a refused root.
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let over = g.over(placement);
    let nothing = g.minus(wall, over);
    let wall2 = g.wall(apart);
    let body = g.minus(nothing, wall2);
    degenerate(g, body, CLIP_EMPTIED);
}

#[test]
fn a_general_boolean_that_leaves_nothing_is_an_empty_tool() {
    // The wall's intersection with a wall far away: the general boolean
    // reports `EmptyResult`, which is an empty operand, not a refusal.
    let placement = building();
    let far = Transform3::from_translation(Vec3::new(100.0, 0.0, 0.0)) * placement;
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let away = g.wall(far);
    let nothing = g.meet(wall, away);
    let body = g.minus(wall, nothing);
    is_the_wall(g, body, |g| g.wall(placement));

    // On its own the intersection is a refused root, as before.
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let away = g.wall(far);
    let nothing = g.meet(wall, away);
    degenerate(g, nothing, "exact boolean intersection leaves nothing");
}

#[test]
fn a_plane_that_crosses_the_envelope_but_not_the_column_leaves_nothing() {
    // A column of radius 0.5 and a plane 0.6 from its axis, facing out: the
    // column's box reaches across the plane, the column does not, so the
    // clip's boolean is what finds nothing.
    let placement = general();
    let r = 0.5;
    let build = |g: &mut Graph| {
        let profile = Profile::Circle(CircleProfile {
            radius: r,
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
    };
    let mut g = Graph::new();
    let nothing = build(&mut g);
    degenerate(g, nothing, CLIP_EMPTIED);

    let mut g = Graph::new();
    let nothing = build(&mut g);
    let wall = g.wall(building());
    let body = g.minus(wall, nothing);
    is_the_wall(g, body, |g| g.wall(building()));
}

#[test]
fn an_instance_of_an_empty_boolean_is_empty() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let local = g.local_wall();
    let miss = g.local_half_space(Point3::new(0.0, 0.0, H + 0.5), Vec3::Z);
    let nothing = g.meet(local, miss);
    let placed = g.instance(nothing, placement);
    let body = g.minus(wall, placed);
    is_the_wall(g, body, |g| g.wall(placement));
}

#[test]
fn an_empty_root_is_a_named_refusal() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let miss = g.miss(placement);
    let nothing = g.meet(wall, miss);
    degenerate(g, nothing, CLIP_EMPTIED);
}
