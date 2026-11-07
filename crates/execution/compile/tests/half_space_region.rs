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
use axiolid_mesh_compile::ReferenceExactCompiler;
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

#[test]
fn a_wall_less_two_half_spaces_is_exact() {
    for placement in [building(), site()] {
        let mut g = Graph::new();
        let wall = g.wall(placement);
        let (h1, h2) = (g.h1(placement), g.h2(placement));
        let both = g.meet(h1, h2);
        let body = g.minus(wall, both);
        // At the site the volume is read about a far origin.
        let (tolerance, expected) = (1e-9, WALL - CORNER);
        let (solid, exact) = if placement == site() {
            let graph = g.finish(vec![body]);
            let (solid, exact) = compile(&graph, body).expect("compiles exactly");
            assert!(
                (volume(&solid) - expected).abs() <= 1e-7,
                "{}",
                volume(&solid)
            );
            (solid, exact)
        } else {
            check(g, body, expected, tolerance)
        };
        assert!(exact, "the report reads exact");
        assert_eq!(solid.topology().shells().len(), 1);
    }
}

#[test]
fn a_wall_meeting_two_half_spaces_keeps_their_corner() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let both = g.meet(h1, h2);
    let body = g.meet(wall, both);
    near(g, body, CORNER);
}

#[test]
fn a_region_under_an_instance_is_placed_with_it() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.local_side(Point3::new(0.0, 0.0, 2.0), h1_normal());
    let h2 = g.local_side(Point3::new(1.0, 0.0, 0.0), -Vec3::X);
    let both = g.meet(h1, h2);
    let placed = g.instance(both, placement);
    let body = g.minus(wall, placed);
    near(g, body, WALL - CORNER);
}

#[test]
fn a_region_of_three_nests_and_keeps_its_leaves_in_order() {
    // (H1 & H2) & {y <= 0}: the half of the corner.
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let (h1, h2) = (g.h1(placement), g.h2(placement));
    let corner = g.meet(h1, h2);
    let back = g.side(Point3::ZERO, -Vec3::Y, placement);
    let region = g.meet(corner, back);
    let body = g.minus(wall, region);
    near(g, body, WALL - CORNER / 2.0);
}

#[test]
fn a_band_cuts_the_wall_in_two() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let above = g.side(Point3::new(0.0, 0.0, 1.0), Vec3::Z, placement);
    let below = g.side(Point3::new(0.0, 0.0, 2.0), -Vec3::Z, placement);
    let band = g.meet(above, below);
    let body = g.minus(wall, band);
    let (solid, exact) = near(g, body, WALL - L * T);
    assert!(exact);
    assert_eq!(solid.topology().shells().len(), 2);
}

#[test]
fn a_leaf_missing_the_wall_leaves_the_region_empty() {
    let placement = building();
    // x >= 10 lies beyond the wall's end.
    let build = |g: &mut Graph, operator: BooleanOperator| {
        let wall = g.wall(placement);
        let h1 = g.h1(placement);
        let far = g.side(Point3::new(10.0, 0.0, 0.0), Vec3::X, placement);
        let region = g.meet(h1, far);
        g.boolean(wall, region, operator)
    };
    let mut g = Graph::new();
    let body = build(&mut g, BooleanOperator::Difference);
    near(g, body, WALL);

    let mut g = Graph::new();
    let body = build(&mut g, BooleanOperator::Intersection);
    let error = refusal(g, body);
    assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
}

#[test]
fn a_bounded_half_space_is_a_leaf_in_either_place() {
    // H1 within x in [-1, 1]: the removed height 1 - 0.2 x over a width 2.
    let placement = building();
    for bounded_first in [true, false] {
        let mut g = Graph::new();
        let wall = g.wall(placement);
        let h1 = g.h1(placement);
        let square = g.square(placement);
        let region = if bounded_first {
            g.meet(square, h1)
        } else {
            g.meet(h1, square)
        };
        let body = g.minus(wall, region);
        near(g, body, WALL - 2.0 * T);
    }
}

#[test]
fn a_leaf_holding_the_wall_adds_nothing() {
    // {z >= -5} holds the whole wall: H1 alone takes 0.3 * 6 = 1.8.
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let floor = g.side(Point3::new(0.0, 0.0, -5.0), Vec3::Z, placement);
    let h1 = g.h1(placement);
    let region = g.meet(floor, h1);
    let body = g.minus(wall, region);
    near(g, body, WALL - 1.8);
}

#[test]
fn two_stand_ins_with_nothing_in_common_leave_the_wall_whole() {
    // z >= 2 + x/2 and z <= 1 + x/2 are parallel and a unit apart: each
    // crosses the wall, and their intersection is empty.
    let placement = building();
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
    let body = g.minus(wall, region);
    near(g, body, WALL);
}

#[test]
fn an_intersection_commutes_with_a_region() {
    let placement = building();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let h1 = g.h1(placement);
    let swapped = g.meet(h1, wall);
    let plain = g.meet(wall, h1);
    let region = {
        let h2 = g.h2(placement);
        g.meet(h1, h2)
    };
    let regional = g.meet(region, wall);
    let graph = g.finish(vec![swapped, plain, regional]);
    let (swapped_body, _) = compile(&graph, swapped).expect("compiles");
    let (plain_body, _) = compile(&graph, plain).expect("compiles");
    assert!((volume(&swapped_body) - 1.8).abs() <= 1e-9 * 1.8);
    assert_eq!(swapped_body, plain_body);
    let (regional_body, _) = compile(&graph, regional).expect("compiles");
    assert!((volume(&regional_body) - CORNER).abs() <= 1e-9 * CORNER);
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
