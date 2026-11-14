//! A difference or an intersection the mesh cut cannot settle is decided by
//! the exact compiler (#423, ADR 0098), as a union is (#391).
//!
//! A door cutter as thick as its wall, placed with the wall under a building
//! placement, meets the wall's faces only up to rounding. The mesh boolean,
//! which has no tolerance, keeps the cut-out slug as a closed shell of its
//! own. The exact compiler reads the faces as one surface within the
//! tolerance and cuts the opening through; its result is tessellated.
//!
//! The referee is never worse than the mesh path: where it refuses, the
//! mesh result stands.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{
    BooleanOperator, Interval, Mat3, Plane3, Point2, Point3, Scalar, Tolerance, Transform2,
    Transform3, Vec2, Vec3,
};
use axiolid_curve::{Curve2, Line2};
use axiolid_inspect::enclosed_volume;
use axiolid_mesh::{component_count, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::{
    REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT, REFEREED_INTERSECTION,
    REFEREED_INTERSECTION_EXACT,
};
use axiolid_mesh_compile::{DeviationBound, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

const EPS: Scalar = 1e-6;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(EPS, 1e-9).expect("a tolerance"))
        .with_chord_error(1e-3)
        .expect("a positive chord error")
}

fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
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
        self.extrude(profile, depth)
    }

    fn extrude(&mut self, profile: NodeId, depth: Scalar) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn placed(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
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

/// The wall `3.65 x 0.25 x 3.67` less a door opening `1.01 x 2.26` through
/// it at `x`, each operand placed by `world` (the construction of #312).
fn wall_with_door(world: Transform3, x: Scalar) -> (GeometryGraph, NodeId) {
    wall_with_doors(world, 3.65, &[x])
}

/// A wall `length x 0.25 x 3.67` less a door opening at each of `xs`.
fn wall_with_doors(world: Transform3, length: Scalar, xs: &[Scalar]) -> (GeometryGraph, NodeId) {
    let mut g = Graph::new();
    let derived = |g: &mut Graph, size: [Scalar; 2], transform: Transform2| {
        let basis = Profile::Rectangle(RectangleProfile {
            x: size[0],
            y: size[1],
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        });
        g.push(GeometryNode::Profile(Profile::Derived {
            basis: Box::new(basis),
            transform,
        }))
    };
    let wall = derived(
        &mut g,
        [length, 0.25],
        Transform2::from_cols(Vec2::NEG_X, Vec2::NEG_Y, Vec2::new(0.5 * length, 0.0)),
    );
    let wall = g.extrude(wall, 3.67);
    let mut root = g.placed(wall, world);
    for x in xs {
        let door = derived(
            &mut g,
            [1.01, 2.26],
            Transform2::from_cols(Vec2::Y, Vec2::NEG_X, Vec2::new(1.13, 0.505)),
        );
        let door = g.extrude(door, 0.25);
        let door = g.placed(
            door,
            world
                * Transform3::from_cols(
                    Vec3::Z,
                    Vec3::NEG_X,
                    Vec3::NEG_Y,
                    Vec3::new(*x, 0.125, 0.2),
                ),
        );
        root = g.boolean(root, door, BooleanOperator::Difference);
    }
    g.finish(root)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> (TriMesh, MeshClosure, DeviationReport) {
    compile_at(graph, root, &options())
}

fn compile_at(
    graph: &GeometryGraph,
    root: NodeId,
    options: &ExecutionOptions,
) -> (TriMesh, MeshClosure, DeviationReport) {
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(graph, root, options)
        .expect("compiles");
    (outcome.mesh, outcome.closure, report)
}

/// The tolerance of a building model, a thousandth.
fn millimetre() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

/// Whether the report names an extra shell as unbounded (#392).
fn names_extra_shell(report: &DeviationReport) -> bool {
    report.contributions.iter().any(
        |c| matches!(c.bound, DeviationBound::Unbounded(reason) if reason.contains("more shells")),
    )
}

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

/// The rectangle from `(x0, y0)` to `(x1, y1)`.
fn corners(x0: Scalar, y0: Scalar, x1: Scalar, y1: Scalar) -> Profile {
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

/// `profile` extruded by `depth` and placed.
fn slab(g: &mut Graph, profile: Profile, depth: Scalar, transform: Transform3) -> NodeId {
    let profile = g.push(GeometryNode::Profile(profile));
    let solid = g.extrude(profile, depth);
    g.placed(solid, transform)
}

fn has(report: &DeviationReport, detail: &str) -> bool {
    report.contributions.iter().any(|c| c.detail == detail)
}

fn refereed(report: &DeviationReport, details: [&str; 2]) -> bool {
    details.iter().any(|detail| has(report, detail))
}

fn volume(mesh: &TriMesh) -> Scalar {
    let volume = enclosed_volume(mesh).expect("a measurable solid");
    0.5 * (volume.lower + volume.upper)
}

/// The door opening through the placed wall leaves one body whatever the
/// door's position: the slug the mesh cut left (#423) is gone.
#[test]
fn a_placed_flush_door_leaves_one_certified_body() {
    for x in [1.35, 3.325] {
        for world in [Transform3::IDENTITY, building()] {
            let (graph, root) = wall_with_door(world, x);
            let (mesh, closure, report) = compile(&graph, root);
            assert_eq!(component_count(&mesh), 1, "x {x}, {world:?}");
            assert_eq!(closure, MeshClosure::Solid);
            assert!(
                report
                    .contributions
                    .iter()
                    .all(|c| matches!(c.bound, DeviationBound::Certified(_))),
                "x {x}, {world:?}: {report:?}"
            );
            let expected = 3.65 * 0.25 * 3.67 - 1.01 * 2.26 * 0.25;
            assert!(
                (volume(&mesh) - expected).abs() < 1e-5,
                "x {x}, {world:?}: {}",
                volume(&mesh)
            );
        }
    }
}

/// The placed opening is decided by the exact compiler, and says so.
#[test]
fn a_placed_flush_door_is_decided_by_the_referee() {
    let (graph, root) = wall_with_door(building(), 1.35);
    let (_, _, report) = compile(&graph, root);
    assert!(
        refereed(&report, [REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT]),
        "{report:?}"
    );
}

/// Two unit boxes, the second placed by the same turn composed another way
/// and `lift` above the first.
fn boxes(operator: BooleanOperator, lift: Scalar) -> (GeometryGraph, NodeId) {
    let a = Transform3::from_rotation_z(0.3);
    let half = Transform3::from_rotation_z(0.15);
    let b = half * half;
    assert_ne!(a.matrix3, b.matrix3, "the placements must differ in bits");
    let mut g = Graph::new();
    let first = g.prism([1.0, 1.0], 1.0);
    let first = g.placed(first, a);
    let second = g.prism([1.0, 1.0], 1.0);
    let second = g.placed(
        second,
        b * Transform3::from_translation(Vec3::new(0.0, 0.0, lift)),
    );
    let root = g.boolean(first, second, operator);
    g.finish(root)
}

/// Near contact alone does not call the referee for an intersection: its
/// mesh has one shell, as the exact result does, so the mesh result stands.
#[test]
fn an_intersection_in_near_contact_with_one_shell_keeps_the_mesh_result() {
    let (graph, root) = boxes(BooleanOperator::Intersection, 0.5 * EPS);
    let (mesh, closure, report) = compile(&graph, root);
    assert_eq!(closure, MeshClosure::Solid);
    assert_eq!(component_count(&mesh), 1);
    assert!((volume(&mesh) - 1.0).abs() < 1e-5, "{}", volume(&mesh));
    assert!(
        !refereed(
            &report,
            [REFEREED_INTERSECTION, REFEREED_INTERSECTION_EXACT]
        ),
        "{report:?}"
    );
}

/// A window `0.25 - 2e` deep in the middle of a `0.25` wall, its caps `e`
/// inside the wall's faces (the construction of `coincident_wall_openings`):
/// two skins `e` thick enclose a void.
fn wall_with_skinned_window(world: Transform3, e: Scalar) -> (GeometryGraph, NodeId) {
    let mut g = Graph::new();
    let wall = slab(&mut g, corners(0.0, 0.0, 5.0, 0.25), 3.0, world);
    let across =
        Transform3::from_cols(Vec3::X, Vec3::Z, Vec3::NEG_Y, Vec3::new(0.0, 0.25 - e, 0.0));
    let window = slab(
        &mut g,
        corners(2.5, 1.0, 3.5, 2.0),
        0.25 - 2.0 * e,
        world * across,
    );
    let root = g.boolean(wall, window, BooleanOperator::Difference);
    g.finish(root)
}

/// Skins thinner than the tolerance that enclose a void are authored (ADR
/// 0096, D7, #291): the mesh holds a shell the exact result, reading the
/// skins away, lacks, but the reading is far above rounding, so the referee
/// does not replace the mesh, and the report is unbounded by name (#392).
#[test]
fn an_authored_void_within_the_tolerance_survives() {
    for world in [Transform3::IDENTITY, building()] {
        for e in [1e-4, 4e-4] {
            let (graph, root) = wall_with_skinned_window(world, e);
            let (mesh, closure, report) = compile_at(&graph, root, &millimetre());
            assert_eq!(closure, MeshClosure::Solid);
            // The wall less the void: the authored volume, not the through
            // window's.
            let authored = 3.75 - (0.25 - 2.0 * e);
            assert!(
                (volume(&mesh) - authored).abs() < 1e-6,
                "e {e:e}: {} against {authored}",
                volume(&mesh)
            );
            assert!(
                !refereed(&report, [REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT]),
                "e {e:e}: {report:?}"
            );
            assert!(names_extra_shell(&report), "e {e:e}: {report:?}");
        }
    }
}

/// The same void, cut in a frame the wall and the cavity share: a cavity
/// `e` inside the wall's faces.
#[test]
fn a_void_in_a_shared_frame_is_unbounded_by_name() {
    for world in [Transform3::IDENTITY, building()] {
        let e = 1e-4;
        let mut g = Graph::new();
        let wall = slab(&mut g, corners(0.0, 0.0, 5.0, 0.25), 3.0, world);
        let cavity = slab(
            &mut g,
            corners(2.5, e, 3.5, 0.25 - e),
            1.0,
            world * Transform3::from_translation(Vec3::new(0.0, 0.0, 1.0)),
        );
        let root = g.boolean(wall, cavity, BooleanOperator::Difference);
        let (graph, root) = g.finish(root);
        let (mesh, _, report) = compile_at(&graph, root, &millimetre());
        let authored = 3.75 - 1.0 * (0.25 - 2.0 * e);
        assert!((volume(&mesh) - authored).abs() < 1e-6, "{}", volume(&mesh));
        assert!(
            !refereed(&report, [REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT]),
            "{report:?}"
        );
        assert!(names_extra_shell(&report), "{report:?}");
    }
}

/// A cutter that leaves a sliver `w` thick detached from the wall: the same
/// shell count, and the same name.
#[test]
fn a_sliver_in_a_shared_frame_is_unbounded_by_name() {
    for world in [Transform3::IDENTITY, building()] {
        let w = 1e-4;
        let mut g = Graph::new();
        let wall = slab(&mut g, corners(0.0, 0.0, 5.0, 0.25), 3.0, world);
        let cutter = slab(
            &mut g,
            corners(w, -1.0, 2.5, 1.0),
            5.0,
            world * Transform3::from_translation(Vec3::new(0.0, 0.0, -1.0)),
        );
        let root = g.boolean(wall, cutter, BooleanOperator::Difference);
        let (graph, root) = g.finish(root);
        let (mesh, _, report) = compile_at(&graph, root, &millimetre());
        assert!(component_count(&mesh) >= 2, "{}", component_count(&mesh));
        assert!(
            !refereed(&report, [REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT]),
            "{report:?}"
        );
        assert!(names_extra_shell(&report), "{report:?}");
    }
}

/// A half-space tool is clipped to the subject and is no node to compile
/// exactly. A notched wall cut by a plane a nanometre above its notch floor
/// leaves two prongs or a slab: near contact within the tolerance and
/// two components, where the referee would be asked, but never for a
/// half-space.
#[test]
fn a_half_space_tool_is_never_refereed() {
    let world = building();
    for operator in [BooleanOperator::Difference, BooleanOperator::Intersection] {
        for sign in [1.0, -1.0] {
            let mut g = Graph::new();
            let wall = slab(&mut g, corners(0.0, 0.0, 5.0, 0.25), 3.0, world);
            let notch = slab(
                &mut g,
                corners(2.0, -1.0, 3.0, 1.0),
                3.0,
                world * Transform3::from_translation(Vec3::new(0.0, 0.0, 1.0)),
            );
            let notched = g.boolean(wall, notch, BooleanOperator::Difference);
            let tool = g.push(GeometryNode::HalfSpace(HalfSpace {
                boundary: Plane3 {
                    origin: world.transform_point3(Point3::new(0.0, 0.0, 1.0 + 1e-9)),
                    normal: world.transform_vector3(Vec3::new(0.0, 0.0, sign)),
                },
                agreement: true,
            }));
            let root = g.boolean(notched, tool, operator);
            let (graph, root) = g.finish(root);
            let (_, _, report) = compile_at(&graph, root, &millimetre());
            assert!(
                !refereed(&report, [REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT])
                    && !refereed(
                        &report,
                        [REFEREED_INTERSECTION, REFEREED_INTERSECTION_EXACT]
                    ),
                "{operator:?} {sign}: {report:?}"
            );
        }
    }
}

/// A sheared placement is no similarity, and the exact compiler refuses it:
/// the mesh cut's result stands, with no referee in its report.
#[test]
fn a_sheared_difference_the_referee_refuses_keeps_the_mesh_result() {
    let shear = Mat3::from_cols(Vec3::X, Vec3::new(0.3, 1.0, 0.0), Vec3::Z);
    let site = Vec3::new(1.0, 2.0, 0.0);
    let placed = Transform3::from_mat3_translation(shear, site);
    let mut g = Graph::new();
    let first = g.prism([1.0, 1.0], 1.0);
    let first = g.placed(first, placed);
    let second = g.prism([1.0, 1.0], 1.0);
    // Half the box along its x axis, and half the tolerance off its y face.
    let second = g.placed(
        second,
        placed * Transform3::from_translation(Vec3::new(0.5, 0.5 * EPS, 0.0)),
    );
    let root = g.boolean(first, second, BooleanOperator::Difference);
    let (graph, root) = g.finish(root);
    let (mesh, closure, report) = compile(&graph, root);
    assert_eq!(closure, MeshClosure::Solid);
    assert!((volume(&mesh) - 0.5).abs() < 1e-5, "{}", volume(&mesh));
    assert!(
        !refereed(&report, [REFEREED_DIFFERENCE, REFEREED_DIFFERENCE_EXACT]),
        "{report:?}"
    );
}
