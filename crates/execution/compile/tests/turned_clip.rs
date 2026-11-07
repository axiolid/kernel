//! A wall turned and placed far from the origin, cut by a half-space clip
//! that reaches through its whole thickness (#400).
//!
//! The wall is `2.5 x 0.125 x 3.6`, its faces at `y = 0` and `y = 0.125`;
//! the clip keeps `z >= 2.7` within the boundary `x` in `1.1..1.2`, `y` from
//! `y_lo` to `0.125`, so it opens a slot `0.1` wide through the wall's
//! thickness. Placed by a turn, the clip's planes in world coordinates lie a
//! rounding residue off the wall's: the mesh cut left a lid of that thickness
//! across the slot and still called the result a Solid. A clip is a
//! constructed stand-in with no coordinates of its own, so it is cut in its
//! subject's frame, where the wall's faces are exactly axis-aligned (ADR 0096
//! D3, extended to half-space tools in ADR 0102).

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{
    BooleanOperator, Plane3, Point2, Point3, Scalar, Tolerance, Transform2, Transform3, Vec2, Vec3,
};
use axiolid_curve::{Curve2, Polyline2};
use axiolid_inspect::{contains, enclosed_volume, thin_sheet};
use axiolid_mesh::{audit_mesh, decompose, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{ReferenceMeshCompiler, SKIN_WITHIN_ROUNDING};
use axiolid_mesh_compile_contract::MeshClosure;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Profile, RectangleProfile};

const WIDTH: Scalar = 2.5;
const THICKNESS: Scalar = 0.125;
const HEIGHT: Scalar = 3.6;
const SLAB_FROM: Scalar = 2.7;
/// The slot the clip opens, in the wall's frame.
const SLOT: (Scalar, Scalar) = (1.1, 1.2);

fn gross() -> Scalar {
    WIDTH * THICKNESS * HEIGHT
}

/// The wall less the slot through its thickness, `0.1` wide, `0.9` tall.
fn net() -> Scalar {
    gross() - (SLOT.1 - SLOT.0) * THICKNESS * (HEIGHT - SLAB_FROM)
}

fn turn(site: Vec3) -> Transform3 {
    Transform3::from_translation(site) * Transform3::from_rotation_z(10.0_f64.to_radians())
}

fn near() -> Transform3 {
    turn(Vec3::new(13.475, -15.95, 4.3))
}

fn far() -> Transform3 {
    turn(Vec3::new(6.0e5, 5.6e6, 200.0))
}

/// How the clip meets the wall.
#[derive(Clone, Copy)]
enum Scene {
    /// `Difference { Instance { subject, wall }, Instance { tool, clip } }`.
    Placed {
        subject: Transform3,
        tool: Transform3,
    },
    /// The clip written in world coordinates: its plane and its frame turned
    /// by `subject`, with no `Instance`.
    Baked { subject: Transform3 },
    /// `Instance { frame, Difference { wall, clip } }`: the cut in the wall's
    /// own frame.
    Whole { frame: Transform3 },
    /// Both turned in world coordinates, with no `Instance`: the wall's
    /// profile carries the turn and so does the clip's frame (#400).
    World { frame: Transform3 },
}

impl Scene {
    /// The frame the wall is placed by.
    fn frame(self) -> Transform3 {
        match self {
            Self::Placed { subject, .. } | Self::Baked { subject } => subject,
            Self::Whole { frame } | Self::World { frame } => frame,
        }
    }
}

fn push(builder: &mut GeometryGraphBuilder, node: GeometryNode) -> NodeId {
    builder.push(node).expect("a valid node")
}

fn wall(builder: &mut GeometryGraphBuilder, placed: Transform2) -> NodeId {
    let profile = push(
        builder,
        GeometryNode::Profile(Profile::Derived {
            basis: Box::new(Profile::Rectangle(RectangleProfile {
                x: WIDTH,
                y: THICKNESS,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            })),
            transform: placed,
        }),
    );
    push(
        builder,
        GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: HEIGHT,
        }),
    )
}

/// The clip: `z >= 2.7` along `z`, within the boundary, its plane through
/// `origin` and its frame `placement`.
fn clip(
    builder: &mut GeometryGraphBuilder,
    (origin, normal): (Vec3, Vec3),
    placement: Transform3,
    y: (Scalar, Scalar),
) -> NodeId {
    let half_space = push(
        builder,
        GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 { origin, normal },
            agreement: true,
        }),
    );
    let boundary = push(
        builder,
        GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::new(SLOT.0, y.0),
                Point2::new(SLOT.1, y.0),
                Point2::new(SLOT.1, y.1),
                Point2::new(SLOT.0, y.1),
            ],
            closed: true,
        })),
    );
    push(
        builder,
        GeometryNode::SolidOperation(SolidOperation::BoundedHalfSpaceAlongZ {
            half_space,
            boundary,
            placement,
        }),
    )
}

fn graph(scene: Scene, y_lo: Scalar) -> (GeometryGraph, NodeId) {
    graph_between(scene, (y_lo, THICKNESS))
}

fn graph_between(scene: Scene, y: (Scalar, Scalar)) -> (GeometryGraph, NodeId) {
    graph_less(scene, y, false)
}

/// [`graph_between`], less an intersection that leaves nothing when
/// `nothing` (#326): a box and a half-space above it.
fn graph_less(scene: Scene, y: (Scalar, Scalar), nothing: bool) -> (GeometryGraph, NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let in_frame = Transform2::from_cols(Vec2::X, Vec2::Y, Vec2::new(WIDTH / 2.0, THICKNESS / 2.0));
    let wall = wall(
        &mut builder,
        match scene {
            // The profile's own frame, turned and placed by the same
            // transform as the clip's, in plan.
            Scene::World { frame } => Transform2::from_cols(
                Vec2::new(frame.matrix3.x_axis.x, frame.matrix3.x_axis.y),
                Vec2::new(frame.matrix3.y_axis.x, frame.matrix3.y_axis.y),
                Vec2::new(
                    frame
                        .transform_point3(Vec3::new(WIDTH / 2.0, THICKNESS / 2.0, 0.0))
                        .x,
                    frame
                        .transform_point3(Vec3::new(WIDTH / 2.0, THICKNESS / 2.0, 0.0))
                        .y,
                ),
            ),
            _ => in_frame,
        },
    );
    let local_plane = (Vec3::new(0.0, 0.0, SLAB_FROM), Vec3::Z);
    let root = match scene {
        Scene::Placed { subject, tool } => {
            let clip = clip(&mut builder, local_plane, Transform3::IDENTITY, y);
            let wall = push(
                &mut builder,
                GeometryNode::Instance(Instance {
                    source: wall,
                    transform: subject,
                }),
            );
            let clip = push(
                &mut builder,
                GeometryNode::Instance(Instance {
                    source: clip,
                    transform: tool,
                }),
            );
            difference(&mut builder, wall, clip)
        }
        Scene::Baked { subject } => {
            let plane = (
                subject.transform_point3(local_plane.0),
                subject.matrix3 * local_plane.1,
            );
            let clip = clip(&mut builder, plane, subject, y);
            let wall = push(
                &mut builder,
                GeometryNode::Instance(Instance {
                    source: wall,
                    transform: subject,
                }),
            );
            difference(&mut builder, wall, clip)
        }
        Scene::Whole { frame } => {
            let clip = clip(&mut builder, local_plane, Transform3::IDENTITY, y);
            let cut = difference(&mut builder, wall, clip);
            push(
                &mut builder,
                GeometryNode::Instance(Instance {
                    source: cut,
                    transform: frame,
                }),
            )
        }
        Scene::World { frame } => {
            let plane = (local_plane.0, local_plane.1);
            let clip = clip(&mut builder, plane, frame, y);
            difference(&mut builder, wall, clip)
        }
    };
    let root = if nothing {
        let profile = push(
            &mut builder,
            GeometryNode::Profile(Profile::Rectangle(RectangleProfile {
                x: 1.0,
                y: 1.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            })),
        );
        let block = push(
            &mut builder,
            GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile,
                direction: Vec3::Z,
                depth: 1.0,
            }),
        );
        let above = push(
            &mut builder,
            GeometryNode::HalfSpace(HalfSpace {
                boundary: Plane3 {
                    origin: Point3::new(0.0, 0.0, 10.0),
                    normal: Vec3::Z,
                },
                agreement: true,
            }),
        );
        let empty = push(
            &mut builder,
            GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: block,
                right: above,
                operator: BooleanOperator::Intersection,
            }),
        );
        difference(&mut builder, root, empty)
    } else {
        root
    };
    (builder.finish(vec![root]).expect("a valid graph"), root)
}

fn difference(builder: &mut GeometryGraphBuilder, left: NodeId, right: NodeId) -> NodeId {
    push(
        builder,
        GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator: BooleanOperator::Difference,
        }),
    )
}

fn try_compile(
    graph: &GeometryGraph,
    root: NodeId,
) -> Result<(TriMesh, MeshClosure), axiolid_contracts::GeomError> {
    let (outcome, _) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(graph, root, &ExecutionOptions::new(Tolerance::MILLIMETRE))?;
    Ok((outcome.mesh, outcome.closure))
}

fn compile(graph: &GeometryGraph, root: NodeId) -> (TriMesh, MeshClosure) {
    try_compile(graph, root).unwrap_or_else(|e| panic!("{e}"))
}

/// Triangles whose centroid, in the wall's frame, lies inside the open slot:
/// across it, where the result must be open. The slot's own walls and floor,
/// and the wall's top, lie on its boundary and are kept out by the margin.
fn lid_triangles(mesh: &TriMesh, frame: Transform3) -> usize {
    let inverse = frame.inverse();
    let margin = 1e-3;
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .filter(|triangle| {
            let centroid = triangle
                .iter()
                .map(|&i| mesh.positions[i as usize])
                .sum::<Point3>()
                / 3.0;
            let local = inverse.transform_point3(centroid);
            local.x > SLOT.0 + margin
                && local.x < SLOT.1 - margin
                && local.z > SLAB_FROM + margin
                && local.z < HEIGHT - margin
        })
        .count()
}

/// The result is the wall with its slot open: one closed solid, measured
/// at `net()` within `slack`, with no thin sheet and no lid.
fn assert_open_slot(scene: Scene, slack: Scalar, what: &str) {
    let (graph, root) = graph(scene, 0.0);
    let (mesh, closure) = compile(&graph, root);
    assert_eq!(closure, MeshClosure::Solid, "{what}");
    assert!(
        audit_mesh(&mesh, Tolerance::MILLIMETRE).is_closed_two_manifold(),
        "{what}"
    );
    assert_eq!(decompose(&mesh).len(), 1, "{what}");
    let volume = enclosed_volume(&mesh).unwrap_or_else(|e| panic!("{what}: {e:?}"));
    assert!(
        volume.lower - slack <= net() && net() <= volume.upper + slack,
        "{what}: [{}, {}] does not hold {} within {slack:e}",
        volume.lower,
        volume.upper,
        net()
    );
    assert!(
        axiolid_heal::self_intersections(&mesh).is_empty(),
        "{what}: the mesh crosses itself"
    );
    assert_eq!(thin_sheet(&mesh, 1e-9), None, "{what}");
    assert_eq!(lid_triangles(&mesh, scene.frame()), 0, "{what}: a lid");
}

#[test]
fn a_clip_through_a_turned_wall_opens_the_slot() {
    let t = near();
    assert_open_slot(
        Scene::Placed {
            subject: t,
            tool: t,
        },
        1e-9,
        "placed",
    );
}

#[test]
fn a_clip_written_in_world_coordinates_opens_the_slot() {
    assert_open_slot(Scene::Baked { subject: near() }, 1e-9, "baked");
}

#[test]
fn a_clip_through_a_wall_at_a_site_opens_the_slot() {
    let t = far();
    // The coordinates round to about 1e-9 out there.
    assert_open_slot(
        Scene::Placed {
            subject: t,
            tool: t,
        },
        1e-6,
        "site",
    );
}

#[test]
fn a_skin_above_the_rounding_stays_authored() {
    let t = near();
    // The mesh cut decides it: the guard reads the rounding, not the
    // tolerance, so the skin never reaches the referee.
    let scene = Scene::Placed {
        subject: t,
        tool: t,
    };
    let (_, _, detail) = compile_reported(scene, (0.0005, THICKNESS)).expect("compiles");
    assert!(!detail.starts_with(REFEREED), "{detail}");
    let (graph, root) = graph(
        Scene::Placed {
            subject: t,
            tool: t,
        },
        0.0005,
    );
    let (mesh, closure) = compile(&graph, root);
    assert_eq!(closure, MeshClosure::Solid);
    assert!(audit_mesh(&mesh, Tolerance::MILLIMETRE).is_closed_two_manifold());
    assert!(
        lid_triangles(&mesh, t) > 0,
        "the half-millimetre skin stays"
    );
    let volume = enclosed_volume(&mesh).expect("measurable");
    let expected = gross() - (SLOT.1 - SLOT.0) * (THICKNESS - 0.0005) * (HEIGHT - SLAB_FROM);
    assert!(
        volume.lower - 1e-9 <= expected && expected <= volume.upper + 1e-9,
        "[{}, {}] does not hold {expected}",
        volume.lower,
        volume.upper
    );
}

#[test]
fn a_tool_offset_in_the_subjects_frame_moves_the_slot() {
    let t = near();
    let shifted = Transform3 {
        matrix3: t.matrix3,
        translation: t.translation + t.matrix3 * Vec3::new(0.5, 0.0, 0.0),
    };
    let scene = Scene::Placed {
        subject: t,
        tool: shifted,
    };
    let (graph, root) = graph(scene, 0.0);
    let (mesh, closure) = compile(&graph, root);
    assert_eq!(closure, MeshClosure::Solid);
    let volume = enclosed_volume(&mesh).expect("measurable");
    assert!(
        volume.lower - 1e-9 <= net() && net() <= volume.upper + 1e-9,
        "[{}, {}]",
        volume.lower,
        volume.upper
    );
    let inside = |local: Vec3| contains(&mesh, t.transform_point3(local)).expect("decided");
    let (moved, home) = (Vec3::new(1.65, 0.06, 3.2), Vec3::new(1.15, 0.06, 3.2));
    assert!(!inside(moved), "the slot is at 1.6 to 1.7");
    assert!(inside(home), "the wall is whole at 1.1 to 1.2");
    assert_eq!(thin_sheet(&mesh, 1e-9), None);
}

#[test]
fn the_clip_is_cut_as_the_wall_cut_in_its_own_frame_and_placed() {
    let t = near();
    let (a, root_a) = graph(
        Scene::Placed {
            subject: t,
            tool: t,
        },
        0.0,
    );
    let (b, root_b) = graph(Scene::Whole { frame: t }, 0.0);
    let (a, _) = compile(&a, root_a);
    let (b, _) = compile(&b, root_b);
    assert_eq!(a.positions, b.positions);
    assert_eq!(a.indices, b.indices);
}

/// The wall and the clip both turned in world coordinates, with no
/// `Instance` to give the cut a frame: the mesh cut is the only one, and
/// whether its lid closes depends on how each turn rounds. Every angle is
/// either a clean solid or refused by name, never a Solid with a lid.
#[test]
fn a_turned_wall_and_clip_are_a_clean_solid_or_refused_by_name() {
    for degrees in 1..=40 {
        let frame = Transform3::from_translation(Vec3::new(13.475, -15.95, 0.0))
            * Transform3::from_rotation_z(f64::from(degrees).to_radians());
        let (graph, root) = graph(Scene::World { frame }, 0.0);
        match try_compile(&graph, root) {
            Ok((mesh, closure)) => {
                assert_eq!(closure, MeshClosure::Solid, "{degrees} degrees");
                assert_eq!(thin_sheet(&mesh, 1e-9), None, "{degrees} degrees");
                assert!(
                    axiolid_heal::self_intersections(&mesh).is_empty(),
                    "{degrees} degrees: the mesh crosses itself"
                );
            }
            Err(GeomError::Degenerate(text)) if text.starts_with(SKIN_WITHIN_ROUNDING) => {}
            Err(other) => panic!("{degrees} degrees: {other}"),
        }
    }
}

/// The wall and the clip turned `degrees` about `z` at `(x, y)`, with no
/// `Instance`: its frame.
fn world(x: Scalar, y: Scalar, degrees: i32) -> Transform3 {
    Transform3::from_translation(Vec3::new(x, y, 0.0))
        * Transform3::from_rotation_z(f64::from(degrees).to_radians())
}

fn compile_reported(
    scene: Scene,
    y: (Scalar, Scalar),
) -> Result<(TriMesh, MeshClosure, String), GeomError> {
    let (graph, root) = graph_between(scene, y);
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &ExecutionOptions::new(Tolerance::MILLIMETRE))?;
    let detail = report.contributions[0].detail.to_owned();
    Ok((outcome.mesh, outcome.closure, detail))
}

const REFEREED: &str = "tessellated from the exact compiler's difference";

/// At 7 degrees the mesh cut leaves the lid (it is refused as a skin within
/// rounding); the exact referee decides the cut instead and opens the slot.
#[test]
fn the_referee_opens_the_slot_where_the_mesh_cut_keeps_a_lid() {
    let frame = world(13.475, -15.95, 7);
    let (mesh, closure, detail) =
        compile_reported(Scene::World { frame }, (0.0, THICKNESS)).expect("decided");
    assert!(detail.starts_with(REFEREED), "{detail}");
    assert_eq!(closure, MeshClosure::Solid);
    assert!(audit_mesh(&mesh, Tolerance::MILLIMETRE).is_closed_two_manifold());
    assert_eq!(decompose(&mesh).len(), 1);
    let volume = enclosed_volume(&mesh).expect("measurable");
    assert!(
        volume.lower - 1e-9 <= net() && net() <= volume.upper + 1e-9,
        "[{}, {}] does not hold {}",
        volume.lower,
        volume.upper,
        net()
    );
    assert_eq!(thin_sheet(&mesh, 1e-9), None);
    assert_eq!(lid_triangles(&mesh, frame), 0, "a lid");
}

/// A skin above the reach stays authored when the referee decides the rest:
/// the boundary stops half a millimetre short of the wall's far face.
#[test]
fn the_referee_keeps_an_authored_skin() {
    let frame = world(13.475, -15.95, 3);
    let (mesh, closure, detail) =
        compile_reported(Scene::World { frame }, (0.0, 0.1245)).expect("decided");
    assert!(detail.starts_with(REFEREED), "{detail}");
    assert_eq!(closure, MeshClosure::Solid);
    assert!(
        lid_triangles(&mesh, frame) > 0,
        "the half-millimetre skin stays"
    );
    let volume = enclosed_volume(&mesh).expect("measurable");
    let expected = gross() - (SLOT.1 - SLOT.0) * 0.1245 * (HEIGHT - SLAB_FROM);
    assert!(
        volume.lower - 1e-9 <= expected && expected <= volume.upper + 1e-9,
        "[{}, {}] does not hold {expected}",
        volume.lower,
        volume.upper
    );
}

/// Pinned by search: at this site the exact compiler cannot classify a
/// point near a face boundary, so there is no referee result and the skin
/// is refused by name with the referee's reason.
#[test]
fn a_skin_the_exact_compiler_cannot_decide_is_refused_by_name() {
    let frame = world(6.0e5, 5.6e6, 1);
    let error = compile_reported(Scene::World { frame }, (0.0, THICKNESS))
        .expect_err("the exact compiler refuses it");
    let GeomError::Degenerate(text) = error else {
        panic!("a degenerate refusal: {error}");
    };
    assert!(text.starts_with(SKIN_WITHIN_ROUNDING), "{text}");
    assert!(text.contains("the exact referee refused"), "{text}");
}

/// Pinned by search: the exact compiler reads features `2e-4` apart as one
/// at this site, far above the rounding reach. Fusing them would erase what
/// is authored (ADR 0096, D7), so the referee's result is not accepted.
#[test]
fn a_referee_that_reads_features_above_the_reach_is_not_accepted() {
    let frame = world(2.0e5, 1.0e5, 26);
    let error = compile_reported(Scene::World { frame }, (0.0, THICKNESS))
        .expect_err("the referee reads above the reach");
    let GeomError::Degenerate(text) = error else {
        panic!("a degenerate refusal: {error}");
    };
    assert!(text.starts_with(SKIN_WITHIN_ROUNDING), "{text}");
    assert!(text.contains("above the rounding reach"), "{text}");
}

/// A skin the mesh cut leaves pending (#400) waits for the root of its chain
/// of clips; a difference by nothing (#326) is the subject, so it is that
/// root's to decide, and the result is the plain case's, bit for bit.
#[test]
fn a_skin_pending_under_an_empty_tool_is_decided_at_the_root() {
    let frame = world(13.475, -15.95, 7);
    let scene = Scene::World { frame };
    let (plain, _, _) = compile_reported(scene, (0.0, THICKNESS)).expect("decided");
    let (graph, root) = graph_less(scene, (0.0, THICKNESS), true);
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .expect("decided");
    let detail = report.contributions[0].detail;
    assert!(detail.starts_with(REFEREED), "{detail}");
    assert_eq!(outcome.closure, MeshClosure::Solid);
    assert_eq!(outcome.mesh, plain);
}
