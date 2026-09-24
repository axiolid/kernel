//! The boundary of a polygonally bounded half-space carries its own frame.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Plane3, Point2, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Polyline2};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, SolidOperation};
use axiolid_primitive::HalfSpace;

fn compile_bounded_half_space(boundary: &[Point2], placement: Transform3) -> TriMesh {
    let mut builder = GeometryGraphBuilder::new();
    let half_space = builder
        .push(GeometryNode::HalfSpace(HalfSpace {
            // A TILTED clip plane. With normal = +Z the internal heuristic
            // happens to pick the world x/y axes, so an authored frame about
            // that same normal is indistinguishable from rotating the result.
            // Tilting makes the guessed basis genuinely differ from the
            // authored one, which is the case IFC files actually hit.
            boundary: Plane3 {
                origin: Vec3::ZERO,
                normal: Vec3::new(0.0, 1.0, 1.0).normalize(),
            },
            agreement: true,
        }))
        .expect("a half-space is a valid node");
    let curve = builder
        .push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: boundary.to_vec(),
            closed: true,
        })))
        .expect("a 2d polyline is a valid node");
    let op = builder
        .push(GeometryNode::SolidOperation(
            SolidOperation::BoundedHalfSpace {
                half_space,
                boundary: curve,
                placement,
            },
        ))
        .expect("a bounded half-space is a valid operation");
    let graph = builder.finish(vec![op]).expect("the graph is valid");
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, op, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .expect("a bounded half-space compiles")
}

/// The graph must carry the boundary's own in-plane orientation.
///
/// `IfcPolygonalBoundedHalfSpace.Position` is schema-defined as independent of
/// `BaseSurface`: it need not share an origin or an orientation with the clip
/// plane. Applying it only to the finished mesh cannot express that, because
/// by then the boundary has already been framed against the plane normal.
#[test]
fn the_boundary_placement_orients_the_profile_not_just_the_result() {
    // An L-shaped boundary, so a rotation is visible in the footprint rather
    // than being absorbed by symmetry.
    let boundary = vec![
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(2.0, 1.0),
        Point2::new(1.0, 1.0),
        Point2::new(1.0, 2.0),
        Point2::new(0.0, 2.0),
    ];

    // Same clip plane, same boundary; only the authored placement differs.
    let quarter_turn = compile_bounded_half_space(
        &boundary,
        // A quarter turn about z, written as the rotated basis directly so the
        // test does not need a quaternion dependency.
        Transform3::from_cols(
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::ZERO,
        ),
    );

    // The expected corner is derived from the authored frame, not hardcoded:
    // the profile point (2,0) lies at origin + x_axis*2 in whatever frame the
    // placement declares. Under the old post-hoc transform the boundary was
    // framed by the plane-normal heuristic instead, so the corner landed
    // somewhere the authored frame never named.
    // The authored x is projected into the clip plane and renormalised, so the
    // frame contributes exactly its rotation about the normal and nothing else.
    // Deriving the expectation the same way keeps this test honest about what
    // the contract promises rather than pinning a magic coordinate.
    let normal = Vec3::new(0.0, 1.0, 1.0).normalize();
    let authored_x = Vec3::new(0.0, 1.0, 0.0);
    let in_plane = (authored_x - normal * authored_x.dot(normal)).normalize();
    let want = in_plane * 2.0;
    assert!(
        has_vertex_at(&quarter_turn, want),
        "the authored placement must orient the PROFILE: expected a vertex at \
         {want:?}, got {:?}",
        quarter_turn.positions
    );
}

/// Whether any mesh vertex sits at `want`.
fn has_vertex_at(mesh: &TriMesh, want: Vec3) -> bool {
    mesh.positions.iter().any(|p| (*p - want).length() < 1e-9)
}

/// The placement's translation places the boundary, not only its rotation
/// (#164). A consumer authoring the boundary frame away from the clip plane's
/// origin would otherwise cut an unrelated region with no error.
#[test]
fn the_boundary_placement_translates_the_profile() {
    let boundary = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ];
    let shift = Vec3::new(4.0, -3.0, 0.0);
    let moved = compile_bounded_half_space(
        &boundary,
        Transform3::from_cols(Vec3::X, Vec3::Y, Vec3::Z, shift),
    );

    // The authored translation projected into the tilted clip plane is where
    // profile point (0, 0) must land; x stays along the projected world x.
    let normal = Vec3::new(0.0, 1.0, 1.0).normalize();
    let anchor = shift - normal * shift.dot(normal);
    assert!(
        has_vertex_at(&moved, anchor),
        "the authored translation must place the PROFILE: expected a vertex at \
         {anchor:?}, got {:?}",
        moved.positions
    );
}
