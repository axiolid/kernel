//! Section requests: faces in the plane become regions (ADR 0091).
//!
//! Every expectation is derived by hand next to its case.

use axiolid_contracts::{CancellationToken, ExecutionOptions, GeomError, Operation};
use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_dispatch::MeshPlaneSectionRegistry;
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    InPlaneSide, MeshPlaneSection, PlaneSection, RegionLimits, SectionLimits, SectionRequest,
    REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES, REFUSED_OVERLAPPING_IN_PLANE_FACES,
    REFUSED_SADDLE_VERTEX,
};
use axiolid_reference::ScalarSection;

/// A box `x` by `y` by `z` as 8 vertices and 12 outward triangles.
fn cuboid(x: Span, y: Span, z: Span) -> (Vec<Point3>, Vec<u32>) {
    let positions = vec![
        Point3::new(x.0, y.0, z.0),
        Point3::new(x.1, y.0, z.0),
        Point3::new(x.1, y.1, z.0),
        Point3::new(x.0, y.1, z.0),
        Point3::new(x.0, y.0, z.1),
        Point3::new(x.1, y.0, z.1),
        Point3::new(x.1, y.1, z.1),
        Point3::new(x.0, y.1, z.1),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    (positions, indices)
}

/// A closed interval.
type Span = (f64, f64);

/// Several boxes as separate shells of one mesh.
fn shells(parts: &[(Span, Span, Span)]) -> TriMesh {
    let mut positions = Vec::new();
    let mut indices = Vec::new();
    for &(x, y, z) in parts {
        let (box_positions, box_indices) = cuboid(x, y, z);
        let offset = u32::try_from(positions.len()).expect("small fixture");
        positions.extend(box_positions);
        indices.extend(box_indices.into_iter().map(|index| index + offset));
    }
    TriMesh::new(positions, indices)
}

fn cube() -> TriMesh {
    shells(&[((0.0, 1.0), (0.0, 1.0), (0.0, 1.0))])
}

fn frame(origin: Point3, x: Vec3, y: Vec3) -> Frame3 {
    Frame3 {
        origin,
        x,
        y,
        z: x.cross(y),
    }
}

fn horizontal(z: f64) -> Frame3 {
    frame(Point3::new(0.0, 0.0, z), Vec3::X, Vec3::Y)
}

fn request(frame: Frame3) -> SectionRequest {
    SectionRequest::new(
        frame,
        SectionLimits::new(1_000, 1_000, 1_000, 100),
        RegionLimits::new(100, 1_000, 1_000),
    )
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn section(mesh: &TriMesh, frame: Frame3) -> Result<PlaneSection, GeomError> {
    ScalarSection::new().section_request(mesh, &request(frame), &options())
}

fn ring(points: &[(f64, f64)]) -> Vec<Point2> {
    points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
}

#[test]
fn cube_top_in_the_plane_is_one_region_below() {
    // Z = 1: the two top triangles face +z, the frame's +z, so the material
    // is Below. The diagonal vertex pair keeps four corners, no chains.
    let result = section(&cube(), horizontal(1.0)).expect("section");
    assert!(result.chains.is_empty());
    assert_eq!(result.regions.len(), 1);
    assert_eq!(result.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        result.regions[0].outer,
        ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    );
    assert_eq!(result.evidence.in_plane_triangles, 2);
    assert_eq!(result.evidence.touching_edges, 0);
    assert_eq!(result.evidence.in_plane_pairs, 0, "the two share an edge");
}

#[test]
fn a_mirrored_frame_flips_the_side_and_the_ring_stays_counter_clockwise() {
    // Frame x = X, y = -Y, so z = X x -Y = -Z. The cube top faces model +Z,
    // which opposes the frame's +z: Above. Local point of model (x, y) is
    // (x, -y): the corners (0,0), (1,0), (1,-1), (0,-1). The drop of axis 2
    // reverses the orientation (s_k < 0). The ring is counter-clockwise from
    // its smallest point (0,-1).
    let mirrored = frame(Point3::new(0.0, 0.0, 1.0), Vec3::X, -Vec3::Y);
    assert_eq!(mirrored.z, -Vec3::Z);
    let result = section(&cube(), mirrored).expect("section");
    assert_eq!(result.regions.len(), 1);
    assert_eq!(result.regions[0].side, InPlaneSide::Above);
    assert_eq!(
        result.regions[0].outer,
        ring(&[(0.0, -1.0), (1.0, -1.0), (1.0, 0.0), (0.0, 0.0)])
    );
}

#[test]
fn a_vertical_frame_drops_another_axis() {
    // Frame x = X, y = Z, so z = X x Z = -Y: the plane y = 0. The front face
    // of the cube (outward -Y) agrees with the frame's +z: Below. Local point
    // of model (x, 0, z) is (x, z). The dominant axis of the normal is 1.
    let vertical = frame(Point3::new(0.0, 0.0, 0.0), Vec3::X, Vec3::Z);
    assert_eq!(vertical.z, -Vec3::Y);
    let result = section(&cube(), vertical).expect("section");
    assert!(result.chains.is_empty());
    assert_eq!(result.regions.len(), 1);
    assert_eq!(result.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        result.regions[0].outer,
        ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    );
}

#[test]
fn a_transversal_plane_still_gives_the_closed_contour() {
    let result = section(&cube(), horizontal(0.5)).expect("section");
    assert!(result.regions.is_empty());
    assert_eq!(result.chains.len(), 1);
    assert!(result.chains[0].is_closed());
    assert_eq!(
        result.chains[0].points,
        ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    );
    assert_eq!(result.evidence.closed_chains, 1);
}

#[test]
fn a_plane_above_the_solid_gives_nothing() {
    let result = section(&cube(), horizontal(2.0)).expect("section");
    assert!(result.chains.is_empty() && result.regions.is_empty());
    assert_eq!(result.evidence.in_plane_triangles, 0);
}

#[test]
fn overlapping_faces_of_one_side_are_refused_by_name() {
    // Boxes [0,2] x [0,1] and [1,3] x [0,1], both to z in [0,1]: the tops
    // overlap in [1,2] x [0,1], both facing +z.
    let mesh = shells(&[
        ((0.0, 2.0), (0.0, 1.0), (0.0, 1.0)),
        ((1.0, 3.0), (0.0, 1.0), (0.0, 1.0)),
    ]);
    assert_unsupported(
        section(&mesh, horizontal(1.0)),
        REFUSED_OVERLAPPING_IN_PLANE_FACES,
    );
}

#[test]
fn faces_of_opposite_sides_that_overlap_are_refused_by_name() {
    // A slab top at z = 1 (+z) under the bottom (-z) of a box standing on it.
    let mesh = shells(&[
        ((0.0, 4.0), (0.0, 4.0), (0.0, 1.0)),
        ((1.0, 2.0), (1.0, 2.0), (1.0, 2.0)),
    ]);
    assert_unsupported(
        section(&mesh, horizontal(1.0)),
        REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES,
    );
}

#[test]
fn faces_of_opposite_sides_that_only_touch_are_not_an_overlap() {
    // A cube top [0,1]^2 facing +z and the bottom of a box [1,2] x [0,1]
    // facing -z meet along the line x = 1 and share no area.
    let mesh = shells(&[
        ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        ((1.0, 2.0), (0.0, 1.0), (1.0, 2.0)),
    ]);
    let result = section(&mesh, horizontal(1.0)).expect("touching is not overlapping");
    assert_eq!(result.regions.len(), 2);
    assert_eq!(result.regions[0].side, InPlaneSide::Below);
    assert_eq!(result.regions[1].side, InPlaneSide::Above);
    assert_eq!(
        result.regions[1].outer,
        ring(&[(1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0)])
    );
}

#[test]
fn a_saddle_vertex_is_refused_by_name() {
    // A bipyramid over the skew ring (1,0,1), (0,1,-1), (-1,0,1), (0,-1,-1)
    // with apexes (0,0,0) and (0,0,-3). The plane z = 0 passes through the
    // upper apex, whose four neighbours alternate above and below the plane:
    // four cut segments meet there.
    let positions = vec![
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(0.0, 1.0, -1.0),
        Point3::new(-1.0, 0.0, 1.0),
        Point3::new(0.0, -1.0, -1.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, -3.0),
    ];
    let mut indices = Vec::new();
    for k in 0..4_u32 {
        let n = (k + 1) % 4;
        indices.extend([4, k, n]);
        indices.extend([5, n, k]);
    }
    let mesh = TriMesh::new(positions, indices);
    assert_unsupported(section(&mesh, horizontal(0.0)), REFUSED_SADDLE_VERTEX);
}

#[test]
fn region_budgets_and_cancellation_are_enforced() {
    let mesh = cube();
    let limits = SectionLimits::new(1_000, 1_000, 1_000, 100);
    for (region_limits, resource) in [
        (
            RegionLimits::new(0, 1_000, 1_000),
            "section in-plane regions",
        ),
        (
            RegionLimits::new(10, 3, 1_000),
            "section in-plane region vertices",
        ),
    ] {
        let request = SectionRequest::new(horizontal(1.0), limits, region_limits);
        assert_eq!(
            ScalarSection::new().section_request(&mesh, &request, &options()),
            Err(GeomError::BudgetExceeded { resource })
        );
    }
    // Two overlapping boxes need pair work; allow none.
    let overlapping = shells(&[
        ((0.0, 2.0), (0.0, 1.0), (0.0, 1.0)),
        ((1.0, 3.0), (0.0, 1.0), (0.0, 1.0)),
    ]);
    let request = SectionRequest::new(horizontal(1.0), limits, RegionLimits::new(10, 100, 0));
    assert_eq!(
        ScalarSection::new().section_request(&overlapping, &request, &options()),
        Err(GeomError::BudgetExceeded {
            resource: "section in-plane pairs"
        })
    );

    let token = CancellationToken::new();
    token.cancel();
    let cancelled = ExecutionOptions::new(Tolerance::METRE).with_cancellation(token);
    assert_eq!(
        ScalarSection::new().section_request(&mesh, &request_for_cube(), &cancelled),
        Err(GeomError::Cancelled)
    );
}

fn request_for_cube() -> SectionRequest {
    request(horizontal(1.0))
}

#[test]
fn the_registry_accepts_the_reference_answer() {
    let mut registry = MeshPlaneSectionRegistry::new();
    registry.register(0, ScalarSection::new());
    let result = registry
        .section_request(&cube(), &request(horizontal(1.0)), &options())
        .expect("registry validates the reference answer");
    assert_eq!(result.regions.len(), 1);
}

#[test]
fn the_legacy_path_still_refuses_an_in_plane_face() {
    let mut registry = MeshPlaneSectionRegistry::new();
    registry.register(0, ScalarSection::new());
    let error = registry
        .section(
            &cube(),
            horizontal(1.0),
            SectionLimits::new(100, 100, 100, 10),
            &options(),
        )
        .expect_err("legacy refuses");
    assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
}

fn assert_unsupported(result: Result<PlaneSection, GeomError>, family: &'static str) {
    match result {
        Err(GeomError::UnsupportedInput {
            operation: Operation::MeshPlaneSection,
            input,
            ..
        }) => assert_eq!(input, family),
        other => panic!("expected a refusal naming {family:?}, got {other:?}"),
    }
}
