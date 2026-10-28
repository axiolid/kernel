//! Section requests: faces in the plane become regions (ADR 0091).
//!
//! Every expectation is derived by hand next to its case.

use axiolid_contracts::{CancellationToken, ExecutionOptions, GeomError, Operation};
use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_dispatch::MeshPlaneSectionRegistry;
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    ChainEnd, InPlaneSide, MeshPlaneSection, PlaneSection, RegionLimits, SectionLimits,
    SectionRequest, SectionRole, REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES,
    REFUSED_OVERLAPPING_IN_PLANE_FACES, REFUSED_SADDLE_VERTEX,
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
    // The two triangles overlap in bounding box: one candidate pair, counted
    // although they share an edge and need no overlap test.
    assert_eq!(result.evidence.in_plane_pairs, 1);
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
fn a_mirrored_frame_still_gets_a_counter_clockwise_contour() {
    // Frame x = X, y = -Y at z = 0.5: local (x, y) is model (x, -y). The
    // contour is oriented exactly from the triangle windings, not from a float
    // area: counter-clockwise from its smallest corner (0,-1), the face
    // diagonals' crossing points dropped as exactly coplanar.
    let mirrored = frame(Point3::new(0.0, 0.0, 0.5), Vec3::X, -Vec3::Y);
    let result = section(&cube(), mirrored).expect("section");
    assert_eq!(result.chains.len(), 1);
    assert_eq!(
        result.chains[0].points,
        ring(&[(0.0, -1.0), (1.0, -1.0), (1.0, 0.0), (0.0, 0.0)])
    );
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
fn registration_runs_the_region_suite_for_a_provider_that_claims_it() {
    // Delegates everything to the reference provider but answers wrongly.
    #[derive(Debug)]
    struct Wrong;
    impl axiolid_contracts::Backend for Wrong {
        fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
            ScalarSection::new().descriptor()
        }
    }
    impl MeshPlaneSection for Wrong {
        fn section(
            &self,
            mesh: &TriMesh,
            frame: Frame3,
            limits: SectionLimits,
            options: &ExecutionOptions,
        ) -> Result<axiolid_mesh_section_contract::SectionOutcome, GeomError> {
            ScalarSection::new().section(mesh, frame, limits, options)
        }
        fn section_request(
            &self,
            mesh: &TriMesh,
            request: &SectionRequest,
            options: &ExecutionOptions,
        ) -> Result<PlaneSection, GeomError> {
            let mut answer = ScalarSection::new().section_request(mesh, request, options)?;
            answer.regions.clear();
            Ok(answer)
        }
    }
    // Delegates only the legacy method: the default section_request refuses.
    #[derive(Debug)]
    struct Legacy;
    impl axiolid_contracts::Backend for Legacy {
        fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
            ScalarSection::new().descriptor()
        }
    }
    impl MeshPlaneSection for Legacy {
        fn section(
            &self,
            mesh: &TriMesh,
            frame: Frame3,
            limits: SectionLimits,
            options: &ExecutionOptions,
        ) -> Result<axiolid_mesh_section_contract::SectionOutcome, GeomError> {
            ScalarSection::new().section(mesh, frame, limits, options)
        }
    }

    let mut registry = MeshPlaneSectionRegistry::new();
    registry
        .register_request_conformant(0, ScalarSection::new())
        .expect("the reference provider is conformant");
    registry
        .register_request_conformant(0, Legacy)
        .expect("a provider that does not claim the method passes");
    let report = registry
        .register_request_conformant(0, Wrong)
        .expect_err("a claimed wrong answer is rejected");
    assert!(report.claimed && !report.is_success());
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

// ---- sheets: open meshes (ADR 0091) --------------------------------------

fn sheet(mesh: &TriMesh, frame: Frame3) -> Result<PlaneSection, GeomError> {
    let request = request(frame).with_role(SectionRole::Sheet);
    ScalarSection::new().section_request(mesh, &request, &options())
}

fn mesh_of(points: &[[f64; 3]], indices: &[u32]) -> TriMesh {
    TriMesh::new(
        points
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
        indices.to_vec(),
    )
}

/// A floor triangle `(0,0,0) (2,0,0) (0,2,0)` in the plane z = 0, and a wall
/// triangle `(2,0,0) (3,-1,1) (3,1,-1)` that meets it in the vertex
/// `(2,0,0)` only. The wall crosses the plane along the line from `(2,0)` to
/// the middle `(3,0)` of its edge `(3,-1,1) (3,1,-1)`. The floor is an
/// in-plane triangle, so the vertex end meets a region; the other end is on
/// an edge used by the wall alone, a boundary edge of the sheet. Six boundary
/// edges: three of each triangle.
fn floor_and_wall() -> TriMesh {
    mesh_of(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [3.0, -1.0, 1.0],
            [3.0, 1.0, -1.0],
        ],
        &[0, 1, 2, 1, 3, 4],
    )
}

#[test]
fn a_chain_ends_on_a_region_at_one_end_and_on_the_mesh_boundary_at_the_other() {
    let found = sheet(&floor_and_wall(), horizontal(0.0)).expect("section");
    assert_eq!(found.chains.len(), 1);
    assert_eq!(found.chains[0].points, ring(&[(2.0, 0.0), (3.0, 0.0)]));
    assert_eq!(
        found.chains[0].ends,
        Some([ChainEnd::InPlaneRegion, ChainEnd::MeshBoundary])
    );
    assert_eq!(found.regions.len(), 1);
    assert_eq!(found.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        found.regions[0].outer,
        ring(&[(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)])
    );
    let evidence = found.evidence;
    assert_eq!(
        (evidence.role, evidence.boundary_edges),
        (SectionRole::Sheet, 6)
    );
    assert_eq!(
        (evidence.in_plane_triangles, evidence.touching_edges),
        (1, 0)
    );

    // Mirrored in x, the frame's z is -Z: the chain runs from (-3,0) to
    // (-2,0), so its ends swap with the points; the floor is now Above.
    let mirrored = frame(Point3::new(0.0, 0.0, 0.0), -Vec3::X, Vec3::Y);
    let found = sheet(&floor_and_wall(), mirrored).expect("section");
    assert_eq!(found.chains[0].points, ring(&[(-3.0, 0.0), (-2.0, 0.0)]));
    assert_eq!(
        found.chains[0].ends,
        Some([ChainEnd::MeshBoundary, ChainEnd::InPlaneRegion])
    );
    assert_eq!(found.regions[0].side, InPlaneSide::Above);
}

#[test]
fn a_solid_cannot_end_a_chain_on_the_mesh_boundary() {
    // A wall quad in the plane y = 0, cut by x = 1: the line ends on its
    // boundary edges. A sheet reports it, a solid cannot have it.
    let wall = mesh_of(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
        ],
        &[0, 1, 2, 0, 2, 3],
    );
    let vertical = frame(Point3::new(1.0, 0.0, 0.0), Vec3::Y, Vec3::Z);
    let found = sheet(&wall, vertical).expect("a sheet may end on the boundary");
    assert_eq!(found.chains.len(), 1);
    assert_eq!(
        found.chains[0].ends,
        Some([ChainEnd::MeshBoundary, ChainEnd::MeshBoundary])
    );
    let error = section(&wall, vertical).expect_err("a solid is closed");
    assert!(
        matches!(error, GeomError::BackendContractViolation { .. }),
        "{error:?}"
    );
}

#[test]
fn an_edge_in_the_plane_with_one_off_plane_triangle_only_touches() {
    // A wall x in [0,2], z in [0,1] at y = 0: its bottom edge lies in the
    // plane z = 0 and has one incident triangle, off the plane.
    let wall = mesh_of(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
        ],
        &[0, 1, 2, 0, 2, 3],
    );
    let found = sheet(&wall, horizontal(0.0)).expect("section");
    assert!(found.chains.is_empty() && found.regions.is_empty());
    assert_eq!(found.evidence.touching_edges, 1);
    assert_eq!(found.evidence.boundary_edges, 4);
}

#[test]
fn a_closed_chain_of_a_sheet_is_counter_clockwise_whatever_its_facing() {
    // The upper half of an octahedron: apex (0,0,1) over the ring (1,0,0),
    // (0,1,0), (-1,0,0), (0,-1,0). The plane z = 0.5 cuts the four edges to
    // the apex at their middles: a diamond (-0.5,0) (0,-0.5) (0.5,0) (0,0.5)
    // counter-clockwise from its smallest point. Four boundary edges.
    let ring_points = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ];
    for outward in [true, false] {
        let mut points = ring_points.to_vec();
        points.push([0.0, 0.0, 1.0]);
        let mut indices = Vec::new();
        for k in 0..4_u32 {
            let n = (k + 1) % 4;
            indices.extend(if outward { [k, n, 4] } else { [n, k, 4] });
        }
        let found = sheet(&mesh_of(&points, &indices), horizontal(0.5)).expect("section");
        assert_eq!(found.chains.len(), 1, "outward {outward}");
        assert_eq!(found.chains[0].ends, None);
        assert_eq!(
            found.chains[0].points,
            ring(&[(-0.5, 0.0), (0.0, -0.5), (0.5, 0.0), (0.0, 0.5)]),
            "outward {outward}"
        );
        assert_eq!(found.evidence.boundary_edges, 4);
    }
}

#[test]
fn a_sheet_folded_flat_onto_itself_is_refused_as_faces_of_both_sides() {
    // (0,0) (1,0) (0,1) face +z; the second triangle b -> a (1,0) (0,0)
    // (0.5,1) winds clockwise, so it faces -z and lies on the same side of
    // the shared edge: the sheet doubles back over itself. The pair shares an
    // edge and must still be tested.
    let folded = mesh_of(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 1.0, 0.0],
        ],
        &[0, 1, 2, 1, 0, 3],
    );
    assert_unsupported(
        sheet(&folded, horizontal(0.0)),
        REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES,
    );
}

#[test]
fn a_registry_refuses_a_sheet_with_inconsistent_winding() {
    let mut registry = MeshPlaneSectionRegistry::new();
    registry.register(0, ScalarSection::new());
    // Both triangles run 0 -> 1.
    let mesh = mesh_of(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
        ],
        &[0, 1, 2, 0, 1, 3],
    );
    // A plane above the sheet: the provider alone would answer with nothing.
    let request = request(horizontal(1.0)).with_role(SectionRole::Sheet);
    assert!(matches!(
        registry.section_request(&mesh, &request, &options()),
        Err(GeomError::NotManifold(_))
    ));
}
