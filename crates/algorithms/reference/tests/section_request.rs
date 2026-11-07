//! Section requests: faces in the plane become regions (ADR 0091).
//!
//! Every expectation is derived by hand next to its case.

use axiolid_contracts::{CancellationToken, ExecutionOptions, GeomError, Operation};
use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_dispatch::MeshPlaneSectionRegistry;
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    ChainEnd, CutRegion, CutRegionLimits, InPlaneSide, MeshPlaneSection, PlaneSection,
    RegionLimits, SectionLimits, SectionRequest, SectionRole, REFUSED_CROSSING_CUT_RINGS,
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
    // Of the five edges (four sides and the diagonal) ten pairs, two of them
    // apart by their boxes (opposite sides): eight segment pairs. Each side
    // also meets the triangle that does not carry it: four more.
    assert_eq!(result.evidence.in_plane_pairs, 8 + 4);
    assert_eq!(result.evidence.in_plane_crossings, 0);
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
fn overlapping_faces_of_one_side_merge() {
    // Boxes [0,2] x [0,1] and [1,3] x [0,1], both to z in [0,1]: the tops
    // overlap in [1,2] x [0,1], both facing +z. Their union is [0,3] x [0,1];
    // the diagonals of the two tops cross the other box's edges at x = 1 and
    // x = 2: two crossings.
    let mesh = shells(&[
        ((0.0, 2.0), (0.0, 1.0), (0.0, 1.0)),
        ((1.0, 3.0), (0.0, 1.0), (0.0, 1.0)),
    ]);
    let result = section(&mesh, horizontal(1.0)).expect("merged");
    assert_eq!(result.regions.len(), 1);
    assert_eq!(result.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        result.regions[0].outer,
        ring(&[(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (0.0, 1.0)])
    );
    assert!(result.regions[0].holes.is_empty());
    assert_eq!(result.evidence.in_plane_triangles, 4);
    assert_eq!(result.evidence.in_plane_crossings, 2);
}

#[test]
fn faces_of_opposite_sides_that_overlap_split_into_below_and_both() {
    // A slab top at z = 1 (+z) under the bottom (-z) of a box standing on it:
    // the slab top minus the footprint is Below, with the footprint as a
    // hole, and the footprint is Both. The footprint's diagonal runs along
    // the slab's, so nothing crosses.
    let mesh = shells(&[
        ((0.0, 4.0), (0.0, 4.0), (0.0, 1.0)),
        ((1.0, 2.0), (1.0, 2.0), (1.0, 2.0)),
    ]);
    let result = section(&mesh, horizontal(1.0)).expect("split");
    assert_eq!(result.regions.len(), 2);
    assert_eq!(result.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        result.regions[0].outer,
        ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)])
    );
    assert_eq!(
        result.regions[0].holes,
        [ring(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])]
    );
    assert_eq!(result.regions[1].side, InPlaneSide::Both);
    assert_eq!(
        result.regions[1].outer,
        ring(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)])
    );
    assert_eq!(result.evidence.in_plane_crossings, 0);
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
fn a_sheet_folded_flat_onto_itself_splits_into_three_regions() {
    // (0,0) (1,0) (0,1) face +z; the second triangle b -> a (1,0) (0,0)
    // (0.5,1) winds clockwise, so it faces -z and lies on the same side of
    // the shared edge: the sheet doubles back over itself. The edges (0,1)-(1,0)
    // and (0,0)-(0.5,1) cross at (1/3, 2/3). The first triangle alone is
    // Below (0,0) (1/3,2/3) (0,1), the second alone is Above (1,0) (0.5,1)
    // (1/3,2/3), and where they overlap, (0,0) (1,0) (1/3,2/3), is Both.
    let folded = mesh_of(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 1.0, 0.0],
        ],
        &[0, 1, 2, 1, 0, 3],
    );
    let result = sheet(&folded, horizontal(0.0)).expect("folded sheet");
    let sides: Vec<_> = result.regions.iter().map(|region| region.side).collect();
    assert_eq!(
        sides,
        [InPlaneSide::Below, InPlaneSide::Above, InPlaneSide::Both]
    );
    // The crossing is rounded once: the same point in all three rings.
    let crossing = result.regions[0].outer[1];
    assert!((crossing.x - 1.0 / 3.0).abs() < 1e-15 && (crossing.y - 2.0 / 3.0).abs() < 1e-15);
    assert_eq!(
        result.regions[0].outer,
        [Point2::new(0.0, 0.0), crossing, Point2::new(0.0, 1.0)]
    );
    assert_eq!(
        result.regions[1].outer,
        [crossing, Point2::new(1.0, 0.0), Point2::new(0.5, 1.0)]
    );
    assert_eq!(
        result.regions[2].outer,
        [Point2::new(0.0, 0.0), Point2::new(1.0, 0.0), crossing]
    );
    assert_eq!(result.evidence.in_plane_crossings, 1);
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

// ---- cut regions (ADR 0092) -----------------------------------------------

/// A prism: the `(x, z)` profile, counter-clockwise, extruded along `y` from
/// `0` to `depth`. `triangles` triangulates the profile.
fn extruded(profile: &[(f64, f64)], triangles: &[[u32; 3]], depth: f64) -> TriMesh {
    let count = u32::try_from(profile.len()).expect("small fixture");
    let mut positions: Vec<Point3> = profile
        .iter()
        .map(|&(x, z)| Point3::new(x, 0.0, z))
        .collect();
    positions.extend(profile.iter().map(|&(x, z)| Point3::new(x, depth, z)));
    let mut indices = Vec::new();
    for &[a, b, c] in triangles {
        indices.extend([a, b, c]);
        indices.extend([a + count, c + count, b + count]);
    }
    for from in 0..count {
        let to = (from + 1) % count;
        indices.extend([from, from + count, to + count, from, to + count, to]);
    }
    TriMesh::new(positions, indices)
}

/// The L-shaped step of the issue: a column `x` in `[0, 1]` up to `z = 2`
/// beside a ledge `x` in `[1, 2]` whose top lies at `z = 1`.
fn l_step() -> TriMesh {
    extruded(
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ],
        &[[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]],
        1.0,
    )
}

fn cut_request(frame: Frame3) -> SectionRequest {
    request(frame).with_cut_regions(CutRegionLimits::new(100, 1_000, 1_000))
}

fn registry() -> MeshPlaneSectionRegistry {
    let mut registry = MeshPlaneSectionRegistry::new();
    registry.register(0, ScalarSection::new());
    registry
}

fn is_refused_as_crossing(result: Result<PlaneSection, GeomError>) -> bool {
    matches!(
        result,
        Err(GeomError::UnsupportedInput {
            input: REFUSED_CROSSING_CUT_RINGS,
            ..
        })
    )
}

#[test]
fn the_cut_of_an_l_step_closes_along_the_region_edge() {
    // The plane z = 1 is the top of the ledge: a Below region [1,2] x [0,1].
    // The column is cut across [0,1] x [0,1]; the ledge's jamb edge x = 1
    // bounds it, so the ring runs along the chain and back along that edge.
    let mesh = l_step();
    let plain = registry()
        .section_request(&mesh, &request(horizontal(1.0)), &options())
        .expect("plain request");
    assert_eq!(plain.cut_regions, None);
    assert_eq!(plain.regions.len(), 1);
    assert_eq!(plain.regions[0].side, InPlaneSide::Below);

    let cut = registry()
        .section_request(&mesh, &cut_request(horizontal(1.0)), &options())
        .expect("cut request");
    assert_eq!(cut.chains, plain.chains);
    assert_eq!(cut.regions, plain.regions);
    assert_eq!(
        cut.cut_regions,
        Some(vec![CutRegion::new(
            ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
            vec![],
        )])
    );
    assert_eq!(
        (cut.evidence.cut_regions, cut.evidence.cut_region_vertices),
        (1, 4)
    );
}

#[test]
fn a_cut_inside_a_face_needs_no_region_edge() {
    let section = registry()
        .section_request(&cube(), &cut_request(horizontal(0.5)), &options())
        .expect("cut request");
    assert_eq!(
        section.cut_regions,
        Some(vec![CutRegion::new(
            ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]),
            vec![],
        )])
    );
}

#[test]
fn the_cube_top_in_the_plane_has_no_cut() {
    let section = registry()
        .section_request(&cube(), &cut_request(horizontal(1.0)), &options())
        .expect("cut request");
    assert_eq!(section.regions.len(), 1);
    assert_eq!(section.cut_regions, Some(vec![]));
}

#[test]
fn crossing_shells_are_refused_for_the_cut_only() {
    let plus = shells(&[
        ((0.0, 3.0), (1.0, 2.0), (0.0, 1.0)),
        ((1.0, 2.0), (0.0, 3.0), (0.0, 1.0)),
    ]);
    registry()
        .section_request(&plus, &request(horizontal(0.5)), &options())
        .expect("the plain request does not look at the cut");
    assert!(is_refused_as_crossing(registry().section_request(
        &plus,
        &cut_request(horizontal(0.5)),
        &options()
    )));
}

#[test]
fn overlapping_boxes_with_collinear_edges_are_refused() {
    let overlapping = shells(&[
        ((0.0, 2.0), (0.0, 1.0), (0.0, 1.0)),
        ((1.0, 3.0), (0.0, 1.0), (0.0, 1.0)),
    ]);
    assert!(is_refused_as_crossing(registry().section_request(
        &overlapping,
        &cut_request(horizontal(0.5)),
        &options()
    )));
}

#[test]
fn a_sheet_has_no_cut() {
    let sheet = cut_request(horizontal(0.5)).with_role(SectionRole::Sheet);
    let error = registry()
        .section_request(&cube(), &sheet, &options())
        .expect_err("a sheet has no inside");
    assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
    let error = ScalarSection::new()
        .section_request(&cube(), &sheet, &options())
        .expect_err("the provider refuses it too");
    assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
}

#[test]
fn cut_region_budgets_are_enforced() {
    let neighbours = shells(&[
        ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        ((1.0, 2.0), (0.0, 1.0), (0.0, 1.0)),
    ]);
    for (limits, resource) in [
        (CutRegionLimits::new(1, 1_000, 1_000), "section cut regions"),
        (
            CutRegionLimits::new(10, 7, 1_000),
            "section cut region vertices",
        ),
        (CutRegionLimits::new(10, 1_000, 0), "section cut pairs"),
    ] {
        let request = request(horizontal(0.5)).with_cut_regions(limits);
        assert_eq!(
            ScalarSection::new().section_request(&neighbours, &request, &options()),
            Err(GeomError::BudgetExceeded { resource })
        );
    }
}

// ---- cut regions decided on the source geometry (#377) ----------------------

/// A horizontal plane at height `z`, turned about the vertical by
/// `pi * k / total`. An even `k` is viewed from below, which mirrors the
/// frame's own sense and gives the projection a negative frame sign.
fn turned(k: u32, total: u32, z: f64) -> Frame3 {
    let (s, c) = (std::f64::consts::PI * f64::from(k) / f64::from(total)).sin_cos();
    let y = if k % 2 == 1 {
        Vec3::new(-s, c, 0.0)
    } else {
        Vec3::new(s, -c, 0.0)
    };
    frame(Point3::new(0.0, 0.0, z), Vec3::new(c, s, 0.0), y)
}

fn cut(mesh: &TriMesh, frame: Frame3) -> Result<PlaneSection, GeomError> {
    let request = request(frame).with_cut_regions(CutRegionLimits::new(100, 10_000, 100_000));
    registry().section_request(mesh, &request, &options())
}

#[test]
fn touching_shells_pass_under_every_rotated_frame() {
    // B stands on the top face of A and touches it along y = 0.1 in the
    // source, whatever the frame.
    let mesh = shells(&[
        ((0.0, 2.0), (0.0, 0.1), (0.0, 1.0)),
        ((0.9, 1.1), (0.1, 1.0), (0.0, 1.0)),
    ]);
    for k in 1..=200 {
        let section = cut(&mesh, turned(k, 201, 0.5)).unwrap_or_else(|error| {
            panic!("frame {k} refused the touching shells: {error:?}");
        });
        let regions = section.cut_regions.expect("cut regions");
        assert_eq!(regions.len(), 2, "frame {k}");
        assert!(regions.iter().all(|region| region.holes.is_empty()));
        assert_eq!(
            regions
                .iter()
                .map(|region| region.outer.len())
                .collect::<Vec<_>>(),
            [4, 4],
            "frame {k}"
        );
    }
}

/// A cable ladder: two rails joined by ten rungs that touch them.
fn cable_ladder() -> TriMesh {
    let mut parts = vec![
        ((0.0, 3.0), (0.0, 0.06), (0.0, 0.1)),
        ((0.0, 3.0), (0.54, 0.6), (0.0, 0.1)),
    ];
    for rung in 0..10 {
        let x = 0.15 + 0.3 * f64::from(rung);
        parts.push(((x, x + 0.03), (0.06, 0.54), (0.02, 0.05)));
    }
    shells(&parts)
}

#[test]
fn a_rotated_cable_ladder_cuts_to_its_bars() {
    let mesh = cable_ladder();
    // The plan at rung height: two rails and ten rungs.
    for k in 1..=100 {
        let section = cut(&mesh, turned(k, 101, 0.035))
            .unwrap_or_else(|error| panic!("plan frame {k} refused: {error:?}"));
        assert_eq!(section.cut_regions.expect("cut regions").len(), 12, "{k}");
    }
    // A slightly oblique vertical section through the first rung: the rail
    // ends and the rung. Its crossing nodes are not representable, so this
    // decides on exact crossings.
    for k in 1..=100 {
        let phi = std::f64::consts::FRAC_PI_2 + 0.05 * (2.0 * f64::from(k) / 101.0 - 1.0);
        let (s, c) = phi.sin_cos();
        let vertical = frame(Point3::new(0.165, 0.3, 0.0), Vec3::new(c, s, 0.0), Vec3::Z);
        let section = cut(&mesh, vertical)
            .unwrap_or_else(|error| panic!("oblique frame {k} refused: {error:?}"));
        assert_eq!(section.cut_regions.expect("cut regions").len(), 3, "{k}");
    }
}

#[test]
fn overlapping_shells_stay_refused_under_every_rotated_frame() {
    let sliver = shells(&[
        ((0.0, 2.0), (0.0, 0.1), (0.0, 1.0)),
        ((0.9, 1.1), (0.05, 1.0), (0.0, 1.0)),
    ]);
    let inside = shells(&[
        ((0.0, 2.0), (0.0, 0.1), (0.0, 1.0)),
        ((0.9, 1.1), (0.02, 0.08), (0.0, 1.0)),
    ]);
    for k in 1..=200 {
        for mesh in [&sliver, &inside] {
            assert!(
                is_refused_as_crossing(cut(mesh, turned(k, 201, 0.5))),
                "frame {k}"
            );
        }
    }
}

#[test]
fn a_sub_ulp_overlap_in_the_source_is_refused() {
    // The mesh itself is rotated in binary64, as a consumer's world transform
    // does. The rounded corners of B may fall strictly inside A; that overlap
    // is in the input, and the cut refuses exactly when it is there.
    let (mut refused, mut accepted) = (0, 0);
    for k in 1..=200 {
        let (s, c) = (std::f64::consts::PI * f64::from(k) / 201.0).sin_cos();
        let turn = |x: f64, y: f64| Point2::new(c * x - s * y, s * x + c * y);
        let mut mesh = shells(&[
            ((0.0, 2.0), (0.0, 0.1), (0.0, 1.0)),
            ((0.9, 1.1), (0.1, 1.0), (0.0, 1.0)),
        ]);
        for position in &mut mesh.positions {
            let flat = turn(position.x, position.y);
            *position = Point3::new(flat.x, flat.y, position.z);
        }
        let (a0, a1) = (turn(0.0, 0.1), turn(2.0, 0.1));
        let overlaps = [turn(0.9, 0.1), turn(1.1, 0.1)].into_iter().any(|corner| {
            axiolid_reference::orient2d(a1, a0, corner).sign()
                == Some(axiolid_contracts::Sign::Positive)
        });
        let result = cut(&mesh, horizontal(0.5));
        if overlaps {
            assert!(is_refused_as_crossing(result), "frame {k}");
            refused += 1;
        } else {
            result.unwrap_or_else(|error| panic!("frame {k} refused: {error:?}"));
            accepted += 1;
        }
    }
    assert!(
        refused > 0 && accepted > 0,
        "{refused} refused, {accepted} accepted"
    );
}

// ---- merging in-plane faces on an exact arrangement (#366) -------------------

/// Concatenate meshes as separate shells.
fn combine(parts: &[(Vec<Point3>, Vec<u32>)]) -> TriMesh {
    let mut positions = Vec::new();
    let mut indices = Vec::new();
    for (part_positions, part_indices) in parts {
        let offset = u32::try_from(positions.len()).expect("small fixture");
        positions.extend(part_positions.iter().copied());
        indices.extend(part_indices.iter().map(|index| index + offset));
    }
    TriMesh::new(positions, indices)
}

/// A prism over a convex counter-clockwise polygon between `z.0` and `z.1`,
/// its bottom and top fanned from the first vertex.
fn prism(polygon: &[(f64, f64)], z: Span) -> (Vec<Point3>, Vec<u32>) {
    let count = u32::try_from(polygon.len()).expect("small fixture");
    let positions: Vec<Point3> = [z.0, z.1]
        .into_iter()
        .flat_map(|z| polygon.iter().map(move |&(x, y)| Point3::new(x, y, z)))
        .collect();
    let mut indices = Vec::new();
    for corner in 1..count - 1 {
        indices.extend([0, corner + 1, corner]);
        indices.extend([count, count + corner, count + corner + 1]);
    }
    for from in 0..count {
        let to = (from + 1) % count;
        indices.extend([from, to, count + to, from, count + to, count + from]);
    }
    (positions, indices)
}

/// `mesh` turned about the vertical by `pi * k / 201` in binary64, as a
/// consumer's world transform does.
fn turned_mesh(mesh: &TriMesh, k: u32) -> TriMesh {
    let (s, c) = (std::f64::consts::PI * f64::from(k) / 201.0).sin_cos();
    let mut turned = mesh.clone();
    for position in &mut turned.positions {
        *position = Point3::new(
            c * position.x - s * position.y,
            s * position.x + c * position.y,
            position.z,
        );
    }
    turned
}

fn sizes(result: &PlaneSection, side: InPlaneSide) -> Vec<(usize, Vec<usize>)> {
    result
        .regions
        .iter()
        .filter(|region| region.side == side)
        .map(|region| {
            (
                region.outer.len(),
                region.holes.iter().map(Vec::len).collect(),
            )
        })
        .collect()
}

#[test]
fn a_footprint_inside_one_slab_triangle_is_a_hole_found_by_containment() {
    // The footprint [1,2] x [2.5,3.5] lies above the slab's diagonal, strictly
    // inside its triangle (0,0) (4,4) (0,4): no edge of the slab meets it.
    let mesh = shells(&[
        ((0.0, 4.0), (0.0, 4.0), (0.0, 1.0)),
        ((1.0, 2.0), (2.5, 3.5), (1.0, 2.0)),
    ]);
    let result = section(&mesh, horizontal(1.0)).expect("section");
    assert_eq!(result.regions.len(), 2);
    assert_eq!(result.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        result.regions[0].holes,
        [ring(&[(1.0, 2.5), (1.0, 3.5), (2.0, 3.5), (2.0, 2.5)])]
    );
    assert_eq!(result.regions[1].side, InPlaneSide::Both);
    assert_eq!(result.evidence.in_plane_crossings, 0);
}

#[test]
fn three_overlapping_tops_and_a_box_on_them_merge() {
    // Tops A [0,2]^2, B [1,3] x [0,2] and C [0.5,2.5] x [1,3] at z = 1, all
    // facing +z; the box D [1.5,2] x [0.5,2.5] stands on them. Below is the
    // union of the tops minus D's footprint; Both is the footprint.
    let mesh = shells(&[
        ((0.0, 2.0), (0.0, 2.0), (0.0, 1.0)),
        ((1.0, 3.0), (0.0, 2.0), (0.0, 1.0)),
        ((0.5, 2.5), (1.0, 3.0), (0.0, 1.0)),
        ((1.5, 2.0), (0.5, 2.5), (1.0, 2.0)),
    ]);
    let result = section(&mesh, horizontal(1.0)).expect("section");
    assert_eq!(result.regions.len(), 2);
    assert_eq!(result.regions[0].side, InPlaneSide::Below);
    assert_eq!(
        result.regions[0].outer,
        ring(&[
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 2.0),
            (2.5, 2.0),
            (2.5, 3.0),
            (0.5, 3.0),
            (0.5, 2.0),
            (0.0, 2.0),
        ])
    );
    assert_eq!(
        result.regions[0].holes,
        [ring(&[(1.5, 0.5), (1.5, 2.5), (2.0, 2.5), (2.0, 0.5)])]
    );
    assert_eq!(result.regions[1].side, InPlaneSide::Both);
    assert_eq!(
        result.regions[1].outer,
        ring(&[(1.5, 0.5), (2.0, 0.5), (2.0, 2.5), (1.5, 2.5)])
    );
}

fn diamond_on_a_square() -> TriMesh {
    // The diamond's bottom is split on y = 1, from (1.25,1) to (2.25,1).
    let diamond = [(1.25, 1.0), (1.75, 0.5), (2.25, 1.0), (1.75, 1.5)];
    combine(&[
        cuboid((0.0, 2.0), (0.0, 2.0), (0.0, 1.0)),
        prism(&diamond, (1.0, 2.0)),
    ])
}

fn assert_diamond(result: &PlaneSection, mirrored: bool, tag: &str) {
    // Below is the square minus the diamond (9 points), Above the diamond
    // beyond x = 2 (3 points), Both their overlap (5 points); the edges x = 2
    // meets the diamond's two edges and its split diagonal: three crossings.
    // A frame that views the plane from below swaps Below and Above.
    let (below, above) = if mirrored {
        (InPlaneSide::Above, InPlaneSide::Below)
    } else {
        (InPlaneSide::Below, InPlaneSide::Above)
    };
    assert_eq!(sizes(result, below), [(9, vec![])], "{tag}");
    assert_eq!(sizes(result, above), [(3, vec![])], "{tag}");
    assert_eq!(sizes(result, InPlaneSide::Both), [(5, vec![])], "{tag}");
    assert_eq!(result.evidence.in_plane_crossings, 3, "{tag}");
}

#[test]
fn a_diamond_across_an_edge_splits_in_every_frame() {
    let mesh = diamond_on_a_square();
    for k in 1..=200 {
        let result = section(&mesh, turned(k, 201, 1.0)).expect("section");
        assert_diamond(&result, k % 2 == 0, &format!("frame {k}"));
    }
}

#[test]
fn a_diamond_across_an_edge_splits_in_a_rotated_mesh() {
    let mesh = diamond_on_a_square();
    for k in 1..=200 {
        let result = section(&turned_mesh(&mesh, k), horizontal(1.0)).expect("section");
        assert_diamond(&result, false, &format!("mesh {k}"));
    }
}

#[test]
fn a_box_on_a_slab_in_a_rotated_mesh_keeps_its_hole() {
    let mesh = shells(&[
        ((0.0, 4.0), (0.0, 4.0), (0.0, 1.0)),
        ((1.0, 2.0), (1.0, 2.0), (1.0, 2.0)),
    ]);
    for k in 1..=200 {
        let result = section(&turned_mesh(&mesh, k), horizontal(1.0)).expect("section");
        assert_eq!(
            sizes(&result, InPlaneSide::Below),
            [(4, vec![4])],
            "mesh {k}"
        );
        assert_eq!(sizes(&result, InPlaneSide::Both), [(4, vec![])], "mesh {k}");
    }
}

#[test]
fn a_gap_of_one_ulp_between_two_tops_is_reported_exactly() {
    // The second top starts one ulp after the first one ends: two regions,
    // not one, and the output keeps the two different x values.
    let edge = 1.0_f64.next_up();
    let mesh = shells(&[
        ((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        ((edge, 2.0), (0.0, 1.0), (0.0, 1.0)),
    ]);
    let result = section(&mesh, horizontal(1.0)).expect("section");
    assert_eq!(result.regions.len(), 2);
    assert_eq!(result.regions[0].outer[1], Point2::new(1.0, 0.0));
    assert_eq!(result.regions[1].outer[0], Point2::new(edge, 0.0));
    assert_eq!(result.evidence.in_plane_crossings, 0);
}

#[test]
fn a_collinear_vertex_where_other_boundary_meets_is_kept() {
    // A top [0,2] x [0,1] (Below) and the bottom of a box [1,2] x [1,2]
    // (Above) meet along y = 1 from x = 1 to 2. The vertex (1,1) is straight
    // on the top's boundary, but the box's corner is there too: it stays.
    let mesh = shells(&[
        ((0.0, 2.0), (0.0, 1.0), (0.0, 1.0)),
        ((1.0, 2.0), (1.0, 2.0), (1.0, 2.0)),
    ]);
    let result = section(&mesh, horizontal(1.0)).expect("section");
    assert_eq!(result.regions.len(), 2);
    assert_eq!(
        result.regions[0].outer,
        ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (0.0, 1.0)])
    );
    assert_eq!(
        result.regions[1].outer,
        ring(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)])
    );
}

// ---- `Both` regions are cut regions (#366) ------------------------------------

fn parts_of(mesh: &TriMesh) -> (Vec<Point3>, Vec<u32>) {
    (mesh.positions.clone(), mesh.indices.clone())
}

#[test]
fn the_cut_of_an_l_step_with_a_block_in_its_notch_has_two_regions_in_every_frame() {
    // The block [1,2] x [0,1] x [1,2] sits on the ledge, its own shell.
    let block = cuboid((1.0, 2.0), (0.0, 1.0), (1.0, 2.0));
    let mesh = combine(&[parts_of(&l_step()), block]);
    for k in 1..=200 {
        let section = cut(&mesh, turned(k, 201, 1.0))
            .unwrap_or_else(|error| panic!("frame {k} refused: {error:?}"));
        assert_eq!(section.cut_regions.expect("cut regions").len(), 2, "{k}");
    }
}

#[test]
fn stacked_boxes_with_crossing_diagonals_cut_to_their_footprint() {
    // The lower top's diagonal runs (0,0)-(1,1), the upper bottom's runs
    // (1,0)-(0,1): they cross at the centre, a binary64 point.
    let flipped = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)];
    let mesh = combine(&[
        cuboid((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)),
        prism(&flipped, (1.0, 2.0)),
    ]);
    let section = cut(&mesh, horizontal(1.0)).expect("section");
    assert_eq!(sizes(&section, InPlaneSide::Both), [(4, vec![])]);
    assert_eq!(section.regions.len(), 1);
    assert_eq!(section.evidence.in_plane_crossings, 1);
    assert_eq!(
        section.cut_regions.expect("cut regions")[0].outer,
        section.regions[0].outer
    );
}

#[test]
fn the_cut_of_a_diamond_across_an_edge_is_its_both_region_bit_for_bit() {
    let mesh = diamond_on_a_square();
    for k in 1..=200 {
        let section = cut(&mesh, turned(k, 201, 1.0))
            .unwrap_or_else(|error| panic!("frame {k} refused: {error:?}"));
        let both: Vec<_> = section
            .regions
            .iter()
            .filter(|region| region.side == InPlaneSide::Both)
            .collect();
        let cut_regions = section.cut_regions.expect("cut regions");
        assert_eq!(both.len(), 1, "{k}");
        assert_eq!(cut_regions.len(), 1, "{k}");
        assert_eq!(cut_regions[0].outer, both[0].outer, "{k}");
        assert_eq!(cut_regions[0].holes, both[0].holes, "{k}");
    }
}

/// A slab `[0,3]^2 x [0,1]` with a triangular pocket A = (0,0), F = (2,1),
/// E = (1,2) sunk to `z = 0.5`; A is also a corner of the slab, so the top
/// is a square with a notch that touches its corner.
fn slab_with_corner_pocket() -> (Vec<Point3>, Vec<u32>) {
    let (mut positions, mut indices) = cuboid((0.0, 3.0), (0.0, 3.0), (0.0, 1.0));
    // Drop the top face; the corners 4 = A, 5 = B, 6 = C, 7 = D stay.
    indices.drain(6..12);
    let mut add = |x: f64, y: f64, z: f64| {
        positions.push(Point3::new(x, y, z));
        u32::try_from(positions.len() - 1).expect("small fixture")
    };
    let (e, f) = (add(1.0, 2.0, 1.0), add(2.0, 1.0, 1.0));
    let (a_low, e_low, f_low) = (add(0.0, 0.0, 0.5), add(1.0, 2.0, 0.5), add(2.0, 1.0, 0.5));
    let (a, b, c, d) = (4, 5, 6, 7);
    for [p, q, r] in [[a, b, f], [b, c, f], [c, e, f], [c, d, e], [d, a, e]] {
        indices.extend([p, q, r]);
    }
    // The walls face into the pocket, to the left of each rim edge, and the
    // floor faces up.
    for (p, q, p_low, q_low) in [
        (a, f, a_low, f_low),
        (f, e, f_low, e_low),
        (e, a, e_low, a_low),
    ] {
        indices.extend([p, q, q_low, p, q_low, p_low]);
    }
    indices.extend([a_low, f_low, e_low]);
    (positions, indices)
}

#[test]
fn a_both_ring_that_touches_itself_is_refused_for_the_cut_only() {
    // A lid [0,3]^2 x [1,2] over the slab: Both is the slab top, whose ring
    // visits the corner (0,0) twice; Above is the pocket's opening.
    let mesh = combine(&[
        slab_with_corner_pocket(),
        cuboid((0.0, 3.0), (0.0, 3.0), (1.0, 2.0)),
    ]);
    let plain = registry()
        .section_request(&mesh, &request(horizontal(1.0)), &options())
        .expect("the plain answer");
    let both: Vec<_> = plain
        .regions
        .iter()
        .filter(|region| region.side == InPlaneSide::Both)
        .collect();
    assert_eq!(both.len(), 1);
    assert_eq!(both[0].outer.len(), 7);
    let corner = Point2::new(0.0, 0.0);
    assert_eq!(both[0].outer.iter().filter(|&&p| p == corner).count(), 2);
    let refused = registry().section_request(&mesh, &cut_request(horizontal(1.0)), &options());
    assert!(matches!(
        refused,
        Err(GeomError::UnsupportedInput {
            input: axiolid_mesh_section_contract::REFUSED_PINCHED_CUT_REGION,
            ..
        })
    ));
}

#[test]
fn a_cut_side_covered_by_in_plane_faces_is_refused_as_crossing() {
    // The block [0.5,2] x [0,1] x [1,2] overlaps the column of the L-step in
    // volume: its bottom covers the plane beside the ledge's jamb, where the
    // cut lies. The plain answer sees nothing wrong.
    let block = cuboid((0.5, 2.0), (0.0, 1.0), (1.0, 2.0));
    let mesh = combine(&[parts_of(&l_step()), block]);
    registry()
        .section_request(&mesh, &request(horizontal(1.0)), &options())
        .expect("the plain answer");
    assert!(is_refused_as_crossing(registry().section_request(
        &mesh,
        &cut_request(horizontal(1.0)),
        &options()
    )));
}

/// A flat square annulus at height `z`, facing up: `outer` minus `inner` as
/// eight triangles, one trapezoid per side.
fn annulus(outer: Span, inner: Span, z: f64) -> (Vec<Point3>, Vec<u32>) {
    let corners = |(low, high): Span| [(low, low), (high, low), (high, high), (low, high)];
    let (o, i) = (corners(outer), corners(inner));
    let mut positions: Vec<Point3> = o.iter().map(|&(x, y)| Point3::new(x, y, z)).collect();
    positions.extend(i.iter().map(|&(x, y)| Point3::new(x, y, z)));
    let mut indices = Vec::new();
    for k in 0..4_u32 {
        let n = (k + 1) % 4;
        indices.extend([k, n, 4 + n, k, 4 + n, 4 + k]);
    }
    (positions, indices)
}

#[test]
fn an_annulus_inside_the_hole_of_another_keeps_its_own_hole() {
    // The hole [3,4]^2 lies in the inner annulus and in the outer one's
    // outline; it belongs to the innermost outer ring, [2,5]^2.
    let flat = combine(&[
        annulus((0.0, 7.0), (1.0, 6.0), 0.0),
        annulus((2.0, 5.0), (3.0, 4.0), 0.0),
    ]);
    let result = sheet(&flat, horizontal(0.0)).expect("sheet");
    assert_eq!(result.regions.len(), 2);
    assert_eq!(
        result.regions[0].outer,
        ring(&[(0.0, 0.0), (7.0, 0.0), (7.0, 7.0), (0.0, 7.0)])
    );
    assert_eq!(
        result.regions[0].holes,
        [ring(&[(1.0, 1.0), (1.0, 6.0), (6.0, 6.0), (6.0, 1.0)])]
    );
    assert_eq!(
        result.regions[1].outer,
        ring(&[(2.0, 2.0), (5.0, 2.0), (5.0, 5.0), (2.0, 5.0)])
    );
    assert_eq!(
        result.regions[1].holes,
        [ring(&[(3.0, 3.0), (3.0, 4.0), (4.0, 4.0), (4.0, 3.0)])]
    );
}

// ---- guards of the merge (#366) ----------------------------------------------

#[test]
fn a_crossing_that_rounds_onto_a_neighbouring_vertex_is_refused_as_collapsed() {
    // The edge (0,1)-(2^60,2) crosses the vertical edge x = 1 at
    // (1, 1 + 2^-60), which is exact in the arrangement and rounds to (1,1),
    // the vertex next to it on the same ring.
    let far = 2.0_f64.powi(60);
    let flat = mesh_of(
        &[
            [0.0, 1.0, 0.0],
            [far, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 1.0, 0.0],
            [5.0, 3.0, 0.0],
            [1.0, 3.0, 0.0],
        ],
        &[0, 1, 2, 3, 4, 5],
    );
    match sheet(&flat, horizontal(0.0)) {
        Err(GeomError::Degenerate(message)) => {
            assert_eq!(message, "in-plane region collapsed under projection");
        }
        other => panic!("expected a collapsed region, got {other:?}"),
    }
}

#[test]
fn a_piece_held_by_a_triangle_whose_box_touches_a_crossing_is_covered() {
    // The slanted edge (0,0)-(3,1) crosses the vertical edge x = 1 at
    // (1, 1/3), which is not a binary64 point. The triangle B starts at
    // v = fl(1/3), below the crossing, and holds the piece above it, so its
    // box must not be pruned by the crossing's enclosure. With B the piece is
    // interior; the outline has no vertex there.
    let low = 1.0_f64 / 3.0;
    let flat = mesh_of(
        &[
            [0.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [3.0, -1.0, 0.0],
            [1.0, 0.9, 0.0],
            [0.5, low, 0.0],
            [2.0, low, 0.0],
            [1.0, 2.0, 0.0],
        ],
        &[0, 1, 2, 3, 4, 5, 6, 7, 8],
    );
    let result = sheet(&flat, horizontal(0.0)).expect("sheet");
    assert_eq!(result.evidence.in_plane_crossings, 8);
    assert_eq!(result.regions.len(), 1, "{:?}", result.regions);
    assert!(result.regions[0].holes.is_empty());
}

#[test]
fn the_hole_containment_pairs_count_against_the_budget() {
    let flat = combine(&[
        annulus((0.0, 7.0), (1.0, 6.0), 0.0),
        annulus((2.0, 5.0), (3.0, 4.0), 0.0),
    ]);
    let full = sheet(&flat, horizontal(0.0)).expect("sheet");
    let pairs = full.evidence.in_plane_pairs;
    assert_eq!(pairs, 133);
    let mut limited = request(horizontal(0.0)).with_role(SectionRole::Sheet);
    limited.region_limits = RegionLimits::new(100, 1_000, pairs);
    assert!(ScalarSection::new()
        .section_request(&flat, &limited, &options())
        .is_ok());
    limited.region_limits = RegionLimits::new(100, 1_000, pairs - 1);
    assert!(matches!(
        ScalarSection::new().section_request(&flat, &limited, &options()),
        Err(GeomError::BudgetExceeded {
            resource: "section in-plane pairs"
        })
    ));
}

// ---- cut boundaries that pinch at a vertex (#373) ----------------------------

/// A mesh whose vertices are shared by exact position.
#[derive(Default)]
struct Welded {
    positions: Vec<Point3>,
    indices: Vec<u32>,
}

impl Welded {
    fn vertex(&mut self, [x, y, z]: [f64; 3]) -> u32 {
        let point = Point3::new(x, y, z);
        let index = self
            .positions
            .iter()
            .position(|&known| known == point)
            .unwrap_or_else(|| {
                self.positions.push(point);
                self.positions.len() - 1
            });
        u32::try_from(index).expect("small fixture")
    }

    /// Triangles fanned from the first point of a convex face,
    /// counter-clockwise seen from outside.
    fn fan(&mut self, face: &[[f64; 3]]) {
        for k in 1..face.len() - 1 {
            for corner in [face[0], face[k], face[k + 1]] {
                let index = self.vertex(corner);
                self.indices.push(index);
            }
        }
    }

    /// A pier from z = 1 to z = 2 over a convex counter-clockwise footprint,
    /// its top shifted by `lean`.
    fn leaning_pier(&mut self, footprint: &[(f64, f64)], lean: (f64, f64)) {
        let foot = |(x, y): (f64, f64)| [x, y, 1.0];
        let top = |(x, y): (f64, f64)| [x + lean.0, y + lean.1, 2.0];
        for k in 0..footprint.len() {
            let (p, q) = (footprint[k], footprint[(k + 1) % footprint.len()]);
            self.fan(&[foot(p), foot(q), top(q), top(p)]);
        }
        let roof: Vec<[f64; 3]> = footprint.iter().map(|&p| top(p)).collect();
        self.fan(&roof);
    }

    fn finish(self) -> TriMesh {
        TriMesh::new(self.positions, self.indices)
    }
}

fn at(p: (f64, f64), z: f64) -> [f64; 3] {
    [p.0, p.1, z]
}

/// A plinth [0,2]^2 x [0,1] with piers on [0,1]^2 and [1,2]^2 that lean
/// apart: at z = 1 they meet at (1,1) only. The plinth top is the other two
/// unit squares; each side carries the midpoint of its top edge.
fn piers_meeting_at_a_corner() -> TriMesh {
    let mut mesh = Welded::default();
    mesh.fan(&[
        [0.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [2.0, 2.0, 0.0],
        [2.0, 0.0, 0.0],
    ]);
    for [p, m, q] in [
        [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)],
        [(2.0, 0.0), (2.0, 1.0), (2.0, 2.0)],
        [(2.0, 2.0), (1.0, 2.0), (0.0, 2.0)],
        [(0.0, 2.0), (0.0, 1.0), (0.0, 0.0)],
    ] {
        mesh.fan(&[at(p, 0.0), at(q, 0.0), at(q, 1.0), at(m, 1.0), at(p, 1.0)]);
    }
    mesh.fan(&[
        [1.0, 0.0, 1.0],
        [2.0, 0.0, 1.0],
        [2.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ]);
    mesh.fan(&[
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
        [1.0, 2.0, 1.0],
        [0.0, 2.0, 1.0],
    ]);
    mesh.leaning_pier(
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        (-0.5, -0.5),
    );
    mesh.leaning_pier(
        &[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)],
        (0.5, 0.5),
    );
    mesh.finish()
}

/// A hexagonal plinth around the origin whose sectors from (2,0), (0,2) and
/// (-2,-2) carry piers leaning away from it; the other three are its top.
fn three_piers_at_one_vertex() -> TriMesh {
    let outline = [
        (2.0, 0.0),
        (2.0, 2.0),
        (0.0, 2.0),
        (-2.0, 0.0),
        (-2.0, -2.0),
        (0.0, -2.0),
    ];
    let leans = [
        Some((0.5, 0.25)),
        None,
        Some((-0.25, 0.5)),
        None,
        Some((-0.25, -0.5)),
        None,
    ];
    let mut mesh = Welded::default();
    let bottom: Vec<[f64; 3]> = outline.iter().rev().map(|&p| at(p, 0.0)).collect();
    mesh.fan(&bottom);
    for k in 0..6 {
        let (p, q) = (outline[k], outline[(k + 1) % 6]);
        mesh.fan(&[at(p, 0.0), at(q, 0.0), at(q, 1.0), at(p, 1.0)]);
        match leans[k] {
            Some(lean) => mesh.leaning_pier(&[(0.0, 0.0), p, q], lean),
            None => mesh.fan(&[at((0.0, 0.0), 1.0), at(p, 1.0), at(q, 1.0)]),
        }
    }
    mesh.finish()
}

/// A block [0,3]^2 x [0,2] with a pyramid void over the floor (1.5,0)
/// (2.5,2) (0.5,2) at z = 1, apex (1.5,1.25,1.5); the wall y = 0 is fanned
/// around (1.5,0,1), where the void touches it.
fn block_with_a_void_touching_its_wall() -> TriMesh {
    let (v, q, r, apex) = (
        [1.5, 0.0, 1.0],
        [2.5, 2.0, 1.0],
        [0.5, 2.0, 1.0],
        [1.5, 1.25, 1.5],
    );
    let wall = [
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [3.0, 0.0, 2.0],
        [0.0, 0.0, 2.0],
    ];
    let mut mesh = Welded::default();
    for k in 0..4 {
        mesh.fan(&[wall[k], wall[(k + 1) % 4], v]);
    }
    mesh.fan(&[
        [3.0, 0.0, 0.0],
        [3.0, 3.0, 0.0],
        [3.0, 3.0, 2.0],
        [3.0, 0.0, 2.0],
    ]);
    mesh.fan(&[
        [3.0, 3.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 3.0, 2.0],
        [3.0, 3.0, 2.0],
    ]);
    mesh.fan(&[
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 2.0],
        [0.0, 3.0, 2.0],
    ]);
    mesh.fan(&[
        [0.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [3.0, 3.0, 0.0],
        [3.0, 0.0, 0.0],
    ]);
    mesh.fan(&[
        [0.0, 0.0, 2.0],
        [3.0, 0.0, 2.0],
        [3.0, 3.0, 2.0],
        [0.0, 3.0, 2.0],
    ]);
    mesh.fan(&[v, q, r]);
    mesh.fan(&[v, apex, q]);
    mesh.fan(&[q, apex, r]);
    mesh.fan(&[r, apex, v]);
    mesh.finish()
}

/// Outer ring and hole lengths of each cut region.
type Shape = Vec<(usize, Vec<usize>)>;

fn cut_shape(section: &PlaneSection) -> Shape {
    section
        .cut_regions
        .as_ref()
        .expect("cut regions")
        .iter()
        .map(|region| {
            (
                region.outer.len(),
                region.holes.iter().map(Vec::len).collect(),
            )
        })
        .collect()
}

#[test]
fn pinched_cuts_split_at_their_vertex_in_every_frame_and_rotated_mesh() {
    // The exact values in the plan frame are conformance fixtures; here the
    // pairing must not depend on the frame or on a binary64 rotation.
    let cases: [(&str, TriMesh, Shape); 3] = [
        (
            "piers",
            piers_meeting_at_a_corner(),
            vec![(4, vec![]), (4, vec![])],
        ),
        (
            "three piers",
            three_piers_at_one_vertex(),
            vec![(3, vec![]), (3, vec![]), (3, vec![])],
        ),
        (
            "void",
            block_with_a_void_touching_its_wall(),
            vec![(5, vec![3])],
        ),
    ];
    for (name, mesh, shape) in cases {
        for k in 1..=200 {
            let section = cut(&mesh, turned(k, 201, 1.0))
                .unwrap_or_else(|error| panic!("{name}, frame {k}: {error:?}"));
            assert_eq!(cut_shape(&section), shape, "{name}, frame {k}");
            let section = cut(&turned_mesh(&mesh, k), horizontal(1.0))
                .unwrap_or_else(|error| panic!("{name}, mesh {k}: {error:?}"));
            assert_eq!(cut_shape(&section), shape, "{name}, mesh {k}");
        }
    }
}
