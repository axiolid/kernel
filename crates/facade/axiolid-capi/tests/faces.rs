use std::ptr::null;

use axiolid::mesh_hidden_lines::VisibleFaceEvidence;
use axiolid_capi::*;

#[path = "faces/seed9.rs"]
mod seed9;

fn context_with(config: AxiolidContextConfig) -> AxiolidContextHandle {
    unsafe {
        let mut handle = AxiolidContextHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_context_create(&config, &mut handle),
            AxiolidStatus::Ok
        );
        handle
    }
}

fn context() -> AxiolidContextHandle {
    context_with(AxiolidContextConfig::default())
}

fn import(context: AxiolidContextHandle, positions: &[f64], indices: &[u32]) -> AxiolidMeshHandle {
    unsafe {
        let mut mesh = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_mesh_import(
                context,
                positions.as_ptr(),
                positions.len() / 3,
                indices.as_ptr(),
                indices.len() / 3,
                &mut mesh,
            ),
            AxiolidStatus::Ok
        );
        mesh
    }
}

fn cube(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    import(
        context,
        &[
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0,
            2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0,
        ],
        &[
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ],
    )
}

fn triangle(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    import(
        context,
        &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        &[0, 1, 2],
    )
}

/// The unit square as two triangles sharing the diagonal 0-2.
fn square(context: AxiolidContextHandle) -> AxiolidMeshHandle {
    import(
        context,
        &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0],
        &[0, 1, 2, 0, 2, 3],
    )
}

fn entry(mesh: AxiolidMeshHandle, id: u64, role: i32) -> AxiolidViewMesh {
    AxiolidViewMesh {
        mesh,
        id,
        groups: null(),
        group_count: 0,
        role,
        group_mode: AXIOLID_GROUP_MODE_SINGLE,
    }
}

fn plan() -> AxiolidView {
    AxiolidView {
        frame: AxiolidViewFrame {
            origin: [0.0; 3],
            x: [1.0, 0.0, 0.0],
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        },
        near: 0.0,
        far: 0.0,
        projection: AXIOLID_PROJECTION_ORTHOGRAPHIC,
        has_near: 0,
        has_far: 0,
        reserved: [0; 2],
    }
}

fn limits() -> AxiolidVisibleFaceLimits {
    AxiolidVisibleFaceLimits {
        max_meshes: 8,
        max_source_vertices: 1000,
        max_source_triangles: 1000,
        max_candidate_pairs: 100_000,
        max_regions: 100,
        max_region_vertices: 1000,
    }
}

const TOLERANCE: AxiolidTolerance = AxiolidTolerance {
    linear: 1.0e-9,
    angular: 1.0e-12,
};

/// A square sheet `[min, max]^2` at height `z`, two triangles.
fn sheet(context: AxiolidContextHandle, min: f64, max: f64, z: f64) -> AxiolidMeshHandle {
    import(
        context,
        &[min, min, z, max, min, z, max, max, z, min, max, z],
        &[0, 1, 2, 0, 2, 3],
    )
}

type Change = Box<dyn Fn(&mut Call)>;

struct Call {
    meshes: Vec<AxiolidViewMesh>,
    view: AxiolidView,
    region: Option<AxiolidImageBox>,
    limits: AxiolidVisibleFaceLimits,
    tolerance: AxiolidTolerance,
}

impl Call {
    fn new(meshes: Vec<AxiolidViewMesh>) -> Self {
        Self {
            meshes,
            view: plan(),
            region: None,
            limits: limits(),
            tolerance: TOLERANCE,
        }
    }

    fn run(&self, context: AxiolidContextHandle) -> (AxiolidStatus, AxiolidResultHandle) {
        let mut result = AxiolidResultHandle::INVALID;
        let region = self.region.as_ref().map_or(null(), |region| region);
        let status = unsafe {
            axiolid_v0_4_visible_faces(
                context,
                self.meshes.as_ptr(),
                self.meshes.len(),
                &self.view,
                region,
                &self.limits,
                self.tolerance,
                &mut result,
            )
        };
        (status, result)
    }
}

fn counts(context: AxiolidContextHandle, result: AxiolidResultHandle) -> AxiolidVisibleFaceCounts {
    let mut counts = AxiolidVisibleFaceCounts {
        regions: 99,
        rings: 99,
        points: 99,
        triangles: 99,
    };
    assert_eq!(
        unsafe { axiolid_v0_4_visible_faces_counts(context, result, &mut counts) },
        AxiolidStatus::Ok
    );
    counts
}

const EMPTY_REGION: AxiolidVisibleRegion = AxiolidVisibleRegion {
    mesh_id: 0,
    first_ring: 0,
    ring_count: 0,
    first_triangle: 0,
    triangle_count: 0,
    group: 0,
    source: AxiolidFaceSource::Unknown,
};

/// The four tables of a result, read with exactly the reported counts.
struct Tables {
    regions: Vec<AxiolidVisibleRegion>,
    rings: Vec<AxiolidRegionRing>,
    points: Vec<[f64; 2]>,
    triangles: Vec<u32>,
}

fn read_tables(context: AxiolidContextHandle, result: AxiolidResultHandle) -> Tables {
    let need = counts(context, result);
    let mut regions = vec![EMPTY_REGION; need.regions];
    let mut rings = vec![
        AxiolidRegionRing {
            first_point: 99,
            point_count: 99
        };
        need.rings
    ];
    let mut xy = vec![f64::NAN; 2 * need.points];
    let mut triangles = vec![u32::MAX; need.triangles];
    assert_eq!(
        unsafe {
            axiolid_v0_4_visible_faces_copy(
                context,
                result,
                regions.as_mut_ptr(),
                need.regions,
                rings.as_mut_ptr(),
                need.rings,
                xy.as_mut_ptr(),
                need.points,
                triangles.as_mut_ptr(),
                need.triangles,
            )
        },
        AxiolidStatus::Ok
    );
    Tables {
        regions,
        rings,
        points: xy.chunks(2).map(|p| [p[0], p[1]]).collect(),
        triangles,
    }
}

fn evidence(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
) -> AxiolidVisibleFaceEvidence {
    let mut out = AxiolidVisibleFaceEvidence {
        source_meshes: 99,
        source_triangles: 99,
        edge_on_triangles: 99,
        removed_triangles: 99,
        boundary_segments: 99,
        intersection_segments: 99,
        crossings: 99,
        pieces: 99,
        tied_sides: 99,
        candidate_pairs: 99,
        exact_escalations: 99,
        regions: 99,
        cap_regions: 99,
        rerouted_edges: 99,
        collapsed_regions: 99,
        collapsed_holes: 99,
        snap_displacement: -1.0,
    };
    assert_eq!(
        unsafe { axiolid_v0_4_visible_faces_evidence(context, result, &mut out) },
        AxiolidStatus::Ok
    );
    out
}

fn last_error_mesh(context: AxiolidContextHandle) -> (AxiolidStatus, u8, u64) {
    let (mut has, mut id) = (9_u8, 9_u64);
    let status = unsafe { axiolid_v0_4_context_last_error_mesh(context, &mut has, &mut id) };
    (status, has, id)
}

const SHEET: i32 = AXIOLID_MESH_ROLE_SHEET;
const SOLID: i32 = AXIOLID_MESH_ROLE_SOLID;

/// Two nested sheets seen from above: the lower owns a ring-shaped region
/// with a hole where the higher one covers it.
fn nested_sheets(context: AxiolidContextHandle) -> Call {
    let low = sheet(context, 0.0, 4.0, 0.0);
    let high = sheet(context, 1.0, 3.0, 1.0);
    Call::new(vec![entry(low, 1, SHEET), entry(high, 2, SHEET)])
}

fn region(
    mesh_id: u64,
    rings: (usize, usize),
    triangles: (usize, usize),
    group: u32,
    source: AxiolidFaceSource,
) -> AxiolidVisibleRegion {
    AxiolidVisibleRegion {
        mesh_id,
        first_ring: rings.0,
        ring_count: rings.1,
        first_triangle: triangles.0,
        triangle_count: triangles.1,
        group,
        source,
    }
}

fn ring(first_point: usize, point_count: usize) -> AxiolidRegionRing {
    AxiolidRegionRing {
        first_point,
        point_count,
    }
}

const SQUARE: [[f64; 2]; 4] = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];

#[test]
fn nested_sheets_give_a_region_with_a_hole() {
    unsafe {
        let context = context();
        let (status, result) = nested_sheets(context).run(context);
        assert_eq!(status, AxiolidStatus::Ok);

        let mut kind = AxiolidGeometryKind::Unknown;
        assert_eq!(
            axiolid_v0_4_result_kind(context, result, &mut kind),
            AxiolidStatus::Ok
        );
        assert_eq!(kind, AxiolidGeometryKind::VisibleFaces);
        let mut taken = AxiolidMeshHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_result_take_mesh(context, result, &mut taken),
            AxiolidStatus::WrongResultKind
        );

        assert_eq!(
            counts(context, result),
            AxiolidVisibleFaceCounts {
                regions: 2,
                rings: 3,
                points: 12,
                triangles: 4,
            }
        );
        let tables = read_tables(context, result);
        // The first ring of a region is its outer ring, the second a hole.
        assert_eq!(
            tables.regions,
            [
                region(1, (0, 2), (0, 2), 0, AxiolidFaceSource::Group),
                region(2, (2, 1), (2, 2), 0, AxiolidFaceSource::Group),
            ]
        );
        assert_eq!(tables.rings, [ring(0, 4), ring(4, 4), ring(8, 4)]);
        assert_eq!(
            tables.points,
            [
                [0.0, 0.0],
                [4.0, 0.0],
                [4.0, 4.0],
                [0.0, 4.0],
                [1.0, 1.0],
                [1.0, 3.0],
                [3.0, 3.0],
                [3.0, 1.0],
                [1.0, 1.0],
                [3.0, 1.0],
                [3.0, 3.0],
                [1.0, 3.0],
            ]
        );
        assert_eq!(tables.triangles, [0, 1, 0, 1]);

        let evidence = evidence(context, result);
        assert_eq!(
            (
                evidence.source_meshes,
                evidence.source_triangles,
                evidence.boundary_segments,
                evidence.pieces,
                evidence.regions,
                evidence.cap_regions
            ),
            (2, 4, 8, 11, 2, 0)
        );
        assert_eq!(evidence.snap_displacement, 0.0);
    }
}

#[test]
fn the_near_plane_cap_of_a_cube_is_a_cap_region() {
    let context = context();
    let cube = cube(context);
    let mut call = Call::new(vec![entry(cube, 7, SOLID)]);
    call.view.has_near = 1;
    call.view.near = 1.0;
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let tables = read_tables(context, result);
    assert_eq!(
        tables.regions,
        [region(7, (0, 1), (0, 0), 0, AxiolidFaceSource::Cap)]
    );
    assert_eq!(tables.rings, [ring(0, 4)]);
    assert_eq!(tables.points, SQUARE);
    assert!(tables.triangles.is_empty());
    let evidence = evidence(context, result);
    assert_eq!((evidence.regions, evidence.cap_regions), (1, 1));
}

#[test]
fn groups_name_the_group_of_a_region() {
    let context = context();
    let cube = cube(context);
    let groups: [u32; 12] = [0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10];
    let mut mesh = entry(cube, 7, SOLID);
    mesh.groups = groups.as_ptr();
    mesh.group_count = groups.len();
    mesh.group_mode = AXIOLID_GROUP_MODE_PER_TRIANGLE;
    let (status, result) = Call::new(vec![mesh]).run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let tables = read_tables(context, result);
    // Seen from above, only the top face (triangles 2 and 3) is frontmost.
    assert_eq!(
        tables.regions,
        [region(7, (0, 1), (0, 2), 2, AxiolidFaceSource::Group)]
    );
    assert_eq!(tables.points, SQUARE);
    assert_eq!(tables.triangles, [2, 3]);

    // Planar groups are named by their smallest triangle.
    mesh.group_mode = AXIOLID_GROUP_MODE_PLANAR;
    mesh.groups = null();
    mesh.group_count = 0;
    let (status, result) = Call::new(vec![mesh]).run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    assert_eq!(
        read_tables(context, result).regions,
        [region(7, (0, 1), (0, 2), 2, AxiolidFaceSource::Group)]
    );
}

/// One scene run through the C ABI and through `Application` directly.
struct Scene {
    /// Positions, indices and role of each mesh; ids are 1, 2, ...
    meshes: Vec<(Vec<f64>, Vec<u32>, i32)>,
    view: AxiolidView,
    tolerance: AxiolidTolerance,
}

impl Scene {
    fn new(view: AxiolidView) -> Self {
        Self {
            meshes: Vec::new(),
            view,
            tolerance: TOLERANCE,
        }
    }

    fn with(mut self, positions: &[f64], indices: &[u32], role: i32) -> Self {
        self.meshes
            .push((positions.to_vec(), indices.to_vec(), role));
        self
    }

    /// The evidence from the C ABI and from the Rust application, and the
    /// number of regions the C ABI returned.
    fn run(&self) -> (AxiolidVisibleFaceEvidence, VisibleFaceEvidence, usize) {
        use axiolid::application::Application;
        use axiolid::contracts::ExecutionOptions;
        use axiolid::core::{Frame3, Point3, Tolerance, Vec3};
        use axiolid::mesh::TriMesh;
        use axiolid::mesh_hidden_lines::{
            HiddenLineMesh, MeshId, MeshRole, View, VisibleFaceLimits, VisibleFaceRequest,
        };

        let context = context();
        let entries = self
            .meshes
            .iter()
            .zip(1_u64..)
            .map(|((positions, indices, role), id)| {
                entry(import(context, positions, indices), id, *role)
            })
            .collect();
        let mut call = Call::new(entries);
        call.view = self.view;
        call.tolerance = self.tolerance;
        call.limits.max_meshes = 64;
        let (status, result) = call.run(context);
        assert_eq!(status, AxiolidStatus::Ok);
        let from_c = evidence(context, result);
        let regions = counts(context, result).regions;

        let meshes: Vec<TriMesh> = self
            .meshes
            .iter()
            .map(|(positions, indices, _)| {
                let points = positions
                    .chunks(3)
                    .map(|p| Point3::new(p[0], p[1], p[2]))
                    .collect();
                TriMesh::new(points, indices.clone())
            })
            .collect();
        let inputs: Vec<HiddenLineMesh<'_>> = meshes
            .iter()
            .zip(&self.meshes)
            .zip(1_u64..)
            .map(|((mesh, (_, _, role)), id)| {
                let role = if *role == SOLID {
                    MeshRole::Solid
                } else {
                    MeshRole::Sheet
                };
                HiddenLineMesh::new(MeshId(id), role, mesh)
            })
            .collect();
        let vec = |v: [f64; 3]| Vec3::new(v[0], v[1], v[2]);
        let frame = &self.view.frame;
        let mut view = View::new(Frame3 {
            origin: vec(frame.origin),
            x: vec(frame.x),
            y: vec(frame.y),
            z: vec(frame.z),
        });
        if self.view.has_near == 1 {
            view = view.with_near(self.view.near);
        }
        let limits = VisibleFaceLimits::new(
            call.limits.max_meshes,
            call.limits.max_source_vertices,
            call.limits.max_source_triangles,
            call.limits.max_candidate_pairs,
            call.limits.max_regions,
            call.limits.max_region_vertices,
        );
        let request = VisibleFaceRequest::new(&inputs, view, limits);
        let options = ExecutionOptions::new(
            Tolerance::new(self.tolerance.linear, self.tolerance.angular).expect("tolerance"),
        );
        let outcome = Application::portable()
            .expect("application")
            .visible_faces(&request, &options)
            .expect("visible faces");
        (from_c, outcome.evidence, regions)
    }
}

/// Every C field equals its Rust counterpart. `exact_escalations` is the one
/// count that may differ between identical requests, so it is only required
/// to be non-zero when the Rust count is.
fn assert_matches_rust(from_c: &AxiolidVisibleFaceEvidence, rust: &VisibleFaceEvidence) {
    assert_eq!(from_c.source_meshes, rust.source_meshes);
    assert_eq!(from_c.source_triangles, rust.source_triangles);
    assert_eq!(from_c.edge_on_triangles, rust.edge_on_triangles);
    assert_eq!(from_c.removed_triangles, rust.removed_triangles);
    assert_eq!(from_c.boundary_segments, rust.boundary_segments);
    assert_eq!(from_c.intersection_segments, rust.intersection_segments);
    assert_eq!(from_c.crossings, rust.crossings);
    assert_eq!(from_c.pieces, rust.pieces);
    assert_eq!(from_c.tied_sides, rust.tied_sides);
    assert_eq!(from_c.candidate_pairs, rust.candidate_pairs);
    assert_eq!(from_c.exact_escalations > 0, rust.exact_escalations > 0);
    assert_eq!(from_c.regions, rust.regions);
    assert_eq!(from_c.cap_regions, rust.cap_regions);
    assert_eq!(from_c.rerouted_edges, rust.rerouted_edges);
    assert_eq!(from_c.collapsed_regions, rust.collapsed_regions);
    assert_eq!(from_c.collapsed_holes, rust.collapsed_holes);
    assert_eq!(from_c.snap_displacement, rust.snap_displacement.get());
}

const CUBE_POSITIONS: [f64; 24] = [
    0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0,
    2.0, 2.0, 0.0, 2.0, 2.0,
];
const CUBE_INDICES: [u32; 36] = [
    0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0,
    4, 3, 4, 7,
];

fn sheet_positions(min: f64, max: f64, z: f64) -> [f64; 12] {
    [min, min, z, max, min, z, max, max, z, min, max, z]
}

/// A cube cut by a near plane over two nested sheets.
fn cube_over_sheets() -> Scene {
    let mut view = plan();
    view.has_near = 1;
    view.near = 1.0;
    Scene::new(view)
        .with(&CUBE_POSITIONS, &CUBE_INDICES, SOLID)
        .with(&sheet_positions(0.0, 4.0, 0.0), &[0, 1, 2, 0, 2, 3], SHEET)
        .with(&sheet_positions(1.0, 3.0, 1.0), &[0, 1, 2, 0, 2, 3], SHEET)
}

#[test]
fn the_evidence_matches_the_application_on_a_scene_with_caps_and_edge_on_faces() {
    let (from_c, rust, regions) = cube_over_sheets().run();
    assert_matches_rust(&from_c, &rust);
    // Properties of the scene, not of the provider's internals.
    assert_eq!(
        (
            from_c.source_meshes,
            from_c.source_triangles,
            from_c.edge_on_triangles,
            from_c.removed_triangles,
            from_c.regions,
            from_c.cap_regions
        ),
        (3, 16, 8, 2, 3, 1)
    );
    assert_eq!(regions, 3);
    assert_eq!(from_c.intersection_segments, 0);
    // Cap, boundaries and ties make these distinct from one another.
    assert!(from_c.boundary_segments > 0 && from_c.pieces > 0 && from_c.tied_sides > 0);
    assert!(from_c.candidate_pairs > 0 && from_c.exact_escalations > 0);
}

#[test]
fn the_evidence_matches_the_application_where_surfaces_intersect() {
    let mut scene = cube_over_sheets();
    // A tilted triangle through both sheets adds intersection curves.
    scene = scene.with(
        &[1.0, 1.0, -1.0, 3.0, 1.0, 1.0, 2.0, 3.0, 1.0],
        &[0, 1, 2],
        SHEET,
    );
    let (from_c, rust, _) = scene.run();
    assert_matches_rust(&from_c, &rust);
    assert_eq!(
        (
            from_c.source_meshes,
            from_c.source_triangles,
            from_c.edge_on_triangles,
            from_c.removed_triangles,
            from_c.regions,
            from_c.cap_regions
        ),
        (4, 17, 8, 2, 3, 1)
    );
    assert!(from_c.intersection_segments > 0 && from_c.crossings > 0);
}

/// The sheet `[-1, 1]^2` behind a triangle a few binary64 steps wide, in an
/// oblique view under a tolerance so small that the audit lets the triangle
/// be one. Its corners round to one cell: the hole it leaves in the sheet and
/// its own region collapse.
fn thin_occluder(second: bool) -> Scene {
    use axiolid::core::{Point3, Vec3};
    let direction = Vec3::new(0.3, -0.5, 0.8).normalize();
    let x = Vec3::Z.cross(direction).normalize();
    let y = direction.cross(x);
    let origin = Point3::new(-0.4, 0.7, 0.1);
    let axis = |v: Vec3| [v.x, v.y, v.z];
    let mut view = plan();
    view.frame = AxiolidViewFrame {
        origin: [origin.x, origin.y, origin.z],
        x: axis(x),
        y: axis(y),
        z: axis(direction),
    };
    let step = |x: f64| f64::from_bits(x.to_bits() + 1);
    let tiny = |x: f64, y: f64| [x, y, 1.0, step(x), y, 1.0, x, step(y), 1.0];
    let mut positions = tiny(0.3, 0.2).to_vec();
    let mut indices = vec![0, 1, 2];
    if second {
        // A second one outside the sheet's outline.
        positions.extend(tiny(-1.7, 0.2));
        indices.extend([3, 4, 5]);
    }
    let mut scene = Scene::new(view)
        .with(&sheet_positions(-1.0, 1.0, 0.0), &[0, 1, 2, 0, 2, 3], SHEET)
        .with(&positions, &indices, SHEET);
    scene.tolerance = AxiolidTolerance {
        linear: 1e-40,
        angular: 1e-40,
    };
    scene
}

#[test]
fn snap_evidence_is_reported_for_a_hole_and_an_occluder_that_collapse() {
    let (from_c, rust, regions) = thin_occluder(false).run();
    assert_matches_rust(&from_c, &rust);
    assert_eq!(
        (
            from_c.rerouted_edges,
            from_c.collapsed_regions,
            from_c.collapsed_holes,
            regions
        ),
        (1, 1, 1, 1)
    );
    assert_eq!(from_c.snap_displacement, 2f64.powi(-53));

    // Two occluders, one outside the sheet: the counts differ from one another.
    let (from_c, rust, _) = thin_occluder(true).run();
    assert_matches_rust(&from_c, &rust);
    assert_eq!(
        (
            from_c.rerouted_edges,
            from_c.collapsed_regions,
            from_c.collapsed_holes
        ),
        (2, 2, 1)
    );
    assert!(from_c.snap_displacement > 0.0);
}

/// The oracle's seeded scene 9 in its isometric view with a near plane at the
/// middle of the depth range: two edges are rerouted and one sliver region
/// collapses, with no collapsed hole.
#[test]
fn snap_evidence_is_reported_for_a_seeded_isometric_scene() {
    use seed9::{BOXES, FRAME, INDICES, NEAR};
    let bits = |bits: &[u64]| bits.iter().map(|b| f64::from_bits(*b)).collect::<Vec<_>>();
    let f = bits(&FRAME);
    let mut view = plan();
    view.frame = AxiolidViewFrame {
        origin: [f[0], f[1], f[2]],
        x: [f[3], f[4], f[5]],
        y: [f[6], f[7], f[8]],
        z: [f[9], f[10], f[11]],
    };
    view.has_near = 1;
    view.near = f64::from_bits(NEAR);
    let mut scene = Scene::new(view);
    for corners in &BOXES {
        let positions: Vec<f64> = corners.iter().flat_map(|p| bits(p)).collect();
        scene = scene.with(&positions, &INDICES, SOLID);
    }
    let (from_c, rust, regions) = scene.run();
    assert_matches_rust(&from_c, &rust);
    assert_eq!(
        (
            from_c.rerouted_edges,
            from_c.collapsed_regions,
            from_c.collapsed_holes,
            regions
        ),
        (2, 1, 0, 23)
    );
    assert_eq!(from_c.snap_displacement, 2f64.powi(-53));
}

#[test]
fn a_region_crops_the_image_and_limits_bound_the_output() {
    let context = context();
    let mut call = nested_sheets(context);
    call.region = Some(AxiolidImageBox {
        min: [0.0, 0.0],
        max: [2.0, 2.0],
    });
    let (status, result) = call.run(context);
    assert_eq!(status, AxiolidStatus::Ok);
    let tables = read_tables(context, result);
    assert_eq!(tables.regions.len(), 2);
    assert!(tables
        .points
        .iter()
        .all(|p| (0.0..=2.0).contains(&p[0]) && (0.0..=2.0).contains(&p[1])));

    // Each limit bounds its own quantity, exactly at the count.
    let ok = AxiolidStatus::Ok;
    let over = AxiolidStatus::LimitExceeded;
    let cases: Vec<(&str, Change, AxiolidStatus)> = vec![
        ("regions", Box::new(|c| c.limits.max_regions = 2), ok),
        ("regions over", Box::new(|c| c.limits.max_regions = 1), over),
        (
            "vertices",
            Box::new(|c| c.limits.max_region_vertices = 12),
            ok,
        ),
        (
            "vertices over",
            Box::new(|c| c.limits.max_region_vertices = 11),
            over,
        ),
        ("meshes", Box::new(|c| c.limits.max_meshes = 2), ok),
        ("meshes over", Box::new(|c| c.limits.max_meshes = 1), over),
        (
            "positions",
            Box::new(|c| c.limits.max_source_vertices = 8),
            ok,
        ),
        (
            "positions over",
            Box::new(|c| c.limits.max_source_vertices = 7),
            over,
        ),
        (
            "triangles",
            Box::new(|c| c.limits.max_source_triangles = 4),
            ok,
        ),
        (
            "triangles over",
            Box::new(|c| c.limits.max_source_triangles = 3),
            over,
        ),
        (
            "pairs over",
            Box::new(|c| c.limits.max_candidate_pairs = 1),
            over,
        ),
    ];
    for (name, change, expected) in cases {
        let mut call = nested_sheets(context);
        change(&mut call);
        assert_eq!(call.run(context).0, expected, "{name}");
    }
}

#[test]
fn copy_writes_nothing_unless_every_capacity_covers_its_count() {
    unsafe {
        let context = context();
        let (_, result) = nested_sheets(context).run(context);
        let mut regions = [EMPTY_REGION; 2];
        let mut rings = [ring(99, 99); 3];
        let mut xy = [f64::NAN; 24];
        let mut triangles = [u32::MAX; 4];
        // Capacities: regions, rings, points, triangles. Each one short in turn.
        for short in 0..4 {
            let mut cap = [2, 3, 12, 4];
            cap[short] -= 1;
            let status = axiolid_v0_4_visible_faces_copy(
                context,
                result,
                regions.as_mut_ptr(),
                cap[0],
                rings.as_mut_ptr(),
                cap[1],
                xy.as_mut_ptr(),
                cap[2],
                triangles.as_mut_ptr(),
                cap[3],
            );
            assert_eq!(status, AxiolidStatus::BufferTooSmall, "table {short}");
            assert_eq!(regions, [EMPTY_REGION; 2], "table {short}");
            assert_eq!(rings, [ring(99, 99); 3], "table {short}");
            assert!(xy.iter().all(|v| v.is_nan()), "table {short}");
            assert_eq!(triangles, [u32::MAX; 4], "table {short}");
        }
        // Exact capacities write the counted extent and no more.
        let status = axiolid_v0_4_visible_faces_copy(
            context,
            result,
            regions.as_mut_ptr(),
            2,
            rings.as_mut_ptr(),
            3,
            xy.as_mut_ptr(),
            12,
            triangles.as_mut_ptr(),
            4,
        );
        assert_eq!(status, AxiolidStatus::Ok);
        assert!(xy[..24].iter().all(|v| !v.is_nan()));

        // A null table is refused only when its capacity is non-zero.
        let null_status = |region_cap, ring_cap, point_cap, triangle_cap| {
            axiolid_v0_4_visible_faces_copy(
                context,
                result,
                std::ptr::null_mut(),
                region_cap,
                std::ptr::null_mut(),
                ring_cap,
                std::ptr::null_mut(),
                point_cap,
                std::ptr::null_mut(),
                triangle_cap,
            )
        };
        assert_eq!(null_status(0, 0, 0, 0), AxiolidStatus::BufferTooSmall);
        assert_eq!(null_status(2, 0, 0, 0), AxiolidStatus::NullPointer);
        assert_eq!(null_status(0, 3, 0, 0), AxiolidStatus::NullPointer);
        assert_eq!(null_status(0, 0, 12, 0), AxiolidStatus::NullPointer);
        assert_eq!(null_status(0, 0, 0, 4), AxiolidStatus::NullPointer);
    }
}

#[test]
fn an_empty_result_copies_with_zero_capacities() {
    unsafe {
        let context = context();
        let sheet = sheet(context, 0.0, 1.0, 0.0);
        let mut call = Call::new(vec![entry(sheet, 1, SHEET)]);
        // The sheet lies at W = 0, in front of the far plane at 5: clipped away.
        call.view.has_far = 1;
        call.view.far = 5.0;
        let (status, result) = call.run(context);
        assert_eq!(status, AxiolidStatus::Ok);
        assert_eq!(
            counts(context, result),
            AxiolidVisibleFaceCounts {
                regions: 0,
                rings: 0,
                points: 0,
                triangles: 0
            }
        );
        let status = axiolid_v0_4_visible_faces_copy(
            context,
            result,
            std::ptr::null_mut(),
            0,
            std::ptr::null_mut(),
            0,
            std::ptr::null_mut(),
            0,
            std::ptr::null_mut(),
            0,
        );
        assert_eq!(status, AxiolidStatus::Ok);
    }
}

#[test]
fn a_refusal_names_the_mesh_it_concerns() {
    let context = context();
    assert_eq!(last_error_mesh(context).0, AxiolidStatus::NoError);
    let triangle = triangle(context);
    // An open triangle is no solid.
    let (status, _) = Call::new(vec![entry(triangle, 9, SOLID)]).run(context);
    assert_eq!(status, AxiolidStatus::OperationFailed);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 9));
    let mut info = AxiolidErrorInfo::default();
    assert_eq!(
        unsafe { axiolid_v0_4_context_last_error_info(context, &mut info) },
        AxiolidStatus::Ok
    );
    assert_eq!(info.operation, AxiolidOperation::MeshVisibleFaces);
    assert_eq!(info.status, AxiolidStatus::OperationFailed);

    // A duplicate id names that id.
    let call = Call::new(vec![entry(triangle, 3, SHEET), entry(triangle, 3, SHEET)]);
    assert_eq!(call.run(context).0, AxiolidStatus::OperationFailed);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 3));

    // A refusal of the whole request names no mesh and keeps its status.
    let mut call = Call::new(vec![entry(triangle, 9, SHEET)]);
    call.limits.max_meshes = 0;
    assert_eq!(call.run(context).0, AxiolidStatus::LimitExceeded);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 0, 0));
    unsafe { axiolid_v0_4_context_last_error_info(context, &mut info) };
    assert_eq!(info.operation, AxiolidOperation::MeshVisibleFaces);
    assert_eq!(info.status, AxiolidStatus::LimitExceeded);
}

#[test]
fn invalid_arguments_are_refused_at_the_boundary() {
    let context = context();
    let triangle = triangle(context);
    let base = || Call::new(vec![entry(triangle, 1, SHEET)]);
    let cases: Vec<(&str, Change)> = vec![
        ("projection", Box::new(|c| c.view.projection = 2)),
        ("view reserved", Box::new(|c| c.view.reserved = [0, 1])),
        ("has_near", Box::new(|c| c.view.has_near = 2)),
        ("has_far", Box::new(|c| c.view.has_far = 2)),
        (
            "origin NaN",
            Box::new(|c| c.view.frame.origin[1] = f64::NAN),
        ),
        (
            "unused far infinite",
            Box::new(|c| c.view.far = f64::INFINITY),
        ),
        (
            "region NaN",
            Box::new(|c| {
                c.region = Some(AxiolidImageBox {
                    min: [f64::NAN, 0.0],
                    max: [1.0, 1.0],
                })
            }),
        ),
        ("role", Box::new(|c| c.meshes[0].role = 0)),
        ("group mode", Box::new(|c| c.meshes[0].group_mode = 4)),
        (
            "single with groups",
            Box::new(|c| c.meshes[0].group_count = 1),
        ),
        ("tolerance", Box::new(|c| c.tolerance.linear = f64::NAN)),
    ];
    for (name, change) in cases {
        let mut call = base();
        change(&mut call);
        let (status, result) = call.run(context);
        assert_eq!(status, AxiolidStatus::InvalidArgument, "{name}");
        assert_eq!(result, AxiolidResultHandle::INVALID, "{name}");
        let mut info = AxiolidErrorInfo::default();
        assert_eq!(
            unsafe { axiolid_v0_4_context_last_error_info(context, &mut info) },
            AxiolidStatus::Ok,
            "{name}"
        );
        assert_eq!(info.status, AxiolidStatus::InvalidArgument, "{name}");
        assert_eq!(info.operation, AxiolidOperation::MeshVisibleFaces, "{name}");
    }
}

#[test]
fn group_counts_are_checked_against_the_mesh_before_they_are_read() {
    let context = context();
    let triangle = triangle(context);
    let groups = [0_u32; 4];
    let mut mesh = entry(triangle, 4, SHEET);
    mesh.group_mode = AXIOLID_GROUP_MODE_PER_TRIANGLE;
    mesh.groups = groups.as_ptr();
    // One triangle, four ids: refused before a slice of four is built.
    mesh.group_count = 4;
    assert_eq!(
        Call::new(vec![mesh]).run(context).0,
        AxiolidStatus::InvalidArgument
    );
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 4));
    // A count whose byte length overflows is refused the same way.
    mesh.group_count = usize::MAX;
    assert_eq!(
        Call::new(vec![mesh]).run(context).0,
        AxiolidStatus::InvalidArgument
    );
    // A null array for a non-empty mesh is refused.
    mesh.groups = null();
    mesh.group_count = 1;
    assert_eq!(
        Call::new(vec![mesh]).run(context).0,
        AxiolidStatus::NullPointer
    );
    // The count is the mesh's own: a two-triangle mesh needs two ids.
    let square = square(context);
    let mut two = entry(square, 5, SHEET);
    two.group_mode = AXIOLID_GROUP_MODE_PER_TRIANGLE;
    two.groups = groups.as_ptr();
    two.group_count = 1;
    assert_eq!(
        Call::new(vec![two]).run(context).0,
        AxiolidStatus::InvalidArgument
    );
    two.group_count = 2;
    assert_eq!(Call::new(vec![two]).run(context).0, AxiolidStatus::Ok);
}

#[test]
fn handles_pointers_and_counts_are_checked() {
    let config = AxiolidContextConfig {
        max_meshes: 2,
        ..AxiolidContextConfig::default()
    };
    let context = context_with(config);
    let triangle = triangle(context);
    let call = Call::new(vec![entry(triangle, 1, SHEET)]);

    let unknown = Call::new(vec![entry(AxiolidMeshHandle(u64::MAX), 5, SHEET)]);
    assert_eq!(unknown.run(context).0, AxiolidStatus::InvalidHandle);
    assert_eq!(last_error_mesh(context), (AxiolidStatus::Ok, 1, 5));
    assert_eq!(
        call.run(AxiolidContextHandle(u64::MAX)).0,
        AxiolidStatus::InvalidHandle
    );
    let many = Call::new(vec![entry(triangle, 1, SHEET); 3]);
    assert_eq!(many.run(context).0, AxiolidStatus::LimitExceeded);

    unsafe {
        let mut result = AxiolidResultHandle::INVALID;
        let (view, limits) = (plan(), limits());
        let entries = [entry(triangle, 1, SHEET)];
        let run = |meshes, count, view: *const _, limits: *const _, out| {
            axiolid_v0_4_visible_faces(context, meshes, count, view, null(), limits, TOLERANCE, out)
        };
        let null_pointer = |status| assert_eq!(status, AxiolidStatus::NullPointer);
        null_pointer(run(null(), 1, &view, &limits, &mut result));
        null_pointer(run(entries.as_ptr(), 1, null(), &limits, &mut result));
        null_pointer(run(entries.as_ptr(), 1, &view, null(), &mut result));
        null_pointer(run(
            entries.as_ptr(),
            1,
            &view,
            &limits,
            std::ptr::null_mut(),
        ));
        // An empty mesh list reaches the contract, which refuses it.
        assert_eq!(
            run(null(), 0, &view, &limits, &mut result),
            AxiolidStatus::OperationFailed
        );

        // Result readers check the pointer, the handle and the kind.
        let (status, drawn) = call.run(context);
        assert_eq!(status, AxiolidStatus::Ok);
        let bad = AxiolidResultHandle(u64::MAX);
        let mut counts = AxiolidVisibleFaceCounts {
            regions: 0,
            rings: 0,
            points: 0,
            triangles: 0,
        };
        assert_eq!(
            axiolid_v0_4_visible_faces_counts(context, drawn, std::ptr::null_mut()),
            AxiolidStatus::NullPointer
        );
        assert_eq!(
            axiolid_v0_4_visible_faces_counts(context, bad, &mut counts),
            AxiolidStatus::InvalidHandle
        );
        assert_eq!(
            axiolid_v0_4_visible_faces_copy(
                context,
                bad,
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                0,
                std::ptr::null_mut(),
                0
            ),
            AxiolidStatus::InvalidHandle
        );
        assert_eq!(
            axiolid_v0_4_visible_faces_evidence(context, drawn, std::ptr::null_mut()),
            AxiolidStatus::NullPointer
        );

        // Another kind of result is not a visible-face result, and vice versa.
        let mut exact = AxiolidResultHandle::INVALID;
        assert_eq!(
            axiolid_v0_4_exact_extrude_rectangle(context, 1.0, 1.0, 1.0, TOLERANCE, &mut exact),
            AxiolidStatus::Ok
        );
        assert_eq!(
            axiolid_v0_4_visible_faces_counts(context, exact, &mut counts),
            AxiolidStatus::WrongResultKind
        );
        let mut count = 0_usize;
        assert_eq!(
            axiolid_v0_4_hidden_lines_count(context, drawn, &mut count),
            AxiolidStatus::WrongResultKind
        );
    }
}

#[test]
fn the_context_lists_the_visible_face_capability() {
    unsafe {
        let context = context();
        let mut count = 0_usize;
        assert_eq!(
            axiolid_v0_4_capability_count(context, &mut count),
            AxiolidStatus::Ok
        );
        assert_eq!(count, 6);
        let found = (0..count)
            .map(|index| {
                let mut descriptor = AxiolidCapabilityDescriptor::default();
                assert_eq!(
                    axiolid_v0_4_capability_get(context, index, &mut descriptor),
                    AxiolidStatus::Ok
                );
                descriptor
            })
            .find(|descriptor| descriptor.operation as i32 == 13)
            .expect("visible faces is listed");
        assert_eq!(found.representation, AxiolidRepresentation::Profile2d);
        let provider = &found.provider[..usize::from(found.provider_len)];
        assert_eq!(provider, b"scalar-visible-faces");
    }
}

#[test]
fn the_documented_layout_is_stable() {
    use std::mem::size_of;
    assert_eq!(size_of::<AxiolidVisibleFaceLimits>(), 48);
    assert_eq!(size_of::<AxiolidVisibleRegion>(), 48);
    assert_eq!(size_of::<AxiolidRegionRing>(), 16);
    assert_eq!(size_of::<AxiolidVisibleFaceCounts>(), 32);
    assert_eq!(size_of::<AxiolidVisibleFaceEvidence>(), 136);
    assert_eq!(AxiolidFaceSource::Group as i32, 1);
    assert_eq!(AxiolidFaceSource::Cap as i32, 2);
    assert_eq!(AxiolidFaceSource::Unknown as i32, 255);
    assert_eq!(AxiolidGeometryKind::VisibleFaces as i32, 4);
}
