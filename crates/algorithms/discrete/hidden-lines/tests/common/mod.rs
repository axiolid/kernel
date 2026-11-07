//! Scene builders shared by the integration tests. Each test binary uses a
//! different subset, hence the `dead_code` allowance on the module.
#![allow(dead_code)]

pub mod topology;

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_hidden_lines::ScalarHiddenLines;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, FaceSource, HiddenLineError, HiddenLineLimits, HiddenLineMesh,
    HiddenLineOutcome, HiddenLineRequest, MeshHiddenLines, MeshId, MeshRole, MeshVisibleFaces,
    OrthographicView, Visibility, VisibleFaceLimits, VisibleFaceOutcome, VisibleFaceRequest,
    VisibleRegion,
};

/// The twelve edges of a box mesh built by [`box_mesh`], as ascending vertex
/// pairs (the six face diagonals are not box edges).
pub const BOX_EDGES: [[u32; 2]; 12] = [
    [0, 1],
    [1, 2],
    [2, 3],
    [0, 3],
    [4, 5],
    [5, 6],
    [6, 7],
    [4, 7],
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],
];

/// An outward-wound box: 0..3 the bottom corners counter-clockwise seen from
/// above, 4..7 the top corners above them.
pub fn box_mesh(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let (x0, y0, z0, x1, y1, z1) = (min[0], min[1], min[2], max[0], max[1], max[2]);
    TriMesh::new(
        vec![
            Point3::new(x0, y0, z0),
            Point3::new(x1, y0, z0),
            Point3::new(x1, y1, z0),
            Point3::new(x0, y1, z0),
            Point3::new(x0, y0, z1),
            Point3::new(x1, y0, z1),
            Point3::new(x1, y1, z1),
            Point3::new(x0, y1, z1),
        ],
        vec![
            0, 3, 2, 0, 2, 1, // bottom
            4, 5, 6, 4, 6, 7, // top
            0, 1, 5, 0, 5, 4, // y = min
            1, 2, 6, 1, 6, 5, // x = max
            2, 3, 7, 2, 7, 6, // y = max
            3, 0, 4, 3, 4, 7, // x = min
        ],
    )
}

pub fn frame(x: Vec3, y: Vec3, z: Vec3, origin: Point3) -> Frame3 {
    Frame3 { origin, x, y, z }
}

/// Looking down the world z axis: image (x, y), nearer is higher.
pub fn plan_frame() -> Frame3 {
    frame(Vec3::X, Vec3::Y, Vec3::Z, Point3::ZERO)
}

/// The isometric frame of the conformance suite, about `origin`.
pub fn isometric_frame(origin: Point3) -> Frame3 {
    let (s2, s6, s3) = (
        std::f64::consts::FRAC_1_SQRT_2,
        1.0 / 6f64.sqrt(),
        1.0 / 3f64.sqrt(),
    );
    frame(
        Vec3::new(s2, -s2, 0.0),
        Vec3::new(s6, s6, -2.0 * s6),
        Vec3::new(s3, s3, s3),
        origin,
    )
}

/// A right-handed orthonormal frame looking along `-direction`.
pub fn looking_along(direction: Vec3, origin: Point3) -> Frame3 {
    let z = direction.normalize();
    let up = if z.z.abs() < 0.9 { Vec3::Z } else { Vec3::X };
    let x = up.cross(z).normalize();
    let y = z.cross(x);
    frame(x, y, z, origin)
}

pub const LIMITS: HiddenLineLimits = HiddenLineLimits::new(64, 4096, 4096, 10_000_000, 1_000_000);

/// The scene's meshes with ids 1, 2, ... in order, all solid.
pub fn solids(meshes: &[TriMesh]) -> Vec<HiddenLineMesh<'_>> {
    meshes
        .iter()
        .enumerate()
        .map(|(i, mesh)| HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Solid, mesh))
        .collect()
}

pub fn try_draw(
    provider: &dyn MeshHiddenLines,
    meshes: &[HiddenLineMesh<'_>],
    view: OrthographicView,
) -> Result<HiddenLineOutcome, HiddenLineError> {
    let request = HiddenLineRequest::new(meshes, view, EdgeSelection::new(Some(0.5)), LIMITS);
    let options = axiolid_contracts::ExecutionOptions::new(
        axiolid_core::Tolerance::new(1e-9, 1e-9).expect("tolerance"),
    );
    provider.hidden_lines(&request, &options)
}

/// `try_draw` with the caller's edge selection.
pub fn try_draw_with(
    provider: &dyn MeshHiddenLines,
    meshes: &[HiddenLineMesh<'_>],
    view: OrthographicView,
    edges: EdgeSelection,
) -> Result<HiddenLineOutcome, HiddenLineError> {
    let request = HiddenLineRequest::new(meshes, view, edges, LIMITS);
    let options = axiolid_contracts::ExecutionOptions::new(
        axiolid_core::Tolerance::new(1e-9, 1e-9).expect("tolerance"),
    );
    provider.hidden_lines(&request, &options)
}

pub fn draw(
    provider: &dyn MeshHiddenLines,
    meshes: &[HiddenLineMesh<'_>],
    view: OrthographicView,
) -> HiddenLineOutcome {
    try_draw(provider, meshes, view).unwrap_or_else(|error| panic!("{error}"))
}

pub const FACE_LIMITS: VisibleFaceLimits =
    VisibleFaceLimits::new(64, 4096, 4096, 10_000_000, 100_000, 1_000_000);

pub fn try_faces(
    provider: &dyn MeshVisibleFaces,
    meshes: &[HiddenLineMesh<'_>],
    view: OrthographicView,
    region: Option<axiolid_core::Aabb2>,
) -> Result<VisibleFaceOutcome, HiddenLineError> {
    let mut request = VisibleFaceRequest::new(meshes, view, FACE_LIMITS);
    if let Some(region) = region {
        request = request.with_region(region);
    }
    let options = axiolid_contracts::ExecutionOptions::new(
        axiolid_core::Tolerance::new(1e-9, 1e-9).expect("tolerance"),
    );
    provider.visible_faces(&request, &options)
}

pub fn faces(
    provider: &dyn MeshVisibleFaces,
    meshes: &[HiddenLineMesh<'_>],
    view: OrthographicView,
) -> VisibleFaceOutcome {
    try_faces(provider, meshes, view, None).unwrap_or_else(|error| panic!("{error}"))
}

/// Shoelace area of a ring (positive counter-clockwise).
pub fn ring_area(ring: &[axiolid_core::Point2]) -> f64 {
    let n = ring.len();
    (0..n)
        .map(|i| ring[i].x * ring[(i + 1) % n].y - ring[(i + 1) % n].x * ring[i].y)
        .sum::<f64>()
        / 2.0
}

/// Area of the regions of one source of one mesh, holes subtracted.
pub fn source_area(outcome: &VisibleFaceOutcome, mesh: u64, source: FaceSource) -> f64 {
    outcome
        .regions
        .iter()
        .filter(|r: &&VisibleRegion| r.mesh == MeshId(mesh) && r.source == source)
        .map(|r| ring_area(&r.outer) + r.holes.iter().map(|h| ring_area(h)).sum::<f64>())
        .sum()
}

pub fn scalar() -> ScalarHiddenLines {
    ScalarHiddenLines::new()
}

/// Total image length of the segments with `visibility`.
pub fn length(outcome: &HiddenLineOutcome, visibility: Visibility) -> f64 {
    outcome
        .segments
        .iter()
        .filter(|s| s.visibility == visibility)
        .map(|s| (s.end - s.start).length())
        .sum()
}

// --- seeded scenes shared by the oracles -------------------------------------

pub struct Lcg(pub u64);

impl Lcg {
    pub fn below(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) % n
    }
}

/// 1 to 12 boxes on a 1/8 grid; a third of them rotated about their centre.
pub fn scene(seed: u64) -> Vec<TriMesh> {
    let mut rng = Lcg(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let count = 1 + rng.below(12) as usize;
    let about_z: [f64; 6] = [0.31, 0.57, 0.93, 1.21, 2.17, 2.6];
    let about_x: [f64; 3] = [0.2, 0.4, 0.7];
    let about_y: [f64; 3] = [0.15, 0.35, 0.55];
    (0..count)
        .map(|_| {
            let min = [
                rng.below(13) as f64 / 8.0,
                rng.below(13) as f64 / 8.0,
                rng.below(9) as f64 / 8.0,
            ];
            let size = [
                (1 + rng.below(8)) as f64 / 8.0,
                (1 + rng.below(8)) as f64 / 8.0,
                (1 + rng.below(8)) as f64 / 8.0,
            ];
            let max = [min[0] + size[0], min[1] + size[1], min[2] + size[2]];
            let mut mesh = box_mesh(min, max);
            if rng.below(3) == 0 {
                let centre = Point3::new(
                    (min[0] + max[0]) / 2.0,
                    (min[1] + max[1]) / 2.0,
                    (min[2] + max[2]) / 2.0,
                );
                let (sz, cz) = about_z[rng.below(6) as usize].sin_cos();
                let (sx, cx) = about_x[rng.below(3) as usize].sin_cos();
                let (sy, cy) = about_y[rng.below(3) as usize].sin_cos();
                for p in &mut mesh.positions {
                    let d = *p - centre;
                    // Generic angles about all three axes: no face is left
                    // parallel to an axis view, so none is edge-on up to noise.
                    let d = Vec3::new(cy * d.x + sy * d.z, d.y, -sy * d.x + cy * d.z);
                    let tilted = Vec3::new(d.x, cx * d.y - sx * d.z, sx * d.y + cx * d.z);
                    *p = centre
                        + Vec3::new(
                            cz * tilted.x - sz * tilted.y,
                            sz * tilted.x + cz * tilted.y,
                            tilted.z,
                        );
                }
            }
            mesh
        })
        .collect()
}

pub fn views() -> Vec<(&'static str, Frame3)> {
    let o = Point3::ZERO;
    vec![
        ("plan", plan_frame()),
        ("+X", frame(Vec3::Y, Vec3::Z, Vec3::X, o)),
        ("-X", frame(Vec3::Z, Vec3::Y, -Vec3::X, o)),
        ("+Y", frame(Vec3::Z, Vec3::X, Vec3::Y, o)),
        ("-Y", frame(Vec3::X, Vec3::Z, -Vec3::Y, o)),
        ("isometric", isometric_frame(Point3::new(0.1, 0.2, 0.3))),
        (
            "oblique",
            looking_along(Vec3::new(0.3, -0.5, 0.8), Point3::new(-0.4, 0.7, 0.1)),
        ),
    ]
}
