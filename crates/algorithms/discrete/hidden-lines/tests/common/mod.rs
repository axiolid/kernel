//! Scene builders shared by the integration tests. Each test binary uses a
//! different subset, hence the `dead_code` allowance on the module.
#![allow(dead_code)]

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_hidden_lines::ScalarHiddenLines;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineError, HiddenLineLimits, HiddenLineMesh, HiddenLineOutcome,
    HiddenLineRequest, MeshHiddenLines, MeshId, MeshRole, OrthographicView, Visibility,
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

pub fn draw(
    provider: &dyn MeshHiddenLines,
    meshes: &[HiddenLineMesh<'_>],
    view: OrthographicView,
) -> HiddenLineOutcome {
    try_draw(provider, meshes, view).unwrap_or_else(|error| panic!("{error}"))
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
