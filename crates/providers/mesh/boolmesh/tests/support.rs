//! Shared fixture builders. Not a test file by itself.
//!
//! Included via `mod support;` in each test binary, so any helper a given
//! binary does not use is dead code there. The allow is on the module, not on
//! individual items, so adding a helper never needs a second annotation.

#![allow(dead_code)]

use axiolid_core::Point3;
use axiolid_mesh::TriMesh;

/// Outward-oriented box. Corner order is `corner*2 + (0 = bottom, 1 = top)`.
///
/// The index list is fixed and verified by `winding.rs`; do not "tidy" it
/// without re-running that gate, because an inverted box is silently accepted
/// by every structural check and only shows up as a wrong volume.
pub fn boxx(cx: f64, cy: f64, z0: f64, sx: f64, sy: f64, sz: f64, angle: f64) -> TriMesh {
    let (hx, hy) = (sx / 2.0, sy / 2.0);
    let (sin, cos) = angle.sin_cos();
    let mut positions = Vec::with_capacity(8);
    for &(dx, dy) in &[(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy)] {
        let (rx, ry) = (dx * cos - dy * sin, dx * sin + dy * cos);
        for &z in &[z0, z0 + sz] {
            positions.push(Point3::new(cx + rx, cy + ry, z));
        }
    }
    let indices = vec![
        0, 4, 2, 0, 6, 4, // bottom
        1, 3, 5, 1, 5, 7, // top
        0, 3, 1, 0, 2, 3, // side 0-1
        2, 5, 3, 2, 4, 5, // side 1-2
        4, 7, 5, 4, 6, 7, // side 2-3
        6, 1, 7, 6, 0, 1, // side 3-0
    ];
    TriMesh::new(positions, indices)
}

/// Absolute enclosed volume via the divergence theorem.
///
/// Independent of the crate's internal helper on purpose: a test that reuses
/// the implementation's own arithmetic cannot detect an error in it.
pub fn volume(mesh: &TriMesh) -> f64 {
    let mut total = 0.0;
    for corner in mesh.indices.as_chunks::<3>().0 {
        let a = mesh.positions[corner[0] as usize];
        let b = mesh.positions[corner[1] as usize];
        let c = mesh.positions[corner[2] as usize];
        total += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (total / 6.0).abs()
}

/// Reverse every triangle, producing a structurally valid inside-out mesh.
pub fn inverted(mesh: &TriMesh) -> TriMesh {
    let mut flipped = mesh.clone();
    for corner in flipped.indices.as_chunks_mut::<3>().0 {
        corner.swap(1, 2);
    }
    flipped
}
