//! Axis-aligned box recognition.
//!
//! The analytic path in [`super::cellular`] is exact only when every operand is
//! an axis-aligned box. This module decides that question and is deliberately
//! conservative: a false positive silently produces wrong geometry, while a
//! false negative only costs the general kernel's time.

use axiolid_core::Point3;
use axiolid_mesh::TriMesh;

/// An axis-aligned bounding box recognised from a mesh.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignedBox {
    /// Minimum corner.
    pub min: [f64; 3],
    /// Maximum corner.
    pub max: [f64; 3],
}

impl AlignedBox {
    /// Extent along each axis.
    fn extent(&self) -> [f64; 3] {
        [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ]
    }

    /// Whether this box has positive volume at the given tolerance.
    fn is_proper(&self, eps: f64) -> bool {
        self.extent().iter().all(|e| *e > eps)
    }
}

/// Recognise `mesh` as an axis-aligned box, or return `None`.
///
/// Verified structurally rather than by bounding box, which is the whole point:
/// every mesh HAS a bounding box, so accepting one on that basis would feed
/// arbitrary geometry to a path that is only correct for boxes. A sphere and
/// the cube around it share a bounding box and must not share an answer.
///
/// Requirements, all necessary:
/// * exactly 36 indices (12 triangles)
/// * every position, referenced or not, lies on the min/max lattice of each
///   axis
/// * every triangle lies in one of the 6 axis-aligned face planes
/// * each face plane carries exactly 2 triangles
///
/// Together these exclude a box with a dent (wrong triangle count), a sheared
/// box (corners off the min/max lattice), and a box with an interior void
/// (extra triangles).
///
/// None of the checks is redundant. `chunks_exact(3)` drops a trailing
/// partial triangle, so a malformed 38-index buffer presents a perfect box
/// to the plane loop and only the length check sees it. The plane loop sees
/// only referenced positions, so an unused off-lattice vertex is caught only
/// by the lattice walk. Each decline case is pinned in
/// `tests/analytic_boxes.rs`.
pub fn recognise(mesh: &TriMesh, eps: f64) -> Option<AlignedBox> {
    if mesh.indices.len() != 36 {
        return None;
    }
    let p = &mesh.positions;
    if p.is_empty() {
        return None;
    }

    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for v in p {
        let c = [v.x, v.y, v.z];
        for a in 0..3 {
            if !c[a].is_finite() {
                return None;
            }
            min[a] = min[a].min(c[a]);
            max[a] = max[a].max(c[a]);
        }
    }

    let bbox = AlignedBox { min, max };
    if !bbox.is_proper(eps) {
        return None;
    }

    // Every position, referenced or not, must sit on a corner of the
    // lattice: each coordinate equals that axis's min or max. This is what rejects a sphere,
    // a sheared box, or any mesh that merely happens to span the same extent.
    for v in p {
        let c = [v.x, v.y, v.z];
        for a in 0..3 {
            let on_min = (c[a] - min[a]).abs() <= eps;
            let on_max = (c[a] - max[a]).abs() <= eps;
            if !on_min && !on_max {
                return None;
            }
        }
    }

    // Each of the 6 face planes must carry exactly 2 triangles. A closed box
    // has 12; any other distribution means the surface is not a plain box even
    // if all corners are lattice points (e.g. a face split into 4 triangles
    // with another face missing).
    let mut plane_counts = [0usize; 6];
    for tri in mesh.indices.chunks_exact(3) {
        let vs: [[f64; 3]; 3] = [
            corner_of(p, tri[0])?,
            corner_of(p, tri[1])?,
            corner_of(p, tri[2])?,
        ];
        let mut matched = None;
        for a in 0..3 {
            if vs.iter().all(|v| (v[a] - min[a]).abs() <= eps) {
                matched = Some(a * 2);
                break;
            }
            if vs.iter().all(|v| (v[a] - max[a]).abs() <= eps) {
                matched = Some(a * 2 + 1);
                break;
            }
        }
        // A triangle spanning two planes is a diagonal: not a box face.
        plane_counts[matched?] += 1;
    }
    if plane_counts.iter().any(|c| *c != 2) {
        return None;
    }

    Some(bbox)
}

/// Fetch a vertex by index as a plain coordinate triple.
fn corner_of(positions: &[Point3], index: u32) -> Option<[f64; 3]> {
    let v = positions.get(index as usize)?;
    Some([v.x, v.y, v.z])
}
