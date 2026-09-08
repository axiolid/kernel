//! `TriMesh` <-> `csg::Manifold` conversion.
//!
//! Orientation is the dangerous part. An inside-out solid is structurally
//! valid: every edge still has exactly two incident faces, so a manifold check
//! passes. `boolmesh` then treats it as a solid whose interior is everywhere
//! outside, and `Difference` behaves like `Union` -- returning a larger mesh,
//! with no error and a plausible-looking triangle count.
//!
//! This was not hypothetical: the ADR 0014 evaluation harness hit exactly this,
//! and the wall's volume *grew* from 2.40 to 2.86 after subtracting three
//! openings. The failure is dangerous precisely because it is quiet.
//!
//! So orientation is checked here, on the way in, using the divergence theorem.

use crate::csg::Manifold;
use axiolid_contracts::GeomError;
use axiolid_core::Point3;
use axiolid_mesh::TriMesh;
use axiolid_mesh_contracts::SolidRequirements;

/// Six times the signed volume of a closed triangle mesh.
///
/// Sums the scalar triple product over every triangle (the divergence theorem
/// applied to `F = (x,y,z)/3`). Positive means outward-facing normals under the
/// right-hand rule. The factor of six is left in: only the sign and a
/// relative-magnitude comparison are ever needed, and dividing would add a
/// rounding step for no benefit.
///
/// This is O(triangles) with no allocation, so it is affordable on every call.
pub(crate) fn six_signed_volume(positions: &[Point3], indices: &[u32]) -> f64 {
    let mut total = 0.0;
    for corner in indices.chunks_exact(3) {
        let a = positions[corner[0] as usize];
        let b = positions[corner[1] as usize];
        let c = positions[corner[2] as usize];
        total += a.dot(b.cross(c));
    }
    total
}

/// Convert a `TriMesh` into a `csg::Manifold`, rejecting inputs whose
/// orientation would silently invert the operation.
///
/// `role` names the argument for the diagnostic, so a caller learns *which*
/// mesh was wrong rather than that some mesh was.
pub(crate) fn to_manifold(mesh: &TriMesh, role: &str) -> Result<Manifold, GeomError> {
    SolidRequirements::Oriented.validate(mesh, role)?;

    let positions: Vec<f64> = mesh
        .positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    let indices: Vec<usize> = mesh.indices.iter().map(|&i| i as usize).collect();

    Manifold::new(&positions, &indices)
        .map_err(|reason| GeomError::NotManifold(format!("{role}: {reason}")))
}

/// Convert a `csg::Manifold` back into a `TriMesh`.
///
/// The result carries no normals: `boolmesh` computes its own face normals, and
/// re-exporting them as *vertex* normals would misrepresent hard edges created
/// by the cut. Downstream code that needs normals should derive them from the
/// topology it actually wants.
pub(crate) fn from_manifold(manifold: &Manifold) -> TriMesh {
    let positions = manifold
        .ps
        .iter()
        .map(|p| Point3::new(p.x, p.y, p.z))
        .collect();
    let indices = manifold
        .get_indices()
        .iter()
        .flat_map(|t| [t.x as u32, t.y as u32, t.z as u32])
        .collect();
    TriMesh::new(positions, indices)
}
