//! Splitting a solid by a plane, either hand-rolled or via a boolean provider.
//!
//! # Why both exist
//!
//! Decomposition needs one primitive: cut a solid with a plane and keep both
//! halves as closed solids. The hand-rolled clipper does this with no
//! dependencies, which matters because `axiolid-decompose` should be usable
//! without pulling in a boolean backend.
//!
//! But deciding what to do with a face lying ON the cut plane is a global
//! question about which side the material is on, and a boolean solver
//! answers it properly. So a caller that already has one can pass it in.
//!
//! # Why a contract and not a direct dependency
//!
//! `axiolid-decompose` is an `algorithms` crate and `boolmesh` is a
//! `providers` crate; that dependency edge is forbidden and the
//! architecture gate enforces it. Both layers may depend on `contracts`,
//! so the provider arrives as a `&dyn MeshBoolean` and any implementation
//! works. The choice is per call, not per build.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::MeshBoolean;

use crate::DecomposeError;

/// How a solid gets cut in two.
///
/// The variants are deliberately not a build-time switch: a caller may use
/// the hand-rolled clipper for one solid and a provider for the next.
pub enum Splitter<'a> {
    /// Clip in-crate, with no boolean backend.
    HandRolled,
    /// Delegate to a mesh boolean provider.
    Provider(&'a dyn MeshBoolean),
}

impl std::fmt::Debug for Splitter<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::HandRolled => f.write_str("Splitter::HandRolled"),
            Self::Provider(_) => f.write_str("Splitter::Provider"),
        }
    }
}

impl Splitter<'_> {
    /// Cut `mesh` by the plane `normal · x = offset`.
    ///
    /// Returns the part behind the plane and the part in front, each closed.
    /// `None` for a side means the plane did not separate any material there.
    pub fn split(
        &self,
        mesh: &TriMesh,
        normal: Vec3,
        offset: Scalar,
        tolerance: Tolerance,
    ) -> Result<(Option<TriMesh>, Option<TriMesh>), DecomposeError> {
        match self {
            Self::HandRolled => Ok((
                crate::clip(mesh, normal, offset, tolerance),
                crate::clip(mesh, -normal, -offset, tolerance),
            )),
            Self::Provider(provider) => {
                let behind = intersect_half_space(*provider, mesh, normal, offset, tolerance)?;
                let front = intersect_half_space(*provider, mesh, -normal, -offset, tolerance)?;
                Ok((behind, front))
            }
        }
    }
}

/// Intersect a solid with the half-space `normal · x <= offset`.
///
/// A boolean provider works on solids, not half-spaces, so the half-space
/// is realised as a box large enough to contain the mesh with room to
/// spare. Sizing it from the mesh's own extent rather than a fixed constant
/// keeps the construction scale-free: a model in millimetres and the same
/// model in kilometres both get a box that comfortably encloses them.
fn intersect_half_space(
    provider: &dyn MeshBoolean,
    mesh: &TriMesh,
    normal: Vec3,
    offset: Scalar,
    tolerance: Tolerance,
) -> Result<Option<TriMesh>, DecomposeError> {
    let Some(reach) = enclosing_reach(mesh) else {
        return Ok(None);
    };

    let tool = half_space_box(normal, offset, reach);
    let options = ExecutionOptions::new(tolerance);
    let outcome = provider
        .boolean(mesh, &tool, BooleanOperator::Intersection, &options)
        .map_err(|error| DecomposeError::SplitFailed(error.to_string()))?;

    // An empty intersection is a legitimate answer: the plane missed this
    // part entirely. It is reported as "no material on that side" rather
    // than as an error, because the caller's next move differs.
    if outcome.mesh.indices.is_empty() {
        return Ok(None);
    }
    Ok(Some(outcome.mesh))
}

/// Half-diagonal of the mesh's bounding box, plus a margin.
///
/// `None` when the mesh has no extent, which means there is nothing to cut.
fn enclosing_reach(mesh: &TriMesh) -> Option<Scalar> {
    if mesh.positions.is_empty() {
        return None;
    }
    let mut min = Point3::new(Scalar::INFINITY, Scalar::INFINITY, Scalar::INFINITY);
    let mut max = Point3::new(
        Scalar::NEG_INFINITY,
        Scalar::NEG_INFINITY,
        Scalar::NEG_INFINITY,
    );
    for p in &mesh.positions {
        min = Point3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
        max = Point3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
    }
    let span = max - min;
    if !span.is_finite() {
        return None;
    }
    let diagonal = span.length();
    if diagonal <= 0.0 {
        return None;
    }
    // Doubling keeps every face of the tool clear of the subject, so the
    // only surface the boolean has to resolve is the cut plane itself.
    Some(diagonal * 2.0)
}

/// A box filling `normal · x <= offset` out to `reach` in every direction.
fn half_space_box(normal: Vec3, offset: Scalar, reach: Scalar) -> TriMesh {
    let unit = normal.normalize();
    // Any orthonormal pair spanning the plane; the box's shape does not
    // depend on which, only its orientation about the normal.
    let seed = if unit.x.abs() <= unit.y.abs() && unit.x.abs() <= unit.z.abs() {
        Vec3::X
    } else if unit.y.abs() <= unit.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    let u = unit.cross(seed).normalize();
    let v = unit.cross(u);

    // Front face sits ON the cut plane; the box extends backwards from it.
    let centre = unit * offset;
    let corner = |du: Scalar, dv: Scalar, dn: Scalar| centre + u * du + v * dv + unit * dn;

    let positions = vec![
        corner(-reach, -reach, 0.0),
        corner(reach, -reach, 0.0),
        corner(reach, reach, 0.0),
        corner(-reach, reach, 0.0),
        corner(-reach, -reach, -reach * 2.0),
        corner(reach, -reach, -reach * 2.0),
        corner(reach, reach, -reach * 2.0),
        corner(-reach, reach, -reach * 2.0),
    ];
    // Wound outward: the front face (on the cut plane) points along `unit`.
    let indices = vec![
        0, 1, 2, 0, 2, 3, // front, on the plane
        4, 6, 5, 4, 7, 6, // back
        0, 4, 5, 0, 5, 1, // sides
        1, 5, 6, 1, 6, 2, //
        2, 6, 7, 2, 7, 3, //
        3, 7, 4, 3, 4, 0, //
    ];
    TriMesh::new(positions, indices)
}
