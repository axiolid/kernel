//! Topological genus from the Euler characteristic.
//!
//! For a closed orientable surface, V - E + F = 2 - 2g. The formula is only
//! meaningful on a closed two-manifold, so this refuses anything else
//! rather than returning a number the caller would have no way to distrust.

use axiolid_mesh::{EdgeAdjacency, TriMesh};
use thiserror::Error;

/// Why a genus could not be computed.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum GenusError {
    /// The mesh is not a closed two-manifold, so Euler's formula does not apply.
    #[error("genus requires a closed two-manifold; found {boundary} boundary and {non_manifold} non-manifold edges")]
    NotClosedManifold {
        /// Edges used by exactly one triangle.
        boundary: usize,
        /// Edges used by more than two triangles.
        non_manifold: usize,
    },
    /// The Euler characteristic was odd, so no integer genus exists.
    #[error("Euler characteristic {characteristic} is odd; the surface is not orientable")]
    NotOrientable {
        /// The computed characteristic.
        characteristic: i64,
    },
}

/// Genus of a closed orientable triangle mesh.
///
/// A sphere or cube is 0, a torus 1, a double torus 2.
///
/// # Errors
///
/// Refuses a mesh with boundary or non-manifold edges: Euler's formula
/// assumes a closed two-manifold, and applying it anyway would produce a
/// plausible-looking integer with no meaning.
pub fn genus(mesh: &TriMesh) -> Result<u32, GenusError> {
    let adjacency = EdgeAdjacency::build(mesh);
    let boundary = adjacency.boundary_edges().count();
    let non_manifold = adjacency.non_manifold_edges().count();
    if boundary > 0 || non_manifold > 0 {
        return Err(GenusError::NotClosedManifold {
            boundary,
            non_manifold,
        });
    }

    // Only vertices actually referenced by a triangle count: an unused
    // position is stray data, not part of the surface. `EdgeAdjacency`
    // already applies that rule, so the two cannot disagree.
    let characteristic = adjacency.euler_characteristic();

    // chi = 2 - 2g, so g = (2 - chi) / 2. An odd characteristic means the
    // input is not a closed orientable surface after all.
    let doubled = 2 - characteristic;
    if doubled % 2 != 0 {
        return Err(GenusError::NotOrientable { characteristic });
    }
    Ok(u32::try_from(doubled / 2).unwrap_or(0))
}
