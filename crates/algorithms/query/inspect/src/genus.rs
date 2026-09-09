//! Topological genus from the Euler characteristic.
//!
//! For a closed orientable surface, V - E + F = 2 - 2g. The formula is only
//! meaningful on a closed two-manifold, so this refuses anything else
//! rather than returning a number the caller would have no way to distrust.

use axiolid_mesh::{component_count, EdgeAdjacency, TriMesh};
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
    /// The mesh has several components, so a single genus does not describe it.
    ///
    /// `chi = 2c - 2g` still holds for the whole surface, but the resulting
    /// `g` is a total across components and is negative whenever `c > g + 1`
    /// -- a cube with eight spherical cavities has `c = 9`, `chi = 18` and a
    /// total genus of `-8`. Reporting that as a `u32` is impossible, and
    /// reporting `0` would be indistinguishable from a genuine sphere, so
    /// this refuses and hands back what it measured.
    #[error(
        "genus requires a single connected component; found {components} \
         (Euler characteristic {characteristic}). Use `decompose` and take \
         the genus of each component separately."
    )]
    MultipleComponents {
        /// How many connected components the mesh has.
        components: usize,
        /// The characteristic of the whole surface, for callers that want it.
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
///
/// Refuses a mesh with more than one connected component for the same
/// reason. `chi = 2 - 2g` is the single-component form of `chi = 2c - 2g`;
/// applying it when `c > 1` yields a negative genus that cannot be
/// represented, and clamping that into `0` would report a solid full of
/// cavities as a sphere. Split with [`axiolid_mesh::decompose`] and call
/// this per component.
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

    // The general form is chi = 2c - 2g. Everything below assumes c == 1,
    // so a multi-component mesh must be refused BEFORE the arithmetic --
    // otherwise `doubled` goes negative and the `u32` conversion silently
    // clamps a cavity-filled solid to the same answer as a sphere.
    let components = component_count(mesh);
    if components != 1 {
        return Err(GenusError::MultipleComponents {
            components,
            characteristic,
        });
    }

    // chi = 2 - 2g, so g = (2 - chi) / 2. An odd characteristic means the
    // input is not a closed orientable surface after all.
    let doubled = 2 - characteristic;
    if doubled % 2 != 0 {
        return Err(GenusError::NotOrientable { characteristic });
    }

    // With c == 1 and chi even, `doubled / 2` is non-negative for every
    // closed orientable surface, so this conversion cannot silently clamp.
    // A failure here would mean one of the guards above was wrong, and
    // saying so is better than inventing a genus.
    u32::try_from(doubled / 2).map_err(|_| GenusError::NotOrientable { characteristic })
}
