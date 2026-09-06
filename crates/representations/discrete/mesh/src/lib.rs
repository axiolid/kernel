#![forbid(unsafe_code)]

//! Mesh representations and zero-copy interoperability views.
//!
//! N-gons remain [`PolygonMesh`] until explicit triangulation. [`TriMesh`] is
//! the compact exchange type for render, spatial, and mesh-kernel algorithms.

pub mod adjacency;
pub mod attribute;
pub mod audit;
pub mod component;
pub mod error;
pub mod polygon;
pub mod triangle;
pub mod view;

pub use adjacency::{EdgeAdjacency, EdgeKey, EdgeUse};
pub use attribute::{AttributeChannel, AttributeFate, Blend, DropReason};
pub use audit::{audit_mesh, audit_mesh_scratch_bytes, try_audit_mesh, MeshAuditError, MeshHealth};
pub use component::{component_count, compose, decompose};
pub use error::MeshValidationError;
pub use polygon::{PolygonFace, PolygonMesh};
pub use triangle::{NormalAttribute, TriMesh};
pub use view::{MeshView, TriangleMeshView};

#[cfg(test)]
mod tests {
    use axiolid_core::Vec3;

    use super::*;

    #[test]
    fn triangle_iteration_preserves_index_order() {
        let mesh = TriMesh::new(
            vec![Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::Z],
            vec![0, 1, 2, 0, 2, 3],
        );
        assert_eq!(mesh.triangles().collect::<Vec<_>>(), [[0, 1, 2], [0, 2, 3]]);
        assert!(mesh.validate_structure().is_ok());
    }

    #[test]
    fn invalid_index_is_a_structured_error() {
        let mesh = TriMesh::new(vec![Vec3::ZERO], vec![0, 1, 2]);
        assert!(matches!(
            mesh.validate_structure(),
            Err(MeshValidationError::PositionIndexOutOfRange { index: 1, .. })
        ));
    }
}
