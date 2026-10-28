//! Refusals that name the mesh at fault.

use core::fmt;

use axiolid_contracts::GeomError;

use crate::MeshId;

/// A refused hidden-line request: the [`GeomError`] and, when one input mesh
/// is at fault, which.
///
/// A consumer that drops elements one at a time reads [`mesh`](Self::mesh)
/// instead of parsing text. `mesh` is `None` for refusals of the whole
/// request: limits, frame, near/far, crease threshold, region, an empty mesh
/// list, memory, cancellation and a missing provider. The inner message does
/// not repeat the mesh id; [`Display`](fmt::Display) prefixes it.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct HiddenLineError {
    /// The input mesh the refusal is about, if it is about one.
    pub mesh: Option<MeshId>,
    /// Why the request was refused.
    pub error: GeomError,
}

impl HiddenLineError {
    /// Refusal of the whole request.
    pub const fn request(error: GeomError) -> Self {
        Self { mesh: None, error }
    }

    /// Refusal caused by one input mesh.
    pub const fn mesh(mesh: MeshId, error: GeomError) -> Self {
        Self {
            mesh: Some(mesh),
            error,
        }
    }
}

impl From<GeomError> for HiddenLineError {
    fn from(error: GeomError) -> Self {
        Self::request(error)
    }
}

impl fmt::Display for HiddenLineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.mesh {
            Some(MeshId(id)) => write!(f, "mesh {id}: {}", self.error),
            None => fmt::Display::fmt(&self.error, f),
        }
    }
}

impl std::error::Error for HiddenLineError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;

    use super::*;

    #[test]
    fn display_prefixes_the_mesh_only_when_there_is_one() {
        let inner = GeomError::NotManifold("edge 0-1 has 3 triangles".into());
        let plain = HiddenLineError::request(inner.clone());
        assert_eq!(plain.to_string(), inner.to_string());
        let named = HiddenLineError::mesh(MeshId(7), inner.clone());
        assert_eq!(named.to_string(), format!("mesh 7: {inner}"));
    }

    #[test]
    fn conversion_and_source_keep_the_geometry_error() {
        let error: HiddenLineError = GeomError::Cancelled.into();
        assert_eq!(error.mesh, None);
        assert_eq!(error.error, GeomError::Cancelled);
        let source = error.source().expect("source").to_string();
        assert_eq!(source, GeomError::Cancelled.to_string());
    }
}
