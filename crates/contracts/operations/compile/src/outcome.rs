//! What a compilation did to the caller's attribute channels.

use axiolid_mesh::{AttributeFate, TriMesh};

/// A compiled mesh and, when the compiler tracks it, each channel's fate.
///
/// Compilation moves and concatenates geometry, and a boolean inside the
/// graph cuts it, so a channel can reach the result unchanged, derived, or
/// not at all. `mesh.attributes` shows what survived; `attribute_fates`
/// says why the rest did not -- the #84 contract, applied to graph
/// compilation (#115).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct CompileOutcome {
    /// The compiled mesh.
    pub mesh: TriMesh,
    /// Each channel seen anywhere in the compiled subgraph, in first-seen
    /// order, with its worst fate on any path to the root.
    ///
    /// `None` when the compiler does not track channels: that is "unknown",
    /// not "nothing was dropped".
    pub attribute_fates: Option<Vec<(String, AttributeFate)>>,
}

impl CompileOutcome {
    /// A mesh from a compiler that does not report channel fates.
    pub const fn untracked(mesh: TriMesh) -> Self {
        Self {
            mesh,
            attribute_fates: None,
        }
    }

    /// A mesh with its channel report.
    pub const fn tracked(mesh: TriMesh, attribute_fates: Vec<(String, AttributeFate)>) -> Self {
        Self {
            mesh,
            attribute_fates: Some(attribute_fates),
        }
    }

    /// The fate of the channel named `name`, if the compiler reported one.
    pub fn fate(&self, name: &str) -> Option<&AttributeFate> {
        self.attribute_fates
            .as_ref()?
            .iter()
            .find(|(seen, _)| seen == name)
            .map(|(_, fate)| fate)
    }
}
