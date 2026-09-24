//! What a compilation did to the caller's attribute channels, and whether
//! the result encloses a volume.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_mesh::{AttributeFate, TriMesh};

/// Whether a compiled mesh bounds a solid (#161).
///
/// A triangle mesh looks the same either way: a closed surface model (say,
/// an IFC `IfcShellBasedSurfaceModel` shaped like a box) is watertight, so
/// its divergence sum is a finite, plausible number. It is still not a
/// volume, because the source never claimed one. This is what tells the
/// two apart.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MeshClosure {
    /// Every part of the mesh is the boundary of a solid; its enclosed
    /// volume is meaningful.
    Solid,
    /// At least one part is a surface with no solid behind it: an open
    /// shell, or a closed shell the source never declared a solid. Area is
    /// meaningful; volume is not.
    Surface,
    /// The compiler did not report it.
    Unknown,
}

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
    /// Whether the mesh bounds a solid. [`MeshClosure::Unknown`] unless the
    /// compiler reports it; read the mesh through [`Self::solid_mesh`]
    /// before taking a volume from it.
    pub closure: MeshClosure,
}

impl CompileOutcome {
    /// A mesh from a compiler that does not report channel fates.
    pub const fn untracked(mesh: TriMesh) -> Self {
        Self {
            mesh,
            attribute_fates: None,
            closure: MeshClosure::Unknown,
        }
    }

    /// A mesh with its channel report.
    pub const fn tracked(mesh: TriMesh, attribute_fates: Vec<(String, AttributeFate)>) -> Self {
        Self {
            mesh,
            attribute_fates: Some(attribute_fates),
            closure: MeshClosure::Unknown,
        }
    }

    /// The same outcome with its closure reported.
    #[must_use]
    pub const fn with_closure(mut self, closure: MeshClosure) -> Self {
        self.closure = closure;
        self
    }

    /// The mesh, only if it is known to bound a solid.
    ///
    /// Use this before measuring volume, centroid or second moments. A
    /// surface model's mesh can be closed and consistently wound, so a
    /// volume computed from it would be a plausible number with no meaning;
    /// this refuses instead of returning it (#161).
    ///
    /// # Errors
    ///
    /// [`GeomError::InvalidInput`] for [`MeshClosure::Surface`] and for
    /// [`MeshClosure::Unknown`]: an unreported closure is not a solid.
    pub fn solid_mesh(&self) -> GeomResult<&TriMesh> {
        match self.closure {
            MeshClosure::Solid => Ok(&self.mesh),
            MeshClosure::Surface => Err(GeomError::InvalidInput(
                "compiled mesh is a surface model: it has area but encloses no volume".to_owned(),
            )),
            MeshClosure::Unknown => Err(GeomError::InvalidInput(
                "compiler did not report whether the mesh bounds a solid".to_owned(),
            )),
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
