//! Named per-vertex data carried alongside positions.
//!
//! A mesh often arrives with more per-vertex information than geometry:
//! material ids, texture coordinates, analysis results, source entity
//! handles. That data is not decoration -- losing it silently is how a
//! quantity takeoff ends up unable to say which wall a triangle came from.
//!
//! Two things make this hard, and both are modelled here rather than
//! wished away:
//!
//! - Not every value can be blended. Interpolating a material id at a new
//!   vertex invents a material that was never authored, so a channel must
//!   declare whether blending is even meaningful.
//! - Some operations genuinely cannot preserve a channel. A boolean cut
//!   creates vertices with no preimage in either operand. Reporting that
//!   honestly beats fabricating a plausible value.

use axiolid_core::Scalar;

/// How a channel's values may be combined when a new vertex appears.
///
/// This is a property of the DATA, not of the operation. An operation asks
/// the channel what is permissible; it does not decide on the channel's
/// behalf.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Blend {
    /// Values are continuous and may be linearly interpolated.
    ///
    /// Appropriate for texture coordinates, temperatures, displacements.
    Linear,
    /// Values are labels. A new vertex takes the value of the nearest
    /// existing one; averaging two ids would invent a third that names
    /// nothing.
    Nearest,
    /// Values cannot be derived for a vertex that did not exist before.
    ///
    /// The channel is dropped rather than guessed at.
    None,
}

/// Per-vertex values under a caller-chosen name.
///
/// `values.len()` must equal the mesh's vertex count. The invariant is
/// checked by [`crate::TriMesh::validate_structure`] rather than enforced
/// at construction, matching how this crate treats dirty imported data:
/// representable, then validated at a trust boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct AttributeChannel {
    /// Caller-chosen identifier, unique within a mesh.
    pub name: String,
    /// One scalar tuple per vertex, `width` entries each.
    pub values: Vec<Scalar>,
    /// Number of scalars per vertex. `2` for a UV, `1` for an id.
    pub width: usize,
    /// How values may be combined when a vertex is created.
    pub blend: Blend,
}

impl AttributeChannel {
    /// Build a channel from flat values.
    pub fn new(name: impl Into<String>, values: Vec<Scalar>, width: usize, blend: Blend) -> Self {
        Self {
            name: name.into(),
            values,
            width,
            blend,
        }
    }

    /// Number of vertices this channel covers.
    ///
    /// Returns `0` for a zero width rather than dividing by it, so a
    /// malformed channel is inspectable instead of panicking.
    pub fn vertex_count(&self) -> usize {
        if self.width == 0 {
            return 0;
        }
        self.values.len() / self.width
    }

    /// The tuple for one vertex, or `None` when out of range.
    pub fn get(&self, vertex: usize) -> Option<&[Scalar]> {
        if self.width == 0 {
            return None;
        }
        let start = vertex.checked_mul(self.width)?;
        self.values.get(start..start + self.width)
    }
}

/// What happened to a channel during an operation.
///
/// Reported rather than rejected, matching how `BooleanEvidence` treats a
/// no-op cut: the caller is told, and decides whether it matters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AttributeFate {
    /// Every vertex kept its original value.
    Preserved,
    /// New vertices received values derived under the channel's blend rule.
    Interpolated,
    /// The channel was not carried through, with the reason why.
    Dropped(DropReason),
}

/// Why a channel could not be carried through an operation.
///
/// Non-exhaustive: new operations find new reasons, and adding one must not
/// break every caller that reports them.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DropReason {
    /// The operation created vertices and the channel forbids derivation.
    ///
    /// Not a failure: [`Blend::None`] is the channel stating that an
    /// invented value would be worse than an absent one.
    NotBlendable,
    /// The provider does not carry attributes through this operation.
    ///
    /// Distinct from [`Self::NotBlendable`]: the data could in principle
    /// have survived, but this implementation does not preserve it. Naming
    /// the provider's limit separately keeps a capability gap from reading
    /// as a property of the data.
    ProviderLimitation,
    /// Vertices the operation merged carried different values.
    ///
    /// A per-vertex channel holds one value per position, so a seam -- two
    /// coincident vertices with different UVs, say -- cannot survive a weld
    /// without keeping one side's value for both. Dropping by name is the
    /// honest answer; a corner-indexed channel is the lossless one.
    ConflictingValues,
}
