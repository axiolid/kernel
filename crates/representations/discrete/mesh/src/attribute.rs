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

/// Per-vertex or per-corner values under a caller-chosen name.
///
/// Per-vertex (the default, `corner_indices: None`): one tuple per position.
///
/// Per-corner (`corner_indices: Some`): one index into `values` for every
/// triangle corner, mirroring [`crate::NormalAttribute::indices`]. This is
/// how source formats store texture coordinates: a box corner shared by
/// three faces carries a different UV in each (#112). Positions stay shared,
/// so the mesh's adjacency and closure do not depend on attribute seams.
///
/// Invariants are checked by [`crate::TriMesh::validate_structure`] rather
/// than enforced at construction, matching how this crate treats dirty
/// imported data: representable, then validated at a trust boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct AttributeChannel {
    /// Caller-chosen identifier, unique within a mesh.
    pub name: String,
    /// Scalar tuples, `width` entries each: one per vertex, or the pool the
    /// corner indices point into.
    pub values: Vec<Scalar>,
    /// Number of scalars per vertex. `2` for a UV, `1` for an id.
    pub width: usize,
    /// How values may be combined when a vertex is created.
    pub blend: Blend,
    /// Optional per-corner tuple indices, one per entry of the mesh's
    /// `indices`. `None` means the channel is per-vertex.
    ///
    /// A triangle whose three entries are all [`Self::UNMAPPED`] carries no
    /// value -- real files texture only some faces. A triangle mixing
    /// mapped and unmapped corners is invalid: it would have a value at
    /// some corners and nothing to interpolate towards at the others.
    pub corner_indices: Option<Vec<u32>>,
}

impl AttributeChannel {
    /// Corner-index marker for a triangle that carries no value.
    pub const UNMAPPED: u32 = u32::MAX;

    /// Build a per-vertex channel from flat values.
    pub fn new(name: impl Into<String>, values: Vec<Scalar>, width: usize, blend: Blend) -> Self {
        Self {
            name: name.into(),
            values,
            width,
            blend,
            corner_indices: None,
        }
    }

    /// Build a per-corner channel: `corner_indices[c]` selects the tuple for
    /// triangle corner `c`, or is [`Self::UNMAPPED`].
    pub fn corner_indexed(
        name: impl Into<String>,
        values: Vec<Scalar>,
        width: usize,
        blend: Blend,
        corner_indices: Vec<u32>,
    ) -> Self {
        Self {
            corner_indices: Some(corner_indices),
            ..Self::new(name, values, width, blend)
        }
    }

    /// Whether values are addressed per triangle corner.
    pub fn is_corner_indexed(&self) -> bool {
        self.corner_indices.is_some()
    }

    /// Number of tuples in `values`.
    ///
    /// Returns `0` for a zero width rather than dividing by it, so a
    /// malformed channel is inspectable instead of panicking.
    pub fn value_count(&self) -> usize {
        if self.width == 0 {
            return 0;
        }
        self.values.len() / self.width
    }

    /// Number of vertices a per-vertex channel covers.
    ///
    /// The same as [`Self::value_count`]; kept for per-vertex callers, where
    /// a tuple IS a vertex.
    pub fn vertex_count(&self) -> usize {
        self.value_count()
    }

    /// The tuple at a `values` index, or `None` when out of range.
    ///
    /// For a per-vertex channel the index is a vertex. For a per-corner one
    /// it is an entry of `corner_indices`; use [`Self::at_corner`] to go
    /// from a mesh corner to its value.
    pub fn get(&self, vertex: usize) -> Option<&[Scalar]> {
        if self.width == 0 {
            return None;
        }
        let start = vertex.checked_mul(self.width)?;
        self.values.get(start..start + self.width)
    }

    /// The tuple at triangle corner `corner` of a mesh whose index buffer is
    /// `mesh_indices`, whichever way the channel is addressed.
    ///
    /// `None` for an unmapped corner or any out-of-range index. That makes
    /// "no value here" a result the caller has to handle, never a zero.
    pub fn at_corner(&self, mesh_indices: &[u32], corner: usize) -> Option<&[Scalar]> {
        let slot = match &self.corner_indices {
            Some(corners) => *corners.get(corner)?,
            None => *mesh_indices.get(corner)?,
        };
        if slot == Self::UNMAPPED {
            return None;
        }
        self.get(slot as usize)
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

impl AttributeFate {
    /// The fate of a channel that went through `self`, then `next`.
    ///
    /// For composed operations (batches, symmetric difference): once dropped
    /// always dropped, and the FIRST reason is kept -- it is the step that
    /// lost the data. Otherwise any interpolation makes the whole
    /// interpolated; only preserved-then-preserved stays preserved.
    #[must_use]
    pub fn then(self, next: Self) -> Self {
        match (self, next) {
            (Self::Dropped(reason), _) | (_, Self::Dropped(reason)) => Self::Dropped(reason),
            (Self::Preserved, Self::Preserved) => Self::Preserved,
            _ => Self::Interpolated,
        }
    }
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
