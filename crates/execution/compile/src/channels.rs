//! Keep attribute channels and normals with the geometry through graph
//! compilation (#115).
//!
//! An `Instance` moves triangles and a `Collection` concatenates them; neither
//! creates a surface point, so every value survives unchanged and nothing
//! here derives one. The only loss is a channel name that two inputs define
//! incompatibly (different width or blend): that is dropped by name and
//! reported, never guessed at.
//!
//! The compiler caches a [`Built`], never a bare `TriMesh`, so every node
//! kind has to state what it did to each channel. A new graph path goes
//! through [`transform()`], [`merge()`] or [`after_boolean`], or wraps a mesh it
//! made itself in [`Built::leaf`] / [`Built::with_closure`]. Rebuilding a
//! mesh from positions and indices on an `Instance` or `Collection` path is
//! the #115 bug: it silently drops every channel. Gate:
//! `tests/graph_channels.rs`.
//!
//! # Closure rides along too (#161)
//!
//! [`Built::leaf`] claims a solid. Anything that may not bound one goes
//! through [`Built::with_closure`]: a B-rep is a solid only when it declares
//! one, never because its shells happen to be watertight, since a surface
//! model's mesh can be closed. A collection is a solid only if every member
//! is ([`combined_closure`]), and a boolean refuses a surface operand.
//! Gate: `tests/surface_models.rs`.

mod boolean;
mod merge;
mod transform;

pub(crate) use boolean::after_boolean;
pub(crate) use merge::merge;
pub(crate) use transform::transform;

use axiolid_mesh::{AttributeFate, TriMesh};
use axiolid_mesh_compile_contract::MeshClosure;

/// A node's mesh, what happened to each channel on the way to it, whether
/// it bounds a solid (#161), and how far the exact surface may lie from it
/// (#232).
#[derive(Debug, Clone)]
pub(crate) struct Built {
    pub mesh: TriMesh,
    pub fates: Fates,
    /// Never [`MeshClosure::Unknown`] inside the compiler: every node
    /// states it.
    pub closure: MeshClosure,
    /// Unreported (unbounded) until the node's path states it with
    /// [`Built::with_deviation`]; instances, collections and booleans
    /// derive it from their operands.
    pub deviation: crate::deviation::Deviation,
}

impl Built {
    /// A solid made here from source data: every channel it carries is
    /// original.
    pub fn leaf(mesh: TriMesh) -> Self {
        Self::with_closure(mesh, MeshClosure::Solid)
    }

    /// A mesh made here from source data, with its closure.
    pub fn with_closure(mesh: TriMesh, closure: MeshClosure) -> Self {
        let mut fates = Fates::default();
        for channel in &mesh.attributes {
            fates.record(&channel.name, AttributeFate::Preserved);
        }
        Self {
            mesh,
            fates,
            closure,
            deviation: crate::deviation::Deviation::default(),
        }
    }

    /// The same node with its deviation stated.
    pub fn with_deviation(mut self, deviation: crate::deviation::Deviation) -> Self {
        self.deviation = deviation;
        self
    }
}

/// The closure of several parts together: a solid only if every part is.
pub(crate) fn combined_closure<'a>(parts: impl IntoIterator<Item = &'a Built>) -> MeshClosure {
    if parts
        .into_iter()
        .all(|part| part.closure == MeshClosure::Solid)
    {
        MeshClosure::Solid
    } else {
        MeshClosure::Surface
    }
}

/// Per-channel fates, in first-seen order, each the worst seen so far.
///
/// Worst means furthest from the source data: `Dropped` over `Interpolated`
/// over `Preserved`. A channel preserved in one member and dropped in
/// another is reported dropped, because the caller cannot trust it whole.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct Fates(Vec<(String, AttributeFate)>);

fn rank(fate: &AttributeFate) -> u8 {
    match fate {
        AttributeFate::Preserved => 0,
        AttributeFate::Interpolated => 1,
        AttributeFate::Dropped(_) => 2,
    }
}

impl Fates {
    /// Record `fate` for `name`, keeping the worse of it and any earlier one.
    pub fn record(&mut self, name: &str, fate: AttributeFate) {
        match self.0.iter_mut().find(|(seen, _)| seen == name) {
            Some((_, current)) if rank(&fate) > rank(current) => *current = fate,
            Some(_) => {}
            None => self.0.push((name.to_owned(), fate)),
        }
    }

    /// Fold another node's fates into this one.
    pub fn absorb(&mut self, other: &Fates) {
        for (name, fate) in &other.0 {
            self.record(name, fate.clone());
        }
    }

    /// The report, in first-seen order.
    pub fn into_vec(self) -> Vec<(String, AttributeFate)> {
        self.0
    }
}

#[cfg(test)]
mod tests;
