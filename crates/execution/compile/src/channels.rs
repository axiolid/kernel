//! Keep attribute channels and normals with the geometry through graph
//! compilation (#115).
//!
//! An `Instance` moves triangles and a `Collection` concatenates them; neither
//! creates a surface point, so every value survives unchanged and nothing
//! here derives one. The only loss is a channel name that two inputs define
//! incompatibly (different width or blend): that is dropped by name and
//! reported, never guessed at.

mod boolean;
mod merge;
mod transform;

pub(crate) use boolean::after_boolean;
pub(crate) use merge::merge;
pub(crate) use transform::transform;

use axiolid_mesh::{AttributeFate, TriMesh};

/// A node's mesh and what happened to each channel on the way to it.
#[derive(Debug, Clone, Default)]
pub(crate) struct Built {
    pub mesh: TriMesh,
    pub fates: Fates,
}

impl Built {
    /// A mesh made here from source data: every channel it carries is
    /// original.
    pub fn leaf(mesh: TriMesh) -> Self {
        let mut fates = Fates::default();
        for channel in &mesh.attributes {
            fates.record(&channel.name, AttributeFate::Preserved);
        }
        Self { mesh, fates }
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
