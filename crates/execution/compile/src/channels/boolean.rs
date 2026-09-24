//! `Boolean`: compose the operands' upstream fates with the provider's.

use axiolid_mesh::{AttributeFate, TriMesh};
use axiolid_mesh_boolean_contract::merge_fates;

use super::{Built, Fates};

/// The boolean's result, with a report covering the whole subgraph.
///
/// The two operands are parallel inputs, so their upstream fates combine
/// worst-of ([`Fates::absorb`]). The boolean then runs AFTER them, so it
/// composes sequentially through [`merge_fates`]: a channel dropped upstream
/// keeps that first reason, and one the provider no longer mentions was lost
/// in the boolean.
///
/// `tool` is `None` for a tool built here from the subject (a bounded half
/// space), which has no history of its own.
pub(crate) fn after_boolean(
    mesh: TriMesh,
    subject: &Fates,
    tool: Option<&Fates>,
    provider: Vec<(String, AttributeFate)>,
) -> Built {
    let mut upstream = subject.clone();
    if let Some(tool) = tool {
        upstream.absorb(tool);
    }
    let upstream = upstream.into_vec();
    let mut composed = merge_fates(&upstream, provider.clone());
    // `merge_fates` walks `upstream`'s names only; a channel the provider
    // reports that nothing upstream tracked is appended as reported.
    for (name, fate) in provider {
        if !composed.iter().any(|(seen, _)| *seen == name) {
            composed.push((name, fate));
        }
    }
    let mut fates = Fates::default();
    for (name, fate) in composed {
        fates.record(&name, fate);
    }
    // Callers refuse surface operands, so a boolean result bounds a solid.
    Built {
        mesh,
        fates,
        closure: axiolid_mesh_compile_contract::MeshClosure::Solid,
    }
}
