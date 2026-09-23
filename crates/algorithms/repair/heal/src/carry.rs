//! Keep per-vertex and per-corner data in step with repaired geometry (#114).
//!
//! A repair that renumbers vertices or removes triangles must apply the same
//! edit to every buffer indexed by vertex or by corner. Doing that in one
//! place is what stops the next buffer from being forgotten -- which is how
//! #114 happened: the repairs rewrote positions and indices only.

use axiolid_mesh::{AttributeChannel, DropReason, TriMesh};

/// Whether any merged vertex disagrees with its representative.
fn conflicts<T: PartialEq>(groups: &[(u32, Vec<u32>)], get: impl Fn(u32) -> Option<T>) -> bool {
    groups
        .iter()
        .any(|(rep, dups)| dups.iter().any(|&d| get(d) != get(*rep)))
}

/// Carry per-vertex data through a weld.
///
/// `groups` pairs each surviving representative with the vertices merged
/// into it, `keep` lists the surviving original vertices in their new order,
/// and `corners` is the index buffer as it was BEFORE the weld.
pub(crate) fn weld(
    mesh: &mut TriMesh,
    groups: &[(u32, Vec<u32>)],
    keep: &[u32],
    corners: &[u32],
    dropped: &mut Vec<(String, DropReason)>,
) {
    let mut kept = Vec::with_capacity(mesh.attributes.len());
    for channel in core::mem::take(&mut mesh.attributes) {
        if conflicts(groups, |v| channel.get(v as usize)) {
            // Coincident vertices carrying different values are a seam. One
            // value per position cannot hold both, and keeping either one
            // smears the other side, so the channel is dropped by name.
            dropped.push((channel.name.clone(), DropReason::ConflictingValues));
            continue;
        }
        let mut values = Vec::with_capacity(keep.len() * channel.width);
        for &v in keep {
            if let Some(tuple) = channel.get(v as usize) {
                values.extend_from_slice(tuple);
            }
        }
        kept.push(AttributeChannel { values, ..channel });
    }
    mesh.attributes = kept;

    let Some(normals) = mesh.normals.as_mut() else {
        return;
    };
    // Corner-indexed normals index normal VALUES, not positions, so a weld
    // leaves them valid as they are.
    if normals.indices.is_some() {
        return;
    }
    if conflicts(groups, |v| normals.values.get(v as usize).copied()) {
        // A hard edge: coincident positions with different normals.
        // `NormalAttribute` already has a corner-indexed form for exactly
        // this, so switching to it is lossless: each corner keeps the normal
        // of the vertex it referenced before the weld.
        normals.indices = Some(corners.to_vec());
    } else {
        normals.values = keep
            .iter()
            .filter_map(|&v| normals.values.get(v as usize).copied())
            .collect();
    }
}

/// Keep only the listed triangles' corner data, in the given order.
pub(crate) fn keep_triangles(mesh: &mut TriMesh, kept: &[usize]) {
    if let Some(indices) = mesh.normals.as_mut().and_then(|n| n.indices.as_mut()) {
        let next: Vec<u32> = kept
            .iter()
            .filter_map(|&t| indices.get(t * 3..t * 3 + 3))
            .flatten()
            .copied()
            .collect();
        *indices = next;
    }
}

/// Reverse one triangle's corner data to match `indices.swap(t*3+1, t*3+2)`.
///
/// Only the corner association moves. Normal VALUES are not negated: a
/// winding repair says the triangle's order was wrong, not that its authored
/// normals were, and the healer has no evidence either way.
pub(crate) fn flip_triangle(mesh: &mut TriMesh, t: usize) {
    if let Some(indices) = mesh.normals.as_mut().and_then(|n| n.indices.as_mut()) {
        if t * 3 + 2 < indices.len() {
            indices.swap(t * 3 + 1, t * 3 + 2);
        }
    }
}
