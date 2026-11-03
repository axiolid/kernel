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
        if channel.is_corner_indexed() {
            // Corner indices address the channel's own values, not
            // positions, so renumbering positions leaves them valid. This
            // is also why a seam never conflicts here: each side of it is
            // already a separate value (#112).
            kept.push(channel);
            continue;
        }
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
    let select = |corners: &mut Vec<u32>| {
        let next: Vec<u32> = kept
            .iter()
            .filter_map(|&t| corners.get(t * 3..t * 3 + 3))
            .flatten()
            .copied()
            .collect();
        *corners = next;
    };
    for_each_corner_buffer(mesh, select);
}

/// Reverse one triangle's corner data to match `indices.swap(t*3+1, t*3+2)`.
///
/// Only the corner association moves. Normal VALUES are not negated: a
/// winding repair says the triangle's order was wrong, not that its authored
/// normals were, and the healer has no evidence either way.
pub(crate) fn flip_triangle(mesh: &mut TriMesh, t: usize) {
    for_each_corner_buffer(mesh, |corners| {
        if t * 3 + 2 < corners.len() {
            corners.swap(t * 3 + 1, t * 3 + 2);
        }
    });
}

/// Apply one edit to every buffer indexed by triangle corner.
///
/// The single place that knows which buffers those are. A repair that
/// edits corners calls this, so adding a new per-corner buffer to the mesh
/// means adding it here once, not in every repair.
fn for_each_corner_buffer(mesh: &mut TriMesh, mut edit: impl FnMut(&mut Vec<u32>)) {
    if let Some(indices) = mesh.normals.as_mut().and_then(|n| n.indices.as_mut()) {
        edit(indices);
    }
    for channel in &mut mesh.attributes {
        if let Some(corners) = channel.corner_indices.as_mut() {
            edit(corners);
        }
    }
}

/// Rebuild every corner buffer so slot `i` holds what slot `source[i]` held.
///
/// A repair that splits a triangle keeps its first piece in place and
/// appends the others; this gives each new corner the value of the corner
/// it grew from. Every source slot must exist in every buffer, which
/// [`corner_buffers_aligned`] decides beforehand.
pub(crate) fn remap_corners(mesh: &mut TriMesh, source: &[usize]) {
    for_each_corner_buffer(mesh, |corners| {
        let next: Vec<u32> = source.iter().map(|&s| corners[s]).collect();
        *corners = next;
    });
}

/// Whether every corner buffer has one entry per index-buffer slot.
///
/// Read-only twin of [`for_each_corner_buffer`]; the two name the same
/// buffers and must be kept in step.
pub(crate) fn corner_buffers_aligned(mesh: &TriMesh) -> bool {
    let len = mesh.indices.len();
    let normals = mesh.normals.as_ref().and_then(|n| n.indices.as_ref());
    normals.into_iter().all(|c| c.len() == len)
        && mesh
            .attributes
            .iter()
            .filter_map(|channel| channel.corner_indices.as_ref())
            .all(|c| c.len() == len)
}

/// Whether corners `a` and `b` hold the same value in every corner-indexed
/// buffer: the same value index, or bit-equal value tuples.
///
/// An unmapped or out-of-range corner equals only a corner with the same
/// index. Read-only twin of [`for_each_corner_buffer`].
pub(crate) fn corners_agree(mesh: &TriMesh, a: usize, b: usize) -> bool {
    let bits = |tuple: &[f64]| tuple.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    if let Some(normals) = mesh.normals.as_ref() {
        if let Some(indices) = normals.indices.as_ref() {
            let (Some(&i), Some(&j)) = (indices.get(a), indices.get(b)) else {
                return false;
            };
            if i != j {
                let value = |k: u32| {
                    normals
                        .values
                        .get(k as usize)
                        .map(|n| bits(&[n.x, n.y, n.z]))
                };
                match (value(i), value(j)) {
                    (Some(x), Some(y)) if x == y => {}
                    _ => return false,
                }
            }
        }
    }
    for channel in &mesh.attributes {
        let Some(corners) = channel.corner_indices.as_ref() else {
            continue;
        };
        let (Some(&i), Some(&j)) = (corners.get(a), corners.get(b)) else {
            return false;
        };
        if i == j {
            continue;
        }
        match (channel.get(i as usize), channel.get(j as usize)) {
            (Some(x), Some(y)) if bits(x) == bits(y) => {}
            _ => return false,
        }
    }
    true
}

/// Whether the per-vertex data of vertex `m`, lying at parameter `t` on the
/// edge `low`-`high`, is what interpolating the endpoints gives.
///
/// Used before a split reuses `m` as a new corner of a triangle that did not
/// have it, so that "preserved" stays true. A `Linear` channel (and a
/// per-vertex normal) must equal `low + (high - low) * t` exactly; any other
/// blend must be bit-equal at `low`, `high` and `m`. A missing value
/// disagrees. Corner-indexed buffers are [`corners_agree`]'s concern.
pub(crate) fn vertex_data_agrees(mesh: &TriMesh, low: u32, high: u32, m: u32, t: f64) -> bool {
    let at = |v: u32| v as usize;
    if let Some(normals) = mesh.normals.as_ref().filter(|n| n.indices.is_none()) {
        let get = |v: u32| normals.values.get(at(v)).copied();
        let (Some(a), Some(b), Some(p)) = (get(low), get(high), get(m)) else {
            return false;
        };
        if a + (b - a) * t != p {
            return false;
        }
    }
    for channel in mesh.attributes.iter().filter(|c| !c.is_corner_indexed()) {
        let (Some(a), Some(b), Some(p)) = (
            channel.get(at(low)),
            channel.get(at(high)),
            channel.get(at(m)),
        ) else {
            return false;
        };
        let agrees = (0..p.len()).all(|k| {
            if channel.blend == axiolid_mesh::Blend::Linear {
                a[k] + (b[k] - a[k]) * t == p[k]
            } else {
                a[k].to_bits() == p[k].to_bits() && b[k].to_bits() == p[k].to_bits()
            }
        });
        if !agrees {
            return false;
        }
    }
    true
}
