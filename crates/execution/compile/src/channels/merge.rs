//! `Collection`: concatenate members, merging channels by name.

use axiolid_mesh::{AttributeChannel, AttributeFate, DropReason, NormalAttribute, TriMesh};

use super::{Built, Fates};

/// Concatenate `members` into one mesh.
///
/// Positions and indices append with a vertex offset. A channel is merged
/// across every member that defines it:
///
/// - all per-vertex and present in every member: stays per-vertex;
/// - otherwise it becomes corner-indexed. A per-vertex member's corners
///   point at its own values (offset into the merged pool); a member without
///   the channel contributes [`AttributeChannel::UNMAPPED`] triangles. This
///   is lossless, so the channel is not dropped just because one item of a
///   product is untextured.
///
/// Members that define the same name with a different width or blend cannot
/// share one channel; it is dropped and reported
/// [`DropReason::IncompatibleChannels`].
///
/// Normals are kept only if every member has them (there is no "unmapped"
/// normal), converted to corner-indexed when any member's are.
pub(crate) fn merge(members: &[&Built]) -> Built {
    let mut mesh = TriMesh::default();
    let mut fates = Fates::default();
    // Vertex and triangle offset of each member in the merged mesh.
    let mut vertex_offsets = Vec::with_capacity(members.len());
    for member in members {
        let offset = mesh.positions.len() as u32;
        vertex_offsets.push(offset);
        mesh.positions.extend_from_slice(&member.mesh.positions);
        mesh.indices
            .extend(member.mesh.indices.iter().map(|&i| i + offset));
        fates.absorb(&member.fates);
    }

    for name in channel_names(members) {
        match merge_channel(members, &name) {
            Ok(channel) => mesh.attributes.push(channel),
            Err(reason) => fates.record(&name, AttributeFate::Dropped(reason)),
        }
    }
    mesh.normals = merge_normals(members, &vertex_offsets);
    Built { mesh, fates }
}

/// Every channel name, in first-seen order.
fn channel_names(members: &[&Built]) -> Vec<String> {
    let mut names: Vec<String> = Vec::new();
    for member in members {
        for channel in &member.mesh.attributes {
            if !names.contains(&channel.name) {
                names.push(channel.name.clone());
            }
        }
    }
    names
}

fn find<'m>(member: &'m Built, name: &str) -> Option<&'m AttributeChannel> {
    member.mesh.attributes.iter().find(|c| c.name == name)
}

/// Merge one named channel across all members.
fn merge_channel(members: &[&Built], name: &str) -> Result<AttributeChannel, DropReason> {
    let present: Vec<&AttributeChannel> = members.iter().filter_map(|m| find(m, name)).collect();
    let first = present[0];
    if present
        .iter()
        .any(|c| c.width != first.width || c.blend != first.blend)
    {
        return Err(DropReason::IncompatibleChannels);
    }
    let width = first.width;

    let everywhere = present.len() == members.len();
    if everywhere && present.iter().all(|c| !c.is_corner_indexed()) {
        let values = present
            .iter()
            .flat_map(|c| c.values.iter().copied())
            .collect();
        return Ok(AttributeChannel::new(name, values, width, first.blend));
    }

    let mut values = Vec::new();
    let mut corners = Vec::new();
    for member in members {
        let corner_count = member.mesh.indices.len();
        let Some(channel) = find(member, name) else {
            corners.extend(std::iter::repeat_n(
                AttributeChannel::UNMAPPED,
                corner_count,
            ));
            continue;
        };
        let base = (values.len() / width.max(1)) as u32;
        values.extend_from_slice(&channel.values);
        // Per-vertex: a corner's value is its vertex's. Corner-indexed: its
        // own entry. Either way, rebase into the merged pool and leave
        // UNMAPPED as it is.
        let own = channel
            .corner_indices
            .as_deref()
            .unwrap_or(&member.mesh.indices);
        corners.extend(own.iter().map(|&slot| {
            if slot == AttributeChannel::UNMAPPED {
                slot
            } else {
                slot + base
            }
        }));
    }
    Ok(AttributeChannel::corner_indexed(
        name,
        values,
        width,
        first.blend,
        corners,
    ))
}

/// Merge normals, or `None` when any member lacks them.
fn merge_normals(members: &[&Built], vertex_offsets: &[u32]) -> Option<NormalAttribute> {
    let all: Vec<&NormalAttribute> = members
        .iter()
        .map(|m| m.mesh.normals.as_ref())
        .collect::<Option<_>>()?;
    if all.is_empty() {
        return None;
    }
    let values = all.iter().flat_map(|n| n.values.iter().copied()).collect();
    if all.iter().all(|n| n.indices.is_none()) {
        // Per-vertex normals line up with positions, which concatenated in
        // the same order.
        return Some(NormalAttribute {
            values,
            indices: None,
        });
    }
    let mut indices = Vec::new();
    let mut base = 0u32;
    for ((member, normals), _) in members.iter().zip(&all).zip(vertex_offsets) {
        let own = normals.indices.as_deref().unwrap_or(&member.mesh.indices);
        indices.extend(own.iter().map(|&slot| slot + base));
        base += normals.values.len() as u32;
    }
    Some(NormalAttribute {
        values,
        indices: Some(indices),
    })
}
