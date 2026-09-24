//! `Instance`: move the triangles, keep every value.

use axiolid_core::Transform3;
use axiolid_mesh::{NormalAttribute, TriMesh};

use super::Built;

/// Apply `transform` to a built mesh.
///
/// Positions take the point transform. Normals take the inverse transpose of
/// its linear part, renormalised: the point transform itself would skew them
/// under non-uniform scale. Attribute VALUES are untouched -- they are
/// texture coordinates, labels, measurements, not geometry -- and so are
/// their fates.
///
/// A mirroring transform reverses winding, so each triangle's corners 1 and
/// 2 are swapped to keep the result outward. Everything addressed per corner
/// (corner-indexed channels and normals) is swapped with them, or corner 1
/// would read corner 2's value.
pub(crate) fn transform(built: &Built, transform: Transform3) -> Built {
    let source = &built.mesh;
    let positions = source
        .positions
        .iter()
        .map(|&p| transform.transform_point3(p))
        .collect();
    let mut mesh = TriMesh::new(positions, source.indices.clone());
    mesh.attributes = source.attributes.clone();

    if let Some(normals) = &source.normals {
        let inverse_transpose = transform.matrix3.inverse().transpose();
        mesh.normals = Some(NormalAttribute {
            values: normals
                .values
                .iter()
                .map(|&n| (inverse_transpose * n).normalize_or_zero())
                .collect(),
            indices: normals.indices.clone(),
        });
    }

    if transform.matrix3.determinant() < 0.0 {
        swap_corners(&mut mesh.indices);
        for channel in &mut mesh.attributes {
            if let Some(corners) = channel.corner_indices.as_mut() {
                swap_corners(corners);
            }
        }
        if let Some(indices) = mesh.normals.as_mut().and_then(|n| n.indices.as_mut()) {
            swap_corners(indices);
        }
    }

    Built {
        mesh,
        fates: built.fates.clone(),
        closure: built.closure,
    }
}

/// Swap corners 1 and 2 of every triangle in a per-corner buffer.
fn swap_corners(buffer: &mut [u32]) {
    for triangle in buffer.chunks_exact_mut(3) {
        triangle.swap(1, 2);
    }
}
