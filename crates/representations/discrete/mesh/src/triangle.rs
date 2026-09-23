//! Indexed triangle mesh representation.

use axiolid_core::{Aabb, Point3, Vec3};

use crate::attribute::AttributeChannel;
use crate::MeshValidationError;

/// Optional independently indexed vertex normals.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NormalAttribute {
    /// Normal values.
    pub values: Vec<Vec3>,
    /// Optional corner indices; absent means normals align with positions.
    pub indices: Option<Vec<u32>>,
}

/// Indexed triangle mesh in local coordinates.
///
/// Dirty source geometry is representable. Call [`TriMesh::validate_structure`]
/// at trust boundaries; manifold validation is a separate, more expensive pass.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TriMesh {
    /// Vertex positions.
    pub positions: Vec<Point3>,
    /// Triangle corner indices, three entries per triangle.
    pub indices: Vec<u32>,
    /// Optional normals preserved from the source.
    pub normals: Option<NormalAttribute>,
    /// Named per-vertex channels carried alongside positions.
    ///
    /// Empty by default: a mesh with no extra data pays nothing for this.
    pub attributes: Vec<AttributeChannel>,
}

impl TriMesh {
    /// Construct a position/index mesh without normals.
    pub fn new(positions: Vec<Point3>, indices: Vec<u32>) -> Self {
        Self {
            positions,
            indices,
            normals: None,
            attributes: Vec::new(),
        }
    }

    /// Number of complete triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Bounds of all positions.
    pub fn bounds(&self) -> Aabb {
        let mut bounds = Aabb::default();
        for &position in &self.positions {
            bounds.extend(position);
        }
        bounds
    }

    /// Cheap index-buffer and attribute-size validation.
    pub fn validate_structure(&self) -> Result<(), MeshValidationError> {
        if self.indices.len() % 3 != 0 {
            return Err(MeshValidationError::IncompleteTriangle {
                index_count: self.indices.len(),
            });
        }
        if let Some(&index) = self
            .indices
            .iter()
            .find(|&&index| index as usize >= self.positions.len())
        {
            return Err(MeshValidationError::PositionIndexOutOfRange {
                index,
                position_count: self.positions.len(),
            });
        }
        if let Some(normals) = &self.normals {
            if let Some(indices) = &normals.indices {
                if indices.len() != self.indices.len() {
                    return Err(MeshValidationError::NormalIndexCount {
                        expected: self.indices.len(),
                        actual: indices.len(),
                    });
                }
                if let Some(&index) = indices
                    .iter()
                    .find(|&&index| index as usize >= normals.values.len())
                {
                    return Err(MeshValidationError::NormalIndexOutOfRange {
                        index,
                        normal_count: normals.values.len(),
                    });
                }
            } else if normals.values.len() != self.positions.len() {
                return Err(MeshValidationError::NormalCount {
                    expected: self.positions.len(),
                    actual: normals.values.len(),
                });
            }
        }
        for (position, channel) in self.attributes.iter().enumerate() {
            if channel.width == 0 {
                return Err(MeshValidationError::AttributeZeroWidth {
                    name: channel.name.clone(),
                });
            }
            // Checked before the count so a duplicate is named as such
            // rather than as whichever length happens to disagree first.
            if self.attributes[..position]
                .iter()
                .any(|earlier| earlier.name == channel.name)
            {
                return Err(MeshValidationError::AttributeDuplicateName {
                    name: channel.name.clone(),
                });
            }
            if let Some(corners) = &channel.corner_indices {
                validate_corner_channel(channel, corners, self.indices.len())?;
            } else if channel.vertex_count() != self.positions.len()
                || channel.values.len() % channel.width != 0
            {
                return Err(MeshValidationError::AttributeCount {
                    name: channel.name.clone(),
                    expected: self.positions.len(),
                    actual: channel.vertex_count(),
                });
            }
        }
        Ok(())
    }

    /// Whether cheap structural validation succeeds.
    pub fn is_structurally_valid(&self) -> bool {
        self.validate_structure().is_ok()
    }

    /// Triangles in deterministic index-buffer order.
    pub fn triangles(&self) -> impl ExactSizeIterator<Item = [u32; 3]> + '_ {
        self.indices
            .chunks_exact(3)
            .map(|triangle| [triangle[0], triangle[1], triangle[2]])
    }
}

/// Check one corner-indexed channel against the mesh's corner count.
///
/// Order matters for the error a caller sees: shape first (whole tuples,
/// one entry per corner), then each entry's range, then the per-triangle
/// all-or-nothing rule. A malformed buffer is named as malformed rather
/// than as whichever triangle happens to trip first.
fn validate_corner_channel(
    channel: &AttributeChannel,
    corners: &[u32],
    corner_count: usize,
) -> Result<(), MeshValidationError> {
    if channel.values.len() % channel.width != 0 {
        return Err(MeshValidationError::AttributeRaggedValues {
            name: channel.name.clone(),
            values: channel.values.len(),
            width: channel.width,
        });
    }
    if corners.len() != corner_count {
        return Err(MeshValidationError::AttributeCornerCount {
            name: channel.name.clone(),
            expected: corner_count,
            actual: corners.len(),
        });
    }
    let value_count = channel.value_count();
    if let Some(&index) = corners
        .iter()
        .find(|&&index| index != AttributeChannel::UNMAPPED && index as usize >= value_count)
    {
        return Err(MeshValidationError::AttributeCornerIndexOutOfRange {
            name: channel.name.clone(),
            index,
            value_count,
        });
    }
    for (triangle, entries) in corners.chunks_exact(3).enumerate() {
        let unmapped = entries
            .iter()
            .filter(|&&index| index == AttributeChannel::UNMAPPED)
            .count();
        if unmapped != 0 && unmapped != 3 {
            return Err(MeshValidationError::AttributePartiallyMapped {
                name: channel.name.clone(),
                triangle,
            });
        }
    }
    Ok(())
}
