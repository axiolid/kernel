//! Laplacian smoothing with an explicitly fixed boundary.
//!
//! Smoothing moves every vertex toward the average of its neighbours. That
//! is unconditionally destructive at a boundary: an open edge has neighbours
//! on one side only, so the average pulls it inward and the mesh shrinks
//! away from its own border.
//!
//! This module therefore treats the boundary as data, not as a special case
//! to be handled later. `fix_boundary` is on by default and a fixed vertex
//! is left BIT-IDENTICAL, not merely close: a caller comparing a smoothed
//! border against the original is asking an exactness question and deserves
//! an exact answer.

use std::collections::BTreeSet;

use axiolid_core::{Point3, Scalar};
use axiolid_mesh::{EdgeAdjacency, TriMesh};

use crate::{RefineError, SmoothReport};

/// How strongly each pass pulls a vertex toward its neighbourhood average.
///
/// A factor of 1 replaces the vertex with the average outright, which is
/// stable for a single pass but oscillates across several on an irregular
/// mesh. Values below 1 damp that.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmoothOptions {
    /// Relaxation factor in `(0, 1]`.
    pub factor: Scalar,
    /// Number of passes.
    pub passes: u32,
    /// Whether boundary vertices stay exactly where they are.
    pub fix_boundary: bool,
}

impl Default for SmoothOptions {
    fn default() -> Self {
        Self {
            factor: 0.5,
            passes: 1,
            fix_boundary: true,
        }
    }
}

/// Laplacian-smooth a triangle mesh.
///
/// # Errors
///
/// Refuses a ragged index buffer, an out-of-range index, and a relaxation
/// factor outside `(0, 1]`.
pub fn smooth(
    mesh: &TriMesh,
    options: SmoothOptions,
) -> Result<(TriMesh, SmoothReport), RefineError> {
    if mesh.indices.len() % 3 != 0 {
        return Err(RefineError::RaggedIndices(mesh.indices.len()));
    }
    let vertex_count = mesh.positions.len();
    for (triangle, chunk) in mesh.indices.chunks_exact(3).enumerate() {
        for &index in chunk {
            if index as usize >= vertex_count {
                return Err(RefineError::IndexOutOfRange(triangle, index));
            }
        }
    }
    if !(options.factor > 0.0 && options.factor <= 1.0) {
        return Err(RefineError::InvalidTarget(options.factor));
    }

    // Adjacency is derived once: connectivity never changes here, only
    // positions, so rebuilding it per pass would be wasted work. Boundary
    // vertices are pinned, since moving them would shrink the surface.
    let adjacency = EdgeAdjacency::build(mesh);
    let boundary: BTreeSet<u32> = adjacency.boundary_vertices().into_iter().collect();
    let neighbours = adjacency.vertex_neighbours();

    let mut positions = mesh.positions.clone();
    let mut moved = 0usize;
    let mut max_movement: Scalar = 0.0;

    for _ in 0..options.passes {
        let source = positions.clone();
        for (index, position) in positions.iter_mut().enumerate() {
            let vertex = index as u32;
            if options.fix_boundary && boundary.contains(&vertex) {
                continue;
            }
            let adjacent = &neighbours[index];
            if adjacent.is_empty() {
                continue;
            }
            let mut sum = Point3::ZERO;
            for &other in adjacent {
                sum += source[other as usize];
            }
            let average = sum / (adjacent.len() as Scalar);
            let target = source[index] + (average - source[index]) * options.factor;
            let movement = (target - source[index]).length();
            if movement > 0.0 {
                moved += 1;
                max_movement = max_movement.max(movement);
            }
            *position = target;
        }
    }

    let attribute_fates = crate::carry_attributes(mesh);
    let out = TriMesh {
        positions,
        indices: mesh.indices.clone(),
        normals: None,
        attributes: Vec::new(),
    };
    let report = SmoothReport {
        vertices_moved: moved,
        boundary_vertices: boundary.len(),
        max_movement,
        attribute_fates,
    };
    Ok((out, report))
}
