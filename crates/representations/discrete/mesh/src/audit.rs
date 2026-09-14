//! Deterministic structural triangle-mesh audit.
//!
//! The audit reports defects instead of rejecting dirty geometry. Open meshes can
//! still support surface distance and intersection; callers that need a watertight
//! solid must check [`MeshHealth::is_closed_two_manifold`].

use std::collections::BTreeMap;
use std::fmt;

use axiolid_core::Tolerance;

use crate::TriangleMeshView;

/// Source-neutral structural facts about one triangle mesh.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MeshHealth {
    /// Number of addressable positions.
    pub positions: usize,
    /// Number of triangle records inspected.
    pub triangles: usize,
    /// Triangles with valid, finite, non-degenerate vertices.
    pub usable_triangles: usize,
    /// Number of invalid triangle-corner indices.
    pub invalid_indices: usize,
    /// Number of non-finite positions.
    pub non_finite_positions: usize,
    /// Number of triangles below the explicit area threshold.
    pub degenerate_triangles: usize,
    /// Undirected edges with exactly one usable incident triangle.
    pub boundary_edges: usize,
    /// Undirected edges with more than two usable incident triangles.
    pub non_manifold_edges: usize,
    /// Two-manifold edges whose incident faces use the same directed edge.
    /// Such a mesh is closed but has inconsistent local winding, so signed
    /// enclosed-volume reduction is not structurally trustworthy.
    pub inconsistent_winding_edges: usize,
    /// First `(triangle, source_index)` which could not address a position.
    pub first_invalid_index: Option<(usize, u64)>,
    /// First non-finite position index.
    pub first_non_finite_position: Option<usize>,
}

impl MeshHealth {
    /// Whether at least one triangle can safely support surface algorithms.
    pub fn is_surface_usable(&self) -> bool {
        self.usable_triangles > 0 && self.invalid_indices == 0 && self.non_finite_positions == 0
    }

    /// Whether the usable mesh is closed and two-manifold.
    pub fn is_closed_two_manifold(&self) -> bool {
        self.is_surface_usable()
            && self.degenerate_triangles == 0
            && self.boundary_edges == 0
            && self.non_manifold_edges == 0
            && self.inconsistent_winding_edges == 0
    }
}

/// A bounded audit could not reserve its declared edge-record scratch.
#[derive(Debug)]
pub enum MeshAuditError {
    /// `3 * triangle_count` or its byte size overflowed `usize`.
    CapacityOverflow,
    /// The allocator refused the exact edge-record reservation.
    Allocation(std::collections::TryReserveError),
}

impl fmt::Display for MeshAuditError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CapacityOverflow => formatter.write_str("mesh audit edge count overflowed"),
            Self::Allocation(error) => write!(formatter, "mesh audit allocation failed: {error}"),
        }
    }
}

impl std::error::Error for MeshAuditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::CapacityOverflow => None,
            Self::Allocation(error) => Some(error),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct EdgeRecord {
    low: u64,
    high: u64,
    direction: i8,
}

impl EdgeRecord {
    /// Filler for the counting sort's scratch buffer. Every slot is
    /// overwritten before it is read; this only avoids `unsafe`.
    const EMPTY: Self = Self {
        low: 0,
        high: 0,
        direction: 0,
    };
}

/// Exact requested scratch bytes for [`try_audit_mesh`].
///
/// The bounded implementation stores at most three fixed-size edge records per
/// source triangle. It may reserve a second buffer of the same size to run a
/// counting sort instead of a comparison sort, so the bound covers two buffers:
/// reporting only one would let a caller admit an audit this function then
/// refuses. `None` means the count overflows.
pub const fn audit_mesh_scratch_bytes(triangle_count: usize) -> Option<usize> {
    match triangle_count.checked_mul(3) {
        Some(edges) => match edges.checked_mul(std::mem::size_of::<EdgeRecord>()) {
            Some(bytes) => bytes.checked_mul(2),
            None => None,
        },
        None => None,
    }
}

#[derive(Debug, Default)]
struct EdgeSummary {
    boundary: usize,
    non_manifold: usize,
    inconsistent_winding: usize,
}

trait EdgeSink {
    fn record(&mut self, low: u64, high: u64, direction: i8);
    fn summarize(&mut self) -> EdgeSummary;
}

#[derive(Debug)]
struct VecEdgeSink {
    edges: Vec<EdgeRecord>,
    /// Second buffer for the counting sort's ping-pong. Reserved up
    /// front so `summarize` cannot fail on allocation half way through.
    scratch: Vec<EdgeRecord>,
    /// Exclusive upper bound on vertex ids, i.e. the counting-sort key
    /// space. Zero disables the counting sort.
    buckets: usize,
}

impl VecEdgeSink {
    fn try_new(triangle_count: usize, positions: usize) -> Result<Self, MeshAuditError> {
        let count = triangle_count
            .checked_mul(3)
            .ok_or(MeshAuditError::CapacityOverflow)?;
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(count)
            .map_err(MeshAuditError::Allocation)?;
        // The counting sort needs a second buffer and a counts array of
        // `positions` entries. It only pays when the key space is
        // comparable to the edge count: a mesh with few triangles over a
        // huge index space would spend more time clearing counts than
        // sorting. Falling back to the comparison sort there keeps the
        // pathological case from regressing.
        let dense = positions <= count.saturating_mul(2).max(1024);
        let fits = u32::try_from(count).is_ok();
        let mut scratch = Vec::new();
        let buckets = if dense && fits {
            match scratch.try_reserve_exact(count) {
                Ok(()) => {
                    scratch.resize(count, EdgeRecord::EMPTY);
                    positions
                }
                // Scratch is an optimisation, not a requirement: losing
                // it costs speed, not correctness.
                Err(_) => 0,
            }
        } else {
            0
        };
        Ok(Self {
            edges,
            scratch,
            buckets,
        })
    }
}

/// Group equal edge keys by counting sort on vertex ids.
///
/// Profiling attributed most of `audit_mesh` to sorting. The keys are
/// vertex indices, bounded by the position count, so a two-pass counting
/// sort replaces the comparison sort. Passes run high then low so the
/// final order is by (low, high) -- the same order the old code produced,
/// which keeps every downstream count identical rather than merely
/// grouped.
fn counting_sort_edges(edges: &mut Vec<EdgeRecord>, scratch: &mut Vec<EdgeRecord>, buckets: usize) {
    debug_assert_eq!(scratch.len(), edges.len());
    let mut counts: Vec<u32> = Vec::new();
    for pass in 0..2 {
        counts.clear();
        counts.resize(buckets + 2, 0);
        for edge in edges.iter() {
            let key = if pass == 0 { edge.high } else { edge.low } as usize;
            counts[key + 1] += 1;
        }
        for index in 0..=buckets {
            counts[index + 1] += counts[index];
        }
        for edge in edges.iter() {
            let key = if pass == 0 { edge.high } else { edge.low } as usize;
            scratch[counts[key] as usize] = *edge;
            counts[key] += 1;
        }
        std::mem::swap(edges, scratch);
    }
}

impl EdgeSink for VecEdgeSink {
    fn record(&mut self, low: u64, high: u64, direction: i8) {
        self.edges.push(EdgeRecord {
            low,
            high,
            direction,
        });
    }

    fn summarize(&mut self) -> EdgeSummary {
        if self.buckets > 0 && self.scratch.len() == self.edges.len() {
            counting_sort_edges(&mut self.edges, &mut self.scratch, self.buckets);
        } else {
            self.edges
                .sort_unstable_by_key(|edge| (edge.low, edge.high));
        }
        let mut summary = EdgeSummary::default();
        let mut start = 0;
        while start < self.edges.len() {
            let key = (self.edges[start].low, self.edges[start].high);
            let mut end = start + 1;
            let mut winding = i128::from(self.edges[start].direction);
            while end < self.edges.len() && (self.edges[end].low, self.edges[end].high) == key {
                winding += i128::from(self.edges[end].direction);
                end += 1;
            }
            match end - start {
                1 => summary.boundary += 1,
                2 if winding != 0 => summary.inconsistent_winding += 1,
                count if count > 2 => summary.non_manifold += 1,
                _ => {}
            }
            start = end;
        }
        summary
    }
}

#[derive(Debug, Default)]
struct MapEdgeSink {
    edges: BTreeMap<(u64, u64), (usize, i128)>,
}

impl EdgeSink for MapEdgeSink {
    fn record(&mut self, low: u64, high: u64, direction: i8) {
        let entry = self.edges.entry((low, high)).or_default();
        entry.0 = entry.0.saturating_add(1);
        entry.1 += i128::from(direction);
    }

    fn summarize(&mut self) -> EdgeSummary {
        EdgeSummary {
            boundary: self
                .edges
                .values()
                .filter(|&&(count, _)| count == 1)
                .count(),
            non_manifold: self.edges.values().filter(|&&(count, _)| count > 2).count(),
            inconsistent_winding: self
                .edges
                .values()
                .filter(|&&(count, winding)| count == 2 && winding != 0)
                .count(),
        }
    }
}

/// Audit a triangle mesh with an explicit source-unit tolerance.
///
/// A triangle is degenerate when its doubled area is at most
/// `tolerance.linear()²`; the implementation compares squared values to avoid a
/// square root. Pass [`Tolerance::ZERO`] for exact-coordinate compatibility.
///
/// This compatibility entry point falls back to the historical map-backed
/// audit if the bounded vector reservation fails. Operations with an explicit
/// memory budget should use [`try_audit_mesh`] and preflight
/// [`audit_mesh_scratch_bytes`] instead.
pub fn audit_mesh<M: TriangleMeshView + ?Sized>(mesh: &M, tolerance: Tolerance) -> MeshHealth {
    match VecEdgeSink::try_new(mesh.triangle_count(), mesh.position_count()) {
        Ok(edges) => audit_with_edges(mesh, tolerance, edges),
        Err(_) => audit_with_edges(mesh, tolerance, MapEdgeSink::default()),
    }
}

/// Audit using one fallible, precomputable edge-record allocation.
///
/// Callers can refuse before allocation by comparing
/// [`audit_mesh_scratch_bytes`] with their memory budget.
pub fn try_audit_mesh<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    tolerance: Tolerance,
) -> Result<MeshHealth, MeshAuditError> {
    let edges = VecEdgeSink::try_new(mesh.triangle_count(), mesh.position_count())?;
    Ok(audit_with_edges(mesh, tolerance, edges))
}

fn audit_with_edges<M: TriangleMeshView + ?Sized, E: EdgeSink>(
    mesh: &M,
    tolerance: Tolerance,
    mut edges: E,
) -> MeshHealth {
    let positions = mesh.position_count();
    let triangles = mesh.triangle_count();
    let first_non_finite_position = (0..positions).find(|&index| !mesh.position(index).is_finite());
    let non_finite_positions = (0..positions)
        .filter(|&index| !mesh.position(index).is_finite())
        .count();
    let mut invalid_indices = 0;
    let mut first_invalid_index = None;
    let mut degenerate_triangles = 0;
    let mut usable_triangles = 0;
    let squared_double_area_limit = tolerance.linear().powi(4);

    for triangle_index in 0..triangles {
        let triangle = mesh.triangle(triangle_index);
        let converted = triangle.map(|source_index| usize::try_from(source_index).ok());
        let [Some(a_index), Some(b_index), Some(c_index)] = converted else {
            for source_index in triangle {
                if usize::try_from(source_index).is_err() {
                    invalid_indices += 1;
                    first_invalid_index.get_or_insert((triangle_index, source_index));
                }
            }
            continue;
        };
        let indices = [a_index, b_index, c_index];
        let mut valid = true;
        for (corner, &index) in indices.iter().enumerate() {
            if index >= positions {
                invalid_indices += 1;
                first_invalid_index.get_or_insert((triangle_index, triangle[corner]));
                valid = false;
            }
        }
        if !valid {
            continue;
        }

        let [a, b, c] = indices.map(|index| mesh.position(index));
        if !a.is_finite() || !b.is_finite() || !c.is_finite() {
            continue;
        }
        let squared_double_area = (b - a).cross(c - a).length_squared();
        if squared_double_area <= squared_double_area_limit {
            degenerate_triangles += 1;
            continue;
        }
        usable_triangles += 1;
        for (left, right) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            let direction = if left < right { 1 } else { -1 };
            edges.record(left.min(right), left.max(right), direction);
        }
    }

    let edge_summary = edges.summarize();
    MeshHealth {
        positions,
        triangles,
        usable_triangles,
        invalid_indices,
        non_finite_positions,
        degenerate_triangles,
        boundary_edges: edge_summary.boundary,
        non_manifold_edges: edge_summary.non_manifold,
        inconsistent_winding_edges: edge_summary.inconsistent_winding,
        first_invalid_index,
        first_non_finite_position,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TriMesh;
    use axiolid_core::Point3;

    /// The bound covers BOTH buffers the counting sort may hold at once.
    /// Charging for one would let a caller admit an audit that then
    /// refuses its own allocation.
    #[test]
    fn scratch_bound_charges_two_buffers_of_three_edge_records_per_triangle() {
        assert_eq!(
            audit_mesh_scratch_bytes(7),
            7usize
                .checked_mul(3)
                .and_then(|count| count.checked_mul(std::mem::size_of::<EdgeRecord>()))
                .and_then(|bytes| bytes.checked_mul(2))
        );
        let sink = VecEdgeSink::try_new(7, 16).expect("small bounded audit allocation");
        assert!(sink.edges.capacity() >= 21);
    }

    /// The counting sort must order records exactly as the comparison
    /// sort did. Grouping alone would be enough for `summarize`, but
    /// proving full order equality is stronger and catches a pass
    /// ordering mistake that grouping would hide.
    #[test]
    fn counting_sort_matches_comparison_sort() {
        let buckets = 64usize;
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut edges: Vec<EdgeRecord> = (0..4096)
            .map(|_| {
                let a = (next() as usize % buckets) as u64;
                let b = (next() as usize % buckets) as u64;
                EdgeRecord {
                    low: a.min(b),
                    high: a.max(b),
                    direction: if next() % 2 == 0 { 1 } else { -1 },
                }
            })
            .collect();
        let mut expected = edges.clone();
        expected.sort_unstable_by_key(|e| (e.low, e.high));

        let mut scratch = vec![EdgeRecord::EMPTY; edges.len()];
        counting_sort_edges(&mut edges, &mut scratch, buckets);

        let keys: Vec<_> = edges.iter().map(|e| (e.low, e.high)).collect();
        let want: Vec<_> = expected.iter().map(|e| (e.low, e.high)).collect();
        assert_eq!(
            keys, want,
            "counting sort must reproduce the comparison order"
        );
    }

    /// Both sinks must agree on a real mesh. The map sink is the
    /// allocation-failure fallback, so a divergence here would mean the
    /// audit silently reports different health under memory pressure.
    #[test]
    fn both_sinks_agree_on_a_defective_mesh() {
        // Two triangles sharing an edge, plus a third fin on that same
        // edge: boundary, non-manifold and winding counts all exercised.
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, -1.0, 0.0),
        ];
        let indices = vec![0, 1, 2, 0, 1, 3, 0, 1, 4];
        let mesh = TriMesh::new(positions, indices);

        let tolerance = Tolerance::MILLIMETRE;
        let mut fast = VecEdgeSink::try_new(mesh.triangle_count(), mesh.position_count())
            .expect("fixture allocation");
        assert!(
            fast.buckets > 0,
            "counting sort must be active for this fixture"
        );
        let via_counting = audit_with_edges(&mesh, tolerance, fast);
        let via_map = audit_with_edges(&mesh, tolerance, MapEdgeSink::default());
        assert_eq!(via_counting, via_map, "sinks must report identical health");
    }
}
