//! Splitting a mesh into connected components, and recombining them (#85).
//!
//! Two providers already walked this graph to count components and returned
//! only the count, discarding the partition they had just built. This owns
//! the partition; the counters become callers of it.
//!
//! Connectivity is by shared vertex index, not by geometric proximity: two
//! triangles that touch in space but reference distinct indices belong to
//! distinct components. Welding is `axiolid-heal`'s job and stays there.

use crate::TriMesh;
use std::collections::BTreeMap;

/// Disjoint-set over vertex indices, unioned across each triangle's corners.
struct Partition {
    parent: Vec<usize>,
}

impl Partition {
    fn new(count: usize) -> Self {
        Self {
            parent: (0..count).collect(),
        }
    }

    fn find(&mut self, mut node: usize) -> usize {
        while self.parent[node] != node {
            self.parent[node] = self.parent[self.parent[node]];
            node = self.parent[node];
        }
        node
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
}

/// Build the vertex partition for a mesh.
fn partition_of(mesh: &TriMesh) -> Partition {
    let mut partition = Partition::new(mesh.positions.len());
    for triangle in mesh.indices.as_chunks::<3>().0 {
        let first = triangle[0] as usize;
        for corner in &triangle[1..] {
            partition.union(first, *corner as usize);
        }
    }
    partition
}

/// Number of connected components, counting only referenced vertices.
///
/// An unreferenced position is unused data, not a component of its own.
pub fn component_count(mesh: &TriMesh) -> usize {
    if mesh.indices.is_empty() {
        return 0;
    }
    let mut partition = partition_of(mesh);
    let mut roots = std::collections::BTreeSet::new();
    for index in &mesh.indices {
        let root = partition.find(*index as usize);
        roots.insert(root);
    }
    roots.len()
}

/// Split a mesh into its connected components.
///
/// Components are ordered by the smallest triangle index they contain, so
/// the same input always produces the same order -- iteration order of a
/// hash map would not, and #85 requires determinism.
///
/// A single-body mesh returns one component whose triangles are in input
/// order, so decomposing a connected mesh is not a reindexing hazard.
pub fn decompose(mesh: &TriMesh) -> Vec<TriMesh> {
    if mesh.indices.is_empty() {
        return Vec::new();
    }
    let mut partition = partition_of(mesh);

    // Group triangles by the root of their first corner. BTreeMap keyed by
    // first appearance keeps the order deterministic and input-shaped.
    let mut order: BTreeMap<usize, usize> = BTreeMap::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (triangle_index, triangle) in mesh.indices.as_chunks::<3>().0.iter().enumerate() {
        let root = partition.find(triangle[0] as usize);
        let slot = *order.entry(root).or_insert_with(|| {
            groups.push(Vec::new());
            groups.len() - 1
        });
        groups[slot].push(triangle_index);
    }

    // A connected mesh comes back exactly as it went in. Reindexing a
    // single-body mesh would make "decompose defensively" cost a rebuild and
    // invalidate any index the caller already holds.
    if groups.len() == 1 && mesh.positions.len() == referenced(mesh) {
        return vec![mesh.clone()];
    }

    groups
        .into_iter()
        .map(|triangles| extract(mesh, &triangles))
        .collect()
}

/// Count of distinct vertices actually referenced by a triangle.
fn referenced(mesh: &TriMesh) -> usize {
    mesh.indices
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len()
}

/// Build one component mesh from a triangle selection.
///
/// Vertices are emitted in first-use order and reindexed, so a component
/// carries only the positions it references.
fn extract(mesh: &TriMesh, triangles: &[usize]) -> TriMesh {
    let mut remap: BTreeMap<u32, u32> = BTreeMap::new();
    let mut positions = Vec::new();
    let mut indices = Vec::with_capacity(triangles.len() * 3);
    for &triangle in triangles {
        for corner in 0..3 {
            let old = mesh.indices[triangle * 3 + corner];
            let new = *remap.entry(old).or_insert_with(|| {
                positions.push(mesh.positions[old as usize]);
                (positions.len() - 1) as u32
            });
            indices.push(new);
        }
    }
    TriMesh::new(positions, indices)
}

/// Combine meshes into one, rebasing each mesh's indices.
///
/// The inverse of [`decompose`] up to vertex ordering: positions are
/// concatenated in argument order and never merged, so composing meshes
/// that share a coordinate leaves them as separate components. Merging
/// coincident vertices is `axiolid-heal`'s weld, deliberately not done here.
pub fn compose(meshes: &[TriMesh]) -> TriMesh {
    let mut positions = Vec::new();
    let mut indices = Vec::new();
    for mesh in meshes {
        let base = positions.len() as u32;
        positions.extend_from_slice(&mesh.positions);
        indices.extend(mesh.indices.iter().map(|index| index + base));
    }
    TriMesh::new(positions, indices)
}
