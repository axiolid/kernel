//! Edge adjacency derived once, so algorithms stop rebuilding it.
//!
//! # Why this exists
//!
//! Healing, genus, smoothing, decimation, and decomposition each needed the
//! same fact: which triangles meet along each edge. Every one of them built
//! its own `BTreeMap<(u32, u32), _>` inline, with its own key convention and
//! its own degenerate-triangle handling. Five implementations of one idea is
//! five places for the invariant to drift.
//!
//! This derives it once. The algorithms ask questions instead of rebuilding
//! the answer.
//!
//! # Not a half-edge structure
//!
//! A true half-edge (or winged-edge) representation stores per-half-edge
//! `next`/`twin`/`face` links and supports *mutation* through them. That is
//! the right structure for algorithms that rewrite connectivity in place.
//!
//! This is deliberately less: an immutable, derived index answering
//! adjacency queries over a `TriMesh` that stays the owner of the data.
//! Every current consumer reads adjacency and writes a *new* mesh, so
//! nothing needs mutable topology, and a mutable structure would add an
//! invariant to maintain for no consumer.
//!
//! The name says what it is. If in-place connectivity editing ever appears,
//! that is a separate type, not a field bolted onto this one.
//!
//! # Degenerate triangles
//!
//! A triangle with a repeated corner (`[4, 4, 7]`) has an edge from a vertex
//! to itself. Counting it as adjacency makes a sound mesh look non-manifold.
//! Such triangles are excluded and counted in
//! [`EdgeAdjacency::degenerate_triangles`], matching what `audit_mesh`
//! already does, so the two never disagree about what is usable.

use crate::TriMesh;
use std::collections::BTreeMap;

/// An undirected edge, canonically ordered so both sides collide.
///
/// Construction is the only place ordering is decided, which is what stops
/// two call sites from disagreeing about whether `(3, 1)` and `(1, 3)` are
/// the same edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeKey {
    lower: u32,
    upper: u32,
}

impl EdgeKey {
    /// Canonical key for an edge between two corners.
    ///
    /// Order-insensitive: `new(a, b) == new(b, a)`.
    #[must_use]
    pub fn new(a: u32, b: u32) -> Self {
        if a <= b {
            Self { lower: a, upper: b }
        } else {
            Self { lower: b, upper: a }
        }
    }

    /// Lower-numbered endpoint.
    #[must_use]
    pub const fn lower(self) -> u32 {
        self.lower
    }

    /// Higher-numbered endpoint.
    #[must_use]
    pub const fn upper(self) -> u32 {
        self.upper
    }

    /// Both endpoints, lower first.
    #[must_use]
    pub const fn endpoints(self) -> (u32, u32) {
        (self.lower, self.upper)
    }
}

/// One triangle's use of an edge, and the direction it traversed.
///
/// Direction is what reveals winding: two triangles sharing an edge are
/// consistently wound when they traverse it in *opposite* directions. A pair
/// agreeing on direction is the classic flipped-face defect.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EdgeUse {
    /// Triangle index into the mesh.
    pub triangle: usize,
    /// Whether this triangle traversed the edge low-to-high.
    pub forward: bool,
}

/// Edge-to-triangle adjacency over a triangle mesh.
///
/// Built once with [`EdgeAdjacency::build`], then queried. Iteration order is
/// by [`EdgeKey`], so any diagnosis derived from it is reproducible -- a
/// report that reorders between runs is useless as an audit record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeAdjacency {
    edges: BTreeMap<EdgeKey, Vec<EdgeUse>>,
    vertex_count: usize,
    degenerate_triangles: usize,
}

impl EdgeAdjacency {
    /// Derive adjacency from a mesh.
    ///
    /// Triangles with a repeated corner are skipped and counted. Corners are
    /// taken as given: an out-of-range index is not adjacency data, and
    /// validating it belongs to `audit_mesh`, not here.
    #[must_use]
    pub fn build(mesh: &TriMesh) -> Self {
        let mut edges: BTreeMap<EdgeKey, Vec<EdgeUse>> = BTreeMap::new();
        let mut degenerate_triangles = 0;

        for (triangle, chunk) in mesh.indices.chunks_exact(3).enumerate() {
            let (a, b, c) = (chunk[0], chunk[1], chunk[2]);
            if a == b || b == c || c == a {
                degenerate_triangles += 1;
                continue;
            }
            for (from, to) in [(a, b), (b, c), (c, a)] {
                let key = EdgeKey::new(from, to);
                edges.entry(key).or_default().push(EdgeUse {
                    triangle,
                    forward: from <= to,
                });
            }
        }

        Self {
            edges,
            vertex_count: mesh.positions.len(),
            degenerate_triangles,
        }
    }

    /// Distinct undirected edges.
    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Triangles skipped for having a repeated corner.
    #[must_use]
    pub const fn degenerate_triangles(&self) -> usize {
        self.degenerate_triangles
    }

    /// Every edge with its uses, ordered by [`EdgeKey`].
    pub fn edges(&self) -> impl Iterator<Item = (EdgeKey, &[EdgeUse])> {
        self.edges.iter().map(|(key, uses)| (*key, uses.as_slice()))
    }

    /// Triangles incident to one edge, or an empty slice if it is absent.
    #[must_use]
    pub fn uses(&self, edge: EdgeKey) -> &[EdgeUse] {
        self.edges.get(&edge).map_or(&[], Vec::as_slice)
    }

    /// Edges used by exactly one triangle: the mesh boundary.
    ///
    /// In a closed shell this is empty, which is what makes it a usable
    /// definition of "the hole" after a clip or a cut.
    pub fn boundary_edges(&self) -> impl Iterator<Item = EdgeKey> + '_ {
        self.edges
            .iter()
            .filter(|(_, uses)| uses.len() == 1)
            .map(|(key, _)| *key)
    }

    /// Edges used by three or more triangles.
    pub fn non_manifold_edges(&self) -> impl Iterator<Item = EdgeKey> + '_ {
        self.edges
            .iter()
            .filter(|(_, uses)| uses.len() > 2)
            .map(|(key, _)| *key)
    }

    /// Edges whose two triangles traverse them the same way.
    ///
    /// Consistently wound neighbours traverse a shared edge in opposite
    /// directions, so agreement means one of the pair is flipped. Reported
    /// only for two-triangle edges: with three or more the pairing is
    /// ambiguous, and that is already a non-manifold defect.
    pub fn inconsistent_edges(&self) -> impl Iterator<Item = EdgeKey> + '_ {
        self.edges
            .iter()
            .filter(|(_, uses)| uses.len() == 2 && uses[0].forward == uses[1].forward)
            .map(|(key, _)| *key)
    }

    /// Whether every edge has exactly two consistently wound triangles.
    #[must_use]
    pub fn is_closed_two_manifold(&self) -> bool {
        self.edges.values().all(|uses| uses.len() == 2)
            && self.inconsistent_edges().next().is_none()
    }

    /// Vertices touched by a boundary edge.
    #[must_use]
    pub fn boundary_vertices(&self) -> Vec<u32> {
        let mut seen = vec![false; self.vertex_count];
        let mut out = Vec::new();
        for key in self.boundary_edges() {
            for corner in [key.lower(), key.upper()] {
                if let Some(slot) = seen.get_mut(corner as usize) {
                    if !*slot {
                        *slot = true;
                        out.push(corner);
                    }
                }
            }
        }
        out.sort_unstable();
        out
    }

    /// Vertex-to-vertex neighbours, indexed by vertex.
    ///
    /// Entry `v` lists the vertices sharing an edge with `v`, ascending.
    /// Vertices used by no usable triangle get an empty list rather than
    /// being omitted, so the result can be indexed directly.
    #[must_use]
    pub fn vertex_neighbours(&self) -> Vec<Vec<u32>> {
        let mut out = vec![Vec::new(); self.vertex_count];
        for key in self.edges.keys() {
            let (a, b) = key.endpoints();
            if let Some(list) = out.get_mut(a as usize) {
                list.push(b);
            }
            if let Some(list) = out.get_mut(b as usize) {
                list.push(a);
            }
        }
        for list in &mut out {
            list.sort_unstable();
            list.dedup();
        }
        out
    }

    /// Triangles sharing an edge with `triangle`, ascending.
    #[must_use]
    pub fn triangle_neighbours(&self, triangle: usize) -> Vec<usize> {
        let mut out = Vec::new();
        for uses in self.edges.values() {
            if uses.iter().any(|use_| use_.triangle == triangle) {
                out.extend(
                    uses.iter()
                        .map(|use_| use_.triangle)
                        .filter(|other| *other != triangle),
                );
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Euler characteristic `V - E + F` over usable triangles.
    ///
    /// Vertices are counted as those actually used by an edge, not the
    /// length of the position array: an unreferenced position is not part of
    /// the surface, and counting it shifts the characteristic silently.
    #[must_use]
    pub fn euler_characteristic(&self) -> i64 {
        let mut used = vec![false; self.vertex_count];
        for key in self.edges.keys() {
            for corner in [key.lower(), key.upper()] {
                if let Some(slot) = used.get_mut(corner as usize) {
                    *slot = true;
                }
            }
        }
        let vertices = used.iter().filter(|seen| **seen).count() as i64;
        let edges = self.edges.len() as i64;
        let faces = self.face_count() as i64;
        vertices - edges + faces
    }

    /// Usable triangles, i.e. those that contributed adjacency.
    #[must_use]
    pub fn face_count(&self) -> usize {
        let mut seen = std::collections::BTreeSet::new();
        for uses in self.edges.values() {
            for use_ in uses {
                seen.insert(use_.triangle);
            }
        }
        seen.len()
    }
}
