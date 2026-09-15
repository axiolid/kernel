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
///
/// # Layout
///
/// Edges are held in compressed-sparse-row form: `keys` ascending, `starts`
/// giving each key's span, and `uses` one flat run per edge. A
/// `BTreeMap<EdgeKey, Vec<EdgeUse>>` costs a node allocation per distinct
/// edge plus a `Vec` allocation per edge that is ever used twice -- roughly
/// 245,000 allocations on an 82k-triangle sphere, against three here.
///
/// The structure is immutable after the build, which is what makes this
/// affordable: CSR cannot accept a late insertion without reflowing, and
/// nothing in the API offers one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeAdjacency {
    /// Distinct edge keys, ascending.
    keys: Vec<EdgeKey>,
    /// `starts[i]..starts[i + 1]` is edge `i`'s span in `uses`.
    ///
    /// Length is `keys.len() + 1`, so the last span needs no special case.
    starts: Vec<u32>,
    /// Uses grouped by edge, ascending by triangle within each edge.
    uses: Vec<EdgeUse>,
    vertex_count: usize,
    degenerate_triangles: usize,
    /// Triangles that contributed adjacency, recorded during the build.
    ///
    /// Counting these by walking the edge map means inserting every
    /// triangle index into a set -- three visits per triangle and an
    /// allocation per node, to recover a number the build already knew.
    face_count: usize,
}

impl EdgeAdjacency {
    /// Derive adjacency from a mesh.
    ///
    /// Triangles with a repeated corner are skipped and counted. Corners are
    /// taken as given: an out-of-range index is not adjacency data, and
    /// validating it belongs to `audit_mesh`, not here.
    #[must_use]
    pub fn build(mesh: &TriMesh) -> Self {
        let triangles = mesh.indices.len() / 3;
        // Every usable triangle contributes exactly three edge records, so
        // the whole build fits in one allocation sized up front.
        let mut records: Vec<(EdgeKey, EdgeUse)> = Vec::with_capacity(triangles * 3);
        let mut degenerate_triangles = 0;
        let mut face_count = 0;

        for (triangle, chunk) in mesh.indices.chunks_exact(3).enumerate() {
            let (a, b, c) = (chunk[0], chunk[1], chunk[2]);
            if a == b || b == c || c == a {
                degenerate_triangles += 1;
                continue;
            }
            // Past the guard this triangle contributes all three of its
            // edges, so it is exactly one face.
            face_count += 1;
            for (from, to) in [(a, b), (b, c), (c, a)] {
                records.push((
                    EdgeKey::new(from, to),
                    EdgeUse {
                        triangle,
                        forward: from <= to,
                    },
                ));
            }
        }

        // Sort by key, then by triangle within a key. `sort_unstable_by_key`
        // would leave equal keys in an unspecified order, and the uses of one
        // edge are documented as ascending by triangle -- callers pair
        // `uses[0]`/`uses[1]` to judge winding, so the order is observable.
        records.sort_unstable_by(|left, right| {
            left.0
                .cmp(&right.0)
                .then_with(|| left.1.triangle.cmp(&right.1.triangle))
        });

        let mut keys: Vec<EdgeKey> = Vec::new();
        let mut starts: Vec<u32> = Vec::new();
        let mut uses: Vec<EdgeUse> = Vec::with_capacity(records.len());
        for (key, use_) in records {
            if keys.last() != Some(&key) {
                keys.push(key);
                // This edge's run begins where the previous one ended.
                starts.push(uses.len() as u32);
            }
            uses.push(use_);
        }
        // Sentinel closing the last run, so `starts[i]..starts[i + 1]` is valid
        // for every edge and `starts.len() == keys.len() + 1` even when empty.
        starts.push(uses.len() as u32);

        Self {
            keys,
            starts,
            uses,
            vertex_count: mesh.positions.len(),
            degenerate_triangles,
            face_count,
        }
    }

    /// Distinct undirected edges.
    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.keys.len()
    }

    /// Triangles skipped for having a repeated corner.
    #[must_use]
    pub const fn degenerate_triangles(&self) -> usize {
        self.degenerate_triangles
    }

    /// Uses of the edge at `index`, which must be in range.
    fn span(&self, index: usize) -> &[EdgeUse] {
        let from = self.starts[index] as usize;
        let to = self.starts[index + 1] as usize;
        &self.uses[from..to]
    }

    /// Every edge with its uses, ordered by [`EdgeKey`].
    pub fn edges(&self) -> impl Iterator<Item = (EdgeKey, &[EdgeUse])> {
        self.keys
            .iter()
            .enumerate()
            .map(|(index, key)| (*key, self.span(index)))
    }

    /// Triangles incident to one edge, or an empty slice if it is absent.
    #[must_use]
    pub fn uses(&self, edge: EdgeKey) -> &[EdgeUse] {
        // Keys are ascending, so this is the CSR equivalent of the map
        // lookup it replaces.
        self.keys
            .binary_search(&edge)
            .map_or(&[], |index| self.span(index))
    }

    /// Edges used by exactly one triangle: the mesh boundary.
    ///
    /// In a closed shell this is empty, which is what makes it a usable
    /// definition of "the hole" after a clip or a cut.
    pub fn boundary_edges(&self) -> impl Iterator<Item = EdgeKey> + '_ {
        self.edges()
            .filter(|(_, uses)| uses.len() == 1)
            .map(|(key, _)| key)
    }

    /// Edges used by three or more triangles.
    pub fn non_manifold_edges(&self) -> impl Iterator<Item = EdgeKey> + '_ {
        self.edges()
            .filter(|(_, uses)| uses.len() > 2)
            .map(|(key, _)| key)
    }

    /// Edges whose two triangles traverse them the same way.
    ///
    /// Consistently wound neighbours traverse a shared edge in opposite
    /// directions, so agreement means one of the pair is flipped. Reported
    /// only for two-triangle edges: with three or more the pairing is
    /// ambiguous, and that is already a non-manifold defect.
    pub fn inconsistent_edges(&self) -> impl Iterator<Item = EdgeKey> + '_ {
        self.edges()
            .filter(|(_, uses)| uses.len() == 2 && uses[0].forward == uses[1].forward)
            .map(|(key, _)| key)
    }

    /// Whether every edge has exactly two consistently wound triangles.
    #[must_use]
    pub fn is_closed_two_manifold(&self) -> bool {
        self.edges().all(|(_, uses)| uses.len() == 2) && self.inconsistent_edges().next().is_none()
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
        for key in &self.keys {
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
        for (_, uses) in self.edges() {
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
        for key in &self.keys {
            for corner in [key.lower(), key.upper()] {
                if let Some(slot) = used.get_mut(corner as usize) {
                    *slot = true;
                }
            }
        }
        let vertices = used.iter().filter(|seen| **seen).count() as i64;
        let edges = self.keys.len() as i64;
        let faces = self.face_count() as i64;
        vertices - edges + faces
    }

    /// Usable triangles, i.e. those that contributed adjacency.
    #[must_use]
    pub const fn face_count(&self) -> usize {
        self.face_count
    }
}
