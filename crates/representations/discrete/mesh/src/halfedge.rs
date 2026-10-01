//! Halfedge surface mesh: O(1) adjacency and in-place connectivity edits.
//!
//! # Why this exists
//!
//! [`EdgeAdjacency`](crate::EdgeAdjacency) answers "which triangles meet
//! along this edge" over a [`TriMesh`] that stays the owner of the data, and
//! every consumer that reads it writes a *new* mesh. Remeshing, hole filling,
//! subdivision, geodesics and parameterisation instead rewrite connectivity
//! in place, many times, and walk the one-ring of a vertex in its cyclic
//! order. That needs a persistent structure whose links survive each edit:
//! this is it.
//!
//! # Layout
//!
//! The layout follows CGAL's `Surface_mesh`: index-based arrays rather than
//! pointers. Each edge `e` owns the two halfedges `2e` and `2e + 1`, so
//! [`HalfedgeMesh::opposite`] is `h ^ 1` and stores nothing. Each halfedge
//! stores its `next`, `prev`, `target` vertex and incident face; a vertex
//! stores one outgoing halfedge and a face one halfedge of its loop. Every
//! navigation query is one array read.
//!
//! A halfedge with no face is a *boundary* halfedge. Boundary halfedges are
//! linked by `next`/`prev` into loops, one per hole, so a hole is walked the
//! same way as a face. An isolated vertex has no halfedge.
//!
//! # Invariants
//!
//! [`HalfedgeMesh::validate`] checks every one of these; every constructor
//! and edit keeps them.
//!
//! - `next` and `prev` are inverse permutations of the live halfedges, and
//!   `face(next(h)) == face(h)`.
//! - `source(h) == target(opposite(h))` and an edge never joins a vertex to
//!   itself.
//! - Every face loop has at least three halfedges and no repeated vertex.
//! - No edge has a boundary on both sides.
//! - The mesh is *simple*: two vertices are joined by at most one edge.
//! - Every vertex is manifold: rotating around it reaches every outgoing
//!   halfedge, at most one of them is a boundary halfedge, and when one is,
//!   it is the vertex's stored halfedge. So [`HalfedgeMesh::is_boundary_vertex`]
//!   is O(1) and a vertex circulation starts at the boundary.
//!
//! Non-manifold input is refused by name ([`HalfedgeBuildError`]) rather
//! than represented: an edge shared by three faces, a vertex pinched between
//! two fans, and two faces disagreeing on winding each have their own
//! variant.
//!
//! # Identity
//!
//! A mesh built from a [`TriMesh`] or face list keeps the input's numbering:
//! vertex `i` is position `i` and face `j` is input face `j`. Edits that
//! remove elements mark them removed instead of renumbering, so ids held by
//! a caller stay valid across edits; [`HalfedgeMesh::compact`] renumbers
//! densely and returns the mapping.

mod build;
mod check;
mod edit;

use core::fmt;

use axiolid_core::Point3;

use crate::TriMesh;

pub use build::HalfedgeBuildError;
pub use check::HalfedgeInvariantError;
pub use edit::HalfedgeEditError;

/// Sentinel for "no element" in the link arrays.
const NONE: u32 = u32::MAX;

macro_rules! element_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            /// Id with the given index.
            #[must_use]
            pub const fn new(index: u32) -> Self {
                Self(index)
            }

            /// Index into the mesh's arrays.
            #[must_use]
            pub const fn index(self) -> usize {
                self.0 as usize
            }

            /// Raw index.
            #[must_use]
            pub const fn get(self) -> u32 {
                self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}", self.0)
            }
        }
    };
}

element_id!(
    /// A vertex of a [`HalfedgeMesh`].
    VertexId
);
element_id!(
    /// A directed halfedge of a [`HalfedgeMesh`].
    HalfedgeId
);
element_id!(
    /// An undirected edge of a [`HalfedgeMesh`]; it owns halfedges `2e` and `2e + 1`.
    EdgeId
);
element_id!(
    /// A face of a [`HalfedgeMesh`].
    FaceId
);

/// One halfedge's links.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Link {
    next: u32,
    prev: u32,
    target: u32,
    /// `NONE` for a boundary halfedge.
    face: u32,
}

const DEAD_LINK: Link = Link {
    next: NONE,
    prev: NONE,
    target: NONE,
    face: NONE,
};

/// Halfedge (doubly-connected edge list) surface mesh of a 2-manifold, with
/// or without boundary, with polygonal faces.
///
/// See the [module documentation](self) for the layout and the invariants.
/// Navigation queries take ids of live elements; passing an id of a removed
/// element returns stale data rather than panicking, and an id past the
/// arrays panics like slice indexing. Edits check liveness and refuse by
/// name.
#[derive(Debug, Clone, PartialEq)]
pub struct HalfedgeMesh {
    positions: Vec<Point3>,
    /// One outgoing halfedge per vertex; `NONE` for an isolated vertex.
    vertex_out: Vec<u32>,
    vertex_live: Vec<bool>,
    links: Vec<Link>,
    edge_live: Vec<bool>,
    face_halfedge: Vec<u32>,
    face_live: Vec<bool>,
    live_vertices: usize,
    live_edges: usize,
    live_faces: usize,
}

/// Old-to-new id mapping returned by [`HalfedgeMesh::compact`].
///
/// `None` marks an element that was removed. Halfedges follow their edge:
/// halfedge `2e + s` becomes `2 * edges[e] + s`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct HalfedgeRemap {
    /// New id of each old vertex.
    pub vertices: Vec<Option<VertexId>>,
    /// New id of each old edge.
    pub edges: Vec<Option<EdgeId>>,
    /// New id of each old face.
    pub faces: Vec<Option<FaceId>>,
}

impl HalfedgeMesh {
    /// Build from an indexed triangle mesh.
    ///
    /// Vertex `i` is `mesh.positions[i]` and face `j` is triangle `j`.
    /// Normals and attribute channels are not carried: the halfedge mesh is
    /// connectivity plus positions.
    ///
    /// # Errors
    ///
    /// A ragged index buffer, an out-of-range or repeated corner, and every
    /// non-manifold or inconsistently wound configuration, each by name.
    pub fn from_tri_mesh(mesh: &TriMesh) -> Result<Self, HalfedgeBuildError> {
        if mesh.indices.len() % 3 != 0 {
            return Err(HalfedgeBuildError::IncompleteTriangle {
                index_count: mesh.indices.len(),
            });
        }
        let triangles: Vec<[u32; 3]> = mesh.triangles().collect();
        Self::from_faces(mesh.positions.clone(), &triangles)
    }

    /// Build from positions and polygonal faces, each a loop of position
    /// indices wound consistently with its neighbours.
    ///
    /// # Errors
    ///
    /// As [`HalfedgeMesh::from_tri_mesh`], plus faces with fewer than three
    /// corners.
    pub fn from_faces<F: AsRef<[u32]>>(
        positions: Vec<Point3>,
        faces: &[F],
    ) -> Result<Self, HalfedgeBuildError> {
        build::build(positions, faces)
    }

    /// Indexed triangle mesh of the live elements.
    ///
    /// Live vertices keep their relative order, so a mesh that has had no
    /// vertex removed keeps its numbering; each triangle starts at its face's
    /// stored halfedge, so building from a [`TriMesh`] and converting back
    /// reproduces its index buffer exactly.
    ///
    /// # Errors
    ///
    /// [`HalfedgeEditError::NotATriangle`] for the first face that is not a
    /// triangle; split it with [`HalfedgeMesh::split_face`] or
    /// [`HalfedgeMesh::split_face_diagonal`] first.
    pub fn to_tri_mesh(&self) -> Result<TriMesh, HalfedgeEditError> {
        let mut remap = vec![NONE; self.positions.len()];
        let mut positions = Vec::with_capacity(self.live_vertices);
        for v in self.vertices() {
            remap[v.index()] = positions.len() as u32;
            positions.push(self.positions[v.index()]);
        }
        let mut indices = Vec::with_capacity(self.live_faces * 3);
        for f in self.faces() {
            let degree = self.face_degree(f);
            if degree != 3 {
                return Err(HalfedgeEditError::NotATriangle { face: f, degree });
            }
            let h = self.face_halfedge(f);
            for corner in [self.source(h), self.target(h), self.target(self.next(h))] {
                indices.push(remap[corner.index()]);
            }
        }
        Ok(TriMesh::new(positions, indices))
    }

    // ---- counts -----------------------------------------------------------

    /// Live vertices, including isolated ones.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.live_vertices
    }

    /// Live edges.
    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.live_edges
    }

    /// Live halfedges: twice the edge count.
    #[must_use]
    pub fn halfedge_count(&self) -> usize {
        self.live_edges * 2
    }

    /// Live faces.
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.live_faces
    }

    /// `V - E + F` over the live elements (isolated vertices count).
    #[must_use]
    pub fn euler_characteristic(&self) -> i64 {
        self.live_vertices as i64 - self.live_edges as i64 + self.live_faces as i64
    }

    /// Whether the vertex id names a live vertex.
    #[must_use]
    pub fn contains_vertex(&self, v: VertexId) -> bool {
        self.vertex_live.get(v.index()).copied().unwrap_or(false)
    }

    /// Whether the edge id names a live edge.
    #[must_use]
    pub fn contains_edge(&self, e: EdgeId) -> bool {
        self.edge_live.get(e.index()).copied().unwrap_or(false)
    }

    /// Whether the halfedge id names a halfedge of a live edge.
    #[must_use]
    pub fn contains_halfedge(&self, h: HalfedgeId) -> bool {
        self.contains_edge(EdgeId(h.0 / 2))
    }

    /// Whether the face id names a live face.
    #[must_use]
    pub fn contains_face(&self, f: FaceId) -> bool {
        self.face_live.get(f.index()).copied().unwrap_or(false)
    }

    // ---- element iteration ------------------------------------------------

    /// Live vertices in id order.
    pub fn vertices(&self) -> impl Iterator<Item = VertexId> + '_ {
        live_ids(&self.vertex_live).map(VertexId)
    }

    /// Live edges in id order.
    pub fn edges(&self) -> impl Iterator<Item = EdgeId> + '_ {
        live_ids(&self.edge_live).map(EdgeId)
    }

    /// Live halfedges in id order.
    pub fn halfedges(&self) -> impl Iterator<Item = HalfedgeId> + '_ {
        self.edges()
            .flat_map(|e| [HalfedgeId(e.0 * 2), HalfedgeId(e.0 * 2 + 1)])
    }

    /// Live faces in id order.
    pub fn faces(&self) -> impl Iterator<Item = FaceId> + '_ {
        live_ids(&self.face_live).map(FaceId)
    }

    /// Live boundary halfedges in id order.
    pub fn boundary_halfedges(&self) -> impl Iterator<Item = HalfedgeId> + '_ {
        self.halfedges().filter(|&h| self.is_boundary_halfedge(h))
    }

    /// One boundary halfedge per hole, each the lowest id on its loop.
    #[must_use]
    pub fn boundary_loops(&self) -> Vec<HalfedgeId> {
        let mut seen = vec![false; self.links.len()];
        let mut loops = Vec::new();
        for h in self.boundary_halfedges() {
            if seen[h.index()] {
                continue;
            }
            loops.push(h);
            for g in self.loop_halfedges(h) {
                seen[g.index()] = true;
            }
        }
        loops
    }

    // ---- O(1) navigation --------------------------------------------------

    /// Position of a vertex.
    #[must_use]
    pub fn position(&self, v: VertexId) -> Point3 {
        self.positions[v.index()]
    }

    /// Move a vertex. Connectivity is untouched.
    pub fn set_position(&mut self, v: VertexId, position: Point3) {
        self.positions[v.index()] = position;
    }

    /// The position array, indexed by vertex id (removed vertices included).
    #[must_use]
    pub fn positions(&self) -> &[Point3] {
        &self.positions
    }

    /// Next halfedge around the same face or hole.
    #[must_use]
    pub fn next(&self, h: HalfedgeId) -> HalfedgeId {
        HalfedgeId(self.links[h.index()].next)
    }

    /// Previous halfedge around the same face or hole.
    #[must_use]
    pub fn prev(&self, h: HalfedgeId) -> HalfedgeId {
        HalfedgeId(self.links[h.index()].prev)
    }

    /// The other halfedge of the same edge.
    #[must_use]
    pub const fn opposite(&self, h: HalfedgeId) -> HalfedgeId {
        HalfedgeId(h.0 ^ 1)
    }

    /// Vertex the halfedge points to.
    #[must_use]
    pub fn target(&self, h: HalfedgeId) -> VertexId {
        VertexId(self.links[h.index()].target)
    }

    /// Vertex the halfedge leaves.
    #[must_use]
    pub fn source(&self, h: HalfedgeId) -> VertexId {
        VertexId(self.links[(h.0 ^ 1) as usize].target)
    }

    /// Face to the left of the halfedge; `None` on a boundary halfedge.
    #[must_use]
    pub fn face(&self, h: HalfedgeId) -> Option<FaceId> {
        let face = self.links[h.index()].face;
        (face != NONE).then_some(FaceId(face))
    }

    /// Edge owning the halfedge.
    #[must_use]
    pub const fn edge(&self, h: HalfedgeId) -> EdgeId {
        EdgeId(h.0 / 2)
    }

    /// One of the edge's two halfedges: `side` 0 or 1.
    ///
    /// # Panics
    ///
    /// When `side` is greater than 1.
    #[must_use]
    pub fn edge_halfedge(&self, e: EdgeId, side: u32) -> HalfedgeId {
        assert!(side < 2, "an edge has two halfedges");
        HalfedgeId(e.0 * 2 + side)
    }

    /// Both endpoints of an edge, in the direction of its halfedge `2e`.
    #[must_use]
    pub fn edge_vertices(&self, e: EdgeId) -> [VertexId; 2] {
        let h = self.edge_halfedge(e, 0);
        [self.source(h), self.target(h)]
    }

    /// Stored outgoing halfedge of a vertex; `None` when isolated.
    ///
    /// For a boundary vertex it is the boundary halfedge leaving it.
    #[must_use]
    pub fn vertex_halfedge(&self, v: VertexId) -> Option<HalfedgeId> {
        let h = self.vertex_out[v.index()];
        (h != NONE).then_some(HalfedgeId(h))
    }

    /// Stored halfedge of a face's loop.
    #[must_use]
    pub fn face_halfedge(&self, f: FaceId) -> HalfedgeId {
        HalfedgeId(self.face_halfedge[f.index()])
    }

    /// Whether the halfedge has no face.
    #[must_use]
    pub fn is_boundary_halfedge(&self, h: HalfedgeId) -> bool {
        self.links[h.index()].face == NONE
    }

    /// Whether either side of the edge is a boundary.
    #[must_use]
    pub fn is_boundary_edge(&self, e: EdgeId) -> bool {
        self.is_boundary_halfedge(HalfedgeId(e.0 * 2))
            || self.is_boundary_halfedge(HalfedgeId(e.0 * 2 + 1))
    }

    /// Whether the vertex lies on a hole. Isolated vertices are not boundary.
    #[must_use]
    pub fn is_boundary_vertex(&self, v: VertexId) -> bool {
        self.vertex_halfedge(v)
            .is_some_and(|h| self.is_boundary_halfedge(h))
    }

    /// Whether the vertex has no incident edge.
    #[must_use]
    pub fn is_isolated(&self, v: VertexId) -> bool {
        self.vertex_out[v.index()] == NONE
    }

    /// Halfedge from `from` to `to`, if the edge exists. O(degree).
    #[must_use]
    pub fn find_halfedge(&self, from: VertexId, to: VertexId) -> Option<HalfedgeId> {
        self.outgoing_halfedges(from)
            .find(|&h| self.target(h) == to)
    }

    // ---- circulators ------------------------------------------------------

    /// Outgoing halfedges of a vertex in counter-clockwise order (seen with
    /// faces wound counter-clockwise), starting at the stored halfedge, so a
    /// boundary vertex starts at its boundary halfedge.
    pub fn outgoing_halfedges(&self, v: VertexId) -> impl Iterator<Item = HalfedgeId> + '_ {
        let start = self.vertex_halfedge(v);
        let mut current = start;
        core::iter::from_fn(move || {
            let h = current?;
            let rotated = self.rotate(h);
            current = (Some(rotated) != start).then_some(rotated);
            Some(h)
        })
    }

    /// Incoming halfedges of a vertex, the opposites of
    /// [`HalfedgeMesh::outgoing_halfedges`] in the same order.
    pub fn incoming_halfedges(&self, v: VertexId) -> impl Iterator<Item = HalfedgeId> + '_ {
        self.outgoing_halfedges(v).map(|h| self.opposite(h))
    }

    /// Neighbouring vertices in counter-clockwise order.
    pub fn vertex_vertices(&self, v: VertexId) -> impl Iterator<Item = VertexId> + '_ {
        self.outgoing_halfedges(v).map(|h| self.target(h))
    }

    /// Faces around a vertex in counter-clockwise order, holes skipped.
    pub fn vertex_faces(&self, v: VertexId) -> impl Iterator<Item = FaceId> + '_ {
        self.outgoing_halfedges(v).filter_map(|h| self.face(h))
    }

    /// Number of edges at a vertex. O(degree).
    #[must_use]
    pub fn degree(&self, v: VertexId) -> usize {
        self.outgoing_halfedges(v).count()
    }

    /// Halfedges of the loop containing `h` (a face or a hole), from `h`.
    pub fn loop_halfedges(&self, h: HalfedgeId) -> impl Iterator<Item = HalfedgeId> + '_ {
        let mut current = Some(h);
        core::iter::from_fn(move || {
            let g = current?;
            let next = self.next(g);
            current = (next != h).then_some(next);
            Some(g)
        })
    }

    /// Halfedges of a face's loop, from its stored halfedge.
    pub fn face_halfedges(&self, f: FaceId) -> impl Iterator<Item = HalfedgeId> + '_ {
        self.loop_halfedges(self.face_halfedge(f))
    }

    /// Corners of a face in winding order, starting at the source of its
    /// stored halfedge.
    pub fn face_vertices(&self, f: FaceId) -> impl Iterator<Item = VertexId> + '_ {
        self.face_halfedges(f).map(|h| self.source(h))
    }

    /// Faces across each edge of a face, `None` across a boundary.
    pub fn face_faces(&self, f: FaceId) -> impl Iterator<Item = Option<FaceId>> + '_ {
        self.face_halfedges(f).map(|h| self.face(self.opposite(h)))
    }

    /// Number of sides of a face. O(degree).
    #[must_use]
    pub fn face_degree(&self, f: FaceId) -> usize {
        self.face_halfedges(f).count()
    }

    // ---- maintenance ------------------------------------------------------

    /// Renumber live elements densely in their current order and drop the
    /// removed ones. Returns the old-to-new mapping so per-element data held
    /// by a caller can follow.
    pub fn compact(&mut self) -> HalfedgeRemap {
        let vertices = dense_map(&self.vertex_live);
        let edges = dense_map(&self.edge_live);
        let faces = dense_map(&self.face_live);
        let map_halfedge = |h: u32| -> u32 {
            if h == NONE {
                return NONE;
            }
            edges[(h / 2) as usize].map_or(NONE, |e| e * 2 + (h & 1))
        };
        let map = |table: &[Option<u32>], i: u32| -> u32 {
            if i == NONE {
                NONE
            } else {
                table[i as usize].unwrap_or(NONE)
            }
        };

        let mut positions = Vec::with_capacity(self.live_vertices);
        let mut vertex_out = Vec::with_capacity(self.live_vertices);
        for (v, slot) in vertices.iter().enumerate() {
            if slot.is_some() {
                positions.push(self.positions[v]);
                vertex_out.push(map_halfedge(self.vertex_out[v]));
            }
        }
        let mut links = Vec::with_capacity(self.live_edges * 2);
        for (h, link) in self.links.iter().enumerate() {
            if edges[h / 2].is_some() {
                links.push(Link {
                    next: map_halfedge(link.next),
                    prev: map_halfedge(link.prev),
                    target: map(&vertices, link.target),
                    face: map(&faces, link.face),
                });
            }
        }
        let mut face_halfedge = Vec::with_capacity(self.live_faces);
        for (f, slot) in faces.iter().enumerate() {
            if slot.is_some() {
                face_halfedge.push(map_halfedge(self.face_halfedge[f]));
            }
        }

        self.vertex_live = vec![true; positions.len()];
        self.edge_live = vec![true; links.len() / 2];
        self.face_live = vec![true; face_halfedge.len()];
        self.positions = positions;
        self.vertex_out = vertex_out;
        self.links = links;
        self.face_halfedge = face_halfedge;

        HalfedgeRemap {
            vertices: vertices.into_iter().map(|v| v.map(VertexId)).collect(),
            edges: edges.into_iter().map(|e| e.map(EdgeId)).collect(),
            faces: faces.into_iter().map(|f| f.map(FaceId)).collect(),
        }
    }

    /// Check every structural invariant listed in the
    /// [module documentation](self). O(halfedges).
    ///
    /// # Errors
    ///
    /// The first violated invariant found.
    pub fn validate(&self) -> Result<(), HalfedgeInvariantError> {
        check::validate(self)
    }

    // ---- internals --------------------------------------------------------

    /// Next outgoing halfedge counter-clockwise around `source(h)`.
    fn rotate(&self, h: HalfedgeId) -> HalfedgeId {
        self.opposite(self.prev(h))
    }

    fn link(&mut self, from: u32, to: u32) {
        self.links[from as usize].next = to;
        self.links[to as usize].prev = from;
    }
}

fn live_ids(live: &[bool]) -> impl Iterator<Item = u32> + '_ {
    live.iter()
        .enumerate()
        .filter(|(_, &alive)| alive)
        .map(|(i, _)| i as u32)
}

fn dense_map(live: &[bool]) -> Vec<Option<u32>> {
    let mut next = 0u32;
    live.iter()
        .map(|&alive| {
            alive.then(|| {
                next += 1;
                next - 1
            })
        })
        .collect()
}
