//! Local connectivity edits. Each one checks its preconditions before
//! touching anything, so a refused edit leaves the mesh exactly as it was,
//! and an accepted one leaves every invariant intact.

use core::fmt;

use axiolid_core::Point3;

use super::{EdgeId, FaceId, HalfedgeId, HalfedgeMesh, Link, VertexId, DEAD_LINK, NONE};

/// Why a local edit or a conversion was refused.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HalfedgeEditError {
    /// The id names an element that was removed or never existed.
    RemovedElement {
        /// `"vertex"`, `"halfedge"`, `"edge"` or `"face"`.
        element: &'static str,
        /// The id.
        index: u32,
    },
    /// The edit needs a face on both sides of the edge.
    BoundaryEdge {
        /// The edge.
        edge: EdgeId,
    },
    /// The edit needs a triangle here.
    NotATriangle {
        /// The face.
        face: FaceId,
        /// Its number of sides.
        degree: usize,
    },
    /// The edit would add an edge between vertices already joined, making
    /// the mesh non-simple.
    EdgeExists {
        /// One endpoint.
        a: VertexId,
        /// The other endpoint.
        b: VertexId,
    },
    /// Collapsing the edge would pinch the surface: its endpoints share a
    /// neighbour that is not opposite the edge, or both lie on a boundary the
    /// edge does not run along (Dey et al.'s link condition).
    LinkCondition {
        /// The edge.
        edge: EdgeId,
    },
    /// The edit would leave a degenerate surface: collapsing an edge of a
    /// tetrahedron, of a lone triangle or of a two-triangle pillow, or
    /// flipping an edge whose two opposite corners coincide.
    WouldDegenerate {
        /// The edge.
        edge: EdgeId,
    },
    /// The halfedge has a face; filling needs a boundary halfedge.
    NotBoundary {
        /// The halfedge.
        halfedge: HalfedgeId,
    },
    /// The two halfedges do not lie on the same face.
    NotSameFace {
        /// First halfedge.
        a: HalfedgeId,
        /// Second halfedge.
        b: HalfedgeId,
    },
    /// The two corners are the same or already joined by a side of the face.
    AdjacentCorners {
        /// First halfedge.
        a: HalfedgeId,
        /// Second halfedge.
        b: HalfedgeId,
    },
    /// The edit needs more elements than 32-bit ids can address.
    TooLarge,
}

impl fmt::Display for HalfedgeEditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RemovedElement { element, index } => {
                write!(f, "{element} {index} is not a live element")
            }
            Self::BoundaryEdge { edge } => write!(f, "edge {edge} lies on a boundary"),
            Self::NotATriangle { face, degree } => {
                write!(f, "face {face} has {degree} sides, not three")
            }
            Self::EdgeExists { a, b } => write!(f, "vertices {a} and {b} are already joined"),
            Self::LinkCondition { edge } => {
                write!(f, "collapsing edge {edge} violates the link condition")
            }
            Self::WouldDegenerate { edge } => {
                write!(f, "editing edge {edge} would leave a degenerate surface")
            }
            Self::NotBoundary { halfedge } => {
                write!(f, "halfedge {halfedge} is not a boundary halfedge")
            }
            Self::NotSameFace { a, b } => {
                write!(f, "halfedges {a} and {b} lie on different faces")
            }
            Self::AdjacentCorners { a, b } => write!(
                f,
                "the targets of halfedges {a} and {b} are already joined by a side"
            ),
            Self::TooLarge => write!(f, "mesh exceeds 32-bit element ids"),
        }
    }
}

impl std::error::Error for HalfedgeEditError {}

type EditResult<T> = Result<T, HalfedgeEditError>;

impl HalfedgeMesh {
    /// Replace an edge between two triangles by the other diagonal of the
    /// quadrilateral they form. Counts, and so the Euler characteristic, are
    /// unchanged; the edge keeps its id.
    ///
    /// # Errors
    ///
    /// [`HalfedgeEditError::BoundaryEdge`], [`HalfedgeEditError::NotATriangle`],
    /// [`HalfedgeEditError::WouldDegenerate`] when the two opposite corners are
    /// one vertex, and [`HalfedgeEditError::EdgeExists`] when they are already
    /// joined.
    pub fn flip_edge(&mut self, e: EdgeId) -> EditResult<()> {
        self.require_edge(e)?;
        let h = e.0 * 2;
        let o = h + 1;
        let (f1, f2) = (self.links[h as usize].face, self.links[o as usize].face);
        if f1 == NONE || f2 == NONE {
            return Err(HalfedgeEditError::BoundaryEdge { edge: e });
        }
        self.require_triangle(FaceId(f1))?;
        self.require_triangle(FaceId(f2))?;
        let (hn, hp) = (self.links[h as usize].next, self.links[h as usize].prev);
        let (on, op) = (self.links[o as usize].next, self.links[o as usize].prev);
        let a = self.links[o as usize].target;
        let b = self.links[h as usize].target;
        let c = self.links[hn as usize].target;
        let d = self.links[on as usize].target;
        if c == d {
            return Err(HalfedgeEditError::WouldDegenerate { edge: e });
        }
        if self.find_halfedge(VertexId(c), VertexId(d)).is_some() {
            return Err(HalfedgeEditError::EdgeExists {
                a: VertexId(c),
                b: VertexId(d),
            });
        }

        // Quad a -> d -> b -> c; new diagonal h: d -> c, o: c -> d.
        self.link(hp, on);
        self.link(on, h);
        self.link(h, hp);
        self.link(op, hn);
        self.link(hn, o);
        self.link(o, op);
        self.links[h as usize].target = c;
        self.links[o as usize].target = d;
        self.links[on as usize].face = f1;
        self.links[hn as usize].face = f2;
        self.face_halfedge[f1 as usize] = h;
        self.face_halfedge[f2 as usize] = o;
        if self.vertex_out[a as usize] == h {
            self.vertex_out[a as usize] = on;
        }
        if self.vertex_out[b as usize] == o {
            self.vertex_out[b as usize] = hn;
        }
        Ok(())
    }

    /// Insert a vertex at `position` on an edge and join it to the opposite
    /// corner of each adjacent triangle, so a triangle mesh stays a triangle
    /// mesh. A non-triangular adjacent face just gains a side. The Euler
    /// characteristic is unchanged.
    ///
    /// The old edge keeps its id and runs from its original source to the
    /// new vertex.
    ///
    /// # Errors
    ///
    /// [`HalfedgeEditError::RemovedElement`] and [`HalfedgeEditError::TooLarge`].
    pub fn split_edge(&mut self, e: EdgeId, position: Point3) -> EditResult<VertexId> {
        self.require_edge(e)?;
        self.require_room(1, 3, 2)?;
        let h = e.0 * 2;
        let o = h + 1;
        let fh = self.links[h as usize].face;
        let fo = self.links[o as usize].face;
        let split_h = fh != NONE && self.face_degree(FaceId(fh)) == 3;
        let split_o = fo != NONE && self.face_degree(FaceId(fo)) == 3;
        let b = self.links[h as usize].target;
        let hn = self.links[h as usize].next;
        let op = self.links[o as usize].prev;

        let m = self.new_vertex(position);
        let g = self.new_edge(m, b);
        // h: a -> m, g: m -> b on h's side; g^1: b -> m, o: m -> a on o's side.
        self.links[h as usize].target = m;
        self.link(h, g);
        self.link(g, hn);
        self.links[g as usize].face = fh;
        self.link(op, g ^ 1);
        self.link(g ^ 1, o);
        self.links[(g ^ 1) as usize].face = fo;
        self.vertex_out[m as usize] = if fo == NONE { o } else { g };
        if self.vertex_out[b as usize] == o {
            self.vertex_out[b as usize] = g ^ 1;
        }
        if split_h {
            self.insert_diagonal(h, hn);
        }
        if split_o {
            let on = self.links[o as usize].next;
            self.insert_diagonal(g ^ 1, on);
        }
        Ok(VertexId(m))
    }

    /// Collapse halfedge `h` by merging its source into its target, which
    /// keeps its id and position (move it with
    /// [`HalfedgeMesh::set_position`]). The incident triangles disappear and
    /// their remaining sides merge pairwise. Removes one vertex, three edges
    /// on an interior edge (two on a boundary edge) and the adjacent faces,
    /// so the Euler characteristic is unchanged.
    ///
    /// # Errors
    ///
    /// [`HalfedgeEditError::NotATriangle`] when a face around the removed
    /// vertex is not a triangle, [`HalfedgeEditError::LinkCondition`] when the
    /// collapse would pinch the surface, and
    /// [`HalfedgeEditError::WouldDegenerate`] for a tetrahedron, a lone
    /// triangle or a two-triangle pillow.
    pub fn collapse_edge(&mut self, h: HalfedgeId) -> EditResult<VertexId> {
        self.require_edge(self.edge(h))?;
        let edge = self.edge(h);
        let h = h.0;
        let o = h ^ 1;
        let a = VertexId(self.links[o as usize].target);
        let b = VertexId(self.links[h as usize].target);
        let fh = self.links[h as usize].face;
        let fo = self.links[o as usize].face;
        for f in self.vertex_faces(a) {
            self.require_triangle(f)?;
        }
        let opposite_corner = |mesh: &Self, side: u32| {
            (mesh.links[side as usize].face != NONE)
                .then(|| VertexId(mesh.links[mesh.links[side as usize].next as usize].target))
        };
        let c = opposite_corner(self, h);
        let d = opposite_corner(self, o);
        if c.is_some() && c == d {
            return Err(HalfedgeEditError::WouldDegenerate { edge });
        }

        // Link condition: every common neighbour is a corner opposite the
        // edge, and the boundary (a virtual vertex joined to every boundary
        // vertex) is a common neighbour only when the edge runs along it.
        let around_a: Vec<VertexId> = self.vertex_vertices(a).collect();
        if self
            .vertex_vertices(b)
            .any(|n| n != a && around_a.contains(&n) && Some(n) != c && Some(n) != d)
        {
            return Err(HalfedgeEditError::LinkCondition { edge });
        }
        if self.is_boundary_vertex(a) && self.is_boundary_vertex(b) && !self.is_boundary_edge(edge)
        {
            return Err(HalfedgeEditError::LinkCondition { edge });
        }
        // The link condition holds on the boundary of a tetrahedron, which
        // still collapses to nothing; with the boundary as a virtual vertex,
        // a lone triangle is that case too.
        match (c, d) {
            (Some(c), Some(d)) => {
                if [a, b, c, d]
                    .iter()
                    .all(|&v| !self.is_boundary_vertex(v) && self.degree(v) == 3)
                {
                    return Err(HalfedgeEditError::WouldDegenerate { edge });
                }
            }
            _ => {
                let hole = if fh == NONE { h } else { o };
                if self.loop_halfedges(HalfedgeId(hole)).count() == 3 {
                    return Err(HalfedgeEditError::WouldDegenerate { edge });
                }
            }
        }

        let incoming: Vec<u32> = self.incoming_halfedges(a).map(|g| g.0).collect();
        let (hn, hp) = (self.links[h as usize].next, self.links[h as usize].prev);
        let (on, op) = (self.links[o as usize].next, self.links[o as usize].prev);
        for g in incoming {
            self.links[g as usize].target = b.0;
        }
        self.link(hp, hn);
        self.link(op, on);
        if self.vertex_out[b.index()] == o {
            self.vertex_out[b.index()] = hn;
        }
        if fh != NONE {
            self.remove_loop(hn);
        }
        if fo != NONE {
            self.remove_loop(on);
        }
        self.remove_edge(edge.0);
        self.remove_vertex(a.0);
        for v in [Some(b), c, d].into_iter().flatten() {
            self.adjust_outgoing(v.0);
        }
        Ok(b)
    }

    /// Insert a vertex at `position` inside a face and fan the face into
    /// triangles around it, one per side. The face keeps its id as the
    /// triangle on its stored halfedge. The Euler characteristic is
    /// unchanged.
    ///
    /// # Errors
    ///
    /// [`HalfedgeEditError::RemovedElement`] and [`HalfedgeEditError::TooLarge`].
    pub fn split_face(&mut self, f: FaceId, position: Point3) -> EditResult<VertexId> {
        self.require_face(f)?;
        let sides: Vec<u32> = self.face_halfedges(f).map(|h| h.0).collect();
        let n = sides.len();
        self.require_room(1, n, n)?;
        let center = self.new_vertex(position);
        // spoke[i]: center -> source(sides[i]); spoke[i] ^ 1 runs back.
        let spokes: Vec<u32> = sides
            .iter()
            .map(|&h| {
                let corner = self.links[(h ^ 1) as usize].target;
                self.new_edge(center, corner)
            })
            .collect();
        for i in 0..n {
            let side = sides[i];
            let up = spokes[(i + 1) % n] ^ 1;
            let down = spokes[i];
            let face = if i == 0 { f.0 } else { self.new_face(side) };
            self.link(side, up);
            self.link(up, down);
            self.link(down, side);
            for g in [side, up, down] {
                self.links[g as usize].face = face;
            }
        }
        self.face_halfedge[f.index()] = sides[0];
        self.vertex_out[center as usize] = spokes[0];
        Ok(VertexId(center))
    }

    /// Split a face by a new edge from `target(a)` to `target(b)`, two
    /// corners of the same face that are not already joined. Returns the new
    /// halfedge running from `target(a)` to `target(b)`; it stays on the
    /// original face, which keeps `a`, and the other part becomes a new
    /// face. The Euler characteristic is unchanged.
    ///
    /// This triangulates a polygon one diagonal at a time, for example a
    /// hole closed with [`HalfedgeMesh::fill_hole`].
    ///
    /// # Errors
    ///
    /// [`HalfedgeEditError::NotSameFace`] (boundary halfedges included),
    /// [`HalfedgeEditError::AdjacentCorners`], and
    /// [`HalfedgeEditError::EdgeExists`] when the corners are joined through
    /// another face.
    pub fn split_face_diagonal(&mut self, a: HalfedgeId, b: HalfedgeId) -> EditResult<HalfedgeId> {
        self.require_edge(self.edge(a))?;
        self.require_edge(self.edge(b))?;
        let face = self.links[a.index()].face;
        if face == NONE || self.links[b.index()].face != face {
            return Err(HalfedgeEditError::NotSameFace { a, b });
        }
        if a == b || self.next(a) == b || self.next(b) == a {
            return Err(HalfedgeEditError::AdjacentCorners { a, b });
        }
        let (u, w) = (self.target(a), self.target(b));
        if self.find_halfedge(u, w).is_some() {
            return Err(HalfedgeEditError::EdgeExists { a: u, b: w });
        }
        self.require_room(0, 1, 1)?;
        Ok(HalfedgeId(self.insert_diagonal(a.0, b.0)))
    }

    /// Close the hole bounded by boundary halfedge `h` with one new face.
    /// The Euler characteristic rises by one.
    ///
    /// # Errors
    ///
    /// [`HalfedgeEditError::NotBoundary`] when `h` already has a face.
    pub fn fill_hole(&mut self, h: HalfedgeId) -> EditResult<FaceId> {
        self.require_edge(self.edge(h))?;
        if !self.is_boundary_halfedge(h) {
            return Err(HalfedgeEditError::NotBoundary { halfedge: h });
        }
        self.require_room(0, 0, 1)?;
        let face = self.new_face(h.0);
        let sides: Vec<HalfedgeId> = self.loop_halfedges(h).collect();
        for g in sides {
            self.links[g.index()].face = face;
        }
        Ok(FaceId(face))
    }

    // ---- preconditions ----------------------------------------------------

    fn require_edge(&self, e: EdgeId) -> EditResult<()> {
        if self.contains_edge(e) {
            Ok(())
        } else {
            Err(HalfedgeEditError::RemovedElement {
                element: "edge",
                index: e.0,
            })
        }
    }

    fn require_face(&self, f: FaceId) -> EditResult<()> {
        if self.contains_face(f) {
            Ok(())
        } else {
            Err(HalfedgeEditError::RemovedElement {
                element: "face",
                index: f.0,
            })
        }
    }

    fn require_triangle(&self, f: FaceId) -> EditResult<()> {
        let degree = self.face_degree(f);
        if degree == 3 {
            Ok(())
        } else {
            Err(HalfedgeEditError::NotATriangle { face: f, degree })
        }
    }

    fn require_room(&self, vertices: usize, edges: usize, faces: usize) -> EditResult<()> {
        let limit = NONE as usize;
        if self.positions.len() + vertices >= limit
            || (self.edge_live.len() + edges) * 2 >= limit
            || self.face_live.len() + faces >= limit
        {
            Err(HalfedgeEditError::TooLarge)
        } else {
            Ok(())
        }
    }

    // ---- primitive mutations ----------------------------------------------

    fn new_vertex(&mut self, position: Point3) -> u32 {
        self.positions.push(position);
        self.vertex_out.push(NONE);
        self.vertex_live.push(true);
        self.live_vertices += 1;
        (self.positions.len() - 1) as u32
    }

    /// New edge with dangling links; returns the halfedge `from -> to`.
    fn new_edge(&mut self, from: u32, to: u32) -> u32 {
        let h = self.links.len() as u32;
        self.links.push(Link {
            target: to,
            ..DEAD_LINK
        });
        self.links.push(Link {
            target: from,
            ..DEAD_LINK
        });
        self.edge_live.push(true);
        self.live_edges += 1;
        h
    }

    fn new_face(&mut self, h: u32) -> u32 {
        self.face_halfedge.push(h);
        self.face_live.push(true);
        self.live_faces += 1;
        (self.face_halfedge.len() - 1) as u32
    }

    fn remove_vertex(&mut self, v: u32) {
        self.vertex_live[v as usize] = false;
        self.vertex_out[v as usize] = NONE;
        self.live_vertices -= 1;
    }

    fn remove_edge(&mut self, e: u32) {
        self.edge_live[e as usize] = false;
        self.links[(e * 2) as usize] = DEAD_LINK;
        self.links[(e * 2 + 1) as usize] = DEAD_LINK;
        self.live_edges -= 1;
    }

    fn remove_face(&mut self, f: u32) {
        self.face_live[f as usize] = false;
        self.face_halfedge[f as usize] = NONE;
        self.live_faces -= 1;
    }

    /// Join `target(h1)` to `target(h2)` across their shared face. The face
    /// keeps `h1` and the new halfedge; the other part becomes a new face.
    fn insert_diagonal(&mut self, h1: u32, h2: u32) -> u32 {
        let face = self.links[h1 as usize].face;
        let u = self.links[h1 as usize].target;
        let w = self.links[h2 as usize].target;
        let n1 = self.links[h1 as usize].next;
        let n2 = self.links[h2 as usize].next;
        let d = self.new_edge(u, w);
        self.link(h1, d);
        self.link(d, n2);
        self.links[d as usize].face = face;
        self.link(h2, d ^ 1);
        self.link(d ^ 1, n1);
        let other = self.new_face(d ^ 1);
        let sides: Vec<HalfedgeId> = self.loop_halfedges(HalfedgeId(d ^ 1)).collect();
        for g in sides {
            self.links[g.index()].face = other;
        }
        self.face_halfedge[face as usize] = h1;
        d
    }

    /// Dissolve a two-sided loop `h0 -> h1 -> h0` left by a collapse: `h0`
    /// takes the place of `opposite(h1)`, whose edge and the loop's face are
    /// removed, so the two sides merge into `h0`'s edge.
    fn remove_loop(&mut self, h0: u32) {
        let h1 = self.links[h0 as usize].next;
        debug_assert_eq!(self.links[h1 as usize].next, h0);
        let o0 = h0 ^ 1;
        let o1 = h1 ^ 1;
        let face = self.links[h0 as usize].face;
        let v0 = self.links[h0 as usize].target;
        let v1 = self.links[h1 as usize].target;
        let Link {
            prev,
            next,
            face: outer,
            ..
        } = self.links[o1 as usize];
        self.link(prev, h0);
        self.link(h0, next);
        self.links[h0 as usize].face = outer;
        if outer != NONE && self.face_halfedge[outer as usize] == o1 {
            self.face_halfedge[outer as usize] = h0;
        }
        if self.vertex_out[v0 as usize] == h1 {
            self.vertex_out[v0 as usize] = o0;
        }
        if self.vertex_out[v1 as usize] == o1 {
            self.vertex_out[v1 as usize] = h0;
        }
        self.remove_face(face);
        self.remove_edge(h1 / 2);
    }

    /// Restore the rule that a boundary vertex stores its boundary halfedge.
    fn adjust_outgoing(&mut self, v: u32) {
        if self.vertex_out[v as usize] == NONE {
            return;
        }
        let boundary = self
            .outgoing_halfedges(VertexId(v))
            .find(|&h| self.is_boundary_halfedge(h));
        if let Some(h) = boundary {
            self.vertex_out[v as usize] = h.0;
        }
    }
}
