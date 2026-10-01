//! Incremental 3D Delaunay triangulation (Bowyer-Watson).
//!
//! # Structure
//!
//! The triangulation covers all of space: besides the finite tetrahedra it
//! keeps one *infinite* cell per convex-hull face, joining that face to a
//! symbolic vertex at infinity. Every cell then has exactly four neighbours,
//! so insertion never special-cases the hull: a point outside it conflicts
//! with the infinite cells over the faces it sees, and those are replaced
//! like any other.
//!
//! A cell stores its vertices and, at index `i`, the neighbour across the
//! face opposite vertex `i`. Finite cells are positively oriented
//! ([`orientation`] > 0). An infinite cell is oriented as if its infinite
//! vertex were a point beyond its hull face.
//!
//! # Insertion
//!
//! 1. Locate: a stochastic visibility walk from the last created cell steps
//!    across any face the point lies strictly beyond, in random face order,
//!    until no face separates them; it ends in the finite cell containing the
//!    point (closed) or in the infinite cell of a hull face the point sees.
//!    A point equal to a vertex of the closed cell is a duplicate.
//! 2. Conflict region: every cell whose perturbed circumsphere contains the
//!    point, grown by breadth-first search from the located cell. Under the
//!    perturbation of [`crate::sos`] no point lies *on* a sphere, so the
//!    region is exactly the set of cells the point invalidates, and it is
//!    star-shaped from the point.
//! 3. Re-triangulate: delete the region and join the point to each of its
//!    boundary faces.

use std::collections::HashMap;
use std::fmt;

use axiolid_core::Point3;
use axiolid_guarantees::Sign;

use crate::hilbert;
use crate::sos::{beyond_face, collinear, in_sphere, orientation};

/// The symbolic vertex at infinity.
const INFINITE: u32 = u32::MAX;

/// Smallest non-zero coordinate magnitude accepted: `2^-100`.
///
/// The in-sphere test is a degree-5 polynomial and the coplanar circle test
/// degree 6. Keeping every non-zero coordinate in `[2^-100, 2^100]` keeps
/// every product those exact evaluations form inside binary64's normal
/// range, which is what makes them exact.
pub const MIN_COORDINATE: f64 = 7.888_609_052_210_118e-31;

/// Largest coordinate magnitude accepted: `2^100`.
pub const MAX_COORDINATE: f64 = 1.267_650_600_228_229_4e30;

/// Why a point or a point set was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Delaunay3Error {
    /// The point at `index` has a NaN or infinite coordinate.
    NonFinite {
        /// Index of the offending point.
        index: usize,
    },
    /// The point at `index` has a non-zero coordinate whose magnitude lies
    /// outside [`MIN_COORDINATE`]`..=`[`MAX_COORDINATE`], where the exact
    /// predicates cannot be evaluated exactly.
    OutOfRange {
        /// Index of the offending point.
        index: usize,
    },
    /// The points span fewer than three dimensions, so there is no
    /// tetrahedron: `None` for no points, else the affine dimension (0 for
    /// one distinct point, 1 collinear, 2 coplanar).
    NotFullDimensional {
        /// Affine dimension of the input, `None` when it is empty.
        dimension: Option<usize>,
    },
}

impl fmt::Display for Delaunay3Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite { index } => write!(f, "point {index} has a non-finite coordinate"),
            Self::OutOfRange { index } => write!(
                f,
                "point {index} has a coordinate outside the exactly decidable range 2^-100..=2^100"
            ),
            Self::NotFullDimensional { dimension: None } => write!(f, "no points"),
            Self::NotFullDimensional {
                dimension: Some(dimension),
            } => write!(
                f,
                "the points span {dimension} dimension(s); a tetrahedralization needs 3"
            ),
        }
    }
}

impl std::error::Error for Delaunay3Error {}

/// What [`Delaunay3::insert`] did with a point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Insertion {
    /// The point became a vertex; its index in [`Delaunay3::points`].
    Vertex(usize),
    /// The point repeats an existing vertex (exactly equal coordinates) and
    /// was not inserted again; the existing vertex's index.
    Duplicate(usize),
}

/// The finite tetrahedra of a triangulation, with their adjacency.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Tetrahedra {
    /// Vertex indices into [`Delaunay3::points`], positively oriented: the
    /// signed volume `(b - a) . ((c - a) x (d - a))` of `[a, b, c, d]` is
    /// positive.
    pub tetrahedra: Vec<[usize; 4]>,
    /// `neighbors[t][i]` is the tetrahedron across the face opposite vertex
    /// `i` of tetrahedron `t`, or `None` on the convex hull.
    pub neighbors: Vec<[Option<usize>; 4]>,
}

#[derive(Debug, Clone, Copy)]
struct Cell {
    vertices: [u32; 4],
    neighbors: [u32; 4],
}

impl Cell {
    fn infinite_index(&self) -> Option<usize> {
        self.vertices.iter().position(|&v| v == INFINITE)
    }
}

/// An incremental 3D Delaunay triangulation of a point set.
///
/// Built with [`Delaunay3::from_points`] or point by point with
/// [`Delaunay3::insert`]. Every decision is an exact sign from
/// `axiolid-predicates`; ties (cospherical points, points in the plane of a
/// hull face) are broken by a symbolic perturbation that depends only on the
/// coordinates, so the result is a Delaunay triangulation -- every
/// circumsphere has no point strictly inside -- and it is the same for every
/// insertion order. Exact duplicates are merged, not inserted twice.
#[derive(Debug, Clone)]
pub struct Delaunay3 {
    points: Vec<Point3>,
    /// Per point, the vertex it is (itself) or duplicates.
    vertex: Vec<u32>,
    cells: Vec<Cell>,
    alive: Vec<bool>,
    free: Vec<u32>,
    marks: Vec<u64>,
    epoch: u64,
    hint: u32,
    rng: u64,
    /// Affinely independent points found so far; four means 3D.
    frame: Vec<u32>,
    /// Distinct points held back until the first tetrahedron exists.
    pending: Vec<u32>,
    seen: HashMap<[u64; 3], u32>,
    scratch: Scratch,
}

/// Buffers one insertion uses, kept to avoid reallocating them.
#[derive(Debug, Clone, Default)]
struct Scratch {
    stack: Vec<u32>,
    region: Vec<u32>,
    /// New cell, slot of the inserted point, outside neighbour, and the
    /// neighbour's slot facing the cavity.
    created: Vec<(Cell, usize, u32, usize)>,
    ids: Vec<u32>,
    /// Open-addressing table: edge key, then new cell and slot.
    edges: Vec<(u64, u64)>,
}

impl Scratch {
    fn clear(&mut self) {
        self.stack.clear();
        self.region.clear();
        self.created.clear();
        self.ids.clear();
        self.edges.clear();
    }
}

impl Default for Delaunay3 {
    fn default() -> Self {
        Self::new()
    }
}

/// Reject non-finite and out-of-range coordinates; fold `-0.0` into `0.0`
/// so that equal points compare equal bit for bit.
fn check_point(p: Point3, index: usize) -> Result<Point3, Delaunay3Error> {
    for value in [p.x, p.y, p.z] {
        if !value.is_finite() {
            return Err(Delaunay3Error::NonFinite { index });
        }
        if value != 0.0 && !(MIN_COORDINATE..=MAX_COORDINATE).contains(&value.abs()) {
            return Err(Delaunay3Error::OutOfRange { index });
        }
    }
    Ok(Point3::new(p.x + 0.0, p.y + 0.0, p.z + 0.0))
}

fn key(p: Point3) -> [u64; 3] {
    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
}

impl Delaunay3 {
    /// An empty triangulation.
    #[must_use]
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            vertex: Vec::new(),
            cells: Vec::new(),
            alive: Vec::new(),
            free: Vec::new(),
            marks: Vec::new(),
            epoch: 0,
            hint: 0,
            rng: 0x9E37_79B9_7F4A_7C15,
            frame: Vec::new(),
            pending: Vec::new(),
            seen: HashMap::new(),
            scratch: Scratch::default(),
        }
    }

    /// The Delaunay triangulation of `points`.
    ///
    /// Points keep their indices: [`Delaunay3::points`] is `points` (with
    /// `-0.0` read as `0.0`). They are inserted in Hilbert-curve order, which
    /// changes only the running time.
    ///
    /// # Errors
    ///
    /// [`Delaunay3Error::NonFinite`] or [`Delaunay3Error::OutOfRange`] for the
    /// first offending point, and [`Delaunay3Error::NotFullDimensional`] when
    /// the points are empty, all equal, collinear or coplanar.
    pub fn from_points(points: &[Point3]) -> Result<Self, Delaunay3Error> {
        let mut triangulation = Self::new();
        triangulation.points.reserve(points.len());
        for (index, &p) in points.iter().enumerate() {
            triangulation.points.push(check_point(p, index)?);
        }
        let count = u32::try_from(points.len()).map_err(|_| Delaunay3Error::OutOfRange {
            index: u32::MAX as usize,
        })?;
        triangulation.vertex = (0..count).collect();
        let all: Vec<u32> = (0..count).collect();
        for index in hilbert::order(&triangulation.points, &all) {
            triangulation.insert_index(index);
        }
        if triangulation.frame.len() < 4 {
            return Err(Delaunay3Error::NotFullDimensional {
                dimension: triangulation.dimension(),
            });
        }
        Ok(triangulation)
    }

    /// Insert one point.
    ///
    /// Before four affinely independent points have arrived there is no
    /// tetrahedron; points are held and inserted as soon as one exists.
    ///
    /// # Errors
    ///
    /// [`Delaunay3Error::NonFinite`] or [`Delaunay3Error::OutOfRange`], with
    /// `index` the index the point would have had; the triangulation is left
    /// unchanged.
    pub fn insert(&mut self, point: Point3) -> Result<Insertion, Delaunay3Error> {
        let index = self.points.len();
        let point = check_point(point, index)?;
        let id = u32::try_from(index)
            .ok()
            .filter(|&id| id != INFINITE)
            .ok_or(Delaunay3Error::OutOfRange { index })?;
        self.points.push(point);
        self.vertex.push(id);
        Ok(self.insert_index(id))
    }

    /// Every point given so far, duplicates included, by index.
    #[must_use]
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// The vertex that point `index` became: `index` itself, or the earlier
    /// point it exactly duplicates.
    ///
    /// # Panics
    ///
    /// If `index` is not a point index.
    #[must_use]
    pub fn vertex_of(&self, index: usize) -> usize {
        self.vertex[index] as usize
    }

    /// Affine dimension of the points so far: `None` when there are none.
    #[must_use]
    pub fn dimension(&self) -> Option<usize> {
        self.frame.len().checked_sub(1)
    }

    /// Number of finite tetrahedra.
    #[must_use]
    pub fn tetrahedron_count(&self) -> usize {
        self.live_cells()
            .filter(|&c| self.cells[c].infinite_index().is_none())
            .count()
    }

    /// The finite tetrahedra and their adjacency.
    #[must_use]
    pub fn tetrahedra(&self) -> Tetrahedra {
        let mut id = vec![usize::MAX; self.cells.len()];
        let mut tetrahedra = Vec::new();
        for c in self.live_cells() {
            let cell = &self.cells[c];
            if cell.infinite_index().is_none() {
                id[c] = tetrahedra.len();
                tetrahedra.push(cell.vertices.map(|v| v as usize));
            }
        }
        let neighbors = self
            .live_cells()
            .filter(|&c| id[c] != usize::MAX)
            .map(|c| {
                self.cells[c].neighbors.map(|n| {
                    let n = id[n as usize];
                    (n != usize::MAX).then_some(n)
                })
            })
            .collect();
        Tetrahedra {
            tetrahedra,
            neighbors,
        }
    }

    /// The convex hull's faces, as vertex triples ordered counter-clockwise
    /// seen from outside: the normal `(b - a) x (c - a)` points out.
    ///
    /// Coplanar hull facets come triangulated (Delaunay within their plane,
    /// under the same perturbation). Empty below three dimensions.
    #[must_use]
    pub fn hull_triangles(&self) -> Vec<[usize; 3]> {
        self.live_cells()
            .filter_map(|c| {
                let cell = &self.cells[c];
                let k = cell.infinite_index()?;
                let mut order = cell.vertices;
                // Move the infinite vertex to the last slot; one transposition
                // flips the orientation, a second among the finite three
                // restores it.
                order.swap(k, 3);
                if k != 3 {
                    order.swap(0, 1);
                }
                Some([order[0], order[1], order[2]].map(|v| v as usize))
            })
            .collect()
    }

    fn live_cells(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.cells.len()).filter(|&c| self.alive[c])
    }

    fn point(&self, v: u32) -> Point3 {
        self.points[v as usize]
    }

    fn next_random(&mut self) -> u64 {
        self.rng ^= self.rng << 13;
        self.rng ^= self.rng >> 7;
        self.rng ^= self.rng << 17;
        self.rng
    }

    /// Insert a stored, checked point by index.
    fn insert_index(&mut self, id: u32) -> Insertion {
        if self.frame.len() == 4 {
            return self.insert_3d(id);
        }
        let p = self.point(id);
        if let Some(&existing) = self.seen.get(&key(p)) {
            self.vertex[id as usize] = existing;
            return Insertion::Duplicate(existing as usize);
        }
        self.seen.insert(key(p), id);
        self.pending.push(id);
        let raises = match self.frame.len() {
            0 | 1 => true,
            2 => !collinear(self.point(self.frame[0]), self.point(self.frame[1]), p),
            _ => {
                let [a, b, c] =
                    [self.frame[0], self.frame[1], self.frame[2]].map(|v| self.point(v));
                orientation(a, b, c, p) != Sign::Zero
            }
        };
        if raises {
            self.frame.push(id);
        }
        if self.frame.len() == 4 {
            self.build_first_tetrahedron();
            let pending = std::mem::take(&mut self.pending);
            self.seen = HashMap::new();
            let rest: Vec<u32> = pending
                .into_iter()
                .filter(|v| !self.frame.contains(v))
                .collect();
            for v in hilbert::order(&self.points, &rest) {
                let inserted = self.insert_3d(v);
                debug_assert_eq!(inserted, Insertion::Vertex(v as usize));
            }
        }
        Insertion::Vertex(id as usize)
    }

    /// One positive tetrahedron and the four infinite cells over its faces.
    fn build_first_tetrahedron(&mut self) {
        let [a, b, mut c, mut d] = [self.frame[0], self.frame[1], self.frame[2], self.frame[3]];
        if orientation(self.point(a), self.point(b), self.point(c), self.point(d)) == Sign::Negative
        {
            std::mem::swap(&mut c, &mut d);
        }
        let finite = [a, b, c, d];
        let mut cells = vec![finite];
        for i in 0..4 {
            let mut ghost = finite;
            ghost[i] = INFINITE;
            // Beyond face i lies opposite vertex i, the negative side: swap
            // two finite slots so the infinite vertex reads as beyond.
            let others: Vec<usize> = (0..4).filter(|&j| j != i).collect();
            ghost.swap(others[0], others[1]);
            cells.push(ghost);
        }
        for vertices in &cells {
            self.allocate(Cell {
                vertices: *vertices,
                neighbors: [INFINITE; 4],
            });
        }
        // Five cells: match faces directly.
        for x in 0..5 {
            for i in 0..4 {
                let face = face_key(cells[x], i);
                for (y, other) in cells.iter().enumerate() {
                    if y == x {
                        continue;
                    }
                    if let Some(j) = (0..4).find(|&j| face_key(*other, j) == face) {
                        self.cells[x].neighbors[i] = y as u32;
                        debug_assert!(j < 4);
                    }
                }
            }
        }
        self.hint = 0;
    }

    fn allocate(&mut self, cell: Cell) -> u32 {
        if let Some(c) = self.free.pop() {
            self.cells[c as usize] = cell;
            self.alive[c as usize] = true;
            c
        } else {
            self.cells.push(cell);
            self.alive.push(true);
            self.marks.push(0);
            (self.cells.len() - 1) as u32
        }
    }

    /// Orientation of cell `c` with its vertex at slot `slot` replaced by
    /// point `p`; the slot may hold the infinite vertex.
    fn orientation_with(&self, c: u32, slot: usize, p: u32) -> Sign {
        let mut vertices = self.cells[c as usize].vertices;
        vertices[slot] = p;
        if vertices.contains(&INFINITE) {
            return Sign::Zero;
        }
        let [a, b, x, y] = vertices.map(|v| self.point(v));
        orientation(a, b, x, y)
    }

    /// The cell containing point `p`: a finite cell whose closure holds it,
    /// or an infinite cell whose hull face it lies strictly beyond.
    fn locate(&mut self, p: u32) -> u32 {
        let mut c = if self.alive[self.hint as usize] {
            self.hint
        } else {
            self.live_cells().next().map_or(0, |c| c as u32)
        };
        let mut previous = INFINITE;
        let budget = 4 * self.cells.len() + 64;
        for _ in 0..budget {
            let cell = self.cells[c as usize];
            if let Some(k) = cell.infinite_index() {
                if self.orientation_with(c, k, p) == Sign::Positive {
                    return c;
                }
                previous = c;
                c = cell.neighbors[k];
                continue;
            }
            let start = (self.next_random() % 4) as usize;
            let mut moved = false;
            for t in 0..4 {
                let i = (start + t) % 4;
                let next = cell.neighbors[i];
                if next == previous {
                    continue;
                }
                if self.orientation_with(c, i, p) == Sign::Negative {
                    previous = c;
                    c = next;
                    moved = true;
                    break;
                }
            }
            if !moved {
                return c;
            }
        }
        self.locate_exhaustively(p)
    }

    /// Exhaustive location: the guarantee behind the walk's step budget.
    fn locate_exhaustively(&self, p: u32) -> u32 {
        let mut beyond = None;
        for c in self.live_cells() {
            let cell = &self.cells[c];
            match cell.infinite_index() {
                Some(k) => {
                    if beyond.is_none() && self.orientation_with(c as u32, k, p) == Sign::Positive {
                        beyond = Some(c as u32);
                    }
                }
                None => {
                    if (0..4).all(|i| self.orientation_with(c as u32, i, p) != Sign::Negative) {
                        return c as u32;
                    }
                }
            }
        }
        beyond.unwrap_or(self.hint)
    }

    fn in_conflict(&self, c: u32, p: u32) -> bool {
        let cell = &self.cells[c as usize];
        match cell.infinite_index() {
            None => in_sphere(&self.points, cell.vertices, p),
            Some(k) => {
                // The face in the orientation that makes "beyond" positive:
                // the infinite slot moved last, as in `hull_triangles`.
                let mut order = cell.vertices;
                order.swap(k, 3);
                if k != 3 {
                    order.swap(0, 1);
                }
                beyond_face(&self.points, [order[0], order[1], order[2]], p)
            }
        }
    }

    fn insert_3d(&mut self, p: u32) -> Insertion {
        let start = self.locate(p);
        let located = self.cells[start as usize];
        if located.infinite_index().is_none() {
            let point = self.point(p);
            if let Some(&v) = located
                .vertices
                .iter()
                .find(|&&v| key(self.point(v)) == key(point))
            {
                self.vertex[p as usize] = v;
                return Insertion::Duplicate(v as usize);
            }
        }

        // Conflict region by depth-first search. Marks: 2*epoch + 1 in
        // conflict, 2*epoch tested and not. Buffers are reused across
        // insertions.
        let mut scratch = std::mem::take(&mut self.scratch);
        scratch.clear();
        self.epoch += 1;
        let conflict = 2 * self.epoch + 1;
        let clear = 2 * self.epoch;
        scratch.stack.push(start);
        scratch.region.push(start);
        self.marks[start as usize] = conflict;
        while let Some(c) = scratch.stack.pop() {
            for i in 0..4 {
                let n = self.cells[c as usize].neighbors[i];
                let mark = self.marks[n as usize];
                if mark == conflict {
                    continue;
                }
                if mark != clear && self.in_conflict(n, p) {
                    self.marks[n as usize] = conflict;
                    scratch.stack.push(n);
                    scratch.region.push(n);
                } else {
                    self.marks[n as usize] = clear;
                    // Boundary face: slot i of conflict cell c, outside n.
                    let j = self.cells[n as usize]
                        .neighbors
                        .iter()
                        .position(|&m| m == c)
                        .expect("adjacency is symmetric");
                    let mut vertices = self.cells[c as usize].vertices;
                    vertices[i] = p;
                    let mut neighbors = [INFINITE; 4];
                    neighbors[i] = n;
                    scratch.created.push((
                        Cell {
                            vertices,
                            neighbors,
                        },
                        i,
                        n,
                        j,
                    ));
                }
            }
        }

        // New cells: each boundary face joined to p, in the slot of the
        // vertex it replaces, which keeps the orientation.
        for &c in &scratch.region {
            self.alive[c as usize] = false;
            self.free.push(c);
        }
        for &(cell, _, outside, j) in &scratch.created {
            let id = self.allocate(cell);
            self.cells[outside as usize].neighbors[j] = id;
            scratch.ids.push(id);
        }

        // Faces between new cells contain p and one edge of the boundary;
        // each such edge is shared by exactly two new cells. Pair them in an
        // open-addressing table keyed by the edge.
        let capacity = (6 * scratch.created.len()).next_power_of_two();
        scratch.edges.clear();
        scratch.edges.resize(capacity, (u64::MAX, 0));
        for (k, &(cell, i, _, _)) in scratch.created.iter().enumerate() {
            for j in (0..4).filter(|&j| j != i) {
                let mut ends = (0..4)
                    .filter(|&slot| slot != i && slot != j)
                    .map(|slot| cell.vertices[slot]);
                let (x, y) = (
                    ends.next().unwrap_or(INFINITE),
                    ends.next().unwrap_or(INFINITE),
                );
                let edge = (u64::from(x.min(y)) << 32) | u64::from(x.max(y));
                let here = (scratch.ids[k], j as u32);
                // Fibonacci hashing; the table is at least 1.5x the entries.
                let mut probe = (edge.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 32) as usize;
                loop {
                    probe &= capacity - 1;
                    let (stored, other) = scratch.edges[probe];
                    if stored == u64::MAX {
                        scratch.edges[probe] = (edge, (u64::from(here.0) << 2) | u64::from(here.1));
                        break;
                    }
                    if stored == edge {
                        let (b, l) = ((other >> 2) as u32, (other & 3) as usize);
                        self.cells[here.0 as usize].neighbors[j] = b;
                        self.cells[b as usize].neighbors[l] = here.0;
                        // Each edge pairs exactly once; free the slot's key
                        // so a stale match is impossible.
                        scratch.edges[probe].0 = u64::MAX - 1;
                        break;
                    }
                    probe += 1;
                }
            }
        }
        self.hint = scratch.ids[0];
        self.scratch = scratch;
        Insertion::Vertex(p as usize)
    }
}

/// The vertices of cell face `i` (opposite slot `i`), sorted.
fn face_key(vertices: [u32; 4], i: usize) -> [u32; 3] {
    let mut face = [0; 3];
    let mut k = 0;
    for (slot, &v) in vertices.iter().enumerate() {
        if slot != i {
            face[k] = v;
            k += 1;
        }
    }
    face.sort_unstable();
    face
}
