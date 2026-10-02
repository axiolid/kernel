//! Subdivision surfaces: Loop on triangle meshes, Catmull-Clark on
//! polygon meshes.
//!
//! Unlike [`crate::refine`], which splits triangles and leaves the surface
//! where it was, a subdivision *scheme* also moves every vertex by a fixed
//! stencil (a mask) of its neighbours. Repeating it converges to a smooth
//! limit surface: C2 except at extraordinary vertices, where it is C1.
//! Each level is a new [`HalfedgeMesh`]; the input is not modified.
//!
//! # Masks
//!
//! **Loop** (Loop 1987), triangles only, 1-to-4 split:
//!
//! - new edge vertex: `3/8 (a + b) + 1/8 (c + d)`, `c`, `d` the corners
//!   opposite the edge;
//! - old vertex of valence `n`: `(1 - n beta) v + beta sum(neighbours)` with
//!   Loop's `beta = (5/8 - (3/8 + 1/4 cos(2 pi / n))^2) / n`.
//!
//! **Catmull-Clark** (Catmull and Clark 1978), any polygons, each `k`-gon
//! split into `k` quads:
//!
//! - face vertex: the centroid of the corners;
//! - new edge vertex: `(a + b + f1 + f2) / 4` with the two face vertices;
//! - old vertex of valence `n`: `(Q + 2 R + (n - 3) v) / n`, `Q` the mean of
//!   the adjacent face vertices and `R` the mean of the incident edges'
//!   midpoints.
//!
//! **Boundaries** (both schemes) follow the cubic B-spline curve, so a
//! boundary converges to a curve that depends on boundary vertices only and
//! two meshes sharing a boundary stay joined: a boundary edge gets its
//! midpoint and a boundary vertex `3/4 v + 1/8 (b1 + b2)` from its two
//! boundary neighbours. Isolated vertices are carried unchanged.
//!
//! # Numbering
//!
//! Each level's vertices are, in order: the input's live vertices (so a
//! compact input keeps its vertex ids), one per live edge in edge-id order,
//! and for Catmull-Clark one per live face in face-id order. Faces follow
//! the input faces in face-id order. The output is deterministic.
//!
//! [`limit_positions`] evaluates where each vertex converges, from the
//! eigen-analysis of the masks above (Halstead, Kass and DeRose 1993).

use core::f64::consts::PI;

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::halfedge::{FaceId, HalfedgeBuildError, HalfedgeId, HalfedgeMesh, VertexId};

/// A subdivision scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SubdivisionScheme {
    /// Loop: triangle meshes, 1-to-4 split, C2 away from extraordinary
    /// vertices.
    Loop,
    /// Catmull-Clark: polygon meshes, all quads after the first level, C2
    /// away from extraordinary vertices.
    CatmullClark,
}

/// Options of [`subdivide`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubdivisionOptions {
    /// The scheme.
    pub scheme: SubdivisionScheme,
    /// Number of levels; zero returns a copy.
    pub levels: u32,
    /// Largest face count any level may produce.
    pub max_faces: usize,
}

impl SubdivisionOptions {
    /// `levels` levels of `scheme` with the default face budget.
    #[must_use]
    pub const fn new(scheme: SubdivisionScheme, levels: u32) -> Self {
        Self {
            scheme,
            levels,
            max_faces: 50_000_000,
        }
    }
}

/// Why subdivision refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum SubdivisionError {
    /// Loop subdivision needs triangles.
    #[error("face {face} has {degree} sides; Loop subdivision needs triangles")]
    NotATriangle {
        /// The face.
        face: FaceId,
        /// Its number of sides.
        degree: usize,
    },
    /// Catmull-Clark limit positions need quads (subdivide once first).
    #[error("face {face} has {degree} sides; Catmull-Clark limit positions need quads")]
    NotAQuad {
        /// The face.
        face: FaceId,
        /// Its number of sides.
        degree: usize,
    },
    /// A level would produce more faces than allowed.
    #[error("subdivision would produce {faces} faces, over the {limit} budget")]
    BudgetExceeded {
        /// Faces the level would have.
        faces: usize,
        /// The budget in force.
        limit: usize,
    },
    /// Building a level's mesh failed. Not expected: every level of a valid
    /// mesh is valid. Reported rather than unwrapped.
    #[error("subdivided mesh could not be built: {0}")]
    Build(#[from] HalfedgeBuildError),
}

/// Subdivide `mesh` `options.levels` times.
///
/// # Errors
///
/// [`SubdivisionError::NotATriangle`] for Loop on a mesh with another face
/// and [`SubdivisionError::BudgetExceeded`] before building a level over
/// the face budget.
pub fn subdivide(
    mesh: &HalfedgeMesh,
    options: &SubdivisionOptions,
) -> Result<HalfedgeMesh, SubdivisionError> {
    let mut current = mesh.clone();
    for _ in 0..options.levels {
        let faces = match options.scheme {
            SubdivisionScheme::Loop => {
                require_triangles(&current)?;
                current.face_count() * 4
            }
            SubdivisionScheme::CatmullClark => current
                .faces()
                .map(|f| current.face_degree(f))
                .sum::<usize>(),
        };
        if faces > options.max_faces {
            return Err(SubdivisionError::BudgetExceeded {
                faces,
                limit: options.max_faces,
            });
        }
        current = match options.scheme {
            SubdivisionScheme::Loop => loop_level(&current)?,
            SubdivisionScheme::CatmullClark => catmull_clark_level(&current)?,
        };
    }
    Ok(current)
}

/// The limit position of every vertex under repeated subdivision by
/// `scheme`, indexed by vertex id (a removed vertex keeps its stored
/// position).
///
/// # Errors
///
/// [`SubdivisionError::NotATriangle`] for Loop and
/// [`SubdivisionError::NotAQuad`] for Catmull-Clark on a face of another
/// degree.
pub fn limit_positions(
    mesh: &HalfedgeMesh,
    scheme: SubdivisionScheme,
) -> Result<Vec<Point3>, SubdivisionError> {
    match scheme {
        SubdivisionScheme::Loop => require_triangles(mesh)?,
        SubdivisionScheme::CatmullClark => {
            for f in mesh.faces() {
                let degree = mesh.face_degree(f);
                if degree != 4 {
                    return Err(SubdivisionError::NotAQuad { face: f, degree });
                }
            }
        }
    }
    let mut out = mesh.positions().to_vec();
    for v in mesh.vertices() {
        let p = mesh.position(v);
        out[v.index()] = if mesh.is_isolated(v) {
            p
        } else if let Some((b1, b2)) = boundary_neighbours(mesh, v) {
            // Limit of the cubic B-spline mask (1, 6, 1) / 8.
            (4.0 * p + b1 + b2) / 6.0
        } else {
            match scheme {
                SubdivisionScheme::Loop => {
                    let n = mesh.degree(v) as Scalar;
                    let gamma = 1.0 / (3.0 / (8.0 * loop_beta(n)) + n);
                    let sum = mesh
                        .vertex_vertices(v)
                        .fold(Vec3::ZERO, |s, u| s + mesh.position(u));
                    (1.0 - n * gamma) * p + gamma * sum
                }
                SubdivisionScheme::CatmullClark => {
                    let n = mesh.degree(v) as Scalar;
                    let mut edges = Vec3::ZERO;
                    let mut diagonals = Vec3::ZERO;
                    for h in mesh.outgoing_halfedges(v) {
                        edges += mesh.position(mesh.target(h));
                        diagonals += mesh.position(mesh.target(mesh.next(h)));
                    }
                    (n * n * p + 4.0 * edges + diagonals) / (n * (n + 5.0))
                }
            }
        };
    }
    Ok(out)
}

fn require_triangles(mesh: &HalfedgeMesh) -> Result<(), SubdivisionError> {
    for f in mesh.faces() {
        let degree = mesh.face_degree(f);
        if degree != 3 {
            return Err(SubdivisionError::NotATriangle { face: f, degree });
        }
    }
    Ok(())
}

/// Loop's `beta` for valence `n`.
fn loop_beta(n: Scalar) -> Scalar {
    let c = 3.0 / 8.0 + 0.25 * (2.0 * PI / n).cos();
    (5.0 / 8.0 - c * c) / n
}

/// The two boundary neighbours of a boundary vertex (next along its
/// boundary halfedge, previous along the incoming one); `None` for an
/// interior vertex.
fn boundary_neighbours(mesh: &HalfedgeMesh, v: VertexId) -> Option<(Point3, Point3)> {
    let h = mesh.vertex_halfedge(v)?;
    if !mesh.is_boundary_halfedge(h) {
        return None;
    }
    let next = mesh.target(h);
    let previous = mesh.source(mesh.prev(h));
    Some((mesh.position(next), mesh.position(previous)))
}

/// Cubic B-spline vertex mask, shared by both schemes on a boundary.
fn boundary_vertex(p: Point3, b1: Point3, b2: Point3) -> Point3 {
    0.75 * p + 0.125 * (b1 + b2)
}

/// Dense numbering of the live input vertices.
fn vertex_numbering(mesh: &HalfedgeMesh) -> (Vec<u32>, Vec<VertexId>) {
    let mut number = vec![u32::MAX; mesh.positions().len()];
    let live: Vec<VertexId> = mesh.vertices().collect();
    for (k, v) in live.iter().enumerate() {
        number[v.index()] = k as u32;
    }
    (number, live)
}

/// Dense numbering of the live edges, offset by `base`.
fn edge_numbering(mesh: &HalfedgeMesh, base: usize) -> Vec<u32> {
    let slots = mesh.edges().last().map_or(0, |e| e.index() + 1);
    let mut number = vec![u32::MAX; slots];
    for (k, e) in mesh.edges().enumerate() {
        number[e.index()] = (base + k) as u32;
    }
    number
}

fn loop_level(mesh: &HalfedgeMesh) -> Result<HalfedgeMesh, SubdivisionError> {
    let (vertex_number, live) = vertex_numbering(mesh);
    let edge_number = edge_numbering(mesh, live.len());
    let mut positions = Vec::with_capacity(live.len() + mesh.edge_count());
    for &v in &live {
        let p = mesh.position(v);
        positions.push(if mesh.is_isolated(v) {
            p
        } else if let Some((b1, b2)) = boundary_neighbours(mesh, v) {
            boundary_vertex(p, b1, b2)
        } else {
            let n = mesh.degree(v) as Scalar;
            let beta = loop_beta(n);
            let sum = mesh
                .vertex_vertices(v)
                .fold(Vec3::ZERO, |s, u| s + mesh.position(u));
            (1.0 - n * beta) * p + beta * sum
        });
    }
    for e in mesh.edges() {
        let h = mesh.edge_halfedge(e, 0);
        let g = mesh.opposite(h);
        let a = mesh.position(mesh.source(h));
        let b = mesh.position(mesh.target(h));
        positions.push(if mesh.is_boundary_edge(e) {
            0.5 * (a + b)
        } else {
            let c = mesh.position(mesh.target(mesh.next(h)));
            let d = mesh.position(mesh.target(mesh.next(g)));
            0.375 * (a + b) + 0.125 * (c + d)
        });
    }
    let mut faces: Vec<[u32; 3]> = Vec::with_capacity(mesh.face_count() * 4);
    for f in mesh.faces() {
        let hs: Vec<HalfedgeId> = mesh.face_halfedges(f).collect();
        let corner = |i: usize| vertex_number[mesh.source(hs[i]).index()];
        let mid = |i: usize| edge_number[mesh.edge(hs[i]).index()];
        for i in 0..3 {
            faces.push([corner(i), mid(i), mid((i + 2) % 3)]);
        }
        faces.push([mid(0), mid(1), mid(2)]);
    }
    Ok(HalfedgeMesh::from_faces(positions, &faces)?)
}

fn catmull_clark_level(mesh: &HalfedgeMesh) -> Result<HalfedgeMesh, SubdivisionError> {
    let (vertex_number, live) = vertex_numbering(mesh);
    let edge_number = edge_numbering(mesh, live.len());
    let face_base = live.len() + mesh.edge_count();
    let mut face_number = vec![u32::MAX; mesh.faces().last().map_or(0, |f| f.index() + 1)];
    let mut face_point = vec![Vec3::ZERO; face_number.len()];
    for (k, f) in mesh.faces().enumerate() {
        face_number[f.index()] = (face_base + k) as u32;
        let (sum, count) = mesh
            .face_vertices(f)
            .fold((Vec3::ZERO, 0usize), |(s, c), v| {
                (s + mesh.position(v), c + 1)
            });
        face_point[f.index()] = sum / count as Scalar;
    }
    let mut positions = Vec::with_capacity(face_base + mesh.face_count());
    for &v in &live {
        let p = mesh.position(v);
        positions.push(if mesh.is_isolated(v) {
            p
        } else if let Some((b1, b2)) = boundary_neighbours(mesh, v) {
            boundary_vertex(p, b1, b2)
        } else {
            let mut q = Vec3::ZERO;
            let mut r = Vec3::ZERO;
            let mut n = 0usize;
            for h in mesh.outgoing_halfedges(v) {
                let f = mesh
                    .face(h)
                    .expect("an interior vertex has faces all round");
                q += face_point[f.index()];
                r += 0.5 * (p + mesh.position(mesh.target(h)));
                n += 1;
            }
            let n = n as Scalar;
            (q / n + 2.0 * r / n + (n - 3.0) * p) / n
        });
    }
    for e in mesh.edges() {
        let h = mesh.edge_halfedge(e, 0);
        let a = mesh.position(mesh.source(h));
        let b = mesh.position(mesh.target(h));
        positions.push(match (mesh.face(h), mesh.face(mesh.opposite(h))) {
            (Some(f1), Some(f2)) => {
                0.25 * (a + b + face_point[f1.index()] + face_point[f2.index()])
            }
            _ => 0.5 * (a + b),
        });
    }
    for f in mesh.faces() {
        positions.push(face_point[f.index()]);
    }
    let mut faces: Vec<[u32; 4]> =
        Vec::with_capacity(mesh.faces().map(|f| mesh.face_degree(f)).sum::<usize>());
    for f in mesh.faces() {
        let hs: Vec<HalfedgeId> = mesh.face_halfedges(f).collect();
        let k = hs.len();
        let centre = face_number[f.index()];
        for i in 0..k {
            let previous = hs[(i + k - 1) % k];
            faces.push([
                vertex_number[mesh.source(hs[i]).index()],
                edge_number[mesh.edge(hs[i]).index()],
                centre,
                edge_number[mesh.edge(previous).index()],
            ]);
        }
    }
    Ok(HalfedgeMesh::from_faces(positions, &faces)?)
}
