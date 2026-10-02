//! Fairing: move chosen vertices to the smoothest surface their fixed
//! surroundings allow.
//!
//! Smoothing ([`crate::smooth`]) relaxes vertices a step at a time and never
//! arrives anywhere in particular. Fairing *solves* for the arrival: the
//! free vertices are placed where the discrete Laplacian (harmonic, a
//! membrane) or bi-Laplacian (biharmonic, a thin plate) of the surface
//! vanishes, with every other vertex held fixed. The biharmonic solution
//! also matches the fixed surface's first derivatives across the border,
//! which is what makes a filled hole continue its surroundings tangentially
//! instead of meeting them at a crease. This is the fairing step of Liepa
//! (2003) and of CGAL's `Polygon_mesh_processing::fair`.
//!
//! # The system
//!
//! With symmetric edge weights `w_ij` the Laplacian `L` (`L_ij = w_ij`,
//! `L_ii = -sum_j w_ij`) is symmetric, and with a positive diagonal mass `M`
//! the biharmonic conditions `(L M^-1 L x)_i = 0` for the free vertices `F`
//! read
//!
//! ```text
//! (L_{:,F})^T M^-1 L_{:,F} x_F = -(L_{:,F})^T M^-1 L_{:,C} x_C
//! ```
//!
//! whose matrix is a Gram matrix: symmetric positive definite whenever each
//! connected group of free vertices touches a fixed one. It is solved for
//! each coordinate by [`axiolid_numeric::conjugate_gradient`]; the harmonic
//! case solves `-L_FF x_F = L_FC x_C` the same way. Only rows of `L` at the
//! free vertices and their neighbours enter, so the cost follows the faired
//! region, not the mesh.
//!
//! Weights are computed once, on the positions before fairing.
//! [`FairingWeights::Cotangent`] uses `(cot a + cot b) / 2` with the
//! barycentric (one-third) vertex area as mass, the standard linear finite
//! element discretisation; [`FairingWeights::Uniform`] uses unit weights and
//! the valence as mass (Kobbelt's umbrella operator), so it ignores the
//! shape of the triangles and also works on polygonal faces.
//!
//! # Exactness
//!
//! Rows of `L` sum to zero, so every affine function of position is
//! reproduced: fixed vertices that all lie in one plane keep the free ones in
//! it, up to the solver tolerance, and in an axis-aligned plane exactly.

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::halfedge::{FaceId, HalfedgeMesh, VertexId};
use axiolid_numeric::{conjugate_gradient, ConjugateGradientOptions, NumericError, SparseMatrix};

/// Edge weights of the discrete Laplacian.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum FairingWeights {
    /// Cotangent weights with barycentric vertex areas as mass: follows the
    /// geometry, needs triangles with non-zero area around the faired
    /// region.
    #[default]
    Cotangent,
    /// Unit weights, valence as mass: connectivity only, any polygon.
    Uniform,
}

/// Which derivative vanishes on the faired region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum FairingOrder {
    /// `L x = 0`: a membrane, continuous with the border (C0).
    Harmonic,
    /// `L M^-1 L x = 0`: a thin plate whose tangent planes continue the
    /// border's (G1).
    #[default]
    Biharmonic,
}

/// Options of [`fair`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FairOptions {
    /// Edge weights.
    pub weights: FairingWeights,
    /// Harmonic or biharmonic.
    pub order: FairingOrder,
    /// Conjugate-gradient stopping rule: relative residual of each
    /// coordinate's system.
    pub relative_tolerance: Scalar,
    /// Conjugate-gradient iteration budget per coordinate; `None` uses the
    /// solver's default (`10 n + 100` for `n` free vertices).
    pub max_iterations: Option<usize>,
}

impl Default for FairOptions {
    fn default() -> Self {
        Self {
            weights: FairingWeights::Cotangent,
            order: FairingOrder::Biharmonic,
            relative_tolerance: 1e-12,
            max_iterations: None,
        }
    }
}

/// What [`fair`] did.
#[derive(Debug, Clone, PartialEq)]
pub struct FairReport {
    /// Vertices moved (distinct).
    pub vertices: usize,
    /// Conjugate-gradient iterations for `x`, `y` and `z`.
    pub iterations: [usize; 3],
    /// Largest relative residual of the three solves, recomputed from the
    /// system matrix.
    pub relative_residual: Scalar,
}

/// Why [`fair`] refused. On every refusal the mesh is unchanged.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum FairError {
    /// The vertex id names no live vertex.
    #[error("vertex {vertex} is not a live vertex")]
    RemovedVertex {
        /// The vertex.
        vertex: VertexId,
    },
    /// The vertex is free and so is every vertex connected to it through
    /// free vertices: nothing fixes where that group goes.
    #[error("vertex {vertex} is not connected to any fixed vertex")]
    Unanchored {
        /// A vertex of the unanchored group.
        vertex: VertexId,
    },
    /// Cotangent weights need triangles here.
    #[error("face {face} has {degree} sides; cotangent weights need triangles")]
    NotATriangle {
        /// The face.
        face: FaceId,
        /// Its number of sides.
        degree: usize,
    },
    /// A triangle near the faired region has (numerically) zero area, so its
    /// cotangent weights are undefined.
    #[error("face {face} is degenerate; its cotangent weights are undefined")]
    DegenerateTriangle {
        /// The face.
        face: FaceId,
    },
    /// The relative tolerance is negative or not finite.
    #[error("relative tolerance {0} is not a finite non-negative number")]
    InvalidTolerance(Scalar),
    /// A coordinate's system did not reach the tolerance within the budget.
    #[error("fairing solve for coordinate {coordinate} stopped at relative residual {relative_residual}")]
    NotConverged {
        /// 0, 1 or 2 for `x`, `y`, `z`.
        coordinate: usize,
        /// Residual reached.
        relative_residual: Scalar,
    },
    /// The linear solve refused its system, which happens only with
    /// cotangent weights on badly shaped triangles that make it indefinite
    /// to working precision.
    #[error("fairing solve refused: {0}")]
    Numeric(NumericError),
}

/// Smallest `|sin|` of a corner angle accepted for a cotangent.
const MIN_SINE: Scalar = 1e-12;

/// Fair the given vertices of a mesh in place; every other vertex is fixed.
///
/// Duplicates in `vertices` are ignored. Fixed vertices are bit-identical
/// afterwards. The output is deterministic: it depends on the mesh and the
/// set of vertices, not on their order.
///
/// # Errors
///
/// [`FairError`]: a removed vertex, a free group with no fixed neighbour,
/// a non-triangle or degenerate triangle near the region with cotangent
/// weights, a bad tolerance, and a solve that does not converge.
pub fn fair(
    mesh: &mut HalfedgeMesh,
    vertices: &[VertexId],
    options: &FairOptions,
) -> Result<FairReport, FairError> {
    if !options.relative_tolerance.is_finite() || options.relative_tolerance < 0.0 {
        return Err(FairError::InvalidTolerance(options.relative_tolerance));
    }
    let mut free: Vec<VertexId> = vertices.to_vec();
    free.sort_unstable();
    free.dedup();
    for &v in &free {
        if !mesh.contains_vertex(v) {
            return Err(FairError::RemovedVertex { vertex: v });
        }
    }
    if free.is_empty() {
        return Ok(FairReport {
            vertices: 0,
            iterations: [0; 3],
            relative_residual: 0.0,
        });
    }
    let slot_count = mesh.positions().len();
    let mut unknown = vec![usize::MAX; slot_count];
    for (k, &v) in free.iter().enumerate() {
        unknown[v.index()] = k;
    }
    require_anchored(mesh, &free, &unknown)?;

    // Rows of L that touch a free vertex: the free vertices and their
    // neighbours, in id order.
    let mut in_rows = vec![false; slot_count];
    for &v in &free {
        in_rows[v.index()] = true;
        for u in mesh.vertex_vertices(v) {
            in_rows[u.index()] = true;
        }
    }
    let rows: Vec<VertexId> = (0..slot_count)
        .filter(|&i| in_rows[i])
        .map(|i| VertexId::new(i as u32))
        .collect();

    // Translate to the centroid of the fixed vertices next to the region so
    // the solve sees coordinates of the region's size, not of its distance
    // from the origin. Fixed vertices only: a coordinate that is zero on all
    // of them then stays exactly zero.
    let anchors: Vec<VertexId> = rows
        .iter()
        .copied()
        .filter(|v| unknown[v.index()] == usize::MAX)
        .collect();
    let origin = anchors
        .iter()
        .fold(Vec3::ZERO, |s, &v| s + mesh.position(v))
        / anchors.len() as Scalar;
    let local = |v: VertexId| mesh.position(v) - origin;

    let n = free.len();
    let mut triplets: Vec<(usize, usize, Scalar)> = Vec::new();
    let mut rhs = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    match options.order {
        FairingOrder::Harmonic => {
            for (a, &v) in free.iter().enumerate() {
                let (diagonal, neighbours) = laplacian_row(mesh, v, options.weights)?;
                triplets.push((a, a, -diagonal));
                for (u, w) in neighbours {
                    let b = unknown[u.index()];
                    if b == usize::MAX {
                        let p = local(u);
                        for (axis, r) in rhs.iter_mut().enumerate() {
                            r[a] += w * p[axis];
                        }
                    } else {
                        triplets.push((a, b, -w));
                    }
                }
            }
        }
        FairingOrder::Biharmonic => {
            for &k in &rows {
                let (diagonal, neighbours) = laplacian_row(mesh, k, options.weights)?;
                let mass = vertex_mass(mesh, k, options.weights)?;
                let mut free_part: Vec<(usize, Scalar)> = Vec::new();
                let mut fixed_part = Vec3::ZERO;
                for (u, w) in core::iter::once((k, diagonal)).chain(neighbours) {
                    let b = unknown[u.index()];
                    if b == usize::MAX {
                        fixed_part += w * local(u);
                    } else {
                        free_part.push((b, w));
                    }
                }
                for &(a, wa) in &free_part {
                    for &(b, wb) in &free_part {
                        triplets.push((a, b, wa * wb / mass));
                    }
                    for (axis, r) in rhs.iter_mut().enumerate() {
                        r[a] -= wa * fixed_part[axis] / mass;
                    }
                }
            }
        }
    }
    let matrix = SparseMatrix::from_triplets(n, n, &triplets).map_err(FairError::Numeric)?;
    let cg = ConjugateGradientOptions {
        relative_tolerance: options.relative_tolerance,
        max_iterations: options.max_iterations,
    };
    let mut solved = [vec![], vec![], vec![]];
    let mut iterations = [0usize; 3];
    let mut relative_residual: Scalar = 0.0;
    for axis in 0..3 {
        let initial: Vec<Scalar> = free.iter().map(|&v| local(v)[axis]).collect();
        let solution = conjugate_gradient(&matrix, &rhs[axis], Some(&initial), cg)
            .map_err(FairError::Numeric)?;
        if solution.status != axiolid_numeric::Status::Converged {
            return Err(FairError::NotConverged {
                coordinate: axis,
                relative_residual: solution.relative_residual,
            });
        }
        iterations[axis] = solution.iterations;
        relative_residual = relative_residual.max(solution.relative_residual);
        solved[axis] = solution.x;
    }
    for (k, &v) in free.iter().enumerate() {
        let p = Point3::new(solved[0][k], solved[1][k], solved[2][k]) + origin;
        mesh.set_position(v, p);
    }
    Ok(FairReport {
        vertices: n,
        iterations,
        relative_residual,
    })
}

/// Refuse a group of free vertices that no fixed vertex touches.
fn require_anchored(
    mesh: &HalfedgeMesh,
    free: &[VertexId],
    unknown: &[usize],
) -> Result<(), FairError> {
    let mut seen = vec![false; free.len()];
    for start in 0..free.len() {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut stack = vec![free[start]];
        let mut anchored = false;
        while let Some(v) = stack.pop() {
            for u in mesh.vertex_vertices(v) {
                let k = unknown[u.index()];
                if k == usize::MAX {
                    anchored = true;
                } else if !seen[k] {
                    seen[k] = true;
                    stack.push(u);
                }
            }
        }
        if !anchored {
            return Err(FairError::Unanchored {
                vertex: free[start],
            });
        }
    }
    Ok(())
}

/// Row `v` of the Laplacian: the diagonal `-sum w` and each neighbour's
/// weight, neighbours in circulation order.
fn laplacian_row(
    mesh: &HalfedgeMesh,
    v: VertexId,
    weights: FairingWeights,
) -> Result<(Scalar, Vec<(VertexId, Scalar)>), FairError> {
    let mut row = Vec::new();
    let mut sum = 0.0;
    for h in mesh.outgoing_halfedges(v) {
        let w = match weights {
            FairingWeights::Uniform => 1.0,
            FairingWeights::Cotangent => {
                let mut w = 0.0;
                for g in [h, mesh.opposite(h)] {
                    if let Some(f) = mesh.face(g) {
                        let opposite = triangle_apex(mesh, g, f)?;
                        w += 0.5 * cotangent(mesh, f, opposite, mesh.source(g), mesh.target(g))?;
                    }
                }
                w
            }
        };
        sum += w;
        row.push((mesh.target(h), w));
    }
    Ok((-sum, row))
}

/// The corner of triangle `f` opposite its halfedge `g`.
fn triangle_apex(
    mesh: &HalfedgeMesh,
    g: axiolid_mesh::halfedge::HalfedgeId,
    f: FaceId,
) -> Result<VertexId, FairError> {
    let degree = mesh.face_degree(f);
    if degree != 3 {
        return Err(FairError::NotATriangle { face: f, degree });
    }
    Ok(mesh.target(mesh.next(g)))
}

/// Cotangent of the angle at `apex` in the triangle `apex, a, b` of face `f`.
fn cotangent(
    mesh: &HalfedgeMesh,
    f: FaceId,
    apex: VertexId,
    a: VertexId,
    b: VertexId,
) -> Result<Scalar, FairError> {
    let p = mesh.position(apex);
    let u = mesh.position(a) - p;
    let v = mesh.position(b) - p;
    let cross = u.cross(v).length();
    if cross.is_nan() || cross <= MIN_SINE * u.length() * v.length() {
        return Err(FairError::DegenerateTriangle { face: f });
    }
    Ok(u.dot(v) / cross)
}

/// Mass of a vertex: its valence for uniform weights, a third of the incident
/// triangles' area for cotangent weights.
fn vertex_mass(
    mesh: &HalfedgeMesh,
    v: VertexId,
    weights: FairingWeights,
) -> Result<Scalar, FairError> {
    match weights {
        FairingWeights::Uniform => Ok(mesh.degree(v) as Scalar),
        FairingWeights::Cotangent => {
            let mut area = 0.0;
            for f in mesh.vertex_faces(v) {
                let corners: Vec<Point3> =
                    mesh.face_vertices(f).map(|c| mesh.position(c)).collect();
                if corners.len() != 3 {
                    return Err(FairError::NotATriangle {
                        face: f,
                        degree: corners.len(),
                    });
                }
                area += 0.5
                    * (corners[1] - corners[0])
                        .cross(corners[2] - corners[0])
                        .length();
            }
            // Every incident triangle passed the cotangent test, so the sum
            // is positive unless the vertex has no face at all.
            if area > 0.0 {
                Ok(area / 3.0)
            } else {
                Err(FairError::Unanchored { vertex: v })
            }
        }
    }
}
