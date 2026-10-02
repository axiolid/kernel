//! Isotropic remeshing: rewrite a triangle mesh so its edges approach one
//! target length and its vertices approach valence six, while the surface
//! stays on the input.
//!
//! # Method
//!
//! The incremental scheme of Botsch and Kobbelt ("A remeshing approach to
//! multiresolution modeling", SGP 2004), as CGAL's
//! `Polygon_mesh_processing::isotropic_remeshing` implements it. With `L`
//! the target edge length, each iteration runs four passes over a
//! [`HalfedgeMesh`]:
//!
//! 1. split every edge longer than `4/3 L` at its midpoint;
//! 2. collapse every edge shorter than `4/5 L` whose collapse creates no
//!    edge longer than `4/3 L`;
//! 3. flip edges where that brings the four corners' valences closer to
//!    optimal (6 inside, 4 on a boundary);
//! 4. move every free vertex towards the area-weighted centroid of its
//!    triangles, within its tangent plane, and project it back onto the
//!    input surface.
//!
//! # What is kept
//!
//! *Protected* edges are the input's boundary edges and, unless
//! [`RemeshOptions::feature_angle`] is `None`, every edge whose two faces'
//! normals differ by more than that angle. A protected edge is never
//! flipped, may be split (its midpoint stays on the edge) and may be
//! collapsed only along itself, removing a vertex its feature line runs
//! straight through. So every protected polyline of the input is kept
//! exactly: its corners stay, bit-identical, and every vertex added on it
//! lies on one of its segments. Vertices on protected edges are not
//! relaxed. The boundary loops, the Euler characteristic and the
//! orientation of the input are therefore preserved.
//!
//! The protected edges cut the input into *patches*. A free vertex is
//! projected onto the closest point of its own patch only, so it never
//! crosses a sharp edge.
//!
//! # What is refused and what is skipped
//!
//! Input the halfedge structure cannot represent -- non-manifold edges or
//! vertices, inconsistent winding, out-of-range or repeated corners -- is
//! refused by name, as are non-finite positions, degenerate (zero-area)
//! input triangles, a non-positive target, a feature angle outside
//! `[0, pi]` and a request whose output would exceed the triangle budget.
//!
//! Individual edits are *skipped*, not forced: a collapse, flip or move that
//! would fold a triangle (turn its normal by a right angle or more) or leave
//! it flatter than a height-to-longest-edge ratio of `1e-3`, a collapse or
//! flip that would lower the smallest angle of the triangles it rewrites
//! below 15 degrees (when it was above), or an edit that breaks the link
//! condition, does not happen. The result is therefore always a
//! valid, consistently oriented mesh, and [`RemeshReport`] states how close
//! it came: the fraction of edges within `[4/5 L, 4/3 L]`, the valence
//! deviation before and after, the smallest angle, and the largest distance
//! of any output vertex from the input surface.
//!
//! # Deviation from the input
//!
//! Protected vertices lie on the input exactly (up to the rounding of a
//! midpoint). Free vertices are projected onto the input, so their distance
//! from it is zero up to rounding unless a projection was skipped to avoid a
//! fold; [`RemeshReport::max_vertex_deviation`] measures it over every
//! output vertex. Between vertices the output is a chord of the input: on a
//! surface of curvature radius `R` its distance grows like `L^2 / (8 R)`,
//! the sagitta of an edge of length `L`. That part is the cost of the
//! requested resolution, not of the algorithm.
//!
//! # Determinism
//!
//! Every pass visits elements in id order and every decision is a function
//! of the current mesh; no hash map is iterated. The same input and options
//! give a bit-identical output on every run.

mod edits;
mod reference;

use axiolid_core::{Point3, Scalar};
use axiolid_mesh::{
    AttributeFate, EdgeId, HalfedgeBuildError, HalfedgeEditError, HalfedgeMesh, TriMesh,
};

use edits::Remesher;
use reference::Reference;

/// Why a remeshing could not be performed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum RemeshError {
    /// The target edge length must be a positive, finite length.
    #[error("target edge length {0} is not a positive finite length")]
    InvalidTarget(Scalar),
    /// The feature angle must lie in `[0, pi]` radians.
    #[error("feature angle {0} is not within [0, pi] radians")]
    InvalidFeatureAngle(Scalar),
    /// The input is not an oriented 2-manifold triangle mesh: the halfedge
    /// structure names the configuration it cannot represent.
    #[error("input is not an oriented manifold triangle mesh: {0}")]
    InvalidMesh(#[from] HalfedgeBuildError),
    /// A vertex position is NaN or infinite.
    #[error("vertex {vertex} has a non-finite position")]
    NonFinitePosition {
        /// The vertex.
        vertex: usize,
    },
    /// An input triangle has no area, so it has no normal to keep.
    #[error("triangle {triangle} is degenerate")]
    DegenerateTriangle {
        /// The triangle.
        triangle: usize,
    },
    /// The output would exceed the triangle budget.
    #[error("remeshing would produce about {produced} triangles, over the {limit} budget")]
    BudgetExceeded {
        /// Triangles the request would produce.
        produced: usize,
        /// The cap.
        limit: usize,
    },
    /// A local edit refused a configuration the remesher had checked; this
    /// indicates a defect, not a property of the input.
    #[error("internal halfedge edit refused: {0}")]
    Internal(HalfedgeEditError),
}

/// What to remesh towards.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct RemeshOptions {
    /// Target edge length `L`, in model units.
    pub target_edge_length: Scalar,
    /// Split-collapse-flip-relax iterations.
    pub iterations: u32,
    /// Dihedral angle, in radians, above which an interior edge is a sharp
    /// feature and protected like a boundary edge. `None` protects only the
    /// boundary.
    pub feature_angle: Option<Scalar>,
}

impl RemeshOptions {
    /// Iterations run by default, as in Botsch and Kobbelt's paper.
    pub const DEFAULT_ITERATIONS: u32 = 5;

    /// Default feature angle: 60 degrees between face normals.
    pub const DEFAULT_FEATURE_ANGLE: Scalar = core::f64::consts::FRAC_PI_3;

    /// Remesh towards `target_edge_length` with the default iterations and
    /// feature angle.
    #[must_use]
    pub fn new(target_edge_length: Scalar) -> Self {
        Self {
            target_edge_length,
            iterations: Self::DEFAULT_ITERATIONS,
            feature_angle: Some(Self::DEFAULT_FEATURE_ANGLE),
        }
    }

    /// The same options with another iteration count.
    #[must_use]
    pub fn with_iterations(mut self, iterations: u32) -> Self {
        self.iterations = iterations;
        self
    }

    /// The same options with another feature angle; `None` protects only
    /// the boundary.
    #[must_use]
    pub fn with_feature_angle(mut self, feature_angle: Option<Scalar>) -> Self {
        self.feature_angle = feature_angle;
        self
    }
}

/// What a remeshing did and how close it came to the target.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct RemeshReport {
    /// Triangles before.
    pub input_triangles: usize,
    /// Triangles after.
    pub output_triangles: usize,
    /// Edge splits performed.
    pub splits: usize,
    /// Edge collapses performed.
    pub collapses: usize,
    /// Edge flips performed.
    pub flips: usize,
    /// Vertex relocations (relaxation and projection) performed.
    pub moves: usize,
    /// Protected (boundary or sharp feature) edges of the output.
    pub protected_edges: usize,
    /// Fraction of output edges with length in `[4/5 L, 4/3 L]`.
    pub edges_in_band: Scalar,
    /// Mean absolute distance of the input's vertex valences from optimal
    /// (6 inside, 4 on a boundary).
    pub valence_deviation_before: Scalar,
    /// The same for the output.
    pub valence_deviation_after: Scalar,
    /// Smallest interior angle of any output triangle, in radians.
    pub min_angle: Scalar,
    /// Largest distance of an output vertex from the input surface, in
    /// model units.
    pub max_vertex_deviation: Scalar,
    /// What happened to each named attribute channel: every one is dropped,
    /// since remeshing creates vertices with no preimage to read from.
    pub attribute_fates: Vec<(String, AttributeFate)>,
}

/// Remesh a triangle mesh towards a uniform edge length.
///
/// See the [module documentation](self) for the method, what is preserved
/// and what the report measures. Normals and attribute channels are not
/// carried; the report lists each channel as dropped.
///
/// # Errors
///
/// [`RemeshError::InvalidTarget`], [`RemeshError::InvalidFeatureAngle`],
/// [`RemeshError::InvalidMesh`] for non-manifold, inconsistently wound,
/// ragged or out-of-range input, [`RemeshError::NonFinitePosition`],
/// [`RemeshError::DegenerateTriangle`] and [`RemeshError::BudgetExceeded`].
pub fn remesh(
    mesh: &TriMesh,
    options: RemeshOptions,
) -> Result<(TriMesh, RemeshReport), RemeshError> {
    let target = options.target_edge_length;
    if !target.is_finite() || target <= 0.0 {
        return Err(RemeshError::InvalidTarget(target));
    }
    if let Some(angle) = options.feature_angle {
        if !(0.0..=core::f64::consts::PI).contains(&angle) {
            return Err(RemeshError::InvalidFeatureAngle(angle));
        }
    }
    let input = HalfedgeMesh::from_tri_mesh(mesh)?;
    if let Some(vertex) = mesh.positions.iter().position(|p| !p.is_finite()) {
        return Err(RemeshError::NonFinitePosition { vertex });
    }
    let triangles: Vec<[Point3; 3]> = mesh
        .triangles()
        .map(|t| t.map(|i| mesh.positions[i as usize]))
        .collect();
    let mut area = 0.0;
    for (index, t) in triangles.iter().enumerate() {
        let n = edits::normal(*t);
        let scale = (t[1] - t[0]).length() * (t[2] - t[0]).length();
        if n.length() <= DEGENERATE_SINE * scale {
            return Err(RemeshError::DegenerateTriangle { triangle: index });
        }
        area += n.length() / 2.0;
    }
    // An equilateral triangle of side L has area sqrt(3)/4 L^2.
    let expected = area / (3.0_f64.sqrt() / 4.0 * target * target);
    if !expected.is_finite() || expected > crate::MAX_TRIANGLES as Scalar {
        return Err(RemeshError::BudgetExceeded {
            produced: if expected.is_finite() {
                expected as usize
            } else {
                usize::MAX
            },
            limit: crate::MAX_TRIANGLES,
        });
    }

    let protected = protected_edges(&input, &triangles, options.feature_angle);
    let patch = patches(&input, &protected);
    let reference = Reference::new(triangles, &patch);
    let valence_deviation_before = valence_deviation(&input);
    let input_triangles = input.face_count();

    let mut remesher = Remesher::new(
        input,
        protected,
        patch,
        &reference,
        target,
        crate::MAX_TRIANGLES,
    );
    for _ in 0..options.iterations {
        remesher.split_long_edges()?;
        remesher.collapse_short_edges();
        remesher.equalize_valences();
        remesher.relax_and_project();
    }

    let result = &remesher.mesh;
    let (low, high) = (target * 4.0 / 5.0, target * 4.0 / 3.0);
    let mut in_band = 0usize;
    let mut protected_count = 0usize;
    for e in result.edges() {
        let [a, b] = result.edge_vertices(e);
        let length = (result.position(a) - result.position(b)).length();
        if (low..=high).contains(&length) {
            in_band += 1;
        }
        if remesher.is_protected(e) {
            protected_count += 1;
        }
    }
    let edges_in_band = if result.edge_count() == 0 {
        1.0
    } else {
        in_band as Scalar / result.edge_count() as Scalar
    };
    let max_vertex_deviation = result
        .vertices()
        .filter(|&v| !result.is_isolated(v))
        .map(|v| reference.distance(result.position(v)))
        .fold(0.0, Scalar::max);
    let min_angle = result
        .faces()
        .map(|f| {
            let mut corners = result.face_vertices(f).map(|v| result.position(v));
            let t = [(); 3].map(|()| corners.next().unwrap_or(Point3::ZERO));
            edits::min_angle(t)
        })
        .fold(core::f64::consts::PI, Scalar::min);
    let valence_deviation_after = valence_deviation(result);
    let out = result.to_tri_mesh().map_err(RemeshError::Internal)?;
    let counts = remesher.counts;

    let attribute_fates = crate::carry_attributes(mesh);
    let report = RemeshReport {
        input_triangles,
        output_triangles: out.triangle_count(),
        splits: counts.splits,
        collapses: counts.collapses,
        flips: counts.flips,
        moves: counts.moves,
        protected_edges: protected_count,
        edges_in_band,
        valence_deviation_before,
        valence_deviation_after,
        min_angle,
        max_vertex_deviation,
        attribute_fates,
    };
    Ok((out, report))
}

/// An input triangle whose corner angle at its first vertex has a sine at
/// most this is refused as degenerate.
const DEGENERATE_SINE: Scalar = 1e-12;

/// Boundary edges, and interior edges whose faces' normals differ by more
/// than the feature angle. Indexed by edge id.
fn protected_edges(
    mesh: &HalfedgeMesh,
    triangles: &[[Point3; 3]],
    feature_angle: Option<Scalar>,
) -> Vec<bool> {
    let mut protected = vec![false; mesh.edges().last().map_or(0, |e| e.index() + 1)];
    for e in mesh.edges() {
        protected[e.index()] = if mesh.is_boundary_edge(e) {
            true
        } else if let Some(angle) = feature_angle {
            let normal = |side| {
                let f = mesh.face(mesh.edge_halfedge(e, side)).map(|f| f.index());
                f.and_then(|f| edits::normal(triangles[f]).try_normalize())
            };
            match (normal(0), normal(1)) {
                (Some(n0), Some(n1)) => n0.dot(n1).clamp(-1.0, 1.0).acos() > angle,
                _ => true,
            }
        } else {
            false
        };
    }
    protected
}

/// Connected components of faces across unprotected edges, numbered in
/// order of their lowest face id. Indexed by face id.
fn patches(mesh: &HalfedgeMesh, protected: &[bool]) -> Vec<u32> {
    let count = mesh.faces().last().map_or(0, |f| f.index() + 1);
    let mut patch = vec![u32::MAX; count];
    let mut next = 0u32;
    let mut stack = Vec::new();
    for seed in mesh.faces() {
        if patch[seed.index()] != u32::MAX {
            continue;
        }
        patch[seed.index()] = next;
        stack.push(seed);
        while let Some(f) = stack.pop() {
            for h in mesh.face_halfedges(f) {
                let e: EdgeId = mesh.edge(h);
                if protected[e.index()] {
                    continue;
                }
                if let Some(g) = mesh.face(mesh.opposite(h)) {
                    if patch[g.index()] == u32::MAX {
                        patch[g.index()] = next;
                        stack.push(g);
                    }
                }
            }
        }
        next += 1;
    }
    patch
}

/// Mean absolute distance of vertex valences from optimal, over vertices
/// with at least one edge.
fn valence_deviation(mesh: &HalfedgeMesh) -> Scalar {
    let mut total = 0i64;
    let mut count = 0usize;
    for v in mesh.vertices().filter(|&v| !mesh.is_isolated(v)) {
        let optimal = if mesh.is_boundary_vertex(v) { 4 } else { 6 };
        total += (mesh.degree(v) as i64 - optimal).abs();
        count += 1;
    }
    if count == 0 {
        0.0
    } else {
        total as Scalar / count as Scalar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protected_edges_cut_a_box_into_its_six_faces() {
        let p = |x: Scalar, y: Scalar, z: Scalar| Point3::new(x, y, z);
        let positions = vec![
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
            p(1.0, 0.0, 1.0),
            p(1.0, 1.0, 1.0),
            p(0.0, 1.0, 1.0),
        ];
        let indices = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ];
        let mesh = TriMesh::new(positions, indices);
        let he = HalfedgeMesh::from_tri_mesh(&mesh).expect("closed box");
        let triangles: Vec<[Point3; 3]> = mesh
            .triangles()
            .map(|t| t.map(|i| mesh.positions[i as usize]))
            .collect();
        let protected =
            protected_edges(&he, &triangles, Some(RemeshOptions::DEFAULT_FEATURE_ANGLE));
        assert_eq!(protected.iter().filter(|&&p| p).count(), 12);
        let patch = patches(&he, &protected);
        assert_eq!(patch, vec![0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5]);
        // Without features, one patch.
        let protected = protected_edges(&he, &triangles, None);
        assert_eq!(patches(&he, &protected), vec![0; 12]);
    }
}
