//! Minkowski sum and difference of planar-faced solids.
//!
//! # The convex case is the tractable one
//!
//! The Minkowski sum of two convex polyhedra is the convex hull of their
//! pairwise vertex sums. That is exact, needs no boolean solver, and is
//! what [`minkowski_sum`] computes.
//!
//! Non-convex operands have no such shortcut. The standard construction
//! decomposes both into convex parts, sums every pair, and unions the
//! results — which is why [`minkowski_sum_with`] takes a boolean provider.
//! Cost is the product of the part counts, so it is budgeted rather than
//! left to run away.
//!
//! # Difference is not the sum run backwards
//!
//! `A ⊖ B` is the set of translations that keep `B` inside `A`:
//! `{ x : x + B ⊆ A }`. It is an erosion, not a hull of pairwise
//! differences, and computing it as one is a common and silent error.
//!
//! When `A` is convex the containment test reduces to the vertices of `B`,
//! because a convex `A` containing every `x + vᵢ` contains their hull and
//! therefore all of `x + B`. That gives
//!
//! ```text
//! A ⊖ B = ⋂ᵢ (A − vᵢ)   over the vertices vᵢ of B
//! ```
//!
//! which is exact and computable with intersections alone. For non-convex
//! `A` the reduction fails — `A` can contain each translated vertex while
//! missing the material between them — so [`minkowski_difference_with`]
//! refuses a non-convex subject by name rather than returning a result that
//! is too large.
//!
//! # Curved operands
//!
//! Both operations are defined here for planar-faced solids. A curved
//! operand is refused, matching the rest of the offset work: the sum of two
//! curved solids is not a polyhedron and approximating it silently would
//! misreport what the result is.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Point3, Scalar, Tolerance, Vec3};
use axiolid_decompose::{convex_decompose_with, split::Splitter, Strategy};
use axiolid_mesh::{audit_mesh, TriMesh};
use axiolid_mesh_boolean_contract::MeshBoolean;
use thiserror::Error;

/// Why a Minkowski operation could not be performed.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum MinkowskiError {
    /// An operand is not a closed two-manifold solid.
    #[error("{operand} is not a closed two-manifold solid: {boundary} boundary and {non_manifold} non-manifold edges")]
    NotASolid {
        /// Which operand failed: `"subject"` or `"tool"`.
        operand: &'static str,
        /// Edges with a single incident triangle.
        boundary: usize,
        /// Edges with more than two incident triangles.
        non_manifold: usize,
    },
    /// An operand has no vertices to sum.
    #[error("{0} is empty")]
    EmptyOperand(&'static str),
    /// [`minkowski_sum`] was given a non-convex operand.
    ///
    /// The convex-only entry point refuses rather than silently returning
    /// the hull, which for a non-convex operand is strictly too large.
    #[error("{operand} is not convex; use minkowski_sum_with to decompose it")]
    NotConvex {
        /// Which operand failed.
        operand: &'static str,
    },
    /// The erosion's subject is not convex.
    ///
    /// Refused rather than approximated: for a non-convex subject the
    /// vertex-wise containment test admits translations that do not
    /// actually fit, so the result would be too large.
    #[error("erosion requires a convex subject, and this one is not")]
    ErosionSubjectNotConvex,
    /// The decomposition of an operand failed.
    #[error("decomposing {operand} failed: {reason}")]
    DecompositionFailed {
        /// Which operand failed.
        operand: &'static str,
        /// What the decomposer reported.
        reason: String,
    },
    /// A hull could not be built from a pairwise vertex sum.
    #[error("the hull of a pairwise sum failed: {0}")]
    HullFailed(String),
    /// A boolean step failed.
    #[error("combining parts failed: {0}")]
    BooleanFailed(String),
    /// The pairwise product exceeds the budget.
    ///
    /// Reported rather than run: the work is the product of the two part
    /// counts, so a modestly non-convex pair can be very expensive, and a
    /// caller deserves to know before it starts rather than after.
    #[error("the decomposition needs {pairs} pairwise sums, over the {limit} budget")]
    BudgetExceeded {
        /// Pairwise sums the request would have taken.
        pairs: usize,
        /// The cap that was not raised.
        limit: usize,
    },
}

/// Largest number of pairwise convex sums before the work is refused.
const MAX_PAIRS: usize = 4096;

/// What was done to produce a Minkowski result.
///
/// Reported rather than implied: a caller that asked for a sum of two
/// solids it believed convex should be able to see whether decomposition
/// was needed, and how much work it cost.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct MinkowskiEvidence {
    /// Convex parts the subject was split into.
    pub subject_parts: usize,
    /// Convex parts the tool was split into.
    pub tool_parts: usize,
    /// Pairwise convex sums performed.
    pub pairwise_sums: usize,
    /// Boolean operations used to combine them.
    pub boolean_operations: usize,
}

impl MinkowskiEvidence {
    /// Whether both operands were already convex.
    pub fn was_convex(&self) -> bool {
        self.subject_parts == 1 && self.tool_parts == 1
    }
}

/// A Minkowski result and what it took to produce.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct MinkowskiOutcome {
    /// The resulting solid.
    pub mesh: TriMesh,
    /// What was done.
    pub evidence: MinkowskiEvidence,
}

/// Minkowski sum of two CONVEX planar-faced solids.
///
/// The result is the convex hull of the pairwise vertex sums, which is
/// exact for convex operands and needs no boolean solver.
///
/// # Errors
///
/// Refuses an operand that is not a closed two-manifold solid, is empty, or
/// is not convex. For non-convex operands use [`minkowski_sum_with`].
pub fn minkowski_sum(
    subject: &TriMesh,
    tool: &TriMesh,
    tolerance: Tolerance,
) -> Result<TriMesh, MinkowskiError> {
    require_solid(subject, "subject", tolerance)?;
    require_solid(tool, "tool", tolerance)?;

    // Convexity is checked by asking the decomposer: a convex solid is its
    // own single part. Reusing that keeps one definition of convex in the
    // codebase rather than a second, subtly different one here.
    if !is_convex(subject, tolerance)? {
        return Err(MinkowskiError::NotConvex { operand: "subject" });
    }
    if !is_convex(tool, tolerance)? {
        return Err(MinkowskiError::NotConvex { operand: "tool" });
    }
    convex_sum(&subject.positions, &tool.positions)
}

/// Minkowski sum of two planar-faced solids, convex or not.
///
/// Non-convex operands are decomposed into convex parts; every pair is
/// summed and the results unioned through `provider`.
///
/// # Errors
///
/// As [`minkowski_sum`] for malformed operands, plus decomposition, hull,
/// boolean, and budget failures.
pub fn minkowski_sum_with(
    subject: &TriMesh,
    tool: &TriMesh,
    tolerance: Tolerance,
    provider: &dyn MeshBoolean,
) -> Result<MinkowskiOutcome, MinkowskiError> {
    require_solid(subject, "subject", tolerance)?;
    require_solid(tool, "tool", tolerance)?;

    let splitter = Splitter::Provider(provider);
    let subject_parts = decompose(subject, "subject", tolerance, &splitter)?;
    let tool_parts = decompose(tool, "tool", tolerance, &splitter)?;

    let pairs = subject_parts.len().saturating_mul(tool_parts.len());
    if pairs > MAX_PAIRS {
        return Err(MinkowskiError::BudgetExceeded {
            pairs,
            limit: MAX_PAIRS,
        });
    }

    let mut summed: Vec<TriMesh> = Vec::with_capacity(pairs);
    for left in &subject_parts {
        for right in &tool_parts {
            summed.push(convex_sum(&left.positions, &right.positions)?);
        }
    }

    let options = ExecutionOptions::new(tolerance);
    let mut boolean_operations = 0usize;
    let mut result = summed
        .first()
        .cloned()
        .ok_or(MinkowskiError::EmptyOperand("subject"))?;
    for piece in summed.iter().skip(1) {
        let outcome = provider
            .boolean(&result, piece, BooleanOperator::Union, &options)
            .map_err(|error| MinkowskiError::BooleanFailed(error.to_string()))?;
        boolean_operations += 1;
        result = outcome.mesh;
    }

    Ok(MinkowskiOutcome {
        mesh: result,
        evidence: MinkowskiEvidence {
            subject_parts: subject_parts.len(),
            tool_parts: tool_parts.len(),
            pairwise_sums: summed.len(),
            boolean_operations,
        },
    })
}

/// Minkowski difference: the translations of `tool` that stay inside
/// `subject`.
///
/// Computed as the intersection of `subject` translated by the negated
/// vertices of `tool`, which is exact when `subject` is convex.
///
/// # Errors
///
/// Refuses malformed operands, and a non-convex subject by name: the
/// vertex-wise containment test is only sufficient for a convex subject, so
/// any other answer would be too large.
pub fn minkowski_difference_with(
    subject: &TriMesh,
    tool: &TriMesh,
    tolerance: Tolerance,
    provider: &dyn MeshBoolean,
) -> Result<MinkowskiOutcome, MinkowskiError> {
    require_solid(subject, "subject", tolerance)?;
    require_solid(tool, "tool", tolerance)?;

    if !is_convex(subject, tolerance)? {
        return Err(MinkowskiError::ErosionSubjectNotConvex);
    }

    // Only the tool's distinct vertices matter: containment of the hull
    // follows from containment of its extreme points, and repeating a
    // vertex would repeat an identical intersection.
    let offsets = distinct_points(&tool.positions, tolerance);
    let options = ExecutionOptions::new(tolerance);

    let mut result = translated(subject, -offsets[0]);
    let mut boolean_operations = 0usize;
    for offset in offsets.iter().skip(1) {
        let shifted = translated(subject, -*offset);
        let outcome = provider
            .boolean(&result, &shifted, BooleanOperator::Intersection, &options)
            .map_err(|error| MinkowskiError::BooleanFailed(error.to_string()))?;
        boolean_operations += 1;
        result = outcome.mesh;
        // An empty intersection is a legitimate answer: the tool does not
        // fit inside the subject in any translation. Stopping early avoids
        // intersecting an empty solid, which many backends reject.
        if result.indices.is_empty() {
            break;
        }
    }

    Ok(MinkowskiOutcome {
        mesh: result,
        evidence: MinkowskiEvidence {
            subject_parts: 1,
            tool_parts: offsets.len(),
            pairwise_sums: 0,
            boolean_operations,
        },
    })
}

/// Convex hull of every pairwise sum of two point sets.
fn convex_sum(left: &[Point3], right: &[Point3]) -> Result<TriMesh, MinkowskiError> {
    let mut sums = Vec::with_capacity(left.len() * right.len());
    for a in left {
        for b in right {
            sums.push(*a + Vec3::new(b.x, b.y, b.z));
        }
    }
    axiolid_construct::hull::convex_hull(&sums)
        .map_err(|error| MinkowskiError::HullFailed(error.to_string()))
}

/// Translate every vertex of a mesh.
fn translated(mesh: &TriMesh, offset: Vec3) -> TriMesh {
    TriMesh::new(
        mesh.positions.iter().map(|p| *p + offset).collect(),
        mesh.indices.clone(),
    )
}

/// Distinct positions of a mesh, as offset vectors.
///
/// Never empty: callers have already established the mesh is a solid, so it
/// has vertices.
fn distinct_points(positions: &[Point3], tolerance: Tolerance) -> Vec<Vec3> {
    let step = tolerance.linear().max(Scalar::EPSILON);
    let mut seen = std::collections::BTreeSet::new();
    let mut out = Vec::new();
    for p in positions {
        let key = (
            (p.x / step).round() as i64,
            (p.y / step).round() as i64,
            (p.z / step).round() as i64,
        );
        if seen.insert(key) {
            out.push(Vec3::new(p.x, p.y, p.z));
        }
    }
    out
}

/// Whether a solid is convex, as the decomposer defines it.
fn is_convex(mesh: &TriMesh, tolerance: Tolerance) -> Result<bool, MinkowskiError> {
    let decomposition = axiolid_decompose::convex_decompose(mesh, Strategy::Exact, tolerance)
        .map_err(|error| MinkowskiError::DecompositionFailed {
            operand: "operand",
            reason: error.to_string(),
        })?;
    Ok(decomposition.is_single_part())
}

/// Decompose an operand into convex parts.
fn decompose(
    mesh: &TriMesh,
    operand: &'static str,
    tolerance: Tolerance,
    splitter: &Splitter<'_>,
) -> Result<Vec<TriMesh>, MinkowskiError> {
    convex_decompose_with(mesh, Strategy::Exact, tolerance, splitter)
        .map(|decomposition| decomposition.parts)
        .map_err(|error| MinkowskiError::DecompositionFailed {
            operand,
            reason: error.to_string(),
        })
}

/// Refuse an operand that is not a closed two-manifold solid.
fn require_solid(
    mesh: &TriMesh,
    operand: &'static str,
    tolerance: Tolerance,
) -> Result<(), MinkowskiError> {
    if mesh.positions.is_empty() || mesh.indices.is_empty() {
        return Err(MinkowskiError::EmptyOperand(operand));
    }
    let health = audit_mesh(mesh, tolerance);
    if !health.is_closed_two_manifold() {
        return Err(MinkowskiError::NotASolid {
            operand,
            boundary: health.boundary_edges,
            non_manifold: health.non_manifold_edges,
        });
    }
    Ok(())
}
