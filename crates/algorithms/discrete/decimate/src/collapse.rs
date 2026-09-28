//! Edge collapse with validity checks and a measured deviation bound.

use axiolid_core::{Point3, Tolerance};
use axiolid_mesh::TriMesh;
use std::collections::{BTreeMap, BTreeSet};
use thiserror::Error;

/// What the caller wants back.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum DecimateTarget {
    /// Reduce until at most this many triangles remain.
    ///
    /// The deviation bound still applies: the budget is honoured only as far
    /// as it can be without exceeding `max_deviation`.
    TriangleBudget(usize),
    /// Collapse every edge whose removal stays within this deviation.
    MaxDeviation(f64),
}

/// Why a decimation could not run.
#[derive(Debug, Clone, PartialEq, Error)]
#[non_exhaustive]
pub enum DecimateError {
    /// The index buffer is not a whole number of triangles.
    #[error("index buffer length {0} is not a multiple of 3")]
    RaggedIndices(usize),
    /// A triangle references a vertex that does not exist.
    #[error("triangle {0} references vertex {1}, which is out of range")]
    IndexOutOfRange(usize, u32),
    /// A deviation bound must be a positive, finite length.
    #[error("deviation bound {0} is not a positive finite length")]
    InvalidBound(f64),
}

/// What a decimation actually did.
///
/// The deviation is measured, not estimated: it is the largest distance any
/// collapsed vertex moved from its original position. A caller that asked
/// for 1mm accuracy can check it got it.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct DecimateReport {
    /// Triangles before.
    pub input_triangles: usize,
    /// Triangles after.
    pub output_triangles: usize,
    /// Edge collapses performed.
    pub collapses: usize,
    /// Collapses rejected because they would have created a defect.
    ///
    /// A non-zero count is not a failure: it is the reason the result is
    /// still a valid manifold rather than a faster, broken one.
    pub rejected_unsafe: usize,
    /// Collapses rejected because they would have exceeded the bound.
    pub rejected_deviation: usize,
    /// Largest distance any vertex moved, in model units.
    pub max_deviation: f64,
}

impl DecimateReport {
    /// Whether the mesh was left untouched.
    pub fn is_noop(&self) -> bool {
        self.collapses == 0
    }
}

/// Decimate a triangle mesh by edge collapse.
///
/// Deterministic: candidate edges are ordered by length then by vertex
/// index, so the same input and target produce the same output on every
/// run, matching the determinism discipline the plan contract established
/// in v0.6.
///
/// # Errors
///
/// Refuses a ragged index buffer, out-of-range indices, and a non-positive
/// deviation bound.
pub fn decimate(
    mesh: &TriMesh,
    target: DecimateTarget,
    tolerance: Tolerance,
) -> Result<(TriMesh, DecimateReport), DecimateError> {
    if mesh.indices.len() % 3 != 0 {
        return Err(DecimateError::RaggedIndices(mesh.indices.len()));
    }
    let vertex_count = mesh.positions.len();
    for (t, chunk) in mesh.indices.chunks_exact(3).enumerate() {
        for &index in chunk {
            if index as usize >= vertex_count {
                return Err(DecimateError::IndexOutOfRange(t, index));
            }
        }
    }

    let bound = match target {
        DecimateTarget::MaxDeviation(d) => {
            if !d.is_finite() || d <= 0.0 {
                return Err(DecimateError::InvalidBound(d));
            }
            d
        }
        // A budget still needs a ceiling, or "reduce to N triangles" would
        // licence arbitrary damage. The caller's tolerance is that ceiling.
        DecimateTarget::TriangleBudget(_) => tolerance.linear(),
    };

    let budget = match target {
        DecimateTarget::TriangleBudget(n) => n,
        DecimateTarget::MaxDeviation(_) => 0,
    };

    run(mesh, budget, bound)
}

/// The collapse loop.
///
/// Each candidate edge is collapsed to its midpoint. The move is accepted
/// only when it neither exceeds the deviation bound nor creates a defect.
fn run(
    mesh: &TriMesh,
    budget: usize,
    bound: f64,
) -> Result<(TriMesh, DecimateReport), DecimateError> {
    let input_triangles = mesh.indices.len() / 3;
    let mut positions = mesh.positions.clone();
    let mut triangles: Vec<[u32; 3]> = mesh
        .indices
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();

    // Deterministic candidate order: shortest edges first, ties broken by
    // vertex index. Sorting by a float alone would leave equal-length edges
    // in hash order and make the output depend on iteration chance.
    let mut candidates: Vec<(u32, u32)> = unique_edges(&triangles).into_iter().collect();
    candidates.sort_by(|a, b| {
        let la = (positions[a.0 as usize] - positions[a.1 as usize]).length();
        let lb = (positions[b.0 as usize] - positions[b.1 as usize]).length();
        la.partial_cmp(&lb)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(b))
    });

    let mut report = DecimateReport {
        input_triangles,
        output_triangles: input_triangles,
        collapses: 0,
        rejected_unsafe: 0,
        rejected_deviation: 0,
        max_deviation: 0.0,
    };
    let mut moved = vec![0.0_f64; positions.len()];
    let mut alive: Vec<bool> = vec![true; positions.len()];

    for (u, v) in candidates {
        if triangles.len() <= budget.max(4) && budget > 0 {
            break;
        }
        if !alive[u as usize] || !alive[v as usize] {
            continue;
        }
        let midpoint = (positions[u as usize] + positions[v as usize]) / 2.0;

        // Deviation is cumulative: a vertex that already moved carries that
        // history, so repeated collapses cannot drift past the bound one
        // small step at a time.
        let deviation = (midpoint - positions[u as usize])
            .length()
            .max((midpoint - positions[v as usize]).length())
            + moved[u as usize].max(moved[v as usize]);
        if deviation > bound {
            report.rejected_deviation += 1;
            continue;
        }

        match try_collapse(&triangles, &positions, u, v, midpoint) {
            Some(next) => {
                triangles = next;
                positions[u as usize] = midpoint;
                moved[u as usize] = deviation;
                alive[v as usize] = false;
                report.collapses += 1;
                report.max_deviation = report.max_deviation.max(deviation);
            }
            None => report.rejected_unsafe += 1,
        }
    }

    report.output_triangles = triangles.len();
    Ok((compact(&triangles, &positions), report))
}

/// Attempt one collapse, returning the new triangle list or nothing.
///
/// Rejects the two ways a collapse damages a mesh:
///
/// - **Inversion.** A triangle whose normal flips has turned inside out. A
///   decimator that permits this produces exactly the defect
///   `axiolid-heal`'s `OrientOutward` exists to repair.
/// - **Non-manifold edges.** Collapsing an edge whose endpoints share
///   neighbours other than the two triangles on it welds unrelated sheets
///   together.
///
/// Boundary vertices are not pinned: a collapse may move a vertex on an open
/// boundary to an edge midpoint like any other vertex, and the deviation
/// bound is what limits how far the silhouette moves.
fn try_collapse(
    triangles: &[[u32; 3]],
    positions: &[Point3],
    u: u32,
    v: u32,
    midpoint: Point3,
) -> Option<Vec<[u32; 3]>> {
    // The link condition: the shared neighbourhood of u and v must be
    // exactly the two triangles on the edge. More than that, and the
    // collapse creates a non-manifold edge.
    let nu = neighbours(triangles, u);
    let nv = neighbours(triangles, v);
    let shared = nu.intersection(&nv).count();
    if shared != 2 {
        return None;
    }

    let mut next = Vec::with_capacity(triangles.len());
    for &t in triangles {
        let touches_u = t.contains(&u);
        let touches_v = t.contains(&v);
        if touches_u && touches_v {
            // Degenerates to a line: this is the triangle being removed.
            continue;
        }
        let mapped = t.map(|c| if c == v { u } else { c });
        if touches_u || touches_v {
            let before = normal(positions, t);
            let after = normal_with(positions, mapped, u, midpoint);
            // A flipped normal means the triangle turned inside out. Zero
            // area after the move is equally unacceptable: it contributes
            // nothing and confuses every downstream audit.
            if after.length_squared() == 0.0 || before.dot(after) <= 0.0 {
                return None;
            }
        }
        next.push(mapped);
    }
    Some(next)
}

/// Vertices sharing a triangle with `vertex`.
fn neighbours(triangles: &[[u32; 3]], vertex: u32) -> BTreeSet<u32> {
    let mut set = BTreeSet::new();
    for t in triangles {
        if t.contains(&vertex) {
            for &c in t {
                if c != vertex {
                    set.insert(c);
                }
            }
        }
    }
    set
}

/// Every undirected edge, each appearing once.
fn unique_edges(triangles: &[[u32; 3]]) -> BTreeSet<(u32, u32)> {
    let mut set = BTreeSet::new();
    for t in triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            set.insert((a.min(b), a.max(b)));
        }
    }
    set
}

/// Unnormalised triangle normal.
fn normal(positions: &[Point3], t: [u32; 3]) -> axiolid_core::Vec3 {
    let a = positions[t[0] as usize];
    let b = positions[t[1] as usize];
    let c = positions[t[2] as usize];
    (b - a).cross(c - a)
}

/// Triangle normal with one vertex moved to a candidate position.
fn normal_with(
    positions: &[Point3],
    t: [u32; 3],
    moved_index: u32,
    moved_to: Point3,
) -> axiolid_core::Vec3 {
    let at = |i: u32| {
        if i == moved_index {
            moved_to
        } else {
            positions[i as usize]
        }
    };
    let (a, b, c) = (at(t[0]), at(t[1]), at(t[2]));
    (b - a).cross(c - a)
}

/// Drop unreferenced vertices and renumber.
///
/// Collapsed vertices leave holes in the position array. Emitting them would
/// leave unreferenced vertices that make the result look defective to the
/// v0.7 diagnosis.
fn compact(triangles: &[[u32; 3]], positions: &[Point3]) -> TriMesh {
    let mut remap: BTreeMap<u32, u32> = BTreeMap::new();
    let mut kept = Vec::new();
    let mut indices = Vec::with_capacity(triangles.len() * 3);
    for t in triangles {
        for &corner in t {
            let next = u32::try_from(kept.len()).unwrap_or(u32::MAX);
            let slot = *remap.entry(corner).or_insert_with(|| {
                kept.push(positions[corner as usize]);
                next
            });
            indices.push(slot);
        }
    }
    TriMesh::new(kept, indices)
}
