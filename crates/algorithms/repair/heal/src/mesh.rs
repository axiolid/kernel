//! Diagnosing and repairing triangle meshes.
//!
//! `MeshHealth` already counts what is wrong with a mesh. A count is not
//! actionable: "4444 inconsistent winding edges" does not say which
//! triangles to fix. This locates each defect against a stable element
//! index so a caller can act on it, and applies only the repairs it was
//! explicitly asked for.

use std::collections::{HashMap, VecDeque};

use axiolid_core::{Scalar, Tolerance};
use axiolid_measure::volume_properties;
use axiolid_mesh::{EdgeAdjacency, TriMesh};

use crate::diagnosis::{Defect, DefectKind, Diagnosis};
use crate::repair::{RepairAction, RepairPlan, RepairReport};
use crate::traits::{Diagnose, Repair};

/// Diagnosis and repair for indexed triangle meshes.
#[derive(Debug, Default, Clone, Copy)]
pub struct MeshHealer;

/// Failure to complete diagnosis or repair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshHealError {
    /// The index buffer is not a whole number of triangles.
    RaggedIndices {
        /// Number of index entries found.
        indices: usize,
    },
}

impl core::fmt::Display for MeshHealError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::RaggedIndices { indices } => {
                write!(f, "index buffer of {indices} is not a multiple of three")
            }
        }
    }
}

impl std::error::Error for MeshHealError {}

impl Diagnose<TriMesh> for MeshHealer {
    type Error = MeshHealError;

    fn diagnose(&self, mesh: &TriMesh, tolerance: Tolerance) -> Result<Diagnosis, Self::Error> {
        if !mesh.indices.len().is_multiple_of(3) {
            return Err(MeshHealError::RaggedIndices {
                indices: mesh.indices.len(),
            });
        }
        let mut defects = Vec::new();
        let count = mesh.indices.len() / 3;
        let limit = tolerance.linear() * tolerance.linear();

        // Degenerate triangles first: they are excluded from the edge
        // topology below, because a sliver's edges are not meaningful
        // adjacency and would produce misleading manifold defects.
        let mut usable = Vec::with_capacity(count);
        for t in 0..count {
            match area(mesh, t) {
                Some(a) if a > limit => usable.push(t),
                Some(_) => defects.push(Defect {
                    kind: DefectKind::DegenerateElement,
                    element: Some(t as u32),
                    detail: None,
                }),
                None => defects.push(Defect {
                    kind: DefectKind::DegenerateElement,
                    element: Some(t as u32),
                    detail: Some("triangle corner does not address a position".to_owned()),
                }),
            }
        }

        // Edge topology over usable triangles only. The key is the
        // undirected edge; the stored direction is what reveals winding.
        //
        // This deliberately does NOT use `axiolid_mesh::EdgeAdjacency`.
        // That type excludes triangles with a repeated corner -- a purely
        // combinatorial test. Healing excludes triangles whose *area* is
        // below tolerance, which is strictly stronger: a sliver with three
        // distinct corners is sound to `EdgeAdjacency` and a defect here.
        // Sharing the structure would mean diagnosing slivers' edges as
        // real adjacency and reporting manifold defects that do not exist.
        //
        // The duplication is the tolerance-aware filter, not the edge map.
        let mut edges: HashMap<(u32, u32), Vec<(usize, bool)>> = HashMap::new();
        for &t in &usable {
            for (a, b) in corners(mesh, t) {
                let forward = a < b;
                let key = if forward { (a, b) } else { (b, a) };
                edges.entry(key).or_default().push((t, forward));
            }
        }

        // Sort for deterministic defect order: HashMap iteration is not
        // stable, and a diagnosis that reorders between runs is unusable
        // as an audit record.
        let mut keys: Vec<_> = edges.keys().copied().collect();
        keys.sort_unstable();

        for key in keys {
            let uses = &edges[&key];
            match uses.len() {
                1 => {
                    defects.push(Defect {
                        kind: DefectKind::OpenShell,
                        element: Some(uses[0].0 as u32),
                        detail: Some(format!("edge {}-{} has one incident face", key.0, key.1)),
                    });
                }
                2 => {
                    // Two faces sharing the SAME directed edge disagree
                    // about which side is outside. Consistent neighbours
                    // traverse a shared edge in opposite directions.
                    if uses[0].1 == uses[1].1 {
                        defects.push(Defect {
                            kind: DefectKind::InconsistentOrientation,
                            element: Some(uses[1].0 as u32),
                            detail: Some(format!(
                                "triangles {} and {} traverse edge {}-{} the same way",
                                uses[0].0, uses[1].0, key.0, key.1
                            )),
                        });
                    }
                }
                n => {
                    defects.push(Defect {
                        kind: DefectKind::NonManifoldEdge,
                        element: Some(uses[0].0 as u32),
                        detail: Some(format!("edge {}-{} has {n} incident faces", key.0, key.1)),
                    });
                }
            }
        }

        // Coincident positions stored separately. Bucketing by a
        // tolerance-sized cell finds them in one pass; comparing every
        // pair would be quadratic on model-scale meshes.
        for (representative, group) in coincident_groups(mesh, tolerance) {
            for duplicate in group {
                defects.push(Defect {
                    kind: DefectKind::DuplicateVertex,
                    element: Some(duplicate),
                    detail: Some(format!("coincident with vertex {representative}")),
                });
            }
        }

        Ok(Diagnosis { defects })
    }
}

/// Group coincident vertices, returning `(representative, duplicates)`.
///
/// Vertices are bucketed into tolerance-sized cells. A vertex is compared
/// only against the 27 cells touching its own, so a pair straddling a cell
/// boundary is still found; comparing all pairs would be quadratic.
fn coincident_groups(mesh: &TriMesh, tolerance: Tolerance) -> Vec<(u32, Vec<u32>)> {
    let eps = tolerance.linear();
    if eps <= 0.0 {
        return Vec::new();
    }
    let cell = |p: axiolid_core::Point3| {
        (
            (p.x / eps).floor() as i64,
            (p.y / eps).floor() as i64,
            (p.z / eps).floor() as i64,
        )
    };
    let mut buckets: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for (i, p) in mesh.positions.iter().enumerate() {
        buckets.entry(cell(*p)).or_default().push(i as u32);
    }

    let eps2 = eps * eps;
    let mut owner: Vec<Option<u32>> = vec![None; mesh.positions.len()];
    let mut groups: Vec<(u32, Vec<u32>)> = Vec::new();
    for i in 0..mesh.positions.len() as u32 {
        if owner[i as usize].is_some() {
            continue;
        }
        let p = mesh.positions[i as usize];
        let (cx, cy, cz) = cell(p);
        let mut mates = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(bucket) = buckets.get(&(cx + dx, cy + dy, cz + dz)) else {
                        continue;
                    };
                    for &j in bucket {
                        if j > i
                            && owner[j as usize].is_none()
                            && (mesh.positions[j as usize] - p).length_squared() <= eps2
                        {
                            mates.push(j);
                        }
                    }
                }
            }
        }
        if !mates.is_empty() {
            mates.sort_unstable();
            owner[i as usize] = Some(i);
            for &j in &mates {
                owner[j as usize] = Some(i);
            }
            groups.push((i, mates));
        }
    }
    groups
}

impl Repair<TriMesh> for MeshHealer {
    type Error = MeshHealError;

    fn repair(
        &self,
        mesh: &TriMesh,
        plan: &RepairPlan,
        tolerance: Tolerance,
    ) -> Result<(TriMesh, RepairReport), Self::Error> {
        if !mesh.indices.len().is_multiple_of(3) {
            return Err(MeshHealError::RaggedIndices {
                indices: mesh.indices.len(),
            });
        }
        let mut out = mesh.clone();
        let mut report = RepairReport::default();
        // Actions run in the caller's order. The plan is ordered on
        // purpose: welding before orientation gives orientation a
        // connected mesh to work with, and the reverse does not.
        for &action in &plan.actions {
            let changed = match action {
                RepairAction::WeldVertices => weld(&mut out, tolerance),
                RepairAction::DropDegenerateElements => drop_degenerate(&mut out, tolerance),
                RepairAction::UnifyOrientation => unify_orientation(&mut out),
                RepairAction::OrientOutward => orient_outward(&mut out, tolerance),
            };
            if changed {
                report.applied.push(action);
            } else {
                report.skipped.push(action);
            }
        }
        Ok((out, report))
    }
}

/// Merge coincident vertices and repoint the index buffer.
///
/// Positions are compacted rather than left orphaned: a welded mesh that
/// still carries unreferenced vertices reports the same duplicate defects
/// on the next diagnosis, which would make the repair look ineffective.
fn weld(mesh: &mut TriMesh, tolerance: Tolerance) -> bool {
    let groups = coincident_groups(mesh, tolerance);
    if groups.is_empty() {
        return false;
    }
    let mut remap: Vec<u32> = (0..mesh.positions.len() as u32).collect();
    for (representative, duplicates) in groups {
        for d in duplicates {
            remap[d as usize] = representative;
        }
    }
    let mut keep: Vec<u32> = Vec::new();
    let mut compact = vec![u32::MAX; mesh.positions.len()];
    for i in 0..mesh.positions.len() as u32 {
        if remap[i as usize] == i {
            compact[i as usize] = keep.len() as u32;
            keep.push(i);
        }
    }
    mesh.positions = keep.iter().map(|&i| mesh.positions[i as usize]).collect();
    for index in &mut mesh.indices {
        *index = compact[remap[*index as usize] as usize];
    }
    // Welding can collapse a triangle to a line; those are degenerate now,
    // but removing them is a different action the caller did not request.
    true
}

/// Remove triangles at or below the tolerance area.
fn drop_degenerate(mesh: &mut TriMesh, tolerance: Tolerance) -> bool {
    let limit = tolerance.linear() * tolerance.linear();
    let count = mesh.indices.len() / 3;
    let mut kept = Vec::with_capacity(mesh.indices.len());
    for t in 0..count {
        if area(mesh, t).is_some_and(|a| a > limit) {
            kept.extend_from_slice(&mesh.indices[t * 3..t * 3 + 3]);
        }
    }
    if kept.len() == mesh.indices.len() {
        return false;
    }
    mesh.indices = kept;
    true
}

/// Make connected-face winding consistent by flood fill.
///
/// Two triangles sharing an edge agree when they traverse that edge in
/// OPPOSITE directions. Starting from a seed, each neighbour that
/// traverses the shared edge the same way is flipped, and the fill
/// continues. This is the defect that a closed manifold mesh can still
/// carry: signed volume comes out negative or partly cancelled while
/// every boundary/manifold count looks perfect.
///
/// Orientation is fixed per connected component; the absolute sense of
/// each component is left as its seed found it, because choosing outward
/// requires a volume convention this function does not own.
fn unify_orientation(mesh: &mut TriMesh) -> bool {
    let count = mesh.indices.len() / 3;
    if count == 0 {
        return false;
    }
    // Adjacency is derived once by the shared structure. Flipping a triangle
    // reverses its corner order but not its undirected edge set, so the map
    // stays valid for the whole traversal.
    let adjacency = EdgeAdjacency::build(mesh);
    let mut visited = vec![false; count];
    let mut flipped = false;
    for seed in 0..count {
        if visited[seed] {
            continue;
        }
        visited[seed] = true;
        let mut queue = VecDeque::from([seed]);
        while let Some(t) = queue.pop_front() {
            for (a, b) in corners(mesh, t) {
                for n in adjacency.triangle_neighbours(t) {
                    if visited[n] {
                        continue;
                    }
                    // Only the neighbour across THIS corner is judged here;
                    // the others are reached on their own corner.
                    let shares = corners(mesh, n)
                        .into_iter()
                        .any(|(c, d)| (c.min(d), c.max(d)) == (a.min(b), a.max(b)));
                    if !shares {
                        continue;
                    }
                    // Same directed edge means the neighbour winds the
                    // same way round a shared edge, which is inconsistent.
                    if corners(mesh, n).into_iter().any(|(c, d)| c == a && d == b) {
                        mesh.indices.swap(n * 3 + 1, n * 3 + 2);
                        flipped = true;
                    }
                    visited[n] = true;
                    queue.push_back(n);
                }
            }
        }
    }
    flipped
}

/// The three directed corner pairs of a triangle.
fn corners(mesh: &TriMesh, t: usize) -> [(u32, u32); 3] {
    let i = &mesh.indices[t * 3..t * 3 + 3];
    [(i[0], i[1]), (i[1], i[2]), (i[2], i[0])]
}

/// Triangle area, or `None` when a corner index is out of range.
fn area(mesh: &TriMesh, t: usize) -> Option<Scalar> {
    let i = &mesh.indices[t * 3..t * 3 + 3];
    let p = |k: usize| mesh.positions.get(i[k] as usize).copied();
    let (a, b, c) = (p(0)?, p(1)?, p(2)?);
    Some((b - a).cross(c - a).length() * 0.5)
}

/// Flip a closed shell whose faces point inward.
///
/// `unify_orientation` makes neighbours agree but leaves the absolute
/// sense as the seed found it, because choosing outward needs a volume
/// convention it does not own. `axiolid-measure` owns that convention
/// now, so this repair completes the pair: unify makes the shell
/// consistent, this makes it consistent the RIGHT WAY ROUND.
///
/// An inside-out shell is structurally perfect -- closed, two-manifold,
/// consistently wound -- so no topological audit finds it. The boolmesh
/// provider records the consequence: `Difference` behaves as `Union` and
/// returns a LARGER mesh with no error.
///
/// Only closed shells are eligible. An open surface has no enclosed
/// volume, so `inward` is not defined for it and it is left untouched.
fn orient_outward(mesh: &mut TriMesh, tolerance: Tolerance) -> bool {
    let Ok(properties) = volume_properties(&*mesh, tolerance) else {
        return false;
    };
    if properties.signed_volume >= 0.0 {
        return false;
    }
    for triangle in mesh.indices.chunks_exact_mut(3) {
        triangle.swap(1, 2);
    }
    true
}
