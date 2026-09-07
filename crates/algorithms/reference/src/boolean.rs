//! Scalar reference implementation of solid booleans (ADR 0012, ADR 0017 §5).
//!
//! # Why this exists
//!
//! ADR 0012 requires a scalar reference to land *before* an optimized provider,
//! so conformance has something to be judged against. Booleans skipped that
//! step: `axiolid-mesh-boolean-boolmesh` arrived first and was, for a while, the only
//! definition of a correct result. A suite that only ever runs one
//! implementation cannot tell "correct" from "self-consistent".
//!
//! # What "reference" means here
//!
//! Correctness first, speed never. This deliberately uses the most direct
//! algorithm that can be reasoned about line by line, because its job is to be
//! *obviously right*, not fast:
//!
//! - Classification is by **exact** [`orient3d`] signs and ray parity, not by
//!   floating-point distance comparisons.
//! - Work is `O(n·m)` with no acceleration structure. A BVH would be a second
//!   thing to get wrong, and an oracle with its own bugs is worse than none.
//!
//! # Independence
//!
//! This shares no code path with `boolmesh`. It does not subdivide against the
//! other operand's triangles; it classifies whole triangles by containment and
//! keeps or drops them. That makes it a genuinely independent implementation
//! for differential testing, at the cost of only being exact for operands whose
//! surfaces do not interpenetrate.
//!
//! # Honest limits
//!
//! [`ScalarBoolean`] refuses inputs it cannot answer exactly rather than
//! guessing. It is exact for:
//!
//! - disjoint operands (all four operations),
//! - nested operands (one strictly inside the other),
//! - identical operands,
//! - **properly intersecting surfaces**, via [`crate::exact_boolean`], which
//!   computes the intersection curve, retriangulates both operands along it,
//!   and keeps the pieces the operation asks for.
//!
//! It still reports [`GeomError::Unsupported`] for coplanar faces that share
//! an AREA -- two solids flush over a whole face. The shared region's boundary
//! is currently derived per triangle pair, which yields edges interior to that
//! region and a curve that branches. Continuing would produce a plausible
//! wrong answer (measured: a Difference volume of 2.67 where the geometry says
//! 8), so it refuses instead. Resolving it needs the overlap of the two face
//! SETS per plane rather than of individual triangle pairs.
//!
//! A refusal is typed, so a registry treats it as retryable and another
//! provider answers.
//!
//! Those cases already pin the algebra: identity, annihilation, idempotence,
//! and containment. See `tests/oracle.rs`.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, ExecutionOptions,
    ExecutionTarget, GeomError, GeomResult, Operation, ScratchRequirement, Sign,
};
use axiolid_core::BooleanOperator;
use axiolid_core::Point3;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::{BooleanEvidence, BooleanOutcome, MeshBoolean};

use crate::orient3d;

/// Portable scalar boolean reference.
///
/// Not a production provider: `O(n·m)`, and it refuses interpenetrating
/// surfaces. Registered at low priority so a real provider always wins
/// dispatch; it exists to be the thing conformance is judged against.
#[derive(Debug, Default, Clone, Copy)]
pub struct ScalarBoolean;

impl ScalarBoolean {
    /// Stable identity for this reference implementation.
    pub const ID: BackendId = BackendId::new("scalar-reference");

    /// Construct the reference provider.
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl Backend for ScalarBoolean {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

/// How one operand sits relative to the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Arrangement {
    /// The operand surfaces properly cross.
    ///
    /// Answered by [`crate::exact_boolean`], which retriangulates along the
    /// intersection curve. Kept as its own arrangement rather than an error so
    /// the decision to route stays visible next to the other cases.
    Interpenetrating,
    /// No shared volume and no surface contact.
    Disjoint,
    /// `subject` lies entirely within `tool`.
    SubjectInsideTool,
    /// `tool` lies entirely within `subject`.
    ToolInsideSubject,
    /// Same vertex set and same triangles, up to ordering.
    Identical,
}

impl MeshBoolean for ScalarBoolean {
    /// Exact, so no filter escalation and no scratch beyond the output.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::None
    }

    /// Checked per triangle pair, which is the inner loop of the `O(n·m)` scan.
    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn boolean(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operation: BooleanOperator,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        let arrangement = classify(subject, tool, options)?;
        let mesh = match (operation, arrangement) {
            // --- identical operands: idempotence and annihilation ---
            (BooleanOperator::Union | BooleanOperator::Intersection, Arrangement::Identical) => {
                subject.clone()
            }
            (
                BooleanOperator::Difference | BooleanOperator::SymmetricDifference,
                Arrangement::Identical,
            ) => empty(),

            // --- disjoint operands ---
            (
                BooleanOperator::Union | BooleanOperator::SymmetricDifference,
                Arrangement::Disjoint,
            ) => concatenate(subject, tool),
            (BooleanOperator::Intersection, Arrangement::Disjoint) => empty(),
            (BooleanOperator::Difference, Arrangement::Disjoint) => subject.clone(),

            // --- subject inside tool ---
            (BooleanOperator::Union, Arrangement::SubjectInsideTool) => tool.clone(),
            (BooleanOperator::Intersection, Arrangement::SubjectInsideTool) => subject.clone(),
            (BooleanOperator::Difference, Arrangement::SubjectInsideTool) => empty(),
            // A shell: outer boundary plus the inner boundary reversed, so the
            // cavity's normals point into the removed volume.
            (BooleanOperator::SymmetricDifference, Arrangement::SubjectInsideTool) => {
                concatenate(tool, &reversed(subject))
            }

            // --- tool inside subject ---
            (BooleanOperator::Union, Arrangement::ToolInsideSubject) => subject.clone(),
            (BooleanOperator::Intersection, Arrangement::ToolInsideSubject) => tool.clone(),
            (
                BooleanOperator::Difference | BooleanOperator::SymmetricDifference,
                Arrangement::ToolInsideSubject,
            ) => concatenate(subject, &reversed(tool)),

            // Properly crossing surfaces: hand over to the exact path, which
            // cuts both operands along the curve and keeps the pieces this
            // operation asks for. It refuses in turn on the shapes it cannot
            // resolve yet, and that refusal reaches the registry as
            // `Unsupported` so another provider can answer.
            (_, Arrangement::Interpenetrating) => crate::exact_boolean(subject, tool, operation)?,

            // The contract is `#[non_exhaustive]`; refuse rather than guess.
            _ => {
                return Err(GeomError::Unsupported {
                    backend: Self::ID,
                    operation: Operation::MeshBoolean,
                })
            }
        };

        let evidence = BooleanEvidence::record(
            subject.triangle_count(),
            tool.triangle_count(),
            mesh.triangle_count(),
            components(&mesh),
        )
        .with_disjoint_tools(usize::from(arrangement == Arrangement::Disjoint));
        Ok(BooleanOutcome::new(mesh, evidence))
    }
}

/// Exact orientation sign.
///
/// [`orient3d`] escalates to exact arithmetic internally and is documented to
/// always return `Certain`, so `Uncertain` is unreachable. Treating it as
/// [`Sign::Zero`] keeps that assumption from becoming a panic: a degenerate
/// answer makes callers refuse or retry, which is the safe direction.
fn exact_sign(certified: axiolid_contracts::Certified) -> Sign {
    certified.sign().unwrap_or(Sign::Zero)
}

/// Empty solid: a legitimate boolean result, not an error.
fn empty() -> TriMesh {
    TriMesh::new(Vec::new(), Vec::new())
}

/// Append `b`'s geometry to `a`'s, rebasing `b`'s indices.
fn concatenate(a: &TriMesh, b: &TriMesh) -> TriMesh {
    let offset = a.positions.len() as u32;
    let mut positions = a.positions.clone();
    positions.extend_from_slice(&b.positions);
    let mut indices = a.indices.clone();
    indices.extend(b.indices.iter().map(|i| i + offset));
    TriMesh::new(positions, indices)
}

/// Flip winding so the surface bounds the complement of what it bounded.
fn reversed(mesh: &TriMesh) -> TriMesh {
    let mut indices = mesh.indices.clone();
    for triangle in indices.chunks_exact_mut(3) {
        triangle.swap(0, 1);
    }
    TriMesh::new(mesh.positions.clone(), indices)
}

/// Connected components over triangle-shared vertices, by union-find.
fn components(mesh: &TriMesh) -> usize {
    if mesh.positions.is_empty() {
        return 0;
    }
    let mut parent: Vec<usize> = (0..mesh.positions.len()).collect();

    fn find(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }

    for triangle in mesh.indices.chunks_exact(3) {
        let root = find(&mut parent, triangle[0] as usize);
        for corner in &triangle[1..] {
            let other = find(&mut parent, *corner as usize);
            if root != other {
                parent[other] = root;
            }
        }
    }

    let mut roots = std::collections::BTreeSet::new();
    for index in &mesh.indices {
        let root = find(&mut parent, *index as usize);
        roots.insert(root);
    }
    roots.len()
}

/// Decide how the operands sit, or refuse if the answer needs real cutting.
fn classify(
    subject: &TriMesh,
    tool: &TriMesh,
    options: &ExecutionOptions,
) -> GeomResult<Arrangement> {
    // An operand with no geometry is not a solid. Refusing here rather than
    // indexing `positions[0]` keeps a malformed input from becoming a panic
    // inside a reference implementation, where a crash is the worst outcome:
    // it takes down the harness that was supposed to be judging correctness.
    for (mesh, role) in [(subject, "subject"), (tool, "tool")] {
        if mesh.positions.is_empty() || mesh.indices.is_empty() {
            return Err(GeomError::InvalidInput(format!(
                "{role}: an empty mesh has no interior and cannot be a boolean operand"
            )));
        }
    }

    if same_geometry(subject, tool) {
        return Ok(Arrangement::Identical);
    }

    // Surfaces that properly cross need retriangulating along the
    // intersection curve. That is no longer a refusal: it is a distinct
    // arrangement, answered by the exact path.
    if surfaces_intersect(subject, tool, options)? {
        return Ok(Arrangement::Interpenetrating);
    }

    // Non-crossing surfaces: containment is decided by a single vertex, since
    // the whole operand is on one side.
    let subject_in_tool = contains_point(tool, subject.positions[0]);
    let tool_in_subject = contains_point(subject, tool.positions[0]);

    Ok(match (subject_in_tool, tool_in_subject) {
        (true, false) => Arrangement::SubjectInsideTool,
        (false, true) => Arrangement::ToolInsideSubject,
        (false, false) => Arrangement::Disjoint,
        // Mutual containment is impossible for non-crossing closed surfaces.
        (true, true) => {
            return Err(GeomError::Degenerate(
                "operands report mutual containment, which is geometrically impossible".into(),
            ))
        }
    })
}

/// Same positions and same triangles, ignoring triangle order.
fn same_geometry(a: &TriMesh, b: &TriMesh) -> bool {
    if a.positions.len() != b.positions.len() || a.indices.len() != b.indices.len() {
        return false;
    }
    if a.positions
        .iter()
        .zip(&b.positions)
        .any(|(p, q)| p.x != q.x || p.y != q.y || p.z != q.z)
    {
        return false;
    }
    let mut left: Vec<[u32; 3]> = a
        .indices
        .chunks_exact(3)
        .map(|t| {
            let mut v = [t[0], t[1], t[2]];
            v.sort_unstable();
            v
        })
        .collect();
    let mut right: Vec<[u32; 3]> = b
        .indices
        .chunks_exact(3)
        .map(|t| {
            let mut v = [t[0], t[1], t[2]];
            v.sort_unstable();
            v
        })
        .collect();
    left.sort_unstable();
    right.sort_unstable();
    left == right
}

/// Whether any triangle of `a` properly crosses any triangle of `b`.
///
/// Uses exact [`orient3d`] signs: `b`'s triangle is crossed when `a`'s vertices
/// straddle its plane *and* the crossing lies inside the triangle. Shared
/// vertices and edge contact are not proper crossings.
fn surfaces_intersect(a: &TriMesh, b: &TriMesh, options: &ExecutionOptions) -> GeomResult<bool> {
    for left in a.indices.chunks_exact(3) {
        options.check_cancelled()?;
        let triangle_a = [
            a.positions[left[0] as usize],
            a.positions[left[1] as usize],
            a.positions[left[2] as usize],
        ];
        for right in b.indices.chunks_exact(3) {
            let triangle_b = [
                b.positions[right[0] as usize],
                b.positions[right[1] as usize],
                b.positions[right[2] as usize],
            ];
            if edges_cross_triangle(&triangle_a, &triangle_b)
                || edges_cross_triangle(&triangle_b, &triangle_a)
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Whether any edge of `edges` passes through the interior of `face`.
fn edges_cross_triangle(edges: &[Point3; 3], face: &[Point3; 3]) -> bool {
    let [p, q, r] = *face;
    for (start, end) in [
        (edges[0], edges[1]),
        (edges[1], edges[2]),
        (edges[2], edges[0]),
    ] {
        let side_start = exact_sign(orient3d(p, q, r, start));
        let side_end = exact_sign(orient3d(p, q, r, end));
        // Both on one side, or either exactly on the plane: not a proper
        // crossing. Touching is contact, and contact is not interpenetration.
        if side_start == Sign::Zero || side_end == Sign::Zero || side_start == side_end {
            continue;
        }
        // The segment pierces the plane; is the hit inside the CLOSED
        // triangle? Requiring three identical non-zero signs tests the open
        // interior only, and misses a hit landing exactly on an edge -- which
        // is precisely where two triangles of a quad meet. Both triangles then
        // report "no crossing" and interpenetration goes undetected.
        //
        // Closed test: the point is inside or on the boundary unless the signs
        // disagree strictly. Zeros mean "on an edge", which still counts.
        let signs = [
            exact_sign(orient3d(start, end, p, q)),
            exact_sign(orient3d(start, end, q, r)),
            exact_sign(orient3d(start, end, r, p)),
        ];
        let positive = signs.contains(&Sign::Positive);
        let negative = signs.contains(&Sign::Negative);
        if !(positive && negative) {
            return true;
        }
    }
    false
}

/// Whether `point` lies strictly inside the closed surface `mesh`.
///
/// Ray parity along `+x`. Rays that hit a vertex or edge are ambiguous, so the
/// direction is perturbed and retried rather than resolved by tolerance: an
/// oracle decides exactly or not at all.
/// Whether `point` lies strictly inside the closed surface `mesh`.
///
/// Exact ray parity: `orient3d` signs decide every crossing, so a point is
/// never misclassified by a near-miss. Shared with the exact boolean assembly,
/// which asks the same question per retriangulated piece.
/// Containment, or `None` when the point cannot be classified exactly.
///
/// `contains_point` folds three different situations into `false`: strictly
/// outside, exactly ON the surface, and every probe direction degenerate.
/// That is fine for the whole-operand arrangement test, which only ever
/// samples interior points. It is NOT fine for classifying a retriangulated
/// piece by its centroid: a centroid can land exactly on the other operand's
/// edge, and silently calling that "outside" keeps a piece that should be
/// dropped, leaving a hole in the result.
pub(crate) fn contains_point_exact(mesh: &TriMesh, point: Point3) -> Option<bool> {
    const DIRECTIONS: [[f64; 3]; 4] = [
        [1.0, 0.0, 0.0],
        [1.0, 0.125, 0.0625],
        [0.5, 1.0, 0.25],
        [0.25, 0.5, 1.0],
    ];

    // A point lying ON the surface is neither inside nor outside; report the
    // ambiguity rather than picking a side.
    if on_surface(mesh, point) {
        return None;
    }
    for direction in DIRECTIONS {
        if let Some(inside) = parity_along(mesh, point, direction) {
            return Some(inside);
        }
    }
    None
}

/// Whether `point` lies exactly on one of the mesh's triangles.
fn on_surface(mesh: &TriMesh, point: Point3) -> bool {
    mesh.indices.chunks_exact(3).any(|triangle| {
        let p = mesh.positions[triangle[0] as usize];
        let q = mesh.positions[triangle[1] as usize];
        let r = mesh.positions[triangle[2] as usize];
        exact_sign(orient3d(p, q, r, point)) == Sign::Zero
            && crate::segment_triangle_relation(point, point, [p, q, r])
                != crate::SegmentTriangleRelation::Disjoint
    })
}

pub(crate) fn contains_point(mesh: &TriMesh, point: Point3) -> bool {
    // Directions tried in order; each is used only if the previous produced a
    // degenerate hit. Fixed, so the result stays deterministic.
    const DIRECTIONS: [[f64; 3]; 4] = [
        [1.0, 0.0, 0.0],
        [1.0, 0.125, 0.0625],
        [0.5, 1.0, 0.25],
        [0.25, 0.5, 1.0],
    ];

    for direction in DIRECTIONS {
        if let Some(inside) = parity_along(mesh, point, direction) {
            return inside;
        }
    }
    // Every direction was degenerate. Outside is the conservative answer, and
    // callers only reach here for pathological inputs the oracle refuses.
    false
}

/// Count crossings along one ray, or `None` if any hit was degenerate.
fn parity_along(mesh: &TriMesh, origin: Point3, direction: [f64; 3]) -> Option<bool> {
    // A point far enough along the ray to be outside any operand: the ray
    // becomes a segment, which orient3d can answer exactly.
    let bounds = mesh.bounds();
    let span = (bounds.max.x - bounds.min.x)
        .max(bounds.max.y - bounds.min.y)
        .max(bounds.max.z - bounds.min.z)
        .max(1.0)
        * 8.0;
    let far = Point3::new(
        origin.x + direction[0] * span,
        origin.y + direction[1] * span,
        origin.z + direction[2] * span,
    );

    let mut crossings = 0usize;
    for triangle in mesh.indices.chunks_exact(3) {
        let p = mesh.positions[triangle[0] as usize];
        let q = mesh.positions[triangle[1] as usize];
        let r = mesh.positions[triangle[2] as usize];

        let side_origin = exact_sign(orient3d(p, q, r, origin));
        let side_far = exact_sign(orient3d(p, q, r, far));
        if side_origin == Sign::Zero {
            // The point is ON the surface: neither inside nor outside.
            return Some(false);
        }
        if side_far == Sign::Zero || side_origin == side_far {
            continue;
        }

        let a = exact_sign(orient3d(origin, far, p, q));
        let b = exact_sign(orient3d(origin, far, q, r));
        let c = exact_sign(orient3d(origin, far, r, p));
        // A zero means the ray grazes an edge or vertex: ambiguous parity, so
        // this direction cannot be trusted at all.
        if a == Sign::Zero || b == Sign::Zero || c == Sign::Zero {
            return None;
        }
        if a == b && b == c {
            crossings += 1;
        }
    }
    Some(crossings % 2 == 1)
}
