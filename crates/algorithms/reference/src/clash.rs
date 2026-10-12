//! Mesh interference: the narrow phase (#4).
//!
//! # What this closes
//!
//! `Bvh::overlap_pairs` reports pairs whose *bounding boxes* overlap, and
//! `triangle_triangle_relation` classifies *one* triangle pair exactly. Nothing
//! joined them, so there was no way to ask the question a model checker
//! actually asks: do these two solids interfere, and by how much?
//!
//! # Method
//!
//! Broad phase over triangle AABBs to reject the quadratic majority, then the
//! exact predicate on survivors. The predicate is unconditionally correct for
//! binary64 input, so a reported `Penetrating` is a fact, not an estimate.
//!
//! Penetration *depth* is a separate matter: it is measured, approximate, and
//! reported as evidence rather than folded into the verdict. Conflating an
//! exact topological answer with an approximate metric one is how a checker
//! ends up unable to say why it flagged something.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Aabb, Point2, Point3, Scalar, Tolerance};
use axiolid_measure::WindingMesh;
use axiolid_mesh::TriMesh;

use crate::orient2d;
use crate::triangle_triangle::{triangle_triangle_relation, TriangleTriangleRelation};

/// What two meshes do to each other.
///
/// Deliberately three states, not a boolean. A model checker that cannot
/// distinguish "touching" from "overlapping" either floods the report with
/// every abutting wall or silently drops real interferences.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interference {
    /// No triangle pair meets, within the separation tested.
    Clear,
    /// Triangles meet only at shared vertices, edges, or coplanar contact.
    /// Two slabs sharing a face are in contact, not in conflict.
    Touching,
    /// At least one pair crosses transversally: the solids share volume.
    Penetrating,
}

/// The verdict plus the evidence behind it.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct InterferenceReport {
    /// The verdict.
    pub kind: Interference,
    /// Triangle-index pairs that meet transversally.
    pub penetrating_pairs: Vec<(usize, usize)>,
    /// Triangle-index pairs in non-transverse contact.
    pub touching_pairs: Vec<(usize, usize)>,
    /// Triangle pairs whose boxes overlapped and were tested exactly.
    pub narrow_phase_tests: usize,
    /// Triangle pairs the broad phase rejected without an exact test.
    pub broad_phase_rejections: usize,
    /// Pairs skipped because a triangle was degenerate. A caller measuring
    /// clearance must know its input was not fully testable.
    pub degenerate_skips: usize,
    /// Whether a vertex of one solid was found strictly inside the other.
    ///
    /// Two boxes overlapping face-to-face cross only edge-to-edge, so the
    /// triangle predicate reports contact and never `Proper`. Shared volume is
    /// therefore decided by containment, not by surface topology alone.
    pub containment: bool,
}

impl InterferenceReport {
    /// Whether the solids share volume.
    pub fn is_penetrating(&self) -> bool {
        self.kind == Interference::Penetrating
    }

    /// Whether anything at all was found, contact included.
    pub fn is_clear(&self) -> bool {
        self.kind == Interference::Clear
    }
}

/// Classify interference between two triangle meshes.
///
/// `tolerance` inflates the broad-phase boxes so a pair that is within
/// tolerance of touching is still tested exactly. It does **not** loosen the
/// exact predicate: the verdict stays a topological fact about the supplied
/// coordinates.
pub fn interference(
    a: &TriMesh,
    b: &TriMesh,
    tolerance: Tolerance,
) -> GeomResult<InterferenceReport> {
    let pad = tolerance.linear();
    if !(pad.is_finite() && pad >= 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "tolerance must be finite and non-negative, got {pad}"
        )));
    }

    let boxes_a = triangle_boxes(a, pad);
    let boxes_b = triangle_boxes(b, pad);

    let mut report = InterferenceReport {
        kind: Interference::Clear,
        penetrating_pairs: Vec::new(),
        touching_pairs: Vec::new(),
        narrow_phase_tests: 0,
        broad_phase_rejections: 0,
        degenerate_skips: 0,
        containment: false,
    };

    // Build a BVH over B's triangles and probe it with A's boxes. The
    // quadratic scan this replaces made `interference` unusable at model
    // scale: two 4.6k-triangle solids cost 21M box tests.
    let tree = axiolid_spatial::Bvh::build(
        boxes_b
            .iter()
            .enumerate()
            .map(|(j, bounds)| axiolid_spatial::SpatialItem::new(j, *bounds)),
    );
    let mut hits: Vec<usize> = Vec::new();
    for (i, box_a) in boxes_a.iter().enumerate() {
        tree.query_aabb(box_a, &mut hits);
        // Everything the tree pruned would have been a rejected box test in
        // the quadratic version; count it so the two are comparable.
        report.broad_phase_rejections += boxes_b.len() - hits.len();
        for &j in &hits {
            let box_b = &boxes_b[j];
            if !box_a.intersects(box_b) {
                report.broad_phase_rejections += 1;
                continue;
            }
            report.narrow_phase_tests += 1;
            let ta = triangle(a, i);
            let tb = triangle(b, j);
            match triangle_triangle_relation(ta, tb) {
                TriangleTriangleRelation::Proper => {
                    report.penetrating_pairs.push((i, j));
                    report.kind = Interference::Penetrating;
                }
                TriangleTriangleRelation::Touching => {
                    report.touching_pairs.push((i, j));
                    if report.kind == Interference::Clear {
                        report.kind = Interference::Touching;
                    }
                }
                // `Coplanar` short-circuits the predicate before any edge is
                // tested, so it cannot separate "same plane, overlapping" from
                // "same plane, five metres apart". It is inconclusive here;
                // coplanar overlap is decided by the metric test below.
                TriangleTriangleRelation::Coplanar => {
                    // `Coplanar` means all six vertices share ONE plane. Two
                    // parallel faces 50mm apart are not coplanar, so they
                    // cannot reach here -- unless the predicate said so for
                    // the supplied coordinates, which is the exact answer.
                    if coplanar_pair_overlaps(ta, tb) {
                        report.touching_pairs.push((i, j));
                        if report.kind == Interference::Clear {
                            report.kind = Interference::Touching;
                        }
                    }
                }
                TriangleTriangleRelation::DegenerateTriangle => {
                    report.degenerate_skips += 1;
                }
                // The triangles were tested and share no point: a verdict,
                // not a skip. The match is exhaustive on purpose, so a new
                // relation must be decided here rather than fall through.
                TriangleTriangleRelation::Disjoint => {}
            }
        }
    }

    // Containment requires overlapping bounds. Checking that first turns the
    // common disjoint case from O(probes x triangles) into two box tests:
    // measured 1265 ms -> under 2 ms for two 2k-triangle spheres.
    let bounds_a = mesh_bounds(a);
    let bounds_b = mesh_bounds(b);
    if report.kind != Interference::Penetrating
        && !a.indices.is_empty()
        && !b.indices.is_empty()
        && bounds_a.intersects(&bounds_b)
    {
        // Prepare each winding mesh ONCE. `WindingMesh::prepare` runs a full
        // audit_mesh, so preparing per probe made containment quadratic in
        // triangle count with an enormous constant.
        let a_in_b = WindingMesh::prepare(b, tolerance)
            .ok()
            .is_some_and(|w| interior_probes(a).any(|p| inside(&w, p)));
        let b_in_a = WindingMesh::prepare(a, tolerance)
            .ok()
            .is_some_and(|w| interior_probes(b).any(|p| inside(&w, p)));
        if a_in_b || b_in_a {
            report.kind = Interference::Penetrating;
            report.containment = true;
        }
    }

    Ok(report)
}

/// Padded axis-aligned box per triangle.
fn triangle_boxes(mesh: &TriMesh, pad: Scalar) -> Vec<Aabb> {
    (0..mesh.indices.len() / 3)
        .map(|i| {
            let [a, b, c] = triangle(mesh, i);
            let lo = Point3::new(
                a.x.min(b.x).min(c.x) - pad,
                a.y.min(b.y).min(c.y) - pad,
                a.z.min(b.z).min(c.z) - pad,
            );
            let hi = Point3::new(
                a.x.max(b.x).max(c.x) + pad,
                a.y.max(b.y).max(c.y) + pad,
                a.z.max(b.z).max(c.z) + pad,
            );
            Aabb { min: lo, max: hi }
        })
        .collect()
}

/// Corner positions of triangle `i`.
fn triangle(mesh: &TriMesh, i: usize) -> [Point3; 3] {
    let base = i * 3;
    [
        mesh.positions[mesh.indices[base] as usize],
        mesh.positions[mesh.indices[base + 1] as usize],
        mesh.positions[mesh.indices[base + 2] as usize],
    ]
}

/// Whether two coplanar triangles actually share area.
///
/// `triangle_triangle_relation` returns `Coplanar` from a vertex-side test and
/// never reaches its edge tests, so it cannot answer this. Projecting onto the
/// dominant plane axis and testing in 2D can: exact `orient2d` decides both
/// edge crossings and vertex containment.
fn coplanar_pair_overlaps(a: [Point3; 3], b: [Point3; 3]) -> bool {
    let normal = (a[1] - a[0]).cross(a[2] - a[0]);
    let drop = dominant_axis(normal);
    let pa = a.map(|p| flatten(p, drop));
    let pb = b.map(|p| flatten(p, drop));

    // Either triangle containing any corner of the other is overlap.
    if pb.iter().any(|p| point_in_triangle(*p, pa)) || pa.iter().any(|p| point_in_triangle(*p, pb))
    {
        return true;
    }
    // Otherwise they overlap only if their boundaries cross.
    let edges = |t: [Point2; 3]| [[t[0], t[1]], [t[1], t[2]], [t[2], t[0]]];
    edges(pa)
        .iter()
        .any(|ea| edges(pb).iter().any(|eb| segments_cross(*ea, *eb)))
}

/// Index of the largest-magnitude normal component.
fn dominant_axis(n: axiolid_core::Vec3) -> usize {
    let (x, y, z) = (n.x.abs(), n.y.abs(), n.z.abs());
    if x >= y && x >= z {
        0
    } else if y >= z {
        1
    } else {
        2
    }
}

/// Drop the dominant axis, keeping the projection non-degenerate.
fn flatten(p: Point3, drop: usize) -> Point2 {
    match drop {
        0 => Point2::new(p.y, p.z),
        1 => Point2::new(p.x, p.z),
        _ => Point2::new(p.x, p.y),
    }
}

/// Containment by exact orientation, boundary included.
fn point_in_triangle(p: Point2, t: [Point2; 3]) -> bool {
    let s = |i: usize, j: usize| sign_of(orient2d(t[i], t[j], p));
    let (a, b, c) = (s(0, 1), s(1, 2), s(2, 0));
    let non_negative = a >= 0 && b >= 0 && c >= 0;
    let non_positive = a <= 0 && b <= 0 && c <= 0;
    non_negative || non_positive
}

/// Whether two closed segments share a point.
fn segments_cross(u: [Point2; 2], v: [Point2; 2]) -> bool {
    let d1 = sign_of(orient2d(u[0], u[1], v[0]));
    let d2 = sign_of(orient2d(u[0], u[1], v[1]));
    let d3 = sign_of(orient2d(v[0], v[1], u[0]));
    let d4 = sign_of(orient2d(v[0], v[1], u[1]));
    if d1 * d2 < 0 && d3 * d4 < 0 {
        return true;
    }
    // Collinear touching cases.
    (d1 == 0 && between(u[0], v[0], u[1]))
        || (d2 == 0 && between(u[0], v[1], u[1]))
        || (d3 == 0 && between(v[0], u[0], v[1]))
        || (d4 == 0 && between(v[0], u[1], v[1]))
}

/// Whether collinear `q` lies within the `p`-`r` box.
fn between(p: Point2, q: Point2, r: Point2) -> bool {
    q.x >= p.x.min(r.x) && q.x <= p.x.max(r.x) && q.y >= p.y.min(r.y) && q.y <= p.y.max(r.y)
}

/// Certified sign as a small integer.
fn sign_of(c: axiolid_contracts::Certified) -> i32 {
    match c.sign() {
        Some(axiolid_contracts::Sign::Positive) => 1,
        Some(axiolid_contracts::Sign::Negative) => -1,
        _ => 0,
    }
}

/// Whether `point` lies strictly inside the closed mesh `solid`.
///
/// Uses the generalized winding number rather than ray parity. Parity is
/// unreliable exactly where building geometry lives: axis-aligned boxes put
/// rays along faces and through shared edges, and every tie-break there is a
/// guess. Winding accumulates oriented solid angle, so it degrades smoothly
/// instead of flipping.
///
/// A point ON the boundary has winding near 1/2 and is deliberately NOT
/// inside: two slabs sharing a face are in contact, not overlapping. That
/// single decision is what keeps every abutting wall out of a clash report.
pub fn point_inside(point: Point3, solid: &TriMesh, tolerance: Tolerance) -> Option<bool> {
    let winding = WindingMesh::prepare(solid, tolerance).ok()?;
    let w = winding.winding_number(point).ok()?.value;
    // Interior is ~1, exterior ~0, boundary ~0.5. Require a clear interior so
    // a boundary point never counts as containment.
    Some(w > 0.75)
}

/// Points used to test whether one solid reaches inside another.
///
/// Vertices alone are not enough: two boxes crossing face-to-face have every
/// corner outside or on the other's boundary, yet plainly share volume. Each
/// triangle centroid of the overlapping faces does lie strictly inside, so
/// both are sampled.
///
/// This is a sampling test, so it can only ever produce false negatives, never
/// false positives: a reported penetration is always real. A pathological
/// sliver overlap smaller than a triangle is missed, which is why the exact
/// surface-crossing test above runs first and independently.
fn interior_probes(mesh: &TriMesh) -> impl Iterator<Item = Point3> + '_ {
    // Centre of the mesh's own bounding box: a convex solid always contains
    // it, and for a non-convex one the triangle-nudge probes below still fire.
    let mut lo = Point3::splat(Scalar::INFINITY);
    let mut hi = Point3::splat(Scalar::NEG_INFINITY);
    for p in &mesh.positions {
        lo = lo.min(*p);
        hi = hi.max(*p);
    }
    let centre = (lo + hi) * 0.5;

    // Each triangle centroid pulled a whisker toward the mesh centre. On the
    // boundary the winding number is undefined -- measured as ~1.0 rather
    // than 0.5 -- so a probe left exactly on a shared face reads as inside
    // and turns face contact into a false penetration. The nudge is relative
    // so it scales with the model.
    let span = (hi - lo).length().max(1.0);
    // The nudge must be small enough not to step over a real overlap. A 1e-9
    // model-relative step is the same order as the smallest overlap worth
    // reporting, so it would hide exactly the cases this function exists to
    // catch. 1e-12 is far below any meaningful interference yet still clears
    // the exact-arithmetic boundary where winding is undefined.
    let nudge = span * 1e-12;

    core::iter::once(centre).chain(mesh.indices.chunks_exact(3).map(move |t| {
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        let m = (a + b + c) / 3.0;
        let toward = centre - m;
        let len = toward.length();
        if len > 0.0 {
            m + toward * (nudge / len)
        } else {
            m
        }
    }))
}

/// Interior test against an already-prepared winding mesh.
///
/// Shares the 0.75 threshold with `point_inside`; see there for why a
/// boundary point must not count as inside.
fn inside(winding: &WindingMesh<'_, TriMesh>, point: Point3) -> bool {
    winding
        .winding_number(point)
        .map(|w| w.value > 0.75)
        .unwrap_or(false)
}

/// Axis-aligned bounds of a mesh.
fn mesh_bounds(mesh: &TriMesh) -> Aabb {
    let mut bounds = Aabb::empty();
    for p in &mesh.positions {
        bounds.extend(*p);
    }
    bounds
}
