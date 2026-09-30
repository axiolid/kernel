//! Certified Hausdorff distance between two triangle meshes (#148, D15).
//!
//! # What is certified
//!
//! [`one_sided_hausdorff`] returns an interval `[lower, upper]` certain to
//! contain `h(A, B) = max over a in A of min over b in B of |a - b|`, the
//! one-sided Hausdorff distance between the surfaces of two meshes, and
//! [`hausdorff_distance`] the two-sided `max(h(A, B), h(B, A))` together
//! with both one-sided intervals. Both ends are guarantees:
//!
//! - `lower` is `d(y, B)` for a sampled point `y` of `A`, bounded below.
//!   `d(y, B)` is the minimum over the triangles `T` of `B` of `d(y, T)`, and
//!   each `d(y, T)` is bounded below through the support function of the
//!   convex `T`: for every direction `u` with `|u| <= 1`,
//!   `d(y, T) >= min over corners c of T of u . (y - c)`. The direction is
//!   taken towards the computed nearest point, where the bound is tight, and
//!   triangles the search never opens are bounded by their boxes.
//! - `upper` bounds every point of `A`. `A`'s triangles are subdivided, and
//!   over one piece `t`, `max over x in t of d(x, B) <= min over T of max
//!   over corners v of t of d(v, T)`: `d(., T)` is convex, so its maximum
//!   over a triangle is at a corner. Each `d(v, T)` is bounded above by the
//!   distance to a point of `T` itself.
//!
//! The witnesses are `point_from`, the sample realising `lower`, and
//! `point_to`, its nearest point found on the other mesh.
//!
//! # Rounding accounted for
//!
//! The input coordinates are taken as exact. Every bound is computed in
//! binary64 and widened by an explicit margin in the standard model
//! (`fl(a op b) = (a op b)(1 + d)`, `|d| <= eps / 2`, `eps = 2^-52`), with
//! no underflow assumed:
//!
//! - A point of `T` is `a + s (b - a) + t (c - a)` with `s, t >= 0` and
//!   `s + t <= 1` up to three units in the last place. Its rounded
//!   evaluation, the parameter overshoot and the rounded length `|v - p|`
//!   together stay below `8.5 eps S`, with `S = |v| + |a| + |b| + |c|`; the
//!   upper bound adds `16 eps S`.
//! - The support bound normalises `u` and shrinks it by `4 eps` so `|u| <= 1`
//!   survives rounding; each `u . (v - c)` is off by at most
//!   `2.5 eps (|v| + |c|)`, and the lower bound subtracts `4 eps (|v| + max |c|)`.
//! - A box gap is shrunk by `4 eps` of itself: coordinate differences keep
//!   their sign and the rounded Euclidean length is within `2 eps`.
//! - Subdivision points are rounded midpoints, so a piece is only close to
//!   the exact piece of `A` it stands for. Each corner carries its drift, the
//!   parent's plus `eps |m|` per midpoint; the distance is 1-Lipschitz, so a
//!   piece's upper bound adds its largest drift and a sample's lower bound
//!   subtracts its own. `point_from` is therefore within that drift (a few
//!   units in the last place of its coordinates) of `A`, not on it.
//! - A flat patch of `B` (below) carries a certified excess, rounded up.
//!
//! No other rounding enters: branch decisions (which Voronoi region, which
//! edge to split, which triangles to try) only choose a point, a direction
//! or a candidate, and any choice keeps the bounds sound; they affect how
//! fast the interval closes, not whether it contains the distance.
//!
//! # Method
//!
//! Branch and bound over `A`'s triangles. Pieces are kept in a queue by
//! upper bound; the largest is split at the midpoint of its longest edge,
//! the midpoint is measured to raise `lower`, and pieces whose bound falls
//! to `lower` are dropped. Nearest queries walk a bounding-volume hierarchy
//! over `B` with a node bound of `max over corners of the box gap`, which is
//! a certified lower bound of what any triangle below the node can offer,
//! seeded with the triangles found nearest to the piece's corners.
//!
//! One triangle at a time, a piece over the seam between two triangles of a
//! flat face stays loose until it is split below the accuracy, along every
//! seam. So flat convex patches of `B` -- two triangles forming a convex
//! quadrilateral, the fan around an interior vertex with a convex link --
//! bound a piece as one convex set, up to a certified excess (see
//! `patches`).
//!
//! The hierarchy is private rather than `axiolid-spatial`'s: that one ranks
//! candidates by box-to-box gap and returns the nearest box, where this
//! search needs a best-first walk under a caller's bound over several points
//! at once, and keeps `axiolid-measure` free of a dependency for it.
//!
//! # Convergence
//!
//! A piece nearest to a single triangle or flat patch of `B` is bounded
//! exactly up to rounding, so only pieces across the medial surfaces of `B`
//! need splitting. Two things end the refinement early, and the interval is
//! sound either way, only wider than asked: the split budget, and a request
//! for an accuracy below the rounding margins, where pieces stop splitting
//! once their edges are within rounding of their coordinates.
//!
//! # Scope
//!
//! This is the distance between SURFACES, like [`crate::mesh_distance`]: a
//! mesh wholly inside a closed mesh is as far from it as from its boundary.
//! Degenerate triangles are measured as their edges. Open, non-manifold and
//! self-intersecting meshes are all fine; the Hausdorff distance is defined
//! for any point sets.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use axiolid_core::Point3;
use axiolid_mesh::TriangleMeshView;

mod patches;

use patches::{flat_patches, Patch};

const EPS: f64 = f64::EPSILON;

/// Pieces split per one-sided query before it reports what it has.
const MAX_SPLITS: usize = 400_000;

/// Triangles per hierarchy leaf.
const LEAF_SIZE: usize = 4;

/// Why a Hausdorff query could not be answered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum HausdorffError {
    /// A mesh had no triangles: the distance to or from nothing is undefined.
    EmptyMesh,
    /// A triangle referenced a position the mesh does not have.
    IndexOutOfRange,
    /// A position used by a triangle is NaN or infinite.
    NonFiniteInput,
    /// The accuracy was negative or NaN.
    InvalidAccuracy,
}

impl core::fmt::Display for HausdorffError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str(match self {
            Self::EmptyMesh => "mesh has no triangles",
            Self::IndexOutOfRange => "triangle index out of range",
            Self::NonFiniteInput => "mesh position is not finite",
            Self::InvalidAccuracy => "accuracy must be non-negative",
        })
    }
}

impl std::error::Error for HausdorffError {}

/// An interval certain to contain a Hausdorff distance, with its witness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HausdorffBounds {
    /// The distance is at least this: `point_from` is certainly this far
    /// from every point of the other mesh (less its drift, see the module
    /// documentation).
    pub lower: f64,
    /// The distance is at most this: every point of the measured mesh is
    /// within it of the other mesh.
    pub upper: f64,
    /// The sample realising `lower`, on the mesh measured FROM (within
    /// rounding of it).
    pub point_from: Point3,
    /// The nearest point found to `point_from` on the other mesh.
    pub point_to: Point3,
}

impl HausdorffBounds {
    /// Width of the interval.
    #[must_use]
    pub fn width(&self) -> f64 {
        self.upper - self.lower
    }

    /// The interval contains `value`.
    #[must_use]
    pub fn contains(&self, value: f64) -> bool {
        self.lower <= value && value <= self.upper
    }
}

/// The two-sided Hausdorff distance and both one-sided ones.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshHausdorff {
    /// `max(h(A, B), h(B, A))`. The witness is that of the side with the
    /// larger lower bound; `point_from` lies on that side's mesh.
    pub distance: HausdorffBounds,
    /// `h(A, B)`: how far `A` strays from `B`. Witness from `A` to `B`.
    pub forward: HausdorffBounds,
    /// `h(B, A)`: how far `B` strays from `A`. Witness from `B` to `A`.
    pub backward: HausdorffBounds,
}

/// One-sided Hausdorff distance `h(from, to)`, to within `accuracy`.
///
/// Refines until `upper - lower <= accuracy` or the split budget runs out;
/// the interval is sound either way.
///
/// # Errors
///
/// An empty mesh, an out-of-range index, a non-finite position used by a
/// triangle, or a negative or NaN accuracy.
pub fn one_sided_hausdorff<A, B>(
    from: &A,
    to: &B,
    accuracy: f64,
) -> Result<HausdorffBounds, HausdorffError>
where
    A: TriangleMeshView + ?Sized,
    B: TriangleMeshView + ?Sized,
{
    let accuracy = check_accuracy(accuracy)?;
    let source = collect(from)?;
    let target = Target::new(collect(to)?, accuracy);
    Ok(search(&source, &target, accuracy))
}

/// Two-sided Hausdorff distance between `a` and `b`, to within `accuracy`.
///
/// Both one-sided distances are refined to `accuracy`, so each of the three
/// intervals returned is at most that wide unless the budget ran out.
///
/// # Errors
///
/// As [`one_sided_hausdorff`].
pub fn hausdorff_distance<A, B>(
    a: &A,
    b: &B,
    accuracy: f64,
) -> Result<MeshHausdorff, HausdorffError>
where
    A: TriangleMeshView + ?Sized,
    B: TriangleMeshView + ?Sized,
{
    let accuracy = check_accuracy(accuracy)?;
    let first = collect(a)?;
    let second = collect(b)?;
    let forward = search(&first, &Target::new(second.clone(), accuracy), accuracy);
    let backward = search(&second, &Target::new(first, accuracy), accuracy);
    let mut distance = if backward.lower > forward.lower {
        backward
    } else {
        forward
    };
    distance.upper = forward.upper.max(backward.upper);
    Ok(MeshHausdorff {
        distance,
        forward,
        backward,
    })
}

fn check_accuracy(accuracy: f64) -> Result<f64, HausdorffError> {
    if accuracy.is_nan() || accuracy < 0.0 {
        Err(HausdorffError::InvalidAccuracy)
    } else {
        Ok(accuracy)
    }
}

/// A mesh's triangles, as points and as position indices.
#[derive(Clone)]
struct Surface {
    triangles: Vec<[Point3; 3]>,
    indices: Vec<[usize; 3]>,
    positions: Vec<Point3>,
    /// Positions some triangle uses.
    used: Vec<bool>,
}

fn collect<M: TriangleMeshView + ?Sized>(mesh: &M) -> Result<Surface, HausdorffError> {
    let positions = mesh.position_count();
    let mut used = vec![false; positions];
    let mut triangles = Vec::with_capacity(mesh.triangle_count());
    let mut indices = Vec::with_capacity(mesh.triangle_count());
    for index in 0..mesh.triangle_count() {
        let mut corners = [Point3::ZERO; 3];
        let mut ids = [0; 3];
        for (slot, corner) in mesh.triangle(index).iter().enumerate() {
            let corner = usize::try_from(*corner)
                .ok()
                .filter(|corner| *corner < positions)
                .ok_or(HausdorffError::IndexOutOfRange)?;
            let point = mesh.position(corner);
            if !point.is_finite() {
                return Err(HausdorffError::NonFiniteInput);
            }
            corners[slot] = point;
            ids[slot] = corner;
            used[corner] = true;
        }
        triangles.push(corners);
        indices.push(ids);
    }
    if triangles.is_empty() {
        return Err(HausdorffError::EmptyMesh);
    }
    let positions = (0..positions).map(|index| mesh.position(index)).collect();
    Ok(Surface {
        triangles,
        indices,
        positions,
        used,
    })
}

/// A triangle of the mesh measured TO, with its rounding scales.
struct Facet {
    corners: [Point3; 3],
    /// `|a| + |b| + |c|`.
    scale: f64,
    /// `max |corner|`.
    reach: f64,
    /// Too thin to project onto: measured as its three edges.
    thin: bool,
    lo: Point3,
    hi: Point3,
}

enum NodeKind {
    Leaf { start: usize, end: usize },
    Branch { left: usize, right: usize },
}

struct Node {
    lo: Point3,
    hi: Point3,
    kind: NodeKind,
}

/// The mesh measured TO, with a median-split hierarchy over its triangles.
struct Target {
    facets: Vec<Facet>,
    /// Facet indices in leaf order.
    order: Vec<usize>,
    nodes: Vec<Node>,
    /// Flat convex patches (see [`patches`]).
    patches: Vec<Patch>,
    /// Per facet, the patches it belongs to.
    patches_of: Vec<Vec<usize>>,
    /// Per patch, the box around its members.
    patch_boxes: Vec<(Point3, Point3)>,
}

impl Target {
    /// Pairs whose hull excess is within a quarter of `accuracy`, or within
    /// rounding of the model, are used.
    fn new(surface: Surface, accuracy: f64) -> Self {
        let facets: Vec<Facet> = surface
            .triangles
            .into_iter()
            .map(|corners| {
                let [a, b, c] = corners;
                let (ab, ac) = (b - a, c - a);
                let area = ab.cross(ac).length_squared();
                Facet {
                    corners,
                    scale: a.length() + b.length() + c.length(),
                    reach: a.length().max(b.length()).max(c.length()),
                    thin: area <= 1e-20 * ab.length_squared() * ac.length_squared(),
                    lo: a.min(b).min(c),
                    hi: a.max(b).max(c),
                }
            })
            .collect();
        let size = facets.iter().map(|f| f.reach).fold(0.0, f64::max);
        let corners: Vec<[Point3; 3]> = facets.iter().map(|f| f.corners).collect();
        let thin: Vec<bool> = facets.iter().map(|f| f.thin).collect();
        let (patches, patches_of) =
            flat_patches(&corners, &thin, (accuracy / 4.0).max(1024.0 * EPS * size));
        let patch_boxes = patches
            .iter()
            .map(|patch| {
                patch.members.iter().fold(
                    (
                        Point3::splat(f64::INFINITY),
                        Point3::splat(f64::NEG_INFINITY),
                    ),
                    |(lo, hi), &m| (lo.min(facets[m].lo), hi.max(facets[m].hi)),
                )
            })
            .collect();
        let mut target = Self {
            order: (0..facets.len()).collect(),
            facets,
            nodes: Vec::new(),
            patches,
            patches_of,
            patch_boxes,
        };
        let count = target.order.len();
        target.build(0, count);
        target
    }

    fn build(&mut self, start: usize, end: usize) -> usize {
        let (mut lo, mut hi) = (
            Point3::splat(f64::INFINITY),
            Point3::splat(f64::NEG_INFINITY),
        );
        let (mut clo, mut chi) = (lo, hi);
        for &index in &self.order[start..end] {
            for corner in self.facets[index].corners {
                lo = lo.min(corner);
                hi = hi.max(corner);
            }
            let centre = centroid(&self.facets[index].corners);
            clo = clo.min(centre);
            chi = chi.max(centre);
        }
        let node = self.nodes.len();
        self.nodes.push(Node {
            lo,
            hi,
            kind: NodeKind::Leaf { start, end },
        });
        if end - start <= LEAF_SIZE {
            return node;
        }
        let extent = chi - clo;
        let axis = if extent.x >= extent.y && extent.x >= extent.z {
            0
        } else if extent.y >= extent.z {
            1
        } else {
            2
        };
        let middle = start + (end - start) / 2;
        let facets = &self.facets;
        self.order[start..end].select_nth_unstable_by(middle - start, |&i, &j| {
            centroid(&facets[i].corners)[axis].total_cmp(&centroid(&facets[j].corners)[axis])
        });
        let left = self.build(start, middle);
        let right = self.build(middle, end);
        self.nodes[node].kind = NodeKind::Branch { left, right };
        node
    }

    /// Certified bounds on `d(v, self)`, and the nearest point found with
    /// its facet.
    fn point(&self, v: Point3) -> Nearest {
        let mut best = f64::INFINITY;
        let mut nearest = (self.facets[0].corners[0], 0);
        let mut lower = f64::INFINITY;
        let root = &self.nodes[0];
        let mut stack = vec![(0usize, box_gap(v, root.lo, root.hi))];
        while let Some((node, gap)) = stack.pop() {
            if gap >= best {
                // Nothing below the node is nearer than `gap`.
                lower = lower.min(gap);
                continue;
            }
            match self.nodes[node].kind {
                NodeKind::Leaf { start, end } => {
                    for &index in &self.order[start..end] {
                        let facet = &self.facets[index];
                        let gap = box_gap(v, facet.lo, facet.hi);
                        if gap >= best {
                            lower = lower.min(gap);
                            continue;
                        }
                        let p = nearest_on(v, facet);
                        let up = upper_at(v, p, facet);
                        lower = lower.min(lower_at(v, p, facet));
                        if up < best {
                            best = up;
                            nearest = (p, index);
                        }
                    }
                }
                NodeKind::Branch { left, right } => {
                    self.push_children(&mut stack, left, right, |n| box_gap(v, n.lo, n.hi));
                }
            }
        }
        Nearest {
            lower: lower.max(0.0).min(best),
            point: nearest.0,
            facet: nearest.1,
        }
    }

    /// A certified upper bound on `max over x in piece of d(x, self)`,
    /// seeded with the facets nearest to the corners.
    fn piece(&self, corners: &[Point3; 3], hints: [usize; 3]) -> f64 {
        let gap = |lo: Point3, hi: Point3| {
            corners
                .iter()
                .map(|&v| box_gap(v, lo, hi))
                .fold(0.0, f64::max)
        };
        let mut best = hints
            .iter()
            .map(|&hint| self.worst(corners, &[hint]))
            .fold(f64::INFINITY, f64::min);
        let mut seen = Vec::new();
        let root = &self.nodes[0];
        let mut stack = vec![(0usize, gap(root.lo, root.hi))];
        while let Some((node, lower)) = stack.pop() {
            if lower >= best {
                continue;
            }
            match self.nodes[node].kind {
                NodeKind::Leaf { start, end } => {
                    for &index in &self.order[start..end] {
                        let facet = &self.facets[index];
                        if gap(facet.lo, facet.hi) < best {
                            best = best.min(self.worst(corners, &[index]));
                        }
                        // A flat convex patch bounds the piece by its hull,
                        // whose distance is convex across the seams.
                        for &id in &self.patches_of[index] {
                            let (lo, hi) = self.patch_boxes[id];
                            if seen.contains(&id) || gap(lo, hi) >= best {
                                continue;
                            }
                            seen.push(id);
                            let patch = &self.patches[id];
                            best = best.min(self.worst(corners, &patch.members) + patch.excess);
                        }
                    }
                }
                NodeKind::Branch { left, right } => {
                    self.push_children(&mut stack, left, right, |n| gap(n.lo, n.hi));
                }
            }
        }
        best
    }

    /// `max over corners v of min over members T of d(v, T)`, bounded above.
    fn worst(&self, corners: &[Point3; 3], members: &[usize]) -> f64 {
        corners
            .iter()
            .map(|&v| {
                members
                    .iter()
                    .map(|&member| {
                        let facet = &self.facets[member];
                        upper_at(v, nearest_on(v, facet), facet)
                    })
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max)
    }

    /// Push both children, the nearer last so it is walked first.
    fn push_children(
        &self,
        stack: &mut Vec<(usize, f64)>,
        left: usize,
        right: usize,
        gap: impl Fn(&Node) -> f64,
    ) {
        let (gl, gr) = (gap(&self.nodes[left]), gap(&self.nodes[right]));
        if gl <= gr {
            stack.push((right, gr));
            stack.push((left, gl));
        } else {
            stack.push((left, gl));
            stack.push((right, gr));
        }
    }
}

fn centroid(corners: &[Point3; 3]) -> Point3 {
    (corners[0] + corners[1] + corners[2]) / 3.0
}

/// What a nearest-point query found.
struct Nearest {
    /// A certified lower bound on the distance.
    lower: f64,
    /// The nearest point found, on the facet below.
    point: Point3,
    facet: usize,
}

/// A certified lower bound on the distance from `v` to the box `[lo, hi]`.
fn box_gap(v: Point3, lo: Point3, hi: Point3) -> f64 {
    let gap = (lo - v).max(v - hi).max(Point3::ZERO);
    gap.length() * (1.0 - 4.0 * EPS)
}

/// A certified upper bound on `d(v, facet)` from `p`, a computed point of it.
fn upper_at(v: Point3, p: Point3, facet: &Facet) -> f64 {
    (v - p).length() + 16.0 * EPS * (v.length() + facet.scale)
}

/// A certified lower bound on `d(v, facet)`, from the support function in
/// the direction of `v - p`.
fn lower_at(v: Point3, p: Point3, facet: &Facet) -> f64 {
    let w = v - p;
    let length = w.length();
    if length <= 0.0 {
        return 0.0;
    }
    let u = w * ((1.0 - 4.0 * EPS) / length);
    let support = facet
        .corners
        .iter()
        .map(|&corner| u.dot(v - corner))
        .fold(f64::INFINITY, f64::min);
    (support - 4.0 * EPS * (v.length() + facet.reach)).max(0.0)
}

/// A point of `facet` near `v`, as `a + s (b - a) + t (c - a)` with
/// `s, t >= 0` and `s + t <= 1` (up to rounding the bounds account for).
fn nearest_on(v: Point3, facet: &Facet) -> Point3 {
    let [a, b, c] = facet.corners;
    let (ab, ac) = (b - a, c - a);
    let (s, t) = if facet.thin {
        thin_params(v, a, b, c)
    } else {
        let (s, t) = region_params(v, a, b, c);
        if s.is_finite() && t.is_finite() {
            (s, t)
        } else {
            thin_params(v, a, b, c)
        }
    };
    let (mut s, mut t) = (s.clamp(0.0, 1.0), t.clamp(0.0, 1.0));
    let sum = s + t;
    if sum > 1.0 {
        s /= sum;
        t /= sum;
    }
    a + ab * s + ac * t
}

/// Parameters of the nearest point by Voronoi region (Ericson, RTCD 5.1.5).
fn region_params(v: Point3, a: Point3, b: Point3, c: Point3) -> (f64, f64) {
    let (ab, ac) = (b - a, c - a);
    let ap = v - a;
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return (0.0, 0.0);
    }
    let bp = v - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return (1.0, 0.0);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return (d1 / (d1 - d3), 0.0);
    }
    let cp = v - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return (0.0, 1.0);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return (0.0, d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (1.0 - w, w);
    }
    let denominator = va + vb + vc;
    (vb / denominator, vc / denominator)
}

/// Parameters of the nearest point on the three edges.
fn thin_params(v: Point3, a: Point3, b: Point3, c: Point3) -> (f64, f64) {
    let along = |from: Point3, to: Point3| {
        let edge = to - from;
        let length = edge.length_squared();
        if length > 0.0 {
            ((v - from).dot(edge) / length).clamp(0.0, 1.0)
        } else {
            0.0
        }
    };
    let (x, y, z) = (along(a, b), along(a, c), along(b, c));
    let candidates = [(x, 0.0), (0.0, y), (1.0 - z, z)];
    let point = |(s, t): (f64, f64)| a + (b - a) * s + (c - a) * t;
    candidates
        .into_iter()
        .min_by(|&p, &q| (v - point(p)).length().total_cmp(&(v - point(q)).length()))
        .unwrap_or((0.0, 0.0))
}

/// A piece of the measured mesh: a triangle with per-corner drift and the
/// facets found nearest to its corners.
struct Piece {
    corners: [Point3; 3],
    drift: [f64; 3],
    hints: [usize; 3],
    upper: f64,
}

impl PartialEq for Piece {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Piece {}

impl PartialOrd for Piece {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Piece {
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper.total_cmp(&other.upper)
    }
}

fn search(source: &Surface, target: &Target, accuracy: f64) -> HausdorffBounds {
    let first = source.triangles[0][0];
    let mut best = HausdorffBounds {
        lower: 0.0,
        upper: f64::INFINITY,
        point_from: first,
        point_to: target.point(first).point,
    };
    // Measure a sample `drift` from the mesh; returns its nearest facet.
    let sample = |v: Point3, drift: f64, best: &mut HausdorffBounds| {
        let nearest = target.point(v);
        let lower = nearest.lower - drift;
        if lower > best.lower {
            best.lower = lower;
            best.point_from = v;
            best.point_to = nearest.point;
        }
        nearest.facet
    };
    let hints: Vec<usize> = (0..source.positions.len())
        .map(|index| {
            if source.used[index] {
                sample(source.positions[index], 0.0, &mut best)
            } else {
                0
            }
        })
        .collect();

    let make = |corners: [Point3; 3], drift: [f64; 3], hints: [usize; 3]| {
        let upper = target.piece(&corners, hints) + drift.iter().copied().fold(0.0, f64::max);
        Piece {
            corners,
            drift,
            hints,
            upper,
        }
    };
    let mut queue: BinaryHeap<Piece> = source
        .triangles
        .iter()
        .zip(&source.indices)
        .map(|(&corners, ids)| make(corners, [0.0; 3], ids.map(|id| hints[id])))
        .collect();

    // Upper bounds of pieces too small to split further.
    let mut settled = 0.0_f64;
    let mut splits = 0;
    loop {
        let open = queue.peek().map_or(0.0, |piece| piece.upper);
        let upper = open.max(settled).max(best.lower);
        if upper - best.lower <= accuracy || splits >= MAX_SPLITS {
            best.upper = upper;
            return best;
        }
        let Some(piece) = queue.pop() else {
            best.upper = upper;
            return best;
        };
        if piece.upper <= best.lower {
            continue;
        }
        let [p0, p1, p2] = piece.corners;
        let lengths = [
            (p1 - p0).length_squared(),
            (p2 - p1).length_squared(),
            (p0 - p2).length_squared(),
        ];
        let edge = if lengths[0] >= lengths[1] && lengths[0] >= lengths[2] {
            0
        } else if lengths[1] >= lengths[2] {
            1
        } else {
            2
        };
        let (i, j, k) = (edge, (edge + 1) % 3, (edge + 2) % 3);
        let (vi, vj, vk) = (piece.corners[i], piece.corners[j], piece.corners[k]);
        let size = vi.length().max(vj.length()).max(vk.length());
        if lengths[edge].sqrt() <= 1024.0 * EPS * size {
            settled = settled.max(piece.upper);
            continue;
        }
        splits += 1;
        let m = (vi + vj) * 0.5;
        let dm = piece.drift[i].max(piece.drift[j]) + EPS * m.length();
        let hm = sample(m, dm, &mut best);
        let [di, dj, dk] = [piece.drift[i], piece.drift[j], piece.drift[k]];
        let [hi, hj, hk] = [piece.hints[i], piece.hints[j], piece.hints[k]];
        for (corners, drift, hints) in [
            ([vi, m, vk], [di, dm, dk], [hi, hm, hk]),
            ([m, vj, vk], [dm, dj, dk], [hm, hj, hk]),
        ] {
            let child = make(corners, drift, hints);
            if child.upper > best.lower {
                queue.push(child);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn facet(corners: [Point3; 3]) -> Facet {
        let [a, b, c] = corners;
        Facet {
            corners,
            scale: a.length() + b.length() + c.length(),
            reach: a.length().max(b.length()).max(c.length()),
            thin: false,
            lo: a.min(b).min(c),
            hi: a.max(b).max(c),
        }
    }

    #[test]
    fn point_bounds_bracket_the_distance_to_a_triangle() {
        let t = facet([
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ]);
        // Above the face, beyond an edge, beyond a corner.
        for (v, d) in [
            (Point3::new(0.5, 0.5, 3.0), 3.0),
            (Point3::new(1.0, -4.0, 3.0), 5.0),
            (Point3::new(-3.0, -4.0, 0.0), 5.0),
            (Point3::new(0.25, 0.25, 0.0), 0.0),
        ] {
            let p = nearest_on(v, &t);
            let (lo, hi) = (lower_at(v, p, &t), upper_at(v, p, &t));
            assert!(lo <= d && d <= hi, "{lo} <= {d} <= {hi}");
            assert!(hi - lo < 1e-13, "{lo} {hi}");
        }
    }

    #[test]
    fn a_thin_triangle_is_measured_along_its_edges() {
        let t = Facet {
            thin: true,
            ..facet([
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ])
        };
        let v = Point3::new(1.5, 2.0, 0.0);
        let p = nearest_on(v, &t);
        assert!((lower_at(v, p, &t) - 2.0).abs() < 1e-13);
        assert!((upper_at(v, p, &t) - 2.0).abs() < 1e-13);
    }

    #[test]
    fn box_gaps_never_exceed_the_distance_to_the_box() {
        let (lo, hi) = (Point3::ZERO, Point3::ONE);
        assert!(box_gap(Point3::new(4.0, 5.0, 0.5), lo, hi) <= 5.0);
        assert!(box_gap(Point3::new(4.0, 5.0, 0.5), lo, hi) > 5.0 - 1e-14);
        assert_eq!(box_gap(Point3::new(0.5, 0.5, 0.5), lo, hi), 0.0);
    }

    #[test]
    fn a_folded_pair_bounds_a_piece_across_its_hull() {
        // A valley: two triangles sharing the y axis, rising to height `h`
        // on either side. A piece spanning the valley's rim touches the
        // surface at its corners, yet its middle hovers over the valley:
        // the bound must pay the fold's excess, not read zero from corners.
        // Shallower than 45 degrees, so the pair projects to a convex
        // quadrilateral along either side's normal.
        let h = 0.3;
        let (a, b) = (Point3::new(0.0, -1.0, 0.0), Point3::new(0.0, 1.0, 0.0));
        let (c, d) = (Point3::new(-1.0, 0.0, h), Point3::new(1.0, 0.0, h));
        let surface = Surface {
            triangles: vec![[a, b, c], [a, b, d]],
            indices: vec![[0, 1, 2], [0, 1, 3]],
            positions: vec![a, b, c, d],
            used: vec![true; 4],
        };
        let target = Target::new(surface, 100.0);
        assert_eq!(
            target.patches.len(),
            1,
            "the fold is a flat-enough pair here"
        );
        // Every corner lies on the surface: the pair read without its
        // excess would bound the piece by rounding alone.
        let piece = [c, d, (b + c) * 0.5];
        let bound = target.piece(&piece, [0, 1, 0]);
        // The middle of `cd` is h / sqrt(1 + h^2) from either side.
        let middle = (c + d) * 0.5;
        let truth = [0, 1]
            .map(|i| (middle - nearest_on(middle, &target.facets[i])).length())
            .into_iter()
            .fold(f64::INFINITY, f64::min);
        assert!((truth - h / (1.0 + h * h).sqrt()).abs() < 1e-12);
        assert!(bound >= truth, "bound {bound} below {truth}");
    }
}
