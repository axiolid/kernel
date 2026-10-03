//! Certify how far a parameterised exact surface lies from a mesh (#232).
//!
//! Used where no construction proof exists and no per-triangle argument
//! applies: a disk swept along a B-spline or an ellipse: the mesh's stations are
//! placed on frames read from the sampled path, so its vertices are not at
//! known surface parameters. The bound is computed against the exact
//! surface alone, by branch and bound over its parameter domain:
//!
//! - A cell `X` with centre `x_c` and half-widths `h` maps into
//!   `T(x_c) + DT(x_c) [-h, h] + B(e2)`, a parallelogram grown by a ball,
//!   by Taylor's theorem with `e2 = 1/2 (A h_x^2 + 2 B h_x h_y + C h_y^2)`
//!   and `A, B, C` certified bounds on `|T_xx|, |T_xy|, |T_yy|` over `X`.
//! - The distance to the mesh is 1-Lipschitz, and the distance to one
//!   triangle is convex, so over a piece of the parallelogram it is at most
//!   the largest of its corners' distances to any one triangle. Taking the
//!   best triangle per piece and the worst piece bounds the distance from
//!   every point of the parallelogram to the mesh; adding `e2` bounds it
//!   over `T(X)`.
//! - The exact distance from `T(x_c)` to the mesh is a sample of the true
//!   maximum, so it bounds that maximum from below.
//!
//! Cells are split, worst upper bound first, until the worst is within 10%
//! of the best lower bound (or a tenth of the target, below which the gap
//! does not matter), or the work budget runs out. The returned value is the
//! worst upper bound over all cells, so stopping early only loosens it.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::TriMesh;

/// A parameterised piece of the exact surface.
pub(crate) trait Patch {
    /// Point and first partials at `(x, y)`; `None` where undefined.
    fn jet(&self, x: Scalar, y: Scalar) -> Option<(Point3, Vec3, Vec3)>;
    /// Certified `sup |T_xx|, |T_xy|, |T_yy|` over the box.
    fn second(&self, x: (Scalar, Scalar), y: (Scalar, Scalar)) -> Option<[Scalar; 3]>;
}

/// One parameter box of one patch.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Cell {
    pub patch: usize,
    pub x: (Scalar, Scalar),
    pub y: (Scalar, Scalar),
}

/// Cell evaluations before the search stops and reports what it has.
const MAX_CELLS: usize = 1 << 19;
/// Stop once the worst upper bound is within this factor of the best
/// lower bound.
const TIGHT: Scalar = 1.1;
/// A cell whose image is this small relative to the target is not split.
const SMALLEST: Scalar = 1e-4;

struct Scored {
    upper: Scalar,
    cell: Cell,
}

impl PartialEq for Scored {
    fn eq(&self, other: &Self) -> bool {
        self.upper.total_cmp(&other.upper) == Ordering::Equal
    }
}
impl Eq for Scored {}
impl PartialOrd for Scored {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Scored {
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper.total_cmp(&other.upper)
    }
}

/// A certified upper bound on the distance from every point of the
/// patches' cells to the mesh; `None` when a patch cannot be bounded
/// somewhere (an undefined point, an unbounded derivative) or the mesh has
/// no triangle.
pub(crate) fn certify(
    patches: &[&dyn Patch],
    cells: Vec<Cell>,
    mesh: &TriMesh,
    target: Scalar,
) -> Option<Scalar> {
    let index = TriangleIndex::new(mesh)?;
    let mut lower: Scalar = 0.0;
    let mut heap = BinaryHeap::with_capacity(cells.len());
    let mut work = 0_usize;
    for cell in cells {
        let (upper, sample) = evaluate(patches[cell.patch], &cell, &index)?;
        lower = lower.max(sample);
        heap.push(Scored { upper, cell });
        work += 1;
    }
    while let Some(top) = heap.pop() {
        let settled = top.upper <= (TIGHT * lower).max(0.1 * target);
        let cell = top.cell;
        let patch = patches[cell.patch];
        let (_, du, dv) = patch.jet(mid(cell.x), mid(cell.y))?;
        let extent_x = du.length() * (cell.x.1 - cell.x.0);
        let extent_y = dv.length() * (cell.y.1 - cell.y.0);
        let tiny = extent_x.max(extent_y) <= SMALLEST * target;
        if settled || tiny || work >= MAX_CELLS {
            return Some(top.upper * (1.0 + 1e-9));
        }
        let halves = if extent_x >= extent_y {
            let m = mid(cell.x);
            [
                Cell {
                    x: (cell.x.0, m),
                    ..cell
                },
                Cell {
                    x: (m, cell.x.1),
                    ..cell
                },
            ]
        } else {
            let m = mid(cell.y);
            [
                Cell {
                    y: (cell.y.0, m),
                    ..cell
                },
                Cell {
                    y: (m, cell.y.1),
                    ..cell
                },
            ]
        };
        for half in halves {
            let (upper, sample) = evaluate(patch, &half, &index)?;
            lower = lower.max(sample);
            heap.push(Scored { upper, cell: half });
            work += 1;
        }
    }
    Some(0.0)
}

fn mid(range: (Scalar, Scalar)) -> Scalar {
    0.5 * (range.0 + range.1)
}

/// Pieces per side the parallelogram is cut into.
const PIECES: usize = 2;

/// `(upper bound over the cell, exact distance at its centre)`.
fn evaluate(patch: &dyn Patch, cell: &Cell, index: &TriangleIndex) -> Option<(Scalar, Scalar)> {
    let (q, du, dv) = patch.jet(mid(cell.x), mid(cell.y))?;
    let [a, b, c] = patch.second(cell.x, cell.y)?;
    let (hx, hy) = (0.5 * (cell.x.1 - cell.x.0), 0.5 * (cell.y.1 - cell.y.0));
    let e2 = 0.5 * (a * hx * hx + 2.0 * b * hx * hy + c * hy * hy);
    let reach = du.length() * hx + dv.length() * hy;
    let centre = index.nearest(q);
    if !(e2.is_finite() && reach.is_finite() && centre.is_finite()) {
        return None;
    }
    let candidates = index.within(q, centre + 2.0 * reach);
    // Distances from the lattice of piece corners to each candidate.
    let side = PIECES + 1;
    let mut lattice = Vec::with_capacity(side * side);
    for i in 0..side {
        for j in 0..side {
            let sx = -1.0 + 2.0 * i as Scalar / PIECES as Scalar;
            let sy = -1.0 + 2.0 * j as Scalar / PIECES as Scalar;
            let p = q + du * (sx * hx) + dv * (sy * hy);
            lattice.push(
                candidates
                    .iter()
                    .map(|&t| index.distance(p, t))
                    .collect::<Vec<_>>(),
            );
        }
    }
    let mut worst: Scalar = 0.0;
    for i in 0..PIECES {
        for j in 0..PIECES {
            let corners = [
                i * side + j,
                (i + 1) * side + j,
                i * side + j + 1,
                (i + 1) * side + j + 1,
            ];
            let best = (0..candidates.len())
                .map(|k| {
                    corners
                        .iter()
                        .map(|&c| lattice[c][k])
                        .fold(0.0, Scalar::max)
                })
                .fold(Scalar::INFINITY, Scalar::min);
            worst = worst.max(best);
        }
    }
    let upper = worst + e2;
    upper.is_finite().then_some((upper, centre))
}

/// Triangles of a mesh in a bounding-volume hierarchy, for exact nearest
/// distances and range queries.
pub(crate) struct TriangleIndex {
    triangles: Vec<[Point3; 3]>,
    nodes: Vec<Node>,
    order: Vec<usize>,
}

struct Node {
    lo: Vec3,
    hi: Vec3,
    /// Leaf: `order[start..end]`; inner: children `start` and `end`.
    start: usize,
    end: usize,
    leaf: bool,
}

const LEAF: usize = 8;

impl TriangleIndex {
    pub(crate) fn new(mesh: &TriMesh) -> Option<Self> {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .chunks_exact(3)
            .filter_map(|t| {
                let corners = [0, 1, 2].map(|k| mesh.positions.get(t[k] as usize).copied());
                match corners {
                    [Some(a), Some(b), Some(c)] => Some([a, b, c]),
                    _ => None,
                }
            })
            .filter(|t| t.iter().all(|p| p.is_finite()))
            .collect();
        if triangles.is_empty() {
            return None;
        }
        let mut index = Self {
            order: (0..triangles.len()).collect(),
            triangles,
            nodes: Vec::new(),
        };
        let count = index.order.len();
        index.build(0, count);
        Some(index)
    }

    fn bounds(&self, range: &[usize]) -> (Vec3, Vec3) {
        let mut lo = Vec3::splat(Scalar::INFINITY);
        let mut hi = Vec3::splat(Scalar::NEG_INFINITY);
        for &t in range {
            for p in self.triangles[t] {
                lo = lo.min(p);
                hi = hi.max(p);
            }
        }
        (lo, hi)
    }

    fn build(&mut self, start: usize, end: usize) -> usize {
        let (lo, hi) = self.bounds(&self.order[start..end]);
        let id = self.nodes.len();
        self.nodes.push(Node {
            lo,
            hi,
            start,
            end,
            leaf: true,
        });
        if end - start <= LEAF {
            return id;
        }
        let extent = hi - lo;
        let axis = if extent.x >= extent.y && extent.x >= extent.z {
            0
        } else if extent.y >= extent.z {
            1
        } else {
            2
        };
        let key = |t: &[Point3; 3]| (t[0][axis] + t[1][axis] + t[2][axis]) / 3.0;
        let triangles = &self.triangles;
        self.order[start..end].sort_by(|&a, &b| key(&triangles[a]).total_cmp(&key(&triangles[b])));
        let middle = start + (end - start) / 2;
        let left = self.build(start, middle);
        let right = self.build(middle, end);
        self.nodes[id] = Node {
            lo,
            hi,
            start: left,
            end: right,
            leaf: false,
        };
        id
    }

    fn box_distance(node: &Node, p: Point3) -> Scalar {
        let d = (node.lo - p).max(p - node.hi).max(Vec3::ZERO);
        d.length()
    }

    /// Exact distance from `p` to triangle `t`.
    pub(crate) fn distance(&self, p: Point3, t: usize) -> Scalar {
        point_triangle_distance(p, self.triangles[t])
    }

    /// Exact distance from `p` to the nearest triangle.
    pub(crate) fn nearest(&self, p: Point3) -> Scalar {
        let mut best = Scalar::INFINITY;
        let mut stack = vec![0_usize];
        while let Some(id) = stack.pop() {
            let node = &self.nodes[id];
            if Self::box_distance(node, p) >= best {
                continue;
            }
            if node.leaf {
                for &t in &self.order[node.start..node.end] {
                    best = best.min(self.distance(p, t));
                }
            } else {
                let (a, b) = (node.start, node.end);
                let (da, db) = (
                    Self::box_distance(&self.nodes[a], p),
                    Self::box_distance(&self.nodes[b], p),
                );
                if da <= db {
                    stack.push(b);
                    stack.push(a);
                } else {
                    stack.push(a);
                    stack.push(b);
                }
            }
        }
        best
    }

    /// Every triangle whose box is within `radius` of `p`.
    pub(crate) fn within(&self, p: Point3, radius: Scalar) -> Vec<usize> {
        let mut out = Vec::new();
        let mut stack = vec![0_usize];
        while let Some(id) = stack.pop() {
            let node = &self.nodes[id];
            if Self::box_distance(node, p) > radius {
                continue;
            }
            if node.leaf {
                out.extend_from_slice(&self.order[node.start..node.end]);
            } else {
                stack.push(node.start);
                stack.push(node.end);
            }
        }
        out
    }
}

/// Distance from `p` to a triangle (Ericson, Real-Time Collision Detection
/// 5.1.5); a triangle with no area is measured as its three edges.
pub(crate) fn point_triangle_distance(p: Point3, [a, b, c]: [Point3; 3]) -> Scalar {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let area = ab.cross(ac).length_squared();
    let scale = ab.length_squared().max(ac.length_squared());
    if area <= 1e-24 * scale * scale || area.is_nan() {
        return segment_distance(p, a, b)
            .min(segment_distance(p, b, c))
            .min(segment_distance(p, c, a));
    }
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.length();
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return bp.length();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (p - (a + ab * v)).length();
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return cp.length();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (p - (a + ac * w)).length();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (b + (c - b) * w)).length();
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    (p - (a + ab * v + ac * w)).length()
}

fn segment_distance(p: Point3, a: Point3, b: Point3) -> Scalar {
    let ab = b - a;
    let len = ab.length_squared();
    let t = if len > 0.0 {
        ((p - a).dot(ab) / len).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p - (a + ab * t)).length()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A cylinder of radius `r` about z, in angle and height.
    struct Cylinder(Scalar);

    impl Patch for Cylinder {
        fn jet(&self, phi: Scalar, z: Scalar) -> Option<(Point3, Vec3, Vec3)> {
            let (s, c) = phi.sin_cos();
            Some((
                Point3::new(self.0 * c, self.0 * s, z),
                Vec3::new(-self.0 * s, self.0 * c, 0.0),
                Vec3::Z,
            ))
        }

        fn second(&self, _x: (Scalar, Scalar), _y: (Scalar, Scalar)) -> Option<[Scalar; 3]> {
            Some([self.0, 0.0, 0.0])
        }
    }

    /// The face of a circumscribed prism tangent at angle 0: the parallelogram
    /// of the cell lies in it, every corner at distance zero, so the whole
    /// bound is the Taylor remainder, and the exact surface curves away
    /// from the face by `r (1 - cos h)` at the cell's edge.
    #[test]
    fn the_taylor_remainder_carries_a_cell_whose_corners_touch_the_mesh() {
        let r = 1.0;
        // One triangle of the face covers the whole cell, so no piece of
        // the parallelogram straddles a diagonal.
        let face = TriMesh::new(
            vec![
                Point3::new(r, -1.0, -1.0),
                Point3::new(r, 1.0, -1.0),
                Point3::new(r, 0.0, 3.0),
            ],
            vec![0, 1, 2],
        );
        let h = 0.2;
        let cell = Cell {
            patch: 0,
            x: (-h, h),
            y: (0.2, 0.8),
        };
        let bound = certify(&[&Cylinder(r)], vec![cell], &face, 1e3).unwrap();
        let measured = r * (1.0 - h.cos());
        assert!(bound >= measured, "{bound} below the edge's {measured}");
        assert!(bound <= 1.1 * 0.5 * r * h * h, "{bound}");
    }

    /// Two faces of a circumscribed prism meeting over angle 0: the exact
    /// surface is furthest from them there, at the cell's centre, and only
    /// the corner furthest from each face covers it.
    #[test]
    fn a_cell_over_a_prism_edge_is_bounded_by_its_far_corners() {
        let r = 1.0;
        let half = core::f64::consts::PI / 12.0;
        let big = r / half.cos();
        let at = |angle: Scalar, z: Scalar| Point3::new(big * angle.cos(), big * angle.sin(), z);
        let prism = TriMesh::new(
            vec![
                at(-2.0 * half, 0.0),
                at(0.0, 0.0),
                at(2.0 * half, 0.0),
                at(-2.0 * half, 1.0),
                at(0.0, 1.0),
                at(2.0 * half, 1.0),
            ],
            vec![0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4],
        );
        let h = 0.1;
        let cell = Cell {
            patch: 0,
            x: (-h, h),
            y: (0.2, 0.8),
        };
        let bound = certify(&[&Cylinder(r)], vec![cell], &prism, 1e3).unwrap();
        let measured = r * (1.0 - half.cos());
        assert!(bound >= measured, "{bound} below the centre's {measured}");
    }

    #[test]
    fn point_triangle_distance_matches_brute_force() {
        let t = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.2, 0.1),
            Point3::new(0.3, 1.1, -0.2),
        ];
        for i in 0..20 {
            for j in 0..20 {
                let p = Point3::new(-0.5 + 0.1 * i as Scalar, -0.5 + 0.1 * j as Scalar, 0.37);
                let mut brute = Scalar::INFINITY;
                for u in 0..=200 {
                    for v in 0..=(200 - u) {
                        let (u, v) = (u as Scalar / 200.0, v as Scalar / 200.0);
                        let q = t[0] + (t[1] - t[0]) * u + (t[2] - t[0]) * v;
                        brute = brute.min((p - q).length());
                    }
                }
                let exact = point_triangle_distance(p, t);
                assert!(
                    exact <= brute + 1e-12 && brute - exact < 0.02,
                    "{exact} vs {brute}"
                );
            }
        }
    }
}
