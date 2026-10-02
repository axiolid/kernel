//! The input surface, kept unchanged while the mesh is rewritten, and the
//! closest-point queries that project vertices back onto it.
//!
//! The input's triangles are grouped into *patches*: the connected pieces
//! left when every protected edge (boundary or sharp feature) is cut. A
//! free vertex of the remeshed surface lies inside one patch and is
//! projected onto that patch only, so a vertex next to a sharp edge never
//! jumps onto the face on the other side of it. Each patch has its own
//! bounding-volume hierarchy.

use axiolid_core::{Point3, Scalar};

/// Triangles per leaf of the hierarchy.
const LEAF_SIZE: usize = 4;

/// The input surface, split into patches, with a closest-point query per
/// patch.
#[derive(Debug, Clone)]
pub(super) struct Reference {
    triangles: Vec<[Point3; 3]>,
    patches: Vec<Bvh>,
}

impl Reference {
    /// `triangles[i]` is input face `i` and `patch[i]` its patch, numbered
    /// densely from zero.
    pub(super) fn new(triangles: Vec<[Point3; 3]>, patch: &[u32]) -> Self {
        let count = patch.iter().map(|&p| p as usize + 1).max().unwrap_or(0);
        let mut members: Vec<Vec<u32>> = vec![Vec::new(); count];
        for (face, &p) in patch.iter().enumerate() {
            members[p as usize].push(face as u32);
        }
        let patches = members
            .into_iter()
            .map(|faces| Bvh::build(&triangles, faces))
            .collect();
        Self { triangles, patches }
    }

    /// Closest point to `point` on the given patch.
    pub(super) fn closest_on_patch(&self, patch: u32, point: Point3) -> Point3 {
        self.patches[patch as usize]
            .closest(&self.triangles, point, Scalar::INFINITY)
            .map_or(point, |(_, q)| q)
    }

    /// Distance from `point` to the whole input surface.
    pub(super) fn distance(&self, point: Point3) -> Scalar {
        let mut best = Scalar::INFINITY;
        for bvh in &self.patches {
            if let Some((d2, _)) = bvh.closest(&self.triangles, point, best) {
                best = best.min(d2);
            }
        }
        best.sqrt()
    }
}

#[derive(Debug, Clone, Copy)]
struct Node {
    min: Point3,
    max: Point3,
    /// Leaf: first index into `order`. Inner: index of the left child; the
    /// right child is `first + 1`.
    first: u32,
    /// Triangles in a leaf; zero for an inner node.
    count: u32,
}

/// Median-split bounding-volume hierarchy over one patch's triangles.
#[derive(Debug, Clone)]
struct Bvh {
    nodes: Vec<Node>,
    order: Vec<u32>,
}

impl Bvh {
    fn build(triangles: &[[Point3; 3]], mut order: Vec<u32>) -> Self {
        let mut nodes = Vec::with_capacity(order.len().max(1) * 2 / LEAF_SIZE + 1);
        if !order.is_empty() {
            nodes.push(bounds(triangles, &order));
            let len = order.len();
            split(triangles, &mut order, &mut nodes, 0, 0, len);
        }
        Self { nodes, order }
    }

    /// Squared distance and closest point, if any triangle is nearer than
    /// `bound` (squared distance).
    fn closest(
        &self,
        triangles: &[[Point3; 3]],
        point: Point3,
        bound: Scalar,
    ) -> Option<(Scalar, Point3)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best: Option<(Scalar, Point3)> = None;
        let mut limit = bound;
        let mut stack = vec![0u32];
        while let Some(index) = stack.pop() {
            let node = self.nodes[index as usize];
            if box_distance_squared(&node, point) >= limit {
                continue;
            }
            if node.count > 0 {
                let start = node.first as usize;
                for &t in &self.order[start..start + node.count as usize] {
                    let q = closest_on_triangle(point, &triangles[t as usize]);
                    let d2 = (q - point).length_squared();
                    if d2 < limit {
                        limit = d2;
                        best = Some((d2, q));
                    }
                }
            } else {
                let left = node.first;
                let right = left + 1;
                let dl = box_distance_squared(&self.nodes[left as usize], point);
                let dr = box_distance_squared(&self.nodes[right as usize], point);
                // Nearer child popped first; ties go left, so the walk is
                // fixed by the tree alone.
                if dl <= dr {
                    stack.push(right);
                    stack.push(left);
                } else {
                    stack.push(left);
                    stack.push(right);
                }
            }
        }
        best
    }
}

fn split(
    triangles: &[[Point3; 3]],
    order: &mut [u32],
    nodes: &mut Vec<Node>,
    node: usize,
    start: usize,
    end: usize,
) {
    if end - start <= LEAF_SIZE {
        nodes[node].first = start as u32;
        nodes[node].count = (end - start) as u32;
        return;
    }
    let centroid = |t: u32| {
        let [a, b, c] = triangles[t as usize];
        (a + b + c) / 3.0
    };
    let mut lo = Point3::splat(Scalar::INFINITY);
    let mut hi = Point3::splat(Scalar::NEG_INFINITY);
    for &t in &order[start..end] {
        let c = centroid(t);
        lo = lo.min(c);
        hi = hi.max(c);
    }
    let extent = hi - lo;
    let axis = if extent.x >= extent.y && extent.x >= extent.z {
        0
    } else if extent.y >= extent.z {
        1
    } else {
        2
    };
    order[start..end].sort_by(|&p, &q| {
        centroid(p)[axis]
            .total_cmp(&centroid(q)[axis])
            .then(p.cmp(&q))
    });
    let mid = start + (end - start) / 2;
    let left = nodes.len();
    nodes.push(bounds(triangles, &order[start..mid]));
    nodes.push(bounds(triangles, &order[mid..end]));
    nodes[node].first = left as u32;
    nodes[node].count = 0;
    split(triangles, order, nodes, left, start, mid);
    split(triangles, order, nodes, left + 1, mid, end);
}

fn bounds(triangles: &[[Point3; 3]], members: &[u32]) -> Node {
    let mut min = Point3::splat(Scalar::INFINITY);
    let mut max = Point3::splat(Scalar::NEG_INFINITY);
    for &t in members {
        for p in triangles[t as usize] {
            min = min.min(p);
            max = max.max(p);
        }
    }
    Node {
        min,
        max,
        first: 0,
        count: 0,
    }
}

fn box_distance_squared(node: &Node, p: Point3) -> Scalar {
    let below = (node.min - p).max(Point3::ZERO);
    let above = (p - node.max).max(Point3::ZERO);
    (below + above).length_squared()
}

/// Closest point of a triangle (Ericson, *Real-Time Collision Detection*,
/// 5.1.5), by Voronoi region of the corners, the edges and the interior.
pub(super) fn closest_on_triangle(p: Point3, [a, b, c]: &[Point3; 3]) -> Point3 {
    let (a, b, c) = (*a, *b, *c);
    let ab = b - a;
    let ac = c - a;
    let ap = p - a;
    let d1 = ab.dot(ap);
    let d2 = ac.dot(ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let d3 = ab.dot(bp);
    let d4 = ac.dot(bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let d5 = ab.dot(cp);
    let d6 = ac.dot(cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    a + ab * (vb * denom) + ac * (vc * denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fixed pseudo-random sequence in `[-1, 1)`.
    fn sequence(seed: u64) -> impl FnMut() -> Scalar {
        let mut state = seed;
        move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as Scalar / (1u64 << 52) as Scalar - 1.0
        }
    }

    #[test]
    fn hierarchy_finds_the_nearest_point_of_brute_force() {
        let mut next = sequence(7);
        let mut point = || Point3::new(next(), next(), next());
        let triangles: Vec<[Point3; 3]> = (0..300)
            .map(|_| {
                let centre = point() * 4.0;
                [
                    centre + point() * 0.3,
                    centre + point() * 0.3,
                    centre + point() * 0.3,
                ]
            })
            .collect();
        // Two patches: even and odd triangles.
        let patch: Vec<u32> = (0..triangles.len() as u32).map(|i| i % 2).collect();
        let reference = Reference::new(triangles.clone(), &patch);
        for _ in 0..500 {
            let q = point() * 5.0;
            let brute = |keep: &dyn Fn(usize) -> bool| {
                triangles
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| keep(*i))
                    .map(|(_, t)| (closest_on_triangle(q, t) - q).length())
                    .fold(Scalar::INFINITY, Scalar::min)
            };
            let all = brute(&|_| true);
            assert!((reference.distance(q) - all).abs() <= 1e-12 * (1.0 + all));
            for p in 0..2u32 {
                let expected = brute(&|i| i as u32 % 2 == p);
                let found = (reference.closest_on_patch(p, q) - q).length();
                assert!((found - expected).abs() <= 1e-12 * (1.0 + expected));
            }
        }
    }

    #[test]
    fn closest_point_on_a_triangle_by_region() {
        let t = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];
        let cases = [
            ((0.5, 0.5, 1.0), (0.5, 0.5, 0.0)),   // face
            ((-1.0, -1.0, 0.0), (0.0, 0.0, 0.0)), // corner a
            ((3.0, -1.0, 0.0), (2.0, 0.0, 0.0)),  // corner b
            ((-1.0, 3.0, 0.0), (0.0, 2.0, 0.0)),  // corner c
            ((1.0, -1.0, 0.0), (1.0, 0.0, 0.0)),  // edge ab
            ((-1.0, 1.0, 0.0), (0.0, 1.0, 0.0)),  // edge ac
            ((2.0, 2.0, 0.0), (1.0, 1.0, 0.0)),   // edge bc
        ];
        for ((px, py, pz), (ex, ey, ez)) in cases {
            let found = closest_on_triangle(Point3::new(px, py, pz), &t);
            assert!(
                (found - Point3::new(ex, ey, ez)).length() < 1e-15,
                "{found:?}"
            );
        }
    }
}
