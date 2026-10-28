//! In-plane regions of a section: faces lying in the cutting plane.
//!
//! Every decision is an exact `orient2d` on the original binary64
//! coordinates with one axis dropped (`Projection`), never on projected or
//! interpolated values. The pipeline, in order:
//!
//! 1. `InPlaneSet`: the in-plane triangles with their facing, and the
//!    directed-edge map that joins them.
//! 2. `check_overlaps`: in-plane triangles must have disjoint interiors; a
//!    sweep over bounding boxes proposes pairs and a separating-axis test with
//!    exact predicates decides them.
//! 3. `build_regions`: same-facing triangles sharing an edge form a
//!    component. Its boundary edges are linked into cycles by fan order at
//!    each vertex, so a vertex where a region pinches is resolved per fan.
//!    Each cycle is an outer ring or a hole by the exact turn at its
//!    lexicographically smallest vertex, collinear vertices are dropped, and
//!    the rings are projected to the plane frame, oriented and made canonical.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation, Sign};
use axiolid_core::{Frame3, Point2, Point3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    InPlaneRegion, InPlaneSide, RegionLimits, REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES,
    REFUSED_OVERLAPPING_IN_PLANE_FACES,
};

use super::{
    internal_topology_error, mesh_index, orient_sign, project_local, sequence_order, sign_times,
    ExactSectionPlane, ScalarSection,
};
use crate::coplanar::{dominant_axis, project};

/// The axis dropped to reach two dimensions, and what that drop does to the
/// frame's own orientation.
#[derive(Debug, Clone, Copy)]
pub(super) struct Projection {
    axis: usize,
    frame_sign: Sign,
}

impl Projection {
    /// Pick the dominant axis of the plane normal, or the first axis that
    /// keeps the frame non-degenerate.
    ///
    /// Dropping an axis is an exact affine bijection on the plane when the
    /// projected frame points `(origin, origin + x, origin + y)` are not
    /// collinear; their exact orientation is `frame_sign`.
    pub(super) fn new(plane: &ExactSectionPlane) -> GeomResult<Self> {
        let dominant = dominant_axis(plane.normal);
        for axis in [dominant, 0, 1, 2] {
            let frame_sign = orient_sign(
                project(plane.origin, axis),
                project(plane.x_point, axis),
                project(plane.y_point, axis),
            )?;
            if frame_sign != Sign::Zero {
                return Ok(Self { axis, frame_sign });
            }
        }
        Err(GeomError::Degenerate(
            "section plane has no non-degenerate axis projection".into(),
        ))
    }

    fn flat(self, point: Point3) -> Point2 {
        project(point, self.axis)
    }

    /// Exact orientation of three points of the plane in the frame's own
    /// `(x, y)` sense: positive when they turn counter-clockwise there.
    fn orient(self, a: Point3, b: Point3, c: Point3) -> GeomResult<Sign> {
        let projected = orient_sign(self.flat(a), self.flat(b), self.flat(c))?;
        Ok(sign_times(projected, self.frame_sign))
    }
}

/// One source triangle lying in the plane.
#[derive(Debug, Clone, Copy)]
pub(super) struct InPlaneTriangle {
    /// Source vertices in winding order.
    pub(super) vertices: [u32; 3],
    /// Positive when the triangle is counter-clockwise in the frame, so its
    /// normal agrees with the frame's `+z`; negative otherwise.
    pub(super) facing: Sign,
}

impl InPlaneTriangle {
    pub(super) fn new(
        mesh: &TriMesh,
        projection: Projection,
        index: usize,
        vertices: [u32; 3],
    ) -> GeomResult<Self> {
        let [a, b, c] = corners(mesh, vertices)?;
        let facing = projection.orient(a, b, c)?;
        if facing == Sign::Zero {
            return Err(GeomError::Degenerate(format!(
                "in-plane source triangle {index} has zero area"
            )));
        }
        Ok(Self { vertices, facing })
    }

    /// The side its material lies on: below the plane when the normal agrees
    /// with the frame's `+z`.
    fn side(self) -> InPlaneSide {
        if self.facing == Sign::Positive {
            InPlaneSide::Below
        } else {
            InPlaneSide::Above
        }
    }
}

fn corners(mesh: &TriMesh, vertices: [u32; 3]) -> GeomResult<[Point3; 3]> {
    let mut points = [Point3::ZERO; 3];
    for (slot, vertex) in points.iter_mut().zip(vertices) {
        *slot = mesh.positions[mesh_index(vertex)?];
    }
    Ok(points)
}

/// The in-plane triangles and their directed edges.
#[derive(Debug)]
pub(super) struct InPlaneSet {
    triangles: Vec<InPlaneTriangle>,
    /// Directed edge `(from, to)` to the position of its triangle.
    directed: BTreeMap<(u32, u32), usize>,
}

impl InPlaneSet {
    pub(super) fn new(triangles: Vec<InPlaneTriangle>) -> GeomResult<Self> {
        let mut directed = BTreeMap::new();
        for (position, triangle) in triangles.iter().enumerate() {
            let [a, b, c] = triangle.vertices;
            for edge in [(a, b), (b, c), (c, a)] {
                if directed.insert(edge, position).is_some() {
                    return Err(GeomError::NotManifold(format!(
                        "in-plane mesh edge {edge:?} is used twice in the same direction"
                    )));
                }
            }
        }
        Ok(Self {
            triangles,
            directed,
        })
    }

    pub(super) fn len(&self) -> usize {
        self.triangles.len()
    }

    /// Whether the directed edge `edge` has an in-plane twin of the same
    /// facing, that is, whether it is interior to a component.
    fn is_interior(&self, position: usize, edge: (u32, u32)) -> bool {
        self.directed
            .get(&(edge.1, edge.0))
            .is_some_and(|&twin| self.triangles[twin].facing == self.triangles[position].facing)
    }
}

fn refuse(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: ScalarSection::ID,
        operation: Operation::MeshPlaneSection,
        input,
    }
}

#[derive(Debug, Clone, Copy)]
struct BoundingBox {
    min_u: f64,
    max_u: f64,
    min_v: f64,
    max_v: f64,
}

impl BoundingBox {
    fn of(points: [Point2; 3]) -> Self {
        let fold = |select: fn(&Point2) -> f64| {
            points
                .iter()
                .map(select)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), value| {
                    (low.min(value), high.max(value))
                })
        };
        let (min_u, max_u) = fold(|point| point.x);
        let (min_v, max_v) = fold(|point| point.y);
        Self {
            min_u,
            max_u,
            min_v,
            max_v,
        }
    }
}

/// Refuse in-plane triangles whose interiors overlap; return how many
/// candidate pairs were examined.
///
/// A sweep over `u` proposes pairs whose bounding boxes overlap. Two
/// triangles of the same facing that share a mesh edge lie on opposite sides
/// of it and are skipped. Every other pair is decided by a separating-axis
/// test: counter-clockwise triangles have disjoint interiors when some edge
/// of one has all three vertices of the other on or to its right.
pub(super) fn check_overlaps(
    mesh: &TriMesh,
    projection: Projection,
    set: &InPlaneSet,
    max_pairs: usize,
    options: &ExecutionOptions,
) -> GeomResult<usize> {
    let mut flat = Vec::with_capacity(set.len());
    let mut boxes = Vec::with_capacity(set.len());
    for triangle in &set.triangles {
        let [a, b, c] = corners(mesh, triangle.vertices)?;
        let points = [projection.flat(a), projection.flat(b), projection.flat(c)];
        // Counter-clockwise in the dropped-axis coordinates.
        let counter_clockwise =
            sign_times(triangle.facing, projection.frame_sign) == Sign::Positive;
        flat.push(if counter_clockwise {
            points
        } else {
            [points[0], points[2], points[1]]
        });
        boxes.push(BoundingBox::of(points));
    }
    let mut order: Vec<usize> = (0..set.len()).collect();
    order.sort_by(|&left, &right| {
        boxes[left]
            .min_u
            .total_cmp(&boxes[right].min_u)
            .then(left.cmp(&right))
    });

    let mut pairs = 0usize;
    for (rank, &i) in order.iter().enumerate() {
        for &j in &order[rank + 1..] {
            if boxes[j].min_u > boxes[i].max_u {
                break;
            }
            options.check_cancelled()?;
            if boxes[j].min_v > boxes[i].max_v || boxes[i].min_v > boxes[j].max_v {
                continue;
            }
            let (first, second) = (set.triangles[i], set.triangles[j]);
            let shared = first
                .vertices
                .iter()
                .filter(|vertex| second.vertices.contains(vertex))
                .count();
            if shared >= 2 && first.facing == second.facing {
                continue;
            }
            pairs += 1;
            if pairs > max_pairs {
                return Err(GeomError::BudgetExceeded {
                    resource: "section in-plane pairs",
                });
            }
            if interiors_overlap(&flat[i], &flat[j])? {
                return Err(refuse(if first.facing == second.facing {
                    REFUSED_OVERLAPPING_IN_PLANE_FACES
                } else {
                    REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES
                }));
            }
        }
    }
    Ok(pairs)
}

/// Whether two counter-clockwise triangles share interior area.
fn interiors_overlap(first: &[Point2; 3], second: &[Point2; 3]) -> GeomResult<bool> {
    for (owner, other) in [(first, second), (second, first)] {
        for corner in 0..3 {
            let (p, q) = (owner[corner], owner[(corner + 1) % 3]);
            let mut separated = true;
            for &r in other {
                if orient_sign(p, q, r)? == Sign::Positive {
                    separated = false;
                    break;
                }
            }
            if separated {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Union-find over triangle positions; the root is the smallest member.
fn components(set: &InPlaneSet, options: &ExecutionOptions) -> GeomResult<Vec<usize>> {
    fn find(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }
    let mut parent: Vec<usize> = (0..set.len()).collect();
    for (position, triangle) in set.triangles.iter().enumerate() {
        options.check_cancelled()?;
        let [a, b, c] = triangle.vertices;
        for edge in [(a, b), (b, c), (c, a)] {
            if set.is_interior(position, edge) {
                let twin = set.directed[&(edge.1, edge.0)];
                let (left, right) = (find(&mut parent, position), find(&mut parent, twin));
                parent[left.max(right)] = left.min(right);
            }
        }
    }
    (0..set.len())
        .map(|node| Ok(find(&mut parent, node)))
        .collect()
}

/// The vertex after `vertex` in the winding of `triangle`.
fn successor(triangle: [u32; 3], vertex: u32) -> GeomResult<u32> {
    triangle
        .iter()
        .position(|&corner| corner == vertex)
        .map(|index| triangle[(index + 1) % 3])
        .ok_or_else(internal_topology_error)
}

/// The boundary edge that leaves `vertex` after the boundary edge that
/// entered it through triangle `entering`.
///
/// Walk the fan around `vertex` from `entering`: the edge leaving `vertex`
/// in the current triangle is either on the boundary, which ends the walk, or
/// shared with a same-facing neighbour, whose own leaving edge is next. The
/// walk stays inside one fan, so several fans meeting at a pinch vertex are
/// linked fan by fan.
fn leaving_edge(
    set: &InPlaneSet,
    boundary: &BTreeSet<(u32, u32)>,
    vertex: u32,
    entering: usize,
) -> GeomResult<(u32, u32)> {
    let mut current = entering;
    for _ in 0..=set.len() {
        let next = successor(set.triangles[current].vertices, vertex)?;
        let edge = (vertex, next);
        if boundary.contains(&edge) {
            return Ok(edge);
        }
        current = *set
            .directed
            .get(&(next, vertex))
            .ok_or_else(internal_topology_error)?;
    }
    Err(internal_topology_error())
}

/// A closed boundary cycle of one component, as source vertices.
struct Cycle {
    component: usize,
    vertices: Vec<u32>,
}

fn trace_cycles(
    set: &InPlaneSet,
    roots: &[usize],
    options: &ExecutionOptions,
) -> GeomResult<Vec<Cycle>> {
    let mut boundary = BTreeSet::new();
    for (position, triangle) in set.triangles.iter().enumerate() {
        let [a, b, c] = triangle.vertices;
        for edge in [(a, b), (b, c), (c, a)] {
            if !set.is_interior(position, edge) {
                boundary.insert(edge);
            }
        }
    }

    let mut used = BTreeSet::new();
    let mut cycles = Vec::new();
    for &start in &boundary {
        if used.contains(&start) {
            continue;
        }
        let mut vertices = Vec::new();
        let mut edge = start;
        loop {
            options.check_cancelled()?;
            used.insert(edge);
            vertices.push(edge.0);
            let entering = set.directed[&edge];
            edge = leaving_edge(set, &boundary, edge.1, entering)?;
            if edge == start {
                break;
            }
            if used.contains(&edge) || vertices.len() > boundary.len() {
                return Err(internal_topology_error());
            }
        }
        cycles.push(Cycle {
            component: roots[set.directed[&start]],
            vertices,
        });
    }
    Ok(cycles)
}

/// Whether `point` lies strictly between `a` and `b` on their common line.
fn strictly_between(a: Point3, point: Point3, b: Point3) -> bool {
    let within = |a: f64, p: f64, b: f64| a.min(b) <= p && p <= a.max(b);
    point != a
        && point != b
        && within(a.x, point.x, b.x)
        && within(a.y, point.y, b.y)
        && within(a.z, point.z, b.z)
}

/// Lexicographic order of model-space points, then by source vertex.
fn corner_order(mesh: &TriMesh, left: u32, right: u32) -> GeomResult<Ordering> {
    let (a, b) = (
        mesh.positions[mesh_index(left)?],
        mesh.positions[mesh_index(right)?],
    );
    Ok(a.x
        .total_cmp(&b.x)
        .then_with(|| a.y.total_cmp(&b.y))
        .then_with(|| a.z.total_cmp(&b.z))
        .then(left.cmp(&right)))
}

/// Whether the cycle is the outer boundary of its component, by the exact
/// turn at its lexicographically smallest vertex.
fn is_outer(
    mesh: &TriMesh,
    projection: Projection,
    facing: Sign,
    ring: &[u32],
) -> GeomResult<bool> {
    let mut smallest = 0;
    for index in 1..ring.len() {
        if corner_order(mesh, ring[index], ring[smallest])?.is_lt() {
            smallest = index;
        }
    }
    let vertex = ring[smallest];
    if ring.iter().filter(|&&corner| corner == vertex).count() != 1 {
        return Err(GeomError::Degenerate(
            "in-plane region boundary visits its extreme vertex twice".into(),
        ));
    }
    let previous = ring[(smallest + ring.len() - 1) % ring.len()];
    let next = ring[(smallest + 1) % ring.len()];
    let [a, b, c] = corners(mesh, [previous, vertex, next])?;
    // The extreme vertex is a convex corner: its turn is counter-clockwise
    // for an outer ring of an up-facing component and for a hole of a
    // down-facing one; the facing product separates the cases.
    match sign_times(projection.orient(a, b, c)?, facing) {
        Sign::Positive => Ok(true),
        Sign::Negative => Ok(false),
        _ => Err(GeomError::Degenerate(
            "in-plane region boundary has a spike at its extreme vertex".into(),
        )),
    }
}

/// Drop vertices that lie strictly between their neighbours on a line,
/// except `protected` ones; a spike is refused.
fn simplify(
    mesh: &TriMesh,
    projection: Projection,
    ring: &mut Vec<u32>,
    protected: &BTreeSet<u32>,
    options: &ExecutionOptions,
) -> GeomResult<()> {
    let collapsed =
        || GeomError::Degenerate("in-plane region collapsed under simplification".into());
    let mut index = 0usize;
    let mut clean = 0usize;
    while clean < ring.len() {
        options.check_cancelled()?;
        if ring.len() < 3 {
            return Err(collapsed());
        }
        let at = index % ring.len();
        let [a, b, c] = corners(
            mesh,
            [
                ring[(at + ring.len() - 1) % ring.len()],
                ring[at],
                ring[(at + 1) % ring.len()],
            ],
        )?;
        if projection.orient(a, b, c)? != Sign::Zero {
            clean += 1;
            index = at + 1;
            continue;
        }
        if !strictly_between(a, b, c) {
            return Err(GeomError::Degenerate(
                "in-plane region boundary has a zero-width spike".into(),
            ));
        }
        if protected.contains(&ring[at]) {
            clean += 1;
            index = at + 1;
            continue;
        }
        ring.remove(at);
        clean = 0;
        index = (at + ring.len().saturating_sub(1)) % ring.len().max(1);
    }
    if ring.len() < 3 {
        return Err(collapsed());
    }
    Ok(())
}

/// Plane-local points of a ring, reversed for a down-facing component,
/// rotated to start at the smallest point.
fn plane_ring(
    mesh: &TriMesh,
    frame: Frame3,
    ring: &[u32],
    reverse: bool,
) -> GeomResult<Vec<Point2>> {
    let mut points = ring
        .iter()
        .map(|&vertex| project_local(frame, mesh.positions[mesh_index(vertex)?]))
        .collect::<GeomResult<Vec<_>>>()?;
    if reverse {
        points.reverse();
    }
    let collapsed = || GeomError::Degenerate("in-plane region collapsed under projection".into());
    if points.len() < 3
        || (0..points.len()).any(|index| points[index] == points[(index + 1) % points.len()])
    {
        return Err(collapsed());
    }
    Ok(canonical_rotation(points))
}

/// Rotate a ring to start at its lexicographically smallest point; of
/// several equal smallest points the rotation that orders lowest wins.
pub(super) fn canonical_rotation(points: Vec<Point2>) -> Vec<Point2> {
    let rotated = |start: usize| -> Vec<Point2> {
        points[start..]
            .iter()
            .chain(&points[..start])
            .copied()
            .collect()
    };
    let smallest = points
        .iter()
        .min_by(|left, right| super::point_order(left, right))
        .copied();
    (0..points.len())
        .filter(|&start| Some(points[start]) == smallest)
        .map(rotated)
        .min_by(|left, right| sequence_order(left, right))
        .unwrap_or(points)
}

/// Regions of the in-plane triangles: one per component, with its holes.
///
/// `protected` vertices (nodes of cut segments) are never simplified away,
/// so a chain end stays a vertex of the ring it meets.
pub(super) fn build_regions(
    mesh: &TriMesh,
    frame: Frame3,
    projection: Projection,
    set: &InPlaneSet,
    protected: &BTreeSet<u32>,
    limits: RegionLimits,
    options: &ExecutionOptions,
) -> GeomResult<Vec<InPlaneRegion>> {
    let roots = components(set, options)?;
    let cycles = trace_cycles(set, &roots, options)?;

    // Component root to its outer ring and holes, in plane-local points.
    let mut outers = BTreeMap::<usize, Vec<Point2>>::new();
    let mut holes = BTreeMap::<usize, Vec<Vec<Point2>>>::new();
    for mut cycle in cycles {
        let triangle = set.triangles[cycle.component];
        let facing = triangle.facing;
        let outer = is_outer(mesh, projection, facing, &cycle.vertices)?;
        simplify(mesh, projection, &mut cycle.vertices, protected, options)?;
        let ring = plane_ring(mesh, frame, &cycle.vertices, facing == Sign::Negative)?;
        if outer {
            if outers.insert(cycle.component, ring).is_some() {
                return Err(GeomError::Degenerate(
                    "in-plane component has two outer boundaries".into(),
                ));
            }
        } else {
            holes.entry(cycle.component).or_default().push(ring);
        }
    }

    let mut regions = Vec::new();
    let mut vertices = 0usize;
    for (component, outer) in outers {
        let mut rings = holes.remove(&component).unwrap_or_default();
        rings.sort_by(|left, right| sequence_order(left, right));
        if regions.len() >= limits.max_regions {
            return Err(GeomError::BudgetExceeded {
                resource: "section in-plane regions",
            });
        }
        vertices += outer.len() + rings.iter().map(Vec::len).sum::<usize>();
        if vertices > limits.max_region_vertices {
            return Err(GeomError::BudgetExceeded {
                resource: "section in-plane region vertices",
            });
        }
        regions.push(InPlaneRegion::new(
            set.triangles[component].side(),
            outer,
            rings,
        ));
    }
    if !holes.is_empty() {
        return Err(GeomError::Degenerate(
            "in-plane component has a hole but no outer boundary".into(),
        ));
    }
    regions.sort_by(|left, right| {
        left.side
            .cmp(&right.side)
            .then_with(|| sequence_order(&left.outer, &right.outer))
    });
    Ok(regions)
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Tolerance, Vec3};

    use super::*;

    /// Eight triangles fan around `c = (1,1,0)` over `[0,2]^2`; the edge
    /// midpoints `(1,0)` and `(1,2)` are lifted to `z = 1`, so the four
    /// triangles touching them leave the plane. The remaining four form two
    /// groups that meet only at `c`: a pinch. Only the in-plane triangles
    /// are given; their missing neighbours make the boundary.
    fn pinched_fan() -> (TriMesh, InPlaneSet, Frame3, Projection) {
        let positions = vec![
            Point3::new(1.0, 1.0, 0.0), // 0: c
            Point3::new(0.0, 0.0, 0.0), // 1: r0
            Point3::new(2.0, 0.0, 0.0), // 2: r2
            Point3::new(2.0, 1.0, 0.0), // 3: r3, collinear between r2 and r4
            Point3::new(2.0, 2.0, 0.0), // 4: r4
            Point3::new(0.0, 2.0, 0.0), // 5: r6
            Point3::new(0.0, 1.0, 0.0), // 6: r7, collinear between r6 and r0
        ];
        let mesh = TriMesh::new(positions, vec![]);
        let frame = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let plane = ExactSectionPlane::new(frame).expect("plane");
        let projection = Projection::new(&plane).expect("projection");
        let triangles = [[0, 2, 3], [0, 3, 4], [0, 5, 6], [0, 6, 1]]
            .into_iter()
            .enumerate()
            .map(|(index, vertices)| {
                InPlaneTriangle::new(&mesh, projection, index, vertices).expect("triangle")
            })
            .collect();
        (
            mesh,
            InPlaneSet::new(triangles).expect("set"),
            frame,
            projection,
        )
    }

    fn regions(protected: &BTreeSet<u32>) -> Vec<InPlaneRegion> {
        let (mesh, set, frame, projection) = pinched_fan();
        build_regions(
            &mesh,
            frame,
            projection,
            &set,
            protected,
            RegionLimits::new(10, 100, 100),
            &ExecutionOptions::new(Tolerance::METRE),
        )
        .expect("regions")
    }

    fn ring(points: &[(f64, f64)]) -> Vec<Point2> {
        points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    #[test]
    fn a_pinch_gives_two_regions_and_drops_collinear_vertices() {
        // The boundary of the left group is c, r6, r7, r0 = (1,1), (0,2),
        // (0,1), (0,0); (0,1) lies between its neighbours and goes. The ring
        // is counter-clockwise from its smallest point (0,0). The right group
        // is (1,1), (2,0), (2,1), (2,2) without (2,1). Sorted by outer ring.
        let found = regions(&BTreeSet::new());
        assert_eq!(found.len(), 2);
        assert_eq!(found[0].side, InPlaneSide::Below);
        assert_eq!(found[0].outer, ring(&[(0.0, 0.0), (1.0, 1.0), (0.0, 2.0)]));
        assert_eq!(found[1].outer, ring(&[(1.0, 1.0), (2.0, 0.0), (2.0, 2.0)]));
        assert!(found.iter().all(|region| region.holes.is_empty()));
    }

    #[test]
    fn a_protected_collinear_vertex_is_kept() {
        let found = regions(&BTreeSet::from([6]));
        assert_eq!(
            found[0].outer,
            ring(&[(0.0, 0.0), (1.0, 1.0), (0.0, 2.0), (0.0, 1.0)])
        );
        assert_eq!(found[1].outer, ring(&[(1.0, 1.0), (2.0, 0.0), (2.0, 2.0)]));
    }
}
