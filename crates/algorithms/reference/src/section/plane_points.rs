//! Exact points of a section plane, for deciding cut rings (ADR 0092, #377).
//!
//! A cut ring node is either a source vertex lying in the plane or the
//! crossing of a source edge with it. Neither has to be a binary64 point, and
//! the coordinates a section returns are rounded: a crossing is interpolated,
//! and the projection to a turned frame multiplies and adds. Rings that touch
//! in the source can cross by an ulp in those coordinates, so every decision
//! about rings (crossing, contact, containment, outer or hole) is made on
//! exact points (`axiolid_exact_arrangement::ExactPoints`, ADR 0094), and the
//! returned coordinates take no part in it.
//!
//! The plane is the one `ExactSectionPlane` classifies against: `O = origin`,
//! `X = x_point`, `Y = y_point`. Let `N = (X - O) x (Y - O)` and
//! `D(Q) = N . (Q - O)`, which is `det[X - O, Y - O, Q - O]` and has degree 3
//! in the source coordinates. `D` vanishes exactly on the plane. A point is
//! stored in homogeneous coordinates `(x, y, w)` with `w > 0`, in the
//! projection that drops one axis (`Projection`), so the projected point is
//! `(x / w, y / w)`:
//!
//! - a source vertex `P` with `D(P) = 0`: `(P_u, P_v, 1)`;
//! - the crossing of the edge `A`-`B`, whose ends have `D` of strictly
//!   opposite signs: `x = D(A) B_u - D(B) A_u`, `y` the same in `v`, and
//!   `w = D(A) - D(B)`, all negated when `w < 0`. The degrees are 4, 4 and 3.
//!   When `x = q w` and `y = r w` hold exactly for binary64 `q` and `r`, the
//!   crossing is stored as `(q, r, 1)`, the same point, so that it takes the
//!   fast path of the predicates.

use std::collections::BTreeMap;

use axiolid_contracts::{GeomResult, Sign};
use axiolid_core::{Point2, Point3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_exact_arrangement::{ExactPoints, PointId};
pub(super) use axiolid_exact_arrangement::{Location, RingBox};
use axiolid_mesh::TriMesh;

use super::regions::Projection;
use super::{
    internal_topology_error, mesh_index, Classification, EdgeKey, ExactSectionPlane, NodeKey,
};
use crate::coplanar::project;

/// The exact points of one section, interned by cut node. The predicates are
/// those of `ExactPoints`, on plain indices.
pub(super) struct PlanePoints {
    origin: [Dyadic; 3],
    normal: [Dyadic; 3],
    axis: usize,
    points: ExactPoints,
    by_node: BTreeMap<NodeKey, usize>,
    /// `D` of the mesh vertices met so far.
    values: BTreeMap<u32, Dyadic>,
}

fn exact(value: f64) -> GeomResult<Dyadic> {
    Dyadic::try_from_f64(value).ok_or_else(|| {
        axiolid_contracts::GeomError::Degenerate("cut ring point is not finite".into())
    })
}

fn exact_point(point: Point3) -> GeomResult<[Dyadic; 3]> {
    Ok([exact(point.x)?, exact(point.y)?, exact(point.z)?])
}

fn dot(left: &[Dyadic; 3], right: &[Dyadic; 3]) -> Dyadic {
    left[0]
        .mul(&right[0])
        .add(&left[1].mul(&right[1]))
        .add(&left[2].mul(&right[2]))
}

/// `p_i q_j - q_i p_j`.
fn minor(p: &[Dyadic; 3], q: &[Dyadic; 3], i: usize, j: usize) -> Dyadic {
    p[i].mul(&q[j]).sub(&q[i].mul(&p[j]))
}

fn sub3(left: &[Dyadic; 3], right: &[Dyadic; 3]) -> [Dyadic; 3] {
    [
        left[0].sub(&right[0]),
        left[1].sub(&right[1]),
        left[2].sub(&right[2]),
    ]
}

impl PlanePoints {
    pub(super) fn new(plane: &ExactSectionPlane, projection: Projection) -> GeomResult<Self> {
        let origin = exact_point(plane.origin)?;
        let to_x = sub3(&exact_point(plane.x_point)?, &origin);
        let to_y = sub3(&exact_point(plane.y_point)?, &origin);
        let normal = [
            minor(&to_x, &to_y, 1, 2),
            minor(&to_x, &to_y, 2, 0),
            minor(&to_x, &to_y, 0, 1),
        ];
        Ok(Self {
            origin,
            normal,
            axis: projection.axis(),
            points: ExactPoints::with_frame_sign(projection.frame_sign()),
            by_node: BTreeMap::new(),
            values: BTreeMap::new(),
        })
    }

    /// The store, for the arrangement.
    pub(super) fn store(&self) -> &ExactPoints {
        &self.points
    }

    /// The store, for the arrangement.
    pub(super) fn store_mut(&mut self) -> &mut ExactPoints {
        &mut self.points
    }

    /// `D(point)`: zero on the plane.
    fn plane_value(&self, point: Point3) -> GeomResult<Dyadic> {
        let offset = sub3(&exact_point(point)?, &self.origin);
        Ok(dot(&self.normal, &offset))
    }

    /// `D` of a mesh vertex, computed once.
    fn vertex_value(&mut self, index: u32, point: Point3) -> GeomResult<Dyadic> {
        if let Some(value) = self.values.get(&index) {
            return Ok(value.clone());
        }
        let value = self.plane_value(point)?;
        self.values.insert(index, value.clone());
        Ok(value)
    }

    /// The exact point of a cut node, interned. `classifications` are the
    /// certified plane sides of the mesh vertices.
    pub(super) fn node(
        &mut self,
        mesh: &TriMesh,
        classifications: &[Classification],
        node: NodeKey,
    ) -> GeomResult<usize> {
        if let Some(&index) = self.by_node.get(&node) {
            return Ok(index);
        }
        let position = |index: u32| {
            mesh.positions
                .get(mesh_index(index)?)
                .copied()
                .ok_or_else(internal_topology_error)
        };
        let index = match node {
            NodeKey::Vertex(vertex) => {
                // The classification is certified exact (`orient3d`) and is
                // the guard here; the plane value is only re-derived in
                // debug builds.
                let certified = classifications
                    .get(mesh_index(vertex)?)
                    .ok_or_else(internal_topology_error)?;
                if certified.sign != Sign::Zero {
                    return Err(internal_topology_error());
                }
                let point = position(vertex)?;
                debug_assert_eq!(self.plane_value(point)?.sign(), Some(Sign::Zero));
                self.push_vertex(point)?
            }
            NodeKey::Edge(EdgeKey(left, right)) => {
                let (a, b) = (position(left)?, position(right)?);
                let (da, db) = (self.vertex_value(left, a)?, self.vertex_value(right, b)?);
                self.crossing_of(&da, &db, a, b)?
            }
        };
        self.by_node.insert(node, index);
        Ok(index)
    }

    /// A source vertex lying in the plane.
    #[cfg(test)]
    pub(super) fn vertex(&mut self, point: Point3) -> GeomResult<usize> {
        if self.plane_value(point)?.sign() != Some(Sign::Zero) {
            return Err(internal_topology_error());
        }
        self.push_vertex(point)
    }

    fn push_vertex(&mut self, point: Point3) -> GeomResult<usize> {
        Ok(self.points.push_flat(project(point, self.axis))?.0)
    }

    /// The crossing of the edge `a`-`b` with the plane; the ends lie strictly
    /// on opposite sides.
    #[cfg(test)]
    pub(super) fn crossing(&mut self, a: Point3, b: Point3) -> GeomResult<usize> {
        let (da, db) = (self.plane_value(a)?, self.plane_value(b)?);
        self.crossing_of(&da, &db, a, b)
    }

    /// The crossing of `a`-`b`, whose plane values are `da` and `db`.
    fn crossing_of(&mut self, da: &Dyadic, db: &Dyadic, a: Point3, b: Point3) -> GeomResult<usize> {
        let opposite = matches!(
            (da.sign(), db.sign()),
            (Some(Sign::Positive), Some(Sign::Negative))
                | (Some(Sign::Negative), Some(Sign::Positive))
        );
        if !opposite {
            return Err(internal_topology_error());
        }
        let (a, b): (Point2, Point2) = (project(a, self.axis), project(b, self.axis));
        let point = [
            da.mul(&exact(b.x)?).sub(&db.mul(&exact(a.x)?)),
            da.mul(&exact(b.y)?).sub(&db.mul(&exact(a.y)?)),
            da.sub(db),
        ];
        // An edge that does not move in a coordinate crosses at that value.
        let fixed = |from: f64, to: f64| (from == to).then_some(from);
        let hints = [fixed(a.x, b.x), fixed(a.y, b.y)];
        Ok(self.points.push_hinted(point, hints)?.0)
    }

    /// The parameter along the first of two crossing source edges at which
    /// they cross (`ExactPoints::segment_parameter`). For output only.
    pub(super) fn segment_parameter(&self, ends: [usize; 4]) -> GeomResult<f64> {
        self.points.segment_parameter(ends.map(PointId))
    }

    /// `[u_lo, v_lo, u_hi, v_hi]`, a sound box of the projected point.
    pub(super) fn bounds(&self, index: usize) -> [f64; 4] {
        self.points.bounds(PointId(index))
    }

    /// The projected position of a point that has one in binary64.
    #[cfg(test)]
    pub(super) fn flat(&self, index: usize) -> Option<Point2> {
        self.points.flat(PointId(index))
    }

    /// The exact homogeneous coordinates of a point.
    pub(super) fn exact(&self, index: usize) -> &[Dyadic; 3] {
        self.points.exact(PointId(index))
    }

    /// Exact orientation of `a`, `b`, `c` in the frame's own `(x, y)` sense:
    /// positive when they turn counter-clockwise there.
    pub(super) fn orient(&self, a: usize, b: usize, c: usize) -> GeomResult<Sign> {
        self.points.orient(PointId(a), PointId(b), PointId(c))
    }

    /// Order of the projected coordinates of `a` and `b` along `axis`
    /// (0 = u, 1 = v). The projection's own axes, not the frame's.
    pub(super) fn compare(
        &self,
        a: usize,
        b: usize,
        axis: usize,
    ) -> GeomResult<std::cmp::Ordering> {
        self.points.compare(PointId(a), PointId(b), axis)
    }

    /// Whether `a` and `b` are the same point.
    pub(super) fn same(&self, a: usize, b: usize) -> GeomResult<bool> {
        self.points.same(PointId(a), PointId(b))
    }

    /// Order by `u`, then `v`.
    pub(super) fn lexicographic(&self, a: usize, b: usize) -> GeomResult<std::cmp::Ordering> {
        self.points.lexicographic(PointId(a), PointId(b))
    }

    /// The axis of the projection along which two distinct points of a line
    /// differ.
    pub(super) fn dominant(&self, a: usize, b: usize) -> GeomResult<usize> {
        self.points.dominant(PointId(a), PointId(b))
    }

    /// `(low, high)` of two points along `axis`.
    pub(super) fn ordered(&self, a: usize, b: usize, axis: usize) -> GeomResult<(usize, usize)> {
        let (low, high) = self.points.ordered(PointId(a), PointId(b), axis)?;
        Ok((low.0, high.0))
    }

    /// Whether collinear `at` lies strictly between `before` and `after`.
    pub(super) fn strictly_between(
        &self,
        before: usize,
        at: usize,
        after: usize,
    ) -> GeomResult<bool> {
        self.points
            .strictly_between(PointId(before), PointId(at), PointId(after))
    }

    /// Whether collinear `point` lies on the closed segment `start`-`end`.
    pub(super) fn within_segment(
        &self,
        point: usize,
        start: usize,
        end: usize,
    ) -> GeomResult<bool> {
        self.points
            .within_segment(PointId(point), PointId(start), PointId(end))
    }

    /// Where `point` lies relative to the closed `ring`.
    pub(super) fn locate(&self, point: usize, ring: &[usize]) -> GeomResult<Location> {
        let ring: Vec<PointId> = ring.iter().copied().map(PointId).collect();
        self.points.locate(PointId(point), &ring)
    }

    /// The boxes of a ring of exact points.
    pub(super) fn ring_box(&self, ring: &[usize]) -> RingBox {
        let ring: Vec<PointId> = ring.iter().copied().map(PointId).collect();
        self.points.ring_box(&ring)
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use axiolid_contracts::GeomError;
    use axiolid_core::{Frame3, Point2, Vec3};

    use super::super::{classify_vertices, node_point, orient_sign, project_local};
    use super::*;

    fn plane_of(frame: Frame3) -> (ExactSectionPlane, Projection) {
        let plane = ExactSectionPlane::new(frame).expect("plane");
        let projection = Projection::new(&plane).expect("projection");
        (plane, projection)
    }

    fn flat_frame() -> Frame3 {
        Frame3 {
            origin: Point3::new(0.0, 0.0, 0.5),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }
    }

    fn turned(k: u32) -> Frame3 {
        let (s, c) = (std::f64::consts::PI * f64::from(k) / 64.0).sin_cos();
        Frame3 {
            origin: Point3::new(0.0, 0.0, 0.5),
            x: Vec3::new(c, s, 0.0),
            y: Vec3::new(-s, c, 0.0),
            z: Vec3::Z,
        }
    }

    /// The vertical edge at `(x, 0.1)` across the plane z = 0.5.
    fn column(x: f64) -> (Point3, Point3) {
        (Point3::new(x, 0.1, 0.0), Point3::new(x, 0.1, 1.0))
    }

    #[test]
    fn crossings_on_one_line_stay_collinear_where_rounded_points_do_not() {
        let xs = [0.3, 0.9, 1.7];
        let mut premise = false;
        for k in 1..64 {
            let frame = turned(k);
            let (plane, projection) = plane_of(frame);
            let mut points = PlanePoints::new(&plane, projection).expect("points");
            let nodes: Vec<usize> = xs
                .iter()
                .map(|&x| {
                    let (a, b) = column(x);
                    points.crossing(a, b).expect("crossing")
                })
                .collect();
            assert_eq!(
                points.orient(nodes[0], nodes[1], nodes[2]).expect("orient"),
                Sign::Zero,
                "frame {k}"
            );
            // The same three points as the section returns them.
            let positions: Vec<Point3> = xs
                .iter()
                .flat_map(|&x| {
                    let (a, b) = column(x);
                    [a, b]
                })
                .collect();
            let mesh = TriMesh::new(positions, vec![]);
            let classifications = classify_vertices(&mesh, plane).expect("classes");
            let rounded: Vec<_> = (0..3u32)
                .map(|i| {
                    let node = NodeKey::Edge(EdgeKey::new(2 * i, 2 * i + 1));
                    let world = node_point(&mesh, &classifications, node).expect("point");
                    project_local(frame, world).expect("local")
                })
                .collect();
            premise |= orient_sign(rounded[0], rounded[1], rounded[2]).expect("sign") != Sign::Zero;
        }
        assert!(premise, "no frame rounds the collinear points off the line");
    }

    #[test]
    fn a_vertex_and_a_crossing_at_the_same_point_are_the_same() {
        let (plane, projection) = plane_of(flat_frame());
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        let crossing = points
            .crossing(Point3::new(0.25, 0.5, 1.0), Point3::new(0.25, 0.5, -1.0))
            .expect("crossing");
        let vertex = points.vertex(Point3::new(0.25, 0.5, 0.5)).expect("vertex");
        assert!(points.same(crossing, vertex).expect("same"));
        let slanted = points
            .crossing(Point3::new(0.0, 0.0, -1.0), Point3::new(1.0, 2.0, 3.0))
            .expect("crossing");
        // The edge meets z = 0.5 at t = 0.375: (0.375, 0.75).
        let other = points
            .vertex(Point3::new(0.375, 0.75, 0.5))
            .expect("vertex");
        assert!(points.same(slanted, other).expect("same"));
        assert!(!points.same(slanted, vertex).expect("same"));
    }

    #[test]
    fn a_crossing_with_a_binary64_position_takes_the_fast_path() {
        let (plane, projection) = plane_of(flat_frame());
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        // z = 0.5 at t = 0.375: (0.375, 0.75), exact in binary64.
        let crossing = points
            .crossing(Point3::new(0.0, 0.0, -1.0), Point3::new(1.0, 2.0, 3.0))
            .expect("crossing");
        assert_eq!(points.flat(crossing), Some(Point2::new(0.375, 0.75)));
        let vertex = points
            .vertex(Point3::new(0.375, 0.75, 0.5))
            .expect("vertex");
        assert!(points.same(crossing, vertex).expect("same"));
        // z = 0.5 at t = 1/3: not representable, so it stays homogeneous.
        let third = points
            .crossing(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.5))
            .expect("crossing");
        assert_eq!(points.flat(third), None);
        assert!(!points.same(third, vertex).expect("same"));
        // A neighbour of the binary64 position is not that position.
        let near = points
            .vertex(Point3::new(0.375_f64.next_up(), 0.75, 0.5))
            .expect("vertex");
        assert!(!points.same(crossing, near).expect("same"));
        assert_eq!(
            points.lexicographic(crossing, near).expect("order"),
            Ordering::Less
        );
    }

    #[test]
    fn a_crossing_does_not_depend_on_the_direction_of_its_edge() {
        let (plane, projection) = plane_of(turned(5));
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        let (a, b) = (Point3::new(0.1, 0.7, -0.3), Point3::new(-0.9, 0.2, 1.4));
        let forward = points.crossing(a, b).expect("crossing");
        let backward = points.crossing(b, a).expect("crossing");
        assert!(points.same(forward, backward).expect("same"));
        assert_eq!(points.exact(forward), points.exact(backward));
    }

    #[test]
    fn the_box_of_a_crossing_contains_its_position() {
        let frame = Frame3 {
            origin: Point3::ZERO,
            ..flat_frame()
        };
        let (plane, projection) = plane_of(frame);
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        // z runs from -1 to 2 along the edge, so z = 0 at t = 1/3.
        let crossing = points
            .crossing(Point3::new(0.0, 0.0, -1.0), Point3::new(1.0, 1.0, 2.0))
            .expect("crossing");
        let [u_lo, v_lo, u_hi, v_hi] = points.bounds(crossing);
        let third = 1.0 / 3.0;
        assert!(u_lo <= third && third <= u_hi && v_lo <= third && third <= v_hi);
        assert!(u_hi - u_lo < 1e-12);
    }

    #[test]
    fn a_mirrored_frame_reverses_the_orientation() {
        let counter = flat_frame();
        let mirrored = Frame3 {
            y: -Vec3::Y,
            z: -Vec3::Z,
            ..counter
        };
        for (frame, expected) in [(counter, Sign::Positive), (mirrored, Sign::Negative)] {
            let (plane, projection) = plane_of(frame);
            let mut points = PlanePoints::new(&plane, projection).expect("points");
            let nodes: Vec<usize> = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
                .into_iter()
                .map(|(x, y)| points.vertex(Point3::new(x, y, 0.5)).expect("vertex"))
                .collect();
            assert_eq!(
                points.orient(nodes[0], nodes[1], nodes[2]).expect("orient"),
                expected
            );
        }
    }

    /// Crossings of slanted edges with the plane z = 0.5, on the line
    /// y = 0.1 at `x = 2 * t` for the given `t`, none of them binary64.
    fn slanted(points: &mut PlanePoints, thirds: &[f64]) -> Vec<usize> {
        thirds
            .iter()
            .map(|&x| {
                points
                    .crossing(Point3::new(0.0, 0.1, 0.0), Point3::new(x, 0.1, 1.5))
                    .expect("crossing")
            })
            .collect()
    }

    #[test]
    fn homogeneous_points_decide_collinearity_exactly_in_either_frame_sense() {
        let counter = flat_frame();
        let mirrored = Frame3 {
            y: -Vec3::Y,
            z: -Vec3::Z,
            ..counter
        };
        for (frame, turn) in [(counter, Sign::Positive), (mirrored, Sign::Negative)] {
            let (plane, projection) = plane_of(frame);
            let mut points = PlanePoints::new(&plane, projection).expect("points");
            // On y = 0.1 at x = 1/3, 2/3 and 4/3.
            let line = slanted(&mut points, &[1.0, 2.0, 4.0]);
            assert!(line.iter().all(|&i| points.flat(i).is_none()));
            assert_eq!(
                points.orient(line[0], line[1], line[2]).expect("orient"),
                Sign::Zero
            );
            // (1/3, 0.1), (2/3, 0.1), (0.5, 0.5) turn left in the frame.
            let above = points.vertex(Point3::new(0.5, 0.5, 0.5)).expect("vertex");
            assert_eq!(
                points.orient(line[0], line[1], above).expect("orient"),
                turn
            );
        }
    }

    #[test]
    fn equal_homogeneous_points_compare_equal() {
        let (plane, projection) = plane_of(flat_frame());
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        // (1/3, 0.1) from two different edges.
        let first = slanted(&mut points, &[1.0])[0];
        let second = points
            .crossing(Point3::new(0.0, 0.1, 0.0), Point3::new(2.0, 0.1, 3.0))
            .expect("crossing");
        assert!(points.same(first, second).expect("same"));
        assert_eq!(
            points.lexicographic(first, second).expect("order"),
            Ordering::Equal
        );
    }

    #[test]
    fn an_edge_that_does_not_move_crosses_at_its_own_coordinates() {
        let (plane, projection) = plane_of(flat_frame());
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        let column = points
            .crossing(Point3::new(0.25, 0.5, 1.0), Point3::new(0.25, 0.5, -1.0))
            .expect("crossing");
        assert_eq!(points.flat(column), Some(Point2::new(0.25, 0.5)));
    }

    #[test]
    fn a_vertex_node_needs_a_certified_plane_side_of_zero() {
        let (plane, projection) = plane_of(flat_frame());
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        let mesh = TriMesh::new(vec![Point3::new(0.0, 0.0, 1.0)], vec![]);
        let classifications = classify_vertices(&mesh, plane).expect("classes");
        let off = points.node(&mesh, &classifications, NodeKey::Vertex(0));
        assert!(matches!(
            off,
            Err(GeomError::BackendContractViolation { .. })
        ));
    }

    #[test]
    fn a_vertex_off_the_plane_is_a_contract_violation() {
        let (plane, projection) = plane_of(flat_frame());
        let mut points = PlanePoints::new(&plane, projection).expect("points");
        let off = points.vertex(Point3::new(0.0, 0.0, 0.5 + f64::EPSILON));
        assert!(matches!(
            off,
            Err(GeomError::BackendContractViolation { .. })
        ));
        let same_side = points.crossing(Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 2.0));
        assert!(matches!(
            same_side,
            Err(GeomError::BackendContractViolation { .. })
        ));
    }
}
