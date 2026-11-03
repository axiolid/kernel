//! In-plane regions of a section: faces lying in the cutting plane.
//!
//! The in-plane triangles may touch, abut and overlap. `arrangement` cuts
//! their source edges into pieces on exact points and says what covers each
//! side of a piece: `Below`, `Above`, `Both` or nothing. A region is a
//! maximal part of the plane of one class, so regions have pairwise disjoint
//! interiors, faces of one side merge whether they overlap or only share a
//! boundary segment, and `Both` is where faces of both sides coincide. The
//! pipeline, in order:
//!
//! 1. Boundary half-edges: a piece whose two sides differ is a half-edge of
//!    each non-empty class, directed with that class on its left.
//! 2. Fans: the pieces at a vertex are ordered counter-clockwise in the
//!    frame's sense by exact signs. The half-edge that follows an incoming
//!    one is the first boundary of its class met turning clockwise from the
//!    way it came, so a vertex where a region pinches is resolved per fan.
//! 3. Cycles: following half-edges closes rings. Each is an outer ring or a
//!    hole by the exact turn at its extreme vertex.
//! 4. Collinear vertices that only two boundary pieces meet at are dropped;
//!    the others, and the vertices of cut segments, stay.
//! 5. Each hole belongs to the innermost outer ring of its class that holds
//!    it, decided with exact point location.
//! 6. The rings are projected to the plane frame once per vertex, so the
//!    regions of all classes share bit-identical points, and made canonical.
//!
//! Only the returned coordinates are rounded. A ring that collapses under
//! that rounding is refused as `Degenerate`.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};
use axiolid_core::{Frame3, Point2, Point3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{InPlaneRegion, InPlaneSide, RegionLimits};

pub(super) use axiolid_exact_arrangement::canonical_rotation;
use axiolid_exact_arrangement::{canonical_start, trace, HalfEdge, PairBudget, PieceId, VertexId};

use super::arrangement::{Class, InPlaneArrangement};
use super::plane_points::PlanePoints;
use super::ring_nesting::TracedRing;
use super::{
    internal_topology_error, mesh_index, orient_sign, sequence_order, sign_times, ExactSectionPlane,
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

    /// The dropped axis: 0 = x, 1 = y, 2 = z.
    pub(super) fn axis(self) -> usize {
        self.axis
    }

    /// The exact orientation of the frame's own axes after the drop.
    pub(super) fn frame_sign(self) -> Sign {
        self.frame_sign
    }

    /// Exact orientation of three points of the plane in the frame's own
    /// `(x, y)` sense: positive when they turn counter-clockwise there.
    pub(super) fn orient(self, a: Point3, b: Point3, c: Point3) -> GeomResult<Sign> {
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
}

pub(super) fn corners(mesh: &TriMesh, vertices: [u32; 3]) -> GeomResult<[Point3; 3]> {
    let mut points = [Point3::ZERO; 3];
    for (slot, vertex) in points.iter_mut().zip(vertices) {
        *slot = mesh.positions[mesh_index(vertex)?];
    }
    Ok(points)
}

/// Whether `point` lies strictly between `a` and `b` on their common line.
pub(super) fn strictly_between(a: Point3, point: Point3, b: Point3) -> bool {
    let within = |a: f64, p: f64, b: f64| a.min(b) <= p && p <= a.max(b);
    point != a
        && point != b
        && within(a.x, point.x, b.x)
        && within(a.y, point.y, b.y)
        && within(a.z, point.z, b.z)
}

/// The regions of a section and what the cut needs to know about them.
#[derive(Debug)]
pub(super) struct Regions {
    pub(super) regions: Vec<InPlaneRegion>,
    /// Source vertices that remain on a ring after collinear ones are dropped.
    pub(super) ring_vertices: BTreeSet<u32>,
    /// The rings of the `Both` regions on their exact points, aligned with
    /// the returned points.
    pub(super) both_rings: Vec<TracedRing>,
    /// Whether a ring of a `Both` region visits a vertex twice.
    pub(super) both_pinched: bool,
    /// Candidate pairs examined: the arrangement's and the hole containment's.
    pub(super) pairs: usize,
}

fn collapsed() -> GeomError {
    GeomError::Degenerate("in-plane region collapsed under projection".into())
}

/// Regions of the arranged in-plane triangles: one per maximal part of a
/// class, with its holes.
///
/// `protected` vertices (nodes of cut segments) are never simplified away,
/// so a chain end stays a vertex of the ring it meets.
pub(super) fn build_regions(
    mesh: &TriMesh,
    frame: Frame3,
    arrangement: &InPlaneArrangement,
    exact: &PlanePoints,
    protected: &BTreeSet<u32>,
    limits: RegionLimits,
    options: &ExecutionOptions,
) -> GeomResult<Regions> {
    let mut budget = PairBudget::new(limits.max_in_plane_pairs, "section in-plane pairs");
    budget.used = arrangement.pairs;

    // Boundary half-edges per class, with the class on their left.
    let mut half_edges = BTreeMap::<Class, Vec<HalfEdge>>::new();
    let mut on_boundary = vec![0usize; arrangement.vertices.len()];
    for (index, piece) in arrangement.pieces.iter().enumerate() {
        if piece.left == piece.right {
            continue;
        }
        let [from, to] = piece.ends;
        on_boundary[from] += 1;
        on_boundary[to] += 1;
        let piece = PieceId(index);
        if arrangement.pieces[index].left != Class::Empty {
            half_edges
                .entry(arrangement.pieces[index].left)
                .or_default()
                .push(HalfEdge {
                    piece,
                    forward: true,
                });
        }
        if arrangement.pieces[index].right != Class::Empty {
            half_edges
                .entry(arrangement.pieces[index].right)
                .or_default()
                .push(HalfEdge {
                    piece,
                    forward: false,
                });
        }
    }
    // A vertex stays on a ring where other boundary meets or a cut ends.
    let keep = |vertex: VertexId| {
        arrangement.vertices[vertex.0]
            .sources
            .iter()
            .any(|source| protected.contains(source))
            || on_boundary[vertex.0] != 2
    };

    // Output points, rounded once per vertex.
    let mut points = BTreeMap::<usize, Point2>::new();
    let mut point_of = |vertex: usize| -> GeomResult<Point2> {
        if let Some(&point) = points.get(&vertex) {
            return Ok(point);
        }
        let point = arrangement.output_point(mesh, frame, exact, vertex)?;
        points.insert(vertex, point);
        Ok(point)
    };
    // A ring as plane-local points and the exact points they come from, in
    // canonical rotation.
    let mut traced = |ring: &[VertexId]| -> GeomResult<TracedRing> {
        let local = ring
            .iter()
            .map(|&vertex| point_of(vertex.0))
            .collect::<GeomResult<Vec<_>>>()?;
        if local.len() < 3
            || (0..local.len()).any(|index| local[index] == local[(index + 1) % local.len()])
        {
            return Err(collapsed());
        }
        let start = canonical_start(&local);
        let nodes: Vec<usize> = ring
            .iter()
            .map(|&vertex| arrangement.vertices[vertex.0].exact)
            .collect();
        Ok(TracedRing {
            nodes: nodes[start..]
                .iter()
                .chain(&nodes[..start])
                .copied()
                .collect(),
            points: canonical_rotation(local),
        })
    };

    let mut ring_vertices = BTreeSet::new();
    let mut regions = Vec::new();
    let mut both_rings = Vec::new();
    let mut both_pinched = false;
    let mut vertices = 0usize;
    for (class, class_half_edges) in &half_edges {
        let Some(side) = class.side() else {
            return Err(internal_topology_error());
        };
        let (rings, class_regions) = trace(
            exact.store(),
            &arrangement.cut,
            class_half_edges,
            &keep,
            &mut budget,
            options,
        )?;
        for ring in &rings {
            for &vertex in &ring.vertices {
                ring_vertices.extend(arrangement.vertices[vertex.0].sources.iter().copied());
            }
        }
        for region in &class_regions {
            if regions.len() >= limits.max_regions {
                return Err(GeomError::BudgetExceeded {
                    resource: "section in-plane regions",
                });
            }
            let outer = traced(&rings[region.outer].vertices)?;
            let mut holes = Vec::new();
            for &hole in &region.holes {
                holes.push(traced(&rings[hole].vertices)?);
            }
            holes.sort_by(|left, right| sequence_order(&left.points, &right.points));
            vertices +=
                outer.points.len() + holes.iter().map(|ring| ring.points.len()).sum::<usize>();
            if vertices > limits.max_region_vertices {
                return Err(GeomError::BudgetExceeded {
                    resource: "section in-plane region vertices",
                });
            }
            if side == InPlaneSide::Both {
                for ring in std::iter::once(&outer).chain(&holes) {
                    let distinct: BTreeSet<usize> = ring.nodes.iter().copied().collect();
                    both_pinched |= distinct.len() != ring.nodes.len();
                }
                both_rings.push(TracedRing {
                    nodes: outer.nodes.clone(),
                    points: outer.points.clone(),
                });
                for ring in &holes {
                    both_rings.push(TracedRing {
                        nodes: ring.nodes.clone(),
                        points: ring.points.clone(),
                    });
                }
            }
            regions.push(InPlaneRegion::new(
                side,
                outer.points,
                holes.into_iter().map(|ring| ring.points).collect(),
            ));
        }
    }
    regions.sort_by(|left, right| {
        left.side
            .cmp(&right.side)
            .then_with(|| sequence_order(&left.outer, &right.outer))
    });
    Ok(Regions {
        regions,
        ring_vertices,
        both_rings,
        both_pinched,
        pairs: budget.used,
    })
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Tolerance, Vec3};

    use super::super::classify_vertices;
    use super::*;

    /// Eight triangles fan around `c = (1,1,0)` over `[0,2]^2`; the edge
    /// midpoints `(1,0)` and `(1,2)` are lifted to `z = 1`, so the four
    /// triangles touching them leave the plane. The remaining four form two
    /// groups that meet only at `c`: a pinch. Only the in-plane triangles
    /// are given; their missing neighbours make the boundary.
    fn built(protected: &BTreeSet<u32>) -> Regions {
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
        let classifications = classify_vertices(&mesh, plane).expect("classes");
        let triangles: Vec<InPlaneTriangle> = [[0, 2, 3], [0, 3, 4], [0, 5, 6], [0, 6, 1]]
            .into_iter()
            .enumerate()
            .map(|(index, vertices)| {
                InPlaneTriangle::new(&mesh, projection, index, vertices).expect("triangle")
            })
            .collect();
        let options = ExecutionOptions::new(Tolerance::METRE);
        let mut exact = PlanePoints::new(&plane, projection).expect("points");
        let arrangement = InPlaneArrangement::build(
            &mesh,
            &classifications,
            &triangles,
            &mut exact,
            1000,
            &options,
        )
        .expect("arrangement");
        build_regions(
            &mesh,
            frame,
            &arrangement,
            &exact,
            protected,
            RegionLimits::new(10, 100, 100),
            &options,
        )
        .expect("regions")
    }

    fn regions(protected: &BTreeSet<u32>) -> Vec<InPlaneRegion> {
        built(protected).regions
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
    fn the_kept_ring_vertices_exclude_the_dropped_collinear_one() {
        // Vertex 6 is (0,1), the collinear one on the left group.
        assert!(!built(&BTreeSet::new()).ring_vertices.contains(&6));
        assert!(built(&BTreeSet::from([6])).ring_vertices.contains(&6));
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
