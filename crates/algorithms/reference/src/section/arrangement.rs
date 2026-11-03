//! The exact arrangement of the in-plane faces of a section (ADR 0093, #366).
//!
//! The in-plane triangles may touch, abut and overlap. Their source edges are
//! cut at every point where they meet into pieces that no other piece
//! crosses, and each piece records which class of material covers each of its
//! sides. Every decision is an exact sign on `PlanePoints`; the returned
//! coordinates only label the output.
//!
//! 1. Cutting: the unique in-plane source edges are the segments of
//!    `axiolid_exact_arrangement::arrange` (ADR 0094), which joins equal
//!    points, splits at crossings, T-junctions and collinear overlaps, and
//!    returns the pieces.
//! 2. Coverage: a triangle that carries a piece covers the side its third
//!    vertex is on; every other triangle either contains the piece in its
//!    interior, covering both sides, or does not. A triangle is convex, so a
//!    piece lies inside it or outside, never partly, and the pieces need no
//!    face tracing.
//!
//! A side is `Below` when only faces whose material lies below the plane
//! cover it, `Above` when only the others do, `Both` when both kinds do and
//! `Empty` when none does.
//!
//! Candidate pairs of segments whose boxes overlap, and of pieces with the
//! triangles whose boxes can hold them, count against
//! `RegionLimits::max_in_plane_pairs`. The crossing points are bounded by the
//! pairs and charged against the memory budget as they are created.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomResult, Sign};
use axiolid_core::{Frame3, Point2};
use axiolid_exact_arrangement::{arrange, union_box, Arrangement, PairBudget, PointId, Segment};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::InPlaneSide;

use super::plane_points::PlanePoints;
use super::regions::InPlaneTriangle;
use super::{
    internal_topology_error, mesh_index, project_local, Classification, EdgeKey, NodeKey,
    REQUEST_SCRATCH_PER_ELEMENT,
};

/// What covers one side of a piece.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Class {
    Empty,
    Below,
    Above,
    Both,
}

impl Class {
    fn of(below: bool, above: bool) -> Self {
        match (below, above) {
            (false, false) => Self::Empty,
            (true, false) => Self::Below,
            (false, true) => Self::Above,
            (true, true) => Self::Both,
        }
    }

    /// The side of a region of this class; `None` for `Empty`.
    pub(super) fn side(self) -> Option<InPlaneSide> {
        match self {
            Self::Empty => None,
            Self::Below => Some(InPlaneSide::Below),
            Self::Above => Some(InPlaneSide::Above),
            Self::Both => Some(InPlaneSide::Both),
        }
    }
}

/// A point of the arrangement: a source vertex, or where source edges cross.
#[derive(Debug, Clone)]
pub(super) struct Vertex {
    /// An exact point of `PlanePoints` at this position.
    pub(super) exact: usize,
    /// The source vertices at this position, ascending.
    pub(super) sources: Vec<u32>,
    /// For a point without source vertices, the crossing that labels it: the
    /// smallest pair of edges, and the exact points of their ends.
    crossing: Option<((EdgeKey, EdgeKey), [usize; 4])>,
}

/// A stretch of source edges no other edge crosses or touches inside.
#[derive(Debug, Clone, Copy)]
pub(super) struct Piece {
    /// Vertex numbers, `ends[0] < ends[1]`.
    pub(super) ends: [usize; 2],
    /// What covers the left of the piece going from `ends[0]` to `ends[1]`.
    pub(super) left: Class,
    /// What covers its right.
    pub(super) right: Class,
}

/// The pieces of the in-plane faces and what covers each side of them.
#[derive(Debug)]
pub(super) struct InPlaneArrangement {
    pub(super) vertices: Vec<Vertex>,
    pub(super) pieces: Vec<Piece>,
    /// The pieces of each source edge from `key.0` to `key.1`, with whether
    /// that direction runs from `ends[0]` to `ends[1]` of the piece.
    edges: BTreeMap<EdgeKey, Vec<(usize, bool)>>,
    /// The arrangement the vertices and pieces were taken from, for tracing.
    pub(super) cut: Arrangement,
    /// Vertices without a source vertex.
    pub(super) crossings: usize,
    /// Candidate pairs examined so far.
    pub(super) pairs: usize,
}

/// One unique in-plane source edge.
struct SourceEdge {
    /// Exact points of the ends, the lower source vertex first.
    ends: [usize; 2],
    /// The triangles with this edge and their third vertices (exact points).
    carriers: Vec<(usize, usize)>,
}

impl InPlaneArrangement {
    /// Arrange `triangles`, the in-plane triangles of `mesh`.
    pub(super) fn build(
        mesh: &TriMesh,
        classifications: &[Classification],
        triangles: &[InPlaneTriangle],
        exact: &mut PlanePoints,
        max_pairs: usize,
        options: &ExecutionOptions,
    ) -> GeomResult<Self> {
        // The exact point of every source vertex, in vertex order.
        let mut sources: Vec<u32> = triangles
            .iter()
            .flat_map(|triangle| triangle.vertices)
            .collect();
        sources.sort_unstable();
        sources.dedup();
        let mut point_of = BTreeMap::<u32, usize>::new();
        let mut source_of = BTreeMap::<usize, u32>::new();
        for &vertex in &sources {
            let point = exact.node(mesh, classifications, NodeKey::Vertex(vertex))?;
            point_of.insert(vertex, point);
            source_of.insert(point, vertex);
        }
        let corners: Vec<[usize; 3]> = triangles
            .iter()
            .map(|triangle| triangle.vertices.map(|vertex| point_of[&vertex]))
            .collect();
        let mut by_key = BTreeMap::<EdgeKey, SourceEdge>::new();
        for (index, triangle) in triangles.iter().enumerate() {
            for (from, to, third) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
                let key = EdgeKey::new(triangle.vertices[from], triangle.vertices[to]);
                let ends = if triangle.vertices[from] < triangle.vertices[to] {
                    [corners[index][from], corners[index][to]]
                } else {
                    [corners[index][to], corners[index][from]]
                };
                by_key
                    .entry(key)
                    .or_insert_with(|| SourceEdge {
                        ends,
                        carriers: Vec::new(),
                    })
                    .carriers
                    .push((index, corners[index][third]));
            }
        }
        let keys: Vec<EdgeKey> = by_key.keys().copied().collect();
        let edges_in_order: Vec<SourceEdge> = by_key.into_values().collect();
        let segments: Vec<Segment> = edges_in_order
            .iter()
            .map(|edge| Segment {
                ends: edge.ends.map(PointId),
            })
            .collect();

        let mut budget = PairBudget::new(max_pairs, "section in-plane pairs");
        let declared = REQUEST_SCRATCH_PER_ELEMENT
            .saturating_mul(mesh.positions.len().max(mesh.triangle_count()));
        let cut = arrange(exact.store_mut(), &segments, &mut budget, declared, options)?;

        let vertices: Vec<Vertex> = cut
            .vertices
            .iter()
            .map(|vertex| {
                let mut sources: Vec<u32> = vertex
                    .inputs
                    .iter()
                    .filter_map(|point| source_of.get(&point.0).copied())
                    .collect();
                sources.sort_unstable();
                Vertex {
                    exact: vertex.point.0,
                    sources,
                    crossing: vertex.crossing.map(|crossing| {
                        (
                            (keys[crossing.segments[0]], keys[crossing.segments[1]]),
                            crossing.ends.map(|end| end.0),
                        )
                    }),
                }
            })
            .collect();
        let ends_of: Vec<(usize, usize)> = cut
            .pieces
            .iter()
            .map(|piece| (piece.ends[0].0, piece.ends[1].0))
            .collect();
        let carriers: Vec<BTreeSet<(usize, usize)>> = cut
            .pieces
            .iter()
            .map(|piece| {
                piece
                    .segments
                    .iter()
                    .flat_map(|&segment| edges_in_order[segment].carriers.iter().copied())
                    .collect()
            })
            .collect();
        let covered = coverage(
            (triangles, &corners),
            &Pieces {
                vertices: &vertices,
                ends_of: &ends_of,
                carriers: &carriers,
            },
            exact,
            &mut budget,
            options,
        )?;
        let pieces = covered
            .iter()
            .zip(&ends_of)
            .map(
                |(&(below_left, above_left, below_right, above_right), &(a, b))| Piece {
                    ends: [a, b],
                    left: Class::of(below_left, above_left),
                    right: Class::of(below_right, above_right),
                },
            )
            .collect();
        let mut edges = BTreeMap::<EdgeKey, Vec<(usize, bool)>>::new();
        for (key, path) in keys.iter().zip(&cut.paths) {
            edges.insert(
                *key,
                path.iter()
                    .map(|&(piece, forward)| (piece.0, forward))
                    .collect(),
            );
        }
        Ok(Self {
            vertices,
            pieces,
            edges,
            crossings: cut.crossings,
            cut,
            pairs: budget.used,
        })
    }

    /// What covers the left of the source edge `from` to `to`, per piece of
    /// it, in order along the edge.
    pub(super) fn left_of(&self, from: u32, to: u32) -> GeomResult<Vec<Class>> {
        let pieces = self
            .edges
            .get(&EdgeKey::new(from, to))
            .ok_or_else(internal_topology_error)?;
        let mut classes: Vec<Class> = pieces
            .iter()
            .map(|&(piece, forward)| {
                // `forward` is the key direction; the edge is walked from
                // `from`, which is the key direction when `from < to`.
                if forward == (from < to) {
                    self.pieces[piece].left
                } else {
                    self.pieces[piece].right
                }
            })
            .collect();
        if from > to {
            classes.reverse();
        }
        Ok(classes)
    }

    /// The returned plane-local position of a vertex: a source vertex or a
    /// crossing interpolated along its first edge, rounded once.
    pub(super) fn output_point(
        &self,
        mesh: &TriMesh,
        frame: Frame3,
        exact: &PlanePoints,
        vertex: usize,
    ) -> GeomResult<Point2> {
        let position = |index: u32| {
            mesh.positions
                .get(mesh_index(index)?)
                .copied()
                .ok_or_else(internal_topology_error)
        };
        let vertex = &self.vertices[vertex];
        match (vertex.sources.first(), vertex.crossing) {
            (Some(&source), _) => project_local(frame, position(source)?),
            (None, Some(((edge, _), ends))) => {
                let (a, b) = (position(edge.0)?, position(edge.1)?);
                let t = exact.segment_parameter(ends)?;
                project_local(frame, a + (b - a) * t)
            }
            (None, None) => Err(internal_topology_error()),
        }
    }
}

/// The pieces found so far: their vertices, ends and carrying triangles.
#[derive(Clone, Copy)]
struct Pieces<'a> {
    vertices: &'a [Vertex],
    ends_of: &'a [(usize, usize)],
    /// Per piece, the triangles over it with their third vertices.
    carriers: &'a [BTreeSet<(usize, usize)>],
}

/// For each piece, what covers its left and right: `(below left, above left,
/// below right, above right)`.
fn coverage(
    (triangles, corners): (&[InPlaneTriangle], &[[usize; 3]]),
    pieces: &Pieces<'_>,
    exact: &PlanePoints,
    budget: &mut PairBudget,
    options: &ExecutionOptions,
) -> GeomResult<Vec<(bool, bool, bool, bool)>> {
    let Pieces {
        vertices,
        ends_of,
        carriers,
    } = *pieces;
    let mut covered = vec![(false, false, false, false); ends_of.len()];
    let mut cover = |piece: usize, left: bool, right: bool, facing: Sign| {
        let slot = &mut covered[piece];
        match (facing == Sign::Positive, left, right) {
            (true, l, r) => {
                slot.0 |= l;
                slot.2 |= r;
            }
            (false, l, r) => {
                slot.1 |= l;
                slot.3 |= r;
            }
        }
    };
    // Carriers cover the side their third vertex is on.
    for (piece, set) in carriers.iter().enumerate() {
        let (from, to) = (
            vertices[ends_of[piece].0].exact,
            vertices[ends_of[piece].1].exact,
        );
        for &(triangle, third) in set {
            match exact.orient(from, to, third)? {
                Sign::Positive => cover(piece, true, false, triangles[triangle].facing),
                Sign::Negative => cover(piece, false, true, triangles[triangle].facing),
                _ => return Err(internal_topology_error()),
            }
        }
    }
    // Every other triangle covers both sides of a piece inside it.
    let piece_boxes: Vec<[f64; 4]> = ends_of
        .iter()
        .map(|&(a, b)| {
            union_box(&[
                exact.bounds(vertices[a].exact),
                exact.bounds(vertices[b].exact),
            ])
        })
        .collect();
    // `[min_u, min_v, max_u, max_v]` that the piece's extremes reach at least:
    // a piece inside a triangle lies in the triangle's box.
    let piece_cores: Vec<[f64; 4]> = ends_of
        .iter()
        .map(|&(a, b)| {
            let (a, b) = (
                exact.bounds(vertices[a].exact),
                exact.bounds(vertices[b].exact),
            );
            [
                a[2].min(b[2]),
                a[3].min(b[3]),
                a[0].max(b[0]),
                a[1].max(b[1]),
            ]
        })
        .collect();
    let triangle_boxes: Vec<[f64; 4]> = corners
        .iter()
        .map(|at| union_box(&at.map(|corner| exact.bounds(corner))))
        .collect();
    let by_min_u = |boxes: &[[f64; 4]]| {
        let mut order: Vec<usize> = (0..boxes.len()).collect();
        order.sort_by(|&left, &right| {
            boxes[left][0]
                .total_cmp(&boxes[right][0])
                .then(left.cmp(&right))
        });
        order
    };
    let (piece_order, triangle_order) = (by_min_u(&piece_boxes), by_min_u(&triangle_boxes));
    // The boxes in sweep order, to scan them without looking pieces and
    // triangles up.
    let in_order = |order: &[usize], boxes: &[[f64; 4]]| -> Vec<[f64; 4]> {
        order.iter().map(|&index| boxes[index]).collect()
    };
    let piece_swept = in_order(&piece_order, &piece_boxes);
    let triangle_swept = in_order(&triangle_order, &triangle_boxes);
    let mut consider = |piece: usize, triangle: usize| -> GeomResult<()> {
        options.check_cancelled()?;
        let (core, held) = (&piece_cores[piece], &triangle_boxes[triangle]);
        if core[0] < held[0]
            || core[1] < held[1]
            || core[2] > held[2]
            || core[3] > held[3]
            || carriers[piece]
                .iter()
                .any(|&(carrier, _)| carrier == triangle)
        {
            return Ok(());
        }
        budget.count()?;
        let (from, to) = (
            vertices[ends_of[piece].0].exact,
            vertices[ends_of[piece].1].exact,
        );
        let facing = triangles[triangle].facing;
        // Counter-clockwise in the frame's sense.
        let at = corners[triangle];
        let ring = if facing == Sign::Positive {
            [at[0], at[1], at[2]]
        } else {
            [at[0], at[2], at[1]]
        };
        if contains_piece(exact, &ring, from, to)? {
            cover(piece, true, true, facing);
        }
        Ok(())
    };
    // Triangles that start inside the span of a piece, then pieces that start
    // strictly inside the span of a triangle; the boxes must overlap in `v`
    // as well.
    for (&piece, held) in piece_order.iter().zip(&piece_swept) {
        options.check_cancelled()?;
        let start = triangle_swept.partition_point(|swept| swept[0] < held[0]);
        for (&triangle, swept) in triangle_order[start..].iter().zip(&triangle_swept[start..]) {
            if swept[0] > held[2] {
                break;
            }
            if swept[1] <= held[3] && held[1] <= swept[3] {
                consider(piece, triangle)?;
            }
        }
    }
    for (&triangle, held) in triangle_order.iter().zip(&triangle_swept) {
        options.check_cancelled()?;
        let start = piece_swept.partition_point(|swept| swept[0] <= held[0]);
        for (&piece, swept) in piece_order[start..].iter().zip(&piece_swept[start..]) {
            if swept[0] > held[2] {
                break;
            }
            if swept[1] <= held[3] && held[1] <= swept[3] {
                consider(piece, triangle)?;
            }
        }
    }
    Ok(covered)
}

/// Whether the piece `from`-`to`, which no edge crosses, lies in the interior
/// of the counter-clockwise triangle `ring`.
fn contains_piece(
    exact: &PlanePoints,
    ring: &[usize; 3],
    from: usize,
    to: usize,
) -> GeomResult<bool> {
    let mut on_line = false;
    for index in 0..3 {
        let (a, b) = (ring[index], ring[(index + 1) % 3]);
        let start = exact.orient(a, b, from)?;
        if start == Sign::Negative {
            return Ok(false);
        }
        let end = exact.orient(a, b, to)?;
        if end == Sign::Negative {
            return Ok(false);
        }
        on_line |= start == Sign::Zero && end == Sign::Zero;
    }
    if on_line {
        // On an edge of the triangle: that triangle carries the piece.
        return Err(internal_topology_error());
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use axiolid_contracts::GeomError;
    use axiolid_core::{Point3, Tolerance, Vec3};

    use super::super::regions::Projection;
    use super::super::{classify_vertices, ExactSectionPlane};
    use super::*;

    fn flat() -> Frame3 {
        Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }
    }

    fn mirrored() -> Frame3 {
        Frame3 {
            y: -Vec3::Y,
            z: -Vec3::Z,
            ..flat()
        }
    }

    fn options() -> ExecutionOptions {
        ExecutionOptions::new(Tolerance::METRE)
    }

    /// Arrange the triangles `faces`, given as corners (counter-clockwise in
    /// `xy`), each with its own three vertices unless the positions repeat
    /// and `share` is set.
    fn arrange_in(
        frame: Frame3,
        faces: &[[(f64, f64); 3]],
        share: bool,
        options: &ExecutionOptions,
        max_pairs: usize,
    ) -> GeomResult<(TriMesh, InPlaneArrangement, PlanePoints)> {
        let mut positions = Vec::<Point3>::new();
        let mut indices = Vec::new();
        for face in faces {
            for &(x, y) in face {
                let point = Point3::new(x, y, 0.0);
                let index = share
                    .then(|| positions.iter().position(|&known| known == point))
                    .flatten()
                    .unwrap_or_else(|| {
                        positions.push(point);
                        positions.len() - 1
                    });
                indices.push(index as u32);
            }
        }
        let mesh = TriMesh::new(positions, indices);
        let plane = ExactSectionPlane::new(frame).expect("plane");
        let projection = Projection::new(&plane).expect("projection");
        let classifications = classify_vertices(&mesh, plane).expect("classes");
        let triangles: Vec<InPlaneTriangle> = mesh
            .triangles()
            .enumerate()
            .map(|(index, vertices)| {
                InPlaneTriangle::new(&mesh, projection, index, vertices).expect("triangle")
            })
            .collect();
        let mut exact = PlanePoints::new(&plane, projection).expect("points");
        let arrangement = InPlaneArrangement::build(
            &mesh,
            &classifications,
            &triangles,
            &mut exact,
            max_pairs,
            options,
        )?;
        Ok((mesh, arrangement, exact))
    }

    fn arrange(faces: &[[(f64, f64); 3]]) -> (TriMesh, InPlaneArrangement) {
        let (mesh, arrangement, _) =
            arrange_in(flat(), faces, true, &options(), 1000).expect("arrangement");
        (mesh, arrangement)
    }

    fn square(x: f64, y: f64, size: f64) -> [[(f64, f64); 3]; 2] {
        let (a, b, c, d) = ((x, y), (x + size, y), (x + size, y + size), (x, y + size));
        [[a, b, c], [a, c, d]]
    }

    fn vertex_at(mesh: &TriMesh, x: f64, y: f64) -> u32 {
        mesh.positions
            .iter()
            .position(|point| point.x == x && point.y == y)
            .expect("vertex") as u32
    }

    #[test]
    fn squares_of_one_side_that_share_an_edge_dissolve_it() {
        let faces: Vec<_> = square(0.0, 0.0, 1.0)
            .into_iter()
            .chain(square(1.0, 0.0, 1.0))
            .collect();
        let (mesh, arrangement) = arrange(&faces);
        let (from, to) = (vertex_at(&mesh, 1.0, 0.0), vertex_at(&mesh, 1.0, 1.0));
        assert_eq!(arrangement.left_of(from, to).expect("left"), [Class::Below]);
        assert_eq!(arrangement.left_of(to, from).expect("left"), [Class::Below]);
        // The outside of the union is empty on the other side.
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 1.0, 0.0));
        assert_eq!(arrangement.left_of(a, b).expect("left"), [Class::Below]);
        assert_eq!(arrangement.left_of(b, a).expect("left"), [Class::Empty]);
        assert_eq!(arrangement.crossings, 0);
    }

    #[test]
    fn a_collinear_partial_overlap_splits_into_three_pieces() {
        // Bottom edges (0,0)-(2,0) and (1,0)-(3,0) overlap on [1, 2].
        let faces = [
            [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)],
            [(1.0, 0.0), (3.0, 0.0), (3.0, 1.0)],
        ];
        let (mesh, arrangement) = arrange(&faces);
        // The slanted edges (2,0)-(0,1) and (1,0)-(3,1) cross at (1.5, 0.25).
        assert_eq!(arrangement.crossings, 1);
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 2.0, 0.0));
        assert_eq!(
            arrangement.left_of(a, b).expect("left"),
            [Class::Below, Class::Below]
        );
        let (c, d) = (vertex_at(&mesh, 1.0, 0.0), vertex_at(&mesh, 3.0, 0.0));
        assert_eq!(
            arrangement.left_of(c, d).expect("left"),
            [Class::Below, Class::Below]
        );
        // Three stretches of the line y = 0, and the pieces of the others.
        let along = arrangement
            .pieces
            .iter()
            .filter(|piece| {
                piece.ends.iter().all(|&end| {
                    arrangement.vertices[end]
                        .sources
                        .first()
                        .is_some_and(|&source| mesh.positions[source as usize].y == 0.0)
                })
            })
            .count();
        assert_eq!(along, 3);
    }

    #[test]
    fn collinear_segments_of_one_vertex_split_the_longer() {
        // Both bottom edges leave (0,0): (0,0)-(2,0) lies on (0,0)-(3,0).
        let faces = [
            [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)],
            [(0.0, 0.0), (3.0, 0.0), (3.0, 1.0)],
        ];
        let (mesh, arrangement) = arrange(&faces);
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 3.0, 0.0));
        assert_eq!(
            arrangement.left_of(a, b).expect("left"),
            [Class::Below, Class::Below]
        );
    }

    #[test]
    fn a_t_junction_splits_the_edge_it_ends_on() {
        let faces = [
            [(0.0, 0.0), (2.0, 0.0), (1.0, 1.0)],
            [(1.0, 0.0), (0.5, -1.0), (1.5, -1.0)],
        ];
        let (mesh, arrangement) = arrange(&faces);
        // Four pieces for the first triangle's three edges (one is split) and
        // three for the second.
        assert_eq!(arrangement.pieces.len(), 7);
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 2.0, 0.0));
        assert_eq!(
            arrangement.left_of(a, b).expect("left"),
            [Class::Below, Class::Below]
        );
        assert_eq!(
            arrangement.left_of(b, a).expect("left"),
            [Class::Empty, Class::Empty]
        );
    }

    #[test]
    fn segments_through_one_point_share_a_vertex() {
        // The diagonal and the anti-diagonal of [0, 2]^2 and a vertical edge
        // all pass through (1, 1).
        let faces = [
            [(0.0, 0.0), (2.0, 2.0), (0.0, 2.0)],
            [(0.0, 2.0), (2.0, 0.0), (2.0, 2.0)],
            [(1.0, 0.0), (5.0, 1.0), (1.0, 2.0)],
        ];
        let (mesh, arrangement) = arrange(&faces);
        // (1,1), (2,.25), (2,1.75), (1.8,1.8) and (1.8,.2).
        assert_eq!(arrangement.crossings, 5);
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 2.0, 2.0));
        assert_eq!(arrangement.left_of(a, b).expect("left").len(), 3);
        let (c, d) = (vertex_at(&mesh, 1.0, 0.0), vertex_at(&mesh, 1.0, 2.0));
        // Cut once, at (1, 1), where the three edges meet: were their three
        // crossing points not one vertex, a zero-length piece would show.
        assert_eq!(arrangement.left_of(c, d).expect("left").len(), 2);
    }

    #[test]
    fn a_crossing_at_a_source_vertex_is_that_vertex() {
        // The diagonals cross at (1, 1), where a third triangle has a vertex.
        let faces = [
            [(0.0, 0.0), (2.0, 2.0), (0.0, 2.0)],
            [(0.0, 2.0), (2.0, 0.0), (2.0, 2.0)],
            [(1.0, 1.0), (1.5, 0.5), (1.5, 1.5)],
        ];
        let (mesh, arrangement) = arrange(&faces);
        assert_eq!(arrangement.crossings, 0);
        let centre = vertex_at(&mesh, 1.0, 1.0);
        let at_centre: Vec<_> = arrangement
            .vertices
            .iter()
            .filter(|vertex| vertex.sources == [centre])
            .collect();
        assert_eq!(at_centre.len(), 1);
        // The diagonal is cut at (1,1) and at (1.5,1.5) by the third triangle.
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 2.0, 2.0));
        assert_eq!(arrangement.left_of(a, b).expect("left").len(), 3);
    }

    #[test]
    fn a_piece_inside_a_triangle_is_covered_on_both_sides() {
        // The small triangle's edge from (0,1), on the big one's boundary, to
        // (1,1), inside it, has the big triangle on both sides.
        let faces = [
            [(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)],
            [(0.0, 1.0), (1.0, 1.0), (0.0, 2.0)],
        ];
        let (mesh, arrangement) = arrange(&faces);
        let (a, b) = (vertex_at(&mesh, 0.0, 1.0), vertex_at(&mesh, 1.0, 1.0));
        assert_eq!(arrangement.left_of(a, b).expect("left"), [Class::Below]);
        assert_eq!(arrangement.left_of(b, a).expect("left"), [Class::Below]);
    }

    #[test]
    fn a_piece_on_the_extension_of_an_edge_is_outside() {
        // The second triangle's bottom edge runs along y = 0 from the first
        // one's vertex (4, 0), outside it; the vertices are not shared.
        let faces = [
            [(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)],
            [(4.0, 0.0), (6.0, 0.0), (5.0, 1.0)],
        ];
        let (_, arrangement, _) =
            arrange_in(flat(), &faces, false, &options(), 1000).expect("arrangement");
        let (a, b) = (3, 4);
        assert_eq!(arrangement.left_of(a, b).expect("left"), [Class::Below]);
        assert_eq!(arrangement.left_of(b, a).expect("left"), [Class::Empty]);
        // The two vertices at (4, 0) are one point.
        assert!(arrangement
            .vertices
            .iter()
            .any(|vertex| vertex.sources.len() == 2));
    }

    #[test]
    fn opposite_facings_overlap_as_both() {
        let up = [(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)];
        let down = [(0.0, 0.0), (0.0, 2.0), (2.0, 0.0)];
        let (mesh, arrangement) = arrange(&[up, down]);
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 2.0, 0.0));
        assert_eq!(arrangement.left_of(a, b).expect("left"), [Class::Both]);
        assert_eq!(arrangement.left_of(b, a).expect("left"), [Class::Empty]);
    }

    #[test]
    fn a_mirrored_frame_mirrors_the_sides() {
        let face = [[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]];
        let (mesh, arrangement, _) =
            arrange_in(mirrored(), &face, true, &options(), 1000).expect("arrangement");
        let (a, b) = (vertex_at(&mesh, 0.0, 0.0), vertex_at(&mesh, 1.0, 0.0));
        // Clockwise in the frame: the material is above, on the frame's right
        // of the edge.
        assert_eq!(arrangement.left_of(a, b).expect("left"), [Class::Empty]);
        assert_eq!(arrangement.left_of(b, a).expect("left"), [Class::Above]);
    }

    #[test]
    fn a_unit_square_counts_its_segment_and_coverage_pairs() {
        let (_, arrangement) = arrange(&square(0.0, 0.0, 1.0));
        // Of five edges, ten pairs, two of them (opposite sides) apart; and
        // each side meets the triangle that does not carry it.
        assert_eq!(arrangement.pairs, 8 + 4);
        assert_eq!(arrangement.crossings, 0);
        let cramped = arrange_in(
            flat(),
            &square(0.0, 0.0, 1.0),
            true,
            &options(),
            arrangement.pairs - 1,
        );
        assert!(matches!(
            cramped,
            Err(GeomError::BudgetExceeded {
                resource: "section in-plane pairs"
            })
        ));
        let exact_fit = arrange_in(
            flat(),
            &square(0.0, 0.0, 1.0),
            true,
            &options(),
            arrangement.pairs,
        );
        assert!(exact_fit.is_ok());
    }

    #[test]
    fn the_output_point_of_a_crossing_is_interpolated_along_its_first_edge() {
        let faces = [
            [(0.0, 0.0), (0.5, 1.0), (0.0, 1.0)],
            [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        ];
        // (0,0)-(0.5,1) crosses (1,0)-(0,1) at (1/3, 2/3), which is not
        // exact in binary64; the other edges share the vertex (0, 1).
        let (mesh, arrangement, exact) =
            arrange_in(flat(), &faces, true, &options(), 1000).expect("arrangement");
        assert_eq!(arrangement.crossings, 1);
        let crossing = arrangement
            .vertices
            .iter()
            .position(|vertex| vertex.sources.is_empty())
            .expect("crossing");
        let at = arrangement
            .output_point(&mesh, flat(), &exact, crossing)
            .expect("point");
        assert!((at.x - 1.0 / 3.0).abs() < 1e-15 && (at.y - 2.0 / 3.0).abs() < 1e-15);
        let source = arrangement
            .output_point(&mesh, flat(), &exact, 0)
            .expect("point");
        assert_eq!(source, Point2::new(0.0, 0.0));
    }

    #[test]
    fn the_crossing_charge_refuses_beyond_the_memory_budget() {
        let faces = [
            [(0.0, 0.0), (0.5, 1.0), (0.0, 1.0)],
            [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        ];
        let fits = options().with_memory_budget(1 << 30);
        assert!(arrange_in(flat(), &faces, true, &fits, 1000).is_ok());
        // Three vertices: the declared scratch alone exceeds the budget.
        let tight = options().with_memory_budget(REQUEST_SCRATCH_PER_ELEMENT * 5);
        assert!(matches!(
            arrange_in(flat(), &faces, true, &tight, 1000),
            Err(GeomError::BudgetExceeded { resource: "memory" })
        ));
        // A request that is not charged for crossings builds without them.
        let none = [square(0.0, 0.0, 1.0)[0]];
        assert!(arrange_in(flat(), &none, true, &tight, 1000).is_ok());
    }

    #[test]
    fn a_cancelled_build_stops() {
        let token = axiolid_contracts::CancellationToken::new();
        token.cancel();
        let cancelled = options().with_cancellation(token);
        let result = arrange_in(flat(), &square(0.0, 0.0, 1.0), true, &cancelled, 1000);
        assert!(matches!(result, Err(GeomError::Cancelled)));
    }
}
