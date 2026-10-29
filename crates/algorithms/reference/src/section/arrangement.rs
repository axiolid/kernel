//! The exact arrangement of the in-plane faces of a section (ADR 0093, #366).
//!
//! The in-plane triangles may touch, abut and overlap. Their source edges are
//! cut at every point where they meet, into pieces that no other piece
//! crosses, and each piece records which class of material covers each of its
//! sides. Every decision is an exact sign on `PlanePoints`; the returned
//! coordinates only label the output.
//!
//! 1. Vertices: each in-plane source vertex is an exact point; vertices at
//!    one position share a point.
//! 2. Segments: the unique in-plane source edges, a sweep over their boxes
//!    proposes pairs. A proper crossing creates a crossing point
//!    (`PlanePoints::segment_crossing`); an end of one segment lying strictly
//!    inside the other, which covers T-junctions and collinear overlaps,
//!    splits it.
//! 3. Equal points: along each segment its points are ordered exactly and
//!    equal neighbours are joined, so three segments through one point, or a
//!    crossing at a source vertex, give one vertex.
//! 4. Pieces: the stretch between two consecutive distinct vertices of a
//!    segment, once however many source edges run over it.
//! 5. Coverage: a triangle that carries a piece covers the side its third
//!    vertex is on; every other triangle either contains the piece in its
//!    interior, covering both sides, or does not. A triangle is convex, so a
//!    piece lies inside it or outside, never partly, and the pieces need no
//!    face tracing.
//!
//! A side is `Below` when only faces whose material lies below the plane
//! cover it, `Above` when only the others do, `Both` when both kinds do and
//! `Empty` when none does.
//!
//! Candidate pairs of segments and of pieces with triangles count against
//! `RegionLimits::max_in_plane_pairs`. The crossing points are bounded by the
//! pairs and charged against the memory budget as they are created.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};
use axiolid_core::{Frame3, Point2};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::InPlaneSide;

use super::plane_points::PlanePoints;
use super::regions::InPlaneTriangle;
use super::{
    internal_topology_error, mesh_index, project_local, Classification, EdgeKey, NodeKey,
    REQUEST_SCRATCH_PER_ELEMENT,
};

/// Bytes charged for one crossing point: the Dyadic coordinates, the
/// enclosure, the label and the union-find entry (ADR 0088).
const CROSSING_SCRATCH_BYTES: usize = 512;

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
    /// Vertices without a source vertex.
    pub(super) crossings: usize,
    /// Candidate pairs examined so far.
    pub(super) pairs: usize,
}

/// A point while the arrangement is built: its exact point and its label.
enum Label {
    Source(u32),
    Crossing((EdgeKey, EdgeKey), [usize; 4]),
}

/// One unique in-plane source edge.
struct Segment {
    key: EdgeKey,
    /// Local point numbers of the ends.
    ends: [usize; 2],
    /// The triangles with this edge and their third vertices (local points).
    carriers: Vec<(usize, usize)>,
    /// `[min_u, min_v, max_u, max_v]`.
    bounds: [f64; 4],
    /// Local points inside the segment, from the pairs examined.
    splits: Vec<usize>,
    /// All local points from the first end to the second, in order.
    path: Vec<usize>,
}

/// Points while building: exact point, label and union-find parent.
struct Points {
    exact: Vec<usize>,
    labels: Vec<Label>,
    parent: Vec<usize>,
}

impl Points {
    fn push(&mut self, exact: usize, label: Label) -> usize {
        self.exact.push(exact);
        self.labels.push(label);
        self.parent.push(self.parent.len());
        self.parent.len() - 1
    }

    fn find(&mut self, mut point: usize) -> usize {
        while self.parent[point] != point {
            self.parent[point] = self.parent[self.parent[point]];
            point = self.parent[point];
        }
        point
    }

    /// Join two points; the root is the smaller number.
    fn union(&mut self, left: usize, right: usize) {
        let (left, right) = (self.find(left), self.find(right));
        self.parent[left.max(right)] = left.min(right);
    }
}

fn overlap(first: &[f64; 4], second: &[f64; 4]) -> bool {
    first[0] <= second[2] && second[0] <= first[2] && first[1] <= second[3] && second[1] <= first[3]
}

fn union_box(boxes: &[[f64; 4]]) -> [f64; 4] {
    boxes.iter().fold(
        [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ],
        |all, one| {
            [
                all[0].min(one[0]),
                all[1].min(one[1]),
                all[2].max(one[2]),
                all[3].max(one[3]),
            ]
        },
    )
}

fn count_pair(pairs: &mut usize, max_pairs: usize) -> GeomResult<()> {
    *pairs += 1;
    if *pairs > max_pairs {
        return Err(GeomError::BudgetExceeded {
            resource: "section in-plane pairs",
        });
    }
    Ok(())
}

/// Sort `items` by an exact comparison that can fail.
fn sort_exact<T: Copy>(
    items: &mut [T],
    mut compare: impl FnMut(T, T) -> GeomResult<std::cmp::Ordering>,
) -> GeomResult<()> {
    let mut failure = None;
    items.sort_by(|&left, &right| match compare(left, right) {
        Ok(order) => order,
        Err(error) => {
            failure.get_or_insert(error);
            std::cmp::Ordering::Equal
        }
    });
    failure.map_or(Ok(()), Err)
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
        let mut points = Points {
            exact: Vec::new(),
            labels: Vec::new(),
            parent: Vec::new(),
        };
        let corners = source_points(mesh, classifications, triangles, exact, &mut points)?;
        let mut segments = segments_of(triangles, &corners, exact, &points);
        let mut pairs = 0usize;
        cut_segments(
            mesh,
            &mut segments,
            exact,
            &mut points,
            (&mut pairs, max_pairs),
            options,
        )?;
        order_segments(&mut segments, exact, &mut points, options)?;

        let (vertices, vertex_of) = number_vertices(&mut points);
        let crossing_vertices = vertices
            .iter()
            .filter(|vertex| vertex.sources.is_empty())
            .count();

        // Pieces, once each however many segments run over them.
        let mut by_ends = BTreeMap::<(usize, usize), usize>::new();
        let mut carriers = Vec::<BTreeSet<(usize, usize)>>::new();
        let mut on_edges = Vec::<(EdgeKey, usize, bool)>::new();
        for segment in &segments {
            options.check_cancelled()?;
            let mut previous: Option<usize> = None;
            for &point in &segment.path {
                let vertex = vertex_of[&points.find(point)];
                if let Some(from) = previous {
                    if from != vertex {
                        let ends = (from.min(vertex), from.max(vertex));
                        let next = carriers.len();
                        let piece = *by_ends.entry(ends).or_insert(next);
                        if piece == next {
                            carriers.push(BTreeSet::new());
                        }
                        carriers[piece].extend(segment.carriers.iter().copied());
                        on_edges.push((segment.key, piece, from < vertex));
                    }
                }
                previous = Some(vertex);
            }
        }
        // Piece numbers in order of their ends.
        let mut ends_of = vec![(0, 0); carriers.len()];
        for (&ends, &piece) in &by_ends {
            ends_of[piece] = ends;
        }
        let mut renumber = vec![0usize; carriers.len()];
        let mut order: Vec<usize> = (0..carriers.len()).collect();
        order.sort_by_key(|&piece| ends_of[piece]);
        for (rank, &piece) in order.iter().enumerate() {
            renumber[piece] = rank;
        }

        let covered = coverage(
            (triangles, &corners),
            &points,
            &Pieces {
                vertices: &vertices,
                ends_of: &ends_of,
                carriers: &carriers,
            },
            exact,
            (&mut pairs, max_pairs),
            options,
        )?;
        let mut pieces = Vec::with_capacity(order.len());
        for &piece in &order {
            let (below_left, above_left, below_right, above_right) = covered[piece];
            pieces.push(Piece {
                ends: [ends_of[piece].0, ends_of[piece].1],
                left: Class::of(below_left, above_left),
                right: Class::of(below_right, above_right),
            });
        }
        let mut edges = BTreeMap::<EdgeKey, Vec<(usize, bool)>>::new();
        for (key, piece, forward) in on_edges {
            edges
                .entry(key)
                .or_default()
                .push((renumber[piece], forward));
        }
        Ok(Self {
            vertices,
            pieces,
            edges,
            crossings: crossing_vertices,
            pairs,
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

/// Intern the source vertices of the triangles, and join those at one
/// position. Returns the local points of each triangle's corners.
fn source_points(
    mesh: &TriMesh,
    classifications: &[Classification],
    triangles: &[InPlaneTriangle],
    exact: &mut PlanePoints,
    points: &mut Points,
) -> GeomResult<Vec<[usize; 3]>> {
    let mut local = BTreeMap::<u32, usize>::new();
    let mut corners = Vec::with_capacity(triangles.len());
    for triangle in triangles {
        let mut at = [0usize; 3];
        for (slot, &vertex) in at.iter_mut().zip(&triangle.vertices) {
            *slot = match local.get(&vertex) {
                Some(&known) => known,
                None => {
                    let point = exact.node(mesh, classifications, NodeKey::Vertex(vertex))?;
                    let new = points.push(point, Label::Source(vertex));
                    local.insert(vertex, new);
                    new
                }
            };
        }
        corners.push(at);
    }
    // Vertices at one position share a point.
    let mut sources: Vec<usize> = local.values().copied().collect();
    sort_exact(&mut sources, |left, right| {
        exact.lexicographic(points.exact[left], points.exact[right])
    })?;
    for pair in sources.windows(2) {
        if exact.same(points.exact[pair[0]], points.exact[pair[1]])? {
            points.union(pair[0], pair[1]);
        }
    }
    Ok(corners)
}

/// The unique in-plane source edges with their carriers and boxes.
fn segments_of(
    triangles: &[InPlaneTriangle],
    corners: &[[usize; 3]],
    exact: &PlanePoints,
    points: &Points,
) -> Vec<Segment> {
    let mut by_key = BTreeMap::<EdgeKey, Segment>::new();
    for (index, triangle) in triangles.iter().enumerate() {
        for (from, to, third) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
            let key = EdgeKey::new(triangle.vertices[from], triangle.vertices[to]);
            let (first, second) = if triangle.vertices[from] < triangle.vertices[to] {
                (corners[index][from], corners[index][to])
            } else {
                (corners[index][to], corners[index][from])
            };
            let segment = by_key.entry(key).or_insert_with(|| {
                let (low, high) = (
                    exact.bounds(points.exact[first]),
                    exact.bounds(points.exact[second]),
                );
                Segment {
                    key,
                    ends: [first, second],
                    carriers: Vec::new(),
                    bounds: [
                        low[0].min(high[0]),
                        low[1].min(high[1]),
                        low[2].max(high[2]),
                        low[3].max(high[3]),
                    ],
                    splits: Vec::new(),
                    path: Vec::new(),
                }
            });
            segment.carriers.push((index, corners[index][third]));
        }
    }
    by_key.into_values().collect()
}

/// Examine the pairs of segments whose boxes overlap: join the proper
/// crossings and the ends that lie inside the other segment.
fn cut_segments(
    mesh: &TriMesh,
    segments: &mut [Segment],
    exact: &mut PlanePoints,
    points: &mut Points,
    (pairs, max_pairs): (&mut usize, usize),
    options: &ExecutionOptions,
) -> GeomResult<()> {
    let mut order: Vec<usize> = (0..segments.len()).collect();
    order.sort_by(|&left, &right| {
        segments[left].bounds[0]
            .total_cmp(&segments[right].bounds[0])
            .then(left.cmp(&right))
    });
    let elements = mesh.positions.len().max(mesh.triangle_count());
    let mut created = 0usize;
    for (rank, &first) in order.iter().enumerate() {
        for &second in &order[rank + 1..] {
            if segments[second].bounds[0] > segments[first].bounds[2] {
                break;
            }
            options.check_cancelled()?;
            if !overlap(&segments[first].bounds, &segments[second].bounds) {
                continue;
            }
            count_pair(pairs, max_pairs)?;
            let (low, high) = (first.min(second), first.max(second));
            let [a, b] = segments[low].ends.map(|end| points.exact[end]);
            let [c, d] = segments[high].ends.map(|end| points.exact[end]);
            let around = [
                exact.orient(a, b, c)?,
                exact.orient(a, b, d)?,
                exact.orient(c, d, a)?,
                exact.orient(c, d, b)?,
            ];
            let proper = around.iter().all(|&sign| sign != Sign::Zero)
                && around[0] != around[1]
                && around[2] != around[3];
            if proper {
                charge_crossing(options, elements, created)?;
                let keys = (segments[low].key, segments[high].key);
                let ends = [a, b, c, d];
                let point = exact.segment_crossing(keys, ends)?;
                let new = points.push(point, Label::Crossing(keys, ends));
                created += 1;
                segments[low].splits.push(new);
                segments[high].splits.push(new);
            } else {
                let (low_ends, high_ends) = (segments[low].ends, segments[high].ends);
                if around[0] == Sign::Zero && exact.strictly_between(a, c, b)? {
                    segments[low].splits.push(high_ends[0]);
                }
                if around[1] == Sign::Zero && exact.strictly_between(a, d, b)? {
                    segments[low].splits.push(high_ends[1]);
                }
                if around[2] == Sign::Zero && exact.strictly_between(c, a, d)? {
                    segments[high].splits.push(low_ends[0]);
                }
                if around[3] == Sign::Zero && exact.strictly_between(c, b, d)? {
                    segments[high].splits.push(low_ends[1]);
                }
            }
        }
    }
    Ok(())
}

/// Refuse the next crossing when the scratch declared per element and the
/// crossings made so far no longer fit the memory budget.
fn charge_crossing(options: &ExecutionOptions, elements: usize, created: usize) -> GeomResult<()> {
    let Some(budget) = options.memory_budget_bytes() else {
        return Ok(());
    };
    let needed = REQUEST_SCRATCH_PER_ELEMENT
        .checked_mul(elements)
        .and_then(|declared| {
            CROSSING_SCRATCH_BYTES
                .checked_mul(created + 1)
                .and_then(|crossings| declared.checked_add(crossings))
        });
    match needed {
        Some(needed) if needed <= budget => Ok(()),
        _ => Err(GeomError::BudgetExceeded { resource: "memory" }),
    }
}

/// Put the points of every segment in order along it, and join equal ones.
fn order_segments(
    segments: &mut [Segment],
    exact: &PlanePoints,
    points: &mut Points,
    options: &ExecutionOptions,
) -> GeomResult<()> {
    for segment in segments {
        options.check_cancelled()?;
        let [first, last] = segment.ends.map(|end| points.exact[end]);
        let axis = exact.dominant(first, last)?;
        let descending = exact.compare(first, last, axis)? == std::cmp::Ordering::Greater;
        let mut inner = std::mem::take(&mut segment.splits);
        inner.sort_unstable();
        inner.dedup();
        sort_exact(&mut inner, |left, right| {
            let order = exact.compare(points.exact[left], points.exact[right], axis)?;
            Ok(if descending { order.reverse() } else { order })
        })?;
        segment.path = std::iter::once(segment.ends[0])
            .chain(inner)
            .chain(std::iter::once(segment.ends[1]))
            .collect();
        for pair in segment.path.windows(2) {
            if exact.same(points.exact[pair[0]], points.exact[pair[1]])? {
                points.union(pair[0], pair[1]);
            }
        }
    }
    Ok(())
}

/// Number the groups of equal points: those with source vertices first, by
/// their smallest vertex, then the crossings, by their smallest pair of edges.
fn number_vertices(points: &mut Points) -> (Vec<Vertex>, BTreeMap<usize, usize>) {
    struct Group {
        sources: Vec<u32>,
        /// The smallest crossing label of the group, with its point.
        crossing: Option<((EdgeKey, EdgeKey), usize)>,
        /// Any point of the group.
        any: usize,
    }
    let mut groups = BTreeMap::<usize, Group>::new();
    for point in 0..points.parent.len() {
        let root = points.find(point);
        let group = groups.entry(root).or_insert(Group {
            sources: Vec::new(),
            crossing: None,
            any: point,
        });
        match &points.labels[point] {
            Label::Source(vertex) => group.sources.push(*vertex),
            Label::Crossing(keys, _) => {
                if group.crossing.is_none_or(|(known, _)| *keys < known) {
                    group.crossing = Some((*keys, point));
                }
            }
        }
    }
    let mut ordered: Vec<(usize, Group)> = groups.into_iter().collect();
    for (_, group) in &mut ordered {
        group.sources.sort_unstable();
    }
    ordered.sort_by(
        |(_, left), (_, right)| match (left.sources.first(), right.sources.first()) {
            (Some(a), Some(b)) => a.cmp(b),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => left.crossing.map(|c| c.0).cmp(&right.crossing.map(|c| c.0)),
        },
    );
    let mut vertex_of = BTreeMap::new();
    let mut vertices = Vec::with_capacity(ordered.len());
    for (number, (root, group)) in ordered.into_iter().enumerate() {
        vertex_of.insert(root, number);
        let crossing = match (group.sources.is_empty(), group.crossing) {
            (true, Some((keys, point))) => match &points.labels[point] {
                Label::Crossing(_, ends) => Some((keys, *ends)),
                Label::Source(_) => None,
            },
            _ => None,
        };
        let point = match (group.sources.is_empty(), group.crossing) {
            (true, Some((_, point))) => point,
            _ => group.any,
        };
        vertices.push(Vertex {
            exact: points.exact[point],
            sources: group.sources,
            crossing,
        });
    }
    (vertices, vertex_of)
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
    points: &Points,
    pieces: &Pieces<'_>,
    exact: &PlanePoints,
    (pairs, max_pairs): (&mut usize, usize),
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
            match exact.orient(from, to, points.exact[third])? {
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
    let triangle_boxes: Vec<[f64; 4]> = corners
        .iter()
        .map(|at| union_box(&at.map(|corner| exact.bounds(points.exact[corner]))))
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
    let mut candidates = Vec::new();
    // Triangles that start inside the span of a piece, then pieces that start
    // strictly inside the span of a triangle.
    for &piece in &piece_order {
        options.check_cancelled()?;
        let start =
            triangle_order.partition_point(|&t| triangle_boxes[t][0] < piece_boxes[piece][0]);
        for &triangle in &triangle_order[start..] {
            if triangle_boxes[triangle][0] > piece_boxes[piece][2] {
                break;
            }
            candidates.push((piece, triangle));
        }
    }
    for &triangle in &triangle_order {
        options.check_cancelled()?;
        let start =
            piece_order.partition_point(|&p| piece_boxes[p][0] <= triangle_boxes[triangle][0]);
        for &piece in &piece_order[start..] {
            if piece_boxes[piece][0] > triangle_boxes[triangle][2] {
                break;
            }
            candidates.push((piece, triangle));
        }
    }
    for (piece, triangle) in candidates {
        options.check_cancelled()?;
        if piece_boxes[piece][1] > triangle_boxes[triangle][3]
            || triangle_boxes[triangle][1] > piece_boxes[piece][3]
            || carriers[piece]
                .iter()
                .any(|&(carrier, _)| carrier == triangle)
        {
            continue;
        }
        count_pair(pairs, max_pairs)?;
        let (from, to) = (
            vertices[ends_of[piece].0].exact,
            vertices[ends_of[piece].1].exact,
        );
        let facing = triangles[triangle].facing;
        // Counter-clockwise in the frame's sense.
        let at = corners[triangle].map(|corner| points.exact[corner]);
        let ring = if facing == Sign::Positive {
            [at[0], at[1], at[2]]
        } else {
            [at[0], at[2], at[1]]
        };
        if contains_piece(exact, &ring, from, to)? {
            cover(piece, true, true, facing);
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
        let (start, end) = (exact.orient(a, b, from)?, exact.orient(a, b, to)?);
        if start == Sign::Negative || end == Sign::Negative {
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
