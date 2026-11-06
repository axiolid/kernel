//! The arrangement of straight segments on exact points.
//!
//! The segments may touch, cross, abut and overlap. They are cut at every
//! point where they meet into pieces that no other piece crosses, and every
//! decision is an exact sign on [`ExactPoints`]:
//!
//! 1. Equal points: the end points of the segments are joined where they are
//!    the same point, whatever their [`PointId`].
//! 2. Cuts: a sweep over the segments' boxes proposes pairs. A proper
//!    crossing creates a point ([`ExactPoints::segment_crossing`]); an end of
//!    one segment lying strictly inside the other, which covers T-junctions
//!    and collinear overlaps, splits it.
//! 3. Equal points again: along each segment its points are ordered exactly
//!    and equal neighbours are joined, so three segments through one point
//!    give one vertex.
//! 4. Pieces: the stretch between two consecutive distinct vertices of a
//!    segment, once however many segments run over it.
//!
//! Vertices are numbered in a deterministic order: those with input points
//! first, by their smallest [`PointId`], then the crossings, by their
//! smallest pair of segment indices. Pieces are numbered by their ends.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};

use crate::points::{inconsistent, ExactPoints, PointId};

/// Bytes charged for one crossing point: the Dyadic coordinates, the
/// enclosures and the bookkeeping of a segment pair.
pub const CROSSING_SCRATCH_BYTES: usize = 512;

/// Bytes charged for one end of a segment lying inside another one (a
/// T-junction split): an entry in two cut lists and its share of the sort.
pub const SPLIT_SCRATCH_BYTES: usize = 64;

/// A vertex of an [`Arrangement`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VertexId(pub usize);

/// A piece of an [`Arrangement`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PieceId(pub usize);

/// A count of candidate pairs against a limit, shared by the steps of one
/// request.
#[derive(Debug, Clone)]
pub struct PairBudget {
    /// Pairs counted so far.
    pub used: usize,
    /// Pairs allowed.
    pub max: usize,
    /// The resource named by `GeomError::BudgetExceeded` when exceeded.
    pub resource: &'static str,
}

impl PairBudget {
    /// A budget of `max` pairs, none used.
    pub const fn new(max: usize, resource: &'static str) -> Self {
        Self {
            used: 0,
            max,
            resource,
        }
    }

    /// Count one candidate pair.
    pub fn count(&mut self) -> GeomResult<()> {
        self.used += 1;
        if self.used > self.max {
            return Err(GeomError::BudgetExceeded {
                resource: self.resource,
            });
        }
        Ok(())
    }
}

/// A segment between two points of the store. The caller keeps its own data
/// per segment index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment {
    /// The ends.
    pub ends: [PointId; 2],
}

/// How the point of a vertex without input points came about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Crossing {
    /// The smallest pair of segments crossing here, ascending.
    pub segments: [usize; 2],
    /// The ends of those segments, in that order.
    pub ends: [PointId; 4],
}

/// A point of the arrangement.
#[derive(Debug, Clone)]
pub struct Vertex {
    /// An exact point at this position.
    pub point: PointId,
    /// The input points joined here, ascending; empty for a crossing.
    pub inputs: Vec<PointId>,
    /// For a vertex without input points, the crossing that created it.
    pub crossing: Option<Crossing>,
}

/// A stretch of segments no other segment crosses or touches inside.
#[derive(Debug, Clone)]
pub struct Piece {
    /// The vertices, `ends[0] < ends[1]`.
    pub ends: [VertexId; 2],
    /// The segments that run over the piece, ascending.
    pub segments: Vec<usize>,
}

/// The pieces of a set of segments.
#[derive(Debug, Clone)]
pub struct Arrangement {
    /// The vertices.
    pub vertices: Vec<Vertex>,
    /// The pieces, ordered by their ends.
    pub pieces: Vec<Piece>,
    /// Per input segment, its pieces in order from its first end to its
    /// second, each with whether that direction runs from `ends[0]` to
    /// `ends[1]` of the piece.
    pub paths: Vec<Vec<(PieceId, bool)>>,
    /// The number of vertices without input points.
    pub crossings: usize,
    /// For each vertex, the pieces at it and whether each leaves it at
    /// `ends[0]`.
    incident: Vec<Vec<(PieceId, bool)>>,
}

/// A point while the arrangement is built: its exact point and its label.
#[derive(Clone, Copy)]
enum Label {
    Input,
    Crossing([usize; 2], [PointId; 4]),
}

struct Cut {
    ends: [usize; 2],
    /// `[min_x, min_y, max_x, max_y]`.
    bounds: [f64; 4],
    /// Local points inside the segment, from the pairs examined.
    splits: Vec<usize>,
    /// All local points from the first end to the second, in order.
    path: Vec<usize>,
}

/// Points while building: exact point, label and union-find parent.
struct Local {
    exact: Vec<PointId>,
    labels: Vec<Label>,
    parent: Vec<usize>,
}

impl Local {
    fn push(&mut self, exact: PointId, label: Label) -> usize {
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

/// Whether two `[min_x, min_y, max_x, max_y]` boxes overlap, edges included.
pub fn boxes_overlap(first: &[f64; 4], second: &[f64; 4]) -> bool {
    first[0] <= second[2] && second[0] <= first[2] && first[1] <= second[3] && second[1] <= first[3]
}

/// The smallest box holding all `boxes`.
pub fn union_box(boxes: &[[f64; 4]]) -> [f64; 4] {
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

/// Sort `items` by an exact comparison that can fail.
pub(crate) fn sort_exact<T: Copy>(
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

/// Cut `segments` into pieces.
///
/// Every examined pair of segments counts against `budget`. Every crossing
/// and every T-junction split is charged: refused with `BudgetExceeded` for
/// `memory` when `declared_bytes`, plus `CROSSING_SCRATCH_BYTES` for each
/// crossing (the next one included), plus `SPLIT_SCRATCH_BYTES` for each
/// split, exceeds `options.memory_budget_bytes()`. `declared_bytes` is what the caller has declared per element of its own
/// input. Cancellation is checked per segment.
///
/// Input contract: the ends are points of `points`, a segment may repeat no
/// end and may not be empty, and segments may cross, touch, share an end or
/// overlap in any way. Two segments that share a `PointId` and lie on one
/// line are cut where one's other end lies strictly inside the other. Two
/// identical segments yield one piece with both as carriers.
pub fn arrange(
    points: &mut ExactPoints,
    segments: &[Segment],
    budget: &mut PairBudget,
    declared_bytes: usize,
    options: &ExecutionOptions,
) -> GeomResult<Arrangement> {
    // Local points: the distinct end points ascending, then the crossings.
    let mut inputs: Vec<PointId> = segments.iter().flat_map(|segment| segment.ends).collect();
    inputs.sort_unstable();
    inputs.dedup();
    let mut local = Local {
        exact: Vec::new(),
        labels: Vec::new(),
        parent: Vec::new(),
    };
    for &point in &inputs {
        local.push(point, Label::Input);
    }
    // Points at one position share a vertex.
    let mut by_position: Vec<usize> = (0..inputs.len()).collect();
    sort_exact(&mut by_position, |left, right| {
        points.lexicographic(local.exact[left], local.exact[right])
    })?;
    for pair in by_position.windows(2) {
        if points.same(local.exact[pair[0]], local.exact[pair[1]])? {
            local.union(pair[0], pair[1]);
        }
    }

    let local_of = |point: PointId| inputs.binary_search(&point).map_err(|_| inconsistent());
    let mut cuts = Vec::with_capacity(segments.len());
    for segment in segments {
        let ends = [local_of(segment.ends[0])?, local_of(segment.ends[1])?];
        let (low, high) = (
            points.bounds(segment.ends[0]),
            points.bounds(segment.ends[1]),
        );
        cuts.push(Cut {
            ends,
            bounds: [
                low[0].min(high[0]),
                low[1].min(high[1]),
                low[2].max(high[2]),
                low[3].max(high[3]),
            ],
            splits: Vec::new(),
            path: Vec::new(),
        });
    }

    cut_segments(
        segments,
        &mut cuts,
        points,
        &mut local,
        (budget, declared_bytes),
        options,
    )?;
    order_cuts(&mut cuts, points, &mut local, options)?;

    let (vertices, vertex_of) = number_vertices(&mut local);
    let crossings = vertices
        .iter()
        .filter(|vertex| vertex.inputs.is_empty())
        .count();

    // Pieces, once each however many segments run over them.
    let mut by_ends = BTreeMap::<(usize, usize), usize>::new();
    let mut carriers = Vec::<BTreeSet<usize>>::new();
    let mut on_segments = Vec::<(usize, usize, bool)>::new();
    for (index, cut) in cuts.iter().enumerate() {
        options.check_cancelled()?;
        let mut previous: Option<usize> = None;
        for &point in &cut.path {
            let vertex = vertex_of[&local.find(point)];
            if let Some(from) = previous {
                if from != vertex {
                    let ends = (from.min(vertex), from.max(vertex));
                    let next = carriers.len();
                    let piece = *by_ends.entry(ends).or_insert(next);
                    if piece == next {
                        carriers.push(BTreeSet::new());
                    }
                    carriers[piece].insert(index);
                    on_segments.push((index, piece, from < vertex));
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
    let pieces: Vec<Piece> = order
        .iter()
        .map(|&piece| Piece {
            ends: [VertexId(ends_of[piece].0), VertexId(ends_of[piece].1)],
            segments: carriers[piece].iter().copied().collect(),
        })
        .collect();
    let mut paths = vec![Vec::new(); segments.len()];
    for (segment, piece, forward) in on_segments {
        paths[segment].push((PieceId(renumber[piece]), forward));
    }
    let mut incident = vec![Vec::new(); vertices.len()];
    for (index, piece) in pieces.iter().enumerate() {
        incident[piece.ends[0].0].push((PieceId(index), true));
        incident[piece.ends[1].0].push((PieceId(index), false));
    }
    Ok(Arrangement {
        vertices,
        pieces,
        paths,
        crossings,
        incident,
    })
}

/// Examine the pairs of segments whose boxes overlap: join the proper
/// crossings and the ends that lie inside the other segment.
fn cut_segments(
    segments: &[Segment],
    cuts: &mut [Cut],
    points: &mut ExactPoints,
    local: &mut Local,
    (budget, declared_bytes): (&mut PairBudget, usize),
    options: &ExecutionOptions,
) -> GeomResult<()> {
    let mut order: Vec<usize> = (0..cuts.len()).collect();
    order.sort_by(|&left, &right| {
        cuts[left].bounds[0]
            .total_cmp(&cuts[right].bounds[0])
            .then(left.cmp(&right))
    });
    // The boxes in sweep order, to scan them without looking segments up.
    let boxes: Vec<[f64; 4]> = order.iter().map(|&index| cuts[index].bounds).collect();
    let mut charged = 0usize;
    for (rank, &first) in order.iter().enumerate() {
        options.check_cancelled()?;
        for (&second, held) in order[rank + 1..].iter().zip(&boxes[rank + 1..]) {
            if held[0] > boxes[rank][2] {
                break;
            }
            if !boxes_overlap(&boxes[rank], held) {
                continue;
            }
            options.check_cancelled()?;
            budget.count()?;
            let (low, high) = (first.min(second), first.max(second));
            let [a, b] = segments[low].ends;
            let [c, d] = segments[high].ends;
            // The same segment twice is one piece and cuts nothing.
            if [a, b] == [c, d] || [a, b] == [d, c] {
                continue;
            }
            // Segments of one point meet only there, unless they lie on one
            // line: then the nearer of the two other ends is inside the
            // longer segment and splits it.
            if let Some((shared, near, far)) = shared_end([a, b], [c, d]) {
                let (low_ends, high_ends) = (cuts[low].ends, cuts[high].ends);
                let local = |end: PointId| {
                    let (own, cut) = if [a, b].contains(&end) {
                        ([a, b], low_ends)
                    } else {
                        ([c, d], high_ends)
                    };
                    cut[usize::from(own[1] == end)]
                };
                if points.orient(shared, near, far)? == Sign::Zero {
                    // `near` and `far` belong to different segments.
                    if points.strictly_between(shared, near, far)? {
                        let longer = if [a, b].contains(&far) { low } else { high };
                        charge(options, declared_bytes, &mut charged, SPLIT_SCRATCH_BYTES)?;
                        cuts[longer].splits.push(local(near));
                    } else if points.strictly_between(shared, far, near)? {
                        let longer = if [a, b].contains(&near) { low } else { high };
                        charge(options, declared_bytes, &mut charged, SPLIT_SCRATCH_BYTES)?;
                        cuts[longer].splits.push(local(far));
                    }
                }
                continue;
            }
            // Both ends of one segment strictly on one side of the other's
            // line: they neither cross nor touch, and most pairs end here.
            let (past_c, past_d) = (points.orient(a, b, c)?, points.orient(a, b, d)?);
            if past_c == past_d && past_c != Sign::Zero {
                continue;
            }
            let (past_a, past_b) = (points.orient(c, d, a)?, points.orient(c, d, b)?);
            if past_a == past_b && past_a != Sign::Zero {
                continue;
            }
            let around = [past_c, past_d, past_a, past_b];
            let proper = around.iter().all(|&sign| sign != Sign::Zero)
                && around[0] != around[1]
                && around[2] != around[3];
            if proper {
                charge(
                    options,
                    declared_bytes,
                    &mut charged,
                    CROSSING_SCRATCH_BYTES,
                )?;
                let ends = [a, b, c, d];
                let point = points.segment_crossing(ends)?;
                let new = local.push(point, Label::Crossing([low, high], ends));
                cuts[low].splits.push(new);
                cuts[high].splits.push(new);
            } else {
                let (low_ends, high_ends) = (cuts[low].ends, cuts[high].ends);
                if around[0] == Sign::Zero && points.strictly_between(a, c, b)? {
                    charge(options, declared_bytes, &mut charged, SPLIT_SCRATCH_BYTES)?;
                    cuts[low].splits.push(high_ends[0]);
                }
                if around[1] == Sign::Zero && points.strictly_between(a, d, b)? {
                    charge(options, declared_bytes, &mut charged, SPLIT_SCRATCH_BYTES)?;
                    cuts[low].splits.push(high_ends[1]);
                }
                if around[2] == Sign::Zero && points.strictly_between(c, a, d)? {
                    charge(options, declared_bytes, &mut charged, SPLIT_SCRATCH_BYTES)?;
                    cuts[high].splits.push(low_ends[0]);
                }
                if around[3] == Sign::Zero && points.strictly_between(c, b, d)? {
                    charge(options, declared_bytes, &mut charged, SPLIT_SCRATCH_BYTES)?;
                    cuts[high].splits.push(low_ends[1]);
                }
            }
        }
    }
    Ok(())
}

/// The end two different segments share, with the other end of each (first
/// segment's first), or `None` when they share none.
fn shared_end(first: [PointId; 2], second: [PointId; 2]) -> Option<(PointId, PointId, PointId)> {
    for (i, &end) in first.iter().enumerate() {
        if let Some(j) = second.iter().position(|&other| other == end) {
            return Some((end, first[1 - i], second[1 - j]));
        }
    }
    None
}

/// Refuse the next charge when the scratch declared per element and the
/// crossings and splits charged so far no longer fit the memory budget.
pub(crate) fn charge(
    options: &ExecutionOptions,
    declared_bytes: usize,
    charged: &mut usize,
    cost: usize,
) -> GeomResult<()> {
    let Some(budget) = options.memory_budget_bytes() else {
        return Ok(());
    };
    let needed = charged
        .checked_add(cost)
        .and_then(|total| declared_bytes.checked_add(total));
    match needed {
        Some(needed) if needed <= budget => {
            *charged += cost;
            Ok(())
        }
        _ => Err(GeomError::BudgetExceeded { resource: "memory" }),
    }
}

/// Put the points of every segment in order along it, and join equal ones.
fn order_cuts(
    cuts: &mut [Cut],
    points: &ExactPoints,
    local: &mut Local,
    options: &ExecutionOptions,
) -> GeomResult<()> {
    for cut in cuts {
        options.check_cancelled()?;
        let [first, last] = cut.ends.map(|end| local.exact[end]);
        let axis = points.dominant(first, last)?;
        let descending = points.compare(first, last, axis)? == std::cmp::Ordering::Greater;
        let mut inner = std::mem::take(&mut cut.splits);
        inner.sort_unstable();
        inner.dedup();
        sort_exact(&mut inner, |left, right| {
            let order = points.compare(local.exact[left], local.exact[right], axis)?;
            Ok(if descending { order.reverse() } else { order })
        })?;
        cut.path = std::iter::once(cut.ends[0])
            .chain(inner)
            .chain(std::iter::once(cut.ends[1]))
            .collect();
        for pair in cut.path.windows(2) {
            if points.same(local.exact[pair[0]], local.exact[pair[1]])? {
                local.union(pair[0], pair[1]);
            }
        }
    }
    Ok(())
}

/// Number the groups of equal points: those with input points first, by their
/// smallest one, then the crossings, by their smallest pair of segments.
fn number_vertices(local: &mut Local) -> (Vec<Vertex>, BTreeMap<usize, usize>) {
    struct Group {
        inputs: Vec<PointId>,
        /// The smallest crossing label of the group, with its point.
        crossing: Option<(Crossing, usize)>,
        /// Any point of the group.
        any: usize,
    }
    let mut groups = BTreeMap::<usize, Group>::new();
    for point in 0..local.parent.len() {
        let root = local.find(point);
        let group = groups.entry(root).or_insert(Group {
            inputs: Vec::new(),
            crossing: None,
            any: point,
        });
        match local.labels[point] {
            Label::Input => group.inputs.push(local.exact[point]),
            Label::Crossing(segments, ends) => {
                if group
                    .crossing
                    .is_none_or(|(known, _)| segments < known.segments)
                {
                    group.crossing = Some((Crossing { segments, ends }, point));
                }
            }
        }
    }
    let mut ordered: Vec<(usize, Group)> = groups.into_iter().collect();
    for (_, group) in &mut ordered {
        group.inputs.sort_unstable();
    }
    ordered.sort_by(
        |(_, left), (_, right)| match (left.inputs.first(), right.inputs.first()) {
            (Some(a), Some(b)) => a.cmp(b),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => left
                .crossing
                .map(|c| c.0.segments)
                .cmp(&right.crossing.map(|c| c.0.segments)),
        },
    );
    let mut vertex_of = BTreeMap::new();
    let mut vertices = Vec::with_capacity(ordered.len());
    for (number, (root, group)) in ordered.into_iter().enumerate() {
        vertex_of.insert(root, number);
        let (crossing, point) = match (group.inputs.is_empty(), group.crossing) {
            (true, Some((crossing, point))) => (Some(crossing), point),
            _ => (None, group.any),
        };
        vertices.push(Vertex {
            point: local.exact[point],
            inputs: group.inputs,
            crossing,
        });
    }
    (vertices, vertex_of)
}

impl Arrangement {
    /// The pieces at `vertex` in counter-clockwise order in the points'
    /// frame, each with whether it leaves the vertex at `ends[0]`.
    ///
    /// Both `orient` and the half-plane split are in the stored axes up to
    /// the frame sign, which reverses the order within a half but keeps the
    /// cyclic order. Two pieces that leave a vertex the same way overlap,
    /// which `arrange` has cut apart, so they are an inconsistency.
    pub fn fan(&self, points: &ExactPoints, vertex: VertexId) -> GeomResult<Vec<(PieceId, bool)>> {
        let here = self.vertices[vertex.0].point;
        let mut spokes = Vec::with_capacity(self.incident[vertex.0].len());
        for &(piece, forward) in &self.incident[vertex.0] {
            let ends = self.pieces[piece.0].ends;
            let to = self.vertices[ends[usize::from(forward)].0].point;
            let upper = points.compare(to, here, 1)? == std::cmp::Ordering::Greater
                || (points.compare(to, here, 1)? == std::cmp::Ordering::Equal
                    && points.compare(to, here, 0)? == std::cmp::Ordering::Greater);
            spokes.push((piece, forward, to, u8::from(!upper)));
        }
        sort_exact(&mut spokes, |left, right| {
            if left.3 != right.3 {
                return Ok(left.3.cmp(&right.3));
            }
            match points.orient(here, left.2, right.2)? {
                Sign::Positive => Ok(std::cmp::Ordering::Less),
                Sign::Negative => Ok(std::cmp::Ordering::Greater),
                _ => Err(inconsistent()),
            }
        })?;
        Ok(spokes
            .into_iter()
            .map(|(piece, forward, _, _)| (piece, forward))
            .collect())
    }
}
