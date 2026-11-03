//! Rings of labelled half-edges, outer boundaries and their holes.
//!
//! The caller says which pieces bound a class (a region of one label) and on
//! which side the class lies: a half-edge is a piece with a direction and the
//! class on its left. [`trace`] closes them into rings and nests them:
//!
//! 1. Cycles: at a vertex the half-edge that follows an incoming one is the
//!    first boundary of the class met turning clockwise from the way it came
//!    (the pieces at a vertex are ordered by [`Arrangement::fan`]), so a
//!    vertex where a region pinches is resolved per fan.
//! 2. Each cycle is an outer ring or a hole by the exact turn at its extreme
//!    vertex.
//! 3. Vertices that the caller does not need are dropped when they lie
//!    strictly between their neighbours on a line.
//! 4. Each hole belongs to the innermost outer ring that holds it, decided
//!    with exact point location.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};

use crate::arrangement::{Arrangement, PairBudget, PieceId, VertexId};
use crate::points::{inconsistent, ExactPoints, Location, PointId, RingBox};

/// A piece with a direction: from `ends[0]` to `ends[1]` when `forward`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HalfEdge {
    /// The piece.
    pub piece: PieceId,
    /// Whether it runs from `ends[0]` to `ends[1]`.
    pub forward: bool,
}

/// A closed boundary of a class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TracedRing {
    /// The half-edges of the cycle, before vertices were dropped.
    pub half_edges: Vec<HalfEdge>,
    /// The vertices kept, in the order of the cycle.
    pub vertices: Vec<VertexId>,
    /// Whether the ring is a hole (clockwise in the frame) rather than an
    /// outer boundary.
    pub hole: bool,
}

/// An outer ring and its holes, as indices into the rings of [`trace`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TracedRegion {
    /// The outer ring.
    pub outer: usize,
    /// The holes, ascending.
    pub holes: Vec<usize>,
}

/// The eight lexicographic orders tried for the extreme vertex of a ring:
/// which axis is major, and whether each runs ascending.
const EXTREME_ORDERS: [(usize, bool, bool); 8] = [
    (0, true, true),
    (0, true, false),
    (0, false, true),
    (0, false, false),
    (1, true, true),
    (1, true, false),
    (1, false, true),
    (1, false, false),
];

/// Whether the cycle is an outer boundary (counter-clockwise in the frame),
/// by the exact turn at its extreme vertex.
///
/// The minimum of a lexicographic order is a convex corner. The first order
/// whose minimum the ring visits once decides; a ring that touches itself at
/// every extreme (a weakly simple ring) is refused.
fn is_outer(arrangement: &Arrangement, exact: &ExactPoints, ring: &[usize]) -> GeomResult<bool> {
    let point = |index: usize| arrangement.vertices[ring[index % ring.len()]].point;
    for (major, major_up, minor_up) in EXTREME_ORDERS {
        let order = |left: usize, right: usize| -> GeomResult<std::cmp::Ordering> {
            let (a, b) = (point(left), point(right));
            let first = exact.compare(a, b, major)?;
            let first = if major_up { first } else { first.reverse() };
            if first.is_ne() {
                return Ok(first);
            }
            let second = exact.compare(a, b, 1 - major)?;
            Ok(if minor_up { second } else { second.reverse() })
        };
        let mut smallest = 0;
        for index in 1..ring.len() {
            if order(index, smallest)?.is_lt() {
                smallest = index;
            }
        }
        let mut visits = 0;
        for index in 0..ring.len() {
            if order(index, smallest)?.is_eq() {
                visits += 1;
            }
        }
        if visits != 1 {
            continue;
        }
        let count = ring.len();
        return match exact.orient(
            point(smallest + count - 1),
            point(smallest),
            point(smallest + 1),
        )? {
            Sign::Positive => Ok(true),
            Sign::Negative => Ok(false),
            _ => Err(GeomError::Degenerate(
                "a ring has a spike at its extreme vertex".into(),
            )),
        };
    }
    Err(GeomError::Degenerate(
        "a ring visits every extreme vertex twice".into(),
    ))
}

/// Drop vertices that lie strictly between their neighbours on a line and
/// that `keep` does not ask for; a spike is refused.
fn simplify(
    arrangement: &Arrangement,
    exact: &ExactPoints,
    ring: &mut Vec<usize>,
    keep: &dyn Fn(VertexId) -> bool,
    options: &ExecutionOptions,
) -> GeomResult<()> {
    let collapsed = || GeomError::Degenerate("a ring collapsed under simplification".into());
    let mut index = 0usize;
    let mut clean = 0usize;
    while clean < ring.len() {
        options.check_cancelled()?;
        if ring.len() < 3 {
            return Err(collapsed());
        }
        let at = index % ring.len();
        let (previous, vertex, next) = (
            ring[(at + ring.len() - 1) % ring.len()],
            ring[at],
            ring[(at + 1) % ring.len()],
        );
        let point = |vertex: usize| arrangement.vertices[vertex].point;
        if exact.orient(point(previous), point(vertex), point(next))? != Sign::Zero {
            clean += 1;
            index = at + 1;
            continue;
        }
        if !exact.strictly_between(point(previous), point(vertex), point(next))? {
            return Err(GeomError::Degenerate(
                "a ring has a zero-width spike".into(),
            ));
        }
        if keep(VertexId(vertex)) {
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

/// Whether `inner` lies in `outer`: its vertices not on the boundary agree.
/// `Degenerate` when they do not, or when none decides.
fn within(exact: &ExactPoints, inner: &[PointId], outer: &[PointId]) -> GeomResult<bool> {
    let (mut inside, mut outside) = (false, false);
    for &node in inner {
        match exact.locate(node, outer)? {
            Location::Inside => inside = true,
            Location::Outside => outside = true,
            Location::Boundary => {}
        }
    }
    match (inside, outside) {
        (true, false) => Ok(true),
        (false, true) => Ok(false),
        _ => Err(GeomError::Degenerate(
            "rings cross or touch without a deciding vertex".into(),
        )),
    }
}

/// Close `half_edges`, which all carry one class on their left, into rings
/// and nest them.
///
/// `keep` says which vertices a ring must keep even where it runs straight
/// through them, typically those where other boundary meets. Each pair of
/// rings whose containment is decided counts against `budget`. Returns the
/// rings in the order of their first half-edge by `(from, to)` vertex, and the
/// regions, one per outer ring in that order.
pub fn trace(
    points: &ExactPoints,
    arrangement: &Arrangement,
    half_edges: &[HalfEdge],
    keep: &dyn Fn(VertexId) -> bool,
    budget: &mut PairBudget,
    options: &ExecutionOptions,
) -> GeomResult<(Vec<TracedRing>, Vec<TracedRegion>)> {
    let mut by_ends = BTreeMap::<(usize, usize), HalfEdge>::new();
    for &half_edge in half_edges {
        let ends = arrangement.pieces[half_edge.piece.0].ends;
        let (from, to) = if half_edge.forward {
            (ends[0].0, ends[1].0)
        } else {
            (ends[1].0, ends[0].0)
        };
        by_ends.insert((from, to), half_edge);
    }
    let mut fans = BTreeMap::<usize, Vec<(PieceId, bool)>>::new();
    let mut seen = BTreeSet::<(usize, usize)>::new();
    let mut cycles = Vec::<(Vec<(usize, usize)>, Vec<usize>, bool)>::new();
    for &start in by_ends.keys() {
        if seen.contains(&start) {
            continue;
        }
        let mut edges = Vec::new();
        let mut ring = Vec::new();
        let mut edge = start;
        loop {
            options.check_cancelled()?;
            seen.insert(edge);
            edges.push(edge);
            ring.push(edge.0);
            let after = next_vertex(points, arrangement, &mut fans, &by_ends, edge)?;
            edge = (edge.1, after);
            if edge == start {
                break;
            }
            if !by_ends.contains_key(&edge) || seen.contains(&edge) || ring.len() > by_ends.len() {
                return Err(inconsistent());
            }
        }
        let outer = is_outer(arrangement, points, &ring)?;
        simplify(arrangement, points, &mut ring, keep, options)?;
        cycles.push((edges, ring, outer));
    }

    let rings: Vec<Vec<PointId>> = cycles
        .iter()
        .map(|(_, ring, _)| {
            ring.iter()
                .map(|&vertex| arrangement.vertices[vertex].point)
                .collect()
        })
        .collect();
    let boxes: Vec<RingBox> = rings.iter().map(|ring| points.ring_box(ring)).collect();
    let mut owner = vec![None::<usize>; cycles.len()];
    for (hole, (_, _, outer)) in cycles.iter().enumerate() {
        if *outer {
            continue;
        }
        let mut containers = Vec::new();
        for (candidate, (_, _, other_outer)) in cycles.iter().enumerate() {
            if !other_outer || !boxes[hole].possibly_within(&boxes[candidate]) {
                continue;
            }
            options.check_cancelled()?;
            budget.count()?;
            if within(points, &rings[hole], &rings[candidate])? {
                containers.push(candidate);
            }
        }
        let mut innermost = Vec::new();
        for &candidate in &containers {
            let mut holds_all = true;
            for &other in &containers {
                if other == candidate {
                    continue;
                }
                budget.count()?;
                if !within(points, &rings[candidate], &rings[other])? {
                    holds_all = false;
                    break;
                }
            }
            if holds_all {
                innermost.push(candidate);
            }
        }
        match innermost.as_slice() {
            [only] => owner[hole] = Some(*only),
            _ => {
                return Err(GeomError::Degenerate(
                    "a hole has no single innermost outer boundary".into(),
                ))
            }
        }
    }
    if owner
        .iter()
        .zip(&cycles)
        .any(|(owned, (_, _, outer))| !outer && owned.is_none())
    {
        return Err(GeomError::Degenerate(
            "a component has a hole but no outer boundary".into(),
        ));
    }
    let regions = cycles
        .iter()
        .enumerate()
        .filter(|(_, (_, _, outer))| *outer)
        .map(|(index, _)| TracedRegion {
            outer: index,
            holes: owner
                .iter()
                .enumerate()
                .filter(|(_, owned)| **owned == Some(index))
                .map(|(hole, _)| hole)
                .collect(),
        })
        .collect();
    let rings = cycles
        .into_iter()
        .map(|(edges, ring, outer)| TracedRing {
            half_edges: edges.iter().map(|ends| by_ends[ends]).collect(),
            vertices: ring.into_iter().map(VertexId).collect(),
            hole: !outer,
        })
        .collect();
    Ok((rings, regions))
}

/// The vertex the boundary leaves `edge.1` for, after it came along `edge`:
/// the first half-edge met turning clockwise from the way it came.
fn next_vertex(
    points: &ExactPoints,
    arrangement: &Arrangement,
    fans: &mut BTreeMap<usize, Vec<(PieceId, bool)>>,
    half_edges: &BTreeMap<(usize, usize), HalfEdge>,
    (from, at): (usize, usize),
) -> GeomResult<usize> {
    let fan = match fans.entry(at) {
        std::collections::btree_map::Entry::Occupied(known) => known.into_mut(),
        std::collections::btree_map::Entry::Vacant(vacant) => {
            vacant.insert(arrangement.fan(points, VertexId(at))?)
        }
    };
    let to_of = |&(piece, forward): &(PieceId, bool)| {
        arrangement.pieces[piece.0].ends[usize::from(forward)].0
    };
    let index = fan
        .iter()
        .position(|spoke| to_of(spoke) == from)
        .ok_or_else(inconsistent)?;
    for step in 1..=fan.len() {
        let to = to_of(&fan[(index + fan.len() - step) % fan.len()]);
        if half_edges.contains_key(&(at, to)) {
            return Ok(to);
        }
    }
    Err(inconsistent())
}
