//! Bridges across gaps narrower than the closing width.
//!
//! On the arrangement of the inputs alone, a candidate is a vertex `V` and a
//! piece `P` it is not on. `C` is the point of `P` nearest to `V`: an end of
//! `P`, or the exact foot of the perpendicular. The pair is bridged when `V`
//! is strictly nearer to `C` than `gap` and the open segment `V`-`C` meets
//! no piece. The minimum distance of two segments that do not cross is
//! reached at an end of one of them, so every narrowest place of a gap is
//! the bridge of some such pair, and nothing is chosen between candidates.
//!
//! Pairs of vertices are found from both sides and kept once.

use std::cmp::Ordering;
use std::collections::BTreeSet;

use axiolid_contracts::{ExecutionOptions, GeomResult, Sign};
use axiolid_core::Point2;
use axiolid_exact_arrangement::{Arrangement, ExactPoints, Nearest, PairBudget, PointId, VertexId};

use super::PieceSource;
use crate::exact_arc::boxes::BoxTree;
use crate::exact_arc::Bounds;

/// A bridge before rounding: its exact ends and what each joins.
#[derive(Debug, Clone)]
pub(super) struct Found {
    pub ends: [PointId; 2],
    pub joins: [Vec<PieceSource>; 2],
}

/// The sorted, deduplicated sources of the segments over a piece.
pub(super) fn piece_sources(
    arrangement: &Arrangement,
    sources: &[PieceSource],
) -> Vec<Vec<PieceSource>> {
    arrangement
        .pieces
        .iter()
        .map(|piece| sorted(piece.segments.iter().map(|&segment| sources[segment])))
        .collect()
}

pub(super) fn sorted(sources: impl IntoIterator<Item = PieceSource>) -> Vec<PieceSource> {
    let mut list: Vec<PieceSource> = sources.into_iter().collect();
    list.sort_unstable();
    list.dedup();
    list
}

/// The corners of the box of an exact point, for a [`Bounds`].
fn corners(points: &ExactPoints, ends: &[PointId]) -> Vec<Point2> {
    ends.iter()
        .flat_map(|&end| {
            let [x_lo, y_lo, x_hi, y_hi] = points.bounds(end);
            [Point2::new(x_lo, y_lo), Point2::new(x_hi, y_hi)]
        })
        .collect()
}

/// Whether the open segment `v`-`c` meets the closed segment `a`-`b`.
///
/// A proper crossing, an end of `a`-`b` strictly inside `v`-`c`, or a
/// stretch of both lying on one line. Touching at `v` or `c` alone does not
/// count: a bridge may start and end on pieces.
pub(super) fn meets_open_segment(
    points: &ExactPoints,
    (v, c): (PointId, PointId),
    (a, b): (PointId, PointId),
) -> GeomResult<bool> {
    let (o1, o2) = (points.orient(v, c, a)?, points.orient(v, c, b)?);
    let (o3, o4) = (points.orient(a, b, v)?, points.orient(a, b, c)?);
    let proper = [o1, o2, o3, o4].iter().all(|&sign| sign != Sign::Zero) && o1 != o2 && o3 != o4;
    if proper {
        return Ok(true);
    }
    if o1 == Sign::Zero && points.strictly_between(v, a, c)? {
        return Ok(true);
    }
    if o2 == Sign::Zero && points.strictly_between(v, b, c)? {
        return Ok(true);
    }
    Ok(o1 == Sign::Zero
        && o2 == Sign::Zero
        && points.within_segment(v, a, b)?
        && points.within_segment(c, a, b)?)
}

/// Where a candidate's nearest point lies.
enum Target {
    Vertex(VertexId),
    Foot,
}

/// Every bridge of `arrangement`, the arrangement of the input segments
/// whose sources are `sources`, for a closing width `gap`. None when `gap`
/// is zero.
pub(super) fn find(
    points: &mut ExactPoints,
    arrangement: &Arrangement,
    sources: &[PieceSource],
    gap: f64,
    budget: &mut PairBudget,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Found>> {
    if gap == 0.0 {
        return Ok(Vec::new());
    }
    let own = piece_sources(arrangement, sources);
    let mut at_vertex: Vec<Vec<usize>> = vec![Vec::new(); arrangement.vertices.len()];
    for (index, piece) in arrangement.pieces.iter().enumerate() {
        at_vertex[piece.ends[0].0].push(index);
        at_vertex[piece.ends[1].0].push(index);
    }
    let joined = |vertex: VertexId| -> Vec<PieceSource> {
        sorted(
            at_vertex[vertex.0]
                .iter()
                .flat_map(|&piece| own[piece].iter().copied()),
        )
    };
    let end_points = |piece: usize| {
        arrangement.pieces[piece]
            .ends
            .map(|end| arrangement.vertices[end.0].point)
    };
    // The pieces' boxes, grown by the gap to find candidates, and plain to
    // find what stands between.
    let boxes = |grow: f64| -> Vec<Bounds> {
        (0..arrangement.pieces.len())
            .map(|piece| Bounds::around(&corners(points, &end_points(piece)), grow))
            .collect()
    };
    let wide = BoxTree::new(boxes(gap));
    let plain = BoxTree::new(boxes(0.0));

    let mut found = Vec::new();
    let mut pairs = BTreeSet::<(VertexId, VertexId)>::new();
    let (mut near, mut between, mut stack) = (Vec::new(), Vec::new(), Vec::new());
    for (index, vertex) in arrangement.vertices.iter().enumerate() {
        options.check_cancelled()?;
        let (here, here_id) = (vertex.point, VertexId(index));
        let own_box = Bounds::around(&corners(points, &[here]), 0.0);
        wide.query_into(|b| b.overlaps(&own_box), &mut near, &mut stack);
        for &piece in &near {
            let ends = arrangement.pieces[piece].ends;
            if ends.contains(&here_id) {
                continue;
            }
            budget.count()?;
            let [start, end] = end_points(piece);
            let (target, closer) = match points.nearest_on_segment(here, start, end)? {
                Nearest::Interior => (
                    Target::Foot,
                    points.compare_line_distance(here, start, end, gap)?,
                ),
                side => {
                    let to = ends[usize::from(side == Nearest::End)];
                    (
                        Target::Vertex(to),
                        points.compare_distance(here, arrangement.vertices[to.0].point, gap)?,
                    )
                }
            };
            if closer != Ordering::Less {
                continue;
            }
            let key = match target {
                Target::Vertex(to) => Some((here_id.min(to), here_id.max(to))),
                Target::Foot => None,
            };
            if key.is_some_and(|key| pairs.contains(&key)) {
                continue;
            }
            let there = match target {
                Target::Vertex(to) => arrangement.vertices[to.0].point,
                Target::Foot => points.foot(here, start, end)?,
            };
            let span = Bounds::around(&corners(points, &[here, there]), 0.0);
            plain.query_into(|b| b.overlaps(&span), &mut between, &mut stack);
            let mut blocked = false;
            for &other in &between {
                let [a, b] = end_points(other);
                if meets_open_segment(points, (here, there), (a, b))? {
                    blocked = true;
                    break;
                }
            }
            if blocked {
                continue;
            }
            if let Some(key) = key {
                pairs.insert(key);
            }
            let there_joins = match target {
                Target::Vertex(to) => joined(to),
                Target::Foot => own[piece].clone(),
            };
            found.push(Found {
                ends: [here, there],
                joins: [joined(here_id), there_joins],
            });
        }
    }
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Whether the open segment `(0, 0)`-`(4, 0)` meets the piece between
    /// the two positions.
    fn meets(a: (f64, f64), b: (f64, f64)) -> bool {
        let mut points = ExactPoints::new();
        let mut at = |(x, y): (f64, f64)| points.push_flat(Point2::new(x, y)).unwrap();
        let ids = [at((0.0, 0.0)), at((4.0, 0.0)), at(a), at(b)];
        meets_open_segment(&points, (ids[0], ids[1]), (ids[2], ids[3])).unwrap()
    }

    #[test]
    fn a_crossing_piece_meets_it() {
        assert!(meets((2.0, -1.0), (2.0, 1.0)));
        assert!(meets((1.0, -1.0), (3.0, 1.0)));
    }

    #[test]
    fn a_piece_touching_an_end_or_standing_off_does_not() {
        // Meeting only the bridge's own ends: a bridge starts and ends on
        // pieces.
        assert!(!meets((0.0, 0.0), (0.0, 3.0)));
        assert!(!meets((4.0, -2.0), (4.0, 0.0)));
        // Ending on the open segment itself blocks: a bridge may not run
        // into another piece.
        assert!(meets((2.0, 0.0), (2.0, 3.0)));
        // Clear of it.
        assert!(!meets((2.0, 1.0), (3.0, 2.0)));
        assert!(!meets((5.0, -1.0), (5.0, 1.0)));
    }

    #[test]
    fn a_collinear_piece_meets_it_when_it_overlaps_and_not_beyond_the_ends() {
        assert!(meets((1.0, 0.0), (2.0, 0.0)));
        assert!(meets((-1.0, 0.0), (1.0, 0.0)));
        // The piece covering the whole bridge.
        assert!(meets((-1.0, 0.0), (5.0, 0.0)));
        // Beyond either end, or meeting it only at an end.
        assert!(!meets((4.0, 0.0), (6.0, 0.0)));
        assert!(!meets((-3.0, 0.0), (0.0, 0.0)));
        assert!(!meets((5.0, 0.0), (6.0, 0.0)));
    }
}
