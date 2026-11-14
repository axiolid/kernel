//! Rings: snap rounding, certification, provenance and order.
//!
//! The exact rings of every bounded free face, closed gaps included, are cut
//! into edges, the straight runs between kept vertices, and each edge goes
//! once into `snap_round`, which routes it through the rounded position of
//! every vertex whose cell it passes through (ADR 0100). A class's boundary
//! is netted on the snapped fragments, so a part of zero width cancels, then
//! traced again into rings by the same rule and certified. A face that
//! rounding pinches apart is returned as its parts; a face, hole or closed
//! gap that vanishes is counted.
//!
//! Every ring edge lists the inputs it lies on. An exact edge carries the
//! sources of its pieces, which are equal along it by the kept-vertex rule;
//! a snapped piece carries the union over the edges routed over it, and a
//! vertex is dropped from a ring only where exactly two such pieces meet in a
//! straight line with equal lists. Bridges are renumbered into the order of
//! the output (rounded ends), and rings start at their smallest point.

use std::collections::BTreeMap;

use axiolid_contracts::{DisplacementBound, ExecutionOptions, GeomError, GeomResult};
use axiolid_core::Point2;
use axiolid_exact_arrangement::{
    canonical_start, check_winding, point_order, sequence_order, snap_round, trace, Arrangement,
    ExactPoints, HalfEdge, PairBudget, Segment, Snapped, TracedRing, VertexId,
};

use super::bridges::{sorted, Found};
use super::classify::{kept_vertices, Classified};
use super::{Bridge, FreeFace, FreeRing, PieceSource, SEGMENT_BYTES};

pub(super) struct Built {
    pub faces: Vec<FreeFace>,
    pub closed_gaps: Vec<FreeFace>,
    pub bridges: Vec<Bridge>,
    pub snap_displacement: DisplacementBound,
    pub rerouted_edges: usize,
    pub collapsed_faces: usize,
    pub collapsed_holes: usize,
    pub collapsed_gaps: usize,
}

/// The elements of a class's boundary over the edges, and the exact ring of
/// each element.
struct Chain {
    elements: Vec<(usize, bool)>,
    rings: Vec<usize>,
}

fn inconsistent(message: &str) -> GeomError {
    GeomError::Degenerate(message.into())
}

/// The bridges of the output, by the rounded ends, and for each bridge found
/// the number it has there (`None` for one that lies in an obstacle).
fn bridge_records(
    points: &ExactPoints,
    found: &[Found],
    alive: &[bool],
) -> GeomResult<(Vec<Bridge>, Vec<Option<usize>>)> {
    let mut rounded: Vec<(usize, [Point2; 2], [Vec<PieceSource>; 2])> = Vec::new();
    for (index, bridge) in found.iter().enumerate() {
        if !alive[index] {
            continue;
        }
        let round = |end| -> GeomResult<Point2> {
            let at = points.round(end)?;
            Ok(Point2::new(at.x + 0.0, at.y + 0.0))
        };
        let mut ends = [round(bridge.ends[0])?, round(bridge.ends[1])?];
        let mut joins = bridge.joins.clone();
        match point_order(&ends[0], &ends[1]) {
            std::cmp::Ordering::Greater => {
                ends.swap(0, 1);
                joins.swap(0, 1);
            }
            std::cmp::Ordering::Equal => {
                let both = sorted(joins[0].iter().chain(&joins[1]).copied());
                joins = [both.clone(), both];
            }
            std::cmp::Ordering::Less => {}
        }
        rounded.push((index, ends, joins));
    }
    rounded.sort_by(|a, b| {
        point_order(&a.1[0], &b.1[0])
            .then_with(|| point_order(&a.1[1], &b.1[1]))
            .then(a.0.cmp(&b.0))
    });
    let mut bridges: Vec<Bridge> = Vec::new();
    let mut number = vec![None; found.len()];
    for (index, ends, joins) in rounded {
        let same = bridges.last().is_some_and(|last| last.ends == ends);
        if same {
            let last = bridges.last_mut().expect("a previous bridge");
            for (own, extra) in last.joins.iter_mut().zip(&joins) {
                *own = sorted(own.iter().chain(extra).copied());
            }
        } else {
            bridges.push(Bridge {
                length: (ends[1] - ends[0]).length(),
                ends,
                joins,
            });
        }
        number[index] = Some(bridges.len() - 1);
    }
    Ok((bridges, number))
}

/// Per kept vertex of a ring, the sources of the edge that leaves it: those
/// of the first half-edge of the cycle that does. The half-edges before the
/// first kept vertex run on to the last edge, and the vertices dropped
/// between kept ones carry equal sources.
fn edge_sources<'a>(
    arrangement: &Arrangement,
    ring: &TracedRing,
    of: impl Fn(HalfEdge) -> &'a [PieceSource],
) -> GeomResult<Vec<&'a [PieceSource]>> {
    let mut found = Vec::with_capacity(ring.vertices.len());
    for &half_edge in &ring.half_edges {
        let ends = arrangement.pieces[half_edge.piece.0].ends;
        let tail = ends[usize::from(!half_edge.forward)];
        if found.len() < ring.vertices.len() && tail == ring.vertices[found.len()] {
            let sources = of(half_edge);
            if sources.is_empty() {
                return Err(inconsistent("a ring edge lies on no input"));
            }
            found.push(sources);
        }
    }
    if found.len() == ring.vertices.len() {
        Ok(found)
    } else {
        Err(inconsistent("a ring's edges do not follow its vertices"))
    }
}

/// The net boundary of each exact ring on its own: `chain` holds the
/// elements of every ring of a class and `element_rings` the ring of each.
fn ring_boundaries(
    snapped: &Snapped,
    chain: &[(usize, bool)],
    element_rings: &[usize],
    rings: usize,
) -> GeomResult<Vec<Vec<HalfEdge>>> {
    (0..rings)
        .map(|ring| {
            let own: Vec<(usize, bool)> = chain
                .iter()
                .zip(element_rings)
                .filter(|&(_, &of)| of == ring)
                .map(|(&element, _)| element)
                .collect();
            Ok(snapped.class_boundary(&own)?.half_edges)
        })
        .collect()
}

/// A snapped ring as the output has it: points and the sources per edge,
/// rotated to start at the smallest point.
fn output_ring(
    snapped: &Snapped,
    ring: &TracedRing,
    sources: &[Vec<PieceSource>],
) -> GeomResult<FreeRing> {
    let at = ring
        .vertices
        .iter()
        .map(|&vertex| snapped.point(vertex))
        .collect::<GeomResult<Vec<Point2>>>()?;
    let per_edge = edge_sources(&snapped.arrangement, ring, |half_edge| {
        &sources[half_edge.piece.0]
    })?;
    let start = canonical_start(&at);
    let count = at.len();
    Ok(FreeRing {
        points: (0..count).map(|i| at[(start + i) % count]).collect(),
        sources: (0..count)
            .map(|i| per_edge[(start + i) % count].to_vec())
            .collect(),
    })
}

pub(super) fn build(
    points: &ExactPoints,
    arrangement: &Arrangement,
    classified: &Classified,
    found: &[Found],
    budget: &mut PairBudget,
    options: &ExecutionOptions,
) -> GeomResult<Built> {
    let (bridges, number) = bridge_records(points, found, &classified.bridge_alive)?;
    let renumber = |sources: &[PieceSource]| -> GeomResult<Vec<PieceSource>> {
        let mut list = Vec::with_capacity(sources.len());
        for &source in sources {
            list.push(match source {
                PieceSource::Bridge { bridge } => PieceSource::Bridge {
                    bridge: number[bridge]
                        .ok_or_else(|| inconsistent("a ring edge lies on a discarded bridge"))?,
                },
                other => other,
            });
        }
        Ok(sorted(list))
    };

    // The edges of all exact rings, once each, and each class's chain over
    // them with the exact ring of every element.
    let vertex_point = |vertex: VertexId| arrangement.vertices[vertex.0].point;
    let mut edge_of: BTreeMap<(VertexId, VertexId), usize> = BTreeMap::new();
    let mut edges: Vec<Segment> = Vec::new();
    let mut edge_lists: Vec<Vec<PieceSource>> = Vec::new();
    let mut chains: Vec<Chain> = Vec::new();
    for class in &classified.classes {
        options.check_cancelled()?;
        let (mut chain, mut element_rings) = (Vec::new(), Vec::new());
        for (index, ring) in class.rings.iter().enumerate() {
            let lists = edge_sources(arrangement, ring, |half_edge| {
                &classified.piece_sources[half_edge.piece.0]
            })?;
            let count = ring.vertices.len();
            for (at, &list) in lists.iter().enumerate() {
                let (from, to) = (ring.vertices[at], ring.vertices[(at + 1) % count]);
                let ends = (from.min(to), from.max(to));
                let edge = match edge_of.get(&ends) {
                    Some(&edge) => edge,
                    None => {
                        edges.push(Segment {
                            ends: [vertex_point(ends.0), vertex_point(ends.1)],
                        });
                        edge_lists.push(renumber(list)?);
                        edge_of.insert(ends, edges.len() - 1);
                        edges.len() - 1
                    }
                };
                chain.push((edge, from < to));
                element_rings.push(index);
            }
        }
        chains.push(Chain {
            elements: chain,
            rings: element_rings,
        });
    }

    let mut built = Built {
        faces: Vec::new(),
        closed_gaps: Vec::new(),
        bridges,
        snap_displacement: DisplacementBound::ZERO,
        rerouted_edges: 0,
        collapsed_faces: 0,
        collapsed_holes: 0,
        collapsed_gaps: 0,
    };
    if edges.is_empty() {
        return Ok(built);
    }
    let snapped = snap_round(points, &edges, budget, edges.len() * SEGMENT_BYTES, options)?;
    built.snap_displacement = snapped.displacement;
    built.rerouted_edges = snapped.rerouted;

    // The inputs each snapped piece lies on, and the vertices to keep.
    let mut lists: Vec<Vec<PieceSource>> = vec![Vec::new(); snapped.arrangement.pieces.len()];
    for (edge, route) in snapped.routes.iter().enumerate() {
        for &(piece, _) in route {
            lists[piece.0].extend(edge_lists[edge].iter().copied());
        }
    }
    for list in &mut lists {
        *list = sorted(list.iter().copied());
    }
    let boundary: Vec<bool> = lists.iter().map(|list| !list.is_empty()).collect();
    let keep_after = kept_vertices(&snapped.points, &snapped.arrangement, &boundary, &lists)?;

    for (
        class,
        Chain {
            elements: chain,
            rings: element_rings,
        },
    ) in classified.classes.iter().zip(&chains)
    {
        options.check_cancelled()?;
        // Which exact rings survive: each is netted alone, so a ring that
        // collapses onto another surviving edge of the class leaves nothing
        // of its own.
        let own = ring_boundaries(&snapped, chain, element_rings, class.rings.len())?;
        if own[class.outer].is_empty() {
            if class.gap {
                built.collapsed_gaps += 1;
            } else {
                built.collapsed_faces += 1;
            }
        } else {
            built.collapsed_holes += class
                .rings
                .iter()
                .enumerate()
                .filter(|&(index, ring)| ring.hole && own[index].is_empty())
                .count();
        }

        let net = snapped.class_boundary(chain)?;
        let (rings, traced) = trace(
            &snapped.points,
            &snapped.arrangement,
            &net.half_edges,
            &|vertex| keep_after[vertex.0],
            budget,
            options,
        )?;
        check_winding(
            &snapped.points,
            &snapped.arrangement,
            &rings,
            &traced,
            budget,
            options,
        )?;
        for region in &traced {
            let outer = output_ring(&snapped, &rings[region.outer], &lists)?;
            let mut holes = region
                .holes
                .iter()
                .map(|&hole| output_ring(&snapped, &rings[hole], &lists))
                .collect::<GeomResult<Vec<_>>>()?;
            holes.sort_by(|a, b| sequence_order(&a.points, &b.points));
            let face = FreeFace { outer, holes };
            if class.gap {
                built.closed_gaps.push(face);
            } else {
                built.faces.push(face);
            }
        }
    }
    for list in [&mut built.faces, &mut built.closed_gaps] {
        list.sort_by(|a, b| sequence_order(&a.outer.points, &b.outer.points));
    }
    Ok(built)
}
