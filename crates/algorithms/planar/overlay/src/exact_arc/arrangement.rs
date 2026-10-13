//! Exact subdivision of the plane by any number of arc rings (#120).
//!
//! The two-operand boolean keeps only the pieces one operation needs. A
//! stepped or stacked solid needs every piece at once: the section changes
//! with height, and the wall between two heights, the ledge where the
//! section steps and each cap are all bounded by pieces of the same
//! boundaries. Computing them from one subdivision is what makes the faces
//! of such a solid agree on every shared vertex, which separate overlays
//! could not promise after rounding.
//!
//! Every decision here is the same exact sign the boolean uses: where
//! boundaries cross, whether two rings share a piece, which rings contain
//! the region on each side of a piece.
//!
//! Soups whose rings share no edges (#218) are cut into many pieces, and
//! two facts spare most of them a search: a segment no other boundary runs
//! along (no other segment meets it in two points) is carried by its own
//! ring only, since another boundary through one of its pieces would cross
//! it there and make a stop; and along such a segment, a stop where
//! exactly one other segment crosses, interior to both, enters or leaves
//! only that segment's ring, which one following from whether the
//! previous piece lay in it. The rings
//! containing each later piece follow from the previous piece's; the first
//! piece, and the piece after any other stop, count windings from a sample.

use axiolid_guarantees::Sign;

use crate::arc::ArcRing;
use crate::OverlayError;

use super::boxes::BoxTree;
use super::edge::Bounds;
use super::edge::{crossings, Edge};
use super::point::{same_point, sign, Pred, XPoint};
use super::{
    assemble, crossing, edges_of, monotone, sample, tangent_of, Carrier, Mono, Piece, WordMap,
};

/// One piece of the subdivision, before rounding.
pub(crate) struct RawEdge {
    /// Index of the exact start vertex in [`Raw::vertices`].
    pub(crate) from: usize,
    /// Index of the exact end vertex.
    pub(crate) to: usize,
    /// Bulge of the piece, for its direction of travel.
    pub(crate) bulge: f64,
    /// Rings that contain the piece, and so both sides of it, without
    /// carrying it; ascending. Kept sparse: a subdivision of thousands of
    /// rings has each piece inside a handful.
    pub(crate) inside: Vec<usize>,
    /// Rings whose boundary carries the piece: `(ring, edge, same way)`.
    pub(crate) sources: Vec<(usize, usize, bool)>,
}

/// The exact subdivision: vertices, pieces, and their memberships.
pub(crate) struct Raw {
    pieces: Vec<Piece>,
    vertices: Vec<XPoint>,
    pub(crate) edges: Vec<RawEdge>,
}

impl Raw {
    /// Rounded vertex positions, one per distinct exact point.
    pub(crate) fn vertex_positions(&self) -> Vec<axiolid_core::Point2> {
        self.vertices.iter().map(XPoint::rounded).collect()
    }

    /// The rounded position of one vertex (see [`Raw::vertex_positions`]).
    pub(crate) fn vertex_position(&self, index: usize) -> axiolid_core::Point2 {
        self.vertices[index].rounded()
    }

    /// Whether a boundary arriving along one use and leaving along the next
    /// goes straight on: both pieces straight, exactly parallel and the
    /// same way. Such a vertex only marks where another ring's boundary
    /// touched and can be dropped without moving the boundary.
    pub(crate) fn straight_on(&self, arriving: (usize, bool), leaving: (usize, bool)) -> bool {
        let (a, b) = (&self.pieces[arriving.0], &self.pieces[leaving.0]);
        if a.arc.is_some() || b.arc.is_some() {
            return false;
        }
        let travel = |piece: &Piece, reversed: bool| {
            if reversed {
                piece.clone().reversed().tangent
            } else {
                piece.tangent.clone()
            }
        };
        let (u, v) = (travel(a, arriving.1), travel(b, leaving.1));
        let at = if leaving.1 { &b.to } else { &b.from };
        let ask = |cross| {
            sign(Pred::Tangents {
                at,
                u: &u,
                v: &v,
                cross,
            })
        };
        ask(true) == Sign::Zero && ask(false) == Sign::Positive
    }

    /// Link the pieces with `keep[i] = Some(reversed)` into regions.
    ///
    /// Each region is its outer ring and its holes, every ring a list of
    /// `(piece, reversed)` uses in travel order.
    pub(crate) fn regions(
        &self,
        keep: &[Option<bool>],
    ) -> Result<Vec<(RingUses, Vec<RingUses>)>, OverlayError> {
        let kept: Vec<Piece> = self
            .pieces
            .iter()
            .zip(keep)
            .filter_map(|(piece, keep)| {
                keep.map(|reversed| {
                    if reversed {
                        piece.clone().reversed()
                    } else {
                        piece.clone()
                    }
                })
            })
            .collect();
        let uses = |ring: &[Piece]| ring.iter().map(|p| (p.tag, p.flipped)).collect();
        Ok(assemble(kept)?
            .into_iter()
            .map(|(outer, holes)| (uses(&outer), holes.iter().map(|h| uses(h)).collect()))
            .collect())
    }
}

/// A ring as `(piece, reversed)` uses.
pub(crate) type RingUses = Vec<(usize, bool)>;

/// One ring's boundary, split into monotone parts for winding numbers.
struct Winder {
    bounds: Bounds,
    parts: Vec<Mono>,
    /// Over the parts' edge boxes, for rings too long to scan.
    tree: Option<BoxTree>,
}

/// Parts up to which a scan beats a tree query.
const SCAN: usize = 32;

impl Winder {
    fn new(own: &[Edge]) -> Self {
        let mut parts = Vec::new();
        let mut boxes = Vec::new();
        for edge in own {
            for part in monotone(edge) {
                parts.push(part);
                boxes.push(edge.bounds);
            }
        }
        let bounds = Bounds::hull(own.iter().map(|e| e.bounds)).expect("a ring has edges");
        let tree = (parts.len() > SCAN).then(|| BoxTree::new(boxes));
        Self {
            bounds,
            parts,
            tree,
        }
    }

    /// The winding number around `s`, which is not on the boundary.
    /// `scratch` holds query buffers reused across calls.
    fn winding(&self, s: &XPoint, scratch: &mut (Vec<usize>, Vec<usize>)) -> i64 {
        match &self.tree {
            None => self.parts.iter().map(|part| crossing(part, s)).sum(),
            Some(tree) => {
                // Only parts that may pass right of `s` count.
                let (sx, sy) = s.enclosures();
                let (hits, stack) = scratch;
                tree.query_into(|b| b.may_hold((sx.0, f64::INFINITY), sy), hits, stack);
                hits.iter()
                    .map(|&index| crossing(&self.parts[index], s))
                    .sum()
            }
        }
    }
}

/// Subdivide the plane by `rings`, each simple and counter-clockwise.
pub(crate) fn build(rings: &[ArcRing]) -> Raw {
    let edges: Vec<Vec<Edge>> = rings.iter().map(edges_of).collect();

    // Broad phase over every edge and every monotone part of every ring
    // (#173): a subdivision of thousands of rings asks, for each edge,
    // which edges of other rings it may meet, and for each piece, which
    // boundaries may carry it or pass right of it. Scanning every ring for
    // each answer made the cost quadratic in the ring count.
    let flat: Vec<(usize, usize)> = edges
        .iter()
        .enumerate()
        .flat_map(|(ring, own)| (0..own.len()).map(move |edge| (ring, edge)))
        .collect();
    let edge_tree = BoxTree::new(
        flat.iter()
            .map(|&(ring, edge)| edges[ring][edge].bounds)
            .collect(),
    );
    // Each ring's monotone parts, for winding numbers. Only rings whose
    // box holds a sample can contain it; a ring with many parts also gets
    // a tree, so the parts right of the sample are found without a scan.
    let windings: Vec<Winder> = edges.iter().map(|own| Winder::new(own)).collect();
    let ring_tree = BoxTree::new(windings.iter().map(|w| w.bounds).collect());

    // Split every edge at every point another ring's boundary meets it.
    // Rings are simple, so a ring never splits itself. The points two
    // edges share are the same seen from either, so each pair is asked
    // once.
    //
    // Each stop carries the edge (flat index) whose boundary crosses this
    // one there when that crossing is simple (#218): two segments, one
    // point in common, interior to both, and no other boundary through
    // it. Only then does passing the stop along the edge change exactly
    // one ring's winding number, by one, so that ring alone is entered or
    // left; any other stop carries `None`.
    let mut cuts: Vec<Vec<(XPoint, Option<usize>)>> = flat
        .iter()
        .map(|&(ring, edge)| {
            let edge = &edges[ring][edge];
            let mut stops = Vec::with_capacity(4);
            stops.extend([(edge.p0.clone(), None), (edge.p1.clone(), None)]);
            stops
        })
        .collect();
    let add = |stops: &mut Vec<(XPoint, Option<usize>)>, x: &XPoint, tag: Option<usize>| {
        match stops.iter_mut().find(|y| same_point(&y.0, x)) {
            // Met again, from another edge or at an end: not simple.
            Some(held) => held.1 = None,
            None => stops.push((x.clone(), tag)),
        }
    };
    // Segments another ring's boundary runs along for a stretch: their
    // pieces may have other carriers.
    let mut overlapped = vec![false; flat.len()];
    // Query buffers, reused for every edge and piece.
    let (mut hits, mut stack) = (Vec::new(), Vec::new());
    let mut ring_hits = Vec::new();
    let mut scratch = (Vec::new(), Vec::new());
    for (index, &(ring, own)) in flat.iter().enumerate() {
        let edge = &edges[ring][own];
        edge_tree.query_into(|b| b.overlaps(&edge.bounds), &mut hits, &mut stack);
        for &near in &hits {
            let (other, their) = flat[near];
            if near <= index || other == ring {
                continue;
            }
            let theirs = &edges[other][their];
            let xs = crossings(edge, theirs);
            let straight = !edge.is_arc() && !theirs.is_arc();
            // Two segments with two points in common share the stretch
            // between them.
            if straight && xs.len() > 1 {
                overlapped[index] = true;
                overlapped[near] = true;
            }
            // A point at an end of either segment is met again: the end
            // is a stop already, or a vertex the next edge of its ring
            // meets too, so `add` clears its tag.
            let simple = straight && xs.len() == 1;
            for x in &xs {
                add(&mut cuts[index], x, simple.then_some(near));
                add(&mut cuts[near], x, simple.then_some(index));
            }
        }
    }

    // At most one piece between consecutive stops.
    let count: usize = cuts.iter().map(|stops| stops.len() - 1).sum();
    let mut pieces = Vec::with_capacity(count);
    let mut out = Vec::with_capacity(count);
    let cut_edges = flat.iter().copied().zip(edges.iter().flatten().zip(cuts));
    for (index, ((ring, own), (edge, mut stops))) in cut_edges.enumerate() {
        stops.sort_by(|a, b| match edge.order(&a.0, &b.0) {
            Sign::Negative => std::cmp::Ordering::Less,
            Sign::Positive => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        });
        let split = stops.len() > 2;
        let arc = match &edge.carrier {
            Carrier::Segment => None,
            Carrier::Arc { circle, turn, .. } => Some((circle.clone(), *turn)),
        };
        let whole = match &edge.carrier {
            Carrier::Arc { bulge, .. } if !split => Some(bulge.to_f64()),
            _ => None,
        };
        // A segment no other boundary runs along (#218). Another ring's
        // boundary through a point inside one of its pieces would cross
        // it there, which would be a stop, so its own ring is the only
        // carrier of every piece; and between two pieces the rings
        // containing them differ only by the ring crossing at a simple
        // stop, entered or left by the direction of the crossing.
        let lone = !edge.is_arc() && !overlapped[index];
        // The rings containing the previous piece of a lone segment.
        let mut before: Option<Vec<usize>> = None;
        for w in stops.windows(2) {
            let piece = Piece {
                sample: sample(edge, &w[0].0, &w[1].0),
                from: w[0].0.clone(),
                to: w[1].0.clone(),
                tangent: tangent_of(edge),
                arc: arc.clone(),
                whole_bulge: whole,
                tag: pieces.len(),
                flipped: false,
            };
            let (sx, sy) = piece.sample.enclosures();
            // The boundaries carrying the piece: `(ring, edge index,
            // twin)`, at most one edge per ring, since the sample lies
            // strictly inside the piece and every ring is simple.
            // An edge with the same two ends as this one (a triangle
            // soup's shared and repeated edges) carries every piece of
            // it, the same way when its ends come in the same order.
            let mut carriers: Vec<(usize, usize, Option<bool>)> = Vec::with_capacity(2);
            if lone && !cfg!(debug_assertions) {
                carriers.push((ring, own, Some(true)));
            } else {
                // Debug builds search for a lone segment too, and check
                // the shortcut against the search below.
                let mut kept_earlier = false;
                edge_tree.query_into(|b| b.may_hold(sx, sy), &mut hits, &mut stack);
                for &index in &hits {
                    let (other, their) = flat[index];
                    if carriers.iter().any(|c| c.0 == other) {
                        continue;
                    }
                    let theirs = &edges[other][their];
                    let twin = edge.same_segment(theirs);
                    if twin.is_some() || theirs.contains(&piece.sample) {
                        // Splitting is mutual, so a boundary shared by two
                        // rings yields identical pieces; the first ring
                        // keeps it. Candidates come ring by ring.
                        if other < ring {
                            kept_earlier = true;
                            break;
                        }
                        carriers.push((other, their, twin));
                    }
                }
                debug_assert!(
                    !lone || (!kept_earlier && carriers == [(ring, own, Some(true))]),
                    "a lone segment is carried by its own ring only"
                );
                if kept_earlier {
                    continue;
                }
            }
            let mut sources = Vec::with_capacity(carriers.len());
            carriers.sort_unstable();
            for &(other, index, twin) in &carriers {
                // A counter-clockwise ring has its inside on the left
                // of its own direction of travel.
                let same = twin.unwrap_or_else(|| {
                    sign(Pred::Tangents {
                        at: &piece.sample,
                        u: &piece.tangent,
                        v: &tangent_of(&edges[other][index]),
                        cross: false,
                    }) == Sign::Positive
                });
                sources.push((other, index, same));
            }
            // Across a simple stop of a lone segment, from the previous
            // piece: the crossing ring's winding number changes by one,
            // between zero and one for a simple ring, so it is left when
            // it held the previous piece and entered otherwise.
            let carried = match (&before, w[0].1) {
                (Some(previous), Some(crossing)) => {
                    let other = flat[crossing].0;
                    let mut inside = Vec::with_capacity(previous.len() + 1);
                    inside.extend_from_slice(previous);
                    match inside.binary_search(&other) {
                        Ok(at) => {
                            inside.remove(at);
                        }
                        Err(at) => inside.insert(at, other),
                    }
                    Some(inside)
                }
                _ => None,
            };
            let mut counted = || -> Vec<usize> {
                // Winding numbers of the other rings, from the parts that
                // may pass right of the sample.
                ring_tree.query_into(|b| b.may_hold(sx, sy), &mut ring_hits, &mut stack);
                // Ascending, as the query returns them.
                ring_hits
                    .iter()
                    .copied()
                    .filter(|other| !carriers.iter().any(|c| c.0 == *other))
                    .filter(|&other| windings[other].winding(&piece.sample, &mut scratch) != 0)
                    .collect()
            };
            let inside = match carried {
                Some(inside) => {
                    debug_assert_eq!(inside, counted(), "rings carried across a stop");
                    inside
                }
                None => counted(),
            };
            before = lone.then(|| inside.clone());
            out.push(RawEdge {
                from: 0,
                to: 0,
                bulge: piece.bulge(),
                inside,
                sources,
            });
            pieces.push(piece);
        }
    }

    // One vertex per distinct exact point, so every face built from this
    // subdivision names a shared vertex by the same index.
    let mut vertices: Vec<XPoint> = Vec::new();
    // Interned points by the low end of their `x` box, in a tree: a point
    // equal to `x` has a box meeting `x`'s, so its low end lies within the
    // widest box seen of `x`'s.
    let mut by_lo: std::collections::BTreeMap<(u64, usize), (f64, usize)> =
        std::collections::BTreeMap::new();
    // Ordered keys for finite floats: flip the sign bit, or every bit of a
    // negative number.
    let key = |v: f64| {
        let bits = v.to_bits();
        if bits >> 63 == 1 {
            !bits
        } else {
            bits | (1 << 63)
        }
    };
    let mut reach = 0.0f64;
    let mut interned = |x: &XPoint| -> usize {
        let ((lo, hi), _) = x.enclosures();
        let floor = lo - reach;
        if floor.is_finite() && hi.is_finite() {
            for (_, &(_, index)) in by_lo.range((key(floor), 0)..=(key(hi), usize::MAX)) {
                if same_point(&vertices[index], x) {
                    return index;
                }
            }
        } else if let Some(index) = vertices.iter().position(|v| same_point(v, x)) {
            return index;
        }
        let index = vertices.len();
        vertices.push(x.clone());
        let width = hi - lo;
        reach = if width.is_nan() {
            f64::INFINITY
        } else {
            reach.max(width)
        };
        if lo.is_finite() {
            by_lo.insert((key(lo), index), (hi, index));
        }
        index
    };
    // Copies of one point share its identity (#218): a stop ends one piece
    // and starts the next, and a crossing splits two edges, so most points
    // are met again as the same copy and need no search.
    let mut by_id: WordMap<usize, usize> = WordMap::default();
    let mut intern = |x: &XPoint| -> usize {
        if let Some(&index) = by_id.get(&x.id()) {
            return index;
        }
        let index = interned(x);
        by_id.insert(x.id(), index);
        index
    };
    for (piece, edge) in pieces.iter().zip(out.iter_mut()) {
        edge.from = intern(&piece.from);
        edge.to = intern(&piece.to);
    }
    Raw {
        pieces,
        vertices,
        edges: out,
    }
}
