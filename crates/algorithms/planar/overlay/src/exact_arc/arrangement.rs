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

use axiolid_guarantees::Sign;

use crate::arc::ArcRing;
use crate::OverlayError;

use super::boxes::BoxTree;
use super::edge::Bounds;
use super::edge::{crossings, Edge};
use super::point::{same_point, sign, Pred, XPoint};
use super::{assemble, crossing, edges_of, monotone, sample, tangent_of, Carrier, Mono, Piece};

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
    fn winding(&self, s: &XPoint) -> i64 {
        match &self.tree {
            None => self.parts.iter().map(|part| crossing(part, s)).sum(),
            Some(tree) => {
                // Only parts that may pass right of `s` count.
                let (sx, sy) = s.enclosures();
                tree.query(|b| b.may_hold((sx.0, f64::INFINITY), sy))
                    .into_iter()
                    .map(|index| crossing(&self.parts[index], s))
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
    let mut cuts: Vec<Vec<XPoint>> = flat
        .iter()
        .map(|&(ring, edge)| vec![edges[ring][edge].p0.clone(), edges[ring][edge].p1.clone()])
        .collect();
    let add = |stops: &mut Vec<XPoint>, x: &XPoint| {
        if !stops.iter().any(|y| same_point(y, x)) {
            stops.push(x.clone());
        }
    };
    for (index, &(ring, own)) in flat.iter().enumerate() {
        let edge = &edges[ring][own];
        for near in edge_tree.query(|b| b.overlaps(&edge.bounds)) {
            let (other, their) = flat[near];
            if near <= index || other == ring {
                continue;
            }
            for x in crossings(edge, &edges[other][their]) {
                add(&mut cuts[index], &x);
                add(&mut cuts[near], &x);
            }
        }
    }

    let mut pieces = Vec::new();
    let mut out = Vec::new();
    for ((ring, _), (edge, mut stops)) in flat.iter().copied().zip(edges.iter().flatten().zip(cuts))
    {
        stops.sort_by(|a, b| match edge.order(a, b) {
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
        for w in stops.windows(2) {
            let piece = Piece {
                sample: sample(edge, &w[0], &w[1]),
                from: w[0].clone(),
                to: w[1].clone(),
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
            let mut carriers: Vec<(usize, usize, Option<bool>)> = Vec::new();
            let mut kept_earlier = false;
            for index in edge_tree.query(|b| b.may_hold(sx, sy)) {
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
            if kept_earlier {
                continue;
            }
            let mut sources = Vec::new();
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
            // Winding numbers of the other rings, from the parts that
            // may pass right of the sample.
            let mut inside: Vec<usize> = ring_tree
                .query(|b| b.may_hold(sx, sy))
                .into_iter()
                .filter(|other| !carriers.iter().any(|c| c.0 == *other))
                .filter(|&other| windings[other].winding(&piece.sample) != 0)
                .collect();
            inside.sort_unstable();
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
    let mut intern = |x: &XPoint| -> usize {
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
