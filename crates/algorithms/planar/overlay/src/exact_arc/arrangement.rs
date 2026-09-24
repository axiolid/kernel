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

use super::edge::{crossings, Edge};
use super::point::{same_point, sign, Pred, XPoint};
use super::{assemble, edges_of, monotone, sample, tangent_of, winding, Carrier, Mono, Piece};

/// One piece of the subdivision, before rounding.
pub(crate) struct RawEdge {
    /// Index of the exact start vertex in [`Raw::vertices`].
    pub(crate) from: usize,
    /// Index of the exact end vertex.
    pub(crate) to: usize,
    /// Bulge of the piece, for its direction of travel.
    pub(crate) bulge: f64,
    /// Whether ring `i` contains the region left of the piece.
    pub(crate) left: Vec<bool>,
    /// Whether ring `i` contains the region right of the piece.
    pub(crate) right: Vec<bool>,
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
        self.vertices.iter().map(XPoint::approx).collect()
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

/// Subdivide the plane by `rings`, each simple and counter-clockwise.
pub(crate) fn build(rings: &[ArcRing]) -> Raw {
    let count = rings.len();
    let edges: Vec<Vec<Edge>> = rings.iter().map(edges_of).collect();
    let parts: Vec<Vec<Mono>> = edges
        .iter()
        .map(|own| own.iter().flat_map(monotone).collect())
        .collect();

    // The ring boundary that carries `x`, if any: `(edge index, edge)`.
    let carrier = |ring: usize, x: &XPoint| -> Option<(usize, &Edge)> {
        let (sx, sy) = x.enclosures();
        edges[ring]
            .iter()
            .enumerate()
            .filter(|(_, edge)| edge.bounds.may_hold(sx, sy))
            .find(|(_, edge)| edge.contains(x))
    };

    let mut pieces = Vec::new();
    let mut out = Vec::new();
    for (ring, own) in edges.iter().enumerate() {
        for edge in own {
            // Split at every point another ring's boundary meets this edge.
            // Rings are simple, so a ring never splits itself.
            let mut stops = vec![edge.p0.clone(), edge.p1.clone()];
            for (other, theirs) in edges.iter().enumerate() {
                if other == ring {
                    continue;
                }
                for their in theirs.iter().filter(|t| t.bounds.overlaps(&edge.bounds)) {
                    for x in crossings(edge, their) {
                        if !stops.iter().any(|y| same_point(y, &x)) {
                            stops.push(x);
                        }
                    }
                }
            }
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
                // Splitting is mutual, so a boundary shared by two rings
                // yields identical pieces; the first ring keeps it.
                if (0..ring).any(|earlier| carrier(earlier, &piece.sample).is_some()) {
                    continue;
                }
                let mut left = vec![false; count];
                let mut right = vec![false; count];
                let mut sources = Vec::new();
                for other in 0..count {
                    if let Some((index, their)) = carrier(other, &piece.sample) {
                        // A counter-clockwise ring has its inside on the left
                        // of its own direction of travel.
                        let same = sign(Pred::Tangents {
                            at: &piece.sample,
                            u: &piece.tangent,
                            v: &tangent_of(their),
                            cross: false,
                        }) == Sign::Positive;
                        if same {
                            left[other] = true;
                        } else {
                            right[other] = true;
                        }
                        sources.push((other, index, same));
                    } else if winding(&piece.sample, &parts[other]) != 0 {
                        left[other] = true;
                        right[other] = true;
                    }
                }
                out.push(RawEdge {
                    from: 0,
                    to: 0,
                    bulge: piece.bulge(),
                    left,
                    right,
                    sources,
                });
                pieces.push(piece);
            }
        }
    }

    // One vertex per distinct exact point, so every face built from this
    // subdivision names a shared vertex by the same index.
    let mut vertices: Vec<XPoint> = Vec::new();
    let mut by_lo: Vec<(f64, f64, usize)> = Vec::new();
    let mut reach = 0.0f64;
    let mut intern = |x: &XPoint| -> usize {
        let ((lo, hi), _) = x.enclosures();
        let floor = lo - reach;
        let first = by_lo.partition_point(|entry| entry.0 < floor);
        for &(_, _, index) in by_lo[first..].iter().take_while(|entry| entry.0 <= hi) {
            if same_point(&vertices[index], x) {
                return index;
            }
        }
        let index = vertices.len();
        vertices.push(x.clone());
        let width = hi - lo;
        reach = if width.is_nan() {
            f64::INFINITY
        } else {
            reach.max(width)
        };
        let at = by_lo.partition_point(|entry| entry.0 < lo);
        by_lo.insert(at, (lo, hi, index));
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
