//! The union of Minkowski terms, read from the term edges that can bound it
//! (#344, ADR 0108).
//!
//! A sum or erosion is built from many convex terms: the hull of the shape
//! moved to each end of each boundary edge, a copy of the set moved by the
//! shape's anchor, and so on (`minkowski`). Put into one arrangement, a
//! shape larger than the region's features makes every hull overlap
//! nearly every other, so the arrangement had `O(n^2)` crossings and every
//! piece lay inside `O(n)` hulls. Most of those edges lie inside another
//! hull and bound nothing. This module keeps the terms and the rule
//! (`Rule::holds`) and changes only how the union is evaluated:
//!
//! 1. **Certify and drop interior term edges.** An edge of any term ring is
//!    left out when some convex piece ring (a ring from `first` on), other
//!    than the edge's own, holds its open segment strictly inside: both
//!    ends in the closed piece, and no side line of the piece holding both
//!    ends. The test is the exact orientation of `f64` points.
//! 2. **Cut only the kept edges, at every term vertex**, with the exact
//!    core's partial mode (`exact_arc::arrangement::build_partial`).
//! 3. **Read each piece's sides with the same rule.** The rings holding a
//!    side are the carriers whose inside it is, then, when a side holds no
//!    piece ring, either a piece ring other than a carrier that holds the
//!    piece's sample strictly (both sides), or the set rings with a nonzero
//!    winding number there.
//! 4. **Link, round, present and settle as `region_of` does**
//!    (`minkowski::region_from`).
//!
//! # Why the result is the same, bit for bit
//!
//! `rings` are oriented as `ArcArrangement::new` orients them. Rings from
//! `first` on are convex, counter-clockwise and of nonzero area. The old
//! route's answer on a side of a piece is `holds(L)`, `L` the ascending
//! list of the rings whose inside holds that side.
//!
//! - **(A) Pieces absorb.** `holds(L)` is the same for every `L` holding a
//!   ring from `first` on: a sum holds such a point, an erosion leaves it
//!   out. So any point strictly inside a piece ring has a fixed answer,
//!   and both sides of anything through it agree.
//! - **(B) Dropped edges never bound the result.** A dropped edge's open
//!   segment lies in the open interior of a piece ring (a half-plane
//!   argument: both ends in each closed half-plane, not both on its line).
//!   By (A) the result is constant across it. A ring vertex of the
//!   result's boundary with only dropped edges through it has a small disc
//!   split into sectors that no ring boundary crosses, so each sector has
//!   one answer, and across each dropped edge it is the absorbing one on
//!   both sides; so every sector has it, and the vertex is not on the
//!   boundary. Every boundary point lies on a kept edge.
//! - **(C) Each piece side is read correctly.** Every ring vertex lying on
//!   a kept edge is a stop, and so is every meeting of two kept edges. A
//!   piece's sample is no ring vertex, and any ring boundary through it is
//!   a carrier or the open segment of a dropped edge. In the second case
//!   the sample lies strictly inside that edge's certifier, a piece that is
//!   no carrier, and (A) decides. Otherwise no non-carrier ring passes
//!   through the sample, so each holds both sides exactly when its winding
//!   number there is nonzero, and no non-carrier piece holds either side.
//!   The lists built differ from the old route's only by piece rings
//!   dropped while a carrier piece is present or none is, which does not
//!   change `holds` (A).
//! - **(D) The kept pieces are the old ones.** A stop inside a piece of the
//!   result's boundary comes, in the old subdivision, from a ring edge
//!   meeting the carrier there. If that edge is dropped, the meeting is
//!   its end, a ring vertex and a stop here too, because an inner point of
//!   a dropped edge is strictly inside a piece and cannot be on the
//!   boundary (B). If it is kept, the meeting is a stop here too.
//!   Conversely every stop here is a stop there. The carriers of a boundary
//!   piece contain it, so none is dropped. So its owner (the first ring
//!   carrying it), direction and tangent agree, and the kept pieces come
//!   out in the old order (rings in order, edges in order, pieces along
//!   each edge).
//! - **(E) Same linking, same rounding.** `assemble` links the same pieces
//!   in the same pool order, so rings, their starts and their order agree.
//!   Holes go to the same outers, since kept pieces share only end
//!   points. Vertices are rounded once (`XPoint::rounded`, correctly
//!   rounded), so equal exact points round equally whatever their recipe;
//!   `presented` and `settle_validated` then get identical input.
//!
//! # Preconditions
//!
//! Only `minkowski.rs` may call [`region`]: it builds the terms, so the rings
//! from `first` on are convex hulls (hulls, chord hulls, trimmed hulls and
//! parallelograms, refused there if of zero area) and `Rule::holds` gives
//! one answer for every list holding one of them. Neither is checked in a
//! release build, since a caller that broke them would get a wrong region
//! and not a refusal; debug builds assert both. A caller outside
//! `minkowski.rs` would need typed refusals first, which is an error-variant
//! decision (`OverlayError` is public).
//!
//! Refusals are unchanged: the shape's validation (`polygon_part`),
//! `validate_arc_ring` on each ring in list order (`oriented_rings`, the
//! loop `ArcArrangement::new` runs) and linking (`SelfIntersection`).
//!
//! # Cost
//!
//! `O((N + I) log N + C)`, `N` the term edges, `I` the crossings among kept
//! edges, `C` the certification and sample tests. Each is an early-exit
//! box-tree search with an `f64`-filtered exact orientation: short for
//! interior edges, where most candidates contain them, and for boundary
//! edges only the few terms whose box holds them. Where little is
//! certified (a disc round a wall of teeth), the subdivision is no smaller
//! than before and the queries add a logarithmic factor.
//!
//! # Verification
//!
//! The tests below pin the certificate on single edges. The tests in
//! `minkowski` compare this route with the arrangement route, kept there as
//! a test-only oracle (`by_arrangement`), by `assert_eq!` on `polygons()`
//! over walls, holes, shape orientations, discs, furniture scenes and
//! seeded random scenes, and count the pieces cut on a furniture scene.
//! `scripts/probe_term_union_mutants.py` lists the faults they must catch;
//! `cargo bench -p axiolid-benchmark --bench minkowski_furniture` times it
//! and hashes the results.

use axiolid_core::{Point2, Tolerance};
use axiolid_guarantees::Sign;

use crate::arc::{ArcRing, ArcVertex};
use crate::arrangement::{oriented_rings, ArrangementRegion, EdgeUse};
use crate::exact_arc::arrangement::{build_partial, Windings};
use crate::exact_arc::boxes::BoxTree;
use crate::exact_arc::{inside_convex, orient_doubles, Bounds, XPoint};
use crate::minkowski::region_from;
use crate::region::Region;
use crate::OverlayError;

/// Whether the open segment `p -> q` lies strictly inside the convex
/// counter-clockwise ring `ccw`: both ends in the closed ring, and no side
/// line holding both ends. Exact.
fn strictly_inside(ccw: &[Point2], p: Point2, q: Point2) -> bool {
    (0..ccw.len()).all(|i| {
        let (u, v) = (ccw[i], ccw[(i + 1) % ccw.len()]);
        let (a, b) = (orient_doubles(u, v, p), orient_doubles(u, v, q));
        a != Sign::Negative && b != Sign::Negative && (a == Sign::Positive || b == Sign::Positive)
    })
}

/// The faces where `holds` is true, read from the term edges that can bound
/// them. `rings[first..]` are convex pieces, and `holds` must give the same
/// answer for every list holding one of them (they absorb). The result is
/// the one `ArcArrangement::new(rings)` read by `region_of(.., holds, true,
/// ..)` gives, bit for bit.
///
/// # Errors
///
/// As `ArcArrangement::new` and `region_of`.
pub(crate) fn region(
    rings: &[ArcRing],
    first: usize,
    holds: impl Fn(&[usize]) -> bool,
    tolerance: Tolerance,
) -> Result<Region, OverlayError> {
    let (oriented, _) = oriented_rings(rings, tolerance)?;
    let points =
        |ring: &ArcRing| -> Vec<Point2> { ring.vertices.iter().map(|v| v.point).collect() };
    let piece_points: Vec<Vec<Point2>> = oriented[first..].iter().map(points).collect();
    debug_assert!(
        piece_points
            .iter()
            .all(|p| (0..p.len()).all(|i| orient_doubles(
                p[i],
                p[(i + 1) % p.len()],
                p[(i + 2) % p.len()]
            ) != Sign::Negative)),
        "a piece ring is convex and counter-clockwise"
    );
    debug_assert!(
        first == 0 || first >= oriented.len() || holds(&[first]) == holds(&[0, first]),
        "a piece absorbs"
    );
    let piece_exact: Vec<Vec<XPoint>> = piece_points
        .iter()
        .map(|p| p.iter().map(|&q| XPoint::from_f64(q)).collect())
        .collect();
    let pieces = BoxTree::new(
        piece_points
            .iter()
            .map(|p| Bounds::around(p, 0.0))
            .collect(),
    );

    // Certify: keep an edge unless a piece other than its own ring holds
    // its open segment strictly inside.
    let keep: Vec<Vec<bool>> = oriented
        .iter()
        .enumerate()
        .map(|(r, ring)| {
            let n = ring.vertices.len();
            (0..n)
                .map(|k| {
                    let (p, q) = (ring.vertices[k].point, ring.vertices[(k + 1) % n].point);
                    pieces
                        .first(
                            |b| {
                                b.may_hold((p.x, p.x), (p.y, p.y))
                                    && b.may_hold((q.x, q.x), (q.y, q.y))
                            },
                            |j| first + j != r && strictly_inside(&piece_points[j], p, q),
                        )
                        .is_none()
                })
                .collect()
        })
        .collect();

    // Cut the kept edges at every term vertex.
    let splits: Vec<Point2> = oriented
        .iter()
        .flat_map(|ring| ring.vertices.iter().map(|v| v.point))
        .collect();
    let raw = build_partial(&oriented, &keep, &splits)?;

    // Read each piece's sides.
    let mut windings = Windings::new(&oriented[..first]);
    let (mut left, mut right, mut held) = (Vec::new(), Vec::new(), Vec::new());
    let mut keep_piece = Vec::with_capacity(raw.edges.len());
    for (i, edge) in raw.edges.iter().enumerate() {
        left.clear();
        right.clear();
        for &(ring, _, same) in &edge.sources {
            if same {
                left.push(ring);
            } else {
                right.push(ring);
            }
        }
        let has_piece = |side: &[usize]| side.iter().any(|&r| r >= first);
        if !has_piece(&left) || !has_piece(&right) {
            let carried = |r: usize| edge.sources.iter().any(|s| s.0 == r);
            let s = raw.sample(i);
            let sx = s.enclosures();
            let sample_at = pieces.first(
                |b| b.may_hold(sx.0, sx.1),
                |j| !carried(first + j) && inside_convex(s, &piece_exact[j]),
            );
            if let Some(j) = sample_at {
                left.push(first + j);
                right.push(first + j);
            } else {
                windings.holding(s, carried, &mut held);
                left.extend_from_slice(&held);
                right.extend_from_slice(&held);
            }
        }
        left.sort_unstable();
        left.dedup();
        right.sort_unstable();
        right.dedup();
        keep_piece.push(match (holds(&left), holds(&right)) {
            (true, false) => Some(false),
            (false, true) => Some(true),
            _ => None,
        });
    }

    // Link, then round, present and settle as the arrangement route does.
    let positions = raw.vertex_positions();
    let regions: Vec<ArrangementRegion> = raw
        .regions(&keep_piece)?
        .into_iter()
        .map(|(outer, holes)| {
            let uses = |ring: Vec<(usize, bool)>| -> Vec<EdgeUse> {
                ring.into_iter()
                    .map(|(edge, reversed)| EdgeUse { edge, reversed })
                    .collect()
            };
            ArrangementRegion {
                outer: uses(outer),
                holes: holes.into_iter().map(uses).collect(),
            }
        })
        .collect();
    let ring = |uses: &[EdgeUse]| {
        ArcRing::new(
            uses.iter()
                .map(|u| {
                    let edge = &raw.edges[u.edge];
                    let (from, bulge) = if u.reversed {
                        (edge.to, -edge.bulge)
                    } else {
                        (edge.from, edge.bulge)
                    };
                    ArcVertex::bulged(positions[from], bulge)
                })
                .collect(),
        )
    };
    region_from(regions, ring, true, tolerance)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square() -> Vec<Point2> {
        [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
            .iter()
            .map(|&(x, y)| Point2::new(x, y))
            .collect()
    }

    fn inside(p: (f64, f64), q: (f64, f64)) -> bool {
        strictly_inside(&square(), Point2::new(p.0, p.1), Point2::new(q.0, q.1))
    }

    #[test]
    fn an_edge_along_a_side_of_a_piece_is_kept() {
        // On the bottom side, whole or in part, and across a corner.
        assert!(!inside((0.0, 0.0), (4.0, 0.0)));
        assert!(!inside((1.0, 0.0), (3.0, 0.0)));
        assert!(!inside((4.0, 1.0), (4.0, 3.0)));
    }

    #[test]
    fn an_edge_with_an_end_outside_the_piece_is_kept() {
        assert!(!inside((1.0, 1.0), (5.0, 1.0)));
        assert!(!inside((-1.0, 2.0), (2.0, 2.0)));
        // Both ends outside, the segment crossing it.
        assert!(!inside((-1.0, 2.0), (5.0, 2.0)));
        // An end a hair outside.
        assert!(!inside((1.0, 1.0), (4.0 + 1e-15, 2.0)));
    }

    #[test]
    fn an_edge_strictly_inside_a_piece_is_dropped() {
        assert!(inside((1.0, 1.0), (3.0, 3.0)));
        assert!(inside((0.5, 0.5), (3.5, 1.0)));
    }

    #[test]
    fn an_edge_ending_on_a_piece_corner_inside_is_dropped() {
        // An end on a corner or on a side, the rest of the segment inside.
        assert!(inside((0.0, 0.0), (2.0, 2.0)));
        assert!(inside((0.0, 2.0), (2.0, 2.0)));
        assert!(inside((0.0, 2.0), (4.0, 2.0)));
    }

    #[test]
    fn a_ring_does_not_certify_its_own_edges() {
        // A convex ring from `first` on: its own edges lie on its side
        // lines, so none is dropped by it, and a second copy of it, being
        // another ring, holds none of them strictly either.
        let ring = ArcRing::from_points(&square());
        let rings = [ring.clone(), ring];
        let (oriented, _) = oriented_rings(&rings, Tolerance::METRE).unwrap();
        let points: Vec<Vec<Point2>> = oriented
            .iter()
            .map(|r| r.vertices.iter().map(|v| v.point).collect())
            .collect();
        for (r, ring) in points.iter().enumerate() {
            for k in 0..ring.len() {
                let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
                assert!(
                    (0..points.len()).all(|j| !strictly_inside(&points[j], p, q)),
                    "ring {r} edge {k}"
                );
            }
        }
    }
}
