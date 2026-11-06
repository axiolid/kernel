//! An independent topology checker for visible-face outcomes.
//!
//! It shares no code with the provider and decides every sign with exact
//! `Dyadic` arithmetic on the binary64 ring points. Rules (ADR 0100):
//!
//! - every ring has at least three points and no two consecutive points are
//!   equal, the wrap included;
//! - no two edges of any rings cross properly;
//! - no ring vertex lies in the relative interior of any edge of any ring,
//!   which also forbids two edges that overlap other than identically;
//! - an outer ring has positive exact signed area and a hole negative.
//!
//! Rings may touch at a vertex, and one ring may revisit a vertex.

use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic};
use axiolid_mesh_hidden_lines_contract::VisibleFaceOutcome;

fn exact(value: f64) -> Dyadic {
    Dyadic::try_from_f64(value).expect("a ring point is finite")
}

/// The sign of the turn `a`, `b`, `c`: `1` counter-clockwise, `-1`
/// clockwise, `0` collinear, in exact arithmetic.
fn turn(a: Point2, b: Point2, c: Point2) -> i8 {
    let (ax, ay) = (exact(a.x), exact(a.y));
    let left = exact(b.x).sub(&ax).mul(&exact(c.y).sub(&ay));
    let right = exact(b.y).sub(&ay).mul(&exact(c.x).sub(&ax));
    match left.sub(&right).sign() {
        Some(axiolid_guarantees::Sign::Positive) => 1,
        Some(axiolid_guarantees::Sign::Negative) => -1,
        _ => 0,
    }
}

/// Twice the exact signed area of a ring: positive counter-clockwise.
fn signed_area_sign(ring: &[Point2]) -> i8 {
    let mut twice = Dyadic::zero();
    for (i, a) in ring.iter().enumerate() {
        let b = ring[(i + 1) % ring.len()];
        twice = twice.add(
            &exact(a.x)
                .mul(&exact(b.y))
                .sub(&exact(b.x).mul(&exact(a.y))),
        );
    }
    match twice.sign() {
        Some(axiolid_guarantees::Sign::Positive) => 1,
        Some(axiolid_guarantees::Sign::Negative) => -1,
        _ => 0,
    }
}

/// Whether `p`, collinear with the edge, lies strictly between its ends.
fn strictly_inside(edge: [Point2; 2], p: Point2) -> bool {
    let [a, b] = edge;
    if p == a || p == b || turn(a, b, p) != 0 {
        return false;
    }
    let (lo_x, hi_x) = (a.x.min(b.x), a.x.max(b.x));
    let (lo_y, hi_y) = (a.y.min(b.y), a.y.max(b.y));
    (lo_x..=hi_x).contains(&p.x) && (lo_y..=hi_y).contains(&p.y)
}

/// Why two edges break the rules, if they do.
pub fn edge_pair_fault(first: [Point2; 2], second: [Point2; 2]) -> Option<&'static str> {
    let [p, q] = first;
    let [r, s] = second;
    let (o1, o2) = (turn(p, q, r), turn(p, q, s));
    let (o3, o4) = (turn(r, s, p), turn(r, s, q));
    if o1 * o2 < 0 && o3 * o4 < 0 {
        return Some("two edges cross");
    }
    let touches = [p, q].iter().any(|&x| strictly_inside(second, x))
        || [r, s].iter().any(|&x| strictly_inside(first, x));
    touches.then_some("a vertex lies inside an edge")
}

fn edges_of(ring: &[Point2]) -> impl Iterator<Item = [Point2; 2]> + '_ {
    (0..ring.len()).map(|i| [ring[i], ring[(i + 1) % ring.len()]])
}

/// Panic unless every ring of the outcome follows the rules above.
pub fn assert_valid_topology(outcome: &VisibleFaceOutcome) {
    let mut edges: Vec<[Point2; 2]> = Vec::new();
    for region in &outcome.regions {
        let rings = std::iter::once((&region.outer, true))
            .chain(region.holes.iter().map(|hole| (hole, false)));
        for (ring, outer) in rings {
            assert!(
                ring.len() >= 3,
                "{region:?}: a ring has {} points",
                ring.len()
            );
            for pair in edges_of(ring) {
                assert!(pair[0] != pair[1], "{region:?}: a point repeats");
            }
            let sign = signed_area_sign(ring);
            assert!(
                sign == if outer { 1 } else { -1 },
                "{region:?}: a {} ring has area sign {sign}",
                if outer { "outer" } else { "hole" }
            );
            edges.extend(edges_of(ring));
        }
    }
    // Sweep over the edges by their smallest x.
    let key = |edge: &[Point2; 2]| (edge[0].x.min(edge[1].x), edge[0].x.max(edge[1].x));
    edges.sort_by(|a, b| key(a).0.total_cmp(&key(b).0));
    for (i, first) in edges.iter().enumerate() {
        let (_, hi) = key(first);
        let (y_lo, y_hi) = (first[0].y.min(first[1].y), first[0].y.max(first[1].y));
        for second in &edges[i + 1..] {
            if key(second).0 > hi {
                break;
            }
            let (lo, up) = (second[0].y.min(second[1].y), second[0].y.max(second[1].y));
            if up < y_lo || lo > y_hi {
                continue;
            }
            if let Some(fault) = edge_pair_fault(*first, *second) {
                panic!("{fault}: {first:?} and {second:?}");
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    #[test]
    fn crossing_edges_and_vertices_inside_edges_are_faults() {
        let diagonal = [at(0.0, 0.0), at(4.0, 4.0)];
        assert_eq!(
            edge_pair_fault(diagonal, [at(0.0, 4.0), at(4.0, 0.0)]),
            Some("two edges cross")
        );
        assert_eq!(
            edge_pair_fault(diagonal, [at(2.0, 2.0), at(2.0, 5.0)]),
            Some("a vertex lies inside an edge")
        );
        assert_eq!(
            edge_pair_fault(diagonal, [at(2.0, 2.0), at(6.0, 6.0)]),
            Some("a vertex lies inside an edge"),
            "an overlap that is not identical"
        );
    }

    #[test]
    fn identical_edges_and_shared_ends_are_not_faults() {
        let diagonal = [at(0.0, 0.0), at(4.0, 4.0)];
        assert_eq!(
            edge_pair_fault(diagonal, [at(4.0, 4.0), at(0.0, 0.0)]),
            None
        );
        assert_eq!(
            edge_pair_fault(diagonal, [at(4.0, 4.0), at(9.0, 0.0)]),
            None
        );
        assert_eq!(
            edge_pair_fault(diagonal, [at(0.0, 4.0), at(5.0, 5.0)]),
            None
        );
    }

    #[test]
    fn the_sign_of_a_ring_is_exact() {
        let up = f64::from_bits(1.0f64.to_bits() + 1);
        let thin = [at(0.0, 0.0), at(1.0, 0.0), at(0.5, up - 1.0)];
        assert_eq!(signed_area_sign(&thin), 1);
        let mut flipped = thin;
        flipped.reverse();
        assert_eq!(signed_area_sign(&flipped), -1);
        assert_eq!(
            signed_area_sign(&[at(0.0, 0.0), at(1.0, 0.0), at(2.0, 0.0)]),
            0
        );
    }

    #[test]
    fn a_ring_that_revisits_a_vertex_may_touch_itself() {
        // A bow tie joined at (2, 2): two triangles, one outer ring.
        let ring = [
            at(0.0, 0.0),
            at(2.0, 2.0),
            at(4.0, 0.0),
            at(4.0, 4.0),
            at(2.0, 2.0),
            at(0.0, 4.0),
        ];
        let edges: Vec<_> = edges_of(&ring).collect();
        for (i, a) in edges.iter().enumerate() {
            for b in &edges[i + 1..] {
                assert_eq!(edge_pair_fault(*a, *b), None, "{a:?} {b:?}");
            }
        }
    }
}
