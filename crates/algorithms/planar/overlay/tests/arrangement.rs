//! The multi-ring subdivision (`ArcArrangement`, #120).
//!
//! Oracle: for two rings, selecting "in A and in B", "in A or in B", "in A
//! and not in B" must give the same area as the two-operand `arc_overlay`,
//! which is tested against closed forms elsewhere. For three rings, the
//! inclusion-exclusion identity ties every selection to pairwise results.
//! Structural checks: every region boundary closes on shared vertex
//! indices, sources name real input edges, and each piece bounds a
//! selection exactly when the selection differs across it.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{
    arc_overlay, arc_ring_area, ArcArrangement, ArcRing, ArcVertex, OverlayOperation,
};

fn disc(cx: f64, cy: f64, r: f64) -> ArcRing {
    ArcRing::circle(Point2::new(cx, cy), r)
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing::new(vec![
        ArcVertex::straight(Point2::new(x0, y0)),
        ArcVertex::straight(Point2::new(x1, y0)),
        ArcVertex::straight(Point2::new(x1, y1)),
        ArcVertex::straight(Point2::new(x0, y1)),
    ])
}

fn tol() -> Tolerance {
    Tolerance::METRE
}

/// Net area of the regions where `inside` holds.
fn area(arr: &ArcArrangement, inside: impl Fn(&[bool]) -> bool) -> f64 {
    let regions = arr.regions(inside).expect("links");
    let mut total = 0.0;
    for region in &regions {
        total += arc_ring_area(&arr.ring(&region.outer));
        for hole in &region.holes {
            total += arc_ring_area(&arr.ring(hole));
        }
    }
    total
}

fn overlay_area(a: &ArcRing, b: &ArcRing, op: OverlayOperation) -> f64 {
    let result = arc_overlay(a, b, op, tol()).expect("overlay");
    result
        .regions
        .iter()
        .map(|r| arc_ring_area(&r.outer) + r.holes.iter().map(arc_ring_area).sum::<f64>())
        .sum()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// Every region boundary is a closed walk over shared vertex indices.
fn assert_closed(arr: &ArcArrangement, inside: impl Fn(&[bool]) -> bool) {
    for region in arr.regions(inside).expect("links") {
        for ring in std::iter::once(&region.outer).chain(&region.holes) {
            for (i, u) in ring.iter().enumerate() {
                let next = ring[(i + 1) % ring.len()];
                let e = &arr.edges()[u.edge];
                let n = &arr.edges()[next.edge];
                let end = if u.reversed { e.from } else { e.to };
                let start = if next.reversed { n.to } else { n.from };
                assert_eq!(end, start, "boundary breaks between uses {i} and next");
            }
        }
    }
}

const PAIRS: [(f64, f64, f64, f64); 6] = [
    // disc centre x, y, radius; rectangle half-width
    (0.0, 0.0, 1.0, 0.5),
    (0.3, -0.2, 1.3, 0.9),
    (1.0, 0.0, 1.0, 1.0),
    (2.5, 0.0, 1.0, 0.7),
    (0.0, 0.0, 0.2, 2.0),
    (0.7, 0.7, 0.9, 0.8),
];

#[test]
fn two_rings_agree_with_the_two_operand_overlay() {
    for &(cx, cy, r, h) in &PAIRS {
        let a = disc(cx, cy, r);
        let b = rect(-h, -h, h, h);
        let arr = ArcArrangement::new(&[a.clone(), b.clone()], tol()).expect("valid");
        for (op, select) in [
            (
                OverlayOperation::Intersection,
                (|m: &[bool]| m[0] && m[1]) as fn(&[bool]) -> bool,
            ),
            (OverlayOperation::Union, |m: &[bool]| m[0] || m[1]),
            (OverlayOperation::Difference, |m: &[bool]| m[0] && !m[1]),
        ] {
            let want = if matches!(op, OverlayOperation::Intersection)
                && overlay_area(&a, &b, op) == 0.0
            {
                0.0
            } else {
                overlay_area(&a, &b, op)
            };
            let got = area(&arr, select);
            assert!(
                close(got, want),
                "{op:?} at {cx},{cy},{r},{h}: {got} vs {want}"
            );
            assert_closed(&arr, select);
        }
    }
}

#[test]
fn three_rings_obey_inclusion_exclusion() {
    let rings = [
        disc(0.0, 0.0, 1.0),
        disc(1.0, 0.2, 1.0),
        rect(-0.4, -1.5, 0.6, 0.3),
    ];
    let arr = ArcArrangement::new(&rings, tol()).expect("valid");
    let single = |i: usize| area(&arr, |m| m[i]);
    let pair = |i: usize, j: usize| area(&arr, |m| m[i] && m[j]);
    let triple = area(&arr, |m| m[0] && m[1] && m[2]);
    let union = area(&arr, |m| m.iter().any(|&x| x));
    let expected =
        single(0) + single(1) + single(2) - pair(0, 1) - pair(0, 2) - pair(1, 2) + triple;
    assert!(close(union, expected), "{union} vs {expected}");
    // Singles and pairs are themselves independently checkable.
    for (i, ring) in rings.iter().enumerate() {
        assert!(close(single(i), arc_ring_area(ring).abs()), "ring {i}");
    }
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        let want = overlay_area(&rings[i], &rings[j], OverlayOperation::Intersection);
        assert!(
            close(pair(i, j), want),
            "pair {i},{j}: {} vs {want}",
            pair(i, j)
        );
    }
    assert!(triple > 0.0, "the three rings share a region");
    // "In exactly one" is the symmetric part; it must close too.
    assert_closed(&arr, |m| m.iter().filter(|&&x| x).count() == 1);
}

#[test]
fn a_ring_inside_another_becomes_a_hole_of_the_difference() {
    let rings = [rect(-2.0, -2.0, 2.0, 2.0), disc(0.0, 0.0, 1.0)];
    let arr = ArcArrangement::new(&rings, tol()).expect("valid");
    let regions = arr.regions(|m| m[0] && !m[1]).expect("links");
    assert_eq!(regions.len(), 1);
    assert_eq!(regions[0].holes.len(), 1, "the disc is a hole");
    let got = area(&arr, |m| m[0] && !m[1]);
    assert!(close(got, 16.0 - std::f64::consts::PI), "{got}");
}

#[test]
fn a_shared_boundary_is_one_piece_with_both_sources() {
    // Two squares sharing the edge x = 1, y in [0, 1].
    let rings = [rect(0.0, 0.0, 1.0, 1.0), rect(1.0, 0.0, 2.0, 1.0)];
    let arr = ArcArrangement::new(&rings, tol()).expect("valid");
    let shared: Vec<_> = arr
        .edges()
        .iter()
        .filter(|e| e.sources.len() == 2)
        .collect();
    assert_eq!(shared.len(), 1, "one piece for the shared edge");
    let e = shared[0];
    // Square 0 has its inside on one side, square 1 on the other.
    assert!(e.inside_left(0) != e.inside_left(1));
    assert!(e.inside_left(0) != e.inside_right(0));
    // The union has no piece there: both sides are in it.
    let union = arr.regions(|m| m[0] || m[1]).expect("links");
    assert_eq!(union.len(), 1);
    assert!(union[0]
        .outer
        .iter()
        .all(|u| arr.edges()[u.edge].sources.len() == 1));
    assert!(close(area(&arr, |m| m[0] || m[1]), 2.0));
}

#[test]
fn sources_name_the_callers_edges_even_for_clockwise_input() {
    let ccw = rect(0.0, 0.0, 2.0, 1.0);
    let mut cw = ccw.clone();
    cw.vertices.reverse();
    for (ring, label) in [(ccw, "ccw"), (cw, "cw")] {
        let arr = ArcArrangement::new(std::slice::from_ref(&ring), tol()).expect("valid");
        assert_eq!(arr.edges().len(), 4);
        for edge in arr.edges() {
            let source = edge.sources[0];
            let n = ring.vertices.len();
            let a = ring.vertices[source.edge].point;
            let b = ring.vertices[(source.edge + 1) % n].point;
            let (from, to) = (arr.vertices()[edge.from], arr.vertices()[edge.to]);
            let (want_from, want_to) = if source.forward { (a, b) } else { (b, a) };
            assert_eq!(
                (from, to),
                (want_from, want_to),
                "{label}: edge {}",
                source.edge
            );
        }
    }
}

#[test]
fn crossings_are_shared_vertices_not_near_duplicates() {
    // A disc through two corners of a square with inexact coordinates.
    let rings = [
        rect(0.0, 0.0, 0.3, 0.3),
        disc(0.15, 0.15, 0.15 * 2f64.sqrt()),
    ];
    let arr = ArcArrangement::new(&rings, tol()).expect("valid");
    let v = arr.vertices();
    for i in 0..v.len() {
        for j in i + 1..v.len() {
            assert!(
                (v[i] - v[j]).length() > 1e-12,
                "vertices {i} and {j} coincide after rounding: {:?}",
                v[i]
            );
        }
    }
    assert_closed(&arr, |m| m[0] ^ m[1]);
}

#[test]
fn a_malformed_ring_is_refused_by_its_reason() {
    let bad = ArcRing::new(vec![
        ArcVertex::straight(Point2::new(0.0, 0.0)),
        ArcVertex::straight(Point2::new(1.0, 0.0)),
    ]);
    assert!(ArcArrangement::new(&[rect(0.0, 0.0, 1.0, 1.0), bad], tol()).is_err());
}

#[test]
fn a_shared_arc_is_one_piece_with_both_sources() {
    // The upper half disc shares the disc's upper arc, the same way: one
    // piece carried by both rings, with both insides on its left.
    let half = ArcRing::new(vec![
        ArcVertex::bulged(Point2::new(1.0, 0.0), 1.0),
        ArcVertex::straight(Point2::new(-1.0, 0.0)),
    ]);
    let arr = ArcArrangement::new(&[disc(0.0, 0.0, 1.0), half], tol()).expect("valid");
    let shared: Vec<_> = arr
        .edges()
        .iter()
        .filter(|e| e.sources.len() == 2)
        .collect();
    assert_eq!(shared.len(), 1, "one piece for the shared arc");
    assert!(shared[0].inside_left(0) == shared[0].inside_left(1));
    assert_eq!(arr.edges().len(), 3, "lower arc, shared arc, diameter");
}

#[test]
fn a_segment_touching_an_arc_stays_inside_its_ring() {
    // The bottom of `dented` is an arc bowed up into the ring, its top at
    // (0, 1); the strip's bottom edge touches it there from inside the
    // ring. The touch is a stop where one other boundary meets the edge,
    // but it crosses nothing: every piece of the strip lies inside
    // `dented` on both sides (#218).
    let dented = ArcRing::new(vec![
        ArcVertex::bulged(Point2::new(-2.0, 0.0), -0.5),
        ArcVertex::straight(Point2::new(2.0, 0.0)),
        ArcVertex::straight(Point2::new(2.0, 4.0)),
        ArcVertex::straight(Point2::new(-2.0, 4.0)),
    ]);
    let strip = rect(-1.0, 1.0, 1.0, 1.5);
    let arr = ArcArrangement::new(&[dented, strip], tol()).expect("valid rings");
    let mut pieces = 0;
    for edge in arr.edges() {
        if edge.sources.iter().all(|s| s.ring == 1) {
            pieces += 1;
            assert!(
                edge.inside_left(0) && edge.inside_right(0),
                "piece {edge:?} left the ring it lies in"
            );
        }
    }
    // The touch splits the strip's bottom edge in two.
    assert_eq!(pieces, 5);
}

/// A ring's vertex inside another ring's edge (#338): the edge is split
/// there, whichever ring comes first, whether the rings share a stretch of
/// the edge or touch at the vertex only, and whether the edge is vertical
/// (every point on it has its ends' `x`) or slanted. Nothing snaps here,
/// so the subdivision itself must record the meeting for the edge it lies
/// inside, though it lies at an end of the other edge.
#[test]
fn a_vertex_inside_another_edge_splits_it() {
    let straight = |p: &[(f64, f64)]| {
        ArcRing::new(
            p.iter()
                .map(|&(x, y)| ArcVertex::straight(Point2::new(x, y)))
                .collect(),
        )
    };
    let square = rect(0.0, 0.0, 2.0, 2.0);
    let triangle = straight(&[(0.0, 2.0), (2.0, 0.0), (4.0, 4.0)]);
    let scenes = [
        // A stretch of the square's right edge shared.
        (&square, rect(2.0, 0.5, 3.0, 1.5)),
        // One vertex on the square's right edge, from outside.
        (
            &square,
            straight(&[(2.0, 1.0), (3.0, 0.5), (3.5, 1.0), (3.0, 1.5)]),
        ),
        // A stretch of the triangle's slanted edge shared, from outside.
        (
            &triangle,
            straight(&[(0.5, 1.5), (0.0, 1.0), (1.0, 0.0), (1.5, 0.5)]),
        ),
        // The same stretches with a vertex where the boundary goes
        // straight on: both its edges run along the other edge, so only
        // the overlap reports it there, as an end of the other edge.
        (
            &square,
            straight(&[(2.0, 0.5), (3.0, 0.5), (3.0, 1.5), (2.0, 1.5), (2.0, 1.0)]),
        ),
        (
            &triangle,
            straight(&[(0.5, 1.5), (0.0, 1.0), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0)]),
        ),
    ];
    for (k, (first, other)) in scenes.iter().enumerate() {
        for rings in [
            vec![(*first).clone(), other.clone()],
            vec![other.clone(), (*first).clone()],
        ] {
            let arr = ArcArrangement::new(&rings, tol()).expect("subdivides");
            let either = area(&arr, |s| s[0] || s[1]);
            let expected = overlay_area(&rings[0], &rings[1], OverlayOperation::Union);
            assert!(
                close(either, expected),
                "scene {k}: {either} against {expected}"
            );
            let both = area(&arr, |s| s[0] && s[1]);
            assert!(close(both, 0.0), "scene {k}: {both} in both");
            assert_closed(&arr, |s| s[0] || s[1]);
            assert_closed(&arr, |s| s[0] != s[1]);
            for ring in &rings {
                for v in &ring.vertices {
                    assert!(
                        arr.vertices().contains(&v.point),
                        "scene {k}: {:?} not a vertex",
                        v.point
                    );
                }
            }
            // Every piece lies on its sources: none spans a vertex of
            // another ring.
            for edge in arr.edges() {
                let (a, b) = (arr.vertices()[edge.from], arr.vertices()[edge.to]);
                for ring in &rings {
                    for v in &ring.vertices {
                        let p = v.point;
                        let inside = (p - a).perp_dot(b - a) == 0.0
                            && (p - a).dot(b - a) > 0.0
                            && (p - b).dot(a - b) > 0.0;
                        assert!(!inside, "scene {k}: {p:?} inside a piece {a:?}-{b:?}");
                    }
                }
            }
        }
    }
}
