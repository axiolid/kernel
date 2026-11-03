//! Exact validity of arc polygons (#352).
//!
//! Dyadic coordinates and bulges of +-1 and +-0.5 make tangencies and
//! shared stretches exact, which is where a tolerance-based check goes
//! wrong. The differential oracle at the end flattens each ring to chords
//! and asks the straight overlay, away from the chords' own error.

use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::{
    check_arc_polygon, overlay, reverse_arc_ring, validate_arc_ring, ArcPolygon, ArcPolygonCheck,
    ArcPolygonDefect, ArcRing, ArcRingRef, ArcVertex, FillRule, OverlayError, OverlayInput,
    OverlayOperation, Polygon, Ring,
};

fn tol() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).expect("valid tolerance")
}

fn v(x: f64, y: f64, b: f64) -> ArcVertex {
    ArcVertex::bulged(Point2::new(x, y), b)
}

fn ring(vs: &[(f64, f64, f64)]) -> ArcRing {
    ArcRing::new(vs.iter().map(|&(x, y, b)| v(x, y, b)).collect())
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ring(&[(x0, y0, 0.0), (x1, y0, 0.0), (x1, y1, 0.0), (x0, y1, 0.0)])
}

fn poly(outer: ArcRing, holes: Vec<ArcRing>) -> ArcPolygon {
    ArcPolygon { outer, holes }
}

fn check(outer: ArcRing, holes: Vec<ArcRing>) -> ArcPolygonCheck {
    check_arc_polygon(&poly(outer, holes), tol())
}

fn near(a: Point2, x: f64, y: f64) -> bool {
    (a.x - x).abs() < 1e-12 && (a.y - y).abs() < 1e-12
}

/// The one self-contact of a check, as (ring, edges, at).
fn only_self_contact(check: &ArcPolygonCheck) -> (ArcRingRef, [usize; 2], Point2) {
    match check.defects.as_slice() {
        [ArcPolygonDefect::SelfContact { ring, edges, at }] => (*ring, *edges, *at),
        other => panic!("expected one self-contact, got {other:?}"),
    }
}

// ---------------------------------------------------------------- valid

#[test]
fn simple_shapes_are_strictly_valid() {
    let cases: Vec<(&str, ArcRing)> = vec![
        ("square", rect(0.0, 0.0, 2.0, 2.0)),
        ("circle", ArcRing::circle(Point2::new(0.0, 0.0), 1.0)),
        ("d shape", ring(&[(-1.0, 0.0, 0.0), (1.0, 0.0, 1.0)])),
        (
            "pillow",
            ring(&[
                (0.0, 0.0, 0.25),
                (2.0, 0.0, 0.25),
                (2.0, 2.0, 0.25),
                (0.0, 2.0, 0.25),
            ]),
        ),
        (
            "minor and major arc",
            ring(&[(0.0, 0.0, 2.0), (2.0, 0.0, 0.0), (1.0, 3.0, 0.125)]),
        ),
    ];
    for (name, outer) in cases {
        let result = check(outer, vec![]);
        assert!(result.is_strict(), "{name}: {result:?}");
    }
}

#[test]
fn nested_and_disjoint_holes_are_strictly_valid() {
    let annulus = check(
        ArcRing::circle(Point2::new(0.0, 0.0), 2.0),
        vec![ArcRing::circle(Point2::new(0.0, 0.0), 1.0)],
    );
    assert!(annulus.is_strict(), "{annulus:?}");
    let three = check(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![
            rect(1.0, 1.0, 2.0, 2.0),
            rect(4.0, 4.0, 5.0, 5.0),
            rect(7.0, 7.0, 8.0, 8.0),
        ],
    );
    assert!(three.is_strict(), "{three:?}");
}

// ---------------------------------------------------------- self-contact

#[test]
fn bow_tie_crosses_itself() {
    let bow = ring(&[
        (0.0, 0.0, 0.0),
        (2.0, 2.0, 0.0),
        (2.0, 0.0, 0.0),
        (0.0, 3.0, 0.0),
    ]);
    let (ring, edges, at) = only_self_contact(&check(bow, vec![]));
    assert_eq!((ring, edges), (ArcRingRef::Outer, [0, 2]));
    assert!(near(at, 1.2, 1.2), "{at:?}");
}

#[test]
fn a_self_contact_skips_the_nesting_check() {
    // The hole lies outside the outer, which nesting would report.
    let bow = ring(&[
        (0.0, 0.0, 0.0),
        (2.0, 2.0, 0.0),
        (2.0, 0.0, 0.0),
        (0.0, 3.0, 0.0),
    ]);
    let result = check(bow, vec![rect(20.0, 20.0, 22.0, 22.0)]);
    let (ring, edges, _) = only_self_contact(&result);
    assert_eq!((ring, edges), (ArcRingRef::Outer, [0, 2]));
    assert!(result.contacts.is_empty());
}

#[test]
fn symmetric_bow_tie_reports_zero_area_and_the_crossing() {
    let bow = ring(&[
        (0.0, 0.0, 0.0),
        (2.0, 2.0, 0.0),
        (2.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
    ]);
    let result = check(bow, vec![rect(20.0, 20.0, 22.0, 22.0)]);
    assert!(
        matches!(
            result.defects.as_slice(),
            [
                ArcPolygonDefect::Ring {
                    ring: ArcRingRef::Outer,
                    error: OverlayError::ZeroArea
                },
                ArcPolygonDefect::SelfContact { edges: [0, 2], .. }
            ]
        ),
        "{result:?}"
    );
}

#[test]
fn arc_bulging_across_the_opposite_edge() {
    let outer = ring(&[
        (0.0, 0.0, -1.0),
        (2.0, 0.0, 0.0),
        (2.0, 0.5, 0.0),
        (0.0, 0.5, 0.0),
    ]);
    let (_, edges, at) = only_self_contact(&check(outer, vec![]));
    assert_eq!(edges, [0, 2]);
    assert!((at.y - 0.5).abs() < 1e-12);
    assert!(
        ((at.x - 1.0).abs() - 0.75_f64.sqrt()).abs() < 1e-12,
        "{at:?}"
    );
}

#[test]
fn consecutive_edges_meeting_again() {
    let outer = ring(&[(0.0, 0.0, -1.0), (2.0, 0.0, 0.0), (0.0, 1.0, 0.0)]);
    let (_, edges, at) = only_self_contact(&check(outer, vec![]));
    assert_eq!(edges, [0, 1]);
    assert!(near(at, 0.4, 0.8), "{at:?}");
}

#[test]
fn fold_back_spikes() {
    let straight = ring(&[
        (0.0, 0.0, 0.0),
        (4.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (2.0, 2.0, 0.0),
    ]);
    let (_, edges, at) = only_self_contact(&check(straight, vec![]));
    assert_eq!(edges, [0, 1]);
    assert!(near(at, 2.0, 0.0), "{at:?}");
    // The second arc runs back over the first on the same circle.
    let arcs = ring(&[(-4.0, 3.0, -0.5), (4.0, 3.0, 1.0), (-4.0, -3.0, 0.0)]);
    let result = check(arcs, vec![]);
    assert!(
        matches!(
            result.defects.first(),
            Some(ArcPolygonDefect::SelfContact { edges: [0, 1], .. })
        ),
        "{result:?}"
    );
}

#[test]
fn vertex_pinched_on_a_non_adjacent_edge() {
    let outer = ring(&[
        (0.0, 0.0, 0.0),
        (4.0, 0.0, 0.0),
        (4.0, 4.0, 0.0),
        (2.0, 0.0, 0.0),
        (0.0, 4.0, 0.0),
    ]);
    let result = check(outer, vec![]);
    assert!(matches!(
        result.defects.first(),
        Some(ArcPolygonDefect::SelfContact { edges: [0, 2], at, .. }) if near(*at, 2.0, 0.0)
    ));
}

fn tangent_ring(height: f64) -> ArcRing {
    // The arc from (1, 0) to (-1, 0) peaks at (0, 1); edge 2 runs on y = height.
    ring(&[
        (1.0, 0.0, 1.0),
        (-1.0, 0.0, 0.0),
        (-3.0, height, 0.0),
        (3.0, height, 0.0),
    ])
}

#[test]
fn tangency_is_a_contact_and_a_near_miss_is_not() {
    let (_, edges, at) = only_self_contact(&check(tangent_ring(1.0), vec![]));
    assert_eq!(edges, [0, 2]);
    assert!(near(at, 0.0, 1.0), "{at:?}");
    // 2^-40 is far below any tolerance and still separates the two.
    let miss = check(tangent_ring(1.0 + 2f64.powi(-40)), vec![]);
    assert!(miss.is_strict(), "{miss:?}");
}

#[test]
fn closing_edge_pairs() {
    // The closing edge shares only vertex 0 with edge 0.
    assert!(check(rect(0.0, 0.0, 2.0, 2.0), vec![]).is_strict());
    // An arc closing edge crossing edge 1.
    let narrow = ring(&[
        (0.0, 0.0, 0.0),
        (0.25, 0.0, 0.0),
        (0.25, 1.0, 0.0),
        (0.0, 1.0, -1.0),
    ]);
    let (_, edges, _) = only_self_contact(&check(narrow, vec![]));
    assert_eq!(edges, [1, 3]);
    // An arc closing edge crossing edge 0.
    let tri = ring(&[(0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (2.0, 2.0, -1.0)]);
    let (_, edges, at) = only_self_contact(&check(tri, vec![]));
    assert_eq!(edges, [0, 2]);
    assert!(near(at, 2.0, 0.0), "{at:?}");
}

// --------------------------------------------------------------- nesting

#[test]
fn hole_outside_the_outer() {
    let result = check(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(20.0, 20.0, 22.0, 22.0)],
    );
    assert_eq!(
        result.defects,
        vec![ArcPolygonDefect::HoleOutsideOuter { hole: 0 }]
    );
    assert!(result.contacts.is_empty());
}

#[test]
fn hole_crossing_the_outer() {
    let result = check(rect(0.0, 0.0, 10.0, 10.0), vec![rect(8.0, 8.0, 12.0, 12.0)]);
    assert_eq!(
        result.defects,
        vec![ArcPolygonDefect::HoleOutsideOuter { hole: 0 }]
    );
    assert_eq!(result.contacts.len(), 2, "{result:?}");
}

#[test]
fn hole_sharing_a_stretch_with_the_outer() {
    let result = check(rect(0.0, 0.0, 10.0, 10.0), vec![rect(2.0, 0.0, 5.0, 5.0)]);
    assert_eq!(
        result.defects,
        vec![ArcPolygonDefect::BoundaryOverlap {
            rings: [ArcRingRef::Outer, ArcRingRef::Hole(0)]
        }]
    );
}

#[test]
fn holes_that_cross_or_nest() {
    for second in [rect(4.0, 4.0, 8.0, 8.0), rect(3.0, 3.0, 5.0, 5.0)] {
        let result = check(
            rect(0.0, 0.0, 10.0, 10.0),
            vec![rect(2.0, 2.0, 6.0, 6.0), second],
        );
        assert_eq!(
            result.defects,
            vec![ArcPolygonDefect::HolesOverlap { holes: [0, 1] }]
        );
    }
}

#[test]
fn contacts_are_valid_but_not_strict() {
    let at_vertex = check(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![ring(&[(8.0, 4.0, 0.0), (10.0, 5.0, 0.0), (8.0, 6.0, 0.0)])],
    );
    assert!(
        at_vertex.is_valid() && !at_vertex.is_strict(),
        "{at_vertex:?}"
    );
    assert_eq!(at_vertex.contacts.len(), 1);
    assert!(near(at_vertex.contacts[0].at, 10.0, 5.0));
    assert_eq!(
        at_vertex.contacts[0].rings,
        [ArcRingRef::Outer, ArcRingRef::Hole(0)]
    );

    let tangent = check(
        ArcRing::circle(Point2::new(0.0, 0.0), 2.0),
        vec![ArcRing::circle(Point2::new(1.0, 0.0), 1.0)],
    );
    assert!(tangent.is_valid(), "{tangent:?}");
    assert_eq!(tangent.contacts.len(), 1);
    assert!(near(tangent.contacts[0].at, 2.0, 0.0));

    let holes = check(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(2.0, 2.0, 4.0, 4.0), rect(4.0, 4.0, 6.0, 6.0)],
    );
    assert!(holes.is_valid() && !holes.is_strict(), "{holes:?}");
    assert_eq!(
        holes.contacts[0].rings,
        [ArcRingRef::Hole(0), ArcRingRef::Hole(1)]
    );
}

// ------------------------------------------------------ degenerate rings

#[test]
fn degenerate_rings_are_reported_and_skip_nesting() {
    let cases = [
        (
            ring(&[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 1.0, 0.0)]),
            OverlayError::RepeatedVertex,
        ),
        (
            ring(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
            OverlayError::ZeroArea,
        ),
        (
            ring(&[(0.0, 0.0, 1.0), (1.5e-9, 0.0, 1.0)]),
            OverlayError::ZeroRadiusArc,
        ),
    ];
    for (bad, expected) in cases {
        assert_eq!(validate_arc_ring(&bad, tol()), Err(expected.clone()));
        // The hole lies outside the outer: nesting would report it.
        let result = check(bad, vec![rect(20.0, 20.0, 22.0, 22.0)]);
        assert_eq!(
            result.defects.first(),
            Some(&ArcPolygonDefect::Ring {
                ring: ArcRingRef::Outer,
                error: expected
            })
        );
        assert!(!result
            .defects
            .iter()
            .any(|d| matches!(d, ArcPolygonDefect::HoleOutsideOuter { .. })));
    }
}

// ------------------------------------------------------------- invariance

fn remap(report: &ArcPolygonCheck, n: usize) -> Vec<(ArcRingRef, [usize; 2])> {
    report
        .defects
        .iter()
        .map(|d| match d {
            ArcPolygonDefect::SelfContact { ring, edges, .. } => {
                let mut e = edges.map(|k| (2 * n - 2 - k) % n);
                e.sort_unstable();
                (*ring, e)
            }
            other => panic!("unexpected {other:?}"),
        })
        .collect()
}

fn fixtures() -> Vec<ArcRing> {
    vec![
        ring(&[
            (0.0, 0.0, 0.0),
            (2.0, 2.0, 0.0),
            (2.0, 0.0, 0.0),
            (0.0, 3.0, 0.0),
        ]),
        tangent_ring(1.0),
        ring(&[
            (0.0, 0.0, -1.0),
            (2.0, 0.0, 0.0),
            (2.0, 0.5, 0.0),
            (0.0, 0.5, 0.0),
        ]),
    ]
}

#[test]
fn reversing_a_ring_remaps_its_edges() {
    for outer in fixtures() {
        let n = outer.vertices.len();
        let forward = check(outer.clone(), vec![]);
        let backward = check(reverse_arc_ring(&outer), vec![]);
        let expected = remap(&forward, n);
        let got: Vec<_> = backward
            .defects
            .iter()
            .map(|d| match d {
                ArcPolygonDefect::SelfContact { ring, edges, .. } => (*ring, *edges),
                other => panic!("unexpected {other:?}"),
            })
            .collect();
        assert_eq!(got, expected);
    }
}

#[test]
fn reversing_a_hole_keeps_the_report() {
    let hole = rect(8.0, 8.0, 12.0, 12.0);
    let outer = rect(0.0, 0.0, 10.0, 10.0);
    let a = check(outer.clone(), vec![hole.clone()]);
    let b = check(outer, vec![reverse_arc_ring(&hole)]);
    assert_eq!(a, b);
}

#[test]
fn permuting_holes_permutes_indices() {
    let (a, b, c) = (
        rect(1.0, 1.0, 5.0, 5.0),
        rect(7.0, 7.0, 8.0, 8.0),
        rect(3.0, 3.0, 6.0, 6.0),
    );
    let outer = rect(0.0, 0.0, 10.0, 10.0);
    let first = check(outer.clone(), vec![a.clone(), b.clone(), c.clone()]);
    assert_eq!(
        first.defects,
        vec![ArcPolygonDefect::HolesOverlap { holes: [0, 2] }]
    );
    let second = check(outer, vec![c, a, b]);
    assert_eq!(
        second.defects,
        vec![ArcPolygonDefect::HolesOverlap { holes: [0, 1] }]
    );
}

#[test]
fn reports_repeat_and_survive_translation() {
    const SHIFT: f64 = 1_048_576.0;
    let shifted = |r: &ArcRing| {
        ArcRing::new(
            r.vertices
                .iter()
                .map(|p| v(p.point.x + SHIFT, p.point.y + SHIFT, p.bulge))
                .collect(),
        )
    };
    for outer in fixtures().into_iter().chain([rect(0.0, 0.0, 2.0, 2.0)]) {
        let base = check(outer.clone(), vec![]);
        assert_eq!(base, check(outer.clone(), vec![]));
        let moved = check(shifted(&outer), vec![]);
        assert_eq!(
            base.defects.len(),
            moved.defects.len(),
            "{base:?} {moved:?}"
        );
        for (x, y) in base.defects.iter().zip(&moved.defects) {
            match (x, y) {
                (
                    ArcPolygonDefect::SelfContact {
                        edges: e1, at: a1, ..
                    },
                    ArcPolygonDefect::SelfContact {
                        edges: e2, at: a2, ..
                    },
                ) => {
                    assert_eq!(e1, e2);
                    assert!(
                        (a1.x + SHIFT - a2.x).abs() < 1e-9 && (a1.y + SHIFT - a2.y).abs() < 1e-9
                    );
                }
                other => panic!("unexpected {other:?}"),
            }
        }
    }
}

// ------------------------------------------------------------ the oracle

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn random_ring(rng: &mut Rng, span: u64) -> ArcRing {
    const BULGES: [f64; 5] = [0.0, 0.5, -0.5, 1.0, -1.0];
    let count = 3 + rng.below(3) as usize;
    let mut vertices: Vec<ArcVertex> = Vec::new();
    while vertices.len() < count {
        let x = (rng.below(span) as f64 - (span / 2) as f64) * 0.25;
        let y = (rng.below(span) as f64 - (span / 2) as f64) * 0.25;
        let bulge = BULGES[rng.below(5) as usize];
        vertices.push(v(x, y, bulge));
    }
    ArcRing::new(vertices)
}

const CHORDS: usize = 64;

/// The ring flattened to chords, with the largest sagitta used.
fn flatten(ring: &ArcRing) -> (Vec<Point2>, f64) {
    let n = ring.vertices.len();
    let mut points = Vec::new();
    let mut sagitta = 0.0f64;
    for i in 0..n {
        let (a, b) = (ring.vertices[i], ring.vertices[(i + 1) % n]);
        points.push(a.point);
        if a.bulge == 0.0 {
            continue;
        }
        let (dx, dy) = (b.point.x - a.point.x, b.point.y - a.point.y);
        let c = dx.hypot(dy);
        if c == 0.0 {
            continue;
        }
        let k = a.bulge;
        let lever = c * (1.0 - k * k) / (4.0 * k);
        let centre = (
            0.5 * (a.point.x + b.point.x) - dy / c * lever,
            0.5 * (a.point.y + b.point.y) + dx / c * lever,
        );
        let radius = c * (1.0 + k * k) / (4.0 * k.abs());
        let theta = 4.0 * k.atan();
        let start = (a.point.y - centre.1).atan2(a.point.x - centre.0);
        for s in 1..CHORDS {
            let angle = start + theta * s as f64 / CHORDS as f64;
            points.push(Point2::new(
                centre.0 + radius * angle.cos(),
                centre.1 + radius * angle.sin(),
            ));
        }
        sagitta = sagitta.max(radius * (1.0 - (theta.abs() / (2.0 * CHORDS as f64)).cos()));
    }
    (points, sagitta)
}

fn dist_point_segment(p: Point2, a: Point2, b: Point2) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    let t = if len2 == 0.0 {
        0.0
    } else {
        (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0)
    };
    (p.x - (a.x + t * dx)).hypot(p.y - (a.y + t * dy))
}

fn cross(o: Point2, a: Point2, b: Point2) -> f64 {
    (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
}

/// Whether two chord edges leave the true arcs' relation undecided: they
/// touch or pass within `margin`, or cross too shallowly or too near an
/// end for the arcs to be sure to cross as well.
fn undecided(a: Point2, b: Point2, c: Point2, d: Point2, margin: f64) -> bool {
    let (d1, d2, d3, d4) = (
        cross(a, b, c),
        cross(a, b, d),
        cross(c, d, a),
        cross(c, d, b),
    );
    if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
        let (u, w) = ((b.x - a.x, b.y - a.y), (d.x - c.x, d.y - c.y));
        let sin = (u.0 * w.1 - u.1 * w.0).abs() / (u.0.hypot(u.1) * w.0.hypot(w.1));
        let t = d3 / (d3 - d4);
        let at = Point2::new(c.x + t * w.0, c.y + t * w.1);
        let near_end = [a, b, c, d]
            .iter()
            .any(|p| (p.x - at.x).hypot(p.y - at.y) < margin / sin.max(1e-9));
        return sin < 0.2 || near_end;
    }
    dist_point_segment(a, c, d)
        .min(dist_point_segment(b, c, d))
        .min(dist_point_segment(c, a, b))
        .min(dist_point_segment(d, a, b))
        <= margin
}

/// Whether any two chord edges that are not neighbours are undecided.
fn any_undecided(rings: &[Vec<Point2>], margin: f64) -> bool {
    let edges: Vec<(usize, usize, Point2, Point2)> = rings
        .iter()
        .enumerate()
        .flat_map(|(r, pts)| (0..pts.len()).map(move |i| (r, i, pts[i], pts[(i + 1) % pts.len()])))
        .collect();
    for (k, &(r1, i1, a, b)) in edges.iter().enumerate() {
        for &(r2, i2, c, d) in &edges[k + 1..] {
            if r1 == r2 {
                let n = rings[r1].len();
                if (i2 - i1).min(n - (i2 - i1)) <= 1 {
                    continue;
                }
            }
            if undecided(a, b, c, d, margin) {
                return true;
            }
        }
    }
    false
}

#[test]
fn agrees_with_the_straight_overlay_on_random_scenes() {
    let mut rng = Rng(0x5eed_0352);
    let (mut compared, mut skipped, mut valid, mut invalid) = (0, 0, 0, 0);
    for _ in 0..1500 {
        let outer = random_ring(&mut rng, 17);
        let holes = if rng.below(2) == 0 {
            vec![]
        } else {
            vec![random_ring(&mut rng, 9)]
        };
        let polygon = poly(outer, holes);
        let flat: Vec<(Vec<Point2>, f64)> = std::iter::once(&polygon.outer)
            .chain(&polygon.holes)
            .map(flatten)
            .collect();
        let sagitta = flat.iter().map(|f| f.1).fold(0.0, f64::max);
        let rings: Vec<Vec<Point2>> = flat.into_iter().map(|f| f.0).collect();
        if any_undecided(&rings, 10.0 * sagitta.max(1e-6)) {
            skipped += 1;
            continue;
        }
        let exact = check_arc_polygon(&polygon, tol()).is_valid();
        // The straight overlay validates every ring and every hole's place.
        let subject = OverlayInput {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Point2::X,
                y: Point2::Y,
            },
            polygons: vec![Polygon {
                outer: Ring {
                    points: rings[0].clone(),
                },
                holes: rings[1..]
                    .iter()
                    .map(|p| Ring { points: p.clone() })
                    .collect(),
            }],
        };
        let empty = OverlayInput {
            polygons: vec![],
            ..subject.clone()
        };
        let chorded = overlay(
            &subject,
            &empty,
            OverlayOperation::Union,
            FillRule::EvenOdd,
            tol(),
        )
        .is_ok();
        assert_eq!(
            exact,
            chorded,
            "{polygon:?} {:?}",
            check_arc_polygon(&polygon, tol())
        );
        compared += 1;
        if exact {
            valid += 1;
        } else {
            invalid += 1;
        }
    }
    eprintln!("compared {compared}, skipped {skipped}, valid {valid}, invalid {invalid}");
    assert!(compared >= 100 && valid >= 10 && invalid >= 10);
}
