//! Opt-in clearance of `check_arc_polygon_with` (#388).
//!
//! Rings that miss each other by less than a distance are valid exactly;
//! with a clearance they are reported, decided by exact signs. The fixtures
//! put lines and arcs a half and twice the tolerance apart, plus a pair at a
//! distance that is a dyadic number, so the verdict at that distance itself
//! shows the comparison is exact and strict.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{
    check_arc_polygon, check_arc_polygon_with, ArcCheckOptions, ArcPolygon, ArcPolygonCheck,
    ArcPolygonDefect, ArcRing, ArcRingRef, ArcVertex,
};

const TOL: f64 = 1e-3;

fn tol() -> Tolerance {
    Tolerance::new(TOL, TOL).expect("valid tolerance")
}

fn ring(vs: &[(f64, f64, f64)]) -> ArcRing {
    ArcRing::new(
        vs.iter()
            .map(|&(x, y, b)| ArcVertex::bulged(Point2::new(x, y), b))
            .collect(),
    )
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ring(&[(x0, y0, 0.0), (x1, y0, 0.0), (x1, y1, 0.0), (x0, y1, 0.0)])
}

fn circle(x: f64, y: f64, r: f64) -> ArcRing {
    ArcRing::circle(Point2::new(x, y), r)
}

fn with_clearance(outer: ArcRing, holes: Vec<ArcRing>, clearance: f64) -> ArcPolygonCheck {
    check_arc_polygon_with(
        &ArcPolygon { outer, holes },
        tol(),
        ArcCheckOptions::new().with_clearance(clearance),
    )
}

/// The clearance defects as ((ring, edge), (ring, edge)), with whether the
/// pair was only undecided.
type Pair = ((ArcRingRef, usize), (ArcRingRef, usize), bool);

fn near(check: &ArcPolygonCheck) -> Vec<Pair> {
    check
        .defects
        .iter()
        .filter_map(|defect| match defect {
            ArcPolygonDefect::Clearance {
                edges, undecided, ..
            } => Some((edges[0], edges[1], *undecided)),
            _ => None,
        })
        .collect()
}

const OUTER: ArcRingRef = ArcRingRef::Outer;
const HOLE: ArcRingRef = ArcRingRef::Hole(0);

#[test]
fn clearance_is_off_by_default() {
    let polygon = ArcPolygon {
        outer: rect(0.0, 0.0, 10.0, 10.0),
        holes: vec![rect(0.0005, 1.0, 5.0, 4.0)],
    };
    let plain = check_arc_polygon(&polygon, tol());
    assert!(plain.is_valid(), "{plain:?}");
    assert_eq!(
        check_arc_polygon_with(&polygon, tol(), ArcCheckOptions::new()),
        plain
    );
    assert_eq!(ArcCheckOptions::default(), ArcCheckOptions::new());
}

#[test]
fn a_hole_near_the_outer_ring_is_reported_at_half_and_not_at_twice_the_tolerance() {
    let close = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(0.5 * TOL, 1.0, 5.0, 4.0)],
        TOL,
    );
    // The hole's left wall and the two walls ending on it are each within
    // the clearance of the outer left side (edge 3).
    assert_eq!(
        near(&close),
        vec![
            ((OUTER, 3), (HOLE, 0), false),
            ((OUTER, 3), (HOLE, 2), false),
            ((OUTER, 3), (HOLE, 3), false),
        ]
    );
    assert!(!close.is_valid());

    let apart = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(2.0 * TOL, 1.0, 5.0, 4.0)],
        TOL,
    );
    assert!(apart.is_valid(), "{apart:?}");
}

#[test]
fn a_slot_in_one_ring_is_reported() {
    let slot = |width: f64| {
        ring(&[
            (0.0, 0.0, 0.0),
            (10.0, 0.0, 0.0),
            (10.0, 10.0, 0.0),
            (5.0 + width, 10.0, 0.0),
            (5.0 + width, 2.0, 0.0),
            (5.0, 2.0, 0.0),
            (5.0, 10.0, 0.0),
            (0.0, 10.0, 0.0),
        ])
    };
    // A slot narrower than the ring tolerance would be a repeated vertex, so
    // the ring tolerance is far below the clearance here.
    let fine = Tolerance::new(1e-9, 1e-9).expect("valid tolerance");
    let check = |width: f64| {
        check_arc_polygon_with(
            &ArcPolygon {
                outer: slot(width),
                holes: vec![],
            },
            fine,
            ArcCheckOptions::new().with_clearance(TOL),
        )
    };
    let narrow = check(0.5 * TOL);
    assert!(
        near(&narrow).contains(&((OUTER, 3), (OUTER, 5), false)),
        "{narrow:?}"
    );
    let wide = check(2.0 * TOL);
    assert!(wide.is_valid(), "{wide:?}");
}

#[test]
fn a_round_hole_near_a_side_is_reported_by_the_arc_edge() {
    // The hole's lower half (edge 0) rises to within half the tolerance of
    // the bottom side (edge 0); its upper half is far away.
    let close = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![circle(5.0, 0.5 + 0.5 * TOL, 0.5)],
        TOL,
    );
    assert_eq!(near(&close), vec![((OUTER, 0), (HOLE, 0), false)]);

    let apart = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![circle(5.0, 0.5 + 2.0 * TOL, 0.5)],
        TOL,
    );
    assert!(apart.is_valid(), "{apart:?}");
}

#[test]
fn a_hole_arc_nearly_tangent_to_the_outer_arc_is_reported() {
    // Internal near-tangency away from the circles' vertices, so only the
    // line-of-centres test sees it: the hole's rim comes within half the
    // tolerance of the outer rim, or twice.
    let hole_at = |gap: f64| {
        let centre = 1.0 - gap;
        circle(0.6 * centre, 0.8 * centre, 1.0)
    };
    let close = with_clearance(circle(0.0, 0.0, 2.0), vec![hole_at(0.5 * TOL)], TOL);
    let pairs = near(&close);
    assert!(!pairs.is_empty(), "{close:?}");
    assert!(pairs
        .iter()
        .all(|&(a, b, undecided)| { a.0 == OUTER && b.0 == HOLE && !undecided }));

    let apart = with_clearance(circle(0.0, 0.0, 2.0), vec![hole_at(2.0 * TOL)], TOL);
    assert!(apart.is_valid(), "{apart:?}");
}

#[test]
fn an_exact_tangency_is_a_contact_and_a_cusp() {
    // The circles touch at their shared vertex (2, 0) and leave it along one
    // direction, down on the lower arcs and up on the upper ones: a zero
    // interior angle, so the arcs stay within the clearance of each other
    // for a stretch. Lower faces upper only back to back, which is no cusp.
    let check = with_clearance(circle(0.0, 0.0, 2.0), vec![circle(1.0, 0.0, 1.0)], TOL);
    assert_eq!(check.contacts.len(), 1);
    assert_eq!(
        near(&check),
        vec![
            ((OUTER, 0), (HOLE, 0), false),
            ((OUTER, 1), (HOLE, 1), false)
        ]
    );
}

#[test]
fn a_cusp_inside_both_edges_is_found() {
    // Circles tangent away from every vertex: the hole at (0, 1) touches the
    // outer circle at (0, 2), inside both upper arcs, decided by the exact
    // line-of-centres test.
    let check = with_clearance(circle(0.0, 0.0, 2.0), vec![circle(0.0, 1.0, 1.0)], TOL);
    assert!(!near(&check).is_empty(), "{check:?}");
    // A line touching a circle along a tangent inside both: the wall of a
    // rectangle hole is tangent to a circle hole at its lowest point.
    let tangent_wall = with_clearance(rect(0.0, 0.0, 10.0, 10.0), vec![circle(5.0, 1.0, 1.0)], TOL);
    assert!(
        near(&tangent_wall).contains(&((OUTER, 0), (HOLE, 0), false)),
        "{tangent_wall:?}"
    );
}

#[test]
fn an_exact_cusp_joint_is_found_but_a_thin_wedge_is_not() {
    // The arc leaves (0, 0) along +x exactly (bulge 0.5 on the chord (6, 8):
    // centre (0, 6.25), radius 6.25), and the closing edge arrives along +x
    // too: a zero interior angle, so the edges stay within 0.5 of each other
    // for about 2.5 units.
    let fine = Tolerance::new(1e-9, 1e-9).expect("valid tolerance");
    let cusp = ArcPolygon {
        outer: ring(&[(0.0, 0.0, 0.5), (6.0, 8.0, 0.0), (11.0, 0.0, 0.0)]),
        holes: vec![],
    };
    for clearance in [0.01, 0.5] {
        let options = ArcCheckOptions::new().with_clearance(clearance);
        let check = check_arc_polygon_with(&cusp, fine, options);
        assert!(
            near(&check).contains(&((OUTER, 0), (OUTER, 2), false)),
            "{clearance}: {check:?}"
        );
    }
    // A wedge of a small but nonzero angle is not reported: the edges are
    // within 0.5 of each other near the apex, but no end point is within
    // 0.5 of the other edge and no common normal is. Its neighbourhood of
    // the meeting point is left out.
    let wedge = ArcPolygon {
        outer: ring(&[(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (10.0, 1.0, 0.0)]),
        holes: vec![],
    };
    let options = ArcCheckOptions::new().with_clearance(0.5);
    let check = check_arc_polygon_with(&wedge, fine, options);
    assert!(check.is_valid(), "{check:?}");
    // It is found as soon as an end comes within the clearance.
    let options = ArcCheckOptions::new().with_clearance(2.0);
    assert!(!check_arc_polygon_with(&wedge, fine, options).is_valid());
}

#[test]
fn a_smooth_joint_is_not_a_cusp() {
    // Two arcs and two lines meeting tangentially in a stadium: the edges
    // leave each joint in opposite directions.
    let stadium = ring(&[
        (0.0, 0.0, 0.0),
        (4.0, 0.0, 1.0),
        (4.0, 2.0, 0.0),
        (0.0, 2.0, 1.0),
    ]);
    let check = with_clearance(stadium, vec![], 0.01);
    assert!(check.is_valid(), "{check:?}");
}

#[test]
fn the_comparison_is_exact_and_strict() {
    // A gap of exactly 2^-10 between two straight walls, and between a
    // wall and a circle: below only for a clearance past it.
    let gap = 2f64.powi(-10);
    let past = gap + 2f64.powi(-40);
    for hole in [rect(gap, 1.0, 5.0, 4.0), circle(5.0, 0.5 + gap, 0.5)] {
        let outer = || rect(0.0, 0.0, 10.0, 10.0);
        let at = with_clearance(outer(), vec![hole.clone()], gap);
        assert!(at.is_valid(), "at the gap: {at:?}");
        let beyond = with_clearance(outer(), vec![hole.clone()], past);
        assert!(!near(&beyond).is_empty(), "past the gap: {beyond:?}");
        assert!(near(&beyond).iter().all(|pair| !pair.2));
    }
}

#[test]
fn a_certified_lower_bound_is_reported_below_the_distance() {
    let close = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(0.5 * TOL, 1.0, 5.0, 4.0)],
        TOL,
    );
    for defect in &close.defects {
        let ArcPolygonDefect::Clearance { lower_bound, .. } = defect else {
            panic!("only clearance defects: {defect:?}");
        };
        assert!(
            (0.4 * TOL..0.5 * TOL).contains(lower_bound),
            "{lower_bound}"
        );
    }
}

#[test]
fn findings_are_ordered_and_repeat() {
    let polygon = ArcPolygon {
        outer: rect(0.0, 0.0, 10.0, 10.0),
        holes: vec![
            rect(0.5 * TOL, 1.0, 5.0, 4.0),
            rect(5.0 + 0.5 * TOL, 1.0, 9.0, 4.0),
        ],
    };
    let options = ArcCheckOptions::new().with_clearance(TOL);
    let first = check_arc_polygon_with(&polygon, tol(), options);
    assert_eq!(first, check_arc_polygon_with(&polygon, tol(), options));
    let pairs = near(&first);
    let mut sorted = pairs.clone();
    sorted.sort_by_key(|&(a, b, _)| (a, b));
    assert_eq!(pairs, sorted);
    // The two holes face each other across 0.5 tol too.
    assert!(pairs
        .iter()
        .any(|&(a, b, _)| a.0 == HOLE && b.0 == ArcRingRef::Hole(1)));
}

#[test]
fn an_invalid_clearance_is_a_defect_and_zero_finds_nothing() {
    let polygon = ArcPolygon {
        outer: rect(0.0, 0.0, 10.0, 10.0),
        holes: vec![rect(0.5 * TOL, 1.0, 5.0, 4.0)],
    };
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        let check =
            check_arc_polygon_with(&polygon, tol(), ArcCheckOptions::new().with_clearance(bad));
        assert_eq!(check.defects.len(), 1, "{check:?}");
        assert!(matches!(
            check.defects[0],
            ArcPolygonDefect::InvalidClearance { .. }
        ));
    }
    let zero = check_arc_polygon_with(&polygon, tol(), ArcCheckOptions::new().with_clearance(0.0));
    assert!(zero.is_valid(), "{zero:?}");
}

#[test]
fn a_crossing_alone_is_not_a_near_miss() {
    // Edges that meet are a contact or a defect already, and the crossing
    // point is at distance zero; the rest of the edges is far apart.
    let check = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(8.0, 4.0, 12.0, 6.0)],
        TOL,
    );
    assert!(
        near(&check).is_empty(),
        "a crossing pair is not a near miss: {check:?}"
    );
}

#[test]
fn a_contact_does_not_hide_a_near_miss_elsewhere_on_the_same_edges() {
    // A crescent hole touches the outer ring's bottom side at (3, 0), a
    // contact, and misses it by 1e-9 at (7, 1e-9), four units away.
    let eps = 1e-9;
    let crescent = ring(&[(3.0, 0.0, -0.5), (7.0, eps, 0.1)]);
    let check = with_clearance(rect(0.0, 0.0, 10.0, 10.0), vec![crescent], 1e-6);
    assert_eq!(check.contacts.len(), 1, "{check:?}");
    assert_eq!(
        near(&check),
        vec![
            ((OUTER, 0), (HOLE, 0), false),
            ((OUTER, 0), (HOLE, 1), false)
        ]
    );
    // Without the near miss only the contact is left, and it is no finding.
    let touching = ring(&[(3.0, 0.0, -0.5), (7.0, 0.0, 0.1)]);
    let check = with_clearance(rect(0.0, 0.0, 10.0, 10.0), vec![touching], 1e-6);
    assert!(near(&check).is_empty(), "{check:?}");
    assert_eq!(check.contacts.len(), 2);
}

#[test]
fn consecutive_edges_that_close_in_on_each_other_are_found() {
    // A sliver of three edges: the two long edges come within 5e-7 of each
    // other at the short edge, away from the joint they share at (0, 0).
    let fine = Tolerance::new(1e-9, 1e-9).expect("valid tolerance");
    let sliver = |width: f64| ArcPolygon {
        outer: ring(&[(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (10.0, width, 0.0)]),
        holes: vec![],
    };
    let options = ArcCheckOptions::new().with_clearance(1e-6);
    let thin = check_arc_polygon_with(&sliver(5e-7), fine, options);
    assert!(
        near(&thin).contains(&((OUTER, 0), (OUTER, 2), false)),
        "{thin:?}"
    );
    let wide = check_arc_polygon_with(&sliver(2e-6), fine, options);
    assert!(wide.is_valid(), "{wide:?}");
}

/// A deterministic generator, so a failing case can be rebuilt.
struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[test]
fn the_verdict_matches_the_analytic_gap_on_random_circles_and_sides() {
    let mut rng = Lcg(0x388);
    let mut cases = 0;
    for _ in 0..300 {
        let angle = rng.unit() * std::f64::consts::TAU;
        let r = 0.3 + rng.unit() * 0.5;
        // A gap within a few clearances either side of the clearance.
        let gap = TOL * (0.2 + rng.unit() * 3.0);
        if (gap - TOL).abs() < 1e-7 {
            continue;
        }
        // Hole circle inside a circle: the gap is R - r - |c|.
        let big = 2.0;
        let c = big - r - gap;
        let hole = circle(c * angle.cos(), c * angle.sin(), r);
        let check = with_clearance(circle(0.0, 0.0, big), vec![hole], TOL);
        assert_eq!(
            near(&check).is_empty(),
            gap > TOL,
            "circles: gap {gap}, angle {angle}"
        );
        // Hole circle above a rectangle's bottom side: the gap is cy - r.
        let (x, y) = (3.0 + rng.unit() * 4.0, r + gap);
        let check = with_clearance(rect(0.0, 0.0, 10.0, 10.0), vec![circle(x, y, r)], TOL);
        assert_eq!(near(&check).is_empty(), gap > TOL, "side: gap {gap}, x {x}");
        cases += 2;
    }
    assert!(cases > 400);
}

fn triangle(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> ArcRing {
    ring(&[(a.0, a.1, 0.0), (b.0, b.1, 0.0), (c.0, c.1, 0.0)])
}

#[test]
fn an_apex_and_a_wall_are_found_whichever_ring_comes_first() {
    // The apex of a triangle comes within half the tolerance of the middle
    // of a wall, so only an end point of the triangle's edges is near.
    let wall = || rect(5.0, 1.0, 9.0, 4.0);
    let apex = || triangle((1.0, 1.0), (1.0, 4.0), (5.0 - 0.5 * TOL, 2.5));
    let outer = || rect(0.0, 0.0, 10.0, 10.0);
    let apex_first = with_clearance(outer(), vec![apex(), wall()], TOL);
    let wall_first = with_clearance(outer(), vec![wall(), apex()], TOL);
    let holes = |check: &ArcPolygonCheck| {
        near(check)
            .into_iter()
            .filter(|&(a, b, _)| a.0 != OUTER && b.0 != OUTER)
            .count()
    };
    assert!(holes(&apex_first) >= 1, "{apex_first:?}");
    assert!(holes(&wall_first) >= 1, "{wall_first:?}");
}

#[test]
fn a_corner_pointing_at_an_arc_is_found_whichever_ring_comes_first() {
    // The apex is within half the tolerance of the circle's lowest point,
    // inside its lower arc: only the radial of an end point is near.
    let apex = || triangle((5.0, 4.5 - 0.5 * TOL), (4.0, 3.0), (6.0, 3.0));
    let round = || circle(5.0, 5.0, 0.5);
    let outer = || rect(0.0, 0.0, 10.0, 10.0);
    for holes in [vec![round(), apex()], vec![apex(), round()]] {
        let check = with_clearance(outer(), holes, TOL);
        assert!(!near(&check).is_empty(), "{check:?}");
    }
}

#[test]
fn only_the_arc_that_faces_a_side_is_reported() {
    // The upper arc is edge 1 of a circle ring, the lower edge 0.
    let top = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![circle(5.0, 10.0 - 0.5 - 0.5 * TOL, 0.5)],
        TOL,
    );
    assert_eq!(near(&top), vec![((OUTER, 2), (HOLE, 1), false)]);
}

#[test]
fn only_the_arc_that_faces_a_rim_is_reported() {
    // A hole in the upper half of a circle, nearly tangent to its upper arc
    // (edge 1) and facing the hole's own upper arc (edge 1).
    let centre = 1.0 - 0.5 * TOL;
    let close = with_clearance(
        circle(0.0, 0.0, 2.0),
        vec![circle(0.6 * centre, 0.8 * centre, 1.0)],
        TOL,
    );
    assert_eq!(near(&close), vec![((OUTER, 1), (HOLE, 1), false)]);
    // Below, the lower arcs face each other.
    let below = with_clearance(
        circle(0.0, 0.0, 2.0),
        vec![circle(0.6 * centre, -0.8 * centre, 1.0)],
        TOL,
    );
    assert_eq!(near(&below), vec![((OUTER, 0), (HOLE, 0), false)]);
}

#[test]
fn two_holes_facing_each_other_on_outer_sides_are_reported_by_their_arcs() {
    // External near-tangency along a diagonal: the centres are 2 + gap
    // apart.
    let gap = 0.5 * TOL;
    let d = (2.0 + gap) / 2f64.sqrt();
    let check = with_clearance(
        rect(-10.0, -10.0, 10.0, 10.0),
        vec![circle(0.0, 0.0, 1.0), circle(d, d, 1.0)],
        TOL,
    );
    assert_eq!(
        near(&check),
        vec![((HOLE, 1), (ArcRingRef::Hole(1), 0), false)]
    );
    let apart = with_clearance(
        rect(-10.0, -10.0, 10.0, 10.0),
        vec![
            circle(0.0, 0.0, 1.0),
            circle(d + 2.0 * TOL, d + 2.0 * TOL, 1.0),
        ],
        TOL,
    );
    assert!(apart.is_valid(), "{apart:?}");
}

#[test]
fn a_short_wall_beside_a_circle_is_not_near_when_its_line_is() {
    // The wall's line passes within half the tolerance of the circle and
    // its box is near, but the circle's nearest point projects below the
    // wall, whose end is further than the clearance away.
    let check = with_clearance(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![circle(5.0, 5.0, 0.5), rect(5.5 + 0.5 * TOL, 5.05, 6.5, 7.0)],
        TOL,
    );
    assert!(check.is_valid(), "{check:?}");
}

#[test]
fn a_wall_beside_a_circle_is_found_whichever_ring_comes_first() {
    // The wall's middle faces the circle's lowest point, which is inside
    // its lower arc, so only the interiors are near.
    let wall = || rect(4.0, 3.0, 6.0, 4.5 - 0.5 * TOL);
    let circle = || circle(5.0, 5.0, 0.5);
    let outer = || rect(0.0, 0.0, 10.0, 10.0);
    for holes in [vec![circle(), wall()], vec![wall(), circle()]] {
        let check = with_clearance(outer(), holes, TOL);
        assert!(!near(&check).is_empty(), "{check:?}");
    }
}

#[test]
fn a_vertex_beside_an_arcs_end_is_near_that_end() {
    // A triangle's apex lies below the diameter of a D shape, where its
    // radial direction misses the arc, within the clearance of the arc's
    // last end point (-1, 0) and of nothing else on it.
    let d_shape = ring(&[(-1.0, 0.0, 0.0), (1.0, 0.0, 1.0)]);
    let apex = triangle((-1.0004, -0.0004), (-3.0, -1.0), (-3.0, -3.0));
    let check = with_clearance(rect(-10.0, -10.0, 10.0, 10.0), vec![d_shape, apex], TOL);
    assert!(
        near(&check)
            .iter()
            .any(|&(a, b, _)| a == (HOLE, 1) && b.0 == ArcRingRef::Hole(1)),
        "{check:?}"
    );
}

#[test]
fn a_ring_that_fails_validation_is_not_tested() {
    // A slot narrower than the ring tolerance is a repeated vertex: the
    // ring is reported for that and its edges are left out.
    let narrow = ring(&[
        (0.0, 0.0, 0.0),
        (10.0, 0.0, 0.0),
        (10.0, 10.0, 0.0),
        (5.0 + 0.5 * TOL, 10.0, 0.0),
        (5.0 + 0.5 * TOL, 2.0, 0.0),
        (5.0, 2.0, 0.0),
        (5.0, 10.0, 0.0),
        (0.0, 10.0, 0.0),
    ]);
    let check = with_clearance(narrow, vec![], TOL);
    assert!(
        matches!(
            check.defects.as_slice(),
            [ArcPolygonDefect::Ring { ring: OUTER, .. }]
        ),
        "{check:?}"
    );
}

#[test]
fn nested_circles_are_found_with_the_inner_one_first() {
    // Holes that nest are a defect of their own; their near rims are still
    // reported, here with the inner circle as the first edge.
    let centre = 1.0 - 0.5 * TOL;
    let check = with_clearance(
        rect(-10.0, -10.0, 10.0, 10.0),
        vec![
            circle(0.6 * centre, 0.8 * centre, 1.0),
            circle(0.0, 0.0, 2.0),
        ],
        TOL,
    );
    assert!(
        near(&check).contains(&((HOLE, 1), (ArcRingRef::Hole(1), 1), false)),
        "{check:?}"
    );
}

#[test]
fn collinear_walls_apart_by_twice_the_tolerance_are_clear() {
    // The gap between them is along their common line, where the distance
    // to the line is zero.
    let fine = Tolerance::new(1e-9, 1e-9).expect("valid tolerance");
    let slot = |width: f64| {
        ring(&[
            (0.0, 0.0, 0.0),
            (10.0, 0.0, 0.0),
            (10.0, 10.0, 0.0),
            (5.0 + width, 10.0, 0.0),
            (5.0 + width, 2.0, 0.0),
            (5.0, 2.0, 0.0),
            (5.0, 10.0, 0.0),
            (0.0, 10.0, 0.0),
        ])
    };
    let options = ArcCheckOptions::new().with_clearance(TOL);
    let wide = check_arc_polygon_with(
        &ArcPolygon {
            outer: slot(2.0 * TOL),
            holes: vec![],
        },
        fine,
        options,
    );
    assert!(wide.is_valid(), "{wide:?}");
}

#[test]
fn a_hole_near_a_joint_of_the_outer_circle_is_found_on_the_arc_it_faces() {
    // Near the circle's two vertices, where a point's radial direction
    // falls in one arc and the nearest end of the other is its joint.
    for (x, expect_edges) in [(-1.0, [0usize, 1]), (1.0, [0, 1])] {
        let r = 0.5;
        let centre = x * (2.0 - r - 0.5 * TOL);
        let close = with_clearance(circle(0.0, 0.0, 2.0), vec![circle(centre, 0.0, r)], TOL);
        let outer_edges: Vec<usize> = near(&close).iter().map(|pair| pair.0 .1).collect();
        assert!(!outer_edges.is_empty(), "{close:?}");
        assert!(outer_edges.iter().all(|e| expect_edges.contains(e)));
        let far = centre * (1.0 - 2.0 * TOL);
        let apart = with_clearance(circle(0.0, 0.0, 2.0), vec![circle(far, 0.0, r)], TOL);
        assert!(near(&apart).is_empty(), "{apart:?}");
    }
}

#[test]
fn a_hole_near_a_major_arc_is_found_whichever_way_the_arc_runs() {
    // Bulge 2 is an arc of about 254 degrees on the circle of radius 1.25
    // about (0, -0.75), closed by its chord. A hole at 200 degrees from the
    // centre is on the arc, but not within half a turn of both its ends.
    let hole = |degrees: f64, gap: f64| {
        let angle = degrees.to_radians();
        let reach = 1.25 - 0.1 - gap;
        circle(reach * angle.cos(), -0.75 + reach * angle.sin(), 0.1)
    };
    let ccw = || ring(&[(-1.0, 0.0, 2.0), (1.0, 0.0, 0.0)]);
    let cw = || ring(&[(1.0, 0.0, -2.0), (-1.0, 0.0, 0.0)]);
    for outer in [ccw, cw] {
        for degrees in [200.0, 270.0, 330.0] {
            let close = with_clearance(outer(), vec![hole(degrees, 0.5 * TOL)], TOL);
            assert!(!near(&close).is_empty(), "{degrees}: {close:?}");
            let apart = with_clearance(outer(), vec![hole(degrees, 2.0 * TOL)], TOL);
            assert!(apart.is_valid(), "{degrees}: {apart:?}");
        }
    }
}

#[test]
fn a_point_at_an_arcs_centre_is_its_radius_away() {
    // A tiny circle around the origin whose radius is above the clearance,
    // and a triangle inside it with a corner at the centre and the other two
    // corners nearer the centre than the rim by more than the clearance.
    let fine = Tolerance::new(1e-9, 1e-9).expect("valid tolerance");
    let options = ArcCheckOptions::new().with_clearance(0.0003);
    let polygon = ArcPolygon {
        outer: circle(0.0, 0.0, 0.0004),
        holes: vec![triangle((0.0, 0.0), (0.00005, 0.0), (0.0, 0.00005))],
    };
    // The triangle is smaller than the clearance, so its own edges are near;
    // only the pairs with the circle matter here.
    let with_circle = |check: &ArcPolygonCheck| near(check).iter().any(|pair| pair.0 .0 == OUTER);
    let check = check_arc_polygon_with(&polygon, fine, options);
    assert!(!with_circle(&check), "{check:?}");
    let options = ArcCheckOptions::new().with_clearance(0.0005);
    let check = check_arc_polygon_with(&polygon, fine, options);
    assert!(with_circle(&check), "{check:?}");
}
