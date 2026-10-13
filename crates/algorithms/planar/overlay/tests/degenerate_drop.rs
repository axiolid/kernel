//! Opt-in dropping of degenerate input rings, every drop reported (#256).
//!
//! By default one degenerate ring refuses the whole overlay. Asked
//! explicitly, `overlay_with` and `union_soup_with` leave out a ring that
//! is too short, has a repeated vertex or has zero area under the exact
//! `ZeroArea` rule, report it with its place, reason and area, and return
//! the exact overlay of the remaining rings. A hole goes with its dropped
//! outer. Self-intersection, non-finite points, holes outside their outer
//! and invalid frames stay refusals.

use axiolid_core::{Frame2, Point2, Tolerance, Vec2};
use axiolid_overlay::{
    overlay, overlay_with, union_soup, union_soup_with, DegenerateRings, DropReason, DroppedRing,
    FillRule, OverlayError, OverlayInput, OverlayOperand, OverlayOperation, OverlayOptions,
    Polygon, ReportedOverlay, Ring, RingRole,
};

/// Linear tolerance 1e-3, so `ZeroArea` is an area of at most 1e-6.
fn tol() -> Tolerance {
    Tolerance::new(1e-3, 1e-9).unwrap()
}

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
    }
}

fn poly(points: &[(f64, f64)]) -> Polygon {
    Polygon {
        outer: ring(points),
        holes: Vec::new(),
    }
}

fn with_holes(outer: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Polygon {
    Polygon {
        outer: ring(outer),
        holes: holes.iter().map(|h| ring(h)).collect(),
    }
}

fn input(polygons: Vec<Polygon>) -> OverlayInput {
    OverlayInput {
        frame: frame(),
        polygons,
    }
}

fn square(x: f64, y: f64, side: f64) -> Vec<(f64, f64)> {
    vec![(x, y), (x + side, y), (x + side, y + side), (x, y + side)]
}

/// A unit square at `(x, y)` with its third corner given twice.
fn doubled_corner(x: f64, y: f64) -> Vec<(f64, f64)> {
    vec![
        (x, y),
        (x + 1.0, y),
        (x + 1.0, y + 1.0),
        (x + 1.0, y + 1.0),
        (x, y + 1.0),
    ]
}

/// A triangle of area `2^-25`, below the `ZeroArea` bound of 1e-6 and not
/// on one line. A triangle's edges are all adjacent, so it never
/// intersects itself.
fn sliver(x: f64, y: f64) -> Vec<(f64, f64)> {
    vec![(x, y), (x + 1.0, y), (x + 0.5, y + 2f64.powi(-24))]
}

fn dropping() -> OverlayOptions {
    OverlayOptions::new().dropping_degenerate_rings()
}

fn union_with(
    subject: Vec<Polygon>,
    clip: Vec<Polygon>,
    options: OverlayOptions,
) -> Result<ReportedOverlay, OverlayError> {
    overlay_with(
        &input(subject),
        &input(clip),
        OverlayOperation::Union,
        FillRule::NonZero,
        tol(),
        options,
    )
}

fn plain_union(subject: Vec<Polygon>, clip: Vec<Polygon>) -> Vec<Polygon> {
    overlay(
        &input(subject),
        &input(clip),
        OverlayOperation::Union,
        FillRule::NonZero,
        tol(),
    )
    .expect("valid operands")
    .polygons
}

fn assert_dropped(
    d: &DroppedRing,
    operand: OverlayOperand,
    polygon: usize,
    role: RingRole,
    reason: DropReason,
) {
    assert_eq!(d.operand, operand, "{d:?}");
    assert_eq!(d.polygon, polygon, "{d:?}");
    assert_eq!(d.ring, role, "{d:?}");
    assert_eq!(d.reason, reason, "{d:?}");
    assert!(d.area <= d.area_bound, "{d:?}");
}

/// The bound is an upper bound, and tight: within a few ulps of the area.
fn assert_area(d: &DroppedRing, expected: f64) {
    assert_eq!(d.area, expected, "{d:?}");
    assert!(d.area_bound >= expected, "{d:?}");
    assert!(
        d.area_bound <= expected * (1.0 + 1e-14) + f64::MIN_POSITIVE,
        "{d:?}"
    );
}

#[test]
fn default_options_refuse_each_degenerate_ring_as_overlay_does() {
    let cases: [(Vec<(f64, f64)>, OverlayError); 3] = [
        (doubled_corner(3.0, 0.0), OverlayError::RepeatedVertex),
        (sliver(3.0, 0.0), OverlayError::ZeroArea),
        (vec![(3.0, 0.0), (4.0, 0.0)], OverlayError::RingTooShort),
    ];
    assert_eq!(OverlayOptions::default(), OverlayOptions::new());
    assert_eq!(
        OverlayOptions::new().degenerate_rings,
        DegenerateRings::Refuse
    );
    for (bad, error) in cases {
        let subject = vec![poly(&square(0.0, 0.0, 1.0)), poly(&bad)];
        let refused = overlay(
            &input(subject.clone()),
            &input(Vec::new()),
            OverlayOperation::Union,
            FillRule::NonZero,
            tol(),
        );
        assert_eq!(refused, Err(error.clone()));
        assert_eq!(
            union_with(subject.clone(), Vec::new(), OverlayOptions::new()),
            Err(error.clone())
        );
        let rings: Vec<Ring> = subject.into_iter().map(|p| p.outer).collect();
        assert_eq!(union_soup(&rings, tol()), Err(error.clone()));
        assert_eq!(
            union_soup_with(&rings, tol(), OverlayOptions::new()),
            Err(error)
        );
    }
}

#[test]
fn default_options_change_nothing_on_valid_operands() {
    let subject = vec![poly(&square(0.0, 0.0, 2.0))];
    let clip = vec![poly(&square(1.0, 1.0, 2.0))];
    let reported = union_with(subject.clone(), clip.clone(), OverlayOptions::new()).unwrap();
    assert!(reported.dropped.is_empty());
    let plain = overlay(
        &input(subject),
        &input(clip),
        OverlayOperation::Union,
        FillRule::NonZero,
        tol(),
    )
    .unwrap();
    assert_eq!(reported.result, plain);
}

#[test]
fn a_repeated_vertex_ring_is_dropped_whole_and_its_area_reported() {
    let valid = poly(&square(0.0, 0.0, 1.0));
    let clip = vec![poly(&square(0.5, 0.5, 1.0))];
    let reported = union_with(
        vec![valid.clone(), poly(&doubled_corner(3.0, 0.0))],
        clip.clone(),
        dropping(),
    )
    .unwrap();
    assert_eq!(reported.result.polygons, plain_union(vec![valid], clip));
    assert_eq!(reported.result.evidence.subject_rings, 1);
    assert_eq!(reported.result.evidence.clip_rings, 1);
    assert_eq!(reported.dropped.len(), 1);
    let d = &reported.dropped[0];
    assert_dropped(
        d,
        OverlayOperand::Subject,
        1,
        RingRole::Outer,
        DropReason::RepeatedVertex,
    );
    // Dropping a repeated-vertex ring can remove real area; the report
    // says how much.
    assert_area(d, 1.0);
}

#[test]
fn a_zero_area_ring_is_dropped_with_its_area_within_the_rule() {
    let valid = poly(&square(0.0, 0.0, 1.0));
    let reported = union_with(
        vec![poly(&sliver(3.0, 0.0)), valid.clone()],
        Vec::new(),
        dropping(),
    )
    .unwrap();
    assert_eq!(
        reported.result.polygons,
        plain_union(vec![valid], Vec::new())
    );
    // The dropped ring took no part: the boolean saw the valid ring only.
    assert_eq!(reported.result.evidence.subject_rings, 1);
    assert_eq!(reported.dropped.len(), 1);
    let d = &reported.dropped[0];
    assert_dropped(
        d,
        OverlayOperand::Subject,
        0,
        RingRole::Outer,
        DropReason::ZeroArea,
    );
    assert_area(d, 2f64.powi(-25));
    assert!(d.area_bound <= tol().linear() * tol().linear());
}

#[test]
fn a_ring_just_above_the_zero_area_rule_is_kept() {
    // Area 1.5e-6 > 1e-6: valid, so it is used, not dropped.
    let thin = vec![(3.0, 0.0), (4.0, 0.0), (3.5, 3e-6)];
    let reported = union_with(vec![poly(&thin)], Vec::new(), dropping()).unwrap();
    assert!(reported.dropped.is_empty());
    assert_eq!(reported.result.polygons.len(), 1);
}

#[test]
fn a_ring_too_short_is_dropped_with_zero_area() {
    let valid = poly(&square(0.0, 0.0, 1.0));
    let reported = union_with(
        vec![valid.clone(), poly(&[(3.0, 0.0), (4.0, 0.0)])],
        Vec::new(),
        dropping(),
    )
    .unwrap();
    assert_eq!(
        reported.result.polygons,
        plain_union(vec![valid], Vec::new())
    );
    // The dropped ring took no part: the boolean saw the valid ring only.
    assert_eq!(reported.result.evidence.subject_rings, 1);
    assert_eq!(reported.dropped.len(), 1);
    let d = &reported.dropped[0];
    assert_dropped(
        d,
        OverlayOperand::Subject,
        1,
        RingRole::Outer,
        DropReason::TooShort,
    );
    assert_eq!((d.area, d.area_bound), (0.0, 0.0));
}

#[test]
fn holes_of_a_dropped_outer_are_dropped_and_reported() {
    let valid = poly(&square(0.0, 0.0, 1.0));
    // Outer: a 4x4 square with a doubled corner. Hole 0 is valid on its
    // own and goes with its outer; hole 1 is a sliver and keeps its own
    // reason.
    let outer = vec![
        (10.0, 0.0),
        (14.0, 0.0),
        (14.0, 4.0),
        (14.0, 4.0),
        (10.0, 4.0),
    ];
    let hole = vec![(11.0, 1.0), (11.0, 2.0), (12.0, 2.0), (12.0, 1.0)];
    let thin = vec![(11.0, 3.0), (11.5, 3.0 + 2f64.powi(-24)), (12.0, 3.0)];
    let reported = union_with(
        vec![valid.clone(), with_holes(&outer, &[&hole, &thin])],
        Vec::new(),
        dropping(),
    )
    .unwrap();
    assert_eq!(
        reported.result.polygons,
        plain_union(vec![valid], Vec::new())
    );
    assert_eq!(reported.result.evidence.subject_rings, 1);
    assert_eq!(reported.dropped.len(), 3);
    let [o, h0, h1] = &reported.dropped[..] else {
        unreachable!()
    };
    assert_dropped(
        o,
        OverlayOperand::Subject,
        1,
        RingRole::Outer,
        DropReason::RepeatedVertex,
    );
    assert_area(o, 16.0);
    assert_dropped(
        h0,
        OverlayOperand::Subject,
        1,
        RingRole::Hole(0),
        DropReason::OuterDropped,
    );
    assert_area(h0, 1.0);
    assert_dropped(
        h1,
        OverlayOperand::Subject,
        1,
        RingRole::Hole(1),
        DropReason::ZeroArea,
    );
    assert_area(h1, 2f64.powi(-25));
}

#[test]
fn a_degenerate_hole_of_a_kept_outer_is_dropped_alone() {
    let outer = square(0.0, 0.0, 4.0);
    let good = vec![(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)];
    let bad = vec![(2.5, 2.5), (2.5, 3.5), (3.5, 3.5), (3.5, 3.5), (3.5, 2.5)];
    let reported = union_with(
        vec![with_holes(&outer, &[&bad, &good])],
        Vec::new(),
        dropping(),
    )
    .unwrap();
    assert_eq!(
        reported.result.polygons,
        plain_union(vec![with_holes(&outer, &[&good])], Vec::new())
    );
    assert_eq!(reported.result.evidence.subject_rings, 2);
    assert_eq!(reported.result.evidence.output_holes, 1);
    assert_eq!(reported.dropped.len(), 1);
    let d = &reported.dropped[0];
    assert_dropped(
        d,
        OverlayOperand::Subject,
        0,
        RingRole::Hole(0),
        DropReason::RepeatedVertex,
    );
    assert_area(d, 1.0);
}

#[test]
fn a_hole_outside_its_kept_outer_is_still_refused() {
    let outer = square(0.0, 0.0, 4.0);
    let outside = vec![(5.0, 1.0), (5.0, 2.0), (6.0, 2.0), (6.0, 1.0)];
    assert_eq!(
        union_with(
            vec![with_holes(&outer, &[&outside])],
            Vec::new(),
            dropping()
        ),
        Err(OverlayError::HoleOutsideOuter)
    );
}

#[test]
fn a_self_intersecting_ring_is_never_dropped() {
    let valid = poly(&square(0.0, 0.0, 1.0));
    let bowtie = vec![(3.0, 0.0), (4.0, 1.0), (4.0, 0.0), (3.0, 1.0)];
    assert_eq!(
        union_with(vec![valid.clone(), poly(&bowtie)], Vec::new(), dropping()),
        Err(OverlayError::SelfIntersection)
    );
    // A repeated vertex does not make a bowtie droppable: merged, it
    // still crosses itself, and its signed area (zero here) is not what
    // leaving it out would remove.
    let doubled = vec![(3.0, 0.0), (4.0, 1.0), (4.0, 1.0), (4.0, 0.0), (3.0, 1.0)];
    let without = OverlayOptions::new();
    assert_eq!(
        union_with(vec![valid.clone(), poly(&doubled)], Vec::new(), without),
        Err(OverlayError::RepeatedVertex)
    );
    assert_eq!(
        union_with(vec![valid.clone(), poly(&doubled)], Vec::new(), dropping()),
        Err(OverlayError::SelfIntersection)
    );
    // The doubled vertex last, so merging must also close the ring.
    let wrapped = vec![(3.0, 0.0), (4.0, 1.0), (4.0, 0.0), (3.0, 1.0), (3.0, 0.0)];
    assert_eq!(
        union_with(vec![valid.clone(), poly(&wrapped)], Vec::new(), dropping()),
        Err(OverlayError::SelfIntersection)
    );
    let rings = vec![valid.outer, ring(&doubled)];
    assert_eq!(
        union_soup_with(&rings, tol(), dropping()),
        Err(OverlayError::SelfIntersection)
    );
}

#[test]
fn a_ring_whose_doubled_vertex_closes_it_is_droppable() {
    // The last vertex repeats the first: merging must not report the
    // ring as crossing itself.
    let closed = vec![(3.0, 0.0), (4.0, 0.0), (4.0, 1.0), (3.0, 1.0), (3.0, 0.0)];
    let reported = union_with(vec![poly(&closed)], Vec::new(), dropping()).unwrap();
    assert_eq!(reported.dropped.len(), 1);
    assert_eq!(reported.dropped[0].reason, DropReason::RepeatedVertex);
    assert_area(&reported.dropped[0], 1.0);
}

#[test]
fn non_finite_points_and_invalid_frames_are_still_refused() {
    let valid = poly(&square(0.0, 0.0, 1.0));
    let nan = vec![(3.0, 0.0), (f64::NAN, 0.0), (3.0, 1.0)];
    assert_eq!(
        union_with(vec![valid.clone(), poly(&nan)], Vec::new(), dropping()),
        Err(OverlayError::NonFinitePoint)
    );
    // Too short and non-finite: no area can be stated, so not droppable.
    let short = vec![(3.0, 0.0), (f64::INFINITY, 0.0)];
    assert_eq!(
        union_with(vec![valid.clone(), poly(&short)], Vec::new(), dropping()),
        Err(OverlayError::NonFinitePoint)
    );
    let mut skewed = input(vec![valid]);
    skewed.frame.y = Vec2::new(1.0, 1.0);
    assert_eq!(
        overlay_with(
            &skewed,
            &skewed,
            OverlayOperation::Union,
            FillRule::NonZero,
            tol(),
            dropping(),
        ),
        Err(OverlayError::InvalidFrame)
    );
}

#[test]
fn clip_rings_are_reported_under_the_clip() {
    let subject = vec![poly(&square(0.0, 0.0, 2.0))];
    let kept = poly(&square(1.0, 1.0, 2.0));
    let clip = vec![
        poly(&sliver(5.0, 5.0)),
        kept.clone(),
        poly(&doubled_corner(8.0, 0.0)),
    ];
    let reported = overlay_with(
        &input(subject.clone()),
        &input(clip),
        OverlayOperation::Difference,
        FillRule::NonZero,
        tol(),
        dropping(),
    )
    .unwrap();
    let plain = overlay(
        &input(subject),
        &input(vec![kept]),
        OverlayOperation::Difference,
        FillRule::NonZero,
        tol(),
    )
    .unwrap();
    assert_eq!(reported.result, plain);
    assert_eq!(reported.dropped.len(), 2);
    assert_dropped(
        &reported.dropped[0],
        OverlayOperand::Clip,
        0,
        RingRole::Outer,
        DropReason::ZeroArea,
    );
    assert_dropped(
        &reported.dropped[1],
        OverlayOperand::Clip,
        2,
        RingRole::Outer,
        DropReason::RepeatedVertex,
    );
}

#[test]
fn subject_drops_are_reported_before_clip_drops() {
    let reported = union_with(
        vec![poly(&square(0.0, 0.0, 1.0)), poly(&sliver(3.0, 0.0))],
        vec![poly(&doubled_corner(6.0, 0.0))],
        dropping(),
    )
    .unwrap();
    let operands: Vec<_> = reported.dropped.iter().map(|d| d.operand).collect();
    assert_eq!(operands, [OverlayOperand::Subject, OverlayOperand::Clip]);
}

#[test]
fn all_rings_degenerate_gives_an_empty_result_and_a_full_report() {
    let reported = union_with(
        vec![poly(&sliver(0.0, 0.0)), poly(&doubled_corner(3.0, 0.0))],
        vec![poly(&[(9.0, 9.0)])],
        dropping(),
    )
    .unwrap();
    assert!(reported.result.polygons.is_empty());
    assert_eq!(reported.result.evidence.subject_rings, 0);
    assert_eq!(reported.result.evidence.clip_rings, 0);
    assert_eq!(reported.result.evidence.output_polygons, 0);
    let reasons: Vec<_> = reported.dropped.iter().map(|d| d.reason).collect();
    assert_eq!(
        reasons,
        [
            DropReason::ZeroArea,
            DropReason::RepeatedVertex,
            DropReason::TooShort
        ]
    );
    let soup = union_soup_with(
        &[ring(&sliver(0.0, 0.0)), ring(&doubled_corner(3.0, 0.0))],
        tol(),
        dropping(),
    )
    .unwrap();
    assert!(soup.polygons.is_empty());
    assert_eq!(soup.dropped.len(), 2);
}

#[test]
fn union_soup_drops_and_reports_by_ring_index() {
    let rings = vec![
        ring(&square(0.0, 0.0, 2.0)),
        ring(&doubled_corner(5.0, 0.0)),
        ring(&square(1.0, 1.0, 2.0)),
        ring(&sliver(9.0, 0.0)),
    ];
    let reported = union_soup_with(&rings, tol(), dropping()).unwrap();
    let valid = [rings[0].clone(), rings[2].clone()];
    assert_eq!(reported.polygons, union_soup(&valid, tol()).unwrap());
    assert_eq!(reported.dropped.len(), 2);
    assert_dropped(
        &reported.dropped[0],
        OverlayOperand::Subject,
        1,
        RingRole::Outer,
        DropReason::RepeatedVertex,
    );
    assert_area(&reported.dropped[0], 1.0);
    assert_dropped(
        &reported.dropped[1],
        OverlayOperand::Subject,
        3,
        RingRole::Outer,
        DropReason::ZeroArea,
    );
    // Default options: union_soup exactly.
    assert_eq!(
        union_soup_with(&valid, tol(), OverlayOptions::new())
            .unwrap()
            .polygons,
        union_soup(&valid, tol()).unwrap()
    );
}

/// Offsets at which every ring here translates exactly: georeferenced
/// plan positions, and 1e7 m.
const OFFSETS: [(f64, f64); 3] = [
    (524_288.0, 5_242_880.0),
    (655_360.0, 5_636_096.0),
    (8_388_608.0, 8_388_608.0),
];

#[test]
fn far_from_origin_rings_drop_and_report_as_near_the_origin() {
    let near = union_with(
        vec![
            poly(&square(0.0, 0.0, 1.0)),
            poly(&sliver(3.0, 0.0)),
            poly(&doubled_corner(6.0, 0.0)),
        ],
        Vec::new(),
        dropping(),
    )
    .unwrap();
    for (dx, dy) in OFFSETS {
        let far = union_with(
            vec![
                poly(&square(dx, dy, 1.0)),
                poly(&sliver(dx + 3.0, dy)),
                poly(&doubled_corner(dx + 6.0, dy)),
            ],
            Vec::new(),
            dropping(),
        )
        .unwrap();
        assert_eq!(far.dropped, near.dropped, "at ({dx}, {dy})");
        assert_eq!(far.result.polygons.len(), 1);
        let shifted: Vec<Point2> = near.result.polygons[0]
            .outer
            .points
            .iter()
            .map(|p| Point2::new(p.x + dx, p.y + dy))
            .collect();
        assert_eq!(far.result.polygons[0].outer.points, shifted);
    }
}
