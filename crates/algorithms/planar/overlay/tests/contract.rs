use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::*;

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::ZERO,
        x: Point2::X,
        y: Point2::Y,
    }
}
fn ring(x: f64, y: f64, w: f64, h: f64) -> Ring {
    Ring {
        points: vec![
            Point2::new(x, y),
            Point2::new(x + w, y),
            Point2::new(x + w, y + h),
            Point2::new(x, y + h),
        ],
    }
}
fn input(x: f64, y: f64, w: f64, h: f64) -> OverlayInput {
    OverlayInput {
        frame: frame(),
        polygons: vec![Polygon {
            outer: ring(x, y, w, h),
            holes: vec![],
        }],
    }
}
fn area(result: &OverlayResult) -> f64 {
    result
        .polygons
        .iter()
        .map(|p| {
            let r = &p.outer.points;
            let a = r
                .iter()
                .zip(r.iter().cycle().skip(1))
                .take(r.len())
                .map(|(a, b)| a.x * b.y - b.x * a.y)
                .sum::<f64>()
                .abs()
                / 2.0;
            a
        })
        .sum()
}
#[test]
fn operations_have_expected_area_and_evidence() {
    let a = input(0.0, 0.0, 2.0, 2.0);
    let b = input(1.0, 0.0, 2.0, 2.0);
    for (op, want) in [
        (OverlayOperation::Intersection, 2.0),
        (OverlayOperation::Union, 6.0),
        (OverlayOperation::Difference, 2.0),
        (OverlayOperation::Xor, 4.0),
    ] {
        let got = overlay(&a, &b, op, FillRule::EvenOdd, Tolerance::ZERO).unwrap();
        assert!((area(&got) - want).abs() < 1e-12, "{op:?}");
        assert_eq!(got.evidence.subject_rings, 1);
        assert_eq!(got.evidence.clip_rings, 1);
        assert_eq!(got.evidence.output_polygons, got.polygons.len());
    }
}

#[test]
fn supports_all_fill_rules_and_is_repeatable() {
    let a = input(0.0, 0.0, 2.0, 2.0);
    let b = input(1.0, 0.0, 2.0, 2.0);
    for fill in [
        FillRule::EvenOdd,
        FillRule::NonZero,
        FillRule::Positive,
        FillRule::Negative,
    ] {
        let first = overlay(&a, &b, OverlayOperation::Union, fill, Tolerance::ZERO).unwrap();
        for _ in 0..8 {
            assert_eq!(
                overlay(&a, &b, OverlayOperation::Union, fill, Tolerance::ZERO).unwrap(),
                first
            );
        }
    }
}

#[test]
fn accepts_explicit_translated_frame_and_rejects_frame_mismatch() {
    let mut a = input(0.0, 0.0, 1.0, 1.0);
    a.frame.origin = Point2::new(20.0, -8.0);
    let mut b = a.clone();
    b.polygons[0].outer = ring(0.5, 0.0, 1.0, 1.0);
    assert!(overlay(
        &a,
        &b,
        OverlayOperation::Intersection,
        FillRule::NonZero,
        Tolerance::ZERO
    )
    .is_ok());
    b.frame.origin.x += 1.0;
    assert_eq!(
        overlay(
            &a,
            &b,
            OverlayOperation::Intersection,
            FillRule::NonZero,
            Tolerance::ZERO
        ),
        Err(OverlayError::InvalidFrame)
    );
}

#[test]
fn tolerance_and_invalid_ring_hole_are_rejected() {
    let mut a = input(0.0, 0.0, 1.0, 1.0);
    a.polygons[0].outer.points.insert(1, Point2::new(1e-4, 0.0));
    assert_eq!(
        overlay(
            &a,
            &input(2.0, 0.0, 1.0, 1.0),
            OverlayOperation::Union,
            FillRule::NonZero,
            Tolerance::new(1e-3, 1e-3).unwrap()
        ),
        Err(OverlayError::RepeatedVertex)
    );
    let mut invalid = input(0.0, 0.0, 1.0, 1.0);
    invalid.polygons[0].holes.push(ring(2.0, 2.0, 1.0, 1.0));
    assert_eq!(
        overlay(
            &invalid,
            &input(4.0, 0.0, 1.0, 1.0),
            OverlayOperation::Union,
            FillRule::NonZero,
            Tolerance::ZERO
        ),
        Err(OverlayError::HoleOutsideOuter)
    );
}

#[test]
fn valid_hole_is_accepted_and_reported() {
    let mut a = input(0.0, 0.0, 4.0, 4.0);
    a.polygons[0].holes.push(ring(1.0, 1.0, 1.0, 1.0));
    let b = input(10.0, 0.0, 1.0, 1.0);
    let got = overlay(
        &a,
        &b,
        OverlayOperation::Union,
        FillRule::EvenOdd,
        Tolerance::ZERO,
    )
    .unwrap();
    assert_eq!(got.evidence.subject_rings, 2);
    assert_eq!(got.evidence.clip_rings, 1);
    assert!(got.polygons.iter().any(|p| p.holes.len() == 1));
    assert_eq!(got.evidence.output_holes, 1);
}

#[test]
fn rejects_self_crossing_ring() {
    let mut a = input(0.0, 0.0, 1.0, 1.0);
    a.polygons[0].outer.points = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
        Point2::new(1.0, 0.0),
    ];
    assert_eq!(
        overlay(
            &a,
            &input(3.0, 0.0, 1.0, 1.0),
            OverlayOperation::Union,
            FillRule::EvenOdd,
            Tolerance::ZERO
        ),
        Err(OverlayError::SelfIntersection)
    );
}

/// A ring with two collinear edges that do not touch (the inner walls of a
/// U, the teeth of a comb) is simple. The self-crossing check once treated
/// "on the same line" as "touching" and refused every such outline, so a
/// U-shaped wall could not take part in any boolean.
#[test]
fn accepts_collinear_edges_that_do_not_touch() {
    let p = Point2::new;
    for (what, points) in [
        (
            "U",
            vec![
                p(0.0, 0.0),
                p(3.0, 0.0),
                p(3.0, 2.0),
                p(2.0, 2.0),
                p(2.0, 1.0),
                p(1.0, 1.0),
                p(1.0, 2.0),
                p(0.0, 2.0),
            ],
        ),
        (
            "comb",
            vec![
                p(0.0, 0.0),
                p(5.0, 0.0),
                p(5.0, 3.0),
                p(4.0, 3.0),
                p(4.0, 1.0),
                p(3.0, 1.0),
                p(3.0, 3.0),
                p(2.0, 3.0),
                p(2.0, 1.0),
                p(1.0, 1.0),
                p(1.0, 3.0),
                p(0.0, 3.0),
            ],
        ),
    ] {
        let mut a = input(0.0, 0.0, 1.0, 1.0);
        a.polygons[0].outer.points = points;
        let got = overlay(
            &a,
            &input(10.0, 0.0, 1.0, 1.0),
            OverlayOperation::Union,
            FillRule::NonZero,
            Tolerance::METRE,
        );
        assert!(got.is_ok(), "{what}: {got:?}");
    }
}

/// The narrowing must not open the door the other way: a vertex that really
/// lands on another edge of its own ring is still a self-intersection.
#[test]
fn rejects_a_vertex_touching_its_own_edge() {
    let p = Point2::new;
    let mut a = input(0.0, 0.0, 1.0, 1.0);
    // The notch tip (2, 0) lies on the bottom edge (0,0)-(4,0).
    a.polygons[0].outer.points = vec![
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 2.0),
        p(3.0, 2.0),
        p(2.0, 0.0),
        p(1.0, 2.0),
        p(0.0, 2.0),
    ];
    assert_eq!(
        overlay(
            &a,
            &input(10.0, 0.0, 1.0, 1.0),
            OverlayOperation::Union,
            FillRule::NonZero,
            Tolerance::METRE,
        ),
        Err(OverlayError::SelfIntersection)
    );
}
