//! Arc ring contract: areas against closed forms, refusals by reason.
//!
//! Every area assertion compares against a value derived independently
//! (pi r^2, a square plus a half-disc, and so on), never against another
//! call of the same function.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{
    arc_edge_radius, arc_ring_area, reverse_arc_ring, validate_arc_ring, ArcRing, ArcVertex,
    OverlayError,
};

/// Bulge of a quarter-circle edge: tan(90deg / 4).
fn quarter_bulge() -> f64 {
    (std::f64::consts::FRAC_PI_8).tan()
}

#[test]
fn a_two_vertex_circle_encloses_pi_r_squared() {
    let ring = ArcRing::circle(Point2::new(3.0, -2.0), 2.5);
    assert_eq!(ring.edge_count(), 2);
    assert!(!ring.is_polygonal());

    let expected = std::f64::consts::PI * 2.5 * 2.5;
    let error = (arc_ring_area(&ring) - expected).abs();
    assert!(error < 1.0e-12, "area off by {error:e}");
}

#[test]
fn a_circle_split_into_quarters_has_the_same_area() {
    // Four quarter arcs, same circle: a different encoding of one shape
    // must not change the measured area.
    let bulge = quarter_bulge();
    let ring = ArcRing::new(vec![
        ArcVertex::bulged(Point2::new(-1.0, 0.0), bulge),
        ArcVertex::bulged(Point2::new(0.0, -1.0), bulge),
        ArcVertex::bulged(Point2::new(1.0, 0.0), bulge),
        ArcVertex::bulged(Point2::new(0.0, 1.0), bulge),
    ]);

    let error = (arc_ring_area(&ring) - std::f64::consts::PI).abs();
    assert!(error < 1.0e-12, "quartered circle off by {error:e}");
}

#[test]
fn a_bulge_sign_distinguishes_a_bump_from_a_bite() {
    // Same 2x2 square, same chord, opposite bulge on the top edge. One
    // adds a half-disc of radius 1, the other removes it. If the segment
    // term were unsigned both would read the same.
    let square = |bulge: f64| {
        ArcRing::new(vec![
            ArcVertex::straight(Point2::new(0.0, 0.0)),
            ArcVertex::straight(Point2::new(2.0, 0.0)),
            ArcVertex::bulged(Point2::new(2.0, 2.0), bulge),
            ArcVertex::straight(Point2::new(0.0, 2.0)),
        ])
    };

    let half_disc = std::f64::consts::PI / 2.0;
    let bump = arc_ring_area(&square(1.0));
    let bite = arc_ring_area(&square(-1.0));

    assert!((bump - (4.0 + half_disc)).abs() < 1.0e-12, "bump {bump}");
    assert!((bite - (4.0 - half_disc)).abs() < 1.0e-12, "bite {bite}");
}

#[test]
fn reversing_a_ring_negates_its_area_and_preserves_geometry() {
    // Bulges are stored per departing edge, so a naive vertex reversal
    // attaches each bulge to the wrong edge. Area must negate EXACTLY.
    let ring = ArcRing::new(vec![
        ArcVertex::straight(Point2::new(0.0, 0.0)),
        ArcVertex::straight(Point2::new(4.0, 0.0)),
        ArcVertex::bulged(Point2::new(4.0, 3.0), 0.5),
        ArcVertex::straight(Point2::new(0.0, 3.0)),
    ]);

    let forward = arc_ring_area(&ring);
    let backward = arc_ring_area(&reverse_arc_ring(&ring));
    assert!(
        (forward + backward).abs() < 1.0e-12,
        "forward {forward} backward {backward} should cancel"
    );

    // Reversing twice returns the original ring exactly.
    assert_eq!(reverse_arc_ring(&reverse_arc_ring(&ring)), ring);
}

#[test]
fn a_valid_circle_and_a_valid_polygon_both_pass() {
    let circle = ArcRing::circle(Point2::new(0.0, 0.0), 1.0);
    assert_eq!(validate_arc_ring(&circle, Tolerance::METRE), Ok(()));

    let square = ArcRing::from_points(&[
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(2.0, 2.0),
        Point2::new(0.0, 2.0),
    ]);
    assert!(square.is_polygonal());
    assert_eq!(validate_arc_ring(&square, Tolerance::METRE), Ok(()));
}

#[test]
fn a_single_vertex_ring_is_refused_as_too_short() {
    // One vertex cannot describe a circle: the chord is degenerate and the
    // centre is unrecoverable.
    let ring = ArcRing::new(vec![ArcVertex::bulged(Point2::new(0.0, 0.0), 1.0)]);
    assert_eq!(
        validate_arc_ring(&ring, Tolerance::METRE),
        Err(OverlayError::RingTooShort)
    );
}

#[test]
fn a_two_vertex_straight_ring_is_refused() {
    // Two straight edges enclose nothing. Allowed vertex counts differ by
    // curvature: two is enough ONLY when the edges are arcs.
    let ring = ArcRing::from_points(&[Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)]);
    assert_eq!(
        validate_arc_ring(&ring, Tolerance::METRE),
        Err(OverlayError::RingTooShort)
    );
}

#[test]
fn a_non_finite_bulge_is_refused() {
    let ring = ArcRing::new(vec![
        ArcVertex::bulged(Point2::new(-1.0, 0.0), f64::INFINITY),
        ArcVertex::bulged(Point2::new(1.0, 0.0), 1.0),
    ]);
    assert_eq!(
        validate_arc_ring(&ring, Tolerance::METRE),
        Err(OverlayError::NonFinitePoint)
    );
}

/// A sub-tolerance arc radius is refused by its own reason.
///
/// The radius of a bulged edge satisfies `R = chord (1 + b^2) / 4|b|`, and
/// `(1 + b^2) / 4|b|` is minimised at `b = 1` with value `1/2`. So
/// `R >= chord / 2` always, and a radius at or below tolerance `t` is
/// reachable only for a chord in `(t, 2t]` -- long enough to clear the
/// repeated-vertex test, short enough to imply an unusable arc. That
/// narrow window is exactly what this variant exists for.
#[test]
fn a_sub_tolerance_arc_radius_is_refused_by_its_own_reason() {
    let tolerance = Tolerance::METRE;
    let chord = 1.5 * tolerance.linear();
    let ring = ArcRing::new(vec![
        ArcVertex::bulged(Point2::new(0.0, 0.0), 1.0),
        ArcVertex::bulged(Point2::new(chord, 0.0), 1.0),
    ]);

    // The chord clears the repeated-vertex check ...
    assert!(chord > tolerance.linear());
    // ... while the implied radius is still not a usable boundary.
    let radius = arc_edge_radius(ring.vertices[0], ring.vertices[1]).expect("an arc");
    assert!(
        radius <= tolerance.linear(),
        "probe needs radius <= tolerance, got {radius:e}"
    );

    assert_eq!(
        validate_arc_ring(&ring, tolerance),
        Err(OverlayError::ZeroRadiusArc)
    );
}
#[test]
fn a_repeated_vertex_is_refused_before_any_radius_test() {
    let ring = ArcRing::new(vec![
        ArcVertex::bulged(Point2::new(0.0, 0.0), 1.0),
        ArcVertex::bulged(Point2::new(0.0, 0.0), 1.0),
    ]);
    assert_eq!(
        validate_arc_ring(&ring, Tolerance::METRE),
        Err(OverlayError::RepeatedVertex)
    );
}

#[test]
fn a_straight_edge_reports_no_radius() {
    let from = ArcVertex::straight(Point2::new(0.0, 0.0));
    let to = ArcVertex::straight(Point2::new(1.0, 0.0));
    assert_eq!(arc_edge_radius(from, to), None);
}

#[test]
fn a_semicircle_edge_reports_half_the_chord_as_radius() {
    // bulge 1 is a semicircle, so the radius is exactly half the chord.
    let from = ArcVertex::bulged(Point2::new(-3.0, 0.0), 1.0);
    let to = ArcVertex::bulged(Point2::new(3.0, 0.0), 1.0);
    let radius = arc_edge_radius(from, to).expect("an arc");
    assert!((radius - 3.0).abs() < 1.0e-12, "radius {radius}");
}
