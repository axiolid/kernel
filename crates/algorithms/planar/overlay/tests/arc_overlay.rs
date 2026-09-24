//! Arc-aware boolean: areas against closed forms (ADR 0050, steps 2-3).
//!
//! Every expected area here is computed from geometry by hand, never
//! from another call into the code under test.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{
    arc_overlay, arc_ring_area, ArcRing, ArcVertex, OverlayError, OverlayOperation,
};

fn disc(cx: f64, cy: f64, r: f64) -> ArcRing {
    ArcRing::circle(Point2::new(cx, cy), r)
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

/// Net area: outers minus holes.
fn net_area(result: &axiolid_overlay::ArcOverlayResult) -> f64 {
    result
        .regions
        .iter()
        .map(|region| {
            let outer = arc_ring_area(&region.outer).abs();
            let holes: f64 = region
                .holes
                .iter()
                .map(|hole| arc_ring_area(hole).abs())
                .sum();
            outer - holes
        })
        .sum()
}

#[test]
fn a_wall_minus_a_round_opening_keeps_the_opening_as_a_hole() {
    // The motivating case: a 10x3 wall with a unit-radius opening.
    // The opening is fully inside, so the backend reports it as a
    // NEGATIVE loop. If the adapter read only positives the opening
    // would vanish and the wall would come back solid.
    let wall = rect(0.0, 0.0, 10.0, 3.0);
    let opening = disc(5.0, 1.5, 1.0);

    let result = arc_overlay(
        &wall,
        &opening,
        OverlayOperation::Difference,
        Tolerance::METRE,
    )
    .expect("a wall with an opening is a valid difference");

    assert_eq!(result.evidence.regions, 1, "one wall");
    assert_eq!(result.evidence.holes, 1, "the opening must survive");
    assert_eq!(result.evidence.arc_edges, 2, "the opening stays curved");

    let expected = 30.0 - std::f64::consts::PI;
    assert!(
        (net_area(&result) - expected).abs() < 1.0e-12,
        "net area {} vs {expected}",
        net_area(&result)
    );
}

#[test]
fn two_overlapping_discs_cut_a_lens_of_the_closed_form_area() {
    // Unit discs at distance 1: the lens area is
    //   2 r^2 acos(d/2r) - (d/2) sqrt(4r^2 - d^2)
    // = 2 acos(0.5) - 0.5 sqrt(3).
    let left = disc(0.0, 0.0, 1.0);
    let right = disc(1.0, 0.0, 1.0);

    let result = arc_overlay(
        &left,
        &right,
        OverlayOperation::Intersection,
        Tolerance::METRE,
    )
    .expect("overlapping discs intersect");

    let expected = 2.0 * (0.5_f64).acos() - 0.5 * 3.0_f64.sqrt();
    assert!(
        (net_area(&result) - expected).abs() < 1.0e-12,
        "lens area {} vs {expected}",
        net_area(&result)
    );
    assert_eq!(result.evidence.line_edges, 0, "a lens has no straight edge");
    assert!(result.evidence.arc_edges >= 2, "both arcs kept");
}

#[test]
fn a_curved_result_is_never_tessellated_into_segments() {
    // The guarantee that separates this path from the polygon one.
    // A tessellating backend would still produce a near-correct area,
    // so area alone cannot catch it: the edge counts can.
    let big = disc(0.0, 0.0, 2.0);
    let small = disc(0.0, 0.0, 1.0);

    let result = arc_overlay(&big, &small, OverlayOperation::Difference, Tolerance::METRE)
        .expect("an annulus is a valid difference");

    assert_eq!(
        result.evidence.line_edges, 0,
        "no straight edges in an annulus"
    );
    assert!(
        result.evidence.arc_edges <= 8,
        "an annulus needs a handful of arcs, got {} -- tessellated?",
        result.evidence.arc_edges
    );

    let expected = std::f64::consts::PI * (4.0 - 1.0);
    assert!(
        (net_area(&result) - expected).abs() < 1.0e-12,
        "annulus area {} vs {expected}",
        net_area(&result)
    );
}

#[test]
fn a_disjoint_union_keeps_both_regions() {
    // Two separate discs: dropping one would halve the area while
    // still returning a structurally valid result.
    let left = disc(0.0, 0.0, 1.0);
    let right = disc(5.0, 0.0, 1.0);

    let result = arc_overlay(&left, &right, OverlayOperation::Union, Tolerance::METRE)
        .expect("disjoint discs union");

    assert_eq!(result.evidence.regions, 2, "both discs must survive");
    let expected = 2.0 * std::f64::consts::PI;
    assert!((net_area(&result) - expected).abs() < 1.0e-12);
}

#[test]
fn a_malformed_operand_is_refused_before_any_backend_work() {
    let good = disc(0.0, 0.0, 1.0);
    let bad = ArcRing {
        vertices: vec![ArcVertex::straight(Point2::new(0.0, 0.0))],
    };

    assert_eq!(
        arc_overlay(&good, &bad, OverlayOperation::Union, Tolerance::METRE),
        Err(OverlayError::RingTooShort)
    );
}

#[test]
fn result_winding_is_outer_ccw_and_holes_cw() {
    // Consumers read containment from winding, so the sign is part of
    // the contract, not a presentation detail. net_area() takes
    // absolute values, so it cannot catch a flipped hole -- this can.
    let wall = rect(0.0, 0.0, 10.0, 3.0);
    let opening = disc(5.0, 1.5, 1.0);

    let result = arc_overlay(
        &wall,
        &opening,
        OverlayOperation::Difference,
        Tolerance::METRE,
    )
    .expect("valid difference");

    for region in &result.regions {
        assert!(
            arc_ring_area(&region.outer) > 0.0,
            "outer must be counter-clockwise, got {}",
            arc_ring_area(&region.outer)
        );
        for hole in &region.holes {
            assert!(
                arc_ring_area(hole) < 0.0,
                "hole must be clockwise, got {}",
                arc_ring_area(hole)
            );
        }
    }
}

#[test]
fn operand_winding_does_not_change_the_answer() {
    // The backend reads subtraction from winding. A caller handing in
    // a clockwise ring means the same region, so the adapter must
    // normalise instead of silently computing a different boolean.
    let wall = rect(0.0, 0.0, 10.0, 3.0);
    let opening = disc(5.0, 1.5, 1.0);
    let reversed_wall = axiolid_overlay::reverse_arc_ring(&wall);
    let reversed_opening = axiolid_overlay::reverse_arc_ring(&opening);

    let expected = 30.0 - std::f64::consts::PI;
    for (subject, clip, label) in [
        (&wall, &opening, "ccw/ccw"),
        (&reversed_wall, &opening, "cw/ccw"),
        (&wall, &reversed_opening, "ccw/cw"),
        (&reversed_wall, &reversed_opening, "cw/cw"),
    ] {
        let result = arc_overlay(
            subject,
            clip,
            OverlayOperation::Difference,
            Tolerance::METRE,
        )
        .unwrap_or_else(|error| panic!("{label} refused: {error:?}"));
        assert!(
            (net_area(&result) - expected).abs() < 1.0e-12,
            "{label}: area {} vs {expected}",
            net_area(&result)
        );
    }
}

/// The same scene drawn in millimetres and in metres, each with its own
/// one-micrometre tolerance, must give the same answer.
///
/// The backend compares positions with a fixed epsilon (1e-5 in drawing
/// units) unless told otherwise: 10 nm in millimetres but 10 um in metres.
/// Features between 1 um and 10 um then vanished in metres only.
fn scaled(ring: &ArcRing, factor: f64) -> ArcRing {
    let mut out = ring.clone();
    for vertex in &mut out.vertices {
        vertex.point = Point2::new(vertex.point.x * factor, vertex.point.y * factor);
    }
    out
}

/// Net area in square millimetres, from millimetre and metre drawings.
fn area_in_both_units(a: &ArcRing, b: &ArcRing, op: OverlayOperation) -> (f64, f64) {
    let mm = arc_overlay(a, b, op, Tolerance::MILLIMETRE).expect("mm overlay");
    let m =
        arc_overlay(&scaled(a, 1e-3), &scaled(b, 1e-3), op, Tolerance::METRE).expect("m overlay");
    (net_area(&mm), net_area(&m) * 1e6)
}

#[test]
fn a_five_micrometre_gap_survives_a_union_in_any_unit() {
    // Two unit squares (mm) 5 um apart: above the 1 um tolerance, so the
    // union keeps the gap and its area is 1 + 0.995.
    let left = rect(0.0, 0.0, 1.0, 1.0);
    let right = rect(1.005, 0.0, 2.0, 1.0);
    let (mm, m) = area_in_both_units(&left, &right, OverlayOperation::Union);
    assert!((mm - 1.995).abs() < 1e-9, "mm drawing: {mm}");
    assert!((m - 1.995).abs() < 1e-9, "m drawing lost the gap: {m}");
}

#[test]
fn a_five_micrometre_overlap_survives_an_intersection_in_any_unit() {
    // Overlap strip 5 um x 1 mm = 0.005 mm^2, above tolerance.
    let left = rect(0.0, 0.0, 1.0, 1.0);
    let right = rect(0.995, 0.0, 2.0, 1.0);
    let (mm, m) = area_in_both_units(&left, &right, OverlayOperation::Intersection);
    assert!((mm - 0.005).abs() < 1e-9, "mm drawing: {mm}");
    assert!((m - 0.005).abs() < 1e-9, "m drawing lost the overlap: {m}");
}

#[test]
fn a_disc_crossing_an_edge_by_micrometres_cuts_the_same_area_in_any_unit() {
    // Disc r = 0.5 mm whose centre sits 4 um outside a 10 mm wall's top
    // edge. The removed area is the circular segment below the edge:
    // r^2 acos(d/r) - d sqrt(r^2 - d^2), with d = 0.004.
    let wall = rect(0.0, 0.0, 10.0, 10.0);
    let disc = disc(5.0, 10.004, 0.5);
    let (r, d) = (0.5_f64, 0.004_f64);
    let segment = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let expected = 100.0 - segment;
    let (mm, m) = area_in_both_units(&wall, &disc, OverlayOperation::Difference);
    assert!(
        (mm - expected).abs() < 1e-9,
        "mm drawing: {mm} vs {expected}"
    );
    assert!((m - expected).abs() < 1e-9, "m drawing: {m} vs {expected}");
}
