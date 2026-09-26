//! Exact boolean over coaxial prisms, differentially tested (#66).
//!
//! Exact boolean previously reached only half-space-bounded difference and
//! intersection. Two prisms sharing an axis are the family where the general
//! polyhedral problem collapses to one the kernel already solves exactly: the
//! 2D boolean of their cross-sections crossed with their height intervals.
//!
//! The oracle is the cross-section area, derived from the input rectangles by
//! hand rather than recorded from a run. A wall with a rectangular opening is
//! exactly this shape, so this is the dominant building-model pattern, not a
//! toy case.

use axiolid_construct::boolean_exact::{boolean_prisms_exact, Prism};
use axiolid_core::{BooleanOperator, Point2, Tolerance};

/// An axis-aligned rectangle ring, counter-clockwise.
fn rect_ring(cx: f64, cy: f64, w: f64, h: f64) -> Vec<Point2> {
    let (hw, hh) = (w / 2.0, h / 2.0);
    vec![
        Point2::new(cx - hw, cy - hh),
        Point2::new(cx + hw, cy - hh),
        Point2::new(cx + hw, cy + hh),
        Point2::new(cx - hw, cy + hh),
    ]
}

fn prism(rings: Vec<Vec<Point2>>, bottom: f64, top: f64) -> Prism {
    Prism { rings, bottom, top }
}

/// Axis-aligned extent of the built solid's base cross-section.
///
/// Deliberately NOT a shoelace over all base vertices: with a hole present
/// those vertices belong to two separate rings, and treating them as one
/// polygon produces a meaningless number. The bounding extent is enough to
/// show the outer boundary was preserved, and the hole is checked separately
/// through the cap's loop count.
fn base_extent(brep: &axiolid_brep::ExactBRep) -> (f64, f64) {
    let points: Vec<(f64, f64)> = brep
        .topology()
        .vertices()
        .iter()
        .filter(|v| v.position.z.abs() < 1e-9)
        .map(|v| (v.position.x, v.position.y))
        .collect();
    let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
    for (x, y) in points {
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
    }
    (max_x - min_x, max_y - min_y)
}

/// Shoelace area over a simple (hole-free) base cross-section.
fn base_area(brep: &axiolid_brep::ExactBRep) -> f64 {
    let mut points: Vec<(f64, f64)> = brep
        .topology()
        .vertices()
        .iter()
        .filter(|v| v.position.z.abs() < 1e-9)
        .map(|v| (v.position.x, v.position.y))
        .collect();
    points.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9);
    (0..points.len())
        .map(|i| {
            let (x0, y0) = points[i];
            let (x1, y1) = points[(i + 1) % points.len()];
            x0 * y1 - x1 * y0
        })
        .sum::<f64>()
        .abs()
        / 2.0
}

/// The case #66 names: a general (non-half-space) exact boolean now succeeds.
///
/// Two overlapping 4x4 columns, offset by 2 in x. Their intersection is a
/// 2x4 region -- previously a typed refusal, now an exact solid.
#[test]
fn a_general_exact_intersection_now_succeeds() {
    let subject = prism(vec![rect_ring(0.0, 0.0, 4.0, 4.0)], 0.0, 3.0);
    let tool = prism(vec![rect_ring(2.0, 0.0, 4.0, 4.0)], 0.0, 3.0);

    let result = boolean_prisms_exact(
        &subject,
        &tool,
        BooleanOperator::Intersection,
        Tolerance::METRE,
    )
    .expect("coaxial prism intersection is exactly constructible");

    // Overlap spans x in [0, 2] and y in [-2, 2]: area 8, not the naive 16.
    assert!(
        (base_area(&result) - 8.0).abs() < 1e-9,
        "expected the 2x4 overlap area 8, got {}",
        base_area(&result)
    );
    assert!(
        result
            .surfaces()
            .iter()
            .all(|s| matches!(s, axiolid_surface::Surface::Plane(_))),
        "a prism boolean must stay planar-faced"
    );
}

/// A wall with an interior opening: the dominant building-model pattern.
///
/// The opening is narrower than the wall in BOTH plan directions, so it
/// leaves a hole rather than cutting the wall in two. The tool spans the full
/// height, so the result is a prism with a hole -- exactly representable.
#[test]
fn a_wall_with_an_interior_opening_differences_exactly() {
    let wall = prism(vec![rect_ring(0.0, 0.0, 10.0, 4.0)], 0.0, 3.0);
    let opening = prism(vec![rect_ring(0.0, 0.0, 2.0, 2.0)], 0.0, 3.0);

    let result = boolean_prisms_exact(
        &wall,
        &opening,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("a full-height interior opening leaves one prism with a hole");

    // The base ring is the wall outline; the hole is a second ring. Area is
    // checked on the outer boundary, which the opening does not touch.
    let (width, height) = base_extent(&result);
    assert!(
        (width - 10.0).abs() < 1e-9 && (height - 4.0).abs() < 1e-9,
        "the outer boundary must be unchanged at 10 x 4, got {width} x {height}"
    );

    // A hole means more than one loop bounds the cap face.
    let cap_bounds = result
        .topology()
        .faces()
        .iter()
        .map(|f| f.bounds.len())
        .max()
        .expect("the solid has faces");
    assert!(
        cap_bounds >= 2,
        "the opening must appear as a hole loop on the cap, got {cap_bounds}"
    );
}

/// Exact volume of a solid, from its planar faces (every face here is
/// planar: rectangle operands).
fn volume(brep: &axiolid_brep::ExactBRep) -> f64 {
    axiolid_measure::exact_properties(brep, Tolerance::METRE)
        .expect("all-planar solid is measurable")
        .signed_volume
}

/// A union of prisms with differing spans is built as the stepped solid,
/// not flattened to either span (#120).
///
/// The oracle is inclusion-exclusion over boxes, from the inputs alone:
/// a 4x4x1 slab plus a 2x4x5 tower sharing a 2x4x1 overlap.
#[test]
fn a_union_with_differing_spans_is_the_stepped_solid() {
    let short = prism(vec![rect_ring(0.0, 0.0, 4.0, 4.0)], 0.0, 1.0);
    let tall = prism(vec![rect_ring(1.0, 0.0, 2.0, 4.0)], 0.0, 5.0);

    let solid = boolean_prisms_exact(&short, &tall, BooleanOperator::Union, Tolerance::METRE)
        .expect("a stepped union is an exact solid");
    let expected = 16.0 * 1.0 + 8.0 * 5.0 - 8.0 * 1.0;
    let got = volume(&solid);
    assert!(
        (got - expected).abs() < 1e-9,
        "volume {got}, expected {expected}"
    );
    let health = axiolid_brep_audit::geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    // Not a prism: the solid reaches both heights.
    let top = solid
        .topology()
        .vertices()
        .iter()
        .map(|v| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(top, 5.0);
}

/// A tool shorter than the subject leaves a pocket: built, not refused.
#[test]
fn a_difference_with_a_short_tool_leaves_a_pocket() {
    let subject = prism(vec![rect_ring(0.0, 0.0, 10.0, 4.0)], 0.0, 3.0);
    // Stops at z = 1.5, halfway up the subject, at a corner.
    let tool = prism(vec![rect_ring(-4.0, -1.0, 2.0, 2.0)], 0.0, 1.5);

    let solid = boolean_prisms_exact(
        &subject,
        &tool,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("a partial-height cut is an exact stepped solid");
    let expected = 10.0 * 4.0 * 3.0 - 2.0 * 2.0 * 1.5;
    let got = volume(&solid);
    assert!(
        (got - expected).abs() < 1e-9,
        "volume {got}, expected {expected}"
    );
}

/// A tool buried inside the subject leaves an enclosed cavity. It used to be
/// refused because the mesh compiler dropped void shells; both now carry it
/// (#120): one solid, one void shell, the subject less the tool.
#[test]
fn a_difference_that_encloses_a_cavity_carries_it_as_a_void() {
    let subject = prism(vec![rect_ring(0.0, 0.0, 10.0, 4.0)], 0.0, 3.0);
    let buried = prism(vec![rect_ring(0.0, 0.0, 2.0, 2.0)], 1.0, 2.0);
    let solid = boolean_prisms_exact(
        &subject,
        &buried,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("a cavity is representable");
    let solids = solid.topology().solids();
    assert_eq!(solids[0].voids.len(), 1, "the cavity is a void shell");
    let volume = axiolid_measure::exact_properties(&solid, Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    assert!((volume - (120.0 - 4.0)).abs() < 1e-9, "volume {volume}");
}
