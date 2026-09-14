//! Does `overlay` ever return a polygon carrying a hole?
//!
//! `polygon_area` subtracts hole areas. If no overlay result can carry a
//! hole, that subtraction is unreachable and the test defending it cannot
//! fail. This pins the answer with a case that produces one.

use axiolid_core::{Frame2, Point2, Tolerance, Vec2};
use axiolid_overlay::{overlay, FillRule, OverlayInput, OverlayOperation, Polygon, Ring};

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

fn ring(x0: f64, y0: f64, x1: f64, y1: f64) -> Ring {
    Ring {
        points: vec![
            Point2::new(x0, y0),
            Point2::new(x1, y0),
            Point2::new(x1, y1),
            Point2::new(x0, y1),
        ],
    }
}

fn plain(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon {
        outer: ring(x0, y0, x1, y1),
        holes: Vec::new(),
    }
}

/// Subtracting an interior square from a larger one must leave a polygon
/// with one hole -- not two disjoint outers, and not a single merged ring.
#[test]
fn subtracting_an_interior_square_returns_a_polygon_with_a_hole() {
    let tolerance = Tolerance::new(1e-9, 1e-9).expect("tolerance");
    let subject = OverlayInput {
        frame: frame(),
        polygons: vec![plain(0.0, 0.0, 4.0, 4.0)],
    };
    let clip = OverlayInput {
        frame: frame(),
        polygons: vec![plain(1.0, 1.0, 3.0, 3.0)],
    };

    let result = overlay(
        &subject,
        &clip,
        OverlayOperation::Difference,
        FillRule::NonZero,
        tolerance,
    )
    .expect("a square minus an interior square overlays");

    assert_eq!(result.polygons.len(), 1, "one enclosing region");
    assert_eq!(
        result.polygons[0].holes.len(),
        1,
        "the interior void must come back as a HOLE, not a second outer ring",
    );
    assert_eq!(result.evidence.output_holes, 1, "the report must count it");
}

/// The reachable hole must actually reach `polygon_area`: 16 - 4 = 12.
#[test]
fn polygon_area_subtracts_a_hole_the_overlay_really_produced() {
    let tolerance = Tolerance::new(1e-9, 1e-9).expect("tolerance");
    let subject = OverlayInput {
        frame: frame(),
        polygons: vec![plain(0.0, 0.0, 4.0, 4.0)],
    };
    let clip = OverlayInput {
        frame: frame(),
        polygons: vec![plain(1.0, 1.0, 3.0, 3.0)],
    };
    let result = overlay(
        &subject,
        &clip,
        OverlayOperation::Difference,
        FillRule::NonZero,
        tolerance,
    )
    .expect("overlays");

    let area = axiolid_overlay::total_area(&result.polygons);
    assert!(
        (area - 12.0).abs() < 1e-9,
        "4x4 minus 2x2 must be 12, got {area}",
    );
}
