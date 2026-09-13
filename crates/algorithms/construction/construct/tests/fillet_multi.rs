//! Multi-corner filleting (ADR 0051).
//!
//! Radii and tangency are checked against the closed forms verified in
//! the ADR, not against another call of the same code.

use axiolid_construct::feature::{fillet_polygon_corner, fillet_polygon_corners};
use axiolid_core::Point2;
use axiolid_surface::Surface;

fn square(half: f64) -> Vec<Point2> {
    vec![
        Point2::new(-half, -half),
        Point2::new(half, -half),
        Point2::new(half, half),
        Point2::new(-half, half),
    ]
}

fn cylinders(solid: &axiolid_brep::ExactBRep) -> Vec<(f64, Point2)> {
    solid
        .surfaces()
        .iter()
        .filter_map(|surface| match surface {
            Surface::Cylinder(cylinder) => Some((
                cylinder.radius,
                Point2::new(cylinder.frame.origin.x, cylinder.frame.origin.y),
            )),
            _ => None,
        })
        .collect()
}

#[test]
fn four_corners_give_four_distinct_blend_surfaces() {
    let ring = square(1.0);
    let radius = 0.25;
    let solid = fillet_polygon_corners(
        &ring,
        &[(0, radius), (1, radius), (2, radius), (3, radius)],
        1.0,
    )
    .expect("four fillets");

    let found = cylinders(&solid);
    assert_eq!(found.len(), 4, "one blend per corner, got {}", found.len());

    // Each centre sits on the interior bisector at r/sin(pi/4) from its
    // corner, which for the axis-aligned square is (1 - r) on both axes.
    let inset = 1.0 - radius;
    for (found_radius, centre) in &found {
        assert!(
            (found_radius - radius).abs() < 1e-12,
            "blend radius {found_radius} should be {radius}"
        );
        assert!(
            (centre.x.abs() - inset).abs() < 1e-12 && (centre.y.abs() - inset).abs() < 1e-12,
            "centre {centre:?} should sit at (+-{inset}, +-{inset})"
        );
    }
}

#[test]
fn each_corner_keeps_its_own_radius() {
    let ring = square(2.0);
    let solid =
        fillet_polygon_corners(&ring, &[(0, 0.3), (2, 0.9)], 1.0).expect("two different radii");

    let mut found: Vec<f64> = cylinders(&solid).into_iter().map(|(r, _)| r).collect();
    found.sort_by(f64::total_cmp);
    assert_eq!(found.len(), 2, "expected two blends, got {}", found.len());
    assert!(
        (found[0] - 0.3).abs() < 1e-12,
        "small radius, got {}",
        found[0]
    );
    assert!(
        (found[1] - 0.9).abs() < 1e-12,
        "large radius, got {}",
        found[1]
    );
}

#[test]
fn radii_that_individually_fit_but_collide_together_are_refused() {
    // Edge length 2. Each setback alone (1.1 and 0.8) is under 2, so a
    // per-corner check passes both. Their sum is 1.9 < 2, which fits;
    // raising the second to 1.0 makes the sum 2.1 and they cross.
    let ring = square(1.0);
    let fits = fillet_polygon_corners(&ring, &[(0, 1.1), (1, 0.8)], 1.0);
    assert!(fits.is_ok(), "1.1 + 0.8 = 1.9 fits on an edge of 2");

    let collides = fillet_polygon_corners(&ring, &[(0, 1.1), (1, 1.0)], 1.0);
    let error = collides.expect_err("1.1 + 1.0 = 2.1 overruns an edge of 2");
    let text = format!("{error:?}");
    assert!(
        text.contains("too large for the edge"),
        "refusal must name the shared edge, got {text}"
    );
}

#[test]
fn one_corner_matches_the_single_corner_path_exactly() {
    // The multi-corner path must not be a second implementation that
    // happens to agree; filleting one corner has to give the same solid.
    let ring = square(1.5);
    let single = fillet_polygon_corner(&ring, 2, 0.4, 1.0).expect("single");
    let multi = fillet_polygon_corners(&ring, &[(2, 0.4)], 1.0).expect("multi");

    assert_eq!(
        single.topology().faces().len(),
        multi.topology().faces().len(),
        "same face count"
    );
    assert_eq!(cylinders(&single), cylinders(&multi), "same blend geometry");
}

#[test]
fn a_duplicate_corner_is_refused_rather_than_filleted_twice() {
    let ring = square(1.0);
    let error = fillet_polygon_corners(&ring, &[(1, 0.2), (1, 0.3)], 1.0)
        .expect_err("the same corner twice is ambiguous");
    let text = format!("{error:?}");
    assert!(text.contains("more than once"), "got {text}");
}
