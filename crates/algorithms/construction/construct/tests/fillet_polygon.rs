//! Fillet beyond rectangles: arbitrary polygon corners.
//!
//! Tangency is the property that matters: the blend must touch both walls
//! exactly, at distance `radius` from the arc centre. That is checked
//! against the geometry, not against another call of the same solver.

use axiolid_construct::feature::fillet_polygon_corner;
use axiolid_core::Point2;
use axiolid_surface::Surface;

/// A regular polygon of `sides` corners, counter-clockwise.
fn regular(sides: usize, radius: f64) -> Vec<Point2> {
    (0..sides)
        .map(|index| {
            let angle = core::f64::consts::TAU * (index as f64) / (sides as f64);
            Point2::new(radius * angle.cos(), radius * angle.sin())
        })
        .collect()
}

/// An L-shaped ring, counter-clockwise, with one reflex corner at index 3.
fn ell() -> Vec<Point2> {
    vec![
        Point2::new(0.0, 0.0),
        Point2::new(3.0, 0.0),
        Point2::new(3.0, 1.0),
        Point2::new(1.0, 1.0),
        Point2::new(1.0, 3.0),
        Point2::new(0.0, 3.0),
    ]
}

#[test]
fn a_hexagon_corner_blends_at_its_own_angle_not_a_right_angle() {
    // A hexagon's interior angle is 120 degrees, so the setback is
    // r / tan(60) = r / sqrt(3) -- the right-angle formula would have used
    // r and put the tangent points in the wrong place.
    let ring = regular(6, 2.0);
    let radius = 0.3;
    let solid = fillet_polygon_corner(&ring, 0, radius, 1.0).expect("hexagon corner fillet");

    let mut cylinders = 0;
    for surface in solid.surfaces() {
        if let Surface::Cylinder(cylinder) = surface {
            cylinders += 1;
            assert!(
                (cylinder.radius - radius).abs() < 1e-12,
                "blend radius {} should be {radius}",
                cylinder.radius
            );
        }
    }
    assert_eq!(cylinders, 1, "exactly one blend face");
}

#[test]
fn the_blend_is_tangent_to_both_walls_at_an_arbitrary_angle() {
    // Tangency: the arc centre is exactly `radius` from each adjacent
    // wall line. Checked against the input ring, independently of the
    // solver that placed the centre.
    for sides in [5_usize, 6, 8, 12] {
        let ring = regular(sides, 2.0);
        let radius = 0.2;
        let solid = fillet_polygon_corner(&ring, 0, radius, 1.0).expect("regular polygon fillet");

        let centre = solid
            .surfaces()
            .iter()
            .find_map(|surface| match surface {
                Surface::Cylinder(cylinder) => Some(cylinder.frame.origin),
                _ => None,
            })
            .expect("a blend face exists");

        // Distance from the centre to each adjacent wall line must be the
        // fillet radius, or the blend is not tangent.
        let here = ring[0];
        let previous = ring[sides - 1];
        let next = ring[1];
        for (a, b) in [(previous, here), (here, next)] {
            let along = (b - a).normalize();
            let offset = Point2::new(centre.x, centre.y) - a;
            let perpendicular = (offset - along * offset.dot(along)).length();
            assert!(
                (perpendicular - radius).abs() < 1e-9,
                "{sides}-gon: wall distance {perpendicular} should equal {radius}"
            );
        }
    }
}

#[test]
fn an_l_shape_convex_corner_is_filleted_and_its_reflex_corner_is_refused() {
    let ring = ell();
    // Index 1 is the outer right-angle corner: a normal convex fillet.
    let solid = fillet_polygon_corner(&ring, 1, 0.25, 1.0).expect("L-shape convex corner");
    assert!(solid
        .surfaces()
        .iter()
        .any(|surface| matches!(surface, Surface::Cylinder(_))));

    // Index 3 is the inner corner, where the arc would bulge into the
    // material. That is a different surface, so it must be refused.
    let error = fillet_polygon_corner(&ring, 3, 0.25, 1.0)
        .expect_err("a reflex corner is not a convex blend");
    assert!(
        format!("{error:?}").contains("reflex"),
        "refusal should name the reflex corner, got {error:?}"
    );
}

#[test]
fn a_radius_that_would_swallow_a_neighbouring_corner_is_refused() {
    // Setback at a right angle equals the radius, so 2.0 on a 3.0 edge
    // fits, but on the 1.0 edge of the L it runs past the next corner.
    let ring = ell();
    let error = fillet_polygon_corner(&ring, 2, 2.0, 1.0)
        .expect_err("an oversized radius changes the profile");
    assert!(
        format!("{error:?}").contains("larger than an adjacent edge"),
        "refusal should name the oversized radius, got {error:?}"
    );
}

#[test]
fn the_tangent_points_sit_at_the_angle_dependent_setback() {
    // The centre alone does not pin the blend: a wrong setback still
    // leaves the centre on the bisector at the right distance while the
    // arc meets the walls in the wrong place. So check the wall corners
    // the extrusion actually produced.
    let sides = 6_usize;
    let ring = regular(sides, 2.0);
    let radius = 0.3;
    // Interior angle of a regular n-gon is (n-2)*pi/n; setback is
    // r / tan(theta/2). For a hexagon that is 0.3 / tan(60) = 0.173...
    let theta = (sides as f64 - 2.0) * core::f64::consts::PI / sides as f64;
    let expected = radius / (theta / 2.0).tan();
    assert!(
        (expected - 0.3 / 3.0_f64.sqrt()).abs() < 1e-12,
        "hexagon setback closed form"
    );

    let solid = fillet_polygon_corner(&ring, 0, radius, 1.0).expect("hexagon fillet");
    // The blend's tangent points become vertices of the solid. Both lie at
    // `expected` from the original corner, along the two adjacent edges.
    let here = ring[0];
    let mut matched = 0;
    for vertex in solid.topology().vertices() {
        if vertex.position.z.abs() > 1e-12 {
            continue;
        }
        let planar = Point2::new(vertex.position.x, vertex.position.y);
        let distance = (planar - here).length();
        if (distance - expected).abs() < 1e-9 {
            matched += 1;
        }
    }
    assert_eq!(
        matched, 2,
        "both tangent points must sit at the angle-dependent setback {expected}"
    );
}
