//! Visibility polygons inside regions with holes (#184). Areas are closed
//! forms from the geometry.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{Polygon, Region, Ring, VisibilityError};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| p(x, y)).collect(),
    }
}

fn tol() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).unwrap()
}

fn region(outer: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Region {
    Region::new(
        vec![Polygon {
            outer: ring(outer),
            holes: holes.iter().map(|h| ring(h)).collect(),
        }],
        tol(),
    )
    .unwrap()
}

const ROOM: [(f64, f64); 4] = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];

/// The visible part, checked to lie in the region, with its area.
fn seen(r: &Region, from: Point2) -> (Region, f64) {
    let v = r.visibility_polygon(from, tol()).unwrap();
    let outside = v.difference(r, tol()).unwrap();
    assert!(
        outside.area() < 1e-9,
        "visible part leaves the region: {outside:?}"
    );
    let area = v.area();
    (v, area)
}

#[test]
fn a_convex_room_is_seen_whole() {
    let r = region(&ROOM, &[]);
    let (_, area) = seen(&r, p(3.0, 7.0));
    assert!((area - 100.0).abs() < 1e-9, "{area}");
}

#[test]
fn an_l_shaped_room_is_seen_round_the_corner_only_so_far() {
    // From (1, 1) the whole long arm, and of the other arm the sliver
    // below the line through the reflex corner (8, 2): area 2/7.
    let l = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (8.0, 10.0),
        (8.0, 2.0),
        (0.0, 2.0),
    ];
    let r = region(&l, &[]);
    let (_, area) = seen(&r, p(1.0, 1.0));
    assert!((area - (20.0 + 2.0 / 7.0)).abs() < 1e-9, "{area}");
}

#[test]
fn a_pillar_casts_a_shadow() {
    // From (1, 5) the pillar's near corners (4, 4) and (4, 6) bound the
    // shadow: a trapezoid from x = 4 to the far wall, heights 2 and 6,
    // less the pillar. Seen: 100 - 4 - 20.
    let pillar = [(4.0, 4.0), (4.0, 6.0), (6.0, 6.0), (6.0, 4.0)];
    let r = region(&ROOM, &[&pillar]);
    let (v, area) = seen(&r, p(1.0, 5.0));
    assert!((area - 76.0).abs() < 1e-9, "{area}");
    // Behind the pillar is hidden; beside it is not.
    let hidden = Region::new(
        vec![Polygon {
            outer: ring(&[(8.0, 4.9), (8.2, 4.9), (8.2, 5.1), (8.0, 5.1)]),
            holes: Vec::new(),
        }],
        tol(),
    )
    .unwrap();
    assert!(v.intersection(&hidden, tol()).unwrap().is_empty());
}

#[test]
fn a_viewpoint_in_line_with_two_corners() {
    // From (2, 2) the pillar's corners (4, 4) and (6, 6) are on one ray.
    // The shadow runs from (4, 6) and (6, 4) out to (6, 10) and (10, 6):
    // area 24. Seen: 100 - 4 - 24.
    let pillar = [(4.0, 4.0), (4.0, 6.0), (6.0, 6.0), (6.0, 4.0)];
    let r = region(&ROOM, &[&pillar]);
    let (v, area) = seen(&r, p(2.0, 2.0));
    assert!((area - 72.0).abs() < 1e-9, "{area}");
    // The shadow's far corners are met exactly.
    let corners: Vec<Point2> = v.polygons()[0].outer.points.clone();
    for c in [p(6.0, 10.0), p(10.0, 6.0)] {
        assert!(corners.contains(&c), "{c:?} not in {corners:?}");
    }
}

#[test]
fn several_holes_and_a_notch() {
    // Two pillars one behind the other and a notch in the outer wall; the
    // result must be star-shaped from the viewpoint and inside the region.
    let outer = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 4.0),
        (7.0, 4.5),
        (10.0, 5.0),
        (10.0, 10.0),
        (0.0, 10.0),
    ];
    let a = [(3.0, 3.0), (3.0, 4.0), (4.0, 4.0), (4.0, 3.0)];
    let b = [(6.0, 6.0), (6.0, 7.5), (7.0, 7.5), (7.0, 6.0)];
    let r = region(&outer, &[&a, &b]);
    let from = p(1.5, 1.0);
    let (v, area) = seen(&r, from);
    assert!(area > 50.0 && area < r.area(), "{area}");
    // Every boundary vertex of the result is seen from the viewpoint: the
    // segment to it stays in the region (sampled finely).
    for q in &v.polygons()[0].outer.points {
        for k in 1..200 {
            let t = k as f64 / 200.0;
            let s = from + (*q - from) * t;
            let dot = Region::new(
                vec![Polygon {
                    outer: ring(&[
                        (s.x - 1e-6, s.y - 1e-6),
                        (s.x + 1e-6, s.y - 1e-6),
                        (s.x + 1e-6, s.y + 1e-6),
                        (s.x - 1e-6, s.y + 1e-6),
                    ]),
                    holes: Vec::new(),
                }],
                tol(),
            )
            .unwrap();
            let inside = r.intersection(&dot, tol()).unwrap().area();
            assert!(inside > 1e-12, "{s:?} on the way to {q:?} is outside");
        }
    }
}

#[test]
fn a_viewpoint_on_the_boundary_or_outside_is_refused() {
    let pillar = [(4.0, 4.0), (4.0, 6.0), (6.0, 6.0), (6.0, 4.0)];
    let r = region(&ROOM, &[&pillar]);
    // (5, 6), on the pillar's top edge, passes a half-open parity count as
    // inside: only the exact on-boundary test refuses it.
    for q in [
        p(0.0, 5.0),
        p(10.0, 10.0),
        p(11.0, 5.0),
        p(5.0, 5.0),
        p(4.0, 5.0),
        p(5.0, 4.0),
        p(5.0, 6.0),
    ] {
        assert_eq!(
            r.visibility_polygon(q, tol()),
            Err(VisibilityError::NotInside),
            "{q:?}"
        );
    }
}

#[test]
fn a_viewpoint_on_a_reflex_corner_is_refused() {
    let l = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (8.0, 10.0),
        (8.0, 2.0),
        (0.0, 2.0),
    ];
    let r = region(&l, &[]);
    for q in [p(8.0, 2.0), p(8.0, 5.0), p(4.0, 2.0)] {
        assert_eq!(
            r.visibility_polygon(q, tol()),
            Err(VisibilityError::NotInside),
            "{q:?}"
        );
    }
}

#[test]
fn corners_in_sight_are_kept_exactly() {
    // Decimal coordinates: v + (w - v) need not round back to w, so the
    // pillar's corners in sight must be taken as they are.
    let room = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let pillar = [(0.45, 0.55), (0.45, 0.9), (0.65, 0.9), (0.65, 0.55)];
    let r = region(&room, &[&pillar]);
    let v = r.visibility_polygon(p(0.1, 0.3), tol()).unwrap();
    let corners = &v.polygons()[0].outer.points;
    // 0.1 + (0.45 - 0.1) is not 0.45 in binary64.
    for c in [p(0.45, 0.55), p(0.45, 0.9), p(0.65, 0.55)] {
        assert!(corners.contains(&c), "{c:?} not in {corners:?}");
    }
}
