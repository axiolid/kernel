//! Certified boundary distance between exact B-reps (#125, ledger row C18).
//!
//! Each expected distance is a closed form from the inputs. The fixtures put
//! the nearest points strictly inside faces -- on a cylinder wall away from
//! its seam and rims, on a disc cap away from its rim -- so the interval can
//! only close through the certified face-domain witnesses, not through edge
//! points.

use axiolid_brep::ExactBRep;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Point2, Tolerance, Vec3};
use axiolid_measure::{boundary_clearance, boundary_distance, Clearance};
use axiolid_overlay::ArcRing;
use axiolid_profile::{CircleProfile, Profile};

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn cylinder(radius: f64, height: f64) -> ExactBRep {
    extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        }),
        Vec3::Z,
        height,
        tol(),
    )
    .expect("a cylinder")
}

/// A prism over `corners` between `bottom` and `top`.
fn block(corners: &[Point2], bottom: f64, top: f64) -> ExactBRep {
    let prism = |scale: f64| {
        let centre = corners.iter().fold(Point2::ZERO, |sum, p| sum + *p) / corners.len() as f64;
        ArcPrism {
            section: ArcRing::from_points(
                &corners
                    .iter()
                    .map(|p| centre + (*p - centre) * scale)
                    .collect::<Vec<_>>(),
            ),
            bottom,
            top,
        }
    };
    boolean_arc_prisms_exact(
        &prism(1.0),
        &prism(2.0),
        BooleanOperator::Intersection,
        tol(),
    )
    .expect("a block")
}

fn contains(bounds: &axiolid_measure::DistanceBounds, expected: f64) {
    assert!(
        bounds.lower <= expected + 1e-12 && expected <= bounds.upper + 1e-12,
        "[{}, {}] must contain {expected}",
        bounds.lower,
        bounds.upper
    );
    // The witnesses are on the boundaries, `upper` apart.
    assert!(((bounds.point_a - bounds.point_b).length() - bounds.upper).abs() < 1e-12);
}

#[test]
fn a_cylinder_wall_and_a_turned_block_close_on_the_wall_interior() {
    // Cylinder r = 1, z in [0, 4]. A 2 x 2 block turned 45 degrees, its
    // near face on the line at distance 3 from the axis along (1, 1)/sqrt2,
    // lifted to z in [1, 2]: the nearest wall points are at u = pi/4, away
    // from the seam and both rims. Distance 3 - 1 = 2.
    let n = Point2::new(1.0, 1.0) / 2.0_f64.sqrt();
    let t = Point2::new(-n.y, n.x);
    let corner = |a: f64, b: f64| n * a + t * b;
    let turned = block(
        &[
            corner(3.0, -1.0),
            corner(5.0, -1.0),
            corner(5.0, 1.0),
            corner(3.0, 1.0),
        ],
        1.0,
        2.0,
    );
    let bounds = boundary_distance(&cylinder(1.0, 4.0), &turned, 1e-9, tol()).expect("bounded");
    contains(&bounds, 2.0);
    assert!(
        bounds.upper - bounds.lower <= 1e-9,
        "interval [{}, {}] did not close",
        bounds.lower,
        bounds.upper
    );
    // The witness on the wall is not on an edge: not at a rim, not the seam.
    // A point up to ~6e-5 past z = 2 is still within 1e-9 of the distance
    // (the gap grows with the square of the overshoot), hence the slack.
    assert!(bounds.point_a.z > 1.0 - 1e-4 && bounds.point_a.z < 2.0 + 1e-4);
    assert!((bounds.point_a.x - n.x).abs() < 1e-4 && (bounds.point_a.y - n.y).abs() < 1e-4);
}

#[test]
fn a_block_over_a_disc_cap_closes_on_the_cap_interior() {
    // A 0.4 x 0.4 block hovering 0.5 over the centre of a radius-1 cap: the
    // nearest cap points are interior to the circular domain.
    let small = block(
        &[
            Point2::new(-0.2, -0.2),
            Point2::new(0.2, -0.2),
            Point2::new(0.2, 0.2),
            Point2::new(-0.2, 0.2),
        ],
        2.5,
        3.0,
    );
    let bounds = boundary_distance(&cylinder(1.0, 2.0), &small, 1e-9, tol()).expect("bounded");
    contains(&bounds, 0.5);
    assert!(bounds.upper - bounds.lower <= 1e-9);
}

#[test]
fn clearance_is_decided_only_when_the_interval_clears_the_limit() {
    let near = block(
        &[
            Point2::new(3.0, -1.0),
            Point2::new(4.0, -1.0),
            Point2::new(4.0, 1.0),
            Point2::new(3.0, 1.0),
        ],
        1.0,
        2.0,
    );
    let wall = cylinder(1.0, 4.0);
    // True distance 2 (seam side, but the answer is what matters here).
    let (_, below) = boundary_clearance(&wall, &near, 2.001, tol()).expect("decided");
    assert_eq!(below, Clearance::Below);
    let (_, above) = boundary_clearance(&wall, &near, 1.999, tol()).expect("decided");
    assert_eq!(above, Clearance::Above);
    // Exactly at the limit, no interval of positive width can clear it.
    let (bounds, at) = boundary_clearance(&wall, &near, 2.0, tol()).expect("bounded");
    assert_eq!(at, Clearance::Indeterminate, "{bounds:?}");
    contains(&bounds, 2.0);
}

#[test]
fn touching_boundaries_measure_zero() {
    let a = block(
        &[
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        0.0,
        1.0,
    );
    let b = block(
        &[
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 1.0),
            Point2::new(1.0, 1.0),
        ],
        0.0,
        1.0,
    );
    let bounds = boundary_distance(&a, &b, 1e-12, tol()).expect("bounded");
    assert_eq!(bounds.lower, 0.0);
    assert!(bounds.upper <= 1e-12, "{bounds:?}");
}

#[test]
fn two_parallel_columns_are_their_axis_gap_less_both_radii() {
    // Round columns r = 0.5 and r = 0.25 with axes 3 apart, both z in
    // [0, 3]: nearest points run along two facing generators.
    let column = |cx: f64, r: f64| {
        let ring = ArcRing::circle(Point2::new(cx, 1.0), r);
        let prism = |scale: f64| ArcPrism {
            section: ArcRing::circle(Point2::new(cx, 1.0), r * scale),
            bottom: 0.0,
            top: 3.0,
        };
        let _ = ring;
        boolean_arc_prisms_exact(
            &prism(1.0),
            &prism(2.0),
            BooleanOperator::Intersection,
            tol(),
        )
        .expect("a column")
    };
    // The nearest points form whole lines, so every slice along them is a
    // near-minimal pair and refinement converges slowly there: 1e-6 is the
    // accuracy asked of this case, not 1e-9 (see the module docs).
    let bounds =
        boundary_distance(&column(0.0, 0.5), &column(3.0, 0.25), 1e-6, tol()).expect("bounded");
    contains(&bounds, 3.0 - 0.75);
    assert!(bounds.upper - bounds.lower <= 1e-6, "{bounds:?}");
}

#[test]
fn a_torus_under_a_slab_closes_on_its_crown() {
    // A circle of radius 0.5 revolved about an axis 2 away (world y): a
    // torus whose crown is the ring y = 0.5. A slab above it, y in
    // [1.25, 1.5], spanning the whole torus in x and z: distance 0.75.
    use axiolid_construct::revolve_exact::revolve_profile_exact;
    use axiolid_core::Point3;
    let torus = revolve_profile_exact(
        &Profile::Circle(CircleProfile {
            radius: 0.5,
            thickness: None,
        }),
        Point3::new(-2.0, 0.0, 0.0),
        Vec3::Y,
        std::f64::consts::TAU,
        tol(),
    )
    .expect("a torus");
    // The slab: a 10 x 0.25 rectangle in (x, y), extruded along z through
    // the torus, which spans z in [-2.5, 2.5].
    let slab = {
        let prism = |scale: f64| ArcPrism {
            section: ArcRing::from_points(&[
                Point2::new(-2.0 - 5.0 * scale, 1.375 - 0.125 * scale),
                Point2::new(-2.0 + 5.0 * scale, 1.375 - 0.125 * scale),
                Point2::new(-2.0 + 5.0 * scale, 1.375 + 0.125 * scale),
                Point2::new(-2.0 - 5.0 * scale, 1.375 + 0.125 * scale),
            ]),
            bottom: -3.0,
            top: 3.0,
        };
        boolean_arc_prisms_exact(
            &prism(1.0),
            &prism(2.0),
            BooleanOperator::Intersection,
            tol(),
        )
        .expect("a slab")
    };
    let bounds = boundary_distance(&torus, &slab, 1e-9, tol()).expect("bounded");
    contains(&bounds, 0.75);
    assert!(bounds.upper - bounds.lower <= 1e-9, "{bounds:?}");
}

#[test]
fn a_block_over_a_hole_measures_to_the_rim_not_the_opening() {
    // A 4 x 4 plate, z in [0, 1], with a round hole of radius 1; a 0.4 x 0.4
    // block hovering 0.5 over the hole's centre. The opening is closer than
    // anything on the plate and faces the block, so a witness taken inside
    // the hole would undercut the true distance: block corner to rim.
    use axiolid_core::{Frame2, Interval, Vec2};
    use axiolid_curve::{Circle2, Curve2, Line2};
    use axiolid_profile::{Contour, ContourProfile, ProfileSegment};

    let quarter = std::f64::consts::FRAC_PI_2;
    let hole = Contour::new(
        (0..4)
            .map(|index| ProfileSegment {
                curve: Curve2::Circle(Circle2 {
                    frame: Frame2 {
                        origin: Point2::ZERO,
                        x: Vec2::X,
                        y: Vec2::Y,
                    },
                    radius: 1.0,
                }),
                domain: Interval::new(quarter * index as f64, quarter * (index + 1) as f64),
                same_sense: true,
            })
            .collect(),
    );
    let corners = [
        Point2::new(-2.0, -2.0),
        Point2::new(2.0, -2.0),
        Point2::new(2.0, 2.0),
        Point2::new(-2.0, 2.0),
    ];
    let outer = Contour::new(
        (0..4)
            .map(|i| ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: corners[i],
                    direction: corners[(i + 1) % 4] - corners[i],
                }),
                domain: Interval::UNIT,
                same_sense: true,
            })
            .collect(),
    );
    let plate = extrude_profile_exact(
        &Profile::Contour(ContourProfile {
            outer,
            holes: vec![hole],
        }),
        Vec3::Z,
        1.0,
        tol(),
    )
    .expect("a plate with a hole");
    let small = block(
        &[
            Point2::new(-0.2, -0.2),
            Point2::new(0.2, -0.2),
            Point2::new(0.2, 0.2),
            Point2::new(-0.2, 0.2),
        ],
        1.5,
        2.0,
    );
    let across = 1.0 - 0.2 * 2.0_f64.sqrt();
    let expected = (across * across + 0.25).sqrt();
    let bounds = boundary_distance(&plate, &small, 1e-6, tol()).expect("bounded");
    contains(&bounds, expected);
    assert!(bounds.upper - bounds.lower <= 1e-6, "{bounds:?}");
}

#[test]
fn a_far_block_is_found_on_a_wall_whose_centre_faces_away() {
    // The whole wall starts as one patch whose centre normal points along
    // -x, while the block sits 10 away along (1, 1)/sqrt2. The nearest wall
    // point (u = pi/4, mid-height) is inside that patch but not where its
    // centre normal points: the normal-cone test must count how far the
    // patch's normals turn, or it drops the patch holding the answer.
    let c = 10.0 / 2.0_f64.sqrt();
    let far = block(
        &[
            Point2::new(c - 0.1, c - 0.1),
            Point2::new(c + 0.1, c - 0.1),
            Point2::new(c + 0.1, c + 0.1),
            Point2::new(c - 0.1, c + 0.1),
        ],
        0.75,
        1.25,
    );
    let expected = (10.0 - 0.1 * 2.0_f64.sqrt()) - 1.0;
    let bounds = boundary_distance(&cylinder(1.0, 2.0), &far, 1e-6, tol()).expect("bounded");
    contains(&bounds, expected);
    assert!(bounds.upper - bounds.lower <= 1e-6, "{bounds:?}");
}
