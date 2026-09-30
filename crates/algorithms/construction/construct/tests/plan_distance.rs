//! Certified distance between plan projections of exact bodies (#217).
//!
//! Every expected value is a closed form in plan, whatever the heights: the
//! bodies are placed so that the distance in space differs, and the plan
//! distance must not see it.

use axiolid_brep::ExactBRep;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Point2, Tolerance};
use axiolid_measure::{
    boundary_distance, plan_boundary_clearance, plan_boundary_distance, plan_overlap, Clearance,
    DistanceBounds, PlanOverlap,
};
use axiolid_overlay::ArcRing;

fn tol() -> Tolerance {
    Tolerance::METRE
}

/// A prism over `section` between `bottom` and `top`, built exactly.
fn prism(section: impl Fn(f64) -> ArcRing, bottom: f64, top: f64) -> ExactBRep {
    let at = |scale: f64| ArcPrism {
        section: section(scale),
        bottom,
        top,
    };
    boolean_arc_prisms_exact(&at(1.0), &at(2.0), BooleanOperator::Intersection, tol())
        .expect("a prism")
}

fn column(centre: Point2, radius: f64, bottom: f64, top: f64) -> ExactBRep {
    prism(|scale| ArcRing::circle(centre, radius * scale), bottom, top)
}

fn block(corners: &[Point2], bottom: f64, top: f64) -> ExactBRep {
    let centre = corners.iter().fold(Point2::ZERO, |sum, p| sum + *p) / corners.len() as f64;
    prism(
        |scale| {
            ArcRing::from_points(
                &corners
                    .iter()
                    .map(|p| centre + (*p - centre) * scale)
                    .collect::<Vec<_>>(),
            )
        },
        bottom,
        top,
    )
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
    vec![
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ]
}

/// The interval holds `expected`, and the witnesses are on the
/// boundaries, their projections `upper` apart.
fn holds(bounds: &DistanceBounds, expected: f64, width: f64) {
    assert!(
        bounds.lower <= expected + 1e-12 && expected <= bounds.upper + 1e-12,
        "[{}, {}] must contain {expected}",
        bounds.lower,
        bounds.upper
    );
    assert!(bounds.upper - bounds.lower <= width, "{bounds:?}");
    let d = bounds.point_a - bounds.point_b;
    assert!((d.x.hypot(d.y) - bounds.upper).abs() < 1e-12, "{bounds:?}");
}

#[test]
fn two_columns_at_different_heights_are_their_axis_gap_less_both_radii() {
    // r = 0.2, centres 1 apart in plan: 0.6, whatever the heights. The
    // second column stands on another floor, so in space they are far.
    let a = column(Point2::new(0.0, 0.0), 0.2, 0.0, 3.0);
    let b = column(Point2::new(0.6, 0.8), 0.2, 5.0, 6.5);
    let plan = plan_boundary_distance(&a, &b, 1e-6, tol()).expect("bounded");
    holds(&plan, 0.6, 1e-6);
    let space = boundary_distance(&a, &b, 1e-3, tol()).expect("bounded");
    assert!(space.lower > 2.0, "{space:?}");
    // Apart, certainly.
    match plan_overlap(&a, &b, tol()).expect("decided") {
        PlanOverlap::Disjoint { gap } => assert!(gap > 0.0 && gap <= 0.6 + 1e-12),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_column_over_a_slab_edge_overlaps_it_in_plan() {
    // A slab [0, 4] x [0, 4], 0.3 thick; a column r = 0.2 on it, its
    // centre 0.1 inside the slab's edge x = 4: the shadows share a patch.
    let slab = block(&rect(0.0, 0.0, 4.0, 4.0), 0.0, 0.3);
    let post = column(Point2::new(3.9, 2.0), 0.2, 0.3, 3.0);
    let bounds = plan_boundary_distance(&slab, &post, 1e-9, tol()).expect("bounded");
    assert_eq!((bounds.lower, bounds.upper), (0.0, 0.0), "{bounds:?}");
    match plan_overlap(&slab, &post, tol()).expect("decided") {
        PlanOverlap::Overlapping { at } => {
            // Inside both shadows.
            assert!(
                at.x > 0.0 && at.x < 4.0 && at.y > 0.0 && at.y < 4.0,
                "{at:?}"
            );
            assert!((at - Point2::new(3.9, 2.0)).length() < 0.2, "{at:?}");
        }
        other => panic!("{other:?}"),
    }
    // Floating 1 m above the slab, the column still overlaps it in plan.
    let high = column(Point2::new(3.9, 2.0), 0.2, 1.3, 3.0);
    assert!(matches!(
        plan_overlap(&slab, &high, tol()).expect("decided"),
        PlanOverlap::Overlapping { .. }
    ));
}

#[test]
fn a_body_inside_another_footprint_measures_zero_in_plan() {
    // In space a block inside a larger one's footprint, but above it,
    // measures the gap between them; in plan its shadow lies inside the
    // other's: zero.
    let big = block(&rect(0.0, 0.0, 10.0, 10.0), 0.0, 1.0);
    let small = block(&rect(4.0, 4.0, 5.0, 5.0), 3.0, 4.0);
    let plan = plan_boundary_distance(&big, &small, 1e-9, tol()).expect("bounded");
    assert_eq!((plan.lower, plan.upper), (0.0, 0.0), "{plan:?}");
    let space = boundary_distance(&big, &small, 1e-9, tol()).expect("bounded");
    assert!(space.lower >= 2.0 - 1e-9, "{space:?}");
}

#[test]
fn a_turned_block_beside_a_cylinder_closes_in_plan() {
    // Cylinder r = 1 about the z axis; a 2 x 2 block turned 45 degrees,
    // its near face on the line at 3 from the axis along (1, 1)/sqrt 2,
    // high above the cylinder's top: plan distance 3 - 1 = 2.
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
        10.0,
        11.0,
    );
    let wall = column(Point2::ZERO, 1.0, 0.0, 4.0);
    let plan = plan_boundary_distance(&wall, &turned, 1e-9, tol()).expect("bounded");
    holds(&plan, 2.0, 1e-9);
}

#[test]
fn clearance_in_plan_is_decided_only_past_the_limit() {
    let a = column(Point2::new(0.0, 0.0), 0.2, 0.0, 3.0);
    let b = block(&rect(1.0, -0.5, 2.0, 0.5), 4.0, 5.0);
    // Plan distance 1 - 0.2 = 0.8.
    let (_, below) = plan_boundary_clearance(&a, &b, 0.801, tol()).expect("decided");
    assert_eq!(below, Clearance::Below);
    let (_, above) = plan_boundary_clearance(&a, &b, 0.799, tol()).expect("decided");
    assert_eq!(above, Clearance::Above);
    let (bounds, at) = plan_boundary_clearance(&a, &b, 0.8, tol()).expect("bounded");
    assert_eq!(at, Clearance::Indeterminate, "{bounds:?}");
}

#[test]
fn touching_footprints_measure_zero_without_an_overlap() {
    // Two blocks side by side at different heights: the shadows touch
    // along x = 1, share no area.
    let a = block(&rect(0.0, 0.0, 1.0, 1.0), 0.0, 1.0);
    let b = block(&rect(1.0, 0.0, 2.0, 1.0), 2.0, 3.0);
    let plan = plan_boundary_distance(&a, &b, 1e-12, tol()).expect("bounded");
    assert_eq!(plan.lower, 0.0);
    assert!(plan.upper <= 1e-12, "{plan:?}");
    assert!(!matches!(
        plan_overlap(&a, &b, tol()).expect("decided"),
        PlanOverlap::Overlapping { .. } | PlanOverlap::Disjoint { .. }
    ));
}
