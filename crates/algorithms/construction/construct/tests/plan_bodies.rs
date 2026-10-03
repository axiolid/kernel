//! Plan distance, clearance and overlap between bodies of several exact
//! solids (#237).
//!
//! Every expected value is a closed form in plan. Body `A` is a block
//! `[-0.5, 0.5]^2 x [0, 1]` and a round column of radius 0.3 on `(3, 0)`
//! from 0 to 2, in its own frame; it is placed turned about `z` (a plan
//! rotation, which plan distances do not see) and, in the tilted cases,
//! first tipped about `x`, which changes its shadow: the block's becomes
//! the rectangle `[-0.5, 0.5] x [y_min, y_max]` of its tipped corners, and
//! the column's the hull of two ellipses of `x` semi-axis 0.3 about its
//! tipped axis' ends, whose sides `x = 3 +- 0.3` run between them. The
//! other bodies stand high above, so distances in space differ.

use axiolid_brep::ExactBRep;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Point2, Point3, Tolerance, Transform3, Vec3};
use axiolid_measure::{
    body_plan_boundary_clearance, body_plan_boundary_distance, body_plan_overlap,
    boundary_clearance, BodyDistance, BodyMeasureError, BodyPlanOverlap, BodySide, Clearance,
    PlacedBody,
};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn moved(brep: ExactBRep, by: Vec3) -> ExactBRep {
    brep.transformed(&Transform3::from_translation(by))
        .expect("rigid")
}

/// A round column of `radius` on `(x, y)`, from `bottom` to `top`.
fn column(radius: f64, x: f64, y: f64, bottom: f64, top: f64) -> ExactBRep {
    let local = extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        }),
        Vec3::Z,
        top - bottom,
        tol(),
    )
    .expect("a column");
    moved(local, Vec3::new(x, y, bottom))
}

/// A block of half-side `half` centred on `(x, y)`, from `bottom` to `top`.
fn block(half: f64, x: f64, y: f64, bottom: f64, top: f64) -> ExactBRep {
    let local = extrude_profile_exact(
        &Profile::Rectangle(RectangleProfile {
            x: 2.0 * half,
            y: 2.0 * half,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Vec3::Z,
        top - bottom,
        tol(),
    )
    .expect("a block");
    moved(local, Vec3::new(x, y, bottom))
}

/// Body `A` in its own frame: the block, then the column.
fn body_a() -> Vec<ExactBRep> {
    vec![
        block(0.5, 0.0, 0.0, 0.0, 1.0),
        column(0.3, 3.0, 0.0, 0.0, 2.0),
    ]
}

/// A turn about `z` and a move: plan distances are blind to it.
fn turn() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.25, 3.0)) * Transform3::from_rotation_z(0.6)
}

/// The tip about `x`, applied before [`turn`].
fn tip() -> Transform3 {
    Transform3::from_rotation_x(0.3)
}

/// Every item placed, for checking witnesses in the world.
fn placed(items: &[ExactBRep], placement: Transform3) -> Vec<ExactBRep> {
    items
        .iter()
        .map(|item| item.transformed(&placement).expect("rigid"))
        .collect()
}

/// The interval holds `expected` and closed to `accuracy`, the witnesses'
/// projections are `upper` apart, the items named are `items`, and each
/// witness lies on the item it names.
fn holds(
    found: &BodyDistance,
    expected: f64,
    accuracy: f64,
    items: (usize, usize),
    a: &[ExactBRep],
    b: &[ExactBRep],
) {
    let bounds = &found.bounds;
    assert!(
        bounds.lower <= expected + 1e-12 && expected <= bounds.upper + 1e-12,
        "[{}, {}] must contain {expected}",
        bounds.lower,
        bounds.upper
    );
    assert!(
        bounds.upper - bounds.lower <= accuracy,
        "[{}, {}] did not close to {accuracy}",
        bounds.lower,
        bounds.upper
    );
    let d = bounds.point_a - bounds.point_b;
    assert!((d.x.hypot(d.y) - bounds.upper).abs() <= 1e-12 * (1.0 + bounds.upper));
    assert_eq!((found.item_a, found.item_b), items, "{found:?}");
    on_item(bounds.point_a, &a[found.item_a]);
    on_item(bounds.point_b, &b[found.item_b]);
}

/// The witness lies on the item it names: a small cube about it comes
/// within a few millimetres of that item's boundary.
fn on_item(point: Point3, item: &ExactBRep) {
    let probe = block(1e-3, point.x, point.y, point.z - 1e-3, point.z + 1e-3);
    let (_, verdict) = boundary_clearance(item, &probe, 3e-3, tol()).expect("bounded");
    assert_eq!(
        verdict,
        Clearance::Below,
        "{point:?} is not on the item it names"
    );
}

/// Round columns from 20 to 21, of radius `r` on `(x, y)` in `frame`
/// coordinates, built in the world.
fn in_frame(frame: Transform3, columns: &[(f64, f64, f64)]) -> Vec<ExactBRep> {
    columns
        .iter()
        .map(|&(x, y, r)| {
            let p = frame.transform_point3(Point3::new(x, y, 0.0));
            column(r, p.x, p.y, 20.0, 21.0)
        })
        .collect()
}

#[test]
fn a_body_turned_about_z_measures_its_plan_turned() {
    // A column r = 0.2 in `A`'s frame on (3, 1.5): 1.5 - 0.3 - 0.2 = 1.0
    // from the column item; from the block, sqrt(2.5^2 + 1^2) - 0.2.
    let a = body_a();
    let world_a = placed(&a, turn());
    let single = in_frame(turn(), &[(3.0, 1.5, 0.2)]);
    let accuracy = 1e-9;
    let found = body_plan_boundary_distance(PlacedBody::new(&a, turn()), &single, accuracy, tol())
        .expect("bounded");
    holds(&found, 1.0, accuracy, (1, 0), &world_a, &single);
    // Either way round.
    let back = body_plan_boundary_distance(&single, PlacedBody::new(&a, turn()), accuracy, tol())
        .expect("bounded");
    holds(&back, 1.0, accuracy, (0, 1), &single, &world_a);

    // Against a two-item body placed alike, its items in `A`'s frame: the
    // same column, and a block of half-side 0.5 on (0, -1.8), whose face
    // y = -1.3 is 0.8 from the block's y = -0.5 along a whole segment.
    let b = vec![
        column(0.2, 3.0, 1.5, 5.0, 6.0),
        block(0.5, 0.0, -1.8, 3.0, 4.0),
    ];
    let world_b = placed(&b, turn());
    let accuracy = 1e-6;
    let found = body_plan_boundary_distance(
        PlacedBody::new(&a, turn()),
        PlacedBody::new(&b, turn()),
        accuracy,
        tol(),
    )
    .expect("bounded");
    holds(&found, 0.8, accuracy, (0, 1), &world_a, &world_b);
    // The clearance stops as soon as the limit is cleared. (At the limit
    // itself the turn's rounding may tip the verdict either way.)
    for (limit, expected) in [(0.801, Clearance::Below), (0.799, Clearance::Above)] {
        let (found, verdict) = body_plan_boundary_clearance(
            PlacedBody::new(&a, turn()),
            PlacedBody::new(&b, turn()),
            limit,
            tol(),
        )
        .expect("bounded");
        assert_eq!(verdict, expected, "{limit}: {found:?}");
        assert!(found.bounds.lower <= 0.8 + 1e-12 && 0.8 <= found.bounds.upper + 1e-12);
    }
}

#[test]
fn a_tilted_body_casts_the_shadow_of_its_tilted_items() {
    // In the tipped frame the block's shadow is [-0.5, 0.5] x [y_min,
    // y_max] and the column's sides are x = 3 +- 0.3 between y = 0 and its
    // top's y.
    let tipped = |p: Point3| tip().transform_point3(p);
    let y_max = [0.0, 1.0]
        .iter()
        .flat_map(|&z| [-0.5, 0.5].map(|y| tipped(Point3::new(0.5, y, z)).y))
        .fold(f64::NEG_INFINITY, f64::max);
    let top = tipped(Point3::new(3.0, 0.0, 2.0)).y;
    assert!(top.abs() > 0.5, "the tip moves the column's top in plan");
    let a = body_a();
    let placement = turn() * tip();
    let world_a = placed(&a, placement);
    let accuracy = 1e-9;

    // A single column r = 0.2 1.2 off the block's far side, on its middle:
    // 1.0, against the tilted block's edge.
    let single = in_frame(turn(), &[(0.0, y_max + 1.2, 0.2)]);
    let found =
        body_plan_boundary_distance(PlacedBody::new(&a, placement), &single, accuracy, tol())
            .expect("bounded");
    holds(&found, 1.0, accuracy, (0, 0), &world_a, &single);

    // A two-item body: that column, and one 0.7 off the column's side,
    // half-way up its tilted axis: 0.7, against the column.
    let b = in_frame(
        turn(),
        &[(0.0, y_max + 1.2, 0.2), (3.3 + 0.2 + 0.7, 0.5 * top, 0.2)],
    );
    let found = body_plan_boundary_distance(PlacedBody::new(&a, placement), &b, accuracy, tol())
        .expect("bounded");
    holds(&found, 0.7, accuracy, (1, 1), &world_a, &b);
    let (_, verdict) =
        body_plan_boundary_clearance(PlacedBody::new(&a, placement), &b, 0.69, tol())
            .expect("bounded");
    assert_eq!(verdict, Clearance::Above);
    let (_, verdict) =
        body_plan_boundary_clearance(PlacedBody::new(&a, placement), &b, 0.71, tol())
            .expect("bounded");
    assert_eq!(verdict, Clearance::Below);
    // Untipped, the column's shadow is its disc on (3, 0): the second
    // column is sqrt(1.2^2 + (top / 2)^2) from its centre, farther.
    let flat = body_plan_boundary_distance(PlacedBody::new(&a, turn()), &b, accuracy, tol())
        .expect("bounded");
    let expected = 1.2f64.hypot(0.5 * top) - 0.3 - 0.2;
    holds(&flat, expected, accuracy, (1, 1), &placed(&a, turn()), &b);
}

#[test]
fn overlap_is_found_between_items_of_different_bodies() {
    // `A`'s block and column do not overlap in plan; `B` is a far block
    // and a block high over `A`'s column, half-side 0.2 on (3.4, 0.1).
    let a = body_a();
    let b = vec![
        block(0.5, 40.0, 40.0, 0.0, 1.0),
        block(0.2, 3.4, 0.1, 5.0, 6.0),
    ];
    let accuracy = 1e-9;
    let found = body_plan_boundary_distance(
        PlacedBody::new(&a, turn()),
        PlacedBody::new(&b, turn()),
        accuracy,
        tol(),
    )
    .expect("bounded");
    assert_eq!((found.bounds.lower, found.bounds.upper), (0.0, 0.0));
    assert_eq!((found.item_a, found.item_b), (1, 1), "{found:?}");
    match body_plan_overlap(
        PlacedBody::new(&a, turn()),
        PlacedBody::new(&b, turn()),
        tol(),
    )
    .expect("decided")
    {
        BodyPlanOverlap::Overlapping { at, item_a, item_b } => {
            assert_eq!((item_a, item_b), (1, 1));
            // In both items' shadows, in `A`'s frame.
            let local = turn()
                .inverse()
                .transform_point3(Point3::new(at.x, at.y, 0.0));
            assert!((Point2::new(local.x, local.y) - Point2::new(3.0, 0.0)).length() < 0.3);
            assert!((local.x - 3.4).abs() < 0.2 && (local.y - 0.1).abs() < 0.2);
        }
        other => panic!("{other:?}"),
    }
    // Tipped, the block's top is a tilted plane, and its shadow still
    // holds (0, 0): a block high over it overlaps item 0.
    let over = vec![
        block(0.5, 40.0, 40.0, 0.0, 1.0),
        block(0.2, 0.0, 0.0, 5.0, 6.0),
    ];
    match body_plan_overlap(
        PlacedBody::new(&a, turn() * tip()),
        PlacedBody::new(&over, turn()),
        tol(),
    )
    .expect("decided")
    {
        BodyPlanOverlap::Overlapping { item_a, item_b, .. } => assert_eq!((item_a, item_b), (0, 1)),
        other => panic!("{other:?}"),
    }
    // Moved 2 clear of it along y, apart: no item of `B` overlaps.
    let clear = vec![
        block(0.5, 40.0, 40.0, 0.0, 1.0),
        block(0.2, 3.0, 2.5, 5.0, 6.0),
    ];
    match body_plan_overlap(
        PlacedBody::new(&a, turn()),
        PlacedBody::new(&clear, turn()),
        tol(),
    )
    .expect("decided")
    {
        // 2.5 - 0.2 - 0.3.
        BodyPlanOverlap::Disjoint { gap } => assert!(gap > 0.0 && gap <= 2.0 + 1e-12),
        other => panic!("{other:?}"),
    }
}

#[test]
fn items_of_one_body_overlapping_in_plan_are_no_overlap() {
    // A footing and a column standing on it overlap in plan and touch in
    // space; against a block 1 off the footing's side, the bodies are
    // apart, and the distance is that gap.
    let footing = vec![
        block(1.0, 0.0, 0.0, 0.0, 0.5),
        column(0.3, 0.0, 0.0, 0.5, 3.0),
    ];
    let other = [block(0.5, 2.5, 0.0, 4.0, 5.0)];
    let found = body_plan_boundary_distance(
        PlacedBody::new(&footing, turn()),
        PlacedBody::new(&other, turn()),
        1e-6,
        tol(),
    )
    .expect("bounded");
    let world = placed(&footing, turn());
    holds(&found, 1.0, 1e-6, (0, 0), &world, &placed(&other, turn()));
    assert!(matches!(
        body_plan_overlap(
            PlacedBody::new(&footing, turn()),
            PlacedBody::new(&other, turn()),
            tol()
        )
        .expect("decided"),
        BodyPlanOverlap::Disjoint { .. }
    ));
}

#[test]
fn an_item_inside_another_items_footprint_measures_zero() {
    // `B`'s second item, a column r = 0.1 high over the axis of `A`'s
    // column item, lies wholly inside its shadow; its first is far away.
    let a = body_a();
    for inner in [
        column(0.1, 3.0, 0.0, 6.0, 7.0),
        block(0.1, 3.0, 0.0, 6.0, 7.0),
    ] {
        let b = vec![block(0.5, -30.0, 10.0, 0.0, 1.0), inner];
        for tilted in [false, true] {
            // Tipped, the column item's shadow still holds the foot of
            // its axis, which the tip about `x` keeps on (3, 0).
            let placement = if tilted { turn() * tip() } else { turn() };
            let found = body_plan_boundary_distance(
                PlacedBody::new(&a, placement),
                PlacedBody::new(&b, turn()),
                1e-9,
                tol(),
            )
            .expect("bounded");
            assert_eq!(found.bounds.lower, 0.0, "{found:?}");
            assert!(found.bounds.upper <= 1e-9, "{found:?}");
            assert_eq!((found.item_a, found.item_b), (1, 1), "{found:?}");
            // Either way round: the outer item on the second side.
            let back = body_plan_boundary_distance(
                PlacedBody::new(&b, turn()),
                PlacedBody::new(&a, placement),
                1e-9,
                tol(),
            )
            .expect("bounded");
            assert!(back.bounds.upper <= 1e-9, "{back:?}");
            assert_eq!((back.item_a, back.item_b), (1, 1), "{back:?}");
            match body_plan_overlap(
                PlacedBody::new(&a, placement),
                PlacedBody::new(&b, turn()),
                tol(),
            )
            .expect("decided")
            {
                BodyPlanOverlap::Overlapping { item_a, item_b, .. } => {
                    assert_eq!((item_a, item_b), (1, 1));
                }
                other => panic!("{other:?}"),
            }
        }
    }
}

#[test]
fn footprints_touching_across_items_measure_zero_without_an_overlap() {
    // `B`'s second item, [0.5, 1.5] x [-0.5, 0.5], touches `A`'s block
    // along x = 0.5 and shares no area with it.
    let a = body_a();
    let b = vec![
        block(0.5, -30.0, 10.0, 0.0, 1.0),
        block(0.5, 1.0, 0.0, 4.0, 5.0),
    ];
    let found = body_plan_boundary_distance(&a, &b, 1e-12, tol()).expect("bounded");
    assert_eq!(found.bounds.lower, 0.0, "{found:?}");
    assert!(found.bounds.upper <= 1e-12, "{found:?}");
    assert_eq!((found.item_a, found.item_b), (0, 1), "{found:?}");
    assert_eq!(
        body_plan_overlap(&a, &b, tol()).expect("decided"),
        BodyPlanOverlap::Undecided
    );
}

#[test]
fn an_empty_body_is_refused_by_name() {
    let a = body_a();
    let empty: [ExactBRep; 0] = [];
    let first = BodyMeasureError::EmptyBody {
        body: BodySide::First,
    };
    let second = BodyMeasureError::EmptyBody {
        body: BodySide::Second,
    };
    assert_eq!(
        body_plan_boundary_distance(&empty, &a, 1e-3, tol()).expect_err("empty"),
        first
    );
    assert_eq!(
        body_plan_boundary_clearance(&a, &empty, 1.0, tol()).expect_err("empty"),
        second
    );
    assert_eq!(
        body_plan_overlap(&empty, &a, tol()).expect_err("empty"),
        first
    );
    assert_eq!(
        body_plan_overlap(&a, &empty, tol()).expect_err("empty"),
        second
    );
}

#[test]
fn a_placement_that_is_not_rigid_is_refused() {
    let a = body_a();
    let stretched = Transform3::from_scale(Vec3::new(1.0, 2.0, 1.0));
    assert!(matches!(
        body_plan_overlap(&a, PlacedBody::new(&a, stretched), tol()).expect_err("not rigid"),
        BodyMeasureError::Placement {
            body: BodySide::Second,
            ..
        }
    ));
}

/// A block over part of a column's disc overlaps it in plan, alone and as
/// an item, however the column is placed. Walls crossing in plan and a
/// plan distance of zero met between two boundary points used to end the
/// search before two level faces showed the overlap.
#[test]
fn a_block_over_part_of_a_disc_overlaps_it() {
    use axiolid_measure::{plan_overlap, PlanOverlap};
    let a = body_a();
    // The disc of radius 0.3 on (3, 0) against [3.0, 3.5] x [-0.15, 0.35]
    // and against [2.9, 3.3] x [-0.1, 0.3].
    for (half, x, y) in [(0.25, 3.25, 0.1), (0.2, 3.1, 0.1)] {
        let b = vec![block(half, x, y, 5.0, 6.0)];
        for placement in [Transform3::IDENTITY, turn()] {
            let world_a = placed(&a, placement);
            let world_b = placed(&b, placement);
            assert!(matches!(
                plan_overlap(&world_a[1], &world_b[0], tol()).expect("decided"),
                PlanOverlap::Overlapping { .. }
            ));
            assert!(matches!(
                body_plan_overlap(
                    PlacedBody::new(&a, placement),
                    PlacedBody::new(&b, placement),
                    tol()
                )
                .expect("decided"),
                BodyPlanOverlap::Overlapping {
                    item_a: 1,
                    item_b: 0,
                    ..
                }
            ));
        }
    }
}
