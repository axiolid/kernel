//! Certified boundary distance and Hausdorff distance between bodies of
//! several exact solids (#229).
//!
//! Every expected value is a closed form. The two-item body is a footing,
//! `[-1, 1]^2 x [0, 0.5]`, and a round column of radius 0.3 on its axis,
//! either standing clear of it (`z` from 0.6 to 3) or on it (from 0.5).
//! Moved by `t` as a whole, the union of disjoint convex items is `|t|`
//! from its copy in the Hausdorff sense: every boundary point has its
//! translate on the other boundary, and the union's support point against
//! `t` is `|t|` from all of the copy.

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Interval, Point2, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_measure::{
    body_boundary_clearance, body_boundary_distance, body_boundary_hausdorff_distance,
    boundary_clearance, boundary_distance, boundary_hausdorff_distance,
    one_sided_body_boundary_hausdorff, one_sided_body_boundary_hausdorff_with_budget,
    BodyHausdorffBounds, BodyMeasureError, BodySide, Clearance, DistanceBounds, HausdorffBounds,
};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
};
use axiolid_topology::Solid;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn moved(brep: ExactBRep, by: Vec3) -> ExactBRep {
    brep.transformed(&Transform3::from_translation(by))
        .expect("rigid")
}

/// A round column of `radius`, from `bottom` to `top`, on the `z` axis.
fn column(radius: f64, bottom: f64, top: f64) -> ExactBRep {
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
    moved(local, Vec3::new(0.0, 0.0, bottom))
}

/// A closed square around the origin, `half` from its centre.
fn square(half: f64) -> Contour {
    let corners = [
        Point2::new(-half, -half),
        Point2::new(half, -half),
        Point2::new(half, half),
        Point2::new(-half, half),
    ];
    Contour::new(
        (0..4)
            .map(|k| ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: corners[k],
                    direction: corners[(k + 1) % 4] - corners[k],
                }),
                domain: Interval::UNIT,
                same_sense: true,
            })
            .collect(),
    )
}

/// A square frame, `outer` and `hole` half-sides, from `bottom` to `top`:
/// one solid with a passage through it.
fn frame(outer: f64, hole: f64, bottom: f64, top: f64) -> ExactBRep {
    let local = extrude_profile_exact(
        &Profile::Contour(ContourProfile {
            outer: square(outer),
            holes: vec![square(hole)],
        }),
        Vec3::Z,
        top - bottom,
        tol(),
    )
    .expect("a frame");
    moved(local, Vec3::new(0.0, 0.0, bottom))
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

fn footing() -> ExactBRep {
    block(1.0, 0.0, 0.0, 0.0, 0.5)
}

/// The footing and a column standing 0.1 clear of it.
fn footing_and_column() -> Vec<ExactBRep> {
    vec![footing(), column(0.3, 0.6, 3.0)]
}

/// The footing and a column standing on it, sharing a disc of its top.
fn column_on_footing() -> Vec<ExactBRep> {
    vec![footing(), column(0.3, 0.5, 3.0)]
}

fn shifted(body: &[ExactBRep], by: Vec3) -> Vec<ExactBRep> {
    body.iter().map(|item| moved(item.clone(), by)).collect()
}

fn holds(bounds: &DistanceBounds, expected: f64, accuracy: f64) {
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
    let gap = (bounds.point_a - bounds.point_b).length();
    assert!((gap - bounds.upper).abs() <= 1e-12 * (1.0 + gap));
}

fn contains(bounds: &HausdorffBounds, expected: f64, accuracy: f64) {
    assert!(
        bounds.lower <= expected + 1e-12 && expected <= bounds.upper + 1e-12,
        "[{}, {}] must contain {expected}",
        bounds.lower,
        bounds.upper
    );
    assert!(
        bounds.width() <= accuracy,
        "[{}, {}] did not close to {accuracy}",
        bounds.lower,
        bounds.upper
    );
    // `point_to` is on the other boundary, so at least `lower` away.
    let gap = (bounds.point_from - bounds.point_to).length();
    assert!(gap + 1e-12 >= bounds.lower, "{gap} < {}", bounds.lower);
}

/// The witness lies on the item it names: a small cube about it comes
/// within a few millimetres of that item's boundary.
fn on_item(point: axiolid_core::Point3, item: &ExactBRep) {
    let probe = block(1e-3, point.x, point.y, point.z - 1e-3, point.z + 1e-3);
    let (_, verdict) = boundary_clearance(item, &probe, 3e-3, tol()).expect("bounded");
    assert_eq!(
        verdict,
        Clearance::Below,
        "{point:?} is not on the item it names"
    );
}

// Distance.

#[test]
fn a_two_item_body_is_as_near_as_its_nearest_item() {
    // A block over `[2, 3] x [-0.5, 0.5] x [2, 3]`: the column wall is 1.7
    // from its face `x = 2` along a whole segment, the footing's top edge
    // `sqrt(1 + 1.5^2)` from its bottom edge.
    let body = footing_and_column();
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    let accuracy = 1e-6;
    let found = body_boundary_distance(&body, &other, accuracy, tol()).expect("bounded");
    holds(&found.bounds, 1.7, accuracy);
    assert_eq!((found.item_a, found.item_b), (1, 0));
    on_item(found.bounds.point_a, &body[1]);
    // The same as the nearest of the pairs.
    let footing = boundary_distance(&body[0], &other[0], accuracy, tol()).expect("bounded");
    holds(&footing, (1.0f64 + 1.5 * 1.5).sqrt(), accuracy);
    // Either way round.
    let back = body_boundary_distance(&other, &body, accuracy, tol()).expect("bounded");
    holds(&back.bounds, 1.7, accuracy);
    assert_eq!((back.item_a, back.item_b), (0, 1));
}

#[test]
fn two_two_item_bodies_meet_at_their_nearest_pair() {
    // Side by side, 5 apart: footing to footing, 3 apart over a whole face.
    // The clearance stops as soon as the limit is cleared.
    let body = footing_and_column();
    let other = shifted(&body, Vec3::new(5.0, 0.0, 0.0));
    let (found, verdict) = body_boundary_clearance(&body, &other, 2.9, tol()).expect("bounded");
    assert_eq!(verdict, Clearance::Above);
    assert!(found.bounds.lower > 2.9 && found.bounds.lower <= 3.0 + 1e-12);
    let (found, verdict) = body_boundary_clearance(&body, &other, 3.1, tol()).expect("bounded");
    assert_eq!(verdict, Clearance::Below);
    assert_eq!((found.item_a, found.item_b), (0, 0));
    // Moved 2 aside and 3 up instead, the moved footing's bottom edge
    // `x = 1, z = 3` is 0.7 from the column's top rim, nearer than any
    // other pair.
    let lifted = shifted(&body, Vec3::new(2.0, 0.0, 3.0));
    let accuracy = 1e-6;
    let found = body_boundary_distance(&body, &lifted, accuracy, tol()).expect("bounded");
    holds(&found.bounds, 0.7, accuracy);
    assert_eq!((found.item_a, found.item_b), (1, 0));
}

#[test]
fn touching_and_overlapping_items_still_measure_the_nearest_pair() {
    // The minimum over pairs is the distance between the unions' boundaries
    // whatever the items' layout: a shared disc or an overlap is interior
    // to the union and never nearer to a body outside.
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    let accuracy = 1e-6;
    for body in [column_on_footing(), vec![footing(), column(0.3, 0.2, 3.0)]] {
        let found = body_boundary_distance(&body, &other, accuracy, tol()).expect("bounded");
        holds(&found.bounds, 1.7, accuracy);
        assert_eq!(found.item_a, 1);
    }
}

#[test]
fn one_item_measures_as_one_solid() {
    let a = [column(0.3, 0.0, 3.0)];
    let b = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    let accuracy = 1e-6;
    let body = body_boundary_distance(&a, &b, accuracy, tol()).expect("bounded");
    let single = boundary_distance(&a[0], &b[0], accuracy, tol()).expect("bounded");
    assert_eq!(body.bounds, single);
    assert_eq!((body.item_a, body.item_b), (0, 0));
}

// Hausdorff.

fn translations() -> Vec<Vec3> {
    vec![
        Vec3::new(1e-4, 0.0, 0.0),
        Vec3::new(0.03, -0.04, 0.0),
        Vec3::new(0.0, 0.0, 0.02),
        Vec3::new(-0.48, 0.6, -0.64) * 0.05,
    ]
}

/// Splits a translate may take, as for one solid.
const TRANSLATE_BUDGET: usize = 200;

#[test]
fn a_moved_two_item_body_is_its_translation_away_and_closes_fast() {
    let body = footing_and_column();
    for t in translations() {
        let copy = shifted(&body, t);
        for accuracy in [1e-6, 1e-9] {
            for (from, to) in [(&body, &copy), (&copy, &body)] {
                let found = one_sided_body_boundary_hausdorff_with_budget(
                    from,
                    to,
                    accuracy,
                    tol(),
                    TRANSLATE_BUDGET,
                )
                .expect("bounded");
                contains(&found.bounds, t.length(), accuracy);
            }
        }
        let both = body_boundary_hausdorff_distance(&body, &copy, 1e-9, tol()).expect("bounded");
        contains(&both.distance.bounds, t.length(), 1e-9);
    }
}

#[test]
fn one_item_moved_alone_is_seeded_from_that_item() {
    // Only the column moves, by `s`: the union's support point against `s`
    // is on the footing, which did not move; the column's own support
    // point, on its rim, is `|s|` from the moved column. Moved obliquely,
    // no face point is that far, and without that seed the rims would have
    // to be bisected down to the accuracy.
    let body = footing_and_column();
    let s = Vec3::new(0.03, 0.04, 0.02);
    let other = vec![footing(), moved(column(0.3, 0.6, 3.0), s)];
    for (from, to) in [(&body, &other), (&other, &body)] {
        let found =
            one_sided_body_boundary_hausdorff_with_budget(from, to, 1e-9, tol(), TRANSLATE_BUDGET)
                .expect("bounded");
        contains(&found.bounds, s.length(), 1e-9);
        assert_eq!((found.item_from, found.item_to), (1, 1));
        on_item(found.bounds.point_from, &from[1]);
    }
}

#[test]
fn a_two_item_body_against_one_of_its_items() {
    // The column's top cap is 2.5 above the footing's top throughout, and
    // no column point is farther from the footing; the footing's boundary
    // is part of the body's, 0 away.
    let body = footing_and_column();
    let single = [footing()];
    let accuracy = 0.01;
    let result = body_boundary_hausdorff_distance(&body, &single, accuracy, tol()).expect("ok");
    contains(&result.forward.bounds, 2.5, accuracy);
    assert_eq!(result.forward.item_from, 1);
    assert!(result.forward.bounds.point_from.z > 3.0 - 2.0 * accuracy);
    contains(&result.backward.bounds, 0.0, accuracy);
    contains(&result.distance.bounds, 2.5, accuracy);
}

#[test]
fn one_item_measures_as_one_solid_in_the_hausdorff_sense() {
    let a = column(0.3, 0.0, 3.0);
    let b = moved(a.clone(), Vec3::new(1e-4, 0.0, 2e-4));
    let body = body_boundary_hausdorff_distance(
        std::slice::from_ref(&a),
        std::slice::from_ref(&b),
        1e-6,
        tol(),
    )
    .expect("bounded");
    let single = boundary_hausdorff_distance(&a, &b, 1e-6, tol()).expect("bounded");
    assert_eq!(body.distance.bounds, single.distance);
    assert_eq!(body.forward.bounds, single.forward);
    assert_eq!(body.backward.bounds, single.backward);
}

#[test]
fn blocks_meeting_at_an_edge_are_one_union_boundary() {
    // Two unit blocks meeting along the vertical edge over `(0.5, 0.5)`:
    // every plane along a face normal holds a face of each, but the
    // diagonal one holds neither, so the contact has no area and the union
    // boundary is the two boundaries.
    let body = vec![
        block(0.5, 0.0, 0.0, 0.0, 1.0),
        block(0.5, 1.0, 1.0, 0.0, 1.0),
    ];
    let t = Vec3::new(0.02, -0.01, 0.03);
    let copy = shifted(&body, t);
    let result = body_boundary_hausdorff_distance(&body, &copy, 1e-9, tol()).expect("touching");
    contains(&result.distance.bounds, t.length(), 1e-9);
}

#[test]
fn a_block_resting_on_an_edge_shares_no_patch() {
    // A block of side 0.5 turned 45 degrees about `y`, its lowest edge on
    // the footing's top: only the footing has a face in the plane between
    // them, so they meet along a segment and the union boundary is the two
    // boundaries.
    let lift = 0.5 + 0.25 * std::f64::consts::SQRT_2;
    let resting = block(0.25, 0.0, 0.0, -0.25, 0.25)
        .transformed(
            &(Transform3::from_translation(Vec3::new(0.0, 0.0, lift))
                * Transform3::from_rotation_y(std::f64::consts::FRAC_PI_4)),
        )
        .expect("rigid");
    let body = vec![footing(), resting];
    let t = Vec3::new(-0.02, 0.01, 0.03);
    let copy = shifted(&body, t);
    let result = body_boundary_hausdorff_distance(&body, &copy, 1e-9, tol()).expect("touching");
    contains(&result.distance.bounds, t.length(), 1e-9);
}

#[test]
fn a_column_through_a_frame_is_apart_without_a_plane_between_them() {
    // No plane separates a column threaded through a frame, but their
    // boundaries are 0.3 apart and neither lies inside the other.
    let body = vec![frame(1.0, 0.5, 1.0, 1.5), column(0.2, 0.0, 3.0)];
    let t = Vec3::new(0.01, 0.02, -0.02);
    let copy = shifted(&body, t);
    let found =
        one_sided_body_boundary_hausdorff_with_budget(&body, &copy, 1e-9, tol(), TRANSLATE_BUDGET)
            .expect("apart");
    contains(&found.bounds, t.length(), 1e-9);
}

fn refused(result: Result<BodyHausdorffBounds, BodyMeasureError>) -> BodyMeasureError {
    result.expect_err("refused")
}

#[test]
fn items_sharing_a_patch_of_face_are_refused_by_name() {
    // The column's base disc is interior to the union, and so is the disc
    // of the footing's top under it: neither is on the union boundary.
    let body = column_on_footing();
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    let error = refused(one_sided_body_boundary_hausdorff(
        &body,
        &other,
        1e-3,
        tol(),
    ));
    assert_eq!(
        error,
        BodyMeasureError::ItemsShareFace {
            body: BodySide::First,
            first: 0,
            second: 1
        }
    );
    let error = body_boundary_hausdorff_distance(&other, &body, 1e-3, tol()).expect_err("refused");
    assert_eq!(
        error,
        BodyMeasureError::ItemsShareFace {
            body: BodySide::Second,
            first: 0,
            second: 1
        }
    );
    // Side by side, sharing a wall.
    let pair = vec![
        block(0.5, 0.0, 0.0, 0.0, 1.0),
        block(0.5, 1.0, 0.0, 0.0, 1.0),
    ];
    assert!(matches!(
        refused(one_sided_body_boundary_hausdorff(
            &pair,
            &other,
            1e-3,
            tol()
        )),
        BodyMeasureError::ItemsShareFace { .. }
    ));
}

#[test]
fn overlapping_items_are_refused_by_name() {
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    // The column sunk 0.3 into the footing.
    let sunk = vec![footing(), column(0.3, 0.2, 3.0)];
    // Two blocks overlapping by half.
    let crossing = vec![
        block(0.5, 0.0, 0.0, 0.0, 1.0),
        block(0.5, 0.5, 0.0, 0.0, 1.0),
    ];
    // A small block wholly inside a large one: their boundaries are apart.
    let nested = vec![
        block(1.0, 0.0, 0.0, 0.0, 1.0),
        block(0.2, 0.0, 0.0, 0.4, 0.6),
    ];
    for body in [sunk, crossing, nested] {
        assert_eq!(
            refused(one_sided_body_boundary_hausdorff(
                &other,
                &body,
                1e-3,
                tol()
            )),
            BodyMeasureError::ItemsOverlap {
                body: BodySide::Second,
                first: 0,
                second: 1
            }
        );
    }
}

/// One B-rep holding both solids, each its own.
fn two_solids(a: &ExactBRep, b: &ExactBRep) -> ExactBRep {
    let mut builder = ExactBRepBuilder::default();
    for item in [a, b] {
        let shells = builder.append(item, false);
        builder.topology_mut().add_solid(Solid {
            outer: shells[0],
            voids: Vec::new(),
        });
    }
    builder.finish().expect("two solids")
}

#[test]
fn an_item_of_several_solids_is_not_shown_apart_by_its_box() {
    // The first item is a small block inside the second plus a block far
    // off: each item has a point outside the other's box, and their
    // boundaries are apart, yet they overlap. Only for one solid per item
    // does that show them apart.
    let inner = block(0.2, 0.0, 0.0, 0.4, 0.6);
    let far = block(0.2, 5.0, 0.0, 0.4, 0.6);
    let body = vec![two_solids(&inner, &far), block(1.0, 0.0, 0.0, 0.0, 1.0)];
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    assert_eq!(
        refused(one_sided_body_boundary_hausdorff(
            &body,
            &other,
            1e-3,
            tol()
        )),
        BodyMeasureError::ItemsOverlap {
            body: BodySide::First,
            first: 0,
            second: 1
        }
    );
}

#[test]
fn an_empty_body_is_refused() {
    let body = footing_and_column();
    let empty: [ExactBRep; 0] = [];
    let first = BodyMeasureError::EmptyBody {
        body: BodySide::First,
    };
    let second = BodyMeasureError::EmptyBody {
        body: BodySide::Second,
    };
    assert_eq!(
        body_boundary_distance(&empty, &body, 1e-3, tol()).expect_err("empty"),
        first
    );
    assert_eq!(
        body_boundary_clearance(&body, &empty, 1.0, tol()).expect_err("empty"),
        second
    );
    assert_eq!(
        body_boundary_hausdorff_distance(&empty, &body, 1e-3, tol()).expect_err("empty"),
        first
    );
    assert_eq!(
        refused(one_sided_body_boundary_hausdorff(
            &body,
            &empty,
            1e-3,
            tol()
        )),
        second
    );
}
