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
    PlacedBody,
};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment,
    RectangleProfile,
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
    // Turned and tilted, their walls are on no axis plane, and no axis
    // lies in the wedge of directions the shared edge is extreme in: only
    // the diagonal plane shows the contact has no area.
    let body = vec![
        block(0.5, 0.0, 0.0, 0.0, 1.0),
        block(0.5, 1.0, 1.0, 0.0, 1.0),
    ];
    let t = Vec3::new(0.02, -0.01, 0.03);
    let turned = Transform3::from_rotation_x(0.3) * Transform3::from_rotation_z(0.6);
    for body in [body.clone(), placed(&body, &turned)] {
        let copy = shifted(&body, t);
        let result = body_boundary_hausdorff_distance(&body, &copy, 1e-9, tol()).expect("touching");
        contains(&result.distance.bounds, t.length(), 1e-9);
    }
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
    // Tilted as a whole, the footing's top is on no axis plane: the block
    // still has no face there, so the contact still has no area.
    let body = vec![footing(), resting];
    let t = Vec3::new(-0.02, 0.01, 0.03);
    for body in [
        body.clone(),
        placed(&body, &Transform3::from_rotation_x(0.3)),
    ] {
        let copy = shifted(&body, t);
        let result = body_boundary_hausdorff_distance(&body, &copy, 1e-9, tol()).expect("touching");
        contains(&result.distance.bounds, t.length(), 1e-9);
    }
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

// Face contact (#229): a patch of face two items share is interior to
// their union, cut away from both, and the union boundary measured.

/// A rectangular slab `2 hx` by `2 hy` centred on `(x, y)`, from `bottom` to
/// `top`.
fn slab(hx: f64, hy: f64, x: f64, y: f64, bottom: f64, top: f64) -> ExactBRep {
    let local = extrude_profile_exact(
        &Profile::Rectangle(RectangleProfile {
            x: 2.0 * hx,
            y: 2.0 * hy,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Vec3::Z,
        top - bottom,
        tol(),
    )
    .expect("a slab");
    moved(local, Vec3::new(x, y, bottom))
}

/// `n` steps, each `going` deep and `rise` high and 1 wide: step `k` spans
/// `x` from `k going` to `n going` and stands on step `k - 1`, which it
/// covers but for that step's tread.
fn stair(n: usize, going: f64, rise: f64) -> Vec<ExactBRep> {
    (0..n)
        .map(|k| {
            let k = k as f64;
            let from = k * going;
            let to = n as f64 * going;
            slab(
                0.5 * (to - from),
                0.5,
                0.5 * (from + to),
                0.0,
                k * rise,
                (k + 1.0) * rise,
            )
        })
        .collect()
}

/// The same stair as one solid: its side profile extruded across it.
fn stair_solid(n: usize, going: f64, rise: f64) -> ExactBRep {
    // Along the bottom, up the back, then down the treads from the top.
    let mut corners = vec![Point2::new(0.0, 0.0)];
    corners.push(Point2::new(n as f64 * going, 0.0));
    corners.push(Point2::new(n as f64 * going, n as f64 * rise));
    for k in (0..n).rev() {
        let (x, z) = (k as f64 * going, (k + 1) as f64 * rise);
        corners.push(Point2::new(x, z));
        if k > 0 {
            corners.push(Point2::new(x, k as f64 * rise));
        }
    }
    let contour = Contour::new(
        (0..corners.len())
            .map(|k| ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: corners[k],
                    direction: corners[(k + 1) % corners.len()] - corners[k],
                }),
                domain: Interval::UNIT,
                same_sense: true,
            })
            .collect(),
    );
    let local = extrude_profile_exact(
        &Profile::Contour(ContourProfile {
            outer: contour,
            holes: Vec::new(),
        }),
        Vec3::Z,
        1.0,
        tol(),
    )
    .expect("a stair");
    // The profile stood up in the plane `y = 0.5`, extruded towards `-y`.
    local
        .transformed(
            &(Transform3::from_translation(Vec3::new(0.0, 0.5, 0.0))
                * Transform3::from_rotation_x(std::f64::consts::FRAC_PI_2)),
        )
        .expect("rigid")
}

fn placed(body: &[ExactBRep], placement: &Transform3) -> Vec<ExactBRep> {
    body.iter()
        .map(|item| item.transformed(placement).expect("rigid"))
        .collect()
}

/// A translate closes within the translate budget at 1e-9, both ways.
fn closes_as_a_translate(body: &[ExactBRep], t: Vec3) {
    let copy = shifted(body, t);
    for (from, to) in [(body, copy.as_slice()), (copy.as_slice(), body)] {
        let found =
            one_sided_body_boundary_hausdorff_with_budget(from, to, 1e-9, tol(), TRANSLATE_BUDGET)
                .expect("measured");
        contains(&found.bounds, t.length(), 1e-9);
    }
}

#[test]
fn a_column_on_its_footing_is_one_union_boundary() {
    let body = column_on_footing();
    let accuracy = 0.01;
    // Against the footing alone: the column's top is 2.5 above it, and the
    // footing's disc under the column is not on the union's boundary: its
    // centre is 0.3 from the column's foot, the nearest point that is.
    let footing = [footing()];
    let result = body_boundary_hausdorff_distance(&body, &footing, accuracy, tol()).expect("cut");
    contains(&result.forward.bounds, 2.5, accuracy);
    contains(&result.backward.bounds, 0.3, accuracy);
    assert!((result.backward.bounds.point_from.z - 0.5).abs() < 1e-9);
    assert!(
        result
            .backward
            .bounds
            .point_from
            .x
            .hypot(result.backward.bounds.point_from.y)
            < 0.3
    );
    // Against the column alone: a bottom corner of the footing is farthest,
    // from the column's foot; the column's base is 0.3 from the union's
    // boundary at its centre.
    let column = [column(0.3, 0.5, 3.0)];
    let result = body_boundary_hausdorff_distance(&body, &column, accuracy, tol()).expect("cut");
    let corner = (2f64.sqrt() - 0.3).hypot(0.5);
    contains(&result.forward.bounds, corner, accuracy);
    contains(&result.backward.bounds, 0.3, accuracy);
    for t in translations() {
        closes_as_a_translate(&body, t);
    }
}

#[test]
fn blocks_sharing_a_wall_are_one_box() {
    // Two unit cubes side by side are the box `[-0.5, 1.5] x [-0.5, 0.5] x
    // [0, 1]`: the shared wall is cut away from both, and every other point
    // of each is on the box. Measured with the wall, its centre would be
    // 0.5 from the box.
    let pair = vec![
        block(0.5, 0.0, 0.0, 0.0, 1.0),
        block(0.5, 1.0, 0.0, 0.0, 1.0),
    ];
    let boxed = [slab(1.0, 0.5, 0.5, 0.0, 0.0, 1.0)];
    let accuracy = 0.05;
    let result = body_boundary_hausdorff_distance(&pair, &boxed, accuracy, tol()).expect("cut");
    contains(&result.forward.bounds, 0.0, accuracy);
    contains(&result.backward.bounds, 0.0, accuracy);
    // A shorter neighbour shares only the lower half of the wall; the
    // upper half stays on the boundary.
    let stepped = vec![
        block(0.5, 0.0, 0.0, 0.0, 1.0),
        block(0.5, 1.0, 0.0, 0.0, 0.5),
    ];
    for t in translations() {
        closes_as_a_translate(&pair, t);
        closes_as_a_translate(&stepped, t);
    }
}

#[test]
fn a_stair_of_stacked_steps_is_its_profile_extruded() {
    // Each step covers the one below but for its tread: the covered part
    // and the step's base are cut away, and what is left is the stair's
    // boundary, the same point set as the stair built as one solid.
    let (n, going, rise) = (4, 0.3, 0.2);
    let steps = stair(n, going, rise);
    let solid = [stair_solid(n, going, rise)];
    let accuracy = 0.05;
    let result = body_boundary_hausdorff_distance(&steps, &solid, accuracy, tol()).expect("cut");
    contains(&result.forward.bounds, 0.0, accuracy);
    contains(&result.backward.bounds, 0.0, accuracy);
    for t in translations() {
        closes_as_a_translate(&steps, t);
    }
}

#[test]
fn a_placed_assembly_keeps_its_contacts_exact() {
    // Turned about the vertical and moved, as an IFC placement does: the
    // contact planes stay horizontal, every height unchanged to the bit,
    // so the contacts are still exact and cut as before.
    let placement = Transform3::from_translation(Vec3::new(3.25, -1.5, 0.75))
        * Transform3::from_rotation_z(0.6);
    let body = placed(&column_on_footing(), &placement);
    let footing = placed(&[footing()], &placement);
    let accuracy = 0.01;
    let result = body_boundary_hausdorff_distance(&body, &footing, accuracy, tol()).expect("cut");
    contains(&result.forward.bounds, 2.5, accuracy);
    contains(&result.backward.bounds, 0.3, accuracy);
    let steps = placed(&stair(4, 0.3, 0.2), &placement);
    for t in translations() {
        closes_as_a_translate(&body, t);
        closes_as_a_translate(&steps, t);
    }
}

#[test]
fn a_frame_on_a_slab_leaves_the_slab_open_inside_it() {
    // A square frame, outer half-side 1 and passage half-side 0.5, 0.5
    // tall, standing on a slab 3 square and 1 thick. Only the ring under
    // the frame's material is cut from the slab's top; inside the passage
    // the top stays on the boundary. The slab point farthest from the
    // union's boundary is on the ring's diagonal, as far from the frame's
    // outer wall as from its inner corner: `(x, x)` with `1 - x =
    // sqrt(2) (x - 0.5)`, at `1 - 1 / sqrt(2)`. Were the passage cut too,
    // its centre would be 0.5 from the inner wall. The frame's top is 0.5
    // above the slab.
    let body = vec![
        slab(1.5, 1.5, 0.0, 0.0, 0.0, 1.0),
        frame(1.0, 0.5, 1.0, 1.5),
    ];
    let alone = [slab(1.5, 1.5, 0.0, 0.0, 0.0, 1.0)];
    let accuracy = 0.01;
    let result = body_boundary_hausdorff_distance(&body, &alone, accuracy, tol()).expect("cut");
    contains(&result.forward.bounds, 0.5, accuracy);
    contains(
        &result.backward.bounds,
        1.0 - std::f64::consts::FRAC_1_SQRT_2,
        accuracy,
    );
    closes_as_a_translate(&body, Vec3::new(0.03, -0.04, 0.02));
}

#[test]
fn a_square_column_and_a_round_one_on_one_footing() {
    // A column of half-side 0.2 and one of radius 0.3, each on the footing.
    // The square's sides are 0.1 inside the round column's wall, at their
    // middles, and so is the footing's top between the square and the
    // circle; the circle is nowhere farther from the square.
    let square = vec![footing(), block(0.2, 0.0, 0.0, 0.5, 3.0)];
    let round = column_on_footing();
    let accuracy = 0.01;
    let result = body_boundary_hausdorff_distance(&square, &round, accuracy, tol()).expect("cut");
    contains(&result.forward.bounds, 0.1, accuracy);
    contains(&result.backward.bounds, 0.1, accuracy);
}

#[test]
fn a_central_block_covers_what_a_corner_block_leaves_open() {
    // On one base 1 thick, a block of half-side 0.8 and height 2 in the
    // middle, or a small one at a corner. The base's top at the centre is
    // open beside the corner block and 0.8 from the open top beside the
    // central one; nothing of the corner body is farther. The central
    // block's top is 2 above the base.
    let base = || block(1.0, 0.0, 0.0, 0.0, 1.0);
    let corner = vec![base(), block(0.05, 0.9, 0.9, 1.0, 1.1)];
    let central = vec![base(), block(0.8, 0.0, 0.0, 1.0, 3.0)];
    let forward = one_sided_body_boundary_hausdorff(&corner, &central, 0.01, tol()).expect("cut");
    contains(&forward.bounds, 0.8, 0.01);
    assert_eq!(forward.item_from, 0);
    let backward = one_sided_body_boundary_hausdorff(&central, &corner, 0.05, tol()).expect("cut");
    contains(&backward.bounds, 2.0, 0.05);
}

#[test]
fn a_lift_below_rounding_is_refused_not_glued() {
    // A column 1e-13 above its footing: neither a certified gap nor exact
    // contact.
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    let lifted = vec![footing(), column(0.3, 0.5 + 1e-13, 3.0)];
    assert!(matches!(
        refused(one_sided_body_boundary_hausdorff(
            &lifted,
            &other,
            1e-3,
            tol()
        )),
        BodyMeasureError::ItemsNearlyShareFace {
            body: BodySide::First,
            first: 0,
            second: 1,
            ..
        }
    ));
}

#[test]
fn contact_off_the_axes_or_on_an_ellipse_is_refused_by_name() {
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    // Tilted item by item, the contact plane is normal to no axis: whether
    // the column's base meets the footing's top, clears it or sinks into
    // it is below rounding, and is not guessed.
    let tilted = placed(&column_on_footing(), &Transform3::from_rotation_x(0.3));
    assert_eq!(
        refused(one_sided_body_boundary_hausdorff(
            &tilted,
            &other,
            1e-3,
            tol()
        )),
        BodyMeasureError::ContactPlaneNotAxisNormal {
            body: BodySide::First,
            first: 0,
            second: 1
        }
    );
    // Turned about the vertical, two blocks sharing a wall have their walls
    // on a plane no coordinate is constant on: refused the same way.
    let walls = placed(
        &[
            block(0.5, 0.0, 0.0, 0.0, 1.0),
            block(0.5, 1.0, 0.0, 0.0, 1.0),
        ],
        &Transform3::from_rotation_z(0.6),
    );
    assert!(matches!(
        refused(one_sided_body_boundary_hausdorff(
            &walls,
            &other,
            1e-3,
            tol()
        )),
        BodyMeasureError::ContactPlaneNotAxisNormal { .. }
    ));
    // An elliptical column's base is bounded by an ellipse, which the
    // plane cannot be cut by exactly.
    let elliptical = moved(
        extrude_profile_exact(
            &Profile::Ellipse(EllipseProfile {
                semi_axis_x: 0.4,
                semi_axis_y: 0.2,
            }),
            Vec3::Z,
            2.5,
            tol(),
        )
        .expect("an elliptical column"),
        Vec3::new(0.0, 0.0, 0.5),
    );
    assert_eq!(
        refused(one_sided_body_boundary_hausdorff(
            &[footing(), elliptical],
            &other,
            1e-3,
            tol()
        )),
        BodyMeasureError::ItemsShareFace {
            body: BodySide::First,
            first: 0,
            second: 1
        }
    );
}

/// One millimetre.
fn millimetre() -> Tolerance {
    Tolerance::new(1e-3, 1e-9).expect("a tolerance")
}

#[test]
fn a_gap_within_the_tolerance_is_measured_as_a_gap() {
    // Lifted 0.3 mm off its footing, the column leaves both discs on the
    // union's boundary, at a 1 mm tolerance as at any other: the footing's
    // disc centre is 0.3 from the boundary of the body whose column stands
    // on the footing, and the other way every point is within 0.3 mm.
    let lift = 3e-4;
    let gapped = vec![footing(), column(0.3, 0.5 + lift, 3.0 + lift)];
    let standing = column_on_footing();
    let accuracy = 0.01;
    let result = body_boundary_hausdorff_distance(&gapped, &standing, accuracy, millimetre())
        .expect("apart");
    contains(&result.forward.bounds, 0.3, accuracy);
    contains(&result.backward.bounds, lift, accuracy);
    closes_as_a_translate(&gapped, Vec3::new(0.03, -0.04, 0.0));
}

#[test]
fn a_penetration_within_the_tolerance_is_refused_with_its_depth() {
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    let depth = 3e-4;
    let sunk = vec![footing(), column(0.3, 0.5 - depth, 3.0)];
    match refused(one_sided_body_boundary_hausdorff(
        &sunk,
        &other,
        1e-3,
        millimetre(),
    )) {
        BodyMeasureError::ItemsNearlyShareFace {
            body: BodySide::First,
            first: 0,
            second: 1,
            gap,
        } => assert!((gap + depth).abs() < 1e-9, "gap {gap}"),
        error => panic!("{error:?}"),
    }
    // Beyond the tolerance it is an overlap.
    assert_eq!(
        refused(one_sided_body_boundary_hausdorff(
            &sunk,
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
fn a_ten_step_stair_closes_as_fast_as_one_solid() {
    // Ten items, nine contacts; also turned and moved as a placement does,
    // where the exact cut leaves slivers along the walls that fall out
    // differently for the moved copy.
    let steps = stair(10, 0.28, 0.175);
    let placement = Transform3::from_translation(Vec3::new(3.25, -1.5, 0.75))
        * Transform3::from_rotation_z(0.6);
    let t = Vec3::new(0.012, -0.007, 0.0);
    for (name, body) in [
        ("axis-aligned", steps.clone()),
        ("placed", placed(&steps, &placement)),
    ] {
        let started = std::time::Instant::now();
        closes_as_a_translate(&body, t);
        let elapsed = started.elapsed();
        println!("ten steps, {name}, both ways at 1e-9: {elapsed:?}");
    }
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

// Placed bodies (#229): items in the body's own frame and one placement.

/// A general placement: turned about all three axes, and moved.
fn general_placement() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.25, 3.0))
        * Transform3::from_rotation_z(0.6)
        * Transform3::from_rotation_y(-0.45)
        * Transform3::from_rotation_x(0.3)
}

#[test]
fn a_placed_body_is_cut_in_its_own_frame_then_turned() {
    let placement = general_placement();
    let column_body = column_on_footing();
    let steps = stair(4, 0.3, 0.2);
    // Placed item by item, the contacts lie on planes no axis is normal
    // to, and are refused; placed as a body, they are cut first.
    let other = [block(0.5, 2.5, 0.0, 2.0, 3.0)];
    for items in [&column_body, &steps] {
        let preplaced = placed(items, &placement);
        assert!(matches!(
            refused(one_sided_body_boundary_hausdorff(
                &preplaced,
                &other,
                1e-3,
                tol()
            )),
            BodyMeasureError::ContactPlaneNotAxisNormal { .. }
        ));
        one_sided_body_boundary_hausdorff(PlacedBody::new(items, placement), &other, 0.05, tol())
            .expect("cut in the body's frame");
    }
    // Against the footing alone, placed alike: the same closed forms as
    // unplaced, and the witnesses in the world -- the footing's point
    // farthest from the union is on its top, under the column.
    let footing_alone = [footing()];
    let accuracy = 0.01;
    let result = body_boundary_hausdorff_distance(
        PlacedBody::new(&column_body, placement),
        PlacedBody::new(&footing_alone, placement),
        accuracy,
        tol(),
    )
    .expect("cut");
    contains(&result.forward.bounds, 2.5, accuracy);
    contains(&result.backward.bounds, 0.3, accuracy);
    let local = placement
        .inverse()
        .transform_point3(result.backward.bounds.point_from);
    assert!(
        (local.z - 0.5).abs() < 1e-9,
        "{local:?} is not on the footing's top"
    );
    assert!(
        local.x.hypot(local.y) < 0.3,
        "{local:?} is not under the column"
    );
    // The stair against the stair built as one solid, placed alike.
    let solid = [stair_solid(4, 0.3, 0.2)];
    let result = body_boundary_hausdorff_distance(
        PlacedBody::new(&steps, placement),
        PlacedBody::new(&solid, placement),
        0.05,
        tol(),
    )
    .expect("cut");
    contains(&result.forward.bounds, 0.0, 0.05);
    contains(&result.backward.bounds, 0.0, 0.05);
    // Against an identically built copy placed a translation away: `|t|`,
    // closed at 1e-9 within the translate budget, both ways.
    for items in [&column_body, &steps] {
        for t in translations() {
            let moved = Transform3::from_translation(t) * placement;
            for (from, to) in [
                (
                    PlacedBody::new(items, placement),
                    PlacedBody::new(items, moved),
                ),
                (
                    PlacedBody::new(items, moved),
                    PlacedBody::new(items, placement),
                ),
            ] {
                let found = one_sided_body_boundary_hausdorff_with_budget(
                    from,
                    to,
                    1e-9,
                    tol(),
                    TRANSLATE_BUDGET,
                )
                .expect("measured");
                contains(&found.bounds, t.length(), 1e-9);
            }
        }
    }
    // Distance takes a placed body as well, and needs no cut: it agrees
    // with the same items placed one by one.
    let far = Transform3::from_translation(Vec3::new(10.0, 0.0, 0.0)) * placement;
    let placed_body = body_boundary_distance(
        PlacedBody::new(&column_body, placement),
        PlacedBody::new(&column_body, far),
        1e-6,
        tol(),
    )
    .expect("bounded");
    let preplaced = body_boundary_distance(
        &placed(&column_body, &placement),
        &placed(&column_body, &far),
        1e-6,
        tol(),
    )
    .expect("bounded");
    assert!((placed_body.bounds.upper - preplaced.bounds.upper).abs() <= 2e-6);
    assert_eq!(
        (placed_body.item_a, placed_body.item_b),
        (preplaced.item_a, preplaced.item_b)
    );
}

#[test]
fn a_placement_that_is_not_rigid_is_refused() {
    let body = column_on_footing();
    let stretched = Transform3::from_scale(Vec3::new(1.0, 2.0, 1.0));
    assert!(matches!(
        body_boundary_hausdorff_distance(PlacedBody::new(&body, stretched), &body, 1e-3, tol())
            .expect_err("not rigid"),
        BodyMeasureError::Placement {
            body: BodySide::First,
            ..
        }
    ));
}

/// A prism over a triangle with corners `(-0.5, -0.3)`, `(0.5, -0.3)` and
/// `(0, 0.5)`, from `bottom` to `top`: two of its walls slanted against
/// the axes.
fn wedge(bottom: f64, top: f64) -> ExactBRep {
    let corners = [
        Point2::new(-0.5, -0.3),
        Point2::new(0.5, -0.3),
        Point2::new(0.0, 0.5),
    ];
    let contour = Contour::new(
        (0..3)
            .map(|k| ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: corners[k],
                    direction: corners[(k + 1) % 3] - corners[k],
                }),
                domain: Interval::UNIT,
                same_sense: true,
            })
            .collect(),
    );
    let local = extrude_profile_exact(
        &Profile::Contour(ContourProfile {
            outer: contour,
            holes: Vec::new(),
        }),
        Vec3::Z,
        top - bottom,
        tol(),
    )
    .expect("a wedge");
    moved(local, Vec3::new(0.0, 0.0, bottom))
}

#[test]
fn an_item_with_slanted_walls_on_an_axis_normal_contact() {
    // A triangular prism with slanted walls stands on the footing: the
    // contact plane is the footing's top, normal to `z`, whatever its
    // walls. Its top is 1 above the footing; the footing's top point
    // farthest from the union's boundary is the triangle's incentre, the
    // inradius `area / semiperimeter` from the walls' feet.
    let body = vec![footing(), wedge(0.5, 1.5)];
    let leg = 0.5f64.hypot(0.8);
    let inradius = 0.4 / (0.5 * (1.0 + 2.0 * leg));
    let accuracy = 0.01;
    for placement in [Transform3::IDENTITY, general_placement()] {
        let result = body_boundary_hausdorff_distance(
            PlacedBody::new(&body, placement),
            PlacedBody::new(&[footing()], placement),
            accuracy,
            tol(),
        )
        .expect("cut");
        contains(&result.forward.bounds, 1.0, accuracy);
        contains(&result.backward.bounds, inradius, accuracy);
    }
    closes_as_a_translate(&body, Vec3::new(0.03, -0.04, 0.02));
}
