//! Certified Hausdorff distance between exact boundaries (#224).
//!
//! Every expected value is a closed form. Two copies of one convex solid
//! offset by a translation `t` are `|t|` apart in the Hausdorff sense, on
//! the boundaries as on the solids: each boundary point has its translate
//! on the other boundary, and the point of `A` extreme against `t` is `|t|`
//! from all of `B`. A square of half-side `a` turned by `theta` about its
//! centre strays `a (cos theta + sin theta - 1)` from the original, at its
//! corners.

use std::f64::consts::FRAC_PI_2;

use axiolid_brep::ExactBRep;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_hausdorff::BoundaryHausdorff;
use axiolid_measure::HausdorffBounds;
use axiolid_measure::{boundary_hausdorff_distance, one_sided_boundary_hausdorff};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
};

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn column(radius: f64, height: f64) -> ExactBRep {
    extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        }),
        Vec3::Z,
        height,
        tol(),
    )
    .expect("a column")
}

fn block(half: f64, height: f64) -> ExactBRep {
    extrude_profile_exact(
        &Profile::Rectangle(RectangleProfile {
            x: 2.0 * half,
            y: 2.0 * half,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Vec3::Z,
        height,
        tol(),
    )
    .expect("a block")
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

/// An arched opening 1 wide, 1.5 to the springing, a semicircular arch of
/// radius 0.5 above, through a wall 0.3 thick: the profile stood upright
/// in the plane `y = 0`, extruded towards `-y`, and placed at `at`.
fn arched_opening(at: Vec3) -> ExactBRep {
    let quarter = |from: f64| ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.0, 1.5),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 0.5,
        }),
        domain: Interval::new(from, from + FRAC_PI_2),
        same_sense: true,
    };
    let profile = Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            line(Point2::new(-0.5, 0.0), Point2::new(0.5, 0.0)),
            line(Point2::new(0.5, 0.0), Point2::new(0.5, 1.5)),
            quarter(0.0),
            quarter(FRAC_PI_2),
            line(Point2::new(-0.5, 1.5), Point2::new(-0.5, 0.0)),
        ]),
        holes: Vec::new(),
    });
    let local = extrude_profile_exact(&profile, Vec3::Z, 0.3, tol()).expect("an opening");
    local
        .transformed(&(Transform3::from_translation(at) * Transform3::from_rotation_x(FRAC_PI_2)))
        .expect("rigid")
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
}

/// `point_to` is a point of the other boundary, so it is at least `lower`
/// from `point_from`.
fn witnessed(bounds: &HausdorffBounds) {
    let gap = (bounds.point_from - bounds.point_to).length();
    assert!(gap + 1e-12 >= bounds.lower, "{gap} < {}", bounds.lower);
}

fn check_all(result: &BoundaryHausdorff, expected: f64, accuracy: f64) {
    contains(&result.distance, expected, accuracy);
    contains(&result.forward, expected, accuracy);
    contains(&result.backward, expected, accuracy);
    witnessed(&result.forward);
    witnessed(&result.backward);
}

#[test]
fn two_round_columns_a_tenth_of_a_millimetre_apart() {
    // Placed by a turn and a shift, then one of them 0.1 mm further.
    let placement =
        Transform3::from_translation(Vec3::new(4.0, -2.0, 0.5)) * Transform3::from_rotation_z(0.6);
    let shift = Vec3::new(0.6e-4, 0.8e-4, 0.0);
    let a = column(0.3, 3.0).transformed(&placement).expect("rigid");
    let b = column(0.3, 3.0)
        .transformed(&(Transform3::from_translation(shift) * placement))
        .expect("rigid");
    let result = boundary_hausdorff_distance(&a, &b, 1e-6, tol()).expect("bounded");
    check_all(&result, 1e-4, 1e-6);
}

#[test]
fn columns_offset_along_and_across_the_axis() {
    // Lifted as well as moved aside: the farthest points are rim points,
    // on edges, and the witnesses reach them along the rims.
    let shift = Vec3::new(1e-4, 0.0, 2e-4);
    let a = column(0.3, 3.0);
    let b = a
        .transformed(&Transform3::from_translation(shift))
        .expect("rigid");
    let result = boundary_hausdorff_distance(&a, &b, 1e-6, tol()).expect("bounded");
    check_all(&result, shift.length(), 1e-6);
}

#[test]
fn an_arched_opening_moved_along_its_wall() {
    let a = arched_opening(Vec3::new(2.0, 1.0, 0.0));
    let b = arched_opening(Vec3::new(2.2, 1.0, 0.0));
    let result = boundary_hausdorff_distance(&a, &b, 1e-6, tol()).expect("bounded");
    check_all(&result, 0.2, 1e-6);
    // The opening's left jamb is 0.2 from all of the moved one.
    let from = result.forward.point_from;
    assert!(from.x <= 2.0 - 0.5 + 1e-6, "{from:?} is not on the left");
}

#[test]
fn an_identical_re_export_is_within_the_accuracy() {
    for accuracy in [1e-6, 1e-9] {
        let a = arched_opening(Vec3::new(10.0, -3.0, 1.0));
        let b = arched_opening(Vec3::new(10.0, -3.0, 1.0));
        let result = boundary_hausdorff_distance(&a, &b, accuracy, tol()).expect("bounded");
        assert!(result.distance.lower >= 0.0);
        assert!(
            result.distance.upper <= accuracy,
            "an identical copy measured up to {}",
            result.distance.upper
        );
        let columns =
            boundary_hausdorff_distance(&column(0.3, 3.0), &column(0.3, 3.0), accuracy, tol())
                .expect("bounded");
        assert!(columns.distance.upper <= accuracy);
    }
}

#[test]
fn a_turned_block_strays_by_its_corners() {
    // Turned about its own axis, no face of one is a translate of a face
    // of the other: the matched bound is the turn's displacement, larger
    // than the distance at the corners, so this closes through the
    // Lipschitz bound (first order: the accuracy sets the patch size along
    // the corners) and witnesses on the corner edges.
    let (half, theta) = (0.5, 0.1);
    let a = block(half, 0.25);
    let b = a
        .transformed(&Transform3::from_rotation_z(theta))
        .expect("rigid");
    let expected = half * (theta.cos() + theta.sin() - 1.0);
    let accuracy = 5e-3;
    let result = boundary_hausdorff_distance(&a, &b, accuracy, tol()).expect("bounded");
    check_all(&result, expected, accuracy);
}

#[test]
fn a_turned_and_lifted_block_strays_by_its_bottom_corners() {
    // Lifted by `lift` as well, a bottom corner of the first block sticks
    // out `s = a (cos theta + sin theta - 1)` past a side of the second and
    // `lift` below its base: it is `sqrt(s^2 + lift^2)` from the nearest
    // bottom edge, the farthest any point strays, and an isolated one. By
    // symmetry the second's top corners stray as far.
    let (half, theta, lift) = (0.5, 0.1, 0.03);
    let a = block(half, 0.25);
    let b = a
        .transformed(
            &(Transform3::from_translation(Vec3::new(0.0, 0.0, lift))
                * Transform3::from_rotation_z(theta)),
        )
        .expect("rigid");
    let s = half * (theta.cos() + theta.sin() - 1.0);
    let accuracy = 1e-3;
    let result = boundary_hausdorff_distance(&a, &b, accuracy, tol()).expect("bounded");
    check_all(&result, s.hypot(lift), accuracy);
}

#[test]
fn nested_columns_of_different_radii() {
    // Radii 0.3 and 0.5 on one base, both 2 tall. A thin wall point at
    // height z is min(0.2, z, 2 - z) from the thick boundary (its wall, or
    // the cap disks it stands over); the thin caps lie on the thick ones.
    // A thick wall point is 0.2 from the thin wall and farther from its
    // caps; a thick cap point at radius r is r - 0.3 from the thin rim.
    // So both one-sided distances are 0.2, over a whole band of the wall.
    // The walls share their trim, so they match with a bound looser than
    // a translation's; the caps do not.
    let thin = column(0.3, 2.0);
    let thick = column(0.5, 2.0);
    let accuracy = 0.02;
    let forward = one_sided_boundary_hausdorff(&thin, &thick, accuracy, tol()).expect("bounded");
    contains(&forward, 0.2, accuracy);
    witnessed(&forward);
    let backward = one_sided_boundary_hausdorff(&thick, &thin, accuracy, tol()).expect("bounded");
    contains(&backward, 0.2, accuracy);
    witnessed(&backward);
}

#[test]
fn a_short_column_and_a_tall_one_stray_differently() {
    // Radius 0.3 on one base, 1 and 2 tall. The short top cap's centre is
    // 0.3 from the tall wall, and every other short point is nearer: 0.3.
    // The tall top cap is 1 above the short one throughout: 1. The
    // two-sided distance is the larger.
    let short = column(0.3, 1.0);
    let tall = column(0.3, 2.0);
    let accuracy = 0.01;
    let result = boundary_hausdorff_distance(&short, &tall, accuracy, tol()).expect("bounded");
    contains(&result.forward, 0.3, accuracy);
    contains(&result.backward, 1.0, accuracy);
    contains(&result.distance, 1.0, accuracy);
    witnessed(&result.forward);
    witnessed(&result.backward);
    // The forward witness is on the short top cap, near its centre.
    let from = result.forward.point_from;
    assert!(
        (from.z - 1.0).abs() < 1e-9,
        "{from:?} is not on the top cap"
    );
}

#[test]
fn a_mirrored_copy_on_itself_is_zero_through_the_lipschitz_bound() {
    // Reflected in a plane through its axis, a column is the same point
    // set, but its curved wall is reparameterised: no face matches, and
    // every patch must shrink to the accuracy.
    let a = column(0.3, 1.0);
    let b = a
        .transformed(&Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0)))
        .expect("rigid");
    let accuracy = 0.05;
    let result = boundary_hausdorff_distance(&a, &b, accuracy, tol()).expect("bounded");
    contains(&result.distance, 0.0, accuracy);
}
