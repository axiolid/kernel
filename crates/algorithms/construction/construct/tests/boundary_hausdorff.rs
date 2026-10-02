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
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{
    BooleanOperator, Frame2, Interval, Point2, Point3, Tolerance, Transform3, Vec2, Vec3,
};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_hausdorff::BoundaryHausdorff;
use axiolid_measure::HausdorffBounds;
use axiolid_measure::{
    boundary_hausdorff_distance, one_sided_boundary_hausdorff,
    one_sided_boundary_hausdorff_with_budget,
};
use axiolid_overlay::{ArcRing, ArcVertex};
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

#[test]
fn a_mirrored_copy_moved_is_witnessed_along_its_rims() {
    // Mirrored and then moved by `t` across and along the axis: the point
    // set is a translate, `|t|` apart, but no wall matches and no face
    // carries the translation, so no support point is seeded. The farthest
    // points are rim points, reached by splitting the rims.
    let a = column(0.3, 1.0);
    let t = Vec3::new(0.03, 0.0, 0.04);
    let b = a
        .transformed(
            &(Transform3::from_translation(t) * Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0))),
        )
        .expect("rigid");
    let accuracy = 2e-3;
    let forward = one_sided_boundary_hausdorff(&a, &b, accuracy, tol()).expect("bounded");
    contains(&forward, t.length(), accuracy);
    witnessed(&forward);
}

// Prisms trimmed in world coordinates (#227). `boolean_arc_prisms_exact`
// builds its caps on the world `xy` frame, so a translate's cap pcurves are
// moved by the translation's plan part, and its wall pcurves are rebuilt
// from moved points: no trim is identical, and every face must be matched
// as a translate.

/// Splits a translate may take: the matched bound holds every face patch
/// at `|t|` at once, so only the witnesses need placing, by halving the
/// edges near the farthest point.
const TRANSLATE_BUDGET: usize = 200;

fn square_section(half: f64) -> ArcRing {
    ArcRing::from_points(&[
        Point2::new(-half, -half),
        Point2::new(half, -half),
        Point2::new(half, half),
        Point2::new(-half, half),
    ])
}

fn round_section(radius: f64) -> ArcRing {
    ArcRing::circle(Point2::new(0.0, 0.0), radius)
}

/// A square of half-side 0.4 with its right side bowed out by a quarter
/// turn's bulge: convex, three straight walls and one round one.
fn arched_section() -> ArcRing {
    ArcRing::new(vec![
        ArcVertex::straight(Point2::new(-0.4, -0.4)),
        ArcVertex::bulged(
            Point2::new(0.4, -0.4),
            (std::f64::consts::FRAC_PI_2 / 4.0).tan(),
        ),
        ArcVertex::straight(Point2::new(0.4, 0.4)),
        ArcVertex::straight(Point2::new(-0.4, 0.4)),
    ])
}

/// The prism over `section` moved by `at`, from `at.z` to `at.z + height`,
/// as the arc prism boolean builds it: intersected with itself.
fn prism_solid(section: &ArcRing, at: Vec3, height: f64) -> ExactBRep {
    let moved = ArcRing::new(
        section
            .vertices
            .iter()
            .map(|v| ArcVertex {
                point: Point2::new(v.point.x + at.x, v.point.y + at.y),
                bulge: v.bulge,
            })
            .collect(),
    );
    let prism = ArcPrism {
        section: moved,
        bottom: at.z,
        top: at.z + height,
    };
    boolean_arc_prisms_exact(&prism, &prism, BooleanOperator::Intersection, tol()).expect("a prism")
}

/// Translations of 1 mm and 0.2 m: across the axis, along it, and oblique.
fn translations() -> Vec<Vec3> {
    let mut all = Vec::new();
    for length in [1e-3, 0.2] {
        for direction in [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 2.0, 2.0) / 3.0,
            Vec3::new(-0.48, 0.6, -0.64),
        ] {
            all.push(direction * length);
        }
    }
    all
}

#[test]
fn prisms_trimmed_in_world_coordinates_match_as_translates() {
    // Copies of one convex solid offset by `t` are `|t|` apart, both ways.
    let at = Vec3::new(3.25, -1.5, 0.4);
    for (name, section) in [
        ("square", square_section(0.5)),
        ("round", round_section(0.3)),
        ("arched", arched_section()),
    ] {
        let a = prism_solid(&section, at, 1.2);
        for t in translations() {
            let b = prism_solid(&section, at + t, 1.2);
            for accuracy in [1e-4, 1e-6] {
                for (from, to) in [(&a, &b), (&b, &a)] {
                    let bounds = one_sided_boundary_hausdorff_with_budget(
                        from,
                        to,
                        accuracy,
                        tol(),
                        TRANSLATE_BUDGET,
                    )
                    .expect("bounded");
                    assert!(
                        bounds.lower <= t.length() + 1e-12 && t.length() <= bounds.upper + 1e-12,
                        "{name} by {t:?}: [{}, {}] must contain {}",
                        bounds.lower,
                        bounds.upper,
                        t.length()
                    );
                    assert!(
                        bounds.width() <= accuracy,
                        "{name} by {t:?}: [{}, {}] did not close to {accuracy} in \
                         {TRANSLATE_BUDGET} splits",
                        bounds.lower,
                        bounds.upper
                    );
                    witnessed(&bounds);
                }
            }
        }
    }
}

#[test]
fn a_translate_closes_before_any_split() {
    // The support point against an oblique `t` lies inside a rim arc, not
    // at its ends: found there, the interval closes without refinement.
    let t = Vec3::new(1.0, 2.0, 2.0) * (1e-3 / 3.0);
    for section in [round_section(0.3), arched_section(), square_section(0.5)] {
        let a = prism_solid(&section, Vec3::ZERO, 1.2);
        let b = prism_solid(&section, t, 1.2);
        let bounds =
            one_sided_boundary_hausdorff_with_budget(&a, &b, 1e-6, tol(), 0).expect("bounded");
        contains(&bounds, t.length(), 1e-6);
    }
}

#[test]
fn a_translated_prism_closes_two_sided() {
    let at = Vec3::new(-7.0, 12.5, 1.0);
    let t = Vec3::new(0.6e-3, -0.8e-3, 0.0);
    let a = prism_solid(&round_section(0.25), at, 2.5);
    let b = prism_solid(&round_section(0.25), at + t, 2.5);
    let result = boundary_hausdorff_distance(&a, &b, 1e-6, tol()).expect("bounded");
    check_all(&result, t.length(), 1e-6);
}

#[test]
fn the_farthest_point_of_a_translated_box_is_its_corner_against_t() {
    let t = Vec3::new(1.0, 2.0, 2.0) * (1e-3 / 3.0);
    let a = prism_solid(&square_section(0.5), Vec3::ZERO, 1.2);
    let b = prism_solid(&square_section(0.5), t, 1.2);
    let bounds = one_sided_boundary_hausdorff_with_budget(&a, &b, 1e-6, tol(), TRANSLATE_BUDGET)
        .expect("bounded");
    let corner = Point3::new(-0.5, -0.5, 0.0);
    assert!(
        (bounds.point_from - corner).length() < 1e-4,
        "{:?} is not the corner extreme against t",
        bounds.point_from
    );
}

/// The interval of a pair that is not a translate: sound, whatever the
/// budget.
fn sound(from: &ExactBRep, to: &ExactBRep, expected: f64) {
    let bounds =
        one_sided_boundary_hausdorff_with_budget(from, to, 1e-7, tol(), 2_000).expect("bounded");
    assert!(
        bounds.lower <= expected + 1e-12 && expected <= bounds.upper + 1e-12,
        "[{}, {}] must contain {expected}",
        bounds.lower,
        bounds.upper
    );
    witnessed(&bounds);
}

#[test]
fn a_slightly_turned_prism_is_not_a_translate() {
    // Turned by `theta` about its axis and lifted by `lift`: a corner of
    // either sticks out `s = a (cos theta + sin theta - 1)` past a side of
    // the other and is `lift` above or below its cap, so the farthest
    // corners are `sqrt(s^2 + lift^2)` from the other boundary, both ways.
    let half = 0.5;
    let at = Vec3::new(1.0, 2.0, 0.0);
    let lift = 1e-3;
    for theta in [1e-6_f64, 1e-3] {
        let (sin, cos) = theta.sin_cos();
        let turned = ArcRing::from_points(
            &square_section(half)
                .vertices
                .iter()
                .map(|v| {
                    Point2::new(
                        cos * v.point.x - sin * v.point.y,
                        sin * v.point.x + cos * v.point.y,
                    )
                })
                .collect::<Vec<_>>(),
        );
        let a = prism_solid(&square_section(half), at, 1.2);
        let b = prism_solid(&turned, at + Vec3::new(0.0, 0.0, lift), 1.2);
        let expected = (half * (cos + sin - 1.0)).hypot(lift);
        sound(&a, &b, expected);
        sound(&b, &a, expected);
    }
}

#[test]
fn a_prism_of_another_radius_is_not_a_translate() {
    // Radii r and r + e, centres d apart across the axis, one height: the
    // wall point of either nearest the other's centre line is `d + e` from
    // the other's wall and farther from its caps; nothing strays farther.
    let (radius, grow) = (0.3, 1e-6);
    let at = Vec3::new(0.5, 0.5, 0.0);
    let t = Vec3::new(1e-3, 0.0, 0.0);
    let a = prism_solid(&round_section(radius), at, 1.2);
    let b = prism_solid(&round_section(radius + grow), at + t, 1.2);
    let expected = t.length() + grow;
    sound(&a, &b, expected);
    sound(&b, &a, expected);
}

#[test]
fn a_prism_taller_by_a_little_is_not_a_translate_of_its_walls() {
    // The second prism is moved by `t` across the axis and is `e` taller:
    // its walls are trimmed at other heights, and only its top cap is a
    // translate, by `(t, e)`. The first strays `max(|t|, e)` (its wall
    // facing away from `t`, or its top cap under the other's), the second
    // `sqrt(|t|^2 + e^2)` at its top corners.
    let at = Vec3::new(-2.0, 0.0, 0.0);
    let t = Vec3::new(1e-3, 0.0, 0.0);
    let taller = 1e-6;
    let a = prism_solid(&square_section(0.5), at, 1.2);
    let b = prism_solid(&square_section(0.5), at + t, 1.2 + taller);
    sound(&a, &b, t.length().max(taller));
    sound(&b, &a, t.length().hypot(taller));
}
