//! Vertical circular arcs and intrinsic (clothoid) profiles, read through
//! the evaluator, composed piecewise and inside an elevated curve, and
//! bounded for flattening (#238).
//!
//! Clothoid references are a 50-digit integration of
//! `(cos t(s), sin t(s))` with `t(s) = atan(g0) + k0 s + c s^2 / 2`, and a
//! 50-digit root of `d(s) = d`, in an independent arbitrary-precision
//! library -- never another run of the kernel's quadrature.

use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw,
    Intrinsic2, Line2,
};
use axiolid_curve_evaluate_contract::{CurveEvaluator, CurveMeasure};
use axiolid_evaluate::banked::{banked_point, banked_tangent};
use axiolid_evaluate::{
    elevated_point, elevated_tangent, elevation_chord_bound, elevation_grade, elevation_height,
    ReferenceCurveEvaluator,
};

fn close(actual: Scalar, expected: Scalar, tolerance: Scalar) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (off by {:e}, allowed {tolerance:e})",
        (actual - expected).abs()
    );
}

/// Sag-to-crest clothoid: curvature 1/4000 -> -1/3000 over 150 m of profile
/// arc length, from grade +3% at height 20.
fn sag_to_crest() -> ElevationLaw {
    ElevationLaw::intrinsic(
        20.0,
        0.03,
        CurvatureLaw::clothoid(1.0 / 4000.0, -1.0 / 3000.0, 150.0),
    )
}

#[test]
fn a_clothoid_between_grades_matches_an_independent_integration() {
    let law = sag_to_crest();
    let cases = [
        (30.0, 20.995_137_632_123_477, 0.035_758_439_286_836_62),
        (75.0, 22.680_253_208_798_398, 0.037_819_512_338_114_1),
        (120.0, 24.280_453_482_341_656, 0.031_985_297_617_067_45),
        (149.5, 25.112_821_605_565_046, 0.023_882_407_211_006_729),
    ];
    for (d, z, g) in cases {
        close(elevation_height(&law, d).unwrap(), z, 1e-13);
        close(elevation_grade(&law, d).unwrap(), g, 1e-15);
    }
    // From straight into a 2000 m sag over 100 m, starting at -2%.
    let transition =
        ElevationLaw::intrinsic(5.0, -0.02, CurvatureLaw::clothoid(0.0, 1.0 / 2000.0, 100.0));
    close(
        elevation_height(&transition, 50.0).unwrap(),
        4.104_237_913_554_467,
        1e-13,
    );
    close(
        elevation_grade(&transition, 50.0).unwrap(),
        -0.013_746_171_991_155_665,
        1e-15,
    );
    close(
        elevation_height(&transition, 99.9).unwrap(),
        3.833_204_072_353_300_3,
        1e-13,
    );
    close(
        elevation_grade(&transition, 99.9).unwrap(),
        0.004_957_516_509_706_674,
        1e-15,
    );
    // The start is the stated height and grade exactly.
    assert_eq!(elevation_height(&transition, 0.0).unwrap(), 5.0);
    assert_eq!(elevation_grade(&transition, 0.0).unwrap(), -0.02);
}

#[test]
fn a_constant_curvature_profile_is_the_closed_form_circular_arc() {
    // Two independent routes to the same circle: quadrature plus inversion,
    // and the closed form. Out to a 6.6 grade, where d(s) is far from s.
    for (grade, radius) in [(-0.02, 2000.0), (0.03, -5000.0), (0.0, 100.0)] {
        let intrinsic = ElevationLaw::intrinsic(10.0, grade, CurvatureLaw::circular(1.0 / radius));
        let arc = ElevationLaw::circular_arc(10.0, grade, radius);
        let reach = radius.abs() * 0.95;
        for step in 1..=20 {
            let d = reach * Scalar::from(step) / 20.0;
            let Some(z) = arc.height_at(d) else { continue };
            close(
                elevation_height(&intrinsic, d).unwrap(),
                z,
                1e-9 * z.abs().max(1.0),
            );
            let g = arc.grade_at(d).unwrap();
            close(
                elevation_grade(&intrinsic, d).unwrap(),
                g,
                1e-9 * g.abs().max(1.0),
            );
        }
    }
}

#[test]
fn an_intrinsic_profile_that_turns_vertical_is_refused_past_that_point() {
    // A 100 m sag from level turns vertical at d = 100.
    let law = ElevationLaw::intrinsic(0.0, 0.0, CurvatureLaw::circular(0.01));
    // sin t = 0.99, a grade of 7.0: still readable, and the circle's.
    close(
        elevation_height(&law, 99.0).unwrap(),
        100.0 * (1.0 - (1.0 - 0.99_f64 * 0.99).sqrt()),
        1e-9,
    );
    for d in [100.5, 150.0, 1000.0] {
        assert!(elevation_height(&law, d).is_err(), "height at {d}");
        assert!(elevation_grade(&law, d).is_err(), "grade at {d}");
    }
    // Before its start and at non-finite distances: refused too.
    // Refused by name, not by an inversion that happens to fail.
    assert!(matches!(
        elevation_height(&law, -1.0),
        Err(GeomError::InvalidInput(_))
    ));
    assert!(elevation_height(&law, Scalar::NAN).is_err());
    // A closed-form arc past its domain refuses through the same route.
    let arc = ElevationLaw::circular_arc(0.0, 0.0, 100.0);
    assert!(elevation_height(&arc, 100.5).is_err());
    assert!(elevation_grade(&arc, 100.5).is_err());
}

#[test]
fn a_profile_that_curls_back_is_not_read_on_its_returning_branch() {
    // Turn 1.8 rad -- past vertical -- then turn back: d(s) rises to 100,
    // falls, and rises again past 100 much later. Plan distances beyond the
    // first vertical have a height on that later branch, but the profile is
    // no longer a function of d there; it must refuse, not pick a branch.
    let law = ElevationLaw::intrinsic(
        0.0,
        0.0,
        CurvatureLaw::piecewise(
            vec![180.0],
            vec![CurvatureLaw::circular(0.01), CurvatureLaw::circular(-0.01)],
        ),
    );
    // Before the vertical it is the circle: sin t = d / 100.
    close(
        elevation_height(&law, 50.0).unwrap(),
        100.0 * (1.0 - (0.75_f64).sqrt()),
        1e-11,
    );
    for d in [100.5, 120.0] {
        assert!(elevation_height(&law, d).is_err(), "height at {d}");
    }
}

/// Arc, clothoid, arc and constant grade, each piece starting at the height
/// and grade where the previous one ends: a continuous rail profile.
fn composed() -> ElevationLaw {
    let arc = ElevationLaw::circular_arc(100.0, -0.02, 3000.0);
    let (z1, g1) = (arc.height_at(90.0).unwrap(), arc.grade_at(90.0).unwrap());
    let spiral = ElevationLaw::intrinsic(
        z1,
        g1,
        CurvatureLaw::clothoid(1.0 / 3000.0, -1.0 / 4000.0, 120.0),
    );
    let (z2, g2) = (
        elevation_height(&spiral, 110.0).unwrap(),
        elevation_grade(&spiral, 110.0).unwrap(),
    );
    let crest = ElevationLaw::circular_arc(z2, g2, -4000.0);
    let (z3, g3) = (
        crest.height_at(60.0).unwrap(),
        crest.grade_at(60.0).unwrap(),
    );
    ElevationLaw::Piecewise {
        breaks: vec![90.0, 200.0, 260.0],
        laws: vec![arc, spiral, crest, ElevationLaw::constant_grade(z3, g3)],
    }
}

#[test]
fn seams_of_a_composed_profile_stay_continuous() {
    let law = composed();
    for seam in [90.0, 200.0, 260.0] {
        let below = seam - 1e-6;
        let (z_left, g_left) = (
            elevation_height(&law, below).unwrap(),
            elevation_grade(&law, below).unwrap(),
        );
        let (z_at, g_at) = (
            elevation_height(&law, seam).unwrap(),
            elevation_grade(&law, seam).unwrap(),
        );
        // Over 1e-6 the height moves by grade * 1e-6 and the grade by
        // well under 1e-9.
        close(z_at - z_left, g_at * 1e-6, 1e-11);
        close(g_at, g_left, 1e-9);
    }
    // The spiral piece is read in its own distance.
    let spiral_start = elevation_height(&law, 90.0).unwrap();
    close(
        spiral_start,
        ElevationLaw::circular_arc(100.0, -0.02, 3000.0)
            .height_at(90.0)
            .unwrap(),
        1e-12,
    );
    // The closed-form pieces still answer height_at; the spiral piece is
    // left to the evaluator.
    assert!(law.height_at(50.0).is_some());
    assert_eq!(law.height_at(150.0), None);
    assert!(law.height_at(230.0).is_some());
}

fn line_plan() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::new(10.0, 20.0),
        direction: Vec2::new(3.0, 4.0),
    })
}

#[test]
fn an_elevated_curve_reads_both_new_laws_through_the_evaluator() {
    let curve = Elevated3::new(line_plan(), composed());
    let evaluator = ReferenceCurveEvaluator::new();
    let wrapped = Curve3::Elevated(curve.clone());
    for d in [0.0, 45.0, 90.0, 150.0, 230.0, 300.0] {
        let z = elevation_height(&curve.elevation, d).unwrap();
        let g = elevation_grade(&curve.elevation, d).unwrap();
        let want = Point3::new(10.0 + 0.6 * d, 20.0 + 0.8 * d, z);
        let point = elevated_point(&curve, d).unwrap();
        assert!((point - want).length() < 1e-12, "{point:?} vs {want:?}");
        let via_contract = evaluator
            .point_at(&wrapped, CurveMeasure::Distance(d))
            .unwrap();
        assert!((via_contract - want).length() < 1e-12);
        let tangent = elevated_tangent(&curve, d).unwrap();
        let scale = (1.0 + g * g).sqrt();
        close(tangent.z, g / scale, 1e-14);
        close(tangent.x, 0.6 / scale, 1e-14);
        close(tangent.length(), 1.0, 1e-14);
    }
    // A sag arc on a circular plan: the vertical turns vertical at its end.
    let steep = Elevated3::new(line_plan(), ElevationLaw::circular_arc(0.0, 0.0, 50.0));
    assert!(elevated_point(&steep, 49.0).is_ok());
    assert!(elevated_point(&steep, 51.0).is_err());
    assert!(elevated_tangent(&steep, 51.0).is_err());
}

/// The largest gap between a law and the chord of its heights over a span,
/// over dense samples.
fn sampled_chord_gap(law: &ElevationLaw, a: Scalar, b: Scalar) -> Scalar {
    let za = elevation_height(law, a).unwrap();
    let zb = elevation_height(law, b).unwrap();
    let n = 400;
    (0..=n)
        .map(|i| {
            let t = Scalar::from(i) / Scalar::from(n);
            let d = a + (b - a) * t;
            (elevation_height(law, d).unwrap() - (za + (zb - za) * t)).abs()
        })
        .fold(0.0, Scalar::max)
}

#[test]
fn chord_bounds_hold_on_dense_samples_and_stay_tight() {
    // The last entry is how loose the bound may be: tight where `z''` is
    // nearly uniform over the span, looser where it varies.
    let cases: Vec<(ElevationLaw, Scalar, Scalar, Scalar)> = vec![
        (
            ElevationLaw::circular_arc(100.0, -0.02, 2000.0),
            0.0,
            120.0,
            1.05,
        ),
        (
            ElevationLaw::circular_arc(100.0, -0.02, 2000.0),
            1500.0,
            1900.0,
            4.0,
        ),
        (
            ElevationLaw::circular_arc(50.0, 0.03, -5000.0),
            100.0,
            400.0,
            1.05,
        ),
        (sag_to_crest(), 0.0, 150.0, 8.0),
        (sag_to_crest(), 100.0, 120.0, 1.5),
        // Across the inflection near d = 64, where z'' changes sign.
        (sag_to_crest(), 60.0, 70.0, 6.0),
        (
            ElevationLaw::intrinsic(0.0, 0.5, CurvatureLaw::circular(0.004)),
            10.0,
            60.0,
            2.0,
        ),
        (
            ElevationLaw::parabolic(3.0, 0.02, -0.03, 200.0),
            0.0,
            200.0,
            1.05,
        ),
        (composed(), 100.0, 190.0, 8.0),
        (composed(), 200.0, 255.0, 1.5),
    ];
    for (law, a, b, looseness) in &cases {
        let bound = elevation_chord_bound(law, *a, *b).expect("bounded");
        let gap = sampled_chord_gap(law, *a, *b);
        assert!(
            gap <= bound,
            "{law:?} over [{a}, {b}]: gap {gap:e} > bound {bound:e}"
        );
        // Sound but not vacuous.
        assert!(
            bound <= looseness * gap + 1e-12,
            "{law:?} over [{a}, {b}]: bound {bound:e} vs gap {gap:e}"
        );
    }
    // A degenerate span strays nowhere; the order of the ends is irrelevant.
    assert_eq!(elevation_chord_bound(&sag_to_crest(), 5.0, 5.0), Some(0.0));
    assert_eq!(
        elevation_chord_bound(&sag_to_crest(), 70.0, 60.0),
        elevation_chord_bound(&sag_to_crest(), 60.0, 70.0)
    );
}

#[test]
fn chord_bounds_refuse_what_they_cannot_bound() {
    let law = composed();
    // A seam strictly inside the span: the grade may jump there.
    assert_eq!(elevation_chord_bound(&law, 80.0, 100.0), None);
    // Also where the piece after the seam would bound its own extension.
    assert_eq!(elevation_chord_bound(&law, 190.0, 210.0), None);
    assert_eq!(elevation_chord_bound(&law, 250.0, 270.0), None);
    // Past the vertical of an arc.
    let arc = ElevationLaw::circular_arc(0.0, 0.0, 100.0);
    assert_eq!(elevation_chord_bound(&arc, 0.0, 100.0), None);
    // Before the start of, or past the vertical of, an intrinsic profile.
    let turning = ElevationLaw::intrinsic(0.0, 0.0, CurvatureLaw::circular(0.01));
    assert_eq!(elevation_chord_bound(&turning, -1.0, 10.0), None);
    assert_eq!(elevation_chord_bound(&turning, 10.0, 120.0), None);
    assert_eq!(elevation_chord_bound(&turning, 0.0, Scalar::INFINITY), None);
}

#[test]
fn a_banked_curve_reads_an_intrinsic_profile_too() {
    // No cant and no pivot: the banked curve is its centreline, so its
    // height and tangent are the profile's.
    let base = Elevated3::new(line_plan(), sag_to_crest());
    let curve = Banked3::new(
        base,
        CantLaw::zero(150.0),
        CantLaw::zero(150.0),
        1.5,
        BankConvention::TangentRotation,
    );
    let (z, g): (Scalar, Scalar) = (22.680_253_208_798_4, 0.037_819_512_338_114_1);
    close(banked_point(&curve, 75.0).unwrap().z, z, 1e-12);
    let tangent = banked_tangent(&curve, 75.0).unwrap();
    let scale = (1.0 + g * g).sqrt();
    close(tangent.z, g / scale, 1e-14);
    close(tangent.y, 0.8 / scale, 1e-14);
}

#[test]
fn a_spiral_plan_carries_an_intrinsic_profile() {
    // The issue's real case: a clothoid in plan and a clothoid in profile.
    let plan = Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::new(1.0, 0.0),
            y: Vec2::new(0.0, 1.0),
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 120.0),
        120.0,
    ));
    let curve = Elevated3::new(plan, sag_to_crest());
    let point = elevated_point(&curve, 75.0).unwrap();
    close(point.z, 22.680_253_208_798_398, 1e-11);
}
