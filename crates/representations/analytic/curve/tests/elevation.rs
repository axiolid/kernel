//! Closed-form vertical circular arcs and the intrinsic profile's storage
//! (#238).
//!
//! Reference values are a 50-digit evaluation of the circle in `(d, z)`,
//! `z = h + R (cos t0 - sqrt(1 - (sin t0 + d/R)^2))`, in an independent
//! arbitrary-precision library -- never another run of this crate.

use axiolid_curve::{CurvatureLaw, ElevationLaw};

fn close(actual: f64, expected: f64, relative: f64) {
    let tolerance = relative * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (off by {:e}, allowed {tolerance:e})",
        (actual - expected).abs()
    );
}

#[test]
fn a_sag_arc_matches_the_circle_through_its_low_point_and_to_its_far_end() {
    // From grade -2% into a 2000 m sag: the low point is at d ~ 40.
    let sag = ElevationLaw::circular_arc(100.0, -0.02, 2000.0);
    let cases = [
        (0.0, 100.0, -0.02),
        (40.0, 99.600_119_976_004_4, 3.998_800_399_892_021e-6),
        (80.0, 100.000_319_968_032, 0.020_008_002_401_760_664),
        (500.0, 153.220_752_482_208_38, 0.236_340_358_558_617_5),
        (-1900.0, 1_613.358_379_025_382_5, -3.989_768_543_256_358_3),
    ];
    for (d, z, g) in cases {
        close(sag.height_at(d).unwrap(), z, 1e-14);
        close(sag.grade_at(d).unwrap(), g, 1e-13);
    }
    // Near the vertical: sin t = 0.99950..., the grade is 31.7.
    close(
        sag.height_at(2039.0).unwrap(),
        2_036.615_792_321_774_5,
        1e-12,
    );
    close(sag.grade_at(2039.0).unwrap(), 31.738_181_108_837_423, 1e-11);
}

#[test]
fn a_crest_arc_has_negative_radius_and_turns_the_grade_down() {
    let crest = ElevationLaw::circular_arc(50.0, 0.03, -5000.0);
    let cases = [
        (60.0, 51.439_630_689_538_21, 0.017_989_419_260_136_18),
        (150.0, 52.248_481_933_155_64, -1.349_089_433_022_498_5e-5),
        (3000.0, -839.574_853_328_575_6, -0.693_754_842_446_360_1),
    ];
    for (d, z, g) in cases {
        close(crest.height_at(d).unwrap(), z, 1e-14);
        close(crest.grade_at(d).unwrap(), g, 1e-13);
    }
    close(
        crest.height_at(5149.0).unwrap(),
        -4_851.187_623_518_919,
        1e-12,
    );
    close(
        crest.grade_at(5149.0).unwrap(),
        -51.769_530_438_172_73,
        1e-10,
    );
}

#[test]
fn beyond_the_vertical_an_arc_has_no_height_and_no_grade() {
    // sin t0 + d/R reaches +1 at d = 2000 (1 + 0.019996) = 2039.992.
    let sag = ElevationLaw::circular_arc(100.0, -0.02, 2000.0);
    for d in [2040.0, 2500.0, -1960.5, f64::INFINITY, f64::NAN] {
        assert_eq!(sag.height_at(d), None, "height at {d}");
        assert_eq!(sag.grade_at(d), None, "grade at {d}");
    }
    let crest = ElevationLaw::circular_arc(50.0, 0.03, -5000.0);
    assert_eq!(crest.height_at(5150.0), None);
    assert_eq!(crest.grade_at(5150.0), None);
    // Exactly vertical at the start is outside the domain too.
    let level_circle = ElevationLaw::circular_arc(0.0, 0.0, 10.0);
    assert_eq!(level_circle.height_at(10.0), None);
    assert!(level_circle.height_at(9.999).is_some());
}

#[test]
fn a_huge_radius_keeps_the_tiny_rise_to_full_relative_precision() {
    // d/R = 1e-9: R (cos t0 - cos t) written naively cancels to ~1e-10
    // absolute against a 1e-5 rise; the stable form keeps every digit.
    let arc = ElevationLaw::circular_arc(0.0, 0.01, 1e6);
    let rise = arc.height_at(1e-3).unwrap();
    let expected = 1.000_000_050_007_500_2e-5;
    assert!(
        ((rise - expected) / expected).abs() < 1e-14,
        "{rise:e} vs {expected:e}"
    );
}

#[test]
fn a_zero_or_non_finite_radius_is_malformed_and_reads_nothing() {
    for radius in [0.0, f64::INFINITY, f64::NAN] {
        let arc = ElevationLaw::circular_arc(1.0, 0.0, radius);
        assert!(!arc.is_well_formed());
        assert_eq!(arc.height_at(1.0), None);
        assert_eq!(arc.grade_at(1.0), None);
    }
    assert!(!ElevationLaw::circular_arc(f64::NAN, 0.0, 10.0).is_well_formed());
    assert!(!ElevationLaw::circular_arc(0.0, f64::INFINITY, 10.0).is_well_formed());
    assert!(ElevationLaw::circular_arc(0.0, 0.0, -10.0).is_well_formed());
}

#[test]
fn an_intrinsic_profile_is_stored_exactly_and_left_to_an_evaluator() {
    let law = ElevationLaw::intrinsic(5.0, -0.02, CurvatureLaw::clothoid(0.0, 1.0 / 2000.0, 100.0));
    assert!(law.is_well_formed());
    // No closed form here: refused, not approximated.
    assert_eq!(law.height_at(50.0), None);
    assert_eq!(law.grade_at(50.0), None);
    let broken = ElevationLaw::intrinsic(
        5.0,
        0.0,
        CurvatureLaw::Piecewise {
            breaks: vec![1.0],
            laws: vec![],
        },
    );
    assert!(!broken.is_well_formed());
    assert!(!ElevationLaw::intrinsic(f64::NAN, 0.0, CurvatureLaw::straight()).is_well_formed());
}

#[test]
fn piece_at_descends_to_the_innermost_law_in_its_own_distance() {
    let inner = ElevationLaw::Piecewise {
        breaks: vec![10.0],
        laws: vec![
            ElevationLaw::level(1.0),
            ElevationLaw::circular_arc(1.0, 0.0, 500.0),
        ],
    };
    let law = ElevationLaw::Piecewise {
        breaks: vec![100.0],
        laws: vec![ElevationLaw::constant_grade(0.0, 0.01), inner],
    };
    let (piece, local) = law.piece_at(115.0).unwrap();
    assert_eq!(piece, &ElevationLaw::circular_arc(1.0, 0.0, 500.0));
    assert_eq!(local, 5.0);
    // A seam belongs to the piece starting there.
    let (piece, local) = law.piece_at(100.0).unwrap();
    assert_eq!(piece, &ElevationLaw::level(1.0));
    assert_eq!(local, 0.0);
    // A leaf is its own piece.
    let leaf = ElevationLaw::level(3.0);
    assert_eq!(leaf.piece_at(7.0), Some((&leaf, 7.0)));
    assert_eq!(leaf.piece_at(f64::NAN), None);
    // The arc composes inside a piecewise law through height_at too.
    close(
        law.height_at(115.0).unwrap(),
        ElevationLaw::circular_arc(1.0, 0.0, 500.0)
            .height_at(5.0)
            .unwrap(),
        0.0,
    );
}
