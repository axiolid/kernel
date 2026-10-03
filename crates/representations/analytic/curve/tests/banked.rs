//! Cant laws, bank angles and the two bank conventions against hand-computed
//! values (ADR 0081).

use axiolid_core::{Point2, Vec2};
use axiolid_curve::{
    bank_angle, BankConvention, BankError, Banked3, CantForm, CantLaw, CantPiece, CantValue,
    Curve2, Elevated3, ElevationLaw, Line2,
};

const EPS: f64 = 1e-12;

fn height(value: Option<CantValue>) -> f64 {
    match value {
        Some(CantValue::Height(v)) => v,
        other => panic!("expected a height, got {other:?}"),
    }
}

fn angle(value: Option<CantValue>) -> f64 {
    match value {
        Some(CantValue::Angle(v)) => v,
        other => panic!("expected an angle, got {other:?}"),
    }
}

fn close(actual: f64, expected: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= EPS,
        "{what}: {actual} != {expected} (off by {:e})",
        actual - expected
    );
}

fn one(piece: CantPiece) -> CantLaw {
    CantLaw::new(vec![piece])
}

#[test]
fn a_constant_cant_holds_its_value_to_both_ends() {
    let law = one(CantPiece::constant(50.0, 0.15));
    for d in [0.0, 25.0, 50.0] {
        close(height(law.value_at(d)), 0.15, "constant");
        close(height(law.rate_at(d)), 0.0, "constant rate");
    }
}

#[test]
fn a_linear_transition_at_thirty_metres_of_a_hundred() {
    let law = one(CantPiece::linear(100.0, 0.0, 0.15));
    close(height(law.value_at(30.0)), 0.045, "linear");
    close(height(law.rate_at(30.0)), 0.0015, "linear rate");
    close(height(law.value_at(100.0)), 0.15, "linear end");
}

#[test]
fn a_bloss_transition_at_xi_three_tenths() {
    // 0.15 (3 - 0.6) 0.09 = 0.0324; rate 0.15 * 6 * 0.3 * 0.7 / 100.
    let law = one(CantPiece::bloss(100.0, 0.0, 0.15));
    close(height(law.value_at(30.0)), 0.0324, "bloss");
    close(height(law.rate_at(30.0)), 0.001_89, "bloss rate");
    close(height(law.rate_at(0.0)), 0.0, "bloss start rate");
    close(height(law.rate_at(100.0)), 0.0, "bloss end rate");
    close(height(law.value_at(100.0)), 0.15, "bloss end");
}

#[test]
fn a_helmert_transition_is_two_parabolas_meeting_at_half_way() {
    let law = CantLaw::new(CantPiece::helmert(100.0, 0.0, 0.15).to_vec());
    assert_eq!(law.seams(), vec![50.0]);
    // 2 dD x^2 at x = 1/4; dD - 2 dD (1 - x)^2 at x = 3/4.
    close(height(law.value_at(25.0)), 0.018_75, "helmert first half");
    close(height(law.value_at(75.0)), 0.131_25, "helmert second half");
    close(height(law.value_at(100.0)), 0.15, "helmert end");
    // Value AND slope continuous at the seam: 4 dD x / L = 0.003 at x = 1/2.
    let first = &law.pieces[0];
    let left = CantLaw::new(vec![first.clone()]);
    close(
        height(left.value_at(50.0)),
        0.075,
        "helmert seam from the left",
    );
    close(
        height(law.value_at(50.0)),
        0.075,
        "helmert seam from the right",
    );
    close(
        height(left.rate_at(50.0)),
        0.003,
        "helmert slope from the left",
    );
    close(
        height(law.rate_at(50.0)),
        0.003,
        "helmert slope from the right",
    );
}

#[test]
fn a_cosine_transition_at_a_quarter() {
    // 0.05 + 0.1 (1 - cos(pi / 4)) / 2; rate 0.1 pi / 2 sin(pi / 4) / 80.
    let law = one(CantPiece::cosine(80.0, 0.05, 0.15));
    close(
        height(law.value_at(20.0)),
        0.064_644_660_940_672_63,
        "cosine",
    );
    close(
        height(law.rate_at(20.0)),
        0.001_388_400_918_174_489,
        "cosine rate",
    );
    close(height(law.value_at(80.0)), 0.15, "cosine end");
}

#[test]
fn a_sine_transition_at_a_quarter() {
    // 0.15 (1/4 - sin(pi / 2) / (2 pi)); rate 0.15 (1 - cos(pi / 2)) / 100.
    let law = one(CantPiece::sine(100.0, 0.0, 0.15));
    close(
        height(law.value_at(25.0)),
        0.013_626_758_536_215_698,
        "sine",
    );
    close(height(law.rate_at(25.0)), 0.0015, "sine rate");
    close(height(law.value_at(100.0)), 0.15, "sine end");
}

#[test]
fn a_viennese_bend_gives_the_angle_itself() {
    let law = one(CantPiece::viennese_bend(100.0, 0.0, 0.1));
    assert!(law.has_angle_pieces());
    // xi = 1/2: xi^4 (35 - 42 + 17.5 - 2.5) = 8 / 16.
    close(angle(law.value_at(50.0)), 0.05, "viennese half way");
    // xi = 1/4: (1/256)(35 - 21 + 4.375 - 0.3125).
    close(
        angle(law.value_at(25.0)),
        0.007_055_664_062_5,
        "viennese quarter",
    );
    close(angle(law.value_at(100.0)), 0.1, "viennese end");
    // Rate 140 xi^3 (1 - xi)^3 dpsi / L: 140 / 64 * 0.1 / 100.
    close(angle(law.rate_at(50.0)), 0.002_187_5, "viennese rate");
    close(angle(law.rate_at(0.0)), 0.0, "viennese start rate");
    close(angle(law.rate_at(100.0)), 0.0, "viennese end rate");
}

#[test]
fn a_law_is_continuous_across_authored_seams_and_owned_by_the_later_piece() {
    let law = CantLaw::new(vec![
        CantPiece::linear(60.0, 0.0, 0.15),
        CantPiece::constant(40.0, 0.15),
        CantPiece::cosine(60.0, 0.15, 0.0),
    ]);
    assert_eq!(law.seams(), vec![60.0, 100.0]);
    close(law.length(), 160.0, "length");
    for seam in [60.0, 100.0] {
        let before = height(law.value_at(seam - 1e-9));
        close(height(law.value_at(seam)), 0.15, "seam value");
        assert!((before - 0.15).abs() < 1e-10, "jump at {seam}");
    }
    // The seam belongs to the piece that starts there: the constant one.
    close(
        height(law.rate_at(60.0)),
        0.0,
        "rate owned by the later piece",
    );
    close(height(law.value_at(160.0)), 0.0, "closed at the far end");
}

#[test]
fn a_law_is_not_extrapolated() {
    let law = one(CantPiece::linear(100.0, 0.0, 0.15));
    assert_eq!(law.value_at(-1e-9), None);
    assert_eq!(law.value_at(100.000_001), None);
    assert_eq!(law.value_at(f64::NAN), None);
    assert_eq!(CantLaw::default().value_at(0.0), None);
    let zero_length = one(CantPiece::constant(0.0, 0.1));
    assert!(!zero_length.is_well_formed());
    assert_eq!(zero_length.value_at(0.0), None);
    let non_finite = one(CantPiece::new(
        10.0,
        CantForm::Polynomial {
            coefficients: vec![f64::NAN],
        },
    ));
    assert!(!non_finite.is_well_formed());
}

#[test]
fn a_cant_beyond_the_rail_head_distance_is_refused_by_name() {
    close(
        bank_angle(0.15, 1.5).unwrap(),
        0.100_167_421_161_559_8,
        "asin 0.1",
    );
    close(
        bank_angle(1.5, 1.5).unwrap(),
        core::f64::consts::FRAC_PI_2,
        "quarter",
    );
    assert_eq!(
        bank_angle(1.6, 1.5),
        Err(BankError::CantExceedsRailHeadDistance {
            cant: 1.6,
            rail_head_distance: 1.5
        })
    );
    assert!(matches!(
        bank_angle(-1.6, 1.5),
        Err(BankError::CantExceedsRailHeadDistance { .. })
    ));
    assert_eq!(bank_angle(0.1, 0.0), Err(BankError::RailHeadDistance(0.0)));
    let message = bank_angle(1.6, 1.5).unwrap_err().to_string();
    assert!(
        message.contains("exceeds the rail-head distance"),
        "{message}"
    );
}

#[test]
fn the_two_conventions_differ_by_one_minus_cos_theta_at_a_two_percent_grade() {
    let grade_cosine = 1.0 / 1.0004_f64.sqrt();
    let (cant, b) = (0.15, 1.5);
    let rotation = BankConvention::TangentRotation
        .roll(cant, b, grade_cosine)
        .unwrap();
    let rise = BankConvention::VerticalRise
        .roll(cant, b, grade_cosine)
        .unwrap();
    close(
        rotation,
        0.100_167_421_161_559_8,
        "tangent rotation is asin(D/b)",
    );
    // asin(D sqrt(1 + g^2) / b).
    close(
        rise,
        0.100_187_519_928_492_4,
        "vertical rise is asin(D/(b cos))",
    );
    // Vertical rise across the heads: b sin(rho) cos(theta).
    let rotated = b * rotation.sin() * grade_cosine;
    let risen = b * rise.sin() * grade_cosine;
    close(risen, cant, "vertical rise is exact");
    close(
        rotated,
        cant * grade_cosine,
        "tangent rotation rises D cos theta",
    );
    // D (1 - cos theta) = 1.9994002e-4 D.
    let difference = risen - rotated;
    close(
        difference / cant,
        1.999_400_199_929_057_8e-4,
        "relative difference",
    );
    // Level track: the two agree.
    assert_eq!(
        BankConvention::TangentRotation.roll(cant, b, 1.0),
        BankConvention::VerticalRise.roll(cant, b, 1.0)
    );
}

#[test]
fn a_vertical_rise_beyond_the_tilted_span_is_refused_where_a_rotation_is_not() {
    // At 10%, b cos(theta) = 1.4925 < 1.4999.
    let grade_cosine = 1.0 / 1.01_f64.sqrt();
    assert!(BankConvention::TangentRotation
        .roll(1.4999, 1.5, grade_cosine)
        .is_ok());
    assert!(matches!(
        BankConvention::VerticalRise.roll(1.4999, 1.5, grade_cosine),
        Err(BankError::CantExceedsVerticalSpan { .. })
    ));
    assert!(matches!(
        BankConvention::TangentRotation.roll_from_angle(1.6, 1.5, 1.0),
        Err(BankError::AngleOutOfRange { .. })
    ));
}

fn straight() -> Elevated3 {
    Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::new(0.0, 0.0),
            direction: Vec2::new(1.0, 0.0),
        }),
        ElevationLaw::constant_grade(100.0, 0.02),
    )
}

#[test]
fn a_banked_curve_reads_its_cant_and_bank_angle_by_station() {
    let curve = Banked3::new(
        straight(),
        CantLaw::new(vec![
            CantPiece::linear(100.0, 0.0, 0.15),
            CantPiece::viennese_bend(100.0, 0.1, 0.0),
        ]),
        CantLaw::zero(200.0),
        1.5,
        BankConvention::TangentRotation,
    );
    close(curve.span(), 200.0, "span");
    close(curve.cant_at(30.0).unwrap(), 0.045, "cant at 30");
    close(
        curve.bank_angle_at(30.0).unwrap(),
        (0.045_f64 / 1.5).asin(),
        "bank angle at 30",
    );
    // The angle piece gives psi; the cant follows as b sin(psi).
    close(curve.bank_angle_at(150.0).unwrap(), 0.05, "angle piece");
    close(
        curve.cant_at(150.0).unwrap(),
        0.074_968_753_906_017_5,
        "b sin(psi)",
    );
    assert!(matches!(
        curve.cant_at(200.5),
        Err(BankError::OutsideLaw { law: "cant", .. })
    ));
}

#[test]
fn a_pivot_must_be_an_elevation() {
    let curve = Banked3::new(
        straight(),
        CantLaw::new(vec![CantPiece::constant(10.0, 0.1)]),
        CantLaw::new(vec![CantPiece::viennese_bend(10.0, 0.0, 0.1)]),
        1.5,
        BankConvention::VerticalRise,
    );
    assert_eq!(curve.pivot_at(5.0), Err(BankError::AngleInPivot));
}
