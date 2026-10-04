//! A round web hole touching an I-beam's flange and its root fillets
//! (#249).
//!
//! The beam `[-B/2, B/2] x [-H/2, H/2] x [0, LEN]` has quarter-cylinder
//! root fillets of radius `RF` between web and flanges, their axes along
//! `z`. The hole, radius `R` with its axis along `x` at height `y` and
//! `z = LEN / 2`, runs past the flange tips or just past the fillets:
//!
//! - **tangent**: the hole's top ruling lies in the top flange's inner
//!   face. It also touches each top fillet where the fillet meets the
//!   flange: the hole/fillet section is a quartic with a double point
//!   there (two loops round the fillet cylinder crossing at the
//!   fillet/flange edge);
//! - **into the fillet**: `10 mm` lower, the hole cuts each fillet
//!   transversally;
//! - **clear**: `50 mm` lower, it misses the fillets.
//!
//! Each result audits clean, is a closed two-manifold and measures the
//! reference volume: the beam's area times its length less the hole's
//! part inside the beam, the fillet part integrated in closed form per
//! slice (a circular segment) and summed by Simpson's rule in the fillet's
//! angle (`removed`). With exact axes at `Tolerance::ZERO` the report is
//! empty; under a general placement at `Tolerance::METRE` the flange is
//! read as touching the hole (#243) and the double point is placed on the
//! contact. A hole a fraction of the tolerance into or short of the
//! flange meets each fillet in two arcs the contact would have to join:
//! refused by name.

use std::f64::consts::PI;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{
    boolean_with_report, BooleanError, BooleanReport, ToleranceDecisionKind,
};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Mat3, Tolerance, Transform3, Vec3};
use axiolid_measure::exact_properties;
use axiolid_profile::{CircleProfile, Profile, SectionProfile};

const H: f64 = 0.5;
const B: f64 = 0.25;
const TW: f64 = 0.125;
const TF: f64 = 0.0625;
const LEN: f64 = 2.0;
const R: f64 = 0.125;
const RF: f64 = 0.03125;
const EPS: f64 = 1e-6;

/// The flange's inner face, the fillets' axes and the hole's axis height
/// when it touches the flange.
const Y_TOP: f64 = H / 2.0 - TF;
const Y_C: f64 = Y_TOP - RF;
const X_C: f64 = TW / 2.0 + RF;
const TOUCHING: f64 = Y_TOP - R;

fn beam(p: Transform3) -> ExactBRep {
    let profile = Profile::Section(SectionProfile::I {
        depth: H,
        width: B,
        web_thickness: TW,
        flange_thickness: TF,
        fillet_radius: Some(RF),
        flange_edge_radius: None,
        flange_slope: None,
    });
    extrude_profile_exact(&profile, Vec3::Z, LEN, Tolerance::ZERO)
        .expect("beam")
        .transformed(&p)
        .expect("rigid")
}

/// Where the hole starts along `x` and how far it runs: past the flange
/// tips, or past the fillets only.
const RUNS: [(f64, f64); 2] = [(-0.75 * B, 1.5 * B), (-TW, 2.0 * TW)];

/// The hole, its axis along `x` at height `y`, turned by an exact matrix.
fn hole(y: f64, (start, depth): (f64, f64), p: Transform3) -> ExactBRep {
    let profile = Profile::Circle(CircleProfile {
        radius: R,
        thickness: None,
    });
    let onto_x = Transform3::from_mat3(Mat3::from_cols(-Vec3::Z, Vec3::Y, Vec3::X));
    extrude_profile_exact(&profile, Vec3::Z, depth, Tolerance::ZERO)
        .expect("hole")
        .transformed(&(p * Transform3::from_translation(Vec3::new(start, y, LEN / 2.0)) * onto_x))
        .expect("rigid")
}

/// Placements with exact entries: none, and a quarter or half turn with a
/// dyadic offset.
fn exact_placements() -> [Transform3; 3] {
    let quarter = Mat3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z);
    let half = Mat3::from_cols(-Vec3::X, -Vec3::Y, Vec3::Z);
    [
        Transform3::IDENTITY,
        Transform3::from_mat3_translation(quarter, Vec3::new(12.5, -4.0, 3.25)),
        Transform3::from_mat3_translation(half, Vec3::new(-7.75, 2.5, 0.0)),
    ]
}

/// A rigid placement with no axis left aligned.
fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

/// The area of the hole's cross-section above the line `d` above its axis.
fn segment(d: f64) -> f64 {
    if d >= R {
        0.0
    } else if d <= -R {
        PI * R * R
    } else {
        R * R * (d / R).acos() - d * (R * R - d * d).sqrt()
    }
}

/// The volume the hole, its axis at height `y`, removes from the beam,
/// for a run covering both top fillets and no more than the flanges.
///
/// Sliced across `x`: through the web the whole disc; through a fillet,
/// at `x = X_C - RF sin t`, the material lies above `Y_C + RF cos t`, so
/// the slice is a segment, integrated by Simpson's rule in `t` (smooth
/// there); under the flange the material lies above `Y_TOP`.
fn removed(y: f64, (start, depth): (f64, f64)) -> f64 {
    assert!(start <= -X_C && start + depth >= X_C);
    let n = 200_000;
    let h = 0.5 * PI / n as f64;
    let mut fillet = 0.0;
    for i in 0..=n {
        let t = i as f64 * h;
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        fillet += w * segment(Y_C + RF * t.cos() - y) * RF * t.cos();
    }
    fillet *= h / 3.0;
    let flange = 2.0 * (start + depth).min(B / 2.0) - 2.0 * X_C;
    PI * R * R * TW + 2.0 * fillet + flange.max(0.0) * segment(Y_TOP - y)
}

fn beam_volume() -> f64 {
    let area = 2.0 * B * TF + (H - 2.0 * TF) * TW + 4.0 * RF * RF * (1.0 - PI / 4.0);
    area * LEN
}

/// The beam minus the hole under `p`: audited and measured against the
/// reference within `slack`.
fn cut(
    y: f64,
    run: (f64, f64),
    p: Transform3,
    tolerance: Tolerance,
    slack: f64,
) -> Result<BooleanReport, BooleanError> {
    let (result, report) = boolean_with_report(
        &beam(p),
        &hole(y, run, p),
        BooleanOperator::Difference,
        tolerance,
    )?;
    let health = geometric_audit(&result, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    let measured = exact_properties(&result, Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    let expected = beam_volume() - removed(y, run);
    assert!(
        (measured - expected).abs() <= slack + 1e-12 * expected,
        "y {y} run {run:?}: volume {measured}, expected {expected}"
    );
    Ok(report)
}

/// Touching, into the fillets, and clear of them.
const HEIGHTS: [f64; 3] = [TOUCHING, TOUCHING - 0.01, TOUCHING - 0.05];

#[test]
fn a_hole_touching_the_flange_and_its_fillets_is_exact_with_exact_axes() {
    // The double point is the exact double root of the fillet/flange edge
    // against the hole; nothing is read within tolerance.
    for p in exact_placements() {
        for y in HEIGHTS {
            for run in RUNS {
                for tolerance in [Tolerance::ZERO, Tolerance::METRE] {
                    let report = cut(y, run, p, tolerance, 0.0)
                        .unwrap_or_else(|e| panic!("y {y} run {run:?}: {e}"));
                    assert!(
                        report.linear() <= report.rounding_floor(),
                        "y {y} run {run:?}: {report:?}"
                    );
                    if tolerance == Tolerance::ZERO {
                        assert!(report.is_exact(), "y {y} run {run:?}: {report:?}");
                    }
                }
            }
        }
    }
}

#[test]
fn a_hole_touching_the_flange_and_its_fillets_under_a_general_placement() {
    // The flange is read as touching the hole, and the trace's double
    // point lies on that contact to rounding.
    for y in HEIGHTS {
        for run in RUNS {
            let report = cut(y, run, general(), Tolerance::METRE, 1e-12)
                .unwrap_or_else(|e| panic!("y {y} run {run:?}: {e}"));
            assert!(report.linear() <= EPS, "y {y} run {run:?}: {report:?}");
            if y == TOUCHING {
                assert!(
                    report.contains(ToleranceDecisionKind::PlaneTouchesCylinder),
                    "{report:?}"
                );
            }
        }
    }
}

#[test]
fn a_hole_a_fraction_of_the_tolerance_off_the_flange_is_refused_by_name() {
    // Read as touching the flange, the hole would have to touch each
    // fillet too, where it meets it in two arcs `sqrt(2 R gap)` from the
    // contact: no single move within the tolerance explains both.
    for fraction in [0.5, -0.5] {
        for p in [Transform3::IDENTITY, general()] {
            for run in RUNS {
                let error = cut(TOUCHING + fraction * EPS, run, p, Tolerance::METRE, 0.0)
                    .expect_err("refused");
                assert_eq!(
                    error,
                    BooleanError::UnsupportedContact,
                    "{fraction} {run:?}"
                );
            }
        }
    }
}

#[test]
fn a_hole_ten_tolerances_into_the_flange_is_decided_exactly() {
    // Beyond the tolerance nothing is read: the hole cuts a groove under
    // the flange and meets each fillet edge in two separate points.
    let y = TOUCHING + 10.0 * EPS;
    let groove = segment(Y_TOP - y) * (B - 2.0 * X_C);
    for run in RUNS {
        let report = cut(
            y,
            run,
            Transform3::IDENTITY,
            Tolerance::METRE,
            1e-6 * groove,
        )
        .unwrap_or_else(|e| panic!("{run:?}: {e}"));
        // (The short run's caps touch the fillet cylinders past the
        // fillet faces, to the construction's rounding: below the floor.)
        assert!(report.linear() <= report.rounding_floor(), "{report:?}");
    }
}
