//! A curved wall less a niche cut from its inner face (#324).
//!
//! The wall is a quarter ring, radii `4.8` to `5` about the origin, `0` to
//! `90` degrees, `3` high. The niche is an annular sector from `40` to `50`
//! degrees, `z` from `1` to `2`, from radius `4.9` in to its inner arc:
//! on the wall's inner cylinder (flush), one unit in the last place inside
//! it (`4.8.next_down()`, so it reaches out of the wall by `1e-15`), or
//! well past it (`4.7`, reaching `0.1` out of the wall). The flush niche's
//! inner face lies on the wall's inner cylinder (the same axis and radius,
//! up to the rounding of building each from its profile) and the 1-ulp
//! niche's within the tolerance of it: both are read as one surface, split
//! along each other's boundaries, the shared patch dropped by the
//! difference.
//!
//! Oracle: the closed-form volume of the annular sectors, measured exactly
//! (ADR 0073) on results that audit clean.

use std::f64::consts::FRAC_PI_2;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{boolean_with_report, BooleanReport, ToleranceDecisionKind};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Frame2, Interval, Point2, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

fn tol() -> Tolerance {
    Tolerance::MILLIMETRE
}

/// An annular sector between radii `inner` and `outer` from angle `from` to
/// `to`, anticlockwise: outer arc, line in, inner arc back, line out, as
/// `ifc-geometry` lowers an `IfcIndexedPolyCurve` of two arcs and two lines.
fn sector(outer: f64, inner: f64, from: f64, to: f64) -> Profile {
    let circle = |radius| {
        Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        })
    };
    let at = |radius: f64, angle: f64| Point2::new(radius * angle.cos(), radius * angle.sin());
    let line = |a: Point2, b: Point2| ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    };
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve: circle(outer),
                domain: Interval::new(from, to),
                same_sense: true,
            },
            line(at(outer, to), at(inner, to)),
            ProfileSegment {
                curve: circle(inner),
                domain: Interval::new(from, to),
                same_sense: false,
            },
            line(at(inner, from), at(outer, from)),
        ]),
        holes: Vec::new(),
    })
}

fn extruded(profile: &Profile, z: f64, depth: f64) -> ExactBRep {
    extrude_profile_exact(profile, Vec3::Z, depth, tol())
        .expect("an extrusion")
        .transformed(&Transform3::from_translation(Vec3::new(0.0, 0.0, z)))
        .expect("rigid")
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

/// The volume of an annular sector extruded by `height`.
fn sector_volume(outer: f64, inner: f64, sweep: f64, height: f64) -> f64 {
    0.5 * sweep * (outer * outer - inner * inner) * height
}

const OUTER: f64 = 5.0;
const INNER: f64 = 4.8;
const HEIGHT: f64 = 3.0;
const NICHE_OUTER: f64 = 4.9;

fn wall() -> ExactBRep {
    extruded(&sector(OUTER, INNER, 0.0, FRAC_PI_2), 0.0, HEIGHT)
}

fn niche(inner: f64) -> ExactBRep {
    let (from, to) = (40.0_f64.to_radians(), 50.0_f64.to_radians());
    extruded(&sector(NICHE_OUTER, inner, from, to), 1.0, 1.0)
}

/// The wall less the niche: audited, closed, and of the closed-form volume
/// of the wall less the part of the niche inside it.
fn cut(niche_inner: f64) -> BooleanReport {
    let (result, report) = boolean_with_report(
        &wall(),
        &niche(niche_inner),
        BooleanOperator::Difference,
        tol(),
    )
    .unwrap_or_else(|e| panic!("niche from {niche_inner}: {e}"));
    let health = geometric_audit(&result, tol());
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    let sweep = 10.0_f64.to_radians();
    let want = sector_volume(OUTER, INNER, FRAC_PI_2, HEIGHT)
        - sector_volume(NICHE_OUTER, niche_inner.max(INNER), sweep, 1.0);
    let got = volume(&result);
    assert!(
        (got - want).abs() <= 1e-9 * want,
        "niche from {niche_inner}: expected {want}, got {got}"
    );
    report
}

/// The reading of the niche's inner cylinder as the wall's: the furthest
/// it moved the niche.
fn coincidence(report: &BooleanReport) -> Option<f64> {
    report
        .decisions()
        .iter()
        .find(|d| d.kind == ToleranceDecisionKind::CoincidentSupports)
        .map(|d| d.linear)
}

#[test]
fn a_niche_flush_with_the_inner_face_is_cut() {
    // Both cylinders are built from rounded profile points (the wall's
    // radius is `4.800000000000001`, the niche's `4.799999999999999`), so
    // even the flush niche is one surface with the wall only within
    // rounding, and the reading is reported.
    let report = cut(INNER);
    let moved = coincidence(&report).unwrap_or_else(|| panic!("{report:?}"));
    assert!(moved <= 1e-14, "{report:?}");
}

#[test]
fn a_niche_one_ulp_inside_the_inner_face_is_read_as_flush() {
    // The niche's inner cylinder is the wall's, moved by the ulp it misses
    // it by and the rounding of both: reported, far within the millimetre.
    let report = cut(INNER.next_down());
    let moved = coincidence(&report).unwrap_or_else(|| panic!("{report:?}"));
    assert!(moved > 0.0 && moved <= 1e-14, "{report:?}");
}

#[test]
fn a_niche_well_inside_the_wall_is_unchanged() {
    // No face of the niche lies on a face of the wall.
    let report = cut(4.7);
    assert_eq!(coincidence(&report), None, "{report:?}");
}
