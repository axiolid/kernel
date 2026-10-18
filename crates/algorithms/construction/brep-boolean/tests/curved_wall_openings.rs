//! A curved wall less a row of openings (#332), and a box tangent to a
//! curved wall's inner face (#327).
//!
//! The wall is an annular sector, radii `6` to `6.25` about the origin,
//! `0` to `60` degrees, `3` high. Each opening is a wedge through it: a
//! quadrilateral from radius `5` to `7` between two angles, its sides
//! radial, extruded from `z = 0.5` to `2.5`, as a delivered model lowers
//! openings whose sides follow the wall's radii. Openings next to each
//! other abut, one's side reaching past its neighbour's as the model's
//! rounding left it: by `DRIFT` (about `15 um` at the wall) or by `NEAR`
//! (about `0.6 mm`, just under the millimetre tolerance). The two side
//! planes then cross on the axis, turned by far more than the angular
//! tolerance, so they are not one surface, but they lie within the linear
//! tolerance of each other in the wall. Read within the tolerance, the
//! sections' cuts there merged at a point on neither face's boundary, and
//! the general boolean refused the cut ("a point lies too close to a face
//! boundary to classify", or faces that do not sew). It is now cut again
//! at a finer tolerance (`FINER_FACTOR` of it), where the sides are apart
//! and decided by the exact predicates. A wider wall (radius `24`) whose
//! openings abut by `7 um` checks that the split keeps such vertices apart
//! in the face's parameters too.
//!
//! The box of #327 lies in the tangent frame of a quarter ring's inner
//! cylinder at 45 degrees, its back face on the tangent plane (or a unit
//! in the last place either side): that face lies in the wall but along
//! the contact ruling, where it touches the inner face. Its points were
//! all within the tolerance of the cylinder near the ruling and refused;
//! read exactly, they are off it, and the box is cut at the caller's
//! tolerance.
//!
//! Oracle: the closed-form volume of the annular sectors and of the box,
//! measured exactly (ADR 0073) on results that audit clean.

use std::f64::consts::PI;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{boolean_with_report, BooleanReport, FINER_FACTOR};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Frame2, Interval, Point2, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

fn tol() -> Tolerance {
    Tolerance::MILLIMETRE
}

fn line(a: Point2, b: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn at(radius: f64, angle: f64) -> Point2 {
    Point2::new(radius * angle.cos(), radius * angle.sin())
}

/// An annular sector between radii `inner` and `outer` from angle `from` to
/// `to`, anticlockwise: outer arc, line in, inner arc back, line out.
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

/// A quadrilateral through `corners`, anticlockwise.
fn quad(corners: [Point2; 4]) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(
            (0..4)
                .map(|i| line(corners[i], corners[(i + 1) % 4]))
                .collect(),
        ),
        holes: Vec::new(),
    })
}

/// A wedge with radial sides from radius `near` to `far` between angles
/// `from` and `to`.
fn wedge(near: f64, far: f64, from: f64, to: f64) -> Profile {
    quad([at(near, from), at(far, from), at(far, to), at(near, to)])
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

/// Audited clean and a closed two-manifold.
fn sound(brep: &ExactBRep, what: &str) {
    let health = geometric_audit(brep, tol());
    assert!(health.is_consistent(), "{what}: {:?}", health.defects());
    let topology = axiolid_topology::audit_brep(brep.topology());
    assert!(topology.is_closed_manifold(), "{what}: {topology:?}");
}

/// The volume of an annular sector extruded by `height`.
fn sector_volume(outer: f64, inner: f64, sweep: f64, height: f64) -> f64 {
    0.5 * sweep * (outer * outer - inner * inner) * height
}

const OUTER: f64 = 6.25;
const INNER: f64 = 6.0;
const HEIGHT: f64 = 3.0;
const SWEEP: f64 = PI / 3.0;
/// How far an opening's side reaches past its neighbour's: more than a
/// thousand times the angular tolerance, `15 um` at the wall.
const DRIFT: f64 = 2.5e-6;
/// The same, `0.6 mm` at the wall: just under the linear tolerance.
const NEAR: f64 = 1e-4;

fn deg(d: f64) -> f64 {
    d.to_radians()
}

/// The openings in the order they are cut, each `(from, to)`: one, two
/// reaching `DRIFT` past the side of the one cut before, one reaching
/// `NEAR` past the first's other side, and one on its own.
fn openings() -> [(f64, f64); 5] {
    [
        (deg(20.0), deg(26.0)),
        (deg(14.0), deg(20.0) + DRIFT),
        (deg(8.0), deg(14.0) + DRIFT),
        (deg(26.0) - NEAR, deg(32.0)),
        (deg(40.0), deg(48.0)),
    ]
}

/// The angle the union of `spans` covers.
fn swept(spans: &[(f64, f64)]) -> f64 {
    let mut sorted = spans.to_vec();
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut total = 0.0;
    let mut reach = f64::NEG_INFINITY;
    for (from, to) in sorted {
        let start = from.max(reach);
        if to > start {
            total += to - start;
        }
        reach = reach.max(to);
    }
    total
}

/// A wall between radii `inner` and `outer` from `0` to `sweep`, less
/// `openings` (wedges a metre past either face) one after another, each
/// step sound and of the closed-form volume; each step's report.
fn cut_all((inner, outer, sweep): (f64, f64, f64), openings: &[(f64, f64)]) -> Vec<BooleanReport> {
    let mut result = extruded(&sector(outer, inner, 0.0, sweep), 0.0, HEIGHT);
    let mut reports = Vec::new();
    let mut removed: Vec<(f64, f64)> = Vec::new();
    for (index, &(from, to)) in openings.iter().enumerate() {
        let opening = extruded(&wedge(inner - 1.0, outer + 1.0, from, to), 0.5, 2.0);
        let (next, report) =
            boolean_with_report(&result, &opening, BooleanOperator::Difference, tol())
                .unwrap_or_else(|e| panic!("opening {index}: {e}"));
        sound(&next, &format!("opening {index}"));
        removed.push((from, to));
        let want = sector_volume(outer, inner, sweep, HEIGHT)
            - sector_volume(outer, inner, swept(&removed), 2.0);
        let got = volume(&next);
        assert!(
            (got - want).abs() <= 1e-9 * want,
            "opening {index}: expected {want}, got {got}"
        );
        result = next;
        reports.push(report);
    }
    reports
}

#[test]
fn a_row_of_abutting_openings_is_cut_at_the_closed_form_volume() {
    let reports = cut_all((INNER, OUTER, SWEEP), &openings());
    let finer = tol().linear() * FINER_FACTOR;
    // The first opening is cut at the millimetre; those abutting another
    // at the finer tolerance. (So is the one on its own: the wall it cuts
    // has the abutments' features, a fraction of a millimetre apart.)
    assert!(!reports[0].finer(), "{:?}", reports[0]);
    for (index, report) in reports.iter().enumerate().skip(1) {
        assert!(report.finer(), "opening {index}: {report:?}");
        assert!(report.linear() <= finer, "opening {index}: {report:?}");
    }
}

#[test]
fn a_wide_wall_whose_openings_abut_by_a_few_micrometres_is_cut() {
    // Radii `24` to `24.1`: abutting sides `3e-7` radians, `7 um`, past
    // each other. That is above the finer tolerance, so the sides are
    // apart there, and below what the face's relative slack in parameters
    // spans along the arc at this radius, which welded the two vertices
    // a split kept apart (two pieces leaving one vertex in one direction).
    let drift = 3e-7;
    let openings = [
        (deg(10.0), deg(12.0)),
        (deg(8.0), deg(10.0) + drift),
        (deg(12.0) - drift, deg(14.0)),
    ];
    let reports = cut_all((24.0, 24.1, deg(20.0)), &openings);
    let finer = tol().linear() * FINER_FACTOR;
    for (index, report) in reports.iter().enumerate().skip(1) {
        assert!(report.finer(), "opening {index}: {report:?}");
        assert!(report.linear() <= finer, "opening {index}: {report:?}");
    }
}

/// A box in the tangent frame at `45` degrees: from `back` to `back + 0.1`
/// along the radius, `0.3` either side of it, `z` from `1` to `2`.
fn tangent_box(back: f64) -> ExactBRep {
    let profile = quad([
        Point2::new(back, -0.3),
        Point2::new(back + 0.1, -0.3),
        Point2::new(back + 0.1, 0.3),
        Point2::new(back, 0.3),
    ]);
    extruded(&profile, 1.0, 1.0)
        .transformed(&Transform3::from_rotation_z(deg(45.0)))
        .expect("rigid")
}

/// The quarter ring of #327: radii `4.8` to `5`, `0` to `90` degrees, `3`
/// high.
fn quarter_ring() -> ExactBRep {
    extruded(&sector(5.0, 4.8, 0.0, PI / 2.0), 0.0, 3.0)
}

#[test]
fn a_box_tangent_to_the_inner_face_is_cut() {
    // The box's back face lies in the wall but along the ruling where it
    // touches the inner cylinder (or a unit in the last place past or
    // short of it), and its front face inside the wall: the whole box,
    // `0.1 x 0.6 x 1`, is removed, the ruling adding nothing. Decided at
    // the caller's tolerance, not the finer one.
    let want = sector_volume(5.0, 4.8, PI / 2.0, 3.0) - 0.1 * 0.6 * 1.0;
    for back in [4.8, 4.8_f64.next_down(), 4.8_f64.next_up()] {
        let (result, report) = boolean_with_report(
            &quarter_ring(),
            &tangent_box(back),
            BooleanOperator::Difference,
            tol(),
        )
        .unwrap_or_else(|e| panic!("back at {back}: {e}"));
        sound(&result, &format!("back at {back}"));
        let got = volume(&result);
        assert!(
            (got - want).abs() <= 1e-12 * want,
            "back at {back}: expected {want}, got {got}"
        );
        assert!(!report.finer(), "back at {back}: {report:?}");
        assert!(report.linear() <= 1e-14, "back at {back}: {report:?}");
    }
}
