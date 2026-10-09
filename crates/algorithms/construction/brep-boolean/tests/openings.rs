//! Openings cut from walls (#228): an opening's caps lying in the faces it
//! cuts, jambs running tangent into an arch, and operands placed by
//! independent rigid motions, so that faces meant to coincide, or to be
//! parallel or tangent, agree only up to rounding.
//!
//! Every result must audit clean, be a closed two-manifold, and measure
//! exactly to the closed form.

use std::f64::consts::{FRAC_PI_2, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Interval, Point2, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn solid(section: ArcRing, bottom: f64, top: f64) -> ExactBRep {
    let p = ArcPrism {
        section,
        bottom,
        top,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a solid")
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

fn difference(a: &ExactBRep, b: &ExactBRep) -> ExactBRep {
    let result = boolean(a, b, BooleanOperator::Difference, tol()).expect("difference");
    let health = geometric_audit(&result, tol());
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    result
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

#[test]
fn a_hole_flush_with_both_faces() {
    let a = solid(square(0.0, 0.0, 4.0, 4.0), 0.0, 1.0);
    let b = solid(square(1.0, 1.0, 2.0, 2.0), 0.0, 1.0);
    close("hole", volume(&difference(&a, &b)), 15.0);
}

#[test]
fn a_hole_flush_with_one_face() {
    let a = solid(square(0.0, 0.0, 4.0, 4.0), 0.0, 1.0);
    let b = solid(square(1.0, 1.0, 2.0, 2.0), 0.0, 2.0);
    close("hole", volume(&difference(&a, &b)), 15.0);
}

#[test]
fn a_recess_flush_with_the_top() {
    let a = solid(square(0.0, 0.0, 4.0, 4.0), 0.0, 1.0);
    let b = solid(square(1.0, 1.0, 2.0, 2.0), 0.5, 1.0);
    close("recess", volume(&difference(&a, &b)), 15.5);
}

#[test]
fn a_notch_flush_with_a_side_and_both_faces() {
    let a = solid(square(0.0, 0.0, 4.0, 4.0), 0.0, 1.0);
    let b = solid(square(3.0, 1.0, 4.0, 2.0), 0.0, 1.0);
    close("notch", volume(&difference(&a, &b)), 15.0);
}

// --- placed walls and openings ----------------------------------------------

/// Wall length, thickness and height.
const L: f64 = 6.0;
const T: f64 = 0.3;
const H: f64 = 3.0;

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn line(a: Point2, b: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: Vec2::new(b.x - a.x, b.y - a.y),
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    }
}

/// A `w x h` rectangle on `y = 0` under a half circle of radius `w / 2`,
/// the half circle split at `split` (radians from its right springing).
/// The jambs run tangent into the arch.
fn arched(w: f64, h: f64, split: f64) -> Profile {
    let r = w / 2.0;
    let p = Point2::new;
    let arc = |from: f64, to: f64| ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: axiolid_core::Frame2 {
                origin: p(0.0, h),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: r,
        }),
        domain: Interval::new(from, to),
        same_sense: true,
    };
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            line(p(-r, 0.0), p(r, 0.0)),
            line(p(r, 0.0), p(r, h)),
            arc(0.0, split),
            arc(split, PI),
            line(p(-r, h), p(-r, 0.0)),
        ]),
        holes: Vec::new(),
    })
}

/// The wall `[-L/2, L/2] x [-T/2, T/2] x [0, H]` under `placement`.
fn wall(placement: Transform3) -> ExactBRep {
    extrude_profile_exact(&rect(L, T), Vec3::Z, H, tol())
        .expect("a wall")
        .transformed(&placement)
        .expect("rigid")
}

/// An opening across the wall: `profile` extruded along the wall's `-y`
/// from `y = start` for `depth`, its profile origin at wall `(cx, cz)`.
fn opening(
    profile: &Profile,
    depth: f64,
    (cx, cz): (f64, f64),
    start: f64,
    placement: Transform3,
) -> ExactBRep {
    extrude_profile_exact(profile, Vec3::Z, depth, tol())
        .expect("an opening")
        .transformed(
            &(placement
                * Transform3::from_translation(Vec3::new(cx, start, cz))
                * Transform3::from_rotation_x(FRAC_PI_2)),
        )
        .expect("rigid")
}

fn placements() -> [(&'static str, Transform3); 3] {
    [
        ("identity", Transform3::IDENTITY),
        ("turned", Transform3::from_rotation_z(0.6)),
        (
            "general",
            Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
                * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7),
        ),
    ]
}

#[test]
fn a_flush_window_under_any_placement() {
    for (name, p) in placements() {
        let cut = difference(
            &wall(p),
            &opening(&rect(0.9, 1.1), T, (-1.5, 1.55), T / 2.0, p),
        );
        close(name, volume(&cut), L * T * H - 0.9 * 1.1 * T);
    }
}

#[test]
fn an_arched_opening_under_any_placement() {
    let (w, h, sill) = (1.0, 1.5, 0.3);
    let area = w * h + PI * w * w / 8.0;
    for (name, p) in placements() {
        for split in [FRAC_PI_2, 0.4 * PI] {
            let profile = arched(w, h, split);
            let w_ = wall(p);
            // Overshooting both faces, and flush with both.
            for (depth, start) in [(T + 0.2, T / 2.0 + 0.1), (T, T / 2.0)] {
                let cut = difference(&w_, &opening(&profile, depth, (1.0, sill), start, p));
                close(name, volume(&cut), L * T * H - area * T);
            }
        }
    }
}

#[test]
fn a_door_on_the_floor_then_a_window_under_any_placement() {
    for (name, p) in placements() {
        let door = opening(&rect(0.9, 2.1), T + 0.2, (-2.0, 1.05), T / 2.0 + 0.1, p);
        let window = opening(&rect(1.2, 1.2), T + 0.2, (0.2, 1.6), T / 2.0 + 0.1, p);
        let one = difference(&wall(p), &door);
        close(name, volume(&one), L * T * H - 0.9 * 2.1 * T);
        let two = difference(&one, &window);
        close(name, volume(&two), L * T * H - (0.9 * 2.1 + 1.44) * T);
    }
}

/// A door standing on a wall's floor face (its bottom exactly coplanar
/// with the wall's) and reaching a rounding error past the wall's face, or
/// stopping that short of it (#276). The door's end face is read as the
/// wall's face; the wall's floor edge along that face then crosses the
/// door's floor edge within the tolerance of the door's corner, and is cut
/// there, not a fraction of the tolerance away, where the piece between
/// lies on the wall's boundary and no point can classify it.
#[test]
fn a_floor_standing_door_a_rounding_error_past_or_short_of_the_face() {
    let wall = solid(square(0.0, 0.0, 4.0, 0.25), 0.0, 3.0);
    for gap in [-1e-9, -1e-12, -4.5e-15, 4.5e-15, 1e-12, 1e-9] {
        let door = solid(square(1.0, gap, 2.0, 0.5), 0.0, 2.0);
        let (cut, report) = axiolid_brep_boolean::boolean_with_report(
            &wall,
            &door,
            BooleanOperator::Difference,
            tol(),
        )
        .unwrap_or_else(|e| panic!("gap {gap}: {e}"));
        let topology = axiolid_topology::audit_brep(cut.topology());
        assert!(topology.is_closed_manifold(), "gap {gap}: {topology:?}");
        // The door's end moved onto the wall's face, and that is reported.
        assert!(!report.is_exact(), "gap {gap}: {report:?}");
        let moved = report
            .decisions()
            .iter()
            .map(|d| d.linear)
            .fold(0.0, f64::max);
        assert!(moved <= tol().linear(), "{report:?}");
        let got = volume(&cut);
        assert!(
            (got - (3.0 - 0.5)).abs() <= 40.0 * moved + 1e-12,
            "gap {gap}: {got}"
        );
    }
    // Ten tolerances short, the skin is kept; ten tolerances past, the
    // wall's floor edge is cut where it crosses the door's, not at the
    // door's corner.
    let door = solid(square(1.0, 1e-5, 2.0, 0.5), 0.0, 2.0);
    let cut = difference(&wall, &door);
    close("skin", volume(&cut), 3.0 - 0.5 + 2.0 * 1e-5);
    let door = solid(square(1.0, -1e-5, 2.0, 0.5), 0.0, 2.0);
    let (cut, report) =
        axiolid_brep_boolean::boolean_with_report(&wall, &door, BooleanOperator::Difference, tol())
            .expect("past by ten tolerances");
    assert!(report.is_exact(), "{report:?}");
    let topology = axiolid_topology::audit_brep(cut.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    close("past", volume(&cut), 3.0 - 0.5);
}

/// The floor-standing door a micrometre to half a millimetre past the
/// wall's face under a millimetre tolerance (#291). Merging the wall's
/// floor edge's cut into the door's corner within the tolerance (#276)
/// left faces that did not sew, and the boolean was refused; it is cut
/// again with ends merged only within the rounding floor, as before #276.
#[test]
fn a_floor_standing_door_an_authored_distance_past_the_face_is_not_refused() {
    let wall = solid(square(0.0, 0.0, 4.0, 0.25), 0.0, 3.0);
    let millimetre = Tolerance::MILLIMETRE;
    for gap in [-1e-6, -1e-4, -5e-4] {
        let door = solid(square(1.0, gap, 2.0, 0.5), 0.0, 2.0);
        let (cut, report) = axiolid_brep_boolean::boolean_with_report(
            &wall,
            &door,
            BooleanOperator::Difference,
            millimetre,
        )
        .unwrap_or_else(|e| panic!("gap {gap}: {e}"));
        let topology = axiolid_topology::audit_brep(cut.topology());
        assert!(topology.is_closed_manifold(), "gap {gap}: {topology:?}");
        // Read within the tolerance, and reported.
        let moved = report
            .decisions()
            .iter()
            .map(|d| d.linear)
            .fold(0.0, f64::max);
        assert!(moved <= millimetre.linear(), "{report:?}");
        let got = volume(&cut);
        assert!(
            (got - (3.0 - 0.5)).abs() <= 40.0 * moved + 1e-12,
            "gap {gap}: {got}"
        );
    }
}
