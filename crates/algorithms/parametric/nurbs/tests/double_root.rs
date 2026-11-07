//! A traced section against a surface that touches it (#397, #427): the
//! double root comes back once, with multiplicity 2, not as a pair of
//! spurious roots either side of it, and a contact of higher order comes
//! back with its order instead of hanging the isolation.
//!
//! The fixture is #381's: the plane `z = 1` meets the torus about `y` with
//! `R = 1`, `r = 0.1` (the top face of that boolean's box) in two traced
//! halves of a loop. The right half reaches its largest `x`, `sqrt(0.21)`,
//! in the middle of its span, where the plane `x = sqrt(0.21)` is tangent to
//! it; a plane `d` lower crosses it twice and one `d` higher misses it.
//! The tip is also where a cylinder of the right radius osculates the
//! half: its curvature there is `11 / sqrt(0.21)`, and the contact of the
//! cylinder about `z` with that radius has order 4 and keeps the curve
//! outside it. The plane `y = 0.1` is a quartic contact too, at the span's
//! end. The boolean of #381 at `Tolerance::ZERO` covers the certificate end
//! to end.

use axiolid_core::{Frame3, Interval, Point3, Vec3};
use axiolid_curve::Curve3;
use axiolid_evaluate::evaluate3;
use axiolid_nurbs::{
    exact_surface_intersection, section_curve_surface_intersection, ExactCurveHit,
    ExactCurveIntersection, ExactCurveRefusal,
};
use axiolid_surface::{Cylinder, Plane, Surface, Torus};

/// The largest `x` of the right half: `sqrt(1.1^2 - 1)`.
const TANGENT: f64 = 0.458_257_569_495_584_2;

fn frame(origin: Point3, axis: Vec3) -> Frame3 {
    let z = axis.normalize();
    let helper = if z.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let x = helper.cross(z).normalize();
    Frame3 {
        origin,
        x,
        y: z.cross(x),
        z,
    }
}

fn plane(origin: Point3, normal: Vec3) -> Surface {
    Surface::Plane(Plane {
        frame: frame(origin, normal),
    })
}

/// The right half of the traced loop of the plane `z = 1` on the torus.
fn half() -> (Curve3, Interval) {
    let torus = Surface::Torus(Torus {
        frame: frame(Point3::ZERO, Vec3::Y),
        major_radius: 1.0,
        minor_radius: 0.1,
    });
    let cut = plane(Point3::new(0.0, 0.0, 1.0), Vec3::Z);
    let section = exact_surface_intersection(&torus, &cut).expect("a section");
    let span = section.spans[0].expect("a graph piece carries its span");
    (section.branches[0].clone(), span)
}

fn hits_on(cutter: &Surface) -> Result<Vec<ExactCurveHit>, ExactCurveRefusal> {
    let (curve, span) = half();
    section_curve_surface_intersection(&curve, span, cutter).map(|found| match found {
        ExactCurveIntersection::Points(hits) => hits,
        other => panic!("{other:?}"),
    })
}

fn hits_at(x: f64) -> Vec<ExactCurveHit> {
    hits_on(&plane(Point3::new(x, 0.0, 0.0), Vec3::X)).expect("hits")
}

/// The curvature radius of the half at its tip.
fn osculating_radius() -> f64 {
    0.21f64.sqrt() / 11.0
}

/// The cylinder about `z` of `radius`, through the tip, shifted by `shift`
/// along `x` (positive: into the half's side).
fn cylinder(radius: f64, shift: f64) -> Surface {
    Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(TANGENT - radius + shift, 0.0, 0.0), Vec3::Z),
        radius,
    })
}

fn hits_y(y: f64) -> Result<Vec<ExactCurveHit>, ExactCurveRefusal> {
    hits_on(&plane(Point3::new(0.0, y, 0.0), Vec3::Y))
}

#[test]
fn the_tangent_plane_meets_the_half_once_doubly() {
    let hits = hits_at(TANGENT);
    assert_eq!(hits.len(), 1, "{hits:?}");
    assert_eq!(hits[0].multiplicity, 2);
    let t = hits[0].parameter.approx();
    assert!(t.abs() < 1e-9, "{t}");
    let p = hits[0].point;
    assert!((p.x - TANGENT).abs() < 1e-9 && p.y.abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9);
    let q = evaluate3(&half().0, t).unwrap();
    assert!((q - p).length() < 1e-9);
}

#[test]
fn a_plane_a_little_lower_crosses_twice_symmetrically() {
    let hits = hits_at(TANGENT - 1e-4);
    assert_eq!(hits.len(), 2, "{hits:?}");
    assert!(hits.iter().all(|h| h.multiplicity == 1));
    let (a, b) = (hits[0].point, hits[1].point);
    assert!((a.y + b.y).abs() < 1e-9, "{a:?} {b:?}");
    assert!(a.y.abs() > 1e-4);
    for p in [a, b] {
        assert!((p.x - (TANGENT - 1e-4)).abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9);
    }
    let (s, t) = (hits[0].parameter.approx(), hits[1].parameter.approx());
    assert!((s + t).abs() < 1e-9, "{s} {t}");
}

#[test]
fn a_plane_a_little_higher_misses() {
    assert!(hits_at(TANGENT + 1e-4).is_empty());
}

#[test]
fn a_plane_a_nanometre_either_side_is_not_a_touch() {
    // 1e-9 is within the old absolute touch rule, but the crossings are
    // 18 micrometres apart and the miss is a miss.
    let hits = hits_at(TANGENT - 1e-9);
    assert_eq!(hits.len(), 2, "{hits:?}");
    assert!(hits
        .iter()
        .all(|h| h.multiplicity == 1 && h.point.y.abs() > 9e-6));
    assert!(hits_at(TANGENT + 1e-9).is_empty());
}

#[test]
fn the_osculating_cylinder_touches_to_the_fourth_order() {
    let hits = hits_on(&cylinder(osculating_radius(), 0.0)).expect("hits");
    assert_eq!(hits.len(), 1, "{hits:?}");
    assert_eq!(hits[0].multiplicity, 4);
    let p = hits[0].point;
    assert!(
        (p.x - TANGENT).abs() < 1e-9 && p.y.abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9,
        "{p:?}"
    );
    assert!(hits[0].parameter.approx().abs() < 1e-9);
}

#[test]
fn a_tighter_cylinder_touches_only_to_the_second_order() {
    let hits = hits_on(&cylinder(0.9 * osculating_radius(), 0.0)).expect("hits");
    assert_eq!(hits.len(), 1, "{hits:?}");
    assert_eq!(hits[0].multiplicity, 2);
}

#[test]
fn the_osculating_cylinder_moved_a_micrometre_either_way_misses_or_crosses_twice() {
    let outward = hits_on(&cylinder(osculating_radius(), -1e-6)).expect("hits");
    assert!(outward.is_empty(), "{outward:?}");
    let inward = hits_on(&cylinder(osculating_radius(), 1e-6)).expect("hits");
    assert_eq!(inward.len(), 2, "{inward:?}");
    assert!(inward.iter().all(|h| h.multiplicity == 1));
    let (a, b) = (inward[0].point, inward[1].point);
    assert!((a.y + b.y).abs() < 1e-9, "{a:?} {b:?}");
    assert!(a.y.abs() > 1e-3 && b.y.abs() > 1e-3, "{a:?} {b:?}");
}

#[test]
fn the_plane_tangent_at_the_end_of_the_span_is_decided() {
    // On the half, y = 0.1 - 1.25 x^4 near x = 0, the end of the span: a
    // contact of order 4 at the curve's end, which reports the order it
    // sees there.
    let hits = hits_y(0.1).expect("decided");
    assert_eq!(hits.len(), 1, "{hits:?}");
    let p = hits[0].point;
    assert!(
        p.x.abs() < 1e-9 && (p.y - 0.1).abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9,
        "{p:?}"
    );
    assert!((2..=4).contains(&hits[0].multiplicity), "{hits:?}");
    let lower = hits_y(0.1 - 1e-4).expect("decided");
    assert_eq!(lower.len(), 1, "{lower:?}");
    assert_eq!(lower[0].multiplicity, 1);
    assert!(hits_y(0.1 + 1e-4).expect("decided").is_empty());
}
