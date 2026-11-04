//! A traced section against a surface that touches it (#397): the double
//! root comes back once, with multiplicity 2, not as a pair of spurious
//! roots either side of it.
//!
//! The fixture is #381's: the plane `z = 1` meets the torus about `y` with
//! `R = 1`, `r = 0.1` (the top face of that boolean's box) in two traced
//! halves of a loop. The right half reaches its largest `x`, `sqrt(0.21)`,
//! in the middle of its span, where the plane `x = sqrt(0.21)` is tangent to
//! it; a plane `d` lower crosses it twice and one `d` higher misses it.
//! (The plane `y = 0.1`, tangent at the span's end, is not used: the old
//! sign-change scan alone did not finish on it either. The boolean of #381
//! at `Tolerance::ZERO` covers the certificate end to end.)

use axiolid_core::{Frame3, Interval, Point3, Vec3};
use axiolid_curve::Curve3;
use axiolid_evaluate::evaluate3;
use axiolid_nurbs::{
    exact_surface_intersection, section_curve_surface_intersection, ExactCurveHit,
    ExactCurveIntersection,
};
use axiolid_surface::{Plane, Surface, Torus};

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

fn hits_at(x: f64) -> Vec<ExactCurveHit> {
    let cutter = plane(Point3::new(x, 0.0, 0.0), Vec3::X);
    let (curve, span) = half();
    match section_curve_surface_intersection(&curve, span, &cutter).expect("hits") {
        ExactCurveIntersection::Points(hits) => hits,
        other => panic!("{other:?}"),
    }
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
