//! Quadratic-graph pcurves and ruled sections (ADR 0076).
//!
//! Derivatives are checked against central differences of the positions,
//! which share no formula with the implicit differentiation under test; the
//! lifted curve is checked to lie on the carrier surface as evaluated by the
//! surface module, which shares no code with the curve's own lift.

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_curve::{Branch, Curve2, Curve3, QuadraticGraph2, RuledCarrier, RuledSection3, Trig2};
use axiolid_evaluate::curve::{derivative2, evaluate2, second_derivative2};
use axiolid_evaluate::surface::evaluate;
use axiolid_evaluate::{derivative3, evaluate3};
use axiolid_surface::{Cone, Surface};

/// A graph whose discriminant stays positive on `[-1, 1]`.
fn graph(branch: Branch) -> QuadraticGraph2 {
    QuadraticGraph2 {
        a: Trig2 {
            constant: 1.0,
            cos: 0.2,
            ..Trig2::default()
        },
        b: Trig2 {
            constant: 0.5,
            sin: -1.5,
            cos2: 0.25,
            ..Trig2::default()
        },
        c: Trig2 {
            constant: -2.0,
            cos: 0.75,
            sin2: 0.5,
            ..Trig2::default()
        },
        branch,
    }
}

fn frame() -> Frame3 {
    Frame3 {
        origin: Point3::new(1.0, -2.0, 0.5),
        x: Vec3::new(0.6, 0.8, 0.0),
        y: Vec3::new(-0.8, 0.6, 0.0),
        z: Vec3::Z,
    }
}

#[test]
fn every_point_is_a_root_of_its_quadratic() {
    for branch in [Branch::Plus, Branch::Minus] {
        let g = graph(branch);
        for i in 0..=20 {
            let t = -1.0 + 0.1 * i as f64;
            let p = evaluate2(&Curve2::QuadraticGraph(g), t).expect("defined");
            let residual = g.a.value(t) * p.y * p.y + g.b.value(t) * p.y + g.c.value(t);
            assert!(residual.abs() < 1e-12, "{branch:?} at {t}: {residual}");
            assert_eq!(p.x, t, "the parameter is the first coordinate");
        }
    }
    // The two branches are the two roots, each the one its name says.
    let (p, m) = (graph(Branch::Plus), graph(Branch::Minus));
    let t = 0.3;
    let (a, b, c) = (p.a.value(t), p.b.value(t), p.c.value(t));
    let root = (b * b - 4.0 * a * c).sqrt();
    assert!((p.height(t).unwrap() - (-b + root) / (2.0 * a)).abs() < 1e-12);
    assert!((m.height(t).unwrap() - (-b - root) / (2.0 * a)).abs() < 1e-12);
}

#[test]
fn derivatives_match_central_differences() {
    let h = 1e-5;
    for branch in [Branch::Plus, Branch::Minus] {
        let curve = Curve2::QuadraticGraph(graph(branch));
        for t in [-0.8, -0.2, 0.0, 0.4, 0.9] {
            let d = derivative2(&curve, t).unwrap();
            let (a, b) = (
                evaluate2(&curve, t - h).unwrap(),
                evaluate2(&curve, t + h).unwrap(),
            );
            let fd = (b - a) / (2.0 * h);
            assert!(
                (d - fd).length() < 1e-6,
                "{branch:?} first at {t}: {d:?} vs {fd:?}"
            );
            let dd = second_derivative2(&curve, t).unwrap();
            let (da, db) = (
                derivative2(&curve, t - h).unwrap(),
                derivative2(&curve, t + h).unwrap(),
            );
            let fdd = (db - da) / (2.0 * h);
            assert!(
                (dd - fdd).length() < 1e-5,
                "{branch:?} second at {t}: {dd:?} vs {fdd:?}"
            );
        }
    }
}

#[test]
fn a_ruled_section_lies_on_its_carrier_and_differentiates_consistently() {
    let cone = Cone {
        frame: frame(),
        radius: 1.5,
        semi_angle: 0.3,
    };
    let carrier = RuledCarrier {
        frame: cone.frame,
        x_radius: cone.radius,
        y_radius: cone.radius,
        slope: cone.semi_angle.tan(),
    };
    let section = Curve3::RuledSection(RuledSection3 {
        carrier,
        graph: graph(Branch::Plus),
    });
    let surface = Surface::Cone(cone);
    let h = 1e-5;
    for t in [-0.7, 0.1, 0.8] {
        let p = evaluate3(&section, t).unwrap();
        let v = graph(Branch::Plus).height(t).unwrap();
        let on = evaluate(&surface, t, v).unwrap();
        assert!(
            (p - on).length() < 1e-12,
            "lift disagrees with the cone at {t}"
        );
        let d = derivative3(&section, t).unwrap();
        let fd =
            (evaluate3(&section, t + h).unwrap() - evaluate3(&section, t - h).unwrap()) / (2.0 * h);
        assert!((d - fd).length() < 1e-6, "tangent at {t}: {d:?} vs {fd:?}");
    }
}

#[test]
fn outside_its_spans_the_graph_is_refused_not_extrapolated() {
    // D = b^2 - 4ac < 0 everywhere: a = 1, b = 0, c = 1.
    let empty = QuadraticGraph2 {
        a: Trig2 {
            constant: 1.0,
            ..Trig2::default()
        },
        b: Trig2::default(),
        c: Trig2 {
            constant: 1.0,
            ..Trig2::default()
        },
        branch: Branch::Plus,
    };
    assert!(evaluate2(&Curve2::QuadraticGraph(empty), 0.3).is_err());
}
