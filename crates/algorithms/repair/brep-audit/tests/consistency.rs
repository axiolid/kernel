//! The geometric audit must reject inconsistent geometry.
//!
//! The solid here is built by hand rather than through a constructor: no
//! constructor produces this defect any more, and depending on the
//! construction crate would create a dependency cycle.

use axiolid_brep_audit::{geometric_audit, GeometricDefect};
use axiolid_core::Tolerance;

#[test]
fn the_audit_rejects_a_chord_pcurve_on_an_arc_edge() {
    // This is the defect the audit exists for, reproduced deliberately:
    // a cap loop that trims an arc edge with a straight chord. It closes,
    // every handle resolves, and the topological audit is clean -- so if
    // the geometric audit does not fail here it is not doing anything.
    //
    // Built by hand because no constructor produces it any more: the fillet
    // path used to, which is how the bug was found.
    use axiolid_brep::ExactBRepBuilder;
    use axiolid_core::{Frame2, Frame3, Interval, Point3, Vec2, Vec3};
    use axiolid_curve::{Circle3, Curve2, Curve3, Line2, Line3};
    use axiolid_surface::{Cylinder, Surface};
    use axiolid_topology::{audit_brep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};

    let radius = 1.0;
    let sweep = std::f64::consts::FRAC_PI_2;
    let mut builder = ExactBRepBuilder::default();

    let start = Point3::new(radius, 0.0, 0.0);
    let end = Point3::new(0.0, radius, 0.0);
    let v0 = builder
        .topology_mut()
        .add_vertex(Vertex { position: start });
    let v1 = builder.topology_mut().add_vertex(Vertex { position: end });

    // The edge genuinely follows a quarter arc.
    let arc = builder.add_curve3(Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    }));
    let edge = builder.topology_mut().add_edge(Edge {
        start: v0,
        end: v1,
        curve: Some(arc),
    });
    builder.set_edge_interval(edge, Interval::new(0.0, sweep));

    // The pcurve is the straight CHORD between the same endpoints: correct
    // at both ends, wrong everywhere between.
    let chord = builder.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(radius, 0.0),
        direction: Vec2::new(-radius, radius),
    }));
    // Close the loop with a straight return edge, so the ONLY defect is the
    // chord pcurve. An open loop would be caught topologically and would
    // prove nothing about geometric auditing.
    let line_back = builder.add_curve3(Curve3::Line(Line3 {
        origin: end,
        direction: start - end,
    }));
    let back = builder.topology_mut().add_edge(Edge {
        start: v1,
        end: v0,
        curve: Some(line_back),
    });
    builder.set_edge_interval(back, Interval::UNIT);
    let back_pcurve = builder.add_curve2(Curve2::Line(Line2 {
        origin: Vec2::new(0.0, radius),
        direction: Vec2::new(radius, -radius),
    }));

    let loop_id = builder.topology_mut().add_loop(Loop {
        edges: vec![
            EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(chord),
            },
            EdgeUse {
                edge: back,
                orientation: Orientation::Forward,
                pcurve: Some(back_pcurve),
            },
        ],
    });
    builder.set_pcurve_interval(loop_id, 0, Interval::UNIT);
    builder.set_pcurve_interval(loop_id, 1, Interval::UNIT);

    // A planar face: parameter space is the xy-plane, so the chord lifts to
    // itself and the disagreement with the arc is pure sagitta.
    let plane = builder.add_surface(Surface::Plane(axiolid_surface::Plane {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
    }));
    let face = builder.topology_mut().add_face(Face {
        surface: Some(plane),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let _ = face;
    let _ = Frame2 {
        origin: Vec2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    };
    let _ = Cylinder {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    };

    let solid = builder.finish().expect("the builder accepts this solid");

    // Topology is clean: this is exactly why the gap existed.
    let topological = audit_brep(solid.topology());
    assert!(
        topological.is_tessellable(),
        "the topological audit must see nothing wrong -- that is the point"
    );

    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        !health.is_consistent(),
        "a chord standing in for an arc must be rejected"
    );
    assert!(matches!(
        health.defects().first(),
        Some(GeometricDefect::PcurveOffCurve { .. })
    ));
    // Sagitta of a unit quarter arc: r(1 - cos(45 deg)) = 0.2929.
    let worst = health.worst_error().expect("a measured error");
    assert!(
        (worst - 0.292_893_218_813_45).abs() < 1e-9,
        "the reported error must be the sagitta, got {worst}"
    );
}
