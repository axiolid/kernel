//! Rigid placement of an exact B-rep (#223): what it keeps, what a
//! reflection reparameterises, and what it refuses.
//!
//! The B-reps here are one-face sheets whose geometry need not agree; the
//! transform is checked support by support. Placement of real solids, with
//! audits, volumes and distances, is tested in `axiolid-construct`.

use std::f64::consts::{FRAC_PI_2, TAU};

use axiolid_brep::{ExactBRep, ExactBRepBuilder, TransformError, RIGID_TOLERANCE};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Transform3, Vec2, Vec3};
use axiolid_curve::{
    Circle2, Circle3, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2, Line2,
    Line3, Sinusoid2,
};
use axiolid_surface::{Cylinder, Plane, Surface};
use axiolid_topology::{Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};

fn frame() -> Frame3 {
    Frame3 {
        origin: Point3::new(1.0, 2.0, 3.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

/// One face on `surface`: a triangle loop whose three uses carry
/// `pcurve`, over `pcurve_interval`, and whose edges carry `curve`.
fn sheet(surface: Surface, curve: Curve3, pcurve: Curve2, pcurve_interval: Interval) -> ExactBRep {
    let mut builder = ExactBRepBuilder::default();
    let curve3 = builder.add_curve3(curve);
    let curve2 = builder.add_curve2(pcurve);
    let surface = builder.add_surface(surface);
    let topology = builder.topology_mut();
    let vertices: Vec<_> = [Point3::X, Point3::Y, Point3::Z]
        .into_iter()
        .map(|position| topology.add_vertex(Vertex { position }))
        .collect();
    let edges: Vec<_> = (0..3)
        .map(|i| {
            topology.add_edge(Edge {
                start: vertices[i],
                end: vertices[(i + 1) % 3],
                curve: Some(curve3),
            })
        })
        .collect();
    let loop_id = topology.add_loop(Loop {
        edges: edges
            .iter()
            .map(|&edge| EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(curve2),
            })
            .collect(),
    });
    let face = topology.add_face(Face {
        surface: Some(surface),
        bounds: vec![FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    builder.set_face_name(face, axiolid_brep::FaceName::Anonymous(7));
    for edge in edges {
        builder.set_edge_interval(edge, Interval::new(0.25, 0.75));
    }
    for use_index in 0..3 {
        builder.set_pcurve_interval(loop_id, use_index, pcurve_interval);
    }
    builder.finish().expect("a structurally valid sheet")
}

fn plane_sheet(pcurve: Curve2) -> ExactBRep {
    sheet(
        Surface::Plane(Plane { frame: frame() }),
        line3(),
        pcurve,
        Interval::new(0.5, 1.5),
    )
}

fn cylinder_sheet(pcurve: Curve2) -> ExactBRep {
    sheet(
        Surface::Cylinder(Cylinder {
            frame: frame(),
            radius: 2.0,
        }),
        Curve3::Circle(Circle3 {
            frame: frame(),
            radius: 2.0,
        }),
        pcurve,
        Interval::new(0.5, 1.5),
    )
}

fn line3() -> Curve3 {
    Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    })
}

fn line2() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::new(0.5, 1.0),
        direction: Vec2::new(1.0, 2.0),
    })
}

/// A quarter turn about `z` after a quarter turn about `x`, then a shift.
fn placement() -> Transform3 {
    Transform3::from_translation(Vec3::new(10.0, -4.0, 2.5))
        * Transform3::from_rotation_z(FRAC_PI_2)
        * Transform3::from_rotation_x(FRAC_PI_2)
}

/// A reflection in the plane `x = 0`, then a shift.
fn mirror() -> Transform3 {
    Transform3::from_translation(Vec3::new(1.0, 0.0, -2.0))
        * Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0))
}

fn close3(a: Vec3, b: Vec3) {
    assert!((a - b).length() < 1e-12, "{a:?} != {b:?}");
}

fn right_handed(frame: &Frame3) {
    close3(frame.x.cross(frame.y), frame.z);
}

#[test]
fn a_scale_a_shear_and_a_non_finite_transform_are_refused() {
    let sheet = plane_sheet(line2());
    let scale = Transform3::from_scale(Vec3::splat(2.0));
    assert_eq!(sheet.transformed(&scale), Err(TransformError::NotRigid));
    let mut shear = Transform3::IDENTITY;
    shear.matrix3.y_axis.x = 1e-6;
    assert_eq!(sheet.transformed(&shear), Err(TransformError::NotRigid));
    // Just outside the bound is refused; well inside it is rigid.
    let mut nearly = Transform3::IDENTITY;
    nearly.matrix3.x_axis.x = 1.0 + 4.0 * RIGID_TOLERANCE;
    assert_eq!(sheet.transformed(&nearly), Err(TransformError::NotRigid));
    nearly.matrix3.x_axis.x = 1.0 + 0.1 * RIGID_TOLERANCE;
    assert!(sheet.transformed(&nearly).is_ok());
    let mut broken = Transform3::IDENTITY;
    broken.translation.y = f64::NAN;
    assert_eq!(sheet.transformed(&broken), Err(TransformError::NonFinite));
    broken = Transform3::IDENTITY;
    broken.matrix3.z_axis.z = f64::INFINITY;
    assert_eq!(sheet.transformed(&broken), Err(TransformError::NonFinite));
}

#[test]
fn a_rotation_moves_points_and_frames_and_keeps_every_parameter() {
    let source = cylinder_sheet(line2());
    let placement = placement();
    let placed = source.transformed(&placement).expect("rigid");

    for (before, after) in source
        .topology()
        .vertices()
        .iter()
        .zip(placed.topology().vertices())
    {
        close3(placement.transform_point3(before.position), after.position);
    }
    let Surface::Cylinder(cylinder) = &placed.surfaces()[0] else {
        panic!("a cylinder stays a cylinder");
    };
    assert_eq!(cylinder.radius, 2.0);
    close3(
        cylinder.frame.origin,
        placement.transform_point3(frame().origin),
    );
    close3(cylinder.frame.z, placement.transform_vector3(Vec3::Z));
    right_handed(&cylinder.frame);
    // Pcurves, intervals, orientation and names are untouched.
    assert_eq!(placed.curves2(), source.curves2());
    let edge = placed.topology().edge_id_at(0).unwrap();
    assert_eq!(placed.edge_interval(edge), source.edge_interval(edge));
    let wire = placed.topology().loop_id_at(0).unwrap();
    assert_eq!(
        placed.pcurve_interval(wire, 1),
        source.pcurve_interval(wire, 1)
    );
    assert_eq!(placed.topology().faces(), source.topology().faces());
    let face = placed.topology().face_id_at(0).unwrap();
    assert_eq!(placed.face_name(face), source.face_name(face));
}

#[test]
fn a_reflection_flips_every_face_and_keeps_frames_right_handed() {
    let placed = plane_sheet(line2()).transformed(&mirror()).expect("rigid");
    assert_eq!(
        placed.topology().faces()[0].orientation,
        Orientation::Reversed
    );
    let Surface::Plane(plane) = &placed.surfaces()[0] else {
        panic!("a plane stays a plane");
    };
    right_handed(&plane.frame);
    // A plane keeps its parameters: its pcurves are untouched.
    assert_eq!(placed.curves2()[0], line2());
    let Curve3::Line(line) = &placed.curves3()[0] else {
        panic!("a line stays a line");
    };
    close3(line.origin, Vec3::new(1.0, 0.0, -2.0));
    close3(line.direction, -Vec3::X);
    // Twice reflected is the original.
    let back = placed
        .transformed(&mirror().inverse())
        .expect("rigid")
        .topology()
        .faces()[0]
        .orientation;
    assert_eq!(back, Orientation::Forward);
}

/// On a curved face a reflection reads the angle backwards: each pcurve
/// point `(u, v)` lands at `(2 pi - u, v)`, the same surface point.
#[test]
fn a_reflection_reflects_pcurves_on_a_curved_face_in_u() {
    let circle = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(1.0, 0.5),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 0.25,
    });
    let wave = Curve2::Sinusoid(Sinusoid2 {
        mean: 1.0,
        cosine: 0.5,
        sine: 0.25,
    });
    let point = |curve: &Curve2, t: f64| -> Point2 {
        match curve {
            Curve2::Line(l) => l.origin + l.direction * t,
            Curve2::Circle(c) => {
                c.frame.origin + (c.frame.x * t.cos() + c.frame.y * t.sin()) * c.radius
            }
            Curve2::Sinusoid(s) => Point2::new(t, s.height(t)),
            _ => unreachable!(),
        }
    };
    for pcurve in [line2(), circle, wave] {
        let source = cylinder_sheet(pcurve.clone());
        let placed = source.transformed(&mirror()).expect("rigid");
        let wire = placed.topology().loop_id_at(0).unwrap();
        let before = source.pcurve_interval(wire, 0).unwrap();
        let after = placed.pcurve_interval(wire, 0).unwrap();
        for (t0, t1) in [(before.start, after.start), (before.end, after.end)] {
            let p = point(&pcurve, t0);
            let q = point(&placed.curves2()[0], t1);
            assert!(
                (q - Point2::new(TAU - p.x, p.y)).length() < 1e-12,
                "{pcurve:?}: {p:?} must land at (2 pi - u, v), got {q:?}"
            );
        }
        let Surface::Cylinder(cylinder) = &placed.surfaces()[0] else {
            panic!("a cylinder stays a cylinder");
        };
        right_handed(&cylinder.frame);
        // The surface point at (2 pi - u, v) is the mirror of the old one.
        let at = |c: &Cylinder, u: f64, v: f64| {
            c.frame.origin + (c.frame.x * u.cos() + c.frame.y * u.sin()) * c.radius + c.frame.z * v
        };
        let Surface::Cylinder(original) = &source.surfaces()[0] else {
            unreachable!()
        };
        close3(
            at(cylinder, TAU - 0.7, 0.3),
            mirror().transform_point3(at(original, 0.7, 0.3)),
        );
        // A circular edge keeps its parameter: its frame keeps x and y.
        let Curve3::Circle(edge) = &placed.curves3()[0] else {
            panic!("a circle stays a circle");
        };
        right_handed(&edge.frame);
        let edge_id = placed.topology().edge_id_at(0).unwrap();
        assert_eq!(placed.edge_interval(edge_id), source.edge_interval(edge_id));
    }
}

#[test]
fn a_pcurve_with_no_closed_form_reflection_is_refused_on_a_curved_face_only() {
    let spiral = Curve2::Intrinsic(Intrinsic2 {
        start: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        curvature: CurvatureLaw::Constant { curvature: 0.5 },
        length: 1.0,
    });
    assert!(matches!(
        cylinder_sheet(spiral.clone()).transformed(&mirror()),
        Err(TransformError::Unsupported(_))
    ));
    // Rotated, every point keeps its parameters; on a plane it does even
    // when reflected.
    assert!(cylinder_sheet(spiral.clone())
        .transformed(&placement())
        .is_ok());
    assert!(plane_sheet(spiral).transformed(&mirror()).is_ok());
}

#[test]
fn an_elevated_alignment_curve_is_refused_by_name() {
    let elevated = Curve3::Elevated(Elevated3::new(
        line2(),
        ElevationLaw::Polynomial {
            coefficients: vec![1.0],
        },
    ));
    let source = sheet(
        Surface::Plane(Plane { frame: frame() }),
        elevated,
        line2(),
        Interval::UNIT,
    );
    let error = source.transformed(&placement()).unwrap_err();
    assert!(matches!(error, TransformError::Unsupported(_)));
    assert!(error.to_string().contains("elevated"));
}

#[test]
fn a_banked_alignment_curve_is_refused_by_name() {
    let base = Elevated3::new(
        line2(),
        ElevationLaw::Polynomial {
            coefficients: vec![1.0],
        },
    );
    let banked = Curve3::Banked(axiolid_curve::Banked3::new(
        base,
        axiolid_curve::CantLaw::new(vec![axiolid_curve::CantPiece::constant(10.0, 0.1)]),
        axiolid_curve::CantLaw::zero(10.0),
        1.5,
        axiolid_curve::BankConvention::TangentRotation,
    ));
    let source = sheet(
        Surface::Plane(Plane { frame: frame() }),
        banked,
        line2(),
        Interval::UNIT,
    );
    let error = source.transformed(&placement()).unwrap_err();
    assert!(matches!(error, TransformError::Unsupported(_)));
    assert!(error.to_string().contains("banked"));
}
