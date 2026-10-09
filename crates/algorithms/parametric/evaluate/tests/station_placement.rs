//! Placing at a station (#264): the rigid placement a section frame stands
//! for maps local `x`, `y`, `z` onto tangent, left lateral and up, a frame
//! carried by a rigid motion moves with it, and only a line's frame is
//! exact. Expected values are closed forms, never the evaluator's own.

use axiolid_core::{Frame2, Point2, Point3, Scalar, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Elevated3, ElevationLaw, Line2, Line3, Polyline2};
use axiolid_evaluate::station::{
    station_frame_is_exact2, station_frame_is_exact3, station_section2, station_section3,
};

const EPS: Scalar = 1e-12;

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?}"
    );
}

#[test]
fn a_placement_maps_x_y_z_onto_tangent_left_lateral_and_up() {
    // A quarter turn into a counter-clockwise R = 5 circle at 0.7 rad:
    // point R (cos, sin, 0), t = (-sin, cos, 0), l = (-cos, -sin, 0)
    // (towards the centre), u = +Z.
    let (radius, angle) = (5.0, 0.7_f64);
    let (s, c) = angle.sin_cos();
    let circle = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    });
    let section = station_section2(&circle, radius * angle).unwrap();
    let placement = section.placement();
    let point = Point3::new(radius * c, radius * s, 0.0);
    let (t, l, u) = (Vec3::new(-s, c, 0.0), Vec3::new(-c, -s, 0.0), Vec3::Z);
    let tight = 1e-9;
    close3(
        placement.transform_point3(Point3::ZERO),
        point,
        tight,
        "origin",
    );
    close3(
        placement.transform_vector3(Vec3::X),
        t,
        tight,
        "x -> tangent",
    );
    close3(
        placement.transform_vector3(Vec3::Y),
        l,
        tight,
        "y -> lateral",
    );
    close3(placement.transform_vector3(Vec3::Z), u, tight, "z -> up");
    close3(
        placement.transform_point3(Point3::new(1.0, 2.0, 3.0)),
        point + t + 2.0 * l + 3.0 * u,
        tight,
        "a local point",
    );
    // Right-handed and rigid.
    assert!((placement.matrix3.determinant() - 1.0).abs() <= EPS);
    // Not the provider layout, whose y is up and z right.
    let frame = section.frame();
    close3(frame.y, u, tight, "provider y is up");
    close3(frame.z, -l, tight, "provider z is right");
}

#[test]
fn a_placement_on_a_grade_tilts_its_up_with_the_section() {
    // 2% grade along +x: t = (1, 0, 0.02) / k, l = +y, u = (-0.02, 0, 1) / k.
    let k = (1.0_f64 + 0.02 * 0.02).sqrt();
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::constant_grade(100.0, 0.02),
    ));
    let placement = station_section3(&curve, 40.0).unwrap().placement();
    close3(
        placement.transform_point3(Point3::ZERO),
        Point3::new(40.0, 0.0, 100.8),
        1e-9,
        "point",
    );
    close3(
        placement.transform_vector3(Vec3::X),
        Vec3::new(1.0, 0.0, 0.02) / k,
        EPS,
        "x",
    );
    close3(placement.transform_vector3(Vec3::Y), Vec3::Y, EPS, "y");
    close3(
        placement.transform_vector3(Vec3::Z),
        Vec3::new(-0.02, 0.0, 1.0) / k,
        EPS,
        "z",
    );
}

#[test]
fn a_moved_frame_follows_the_rigid_motion() {
    let line = Curve2::Line(Line2 {
        origin: Point2::new(1.0, 0.0),
        direction: Vec2::new(0.0, 3.0),
    });
    let section = station_section2(&line, 2.0).unwrap();
    // A quarter turn about +Z, then up by 4.
    let rigid = Transform3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z, Vec3::new(0.0, 0.0, 4.0));
    let moved = section.moved(rigid);
    close3(moved.point, Point3::new(-2.0, 1.0, 4.0), EPS, "point");
    close3(moved.tangent, -Vec3::X, EPS, "tangent");
    close3(moved.lateral, -Vec3::Y, EPS, "lateral");
    close3(moved.up, Vec3::Z, EPS, "up");
}

#[test]
fn only_a_line_frame_is_exact() {
    assert!(station_frame_is_exact2(&Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })));
    assert!(station_frame_is_exact3(&Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    })));
    assert!(!station_frame_is_exact2(&Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 1.0,
    })));
    assert!(!station_frame_is_exact2(&Curve2::Polyline(Polyline2 {
        points: vec![Point2::ZERO, Point2::new(1.0, 0.0)],
        closed: false,
    })));
    assert!(!station_frame_is_exact3(&Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::level(0.0),
    ))));
}
