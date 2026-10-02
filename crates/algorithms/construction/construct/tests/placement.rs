//! Exact B-reps under rigid placements, and exact swept disks (#223).
//!
//! Every expected value is a closed form from the inputs. A body is built
//! in its own frame, then tilted, turned and shifted; its certified
//! boundary distance to a slab whose top is the plane `z = 0` is then its
//! lowest point's height, which the tilt gives directly:
//!
//! - a disk of radius `r` in a plane tilted by `a` from horizontal reaches
//!   `r sin a` below its centre;
//! - an ellipse with semi-axis `b` along the tilt axis' normal reaches
//!   `b sin a` below;
//! - a torus whose axis is `a` from vertical has its tube circle reach
//!   `R sin a` below the centre, and the tube `r` further.
//!
//! Each tilt is chosen so the lowest point is isolated, and every solid is
//! audited and measured: a placement that broke a pcurve or an orientation
//! fails the audit or flips the volume's sign.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_brep::{ExactBRep, TransformError};
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_construct::swept_disk_exact::{
    swept_disk_along_arc_exact, swept_disk_along_line_exact,
};
use axiolid_contracts::GeomError;
use axiolid_core::{Frame3, Interval, Point3, Tolerance, Transform3, Vec3};
use axiolid_curve::Circle3;
use axiolid_measure::{boundary_distance, exact_properties, DistanceBounds};
use axiolid_profile::{CircleProfile, EllipseProfile, Profile, RectangleProfile};

fn tol() -> Tolerance {
    Tolerance::METRE
}

/// Audit clean, closed, and return the volume.
fn checked(solid: &ExactBRep) -> f64 {
    let health = geometric_audit(solid, tol());
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(solid.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    exact_properties(solid, tol())
        .expect("measurable")
        .signed_volume
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

/// A 40 x 40 slab whose top face is the plane `z = 0`.
fn slab() -> ExactBRep {
    let block = extrude_profile_exact(
        &Profile::Rectangle(RectangleProfile {
            x: 40.0,
            y: 40.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Vec3::Z,
        1.0,
        tol(),
    )
    .expect("a block");
    block
        .transformed(&Transform3::from_translation(Vec3::new(0.0, 0.0, -1.0)))
        .expect("a translation is rigid")
}

/// The certified distance to the slab, which must contain `expected` and
/// close to `accuracy`.
fn distance_to_slab(solid: &ExactBRep, expected: f64, accuracy: f64) -> DistanceBounds {
    let bounds = boundary_distance(solid, &slab(), accuracy, tol()).expect("bounded");
    assert!(
        bounds.lower <= expected + 1e-12 && expected <= bounds.upper + 1e-12,
        "[{}, {}] must contain {expected}",
        bounds.lower,
        bounds.upper
    );
    assert!(
        bounds.upper - bounds.lower <= accuracy,
        "[{}, {}] did not close to {accuracy}",
        bounds.lower,
        bounds.upper
    );
    bounds
}

/// Tilt by `tilt` about x, turn by `turn` about z, then shift to `at`.
fn placement(tilt: f64, turn: f64, at: Vec3) -> Transform3 {
    Transform3::from_translation(at)
        * Transform3::from_rotation_z(turn)
        * Transform3::from_rotation_x(tilt)
}

/// A reflection in the plane `x = 0`.
fn mirror() -> Transform3 {
    Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0))
}

fn circle_column(radius: f64, height: f64) -> ExactBRep {
    extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        }),
        Vec3::Z,
        height,
        tol(),
    )
    .expect("a column")
}

fn ellipse_column(a: f64, b: f64, height: f64) -> ExactBRep {
    extrude_profile_exact(
        &Profile::Ellipse(EllipseProfile {
            semi_axis_x: a,
            semi_axis_y: b,
        }),
        Vec3::Z,
        height,
        tol(),
    )
    .expect("an elliptical column")
}

/// A torus of major radius `major` and tube `minor`, centred at the origin
/// with its axis along y: a full turn of a disk.
fn torus(major: f64, minor: f64) -> ExactBRep {
    let local = revolve_profile_exact(
        &Profile::Circle(CircleProfile {
            radius: minor,
            thickness: None,
        }),
        Point3::new(-major, 0.0, 0.0),
        Vec3::Y,
        TAU,
        tol(),
    )
    .expect("a torus");
    local
        .transformed(&Transform3::from_translation(Vec3::new(major, 0.0, 0.0)))
        .expect("rigid")
}

#[test]
fn a_tilted_circular_column_is_its_lowest_rim_point_above_the_slab() {
    let (r, h, tilt) = (0.3, 3.0, 0.4);
    let column = circle_column(r, h)
        .transformed(&placement(tilt, 0.7, Vec3::new(2.0, -1.0, 1.5)))
        .expect("rigid");
    close("volume", checked(&column), PI * r * r * h);
    let bounds = distance_to_slab(&column, 1.5 - r * tilt.sin(), 1e-9);
    assert!(
        bounds.point_b.z.abs() < 1e-9,
        "the slab witness is on its top"
    );
}

#[test]
fn a_tilted_elliptical_column_is_its_lowest_rim_point_above_the_slab() {
    // Semi-axis b lies along local y, which the tilt about x lifts; then
    // the same column tilted about y lowers by its semi-axis a instead.
    let (a, b, h, tilt) = (0.5, 0.2, 2.0, 0.6);
    let about_x = ellipse_column(a, b, h)
        .transformed(&placement(tilt, -1.1, Vec3::new(-3.0, 0.5, 1.2)))
        .expect("rigid");
    close("volume", checked(&about_x), PI * a * b * h);
    distance_to_slab(&about_x, 1.2 - b * tilt.sin(), 1e-9);

    let about_y = ellipse_column(a, b, h)
        .transformed(
            &(Transform3::from_translation(Vec3::new(1.0, 1.0, 0.9))
                * Transform3::from_rotation_y(tilt)),
        )
        .expect("rigid");
    distance_to_slab(&about_y, 0.9 - a * tilt.sin(), 1e-9);
}

#[test]
fn a_lying_column_rests_its_whole_side_at_its_axis_height_less_the_radius() {
    // A horizontal extrusion: every slice along the column is a nearest
    // pair, the slow case, so the accuracy asked is looser.
    let (r, h) = (0.25, 2.0);
    let column = circle_column(r, h)
        .transformed(&placement(FRAC_PI_2, 0.3, Vec3::new(0.0, 0.0, 1.0)))
        .expect("rigid");
    close("volume", checked(&column), PI * r * r * h);
    distance_to_slab(&column, 1.0 - r, 1e-6);
}

#[test]
fn a_full_turn_revolution_under_a_tilt_is_its_tube_bottom_above_the_slab() {
    // The axis starts along y; turning it by `tilt` about x leaves it
    // `pi/2 - tilt` from vertical, so the tube circle dips `R cos(tilt)`.
    let (major, minor, tilt) = (1.0, 0.2, 0.5);
    let solid = torus(major, minor)
        .transformed(&placement(tilt, 0.9, Vec3::new(0.5, -0.5, 2.0)))
        .expect("rigid");
    close(
        "volume",
        checked(&solid),
        2.0 * PI * PI * major * minor * minor,
    );
    distance_to_slab(&solid, 2.0 - major * tilt.cos() - minor, 1e-9);
}

#[test]
fn a_reflected_solid_stays_outward_and_measures_the_same() {
    let (r, h, tilt) = (0.3, 3.0, 0.4);
    let placed = placement(tilt, 0.7, Vec3::new(2.0, -1.0, 1.5));
    let column = circle_column(r, h)
        .transformed(&(mirror() * placed))
        .expect("rigid");
    close("volume", checked(&column), PI * r * r * h);
    distance_to_slab(&column, 1.5 - r * tilt.sin(), 1e-9);

    let (a, b) = (0.5, 0.2);
    let ellipse = ellipse_column(a, b, h)
        .transformed(&(placed * mirror()))
        .expect("rigid");
    close("volume", checked(&ellipse), PI * a * b * h);
    distance_to_slab(&ellipse, 1.5 - b * tilt.sin(), 1e-9);

    let (major, minor) = (1.0, 0.2);
    let ring = torus(major, minor)
        .transformed(&(mirror() * placement(0.5, 0.9, Vec3::new(0.5, -0.5, 2.0))))
        .expect("rigid");
    close(
        "volume",
        checked(&ring),
        2.0 * PI * PI * major * minor * minor,
    );
    distance_to_slab(&ring, 2.0 - major * 0.5_f64.cos() - minor, 1e-9);

    // A partial turn adds planar walls and arc edges.
    let wedge = revolve_profile_exact(
        &Profile::Rectangle(RectangleProfile {
            x: 0.4,
            y: 0.6,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Point3::new(-1.0, 0.0, 0.0),
        Vec3::Y,
        2.0,
        tol(),
    )
    .expect("a wedge");
    let reflected = wedge
        .transformed(&(mirror() * placement(0.3, 0.2, Vec3::new(0.0, 0.0, 3.0))))
        .expect("rigid");
    close("volume", checked(&reflected), 2.0 * 1.0 * 0.4 * 0.6);
}

#[test]
fn a_scaled_or_sheared_placement_is_refused_not_approximated() {
    let column = circle_column(0.3, 1.0);
    assert_eq!(
        column.transformed(&Transform3::from_scale(Vec3::new(1.0, 2.0, 1.0))),
        Err(TransformError::NotRigid)
    );
    assert_eq!(
        column.transformed(&Transform3::from_scale(Vec3::splat(1.001))),
        Err(TransformError::NotRigid)
    );
}

#[test]
fn a_disk_swept_along_a_tilted_segment_is_a_capped_cylinder() {
    let (start, end) = (Point3::new(0.0, 0.0, 3.0), Point3::new(1.0, 2.0, 1.5));
    let r = 0.2;
    let length = (end - start).length();
    let solid = swept_disk_along_line_exact(start, end, r, None, tol()).expect("exact");
    close("volume", checked(&solid), PI * r * r * length);
    // The lowest point is on the end cap's rim, `r sqrt(1 - d_z^2)` below
    // the end point.
    let d = (end - start) / length;
    distance_to_slab(&solid, 1.5 - r * (1.0 - d.z * d.z).sqrt(), 1e-9);

    let ri = 0.12;
    let pipe = swept_disk_along_line_exact(start, end, r, Some(ri), tol()).expect("exact");
    close(
        "hollow volume",
        checked(&pipe),
        PI * (r * r - ri * ri) * length,
    );
}

#[test]
fn a_disk_swept_along_an_arc_is_a_capped_torus_wedge() {
    // An arc in the vertical plane y = 0.5, centre at height 3, running
    // through its lowest point at t = 3 pi / 2.
    let arc = Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.2, 0.5, 3.0),
            x: Vec3::X,
            y: Vec3::Z,
            z: -Vec3::Y,
        },
        radius: 2.0,
    };
    let (r, span) = (0.25, Interval::new(PI, 2.0 * PI - 0.2));
    let angle = span.end - span.start;
    let solid = swept_disk_along_arc_exact(&arc, span, r, None, tol()).expect("exact");
    close("volume", checked(&solid), PI * r * r * 2.0 * angle);
    distance_to_slab(&solid, 3.0 - 2.0 - r, 1e-9);

    // Run backwards over the same arc: the same solid.
    let reversed =
        swept_disk_along_arc_exact(&arc, Interval::new(span.end, span.start), r, None, tol())
            .expect("exact");
    close(
        "reversed volume",
        checked(&reversed),
        PI * r * r * 2.0 * angle,
    );
    distance_to_slab(&reversed, 3.0 - 2.0 - r, 1e-9);

    // A hollow disk round a whole circle: a torus with a toroidal bore.
    let ri = 0.1;
    let ring = swept_disk_along_arc_exact(&arc, Interval::new(0.0, TAU), r, Some(ri), tol())
        .expect("exact");
    close(
        "ring volume",
        checked(&ring),
        PI * (r * r - ri * ri) * TAU * 2.0,
    );
    distance_to_slab(&ring, 3.0 - 2.0 - r, 1e-9);
}

#[test]
fn the_swept_disk_endpoints_follow_the_directrix() {
    // The caps are centred on the directrix ends: the solid's extent
    // along the arc's start radial reaches exactly R + r there.
    let arc = Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: 3.0,
    };
    // A quarter from t = pi/2: the solid lies in x <= 0, y >= 0.
    let solid = swept_disk_along_arc_exact(&arc, Interval::new(FRAC_PI_2, PI), 0.5, None, tol())
        .expect("exact");
    for vertex in solid.topology().vertices() {
        let p = vertex.position;
        assert!(p.x <= 1e-12 && p.y >= -1e-12, "{p:?} outside the quadrant");
        assert!(p.z.abs() <= 0.5 + 1e-12);
    }
}

#[test]
fn a_swept_disk_that_cannot_be_exact_is_refused() {
    let arc = Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: 1.0,
    };
    let unsupported = |result: Result<ExactBRep, GeomError>| {
        assert!(
            matches!(result, Err(GeomError::UnsupportedInput { .. })),
            "{result:?}"
        );
    };
    let invalid = |result: Result<ExactBRep, GeomError>| {
        assert!(
            matches!(result, Err(GeomError::InvalidInput(_))),
            "{result:?}"
        );
    };
    // The disk reaches the axis.
    unsupported(swept_disk_along_arc_exact(
        &arc,
        Interval::new(0.0, 1.0),
        1.0,
        None,
        tol(),
    ));
    // Beyond a full turn.
    unsupported(swept_disk_along_arc_exact(
        &arc,
        Interval::new(0.0, 7.0),
        0.2,
        None,
        tol(),
    ));
    // A bore as wide as the disk, a zero radius, a zero span, a point.
    invalid(swept_disk_along_arc_exact(
        &arc,
        Interval::new(0.0, 1.0),
        0.2,
        Some(0.2),
        tol(),
    ));
    invalid(swept_disk_along_arc_exact(
        &arc,
        Interval::new(1.0, 1.0),
        0.2,
        None,
        tol(),
    ));
    invalid(swept_disk_along_line_exact(
        Point3::ZERO,
        Point3::X,
        0.0,
        None,
        tol(),
    ));
    invalid(swept_disk_along_line_exact(
        Point3::X,
        Point3::X,
        0.1,
        None,
        tol(),
    ));
    // A frame that is not orthonormal.
    let skewed = Circle3 {
        frame: Frame3 {
            x: Vec3::new(1.0, 0.1, 0.0),
            ..arc.frame
        },
        radius: 1.0,
    };
    invalid(swept_disk_along_arc_exact(
        &skewed,
        Interval::new(0.0, 1.0),
        0.2,
        None,
        tol(),
    ));
}
