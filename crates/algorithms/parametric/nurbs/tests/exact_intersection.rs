//! Exact elementary surface intersections.
//!
//! Every derivation is checked by substituting sampled points of the
//! returned curve back into BOTH surface equations. A test that merely
//! asserted "an ellipse came back" would pass on a wrong ellipse; residual
//! checks against the operands catch a wrong frame, radius, or semi-axis.

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_curve::Curve3;
use axiolid_nurbs::{exact_surface_intersection, Derivation, ExactIntersectionRefusal};
use axiolid_surface::{Cylinder, Plane, Sphere, Surface};

const TAU: f64 = std::f64::consts::TAU;

/// A frame with the given origin and z axis, x and y completed arbitrarily.
fn frame(origin: Point3, z: Vec3) -> Frame3 {
    let z = z.normalize();
    let seed = if z.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let x = z.cross(seed).normalize();
    Frame3 {
        origin,
        x,
        y: z.cross(x),
        z,
    }
}

/// Sample points along a returned curve, in model space.
fn sample(curve: &Curve3, count: usize) -> Vec<Point3> {
    (0..count)
        .map(|index| {
            let t = TAU * (index as f64) / (count as f64);
            match curve {
                Curve3::Circle(circle) => {
                    circle.frame.origin
                        + circle.frame.x * (circle.radius * t.cos())
                        + circle.frame.y * (circle.radius * t.sin())
                }
                Curve3::Ellipse(ellipse) => {
                    ellipse.frame.origin
                        + ellipse.frame.x * (ellipse.semi_axis_x * t.cos())
                        + ellipse.frame.y * (ellipse.semi_axis_y * t.sin())
                }
                Curve3::Line(line) => line.origin + line.direction * (t - TAU / 2.0),
                other => panic!("unexpected curve kind: {other:?}"),
            }
        })
        .collect()
}

/// Distance from a point to a surface, zero exactly on the surface.
///
/// Independent of the derivation under test: written straight from each
/// surface's defining equation so it can contradict a wrong answer.
fn residual(surface: &Surface, point: Point3) -> f64 {
    match surface {
        Surface::Plane(plane) => plane.frame.z.normalize().dot(point - plane.frame.origin),
        Surface::Sphere(sphere) => (point - sphere.frame.origin).length() - sphere.radius,
        Surface::Cylinder(cylinder) => {
            let axis = cylinder.frame.z.normalize();
            let offset = point - cylinder.frame.origin;
            let radial = offset - axis * axis.dot(offset);
            radial.length() - cylinder.radius
        }
        other => panic!("no residual for {other:?}"),
    }
}

/// Assert every sampled point lies on both operands.
fn assert_on_both(first: &Surface, second: &Surface, curve: &Curve3) {
    for point in sample(curve, 24) {
        let first_residual = residual(first, point).abs();
        let second_residual = residual(second, point).abs();
        assert!(
            first_residual < 1.0e-12,
            "point {point:?} off first surface by {first_residual:e}"
        );
        assert!(
            second_residual < 1.0e-12,
            "point {point:?} off second surface by {second_residual:e}"
        );
    }
}

#[test]
fn sphere_plane_gives_a_circle_on_both_surfaces() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 5.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 3.0), Vec3::new(0.0, 0.0, 1.0)),
    });

    let result = exact_surface_intersection(&sphere, &plane).expect("derivable");
    assert_eq!(result.derivation, Derivation::SpherePlaneCircle);

    // r^2 - d^2 = 25 - 9 = 16, so the section radius is exactly 4.
    match result.single() {
        Curve3::Circle(circle) => {
            assert!((circle.radius - 4.0).abs() < 1.0e-12, "{}", circle.radius);
        }
        other => panic!("expected a circle, got {other:?}"),
    }
    assert_on_both(&sphere, &plane, result.single());
}

#[test]
fn a_tangent_plane_is_refused_rather_than_returning_a_degenerate_circle() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    // Plane exactly at the pole: touches at one point.
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 2.0), Vec3::new(0.0, 0.0, 1.0)),
    });

    assert_eq!(
        exact_surface_intersection(&sphere, &plane),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}

#[test]
fn a_missing_plane_is_refused_as_disjoint() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 1.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 4.0), Vec3::new(0.0, 0.0, 1.0)),
    });

    assert_eq!(
        exact_surface_intersection(&sphere, &plane),
        Err(ExactIntersectionRefusal::Disjoint)
    );
}

#[test]
fn a_perpendicular_plane_cuts_a_cylinder_in_a_circle_of_its_own_radius() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 3.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 7.0), Vec3::new(0.0, 0.0, 1.0)),
    });

    let result = exact_surface_intersection(&cylinder, &plane).expect("derivable");
    assert_eq!(
        result.derivation,
        Derivation::CylinderPlanePerpendicularCircle
    );
    match result.single() {
        Curve3::Circle(circle) => {
            assert!((circle.radius - 3.0).abs() < 1.0e-12);
            // The section sits at the plane, not at the cylinder origin.
            assert!((circle.frame.origin.z - 7.0).abs() < 1.0e-12);
        }
        other => panic!("expected a circle, got {other:?}"),
    }
    assert_on_both(&cylinder, &plane, result.single());
}

#[test]
fn a_tilted_plane_cuts_a_cylinder_in_an_ellipse_stretched_by_one_over_cosine() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    // 45 degrees: cos(theta) = 1/sqrt(2), so the major semi-axis is 2*sqrt(2).
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 1.0)),
    });

    let result = exact_surface_intersection(&cylinder, &plane).expect("derivable");
    assert_eq!(result.derivation, Derivation::CylinderPlaneObliqueEllipse);
    match result.single() {
        Curve3::Ellipse(ellipse) => {
            assert!(
                (ellipse.semi_axis_x - 2.0).abs() < 1.0e-12,
                "minor {}",
                ellipse.semi_axis_x
            );
            let expected = 2.0 * std::f64::consts::SQRT_2;
            assert!(
                (ellipse.semi_axis_y - expected).abs() < 1.0e-12,
                "major {} expected {expected}",
                ellipse.semi_axis_y
            );
        }
        other => panic!("expected an ellipse, got {other:?}"),
    }
    assert_on_both(&cylinder, &plane, result.single());
}

#[test]
fn a_plane_parallel_to_the_cylinder_axis_is_refused() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 1.0,
    });
    // Normal perpendicular to the axis: the section is two lines, not one curve.
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
    });

    assert_eq!(
        exact_surface_intersection(&cylinder, &plane),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}

#[test]
fn two_planes_meet_in_a_line_on_both_surfaces() {
    let first = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
    });
    let second = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
    });

    let result = exact_surface_intersection(&first, &second).expect("derivable");
    assert_eq!(result.derivation, Derivation::PlanePlaneLine);
    assert_on_both(&first, &second, result.single());
}

#[test]
fn an_offset_plane_pair_meets_on_a_line_through_both_offsets() {
    // z = 2 and x = -1 meet along the line x = -1, z = 2.
    let first = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 2.0), Vec3::new(0.0, 0.0, 1.0)),
    });
    let second = Surface::Plane(Plane {
        frame: frame(Point3::new(-1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
    });

    let result = exact_surface_intersection(&first, &second).expect("derivable");
    match result.single() {
        Curve3::Line(line) => {
            assert!((line.origin.x + 1.0).abs() < 1.0e-12, "{:?}", line.origin);
            assert!((line.origin.z - 2.0).abs() < 1.0e-12, "{:?}", line.origin);
        }
        other => panic!("expected a line, got {other:?}"),
    }
    assert_on_both(&first, &second, result.single());
}

#[test]
fn parallel_planes_are_refused() {
    let first = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
    });
    let second = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 5.0), Vec3::new(0.0, 0.0, 1.0)),
    });

    assert_eq!(
        exact_surface_intersection(&first, &second),
        Err(ExactIntersectionRefusal::Disjoint)
    );
}

#[test]
fn an_unsupported_pair_is_refused_explicitly_rather_than_approximated() {
    let torus = Surface::Torus(axiolid_surface::Torus {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.5), Vec3::new(0.0, 0.0, 1.0)),
    });

    // A torus section is a real curve, but not a conic this module
    // derives. The refusal names that gap instead of fitting something.
    assert_eq!(
        exact_surface_intersection(&torus, &plane),
        Err(ExactIntersectionRefusal::UnsupportedPair)
    );
}

#[test]
fn argument_order_does_not_change_the_result() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 1.0)),
    });

    let forward = exact_surface_intersection(&cylinder, &plane).expect("derivable");
    let reversed = exact_surface_intersection(&plane, &cylinder).expect("derivable");
    assert_eq!(forward, reversed);
}

#[test]
fn two_spheres_meet_in_a_circle_on_both_surfaces() {
    let first = Surface::Sphere(Sphere {
        frame: frame(Point3::new(1.0, 2.0, -1.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 5.0,
    });
    let second = Surface::Sphere(Sphere {
        frame: frame(Point3::new(4.0, -1.0, 3.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 4.0,
    });

    let result = exact_surface_intersection(&first, &second).expect("derivable");
    assert_eq!(result.derivation, Derivation::SphereSphereCircle);
    assert_on_both(&first, &second, result.single());
}

#[test]
fn a_tangent_sphere_pair_refuses_rather_than_returning_a_point() {
    let first = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 1.0,
    });
    let second = Surface::Sphere(Sphere {
        frame: frame(Point3::new(3.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });

    assert_eq!(
        exact_surface_intersection(&first, &second),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}

#[test]
fn equal_radius_cylinders_on_crossing_axes_cut_two_ellipses() {
    let first = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.3, -0.2, 0.5), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let second = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.3, -0.2, 0.5), Vec3::new(0.0, 1.0, 0.4)),
        radius: 2.0,
    });

    let result = exact_surface_intersection(&first, &second).expect("derivable");
    assert_eq!(
        result.derivation,
        Derivation::CylinderCylinderSteinmetzEllipses
    );
    // BOTH components must come back. Returning one would silently lose
    // half the intersection while still looking like a valid answer.
    assert_eq!(result.branches.len(), 2);
    for branch in &result.branches {
        assert_on_both(&first, &second, branch);
    }

    // The two ellipses must be DISTINCT. Both lie on both cylinders, so an
    // on-surface check alone passes even if the same ellipse is returned
    // twice -- which would silently lose half the intersection while
    // reporting the right branch count.
    let (first_ellipse, second_ellipse) = match (&result.branches[0], &result.branches[1]) {
        (Curve3::Ellipse(a), Curve3::Ellipse(b)) => (a, b),
        other => panic!("expected two ellipses, got {other:?}"),
    };
    let normal_alignment = first_ellipse.frame.z.dot(second_ellipse.frame.z).abs();
    assert!(
        normal_alignment < 0.99,
        "the two Steinmetz planes are nearly identical (|n1.n2| = {normal_alignment})"
    );
}

#[test]
fn unequal_radius_crossing_cylinders_refuse_because_the_curve_is_not_planar() {
    let first = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let second = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
        radius: 1.2,
    });

    // This is a genuine space quartic: sampled points have a third
    // singular value of 2.83, so no plane contains them. There is no
    // exact conic to return, and fitting one would be a lie.
    assert_eq!(
        exact_surface_intersection(&first, &second),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}

#[test]
fn parallel_cylinders_meet_in_two_axis_parallel_lines() {
    let first = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let second = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(1.5, 0.0, 4.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });

    let result = exact_surface_intersection(&first, &second).expect("derivable");
    assert_eq!(result.derivation, Derivation::ParallelCylinderLines);
    assert_eq!(result.branches.len(), 2);
    for branch in &result.branches {
        assert_on_both(&first, &second, branch);
    }
}

#[test]
fn tangent_parallel_cylinders_share_exactly_one_line() {
    let first = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let second = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(4.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });

    let result = exact_surface_intersection(&first, &second).expect("derivable");
    assert_eq!(result.branches.len(), 1);
    assert_on_both(&first, &second, result.single());
}

#[test]
fn identical_cylinders_refuse_rather_than_naming_a_curve() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });

    assert_eq!(
        exact_surface_intersection(&cylinder, &cylinder),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}
