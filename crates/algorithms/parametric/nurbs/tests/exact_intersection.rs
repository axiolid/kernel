//! Exact elementary surface intersections.
//!
//! Every derivation is checked by substituting sampled points of the
//! returned curve back into BOTH surface equations. A test that merely
//! asserted "an ellipse came back" would pass on a wrong ellipse; residual
//! checks against the operands catch a wrong frame, radius, or semi-axis.

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_curve::Curve3;
use axiolid_nurbs::{exact_surface_intersection, Derivation, ExactIntersectionRefusal};
use axiolid_surface::{Cone, Cylinder, Plane, Sphere, Surface, Torus};

const TAU: f64 = std::f64::consts::TAU;

/// A frame with the given origin and z axis, x and y completed arbitrarily.
/// Build a frame whose `z` keeps the caller's length.
///
/// The ordinary `frame` helper normalises, which would make a test about
/// non-unit axes silently vacuous.
fn frame_keeping_axis_length(origin: Point3, z: Vec3) -> Frame3 {
    let unit = z.normalize();
    let seed = if unit.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let x = unit.cross(seed).normalize();
    Frame3 {
        origin,
        x,
        y: unit.cross(x),
        z,
    }
}

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
        Surface::Cone(cone) => {
            // rho = radius + z * tan(semi_angle), in cone-local coordinates.
            let axis = cone.frame.z.normalize();
            let offset = point - cone.frame.origin;
            let height = axis.dot(offset);
            let radial = (offset - axis * height).length();
            radial - (cone.radius + height * cone.semi_angle.tan())
        }
        Surface::Torus(torus) => {
            // (rho - major)^2 + z^2 = minor^2, in torus-local coordinates.
            let axis = torus.frame.z.normalize();
            let offset = point - torus.frame.origin;
            let height = axis.dot(offset);
            let radial = (offset - axis * height).length();
            let planar = radial - torus.major_radius;
            planar * planar + height * height - torus.minor_radius * torus.minor_radius
        }
        other => panic!("no residual for {other:?}"),
    }
}

/// Assert every sampled point lies on both operands.
///
/// The bound is relative to how far the sampled point sits from the origin.
/// An absolute bound is itself unit-dependent: the same shape modelled in
/// millimetres carries a thousand times more magnitude, so its rounding is
/// a thousand times coarser and an absolute threshold fails for reasons
/// that have nothing to do with the derivation being tested.
fn assert_on_both(first: &Surface, second: &Surface, curve: &Curve3) {
    for point in sample(curve, 24) {
        let magnitude = point.length().max(1.0);
        let bound = 1.0e-12 * magnitude;
        let first_residual = residual(first, point).abs();
        let second_residual = residual(second, point).abs();
        assert!(
            first_residual < bound,
            "point {point:?} off first surface by {first_residual:e}"
        );
        assert!(
            second_residual < bound,
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
fn a_plane_parallel_to_the_cylinder_axis_cuts_two_rulings() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    // Normal perpendicular to the axis, offset 1.2 from it: the section is
    // the two rulings at half-chord sqrt(4 - 1.44) = 1.6.
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(1.2, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
    });

    let result = exact_surface_intersection(&cylinder, &plane).expect("derivable");
    assert_eq!(result.derivation, Derivation::CylinderPlaneParallelRulings);
    assert_eq!(result.branches.len(), 2);
    for branch in &result.branches {
        assert_on_both(&cylinder, &plane, branch);
    }
    // The two rulings must be DISTINCT: returning one line twice would
    // satisfy every on-surface check while losing half the section.
    let origins: Vec<Point3> = result
        .branches
        .iter()
        .map(|branch| match branch {
            Curve3::Line(line) => line.origin,
            other => panic!("expected lines, got {other:?}"),
        })
        .collect();
    assert!((origins[0] - origins[1]).length() > 1.0e-9);
}

#[test]
fn a_plane_tangent_to_the_cylinder_shares_exactly_one_ruling() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(2.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
    });

    let result = exact_surface_intersection(&cylinder, &plane).expect("derivable");
    assert_eq!(result.branches.len(), 1);
    assert_on_both(&cylinder, &plane, result.single());
}

#[test]
fn a_plane_clear_of_the_cylinder_is_disjoint() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(3.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
    });

    assert_eq!(
        exact_surface_intersection(&cylinder, &plane),
        Err(ExactIntersectionRefusal::Disjoint)
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
fn a_spline_surface_pair_is_refused_explicitly_rather_than_approximated() {
    // Two offset tori: both are surfaces of revolution, but the axes do
    // not coincide, so no closed-form circle family exists. Previously
    // coaxial torus/plane stood here; that case is now derived, so this
    // guards a pair that genuinely remains outside the closed forms.
    let first = Surface::Torus(Torus {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let second = Surface::Torus(Torus {
        frame: frame(Point3::new(2.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        major_radius: 3.0,
        minor_radius: 1.0,
    });

    assert_eq!(
        exact_surface_intersection(&first, &second),
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
fn unequal_radius_crossing_cylinders_meet_in_a_ruled_section() {
    let first = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 2.0,
    });
    let second = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
        radius: 1.2,
    });

    // A genuine space quartic: no plane contains it, so no conic is
    // returned. It is exact as a root branch of a quadratic over the
    // carrier's angle instead (ADR 0076); `ruled_section.rs` checks the
    // points. The thinner pipe pierces the wider: two closed loops.
    let curve = exact_surface_intersection(&first, &second).expect("a ruled section");
    assert_eq!(curve.derivation, Derivation::RuledQuadricSection);
    assert_eq!(
        curve.branches.len(),
        4,
        "two loops, each a plus and a minus piece"
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

#[test]
fn a_coaxial_sphere_and_cylinder_meet_in_two_circles() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 5.0,
    });
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 3.0,
    });

    let result = exact_surface_intersection(&sphere, &cylinder).expect("derivable");
    assert_eq!(result.derivation, Derivation::CoaxialRevolutionCircles);
    assert_eq!(
        result.branches.len(),
        2,
        "a sphere is cut twice by a cylinder"
    );

    // 3^2 + z^2 = 5^2 gives z = -4 and z = +4, radius 3 on both.
    let mut heights = Vec::new();
    for branch in &result.branches {
        match branch {
            Curve3::Circle(circle) => {
                assert!((circle.radius - 3.0).abs() < 1.0e-12);
                heights.push(circle.frame.origin.z);
            }
            other => panic!("expected circles, got {other:?}"),
        }
        assert_on_both(&sphere, &cylinder, branch);
    }
    heights.sort_by(f64::total_cmp);
    assert!((heights[0] + 4.0).abs() < 1.0e-12, "lower circle at z=-4");
    assert!((heights[1] - 4.0).abs() < 1.0e-12, "upper circle at z=+4");
}

#[test]
fn a_coaxial_torus_and_plane_meet_in_two_circles() {
    let torus = Surface::Torus(Torus {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        major_radius: 4.0,
        minor_radius: 1.5,
    });
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.9), Vec3::new(0.0, 0.0, 1.0)),
    });

    let result = exact_surface_intersection(&torus, &plane).expect("derivable");
    assert_eq!(result.derivation, Derivation::CoaxialRevolutionCircles);
    assert_eq!(
        result.branches.len(),
        2,
        "a plane cuts the tube inner and outer"
    );

    // rho = 4 +- sqrt(1.5^2 - 0.9^2) = 4 +- 1.2
    let mut radii = Vec::new();
    for branch in &result.branches {
        match branch {
            Curve3::Circle(circle) => {
                assert!((circle.frame.origin.z - 0.9).abs() < 1.0e-12);
                radii.push(circle.radius);
            }
            other => panic!("expected circles, got {other:?}"),
        }
        assert_on_both(&torus, &plane, branch);
    }
    radii.sort_by(f64::total_cmp);
    assert!((radii[0] - 2.8).abs() < 1.0e-12, "inner radius 4 - 1.2");
    assert!((radii[1] - 5.2).abs() < 1.0e-12, "outer radius 4 + 1.2");
}

#[test]
fn a_coaxial_cone_and_sphere_meet_in_a_circle() {
    let cone = Surface::Cone(Cone {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 1.0,
        semi_angle: std::f64::consts::FRAC_PI_6,
    });
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 3.0,
    });

    let result = exact_surface_intersection(&cone, &sphere).expect("derivable");
    assert_eq!(result.derivation, Derivation::CoaxialRevolutionCircles);
    // The other profile root has rho < 0, which generates no circle.
    assert_eq!(
        result.branches.len(),
        1,
        "only the positive-radius root lifts"
    );
    assert_on_both(&cone, &sphere, &result.branches[0]);
}

#[test]
fn a_coaxial_torus_and_sphere_meet_in_two_circles() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 5.0,
    });
    let torus = Surface::Torus(Torus {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        major_radius: 4.0,
        minor_radius: 1.5,
    });

    let result = exact_surface_intersection(&sphere, &torus).expect("derivable");
    assert_eq!(result.branches.len(), 2);
    for branch in &result.branches {
        assert_on_both(&sphere, &torus, branch);
    }
}

#[test]
fn an_offset_sphere_and_cylinder_meet_in_a_ruled_section_not_a_circle() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 5.0,
    });
    // Axis shifted sideways: the shared rotational symmetry is gone and
    // the intersection is a space quartic, not a circle. It is returned
    // as a ruled section on the cylinder, never as a plausible circle.
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 3.0,
    });
    let curve = exact_surface_intersection(&sphere, &cylinder).expect("a ruled section");
    assert_eq!(curve.derivation, Derivation::RuledQuadricSection);
    assert!(curve
        .branches
        .iter()
        .all(|branch| matches!(branch, Curve3::RuledSection(_))));
}

#[test]
fn a_hyperbolic_cone_section_is_refused_as_unrepresentable() {
    let cone = Surface::Cone(Cone {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 1.0,
        semi_angle: std::f64::consts::FRAC_PI_6,
    });
    // Plane containing the axis direction: the section is a hyperbola,
    // and Curve3 has no hyperbola variant.
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)),
    });

    assert_eq!(
        exact_surface_intersection(&cone, &plane),
        Err(ExactIntersectionRefusal::UnrepresentableConic),
        "a parabola or hyperbola must be named, not swapped for an ellipse"
    );
}

/// The same shape gets the same verdict in metres, millimetres and
/// kilometres.
///
/// An absolute coaxiality threshold silently changes meaning with the
/// modelling unit: a pair accepted in metres was refused in millimetres,
/// which is the unit most building models are authored in. The failure
/// was invisible -- a refusal, not a wrong curve -- so this pins the
/// verdict rather than any particular coordinate.
#[test]
fn the_coaxiality_verdict_does_not_depend_on_modelling_units() {
    for scale in [1.0_f64, 1000.0, 0.001] {
        let radius = 3.0 * scale;
        // A lateral offset fixed at 1e-13 OF THE RADIUS: the same shape
        // every time, only the stored numbers differ.
        let offset = radius * 1.0e-13;
        let sphere = Surface::Sphere(Sphere {
            frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
            radius: 5.0 * scale,
        });
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame(Point3::new(offset, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
            radius,
        });

        let result = exact_surface_intersection(&sphere, &cylinder)
            .unwrap_or_else(|error| panic!("scale {scale} refused: {error:?}"));
        assert_eq!(result.branches.len(), 2, "scale {scale}");
        for branch in &result.branches {
            assert_on_both(&sphere, &cylinder, branch);
        }
    }
}

/// A frame whose `z` is not unit length describes the same surface.
///
/// Nothing in the type system forces `Frame3::z` to be normalised, and
/// the coaxial derivation assumed it was: a frame storing a doubled axis
/// was refused outright, though it denotes exactly the same geometry.
#[test]
fn a_non_unit_frame_axis_describes_the_same_surface() {
    let mut results = Vec::new();
    // The two operands carry DIFFERENTLY scaled axes. Using one length for
    // both would leave the parallel test comparing two equally scaled
    // vectors, which hides a missing normalisation on either side.
    for (first_length, second_length) in [(1.0_f64, 1.0_f64), (2.0, 1.0), (1.0, 0.25)] {
        let sphere = Surface::Sphere(Sphere {
            frame: frame_keeping_axis_length(
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, first_length),
            ),
            radius: 5.0,
        });
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame_keeping_axis_length(
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, second_length),
            ),
            radius: 3.0,
        });

        let result = exact_surface_intersection(&sphere, &cylinder).unwrap_or_else(|error| {
            panic!("axis lengths {first_length}/{second_length} refused: {error:?}")
        });
        for branch in &result.branches {
            assert_on_both(&sphere, &cylinder, branch);
        }

        // Only the FIRST operand's axis drives the derivation, so the
        // reversed order must be checked too: otherwise a missing
        // normalisation on the second operand never shows up.
        let reversed = exact_surface_intersection(&cylinder, &sphere).unwrap_or_else(|error| {
            panic!("reversed {first_length}/{second_length} refused: {error:?}")
        });
        for branch in &reversed.branches {
            assert_on_both(&sphere, &cylinder, branch);
        }

        results.push(result);
    }

    // Same geometry in, same curves out -- not merely 'also derivable'.
    assert_eq!(results[0], results[1]);
    assert_eq!(results[0], results[2]);
}

/// A zero-length frame axis is refused, not divided by.
///
/// `Frame3` cannot express 'this axis is valid', so a degenerate frame
/// reaches the derivation like any other. Normalising it would divide by
/// zero and hand back a curve built from NaNs, which is far worse than a
/// refusal because it looks like an answer.
#[test]
fn a_degenerate_frame_axis_is_refused() {
    let zero = Vec3::new(0.0, 0.0, 0.0);
    let sphere = Surface::Sphere(Sphere {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::new(1.0, 0.0, 0.0),
            y: Vec3::new(0.0, 1.0, 0.0),
            z: zero,
        },
        radius: 5.0,
    });
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
        radius: 3.0,
    });

    let result = exact_surface_intersection(&sphere, &cylinder);
    assert!(result.is_err(), "degenerate axis produced {result:?}");

    // And nothing NaN-shaped escaped in the other order either.
    let reversed = exact_surface_intersection(&cylinder, &sphere);
    if let Ok(curve) = reversed {
        for branch in &curve.branches {
            if let Curve3::Circle(circle) = branch {
                assert!(circle.radius.is_finite(), "NaN radius escaped");
                assert!(circle.frame.origin.is_finite(), "NaN centre escaped");
            }
        }
    }
}

// --- exact decisions (#119) ---------------------------------------------------
//
// Each case below holds EXACTLY for the given doubles, and each was
// misjudged by the earlier `f64` comparisons (measured before the fix).
// Axes are kept at their integer lengths so the exact claim is checkable.

/// A plane exactly tangent to a sphere along an oblique normal.
///
/// Normal (3, 2, 6) has length 7; the plane through 7 * (3, 2, 6) / 7 * r
/// ... is built directly: origin (3, 2, 6), normal (3, 2, 6), radius 7.
/// `r^2 |n|^2 = (n . o)^2` holds exactly. In `f64` the unit-normal
/// distance came out below 7 and the pair was reported as a circle of
/// radius about 1e-7.
#[test]
fn an_oblique_tangent_plane_is_a_touch_not_a_tiny_circle() {
    for (n, len) in [
        ([3.0, 2.0, 6.0], 7.0),
        ([2.0, 3.0, 6.0], 7.0),
        ([1.0, 4.0, 8.0], 9.0),
    ] {
        for k in [1.0, 0.5, 3.0] {
            let normal = Vec3::new(n[0], n[1], n[2]);
            let sphere = Surface::Sphere(Sphere {
                frame: frame(Point3::ZERO, Vec3::Z),
                radius: k * len,
            });
            let plane = Surface::Plane(Plane {
                frame: frame_keeping_axis_length(Point3::ZERO + normal * k, normal),
            });
            assert_eq!(
                exact_surface_intersection(&sphere, &plane),
                Err(ExactIntersectionRefusal::NotRegularCurve),
                "normal {n:?} scale {k}"
            );
        }
    }
}

/// A plane exactly parallel to a cylinder axis, both given obliquely.
///
/// Axis (-3, -3, -3) and normal (-3, 1, 2) have dot product exactly zero.
/// In `f64` the normalised cosine was 2.8e-17, so the pair was treated as
/// an oblique cut and returned an ellipse with a major semi-axis of about
/// 3.6e16 instead of the two rulings.
#[test]
fn an_exactly_parallel_oblique_plane_gives_rulings_not_a_huge_ellipse() {
    for (axis, normal) in [
        (Vec3::new(-3.0, -3.0, -3.0), Vec3::new(-3.0, 1.0, 2.0)),
        (Vec3::new(1.0, 2.0, 2.0), Vec3::new(2.0, 1.0, -2.0)),
        (Vec3::new(3.0, 1.0, 1.0), Vec3::new(1.0, -2.0, -1.0)),
    ] {
        assert_eq!(axis.dot(normal), 0.0, "fixture must be exactly parallel");
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame_keeping_axis_length(Point3::ZERO, axis),
            radius: 1.0,
        });
        let plane = Surface::Plane(Plane {
            frame: frame_keeping_axis_length(Point3::ZERO, normal),
        });
        let curve = exact_surface_intersection(&cylinder, &plane)
            .unwrap_or_else(|e| panic!("axis {axis:?}: {e:?}"));
        assert_eq!(curve.derivation, Derivation::CylinderPlaneParallelRulings);
        assert_eq!(
            curve.branches.len(),
            2,
            "a plane through the axis cuts two rulings"
        );
    }
}

/// A plane exactly perpendicular to an oblique cylinder axis.
///
/// Normal = 2 * axis. In `f64` the normalised cosine came out as
/// 0.9999999999999996, so the section was classed as an oblique ellipse
/// rather than the circle it is.
#[test]
fn an_exactly_perpendicular_oblique_plane_gives_a_circle() {
    for axis in [
        Vec3::new(-3.0, -3.0, -1.0),
        Vec3::new(1.0, 2.0, 2.0),
        Vec3::new(2.0, 3.0, 6.0),
    ] {
        for k in [2.0, 3.0, -1.0] {
            let cylinder = Surface::Cylinder(Cylinder {
                frame: frame_keeping_axis_length(Point3::ZERO, axis),
                radius: 1.5,
            });
            let plane = Surface::Plane(Plane {
                frame: frame_keeping_axis_length(Point3::ZERO, axis * k),
            });
            let curve = exact_surface_intersection(&cylinder, &plane)
                .unwrap_or_else(|e| panic!("axis {axis:?} k {k}: {e:?}"));
            assert_eq!(
                curve.derivation,
                Derivation::CylinderPlanePerpendicularCircle,
                "axis {axis:?} k {k}"
            );
            match &curve.branches[0] {
                Curve3::Circle(c) => assert!((c.radius - 1.5).abs() < 1e-12),
                other => panic!("expected a circle, got {other:?}"),
            }
        }
    }
}

/// A plane tilted from perpendicular by less than `f64` can show in a
/// normalised cosine is still an ellipse, not a circle.
#[test]
fn a_barely_tilted_plane_is_still_an_ellipse() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame_keeping_axis_length(Point3::ZERO, Vec3::Z),
        radius: 1.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame_keeping_axis_length(Point3::ZERO, Vec3::new(1e-9, 0.0, 1.0)),
    });
    let curve = exact_surface_intersection(&cylinder, &plane).unwrap();
    assert_eq!(curve.derivation, Derivation::CylinderPlaneObliqueEllipse);
}

/// A plane tangent to a cylinder along an oblique normal shares one ruling.
#[test]
fn an_oblique_tangent_plane_shares_one_ruling() {
    // Axis z, normal (3, 4, 0) (length 5), plane through (3, 4, 0) * k:
    // distance 5k, radius 5k, exactly.
    // Scales whose products with 3, 4 and 5 are exact doubles.
    for k in [1.0, 0.25, 7.0] {
        let cylinder = Surface::Cylinder(Cylinder {
            frame: frame_keeping_axis_length(Point3::ZERO, Vec3::Z),
            radius: 5.0 * k,
        });
        let normal = Vec3::new(3.0, 4.0, 0.0);
        let plane = Surface::Plane(Plane {
            frame: frame_keeping_axis_length(Point3::ZERO + normal * k, normal),
        });
        let curve = exact_surface_intersection(&cylinder, &plane).unwrap();
        assert_eq!(curve.branches.len(), 1, "scale {k}");
    }
}

/// Tangency is judged on the doubles actually given. `0.2 * 3` rounds up
/// to 0.6000000000000001, so a plane through (0.2 * 3, 0.2 * 4, 0) with
/// radius 0.2 * 5 lies about 2e-16 outside the cylinder: those numbers
/// describe a plane that misses, and saying so is the exact answer.
#[test]
fn a_nearly_tangent_plane_is_judged_on_its_actual_doubles() {
    let k = 0.2;
    let normal = Vec3::new(3.0, 4.0, 0.0);
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame_keeping_axis_length(Point3::ZERO, Vec3::Z),
        radius: 5.0 * k,
    });
    let plane = Surface::Plane(Plane {
        frame: frame_keeping_axis_length(Point3::ZERO + normal * k, normal),
    });
    assert_eq!(
        exact_surface_intersection(&cylinder, &plane),
        Err(ExactIntersectionRefusal::Disjoint)
    );
}

/// A plane off parallel to the axis by far less than any float threshold
/// still cuts an ellipse (a very long one), not rulings. Its normal is
/// (1, 0, 2^-45): the dot product with the axis is exactly 2^-45.
#[test]
fn a_plane_barely_off_parallel_cuts_an_ellipse_not_rulings() {
    let tilt = 2f64.powi(-45);
    let cylinder = Surface::Cylinder(Cylinder {
        frame: frame_keeping_axis_length(Point3::ZERO, Vec3::Z),
        radius: 1.0,
    });
    let plane = Surface::Plane(Plane {
        frame: frame_keeping_axis_length(Point3::ZERO, Vec3::new(1.0, 0.0, tilt)),
    });
    let result = exact_surface_intersection(&cylinder, &plane).expect("an ellipse");
    assert_eq!(result.derivation, Derivation::CylinderPlaneObliqueEllipse);
    match result.single() {
        Curve3::Ellipse(e) => {
            let (minor, major) = (
                e.semi_axis_x.min(e.semi_axis_y),
                e.semi_axis_x.max(e.semi_axis_y),
            );
            assert!((minor - 1.0).abs() < 1e-12, "minor {minor}");
            assert!(major > 1e12, "major {major} must reflect the tiny tilt");
        }
        other => panic!("expected an ellipse, got {other:?}"),
    }
}

/// A plane perpendicular to a cone's axis cuts a circle, decided exactly
/// even when the axis is an unnormalised oblique vector like (1, 2, 2).
#[test]
fn a_plane_perpendicular_to_an_oblique_cone_axis_cuts_a_circle() {
    let axis = Vec3::new(1.0, 2.0, 2.0);
    let cone = Surface::Cone(Cone {
        frame: frame_keeping_axis_length(Point3::ZERO, axis),
        radius: 1.0,
        semi_angle: std::f64::consts::FRAC_PI_6,
    });
    let plane = Surface::Plane(Plane {
        frame: frame_keeping_axis_length(Point3::ZERO + axis, axis * 3.0),
    });
    let result = exact_surface_intersection(&cone, &plane).expect("a circle");
    assert!(
        matches!(result.single(), Curve3::Circle(_)),
        "got {:?}",
        result.single()
    );
}
