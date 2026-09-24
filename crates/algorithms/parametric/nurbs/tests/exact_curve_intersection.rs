//! Exact curve/surface and curve/curve intersection (#119).
//!
//! Hand cases have closed-form answers (stated inline). The random test
//! checks every reported point against the implicit equation of the other
//! operand and counts sign changes of that equation along the curve by
//! dense sampling, an oracle that shares no code with the solver.

use axiolid_core::{Frame2, Frame3, Point2, Point3, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Ellipse3, Line2, Line3};
use axiolid_evaluate::evaluate3;
use axiolid_nurbs::{
    exact_curve_curve_intersection2, exact_curve_curve_intersection3,
    exact_curve_surface_intersection, ExactCurveIntersection, ExactCurveParameter,
    ExactCurveRefusal,
};
use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface, Torus};

fn world() -> Frame3 {
    frame(Point3::ZERO)
}

fn frame(origin: Point3) -> Frame3 {
    Frame3 {
        origin,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn line(o: [f64; 3], d: [f64; 3]) -> Curve3 {
    Curve3::Line(Line3 {
        origin: Point3::new(o[0], o[1], o[2]),
        direction: Vec3::new(d[0], d[1], d[2]),
    })
}

fn circle(f: Frame3, r: f64) -> Curve3 {
    Curve3::Circle(Circle3 {
        frame: f,
        radius: r,
    })
}

fn points(result: ExactCurveIntersection) -> Vec<(f64, usize)> {
    match result {
        ExactCurveIntersection::Points(hits) => hits
            .into_iter()
            .map(|h| (h.parameter.approx(), h.multiplicity))
            .collect(),
        other => panic!("expected points, got {other:?}"),
    }
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-12 * (1.0 + a.abs().max(b.abs()))
}

fn assert_params(got: &[(f64, usize)], want: &[(f64, usize)]) {
    assert_eq!(got.len(), want.len(), "got {got:?}, want {want:?}");
    for (g, w) in got.iter().zip(want) {
        assert!(near(g.0, w.0) && g.1 == w.1, "got {got:?}, want {want:?}");
    }
}

#[test]
fn a_line_through_a_sphere_crosses_twice_and_a_touching_line_once_doubly() {
    let sphere = Surface::Sphere(Sphere {
        frame: world(),
        radius: 2.0,
    });
    // x from -3 along +x: enters at x = -2 (t = 1), leaves at x = 2 (t = 5).
    let through =
        exact_curve_surface_intersection(&line([-3.0, 0.0, 0.0], [1.0, 0.0, 0.0]), &sphere);
    assert_params(&points(through.unwrap()), &[(1.0, 1), (5.0, 1)]);
    // y = 2 grazes the top: one double point at t = 0.
    let touch = exact_curve_surface_intersection(&line([0.0, 2.0, 0.0], [1.0, 0.0, 0.0]), &sphere);
    assert_params(&points(touch.unwrap()), &[(0.0, 2)]);
    // y = 2 + 2^-40 misses: exactness, not a tolerance, decides it.
    let miss = exact_curve_surface_intersection(
        &line([0.0, 2.0 + 2f64.powi(-40), 0.0], [1.0, 0.0, 0.0]),
        &sphere,
    );
    assert_params(&points(miss.unwrap()), &[]);
}

#[test]
fn irrational_crossings_are_isolated_not_rounded() {
    // x^2 + y^2 = 2 on y = 0: x = +- sqrt(2), never a dyadic.
    let cylinder = Surface::Cylinder(Cylinder {
        frame: world(),
        radius: 2f64.sqrt(),
    });
    // sqrt(2) as f64 is not sqrt(2); the cylinder radius is that double r,
    // so the roots are +- r exactly, i.e. the dyadic r itself.
    let r = 2f64.sqrt();
    let hits = points(
        exact_curve_surface_intersection(&line([0.0, 0.0, 5.0], [1.0, 0.0, 0.0]), &cylinder)
            .unwrap(),
    );
    assert_params(&hits, &[(-r, 1), (r, 1)]);
    // Radius 1, line y = 1/2: x = +- sqrt(3)/2, a real irrational root.
    let unit = Surface::Cylinder(Cylinder {
        frame: world(),
        radius: 1.0,
    });
    let hits = points(
        exact_curve_surface_intersection(&line([0.0, 0.5, 0.0], [1.0, 0.0, 0.0]), &unit).unwrap(),
    );
    let s = 3f64.sqrt() / 2.0;
    assert_params(&hits, &[(-s, 1), (s, 1)]);
}

#[test]
fn a_line_on_a_cylinder_wall_is_contained_and_a_parallel_one_misses() {
    let cylinder = Surface::Cylinder(Cylinder {
        frame: world(),
        radius: 1.0,
    });
    let on_wall =
        exact_curve_surface_intersection(&line([1.0, 0.0, 0.0], [0.0, 0.0, 3.0]), &cylinder);
    assert_eq!(on_wall.unwrap(), ExactCurveIntersection::Contained);
    let inside =
        exact_curve_surface_intersection(&line([0.5, 0.0, 0.0], [0.0, 0.0, 1.0]), &cylinder);
    assert_params(&points(inside.unwrap()), &[]);
}

#[test]
fn a_circle_meets_a_plane_including_at_its_antipode() {
    // Circle radius 1 in the xz-plane (x = X, y = Z): theta measured from
    // +x toward +z. Plane x = -1 touches it at theta = pi only, the point
    // where the half-angle parameter is infinite.
    let f = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Z,
        z: -Vec3::Y,
    };
    let plane_at = |x: f64| {
        Surface::Plane(Plane {
            frame: Frame3 {
                origin: Point3::new(x, 0.0, 0.0),
                x: Vec3::Y,
                y: Vec3::Z,
                z: Vec3::X,
            },
        })
    };
    let hits = exact_curve_surface_intersection(&circle(f, 1.0), &plane_at(-1.0)).unwrap();
    match &hits {
        ExactCurveIntersection::Points(h) => {
            assert_eq!(h.len(), 1);
            assert_eq!(h[0].parameter, ExactCurveParameter::Antipode);
            assert_eq!(h[0].multiplicity, 2, "touching, not crossing");
        }
        other => panic!("{other:?}"),
    }
    // x = 0 crosses at theta = pi/2 and 3 pi/2, in that order.
    let hits = points(exact_curve_surface_intersection(&circle(f, 1.0), &plane_at(0.0)).unwrap());
    let half = std::f64::consts::FRAC_PI_2;
    assert_params(&hits, &[(half, 1), (3.0 * half, 1)]);
    // x = 1/2 crosses at theta = pi/3 and 5 pi/3; x = -1/2 at 2pi/3, 4pi/3.
    let third = std::f64::consts::PI / 3.0;
    let hits = points(exact_curve_surface_intersection(&circle(f, 1.0), &plane_at(0.5)).unwrap());
    assert_params(&hits, &[(third, 1), (5.0 * third, 1)]);
    let hits = points(exact_curve_surface_intersection(&circle(f, 1.0), &plane_at(-0.5)).unwrap());
    assert_params(&hits, &[(2.0 * third, 1), (4.0 * third, 1)]);
    // A circle lying in the plane is contained.
    let flat = Surface::Plane(Plane { frame: world() });
    assert_eq!(
        exact_curve_surface_intersection(&circle(world(), 1.0), &flat).unwrap(),
        ExactCurveIntersection::Contained
    );
}

#[test]
fn a_plane_crossing_a_circle_at_its_antipode_counts_it_once() {
    // Unit circle in the xy-plane, theta from +x. The plane y = 0 crosses
    // it at theta = 0 and at theta = pi, both transversally: the antipode
    // (where the half-angle parameter is infinite) with multiplicity ONE.
    let plane = Surface::Plane(Plane {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::Z,
            y: Vec3::X,
            z: Vec3::Y,
        },
    });
    let hits = exact_curve_surface_intersection(&circle(world(), 1.0), &plane).unwrap();
    let ExactCurveIntersection::Points(h) = &hits else {
        panic!("{hits:?}");
    };
    assert_eq!(h.len(), 2, "{h:?}");
    assert!(near(h[0].parameter.approx(), 0.0), "{h:?}");
    assert_eq!(h[0].multiplicity, 1);
    assert_eq!(h[1].parameter, ExactCurveParameter::Antipode);
    assert_eq!(h[1].multiplicity, 1, "crossing, not touching");
}

#[test]
fn a_cone_counts_only_its_own_nappe() {
    // Radius 1 at z = 0, semi-angle -45 deg: radius 1 - z, apex at z = 1,
    // the modelled nappe is z <= 1. The vertical line x = 0.5 meets the
    // nappe at z = 0.5 and the mirror nappe (radius z - 1) at z = 1.5.
    let cone = Surface::Cone(Cone {
        frame: world(),
        radius: 1.0,
        semi_angle: -std::f64::consts::FRAC_PI_4,
    });
    let slope = cone_slope(&cone);
    let hits = points(
        exact_curve_surface_intersection(&line([0.5, 0.0, 0.0], [0.0, 0.0, 1.0]), &cone).unwrap(),
    );
    assert_eq!(hits.len(), 1, "the mirror nappe must not count: {hits:?}");
    // radius + z * slope = 0.5  =>  z = -0.5 / slope.
    assert!(near(hits[0].0, -0.5 / slope), "{hits:?}");
    // A generator through the apex lies on both nappes: only the part with
    // z <= 1 is on the surface, a ray, which is refused by name.
    // Built from the double slope itself (x = 1 + z * slope), so it lies
    // on the cone exactly; tan(-pi/4) as a double is not -1.
    let generator = line([1.0, 0.0, 0.0], [slope, 0.0, 1.0]);
    assert_eq!(
        exact_curve_surface_intersection(&generator, &cone),
        Err(ExactCurveRefusal::PartialOverlap)
    );
    // The same generator run backwards: the on-nappe condition now has a
    // positive leading coefficient, so only its sign change (the apex)
    // shows that part of the line is off the surface.
    let reversed = line([1.0, 0.0, 0.0], [-slope, 0.0, -1.0]);
    assert_eq!(
        exact_curve_surface_intersection(&reversed, &cone),
        Err(ExactCurveRefusal::PartialOverlap)
    );
}

fn cone_slope(cone: &Surface) -> f64 {
    match cone {
        Surface::Cone(c) => c.semi_angle.tan(),
        _ => unreachable!(),
    }
}

#[test]
fn a_line_along_a_torus_axis_misses_and_one_across_it_meets_four_times() {
    let torus = Surface::Torus(Torus {
        frame: world(),
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let axis = exact_curve_surface_intersection(&line([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]), &torus);
    assert_params(&points(axis.unwrap()), &[]);
    // Along x at z = 0: |x| in {2, 4}.
    let across = exact_curve_surface_intersection(&line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), &torus);
    assert_params(
        &points(across.unwrap()),
        &[(-4.0, 1), (-2.0, 1), (2.0, 1), (4.0, 1)],
    );
    // Along x at z = 1 (the tube's top): touches at |x| = 3, doubly.
    let top = exact_curve_surface_intersection(&line([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]), &torus);
    assert_params(&points(top.unwrap()), &[(-3.0, 2), (3.0, 2)]);
}

#[test]
fn an_ellipse_around_an_elliptical_cylinder() {
    let cyl = Surface::EllipticalCylinder(EllipticalCylinder {
        frame: world(),
        semi_axis_x: 2.0,
        semi_axis_y: 1.0,
    });
    // The same ellipse, lifted to z = 3: contained.
    let same = Curve3::Ellipse(Ellipse3 {
        frame: frame(Point3::new(0.0, 0.0, 3.0)),
        semi_axis_x: 2.0,
        semi_axis_y: 1.0,
    });
    assert_eq!(
        exact_curve_surface_intersection(&same, &cyl).unwrap(),
        ExactCurveIntersection::Contained
    );
    // A circle of radius 1 touches it at (0, +-1): theta = pi/2, 3 pi/2.
    let small = circle(world(), 1.0);
    let half = std::f64::consts::FRAC_PI_2;
    assert_params(
        &points(exact_curve_surface_intersection(&small, &cyl).unwrap()),
        &[(half, 2), (3.0 * half, 2)],
    );
}

#[test]
fn plane_curves_meet_exactly() {
    let c = |x: f64, y: f64, r: f64| {
        Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(x, y),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: r,
        })
    };
    // Two unit circles 2 apart touch at theta = 0 on the first.
    let hits =
        points(exact_curve_curve_intersection2(&c(0.0, 0.0, 1.0), &c(2.0, 0.0, 1.0)).unwrap());
    assert_params(&hits, &[(0.0, 2)]);
    // Distance 1: crossings at theta = +- pi/3.
    let third = std::f64::consts::PI / 3.0;
    let hits =
        points(exact_curve_curve_intersection2(&c(0.0, 0.0, 1.0), &c(1.0, 0.0, 1.0)).unwrap());
    assert_params(&hits, &[(third, 1), (5.0 * third, 1)]);
    // Concentric: none. Identical: contained.
    let hits =
        points(exact_curve_curve_intersection2(&c(0.0, 0.0, 1.0), &c(0.0, 0.0, 2.0)).unwrap());
    assert_params(&hits, &[]);
    assert_eq!(
        exact_curve_curve_intersection2(&c(0.0, 0.0, 1.0), &c(0.0, 0.0, 1.0)).unwrap(),
        ExactCurveIntersection::Contained
    );
    // Line y = 1 touches the unit circle at theta = pi/2.
    let l = Curve2::Line(Line2 {
        origin: Point2::new(0.0, 1.0),
        direction: Vec2::X,
    });
    let hits = points(exact_curve_curve_intersection2(&c(0.0, 0.0, 1.0), &l).unwrap());
    assert_params(&hits, &[(std::f64::consts::FRAC_PI_2, 2)]);
}

#[test]
fn skew_space_curves_miss_and_crossing_ones_meet_once() {
    let a = line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    let skew = line([0.0, 1.0, 1.0], [0.0, 0.0, 1.0]);
    assert_params(
        &points(exact_curve_curve_intersection3(&a, &skew).unwrap()),
        &[],
    );
    // Skew again, but with the first line along y: now the first of the
    // three line equations has a root (t = 0) that the others do not
    // share. Only their COMMON root is an intersection.
    let along_y = line([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let above = line([1.0, 0.0, 1.0], [0.0, 0.0, 1.0]);
    assert_params(
        &points(exact_curve_curve_intersection3(&along_y, &above).unwrap()),
        &[],
    );
    let crossing = line([3.0, -1.0, 0.0], [0.0, 1.0, 0.0]);
    assert_params(
        &points(exact_curve_curve_intersection3(&a, &crossing).unwrap()),
        &[(3.0, 1)],
    );
    // A tilted circle through the x axis at x = +-1.
    let tilted = circle(
        Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::new(0.0, 1.0, 1.0),
            z: Vec3::new(0.0, -1.0, 1.0),
        },
        1.0,
    );
    assert_params(
        &points(exact_curve_curve_intersection3(&a, &tilted).unwrap()),
        &[(-1.0, 1), (1.0, 1)],
    );
}

#[test]
fn malformed_and_unsupported_inputs_are_refused_by_name() {
    let sphere = Surface::Sphere(Sphere {
        frame: world(),
        radius: 1.0,
    });
    assert_eq!(
        exact_curve_surface_intersection(&line([f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0]), &sphere),
        Err(ExactCurveRefusal::NonFinite)
    );
    assert_eq!(
        exact_curve_surface_intersection(&line([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]), &sphere),
        Err(ExactCurveRefusal::Degenerate)
    );
    assert_eq!(
        exact_curve_surface_intersection(&circle(world(), 0.0), &sphere),
        Err(ExactCurveRefusal::Degenerate)
    );
}

// --- random oracle -------------------------------------------------------

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next()
    }
    /// Snapped to 1/8 so that tangencies and containment happen often.
    fn grid(&mut self, lo: f64, hi: f64) -> f64 {
        (self.range(lo, hi) * 8.0).round() / 8.0
    }
    fn unit(&mut self) -> Vec3 {
        loop {
            let v = Vec3::new(
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
            );
            if v.length() > 0.2 {
                return v.normalize();
            }
        }
    }
}

fn random_frame(rng: &mut Rng, axis_aligned: bool) -> Frame3 {
    let origin = Point3::new(
        rng.grid(-2.0, 2.0),
        rng.grid(-2.0, 2.0),
        rng.grid(-2.0, 2.0),
    );
    if axis_aligned {
        return frame(origin);
    }
    let z = rng.unit();
    let x = z.any_orthonormal_vector();
    Frame3 {
        origin,
        x,
        y: z.cross(x),
        z,
    }
}

fn radius_of(surface: &Surface) -> Option<f64> {
    match surface {
        Surface::Cylinder(s) => Some(s.radius),
        Surface::Sphere(s) => Some(s.radius),
        _ => None,
    }
}

fn surface_origin(surface: &Surface) -> Point3 {
    match surface {
        Surface::Plane(s) => s.frame.origin,
        Surface::Cylinder(s) => s.frame.origin,
        Surface::Sphere(s) => s.frame.origin,
        Surface::Torus(s) => s.frame.origin,
        _ => unreachable!(),
    }
}

/// Implicit function of a surface, in world coordinates (orthonormal
/// frames only, as the oracle generates).
fn implicit(surface: &Surface, p: Point3) -> f64 {
    let local = |f: &Frame3| {
        let d = p - f.origin;
        (d.dot(f.x), d.dot(f.y), d.dot(f.z))
    };
    match surface {
        Surface::Plane(s) => local(&s.frame).2,
        Surface::Cylinder(s) => {
            let (x, y, _) = local(&s.frame);
            x * x + y * y - s.radius * s.radius
        }
        Surface::Sphere(s) => {
            let (x, y, z) = local(&s.frame);
            x * x + y * y + z * z - s.radius * s.radius
        }
        Surface::Torus(s) => {
            let (x, y, z) = local(&s.frame);
            let q = (x * x + y * y).sqrt() - s.major_radius;
            q * q + z * z - s.minor_radius * s.minor_radius
        }
        _ => unreachable!(),
    }
}

#[test]
fn random_curves_and_surfaces_match_a_sampling_oracle() {
    let mut rng = Rng(0x0119_5eed);
    let (mut scenes, mut hits_total, mut tangents) = (0, 0, 0);
    for scene in 0..600 {
        let aligned = scene % 2 == 0;
        let surface = match scene % 4 {
            0 => Surface::Plane(Plane {
                frame: random_frame(&mut rng, aligned),
            }),
            1 => Surface::Cylinder(Cylinder {
                frame: random_frame(&mut rng, aligned),
                radius: rng.grid(0.5, 2.0).max(0.25),
            }),
            2 => Surface::Sphere(Sphere {
                frame: random_frame(&mut rng, aligned),
                radius: rng.grid(0.5, 2.0).max(0.25),
            }),
            _ => Surface::Torus(Torus {
                frame: random_frame(&mut rng, aligned),
                major_radius: rng.grid(1.5, 2.5),
                minor_radius: rng.grid(0.25, 1.0).max(0.25),
            }),
        };
        // Anchor the curve near the surface so most scenes meet it; on
        // axis-aligned scenes use axis directions and 1/8 offsets, so
        // grazing contact (a tangency) happens often.
        let base = surface_origin(&surface);
        let offset = |rng: &mut Rng| rng.grid(-2.5, 2.5);
        let origin = Point3::new(
            base.x + offset(&mut rng),
            base.y + offset(&mut rng),
            base.z + offset(&mut rng),
        );
        let graze = radius_of(&surface).filter(|_| aligned && scene % 5 == 0);
        let curve = if let Some(r) = graze {
            // A line exactly one radius from an axis-aligned cylinder's axis
            // or a sphere's centre, along x: it touches, doubly.
            line(
                [base.x + offset(&mut rng), base.y + r, base.z],
                [1.0, 0.0, 0.0],
            )
        } else if scene % 3 == 0 {
            let d = if aligned {
                [Vec3::X, Vec3::Y, Vec3::Z][scene % 7 % 3]
            } else {
                rng.unit()
            };
            line([origin.x, origin.y, origin.z], [d.x, d.y, d.z])
        } else {
            let mut f = random_frame(&mut rng, aligned);
            f.origin = origin;
            circle(f, rng.grid(0.5, 2.5).max(0.25))
        };
        let result = exact_curve_surface_intersection(&curve, &surface).unwrap();
        let hits = if let ExactCurveIntersection::Points(hits) = result {
            hits
        } else {
            // Contained: the implicit function must vanish along the curve.
            for k in 0..64 {
                let t = f64::from(k) / 64.0 * std::f64::consts::TAU;
                let p = evaluate3(&curve, t).unwrap();
                assert!(
                    implicit(&surface, p).abs() < 1e-9,
                    "scene {scene}: not contained"
                );
            }
            scenes += 1;
            continue;
        };
        // Every reported point lies on the surface.
        for h in &hits {
            assert!(
                implicit(&surface, h.point).abs() < 1e-8,
                "scene {scene}: hit off the surface: {h:?}"
            );
        }
        // Transverse hits equal the sign changes of the implicit function
        // sampled along the curve, away from tangencies.
        let is_line = matches!(curve, Curve3::Line(_));
        let (lo, hi, n) = if is_line {
            // Past every reported hit, so the window cannot hide one: a line
            // almost parallel to a cylinder can meet it far away. Beyond the
            // outermost root the sign no longer changes.
            let reach = hits
                .iter()
                .map(|h| h.parameter.approx().abs())
                .fold(40.0, f64::max)
                + 1.0;
            (-reach, reach, 40_000)
        } else {
            (0.0, std::f64::consts::TAU, 20_000)
        };
        let sample = |k: usize| {
            let t = lo + (hi - lo) * k as f64 / n as f64;
            implicit(&surface, evaluate3(&curve, t).unwrap())
        };
        let mut changes = 0;
        let mut last = sample(0);
        for k in 1..=n {
            let v = sample(k);
            if (v > 0.0) != (last > 0.0) && v != 0.0 && last != 0.0 {
                changes += 1;
            }
            if v != 0.0 {
                last = v;
            }
        }
        if !is_line && (sample(0) > 0.0) != (last > 0.0) {
            // A closed curve: the wrap-around is not a change.
        }
        let odd = hits.iter().filter(|h| h.multiplicity % 2 == 1).count();
        assert_eq!(
            odd, changes,
            "scene {scene}: {odd} odd-multiplicity hits vs {changes} sampled sign changes; {hits:?}"
        );
        // Parameters strictly increase along the curve.
        for w in hits.windows(2) {
            assert!(
                w[0].parameter.approx() < w[1].parameter.approx(),
                "scene {scene}: out of order {hits:?}"
            );
        }
        hits_total += hits.len();
        tangents += hits.iter().filter(|h| h.is_tangent()).count();
        scenes += 1;
    }
    eprintln!("scenes={scenes} hits={hits_total} tangent={tangents}");
    // Measured 434 hits, 30 of them tangent (mostly from the grazing
    // family). The floors catch a generator that stops meeting or grazing.
    assert!(hits_total > 400, "too few hits: {hits_total}");
    assert!(tangents >= 20, "too few tangencies exercised: {tangents}");
}
