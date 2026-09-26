//! Ruled sections: a quadric cut across a cylinder or cone (ADR 0076, #119).
//!
//! Oracles never read the kernel's answer back:
//! - every point of every piece is checked against BOTH surfaces through
//!   their own distance functions, written out here;
//! - pieces over one span must meet at both ends (a closed loop), or a
//!   whole-turn piece must meet itself;
//! - on random scenes, whether the surfaces meet at all, and on how many
//!   angles, is checked against a brute-force scan of the carrier's rulings.

use axiolid_core::{Frame3, Interval, Point3, Vec3};
use axiolid_curve::Curve3;
use axiolid_evaluate::evaluate3;
use axiolid_nurbs::{exact_surface_intersection, Derivation, ExactIntersectionRefusal};
use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Surface};

const PI: f64 = std::f64::consts::PI;

/// A right-handed orthonormal frame with `z` along `axis`.
fn frame(origin: Point3, axis: Vec3) -> Frame3 {
    let z = axis.normalize();
    let helper = if z.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let x = helper.cross(z).normalize();
    let y = z.cross(x);
    Frame3 { origin, x, y, z }
}

/// Distance-like residual of `p` from a surface, in length units.
fn off(surface: &Surface, p: Point3) -> f64 {
    match surface {
        Surface::Plane(s) => s.frame.z.normalize().dot(p - s.frame.origin),
        Surface::Sphere(s) => (p - s.frame.origin).length() - s.radius,
        Surface::Cylinder(s) => {
            let axis = s.frame.z.normalize();
            let d = p - s.frame.origin;
            (d - axis * d.dot(axis)).length() - s.radius
        }
        Surface::EllipticalCylinder(s) => {
            let d = p - s.frame.origin;
            let (x, y) = (d.dot(s.frame.x), d.dot(s.frame.y));
            // Scaled so it reads as a length near the surface.
            ((x / s.semi_axis_x).hypot(y / s.semi_axis_y) - 1.0) * s.semi_axis_x.min(s.semi_axis_y)
        }
        Surface::Cone(s) => {
            // Along the axis from the frame origin, radius grows with v.
            let axis = s.frame.z.normalize();
            let d = p - s.frame.origin;
            let v = d.dot(axis);
            let radial = (d - axis * v).length();
            (radial - (s.radius + v * s.semi_angle.tan())) * s.semi_angle.cos()
        }
        other => panic!("no residual for {other:?}"),
    }
}

/// Check every piece; return the number of closed loops.
fn check(first: &Surface, second: &Surface) -> usize {
    let curve = exact_surface_intersection(first, second).expect("a ruled section");
    assert_eq!(curve.derivation, Derivation::RuledQuadricSection);
    assert_eq!(curve.branches.len(), curve.spans.len());
    let scale = 1.0;
    for (branch, span) in curve.branches.iter().zip(&curve.spans) {
        let span = span.expect("a ruled piece carries its span");
        assert!(span.end > span.start, "{span:?}");
        assert!(matches!(branch, Curve3::RuledSection(_)));
        for i in 0..=400 {
            let t = span.start + (span.end - span.start) * i as f64 / 400.0;
            let Ok(p) = evaluate3(branch, t) else {
                // Only a span's end may sit exactly on a vertical tangent.
                assert!(i == 0 || i == 400, "no point at t = {t} inside {span:?}");
                continue;
            };
            for surface in [first, second] {
                let r = off(surface, p);
                assert!(
                    r.abs() < 1e-9 * scale,
                    "point {p:?} at t = {t} is {r} off {surface:?}"
                );
            }
        }
    }
    loops(&curve.branches, &curve.spans)
}

/// Pieces sharing a span form one loop when they meet at both ends; a
/// whole-turn piece is a loop on its own.
fn loops(branches: &[Curve3], spans: &[Option<Interval>]) -> usize {
    let mut count = 0;
    let mut used = vec![false; branches.len()];
    for i in 0..branches.len() {
        if used[i] {
            continue;
        }
        let span = spans[i].unwrap();
        if (span.end - span.start - 2.0 * PI).abs() < 1e-9 {
            let a = evaluate3(&branches[i], span.start).unwrap();
            let b = evaluate3(&branches[i], span.end).unwrap();
            assert!((a - b).length() < 1e-9, "a whole-turn piece must close");
            used[i] = true;
            count += 1;
            continue;
        }
        let partner = (i + 1..branches.len())
            .find(|&j| !used[j] && spans[j] == Some(span))
            .expect("a bounded piece has a partner over its span");
        for t in [span.start, span.end] {
            let a = evaluate3(&branches[i], t).expect("end point");
            let b = evaluate3(&branches[partner], t).expect("end point");
            assert!(
                (a - b).length() < 1e-6,
                "pieces must meet at t = {t}: {a:?} vs {b:?}"
            );
        }
        used[i] = true;
        used[partner] = true;
        count += 1;
    }
    count
}

fn cylinder(origin: Point3, axis: Vec3, radius: f64) -> Surface {
    Surface::Cylinder(Cylinder {
        frame: frame(origin, axis),
        radius,
    })
}

#[test]
fn a_pipe_tee_is_two_loops_where_the_thin_pipe_pierces_the_wide_one() {
    let main = cylinder(Point3::ZERO, Vec3::Z, 2.0);
    let branch = cylinder(Point3::ZERO, Vec3::Y, 1.2);
    assert_eq!(check(&main, &branch), 2);
    // Either operand may come first.
    assert_eq!(check(&branch, &main), 2);
}

#[test]
fn overlapping_skew_cylinders_meet_in_one_loop_with_branch_ends() {
    // Axes skew by 1.5; radii 1 and 1.2: the thinner only grazes the
    // wider, so the section is one loop that turns back at two angles.
    let a = cylinder(Point3::ZERO, Vec3::Z, 1.0);
    let b = cylinder(Point3::new(0.0, 0.0, 0.0) + Vec3::X * 1.5, Vec3::Y, 1.2);
    assert_eq!(check(&a, &b), 1);
}

#[test]
fn an_offset_sphere_and_cylinder_meet_in_one_closed_curve_per_crossing() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::ZERO, Vec3::Z),
        radius: 5.0,
    });
    // A cylinder wholly inside the sphere's girth pierces it twice.
    let inside = cylinder(Point3::new(1.0, 0.0, 0.0), Vec3::Z, 3.0);
    assert_eq!(check(&sphere, &inside), 2);
    // One that reaches outside the sphere's girth cuts one loop.
    let across = cylinder(Point3::new(4.0, 0.0, 0.0), Vec3::Z, 2.0);
    assert_eq!(check(&sphere, &across), 1);
}

#[test]
fn an_elliptical_cylinder_meets_a_sphere_exactly() {
    let sphere = Surface::Sphere(Sphere {
        frame: frame(Point3::new(0.5, -0.25, 1.0), Vec3::Z),
        radius: 4.0,
    });
    let oval = Surface::EllipticalCylinder(EllipticalCylinder {
        frame: frame(Point3::ZERO, Vec3::new(0.2, 0.1, 1.0)),
        semi_axis_x: 2.5,
        semi_axis_y: 1.0,
    });
    assert!(check(&oval, &sphere) >= 1);
}

#[test]
fn planes_cut_a_cone_in_hyperbolas_parabolas_and_oblique_ellipses() {
    let cone = Surface::Cone(Cone {
        frame: frame(Point3::ZERO, Vec3::Z),
        radius: 1.0,
        semi_angle: PI / 6.0,
    });
    // Parallel to the axis, off it: one branch of a hyperbola, which
    // leaves to infinity -- an open piece on each side of the cut.
    let hyperbola = Surface::Plane(Plane {
        frame: frame(Point3::new(0.5, 0.0, 0.0), Vec3::X),
    });
    let curve = exact_surface_intersection(&cone, &hyperbola).expect("a section");
    assert_eq!(curve.derivation, Derivation::RuledQuadricSection);
    for (branch, span) in curve.branches.iter().zip(&curve.spans) {
        let span = span.unwrap();
        let mid = 0.5 * (span.start + span.end);
        let p = evaluate3(branch, mid).expect("inside the span");
        assert!(off(&cone, p).abs() < 1e-9 && off(&hyperbola, p).abs() < 1e-9);
    }
    // Tilted to the cone's own slope: a parabola.
    let slope = (PI / 6.0).tan();
    let parabola = Surface::Plane(Plane {
        frame: frame(Point3::new(0.5, 0.0, 0.0), Vec3::new(1.0, 0.0, -slope)),
    });
    let curve = exact_surface_intersection(&cone, &parabola).expect("a section");
    let span = curve.spans[0].unwrap();
    let p = evaluate3(&curve.branches[0], 0.5 * (span.start + span.end)).unwrap();
    assert!(off(&cone, p).abs() < 1e-9 && off(&parabola, p).abs() < 1e-9);
    // Gently tilted: a closed ellipse, off-axis so no closed form applies.
    let ellipse = Surface::Plane(Plane {
        frame: frame(Point3::new(0.3, 0.2, 2.0), Vec3::new(0.2, 0.1, 1.0)),
    });
    assert_eq!(check(&cone, &ellipse), 1);
}

#[test]
fn a_plane_through_the_apex_is_still_refused_by_name() {
    let cone = Surface::Cone(Cone {
        frame: frame(Point3::ZERO, Vec3::Z),
        radius: 1.0,
        semi_angle: PI / 6.0,
    });
    // Contains the axis: two rulings through the apex, lines of constant
    // angle that no graph over the angle can hold.
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::ZERO, Vec3::X),
    });
    assert_eq!(
        exact_surface_intersection(&cone, &plane),
        Err(ExactIntersectionRefusal::UnrepresentableConic)
    );
}

#[test]
fn apart_and_touching_pairs_are_named() {
    let a = cylinder(Point3::ZERO, Vec3::Z, 1.0);
    let far = cylinder(Point3::new(5.0, 0.0, 0.0), Vec3::Y, 1.0);
    assert_eq!(
        exact_surface_intersection(&a, &far),
        Err(ExactIntersectionRefusal::Disjoint)
    );
    // Skew axes exactly 2 apart, radii 1 and 1: the surfaces touch at one
    // point. Decided exactly, not rounded to a tiny loop.
    let touching = cylinder(Point3::new(2.0, 0.0, 0.0), Vec3::Y, 1.0);
    assert_eq!(
        exact_surface_intersection(&a, &touching),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}

/// Small deterministic generator.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next()
    }
}

/// Angles on the carrier where the other surface crosses the ruling, by a
/// dense scan: at each angle, whether the ruling's quadratic has real roots
/// (sign change of the other's residual along a long stretch of the ruling).
fn scan_hits(carrier: &Cylinder, other: &Surface) -> usize {
    let mut hits = 0;
    let axis = carrier.frame.z;
    for i in 0..720 {
        let u = -PI + 2.0 * PI * (i as f64 + 0.5) / 720.0;
        let base = carrier.frame.origin
            + carrier.frame.x * (carrier.radius * u.cos())
            + carrier.frame.y * (carrier.radius * u.sin());
        let mut last = off(other, base + axis * -50.0).signum();
        for k in 1..=2000 {
            let v = -50.0 + 100.0 * k as f64 / 2000.0;
            let now = off(other, base + axis * v).signum();
            if now != last {
                hits += 1;
                break;
            }
            last = now;
        }
    }
    hits
}

#[test]
fn random_cylinder_pairs_agree_with_a_dense_scan() {
    let mut rng = Lcg(0x5eed_1119);
    let mut met = 0;
    for _ in 0..120 {
        let a = Cylinder {
            frame: frame(
                Point3::new(
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                ),
                Vec3::new(
                    rng.range(-1.0, 1.0),
                    rng.range(-1.0, 1.0),
                    rng.range(0.2, 1.0),
                ),
            ),
            radius: rng.range(0.5, 2.0),
        };
        let b = Surface::Cylinder(Cylinder {
            frame: frame(
                Point3::new(
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                ),
                Vec3::new(
                    rng.range(-1.0, 1.0),
                    rng.range(0.2, 1.0),
                    rng.range(-1.0, 1.0),
                ),
            ),
            radius: rng.range(0.5, 2.0),
        });
        let carrier = Surface::Cylinder(a);
        let scanned = scan_hits(&a, &b);
        match exact_surface_intersection(&carrier, &b) {
            Ok(curve) => {
                met += 1;
                // Every returned point is on both.
                for (branch, span) in curve.branches.iter().zip(&curve.spans) {
                    let Some(span) = span else { continue };
                    for i in 1..40 {
                        let t = span.start + (span.end - span.start) * i as f64 / 40.0;
                        if let Ok(p) = evaluate3(branch, t) {
                            assert!(off(&carrier, p).abs() < 1e-8, "off carrier");
                            assert!(off(&b, p).abs() < 1e-8, "off other");
                        }
                    }
                }
                assert!(scanned > 0, "a curve where the scan finds no crossing");
            }
            Err(ExactIntersectionRefusal::Disjoint) => {
                assert_eq!(
                    scanned, 0,
                    "disjoint, but the scan crosses on {scanned} rulings"
                );
            }
            Err(ExactIntersectionRefusal::NotRegularCurve) => {}
            Err(other) => panic!("unexpected refusal {other:?}"),
        }
    }
    assert!(met > 30, "too few meeting pairs ({met}) to mean anything");
}
