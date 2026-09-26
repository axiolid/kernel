//! Sections with no closed form, traced as implicit curves (#119,
//! ADR 0077): a torus against a cylinder, cone or torus off its axis.
//!
//! Oracles written independently of the tracer:
//! - every point of every curve lies on both surfaces, by distance
//!   functions of the surfaces' own definitions;
//! - each curve is continuous and a closed loop closes in space;
//! - completeness: on every tube circle of the torus, the number of times
//!   the traced curves cross it equals the number of sign changes of the
//!   other surface's distance around that circle, scanned densely.

use axiolid_core::{Frame3, Point2, Point3, Vec3};
use axiolid_curve::Curve3;
use axiolid_evaluate::evaluate3;
use axiolid_nurbs::{
    exact_surface_intersection, implicit_surface_intersection, Derivation, ExactIntersectionRefusal,
};
use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Surface, Torus};

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

fn frame(origin: Point3, axis: Vec3) -> Frame3 {
    let z = axis.normalize();
    let helper = if z.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let x = helper.cross(z).normalize();
    let y = z.cross(x);
    Frame3 { origin, x, y, z }
}

fn torus(origin: Point3, axis: Vec3, major: f64, minor: f64) -> Surface {
    Surface::Torus(Torus {
        frame: frame(origin, axis),
        major_radius: major,
        minor_radius: minor,
    })
}

fn cylinder(origin: Point3, axis: Vec3, radius: f64) -> Surface {
    Surface::Cylinder(Cylinder {
        frame: frame(origin, axis),
        radius,
    })
}

/// Signed distance-like residual from a surface.
fn off(surface: &Surface, p: Point3) -> f64 {
    match surface {
        Surface::Cylinder(s) => {
            let axis = s.frame.z.normalize();
            let d = p - s.frame.origin;
            (d - axis * d.dot(axis)).length() - s.radius
        }
        Surface::EllipticalCylinder(s) => {
            let d = p - s.frame.origin;
            let (x, y) = (d.dot(s.frame.x), d.dot(s.frame.y));
            ((x / s.semi_axis_x).hypot(y / s.semi_axis_y) - 1.0) * s.semi_axis_x.min(s.semi_axis_y)
        }
        Surface::Cone(s) => {
            let axis = s.frame.z.normalize();
            let d = p - s.frame.origin;
            let v = d.dot(axis);
            let radial = (d - axis * v).length();
            (radial - (s.radius + v * s.semi_angle.tan())) * s.semi_angle.cos()
        }
        Surface::Torus(t) => {
            let axis = t.frame.z.normalize();
            let d = p - t.frame.origin;
            let h = d.dot(axis);
            let radial = (d - axis * h).length();
            (radial - t.major_radius).hypot(h) - t.minor_radius
        }
        other => panic!("no residual for {other:?}"),
    }
}

/// Check every curve; return how many there are.
fn check(first: &Surface, second: &Surface) -> usize {
    let curve = exact_surface_intersection(first, second).expect("a traced section");
    assert_eq!(curve.derivation, Derivation::ImplicitTrace);
    for (branch, span) in curve.branches.iter().zip(&curve.spans) {
        let span = span.expect("a traced curve carries its span");
        let n = 2000;
        let mut previous: Option<Point3> = None;
        for i in 0..=n {
            let t = span.start + (span.end - span.start) * i as f64 / n as f64;
            let p = evaluate3(branch, t).unwrap_or_else(|e| panic!("no point at {t}: {e}"));
            for surface in [first, second] {
                let r = off(surface, p);
                assert!(r.abs() < 1e-9, "point at {t} is {r} off {surface:?}");
            }
            if let Some(q) = previous {
                assert!((p - q).length() < 0.2, "jump at {t}: {q:?} -> {p:?}");
            }
            previous = Some(p);
        }
    }
    curve.branches.len()
}

/// The torus point at `(u, v)`.
fn torus_point(t: &Torus, u: f64, v: f64) -> Point3 {
    let ring = t.major_radius + t.minor_radius * v.cos();
    t.frame.origin
        + t.frame.x * (ring * u.cos())
        + t.frame.y * (ring * u.sin())
        + t.frame.z * (t.minor_radius * v.sin())
}

/// On each of `k` tube circles, the scan's sign changes and the traced
/// curves' crossings; returns the pairs that disagree.
fn completeness(t: &Torus, other: &Surface, branches: &[Curve3]) -> Vec<(f64, usize, usize)> {
    let mut bad = Vec::new();
    let pcurves: Vec<Vec<(f64, f64)>> = branches
        .iter()
        .map(|b| {
            let Curve3::ImplicitSection(s) = b else {
                panic!("expected an implicit section");
            };
            let n = 20_000;
            (0..=n)
                .map(|i| {
                    let p = s.curve.point(s.curve.end() * i as f64 / n as f64).unwrap();
                    (p.x, p.y)
                })
                .collect()
        })
        .collect();
    for i in 0..90 {
        let v = -PI + TAU * (i as f64 + 0.5) / 90.0;
        let mut scan = 0;
        let steps = 4000;
        let mut last = off(other, torus_point(t, -PI, v)).signum();
        for k in 1..=steps {
            let now = off(
                other,
                torus_point(t, -PI + TAU * k as f64 / steps as f64, v),
            )
            .signum();
            if now != last {
                scan += 1;
            }
            last = now;
        }
        let mut traced = 0;
        for pc in &pcurves {
            for w in pc.windows(2) {
                // Crossings of v in the traced curve's unwrapped angle.
                let (a, b) = (w[0].1, w[1].1);
                let k = ((a.min(b) - v) / TAU).ceil();
                let level = v + k * TAU;
                if level >= a.min(b) && level < a.max(b) {
                    traced += 1;
                }
            }
        }
        if scan != traced {
            bad.push((v, scan, traced));
        }
    }
    bad
}

#[test]
fn a_pipe_through_a_pipe_bend_is_two_loops() {
    // A straight pipe along x runs through the bend's tube where it is
    // tangent to x, entering and leaving the torus once each.
    let bend = torus(Point3::ZERO, Vec3::Z, 3.0, 1.0);
    let pipe = cylinder(Point3::new(0.0, 3.0, 0.1), Vec3::X, 0.5);
    assert_eq!(check(&bend, &pipe), 2);
    assert_eq!(check(&pipe, &bend), 2);
}

#[test]
fn a_skew_pipe_and_a_torus_agree_with_a_dense_scan() {
    let ring = torus(
        Point3::new(0.2, -0.1, 0.3),
        Vec3::new(0.1, 0.2, 1.0),
        3.0,
        1.0,
    );
    let pipe = cylinder(Point3::new(0.4, 2.0, 0.0), Vec3::new(1.0, 0.3, 0.4), 0.7);
    let n = check(&ring, &pipe);
    assert!(n >= 1);
    let Surface::Torus(t) = ring else {
        unreachable!()
    };
    let curve = exact_surface_intersection(&ring, &pipe).unwrap();
    let bad = completeness(&t, &pipe, &curve.branches);
    assert!(bad.is_empty(), "{bad:?}");
}

#[test]
fn a_cone_through_a_torus_off_its_axis() {
    let ring = torus(Point3::ZERO, Vec3::Z, 3.0, 1.0);
    let hopper = Surface::Cone(Cone {
        frame: frame(Point3::new(2.8, 0.3, -3.0), Vec3::new(0.2, 0.1, 1.0)),
        radius: 0.2,
        semi_angle: 0.15,
    });
    assert!(check(&ring, &hopper) >= 1);
    let Surface::Torus(t) = ring else {
        unreachable!()
    };
    let curve = exact_surface_intersection(&ring, &hopper).unwrap();
    let bad = completeness(&t, &hopper, &curve.branches);
    assert!(bad.is_empty(), "{bad:?}");
}

#[test]
fn an_elliptical_duct_through_a_torus() {
    let ring = torus(Point3::ZERO, Vec3::Z, 3.0, 1.0);
    let duct = Surface::EllipticalCylinder(EllipticalCylinder {
        frame: frame(Point3::new(0.0, -3.0, 0.0), Vec3::new(1.0, 0.1, 0.2)),
        semi_axis_x: 0.6,
        semi_axis_y: 0.3,
    });
    assert!(check(&ring, &duct) >= 2);
}

#[test]
fn two_linked_tori_meet_in_closed_curves() {
    // A chain whose links are too tight: each torus's centre circle runs
    // through the other's hole, 1.5 from the other's centre circle, so the
    // tubes (radius 1) cut into each other.
    let a = torus(Point3::ZERO, Vec3::Z, 3.0, 1.0);
    let b = torus(Point3::new(3.0, 0.0, 0.0), Vec3::Y, 1.5, 1.0);
    let n = check(&a, &b);
    assert!(n >= 2, "{n}");
    let Surface::Torus(t) = a else { unreachable!() };
    let curve = exact_surface_intersection(&a, &b).unwrap();
    let bad = completeness(&t, &b, &curve.branches);
    assert!(bad.is_empty(), "{bad:?}");
}

#[test]
fn apart_and_touching_tori_are_named() {
    let a = torus(Point3::ZERO, Vec3::Z, 3.0, 1.0);
    let far = torus(Point3::new(20.0, 0.0, 0.0), Vec3::X, 3.0, 1.0);
    assert_eq!(
        exact_surface_intersection(&a, &far),
        Err(ExactIntersectionRefusal::Disjoint)
    );
    // Two tori in one plane, their outer equators touching at one point.
    let beside = torus(Point3::new(8.0, 0.0, 0.0), Vec3::Z, 3.0, 1.0);
    assert_eq!(
        exact_surface_intersection(&a, &beside),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}

struct Lcg(u64);

impl Lcg {
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
}

#[test]
fn random_torus_and_cylinder_pairs_agree_with_a_dense_scan() {
    let mut rng = Lcg(0x0077_0119);
    let mut met = 0;
    for _ in 0..40 {
        let t = Torus {
            frame: frame(
                Point3::new(
                    rng.range(-0.5, 0.5),
                    rng.range(-0.5, 0.5),
                    rng.range(-0.5, 0.5),
                ),
                Vec3::new(rng.range(-0.4, 0.4), rng.range(-0.4, 0.4), 1.0),
            ),
            major_radius: rng.range(2.5, 4.0),
            minor_radius: rng.range(0.5, 1.2),
        };
        let pipe = cylinder(
            Point3::new(
                rng.range(-3.0, 3.0),
                rng.range(-3.0, 3.0),
                rng.range(-0.5, 0.5),
            ),
            Vec3::new(
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ),
            rng.range(0.2, 1.0),
        );
        match exact_surface_intersection(&Surface::Torus(t), &pipe) {
            Ok(curve) => {
                met += 1;
                let bad = completeness(&t, &pipe, &curve.branches);
                // A tube circle grazing the section right at a turning point
                // can disagree by the scan's resolution; none may disagree by
                // more than that.
                assert!(bad.len() <= 1, "{bad:?}");
            }
            Err(ExactIntersectionRefusal::Disjoint) => {
                let bad = completeness(&t, &pipe, &[]);
                assert!(bad.is_empty(), "missed: {bad:?}");
            }
            Err(other) => panic!("{other:?}"),
        }
    }
    assert!(met >= 10, "{met}");
}

// --- B9 for the section families ------------------------------------------

use axiolid_core::Interval;
use axiolid_nurbs::{
    exact_curve_surface_intersection, section_curve_surface_intersection, ExactCurveIntersection,
};
use axiolid_surface::{Plane, Sphere};

/// Sign changes of `off(surface)` along a curve, sampled densely.
fn scan_crossings(curve: &Curve3, span: Interval, surface: &Surface) -> usize {
    let n = 20_000;
    let mut count = 0;
    let mut last = None;
    for i in 0..=n {
        let t = span.start + (span.end - span.start) * i as f64 / n as f64;
        let Ok(p) = evaluate3(curve, t) else { continue };
        let s = match surface {
            Surface::Plane(q) => q.frame.z.normalize().dot(p - q.frame.origin),
            Surface::Sphere(q) => (p - q.frame.origin).length() - q.radius,
            other => off(other, p),
        }
        .signum();
        if let Some(l) = last {
            if l != s {
                count += 1;
            }
        }
        last = Some(s);
    }
    count
}

#[test]
fn a_traced_section_against_planes_and_spheres_agrees_with_sampling() {
    let bend = torus(Point3::ZERO, Vec3::Z, 3.0, 1.0);
    let pipe = cylinder(Point3::new(0.0, 3.0, 0.1), Vec3::X, 0.5);
    let curve = exact_surface_intersection(&bend, &pipe).unwrap();
    let cutters = [
        Surface::Plane(Plane {
            frame: frame(Point3::new(1.9, 3.0, 0.0), Vec3::new(1.0, 0.2, 0.1)),
        }),
        Surface::Plane(Plane {
            frame: frame(Point3::new(0.0, 3.0, 0.2), Vec3::Z),
        }),
        Surface::Sphere(Sphere {
            frame: frame(Point3::new(2.2, 3.1, 0.0), Vec3::Z),
            radius: 0.8,
        }),
    ];
    let mut total = 0;
    for branch in &curve.branches {
        let Curve3::ImplicitSection(s) = branch else {
            unreachable!()
        };
        let span = Interval::new(0.0, s.curve.end());
        for cutter in &cutters {
            let ExactCurveIntersection::Points(hits) =
                exact_curve_surface_intersection(branch, cutter).expect("hits")
            else {
                panic!("not contained");
            };
            for hit in &hits {
                let p = evaluate3(branch, hit.parameter.approx()).unwrap();
                assert!((p - hit.point).length() < 1e-9);
                let r = match cutter {
                    Surface::Plane(q) => q.frame.z.normalize().dot(p - q.frame.origin),
                    Surface::Sphere(q) => (p - q.frame.origin).length() - q.radius,
                    _ => unreachable!(),
                };
                assert!(r.abs() < 1e-9, "{r}");
            }
            let crossings = hits.iter().filter(|h| h.multiplicity == 1).count();
            assert_eq!(crossings, scan_crossings(branch, span, cutter));
            total += crossings;
        }
        // It lies on both of its own surfaces.
        assert_eq!(
            exact_curve_surface_intersection(branch, &pipe),
            Ok(ExactCurveIntersection::Contained)
        );
    }
    assert!(total >= 4, "{total}");
}

#[test]
fn a_ruled_section_against_a_plane_by_its_span() {
    // A pipe tee: the thin pipe's loops on the wide one.
    let wide = cylinder(Point3::ZERO, Vec3::Z, 2.0);
    let thin = cylinder(Point3::new(0.0, 0.0, 1.0), Vec3::X, 0.7);
    let curve = exact_surface_intersection(&wide, &thin).unwrap();
    let cutter = Surface::Plane(Plane {
        frame: frame(Point3::new(0.0, 0.0, 1.2), Vec3::new(0.1, 0.3, 1.0)),
    });
    let mut total = 0;
    for (branch, span) in curve.branches.iter().zip(&curve.spans) {
        let span = span.expect("a ruled piece carries its span");
        let ExactCurveIntersection::Points(hits) =
            section_curve_surface_intersection(branch, span, &cutter).expect("hits")
        else {
            panic!("not contained");
        };
        for hit in &hits {
            let t = hit.parameter.approx();
            assert!(t >= span.start.min(span.end) - 1e-9 && t <= span.start.max(span.end) + 1e-9);
            let p = evaluate3(branch, t).unwrap();
            assert!((p - hit.point).length() < 1e-8, "{p:?} vs {:?}", hit.point);
        }
        assert_eq!(hits.len(), scan_crossings(branch, span, &cutter));
        total += hits.len();
    }
    assert!(total >= 2, "{total}");
}

#[test]
fn branches_crossing_where_a_pipe_touches_a_torus_meet_at_a_vertex() {
    // A pipe of radius 1/2 along y, inside the torus's tube and touching
    // it at the outer equator (3, 0, 0). Across the tube the pipe curves
    // more than the torus, along it less, so there the two branches of the
    // section cross: a saddle of the field on its zero set.
    let t = Torus {
        frame: frame(Point3::ZERO, Vec3::Z),
        major_radius: 2.0,
        minor_radius: 1.0,
    };
    let pipe = cylinder(Point3::new(2.5, 0.0, 0.0), Vec3::Y, 0.5);
    let torus_surface = Surface::Torus(t);
    check(&torus_surface, &pipe);
    let curve = exact_surface_intersection(&torus_surface, &pipe).expect("a traced section");
    let touch = Point3::new(3.0, 0.0, 0.0);
    // Every branch into the crossing ends there.
    let mut at_touch = 0;
    for (branch, span) in curve.branches.iter().zip(&curve.spans) {
        let span = span.unwrap();
        for t in [span.start, span.end] {
            if (evaluate3(branch, t).unwrap() - touch).length() < 1e-12 {
                at_touch += 1;
            }
        }
    }
    assert_eq!(at_touch, 4, "{} branches", curve.branches.len());
    // The straight bridges into the crossing stay on both surfaces.
    let mut bridges = 0;
    for branch in &curve.branches {
        let Curve3::ImplicitSection(s) = branch else {
            panic!("expected an implicit section");
        };
        for (i, cell) in s.curve.cells.iter().enumerate() {
            if cell.bridge.is_none() {
                continue;
            }
            bridges += 1;
            for k in 0..=50 {
                let p = evaluate3(branch, i as f64 + k as f64 / 50.0).unwrap();
                for surface in [&torus_surface, &pipe] {
                    assert!(off(surface, p).abs() < 1e-9, "bridge point {p:?}");
                }
            }
        }
    }
    assert_eq!(bridges, 4);
    let bad = completeness(&t, &pipe, &curve.branches);
    assert!(bad.is_empty(), "{bad:?}");
}

#[test]
fn bridges_into_a_crossing_stay_on_both_surfaces() {
    // The crossing of a pipe touching a torus's tube from inside, traced on
    // the pipe: there the torus's equation cancels large terms, rounding
    // hides more, and the bridges into the crossing are long enough for
    // their shape to matter.
    let t = Surface::Torus(Torus {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::Z,
            y: Vec3::X,
            z: Vec3::Y,
        },
        major_radius: 4.0,
        minor_radius: 1.0,
    });
    let c = 0.5f64.sqrt();
    let pipe = Surface::Cylinder(Cylinder {
        frame: Frame3 {
            // Its heights start 8 below the torus, as a pipe from z = -8
            // up is read.
            origin: Point3::new(-(4.0 + c) + 0.5 * c, 0.5 * c, -8.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: 0.5,
    });
    let touch = Point3::new(-(4.0 + c), c, 0.0);
    let window = (Point2::new(0.0, 0.0), Point2::new(PI, 16.0));
    let curves = implicit_surface_intersection(&pipe, &t, Some(window)).expect("a traced section");
    let (mut bridges, mut at_touch, mut longest) = (0, 0, 0.0f64);
    for s in &curves {
        let branch = Curve3::ImplicitSection(s.clone());
        for t_end in [0.0, s.curve.end()] {
            if (evaluate3(&branch, t_end).unwrap() - touch).length() < 1e-12 {
                at_touch += 1;
            }
        }
        for (i, cell) in s.curve.cells.iter().enumerate() {
            if cell.bridge.is_none() {
                continue;
            }
            bridges += 1;
            longest = longest.max((cell.to - cell.from).abs());
            for k in 0..=100 {
                let p = evaluate3(&branch, i as f64 + k as f64 / 100.0).unwrap();
                for surface in [&t, &pipe] {
                    assert!(
                        off(surface, p).abs() < 1e-10,
                        "bridge point {p:?} is {} off",
                        off(surface, p)
                    );
                }
            }
        }
    }
    assert_eq!((bridges, at_touch), (4, 4));
    assert!(longest > 1e-6, "bridges of {longest} test nothing");
}

#[test]
fn a_pipe_touching_the_top_of_a_torus_tube_meets_it_in_a_tacnode() {
    // A pipe of radius 1/2 along z inside a torus's tube, touching it at
    // the tube's top (-4, 1, 0): across the tube the pipe curves more, and
    // along it the torus parts from the shared tangent plane only as
    // z^4 / 128. The Hessian is singular there; the two branches touch
    // each other (x ~ z^2) and all four halves end at the contact.
    let t = Surface::Torus(Torus {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::Z,
            y: Vec3::X,
            z: Vec3::Y,
        },
        major_radius: 4.0,
        minor_radius: 1.0,
    });
    let pipe = Surface::Cylinder(Cylinder {
        frame: Frame3 {
            origin: Point3::new(-4.0, 0.5, -8.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: 0.5,
    });
    let curves = implicit_surface_intersection(&t, &pipe, None).expect("a traced section");
    let mut meeting: Vec<Point3> = Vec::new();
    let mut bridges = 0;
    for s in &curves {
        let branch = Curve3::ImplicitSection(s.clone());
        for t_end in [0.0, s.curve.end()] {
            let p = evaluate3(&branch, t_end).unwrap();
            if (p - Point3::new(-4.0, 1.0, 0.0)).length() < 0.05 {
                meeting.push(p);
            }
        }
        for (i, cell) in s.curve.cells.iter().enumerate() {
            let steps = if cell.bridge.is_some() {
                bridges += 1;
                200
            } else {
                4
            };
            for k in 0..=steps {
                let p = evaluate3(&branch, i as f64 + k as f64 / steps as f64).unwrap();
                for surface in [&t, &pipe] {
                    assert!(
                        off(surface, p).abs() < 1e-9,
                        "{p:?} is {} off",
                        off(surface, p)
                    );
                }
            }
        }
    }
    assert_eq!((meeting.len(), bridges), (4, 4));
    assert!(meeting.iter().all(|p| (*p - meeting[0]).length() < 1e-12));
}
