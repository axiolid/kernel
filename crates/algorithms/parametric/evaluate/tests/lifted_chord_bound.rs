//! Certified chord bounds of lifted pcurves on analytic carriers (#414).
//!
//! Each bound is checked against a dense sample: no sampled point may lie
//! further from the chord of the evaluated ends than the bound, and the
//! bound must stay within ten times the sampled maximum, so it is a usable
//! certificate and not merely a safe one. Where the curve has an independent
//! closed form, the points are computed from it; the lifted evaluator's own
//! points are the sample otherwise.

mod common;

use axiolid_core::{Frame3, Point2, Point3, Scalar, Vec3};
use axiolid_curve::{
    AngleGraph2, Axis, Basis, Branch, Carrier, Circle3, Curve2, Curve3, Ellipse3, Field2,
    ImplicitCell, ImplicitCurve2, ImplicitSection3, LiftedCurve2, Polyline3, QuadraticGraph2,
    RuledCarrier, RuledSection3, SeriesField2, TorusCarrier, TorusSection3, Trig2,
};
use axiolid_evaluate::bound::{chord_bound2, continuity_breaks2};
use axiolid_evaluate::curve::evaluate2;
use core::f64::consts::{PI, TAU};

const SAMPLES: usize = 4000;

/// How far above the sampled maximum a bound may sit and stay useful.
const USEFUL: Scalar = 10.0;

/// The largest distance of the sampled points from the chord of the ends.
fn sampled(point: &impl Fn(Scalar) -> Point2, lo: Scalar, hi: Scalar) -> Scalar {
    let (p, q) = (point(lo), point(hi));
    let chord = q - p;
    (1..SAMPLES)
        .map(|i| {
            let t = lo + (hi - lo) * i as Scalar / SAMPLES as Scalar;
            let r = point(t) - p;
            let s = (r.dot(chord) / chord.length_squared()).clamp(0.0, 1.0);
            (r - chord * s).length()
        })
        .fold(0.0, Scalar::max)
}

/// A lifted pcurve whose guide is `unwrapped` at 129 evenly spaced `t`.
fn lifted(
    curve: Curve3,
    carrier: Carrier,
    (start, end): (Scalar, Scalar),
    unwrapped: impl Fn(Scalar) -> Point2,
) -> Curve2 {
    let n = 129;
    Curve2::Lifted(LiftedCurve2 {
        curve: Box::new(curve),
        carrier,
        start,
        end,
        guide: (0..n)
            .map(|i| unwrapped(start + (end - start) * i as Scalar / (n - 1) as Scalar))
            .collect(),
    })
}

/// Check every span: bounded, sound against the samples of `point`, useful.
///
/// On a flat carrier the bound stays within `USEFUL` times the sample. On a
/// curved one the pcurve's second derivative is the difference of two
/// vectors the bound adds in magnitude (`c''` and the carrier's own
/// curvature along `p'`, nearly cancelling for a curve that hugs it), so
/// there it only has to stay within `USEFUL` times the deviation of a bend
/// of 1 where the sample is smaller, as the graph bounds' tests do.
fn check(
    curve: &Curve2,
    point: impl Fn(Scalar) -> Point2,
    spans: &[(Scalar, Scalar)],
    what: &str,
    flat: bool,
) {
    for &(lo, hi) in spans {
        let sample = sampled(&point, lo, hi);
        let bound =
            chord_bound2(curve, lo, hi).unwrap_or_else(|| panic!("{what} [{lo}, {hi}]: no bound"));
        assert!(
            sample <= bound,
            "{what} [{lo}, {hi}]: sampled {sample} above the bound {bound}"
        );
        let reference = if flat {
            sample
        } else {
            sample.max((hi - lo) * (hi - lo) / 8.0)
        };
        assert!(
            bound <= USEFUL * reference + 1e-12,
            "{what} [{lo}, {hi}]: bound {bound} is over {USEFUL}x the reference {reference}"
        );
    }
}

fn evaluated(curve: &Curve2) -> impl Fn(Scalar) -> Point2 + '_ {
    |t| evaluate2(curve, t).expect("evaluates")
}

/// `raw` carried to the turn nearest `reference`.
fn near(raw: Scalar, reference: Scalar) -> Scalar {
    raw + ((reference - raw) / TAU).round() * TAU
}

#[test]
fn a_circle_on_a_rotated_plane() {
    let carrier = common::rotated();
    let frame = Frame3 {
        origin: carrier.origin + carrier.x * 0.2 + carrier.y * 0.1,
        ..carrier
    };
    let curve = lifted(
        Curve3::Circle(Circle3 { frame, radius: 0.6 }),
        Carrier::Plane(carrier),
        (0.0, TAU),
        |t| Point2::new(0.2 + 0.6 * t.cos(), 0.1 + 0.6 * t.sin()),
    );
    let point = |t: Scalar| Point2::new(0.2 + 0.6 * t.cos(), 0.1 + 0.6 * t.sin());
    check(
        &curve,
        point,
        &[(0.0, 0.5), (1.0, 2.0), (0.0, 3.0), (4.0, 6.2), (5.0, 5.001)],
        "circle on a plane",
        true,
    );
    // The same points as the evaluator reads them.
    check(
        &curve,
        evaluated(&curve),
        &[(0.3, 1.1)],
        "evaluated circle",
        true,
    );
}

#[test]
fn an_oblique_ellipse_on_a_cylinder_across_the_seam() {
    let c = common::rotated();
    let (r, tan) = (0.7, 0.5);
    let norm = (1.0 + tan * tan as Scalar).sqrt();
    let y = (c.y + c.z * tan) / norm;
    let ellipse = Curve3::Ellipse(Ellipse3 {
        frame: Frame3 {
            origin: c.origin,
            x: c.x,
            y,
            z: c.x.cross(y),
        },
        semi_axis_x: r,
        semi_axis_y: r * norm,
    });
    let carrier = Carrier::Ruled(RuledCarrier {
        frame: c,
        x_radius: r,
        y_radius: r,
        slope: 0.0,
    });
    let point = |t: Scalar| Point2::new(t, r * tan * t.sin());
    let curve = lifted(ellipse, carrier, (0.0, TAU), point);
    // `(2.9, 3.4)` crosses `u = pi`, where the principal angle wraps.
    check(
        &curve,
        point,
        &[(0.3, 0.9), (2.9, 3.4), (5.5, 6.0), (0.0, 3.0)],
        "ellipse on a cylinder",
        false,
    );
}

#[test]
fn a_tilted_great_circle_on_a_sphere() {
    let c = common::rotated();
    let (radius, tilt) = (0.8, 0.5_f64);
    let y = c.y * tilt.cos() + c.z * tilt.sin();
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: c.origin,
            x: c.x,
            y,
            z: c.x.cross(y),
        },
        radius,
    });
    let carrier = Carrier::Sphere { frame: c, radius };
    let point = move |t: Scalar| {
        let raw = (t.sin() * tilt.cos()).atan2(t.cos());
        Point2::new(near(raw, t), (t.sin() * tilt.sin()).asin())
    };
    let curve = lifted(circle, carrier, (0.0, TAU), point);
    check(
        &curve,
        point,
        &[(0.2, 0.7), (1.3, 1.8), (2.9, 3.4), (4.0, 5.5), (0.0, 3.0)],
        "great circle on a sphere",
        false,
    );
}

fn trig(constant: Scalar, cos: Scalar, sin: Scalar) -> Trig2 {
    Trig2 {
        constant,
        cos,
        sin,
        cos2: 0.0,
        sin2: 0.0,
    }
}

/// A plane's cut of a torus, read on that plane: the torus point at
/// `(angle(t), t)` projected on the plane's axes.
#[test]
fn a_plane_section_of_a_torus_on_its_plane() {
    let frame = common::rotated();
    let (major, minor) = (0.9, 0.3);
    // The plane `n . (P - O) = d` in the torus frame.
    let raw = Vec3::new(0.3, 0.2, 0.9);
    let n = raw / raw.length();
    let d = 0.05;
    let world = frame.x * n.x + frame.y * n.y + frame.z * n.z;
    let first = if world.x.abs() < 0.8 {
        Vec3::X
    } else {
        Vec3::Y
    };
    let px = world.cross(first);
    let px = px / px.length();
    let py = world.cross(px);
    let plane = Frame3 {
        origin: frame.origin + world * d,
        x: px,
        y: py,
        z: world,
    };
    // `(R + r cos t)(n_x cos u + n_y sin u) + n_z r sin t = d`.
    let graph = |branch| AngleGraph2 {
        a: trig(major * n.x, minor * n.x, 0.0),
        b: trig(major * n.y, minor * n.y, 0.0),
        c: trig(d, 0.0, -n.z * minor),
        branch,
    };
    for branch in [Branch::Plus, Branch::Minus] {
        let g = graph(branch);
        let section = Curve3::TorusSection(TorusSection3 {
            torus: TorusCarrier {
                frame,
                major_radius: major,
                minor_radius: minor,
            },
            graph: g,
        });
        let point = move |t: Scalar| {
            let u = g.angle(t).expect("the plane cuts the tube at t");
            let ring = major + minor * t.cos();
            let p = frame.origin
                + (frame.x * u.cos() + frame.y * u.sin()) * ring
                + frame.z * (minor * t.sin());
            let r = p - plane.origin;
            // The point is on the plane.
            assert!(r.dot(world).abs() < 1e-12);
            Point2::new(r.dot(px), r.dot(py))
        };
        let curve = lifted(section, Carrier::Plane(plane), (-0.8, 0.8), point);
        check(
            &curve,
            point,
            &[(-0.6, -0.5), (-0.1, 0.1), (0.3, 0.6), (-0.7, 0.7)],
            &format!("torus section {branch:?}"),
            true,
        );
    }
}

/// A sphere's cut of a cone, read on that sphere.
#[test]
fn a_sphere_section_of_a_cone_on_its_sphere() {
    let frame = common::rotated();
    let (r0, s) = (0.5, 0.3);
    let (cx, cy, v0, sphere_radius) = (0.2, 0.1, 0.3, 0.9);
    // `|P - C|^2 = R^2` on `P = O + rho (cos t X + sin t Y) + v Z`, `rho = r0
    // + s v`, a quadratic in `v`.
    let graph = |branch| QuadraticGraph2 {
        a: trig(s * s + 1.0, 0.0, 0.0),
        b: trig(2.0 * r0 * s - 2.0 * v0, -2.0 * s * cx, -2.0 * s * cy),
        c: trig(
            r0 * r0 + cx * cx + cy * cy + v0 * v0 - sphere_radius * sphere_radius,
            -2.0 * r0 * cx,
            -2.0 * r0 * cy,
        ),
        branch,
    };
    let centre = frame.origin + frame.x * cx + frame.y * cy + frame.z * v0;
    let sphere = Frame3 {
        origin: centre,
        ..frame
    };
    // The other root runs within a tenth of a radian of the sphere's pole.
    {
        let branch = Branch::Plus;
        let g = graph(branch);
        let section = Curve3::RuledSection(RuledSection3 {
            carrier: RuledCarrier {
                frame,
                x_radius: r0,
                y_radius: r0,
                slope: s,
            },
            graph: g,
        });
        let point = move |t: Scalar| {
            // The textbook root of the quadratic.
            let (a, b, c) = (g.a.value(t), g.b.value(t), g.c.value(t));
            let v = (-b + branch.sign() * (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
            let rho = r0 + s * v;
            let p = frame.origin + (frame.x * t.cos() + frame.y * t.sin()) * rho + frame.z * v;
            let l = p - centre;
            let (x, y, z) = (l.dot(frame.x), l.dot(frame.y), l.dot(frame.z));
            assert!(
                (l.length() - sphere_radius).abs() < 1e-12,
                "{t}: {} vs {sphere_radius}, v {v}",
                l.length()
            );
            Point2::new(near(y.atan2(x), t), z.atan2(x.hypot(y)))
        };
        let curve = lifted(
            section,
            Carrier::Sphere {
                frame: sphere,
                radius: sphere_radius,
            },
            (0.0, TAU),
            point,
        );
        check(
            &curve,
            point,
            &[(0.2, 0.5), (1.3, 1.7), (2.9, 3.4), (5.0, 6.0)],
            &format!("cone section {branch:?}"),
            false,
        );
    }
}

/// `v - u^2 = 0`, one regular cell over `u` in `[-0.5, 0.5]`.
fn parabola() -> ImplicitCurve2 {
    ImplicitCurve2 {
        field: Field2::Series(SeriesField2 {
            u: Basis::Power,
            v: Basis::Power,
            coefficients: vec![vec![0.0, 1.0], vec![0.0], vec![-1.0]],
        }),
        cells: vec![ImplicitCell {
            axis: Axis::U,
            from: -0.5,
            to: 0.5,
            low: -0.1,
            high: 0.4,
            bridge: None,
        }],
    }
}

#[test]
fn an_implicit_section_on_a_tilted_plane() {
    let plane = common::rotated();
    let section = Curve3::ImplicitSection(ImplicitSection3 {
        carrier: Carrier::Plane(plane),
        curve: parabola(),
    });
    // The parameter is `u = -0.5 + t` over the cell.
    let point = |t: Scalar| Point2::new(t - 0.5, (t - 0.5) * (t - 0.5));
    let curve = lifted(section, Carrier::Plane(plane), (0.0, 1.0), point);
    check(
        &curve,
        point,
        &[(0.0, 0.5), (0.1, 0.2), (0.25, 0.75), (0.0, 1.0)],
        "parabola on a plane",
        true,
    );
    // A stretch across a cell join is the caller's to split.
    assert!(continuity_breaks2(&curve, 1).is_empty());
}

#[test]
fn a_span_across_an_implicit_cell_join_is_refused() {
    let mut two = parabola();
    two.cells.push(ImplicitCell {
        axis: Axis::U,
        from: 0.5,
        to: 0.9,
        low: 0.1,
        high: 1.0,
        bridge: None,
    });
    let plane = common::rotated();
    let curve = lifted(
        Curve3::ImplicitSection(ImplicitSection3 {
            carrier: Carrier::Plane(plane),
            curve: two,
        }),
        Carrier::Plane(plane),
        (0.0, 2.0),
        |t| Point2::new(t, t),
    );
    assert_eq!(continuity_breaks2(&curve, 1), vec![1.0]);
    assert!(chord_bound2(&curve, 0.5, 1.5).is_none());
    assert!(chord_bound2(&curve, 0.2, 0.8).is_some());
}

/// The oblique ellipse of the seam test with a chosen guide.
fn ellipse_with_guide(guide: Vec<Point2>, start: Scalar, end: Scalar) -> Curve2 {
    let c = common::rotated();
    let (r, tan) = (0.7, 0.5);
    let norm = (1.0 + tan * tan as Scalar).sqrt();
    let y = (c.y + c.z * tan) / norm;
    Curve2::Lifted(LiftedCurve2 {
        curve: Box::new(Curve3::Ellipse(Ellipse3 {
            frame: Frame3 {
                origin: c.origin,
                x: c.x,
                y,
                z: c.x.cross(y),
            },
            semi_axis_x: r,
            semi_axis_y: r * norm,
        })),
        carrier: Carrier::Ruled(RuledCarrier {
            frame: c,
            x_radius: r,
            y_radius: r,
            slope: 0.0,
        }),
        start,
        end,
        guide,
    })
}

#[test]
fn a_guide_that_changes_the_turn_inside_a_span_is_refused() {
    let v = |t: Scalar| 0.35 * t.sin();
    let honest = vec![
        Point2::new(0.0, v(0.0)),
        Point2::new(1.5, v(1.5)),
        Point2::new(3.0, v(3.0)),
    ];
    let curve = ellipse_with_guide(honest.clone(), 0.0, 3.0);
    assert!(chord_bound2(&curve, 0.5, 2.5).is_some());
    // The same ends, the middle node a turn's worth away: the evaluator's
    // nearest turn flips inside the span.
    let mut skewed = honest;
    skewed[1].x += 4.0;
    let curve = ellipse_with_guide(skewed, 0.0, 3.0);
    assert!(chord_bound2(&curve, 0.5, 2.5).is_none());
    // An empty guide cannot place any point.
    assert!(chord_bound2(&ellipse_with_guide(Vec::new(), 0.0, 3.0), 0.5, 2.5).is_none());
}

#[test]
fn a_span_through_a_sphere_pole_is_refused() {
    let c = common::rotated();
    let radius = 0.8;
    // A meridian: the circle through the poles.
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: c.origin,
            x: c.x,
            y: c.z,
            z: c.x.cross(c.z),
        },
        radius,
    });
    let curve = lifted(
        circle,
        Carrier::Sphere { frame: c, radius },
        (0.0, TAU),
        |t| {
            Point2::new(
                if t < PI / 2.0 || t > 1.5 * PI {
                    0.0
                } else {
                    PI
                },
                t,
            )
        },
    );
    // The pole is at `t = pi / 2`.
    assert!(chord_bound2(&curve, 1.0, 2.0).is_none());
    assert!(chord_bound2(&curve, 1.4, 1.7).is_none());
}

#[test]
fn a_spline_carrier_without_a_pair_and_a_curve_without_bounds_have_none() {
    let c = common::rotated();
    let circle = Curve3::Circle(Circle3 {
        frame: c,
        radius: 0.5,
    });
    let on_spline = lifted(
        circle.clone(),
        Carrier::Spline(Box::new(common::patch())),
        (0.0, 1.0),
        |t| Point2::new(t, t),
    );
    assert!(chord_bound2(&on_spline, 0.1, 0.2).is_none());
    // A polyline has no derivative bound for a lifted curve to carry.
    let polyline = Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        closed: false,
    });
    let on_plane = lifted(
        polyline,
        Carrier::Plane(Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        }),
        (0.0, 2.0),
        |t| Point2::new(t.min(1.0), (t - 1.0).max(0.0)),
    );
    assert!(chord_bound2(&on_plane, 0.1, 0.9).is_none());
    // Non-finite spans and a point span.
    assert!(chord_bound2(&on_plane, 0.1, Scalar::NAN).is_none());
    assert_eq!(chord_bound2(&on_plane, 0.5, 0.5), Some(0.0));
}
