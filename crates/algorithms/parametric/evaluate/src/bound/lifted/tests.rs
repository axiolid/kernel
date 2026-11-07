//! The piece certificate of [`super`] below the split into pieces, against
//! the sampled derivatives of the lifted pcurve itself: the slope and bend
//! bounds must hold them, and on a flat carrier, where nothing cancels, be
//! close to them.

use super::*;
use crate::curve::{derivative2, second_derivative2};
use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_curve::{Circle3, Curve2, Curve3};
use core::f64::consts::TAU;

const SAMPLES: usize = 400;

fn frame() -> Frame3 {
    let (s, c) = 0.7_f64.sin_cos();
    let (s2, c2) = 0.4_f64.sin_cos();
    let x = Vec3::new(c, s, 0.0);
    let y = Vec3::new(-s, c, 0.0) * c2 + Vec3::new(0.0, 0.0, s2);
    Frame3 {
        origin: Point3::new(0.4, -0.3, 0.2),
        x,
        y,
        z: x.cross(y),
    }
}

/// A circle of `radius` about the carrier's origin tilted by `tilt` about
/// its x axis, read on the sphere (or plane) of the same frame, its guide
/// the unwrapped pcurve.
fn lifted(carrier: Carrier, tilt: Scalar, radius: Scalar) -> LiftedCurve2 {
    let f = frame();
    let y = f.y * tilt.cos() + f.z * tilt.sin();
    let sphere = matches!(carrier, Carrier::Sphere { .. });
    let point = move |t: Scalar| {
        if sphere {
            let raw = (t.sin() * tilt.cos()).atan2(t.cos());
            let u = raw + ((t - raw) / TAU).round() * TAU;
            Point2::new(u, (t.sin() * tilt.sin()).asin())
        } else {
            Point2::new(radius * t.cos(), radius * t.sin())
        }
    };
    let n = 129;
    LiftedCurve2 {
        curve: Box::new(Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: f.origin,
                x: f.x,
                y: if sphere { y } else { f.y },
                z: f.x.cross(if sphere { y } else { f.y }),
            },
            radius,
        })),
        carrier,
        start: 0.0,
        end: TAU,
        guide: (0..n)
            .map(|i| point(TAU * i as Scalar / (n - 1) as Scalar))
            .collect(),
    }
}

/// `sup |p'|` and `sup |p''|` over `[a, b]`.
fn sampled(l: &LiftedCurve2, a: Scalar, b: Scalar) -> (Scalar, Scalar) {
    let curve = Curve2::Lifted(l.clone());
    let mut out = (0.0_f64, 0.0_f64);
    for i in 0..=SAMPLES {
        let t = a + (b - a) * i as Scalar / SAMPLES as Scalar;
        out.0 = out.0.max(derivative2(&curve, t).unwrap().length());
        out.1 = out.1.max(second_derivative2(&curve, t).unwrap().length());
    }
    out
}

fn piece_of(l: &LiftedCurve2, a: Scalar, b: Scalar) -> PieceBound {
    let oracle = SurfaceBoundOracle::for_carrier(&l.carrier).unwrap();
    piece(l, &oracle, a, b).expect("certified")
}

#[test]
fn on_a_flat_carrier_the_piece_bounds_are_the_sampled_derivatives() {
    let l = lifted(Carrier::Plane(frame()), 0.0, 0.6);
    for (a, b) in [(0.2, 0.5), (1.0, 1.1), (4.0, 4.4)] {
        let (slope, bend) = sampled(&l, a, b);
        let p = piece_of(&l, a, b);
        assert!(
            slope <= p.slope && p.slope <= 1.01 * slope,
            "{a}: {slope} {p:?}"
        );
        let bound = p.bend.unwrap();
        assert!(bend <= bound && bound <= 1.01 * bend, "{a}: {bend} {p:?}");
        // The ends of a circle on its own plane are within rounding.
        assert!(p.start < 1e-12 && p.end < 1e-12, "{p:?}");
    }
}

#[test]
fn on_a_sphere_the_piece_bounds_hold_the_sampled_derivatives() {
    let radius = 0.8;
    let mut checked = 0;
    for tilt in [0.5, 1.0, 1.25] {
        let l = lifted(
            Carrier::Sphere {
                frame: frame(),
                radius,
            },
            tilt,
            radius,
        );
        for (a, b) in [(0.3, 0.4), (1.45, 1.55), (1.6, 1.7), (3.0, 3.1), (5.0, 5.2)] {
            let Some(p) =
                SurfaceBoundOracle::for_carrier(&l.carrier).and_then(|o| piece(&l, &o, a, b))
            else {
                // Near the pole the piece is refused; never silently held.
                continue;
            };
            checked += 1;
            let (slope, bend) = sampled(&l, a, b);
            assert!(
                slope <= p.slope,
                "tilt {tilt} [{a}, {b}]: {slope} over {p:?}"
            );
            let bound = p.bend.unwrap();
            assert!(bend <= bound, "tilt {tilt} [{a}, {b}]: {bend} over {p:?}");
        }
    }
    assert!(checked >= 6, "only {checked} pieces certified");
}

/// `|D p''| <= (C2 + |D p'|^2 (S_uu / a^2 + 2 S_uv / (a b) + S_vv / b^2)) /
/// sigma`, unscaled through `m`: each of the two terms alone.
#[test]
fn the_bend_adds_the_curves_bend_and_the_carriers_curvature_along_the_speed() {
    let carrier = Carrier::Sphere {
        frame: frame(),
        radius: 0.8,
    };
    let oracle = SurfaceBoundOracle::for_carrier(&carrier).unwrap();
    let boxed = [(0.3, 0.32), (0.4, 0.42)];
    let sep = oracle.separation(boxed[0], boxed[1]).unwrap();
    let b = oracle.bounds(boxed[0], boxed[1]).unwrap();
    let [x, y] = sep.scale;
    let m = x.min(y);
    let sigma = Iv::point(sep.sigma);
    let coefficient = b.duu / (x * x) + 2.0 * b.duv / (x * y) + b.dvv / (y * y);
    let close = |got: Scalar, want: Scalar| (got / want - 1.0).abs() < 1e-9;
    let curve = bend(&oracle, boxed, sep.scale, 0.0, 2.0, sigma, m).unwrap();
    assert!(close(curve, 2.0 / (sep.sigma * m)), "{curve}");
    let carrier = bend(&oracle, boxed, sep.scale, 3.0, 0.0, sigma, m).unwrap();
    assert!(
        close(carrier, 9.0 * coefficient / (sep.sigma * m)),
        "{carrier}"
    );
}
