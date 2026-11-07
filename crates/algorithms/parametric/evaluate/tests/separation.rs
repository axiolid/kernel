//! The certified separation of a carrier's parameters over a box (#414).
//!
//! `sigma` must stay below the smallest singular value of the Jacobian in
//! scaled coordinates at every sampled point, and below the secant ratio
//! `|S(x) - S(y)| / |D (x - y)|` of every sampled pair; and a box that
//! touches a pole, an apex, a crease or is too wide has none.

mod common;

use axiolid_core::{Point3, Scalar};
use axiolid_curve::{BSplineSurface, Carrier, RuledCarrier};
use axiolid_evaluate::bound::SurfaceBoundOracle;

/// A deterministic pseudo-random number in `[0, 1)`.
fn next(state: &mut u64) -> Scalar {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as Scalar / (1_u64 << 53) as Scalar
}

/// A box of half widths `(hu, hv)` about the centre of a case's box.
fn about(case: &common::Case, hu: Scalar, hv: Scalar) -> ((Scalar, Scalar), (Scalar, Scalar)) {
    let (um, vm) = (0.5 * (case.u.0 + case.u.1), 0.5 * (case.v.0 + case.v.1));
    ((um - hu, um + hu), (vm - hv, vm + hv))
}

/// The smallest singular value of `[S_u / a, S_v / b]`.
fn smallest_singular(carrier: &Carrier, u: Scalar, v: Scalar, scale: [Scalar; 2]) -> Scalar {
    let jet = carrier.jet(u, v);
    let (su, sv) = (jet.u / scale[0], jet.v / scale[1]);
    let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let mean = 0.5 * (a + c);
    (mean - (0.25 * (a - c) * (a - c) + b * b).sqrt()).sqrt()
}

#[test]
fn sigma_is_below_the_scaled_singular_values_and_the_secant_ratios() {
    let mut rng = 7_u64;
    for case in common::cases() {
        // A rational spline's quotient-rule bounds are loose; its box is
        // narrower.
        let width = if case.name == "spline" {
            1.0 / 12.0
        } else {
            1.0
        };
        let (u, v) = about(&case, 0.12 * width, 0.1 * width);
        let oracle = SurfaceBoundOracle::for_carrier(&case.carrier).unwrap();
        let sep = oracle
            .separation(u, v)
            .unwrap_or_else(|| panic!("{} has no separation", case.name));
        assert!(sep.sigma > 0.0 && sep.sigma <= 1.0, "{}", case.name);
        let point = |rng: &mut u64| (u.0 + (u.1 - u.0) * next(rng), v.0 + (v.1 - v.0) * next(rng));
        for _ in 0..2000 {
            let (pu, pv) = point(&mut rng);
            let s = smallest_singular(&case.carrier, pu, pv, sep.scale);
            assert!(
                sep.sigma <= s * (1.0 + 1e-12),
                "{}: sigma {} above the singular value {s}",
                case.name,
                sep.sigma
            );
            let (qu, qv) = point(&mut rng);
            let (p, q) = (
                case.carrier.jet(pu, pv).point,
                case.carrier.jet(qu, qv).point,
            );
            let dx = ((qu - pu) * sep.scale[0]).hypot((qv - pv) * sep.scale[1]);
            if dx > 1e-9 {
                let ratio = (q - p).length() / dx;
                assert!(
                    sep.sigma <= ratio * (1.0 + 1e-12),
                    "{}: sigma {} above the secant ratio {ratio}",
                    case.name,
                    sep.sigma
                );
            }
        }
    }
}

#[test]
fn sigma_is_not_slack_by_more_than_the_box_allows() {
    // On a tiny box the separation is the centre's own singular value.
    for case in common::cases() {
        let (u, v) = about(&case, 1e-6, 1e-6);
        let oracle = SurfaceBoundOracle::for_carrier(&case.carrier).unwrap();
        let sep = oracle.separation(u, v).unwrap();
        let mid = (0.5 * (u.0 + u.1), 0.5 * (v.0 + v.1));
        let s = smallest_singular(&case.carrier, mid.0, mid.1, sep.scale);
        assert!(sep.sigma <= s && sep.sigma >= s - 1e-4, "{}", case.name);
    }
}

fn refused(carrier: &Carrier, u: (Scalar, Scalar), v: (Scalar, Scalar)) -> bool {
    SurfaceBoundOracle::for_carrier(carrier)
        .unwrap()
        .separation(u, v)
        .is_none()
}

#[test]
fn a_pole_an_apex_and_a_wide_cylinder_box_have_no_separation() {
    let sphere = Carrier::Sphere {
        frame: common::skewed(),
        radius: 0.7,
    };
    assert!(!refused(&sphere, (0.2, 0.4), (0.2, 0.4)));
    // A box on the pole.
    assert!(refused(
        &sphere,
        (0.2, 0.4),
        (1.2, core::f64::consts::FRAC_PI_2)
    ));
    // A cone whose box holds the apex, `radius + slope v = 0` at `v = -1`.
    let cone = Carrier::Ruled(RuledCarrier {
        frame: common::skewed(),
        x_radius: 0.4,
        y_radius: 0.4,
        slope: 0.4,
    });
    assert!(!refused(&cone, (0.2, 0.4), (0.0, 0.2)));
    assert!(refused(&cone, (0.2, 0.4), (-1.1, -0.9)));
    // A cylinder box at least a half turn wide.
    let cylinder = Carrier::Ruled(RuledCarrier {
        frame: common::rotated(),
        x_radius: 0.5,
        y_radius: 0.5,
        slope: 0.0,
    });
    assert!(!refused(&cylinder, (0.0, 0.6), (0.0, 1.0)));
    assert!(refused(&cylinder, (0.0, core::f64::consts::PI), (0.0, 1.0)));
}

/// A quadratic patch with a double knot at `u = 0.5`, a crease there.
fn creased() -> BSplineSurface {
    let mut s = common::patch();
    s.control_points = (0..5)
        .map(|i| {
            (0..3)
                .map(|j| {
                    Point3::new(
                        i as Scalar * 0.25,
                        j as Scalar * 0.5,
                        (i % 2) as Scalar * 0.2,
                    )
                })
                .collect()
        })
        .collect();
    s.u_knots = vec![0.0, 0.5, 1.0];
    s.u_multiplicities = vec![3, 2, 3];
    s.weights = None;
    s
}

#[test]
fn a_crease_or_the_edge_of_a_spline_has_no_separation() {
    let carrier = Carrier::Spline(Box::new(creased()));
    assert!(!refused(&carrier, (0.1, 0.3), (0.3, 0.6)));
    assert!(refused(&carrier, (0.4, 0.6), (0.3, 0.6)));
    // Past the knot domain.
    assert!(refused(&carrier, (0.9, 1.1), (0.3, 0.6)));
}

#[test]
fn a_degenerate_or_reversed_box_is_read_in_order_and_non_finite_has_none() {
    let case = &common::cases()[0];
    let oracle = SurfaceBoundOracle::for_carrier(&case.carrier).unwrap();
    let a = oracle.separation((0.0, 0.1), (0.0, 0.1)).unwrap();
    let b = oracle.separation((0.1, 0.0), (0.1, 0.0)).unwrap();
    assert_eq!(a, b);
    assert!(oracle.separation((0.0, Scalar::NAN), (0.0, 0.1)).is_none());
}
