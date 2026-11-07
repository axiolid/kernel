//! The derivative bounds of a curve's carrier (#414): dense jets read with
//! `Carrier::jet` never exceed the bounds of `SurfaceBoundOracle::for_carrier`,
//! and the bounds are not slack beyond the frames' own stretch.

mod common;

use axiolid_core::Scalar;
use axiolid_curve::{Carrier, RuledCarrier};
use axiolid_evaluate::bound::SurfaceBoundOracle;

const GRID: usize = 41;

#[test]
fn carrier_bounds_cover_the_sampled_jets_and_stay_close() {
    for case in common::cases() {
        let oracle = SurfaceBoundOracle::for_carrier(&case.carrier).unwrap();
        let b = oracle.bounds(case.u, case.v).unwrap();
        assert!(b.smooth, "{}", case.name);
        let mut seen = [0.0_f64; 5];
        for i in 0..GRID {
            for j in 0..GRID {
                let f = |lo: Scalar, hi: Scalar, k: usize| {
                    lo + (hi - lo) * k as Scalar / (GRID - 1) as Scalar
                };
                let jet = case
                    .carrier
                    .jet(f(case.u.0, case.u.1, i), f(case.v.0, case.v.1, j));
                let norms = [
                    jet.u.length(),
                    jet.v.length(),
                    jet.uu.length(),
                    jet.uv.length(),
                    jet.vv.length(),
                ];
                for (s, n) in seen.iter_mut().zip(norms) {
                    *s = s.max(n);
                }
            }
        }
        let bounds = [b.du, b.dv, b.duu, b.duv, b.dvv];
        for ((s, bound), what) in seen.iter().zip(bounds).zip(["u", "v", "uu", "uv", "vv"]) {
            assert!(*s <= bound, "{}: |S_{what}| {s} exceeds {bound}", case.name);
            // A skewed frame's stretch loosens the bound, a rational
            // spline's quotient rule much more; an analytic carrier's is
            // within a few times the truth.
            assert!(
                case.name == "spline" || bound <= 6.0 * s + 1e-12,
                "{}: |S_{what}| bound {bound} is slack against {s}",
                case.name
            );
        }
    }
}

#[test]
fn a_non_finite_carrier_has_no_oracle() {
    let mut carrier = RuledCarrier {
        frame: common::skewed(),
        x_radius: f64::NAN,
        y_radius: 1.0,
        slope: 0.0,
    };
    assert!(SurfaceBoundOracle::for_carrier(&Carrier::Ruled(carrier)).is_none());
    carrier.x_radius = 1.0;
    assert!(SurfaceBoundOracle::for_carrier(&Carrier::Ruled(carrier)).is_some());
}
