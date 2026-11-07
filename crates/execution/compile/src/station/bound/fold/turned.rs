//! The turned section's frame in closed form (#406).
//!
//! Between two stations the section's orientation is interpolated as the
//! kernel does ([`oriented`](crate::station::oriented), `SectionFrame::oriented`):
//! `a(u) = lerp(axis_a, axis_b, u)` and `r(u) = lerp(ref_a, ref_b, u)`, the
//! unit inputs of `StationOrientation::unit_axes`. The turned up is `e = a /
//! |a|`; with `g = r - (r.e) e` the turned tangent, the section plane's
//! normal, is `n = g / |g|`, and the turned lateral is `y^ = e x n`. All are
//! components in the BASE frame `E = [T, L, U]`.
//!
//! The vertex `p = q + x y^ + y e` of a section advances along the plane's
//! normal at a rate that holds the turn of `y^` and `e` along `n`:
//!
//! - `c_x = y^'.n = -(r'.y^ - (r.e)(a'.y^) / |a|) / |g|`;
//! - `c_y = e'.n = (a'.n) / |a|`,
//!
//! with `a' = (axis_b - axis_a) / delta` and `r' = (ref_b - ref_a) / delta`
//! per unit of station distance (`x'` and `y'` along `y^` and `e` are
//! perpendicular to `n` and drop out).

use axiolid_core::{Scalar, Vec3};

/// Below this a vector has no direction to speak of.
const DEGENERATE: Scalar = 1e-12;

/// The turned frame at one fraction of a station interval.
#[derive(Debug, Clone, Copy)]
pub(super) struct TurnedFrame {
    /// The section plane's normal `n`, the turned tangent.
    pub(super) normal: Vec3,
    /// The turned lateral `y^`.
    pub(super) lateral: Vec3,
    /// The turned up `e`.
    pub(super) up: Vec3,
    /// `c_x = y^'.n`.
    pub(super) lateral_rate: Scalar,
    /// `c_y = e'.n`.
    pub(super) up_rate: Scalar,
    /// The rounding of `c_x` and `c_y`.
    pub(super) rate_error: [Scalar; 2],
}

impl TurnedFrame {
    /// The base frame itself: `n = X`, `y^ = Y`, `e = Z`, no turn.
    pub(super) const BASE: Self = Self {
        normal: Vec3::X,
        lateral: Vec3::Y,
        up: Vec3::Z,
        lateral_rate: 0.0,
        up_rate: 0.0,
        rate_error: [0.0; 2],
    };
}

/// The frame the orientations `a` and `b` (unit axis and unit reference
/// direction each) interpolate to at the fraction `u` of an interval `delta`
/// long, with its rates; `None` where `|a|` or `|g|` is not positive.
pub(super) fn turned_frame(
    a: (Vec3, Vec3),
    b: (Vec3, Vec3),
    u: Scalar,
    delta: Scalar,
) -> Option<TurnedFrame> {
    let axis = a.0.lerp(b.0, u);
    let reference = a.1.lerp(b.1, u);
    let (axis_rate, reference_rate) = ((b.0 - a.0) / delta, (b.1 - a.1) / delta);
    let size = axis.length();
    if size.is_nan() || size <= DEGENERATE {
        return None;
    }
    let up = axis / size;
    let g = reference - up * reference.dot(up);
    let g_size = g.length();
    if g_size.is_nan() || g_size <= DEGENERATE {
        return None;
    }
    let normal = g / g_size;
    let lateral = up.cross(normal);
    let up_rate = axis_rate.dot(normal) / size;
    let lateral_rate =
        -(reference_rate.dot(lateral) - reference.dot(up) * axis_rate.dot(lateral) / size) / g_size;
    let rounding = 64.0 * Scalar::EPSILON;
    let out = TurnedFrame {
        normal,
        lateral,
        up,
        lateral_rate,
        up_rate,
        rate_error: [
            rounding * (reference_rate.length() + reference.length() * axis_rate.length() / size)
                / g_size,
            rounding * axis_rate.length() / size,
        ],
    };
    (out.normal.is_finite()
        && out.lateral.is_finite()
        && out.up.is_finite()
        && out.lateral_rate.is_finite()
        && out.up_rate.is_finite()
        && out.rate_error.iter().all(|e| e.is_finite()))
    .then_some(out)
}

#[cfg(test)]
mod tests {
    use axiolid_model::StationOrientation;

    use super::*;
    use crate::station::oriented;

    fn orientation(axis: Vec3, reference: Vec3) -> (Vec3, Vec3) {
        StationOrientation::new(Some(axis), Some(reference))
            .unit_axes()
            .unwrap()
    }

    /// The kernel's own interpolated frame in the identity frame, and the
    /// rates by central differences of it.
    #[test]
    fn the_turned_frame_and_its_rates_match_the_kernels_orientation() {
        let identity = axiolid_reference::station::station_section3(
            &axiolid_curve::Curve3::Line(axiolid_curve::Line3 {
                origin: axiolid_core::Point3::ZERO,
                direction: Vec3::X,
            }),
            0.0,
        )
        .unwrap();
        let a = orientation(Vec3::new(0.2, -0.3, 1.0), Vec3::new(1.0, 0.4, 0.1));
        let b = orientation(Vec3::new(-0.5, 0.6, 0.7), Vec3::new(0.3, 1.0, -0.2));
        let (sa, sb) = (
            StationOrientation::new(Some(a.0), Some(a.1)),
            StationOrientation::new(Some(b.0), Some(b.1)),
        );
        let delta = 8.0;
        let at = |s: Scalar| oriented(&identity, &sa, &sb, s / delta).unwrap();
        let h = 1e-5;
        for s in [1.0, 3.5, 7.0] {
            let turned = turned_frame(a, b, s / delta, delta).unwrap();
            let kernel = at(s);
            assert!((turned.normal - kernel.tangent).length() < 1e-12);
            assert!((turned.lateral - kernel.lateral).length() < 1e-12);
            assert!((turned.up - kernel.up).length() < 1e-12);
            let rate = |axis: fn(&axiolid_reference::station::SectionFrame) -> Vec3| {
                (axis(&at(s + h)) - axis(&at(s - h))) / (2.0 * h)
            };
            let c_x = rate(|f| f.lateral).dot(kernel.tangent);
            let c_y = rate(|f| f.up).dot(kernel.tangent);
            assert!((turned.lateral_rate - c_x).abs() < 1e-8, "{s}");
            assert!((turned.up_rate - c_y).abs() < 1e-8, "{s}");
            assert!(c_x.abs() > 1e-3 && c_y.abs() > 1e-3);
        }
    }

    #[test]
    fn the_base_orientation_does_not_turn() {
        let base = orientation(Vec3::Z, Vec3::X);
        let turned = turned_frame(base, base, 0.3, 5.0).unwrap();
        assert!((turned.normal - Vec3::X).length() < 1e-15);
        assert!((turned.lateral - Vec3::Y).length() < 1e-15);
        assert!((turned.up - Vec3::Z).length() < 1e-15);
        assert_eq!((turned.lateral_rate, turned.up_rate), (0.0, 0.0));
    }

    #[test]
    fn an_axis_passing_through_zero_has_no_frame() {
        let a = orientation(Vec3::Z, Vec3::X);
        let b = orientation(-Vec3::Z, Vec3::X);
        assert!(turned_frame(a, b, 0.5, 5.0).is_none());
    }
}
