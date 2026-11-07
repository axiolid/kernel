//! Derivative bounds of the elementary surfaces over a parameter box.
//!
//! The formulas are shared by [`axiolid_surface::Surface`] and by a curve's
//! [`axiolid_curve::Carrier`], which name the same surfaces in the same
//! parameters (ADR 0077). A general ruled carrier `O + rx cos u X + ry sin u
//! Y + v Z`, `rx = x_radius + slope v` (a cylinder, an elliptical cylinder
//! or a cone), is bounded from the radii at the ends of the `v` range:
//! `|S_u|, |S_uu| <= s max(|rx|, |ry|)`, `|S_v| <= s sqrt(1 + slope^2)`,
//! `|S_uv| <= s |slope|` and `S_vv = 0`, `s` being the frame's stretch.

use axiolid_core::{Frame3, Scalar};
use axiolid_curve::RuledCarrier;

use super::{cos_range, frame_stretch3, max_abs_cos, max_abs_sin, SurfaceBounds};

/// A smooth analytic surface's bounds, the interpolation coefficients being
/// filled in by the caller.
fn smooth(du: Scalar, dv: Scalar, duu: Scalar, duv: Scalar, dvv: Scalar) -> SurfaceBounds {
    SurfaceBounds {
        du,
        dv,
        duu,
        duv,
        dvv,
        smooth: true,
        interpolation: [0.0; 3],
    }
}

/// `O + u X + v Y`.
pub(super) fn plane(frame: &Frame3) -> SurfaceBounds {
    smooth(frame.x.length(), frame.y.length(), 0.0, 0.0, 0.0)
}

/// `O + r cos v (cos u X + sin u Y) + r sin v Z`.
pub(super) fn sphere(frame: &Frame3, radius: Scalar, v0: Scalar, v1: Scalar) -> SurfaceBounds {
    let s = frame_stretch3(frame) * radius.abs();
    let cos = max_abs_cos(v0, v1);
    smooth(s * cos, s, s * cos, s * max_abs_sin(v0, v1), s)
}

/// A torus of major radius `major` and minor radius `minor`.
pub(super) fn torus(
    frame: &Frame3,
    major: Scalar,
    minor: Scalar,
    v0: Scalar,
    v1: Scalar,
) -> SurfaceBounds {
    let s = frame_stretch3(frame);
    let (cmin, cmax) = cos_range(v0, v1);
    let r = minor.abs();
    let ring = (major + r * cmin).abs().max((major + r * cmax).abs());
    smooth(
        s * ring,
        s * r,
        s * ring,
        s * r * max_abs_sin(v0, v1),
        s * r,
    )
}

/// A cylinder, elliptical cylinder or cone over `v` in `[v0, v1]`.
pub(super) fn ruled(k: &RuledCarrier, v0: Scalar, v1: Scalar) -> SurfaceBounds {
    let s = frame_stretch3(&k.frame);
    let radius = |v: Scalar| {
        (k.x_radius + k.slope * v)
            .abs()
            .max((k.y_radius + k.slope * v).abs())
    };
    let rho = radius(v0).max(radius(v1));
    smooth(
        s * rho,
        s * (1.0 + k.slope * k.slope).sqrt(),
        s * rho,
        s * k.slope.abs(),
        0.0,
    )
}
