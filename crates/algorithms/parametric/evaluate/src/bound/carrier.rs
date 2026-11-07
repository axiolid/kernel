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
use axiolid_curve::{Carrier, RuledCarrier};

use super::interval::{cos_iv, sin_iv, Iv};

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
    // Linear in the cosine, so right for either sign of `minor`.
    let ring = (major + minor * cmin)
        .abs()
        .max((major + minor * cmax).abs());
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

// --- enclosures of a point --------------------------------------------------

/// `a_0 + sum a_i b_i` of enclosures.
fn sum(first: Iv, terms: &[(Iv, Iv)]) -> Iv {
    terms.iter().fold(first, |acc, &(a, b)| acc.add(a.mul(b)))
}

/// The enclosure of an analytic carrier's point at `(u, v)`, the one a
/// closed form gives in outward-rounded arithmetic; `None` for a B-spline
/// carrier.
pub(super) fn point_enclosure(carrier: &Carrier, u: Scalar, v: Scalar) -> Option<[Iv; 3]> {
    let p = Iv::point;
    let (cu, su) = (cos_iv(u, u), sin_iv(u, u));
    let axes = |f: &Frame3| [f.x, f.y, f.z].map(|a| [a.x, a.y, a.z].map(p));
    let origin = |f: &Frame3| [f.origin.x, f.origin.y, f.origin.z].map(p);
    Some(match carrier {
        Carrier::Plane(f) => {
            let ([x, y, _], o) = (axes(f), origin(f));
            core::array::from_fn(|k| sum(o[k], &[(p(u), x[k]), (p(v), y[k])]))
        }
        Carrier::Ruled(c) => {
            let (ax, o) = (axes(&c.frame), origin(&c.frame));
            let rx = p(c.x_radius).add(p(c.slope).mul(p(v)));
            let ry = p(c.y_radius).add(p(c.slope).mul(p(v)));
            core::array::from_fn(|k| {
                sum(
                    o[k],
                    &[
                        (ax[0][k], rx.mul(cu)),
                        (ax[1][k], ry.mul(su)),
                        (ax[2][k], p(v)),
                    ],
                )
            })
        }
        Carrier::Sphere { frame, radius } => {
            let (ax, o) = (axes(frame), origin(frame));
            let (cv, sv) = (cos_iv(v, v), sin_iv(v, v));
            let (rc, rs) = (p(*radius).mul(cv), p(*radius).mul(sv));
            core::array::from_fn(|k| {
                let ring = sum(Iv::point(0.0), &[(ax[0][k], cu), (ax[1][k], su)]);
                sum(o[k], &[(ring, rc), (ax[2][k], rs)])
            })
        }
        Carrier::Torus(t) => {
            let (ax, o) = (axes(&t.frame), origin(&t.frame));
            let (cv, sv) = (cos_iv(v, v), sin_iv(v, v));
            let ring = p(t.major_radius).add(p(t.minor_radius).mul(cv));
            let rs = p(t.minor_radius).mul(sv);
            core::array::from_fn(|k| {
                let dir = sum(Iv::point(0.0), &[(ax[0][k], cu), (ax[1][k], su)]);
                sum(o[k], &[(dir, ring), (ax[2][k], rs)])
            })
        }
        Carrier::Spline(_) => return None,
    })
}
