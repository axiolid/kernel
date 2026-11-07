//! The separation of a surface's parameters over a box (#414).
//!
//! A parameter box maps injectively, and by a certified rate, when the
//! surface's Jacobian `J = [S_u S_v]` is bounded below over it. Take the
//! box's centre `m` and the scaling `D = diag(|S_u(m)|, |S_v(m)|)`, so that
//! `J(m) D^-1` has unit columns, the smallest singular value of which is
//! `sqrt(1 - |cos theta|)` (`theta` the angle between the columns; the
//! entries of the 2x2 Gram matrix, bounded by Gershgorin). Over the box the
//! columns move by at most the second partials' suprema times the box's
//! half widths, and Weyl's inequality carries the singular value across: a
//! matrix within `delta` (in the spectral norm) of one with smallest
//! singular value `s` has smallest singular value at least `s - delta`,
//! and `delta` is at most the Frobenius norm of the scaled column spread.
//! Every matrix of the convex hull of `J` over the box is within `delta`
//! too, which is what the injectivity needs: for `x`, `y` in the box,
//! `S(x) - S(y)` is the mean of `J` along the segment applied to `x - y`.
//!
//! The computed jet at `m` is padded by [`JET_PAD`] per scaled column. A
//! box where the separation cannot be shown positive (a pole or apex, a
//! fold, a box too wide, a crease of a B-spline) has none.

use axiolid_core::{Point3, Scalar, Vec3};

use super::interval::Iv;
use super::{ordered, SurfaceBoundOracle};

/// Relative rounding allowance of the computed jet at a box's centre, per
/// scaled column.
const JET_PAD: Scalar = 256.0 * Scalar::EPSILON;

/// A certified lower bound of a surface's Jacobian over a parameter box.
///
/// With `D = diag(scale)`, `|S(x) - S(y)| >= sigma |D (x - y)|` for `x`,
/// `y` in the box, and `|J(x) w| >= sigma |D w|` at every point of it.
/// `scale` is `|S_u|, |S_v|` at the box's centre as computed; it depends
/// only on the centre.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParameterSeparation {
    /// The scaling `|S_u|, |S_v|` at the box's centre.
    pub scale: [Scalar; 2],
    /// The separation, positive.
    pub sigma: Scalar,
}

impl SurfaceBoundOracle<'_> {
    /// The point and the first partials at `(u, v)`, or `None` where the
    /// surface does not evaluate to finite numbers there.
    pub(crate) fn jet(&self, u: Scalar, v: Scalar) -> Option<(Point3, Vec3, Vec3)> {
        use super::Source;
        let (point, du, dv) = match self.source {
            Source::Carrier(c) => {
                let j = c.jet(u, v);
                (j.point, j.u, j.v)
            }
            Source::Surface(s) => {
                let j = crate::surface::jet(s, u, v).ok()?;
                (j.point, j.du, j.dv)
            }
        };
        let finite = point.is_finite() && du.is_finite() && dv.is_finite();
        finite.then_some((point, du, dv))
    }

    /// The certified separation over the box `[u.0, u.1] x [v.0, v.1]`, or
    /// `None` where it cannot be shown positive: a pole, an apex, a fold, a
    /// box too wide, or a B-spline crease or edge inside it (`!smooth`).
    #[must_use]
    pub fn separation(
        &self,
        u: (Scalar, Scalar),
        v: (Scalar, Scalar),
    ) -> Option<ParameterSeparation> {
        let (u0, u1) = ordered(u.0, u.1);
        let (v0, v1) = ordered(v.0, v.1);
        let bounds = self.bounds((u0, u1), (v0, v1))?;
        if !bounds.smooth {
            return None;
        }
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let (_, su, sv) = self.jet(um, vm)?;
        let comp = |a: Vec3, b: Vec3| {
            Iv::point(a.x)
                .mul(Iv::point(b.x))
                .add(Iv::point(a.y).mul(Iv::point(b.y)))
                .add(Iv::point(a.z).mul(Iv::point(b.z)))
        };
        let (a, b) = (su.length(), sv.length());
        if !(a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0) {
            return None;
        }
        // The Gram matrix of the scaled columns.
        let (sa, sb) = (Iv::point(a), Iv::point(b));
        let (guu, gvv) = (comp(su, su).div(sa.mul(sa)), comp(sv, sv).div(sb.mul(sb)));
        let guv = comp(su, sv).div(sa.mul(sb));
        if !(guu.finite() && gvv.finite() && guv.finite()) {
            return None;
        }
        let lambda = Iv::point(guu.lo.min(gvv.lo)).add(Iv::point(-guv.mag()));
        let centre = lambda.sqrt().lo;
        // The columns' spread over the box, in scaled parameters.
        let half = |lo: Scalar, hi: Scalar, m: Scalar| {
            Iv::point(hi)
                .add(Iv::point(-m))
                .hi
                .max(Iv::point(m).add(Iv::point(-lo)).hi)
        };
        let (ru, rv) = (Iv::point(half(u0, u1, um)), Iv::point(half(v0, v1, vm)));
        let column = |x: Scalar, y: Scalar, scale: Scalar| {
            Iv::point(x)
                .mul(ru)
                .add(Iv::point(y).mul(rv))
                .div(Iv::point(scale))
        };
        let (cu, cv) = (
            column(bounds.duu, bounds.duv, a),
            column(bounds.duv, bounds.dvv, b),
        );
        let spread = cu.sq().add(cv.sq()).sqrt().hi;
        let sigma = Iv::point(centre)
            .add(Iv::point(-spread))
            .add(Iv::point(-2.0 * JET_PAD))
            .lo;
        (sigma.is_finite() && sigma > 0.0).then_some(ParameterSeparation {
            scale: [a, b],
            sigma,
        })
    }
}
