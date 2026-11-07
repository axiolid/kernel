//! Enclosures of a surface's point at given parameters (#414).
//!
//! The certified chord bounds of lifted pcurves and pair sections measure
//! how far an evaluated parameter pair is from the curve it stands for by
//! the residual of its image: the point the surface gives there against the
//! space curve's. The point is enclosed in outward-rounded arithmetic, so
//! the residual is certified rather than estimated: the closed forms of the
//! analytic carriers (`carrier`), and for a B-spline surface the rational
//! tensor sum of its basis functions, each by the recurrence of the NURBS
//! book (algorithm A2.2) over intervals, in the same span as the evaluator
//! and clamped to the knot domain as it is.

use axiolid_core::Scalar;
use axiolid_curve::Carrier;
use axiolid_surface::Surface;

use super::interval::Iv;
use super::{carrier, Source, SplineNets, SurfaceBoundOracle};
use crate::curve::span_in;

impl SurfaceBoundOracle<'_> {
    /// The enclosure of the surface's point at `(u, v)`, or `None` for a
    /// surface family without one.
    pub(super) fn point_enclosure(&self, u: Scalar, v: Scalar) -> Option<[Iv; 3]> {
        match self.source {
            Source::Carrier(Carrier::Spline(_)) | Source::Surface(Surface::BSpline(_)) => {
                self.spline.as_ref()?.point_enclosure(u, v)
            }
            Source::Carrier(c) => carrier::point_enclosure(c, u, v),
            Source::Surface(_) => None,
        }
    }
}

/// The `degree + 1` non-zero basis functions at `u` in `span`, enclosed.
fn basis(knots: &[Scalar], degree: usize, span: usize, u: Scalar) -> Vec<Iv> {
    let p = Iv::point;
    let mut n = vec![p(1.0)];
    let mut left = vec![p(0.0); degree + 1];
    let mut right = vec![p(0.0); degree + 1];
    for j in 1..=degree {
        left[j] = p(u).add(p(-knots[span + 1 - j]));
        right[j] = p(knots[span + j]).add(p(-u));
        let mut saved = p(0.0);
        for r in 0..j {
            let temp = n[r].div(right[r + 1].add(left[j - r]));
            n[r] = saved.add(right[r + 1].mul(temp));
            saved = left[j - r].mul(temp);
        }
        n.push(saved);
    }
    n
}

impl SplineNets {
    /// The enclosure of the B-spline surface's point at `(u, v)`; `None`
    /// where the knot vectors or weights do not give a positive
    /// denominator.
    fn point_enclosure(&self, u: Scalar, v: Scalar) -> Option<[Iv; 3]> {
        let clamp = |x: Scalar, knots: &[Scalar], d: usize, n: usize| x.clamp(knots[d], knots[n]);
        let (u, v) = (
            clamp(u, &self.u_knots, self.p, self.rows),
            clamp(v, &self.v_knots, self.q, self.cols),
        );
        let su = span_in(&self.u_knots, self.rows, self.p, u);
        let sv = span_in(&self.v_knots, self.cols, self.q, v);
        let (nu, nv) = (
            basis(&self.u_knots, self.p, su, u),
            basis(&self.v_knots, self.q, sv, v),
        );
        let mut numerator = [Iv::point(0.0); 3];
        let mut denominator = Iv::point(0.0);
        for (i, bu) in nu.iter().enumerate() {
            for (j, bv) in nv.iter().enumerate() {
                let (row, col) = (su - self.p + i, sv - self.q + j);
                let weighted = bu.mul(*bv).scale(self.weights[row][col]);
                denominator = denominator.add(weighted);
                for (k, n) in numerator.iter_mut().enumerate() {
                    *n = n.add(weighted.scale(self.points[row][col][k]));
                }
            }
        }
        if !(denominator.finite() && denominator.lo > 0.0) {
            return None;
        }
        let point = numerator.map(|n| n.div(denominator));
        point.iter().all(|x| x.finite()).then_some(point)
    }
}
