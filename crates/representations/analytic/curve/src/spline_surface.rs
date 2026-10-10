//! Polynomial and rational B-spline surfaces.
//!
//! The data type lives here, next to the B-spline curve, so that a curve
//! lying on a B-spline surface (an [`crate::ImplicitSection3`] traced in the
//! surface's parameters, ADR 0077) can carry its carrier.
//! `axiolid_surface` re-exports it as `axiolid_surface::BSplineSurface`,
//! unchanged.

use axiolid_core::{Point3, Scalar, Vec3};

use crate::implicit::SurfaceJet;
use crate::spline::KnotSpec;

/// Tensor-product B-spline surface preserving exact knot and weight data.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct BSplineSurface {
    /// Degree along the first parameter axis.
    pub u_degree: u16,
    /// Degree along the second parameter axis.
    pub v_degree: u16,
    /// Rectangular control net, row-major in `u` then `v`.
    pub control_points: Vec<Vec<Point3>>,
    /// Distinct knots along `u`.
    pub u_knots: Vec<Scalar>,
    /// Multiplicities matching `u_knots`.
    pub u_multiplicities: Vec<u32>,
    /// Distinct knots along `v`.
    pub v_knots: Vec<Scalar>,
    /// Multiplicities matching `v_knots`.
    pub v_multiplicities: Vec<u32>,
    /// Optional rational weight net matching the control net shape.
    pub weights: Option<Vec<Vec<Scalar>>>,
    /// Whether the surface closes along `u`.
    pub u_closed: bool,
    /// Whether the surface closes along `v`.
    pub v_closed: bool,
    /// Source knot convention.
    pub knot_spec: KnotSpec,
    /// Whether the source declares self intersection.
    pub self_intersect: Option<bool>,
}

/// The full knot vector from distinct knots and their multiplicities.
fn expand(knots: &[Scalar], multiplicities: &[u32]) -> Vec<Scalar> {
    let mut out = Vec::new();
    for (&k, &m) in knots.iter().zip(multiplicities) {
        out.extend(core::iter::repeat_n(k, m as usize));
    }
    out
}

/// The span index `i` with `knots[i] <= t < knots[i + 1]`, clamped to the
/// valid range `[degree, count - 1]`.
fn span(knots: &[Scalar], degree: usize, count: usize, t: Scalar) -> usize {
    let mut s = degree;
    for (k, knot) in knots.iter().enumerate().take(count).skip(degree) {
        if *knot <= t {
            s = k;
        } else {
            break;
        }
    }
    s
}

/// Basis functions and their first two derivatives at `t` in span `s`
/// (Piegl and Tiller, A2.3): `out[k][j]` is the `k`-th derivative of
/// `N_{s - p + j}`.
#[allow(clippy::needless_range_loop)] // Piegl and Tiller's indices, kept
fn basis(knots: &[Scalar], p: usize, s: usize, t: Scalar) -> [Vec<Scalar>; 3] {
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots[s + 1 - j];
        right[j] = knots[s + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = if ndu[j][r] == 0.0 {
                0.0
            } else {
                ndu[r][j - 1] / ndu[j][r]
            };
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let orders = 2.min(p);
    let mut ders = [vec![0.0; p + 1], vec![0.0; p + 1], vec![0.0; p + 1]];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    let mut a = [vec![0.0; p + 1], vec![0.0; p + 1]];
    for r in 0..=p {
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=orders {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p as isize - k as isize;
            if r >= k {
                let denom = ndu[(pk + 1) as usize][rk as usize];
                a[s2][0] = if denom == 0.0 { 0.0 } else { a[s1][0] / denom };
                d = a[s2][0] * ndu[rk as usize][pk as usize];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk { k - 1 } else { p - r };
            for j in j1..=j2 {
                let denom = ndu[(pk + 1) as usize][(rk + j as isize) as usize];
                a[s2][j] = if denom == 0.0 {
                    0.0
                } else {
                    (a[s1][j] - a[s1][j - 1]) / denom
                };
                d += a[s2][j] * ndu[(rk + j as isize) as usize][pk as usize];
            }
            if r as isize <= pk {
                let denom = ndu[(pk + 1) as usize][r];
                a[s2][k] = if denom == 0.0 {
                    0.0
                } else {
                    -a[s1][k - 1] / denom
                };
                d += a[s2][k] * ndu[r][pk as usize];
            }
            ders[k][r] = d;
            core::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as Scalar;
    for k in 1..=orders {
        for value in &mut ders[k] {
            *value *= factor;
        }
        factor *= (p - k) as Scalar;
    }
    ders
}

impl BSplineSurface {
    /// The parameter domain `((u0, u1), (v0, v1))`, or `None` for a net or
    /// knot vector that does not describe a surface.
    #[must_use]
    pub fn domain(&self) -> Option<((Scalar, Scalar), (Scalar, Scalar))> {
        let (p, q) = (usize::from(self.u_degree), usize::from(self.v_degree));
        let rows = self.control_points.len();
        let cols = self.control_points.first()?.len();
        let (ku, kv) = (
            expand(&self.u_knots, &self.u_multiplicities),
            expand(&self.v_knots, &self.v_multiplicities),
        );
        if p == 0 || q == 0 || ku.len() != rows + p + 1 || kv.len() != cols + q + 1 {
            return None;
        }
        Some(((ku[p], ku[rows]), (kv[q], kv[cols])))
    }

    /// The point at `(u, v)` with its first and second partials, clamped
    /// to the domain; `None` for a malformed net or a zero weight.
    #[must_use]
    #[allow(clippy::needless_range_loop)]
    pub fn jet(&self, u: Scalar, v: Scalar) -> Option<SurfaceJet> {
        let (p, q) = (usize::from(self.u_degree), usize::from(self.v_degree));
        let ((u0, u1), (v0, v1)) = self.domain()?;
        let (rows, cols) = (self.control_points.len(), self.control_points[0].len());
        let (ku, kv) = (
            expand(&self.u_knots, &self.u_multiplicities),
            expand(&self.v_knots, &self.v_multiplicities),
        );
        let (u, v) = (u.clamp(u0, u1), v.clamp(v0, v1));
        let (su, sv) = (span(&ku, p, rows, u), span(&kv, q, cols, v));
        let bu = basis(&ku, p, su, u);
        let bv = basis(&kv, q, sv, v);
        // Homogeneous sums and their partials: A (point times weight), W.
        let mut a = [[Vec3::ZERO; 3]; 3];
        let mut w = [[0.0; 3]; 3];
        for i in 0..=p {
            for j in 0..=q {
                let (row, col) = (su - p + i, sv - q + j);
                let weight = self.weights.as_ref().map_or(1.0, |net| net[row][col]);
                let point = self.control_points[row][col] * weight;
                for k in 0..3 {
                    for l in 0..3 - k {
                        let n = bu[k][i] * bv[l][j];
                        a[k][l] += point * n;
                        w[k][l] += weight * n;
                    }
                }
            }
        }
        if w[0][0] == 0.0 || !w[0][0].is_finite() {
            return None;
        }
        let inv = 1.0 / w[0][0];
        let s = a[0][0] * inv;
        let s_u = (a[1][0] - s * w[1][0]) * inv;
        let s_v = (a[0][1] - s * w[0][1]) * inv;
        let s_uu = (a[2][0] - s_u * (2.0 * w[1][0]) - s * w[2][0]) * inv;
        let s_vv = (a[0][2] - s_v * (2.0 * w[0][1]) - s * w[0][2]) * inv;
        let s_uv = (a[1][1] - s_u * w[0][1] - s_v * w[1][0] - s * w[1][1]) * inv;
        let jet = SurfaceJet {
            point: s,
            u: s_u,
            v: s_v,
            uu: s_uu,
            uv: s_uv,
            vv: s_vv,
        };
        (jet.point.is_finite() && jet.u.is_finite() && jet.v.is_finite()).then_some(jet)
    }
}
