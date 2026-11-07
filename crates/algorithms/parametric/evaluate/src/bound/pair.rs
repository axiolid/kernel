//! Chord bounds of the pcurves of a pair section (#414).
//!
//! Between two nodes a [`PairSection3`] is defined by four equations in the
//! parameters `x = (a, b)` on its two surfaces: `S1(a) = S2(b)` and the
//! plane across the chord `P0 -> P1` at `P0 + s d`, `d = P1 - P0`, `s` the
//! chord's local parameter (so `t = i + s`). With `n = d / |d|`,
//!
//! `G(s, x) = (S1(a) - S2(b), n . (S1(a) - P0) - s |d|) = 0`,
//!
//! whose matrix in `x` is `M = [[J1, -J2], [n^T J1, 0]]` and whose rate in
//! `s` is `(0, 0, 0, -|d|)`. This is the continuation of `lifted` on a
//! four-dimensional system.
//!
//! - **Separation.** In parameters scaled by `D = diag(|S1_u|, |S1_v|,
//!   |S2_u|, |S2_v|)` at the boxes' centres, `M' = M D^-1` has the smallest
//!   singular value `mu` over the product of the two boxes (a convex 4D
//!   box, so its hull argument gives injectivity of `x -> G(s, x)` for
//!   each `s`). Take `R`, the floating-point inverse of `M'` at the
//!   centres. If `eps = |I - R M'_c|_F < 1` then `sigma_min(M'_c) >= (1 -
//!   eps) / |R|_F`, and each surface's scaled columns move over its box by
//!   at most the spread of `separation`, so by Weyl's inequality `mu >= (1 -
//!   eps) / |R|_F - |dM'|_F`, `|dM'|_F^2 <= 2 d1^2 + d2^2` (the first
//!   surface's columns enter twice, in `J1` and in `n^T J1`) plus an
//!   allowance for the rounding of the matrix.
//! - **Speed and bend.** `M' y' = (0, 0, 0, |d|)` gives `|y'| <= |d| / mu`
//!   for `y = D x`; and `M x'' = -(Q1[a', a'] - Q2[b', b'], n . Q1[a', a'])`
//!   (`Q_k` the second derivative of surface `k`) gives `|y''| <= |y'|^2
//!   (sqrt 2 q1 + q2) / mu`, `q_k = S_k,uu / |S_k,u|^2 + 2 S_k,uv /
//!   (|S_k,u| |S_k,v|) + S_k,vv / |S_k,v|^2`. A side's pcurve is `y`
//!   unscaled through the smaller of its own two scale factors.
//! - **Reach and ends.** As in `lifted`, over the 4D box of the ends'
//!   hulls widened per surface: the lift of an exact end moves at most
//!   `|d| h / (2 mu)` to the middle of the piece and starts within
//!   `|G(s_i, x_i)| / mu` of the evaluated end, the residual enclosed in
//!   outward-rounded arithmetic (`SurfaceBoundOracle::point_enclosure`).
//!
//! The boxes are clamped to a B-spline's knot domain: the section lies in
//! it over the span, as `FaceBound::lens` assumes, and an evaluated end
//! outside it is refused.

use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_curve::PairSection3;

use super::continuation::{split_chord, PieceBound};
use super::interval::Iv;
use super::separation::BoxJet;
use super::SurfaceBoundOracle;

/// Growth steps of the boxes before the bound is refused.
const RHO_STEPS: usize = 64;

/// Allowance for the rounding of the scaled matrix `M'` (16 entries, each
/// within a few hundred ulps of its true value).
const MATRIX_PAD: Scalar = 16.0 * 256.0 * Scalar::EPSILON;

/// A certified bound on the distance from a pair section's pcurve on one of
/// its surfaces (`first` or the second) over `[lo, hi]`, within one chord,
/// to the chord joining its evaluated ends, or `None` where a piece cannot
/// be certified.
pub(super) fn chord_bound(
    section: &PairSection3,
    first: bool,
    lo: Scalar,
    hi: Scalar,
) -> Option<Scalar> {
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return None;
    }
    let chord = Chord::new(section, first, lo, hi)?;
    split_chord(lo, hi, |a, b| chord.piece(a, b))
}

/// A chord of the section, with the oracles of its two surfaces.
struct Chord<'a> {
    section: &'a PairSection3,
    /// Whether the pcurve is the first surface's.
    first: bool,
    index: usize,
    origin: Point3,
    d: Vec3,
    /// `|d|`.
    length: Iv,
    oracles: [SurfaceBoundOracle<'a>; 2],
}

/// An evaluated end: the parameters on both surfaces and `|G|` there.
struct End {
    at: [Point2; 2],
    residual: Scalar,
}

impl<'a> Chord<'a> {
    /// The chord of the section holding `[lo, hi]`, or `None` over more than
    /// one chord, for a degenerate chord, or for non-finite nodes.
    fn new(section: &'a PairSection3, first: bool, lo: Scalar, hi: Scalar) -> Option<Self> {
        let n = section.nodes.len();
        if n < 2 {
            return None;
        }
        if section
            .nodes
            .iter()
            .any(|k| !(k.point.is_finite() && k.first.is_finite() && k.second.is_finite()))
        {
            return None;
        }
        let index = (lo.floor().max(0.0) as usize).min(n - 2);
        if lo < index as Scalar || hi > (index + 1) as Scalar {
            return None;
        }
        let (p0, p1) = (section.nodes[index].point, section.nodes[index + 1].point);
        let d = p1 - p0;
        let squares = [d.x, d.y, d.z]
            .map(|c| Iv::point(c).sq())
            .into_iter()
            .fold(Iv::point(0.0), Iv::add);
        let length = squares.sqrt();
        if !(length.finite() && length.lo > 0.0) {
            return None;
        }
        Some(Self {
            section,
            first,
            index,
            origin: p0,
            d,
            length,
            oracles: [
                SurfaceBoundOracle::for_carrier(&section.first)?,
                SurfaceBoundOracle::for_carrier(&section.second)?,
            ],
        })
    }
}

impl Chord<'_> {
    fn end(&self, t: Scalar) -> Option<End> {
        let (a, b, _) = self.section.solve(t)?;
        let s = t - self.index as Scalar;
        for (oracle, p) in self.oracles.iter().zip([a, b]) {
            if let Some((u, v)) = oracle.domain() {
                if p.x < u.0 || p.x > u.1 || p.y < v.0 || p.y > v.1 {
                    return None;
                }
            }
        }
        let (s1, s2) = (
            self.oracles[0].point_enclosure(a.x, a.y)?,
            self.oracles[1].point_enclosure(b.x, b.y)?,
        );
        let origin = [self.origin.x, self.origin.y, self.origin.z];
        let d = [self.d.x, self.d.y, self.d.z];
        let mut squares = Iv::point(0.0);
        let mut plane = Iv::point(s).mul(self.length).neg();
        for k in 0..3 {
            squares = squares.add(s1[k].add(s2[k].neg()).sq());
            let n = Iv::point(d[k]).div(self.length);
            plane = plane.add(n.mul(s1[k].add(Iv::point(-origin[k]))));
        }
        let residual = squares.add(plane.sq()).sqrt().hi;
        residual.is_finite().then_some(End {
            at: [a, b],
            residual,
        })
    }

    /// One piece `[a, b]`; see the module docs.
    fn piece(&self, a: Scalar, b: Scalar) -> Option<PieceBound> {
        let (start, finish) = (self.end(a)?, self.end(b)?);
        let h = Iv::point(b).add(Iv::point(-a)).hi;
        let residual = start.residual.max(finish.residual);
        let side = if self.first { 0 } else { 1 };
        let hull = |x: Scalar, y: Scalar| (x.min(y), x.max(y));
        let hulls = [0, 1].map(|k| {
            (
                hull(start.at[k].x, finish.at[k].x),
                hull(start.at[k].y, finish.at[k].y),
            )
        });
        let mut rho = [0.0_f64; 2];
        for _ in 0..RHO_STEPS {
            let boxes = [0, 1].map(|k| self.grown(k, hulls[k], rho[k]));
            let jets = [
                self.oracles[0].box_jet(boxes[0].0, boxes[0].1)?,
                self.oracles[1].box_jet(boxes[1].0, boxes[1].1)?,
            ];
            let mu = self.separation(&jets)?;
            let m = jets.map(|j| j.scale[0].min(j.scale[1]));
            // `|y'| <= |d| / mu`, and the end error `|G| / mu`, in scaled
            // parameters.
            let speed = self.length.div(Iv::point(mu)).hi;
            let error = Iv::point(residual).div(Iv::point(mu)).hi;
            let reach = Iv::point(speed)
                .mul(Iv::point(h))
                .scale(0.5)
                .add(Iv::point(error))
                .hi;
            let reaches = m.map(|m| Iv::point(reach).div(Iv::point(m)).hi);
            if reaches[0] < rho[0] && reaches[1] < rho[1] {
                let floor = Iv::point(mu).mul(Iv::point(m[side]));
                let bend = bend(&jets, speed, mu)?;
                return Some(PieceBound {
                    slope: Iv::point(speed).div(Iv::point(m[side])).hi,
                    bend: Some(Iv::point(bend).div(Iv::point(m[side])).hi),
                    start: Iv::point(start.residual).div(floor).hi,
                    end: Iv::point(finish.residual).div(floor).hi,
                });
            }
            for k in 0..2 {
                rho[k] = (1.1 * reaches[k])
                    .max(1.1 * rho[k])
                    .max(Scalar::MIN_POSITIVE);
            }
        }
        None
    }

    /// The box about the hull of one surface's parameters, widened by `rho`
    /// and clamped to a spline's knot domain.
    fn grown(
        &self,
        k: usize,
        hull: ((Scalar, Scalar), (Scalar, Scalar)),
        rho: Scalar,
    ) -> ((Scalar, Scalar), (Scalar, Scalar)) {
        let widen = |(lo, hi): (Scalar, Scalar)| {
            (
                Iv::point(lo).add(Iv::point(-rho)).lo,
                Iv::point(hi).add(Iv::point(rho)).hi,
            )
        };
        let (u, v) = (widen(hull.0), widen(hull.1));
        match self.oracles[k].domain() {
            Some((du, dv)) => (
                (u.0.clamp(du.0, du.1), u.1.clamp(du.0, du.1)),
                (v.0.clamp(dv.0, dv.1), v.1.clamp(dv.0, dv.1)),
            ),
            None => (u, v),
        }
    }

    /// The certified smallest singular value `mu` of the scaled matrix over
    /// the two boxes, or `None` where it cannot be shown positive.
    fn separation(&self, jets: &[BoxJet; 2]) -> Option<Scalar> {
        let n = self.d / self.length.hi;
        let [j1, j2] = jets;
        let scaled = |j: &BoxJet, k: usize| j.columns[k] / j.scale[k];
        let columns = [scaled(j1, 0), scaled(j1, 1), -scaled(j2, 0), -scaled(j2, 1)];
        let mut m = [[0.0; 4]; 4];
        for (j, c) in columns.iter().enumerate() {
            m[0][j] = c.x;
            m[1][j] = c.y;
            m[2][j] = c.z;
            m[3][j] = if j < 2 { n.dot(*c) } else { 0.0 };
        }
        let r = inverse(m)?;
        // `E = I - R M`, in enclosures.
        let mut e_squares = Iv::point(0.0);
        let mut r_squares = Iv::point(0.0);
        for (i, row) in r.iter().enumerate() {
            for j in 0..4 {
                let product = row
                    .iter()
                    .zip(&m)
                    .map(|(&rik, mk)| Iv::point(rik).mul(Iv::point(mk[j])))
                    .fold(Iv::point(0.0), Iv::add);
                let target = if i == j { 1.0 } else { 0.0 };
                let entry = Iv::point(target).add(product.neg());
                e_squares = e_squares.add(entry.sq());
                r_squares = r_squares.add(Iv::point(row[j]).sq());
            }
        }
        let (eps, norm) = (e_squares.sqrt().hi, r_squares.sqrt().hi);
        if !(eps < 1.0 && norm.is_finite() && norm > 0.0) {
            return None;
        }
        let delta = Iv::point(j1.spread)
            .sq()
            .scale(2.0)
            .add(Iv::point(j2.spread).sq())
            .sqrt()
            .add(Iv::point(MATRIX_PAD))
            .hi;
        let centre = Iv::point(Iv::point(1.0).add(Iv::point(-eps)).lo)
            .div(Iv::point(norm))
            .lo;
        let mu = Iv::point(centre).add(Iv::point(-delta)).lo;
        (mu.is_finite() && mu > 0.0).then_some(mu)
    }
}

/// `sup |y''|` from `|y''| <= |y'|^2 (sqrt 2 q1 + q2) / mu`.
fn bend(jets: &[BoxJet; 2], speed: Scalar, mu: Scalar) -> Option<Scalar> {
    let q = jets.map(|j| {
        let (a, b) = (Iv::point(j.scale[0]), Iv::point(j.scale[1]));
        Iv::point(j.bounds.duu)
            .div(a.mul(a))
            .add(Iv::point(2.0 * j.bounds.duv).div(a.mul(b)))
            .add(Iv::point(j.bounds.dvv).div(b.mul(b)))
    });
    let bend = Iv::point(core::f64::consts::SQRT_2.next_up())
        .mul(q[0])
        .add(q[1])
        .mul(Iv::point(speed).sq())
        .div(Iv::point(mu))
        .hi;
    bend.is_finite().then_some(bend)
}

/// The inverse of a 4x4 matrix by Gauss-Jordan elimination with partial
/// pivoting, or `None` if it is singular. Rounding is the caller's to
/// certify: the result is only an approximation of the inverse.
fn inverse(mut m: [[Scalar; 4]; 4]) -> Option<[[Scalar; 4]; 4]> {
    let mut r = [[0.0; 4]; 4];
    for (i, row) in r.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for col in 0..4 {
        let pivot = (col..4).max_by(|&a, &b| m[a][col].abs().total_cmp(&m[b][col].abs()))?;
        if !(m[pivot][col].is_finite() && m[pivot][col] != 0.0) {
            return None;
        }
        m.swap(col, pivot);
        r.swap(col, pivot);
        let p = m[col][col];
        for k in 0..4 {
            m[col][k] /= p;
            r[col][k] /= p;
        }
        for row in 0..4 {
            if row != col {
                let f = m[row][col];
                for k in 0..4 {
                    m[row][k] -= f * m[col][k];
                    r[row][k] -= f * r[col][k];
                }
            }
        }
    }
    r.iter().flatten().all(|x| x.is_finite()).then_some(r)
}

#[cfg(test)]
mod tests;
