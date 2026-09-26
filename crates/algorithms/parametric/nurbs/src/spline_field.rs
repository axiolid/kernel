//! An analytic surface's equation read on a B-spline surface, as
//! piecewise Bernstein polynomials (ADR 0077).
//!
//! The B-spline surface is split into its rational Bezier patches (knot
//! insertion to full multiplicity). On each patch the homogeneous point
//! `(X, W)` is a Bernstein polynomial of degree `(p, q)` in the patch's local
//! coordinates, and the other surface's implicit equation, homogenised to
//! degree `d` in `(X, W)`, is then a Bernstein polynomial of degree
//! `(d p, d q)`: products of Bernstein polynomials stay Bernstein. Its zero
//! set is the section wherever the weight is positive, which a B-spline's
//! weights always are. The coefficients bound the field on each patch, so
//! the trace over it is certified exactly as over a series field.

use axiolid_core::{Scalar, Vec3};
use axiolid_curve::{BSplineSurface, Curve3, Field2, PatchField2};
use axiolid_surface::Surface;

/// A homogeneous control point `(x w, y w, z w, w)`.
type H = [Scalar; 4];

/// The full knot vector from distinct knots and multiplicities.
fn expand(knots: &[Scalar], multiplicities: &[u32]) -> Vec<Scalar> {
    let mut out = Vec::new();
    for (&k, &m) in knots.iter().zip(multiplicities) {
        out.extend(core::iter::repeat_n(k, m as usize));
    }
    out
}

/// A clamped B-spline's Bezier segments (Piegl and Tiller, A5.6): the
/// breaks between them and each segment's `p + 1` control points. `None`
/// for an unclamped knot vector.
fn decompose(knots: &[Scalar], p: usize, control: &[H]) -> Option<(Vec<Scalar>, Vec<Vec<H>>)> {
    let m = knots.len() - 1;
    let n = control.len();
    if knots.len() != n + p + 1 {
        return None;
    }
    // Clamped at both ends.
    if knots[..=p].iter().any(|k| *k != knots[0]) || knots[m - p..].iter().any(|k| *k != knots[m]) {
        return None;
    }
    let mut segments: Vec<Vec<H>> = vec![control[..=p].to_vec()];
    let mut breaks = vec![knots[p]];
    let (mut a, mut b) = (p, p + 1);
    let mut alphas = vec![0.0; p + 1];
    while b < m {
        let i = b;
        while b < m && knots[b + 1] == knots[b] {
            b += 1;
        }
        let mult = b - i + 1;
        let last = segments.len() - 1;
        if mult < p {
            let numer = knots[b] - knots[a];
            for j in (mult + 1..=p).rev() {
                alphas[j - mult - 1] = numer / (knots[a + j] - knots[a]);
            }
            let r = p - mult;
            let mut next = vec![[0.0; 4]; p + 1];
            for j in 1..=r {
                let save = r - j;
                let s = mult + j;
                for k in (s..=p).rev() {
                    let alpha = alphas[k - s];
                    let (hi, lo) = (segments[last][k], segments[last][k - 1]);
                    for c in 0..4 {
                        segments[last][k][c] = alpha * hi[c] + (1.0 - alpha) * lo[c];
                    }
                }
                if b < m {
                    next[save] = segments[last][p];
                }
            }
            if b < m {
                for idx in (p - mult)..=p {
                    next[idx] = control[b - p + idx];
                }
                breaks.push(knots[b]);
                segments.push(next);
            }
        } else if b < m {
            let mut next = vec![[0.0; 4]; p + 1];
            for idx in (p - mult)..=p {
                next[idx] = control[b - p + idx];
            }
            breaks.push(knots[b]);
            segments.push(next);
        }
        if b < m {
            a = b;
            b += 1;
        }
    }
    breaks.push(knots[m - p]);
    // Degenerate spans (repeated interior breaks) carry no segment.
    Some((breaks, segments))
}

/// A tensor-product Bernstein polynomial.
#[derive(Clone, Debug)]
struct Bern {
    p: usize,
    q: usize,
    c: Vec<Scalar>,
}

fn binomial(n: usize, k: usize) -> Scalar {
    let mut r = 1.0;
    for i in 0..k {
        r = r * (n - i) as Scalar / (i + 1) as Scalar;
    }
    r
}

impl Bern {
    fn at(&self, a: usize, b: usize) -> Scalar {
        self.c[a * (self.q + 1) + b]
    }

    fn mul(&self, o: &Bern) -> Bern {
        let (p, q) = (self.p + o.p, self.q + o.q);
        let mut c = vec![0.0; (p + 1) * (q + 1)];
        for a1 in 0..=self.p {
            for b1 in 0..=self.q {
                let x = self.at(a1, b1) * binomial(self.p, a1) * binomial(self.q, b1);
                if x == 0.0 {
                    continue;
                }
                for a2 in 0..=o.p {
                    for b2 in 0..=o.q {
                        let y = o.at(a2, b2) * binomial(o.p, a2) * binomial(o.q, b2);
                        c[(a1 + a2) * (q + 1) + b1 + b2] += x * y;
                    }
                }
            }
        }
        for a in 0..=p {
            for b in 0..=q {
                c[a * (q + 1) + b] /= binomial(p, a) * binomial(q, b);
            }
        }
        Bern { p, q, c }
    }

    fn add(&self, o: &Bern, scale: Scalar) -> Bern {
        debug_assert_eq!((self.p, self.q), (o.p, o.q));
        Bern {
            p: self.p,
            q: self.q,
            c: self
                .c
                .iter()
                .zip(&o.c)
                .map(|(x, y)| x + scale * y)
                .collect(),
        }
    }

    fn scaled(&self, s: Scalar) -> Bern {
        Bern {
            p: self.p,
            q: self.q,
            c: self.c.iter().map(|x| x * s).collect(),
        }
    }
}

/// The other surface's equation homogenised in `(X, W)`: every term of
/// degree `d`, so products and sums stay in one Bernstein degree.
fn equation(point: &[Bern; 4], other: &Surface) -> Option<Bern> {
    let w = &point[3];
    // Coordinate along `axis` from `origin`, homogeneous: `axis . X - (axis . O) W`.
    let local = |origin: Vec3, axis: Vec3| -> Bern {
        let axis = axis.normalize();
        let mut out = w.scaled(-axis.dot(origin));
        for k in 0..3 {
            out = out.add(&point[k], axis[k]);
        }
        out
    };
    let frame = |f: &axiolid_core::Frame3| {
        [
            local(f.origin, f.x),
            local(f.origin, f.y),
            local(f.origin, f.z),
        ]
    };
    let ww = w.mul(w);
    Some(match other {
        Surface::Plane(p) => local(p.frame.origin, p.frame.z),
        Surface::Cylinder(c) => {
            let [x, y, _] = frame(&c.frame);
            x.mul(&x)
                .add(&y.mul(&y), 1.0)
                .add(&ww, -c.radius * c.radius)
        }
        Surface::EllipticalCylinder(c) => {
            let [x, y, _] = frame(&c.frame);
            let (a2, b2) = (c.semi_axis_x.powi(2), c.semi_axis_y.powi(2));
            x.mul(&x).scaled(b2).add(&y.mul(&y), a2).add(&ww, -a2 * b2)
        }
        Surface::Cone(c) => {
            let [x, y, z] = frame(&c.frame);
            let radius = w.scaled(c.radius).add(&z, c.semi_angle.tan());
            x.mul(&x)
                .add(&y.mul(&y), 1.0)
                .add(&radius.mul(&radius), -1.0)
        }
        Surface::Sphere(s) => {
            let [x, y, z] = frame(&s.frame);
            x.mul(&x)
                .add(&y.mul(&y), 1.0)
                .add(&z.mul(&z), 1.0)
                .add(&ww, -s.radius * s.radius)
        }
        Surface::Torus(t) => {
            let [x, y, z] = frame(&t.frame);
            let (big, small) = (t.major_radius, t.minor_radius);
            let planar = x.mul(&x).add(&y.mul(&y), 1.0);
            let inner = planar
                .add(&z.mul(&z), 1.0)
                .add(&ww, big * big - small * small);
            inner.mul(&inner).add(&ww.mul(&planar), -4.0 * big * big)
        }
        _ => return None,
    })
}

/// Degrees, breaks and homogeneous Bezier control nets (`[iu][jv][a][b]`)
/// of a clamped B-spline surface.
#[allow(clippy::type_complexity)]
fn bezier_net(
    b: &BSplineSurface,
) -> Option<(
    usize,
    usize,
    Vec<Scalar>,
    Vec<Scalar>,
    Vec<Vec<Vec<Vec<H>>>>,
)> {
    let (p, q) = (usize::from(b.u_degree), usize::from(b.v_degree));
    let rows = b.control_points.len();
    let cols = b.control_points.first()?.len();
    let ku = expand(&b.u_knots, &b.u_multiplicities);
    let kv = expand(&b.v_knots, &b.v_multiplicities);
    let homogeneous = |i: usize, j: usize| -> H {
        let w = b.weights.as_ref().map_or(1.0, |net| net[i][j]);
        let pt = b.control_points[i][j];
        [pt.x * w, pt.y * w, pt.z * w, w]
    };
    // Along v, row by row.
    let mut rows_split: Vec<Vec<Vec<H>>> = Vec::with_capacity(rows);
    let mut v_breaks = Vec::new();
    for i in 0..rows {
        let row: Vec<H> = (0..cols).map(|j| homogeneous(i, j)).collect();
        let (breaks, segments) = decompose(&kv, q, &row)?;
        v_breaks = breaks;
        rows_split.push(segments);
    }
    let nv = rows_split[0].len();
    // Along u, per v segment and local column.
    let mut u_breaks = Vec::new();
    // patches[iu][jv][a][b]
    let mut net: Vec<Vec<Vec<Vec<H>>>> = Vec::new();
    for jv in 0..nv {
        for bq in 0..=q {
            let column: Vec<H> = (0..rows).map(|i| rows_split[i][jv][bq]).collect();
            let (breaks, segments) = decompose(&ku, p, &column)?;
            u_breaks = breaks;
            if net.is_empty() {
                net = vec![vec![vec![vec![[0.0; 4]; q + 1]; p + 1]; nv]; segments.len()];
            }
            for (iu, seg) in segments.iter().enumerate() {
                for (a, h) in seg.iter().enumerate() {
                    net[iu][jv][a][bq] = *h;
                }
            }
        }
    }
    let nu = net.len();
    if u_breaks.len() != nu + 1 || v_breaks.len() != nv + 1 {
        return None;
    }
    Some((p, q, u_breaks, v_breaks, net))
}

/// `other`'s equation on the B-spline surface `b`, cell by Bezier cell;
/// `None` for a B-spline `other`, or an unclamped or malformed `b`.
#[allow(clippy::needless_range_loop)] // patch grid by index
pub(crate) fn spline_section_field(b: &BSplineSurface, other: &Surface) -> Option<Field2> {
    let (p, q, u_breaks, v_breaks, net) = bezier_net(b)?;
    let (nu, nv) = (u_breaks.len() - 1, v_breaks.len() - 1);
    let mut patches = Vec::with_capacity(nu * nv);
    let (mut du, mut dv) = (0, 0);
    for iu in 0..nu {
        for jv in 0..nv {
            let comp = |k: usize| Bern {
                p,
                q,
                c: (0..=p)
                    .flat_map(|a| (0..=q).map(move |bb| (a, bb)))
                    .map(|(a, bb)| net[iu][jv][a][bb][k])
                    .collect(),
            };
            let point = [comp(0), comp(1), comp(2), comp(3)];
            let f = equation(&point, other)?;
            du = f.p;
            dv = f.q;
            patches.push(f.c);
        }
    }
    Some(Field2::Patches(PatchField2 {
        u_breaks,
        v_breaks,
        u_degree: du,
        v_degree: dv,
        patches,
    }))
}

/// The four homogeneous coordinates `(x w, y w, z w, w)` of a B-spline
/// surface as patch fields over its Bezier cells, for bounding its image
/// over parameter boxes.
#[allow(clippy::needless_range_loop)]
pub(crate) fn homogeneous_fields(b: &BSplineSurface) -> Option<[Field2; 4]> {
    let (p, q, u_breaks, v_breaks, net) = bezier_net(b)?;
    let (nu, nv) = (u_breaks.len() - 1, v_breaks.len() - 1);
    let make = |k: usize| {
        let mut patches = Vec::with_capacity(nu * nv);
        for iu in 0..nu {
            for jv in 0..nv {
                let mut c = Vec::with_capacity((p + 1) * (q + 1));
                for a in 0..=p {
                    for bb in 0..=q {
                        c.push(net[iu][jv][a][bb][k]);
                    }
                }
                patches.push(c);
            }
        }
        Field2::Patches(PatchField2 {
            u_breaks: u_breaks.clone(),
            v_breaks: v_breaks.clone(),
            u_degree: p,
            v_degree: q,
            patches,
        })
    };
    Some([make(0), make(1), make(2), make(3)])
}

/// A box `(lo, hi)` holding a B-spline curve over `[a, b]`: the control
/// hull of each overlapping Bezier piece restricted to the range, in
/// homogeneous coordinates divided by the weight's bounds.
pub(crate) fn curve_hull(curve: &Curve3, a: Scalar, b: Scalar) -> Option<(Vec3, Vec3)> {
    let Curve3::BSpline(c) = curve else {
        return None;
    };
    let p = usize::from(c.degree);
    let knots = expand(&c.knots, &c.multiplicities);
    let control: Vec<H> = c
        .control_points
        .iter()
        .enumerate()
        .map(|(i, pt)| {
            let w = c.weights.as_ref().map_or(1.0, |ws| ws[i]);
            [pt.x * w, pt.y * w, pt.z * w, w]
        })
        .collect();
    let (breaks, segments) = decompose(&knots, p, &control)?;
    let mut lo = Vec3::splat(Scalar::INFINITY);
    let mut hi = Vec3::splat(Scalar::NEG_INFINITY);
    for (i, seg) in segments.iter().enumerate() {
        let (k0, k1) = (breaks[i], breaks[i + 1]);
        if k1 < a || k0 > b || k1 <= k0 {
            continue;
        }
        let (s0, s1) = (
            ((a - k0) / (k1 - k0)).clamp(0.0, 1.0),
            ((b - k0) / (k1 - k0)).clamp(0.0, 1.0),
        );
        let restricted: Vec<Vec<Scalar>> = (0..4)
            .map(|k| restrict1(&seg.iter().map(|h| h[k]).collect::<Vec<_>>(), s0, s1))
            .collect();
        let range = |v: &[Scalar]| {
            let size = v.iter().fold(0.0 as Scalar, |m, x| m.max(x.abs()));
            let pad = 32.0 * Scalar::EPSILON * size;
            (
                v.iter().copied().fold(Scalar::INFINITY, Scalar::min) - pad,
                v.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max) + pad,
            )
        };
        let (w_lo, w_hi) = range(&restricted[3]);
        if w_lo <= 0.0 {
            return None;
        }
        for k in 0..3 {
            let (x_lo, x_hi) = range(&restricted[k]);
            let q = [x_lo / w_lo, x_lo / w_hi, x_hi / w_lo, x_hi / w_hi];
            lo[k] = lo[k].min(q.iter().copied().fold(Scalar::INFINITY, Scalar::min));
            hi[k] = hi[k].max(q.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max));
        }
    }
    lo.is_finite().then_some((lo, hi))
}

/// A Bernstein polynomial's coefficients over `[s0, s1]` of its `[0, 1]`.
fn restrict1(c: &[Scalar], s0: Scalar, s1: Scalar) -> Vec<Scalar> {
    let n = c.len();
    let casteljau = |c: &[Scalar], s: Scalar, keep_left: bool| -> Vec<Scalar> {
        let mut w = c.to_vec();
        let mut out = vec![0.0; n];
        if keep_left {
            out[0] = w[0];
        } else {
            out[n - 1] = w[n - 1];
        }
        for k in 1..n {
            for i in 0..n - k {
                w[i] = w[i] * (1.0 - s) + w[i + 1] * s;
            }
            if keep_left {
                out[k] = w[0];
            } else {
                out[n - 1 - k] = w[n - 1 - k];
            }
        }
        out
    };
    if s1 <= s0 {
        let mut w = c.to_vec();
        for k in 1..n {
            for i in 0..n - k {
                w[i] = w[i] * (1.0 - s0) + w[i + 1] * s0;
            }
        }
        return vec![w[0]; n];
    }
    let left = if s1 < 1.0 {
        casteljau(c, s1, true)
    } else {
        c.to_vec()
    };
    if s0 <= 0.0 {
        return left;
    }
    casteljau(&left, s0 / s1, false)
}
