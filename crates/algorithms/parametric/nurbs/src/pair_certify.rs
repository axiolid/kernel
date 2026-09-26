//! Certificates for sections of two B-splines (#119, ADR 0077).
//!
//! Interval enclosures of a B-spline surface (or curve) and its first
//! partials over any parameter box come from Bernstein coefficients: the
//! point from its rational control points (positive weights keep the curve
//! in their hull), each partial as `(W X_u - W_u X) / W^2` with the
//! numerator's and `W^2`'s coefficients bounding them. Boxes past the
//! domain continue its edge polynomials. Every enclosure is widened by a
//! margin for the rounding of the coefficients.
//!
//! Krawczyk's test on those enclosures then proves two things:
//! - [`certify_chord`]: over a chord of a traced pair section, for every
//!   level of the plane sweeping across it, exactly one point of the
//!   section lies in a box, so the chord's two nodes lie on one arc and the
//!   curve between them is that arc;
//! - [`isolate`]: every crossing of a curve (a sub-patch edge, or a
//!   B-spline curve) with a surface, each in a box proven to hold exactly
//!   one; boxes proven to hold none are dropped, and the rest refused.

use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_curve::{BSplineSurface, PairNode};

use crate::spline_field::{bezier_net, decompose};

/// A homogeneous control point.
type H = [Scalar; 4];

/// A closed interval, rounded outward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct I {
    pub(crate) lo: Scalar,
    pub(crate) hi: Scalar,
}

fn down(x: Scalar) -> Scalar {
    x - x.abs() * 4.0 * Scalar::EPSILON - Scalar::MIN_POSITIVE
}

fn up(x: Scalar) -> Scalar {
    x + x.abs() * 4.0 * Scalar::EPSILON + Scalar::MIN_POSITIVE
}

impl I {
    pub(crate) fn new(a: Scalar, b: Scalar) -> Self {
        Self {
            lo: down(a.min(b)),
            hi: up(a.max(b)),
        }
    }

    pub(crate) fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    /// `x` give or take `r`.
    pub(crate) fn around(x: Scalar, r: Scalar) -> Self {
        Self::new(x - r.abs(), x + r.abs())
    }

    fn add(self, o: Self) -> Self {
        Self::new(self.lo + o.lo, self.hi + o.hi)
    }

    fn sub(self, o: Self) -> Self {
        Self::new(self.lo - o.hi, self.hi - o.lo)
    }

    fn mul(self, o: Self) -> Self {
        let c = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self::new(
            c.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            c.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        )
    }

    fn scale(self, k: Scalar) -> Self {
        Self::new(self.lo * k, self.hi * k)
    }

    /// `self / o` for `o` bounded away from zero.
    fn div(self, o: Self) -> Option<Self> {
        if o.lo <= 0.0 && o.hi >= 0.0 {
            return None;
        }
        let c = [
            self.lo / o.lo,
            self.lo / o.hi,
            self.hi / o.lo,
            self.hi / o.hi,
        ];
        Some(Self::new(
            c.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            c.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        ))
    }

    fn hull(self, o: Self) -> Self {
        Self {
            lo: self.lo.min(o.lo),
            hi: self.hi.max(o.hi),
        }
    }

    fn meets(self, o: Self) -> bool {
        self.lo <= o.hi && o.lo <= self.hi
    }

    /// Strictly inside `o`.
    fn within(self, o: Self) -> bool {
        self.lo > o.lo && self.hi < o.hi
    }

    fn mid(self) -> Scalar {
        0.5 * (self.lo + self.hi)
    }

    fn is_finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }
}

/// The interval hull of values, widened by a margin for their rounding.
fn range(values: impl Iterator<Item = Scalar>) -> I {
    let (mut lo, mut hi, mut size) = (Scalar::INFINITY, Scalar::NEG_INFINITY, 0.0 as Scalar);
    for x in values {
        lo = lo.min(x);
        hi = hi.max(x);
        size = size.max(x.abs());
    }
    let margin = 1e-12 * size;
    I::new(lo - margin, hi + margin)
}

/// A tensor Bernstein polynomial's coefficients, `c[a][b]`.
type Poly = Vec<Vec<Scalar>>;

/// The Bernstein coefficients of a polynomial of degree `c.len() - 1` over
/// `[s0, s1]` (any reals: past `[0, 1]` it continues), by blossoming.
fn restrict1(c: &[Scalar], s0: Scalar, s1: Scalar) -> Vec<Scalar> {
    let n = c.len() - 1;
    (0..=n)
        .map(|i| {
            let mut w = c.to_vec();
            for k in 1..=n {
                let t = if k <= i { s1 } else { s0 };
                for j in 0..=n - k {
                    w[j] = w[j] * (1.0 - t) + w[j + 1] * t;
                }
            }
            w[0]
        })
        .collect()
}

pub(crate) fn restrict2(c: &Poly, u: (Scalar, Scalar), v: (Scalar, Scalar)) -> Poly {
    let rows: Poly = c.iter().map(|row| restrict1(row, v.0, v.1)).collect();
    let (p, q) = (rows.len(), rows[0].len());
    let mut out = vec![vec![0.0; q]; p];
    for b in 0..q {
        let column: Vec<Scalar> = rows.iter().map(|r| r[b]).collect();
        for (a, x) in restrict1(&column, u.0, u.1).into_iter().enumerate() {
            out[a][b] = x;
        }
    }
    out
}

fn binomial(n: usize, k: usize) -> Scalar {
    let mut r = 1.0;
    for i in 0..k {
        r = r * (n - i) as Scalar / (i + 1) as Scalar;
    }
    r
}

fn du(c: &Poly) -> Poly {
    let p = c.len() - 1;
    if p == 0 {
        return vec![vec![0.0; c[0].len()]];
    }
    (0..p)
        .map(|a| {
            (0..c[0].len())
                .map(|b| p as Scalar * (c[a + 1][b] - c[a][b]))
                .collect()
        })
        .collect()
}

fn dv(c: &Poly) -> Poly {
    let q = c[0].len() - 1;
    if q == 0 {
        return vec![vec![0.0]; c.len()];
    }
    c.iter()
        .map(|row| {
            (0..q)
                .map(|b| q as Scalar * (row[b + 1] - row[b]))
                .collect()
        })
        .collect()
}

fn mul(x: &Poly, y: &Poly) -> Poly {
    let (p1, q1) = (x.len() - 1, x[0].len() - 1);
    let (p2, q2) = (y.len() - 1, y[0].len() - 1);
    let (p, q) = (p1 + p2, q1 + q2);
    let mut out = vec![vec![0.0; q + 1]; p + 1];
    for a1 in 0..=p1 {
        for b1 in 0..=q1 {
            let xv = x[a1][b1] * binomial(p1, a1) * binomial(q1, b1);
            if xv == 0.0 {
                continue;
            }
            for a2 in 0..=p2 {
                for b2 in 0..=q2 {
                    out[a1 + a2][b1 + b2] += xv * y[a2][b2] * binomial(p2, a2) * binomial(q2, b2);
                }
            }
        }
    }
    for (a, row) in out.iter_mut().enumerate() {
        for (b, v) in row.iter_mut().enumerate() {
            *v /= binomial(p, a) * binomial(q, b);
        }
    }
    out
}

fn sub(x: &Poly, y: &Poly) -> Poly {
    x.iter()
        .zip(y)
        .map(|(r, s)| r.iter().zip(s).map(|(a, b)| a - b).collect())
        .collect()
}

/// One Bezier cell of a B-spline and the polynomials its enclosures read.
#[derive(Debug, Clone)]
struct Piece {
    lo: Point2,
    hi: Point2,
    /// Whether the cell's sides are the domain's own, past which its
    /// polynomials continue: `[u low, u high, v low, v high]`.
    open: [bool; 4],
    /// Homogeneous components `X, Y, Z, W`.
    x: [Poly; 4],
    /// `W X_u - W_u X` per coordinate, in local parameters.
    nu: [Poly; 3],
    /// `W X_v - W_v X`.
    nv: [Poly; 3],
    /// `W^2`.
    w2: Poly,
}

/// A point's and its partials' enclosures.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Jet {
    pub(crate) p: [I; 3],
    pub(crate) u: [I; 3],
    pub(crate) v: [I; 3],
}

impl Jet {
    fn hull(self, o: Self) -> Self {
        let h = |a: [I; 3], b: [I; 3]| [a[0].hull(b[0]), a[1].hull(b[1]), a[2].hull(b[2])];
        Self {
            p: h(self.p, o.p),
            u: h(self.u, o.u),
            v: h(self.v, o.v),
        }
    }
}

/// Enclosures of one B-spline surface (or curve, with `v` fixed at 0).
#[derive(Debug, Clone)]
pub(crate) struct Enclosure {
    pieces: Vec<Piece>,
}

fn piece(net: &[Vec<H>], lo: Point2, hi: Point2, open: [bool; 4]) -> Option<Piece> {
    if !net.iter().flatten().all(|h| h[3] > 0.0) {
        return None;
    }
    let comp = |k: usize| -> Poly {
        net.iter()
            .map(|row| row.iter().map(|h| h[k]).collect())
            .collect()
    };
    let x = [comp(0), comp(1), comp(2), comp(3)];
    let w = &x[3];
    let (wu, wv) = (du(w), dv(w));
    let nu = [0, 1, 2].map(|k| sub(&mul(w, &du(&x[k])), &mul(&wu, &x[k])));
    let nv = [0, 1, 2].map(|k| sub(&mul(w, &dv(&x[k])), &mul(&wv, &x[k])));
    let w2 = mul(w, w);
    Some(Piece {
        lo,
        hi,
        open,
        x,
        nu,
        nv,
        w2,
    })
}

impl Enclosure {
    /// A clamped B-spline surface's enclosures; `None` for an unclamped one
    /// or a weight that is not positive.
    pub(crate) fn surface(b: &BSplineSurface) -> Option<Self> {
        let (_, _, ub, vb, net) = bezier_net(b)?;
        let (nu, nv) = (ub.len() - 1, vb.len() - 1);
        let mut pieces = Vec::with_capacity(nu * nv);
        for (iu, row) in net.iter().enumerate() {
            for (jv, cell) in row.iter().enumerate() {
                pieces.push(piece(
                    cell,
                    Point2::new(ub[iu], vb[jv]),
                    Point2::new(ub[iu + 1], vb[jv + 1]),
                    [iu == 0, iu + 1 == nu, jv == 0, jv + 1 == nv],
                )?);
            }
        }
        Some(Self { pieces })
    }

    /// A clamped B-spline curve's enclosures, read as a surface constant in
    /// `v`: its parameter is `u`, and `v` is `0`.
    pub(crate) fn curve(knots: &[Scalar], degree: usize, control: &[H]) -> Option<Self> {
        let (breaks, segments) = decompose(knots, degree, control)?;
        let n = segments.len();
        let mut pieces = Vec::with_capacity(n);
        for (i, seg) in segments.iter().enumerate() {
            let net: Vec<Vec<H>> = seg.iter().map(|h| vec![*h]).collect();
            pieces.push(piece(
                &net,
                Point2::new(breaks[i], 0.0),
                Point2::new(breaks[i + 1], 0.0),
                [i == 0, i + 1 == n, true, true],
            )?);
        }
        Some(Self { pieces })
    }

    /// The point and its first partials at `p`, from the same polynomials
    /// the enclosures bound (past the domain, its edge cells continue).
    pub(crate) fn at(&self, p: Point2) -> Option<(Point3, Vec3, Vec3)> {
        let c = self.pieces.iter().find(|c| {
            (c.open[0] || p.x >= c.lo.x)
                && (c.open[1] || p.x <= c.hi.x)
                && (c.open[2] || p.y >= c.lo.y)
                && (c.open[3] || p.y <= c.hi.y)
        })?;
        let (wu, wv) = (c.hi.x - c.lo.x, c.hi.y - c.lo.y);
        let local = |x: Scalar, o: Scalar, w: Scalar| if w == 0.0 { 0.0 } else { (x - o) / w };
        let (su, sv) = (local(p.x, c.lo.x, wu), local(p.y, c.lo.y, wv));
        let value = |poly: &Poly| restrict2(poly, (su, su), (sv, sv))[0][0];
        let w = value(&c.x[3]);
        let w2 = value(&c.w2);
        if !(w > 0.0 && w2 > 0.0) {
            return None;
        }
        let point = Point3::new(value(&c.x[0]) / w, value(&c.x[1]) / w, value(&c.x[2]) / w);
        let partial = |n: &[Poly; 3], width: Scalar| {
            if width == 0.0 {
                Vec3::ZERO
            } else {
                Vec3::new(value(&n[0]), value(&n[1]), value(&n[2])) / (w2 * width)
            }
        };
        Some((point, partial(&c.nu, wu), partial(&c.nv, wv)))
    }

    /// Enclosures of the point and its partials over `[lo, hi]`.
    pub(crate) fn jet(&self, lo: Point2, hi: Point2) -> Option<Jet> {
        let mut out: Option<Jet> = None;
        for c in &self.pieces {
            // The part of the box this cell answers for.
            let (a0, a1) = (
                if c.open[0] { lo.x } else { lo.x.max(c.lo.x) },
                if c.open[1] { hi.x } else { hi.x.min(c.hi.x) },
            );
            let (b0, b1) = (
                if c.open[2] { lo.y } else { lo.y.max(c.lo.y) },
                if c.open[3] { hi.y } else { hi.y.min(c.hi.y) },
            );
            if a0 > a1 || b0 > b1 {
                continue;
            }
            let (wu, wv) = (c.hi.x - c.lo.x, c.hi.y - c.lo.y);
            let local = |x: Scalar, o: Scalar, w: Scalar| if w == 0.0 { 0.0 } else { (x - o) / w };
            let su = (local(a0, c.lo.x, wu), local(a1, c.lo.x, wu));
            let sv = (local(b0, c.lo.y, wv), local(b1, c.lo.y, wv));
            let x: Vec<Poly> = c.x.iter().map(|p| restrict2(p, su, sv)).collect();
            let w = &x[3];
            if !w.iter().flatten().all(|v| *v > 0.0) {
                return None;
            }
            let p = [0, 1, 2].map(|k| {
                range(
                    x[k].iter()
                        .flatten()
                        .zip(w.iter().flatten())
                        .map(|(a, b)| a / b),
                )
            });
            let w2 = range(restrict2(&c.w2, su, sv).into_iter().flatten());
            let partial = |n: &Poly, width: Scalar| -> Option<I> {
                if width == 0.0 {
                    return Some(I::point(0.0));
                }
                range(restrict2(n, su, sv).into_iter().flatten())
                    .div(w2)
                    .map(|r| r.scale(1.0 / width))
            };
            let u = [
                partial(&c.nu[0], wu)?,
                partial(&c.nu[1], wu)?,
                partial(&c.nu[2], wu)?,
            ];
            let v = [
                partial(&c.nv[0], wv)?,
                partial(&c.nv[1], wv)?,
                partial(&c.nv[2], wv)?,
            ];
            let jet = Jet { p, u, v };
            out = Some(match out {
                Some(o) => o.hull(jet),
                None => jet,
            });
        }
        out.filter(|j| j.p.iter().chain(&j.u).chain(&j.v).all(|i| i.is_finite()))
    }
}

/// A 4x4 interval matrix times a 4-vector of intervals.
fn apply4(m: &[[I; 4]; 4], x: &[I; 4]) -> [I; 4] {
    [0, 1, 2, 3].map(|i| (0..4).fold(I::point(0.0), |acc, j| acc.add(m[i][j].mul(x[j]))))
}

/// The inverse of a point matrix, column by column.
fn inverse<const N: usize>(
    m: [[Scalar; N]; N],
    solve: impl Fn([[Scalar; N]; N], [Scalar; N]) -> Option<[Scalar; N]>,
) -> Option<[[Scalar; N]; N]> {
    let mut out = [[0.0; N]; N];
    for j in 0..N {
        let mut e = [0.0; N];
        e[j] = 1.0;
        let col = solve(m, e)?;
        for i in 0..N {
            out[i][j] = col[i];
        }
    }
    Some(out)
}

/// A margin for a point evaluated from coefficients of this size.
fn rounding(p: Point3) -> Scalar {
    1e-12 * (1.0 + p.x.abs().max(p.y.abs()).max(p.z.abs()))
}

/// The box `X` about the chord from `n0` to `n1` in both surfaces'
/// parameters `(a, b)` within which, for every level `s` in `[0, 1]` of the
/// plane `d . (P - n0) = s |d|^2` (`d = n1 - n0`), exactly one point lies on
/// both surfaces: Krawczyk's test on `S1(a) - S2(b)` and the plane, with
/// `s` an interval. `None` where no box tried passes.
#[allow(clippy::needless_range_loop)]
pub(crate) fn certify_chord(
    e1: &Enclosure,
    e2: &Enclosure,
    n0: &PairNode,
    n1: &PairNode,
    centre: (Point2, Point2),
) -> Option<[I; 4]> {
    let d = n1.point - n0.point;
    let dd = d.dot(d);
    if dd == 0.0 || !dd.is_finite() {
        return None;
    }
    let x0 = [centre.0.x, centre.0.y, centre.1.x, centre.1.y];
    let ends = [
        [n0.first.x, n0.first.y, n0.second.x, n0.second.y],
        [n1.first.x, n1.first.y, n1.second.x, n1.second.y],
    ];
    // Radii that reach both nodes from the centre, with room to contract.
    let reach: [Scalar; 4] =
        [0, 1, 2, 3].map(|k| (ends[0][k] - x0[k]).abs().max((ends[1][k] - x0[k]).abs()));
    let widest = reach.iter().copied().fold(0.0, Scalar::max);
    // Point jets and the preconditioner at the centre.
    let (p1, u1, v1) = e1.at(Point2::new(x0[0], x0[1]))?;
    let (p2, u2, v2) = e2.at(Point2::new(x0[2], x0[3]))?;
    let m = [
        [u1.x, v1.x, -u2.x, -v2.x],
        [u1.y, v1.y, -u2.y, -v2.y],
        [u1.z, v1.z, -u2.z, -v2.z],
        [d.dot(u1), d.dot(v1), 0.0, 0.0],
    ];
    let y = inverse(m, axiolid_curve::pair_section::solve4)?;
    // f at the centre, for every level: the gap and the plane.
    let gap = p1 - p2;
    let margin = rounding(p1).max(rounding(p2));
    let lift = d.dot(p1 - n0.point);
    let f = [
        I::around(gap.x, margin),
        I::around(gap.y, margin),
        I::around(gap.z, margin),
        I::new(lift - dd, lift).add(I::around(0.0, margin * d.length())),
    ];
    let yf: [I; 4] =
        [0, 1, 2, 3].map(|i| (0..4).fold(I::point(0.0), |acc, j| acc.add(f[j].scale(y[i][j]))));
    let dv = [d.x, d.y, d.z];
    for grow in [1.25, 2.0, 4.0] {
        let r: [Scalar; 4] = [0, 1, 2, 3]
            .map(|k| grow * reach[k] + 0.1 * grow * widest + 1e-12 * (1.0 + x0[k].abs()));
        let x: [I; 4] = [0, 1, 2, 3].map(|k| I::around(x0[k], r[k]));
        let (Some(a), Some(b)) = (
            e1.jet(Point2::new(x[0].lo, x[1].lo), Point2::new(x[0].hi, x[1].hi)),
            e2.jet(Point2::new(x[2].lo, x[3].lo), Point2::new(x[2].hi, x[3].hi)),
        ) else {
            continue;
        };
        // J over the box.
        let mut jm = [[I::point(0.0); 4]; 4];
        for k in 0..3 {
            jm[k] = [a.u[k], a.v[k], b.u[k].scale(-1.0), b.v[k].scale(-1.0)];
        }
        let dot = |g: [I; 3]| (0..3).fold(I::point(0.0), |acc, k| acc.add(g[k].scale(dv[k])));
        jm[3] = [dot(a.u), dot(a.v), I::point(0.0), I::point(0.0)];
        // I - Y J.
        let mut c = [[I::point(0.0); 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                let yj = (0..4).fold(I::point(0.0), |acc, k| acc.add(jm[k][j].scale(y[i][k])));
                c[i][j] = I::point(if i == j { 1.0 } else { 0.0 }).sub(yj);
            }
        }
        let offset: [I; 4] = [0, 1, 2, 3].map(|k| x[k].sub(I::point(x0[k])));
        let spread = apply4(&c, &offset);
        let k: [I; 4] = [0, 1, 2, 3].map(|i| I::point(x0[i]).sub(yf[i]).add(spread[i]));
        if (0..4).all(|i| k[i].within(x[i])) {
            // Both nodes must lie in the box: they are its solutions at
            // levels 0 and 1.
            let holds = |e: &[Scalar; 4]| (0..4).all(|i| e[i] > x[i].lo && e[i] < x[i].hi);
            if holds(&ends[0]) && holds(&ends[1]) {
                return Some(x);
            }
        }
    }
    None
}

/// Whether a node lies in a certified chord box, at a level of its chord.
pub(crate) fn in_chord(x: &[I; 4], n0: &PairNode, n1: &PairNode, node: &PairNode) -> bool {
    let v = [node.first.x, node.first.y, node.second.x, node.second.y];
    if !(0..4).all(|i| v[i] >= x[i].lo && v[i] <= x[i].hi) {
        return false;
    }
    let d = n1.point - n0.point;
    let dd = d.dot(d);
    let s = d.dot(node.point - n0.point) / dd;
    (-1e-9..=1.0 + 1e-9).contains(&s)
}

/// A curve crossing a surface: the curve's parameter, the surface's, and
/// the point.
pub(crate) type Crossing = (Scalar, Point2, Point3);

/// Every crossing of a curve with a surface, for curve parameters in `t`
/// and surface parameters in `[lo, hi]`, each proven unique in a box.
///
/// `curve(t)` gives the curve's point and derivative, `curve_box` their
/// enclosures over an interval. Boxes proven to hold no crossing are
/// dropped; one that can be neither proven empty nor proven to hold one
/// crossing within `depth` halvings -- the curve touching the surface --
/// makes the whole answer `None`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn isolate(
    t: (Scalar, Scalar),
    lo: Point2,
    hi: Point2,
    curve: &dyn Fn(Scalar) -> Option<(Point3, Vec3)>,
    curve_box: &dyn Fn(Scalar, Scalar) -> Option<([I; 3], [I; 3])>,
    enclosure: &Enclosure,
) -> Option<Vec<Crossing>> {
    let mut out: Vec<Crossing> = Vec::new();
    let mut queue = vec![([t.0, t.1], [lo.x, hi.x], [lo.y, hi.y], 0u32)];
    let mut work = 0usize;
    while let Some((tb, ub, vb, depth)) = queue.pop() {
        work += 1;
        if work > 50_000 {
            return None;
        }
        let x = [
            I::new(tb[0], tb[1]),
            I::new(ub[0], ub[1]),
            I::new(vb[0], vb[1]),
        ];
        let (cp, cd) = curve_box(tb[0], tb[1])?;
        let sj = enclosure.jet(Point2::new(ub[0], vb[0]), Point2::new(ub[1], vb[1]))?;
        // Apart: no crossing here.
        if !(0..3).all(|k| cp[k].meets(sj.p[k])) {
            continue;
        }
        match krawczyk3(x, curve, cd, sj, enclosure) {
            Verdict::None => continue,
            Verdict::One(root) => {
                push_root(&mut out, root);
                continue;
            }
            Verdict::Unknown => {}
        }
        if depth > 40 {
            let centre = [
                0.5 * (tb[0] + tb[1]),
                0.5 * (ub[0] + ub[1]),
                0.5 * (vb[0] + vb[1]),
            ];
            let inside_asked = |r: &[Scalar; 3]| {
                r[0] >= t.0
                    && r[0] <= t.1
                    && r[1] >= lo.x
                    && r[1] <= hi.x
                    && r[2] >= lo.y
                    && r[2] <= hi.y
            };
            // A box this small that can be neither proven empty nor
            // proven to hold one crossing: the curve touches the surface
            // here (its crossing is double, the Jacobian singular), or
            // passes within rounding of it. The touching point, where a
            // damped Newton finds one, is kept; none found means no
            // crossing above rounding.
            let Some(root) = newton3(centre, curve, enclosure) else {
                match touch3(centre, curve, enclosure) {
                    Some(r) if !tangent(&r, curve, enclosure) => return None,
                    Some(r) => {
                        if inside_asked(&r.0) {
                            push_root(&mut out, r);
                        }
                    }
                    None => {}
                }
                continue;
            };
            let w = [tb[1] - tb[0], ub[1] - ub[0], vb[1] - vb[0]];
            let grown = [0, 1, 2]
                .map(|k| I::around(root.0[k], 3.0 * w[k] + 1e-13 * (1.0 + root.0[k].abs())));
            let covers = grown[0].lo <= tb[0]
                && grown[0].hi >= tb[1]
                && grown[1].lo <= ub[0]
                && grown[1].hi >= ub[1]
                && grown[2].lo <= vb[0]
                && grown[2].hi >= vb[1];
            let (_, gd) = curve_box(grown[0].lo, grown[0].hi)?;
            let gj = enclosure.jet(
                Point2::new(grown[1].lo, grown[2].lo),
                Point2::new(grown[1].hi, grown[2].hi),
            )?;
            match krawczyk3(grown, curve, gd, gj, enclosure) {
                Verdict::One(r) if covers => {
                    // Keep it only where it lies in the box asked about.
                    if inside_asked(&r.0) {
                        push_root(&mut out, r);
                    }
                    continue;
                }
                Verdict::None if covers => continue,
                // Newton reached a point it cannot prove alone: kept where
                // the curve touches the surface there (the Jacobian is
                // singular, so no proof can exist); a regular crossing
                // that cannot be proven is refused.
                _ if tangent(&root, curve, enclosure) => {
                    if inside_asked(&root.0) {
                        push_root(&mut out, root);
                    }
                    continue;
                }
                _ => return None,
            }
        }
        // Split the parameter whose image is largest.
        let size = |i: [I; 3]| (0..3).map(|k| i[k].hi - i[k].lo).fold(0.0, Scalar::max);
        let along_curve = size(cp) >= size(sj.p);
        let f = 0.5;
        if along_curve {
            let m = tb[0] + (tb[1] - tb[0]) * f;
            queue.push(([tb[0], m], ub, vb, depth + 1));
            queue.push(([m, tb[1]], ub, vb, depth + 1));
        } else {
            let mu = ub[0] + (ub[1] - ub[0]) * f;
            let mv = vb[0] + (vb[1] - vb[0]) * f;
            for (a, b) in [(ub[0], mu), (mu, ub[1])] {
                for (c, e) in [(vb[0], mv), (mv, vb[1])] {
                    queue.push((tb, [a, b], [c, e], depth + 1));
                }
            }
        }
    }
    Some(out)
}

fn push_root(out: &mut Vec<Crossing>, root: ([Scalar; 3], Point3)) {
    let (x, p) = root;
    if !out.iter().any(|(t, uv, _)| {
        (t - x[0]).abs() <= 1e-10 * (1.0 + t.abs())
            && (uv.x - x[1]).abs() <= 1e-10 * (1.0 + uv.x.abs())
            && (uv.y - x[2]).abs() <= 1e-10 * (1.0 + uv.y.abs())
    }) {
        out.push((x[0], Point2::new(x[1], x[2]), p));
    }
}

/// What Krawczyk's test says of a box.
enum Verdict {
    None,
    One(([Scalar; 3], Point3)),
    Unknown,
}

/// Newton on `C(t) = S(u, v)` from `x`.
fn newton3(
    mut x: [Scalar; 3],
    curve: &dyn Fn(Scalar) -> Option<(Point3, Vec3)>,
    surface: &Enclosure,
) -> Option<([Scalar; 3], Point3)> {
    for _ in 0..60 {
        let (c, dc) = curve(x[0])?;
        let (sp, su, sv) = surface.at(Point2::new(x[1], x[2]))?;
        let r = c - sp;
        let m = [
            [dc.x, -su.x, -sv.x],
            [dc.y, -su.y, -sv.y],
            [dc.z, -su.z, -sv.z],
        ];
        let step = crate::pair_trace::solve3(m, [-r.x, -r.y, -r.z])?;
        for k in 0..3 {
            x[k] += step[k];
        }
        if step.iter().map(|v| v.abs()).fold(0.0, Scalar::max)
            <= 4.0 * Scalar::EPSILON * (1.0 + x.iter().map(|v| v.abs()).fold(0.0, Scalar::max))
        {
            break;
        }
    }
    let (c, _) = curve(x[0])?;
    let miss = (c - surface.at(Point2::new(x[1], x[2]))?.0).length();
    (miss <= 1e-9 * (1.0 + c.length())).then_some((x, c))
}

/// Whether the curve touches the surface at `root`: its tangent lies in
/// the surface's tangent plane, to within a millionth of a radian.
fn tangent(
    root: &([Scalar; 3], Point3),
    curve: &dyn Fn(Scalar) -> Option<(Point3, Vec3)>,
    surface: &Enclosure,
) -> bool {
    let (Some((_, dc)), Some((_, su, sv))) = (
        curve(root.0[0]),
        surface.at(Point2::new(root.0[1], root.0[2])),
    ) else {
        return false;
    };
    let n = su.cross(sv);
    let (ln, lc) = (n.length(), dc.length());
    ln > 0.0 && lc > 0.0 && (n.dot(dc) / (ln * lc)).abs() <= 1e-6
}

/// Levenberg-Marquardt on `|C(t) - S(u, v)|^2` from `x`: converges, if
/// only linearly, where the curve touches the surface and Newton's
/// Jacobian is singular. The point where the gap vanishes to rounding.
#[allow(clippy::needless_range_loop)]
fn touch3(
    mut x: [Scalar; 3],
    curve: &dyn Fn(Scalar) -> Option<(Point3, Vec3)>,
    surface: &Enclosure,
) -> Option<([Scalar; 3], Point3)> {
    for _ in 0..400 {
        let (c, dc) = curve(x[0])?;
        let (sp, su, sv) = surface.at(Point2::new(x[1], x[2]))?;
        let r = c - sp;
        let cols = [dc, -su, -sv];
        // (J^T J + mu I) step = -J^T r.
        let mut m = [[0.0; 3]; 3];
        let mut g = [0.0; 3];
        let mut scale = 0.0 as Scalar;
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] = cols[i].dot(cols[j]);
            }
            g[i] = -cols[i].dot(r);
            scale = scale.max(m[i][i]);
        }
        for i in 0..3 {
            m[i][i] += 1e-12 * scale;
        }
        let step = crate::pair_trace::solve3(m, g)?;
        for k in 0..3 {
            x[k] += step[k];
        }
        if step.iter().map(|v| v.abs()).fold(0.0, Scalar::max)
            <= 4.0 * Scalar::EPSILON * (1.0 + x.iter().map(|v| v.abs()).fold(0.0, Scalar::max))
        {
            break;
        }
    }
    let (c, _) = curve(x[0])?;
    let miss = (c - surface.at(Point2::new(x[1], x[2]))?.0).length();
    (miss <= 1e-12 * (1.0 + c.length())).then_some((x, c))
}

/// Krawczyk's test for `C(t) - S(u, v) = 0` over the box `x`, given
/// enclosures of `C'` over its `t` and of `S`'s partials over its `(u, v)`.
#[allow(clippy::needless_range_loop)]
fn krawczyk3(
    x: [I; 3],
    curve: &dyn Fn(Scalar) -> Option<(Point3, Vec3)>,
    cd: [I; 3],
    sj: Jet,
    surface: &Enclosure,
) -> Verdict {
    let x0 = [x[0].mid(), x[1].mid(), x[2].mid()];
    let (Some((c, dc)), Some((sp, su, sv))) = (curve(x0[0]), surface.at(Point2::new(x0[1], x0[2])))
    else {
        return Verdict::Unknown;
    };
    let m = [
        [dc.x, -su.x, -sv.x],
        [dc.y, -su.y, -sv.y],
        [dc.z, -su.z, -sv.z],
    ];
    let Some(y) = inverse(m, crate::pair_trace::solve3) else {
        return Verdict::Unknown;
    };
    let gap = c - sp;
    let margin = rounding(c).max(rounding(sp));
    let f = [
        I::around(gap.x, margin),
        I::around(gap.y, margin),
        I::around(gap.z, margin),
    ];
    let mut k = [I::point(0.0); 3];
    for i in 0..3 {
        let yf = (0..3).fold(I::point(0.0), |acc, jj| acc.add(f[jj].scale(y[i][jj])));
        let mut spread = I::point(0.0);
        for jj in 0..3 {
            // (I - Y J)_{i jj}
            let yj = (0..3).fold(I::point(0.0), |acc, kk| {
                let entry = match jj {
                    0 => cd[kk],
                    1 => sj.u[kk].scale(-1.0),
                    _ => sj.v[kk].scale(-1.0),
                };
                acc.add(entry.scale(y[i][kk]))
            });
            let cij = I::point(if i == jj { 1.0 } else { 0.0 }).sub(yj);
            spread = spread.add(cij.mul(x[jj].sub(I::point(x0[jj]))));
        }
        k[i] = I::point(x0[i]).sub(yf).add(spread);
    }
    if (0..3).any(|i| !k[i].meets(x[i])) {
        return Verdict::None;
    }
    if (0..3).all(|i| k[i].within(x[i])) {
        // One crossing, and Newton from the centre stays with it.
        if let Some(root) = newton3(x0, curve, surface) {
            if (0..3).all(|i| root.0[i] >= x[i].lo && root.0[i] <= x[i].hi) {
                return Verdict::One(root);
            }
        }
    }
    Verdict::Unknown
}
