//! Curves given implicitly: where a field over a surface's parameters is
//! zero (ADR 0077).
//!
//! Where two analytic surfaces meet, the section read in one surface's
//! parameters `(u, v)` is the zero set of the other surface's implicit
//! equation composed with the first one's parameterisation. For planes,
//! quadrics and tori that composition is a [`Field2`]: a finite sum of
//! products of powers or harmonics of `u` and of `v`, known exactly from
//! the two surfaces. Most such sections have no closed form -- a torus
//! against a cylinder is a quartic in space -- but the field always has.
//!
//! An [`ImplicitCurve2`] is one connected stretch of such a zero set,
//! carried as a chain of [`ImplicitCell`]s. In each cell the field is
//! strictly monotone along one parameter, so at every value of the other
//! parameter the curve is the field's *unique* zero in the cell's bracket.
//! A point on the curve is therefore defined, not approximated: the root
//! is found to full precision by a safeguarded Newton iteration that
//! cannot leave the bracket or pick a different branch, and derivatives
//! follow from the implicit function theorem.
//!
//! [`ImplicitSection3`] carries the same curve in space, on its
//! [`Carrier`] surface.

use axiolid_core::{Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use core::f64::consts::{PI, TAU};

use crate::quadric_section::RuledCarrier;
use crate::torus_section::TorusCarrier;

/// How a [`SeriesField2`] varies along one parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Basis {
    /// Powers: term `k` is `x^k`.
    Power,
    /// Harmonics: term `0` is `1`, term `2k - 1` is `cos(k x)` and term
    /// `2k` is `sin(k x)`.
    Fourier,
}

impl Basis {
    /// The first `n` terms' values at `x` into `out`, exactly as
    /// [`Basis::terms`] computes them.
    fn values_into(self, x: Scalar, out: &mut [Scalar]) {
        let n = out.len();
        match self {
            Basis::Power => {
                let mut p = 1.0;
                for slot in out.iter_mut() {
                    *slot = p;
                    p *= x;
                }
            }
            Basis::Fourier => {
                if n > 0 {
                    out[0] = 1.0;
                }
                let mut k = 1;
                while 2 * k - 1 < n {
                    let w = k as Scalar;
                    let (s, c) = (w * x).sin_cos();
                    out[2 * k - 1] = c;
                    if 2 * k < n {
                        out[2 * k] = s;
                    }
                    k += 1;
                }
            }
        }
    }

    /// Values, first and second derivatives of the first `n` terms at `x`.
    fn terms(self, x: Scalar, n: usize) -> (Vec<Scalar>, Vec<Scalar>, Vec<Scalar>) {
        let (mut f, mut d, mut dd) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        match self {
            Basis::Power => {
                let mut p = 1.0;
                for slot in f.iter_mut() {
                    *slot = p;
                    p *= x;
                }
                for k in 1..n {
                    d[k] = k as Scalar * f[k - 1];
                }
                for k in 2..n {
                    dd[k] = (k * (k - 1)) as Scalar * f[k - 2];
                }
            }
            Basis::Fourier => {
                if n > 0 {
                    f[0] = 1.0;
                }
                let mut k = 1;
                while 2 * k - 1 < n {
                    let w = k as Scalar;
                    let (s, c) = (w * x).sin_cos();
                    f[2 * k - 1] = c;
                    d[2 * k - 1] = -w * s;
                    dd[2 * k - 1] = -w * w * c;
                    if 2 * k < n {
                        f[2 * k] = s;
                        d[2 * k] = w * c;
                        dd[2 * k] = -w * w * s;
                    }
                    k += 1;
                }
            }
        }
        (f, d, dd)
    }
}

/// A field over the parameter plane: `sum c[i][j] B_i(u) B_j(v)`, with the
/// bases of [`Basis`] along each parameter. The form a plane's, quadric's or
/// torus's equation takes on an analytic surface.
#[derive(Debug, Clone, PartialEq)]
pub struct SeriesField2 {
    /// The basis along `u`.
    pub u: Basis,
    /// The basis along `v`.
    pub v: Basis,
    /// `coefficients[i][j]` multiplies term `i` in `u` and term `j` in `v`.
    pub coefficients: Vec<Vec<Scalar>>,
}

/// A field over a surface's parameters whose zero set is a section curve:
/// another surface's equation read in this surface's parameters.
#[derive(Debug, Clone, PartialEq)]
pub enum Field2 {
    /// Powers and harmonics, on an analytic surface.
    Series(SeriesField2),
    /// Piecewise Bernstein polynomials, on a B-spline surface.
    Patches(PatchField2),
}

impl Field2 {
    /// The field's value at `p`.
    #[must_use]
    pub fn value(&self, p: Point2) -> Scalar {
        match self {
            Self::Series(f) => f.value(p),
            Self::Patches(f) => f.jet(p).value,
        }
    }

    /// Value, gradient and Hessian at `p`.
    #[must_use]
    pub fn jet(&self, p: Point2) -> Jet2 {
        match self {
            Self::Series(f) => f.jet(p),
            Self::Patches(f) => f.jet(p),
        }
    }

    /// A bound on the field's size over its whole domain, and the scale
    /// its rounding is measured in.
    #[must_use]
    pub fn magnitude(&self) -> Scalar {
        match self {
            Self::Series(f) => f.magnitude(),
            Self::Patches(f) => f.magnitude(),
        }
    }

    /// The sum of the magnitudes of the terms that make up the value at
    /// `p`: the scale the value's rounding is measured in there. A power
    /// series read far from its origin has terms much larger than its
    /// coefficients, and a value that cancels them to near zero carries
    /// their rounding.
    #[must_use]
    pub fn scale_at(&self, p: Point2) -> Scalar {
        match self {
            Self::Series(f) => {
                let (n, m) = f.size();
                let (fu, _, _) = f.u.terms(p.x, n);
                let (fv, _, _) = f.v.terms(p.y, m);
                let mut sum = 0.0;
                for (i, row) in f.coefficients.iter().enumerate() {
                    for (j, &c) in row.iter().enumerate() {
                        sum += (c * fu[i] * fv[j]).abs();
                    }
                }
                sum
            }
            // Bernstein weights sum to one: the coefficients bound it.
            Self::Patches(f) => f.magnitude(),
        }
    }

    /// Whether every coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        match self {
            Self::Series(f) => f.is_finite(),
            Self::Patches(f) => f.is_finite(),
        }
    }
}

/// A piecewise polynomial field: on each cell of the grid of `u_breaks` by
/// `v_breaks`, a tensor-product Bernstein polynomial of degree
/// `(u_degree, v_degree)` in the cell's local coordinates `s`, `t` in
/// `[0, 1]`. Its coefficients bound it (the convex hull property), which is
/// what makes a trace on it certified.
#[derive(Debug, Clone, PartialEq)]
pub struct PatchField2 {
    /// Cell boundaries along `u`, increasing.
    pub u_breaks: Vec<Scalar>,
    /// Cell boundaries along `v`, increasing.
    pub v_breaks: Vec<Scalar>,
    /// Degree along `u`.
    pub u_degree: usize,
    /// Degree along `v`.
    pub v_degree: usize,
    /// Per cell (index `i * (v_breaks.len() - 1) + j` for cell `i` along `u`
    /// and `j` along `v`), the coefficients, index `a * (v_degree + 1) + b`.
    pub patches: Vec<Vec<Scalar>>,
}

/// Bernstein basis values and first two derivatives of degree `n` at `s`.
fn bernstein(n: usize, s: Scalar) -> [Vec<Scalar>; 3] {
    let basis = |n: usize| -> Vec<Scalar> {
        // de Casteljau-style build-up, stable on [0, 1].
        let mut b = vec![0.0; n + 1];
        b[0] = 1.0;
        for k in 1..=n {
            let mut prev = 0.0;
            for slot in b.iter_mut().take(k + 1) {
                let here = *slot;
                *slot = here * (1.0 - s) + prev * s;
                prev = here;
            }
        }
        b
    };
    let b0 = basis(n);
    let mut b1 = vec![0.0; n + 1];
    let mut b2 = vec![0.0; n + 1];
    if n >= 1 {
        let lower = basis(n - 1);
        for a in 0..=n {
            let left = if a >= 1 { lower[a - 1] } else { 0.0 };
            let right = if a < n { lower[a] } else { 0.0 };
            b1[a] = n as Scalar * (left - right);
        }
    }
    if n >= 2 {
        let lower = basis(n - 2);
        let at = |k: isize| {
            if k >= 0 && (k as usize) <= n - 2 {
                lower[k as usize]
            } else {
                0.0
            }
        };
        for a in 0..=n {
            let a = a as isize;
            b2[a as usize] = (n * (n - 1)) as Scalar * (at(a - 2) - 2.0 * at(a - 1) + at(a));
        }
    }
    [b0, b1, b2]
}

/// The Bernstein coefficients, over `[s0, s1]`, of the polynomial whose
/// coefficients over `[0, 1]` are `c` -- by de Casteljau, which holds for
/// parameters outside `[0, 1]` as well (extrapolation).
#[allow(clippy::needless_range_loop)] // de Casteljau's triangle, by index
fn restrict(c: &[Scalar], s0: Scalar, s1: Scalar) -> Vec<Scalar> {
    let n = c.len();
    // The polynomial's coefficients over [0, s]: the left points of de
    // Casteljau at s.
    let left = |c: &[Scalar], s: Scalar| -> Vec<Scalar> {
        let mut w = c.to_vec();
        let mut out = vec![0.0; n];
        out[0] = w[0];
        for k in 1..n {
            for i in 0..n - k {
                w[i] = w[i] * (1.0 - s) + w[i + 1] * s;
            }
            out[k] = w[0];
        }
        out
    };
    // Coefficients over [s, 1]: the right points of de Casteljau at s.
    let right = |c: &[Scalar], s: Scalar| -> Vec<Scalar> {
        let mut w = c.to_vec();
        let mut out = vec![0.0; n];
        out[n - 1] = w[n - 1];
        for k in 1..n {
            for i in 0..n - k {
                w[i] = w[i] * (1.0 - s) + w[i + 1] * s;
            }
            out[n - 1 - k] = w[n - 1 - k];
        }
        out
    };
    if s1 <= s0 {
        // A single parameter: every coefficient is the value there.
        let mut w = c.to_vec();
        for k in 1..n {
            for i in 0..n - k {
                w[i] = w[i] * (1.0 - s0) + w[i + 1] * s0;
            }
        }
        return vec![w[0]; n];
    }
    if s0 == 0.0 && s1 == 1.0 {
        return c.to_vec();
    }
    // Over [0, s1], then the part [s0 / s1, 1] of that; when s1 is zero or
    // close to it, over [s0, 1] first instead.
    if s1.abs() >= (1.0 - s0).abs() {
        let over = left(c, s1);
        right(&over, s0 / s1)
    } else {
        let over = right(c, s0);
        left(&over, (s1 - s0) / (1.0 - s0))
    }
}

impl PatchField2 {
    fn cells(&self) -> (usize, usize) {
        (self.u_breaks.len() - 1, self.v_breaks.len() - 1)
    }

    /// The cell index along one axis holding `x`, clamped.
    fn cell_of(breaks: &[Scalar], x: Scalar) -> usize {
        let n = breaks.len() - 1;
        let mut i = 0;
        while i + 1 < n && x >= breaks[i + 1] {
            i += 1;
        }
        i
    }

    /// Value, gradient and Hessian at `p`; beyond the grid, the nearest
    /// edge cell's polynomial continued.
    #[must_use]
    pub fn jet(&self, p: Point2) -> Jet2 {
        let (_, m) = self.cells();
        let i = Self::cell_of(&self.u_breaks, p.x);
        let j = Self::cell_of(&self.v_breaks, p.y);
        let (hu, hv) = (
            self.u_breaks[i + 1] - self.u_breaks[i],
            self.v_breaks[j + 1] - self.v_breaks[j],
        );
        // Past the grid's first or last cell the edge patch's polynomial
        // continues: a trace may look a little beyond a spline's domain.
        let s = (p.x - self.u_breaks[i]) / hu;
        let t = (p.y - self.v_breaks[j]) / hv;
        let bu = bernstein(self.u_degree, s);
        let bv = bernstein(self.v_degree, t);
        let c = &self.patches[i * m + j];
        let q = self.v_degree + 1;
        let mut jet = Jet2 {
            value: 0.0,
            gradient: Vec2::ZERO,
            uu: 0.0,
            uv: 0.0,
            vv: 0.0,
        };
        for a in 0..=self.u_degree {
            for b in 0..=self.v_degree {
                let k = c[a * q + b];
                jet.value += k * bu[0][a] * bv[0][b];
                jet.gradient.x += k * bu[1][a] * bv[0][b];
                jet.gradient.y += k * bu[0][a] * bv[1][b];
                jet.uu += k * bu[2][a] * bv[0][b];
                jet.uv += k * bu[1][a] * bv[1][b];
                jet.vv += k * bu[0][a] * bv[2][b];
            }
        }
        jet.gradient.x /= hu;
        jet.gradient.y /= hv;
        jet.uu /= hu * hu;
        jet.uv /= hu * hv;
        jet.vv /= hv * hv;
        jet
    }

    /// The largest coefficient: a bound on the field's size (convex hull).
    #[must_use]
    pub fn magnitude(&self) -> Scalar {
        self.patches
            .iter()
            .flatten()
            .fold(0.0, |m: Scalar, c| m.max(c.abs()))
    }

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.patches.iter().flatten().all(|c| c.is_finite())
            && self
                .u_breaks
                .iter()
                .chain(&self.v_breaks)
                .all(|b| b.is_finite())
    }

    /// A bound over the box: the hull of every overlapping cell's
    /// coefficients restricted to the box.
    #[must_use]
    pub fn bound(&self, cell: &Cell) -> Range {
        let (n, m) = self.cells();
        let q = self.v_degree + 1;
        let mut out: Option<Range> = None;
        for i in 0..n {
            let (a0, a1) = (self.u_breaks[i], self.u_breaks[i + 1]);
            // The first and last cells reach past the grid's ends.
            let reach_lo = if i == 0 { Scalar::NEG_INFINITY } else { a0 };
            let reach_hi = if i + 1 == n { Scalar::INFINITY } else { a1 };
            if reach_hi < cell.lo.x || reach_lo > cell.hi.x {
                continue;
            }
            let (s0, s1) = (
                (cell.lo.x.max(reach_lo) - a0) / (a1 - a0),
                (cell.hi.x.min(reach_hi) - a0) / (a1 - a0),
            );
            for j in 0..m {
                let (b0, b1) = (self.v_breaks[j], self.v_breaks[j + 1]);
                let reach_lo = if j == 0 { Scalar::NEG_INFINITY } else { b0 };
                let reach_hi = if j + 1 == m { Scalar::INFINITY } else { b1 };
                if reach_hi < cell.lo.y || reach_lo > cell.hi.y {
                    continue;
                }
                let (t0, t1) = (
                    (cell.lo.y.max(reach_lo) - b0) / (b1 - b0),
                    (cell.hi.y.min(reach_hi) - b0) / (b1 - b0),
                );
                let c = &self.patches[i * m + j];
                // Restrict every row along v, then every column along u.
                let mut rows: Vec<Vec<Scalar>> = (0..=self.u_degree)
                    .map(|a| restrict(&c[a * q..(a + 1) * q], t0, t1))
                    .collect();
                for b in 0..q {
                    let column: Vec<Scalar> = rows.iter().map(|r| r[b]).collect();
                    let restricted = restrict(&column, s0, s1);
                    for (a, value) in restricted.into_iter().enumerate() {
                        rows[a][b] = value;
                    }
                }
                let (mut lo, mut hi) = (Scalar::INFINITY, Scalar::NEG_INFINITY);
                let mut size: Scalar = 0.0;
                for value in rows.iter().flatten() {
                    lo = lo.min(*value);
                    hi = hi.max(*value);
                    size = size.max(value.abs());
                }
                let r = Range { lo, hi }
                    .widen(64.0 * Scalar::EPSILON * size * (q + self.u_degree + 1) as Scalar);
                out = Some(match out {
                    None => r,
                    Some(o) => Range {
                        lo: o.lo.min(r.lo),
                        hi: o.hi.max(r.hi),
                    },
                });
            }
        }
        out.unwrap_or(Range::point(0.0))
    }

    /// The partial derivative along `u` (`along_u`) or `v`, in the same
    /// cells, one degree lower along that axis.
    #[must_use]
    pub fn partial(&self, along_u: bool) -> Self {
        let (n, m) = self.cells();
        let (p, q) = (self.u_degree, self.v_degree);
        let mut patches = Vec::with_capacity(self.patches.len());
        for i in 0..n {
            for j in 0..m {
                let c = &self.patches[i * m + j];
                let at = |a: usize, b: usize| c[a * (q + 1) + b];
                if along_u {
                    let h = self.u_breaks[i + 1] - self.u_breaks[i];
                    if p == 0 {
                        patches.push(vec![0.0; q + 1]);
                        continue;
                    }
                    let mut d = vec![0.0; p * (q + 1)];
                    for a in 0..p {
                        for b in 0..=q {
                            d[a * (q + 1) + b] = p as Scalar * (at(a + 1, b) - at(a, b)) / h;
                        }
                    }
                    patches.push(d);
                } else {
                    let h = self.v_breaks[j + 1] - self.v_breaks[j];
                    if q == 0 {
                        patches.push(vec![0.0; p + 1]);
                        continue;
                    }
                    let mut d = vec![0.0; (p + 1) * q];
                    for a in 0..=p {
                        for b in 0..q {
                            d[a * q + b] = q as Scalar * (at(a, b + 1) - at(a, b)) / h;
                        }
                    }
                    patches.push(d);
                }
            }
        }
        let (u_degree, v_degree) = if along_u {
            (p.saturating_sub(1), q)
        } else {
            (p, q.saturating_sub(1))
        };
        Self {
            u_breaks: self.u_breaks.clone(),
            v_breaks: self.v_breaks.clone(),
            u_degree,
            v_degree,
            patches,
        }
    }
}

/// A field's value, gradient and Hessian at one point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Jet2 {
    /// The value.
    pub value: Scalar,
    /// `(dF/du, dF/dv)`.
    pub gradient: Vec2,
    /// `d2F/du2`.
    pub uu: Scalar,
    /// `d2F/dudv`.
    pub uv: Scalar,
    /// `d2F/dv2`.
    pub vv: Scalar,
}

impl SeriesField2 {
    fn size(&self) -> (usize, usize) {
        let n = self.coefficients.len();
        let m = self.coefficients.iter().map(Vec::len).max().unwrap_or(0);
        (n, m)
    }

    /// The field's value at `p`.
    #[must_use]
    pub fn value(&self, p: Point2) -> Scalar {
        // The same terms, summed in the same order, as `jet`'s value: the
        // two agree to the last bit. Small series stay on the stack.
        let (n, m) = self.size();
        let (mut su, mut sv) = ([0.0; 16], [0.0; 16]);
        let (mut hu, mut hv) = (Vec::new(), Vec::new());
        let fu: &mut [Scalar] = if n <= 16 {
            &mut su[..n]
        } else {
            hu.resize(n, 0.0);
            &mut hu
        };
        let fv: &mut [Scalar] = if m <= 16 {
            &mut sv[..m]
        } else {
            hv.resize(m, 0.0);
            &mut hv
        };
        self.u.values_into(p.x, fu);
        self.v.values_into(p.y, fv);
        let mut value = 0.0;
        for (i, row) in self.coefficients.iter().enumerate() {
            for (j, &c) in row.iter().enumerate() {
                if c == 0.0 {
                    continue;
                }
                value += c * fu[i] * fv[j];
            }
        }
        value
    }

    /// Value, gradient and Hessian at `p`.
    #[must_use]
    pub fn jet(&self, p: Point2) -> Jet2 {
        let (n, m) = self.size();
        let (fu, du, ddu) = self.u.terms(p.x, n);
        let (fv, dv, ddv) = self.v.terms(p.y, m);
        let mut jet = Jet2 {
            value: 0.0,
            gradient: Vec2::ZERO,
            uu: 0.0,
            uv: 0.0,
            vv: 0.0,
        };
        for (i, row) in self.coefficients.iter().enumerate() {
            for (j, &c) in row.iter().enumerate() {
                if c == 0.0 {
                    continue;
                }
                jet.value += c * fu[i] * fv[j];
                jet.gradient.x += c * du[i] * fv[j];
                jet.gradient.y += c * fu[i] * dv[j];
                jet.uu += c * ddu[i] * fv[j];
                jet.uv += c * du[i] * dv[j];
                jet.vv += c * fu[i] * ddv[j];
            }
        }
        jet
    }

    /// The sum of the coefficients' magnitudes: a bound on the field over
    /// the whole Fourier range, and the scale its rounding is measured in.
    #[must_use]
    pub fn magnitude(&self) -> Scalar {
        self.coefficients.iter().flatten().map(|c| c.abs()).sum()
    }

    /// Whether every coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.coefficients.iter().flatten().all(|c| c.is_finite())
    }
}

// --- Bounds over boxes ------------------------------------------------------
//
// Interval arithmetic over the terms (a power's or a harmonic's exact range
// over an interval), tightened by the mean-value form, and widened by a
// margin that covers the rounding of the sums.

/// A closed interval of reals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Range {
    /// Lower end.
    pub lo: Scalar,
    /// Upper end.
    pub hi: Scalar,
}

impl Range {
    /// The interval `[x, x]`.
    pub fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    /// The interval between `a` and `b`, in either order.
    pub fn new(a: Scalar, b: Scalar) -> Self {
        Self {
            lo: a.min(b),
            hi: a.max(b),
        }
    }

    fn add(self, o: Self) -> Self {
        Self {
            lo: self.lo + o.lo,
            hi: self.hi + o.hi,
        }
    }

    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self {
            lo: p.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            hi: p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        }
    }

    fn scale(self, c: Scalar) -> Self {
        Self::new(self.lo * c, self.hi * c)
    }

    fn intersect(self, o: Self) -> Self {
        Self {
            lo: self.lo.max(o.lo),
            hi: self.hi.min(o.hi),
        }
    }

    fn widen(self, by: Scalar) -> Self {
        Self {
            lo: self.lo - by,
            hi: self.hi + by,
        }
    }

    /// Whether zero lies in the interval.
    pub fn straddles_zero(self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }
}

/// The range of `cos(x)` over `[a, b]`.
fn cos_range(a: Scalar, b: Scalar) -> Range {
    if b - a >= TAU {
        return Range { lo: -1.0, hi: 1.0 };
    }
    let (ca, cb) = (a.cos(), b.cos());
    let mut r = Range::new(ca, cb);
    // A maximum at 2 pi m, a minimum at pi + 2 pi m.
    let m = (a / TAU).ceil();
    if m * TAU <= b {
        r.hi = 1.0;
    }
    let m = ((a - PI) / TAU).ceil();
    if PI + m * TAU <= b {
        r.lo = -1.0;
    }
    r
}

/// The range of basis term `k` over `[a, b]`.
fn term_range(basis: Basis, k: usize, a: Scalar, b: Scalar) -> Range {
    match basis {
        Basis::Power => {
            if k == 0 {
                return Range::point(1.0);
            }
            let (pa, pb) = (a.powi(k as i32), b.powi(k as i32));
            // Monotone unless an even power straddles zero.
            if k % 2 == 1 || a >= 0.0 || b <= 0.0 {
                Range::new(pa, pb)
            } else {
                Range {
                    lo: 0.0,
                    hi: pa.max(pb),
                }
            }
        }
        Basis::Fourier => {
            if k == 0 {
                return Range::point(1.0);
            }
            let w = k.div_ceil(2) as Scalar;
            if k % 2 == 1 {
                cos_range(w * a, w * b)
            } else {
                // sin(y) = cos(y - pi/2).
                cos_range(w * a - 0.5 * PI, w * b - 0.5 * PI)
            }
        }
    }
}

/// A box in the parameter plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cell {
    /// Lower corner.
    pub lo: Point2,
    /// Upper corner.
    pub hi: Point2,
}

impl Cell {
    /// The box's centre.
    pub fn centre(&self) -> Point2 {
        (self.lo + self.hi) * 0.5
    }
}

/// Naive interval bound of the field over the box.
fn naive(field: &SeriesField2, cell: &Cell) -> Range {
    let mut total = Range::point(0.0);
    let (n, m) = size(field);
    let ru: Vec<Range> = (0..n)
        .map(|k| term_range(field.u, k, cell.lo.x, cell.hi.x))
        .collect();
    let rv: Vec<Range> = (0..m)
        .map(|k| term_range(field.v, k, cell.lo.y, cell.hi.y))
        .collect();
    for (i, row) in field.coefficients.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            if c != 0.0 {
                total = total.add(ru[i].mul(rv[j]).scale(c));
            }
        }
    }
    total
}

/// The margin covering rounding in a sum over the field's terms at the
/// given parameter magnitudes.
fn margin(field: &SeriesField2, cell: &Cell) -> Scalar {
    let (n, m) = size(field);
    // Harmonics never exceed one; powers grow with the parameter.
    let reach = |basis: Basis, x: Scalar, terms: usize| match basis {
        Basis::Fourier => 1.0,
        Basis::Power => (1.0 + x.abs()).powi(terms.saturating_sub(1) as i32),
    };
    let scale = field.magnitude()
        * reach(field.u, cell.lo.x.abs().max(cell.hi.x.abs()), n)
        * reach(field.v, cell.lo.y.abs().max(cell.hi.y.abs()), m);
    64.0 * Scalar::EPSILON * scale
}

/// A bound certain to hold the field's values over the box (up to the
/// rounding margin included in it), given the field's partials.
pub fn bound(field: &Field2, du: &Field2, dv: &Field2, cell: &Cell) -> Range {
    match (field, du, dv) {
        (Field2::Series(f), Field2::Series(fu), Field2::Series(fv)) => {
            let direct = naive(f, cell);
            let c = cell.centre();
            let (hu, hv) = (0.5 * (cell.hi.x - cell.lo.x), 0.5 * (cell.hi.y - cell.lo.y));
            let mean = Range::point(f.value(c))
                .add(naive(fu, cell).mul(Range { lo: -hu, hi: hu }))
                .add(naive(fv, cell).mul(Range { lo: -hv, hi: hv }));
            direct.intersect(mean).widen(margin(f, cell))
        }
        _ => bound_simple(field, cell),
    }
}

/// A bound of the field over the box without the mean-value tightening.
pub fn bound_simple(field: &Field2, cell: &Cell) -> Range {
    match field {
        Field2::Series(f) => naive(f, cell).widen(margin(f, cell)),
        Field2::Patches(f) => f.bound(cell),
    }
}

/// The coefficient table's extent in `u` and `v`.
pub fn size(field: &SeriesField2) -> (usize, usize) {
    (
        field.coefficients.len(),
        field.coefficients.iter().map(Vec::len).max().unwrap_or(0),
    )
}

/// The partial derivative of a field along `u` (`along_u`) or `v`.
pub fn partial(field: &Field2, along_u: bool) -> Field2 {
    match field {
        Field2::Series(f) => Field2::Series(series_partial(f, along_u)),
        Field2::Patches(f) => Field2::Patches(f.partial(along_u)),
    }
}

fn series_partial(field: &SeriesField2, along_u: bool) -> SeriesField2 {
    let (n, m) = size(field);
    let mut out = vec![vec![0.0; m]; n];
    let basis = if along_u { field.u } else { field.v };
    let map = |k: usize| -> Option<(usize, Scalar)> {
        match basis {
            Basis::Power => (k > 0).then(|| (k - 1, k as Scalar)),
            Basis::Fourier => {
                if k == 0 {
                    None
                } else {
                    let w = k.div_ceil(2) as Scalar;
                    if k % 2 == 1 {
                        // cos(w x) -> -w sin(w x)
                        Some((k + 1, -w))
                    } else {
                        // sin(w x) -> w cos(w x)
                        Some((k - 1, w))
                    }
                }
            }
        }
    };
    for (i, row) in field.coefficients.iter().enumerate() {
        for (j, &c) in row.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            if along_u {
                if let Some((k, f)) = map(i) {
                    grow(&mut out, k, j);
                    out[k][j] += c * f;
                }
            } else if let Some((k, f)) = map(j) {
                grow(&mut out, i, k);
                out[i][k] += c * f;
            }
        }
    }
    SeriesField2 {
        u: field.u,
        v: field.v,
        coefficients: out,
    }
}

/// Grow a coefficient table to hold index `(i, j)`.
pub fn grow(c: &mut Vec<Vec<Scalar>>, i: usize, j: usize) {
    if c.len() <= i {
        let m = c.first().map_or(0, Vec::len);
        c.resize(i + 1, vec![0.0; m]);
    }
    if c[0].len() <= j {
        for row in c.iter_mut() {
            row.resize(j + 1, 0.0);
        }
    }
    for row in c.iter_mut() {
        if row.len() <= j {
            row.resize(j + 1, 0.0);
        }
    }
}

/// Which parameter a cell runs along; the other is solved for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// `u` runs freely; `v` is the field's zero.
    U,
    /// `v` runs freely; `u` is the field's zero.
    V,
}

/// One stretch of an [`ImplicitCurve2`]: as the free parameter runs from
/// `from` to `to`, the curve is the unique zero of the field for the other
/// parameter in `[low, high]`, where the field is strictly monotone in it.
///
/// A *bridge* is the last stretch into a point where two branches cross
/// (a saddle of the field on its zero set, where the surfaces touch).
/// Near there the zero cannot be isolated with certainty, so the bridge
/// carries the solved parameter as the cubic that matches the branch's
/// value and slope at both ends: at the certified end the field's own, at
/// the crossing the direction where the field's Hessian vanishes, which is
/// the branch's tangent there. It leaves the branch by about its length to
/// the fourth power (ADR 0077).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImplicitCell {
    /// The free parameter.
    pub axis: Axis,
    /// Where the free parameter starts.
    pub from: Scalar,
    /// Where it ends (it may run either way).
    pub to: Scalar,
    /// Lower end of the bracket of the solved parameter.
    pub low: Scalar,
    /// Upper end of the bracket.
    pub high: Scalar,
    /// For a bridge into a crossing, the slopes `d solved / d free` at
    /// `from` and at `to`.
    pub bridge: Option<(Scalar, Scalar)>,
}

impl ImplicitCell {
    /// A bridge from `from` to `to` in `(u, v)`, leaving along `start` and
    /// arriving along `end` (directions in `(u, v)`), free along the
    /// parameter both directions move in most.
    #[must_use]
    pub fn bridge(from: Point2, to: Point2, start: Vec2, end: Vec2) -> Self {
        let share = |d: Vec2, along_u: bool| {
            let l = d.length();
            if l == 0.0 {
                0.0
            } else if along_u {
                d.x.abs() / l
            } else {
                d.y.abs() / l
            }
        };
        let chord = to - from;
        let along_u = share(start, true)
            .min(share(end, true))
            .min(share(chord, true))
            >= share(start, false)
                .min(share(end, false))
                .min(share(chord, false));
        let slope = |d: Vec2| {
            let (free, solved) = if along_u { (d.x, d.y) } else { (d.y, d.x) };
            if free == 0.0 {
                0.0
            } else {
                solved / free
            }
        };
        let (axis, f0, f1, w0, w1) = if along_u {
            (Axis::U, from.x, to.x, from.y, to.y)
        } else {
            (Axis::V, from.y, to.y, from.x, to.x)
        };
        Self {
            axis,
            from: f0,
            to: f1,
            low: w0,
            high: w1,
            bridge: Some((slope(start), slope(end))),
        }
    }

    /// A bridge's solved value and its first and second derivatives in the
    /// local parameter `s`.
    fn hermite(&self, s: Scalar) -> (Scalar, Scalar, Scalar) {
        let (m0, m1) = self.bridge.unwrap_or((0.0, 0.0));
        let span = self.to - self.from;
        let (a, b) = (m0 * span, m1 * span);
        let (w0, w1) = (self.low, self.high);
        let (s2, s3) = (s * s, s * s * s);
        let value = (2.0 * s3 - 3.0 * s2 + 1.0) * w0
            + (s3 - 2.0 * s2 + s) * a
            + (-2.0 * s3 + 3.0 * s2) * w1
            + (s3 - s2) * b;
        let first = (6.0 * s2 - 6.0 * s) * w0
            + (3.0 * s2 - 4.0 * s + 1.0) * a
            + (-6.0 * s2 + 6.0 * s) * w1
            + (3.0 * s2 - 2.0 * s) * b;
        let second = (12.0 * s - 6.0) * w0
            + (6.0 * s - 4.0) * a
            + (-12.0 * s + 6.0) * w1
            + (6.0 * s - 2.0) * b;
        (value, first, second)
    }

    /// The local parameter of a free value.
    fn local(&self, free: Scalar) -> Scalar {
        let span = self.to - self.from;
        if span == 0.0 {
            0.0
        } else {
            (free - self.from) / span
        }
    }

    /// The part of the cell over local parameters `[s0, s1]`.
    #[must_use]
    pub fn part(&self, s0: Scalar, s1: Scalar) -> Self {
        let (f0, f1) = (self.free(s0), self.free(s1));
        match self.bridge {
            Some(_) => {
                let span = self.to - self.from;
                let slope = |s: Scalar| {
                    if span == 0.0 {
                        0.0
                    } else {
                        self.hermite(s).1 / span
                    }
                };
                Self {
                    from: f0,
                    to: f1,
                    low: self.hermite(s0).0,
                    high: self.hermite(s1).0,
                    bridge: Some((slope(s0), slope(s1))),
                    ..*self
                }
            }
            None => Self {
                from: f0,
                to: f1,
                ..*self
            },
        }
    }

    /// The cell run backwards.
    #[must_use]
    pub fn reversed(&self) -> Self {
        match self.bridge {
            Some((m0, m1)) => Self {
                from: self.to,
                to: self.from,
                low: self.high,
                high: self.low,
                bridge: Some((m1, m0)),
                ..*self
            },
            None => Self {
                from: self.to,
                to: self.from,
                ..*self
            },
        }
    }

    /// The smallest and largest solved value the cell can take: a
    /// bridge's range is bounded by its Bezier control values.
    #[must_use]
    pub fn solved_range(&self) -> (Scalar, Scalar) {
        match self.bridge {
            Some((m0, m1)) => {
                let span = self.to - self.from;
                let c = [
                    self.low,
                    self.low + m0 * span / 3.0,
                    self.high - m1 * span / 3.0,
                    self.high,
                ];
                (
                    c.iter().copied().fold(Scalar::INFINITY, Scalar::min),
                    c.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
                )
            }
            None => (self.low.min(self.high), self.low.max(self.high)),
        }
    }

    /// `(u, v)` from the free and solved values.
    fn place(&self, free: Scalar, solved: Scalar) -> Point2 {
        match self.axis {
            Axis::U => Point2::new(free, solved),
            Axis::V => Point2::new(solved, free),
        }
    }

    /// The free value at local parameter `s` in `[0, 1]`.
    fn free(&self, s: Scalar) -> Scalar {
        self.from + (self.to - self.from) * s
    }
}

/// A stretch of a field's zero set, in the parameters of the surface the
/// field lives on. The parameter `t` runs over `[0, cells.len()]`: cell `i`
/// covers `[i, i + 1]`, its free parameter moving linearly from `from` to
/// `to`.
#[derive(Debug, Clone, PartialEq)]
pub struct ImplicitCurve2 {
    /// The field whose zero the curve is.
    pub field: Field2,
    /// The cells, each starting where the previous ends.
    pub cells: Vec<ImplicitCell>,
}

impl ImplicitCurve2 {
    /// The parameter range, `[0, cells.len()]`.
    #[must_use]
    pub fn end(&self) -> Scalar {
        self.cells.len() as Scalar
    }

    /// The solved value of `cell` at free value `free`.
    #[must_use]
    pub fn solve_cell(&self, cell: &ImplicitCell, free: Scalar) -> Option<Scalar> {
        self.solve(cell, free)
    }

    /// The cell holding `t` and the local parameter in it.
    fn locate(&self, t: Scalar) -> Option<(&ImplicitCell, Scalar)> {
        if !t.is_finite() || self.cells.is_empty() {
            return None;
        }
        let last = self.cells.len() - 1;
        let slack = 1e-12 * (1.0 + self.end());
        if t < -slack || t > self.end() + slack {
            return None;
        }
        let index = (t.floor().max(0.0) as usize).min(last);
        let s = (t - index as Scalar).clamp(0.0, 1.0);
        Some((&self.cells[index], s))
    }

    /// The solved value in `cell` at free value `free`: the field's unique
    /// zero in the bracket.
    fn solve(&self, cell: &ImplicitCell, free: Scalar) -> Option<Scalar> {
        if cell.bridge.is_some() {
            return Some(cell.hermite(cell.local(free)).0);
        }
        let at = |w: Scalar| self.field.jet(cell.place(free, w));
        let along = |jet: &Jet2| match cell.axis {
            Axis::U => jet.gradient.y,
            Axis::V => jet.gradient.x,
        };
        let (mut lo, mut hi) = (cell.low, cell.high);
        let (f_lo, f_hi) = (at(lo).value, at(hi).value);
        if f_lo == 0.0 {
            return Some(lo);
        }
        if f_hi == 0.0 {
            return Some(hi);
        }
        // A zero that sits exactly on the bracket's end may round to the
        // wrong side there; within rounding of the field's scale the end is
        // the root.
        let scale = 1e-13 * self.field.magnitude().max(1.0);
        if f_lo.signum() == f_hi.signum() {
            return if f_lo.abs() <= scale && f_lo.abs() <= f_hi.abs() {
                Some(lo)
            } else if f_hi.abs() <= scale {
                Some(hi)
            } else {
                None
            };
        }
        let rising = f_hi > f_lo;
        let mut w = 0.5 * (lo + hi);
        for _ in 0..200 {
            let jet = at(w);
            let f = jet.value;
            if f == 0.0 {
                return Some(w);
            }
            if (f > 0.0) == rising {
                hi = w;
            } else {
                lo = w;
            }
            let slope = along(&jet);
            let newton = w - f / slope;
            let next = if slope != 0.0 && newton > lo && newton < hi {
                newton
            } else {
                0.5 * (lo + hi)
            };
            if (next - w).abs() <= 4.0 * Scalar::EPSILON * (1.0 + w.abs())
                || hi - lo <= 4.0 * Scalar::EPSILON * (1.0 + w.abs())
            {
                return Some(next);
            }
            w = next;
        }
        Some(w)
    }

    /// The point at `t`, or `None` outside `[0, cells.len()]`.
    #[must_use]
    pub fn point(&self, t: Scalar) -> Option<Point2> {
        let (cell, s) = self.locate(t)?;
        let free = cell.free(s);
        Some(cell.place(free, self.solve(cell, free)?))
    }

    /// `(solved', solved'')` against the free parameter, by the implicit
    /// function theorem, and the cell's rate `d free / dt`.
    fn slopes(&self, t: Scalar) -> Option<(&ImplicitCell, Scalar, Scalar, Scalar)> {
        let (cell, s) = self.locate(t)?;
        if cell.bridge.is_some() {
            let span = cell.to - cell.from;
            if span == 0.0 {
                return Some((cell, 0.0, 0.0, span));
            }
            let (_, d1, d2) = cell.hermite(s);
            return Some((cell, d1 / span, d2 / (span * span), span));
        }
        let free = cell.free(s);
        let solved = self.solve(cell, free)?;
        let jet = self.field.jet(cell.place(free, solved));
        let (f_free, f_solved, f_ff, f_fs, f_ss) = match cell.axis {
            Axis::U => (jet.gradient.x, jet.gradient.y, jet.uu, jet.uv, jet.vv),
            Axis::V => (jet.gradient.y, jet.gradient.x, jet.vv, jet.uv, jet.uu),
        };
        if f_solved == 0.0 {
            return None;
        }
        let first = -f_free / f_solved;
        let second = -(f_ff + 2.0 * f_fs * first + f_ss * first * first) / f_solved;
        Some((cell, first, second, cell.to - cell.from))
    }

    /// `dP/dt`.
    #[must_use]
    pub fn derivative(&self, t: Scalar) -> Option<Vec2> {
        let (cell, first, _, rate) = self.slopes(t)?;
        let d = match cell.axis {
            Axis::U => Vec2::new(1.0, first),
            Axis::V => Vec2::new(first, 1.0),
        };
        Some(d * rate)
    }

    /// `d2P/dt2`.
    #[must_use]
    pub fn second_derivative(&self, t: Scalar) -> Option<Vec2> {
        let (cell, _, second, rate) = self.slopes(t)?;
        let d = match cell.axis {
            Axis::U => Vec2::new(0.0, second),
            Axis::V => Vec2::new(second, 0.0),
        };
        Some(d * (rate * rate))
    }

    /// The parameter of a point on the curve: the cell whose box holds it,
    /// and where its free value falls in the cell.
    #[must_use]
    pub fn parameter_of(&self, p: Point2) -> Option<Scalar> {
        let mut best: Option<(Scalar, Scalar)> = None;
        for (index, cell) in self.cells.iter().enumerate() {
            let (free, solved) = match cell.axis {
                Axis::U => (p.x, p.y),
                Axis::V => (p.y, p.x),
            };
            let (lo, hi) = (cell.from.min(cell.to), cell.from.max(cell.to));
            let slack = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
            if free < lo - slack || free > hi + slack {
                continue;
            }
            let span = cell.to - cell.from;
            let s = if span == 0.0 {
                0.0
            } else {
                ((free - cell.from) / span).clamp(0.0, 1.0)
            };
            let Some(on) = self.solve(cell, cell.free(s)) else {
                continue;
            };
            let miss = (on - solved).abs();
            if best.is_none_or(|(m, _)| miss < m) {
                best = Some((miss, index as Scalar + s));
            }
        }
        best.map(|(_, t)| t)
    }

    /// The same curve moved by whole periods `(du, dv)` in parameters.
    #[must_use]
    pub fn shifted(&self, du: Scalar, dv: Scalar) -> Self {
        let cells = self
            .cells
            .iter()
            .map(|cell| {
                let (df, ds) = match cell.axis {
                    Axis::U => (du, dv),
                    Axis::V => (dv, du),
                };
                ImplicitCell {
                    from: cell.from + df,
                    to: cell.to + df,
                    low: cell.low + ds,
                    high: cell.high + ds,
                    ..*cell
                }
            })
            .collect();
        Self {
            field: self.field.clone(),
            cells,
        }
    }

    /// The offset from the curve's start to its end in parameters when it
    /// closes up to whole periods of `2 pi` in the periodic parameters
    /// (`(0, 0)` for a loop that does not wind), or `None` when it is open.
    #[must_use]
    pub fn closure(&self, periodic_u: bool, periodic_v: bool) -> Option<Vec2> {
        let (a, b) = (self.point(0.0)?, self.point(self.end())?);
        let d = b - a;
        let snap = |x: Scalar, periodic: bool| {
            if periodic {
                (x / TAU).round() * TAU
            } else {
                0.0
            }
        };
        let offset = Vec2::new(snap(d.x, periodic_u), snap(d.y, periodic_v));
        let miss = (d - offset).length();
        (miss <= 1e-8 * (1.0 + a.x.abs().max(a.y.abs()))).then_some(offset)
    }

    /// The stretch from `t0` to `t1`, as a curve of its own over
    /// `[0, cells]`. With `t0 > t1` the stretch runs on past the end and
    /// round from the start, which needs the curve's `closure` offset.
    #[must_use]
    pub fn sub(&self, t0: Scalar, t1: Scalar, closure: Option<Vec2>) -> Option<Self> {
        let n = self.end();
        if t0 > t1 {
            let offset = closure?;
            // Either part may be empty when a cut sits at the loop's start.
            let first = self.sub(t0, n, None);
            let second = self
                .sub(0.0, t1, None)
                .map(|c| c.shifted(offset.x, offset.y));
            return match (first, second) {
                (Some(mut a), Some(b)) => {
                    a.cells.extend(b.cells);
                    Some(a)
                }
                (a, b) => a.or(b),
            };
        }
        let (t0, t1) = (t0.clamp(0.0, n), t1.clamp(0.0, n));
        let mut cells = Vec::new();
        for (index, cell) in self.cells.iter().enumerate() {
            let (c0, c1) = (index as Scalar, index as Scalar + 1.0);
            let (a, b) = (t0.max(c0), t1.min(c1));
            if b - a <= 1e-12 {
                continue;
            }
            cells.push(cell.part(a - c0, b - c0));
        }
        (!cells.is_empty()).then(|| Self {
            field: self.field.clone(),
            cells,
        })
    }

    /// A closed curve's whole loop, starting at `t` instead of at `0`.
    #[must_use]
    pub fn rotated(&self, t: Scalar, closure: Vec2) -> Option<Self> {
        let n = self.end();
        if t <= 1e-12 || t >= n - 1e-12 {
            return Some(self.clone());
        }
        let index = (t.floor() as usize).min(self.cells.len() - 1);
        let s = t - index as Scalar;
        let cell = self.cells[index];
        let shift = |c: ImplicitCell| {
            let (df, ds) = match c.axis {
                Axis::U => (closure.x, closure.y),
                Axis::V => (closure.y, closure.x),
            };
            ImplicitCell {
                from: c.from + df,
                to: c.to + df,
                low: c.low + ds,
                high: c.high + ds,
                ..c
            }
        };
        let mut cells = Vec::with_capacity(self.cells.len() + 1);
        if s < 1.0 - 1e-12 {
            cells.push(cell.part(s, 1.0));
        }
        cells.extend(self.cells[index + 1..].iter().copied());
        cells.extend(self.cells[..index].iter().copied().map(shift));
        if s > 1e-12 {
            cells.push(shift(cell.part(0.0, s)));
        }
        Some(Self {
            field: self.field.clone(),
            cells,
        })
    }

    /// The stretches of the curve inside the box `[lo, hi]`, each as a curve
    /// of its own. Crossings of the box's sides are found by a scan of each
    /// cell and bisection on the distance to the box.
    #[must_use]
    pub fn clipped(&self, lo: Point2, hi: Point2) -> Vec<Self> {
        let outside = |t: Scalar| -> Scalar {
            self.point(t).map_or(Scalar::INFINITY, |p| {
                (lo.x - p.x).max(p.x - hi.x).max(lo.y - p.y).max(p.y - hi.y)
            })
        };
        let n = self.end();
        let steps = 16 * self.cells.len().max(1);
        let mut out = Vec::new();
        let mut start: Option<Scalar> = None;
        let mut previous = (0.0, outside(0.0) <= 0.0);
        if previous.1 {
            start = Some(0.0);
        }
        for k in 1..=steps {
            let t = n * k as Scalar / steps as Scalar;
            let inside = outside(t) <= 0.0;
            if inside != previous.1 {
                // Bisect the change.
                let (mut a, mut b) = (previous.0, t);
                for _ in 0..80 {
                    let m = 0.5 * (a + b);
                    if (outside(m) <= 0.0) == previous.1 {
                        a = m;
                    } else {
                        b = m;
                    }
                }
                let cross = 0.5 * (a + b);
                if inside {
                    start = Some(cross);
                } else if let Some(s) = start.take() {
                    out.extend(self.sub(s, cross, None));
                }
            }
            previous = (t, inside);
        }
        if let Some(s) = start {
            out.extend(self.sub(s, n, None));
        }
        out
    }

    /// The same curve run backwards: `t` becomes `cells.len() - t`.
    #[must_use]
    pub fn reversed(&self) -> Self {
        Self {
            field: self.field.clone(),
            cells: self
                .cells
                .iter()
                .rev()
                .map(ImplicitCell::reversed)
                .collect(),
        }
    }

    /// Parameters in `(lo, hi)` where the curve may stop being monotone in
    /// `u` or `v`: every cell boundary, and inside each cell every point
    /// where the solved parameter turns (the field's partial along the free
    /// parameter changes sign on the curve). Between consecutive values the
    /// curve is monotone in both parameters.
    ///
    /// Each turning point is isolated with interval bounds of that partial
    /// over boxes that certainly hold the curve (the solved parameter moves
    /// at most `max |F_free| / min |F_solved|` per unit of the free one), so
    /// none is missed; a double root, where the sign does not change, is
    /// not a turn and is not reported.
    #[must_use]
    pub fn turning_points(&self, lo: Scalar, hi: Scalar) -> Vec<Scalar> {
        let (lo, hi) = (lo.min(hi), lo.max(hi));
        let du = partial(&self.field, true);
        let dv = partial(&self.field, false);
        let mut out = Vec::new();
        for (index, cell) in self.cells.iter().enumerate() {
            let (t0, t1) = (index as Scalar, index as Scalar + 1.0);
            if t0 > lo && t0 < hi {
                out.push(t0);
            }
            // A bridge is straight: it never turns.
            if t1 <= lo || t0 >= hi || cell.bridge.is_some() {
                continue;
            }
            let (d_free, d_solved) = match cell.axis {
                Axis::U => (&du, &dv),
                Axis::V => (&dv, &du),
            };
            let s0 = (lo - t0).max(0.0);
            let s1 = (hi - t0).min(1.0);
            self.turns_in(cell, d_free, d_solved, s0, s1, 0, &mut |s| out.push(t0 + s));
        }
        out.sort_by(Scalar::total_cmp);
        out.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
        out
    }

    /// The box `[free(s0), free(s1)] x [...]` certain to hold the cell's
    /// curve for `s` in `[s0, s1]`.
    fn hull(
        &self,
        cell: &ImplicitCell,
        d_free: &Field2,
        d_solved: &Field2,
        s0: Scalar,
        s1: Scalar,
    ) -> Option<Cell> {
        let (f0, f1) = (cell.free(s0), cell.free(s1));
        let w0 = self.solve(cell, f0)?;
        if cell.bridge.is_some() {
            // Its part over `[s0, s1]` lies in its control values' range.
            let (w_lo, w_hi) = cell.part(s0, s1).solved_range();
            let (a, b) = (cell.place(f0.min(f1), w_lo), cell.place(f0.max(f1), w_hi));
            return Some(Cell {
                lo: a.min(b),
                hi: a.max(b),
            });
        }
        let bracket = |f_lo: Scalar, f_hi: Scalar, w_lo: Scalar, w_hi: Scalar| {
            let (a, b) = match cell.axis {
                Axis::U => (Point2::new(f_lo, w_lo), Point2::new(f_hi, w_hi)),
                Axis::V => (Point2::new(w_lo, f_lo), Point2::new(w_hi, f_hi)),
            };
            Cell { lo: a, hi: b }
        };
        let (f_lo, f_hi) = (f0.min(f1), f0.max(f1));
        let whole = bracket(f_lo, f_hi, cell.low, cell.high);
        let free = bound_simple(d_free, &whole);
        let solved = bound_simple(d_solved, &whole);
        let floor = solved.lo.abs().min(solved.hi.abs());
        let reach = if solved.straddles_zero() || floor == 0.0 {
            Scalar::INFINITY
        } else {
            free.lo.abs().max(free.hi.abs()) / floor * (f_hi - f_lo)
        };
        let (w_lo, w_hi) = ((w0 - reach).max(cell.low), (w0 + reach).min(cell.high));
        Some(bracket(f_lo, f_hi, w_lo, w_hi))
    }

    #[allow(clippy::too_many_arguments)]
    fn turns_in(
        &self,
        cell: &ImplicitCell,
        d_free: &Field2,
        d_solved: &Field2,
        s0: Scalar,
        s1: Scalar,
        depth: u32,
        found: &mut dyn FnMut(Scalar),
    ) {
        let Some(hull) = self.hull(cell, d_free, d_solved, s0, s1) else {
            return;
        };
        if !bound_simple(d_free, &hull).straddles_zero() {
            return;
        }
        let g = |s: Scalar| -> Option<Scalar> {
            let free = cell.free(s);
            let solved = self.solve(cell, free)?;
            Some(d_free.value(cell.place(free, solved)))
        };
        if depth >= 40 || s1 - s0 <= 1e-12 {
            // Isolated to the last bits: a turn only where the sign changes.
            if let (Some(a), Some(b)) = (g(s0), g(s1)) {
                if (a < 0.0) != (b < 0.0) && a != 0.0 {
                    found(0.5 * (s0 + s1));
                }
            }
            return;
        }
        let m = 0.5 * (s0 + s1);
        self.turns_in(cell, d_free, d_solved, s0, m, depth + 1, found);
        self.turns_in(cell, d_free, d_solved, m, s1, depth + 1, found);
    }

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.field.is_finite()
            && self.cells.iter().all(|c| {
                c.from.is_finite() && c.to.is_finite() && c.low.is_finite() && c.high.is_finite()
            })
    }
}

/// The surface a traced curve lies on, as the curve needs it, in the same
/// parameterisation as the matching `axiolid_surface` family.
#[derive(Debug, Clone, PartialEq)]
pub enum Carrier {
    /// `O + u X + v Y`.
    Plane(Frame3),
    /// A cylinder, elliptical cylinder or cone.
    Ruled(RuledCarrier),
    /// `O + r cos v (cos u X + sin u Y) + r sin v Z`.
    Sphere {
        /// Centre and axes.
        frame: Frame3,
        /// Radius.
        radius: Scalar,
    },
    /// A torus.
    Torus(TorusCarrier),
    /// A B-spline surface.
    Spline(Box<crate::spline_surface::BSplineSurface>),
}

/// A point's partial derivatives on a [`Carrier`], to second order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceJet {
    /// The point.
    pub point: Point3,
    /// `dP/du`.
    pub u: Vec3,
    /// `dP/dv`.
    pub v: Vec3,
    /// `d2P/du2`.
    pub uu: Vec3,
    /// `d2P/dudv`.
    pub uv: Vec3,
    /// `d2P/dv2`.
    pub vv: Vec3,
}

impl Carrier {
    /// The point at `(u, v)` with its partials.
    #[must_use]
    pub fn jet(&self, u: Scalar, v: Scalar) -> SurfaceJet {
        let (su, cu) = u.sin_cos();
        match self {
            Carrier::Plane(f) => SurfaceJet {
                point: f.origin + f.x * u + f.y * v,
                u: f.x,
                v: f.y,
                uu: Vec3::ZERO,
                uv: Vec3::ZERO,
                vv: Vec3::ZERO,
            },
            Carrier::Ruled(k) => {
                let f = &k.frame;
                let (rx, ry) = (k.x_radius + k.slope * v, k.y_radius + k.slope * v);
                SurfaceJet {
                    point: f.origin + f.x * (rx * cu) + f.y * (ry * su) + f.z * v,
                    u: f.x * (-rx * su) + f.y * (ry * cu),
                    v: f.x * (k.slope * cu) + f.y * (k.slope * su) + f.z,
                    uu: f.x * (-rx * cu) + f.y * (-ry * su),
                    uv: f.x * (-k.slope * su) + f.y * (k.slope * cu),
                    vv: Vec3::ZERO,
                }
            }
            Carrier::Sphere {
                frame: f,
                radius: r,
            } => {
                let (sv, cv) = v.sin_cos();
                let ring = f.x * cu + f.y * su;
                let ring_u = f.x * (-su) + f.y * cu;
                SurfaceJet {
                    point: f.origin + ring * (r * cv) + f.z * (r * sv),
                    u: ring_u * (r * cv),
                    v: ring * (-r * sv) + f.z * (r * cv),
                    uu: ring * (-r * cv),
                    uv: ring_u * (-r * sv),
                    vv: ring * (-r * cv) + f.z * (-r * sv),
                }
            }
            Carrier::Torus(t) => {
                let f = &t.frame;
                let (sv, cv) = v.sin_cos();
                let r = t.minor_radius;
                let ring = t.major_radius + r * cv;
                let dir = f.x * cu + f.y * su;
                let dir_u = f.x * (-su) + f.y * cu;
                SurfaceJet {
                    point: f.origin + dir * ring + f.z * (r * sv),
                    u: dir_u * ring,
                    v: dir * (-r * sv) + f.z * (r * cv),
                    uu: dir * (-ring),
                    uv: dir_u * (-r * sv),
                    vv: dir * (-r * cv) + f.z * (-r * sv),
                }
            }
            Carrier::Spline(b) => b.jet(u, v).unwrap_or(SurfaceJet {
                point: Point3::splat(Scalar::NAN),
                u: Vec3::splat(Scalar::NAN),
                v: Vec3::splat(Scalar::NAN),
                uu: Vec3::splat(Scalar::NAN),
                uv: Vec3::splat(Scalar::NAN),
                vv: Vec3::splat(Scalar::NAN),
            }),
        }
    }

    /// Principal parameters of a point on the carrier: angles in
    /// `(-pi, pi]` (a sphere's latitude in `[-pi/2, pi/2]`); a caller
    /// reading a curve that runs past them adds whole turns. A B-spline
    /// carrier has no closed-form inverse: `NaN`, and the caller inverts
    /// the surface itself.
    #[must_use]
    pub fn parameters(&self, p: Point3) -> (Scalar, Scalar) {
        let local = |f: &Frame3| {
            let d = p - f.origin;
            (d.dot(f.x), d.dot(f.y), d.dot(f.z))
        };
        match self {
            Carrier::Plane(f) => {
                let (x, y, _) = local(f);
                (x, y)
            }
            Carrier::Ruled(k) => {
                let (x, y, z) = local(&k.frame);
                let (rx, ry) = (k.x_radius + k.slope * z, k.y_radius + k.slope * z);
                ((y / ry).atan2(x / rx), z)
            }
            Carrier::Sphere { frame, .. } => {
                let (x, y, z) = local(frame);
                (y.atan2(x), z.atan2(x.hypot(y)))
            }
            Carrier::Torus(t) => {
                let (x, y, z) = local(&t.frame);
                (y.atan2(x), z.atan2(x.hypot(y) - t.major_radius))
            }
            Carrier::Spline(_) => (Scalar::NAN, Scalar::NAN),
        }
    }

    /// Whether each parameter is an angle, periodic with `2 pi`.
    #[must_use]
    pub fn periodic(&self) -> (bool, bool) {
        match self {
            Carrier::Plane(_) | Carrier::Spline(_) => (false, false),
            Carrier::Ruled(_) | Carrier::Sphere { .. } => (true, false),
            Carrier::Torus(_) => (true, true),
        }
    }

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        let frame = |f: &Frame3| {
            f.origin.is_finite() && f.x.is_finite() && f.y.is_finite() && f.z.is_finite()
        };
        match self {
            Carrier::Plane(f) => frame(f),
            Carrier::Ruled(k) => k.is_finite(),
            Carrier::Sphere { frame: f, radius } => frame(f) && radius.is_finite(),
            Carrier::Torus(t) => {
                frame(&t.frame) && t.major_radius.is_finite() && t.minor_radius.is_finite()
            }
            Carrier::Spline(b) => {
                b.control_points.iter().flatten().all(|p| p.is_finite())
                    && b.weights
                        .as_ref()
                        .is_none_or(|w| w.iter().flatten().all(|x| x.is_finite()))
            }
        }
    }
}

/// An [`ImplicitCurve2`] on its carrier, in space: the point at `t` is the
/// carrier's point at `curve.point(t)`.
#[derive(Debug, Clone, PartialEq)]
pub struct ImplicitSection3 {
    /// The surface the curve lies on.
    pub carrier: Carrier,
    /// The curve in the carrier's parameters.
    pub curve: ImplicitCurve2,
}

impl ImplicitSection3 {
    /// The point at `t`.
    #[must_use]
    pub fn point(&self, t: Scalar) -> Option<Point3> {
        let p = self.curve.point(t)?;
        Some(self.carrier.jet(p.x, p.y).point)
    }

    /// `dP/dt = P_u u' + P_v v'`.
    #[must_use]
    pub fn tangent(&self, t: Scalar) -> Option<Vec3> {
        let p = self.curve.point(t)?;
        let d = self.curve.derivative(t)?;
        let jet = self.carrier.jet(p.x, p.y);
        Some(jet.u * d.x + jet.v * d.y)
    }

    /// `d2P/dt2 = P_uu u'^2 + 2 P_uv u' v' + P_vv v'^2 + P_u u'' + P_v v''`.
    #[must_use]
    pub fn bend(&self, t: Scalar) -> Option<Vec3> {
        let p = self.curve.point(t)?;
        let d = self.curve.derivative(t)?;
        let dd = self.curve.second_derivative(t)?;
        let jet = self.carrier.jet(p.x, p.y);
        Some(
            jet.uu * (d.x * d.x)
                + jet.uv * (2.0 * d.x * d.y)
                + jet.vv * (d.y * d.y)
                + jet.u * dd.x
                + jet.v * dd.y,
        )
    }

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.carrier.is_finite() && self.curve.is_finite()
    }
}

/// A curve in space read in an analytic surface's parameters: the pcurve
/// at `t` is the carrier's parameters of `curve`'s point at `t`, so it
/// shares the edge's parameter exactly. This is the pcurve, on the analytic
/// face, of a section that only the other face's surface can carry (a
/// B-spline's section, ADR 0077): the inverse is in closed form for planes,
/// ruled surfaces, spheres and tori.
///
/// Angles are defined up to whole turns; `guide` holds the parameters at
/// evenly spaced `t` over `[start, end]`, unwrapped along the curve, and a
/// point is read at the turn nearest the guide there. Evaluation lives in
/// `axiolid-evaluate`, which evaluates the space curve.
#[derive(Debug, Clone, PartialEq)]
pub struct LiftedCurve2 {
    /// The curve in space.
    pub curve: Box<crate::Curve3>,
    /// The surface it is read on.
    pub carrier: Carrier,
    /// Parameter range the guide covers.
    pub start: Scalar,
    /// End of that range.
    pub end: Scalar,
    /// Unwrapped parameters at evenly spaced `t` from `start` to `end`.
    pub guide: Vec<Point2>,
}

impl LiftedCurve2 {
    /// The guide's parameters at `t`, interpolated.
    #[must_use]
    pub fn guide_at(&self, t: Scalar) -> Option<Point2> {
        let n = self.guide.len();
        if n == 0 {
            return None;
        }
        if n == 1 || self.end == self.start {
            return Some(self.guide[0]);
        }
        let x = ((t - self.start) / (self.end - self.start) * (n - 1) as Scalar)
            .clamp(0.0, (n - 1) as Scalar);
        let i = (x.floor() as usize).min(n - 2);
        let f = x - i as Scalar;
        Some(self.guide[i] + (self.guide[i + 1] - self.guide[i]) * f)
    }

    /// The carrier's parameters of a space point, at the turns nearest the
    /// guide at `t`.
    #[must_use]
    pub fn unwrap_at(&self, t: Scalar, point: Point3) -> Option<Point2> {
        let (u, v) = self.carrier.parameters(point);
        if !u.is_finite() || !v.is_finite() {
            return None;
        }
        let near = self.guide_at(t)?;
        let (pu, pv) = self.carrier.periodic();
        let snap = |x: Scalar, g: Scalar, periodic: bool| {
            if periodic {
                x + ((g - x) / TAU).round() * TAU
            } else {
                x
            }
        };
        Some(Point2::new(snap(u, near.x, pu), snap(v, near.y, pv)))
    }

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.carrier.is_finite()
            && self.start.is_finite()
            && self.end.is_finite()
            && self.guide.iter().all(|p| p.is_finite())
    }
}
