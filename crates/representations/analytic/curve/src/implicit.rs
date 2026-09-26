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

use crate::quadric_section::RuledCarrier;
use crate::torus_section::TorusCarrier;

/// How a [`Field2`] varies along one parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Basis {
    /// Powers: term `k` is `x^k`.
    Power,
    /// Harmonics: term `0` is `1`, term `2k - 1` is `cos(k x)` and term
    /// `2k` is `sin(k x)`.
    Fourier,
}

impl Basis {
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
/// bases of [`Basis`] along each parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct Field2 {
    /// The basis along `u`.
    pub u: Basis,
    /// The basis along `v`.
    pub v: Basis,
    /// `coefficients[i][j]` multiplies term `i` in `u` and term `j` in `v`.
    pub coefficients: Vec<Vec<Scalar>>,
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

impl Field2 {
    fn size(&self) -> (usize, usize) {
        let n = self.coefficients.len();
        let m = self.coefficients.iter().map(Vec::len).max().unwrap_or(0);
        (n, m)
    }

    /// The field's value at `p`.
    #[must_use]
    pub fn value(&self, p: Point2) -> Scalar {
        self.jet(p).value
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
}

impl ImplicitCell {
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

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.field.is_finite()
            && self.cells.iter().all(|c| {
                c.from.is_finite() && c.to.is_finite() && c.low.is_finite() && c.high.is_finite()
            })
    }
}

/// An analytic surface, as a curve on it needs it, in the same
/// parameterisation as the matching `axiolid_surface` family.
#[derive(Debug, Clone, Copy, PartialEq)]
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
        }
    }

    /// Principal parameters of a point on the carrier: angles in
    /// `(-pi, pi]` (a sphere's latitude in `[-pi/2, pi/2]`); a caller
    /// reading a curve that runs past them adds whole turns.
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
        }
    }

    /// Whether each parameter is an angle, periodic with `2 pi`.
    #[must_use]
    pub fn periodic(&self) -> (bool, bool) {
        match self {
            Carrier::Plane(_) => (false, false),
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
