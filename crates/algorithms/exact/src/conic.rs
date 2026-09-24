//! Conics: exact line/conic and conic/conic intersection.
//!
//! A conic is `A x^2 + B xy + C y^2 + D x + E y + F = 0` with exact dyadic
//! coefficients; circles and ellipses built from `f64` data are exact.
//!
//! # Line meets conic
//!
//! Along `p(t) = from + t*d` the conic is a quadratic in `t`, so a hit is
//! a [`Root2`] and everything [`crate::root`] decides applies: exact
//! tangency, exact ordering along the line, across different conics.
//!
//! # Conic meets conic
//!
//! Up to four points whose coordinates are, in general, roots of a quartic
//! that radicals cannot express usefully. They are computed the way CGAL's
//! algebraic kernel does it:
//!
//! 1. Shear `u = x + k*y` for a small integer `k`, chosen so no two
//!    intersection points share `u` (finitely many `k` are bad).
//! 2. Eliminate `y`: the resultant `R(u)` has the points' `u` as roots.
//! 3. The common root in `y` is rational in `u`: `y = N(u) / D(u)` from
//!    the first subresultant, valid because the shear made `D(u) != 0` at
//!    every real root.
//! 4. Each point is a [`RealRoot`] `u0` plus those polynomials. Exact
//!    coordinates as [`RealRoot`]s follow from a second resultant, and the
//!    sign of any conic or line at the point is the sign of one polynomial
//!    at `u0`.
//!
//! Tangency is exact: in sheared coordinates the intersection multiplicity
//! at a point is the multiplicity of `u0` as a root of `R`.

use axiolid_core::Point2;
use axiolid_guarantees::Sign;
use num_bigint::BigInt;

use crate::arith::{sign_product, Arith};
use crate::certify::{require_finite, ExactError};
use crate::construct::{Branch, Line};
use crate::dyadic::Dyadic;
use crate::poly::{IntPoly, RealRoot};
use crate::root::Root2;

/// `A x^2 + B xy + C y^2 + D x + E y + F = 0`, coefficients exact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Conic {
    coeffs: [Dyadic; 6],
}

fn d(value: f64) -> Dyadic {
    Dyadic::from_f64(value)
}

impl Conic {
    /// From the six coefficients `[A, B, C, D, E, F]`.
    ///
    /// Refuses non-finite input, and a quadratic part that is identically
    /// zero (that is a line; use the line APIs).
    pub fn from_coefficients(coeffs: [f64; 6]) -> Result<Self, ExactError> {
        require_finite(&coeffs)?;
        Self::from_dyadic(coeffs.map(d))
    }

    fn from_dyadic(coeffs: [Dyadic; 6]) -> Result<Self, ExactError> {
        if coeffs[..3].iter().all(|c| c.sign() == Some(Sign::Zero)) {
            return Err(ExactError::DegenerateConic);
        }
        Ok(Self { coeffs })
    }

    /// The circle about `centre` with `radius > 0`.
    pub fn circle(centre: Point2, radius: f64) -> Result<Self, ExactError> {
        require_finite(&[centre.x, centre.y, radius])?;
        if radius <= 0.0 {
            return Err(ExactError::NegativeRadius);
        }
        let (cx, cy, r) = (d(centre.x), d(centre.y), d(radius));
        let two = d(2.0);
        Self::from_dyadic([
            d(1.0),
            Dyadic::zero(),
            d(1.0),
            two.mul(&cx).neg(),
            two.mul(&cy).neg(),
            cx.square().add(&cy.square()).sub(&r.square()),
        ])
    }

    /// The ellipse about `centre` with semi-axis `a` along `axis` and `b`
    /// across it. `axis` need not be unit length: the implicit form is
    /// scaled by `|axis|^2`, so no square root or division is needed.
    pub fn ellipse(centre: Point2, axis: Point2, a: f64, b: f64) -> Result<Self, ExactError> {
        require_finite(&[centre.x, centre.y, axis.x, axis.y, a, b])?;
        if a <= 0.0 || b <= 0.0 {
            return Err(ExactError::NegativeRadius);
        }
        if axis.x == 0.0 && axis.y == 0.0 {
            return Err(ExactError::DegenerateLine);
        }
        // b^2 ((p-c).u)^2 + a^2 ((p-c).u_perp)^2 = a^2 b^2 |u|^2,
        // u = (ux, uy), u_perp = (-uy, ux). With q = p - c:
        //   (q.u)^2    = ux^2 qx^2 + 2 ux uy qx qy + uy^2 qy^2
        //   (q.uperp)^2 = uy^2 qx^2 - 2 ux uy qx qy + ux^2 qy^2
        let (ux, uy) = (d(axis.x), d(axis.y));
        let (a2, b2) = (d(a).square(), d(b).square());
        let two = d(2.0);
        let qa = b2.mul(&ux.square()).add(&a2.mul(&uy.square()));
        let qb = two.mul(&ux).mul(&uy).mul(&b2.sub(&a2));
        let qc = b2.mul(&uy.square()).add(&a2.mul(&ux.square()));
        let norm2 = ux.square().add(&uy.square());
        let rhs = a2.mul(&b2).mul(&norm2);
        Self::from_centred(qa, qb, qc, rhs, centre)
    }

    /// `qa X^2 + qb X Y + qc Y^2 = rhs` in `X = x - cx`, `Y = y - cy`.
    fn from_centred(
        qa: Dyadic,
        qb: Dyadic,
        qc: Dyadic,
        rhs: Dyadic,
        centre: Point2,
    ) -> Result<Self, ExactError> {
        let (cx, cy) = (d(centre.x), d(centre.y));
        let two = d(2.0);
        let dd = two.mul(&qa).mul(&cx).add(&qb.mul(&cy)).neg();
        let ee = two.mul(&qc).mul(&cy).add(&qb.mul(&cx)).neg();
        let ff = qa
            .mul(&cx.square())
            .add(&qb.mul(&cx).mul(&cy))
            .add(&qc.mul(&cy.square()))
            .sub(&rhs);
        Self::from_dyadic([qa, qb, qc, dd, ee, ff])
    }

    /// The coefficients `[A, B, C, D, E, F]`.
    #[must_use]
    pub fn coefficients(&self) -> &[Dyadic; 6] {
        &self.coeffs
    }

    /// Exact value at an `f64` point.
    #[must_use]
    pub fn eval(&self, p: Point2) -> Dyadic {
        let (x, y) = (d(p.x), d(p.y));
        let [a, b, c, dd, e, f] = &self.coeffs;
        a.mul(&x.square())
            .add(&b.mul(&x).mul(&y))
            .add(&c.mul(&y.square()))
            .add(&dd.mul(&x))
            .add(&e.mul(&y))
            .add(f)
    }

    /// Integer coefficients with the same zero set.
    fn integer(&self) -> [BigInt; 6] {
        let poly = IntPoly::from_dyadic(&self.coeffs);
        let mut out: [BigInt; 6] = Default::default();
        // `from_dyadic` scales all six by one power of two; it trims
        // trailing zeros, which the padding restores.
        for (slot, c) in out.iter_mut().zip(poly.coeffs()) {
            *slot = c.clone();
        }
        out
    }
}

// ------------------------------------------------------------- line x conic

/// Where a line meets a conic.
#[derive(Debug, Clone, PartialEq)]
pub enum ConicLineHits {
    /// The line misses the conic.
    None,
    /// The line is tangent: one hit of multiplicity two.
    Tangent(ConicLineHit),
    /// Two distinct hits, first then second along the line.
    Secant(ConicLineHit, ConicLineHit),
    /// Exactly one hit because the line is parallel to an asymptote (or
    /// the axis of a parabola): the quadratic in `t` degenerates to linear.
    Single(ConicLineHit),
    /// The whole line lies on the conic (a degenerate conic).
    OnConic,
}

/// One exact hit of a line with a conic: `t = (-beta +- sqrt(disc)) / alpha`
/// on the line, or `t = -gamma / (2 beta)` in the linear case.
#[derive(Debug, Clone, PartialEq)]
pub struct ConicLineHit {
    line: Line,
    t: Root2<Dyadic>,
}

/// Quadratic `alpha t^2 + 2 beta t + gamma` of the conic along the line.
fn along(line: Line, conic: &Conic) -> (Dyadic, Dyadic, Dyadic) {
    let (fx, fy) = (d(line.from().x), d(line.from().y));
    let (dx, dy) = (d(line.to().x).sub(&fx), d(line.to().y).sub(&fy));
    let [a, b, c, dd, e, _] = &conic.coeffs;
    let two = d(2.0);
    let alpha = a
        .mul(&dx.square())
        .add(&b.mul(&dx).mul(&dy))
        .add(&c.mul(&dy.square()));
    // 2*beta = 2A fx dx + B (fx dy + fy dx) + 2C fy dy + D dx + E dy
    let two_beta = two
        .mul(a)
        .mul(&fx)
        .mul(&dx)
        .add(&b.mul(&fx.mul(&dy).add(&fy.mul(&dx))))
        .add(&two.mul(c).mul(&fy).mul(&dy))
        .add(&dd.mul(&dx))
        .add(&e.mul(&dy));
    let gamma = conic.eval(line.from());
    // Along the line the conic is alpha t^2 + two_beta t + gamma. Return
    // (2 alpha, two_beta, 2 gamma) = (al, be, ga): the roots are then
    // t = (-be +- sqrt(be^2 - al*ga)) / al with no halving anywhere.
    (two.mul(&alpha), two_beta, two.mul(&gamma))
}

/// Exact intersection of a line with a conic.
pub fn line_conic_hits(line: Line, conic: &Conic) -> Result<ConicLineHits, ExactError> {
    let (al, be, ga) = along(line, conic);
    let exact = |x: &Dyadic| x.sign().expect("exact");
    if exact(&al) == Sign::Zero {
        return Ok(match (exact(&be), exact(&ga)) {
            (Sign::Zero, Sign::Zero) => ConicLineHits::OnConic,
            (Sign::Zero, _) => ConicLineHits::None,
            // Linear: be t + gamma = 0, i.e. t = -ga / (2 be).
            _ => ConicLineHits::Single(ConicLineHit {
                line,
                t: Root2 {
                    a: ga.neg(),
                    b: Dyadic::zero(),
                    c: Dyadic::zero(),
                    d: d(2.0).mul(&be),
                },
            }),
        });
    }
    let disc = be.square().sub(&al.mul(&ga));
    let hit = |branch: Branch| ConicLineHit {
        line,
        t: Root2 {
            a: be.neg(),
            b: match branch {
                Branch::Minus => d(-1.0),
                Branch::Plus => d(1.0),
            },
            c: disc.clone(),
            d: al.clone(),
        },
    };
    Ok(match exact(&disc) {
        Sign::Negative => ConicLineHits::None,
        Sign::Zero => ConicLineHits::Tangent(hit(Branch::Minus)),
        _ => {
            // Order along the line: t = (-be -+ sqrt)/al, so the Minus
            // branch comes first only when al > 0.
            let (first, second) = if exact(&al) == Sign::Positive {
                (Branch::Minus, Branch::Plus)
            } else {
                (Branch::Plus, Branch::Minus)
            };
            ConicLineHits::Secant(hit(first), hit(second))
        }
    })
}

impl ConicLineHit {
    /// The line this hit lies on.
    #[must_use]
    pub const fn line(&self) -> Line {
        self.line
    }

    /// Exact sign of `t - value`.
    pub fn cmp_param(&self, value: f64) -> Result<Sign, ExactError> {
        require_finite(&[value])?;
        let point = Root2 {
            a: d(value),
            b: Dyadic::zero(),
            c: Dyadic::zero(),
            d: d(1.0),
        };
        self.t.cmp_sign(&point).ok_or(ExactError::Undefined)
    }

    /// Exact sign of `self - other` along their common line.
    pub fn compare_along(&self, other: &Self) -> Result<Sign, ExactError> {
        if self.line != other.line {
            return Err(ExactError::DifferentLines);
        }
        self.t.cmp_sign(&other.t).ok_or(ExactError::Undefined)
    }

    /// An approximate point, for output only.
    #[must_use]
    pub fn approx_point(&self) -> Point2 {
        let t = approx_root2(&self.t);
        let (from, to) = (self.line.from(), self.line.to());
        Point2::new(from.x + t * (to.x - from.x), from.y + t * (to.y - from.y))
    }
}

fn approx_root2(r: &Root2<Dyadic>) -> f64 {
    let root = r.c.to_f64().max(0.0).sqrt();
    (r.a.to_f64() + r.b.to_f64() * root) / r.d.to_f64()
}

// ------------------------------------------------------------ conic x conic

/// How two conics meet.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConicIntersection {
    /// Finitely many points (possibly none), in increasing order of the
    /// sheared coordinate; use [`ConicPoint::x`] to sort by `x`.
    Points(Vec<ConicPoint>),
    /// The conics share a component (a whole curve), so the intersection
    /// is not a finite point set.
    Overlapping,
}

/// One exact intersection point of two conics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConicPoint {
    u: RealRoot,
    shear: i64,
    /// y = y_num(u) / den(u), x = x_num(u) / den(u).
    y_num: IntPoly,
    x_num: IntPoly,
    den: IntPoly,
    tangent: bool,
}

// Small integer-polynomial helpers (coefficients lowest degree first).
fn ip(c: &[i64]) -> IntPoly {
    IntPoly::new(c.iter().map(|&v| BigInt::from(v)).collect())
}

fn pconst(c: &BigInt) -> IntPoly {
    IntPoly::new(vec![c.clone()])
}

fn padd(a: &IntPoly, b: &IntPoly) -> IntPoly {
    let n = a.coeffs().len().max(b.coeffs().len());
    let get = |p: &IntPoly, i: usize| p.coeffs().get(i).cloned().unwrap_or_default();
    IntPoly::new((0..n).map(|i| get(a, i) + get(b, i)).collect())
}

fn pneg(a: &IntPoly) -> IntPoly {
    IntPoly::new(a.coeffs().iter().map(|c| -c).collect())
}

fn psub(a: &IntPoly, b: &IntPoly) -> IntPoly {
    padd(a, &pneg(b))
}

fn pmul(a: &IntPoly, b: &IntPoly) -> IntPoly {
    if a.is_zero() || b.is_zero() {
        return IntPoly::new(vec![]);
    }
    let mut out = vec![BigInt::from(0); a.coeffs().len() + b.coeffs().len() - 1];
    for (i, x) in a.coeffs().iter().enumerate() {
        for (j, y) in b.coeffs().iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    IntPoly::new(out)
}

/// The conic in sheared coordinates `u = x + k y`, as a quadratic in `y`
/// with coefficients in `Z[u]`: `q2 y^2 + q1(u) y + q0(u)`.
fn sheared(c: &[BigInt; 6], k: i64) -> [IntPoly; 3] {
    let [a, b, cc, dd, e, f] = c;
    let k = BigInt::from(k);
    // x = u - k y:
    // A(u - ky)^2 + B(u - ky)y + C y^2 + D(u - ky) + E y + F
    // y^2: A k^2 - B k + C
    // y^1: (B - 2Ak) u + (E - Dk)
    // y^0: A u^2 + D u + F
    let q2 = pconst(&(a * &k * &k - b * &k + cc));
    let q1 = IntPoly::new(vec![e - dd * &k, b - BigInt::from(2) * a * &k]);
    let q0 = IntPoly::new(vec![f.clone(), dd.clone(), a.clone()]);
    [q2, q1, q0]
}

/// Integer determinant by Bareiss fraction-free elimination.
fn determinant(mut m: Vec<Vec<BigInt>>) -> BigInt {
    let n = m.len();
    let mut sign = BigInt::from(1);
    let mut prev = BigInt::from(1);
    for k in 0..n {
        if m[k][k].sign() == num_bigint::Sign::NoSign {
            let Some(swap) = (k + 1..n).find(|&r| m[r][k].sign() != num_bigint::Sign::NoSign)
            else {
                return BigInt::from(0);
            };
            m.swap(k, swap);
            sign = -sign;
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let v = &m[i][j] * &m[k][k] - &m[i][k] * &m[k][j];
                m[i][j] = v / &prev;
            }
        }
        prev = m[k][k].clone();
    }
    sign * &m[n - 1][n - 1]
}

/// Resultant of two integer polynomials with formal degrees `n`, `m`.
fn resultant(p: &IntPoly, n: usize, q: &IntPoly, m: usize) -> BigInt {
    let size = n + m;
    if size == 0 {
        return BigInt::from(1);
    }
    let coef = |poly: &IntPoly, i: usize| poly.coeffs().get(i).cloned().unwrap_or_default();
    let mut rows = Vec::with_capacity(size);
    for r in 0..m {
        let mut row = vec![BigInt::from(0); size];
        for i in 0..=n {
            row[r + i] = coef(p, n - i);
        }
        rows.push(row);
    }
    for r in 0..n {
        let mut row = vec![BigInt::from(0); size];
        for i in 0..=m {
            row[r + i] = coef(q, m - i);
        }
        rows.push(row);
    }
    determinant(rows)
}

/// `Res_u(r(u), den(u) * Y - num(u))` as a polynomial in `Y`, by
/// evaluation at `Y = 0..=deg r` and exact interpolation (scaled by
/// `deg r!`, which does not move roots).
fn eliminate(r: &IntPoly, num: &IntPoly, den: &IntPoly) -> IntPoly {
    let n = r.degree().unwrap_or(0);
    let m = num.degree().unwrap_or(0).max(den.degree().unwrap_or(0));
    let values: Vec<BigInt> = (0..=n as i64)
        .map(|y| {
            let line = psub(&pmul(den, &ip(&[y])), num);
            resultant(r, n, &line, m)
        })
        .collect();
    // n! * M(Y) = sum_i (-1)^(n-i) C(n,i) v_i prod_{j != i} (Y - j)
    let mut out = IntPoly::new(vec![]);
    for (i, v) in values.iter().enumerate() {
        let mut term = pconst(&(v * binomial(n, i)));
        if (n - i) % 2 == 1 {
            term = pneg(&term);
        }
        for j in 0..=n {
            if j != i {
                term = pmul(&term, &ip(&[-(j as i64), 1]));
            }
        }
        out = padd(&out, &term);
    }
    out
}

fn binomial(n: usize, k: usize) -> BigInt {
    let mut out = BigInt::from(1);
    for i in 0..k {
        out = out * BigInt::from(n - i) / BigInt::from(i + 1);
    }
    out
}

/// `sum c_i * p_i` for exact dyadic scalars, scaled to integers by one
/// positive power of two (roots and signs unchanged).
fn combine(terms: &[(Dyadic, &IntPoly)]) -> IntPoly {
    let len = terms
        .iter()
        .map(|(_, p)| p.coeffs().len())
        .max()
        .unwrap_or(0);
    let coeffs: Vec<Dyadic> = (0..len)
        .map(|i| {
            terms.iter().fold(Dyadic::zero(), |acc, (c, p)| {
                let pi = p.coeffs().get(i).cloned().unwrap_or_default();
                acc.add(&c.mul(&Dyadic::from_parts(pi, 0)))
            })
        })
        .collect();
    IntPoly::from_dyadic(&coeffs)
}

/// Exact intersection of two conics.
pub fn conic_intersections(first: &Conic, second: &Conic) -> Result<ConicIntersection, ExactError> {
    let (c1, c2) = (first.integer(), second.integer());
    for k in 0..=MAX_SHEAR {
        let [a2, a1, a0] = sheared(&c1, k);
        let [b2, b1, b0] = sheared(&c2, k);
        // Both must stay genuinely quadratic in y, or the resultant formula
        // below is not the resultant.
        if a2.is_zero() || b2.is_zero() {
            continue;
        }
        // y-resultant of two quadratics:
        // (a2 b0 - a0 b2)^2 - (a2 b1 - a1 b2)(a1 b0 - a0 b1)
        let n = psub(&pmul(&a2, &b0), &pmul(&a0, &b2));
        let den = psub(&pmul(&b2, &a1), &pmul(&a2, &b1));
        let res = psub(
            &pmul(&n, &n),
            &pmul(
                &psub(&pmul(&a2, &b1), &pmul(&a1, &b2)),
                &psub(&pmul(&a1, &b0), &pmul(&a0, &b1)),
            ),
        );
        if res.is_zero() {
            return Ok(ConicIntersection::Overlapping);
        }
        if res.degree() == Some(0) {
            return Ok(ConicIntersection::Points(Vec::new()));
        }
        // y = n(u) / den(u) needs den(u0) != 0 at every real root u0; a
        // shared real root means two points share u (or a vertical
        // tangent in these coordinates): try the next shear.
        let sf = res.square_free();
        let shared = sf.gcd(&den);
        if shared.degree().unwrap_or(0) >= 1 && !shared.real_roots().is_empty() {
            continue;
        }
        let doubled = res.gcd(&res.derivative());
        // Complex roots shared with den must go before coordinates are
        // eliminated: there num vanishes too, so den*Y - num would share
        // that root for every Y and the eliminant would be identically
        // zero. Real roots are unaffected (none is shared, checked above).
        let sf = if shared.degree().unwrap_or(0) >= 1 {
            sf.exact_div(&shared)
        } else {
            sf
        };
        let x_num = psub(&pmul(&ip(&[0, 1]), &den), &pmul(&ip(&[k]), &n));
        let points = sf
            .real_roots()
            .into_iter()
            .map(|u| {
                let tangent =
                    doubled.degree().unwrap_or(0) >= 1 && u.sign_of(&doubled) == Sign::Zero;
                ConicPoint {
                    u,
                    shear: k,
                    y_num: n.clone(),
                    x_num: x_num.clone(),
                    den: den.clone(),
                    tangent,
                }
            })
            .collect();
        return Ok(ConicIntersection::Points(points));
    }
    Err(ExactError::DegenerateConic)
}

/// Shears tried before giving up. Each pair of points and each direction
/// of asymptote rules out at most one or two, so 0..=16 always suffices
/// for two proper conics.
const MAX_SHEAR: i64 = 16;

impl ConicPoint {
    /// Whether the conics are tangent here (intersection multiplicity at
    /// least two), decided exactly.
    #[must_use]
    pub fn is_tangent(&self) -> bool {
        self.tangent
    }

    /// Exact sign of the conic `other` at this point: zero on it, and for
    /// a circle or ellipse negative inside, positive outside.
    #[must_use]
    pub fn sign_of_conic(&self, other: &Conic) -> Sign {
        // Multiply through by den^2 > 0: every term is a polynomial in u.
        let [a, b, c, dd, e, f] = other.integer();
        let (xn, yn, den) = (&self.x_num, &self.y_num, &self.den);
        let terms = [
            pmul(&pconst(&a), &pmul(xn, xn)),
            pmul(&pconst(&b), &pmul(xn, yn)),
            pmul(&pconst(&c), &pmul(yn, yn)),
            pmul(&pconst(&dd), &pmul(xn, den)),
            pmul(&pconst(&e), &pmul(yn, den)),
            pmul(&pconst(&f), &pmul(den, den)),
        ];
        let total = terms
            .iter()
            .fold(IntPoly::new(vec![]), |acc, t| padd(&acc, t));
        if total.is_zero() {
            return Sign::Zero;
        }
        self.u.sign_of(&total)
    }

    /// Exact side of the directed line `a -> b`: positive to the left.
    #[must_use]
    pub fn side_of_line(&self, a: Point2, b: Point2) -> Sign {
        // Times den(u0):  dx*(y_num - ay*den) - dy*(x_num - ax*den),
        // assembled with exact dyadic coefficients so the scaling to
        // integers is one common positive factor.
        let (dx, dy) = (d(b.x).sub(&d(a.x)), d(b.y).sub(&d(a.y)));
        let shift = dx.mul(&d(a.y)).sub(&dy.mul(&d(a.x)));
        let value = combine(&[
            (dx, &self.y_num),
            (dy.neg(), &self.x_num),
            (shift.neg(), &self.den),
        ]);
        if value.is_zero() {
            return Sign::Zero;
        }
        sign_product(self.u.sign_of(&value), self.u.sign_of(&self.den))
    }

    /// The exact `x` coordinate, as a root of an integer polynomial.
    #[must_use]
    pub fn x(&self) -> RealRoot {
        self.coordinate(&self.x_num)
    }

    /// The exact `y` coordinate, as a root of an integer polynomial.
    #[must_use]
    pub fn y(&self) -> RealRoot {
        self.coordinate(&self.y_num)
    }

    /// The root `num(u0) / den(u0)` among the roots of the eliminant.
    fn coordinate(&self, num: &IntPoly) -> RealRoot {
        let eliminant = eliminate(self.u.poly(), num, &self.den);
        let sd = self.u.sign_of(&self.den);
        // value - y0 has the sign of (num - y0 den)(u0) * sign(den(u0)).
        let side = |y0: &Dyadic| {
            let shifted = IntPoly::from_dyadic(
                &(0..num.coeffs().len().max(self.den.coeffs().len()))
                    .map(|i| {
                        let n = num.coeffs().get(i).cloned().unwrap_or_default();
                        let dd = self.den.coeffs().get(i).cloned().unwrap_or_default();
                        Dyadic::from_parts(n, 0).sub(&y0.mul(&Dyadic::from_parts(dd, 0)))
                    })
                    .collect::<Vec<_>>(),
            );
            if shifted.is_zero() {
                return Sign::Zero;
            }
            sign_product(self.u.sign_of(&shifted), sd)
        };
        eliminant
            .real_roots()
            .into_iter()
            .find(|candidate| {
                let (lo, hi) = candidate.bounds();
                if candidate.is_exact() {
                    side(lo) == Sign::Zero
                } else {
                    side(lo) == Sign::Positive && side(hi) == Sign::Negative
                }
            })
            .expect("the coordinate is a root of its eliminant")
    }

    /// An approximate point, for output only.
    #[must_use]
    pub fn approx(&self) -> Point2 {
        Point2::new(self.x().approx(), self.y().approx())
    }

    /// The shear used internally (diagnostics).
    #[must_use]
    pub fn shear(&self) -> i64 {
        self.shear
    }
}
