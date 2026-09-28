//! Certified signs of series fields (ADR 0077, #181), where `f64` cannot
//! decide.
//!
//! A [`SeriesField2`] is `sum c[a][b] B_a(u) B_b(v)` with powers or
//! harmonics `cos(k x)`, `sin(k x)` along each parameter. Its derivatives
//! are series of the same terms, each an integer times a basis value, so
//! every question below is asked of the coefficients directly.
//!
//! - **Zero at a point.** At `f64` parameters every harmonic angle is an
//!   integer multiple of one dyadic `w`, so with `z = e^(i w)` the value is
//!   a Laurent polynomial in `z` with dyadic complex coefficients. `z` is
//!   transcendental (Lindemann: `e^a` is, for algebraic `a != 0`), so the
//!   value is zero exactly when every collected coefficient is: an exact
//!   test in dyadic arithmetic.
//! - **Sign at a point.** Otherwise the value is nonzero, and
//!   [`FixedInterval`] evaluation at rising precision shows its sign.
//! - **Sign over a box.** A Taylor form about a corner: the value there
//!   and the derivatives up to order `K - 1`, all certified, make a
//!   polynomial over the box, and the order-`K` derivatives' bounds over
//!   the box bound the remainder. The polynomial's Bernstein coefficients
//!   all clearing the remainder with one sign prove that sign over the
//!   box. `K` rises until they do, or cannot.
//!
//! The trace consults this tier only where the field is flat
//! ([`SeriesTier::flat`]): along a line of contact or at a singular point,
//! where `f64` bounds fail at every size. Elsewhere its answers would
//! change nothing but the cost.

use std::collections::BTreeMap;

use axiolid_core::{Point2, Scalar};
use axiolid_curve::implicit::{size, Basis, SeriesField2};
use axiolid_exact::{Arith, Dyadic, FixedInterval};
use axiolid_guarantees::Sign;

/// A series field's certified tier.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SeriesTier<'a> {
    field: &'a SeriesField2,
}

/// A basis value: what a term's derivative is, up to an integer factor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Term {
    One,
    Power(usize),
    Cos(usize),
    Sin(usize),
}

/// The `i`-th derivative of basis term `a` along one parameter: an integer
/// factor and a basis value, `None` for a zero derivative or a factor past
/// `i128`.
fn derivative(basis: Basis, a: usize, i: usize) -> Option<Option<(i128, Term)>> {
    match basis {
        Basis::Power => {
            if a < i {
                return Some(None);
            }
            let mut factor: i128 = 1;
            for t in (a - i + 1)..=a {
                factor = factor.checked_mul(t as i128)?;
            }
            let e = a - i;
            Some(Some((
                factor,
                if e == 0 { Term::One } else { Term::Power(e) },
            )))
        }
        Basis::Fourier => {
            if a == 0 {
                return Some((i == 0).then_some((1, Term::One)));
            }
            let k = a.div_ceil(2);
            let power = (k as i128).checked_pow(u32::try_from(i).ok()?)?;
            let (cos, sin) = (Term::Cos(k), Term::Sin(k));
            // cos -> -sin -> -cos -> sin; sin -> cos -> -sin -> -cos.
            let (sign, term) = if a % 2 == 1 {
                [(1, cos), (-1, sin), (-1, cos), (1, sin)][i % 4]
            } else {
                [(1, sin), (1, cos), (-1, sin), (-1, cos)][i % 4]
            };
            Some(Some((sign * power, term)))
        }
    }
}

/// Precisions tried for a point's sign, in fractional bits.
const PRECISIONS: [u32; 6] = [128, 256, 512, 1024, 2048, 4096];

/// The highest precision a box's corner value is taken to.
const MAX_BITS: u32 = 8192;

/// How small a gradient is against the sum of its terms' magnitudes for
/// the field to count as flat.
const FLAT: Scalar = 1e-4;

/// The highest Taylor order of a box's form.
const MAX_ORDER: usize = 32;

/// Basis values at one parameter, certified.
struct Table {
    powers: Vec<FixedInterval>,
    cos: Vec<FixedInterval>,
    sin: Vec<FixedInterval>,
}

impl Table {
    fn at(basis: Basis, x: Scalar, terms: usize, bits: u32) -> Option<Self> {
        let d = Dyadic::try_from_f64(x)?;
        let mut table = Table {
            powers: Vec::new(),
            cos: Vec::new(),
            sin: Vec::new(),
        };
        match basis {
            Basis::Power => {
                let mut p = Dyadic::from_f64(1.0);
                for _ in 0..terms.max(1) {
                    table.powers.push(FixedInterval::from_dyadic(&p, bits));
                    p = p.mul(&d);
                }
            }
            Basis::Fourier => {
                // cos((k + 1) x) = 2 cos(x) cos(k x) - cos((k - 1) x), and
                // likewise for sin: one certified pair, then exact integer
                // steps and rounded products, each outward.
                // The intervals widen by at most about `(1 + sqrt 2)^k`:
                // two guard bits per harmonic cover it.
                let n = terms.div_ceil(2);
                let w = bits + 2 * n as u32 + 8;
                let (s1, c1) = FixedInterval::sin_cos(&d, w)?;
                let twice = c1.mul_int(2);
                let mut sin = vec![FixedInterval::integer(0, w), s1];
                let mut cos = vec![FixedInterval::integer(1, w), c1];
                for k in 2..=n {
                    sin.push(twice.mul(&sin[k - 1]).sub(&sin[k - 2]));
                    cos.push(twice.mul(&cos[k - 1]).sub(&cos[k - 2]));
                }
                table.sin = sin.iter().map(|x| x.with_bits(bits)).collect();
                table.cos = cos.iter().map(|x| x.with_bits(bits)).collect();
            }
        }
        Some(table)
    }

    fn value(&self, term: Term, bits: u32) -> FixedInterval {
        match term {
            Term::One => FixedInterval::integer(1, bits),
            Term::Power(e) => self.powers[e].clone(),
            Term::Cos(k) => self.cos[k].clone(),
            Term::Sin(k) => self.sin[k].clone(),
        }
    }
}

/// A basis term's `f64` value at `x`.
fn term_f64(term: Term, x: Scalar) -> Scalar {
    match term {
        Term::One => 1.0,
        Term::Power(e) => x.powi(e as i32),
        Term::Cos(k) => (k as Scalar * x).cos(),
        Term::Sin(k) => (k as Scalar * x).sin(),
    }
}

/// A bound on a basis term's magnitude over `[lo, hi]`.
fn term_bound(term: Term, lo: Scalar, hi: Scalar) -> Scalar {
    match term {
        Term::Power(e) => lo.abs().max(hi.abs()).powi(e as i32),
        _ => 1.0,
    }
}

/// An exact dyadic integer.
fn dyadic_int(n: i128) -> Dyadic {
    // In 32-bit chunks, each exact as an `f64`.
    let negative = n < 0;
    let mut m = n.unsigned_abs();
    let mut out = Dyadic::from_f64(0.0);
    let mut scale = Dyadic::from_f64(1.0);
    let step = Dyadic::from_f64(4_294_967_296.0);
    while m > 0 {
        let chunk = (m & 0xffff_ffff) as Scalar;
        out = out.add(&Dyadic::from_f64(chunk).mul(&scale));
        scale = scale.mul(&step);
        m >>= 32;
    }
    if negative {
        out.neg()
    } else {
        out
    }
}

/// `x = m * 2^e` with an integer `m`.
fn parts(x: Scalar) -> (i128, i32) {
    if x == 0.0 {
        return (0, 0);
    }
    let bits = x.to_bits();
    let biased = ((bits >> 52) & 0x7ff) as i32;
    let fraction = (bits & ((1u64 << 52) - 1)) as i128;
    let (m, e) = if biased == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1i128 << 52), biased - 1075)
    };
    (if x < 0.0 { -m } else { m }, e)
}

/// A complex dyadic coefficient of a Laurent polynomial.
type Complex = (Dyadic, Dyadic);

impl<'a> SeriesTier<'a> {
    pub(crate) fn new(field: &'a SeriesField2) -> Self {
        Self { field }
    }

    /// The `(i, j)` derivative's terms: coefficient, factors and values.
    fn terms(&self, d: (usize, usize)) -> Option<Vec<(Scalar, i128, Term, Term)>> {
        let f = self.field;
        let (_, m) = size(f);
        let mut along_v = Vec::with_capacity(m);
        for b in 0..m {
            along_v.push(derivative(f.v, b, d.1)?);
        }
        let mut out = Vec::new();
        for (a, row) in f.coefficients.iter().enumerate() {
            if row.iter().all(|&c| c == 0.0) {
                continue;
            }
            let Some(u) = derivative(f.u, a, d.0)? else {
                continue;
            };
            for (&c, v) in row.iter().zip(&along_v) {
                if let (true, Some(v)) = (c != 0.0, v) {
                    out.push((c, u.0.checked_mul(v.0)?, u.1, v.1));
                }
            }
        }
        Some(out)
    }

    fn tables(&self, p: Point2, bits: u32) -> Option<(Table, Table)> {
        let (n, m) = size(self.field);
        Some((
            Table::at(self.field.u, p.x, n, bits)?,
            Table::at(self.field.v, p.y, m, bits)?,
        ))
    }

    /// The derivative `d` at the tables' point, certified.
    fn certified(&self, t: &(Table, Table), d: (usize, usize), bits: u32) -> Option<FixedInterval> {
        let mut sum = FixedInterval::integer(0, bits);
        for (c, factor, tu, tv) in self.terms(d)? {
            let term = FixedInterval::from_f64(c, bits)?
                .mul(&t.0.value(tu, bits))
                .mul(&t.1.value(tv, bits))
                .mul_int(factor);
            sum = sum.add(&term);
        }
        Some(sum)
    }

    /// The derivative `d` at `p` in `f64`, and the scale of its rounding.
    fn approximate(&self, p: Point2, d: (usize, usize)) -> Option<(Scalar, Scalar)> {
        let (mut value, mut scale) = (0.0, 0.0);
        for (c, factor, tu, tv) in self.terms(d)? {
            let term = c * factor as Scalar * term_f64(tu, p.x) * term_f64(tv, p.y);
            value += term;
            scale += term.abs();
        }
        Some((value, 64.0 * Scalar::EPSILON * scale))
    }

    /// A bound on the derivative `d`'s magnitude over the box, `None` where
    /// it vanishes identically.
    fn bound(&self, lo: Point2, hi: Point2, d: (usize, usize)) -> Option<Option<Scalar>> {
        let terms = self.terms(d)?;
        if terms.is_empty() {
            return Some(None);
        }
        let mut sum = 0.0;
        for (c, factor, tu, tv) in terms {
            sum += (c * factor as Scalar).abs()
                * term_bound(tu, lo.x, hi.x)
                * term_bound(tv, lo.y, hi.y);
        }
        // Past the sum's own rounding, and never below the smallest normal
        // (a bound underflowing to zero would claim too much).
        Some(Some(sum * (1.0 + 1e-10) + Scalar::MIN_POSITIVE))
    }

    /// Whether the derivative `d` is exactly zero at `p`; `None` where the
    /// harmonics' common angle does not fit the test.
    fn vanishes(&self, p: Point2, d: (usize, usize)) -> Option<bool> {
        let f = self.field;
        // The common angle `w = 2^e` of the harmonics: each parameter
        // along a Fourier basis is an integer multiple of it.
        let (mu, eu) = parts(p.x);
        let (mv, ev) = parts(p.y);
        let mut e = i32::MAX;
        if f.u == Basis::Fourier && mu != 0 {
            e = e.min(eu);
        }
        if f.v == Basis::Fourier && mv != 0 {
            e = e.min(ev);
        }
        let multiple = |basis: Basis, m: i128, ex: i32| -> Option<i128> {
            if basis == Basis::Power || m == 0 {
                return Some(0);
            }
            let shift = u32::try_from(ex - e).ok()?;
            // Room for the harmonic's number too.
            (shift < 60).then_some(m << shift)
        };
        let (pu, pv) = (multiple(f.u, mu, eu)?, multiple(f.v, mv, ev)?);
        let half = Dyadic::from_f64(0.5);
        let zero = Dyadic::from_f64(0.0);
        // A basis value as `sum coefficient * z^exponent`.
        let laurent = |term: Term, x: Scalar, turns: i128| -> Option<Vec<(i128, Complex)>> {
            Some(match term {
                Term::One => vec![(0, (Dyadic::from_f64(1.0), zero.clone()))],
                Term::Power(k) => {
                    let x = Dyadic::from_f64(x);
                    let mut value = Dyadic::from_f64(1.0);
                    for _ in 0..k {
                        value = value.mul(&x);
                    }
                    vec![(0, (value, zero.clone()))]
                }
                Term::Cos(k) => {
                    let n = (k as i128).checked_mul(turns)?;
                    vec![
                        (n, (half.clone(), zero.clone())),
                        (-n, (half.clone(), zero.clone())),
                    ]
                }
                // sin = (z^n - z^-n) / 2i.
                Term::Sin(k) => {
                    let n = (k as i128).checked_mul(turns)?;
                    vec![
                        (n, (zero.clone(), half.neg())),
                        (-n, (zero.clone(), half.clone())),
                    ]
                }
            })
        };
        let mut sum: BTreeMap<i128, Complex> = BTreeMap::new();
        for (c, factor, tu, tv) in self.terms(d)? {
            let scale = Dyadic::from_f64(c).mul(&dyadic_int(factor));
            for (nu, (ar, ai)) in laurent(tu, p.x, pu)? {
                for (nv, (br, bi)) in laurent(tv, p.y, pv)? {
                    let re = ar.mul(&br).sub(&ai.mul(&bi)).mul(&scale);
                    let im = ar.mul(&bi).add(&ai.mul(&br)).mul(&scale);
                    let slot = sum
                        .entry(nu.checked_add(nv)?)
                        .or_insert_with(|| (zero.clone(), zero.clone()));
                    slot.0 = slot.0.add(&re);
                    slot.1 = slot.1.add(&im);
                }
            }
        }
        Some(
            sum.values()
                .all(|(re, im)| re.sign() == Some(Sign::Zero) && im.sign() == Some(Sign::Zero)),
        )
    }

    /// Whether the derivative `d` is flat at `p`: its gradient small
    /// against the sum of its terms' magnitudes, as along a line of
    /// contact or at a singular point, where `f64` bounds fail at every
    /// size. Elsewhere subdivision settles a sign as it always has.
    pub(crate) fn flat(&self, p: Point2, d: (usize, usize)) -> bool {
        [(1, 0), (0, 1)].iter().all(|step| {
            self.approximate(p, (d.0 + step.0, d.1 + step.1))
                .is_some_and(|(a, rounding)| a.abs() <= FLAT * rounding / (64.0 * Scalar::EPSILON))
        })
    }

    /// The sign of the derivative `d` at `p`: `Zero` only where it is
    /// exactly zero, `None` where no precision tried decides.
    pub(crate) fn sign_at(&self, p: Point2, d: (usize, usize)) -> Option<Sign> {
        if !p.is_finite() {
            return None;
        }
        if self.vanishes(p, d) == Some(true) {
            return Some(Sign::Zero);
        }
        for bits in PRECISIONS {
            let t = self.tables(p, bits)?;
            match self.certified(&t, d, bits)?.sign() {
                Some(Sign::Zero) | None => {}
                Some(s) => return Some(s),
            }
        }
        None
    }

    /// The strict sign the derivative `d` keeps over the box `[lo, hi]`,
    /// or `None` where it is not proven to keep one.
    pub(crate) fn keeps_sign(&self, lo: Point2, hi: Point2, d: (usize, usize)) -> Option<Sign> {
        if !(lo.is_finite() && hi.is_finite()) || lo.x > hi.x || lo.y > hi.y {
            return None;
        }
        if lo == hi {
            return self.sign_at(lo, d).filter(|s| *s != Sign::Zero);
        }
        // The Taylor form about the corner `lo`, over `[lo, lo + span]`
        // (holding the box): `span` rounded up.
        let span = (
            if hi.x > lo.x {
                (hi.x - lo.x).next_up()
            } else {
                0.0
            },
            if hi.y > lo.y {
                (hi.y - lo.y).next_up()
            } else {
                0.0
            },
        );
        let form = Form {
            tier: self,
            lo,
            hi,
            d,
            span,
        };
        // Screen in `f64` where the value clears its rounding: a regular
        // zero in the box (the common case) shows at once, uncertified.
        let (g, rounding) = self.approximate(lo, d)?;
        let hidden = g.abs() <= 64.0 * rounding;
        let approximate = |i: usize, j: usize| {
            let (a, r) = self.approximate(lo, (d.0 + i, d.1 + j))?;
            // A value hidden by rounding may be of either sign.
            if hidden && (i, j) == (0, 0) {
                let e = a.abs() + 64.0 * r;
                return Some(Bound { lo: -e, hi: e });
            }
            Some(Bound {
                lo: a - r,
                hi: a + r,
            })
        };
        match form.verdict(approximate)? {
            Verdict::Sign(_) => {}
            // Visible terms of both signs (about a saddle, say): no
            // precision helps.
            Verdict::Mixed => return None,
            Verdict::Unclear if !hidden => return None,
            Verdict::Unclear => {}
        }
        // Zero at the corner: no strict sign over the box.
        if self.vanishes(lo, d) != Some(false) {
            return None;
        }
        // The corner's value, to a precision well past its magnitude.
        let mut bits = PRECISIONS[0];
        let t = loop {
            let t = self.tables(lo, bits)?;
            let (lower, _) = self.certified(&t, d, bits)?.magnitude();
            if lower > 0.0 {
                let needed = (-lower.log2()).ceil().max(0.0) as u32 + 64;
                if needed <= bits {
                    break t;
                }
                bits = needed.next_multiple_of(64);
            } else {
                bits *= 2;
            }
            if bits > MAX_BITS {
                return None;
            }
        };
        let certified = |i: usize, j: usize| {
            let e = self.certified(&t, (d.0 + i, d.1 + j), bits)?.enclosure();
            Some(Bound {
                lo: e.lo(),
                hi: e.hi(),
            })
        };
        match form.verdict(certified)? {
            Verdict::Sign(s) => Some(s),
            Verdict::Mixed | Verdict::Unclear => None,
        }
    }
}

/// An `f64` interval with every operation rounded outward.
#[derive(Debug, Clone, Copy)]
struct Bound {
    lo: Scalar,
    hi: Scalar,
}

impl Bound {
    fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    /// `x`, known to a relative `1e-14`.
    fn near(x: Scalar) -> Self {
        let e = x.abs() * 1e-14;
        Self {
            lo: (x - e).next_down(),
            hi: (x + e).next_up(),
        }
    }

    fn add(self, o: Self) -> Self {
        Self {
            lo: (self.lo + o.lo).next_down(),
            hi: (self.hi + o.hi).next_up(),
        }
    }

    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        let lo = p.iter().copied().fold(Scalar::INFINITY, Scalar::min);
        let hi = p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max);
        Self {
            lo: lo.next_down(),
            hi: hi.next_up(),
        }
    }
}

/// What a Taylor form's Bernstein coefficients say.
enum Verdict {
    /// All strictly of one sign past the remainder: proven.
    Sign(Sign),
    /// Some certainly of each sign.
    Mixed,
    /// Neither, at every order tried.
    Unclear,
}

/// The Taylor form of a derivative about a box's corner.
struct Form<'t, 'a> {
    tier: &'t SeriesTier<'a>,
    lo: Point2,
    hi: Point2,
    d: (usize, usize),
    span: (Scalar, Scalar),
}

impl Form<'_, '_> {
    /// Raise the order until the form's Bernstein coefficients decide,
    /// with `coefficient(i, j)` the derivative `(i, j)` at the corner.
    ///
    /// Over `[lo, lo + span]`, `g(lo + t span) = sum a_ij t_u^i t_v^j +
    /// remainder`: `a_ij` the derivative over `i! j!` times the spans'
    /// powers, `t` in the unit square, and the remainder at most the
    /// order-`K` derivatives' bounds times the same weights. The
    /// polynomial's Bernstein coefficients on the unit square bound it, so
    /// all of them clearing the remainder with one sign prove that sign.
    fn verdict(&self, coefficient: impl Fn(usize, usize) -> Option<Bound>) -> Option<Verdict> {
        let tables = tables();
        let (along_u, along_v) = (self.span.0 > 0.0, self.span.1 > 0.0);
        // `span_u^i span_v^j / (i! j!)`.
        let weight = |i: usize, j: usize| {
            let mut w = tables.inverse_factorial[i].mul(tables.inverse_factorial[j]);
            for _ in 0..i {
                w = w.mul(Bound::point(self.span.0));
            }
            for _ in 0..j {
                w = w.mul(Bound::point(self.span.1));
            }
            w
        };
        let present = |i: usize, j: usize| (along_u || i == 0) && (along_v || j == 0);
        // The derivatives' bounds over the box, each taken once.
        let memo = core::cell::RefCell::new(vec![vec![None; MAX_ORDER + 1]; MAX_ORDER + 1]);
        let bound = |i: usize, j: usize| -> Option<Option<Scalar>> {
            if let Some(b) = memo.borrow()[i][j] {
                return Some(b);
            }
            let b = self
                .tier
                .bound(self.lo, self.hi, (self.d.0 + i, self.d.1 + j))?;
            memo.borrow_mut()[i][j] = Some(b);
            Some(b)
        };
        // The highest power present along each parameter.
        let (mut top_u, mut top_v) = (0, 0);
        let zero = Bound::point(0.0);
        let mut a = vec![vec![zero; MAX_ORDER]; MAX_ORDER];
        let mut filled = 0;
        for &order in &ORDERS {
            // The coefficients of total order below `order`.
            for k in filled..order {
                for (i, row) in a.iter_mut().enumerate().take(k + 1) {
                    let j = k - i;
                    if !present(i, j) || bound(i, j)?.is_none() {
                        continue;
                    }
                    row[j] = coefficient(i, j)?.mul(weight(i, j));
                    (top_u, top_v) = (top_u.max(i), top_v.max(j));
                }
            }
            filled = order;
            // The remainder's bound at this order.
            let mut remainder: Scalar = 0.0;
            for i in 0..=order {
                let j = order - i;
                if !present(i, j) {
                    continue;
                }
                if let Some(b) = bound(i, j)? {
                    remainder = (remainder + b * weight(i, j).hi).next_up();
                }
            }
            // Bernstein coefficients of the polynomial's degree along each
            // parameter, converted one parameter at a time.
            let (nu, nv) = (top_u, top_v);
            let ratio = &tables.ratio;
            let mut c = vec![vec![zero; nv + 1]; nu + 1];
            for (p, row) in c.iter_mut().enumerate() {
                for (j, slot) in row.iter_mut().enumerate() {
                    for (r, row) in ratio[nu][p].iter().zip(&a) {
                        *slot = slot.add(r.mul(row[j]));
                    }
                }
            }
            let (mut positive, mut negative, mut unclear) = (false, false, false);
            let mut largest: Scalar = 0.0;
            for row in &c {
                for q in 0..=nv {
                    let mut b = zero;
                    for (j, value) in row.iter().enumerate().take(q + 1) {
                        b = b.add(ratio[nv][q][j].mul(*value));
                    }
                    largest = largest.max(b.lo.abs()).max(b.hi.abs());
                    if b.lo > remainder {
                        positive = true;
                    } else if b.hi < -remainder {
                        negative = true;
                    } else {
                        unclear = true;
                    }
                }
            }
            match (positive, negative, unclear) {
                (true, false, false) => return Some(Verdict::Sign(Sign::Positive)),
                (false, true, false) => return Some(Verdict::Sign(Sign::Negative)),
                (true, true, _) => return Some(Verdict::Mixed),
                // A coefficient near zero that a smaller remainder would
                // not settle: higher orders barely move it.
                _ if remainder <= 1e-6 * largest => return Some(Verdict::Unclear),
                _ => {}
            }
        }
        Some(Verdict::Unclear)
    }
}

/// The Taylor orders tried, rising.
const ORDERS: [usize; 11] = [1, 2, 3, 4, 6, 8, 11, 14, 18, 24, MAX_ORDER];

/// Constants of the Bernstein conversion, as intervals.
struct Tables {
    /// `1 / n!`.
    inverse_factorial: Vec<Bound>,
    /// `ratio[n][p][i] = C(p, i) / C(n, i)`: the monomial `t^i`'s weight
    /// in the degree-`n` Bernstein coefficient `p`.
    ratio: Vec<Vec<Vec<Bound>>>,
}

fn tables() -> &'static Tables {
    static TABLES: std::sync::OnceLock<Tables> = std::sync::OnceLock::new();
    TABLES.get_or_init(|| {
        let factorial = |n: usize| (1..=n).map(|k| k as Scalar).product::<Scalar>();
        let binomial = |n: usize, k: usize| factorial(n) / (factorial(k) * factorial(n - k));
        Tables {
            inverse_factorial: (0..=MAX_ORDER)
                .map(|n| Bound::near(1.0 / factorial(n)))
                .collect(),
            ratio: (0..MAX_ORDER)
                .map(|n| {
                    (0..=n)
                        .map(|p| {
                            (0..=p)
                                .map(|i| Bound::near(binomial(p, i) / binomial(n, i)))
                                .collect()
                        })
                        .collect()
                })
                .collect(),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn series(u: Basis, v: Basis, coefficients: Vec<Vec<Scalar>>) -> SeriesField2 {
        SeriesField2 { u, v, coefficients }
    }

    #[test]
    fn exact_zeros_are_found_by_the_laurent_test() {
        // sin(u) is zero at u = 0 only among f64s.
        let f = series(
            Basis::Fourier,
            Basis::Power,
            vec![vec![0.0], vec![0.0], vec![1.0]],
        );
        let t = SeriesTier::new(&f);
        assert_eq!(t.sign_at(Point2::new(0.0, 0.3), (0, 0)), Some(Sign::Zero));
        assert_eq!(
            t.sign_at(Point2::new(1e-300, 0.3), (0, 0)),
            Some(Sign::Positive)
        );
        assert_eq!(
            t.sign_at(Point2::new(-0.5, 0.3), (0, 0)),
            Some(Sign::Negative)
        );
        // Its derivative cos(u) at u = 0 is one.
        assert_eq!(
            t.sign_at(Point2::new(0.0, 0.3), (1, 0)),
            Some(Sign::Positive)
        );
        // cos(u) cos(v) - sin(u) sin(v) - cos(2u) vanishes on the diagonal
        // v = u, where two harmonics of different angles cancel exactly.
        let f = series(
            Basis::Fourier,
            Basis::Fourier,
            vec![
                vec![0.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, -1.0],
                vec![-1.0, 0.0, 0.0],
            ],
        );
        let t = SeriesTier::new(&f);
        for x in [0.0, 0.3, -1.25, 2.0, 1e-7] {
            assert_eq!(
                t.sign_at(Point2::new(x, x), (0, 0)),
                Some(Sign::Zero),
                "{x}"
            );
        }
        // Off the diagonal it is not zero.
        assert_ne!(
            t.sign_at(Point2::new(0.3, 0.30000000000000004), (0, 0)),
            Some(Sign::Zero)
        );
        // cos(u) - cos(2v) at u = 2v: angles of different binary exponents,
        // brought to one common angle.
        let f = series(
            Basis::Fourier,
            Basis::Fourier,
            vec![vec![0.0, 0.0, 0.0, -1.0], vec![1.0, 0.0, 0.0, 0.0]],
        );
        let t = SeriesTier::new(&f);
        for v in [0.3, 0.7, -1.1, 3e-5] {
            assert_eq!(
                t.sign_at(Point2::new(2.0 * v, v), (0, 0)),
                Some(Sign::Zero),
                "{v}"
            );
        }
    }

    #[test]
    fn signs_below_rounding_are_certified() {
        // cos(u) - cos(v) at v one ulp past u: far below what `f64`
        // evaluation of the two harmonics resolves.
        let f = series(
            Basis::Fourier,
            Basis::Fourier,
            vec![vec![0.0, -1.0], vec![1.0, 0.0]],
        );
        let t = SeriesTier::new(&f);
        let u: f64 = 0.7;
        let v = f64::from_bits(u.to_bits() + 1);
        // cos is decreasing there: cos(u) > cos(v).
        assert_eq!(t.sign_at(Point2::new(u, v), (0, 0)), Some(Sign::Positive));
        assert_eq!(t.sign_at(Point2::new(v, u), (0, 0)), Some(Sign::Negative));
        assert_eq!(t.sign_at(Point2::new(u, u), (0, 0)), Some(Sign::Zero));
    }

    #[test]
    fn a_box_off_a_high_order_zero_keeps_its_sign() {
        // sin(u)^7 in harmonics: at u in [1e-3, 1.1e-3] it is about 1e-21,
        // far below the harmonics' rounding, and positive throughout.
        let mut coefficients = vec![vec![0.0]; 15];
        // sin^7 u = (35 sin u - 21 sin 3u + 7 sin 5u - sin 7u) / 64.
        coefficients[2][0] = 35.0 / 64.0;
        coefficients[6][0] = -21.0 / 64.0;
        coefficients[10][0] = 7.0 / 64.0;
        coefficients[14][0] = -1.0 / 64.0;
        let f = series(Basis::Fourier, Basis::Power, coefficients);
        let t = SeriesTier::new(&f);
        let (lo, hi) = (Point2::new(1e-3, -1.0), Point2::new(1.1e-3, 1.0));
        assert_eq!(t.keeps_sign(lo, hi, (0, 0)), Some(Sign::Positive));
        let (lo, hi) = (Point2::new(-1.1e-3, -1.0), Point2::new(-1e-3, 1.0));
        assert_eq!(t.keeps_sign(lo, hi, (0, 0)), Some(Sign::Negative));
        // Across the zero no sign is kept.
        let (lo, hi) = (Point2::new(-1e-3, -1.0), Point2::new(1e-3, 1.0));
        assert_eq!(t.keeps_sign(lo, hi, (0, 0)), None);
        // The derivative 7 sin^6 cos is positive on both sides.
        let (lo, hi) = (Point2::new(-1.1e-3, 0.0), Point2::new(-1e-3, 0.0));
        assert_eq!(t.keeps_sign(lo, hi, (1, 0)), Some(Sign::Positive));
    }
}
