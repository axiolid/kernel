//! Exact signs of B-spline fields (ADR 0077), where `f64` cannot decide,
//! and the tier that dispatches to the series fields' own (`exact_series`).
//!
//! On one cell `[a, b] x [c, d]` of a [`PatchField2`] of degree `(p, q)`,
//! the field times `(b - a)^p (d - c)^q` is
//!
//! `N(u, v) = sum c[k][l] C(p, k) C(q, l) (u - a)^k (b - u)^(p - k) (v - c)^l (d - v)^(q - l)`,
//!
//! a polynomial with the field's sign and no division, so its value at any
//! `f64` point is exact in dyadic arithmetic. Over a box `[x0, x1]` in `u`,
//! `(u - a)` and `(b - u)` are combinations of `t = u - x0` and
//! `w = x1 - u`, both nonnegative there; expanding gives `N` as a form in
//! `(w, t)` whose coefficients all sharing one strict sign prove that sign
//! over the box. Derivatives keep this form: their coefficients are the
//! coefficients' finite differences, up to a positive factor.
//!
//! Every question is asked in outward-rounded intervals first and in
//! dyadic big integers only where the interval cannot decide
//! (`axiolid_exact`'s two tiers). Past the grid the edge cells' polynomials
//! continue, as the field's own evaluation does.

use axiolid_core::{Point2, Scalar};
use axiolid_curve::{Field2, PatchField2};
use axiolid_exact::{Arith, Dyadic, Interval};
use axiolid_guarantees::Sign;

use crate::exact_series::SeriesTier;

/// A field's certified tier: exact dyadic Bernstein arithmetic for a
/// B-spline field, certified harmonics for a series field
/// ([`SeriesTier`]). A series field's tier answers only where the field is
/// flat ([`SeriesTier::flat`]): elsewhere `f64` and subdivision decide, as
/// they always have, at a fraction of the cost.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Exact<'a> {
    Patches(PatchTier<'a>),
    Series(SeriesTier<'a>),
}

impl<'a> Exact<'a> {
    /// The certified tier of a field with finite coefficients.
    pub(crate) fn of(field: &'a Field2) -> Option<Self> {
        match field {
            Field2::Patches(f) if f.is_finite() => Some(Self::Patches(PatchTier { field: f })),
            Field2::Series(f) if f.is_finite() => Some(Self::Series(SeriesTier::new(f))),
            _ => None,
        }
    }

    /// The sign of the field's derivative of orders `d` (along `u`, `v`) at
    /// `p`: `Zero` only where it is exactly zero.
    pub(crate) fn sign_at(&self, p: Point2, d: (usize, usize)) -> Option<Sign> {
        match self {
            Self::Patches(t) => t.sign_at(p, d),
            Self::Series(t) => t.flat(p, d).then(|| t.sign_at(p, d)).flatten(),
        }
    }

    /// The strict sign the derivative of orders `d` keeps over the box
    /// `[lo, hi]` (either side may be a point), or `None` where it is not
    /// proven to keep one.
    pub(crate) fn keeps_sign(&self, lo: Point2, hi: Point2, d: (usize, usize)) -> Option<Sign> {
        match self {
            Self::Patches(t) => t.keeps_sign(lo, hi, d),
            Self::Series(t) => {
                let centre = lo + (hi - lo) * 0.5;
                t.flat(centre, d).then(|| t.keeps_sign(lo, hi, d)).flatten()
            }
        }
    }
}

/// A B-spline field's exact tier.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PatchTier<'a> {
    field: &'a PatchField2,
}

/// One axis of a box within a cell: a stretch `[x0, x1]`, or a point.
#[derive(Debug, Clone, Copy)]
struct Span {
    x0: Scalar,
    x1: Scalar,
}

impl PatchTier<'_> {
    fn cells(&self) -> (usize, usize) {
        (self.field.u_breaks.len() - 1, self.field.v_breaks.len() - 1)
    }

    /// The cell along one axis holding `x`, as the field's evaluation picks
    /// it (a break belongs to the cell after it; past the ends, the edge).
    fn cell_of(breaks: &[Scalar], x: Scalar) -> usize {
        let n = breaks.len() - 1;
        let mut i = 0;
        while i + 1 < n && x >= breaks[i + 1] {
            i += 1;
        }
        i
    }

    /// The sign of the field's derivative of orders `d` (along `u`, `v`) at
    /// `p`, exactly.
    pub(crate) fn sign_at(&self, p: Point2, d: (usize, usize)) -> Option<Sign> {
        if !p.is_finite() {
            return None;
        }
        let i = Self::cell_of(&self.field.u_breaks, p.x);
        let j = Self::cell_of(&self.field.v_breaks, p.y);
        let (su, sv) = (Span { x0: p.x, x1: p.x }, Span { x0: p.y, x1: p.y });
        // The interval tier's answers are certain; only its silence goes
        // to the exact one.
        let verdict = match self.signs::<Interval>(i, j, su, sv, d) {
            Some(v) => v,
            None => self.signs::<Dyadic>(i, j, su, sv, d)?,
        };
        match verdict {
            Verdict::Constant(s) => Some(s),
            Verdict::Varies => None,
        }
    }

    /// The strict sign the derivative of orders `d` keeps over the box
    /// `[lo, hi]` (either side may be a point), or `None` where it is not
    /// proven to keep one.
    pub(crate) fn keeps_sign(&self, lo: Point2, hi: Point2, d: (usize, usize)) -> Option<Sign> {
        if !(lo.is_finite() && hi.is_finite()) || lo.x > hi.x || lo.y > hi.y {
            return None;
        }
        let (n, m) = self.cells();
        let (ub, vb) = (&self.field.u_breaks, &self.field.v_breaks);
        let mut sign: Option<Sign> = None;
        for i in Self::cell_of(ub, lo.x)..=Self::cell_of(ub, hi.x) {
            for j in Self::cell_of(vb, lo.y)..=Self::cell_of(vb, hi.y) {
                // The part of the box this cell answers for (edge cells
                // continue past the grid).
                let clip = |lo: Scalar, hi: Scalar, k: usize, breaks: &[Scalar], count: usize| {
                    let a = if k == 0 { lo } else { lo.max(breaks[k]) };
                    let b = if k + 1 == count {
                        hi
                    } else {
                        hi.min(breaks[k + 1])
                    };
                    Span { x0: a, x1: b }
                };
                let su = clip(lo.x, hi.x, i, ub, n);
                let sv = clip(lo.y, hi.y, j, vb, m);
                if su.x0 > su.x1 || sv.x0 > sv.x1 {
                    continue;
                }
                let verdict = match self.signs::<Interval>(i, j, su, sv, d) {
                    Some(v) => v,
                    None => self.signs::<Dyadic>(i, j, su, sv, d)?,
                };
                let here = match verdict {
                    Verdict::Constant(s) => s,
                    Verdict::Varies => return None,
                };
                if here == Sign::Zero || sign.is_some_and(|s| s != here) {
                    return None;
                }
                sign = Some(here);
            }
        }
        sign
    }

    /// The sign of the derivative `d` on cell `(i, j)` over the spans, in
    /// arithmetic `T`: `None` where `T` cannot decide.
    fn signs<T: Arith>(
        &self,
        i: usize,
        j: usize,
        su: Span,
        sv: Span,
        d: (usize, usize),
    ) -> Option<Verdict> {
        let f = self.field;
        let (p, q) = (f.u_degree, f.v_degree);
        if d.0 > p || d.1 > q {
            return Some(Verdict::Constant(Sign::Zero));
        }
        let (_, m) = self.cells();
        let raw = &f.patches[i * m + j];
        // Coefficients, differenced `d` times along each axis.
        let mut c: Vec<Vec<T>> = (0..=p)
            .map(|a| (0..=q).map(|b| T::from_f64(raw[a * (q + 1) + b])).collect())
            .collect();
        for _ in 0..d.0 {
            c = (0..c.len() - 1)
                .map(|a| (0..c[0].len()).map(|b| c[a + 1][b].sub(&c[a][b])).collect())
                .collect();
        }
        for _ in 0..d.1 {
            c = c
                .iter()
                .map(|row| {
                    (0..row.len() - 1)
                        .map(|b| row[b + 1].sub(&row[b]))
                        .collect()
                })
                .collect();
        }
        // Restricted to the spans: along `u` column by column, then along
        // `v` row by row (a point collapses its axis to one value).
        let restrict = |line: &[T], a: Scalar, b: Scalar, span: Span| -> Vec<T> {
            if span.x0 == span.x1 {
                vec![casteljau_point(line, a, b, span.x0)]
            } else {
                let left = casteljau_left(line, a, b, span.x1);
                casteljau_right(&left, a, span.x1, span.x0)
            }
        };
        let (ua, ub) = (f.u_breaks[i], f.u_breaks[i + 1]);
        let (va, vb) = (f.v_breaks[j], f.v_breaks[j + 1]);
        let columns: Vec<Vec<T>> = (0..c[0].len())
            .map(|b| {
                let column: Vec<T> = c.iter().map(|row| row[b].clone()).collect();
                restrict(&column, ua, ub, su)
            })
            .collect();
        let rows: Vec<Vec<T>> = (0..columns[0].len())
            .map(|a| {
                let row: Vec<T> = columns.iter().map(|col| col[a].clone()).collect();
                restrict(&row, va, vb, sv)
            })
            .collect();
        let mut sign: Option<Sign> = None;
        let mut undecided = false;
        for row in &rows {
            for e in row {
                // Two coefficients of different certain signs settle it
                // whatever the undecided ones are.
                let Some(s) = e.sign() else {
                    undecided = true;
                    continue;
                };
                match sign {
                    None => sign = Some(s),
                    Some(t) if t == s => {}
                    Some(_) => return Some(Verdict::Varies),
                }
            }
        }
        if undecided {
            return None;
        }
        // A point is decided by its one value, zero included; a stretch
        // keeps a sign only where every coefficient is strictly of it.
        let stretch = su.x0 < su.x1 || sv.x0 < sv.x1;
        match sign {
            Some(Sign::Zero) if stretch => Some(Verdict::Varies),
            Some(s) => Some(Verdict::Constant(s)),
            None => Some(Verdict::Varies),
        }
    }
}

/// What a cell's coefficients say.
enum Verdict {
    /// All of one sign.
    Constant(Sign),
    /// Of different signs: no sign proven over the box.
    Varies,
}

/// The value at `x` of the Bernstein polynomial `c` on `[a, b]`, times
/// `(b - a)^n`.
fn casteljau_point<T: Arith>(c: &[T], a: Scalar, b: Scalar, x: Scalar) -> T {
    let (wl, wr) = (
        T::from_f64(b).sub(&T::from_f64(x)),
        T::from_f64(x).sub(&T::from_f64(a)),
    );
    let mut w = c.to_vec();
    for r in 1..w.len() {
        for i in 0..w.len() - r {
            w[i] = wl.mul(&w[i]).add(&wr.mul(&w[i + 1]));
        }
    }
    w[0].clone()
}

/// The coefficients of `c` (on `[a, b]`) over `[a, x]`, all times
/// `(b - a)^n`.
fn casteljau_left<T: Arith>(c: &[T], a: Scalar, b: Scalar, x: Scalar) -> Vec<T> {
    let n = c.len() - 1;
    let (wl, wr) = (
        T::from_f64(b).sub(&T::from_f64(x)),
        T::from_f64(x).sub(&T::from_f64(a)),
    );
    let scale = T::from_f64(b).sub(&T::from_f64(a));
    let mut w = c.to_vec();
    let mut out = Vec::with_capacity(n + 1);
    out.push(w[0].clone());
    for r in 1..=n {
        for i in 0..=n - r {
            w[i] = wl.mul(&w[i]).add(&wr.mul(&w[i + 1]));
        }
        out.push(w[0].clone());
    }
    // Level `r` carries `(b - a)^r`: bring every one to `(b - a)^n`.
    let mut factor = T::from_f64(1.0);
    for k in (0..=n).rev() {
        out[k] = out[k].mul(&factor);
        factor = factor.mul(&scale);
    }
    out
}

/// The coefficients of `c` (on `[a, b]`) over `[x, b]`, all times
/// `(b - a)^n`.
fn casteljau_right<T: Arith>(c: &[T], a: Scalar, b: Scalar, x: Scalar) -> Vec<T> {
    let n = c.len() - 1;
    let (wl, wr) = (
        T::from_f64(b).sub(&T::from_f64(x)),
        T::from_f64(x).sub(&T::from_f64(a)),
    );
    let scale = T::from_f64(b).sub(&T::from_f64(a));
    let mut w = c.to_vec();
    let mut out = vec![T::from_f64(0.0); n + 1];
    out[n] = w[n].clone();
    for r in 1..=n {
        for i in 0..=n - r {
            w[i] = wl.mul(&w[i]).add(&wr.mul(&w[i + 1]));
        }
        out[n - r] = w[n - r].clone();
    }
    let mut factor = T::from_f64(1.0);
    for value in &mut out {
        *value = value.mul(&factor);
        factor = factor.mul(&scale);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `x^k` on `[-1, 1]` as a patch field on the unit square in `(u, v)`
    /// (`x = 2u - 1`), constant in `v`, split into two cells at `u = 1/2`
    /// (the second cell's coefficients by the polynomial's values).
    fn power(k: usize) -> Field2 {
        // Bernstein coefficients of s^k over [lo, hi] (s = 2u - 1): the
        // blossom, the product of k endpoints.
        let over = |lo: f64, hi: f64| -> Vec<f64> {
            (0..=k)
                .flat_map(|i| {
                    let c = lo.powi((k - i) as i32) * hi.powi(i as i32);
                    [c, c]
                })
                .collect()
        };
        Field2::Patches(PatchField2 {
            u_breaks: vec![0.0, 0.5, 1.0],
            v_breaks: vec![0.0, 1.0],
            u_degree: k,
            v_degree: 1,
            patches: vec![over(-1.0, 0.0), over(0.0, 1.0)],
        })
    }

    #[test]
    fn signs_below_rounding_are_exact() {
        let f = power(9);
        let e = Exact::of(&f).unwrap();
        // x = 2u - 1: u a hair either side of 1/2 gives x^9 around 1e-135,
        // far below any rounding of coefficients near 1.
        for (u, want) in [
            (0.5 + 1e-15, Sign::Positive),
            (0.5 - 1e-15, Sign::Negative),
            (0.5, Sign::Zero),
            (0.3, Sign::Negative),
            (0.9, Sign::Positive),
        ] {
            assert_eq!(
                e.sign_at(Point2::new(u, 0.4), (0, 0)),
                Some(want),
                "u = {u}"
            );
        }
        // Where f64 decides, it agrees.
        for i in 0..=20 {
            let u = 0.05 * i as f64;
            let v = f.value(Point2::new(u, 0.3));
            if v.abs() > 1e-12 {
                let want = if v > 0.0 {
                    Sign::Positive
                } else {
                    Sign::Negative
                };
                assert_eq!(
                    e.sign_at(Point2::new(u, 0.3), (0, 0)),
                    Some(want),
                    "u = {u}"
                );
            }
        }
    }

    #[test]
    fn boxes_keep_a_sign_only_where_it_is_proven() {
        let f = power(9);
        let e = Exact::of(&f).unwrap();
        let at = |u: f64, v: f64| Point2::new(u, v);
        assert_eq!(
            e.keeps_sign(at(0.5 + 1e-12, 0.0), at(0.51, 1.0), (0, 0)),
            Some(Sign::Positive)
        );
        assert_eq!(
            e.keeps_sign(at(0.2, 0.0), at(0.5 - 1e-12, 1.0), (0, 0)),
            Some(Sign::Negative)
        );
        // Across the zero, and across the cell boundary: no sign.
        assert_eq!(e.keeps_sign(at(0.49, 0.0), at(0.51, 1.0), (0, 0)), None);
        // Touching zero at its edge: not strictly of one sign.
        assert_eq!(e.keeps_sign(at(0.5, 0.0), at(0.6, 1.0), (0, 0)), None);
        // Spanning both cells, away from zero on each side separately.
        assert_eq!(
            e.keeps_sign(at(0.55, 0.2), at(0.95, 0.4), (0, 0)),
            Some(Sign::Positive)
        );
        // The derivative 9 x^8 >= 0 vanishes only at x = 0: positive
        // on a box off it, even one reaching into the other cell.
        assert_eq!(
            e.keeps_sign(at(0.3, 0.0), at(0.5 - 1e-9, 1.0), (1, 0)),
            Some(Sign::Positive)
        );
        // Constant along v: that derivative is zero, never strictly signed.
        assert_eq!(e.keeps_sign(at(0.6, 0.0), at(0.7, 1.0), (0, 1)), None);
        // Past the grid the edge cell's polynomial continues.
        assert_eq!(e.sign_at(at(1.2, 0.5), (0, 0)), Some(Sign::Positive));
        assert_eq!(e.sign_at(at(-0.2, 0.5), (0, 0)), Some(Sign::Negative));
    }
}
