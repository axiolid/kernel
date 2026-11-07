//! The derivative tower of a second field along an implicit curve (#427,
//! ADR 0077), for the Rolle chain of [`crate::root_chain`].
//!
//! With `F` the field whose zero set is the curve and `H` the other, put
//! `H_0 = H` and `H_{j+1} = det(grad F, grad H_j) = F_u ∂_v H_j - F_v ∂_u H_j`.
//! Along the curve `d(G|C)/ds = ±det(grad F, grad G) / |grad F|`, so
//! `H_{j+1}|C` is `H_j|C`'s derivative up to a factor that never vanishes,
//! with a sign that is constant along a regular stretch. These values do
//! not depend on which parameter a cell leaves free, so a stretch can cross
//! a cell join.
//!
//! In a cell the curve is `(f, w(f))` with `f` free and `w` solved, and
//! `H_{j+1}|C = σ F_w|C d/df (H_j|C)` with `σ = -1` for `Axis::U` and `+1`
//! for `Axis::V`. So the tower follows from univariate Taylor coefficients
//! alone. Those of `w` come from the implicit equation, and the coefficient
//! of order `n` of a composite is a polynomial of the partials of the
//! fields up to order `n` and of the coefficients of `w`. Each partial is
//! enclosed over one box certain to hold the stretch's points (the
//! interval Taylor coefficients over the hull, the order-0 model), so the
//! results enclose every point of the stretch at once.

use std::cell::OnceCell;

use axiolid_core::{Point2, Scalar};
use axiolid_curve::implicit::{bound_simple, partial, Cell};
use axiolid_curve::{Axis, Field2, ImplicitCurve2};

use crate::interval::I;
use crate::root_chain::{Tower, MAX_CONTACT_ORDER};

/// Partials stored per field: `(i, j)` with `i + j <= MAX_CONTACT_ORDER`.
const STRIDE: usize = MAX_CONTACT_ORDER + 1;

/// A curve's field and a second field, with their partials built on demand.
pub(crate) struct CurveTower<'a> {
    curve: &'a ImplicitCurve2,
    other: &'a Field2,
    /// `∂_v^j ∂_u^i` of the curve's field at `i * STRIDE + j`.
    own: Vec<OnceCell<Field2>>,
    /// The same for the second field.
    second: Vec<OnceCell<Field2>>,
}

fn part<'f>(root: &'f Field2, cells: &'f [OnceCell<Field2>], i: usize, j: usize) -> &'f Field2 {
    if i == 0 && j == 0 {
        return root;
    }
    cells[i * STRIDE + j].get_or_init(|| {
        if j > 0 {
            partial(part(root, cells, i, j - 1), false)
        } else {
            partial(part(root, cells, i - 1, 0), true)
        }
    })
}

fn factorial(n: usize) -> Scalar {
    (1..=n).map(|k| k as Scalar).product()
}

impl<'a> CurveTower<'a> {
    pub(crate) fn new(curve: &'a ImplicitCurve2, other: &'a Field2) -> Self {
        let empty = || (0..STRIDE * STRIDE).map(|_| OnceCell::new()).collect();
        Self {
            curve,
            other,
            own: empty(),
            second: empty(),
        }
    }

    fn own_part(&self, i: usize, j: usize) -> &Field2 {
        part(&self.curve.field, &self.own, i, j)
    }

    fn other_part(&self, i: usize, j: usize) -> &Field2 {
        part(self.other, &self.second, i, j)
    }

    /// The tower over local parameters `[s0, s1]` of one cell: enclosures of
    /// `H_0..=H_order`. `None` for a bridge, or where the solved parameter's
    /// partial is not bounded away from zero.
    fn cell_enclose(&self, index: usize, s0: Scalar, s1: Scalar, order: usize) -> Option<Vec<I>> {
        let cell = &self.curve.cells[index];
        if cell.bridge.is_some() {
            return None;
        }
        let n = order;
        let axis = cell.axis;
        let free = |s: Scalar| cell.from + (cell.to - cell.from) * s;
        let place = |f: Scalar, w: Scalar| match axis {
            Axis::U => Point2::new(f, w),
            Axis::V => Point2::new(w, f),
        };
        let boxed = |f: (Scalar, Scalar), w: (Scalar, Scalar)| {
            let (a, b) = (place(f.0, w.0), place(f.1, w.1));
            Cell {
                lo: a.min(b),
                hi: a.max(b),
            }
        };
        let (f0, f1) = (free(s0), free(s1));
        let (f_lo, f_hi) = (f0.min(f1), f0.max(f1));
        let w0 = self.curve.solve_cell(cell, f0)?;
        // The partial along the free parameter and along the solved one.
        let (d_free, d_solved) = match axis {
            Axis::U => (self.own_part(1, 0), self.own_part(0, 1)),
            Axis::V => (self.own_part(0, 1), self.own_part(1, 0)),
        };
        let whole = boxed((f_lo, f_hi), (cell.low, cell.high));
        let solved = bound_simple(d_solved, &whole);
        if solved.lo <= 0.0 && solved.hi >= 0.0 {
            return None;
        }
        let floor = solved.lo.abs().min(solved.hi.abs());
        let big = {
            let r = bound_simple(d_free, &whole);
            r.lo.abs().max(r.hi.abs())
        };
        let reach = big / floor * (f_hi - f_lo);
        // The solved value itself is a rounded root: its residual bounds
        // its distance from the true one.
        let residual = bound_simple(&self.curve.field, &boxed((f0, f0), (w0, w0)));
        let slack = residual.lo.abs().max(residual.hi.abs()) / floor;
        let span = I::around(w0, reach + slack);
        let hull = boxed(
            (f_lo, f_hi),
            (span.lo.max(cell.low), span.hi.min(cell.high)),
        );
        let mut p = vec![vec![I::point(0.0); n + 1]; n + 1];
        let mut q = vec![vec![I::point(0.0); n + 1]; n + 1];
        for a in 0..=n {
            for b in 0..=n - a {
                p[a][b] = self.taylor(true, axis, &hull, a, b)?;
                q[a][b] = self.taylor(false, axis, &hull, a, b)?;
            }
        }
        // D[b][m]: the coefficient of order m of (w - w0)^b.
        let mut d = vec![vec![I::point(0.0); n + 1]; n + 1];
        d[0][0] = I::point(1.0);
        let mut w = vec![I::point(0.0); n + 1];
        for m in 1..=n {
            for b in 2..=m {
                let mut acc = I::point(0.0);
                for i in 1..=m + 1 - b {
                    acc = acc.add(w[i].mul(d[b - 1][m - i]));
                }
                d[b][m] = acc;
            }
            let mut acc = I::point(0.0);
            for a in 0..=m {
                for b in 0..=m - a {
                    if (a, b) == (0, 1) {
                        continue;
                    }
                    acc = acc.add(p[a][b].mul(d[b][m - a]));
                }
            }
            let zero = I::point(0.0);
            w[m] = zero.sub(acc).div(p[0][1])?;
            d[1][m] = w[m];
        }
        let sum = |c: &dyn Fn(usize, usize) -> I, m: usize| {
            let mut acc = I::point(0.0);
            for a in 0..=m {
                for (b, row) in d.iter().enumerate().take(m - a + 1) {
                    let x = row[m - a];
                    if x.lo == 0.0 && x.hi == 0.0 {
                        continue;
                    }
                    acc = acc.add(c(a, b).mul(x));
                }
            }
            acc
        };
        let h: Vec<I> = (0..=n).map(|m| sum(&|a, b| q[a][b], m)).collect();
        // F_w along the curve.
        let g: Vec<I> = (0..n)
            .map(|m| sum(&|a, b| p[a][b + 1].scale((b + 1) as Scalar), m))
            .collect();
        let sigma: Scalar = match axis {
            Axis::U => -1.0,
            Axis::V => 1.0,
        };
        let mut level = h;
        let mut out = vec![level[0]];
        for j in 0..n {
            let mut next = Vec::with_capacity(n - j);
            for m in 0..n - j {
                let mut acc = I::point(0.0);
                for i in 0..=m {
                    let derived = level[m - i + 1].scale((m - i + 1) as Scalar);
                    acc = acc.add(g[i].mul(derived));
                }
                next.push(acc.scale(sigma));
            }
            out.push(next[0]);
            level = next;
        }
        out.iter().all(|x| x.is_finite()).then_some(out)
    }

    /// The Taylor coefficient `∂_f^a ∂_w^b / (a! b!)` of the curve's field
    /// (`own`) or of the second one, enclosed over `hull`.
    fn taylor(&self, own: bool, axis: Axis, hull: &Cell, a: usize, b: usize) -> Option<I> {
        let (i, j) = match axis {
            Axis::U => (a, b),
            Axis::V => (b, a),
        };
        let field = if own {
            self.own_part(i, j)
        } else {
            self.other_part(i, j)
        };
        I::from_range(bound_simple(field, hull)).div(I::point(factorial(a) * factorial(b)))
    }

    fn cell_of(&self, t: Scalar) -> usize {
        (t.floor().max(0.0) as usize).min(self.curve.cells.len() - 1)
    }
}

impl Tower for CurveTower<'_> {
    fn enclose(&self, t0: Scalar, t1: Scalar, order: usize) -> Option<Vec<I>> {
        let n = self.curve.cells.len();
        if order > MAX_CONTACT_ORDER
            || n == 0
            || t0.is_nan()
            || t1.is_nan()
            || t0 > t1
            || t0 < 0.0
            || t1 > n as Scalar
        {
            return None;
        }
        if t0 == t1 {
            // A join belongs to the cell on its right, or to the one on its
            // left when that is a bridge.
            let c = self.cell_of(t0);
            let s = |c: usize| (t0 - c as Scalar).clamp(0.0, 1.0);
            return self.cell_enclose(c, s(c), s(c), order).or_else(|| {
                (c > 0 && t0 == c as Scalar)
                    .then(|| self.cell_enclose(c - 1, 1.0, 1.0, order))
                    .flatten()
            });
        }
        let first = self.cell_of(t0);
        let last = ((t1.ceil() as usize).max(1) - 1).min(n - 1);
        let mut out: Option<Vec<I>> = None;
        for c in first..=last {
            let (a, b) = (
                (t0 - c as Scalar).clamp(0.0, 1.0),
                (t1 - c as Scalar).clamp(0.0, 1.0),
            );
            let e = self.cell_enclose(c, a, b, order)?;
            out = Some(match out {
                None => e,
                Some(o) => o.iter().zip(&e).map(|(x, y)| x.hull(*y)).collect(),
            });
        }
        out
    }

    fn joins(&self, t0: Scalar, t1: Scalar) -> Vec<Scalar> {
        let n = self.curve.cells.len();
        let first = t0.floor() as i64 + 1;
        let last = (t1.ceil() as i64 - 1).min(n as i64 - 1);
        (first..=last)
            .map(|k| k as Scalar)
            .filter(|&k| k > t0 && k < t1)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::implicit_trace::{trace_within, Periodic};
    use axiolid_curve::implicit::{Basis, SeriesField2};

    fn power(terms: &[(usize, usize, Scalar)]) -> Field2 {
        let mut c = vec![vec![0.0; 5]; 5];
        for &(i, j, x) in terms {
            c[i][j] += x;
        }
        Field2::Series(SeriesField2 {
            u: Basis::Power,
            v: Basis::Power,
            coefficients: c,
        })
    }

    fn traced(field: &Field2) -> ImplicitCurve2 {
        let window = Cell {
            lo: Point2::new(-1.0, -1.0),
            hi: Point2::new(1.0, 1.0),
        };
        let mut curves = trace_within(field, window, Periodic { u: false, v: false }, 1 << 20)
            .expect("a traced curve");
        assert_eq!(curves.len(), 1);
        curves.remove(0)
    }

    fn contains(e: I, x: Scalar) -> bool {
        e.lo <= x && x <= e.hi
    }

    #[test]
    fn the_tower_of_a_graph_is_the_signed_derivatives() {
        // F = v - 0.3 u^2, H = F + u^4 + 0.5 u: on the curve H is
        // h(u) = u^4 + 0.5 u, and H_j|C = (-1)^j h^(j)(u).
        let f = power(&[(0, 1, 1.0), (2, 0, -0.3)]);
        let h = power(&[(0, 1, 1.0), (2, 0, -0.3), (4, 0, 1.0), (1, 0, 0.5)]);
        let curve = traced(&f);
        let tower = CurveTower::new(&curve, &h);
        for i in 0..=40 {
            let t = curve.end() * i as Scalar / 40.0;
            let u = curve.point(t).unwrap().x;
            let derivatives = [
                u.powi(4) + 0.5 * u,
                -(4.0 * u.powi(3) + 0.5),
                12.0 * u * u,
                -24.0 * u,
                24.0,
            ];
            let e = tower.enclose(t, t, 4).expect("an enclosure");
            for (j, x) in derivatives.iter().enumerate() {
                assert!(contains(e[j], *x), "H_{j} at t = {t}: {:?} vs {x}", e[j]);
                assert!(e[j].hi - e[j].lo < 1e-9, "H_{j} at t = {t}: {:?}", e[j]);
            }
        }
    }

    #[test]
    fn the_first_level_is_the_determinant_of_the_gradients_on_both_axes() {
        let f = power(&[(2, 0, 1.0), (0, 2, 1.0), (0, 0, -0.25)]);
        let h = power(&[(1, 0, 1.0), (0, 1, 2.0), (0, 0, 0.1)]);
        let curve = traced(&f);
        assert!(curve.cells.iter().any(|c| c.axis == Axis::U));
        assert!(curve.cells.iter().any(|c| c.axis == Axis::V));
        let tower = CurveTower::new(&curve, &h);
        for i in 0..=64 {
            let t = curve.end() * i as Scalar / 64.0;
            let p = curve.point(t).unwrap();
            let (jf, jh) = (f.jet(p), h.jet(p));
            let det = jf.gradient.x * jh.gradient.y - jf.gradient.y * jh.gradient.x;
            let e = tower.enclose(t, t, 1).expect("an enclosure");
            assert!(contains(e[1], det), "t = {t}: {:?} vs {det}", e[1]);
        }
    }

    #[test]
    fn a_stretch_across_joins_holds_the_values_inside_it() {
        let f = power(&[(2, 0, 1.0), (0, 2, 1.0), (0, 0, -0.25)]);
        let h = power(&[(1, 0, 1.0), (0, 1, 2.0), (0, 0, 0.1), (3, 0, 0.5)]);
        let curve = traced(&f);
        assert!(curve.cells.len() >= 3);
        let tower = CurveTower::new(&curve, &h);
        let (t0, t1) = (0.4, curve.end() - 0.4);
        assert!(!tower.joins(t0, t1).is_empty());
        let range = tower.enclose(t0, t1, 3).expect("a range enclosure");
        for i in 0..=50 {
            let t = t0 + (t1 - t0) * i as Scalar / 50.0;
            let point = tower.enclose(t, t, 3).expect("a point enclosure");
            for j in 0..=3 {
                let mid = point[j].mid();
                assert!(contains(range[j], mid), "H_{j} at t = {t}");
            }
        }
    }
}
