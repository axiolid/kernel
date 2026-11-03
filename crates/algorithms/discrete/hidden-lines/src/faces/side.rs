//! The point just beside a piece, and signs of affine forms at it.
//!
//! For a piece from `A` to `B` the side point is
//! `P(ε) = A + ε·d + ε²·n`, with `d` along the piece and `n` its left normal
//! (the right side uses `-n`). For an infinitesimal `ε` it lies inside the
//! open arrangement face beside the piece near `A`. An affine form
//! `f(P) = α·U + β·V + γ` then has the sign of the first nonzero of
//!
//! ```text
//! α·x₁ + β·y₁ + γ·w₁      (the value at A, times w₁ > 0)
//! α·d_u + β·d_v           (the slope along the piece)
//! α·n_u + β·n_v           (the slope across it)
//! ```
//!
//! and is `Equal` only when `f` vanishes identically on the path, which for a
//! difference of two depths means a true tie. Each of the three terms is a
//! [`SignExpr`], decided by the interval filter and then dyadic arithmetic.

use std::cmp::Ordering;

use axiolid_contracts::GeomResult;
use axiolid_exact::{Arith, Dyadic, SignExpr};
use axiolid_exact_arrangement::{ExactPoints, PointId};
use axiolid_guarantees::Sign;

use crate::decide::decide;
use crate::view::Vert;

/// An affine form `α·U + β·V + γ` in the image coordinates.
#[derive(Debug, Clone)]
pub(crate) struct Form {
    pub coeffs: [Dyadic; 3],
}

impl Form {
    pub(crate) fn new(alpha: Dyadic, beta: Dyadic, gamma: Dyadic) -> Self {
        Self {
            coeffs: [alpha, beta, gamma],
        }
    }

    fn negated(self) -> Self {
        Self {
            coeffs: self.coeffs.map(|c| c.neg()),
        }
    }

    /// `±orient2d(p, q, P)` for the image line `p'q'`.
    pub(crate) fn orient(verts: &[Vert], p: u32, q: u32) -> Self {
        let c = |v: u32, axis: usize| &verts[v as usize].exact[axis];
        let alpha = c(p, 1).sub(c(q, 1));
        let beta = c(q, 0).sub(c(p, 0));
        let gamma = alpha.mul(c(p, 0)).add(&beta.mul(c(p, 1))).neg();
        Self::new(alpha, beta, gamma)
    }

    /// `coordinate - bound` (lower) or `bound - coordinate` (upper).
    pub(crate) fn bound(axis: usize, upper: bool, bound: f64) -> Self {
        let (one, zero) = (Dyadic::from_f64(1.0), Dyadic::zero());
        let bound = Dyadic::from_f64(bound);
        let (mut alpha, mut beta) = (zero.clone(), zero);
        let coefficient = if upper { one.neg() } else { one };
        if axis == 0 {
            alpha = coefficient;
        } else {
            beta = coefficient;
        }
        let gamma = if upper { bound } else { bound.neg() };
        Self::new(alpha, beta, gamma)
    }
}

/// The plane of a triangle with nonzero image area, `W = N(U, V) / n_w` with
/// `N = c - n_u·U - n_v·V` and `n = (t1 - t0) × (t2 - t0)` in view space.
///
/// `n_w` is the image orientation of the triangle, so its sign is the
/// triangle's `sigma`.
#[derive(Debug, Clone)]
pub(crate) struct Plane {
    n: [Dyadic; 3],
    c: Dyadic,
    /// The three edge forms, `orient2d(t_i, t_{i+1}, P)`.
    pub sides: [Form; 3],
    /// Whether `n_w` is negative (the triangle is clockwise in the image).
    clockwise: bool,
}

impl Plane {
    /// The plane of the triangle `t`; `clockwise` is its image orientation.
    pub(crate) fn new(verts: &[Vert], t: [u32; 3], clockwise: bool) -> Self {
        let c = |v: u32, axis: usize| &verts[v as usize].exact[axis];
        let d = |to: u32, axis: usize| c(to, axis).sub(c(t[0], axis));
        let (a, b) = (
            [d(t[1], 0), d(t[1], 1), d(t[1], 2)],
            [d(t[2], 0), d(t[2], 1), d(t[2], 2)],
        );
        let n = [
            a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
            a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
            a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
        ];
        let at_t0 = n[0]
            .mul(c(t[0], 0))
            .add(&n[1].mul(c(t[0], 1)))
            .add(&n[2].mul(c(t[0], 2)));
        let sides = [(0, 1), (1, 2), (2, 0)].map(|(i, j)| Form::orient(verts, t[i], t[j]));
        Self {
            n,
            c: at_t0,
            sides,
            clockwise,
        }
    }

    /// `sign(n_w) · f`, so that the sign of the form is that of the quantity
    /// `f / n_w`.
    fn oriented(&self, form: Form) -> Form {
        if self.clockwise {
            form.negated()
        } else {
            form
        }
    }

    /// `bound - W` as a form (positive: the plane is below `bound`).
    pub(crate) fn below(&self, bound: f64) -> Form {
        let gamma = Dyadic::from_f64(bound).mul(&self.n[2]).sub(&self.c);
        self.oriented(Form::new(self.n[0].clone(), self.n[1].clone(), gamma))
    }

    /// `W - bound` as a form.
    pub(crate) fn above(&self, bound: f64) -> Form {
        self.below(bound).negated()
    }

    /// `W_self - W_other` as a form (positive: `self` is in front).
    pub(crate) fn in_front_of(&self, other: &Self) -> Form {
        let (x, y) = (self, other);
        let alpha = x.n[0].mul(&y.n[2]).neg().add(&y.n[0].mul(&x.n[2]));
        let beta = x.n[1].mul(&y.n[2]).neg().add(&y.n[1].mul(&x.n[2]));
        let gamma = x.c.mul(&y.n[2]).sub(&y.c.mul(&x.n[2]));
        let form = Form::new(alpha, beta, gamma);
        if x.clockwise == y.clockwise {
            form
        } else {
            form.negated()
        }
    }
}

/// The side point of one side of one piece.
#[derive(Debug, Clone)]
pub(crate) struct SidePoint {
    a: [Dyadic; 3],
    d: [Dyadic; 2],
    n: [Dyadic; 2],
}

impl SidePoint {
    /// The point beside the piece `from -> to`, on its left or its right.
    pub(crate) fn new(points: &ExactPoints, from: PointId, to: PointId, left: bool) -> Self {
        let (a, b) = (points.exact(from), points.exact(to));
        let d = [0, 1].map(|axis| b[axis].mul(&a[2]).sub(&a[axis].mul(&b[2])));
        let n = if left {
            [d[1].neg(), d[0].clone()]
        } else {
            [d[1].clone(), d[0].neg()]
        };
        Self { a: a.clone(), d, n }
    }
}

/// One term of the lexicographic sign.
struct Term<'a> {
    side: &'a SidePoint,
    form: &'a Form,
    order: usize,
}

impl SignExpr for Term<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let [alpha, beta, gamma] = self.form.coeffs.each_ref().map(T::from_dyadic);
        let value = |u: &Dyadic, v: &Dyadic| {
            alpha
                .mul(&T::from_dyadic(u))
                .add(&beta.mul(&T::from_dyadic(v)))
        };
        match self.order {
            0 => value(&self.side.a[0], &self.side.a[1])
                .add(&gamma.mul(&T::from_dyadic(&self.side.a[2])))
                .sign(),
            1 => value(&self.side.d[0], &self.side.d[1]).sign(),
            _ => value(&self.side.n[0], &self.side.n[1]).sign(),
        }
    }
}

/// The sign of `form` at the side point, as an ordering against zero.
pub(crate) fn sign_at(
    side: &SidePoint,
    form: &Form,
    escalations: &mut usize,
) -> GeomResult<Ordering> {
    for order in 0..3 {
        let term = Term { side, form, order };
        let sign = decide(&term, escalations)?;
        if sign != Ordering::Equal {
            return Ok(sign);
        }
    }
    Ok(Ordering::Equal)
}

#[cfg(test)]
mod tests {
    use axiolid_core::Point2;
    use axiolid_mesh_hidden_lines_contract::MeshRole::Sheet;

    use super::*;
    use crate::faces::testing::{p, plan, sheet, Setup};
    use axiolid_mesh::TriMesh;

    #[test]
    fn depth_ties_are_equal_only_for_coplanar_triangles() {
        // Triangles 0 and 2 are coplanar (z = 0); triangle 4 is a ramp that
        // meets z = 0 along y = 2.
        let ramp = TriMesh::new(
            vec![
                p(0.0, 0.0, -1.0),
                p(4.0, 0.0, -1.0),
                p(4.0, 4.0, 1.0),
                p(0.0, 4.0, 1.0),
            ],
            vec![0, 1, 2, 0, 2, 3],
        );
        let setup = Setup::new(
            &[
                (Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0)),
                (Sheet, sheet(1.0, 1.0, 3.0, 3.0, 0.0)),
                (Sheet, ramp),
            ],
            plan(),
            None,
        );
        let mut work = setup.work(usize::MAX);
        let faces = setup.faces(&mut work);
        let plane = |t: usize| faces.planes[t].as_ref().expect("plane");
        let mut points = setup.points();
        let a = points.store.push_flat(Point2::new(1.5, 2.0)).expect("a");
        let b = points.store.push_flat(Point2::new(2.5, 2.0)).expect("b");
        let side = SidePoint::new(&points.store, a, b, true);
        let mut escalations = 0;
        // Meshes: 0 -> triangles 0, 1; 1 -> 2, 3; 2 -> 4, 5.
        let coplanar = sign_at(&side, &plane(0).in_front_of(plane(2)), &mut escalations);
        assert_eq!(coplanar.expect("sign"), Ordering::Equal);
        // Along y = 2 the ramp and the sheet agree on the line but not
        // beside it.
        let crossing =
            sign_at(&side, &plane(0).in_front_of(plane(4)), &mut escalations).expect("sign");
        assert_ne!(crossing, Ordering::Equal);
        let reverse = sign_at(&side, &plane(4).in_front_of(plane(0)), &mut escalations);
        assert_eq!(reverse.expect("sign"), crossing.reverse());
    }
}
