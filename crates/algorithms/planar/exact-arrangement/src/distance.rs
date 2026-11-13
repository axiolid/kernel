//! Exact distance predicates on [`ExactPoints`].
//!
//! Distances are Euclidean in the stored coordinates and never leave exact
//! arithmetic: a square root never appears, because every question is a
//! comparison of squares. For homogeneous points `P = (px, py, pw)` with
//! `pw > 0`, put `D(P, Q) = (qx pw - px qw, qy pw - py qw)`, which is
//! `pw qw (Q - P)`. For a point `V` and a segment `A`-`B`, let `u = D(A, B)`,
//! `v = D(A, V)`, `dot = u . v`, `uu = u . u` and `cross = ux vy - uy vx`.
//! Every weight is positive, so:
//!
//! - `V` projects at or before `A` when `dot <= 0`, and at or beyond `B` when
//!   `dot wb >= uu wv`;
//! - the foot of the perpendicular is
//!   `(ax uu wv + dot ux, ay uu wv + dot uy, wa uu wv)`;
//! - the distance from `V` to the line `AB` is compared with `g` by the sign
//!   of `cross^2 - g^2 wa^2 wv^2 uu`;
//! - the distance between two points is compared with `g` by the sign of
//!   `v . v - g^2 (wa wv)^2` with `v = D(a, b)`.
//!
//! Each sign is the interval tier first, then `Dyadic`, through `certify`.

use std::cmp::Ordering;

use axiolid_contracts::{GeomError, GeomResult, Sign};
use axiolid_exact::{certify, Arith, Dyadic, SignExpr};

use crate::points::{ordering, undecided, ExactPoints, PointId};

/// Where on a segment the point nearest to a given point lies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Nearest {
    /// At the first end, which includes the projection falling on it.
    Start,
    /// Strictly between the ends: the foot of the perpendicular.
    Interior,
    /// At the second end, which includes the projection falling on it.
    End,
}

type Homogeneous = [Dyadic; 3];

fn lift<T: Arith>(point: &Homogeneous) -> [T; 3] {
    [
        T::from_dyadic(&point[0]),
        T::from_dyadic(&point[1]),
        T::from_dyadic(&point[2]),
    ]
}

/// `D(p, q) = (qx pw - px qw, qy pw - py qw)`.
fn delta<T: Arith>(p: &[T; 3], q: &[T; 3]) -> [T; 2] {
    [
        q[0].mul(&p[2]).sub(&p[0].mul(&q[2])),
        q[1].mul(&p[2]).sub(&p[1].mul(&q[2])),
    ]
}

fn dot<T: Arith>(u: &[T; 2], v: &[T; 2]) -> T {
    u[0].mul(&v[0]).add(&u[1].mul(&v[1]))
}

#[derive(Clone, Copy)]
enum Question {
    /// The sign of `dot`.
    Behind,
    /// The sign of `dot wb - uu wv`.
    Beyond,
    /// The sign of `cross^2 - g^2 wa^2 wv^2 uu`.
    Line,
}

struct Projection<'a> {
    question: Question,
    start: &'a Homogeneous,
    end: &'a Homogeneous,
    point: &'a Homogeneous,
    length: f64,
}

impl SignExpr for Projection<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let (a, b, v) = (
            lift::<T>(self.start),
            lift::<T>(self.end),
            lift::<T>(self.point),
        );
        let (u, w) = (delta(&a, &b), delta(&a, &v));
        let along = dot(&u, &w);
        let uu = dot(&u, &u);
        match self.question {
            Question::Behind => along.sign(),
            Question::Beyond => along.mul(&b[2]).sub(&uu.mul(&v[2])).sign(),
            Question::Line => {
                let cross = u[0].mul(&w[1]).sub(&u[1].mul(&w[0]));
                let g = T::from_f64(self.length);
                let scale = g.mul(&g).mul(&a[2]).mul(&a[2]).mul(&v[2]).mul(&v[2]);
                cross.mul(&cross).sub(&scale.mul(&uu)).sign()
            }
        }
    }
}

struct PointPair<'a> {
    first: &'a Homogeneous,
    second: &'a Homogeneous,
    length: f64,
}

impl SignExpr for PointPair<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let (a, b) = (lift::<T>(self.first), lift::<T>(self.second));
        let v = delta(&a, &b);
        let g = T::from_f64(self.length);
        let weights = a[2].mul(&b[2]);
        dot(&v, &v)
            .sub(&g.mul(&g).mul(&weights).mul(&weights))
            .sign()
    }
}

fn certified<E: SignExpr>(expr: &E) -> GeomResult<Sign> {
    certify(expr).map_err(|_| undecided())
}

fn check_length(length: f64) -> GeomResult<()> {
    if length.is_finite() && length >= 0.0 {
        Ok(())
    } else {
        Err(GeomError::InvalidInput(
            "a distance bound is not finite and non-negative".into(),
        ))
    }
}

impl ExactPoints {
    fn check_segment(&self, start: PointId, end: PointId) -> GeomResult<()> {
        if self.same(start, end)? {
            return Err(GeomError::InvalidInput("a segment has equal ends".into()));
        }
        Ok(())
    }

    /// Where the point of the segment `start`-`end` nearest to `point` lies:
    /// an end when the projection of `point` falls on or beyond it, else the
    /// foot of the perpendicular. `InvalidInput` when the ends are equal.
    pub fn nearest_on_segment(
        &self,
        point: PointId,
        start: PointId,
        end: PointId,
    ) -> GeomResult<Nearest> {
        self.check_segment(start, end)?;
        let ask = |question| {
            certified(&Projection {
                question,
                start: self.exact(start),
                end: self.exact(end),
                point: self.exact(point),
                length: 0.0,
            })
        };
        if ask(Question::Behind)? != Sign::Positive {
            return Ok(Nearest::Start);
        }
        if ask(Question::Beyond)? != Sign::Negative {
            return Ok(Nearest::End);
        }
        Ok(Nearest::Interior)
    }

    /// The foot of the perpendicular from `point` to the line through `start`
    /// and `end`, an exact rational point. It lies on the segment exactly
    /// when [`ExactPoints::nearest_on_segment`] says `Interior`.
    /// `InvalidInput` when the ends are equal.
    pub fn foot(&mut self, point: PointId, start: PointId, end: PointId) -> GeomResult<PointId> {
        self.check_segment(start, end)?;
        let [a, b, v] = [start, end, point].map(|id| self.exact(id).clone());
        let (u, w) = (delta(&a, &b), delta(&a, &v));
        let (along, uu) = (dot(&u, &w), dot(&u, &u));
        let scale = uu.mul(&v[2]);
        let foot = [
            a[0].mul(&scale).add(&along.mul(&u[0])),
            a[1].mul(&scale).add(&along.mul(&u[1])),
            a[2].mul(&scale),
        ];
        self.push(foot)
    }

    /// The order of the distance between `a` and `b` against `length`:
    /// `Less` when they are closer, `Equal` at exactly `length`.
    /// `InvalidInput` unless `length` is finite and not negative.
    pub fn compare_distance(&self, a: PointId, b: PointId, length: f64) -> GeomResult<Ordering> {
        check_length(length)?;
        certified(&PointPair {
            first: self.exact(a),
            second: self.exact(b),
            length,
        })
        .map(ordering)
    }

    /// The order of the distance from `point` to the line through `start` and
    /// `end` against `length`: `Less` when it is closer, `Equal` at exactly
    /// `length`. `InvalidInput` unless `length` is finite and not negative,
    /// or when the ends are equal.
    pub fn compare_line_distance(
        &self,
        point: PointId,
        start: PointId,
        end: PointId,
        length: f64,
    ) -> GeomResult<Ordering> {
        check_length(length)?;
        self.check_segment(start, end)?;
        certified(&Projection {
            question: Question::Line,
            start: self.exact(start),
            end: self.exact(end),
            point: self.exact(point),
            length,
        })
        .map(ordering)
    }
}
