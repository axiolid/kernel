//! The sign questions of visibility, each written once for both tiers.
//!
//! A [`Constraint`] is an affine function of the point `P` on an edge, in
//! view space. It is evaluated at the edge's two end vertices; the sign pair
//! decides where along the edge it holds (see `param`). Constants such as the
//! occluder's vertices enter as exact view-space values.

use axiolid_exact::{Arith, SignExpr};
use axiolid_guarantees::Sign;

use crate::view::Vert;

/// The depth an occluder's plane is compared against.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Lift {
    /// The depth of the edge point itself.
    Own,
    /// The near plane.
    Near(f64),
}

/// What an affine constraint measures. `P` is the edge point, `'` the image
/// (`U`, `V`) of a point and `W` its depth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Kind {
    /// `±orient2d(p, q, P')`: the side of the image line `p'q'`.
    ImageSide { p: u32, q: u32, neg: bool },
    /// `±det3(t0, t1, t2, (P', lift))`: the position of `P` (or the near
    /// plane) against the plane of the triangle `t`, at the image of `P`.
    PlaneAbove { t: [u32; 3], neg: bool, lift: Lift },
    /// `(P' - C')·(D' - C')`: `P'` is not before the start of `C'D'`.
    ExtentFrom { from: u32, to: u32 },
    /// `(D' - P')·(D' - C')`: `P'` is not past the end of `C'D'`.
    ExtentTo { from: u32, to: u32 },
    /// `|F'|²·(W_F(P') - W_P)`: the segment `CD` is in front of `P`.
    SegmentAbove { from: u32, to: u32 },
    /// `|F'|²·(near - W_F(P'))`: the segment `CD` is not past the near plane.
    SegmentAtOrBelow { from: u32, to: u32, near: f64 },
    /// `near - W_P`.
    BelowNear(f64),
    /// `W_P - far`.
    AboveFar(f64),
    /// `bound - coordinate` (upper) or `coordinate - bound` (lower).
    Region {
        axis: usize,
        upper: bool,
        bound: f64,
    },
}

/// An affine constraint `L(P) >= 0`, or `L(P) > 0` when `strict`.
///
/// Strictness only matters where `L` is identically zero along the edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Constraint {
    pub kind: Kind,
    pub strict: bool,
}

impl Constraint {
    pub(crate) const fn new(kind: Kind, strict: bool) -> Self {
        Self { kind, strict }
    }

    /// `L(P)` for the vertex `at`.
    pub(crate) fn value<T: Arith>(&self, verts: &[Vert], at: u32) -> T {
        let c = |vertex: u32, axis: usize| T::from_dyadic(&verts[vertex as usize].exact[axis]);
        let signed = |neg: bool, x: T| if neg { x.neg() } else { x };
        match self.kind {
            Kind::ImageSide { p, q, neg } => {
                let side = c(q, 0)
                    .sub(&c(p, 0))
                    .mul(&c(at, 1).sub(&c(p, 1)))
                    .sub(&c(q, 1).sub(&c(p, 1)).mul(&c(at, 0).sub(&c(p, 0))));
                signed(neg, side)
            }
            Kind::PlaneAbove {
                t: [t0, t1, t2],
                neg,
                lift,
            } => {
                let (au, av, aw) = (
                    c(t1, 0).sub(&c(t0, 0)),
                    c(t1, 1).sub(&c(t0, 1)),
                    c(t1, 2).sub(&c(t0, 2)),
                );
                let (bu, bv, bw) = (
                    c(t2, 0).sub(&c(t0, 0)),
                    c(t2, 1).sub(&c(t0, 1)),
                    c(t2, 2).sub(&c(t0, 2)),
                );
                let depth = match lift {
                    Lift::Own => c(at, 2),
                    Lift::Near(near) => T::from_f64(near),
                };
                let (cu, cv, cw) = (
                    c(at, 0).sub(&c(t0, 0)),
                    c(at, 1).sub(&c(t0, 1)),
                    depth.sub(&c(t0, 2)),
                );
                let det = au
                    .mul(&bv.mul(&cw).sub(&bw.mul(&cv)))
                    .sub(&av.mul(&bu.mul(&cw).sub(&bw.mul(&cu))))
                    .add(&aw.mul(&bu.mul(&cv).sub(&bv.mul(&cu))));
                signed(neg, det)
            }
            Kind::ExtentFrom { from, to } => {
                let (fu, fv) = (c(to, 0).sub(&c(from, 0)), c(to, 1).sub(&c(from, 1)));
                c(at, 0)
                    .sub(&c(from, 0))
                    .mul(&fu)
                    .add(&c(at, 1).sub(&c(from, 1)).mul(&fv))
            }
            Kind::ExtentTo { from, to } => {
                let (fu, fv) = (c(to, 0).sub(&c(from, 0)), c(to, 1).sub(&c(from, 1)));
                c(to, 0)
                    .sub(&c(at, 0))
                    .mul(&fu)
                    .add(&c(to, 1).sub(&c(at, 1)).mul(&fv))
            }
            Kind::SegmentAbove { from, to } => {
                let (len2, along) = segment_terms(&c, from, to, at);
                len2.mul(&c(from, 2).sub(&c(at, 2)))
                    .add(&along.mul(&c(to, 2).sub(&c(from, 2))))
            }
            Kind::SegmentAtOrBelow { from, to, near } => {
                let (len2, along) = segment_terms(&c, from, to, at);
                T::from_f64(near).mul(&len2).sub(
                    &len2
                        .mul(&c(from, 2))
                        .add(&along.mul(&c(to, 2).sub(&c(from, 2)))),
                )
            }
            Kind::BelowNear(near) => T::from_f64(near).sub(&c(at, 2)),
            Kind::AboveFar(far) => c(at, 2).sub(&T::from_f64(far)),
            Kind::Region { axis, upper, bound } => {
                if upper {
                    T::from_f64(bound).sub(&c(at, axis))
                } else {
                    c(at, axis).sub(&T::from_f64(bound))
                }
            }
        }
    }
}

/// `|F'|²` and `(P' - C')·F'` for the segment `CD` and the point `at`.
fn segment_terms<T: Arith>(c: &impl Fn(u32, usize) -> T, from: u32, to: u32, at: u32) -> (T, T) {
    let (fu, fv) = (c(to, 0).sub(&c(from, 0)), c(to, 1).sub(&c(from, 1)));
    let len2 = fu.mul(&fu).add(&fv.mul(&fv));
    let along = c(at, 0)
        .sub(&c(from, 0))
        .mul(&fu)
        .add(&c(at, 1).sub(&c(from, 1)).mul(&fv));
    (len2, along)
}

/// The sign of a constraint at one vertex.
pub(crate) struct ValueSign<'a> {
    pub verts: &'a [Vert],
    pub constraint: &'a Constraint,
    pub at: u32,
}

impl SignExpr for ValueSign<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        self.constraint.value::<T>(self.verts, self.at).sign()
    }
}

/// The sign of `a_y·b_x - a_x·b_y`, where `a_i`, `b_i` are the values of
/// constraint `i` at the edge's start and end: with the signs of
/// `a_i - b_i` it orders the two roots `a_i / (a_i - b_i)`.
pub(crate) struct RootOrder<'a> {
    pub verts: &'a [Vert],
    pub x: &'a Constraint,
    pub y: &'a Constraint,
    pub start: u32,
    pub end: u32,
}

impl SignExpr for RootOrder<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let (ax, bx) = (
            self.x.value::<T>(self.verts, self.start),
            self.x.value::<T>(self.verts, self.end),
        );
        let (ay, by) = (
            self.y.value::<T>(self.verts, self.start),
            self.y.value::<T>(self.verts, self.end),
        );
        ay.mul(&bx).sub(&ax.mul(&by)).sign()
    }
}

/// `orient2d` of three vertices' images.
pub(crate) struct Orient2<'a> {
    pub verts: &'a [Vert],
    pub a: u32,
    pub b: u32,
    pub c: u32,
}

impl SignExpr for Orient2<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let c = |vertex: u32, axis: usize| T::from_dyadic(&self.verts[vertex as usize].exact[axis]);
        c(self.b, 0)
            .sub(&c(self.a, 0))
            .mul(&c(self.c, 1).sub(&c(self.a, 1)))
            .sub(
                &c(self.b, 1)
                    .sub(&c(self.a, 1))
                    .mul(&c(self.c, 0).sub(&c(self.a, 0))),
            )
            .sign()
    }
}

/// The difference of one view coordinate between two vertices.
pub(crate) struct Delta<'a> {
    pub verts: &'a [Vert],
    pub from: u32,
    pub to: u32,
    pub axis: usize,
}

impl SignExpr for Delta<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        T::from_dyadic(&self.verts[self.to as usize].exact[self.axis])
            .sub(&T::from_dyadic(
                &self.verts[self.from as usize].exact[self.axis],
            ))
            .sign()
    }
}
