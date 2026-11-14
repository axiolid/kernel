//! The sign questions of visibility, each written once for both tiers.
//!
//! A [`Constraint`] is an affine function of the point `P` on an edge, in
//! view space. It is evaluated at the edge's two end vertices; the sign pair
//! decides where along the edge it holds (see `param`). Constants such as the
//! occluder's vertices enter as exact view-space values.
//!
//! A perspective view (ADR 0107) keeps every question affine in the edge
//! parameter by asking it of vectors from the eye, which is the origin of
//! view space. The image of `X` is `(U, V) / h` with `h = -W`, so the side of
//! an image line is the determinant of the rows `(h, U, V)`. The eye ray
//! through `P` hits a triangle exactly when `P = a·t0 + b·t1 + c·t2` with
//! `a, b, c >= 0`, and then reaches the triangle's plane at `P / s` with
//! `s = a + b + c`: the plane is in front of `P` when `s > 1`, and the hit is
//! not past the near plane when `near·s - W_P >= 0`. No determinant divides
//! by a vertex weight, so a triangle behind the eye needs no special case.

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
/// `perspective` selects the eye-ray reading of the same question.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Constraint {
    pub kind: Kind,
    pub strict: bool,
    pub perspective: bool,
}

impl Constraint {
    pub(crate) const fn new(kind: Kind, strict: bool) -> Self {
        Self::in_view(kind, strict, false)
    }

    /// The constraint in a view that is perspective or not.
    pub(crate) const fn in_view(kind: Kind, strict: bool, perspective: bool) -> Self {
        Self {
            kind,
            strict,
            perspective,
        }
    }

    /// Whether `L` is zero at `at` by vertex identity alone: the image point
    /// is a vertex of the line, or the point is a vertex of the plane it is
    /// lifted to (its own depth). The same vertex has equal coordinates, so
    /// the determinant has a repeated row. Exact, never a tolerance.
    pub(crate) fn zero_by_identity(&self, at: u32) -> bool {
        match self.kind {
            Kind::ImageSide { p, q, .. } => at == p || at == q,
            Kind::PlaneAbove {
                t, lift: Lift::Own, ..
            } => t.contains(&at),
            _ => false,
        }
    }

    /// `L(P)` for the vertex `at`.
    pub(crate) fn value<T: Arith>(&self, verts: &[Vert], at: u32) -> T {
        if self.perspective {
            return self.eye_value(verts, at);
        }
        let c = |vertex: u32, axis: usize| T::from_dyadic(verts[vertex as usize].affine(axis));
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

type Vec3x<T> = [T; 3];

fn sub3<T: Arith>(a: &Vec3x<T>, b: &Vec3x<T>) -> Vec3x<T> {
    [a[0].sub(&b[0]), a[1].sub(&b[1]), a[2].sub(&b[2])]
}

fn scale3<T: Arith>(k: &T, a: &Vec3x<T>) -> Vec3x<T> {
    [k.mul(&a[0]), k.mul(&a[1]), k.mul(&a[2])]
}

fn dot3<T: Arith>(a: &Vec3x<T>, b: &Vec3x<T>) -> T {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn cross3<T: Arith>(a: &Vec3x<T>, b: &Vec3x<T>) -> Vec3x<T> {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

/// The view-space vector from the eye to the vertex.
fn from_eye<T: Arith>(verts: &[Vert], at: u32) -> Vec3x<T> {
    let exact = &verts[at as usize].exact;
    [
        T::from_dyadic(&exact[0]),
        T::from_dyadic(&exact[1]),
        T::from_dyadic(&exact[2]),
    ]
}

/// The determinant of the rows `(h, U, V)` of three vertices: the image
/// orientation of a perspective view, up to the positive weights `h`.
pub(crate) fn image_orient<T: Arith>(verts: &[Vert], a: u32, b: u32, c: u32) -> T {
    let row = |v: u32| {
        let exact = &verts[v as usize].exact;
        [
            T::from_dyadic(&exact[3]),
            T::from_dyadic(&exact[0]),
            T::from_dyadic(&exact[1]),
        ]
    };
    dot3(&row(a), &cross3(&row(b), &row(c)))
}

/// `h_from·X_to - h_to·X_from`: the difference of one image coordinate
/// between two perspective vertices, up to the positive weights `h`.
pub(crate) fn image_delta<T: Arith>(verts: &[Vert], from: u32, to: u32, axis: usize) -> T {
    let (f, t) = (&verts[from as usize].exact, &verts[to as usize].exact);
    T::from_dyadic(&f[3])
        .mul(&T::from_dyadic(&t[axis]))
        .sub(&T::from_dyadic(&t[3]).mul(&T::from_dyadic(&f[axis])))
}

impl Constraint {
    /// `L(P)` of a perspective view, for the vertex `at`.
    fn eye_value<T: Arith>(&self, verts: &[Vert], at: u32) -> T {
        let v = |vertex: u32| from_eye::<T>(verts, vertex);
        let signed = |neg: bool, x: T| if neg { x.neg() } else { x };
        match self.kind {
            Kind::ImageSide { p, q, neg } => signed(neg, image_orient(verts, p, q, at)),
            Kind::PlaneAbove {
                t: [t0, t1, t2],
                neg,
                lift,
            } => {
                let (a, b) = (sub3(&v(t1), &v(t0)), sub3(&v(t2), &v(t0)));
                // `P - t0`; or, at the near plane, `near·P - W_P·t0`, the
                // eye ray's point at depth `near` scaled by `W_P` to stay
                // affine.
                let c = match lift {
                    Lift::Own => sub3(&v(at), &v(t0)),
                    Lift::Near(near) => {
                        let p = v(at);
                        sub3(&scale3(&T::from_f64(near), &p), &scale3(&p[2], &v(t0)))
                    }
                };
                signed(neg, dot3(&a, &cross3(&b, &c)))
            }
            Kind::ExtentFrom { from, to } => {
                let (c, d) = (v(from), v(to));
                dot3(&cross3(&c, &v(at)), &cross3(&c, &d))
            }
            Kind::ExtentTo { from, to } => {
                let (c, d) = (v(from), v(to));
                dot3(&cross3(&v(at), &d), &cross3(&c, &d))
            }
            Kind::SegmentAbove { from, to } => {
                let (s, normal) = segment_scale(&v(from), &v(to), &v(at));
                s.sub(&dot3(&normal, &normal))
            }
            Kind::SegmentAtOrBelow { from, to, near } => {
                let (s, normal) = segment_scale(&v(from), &v(to), &v(at));
                T::from_f64(near)
                    .mul(&s)
                    .sub(&dot3(&normal, &normal).mul(&v(at)[2]))
            }
            Kind::BelowNear(near) => T::from_f64(near).sub(&v(at)[2]),
            Kind::AboveFar(far) => v(at)[2].sub(&T::from_f64(far)),
            Kind::Region { axis, upper, bound } => {
                let exact = &verts[at as usize].exact;
                let coordinate = T::from_dyadic(&exact[axis]);
                let scaled = T::from_f64(bound).mul(&T::from_dyadic(&exact[3]));
                if upper {
                    scaled.sub(&coordinate)
                } else {
                    coordinate.sub(&scaled)
                }
            }
        }
    }
}

/// For the segment `CD` and the point `P` in one plane through the eye:
/// `(P × (D - C))·(C × D)`, which is `s·|C × D|²` when `P / s` is on `CD`,
/// and the normal `C × D`.
fn segment_scale<T: Arith>(c: &Vec3x<T>, d: &Vec3x<T>, p: &Vec3x<T>) -> (T, Vec3x<T>) {
    let normal = cross3(c, d);
    (dot3(&cross3(p, &sub3(d, c)), &normal), normal)
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
        if self.constraint.zero_by_identity(self.at) {
            return Some(Sign::Zero);
        }
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
    pub perspective: bool,
}

impl SignExpr for Orient2<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        if self.c == self.a || self.c == self.b {
            return Some(Sign::Zero);
        }
        if self.perspective {
            return image_orient::<T>(self.verts, self.a, self.b, self.c).sign();
        }
        let c = |vertex: u32, axis: usize| T::from_dyadic(self.verts[vertex as usize].affine(axis));
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

/// The side of the eye (the origin of view space) of the line through two
/// corners of a triangle whose plane contains the eye: the sign of
/// `(v_from × v_to)·n` with `n` the triangle's normal. The three cross
/// products of such a triangle are all parallel to `n`, so their signs are
/// comparable, and the eye is in the closed triangle exactly when no two of
/// the three are strictly opposite.
pub(crate) struct EyeSide<'a> {
    pub verts: &'a [Vert],
    pub t: [u32; 3],
    /// The corner `t[k]` and the next one.
    pub k: usize,
}

impl SignExpr for EyeSide<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let v = |vertex: u32| from_eye::<T>(self.verts, vertex);
        let [t0, t1, t2] = self.t;
        let normal = cross3(&sub3(&v(t1), &v(t0)), &sub3(&v(t2), &v(t0)));
        let (from, to) = (self.t[self.k], self.t[(self.k + 1) % 3]);
        dot3(&cross3(&v(from), &v(to)), &normal).sign()
    }
}

/// The difference of one view coordinate between two vertices.
pub(crate) struct Delta<'a> {
    pub verts: &'a [Vert],
    pub from: u32,
    pub to: u32,
    pub axis: usize,
    pub perspective: bool,
}

impl SignExpr for Delta<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        if self.from == self.to {
            return Some(Sign::Zero);
        }
        if self.perspective {
            return image_delta::<T>(self.verts, self.from, self.to, self.axis).sign();
        }
        T::from_dyadic(self.verts[self.to as usize].affine(self.axis))
            .sub(&T::from_dyadic(
                self.verts[self.from as usize].affine(self.axis),
            ))
            .sign()
    }
}
