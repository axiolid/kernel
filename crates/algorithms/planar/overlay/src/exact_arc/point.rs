//! Exact planar points and the sign questions asked about them.
//!
//! Every point the arc overlay handles is an input vertex (a pair of
//! doubles), a segment/segment crossing (rational), or a crossing involving
//! a circle, which carries one square root. One representation covers all
//! three:
//!
//! ```text
//! x = (xa + xb * sqrt(d)) / w      y = (ya + yb * sqrt(d)) / w
//! ```
//!
//! with dyadic `xa, xb, ya, yb, w, d`, `w > 0` and `d >= 0`. A question
//! about up to three points embeds them in one [`Tower`] (one radical per
//! distinct radicand) and runs the interval filter first, exact arithmetic
//! only when the filter cannot decide.

use axiolid_core::Point2;
use axiolid_exact::{certify, Arith, Dyadic, Interval, Nested, SignExpr, Tower};
use axiolid_guarantees::Sign;

/// The exact sign of a dyadic value.
pub(crate) fn sgn(value: &Dyadic) -> Sign {
    value.sign().expect("dyadic signs are always decided")
}

/// A dyadic value from a finite double. Callers validated finiteness.
pub(crate) fn dy(value: f64) -> Dyadic {
    Dyadic::from_f64(value)
}

/// `x * 2^n` without intermediate overflow; saturates like `f64` would.
fn scale2(mut x: f64, mut n: i64) -> f64 {
    while n > 1000 && x.is_finite() && x != 0.0 {
        x *= 2f64.powi(1000);
        n -= 1000;
    }
    while n < -1000 && x != 0.0 {
        x *= 2f64.powi(-1000);
        n += 1000;
    }
    x * 2f64.powi(n as i32)
}

/// `(a + b * sqrt(d)) / w` to about 50 bits, for output only.
///
/// Mantissas and exponents are combined separately, so coordinates whose
/// homogeneous parts are far outside the `f64` range still come out right.
/// A sum that cancels loses relative precision, which only matters for
/// coordinates tiny next to their own parts; output rounding is documented
/// as approximate.
fn approx(a: &Dyadic, b: &Dyadic, d: &Dyadic, w: &Dyadic) -> f64 {
    let mut terms: Vec<(f64, i64)> = Vec::with_capacity(2);
    if sgn(a) != Sign::Zero {
        terms.push(a.approx_parts());
    }
    if sgn(b) != Sign::Zero && sgn(d) != Sign::Zero {
        let (mb, eb) = b.approx_parts();
        let (mut md, mut ed) = d.approx_parts();
        if ed % 2 != 0 {
            md *= 2.0;
            ed -= 1;
        }
        terms.push((mb * md.sqrt(), eb + ed / 2));
    }
    let Some(top) = terms.iter().map(|t| t.1).max() else {
        return 0.0;
    };
    let sum: f64 = terms.iter().map(|&(m, e)| scale2(m, e - top)).sum();
    let (mw, ew) = w.approx_parts();
    scale2(sum / mw, top - ew)
}

/// A circle `alpha * |X|^2 + bx * x + by * y + gamma = 0`, `alpha > 0`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Circle {
    pub(crate) alpha: Dyadic,
    pub(crate) bx: Dyadic,
    pub(crate) by: Dyadic,
    pub(crate) gamma: Dyadic,
}

impl Circle {
    /// Whether two circle equations describe the same circle.
    ///
    /// Both have `alpha > 0`, so they agree iff they are proportional.
    pub(crate) fn same_as(&self, other: &Self) -> bool {
        let cross =
            |a: &Dyadic, b: &Dyadic| sgn(&other.alpha.mul(a).sub(&self.alpha.mul(b))) == Sign::Zero;
        cross(&self.bx, &other.bx) && cross(&self.by, &other.by) && cross(&self.gamma, &other.gamma)
    }

    /// `bx^2 + by^2 - 4 alpha gamma`, which is `(2 alpha r)^2`.
    pub(crate) fn delta(&self) -> Dyadic {
        self.bx
            .square()
            .add(&self.by.square())
            .sub(&dy(4.0).mul(&self.alpha).mul(&self.gamma))
    }

    /// The highest (`Positive`) or lowest (`Negative`) point, exactly:
    /// `(-bx, -by +- sqrt(delta)) / (2 alpha)`.
    pub(crate) fn extreme(&self, which: Sign) -> XPoint {
        let up = if which == Sign::Negative {
            dy(-1.0)
        } else {
            dy(1.0)
        };
        XPoint::new(
            self.bx.neg(),
            Dyadic::zero(),
            self.by.neg(),
            up,
            dy(2.0).mul(&self.alpha),
            self.delta(),
        )
    }

    /// Approximate centre, for output only.
    pub(crate) fn approx_centre(&self) -> Point2 {
        let two_alpha = self.alpha.to_f64() * 2.0;
        Point2::new(-self.bx.to_f64() / two_alpha, -self.by.to_f64() / two_alpha)
    }
}

/// An exact point; see the module documentation.
#[derive(Debug, Clone)]
pub(crate) struct XPoint {
    xa: Dyadic,
    xb: Dyadic,
    ya: Dyadic,
    yb: Dyadic,
    w: Dyadic,
    /// Zero for rational points, which then embed without a radical.
    d: Dyadic,
    approx: Point2,
    /// Sound boxes around `x` and `y`: cheap early answers for coordinate
    /// comparisons and point identity, which dominate the overlay's work.
    bx: Interval,
    by: Interval,
}

/// A sound enclosure of `(a + b * sqrt(d)) / w`.
fn enclose(a: &Dyadic, b: &Dyadic, d: &Dyadic, w: &Dyadic) -> Interval {
    let mut num = a.enclosure();
    if sgn(b) != Sign::Zero {
        let root = d.enclosure().sqrt_enclosure().unwrap_or(Interval::WHOLE);
        num = num.add(&b.enclosure().mul(&root));
    }
    num.quotient(w.enclosure())
}

impl XPoint {
    /// A general point. `w` must be non-zero and `d` non-negative; the
    /// sign of `w` is normalised to positive here.
    pub(crate) fn new(
        xa: Dyadic,
        xb: Dyadic,
        ya: Dyadic,
        yb: Dyadic,
        w: Dyadic,
        d: Dyadic,
    ) -> Self {
        debug_assert_ne!(sgn(&w), Sign::Zero, "a point needs a non-zero weight");
        debug_assert_ne!(sgn(&d), Sign::Negative, "a radicand must not be negative");
        let flip = sgn(&w) == Sign::Negative;
        let fix = |v: Dyadic| if flip { v.neg() } else { v };
        let (xa, xb, ya, yb, w) = (fix(xa), fix(xb), fix(ya), fix(yb), fix(w));
        let radical = sgn(&d) != Sign::Zero && (sgn(&xb) != Sign::Zero || sgn(&yb) != Sign::Zero);
        let (xb, yb, d) = if radical {
            (xb, yb, d)
        } else {
            (Dyadic::zero(), Dyadic::zero(), Dyadic::zero())
        };
        let approx = Point2::new(approx(&xa, &xb, &d, &w), approx(&ya, &yb, &d, &w));
        let bx = enclose(&xa, &xb, &d, &w);
        let by = enclose(&ya, &yb, &d, &w);
        Self {
            bx,
            by,
            xa,
            xb,
            ya,
            yb,
            w,
            d,
            approx,
        }
    }

    /// `(x / w, y / w)`.
    pub(crate) fn rational(x: Dyadic, y: Dyadic, w: Dyadic) -> Self {
        Self::new(x, Dyadic::zero(), y, Dyadic::zero(), w, Dyadic::zero())
    }

    /// An input vertex, exactly.
    pub(crate) fn from_f64(p: Point2) -> Self {
        Self::rational(dy(p.x), dy(p.y), dy(1.0))
    }

    /// True when the point carries no square root.
    pub(crate) fn is_rational(&self) -> bool {
        sgn(&self.d) == Sign::Zero
    }

    /// A nearby double pair, for output only; never for decisions.
    pub(crate) fn approx(&self) -> Point2 {
        self.approx
    }

    /// Sound enclosures of `x` and `y`, as `(lo, hi)` pairs.
    pub(crate) fn enclosures(&self) -> ((f64, f64), (f64, f64)) {
        ((self.bx.lo(), self.bx.hi()), (self.by.lo(), self.by.hi()))
    }
}

/// A sign question about exact points.
pub(crate) enum Pred<'a> {
    /// Positive when `c` lies left of the directed line `a -> b`.
    Orient(&'a XPoint, &'a XPoint, &'a XPoint),
    /// Sign of `dot(b - a, dir)`.
    Dot(&'a XPoint, &'a XPoint, &'a (Dyadic, Dyadic)),
    /// Sign of the circle's equation at the point: zero on it, negative
    /// inside.
    OnCircle(&'a XPoint, &'a Circle),
    /// Sign of `x(a) - x(b)`.
    DiffX(&'a XPoint, &'a XPoint),
    /// Sign of `y(a) - y(b)`.
    DiffY(&'a XPoint, &'a XPoint),
    /// Sign of `cross(u, v)` (or `dot(u, v)` when `cross` is false) for
    /// the travel tangents `u`, `v` of two curves at the point.
    Tangents {
        at: &'a XPoint,
        u: &'a Tangent,
        v: &'a Tangent,
        cross: bool,
    },
}

/// The travel direction of a curve, as a function of the point.
#[derive(Debug, Clone)]
pub(crate) enum Tangent {
    /// A constant direction (a segment traversed some way).
    Fixed(Dyadic, Dyadic),
    /// `factor * perp(grad F)` for a circle `F`: `factor` positive is
    /// counter-clockwise travel.
    Circle(Circle, Sign),
}

/// Points embedded in one tower, sharing radicals with equal radicands.
struct Embed<T> {
    tower: Tower<T>,
    roots: Vec<(Dyadic, Nested<T>)>,
}

impl<T: Arith> Embed<T> {
    fn new() -> Self {
        Self {
            tower: Tower::new(),
            roots: Vec::new(),
        }
    }

    fn k(&self, value: &Dyadic) -> Nested<T> {
        self.tower.value(T::from_dyadic(value))
    }

    /// `[X, Y, W]` with `x = X / W`, `y = Y / W`, `W > 0`.
    fn point(&mut self, p: &XPoint) -> Option<[Nested<T>; 3]> {
        let w = self.k(&p.w);
        if p.is_rational() {
            return Some([self.k(&p.xa), self.k(&p.ya), w]);
        }
        let known = self
            .roots
            .iter()
            .find(|(d, _)| *d == p.d)
            .map(|(_, root)| root.clone());
        let root = match known {
            Some(root) => root,
            None => {
                let radicand = self.k(&p.d);
                // At most three points per question: depth <= 3.
                let root = self.tower.sqrt(&radicand).ok()?;
                self.roots.push((p.d.clone(), root.clone()));
                root
            }
        };
        let t = &self.tower;
        let x = t.add(&self.k(&p.xa), &t.mul(&self.k(&p.xb), &root));
        let y = t.add(&self.k(&p.ya), &t.mul(&self.k(&p.yb), &root));
        Some([x, y, w])
    }
}

impl SignExpr for Pred<'_> {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let mut e = Embed::<T>::new();
        let value = match self {
            Pred::Orient(a, b, c) => {
                let [ax, ay, aw] = e.point(a)?;
                let [bx, by, bw] = e.point(b)?;
                let [cx, cy, cw] = e.point(c)?;
                let t = &e.tower;
                // det [[ax ay aw] [bx by bw] [cx cy cw]]; all w > 0.
                let m1 = t.sub(&t.mul(&by, &cw), &t.mul(&bw, &cy));
                let m2 = t.sub(&t.mul(&bx, &cw), &t.mul(&bw, &cx));
                let m3 = t.sub(&t.mul(&bx, &cy), &t.mul(&by, &cx));
                t.add(&t.sub(&t.mul(&ax, &m1), &t.mul(&ay, &m2)), &t.mul(&aw, &m3))
            }
            Pred::Dot(a, b, dir) => {
                let [ax, ay, aw] = e.point(a)?;
                let [bx, by, bw] = e.point(b)?;
                let (dx, dy) = (e.k(&dir.0), e.k(&dir.1));
                let t = &e.tower;
                // (b - a) * aw * bw, dotted with dir.
                let ex = t.sub(&t.mul(&bx, &aw), &t.mul(&ax, &bw));
                let ey = t.sub(&t.mul(&by, &aw), &t.mul(&ay, &bw));
                t.add(&t.mul(&ex, &dx), &t.mul(&ey, &dy))
            }
            Pred::OnCircle(p, c) => {
                let [x, y, w] = e.point(p)?;
                let (alpha, bx, by, gamma) = (e.k(&c.alpha), e.k(&c.bx), e.k(&c.by), e.k(&c.gamma));
                let t = &e.tower;
                // w^2 * F(x / w, y / w)
                let square = t.add(&t.mul(&x, &x), &t.mul(&y, &y));
                let linear = t.add(&t.mul(&bx, &x), &t.mul(&by, &y));
                t.add(
                    &t.add(&t.mul(&alpha, &square), &t.mul(&linear, &w)),
                    &t.mul(&gamma, &t.mul(&w, &w)),
                )
            }
            Pred::DiffX(a, b) | Pred::DiffY(a, b) => {
                let [ax, ay, aw] = e.point(a)?;
                let [bx, by, bw] = e.point(b)?;
                let t = &e.tower;
                if matches!(self, Pred::DiffX(..)) {
                    t.sub(&t.mul(&ax, &bw), &t.mul(&bx, &aw))
                } else {
                    t.sub(&t.mul(&ay, &bw), &t.mul(&by, &aw))
                }
            }
            Pred::Tangents { at, u, v, cross } => {
                let [x, y, w] = e.point(at)?;
                let vector = |tangent: &Tangent| -> [Nested<T>; 2] {
                    match tangent {
                        Tangent::Fixed(dx, dy) => [e.k(dx), e.k(dy)],
                        Tangent::Circle(c, factor) => {
                            let t = &e.tower;
                            let two_alpha = t.mul(&e.k(&dy(2.0)), &e.k(&c.alpha));
                            // grad * w = (2 alpha x + bx w, 2 alpha y + by w)
                            let gx = t.add(&t.mul(&two_alpha, &x), &t.mul(&e.k(&c.bx), &w));
                            let gy = t.add(&t.mul(&two_alpha, &y), &t.mul(&e.k(&c.by), &w));
                            // perp(g) = (-gy, gx), times factor; w > 0.
                            if *factor == Sign::Negative {
                                [gy, t.neg(&gx)]
                            } else {
                                [t.neg(&gy), gx]
                            }
                        }
                    }
                };
                let [ux, uy] = vector(u);
                let [vx, vy] = vector(v);
                let t = &e.tower;
                if *cross {
                    t.sub(&t.mul(&ux, &vy), &t.mul(&uy, &vx))
                } else {
                    t.add(&t.mul(&ux, &vx), &t.mul(&uy, &vy))
                }
            }
        };
        e.tower.sign(&value)
    }
}

/// The sign from the points' cached boxes in plain interval arithmetic,
/// or `None`. Sound for the same reason the tower's interval tier is (every
/// box encloses its coordinate, every operation rounds outward), and far
/// cheaper: no tower, no allocation. Decides almost every question whose
/// answer is not zero.
fn quick(pred: &Pred<'_>) -> Option<Sign> {
    let value = match pred {
        Pred::Orient(a, b, c) => {
            let (ux, uy) = (b.bx.sub(&a.bx), b.by.sub(&a.by));
            let (vx, vy) = (c.bx.sub(&a.bx), c.by.sub(&a.by));
            ux.mul(&vy).sub(&uy.mul(&vx))
        }
        Pred::Dot(a, b, dir) => {
            let (dx, dy) = (dir.0.enclosure(), dir.1.enclosure());
            b.bx.sub(&a.bx).mul(&dx).add(&b.by.sub(&a.by).mul(&dy))
        }
        Pred::OnCircle(p, c) => {
            let square = p.bx.mul(&p.bx).add(&p.by.mul(&p.by));
            c.alpha
                .enclosure()
                .mul(&square)
                .add(&c.bx.enclosure().mul(&p.bx))
                .add(&c.by.enclosure().mul(&p.by))
                .add(&c.gamma.enclosure())
        }
        Pred::DiffX(a, b) => a.bx.sub(&b.bx),
        Pred::DiffY(a, b) => a.by.sub(&b.by),
        Pred::Tangents { .. } => return None,
    };
    // A zero from boxes is only a point [0, 0], i.e. two exact equal
    // inputs; still leave zeros to the exact path, which is authoritative.
    match value.sign() {
        Some(Sign::Zero) | None => None,
        decided => decided,
    }
}

/// The proven sign of a predicate.
///
/// # Panics
///
/// Never for points built by this module: every radicand is checked
/// non-negative before a point carrying it is made, so the exact tier
/// always has a real value to decide.
pub(crate) fn sign(pred: Pred<'_>) -> Sign {
    if let Some(sign) = quick(&pred) {
        return sign;
    }
    certify(&pred).expect("predicates over constructed points are always defined")
}

/// Positive when `c` lies left of `a -> b`.
pub(crate) fn orient(a: &XPoint, b: &XPoint, c: &XPoint) -> Sign {
    sign(Pred::Orient(a, b, c))
}

/// Sign of `y(a) - y(b)`, from the boxes when they are apart.
pub(crate) fn cmp_y(a: &XPoint, b: &XPoint) -> Sign {
    box_cmp(a.by, b.by).unwrap_or_else(|| sign(Pred::DiffY(a, b)))
}

fn box_cmp(a: Interval, b: Interval) -> Option<Sign> {
    if a.hi() < b.lo() {
        Some(Sign::Negative)
    } else if b.hi() < a.lo() {
        Some(Sign::Positive)
    } else {
        None
    }
}

/// Exact equality.
pub(crate) fn same_point(a: &XPoint, b: &XPoint) -> bool {
    if a.bx.disjoint(b.bx) || a.by.disjoint(b.by) {
        return false;
    }
    // Same radicand (or none) and proportional coefficients: equal, by
    // plain dyadic cross-multiplication. The common case (a crossing met
    // again from the other edge, a vertex shared by both operands) ends
    // here without building a tower. Not proportional does not prove
    // unequal (a square radicand, say), so that falls through.
    if a.d == b.d {
        let prop = |p: &Dyadic, q: &Dyadic| sgn(&p.mul(&b.w).sub(&q.mul(&a.w))) == Sign::Zero;
        if prop(&a.xa, &b.xa) && prop(&a.xb, &b.xb) && prop(&a.ya, &b.ya) && prop(&a.yb, &b.yb) {
            return true;
        }
        if a.is_rational() {
            // Both rational: proportionality is also necessary.
            return false;
        }
    }
    sign(Pred::DiffX(a, b)) == Sign::Zero && sign(Pred::DiffY(a, b)) == Sign::Zero
}
