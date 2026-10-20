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

use std::sync::{Arc, OnceLock};

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

/// `a / w` (with `w > 0`) correctly rounded to the nearest double, ties to
/// even, starting from a guess within a few units in the last place.
///
/// Each step compares the exact value with the midpoint between the guess
/// and its neighbour, which is a dyadic and so exact: the guess moves
/// towards the value until neither midpoint lies between them.
///
/// The comparison is filtered (#218): `a - r w` is exact, so the box of
/// `(a - r w) / w`, the value's offset from `r`, is tight to a few parts in
/// `2^50` of the offset itself, and it lies clear of the half gap to the
/// neighbour unless the value is that close to a midpoint. Only then is
/// the midpoint built and compared exactly.
pub(crate) fn round_ratio(a: &Dyadic, w: &Dyadic, guess: f64) -> f64 {
    if !guess.is_finite() {
        return guess;
    }
    // Sign of `a / w - (lo + hi) / 2`, exactly.
    let exact = |lo: f64, hi: f64| sgn(&a.sub(&dy(lo).add(&dy(hi)).mul(&dy(0.5)).mul(w)));
    let w_box = w.enclosure();
    // The box of `a / w - r`, for the `r` it was taken at.
    let mut offset: Option<(f64, Interval)> = None;
    let mut beyond = |r: f64, lo: f64, hi: f64| -> Sign {
        // Half gaps are exact doubles away from the subnormal range.
        if r.abs() >= 2f64.powi(-1000) {
            let at = match offset {
                Some((at, b)) if at.to_bits() == r.to_bits() => b,
                _ => {
                    let b = a.sub(&dy(r).mul(w)).enclosure().quotient(w_box);
                    offset = Some((r, b));
                    b
                }
            };
            // `(lo + hi) / 2 - r`, exactly: one of them is `r`.
            let mid = if lo == r {
                0.5 * (hi - r)
            } else {
                -0.5 * (r - lo)
            };
            if let Some(decided @ (Sign::Positive | Sign::Negative)) =
                at.sub(&Interval::point(mid)).sign()
            {
                return decided;
            }
        }
        exact(lo, hi)
    };
    let odd = |r: f64| r.to_bits() & 1 == 1;
    let mut r = guess;
    loop {
        let up = r.next_up();
        if !up.is_finite() {
            break;
        }
        match beyond(r, r, up) {
            Sign::Positive => r = up,
            Sign::Zero if odd(r) => r = up,
            _ => break,
        }
    }
    loop {
        let down = r.next_down();
        if !down.is_finite() {
            break;
        }
        match beyond(r, down, r) {
            Sign::Negative => r = down,
            Sign::Zero if odd(r) => r = down,
            _ => break,
        }
    }
    r
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

/// The exact coefficients of a point; see the module documentation.
#[derive(Debug)]
struct Coeffs {
    xa: Dyadic,
    xb: Dyadic,
    ya: Dyadic,
    yb: Dyadic,
    w: Dyadic,
    /// Zero for rational points, which then embed without a radical.
    d: Dyadic,
}

/// A point whose coefficients wait until a question needs them.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Recipe {
    /// Built eagerly.
    Built,
    /// An input vertex: `(x, y) / 1`.
    Input(Point2),
    /// `from + u * (to - from)` on a straight input edge, `u` a double:
    /// a piece's sample point (#218).
    Along { from: Point2, to: Point2, u: f64 },
    /// The crossing of the lines through two input segments `p0 -> p1`
    /// and `q0 -> q1`, which are not parallel (#218).
    Crossing {
        p0: Point2,
        p1: Point2,
        q0: Point2,
        q1: Point2,
    },
}

/// The coefficients with `w > 0`, and none of a radical a point does not
/// need.
fn normalise(xa: Dyadic, xb: Dyadic, ya: Dyadic, yb: Dyadic, w: Dyadic, d: Dyadic) -> Coeffs {
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
    Coeffs {
        xa,
        xb,
        ya,
        yb,
        w,
        d,
    }
}

/// `(x / w, y / w)`, rational.
fn rational_coeffs(x: Dyadic, y: Dyadic, w: Dyadic) -> Coeffs {
    normalise(x, Dyadic::zero(), y, Dyadic::zero(), w, Dyadic::zero())
}

/// `from + u * (to - from)` in homogeneous dyadics, weight one.
fn along_coeffs(from: Point2, to: Point2, u: f64) -> (Dyadic, Dyadic, Dyadic) {
    let (x0, y0) = (dy(from.x), dy(from.y));
    let u = dy(u);
    (
        x0.add(&u.mul(&dy(to.x).sub(&x0))),
        y0.add(&u.mul(&dy(to.y).sub(&y0))),
        dy(1.0),
    )
}

/// The homogeneous crossing of the lines `p0 + t (p1 - p0)` and
/// `q0 + s (q1 - q0)`: `t = cross(q0 - p0, e) / cross(d, e)`, so the point
/// is `(p0 * den + num * d) / den`. `None` when the lines are parallel.
fn line_crossing(
    p0: Point2,
    p1: Point2,
    q0: Point2,
    q1: Point2,
) -> Option<(Dyadic, Dyadic, Dyadic)> {
    let (p0x, p0y) = (dy(p0.x), dy(p0.y));
    let (q0x, q0y) = (dy(q0.x), dy(q0.y));
    let d = (dy(p1.x).sub(&p0x), dy(p1.y).sub(&p0y));
    let e = (dy(q1.x).sub(&q0x), dy(q1.y).sub(&q0y));
    let den = d.0.mul(&e.1).sub(&d.1.mul(&e.0));
    if sgn(&den) == Sign::Zero {
        return None;
    }
    let (rx, ry) = (q0x.sub(&p0x), q0y.sub(&p0y));
    let num = rx.mul(&e.1).sub(&ry.mul(&e.0));
    Some((
        p0x.mul(&den).add(&num.mul(&d.0)),
        p0y.mul(&den).add(&num.mul(&d.1)),
        den,
    ))
}

/// The coefficients behind a recipe.
fn coeffs_of(recipe: Recipe) -> Coeffs {
    match recipe {
        Recipe::Input(p) => rational_coeffs(dy(p.x), dy(p.y), dy(1.0)),
        Recipe::Along { from, to, u } => {
            let (x, y, w) = along_coeffs(from, to, u);
            rational_coeffs(x, y, w)
        }
        Recipe::Crossing { p0, p1, q0, q1 } => {
            let (x, y, w) =
                line_crossing(p0, p1, q0, q1).expect("a lazy crossing has crossing lines");
            rational_coeffs(x, y, w)
        }
        Recipe::Built => unreachable!("a built point has its coefficients"),
    }
}

/// A point's coefficients and, while they wait, how to build them.
#[derive(Debug)]
struct Shared {
    recipe: Recipe,
    coeffs: OnceLock<Coeffs>,
}

/// An exact point; see the module documentation.
///
/// The coefficients are shared, so a copy of a point costs no arithmetic.
/// An input vertex, a piece's sample on a straight edge and a crossing of
/// two input segments build them only when a question about the point gets
/// past its boxes (#218): most never do.
#[derive(Debug, Clone)]
pub(crate) struct XPoint {
    exact: Arc<Shared>,
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

/// Whether a box from interval arithmetic on doubles can stand in for the
/// one exact coefficients give: finite, and at most `2^-46` of its
/// magnitude wide (a few dozen units in the last place), so most questions
/// end at it and its centre is a close rounding guess.
fn narrow(b: Interval) -> bool {
    let (lo, hi) = (b.lo(), b.hi());
    lo.is_finite() && hi.is_finite() && hi - lo <= lo.abs().max(hi.abs()) * 2f64.powi(-46)
}

/// The centre of a finite box, inside it.
fn centre(b: Interval) -> f64 {
    (0.5 * b.lo() + 0.5 * b.hi()).clamp(b.lo(), b.hi())
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
        let c = normalise(xa, xb, ya, yb, w, d);
        let approx = Point2::new(
            approx(&c.xa, &c.xb, &c.d, &c.w),
            approx(&c.ya, &c.yb, &c.d, &c.w),
        );
        let bx = enclose(&c.xa, &c.xb, &c.d, &c.w);
        let by = enclose(&c.ya, &c.yb, &c.d, &c.w);
        Self {
            bx,
            by,
            exact: Arc::new(Shared {
                recipe: Recipe::Built,
                coeffs: OnceLock::from(c),
            }),
            approx,
        }
    }

    /// A point whose coefficients wait, with sound boxes around it.
    fn lazy(recipe: Recipe, approx: Point2, bx: Interval, by: Interval) -> Self {
        Self {
            exact: Arc::new(Shared {
                recipe,
                coeffs: OnceLock::new(),
            }),
            approx,
            bx,
            by,
        }
    }

    /// The exact coefficients, built on first use for a lazy point.
    fn c(&self) -> &Coeffs {
        self.exact
            .coeffs
            .get_or_init(|| coeffs_of(self.exact.recipe))
    }

    /// `(x / w, y / w)`.
    pub(crate) fn rational(x: Dyadic, y: Dyadic, w: Dyadic) -> Self {
        Self::new(x, Dyadic::zero(), y, Dyadic::zero(), w, Dyadic::zero())
    }

    /// An input vertex, exactly.
    pub(crate) fn from_f64(p: Point2) -> Self {
        // Exact boxes, and single-point boxes let equality of input
        // vertices be read off without arithmetic (see `same_point`). The
        // coefficients `(x, y) / 1` wait until a question needs them; the
        // approximation is the vertex itself, `-0.0` read as `0.0` as the
        // coefficients would give it.
        Self::lazy(
            Recipe::Input(p),
            Point2::new(p.x + 0.0, p.y + 0.0),
            Interval::point(p.x),
            Interval::point(p.y),
        )
    }

    /// `from + u * (to - from)`, exactly, for doubles (#218): the point
    /// [`XPoint::rational`] gives for that expression, with its
    /// coefficients waiting behind a box from interval arithmetic. A box
    /// too wide to stand in builds them at once.
    pub(crate) fn along(from: Point2, to: Point2, u: f64) -> Self {
        let at = |a: f64, b: f64| {
            let a_ = Interval::point(a);
            a_.add(&Interval::point(u).mul(&Interval::point(b).sub(&a_)))
        };
        let (bx, by) = (at(from.x, to.x), at(from.y, to.y));
        if !(narrow(bx) && narrow(by)) {
            let (x, y, w) = along_coeffs(from, to, u);
            return Self::rational(x, y, w);
        }
        Self::lazy(
            Recipe::Along { from, to, u },
            Point2::new(centre(bx), centre(by)),
            bx,
            by,
        )
    }

    /// The crossing of the lines through the input segments `p0 -> p1`
    /// and `q0 -> q1`, or `None` when they are parallel.
    ///
    /// Its boxes come from interval arithmetic on the doubles, and the
    /// exact coefficients wait until a question gets past them (#218). A
    /// box too wide to stand in (nearly parallel lines) builds them at
    /// once, so the point never decides less from its boxes than an eager
    /// one would by much, and its approximation, the boxes' centre, is
    /// always a close guess for rounding.
    pub(crate) fn crossing(p0: Point2, p1: Point2, q0: Point2, q1: Point2) -> Option<Self> {
        let i = Interval::point;
        let (dx, dy_) = (i(p1.x).sub(&i(p0.x)), i(p1.y).sub(&i(p0.y)));
        let (ex, ey) = (i(q1.x).sub(&i(q0.x)), i(q1.y).sub(&i(q0.y)));
        let den = dx.mul(&ey).sub(&dy_.mul(&ex));
        let (rx, ry) = (i(q0.x).sub(&i(p0.x)), i(q0.y).sub(&i(p0.y)));
        let num = rx.mul(&ey).sub(&ry.mul(&ex));
        let t = num.quotient(den);
        let (bx, by) = (i(p0.x).add(&t.mul(&dx)), i(p0.y).add(&t.mul(&dy_)));
        if narrow(bx) && narrow(by) {
            return Some(Self::lazy(
                Recipe::Crossing { p0, p1, q0, q1 },
                Point2::new(centre(bx), centre(by)),
                bx,
                by,
            ));
        }
        let (x, y, w) = line_crossing(p0, p1, q0, q1)?;
        Some(Self::rational(x, y, w))
    }

    /// True when the point carries no square root.
    pub(crate) fn is_rational(&self) -> bool {
        self.exact.recipe != Recipe::Built || sgn(&self.c().d) == Sign::Zero
    }

    /// A nearby double pair, for output only; never for decisions.
    pub(crate) fn approx(&self) -> Point2 {
        self.approx
    }

    /// The output position: each coordinate of a rational point correctly
    /// rounded to nearest (ties to even), so an input vertex comes back
    /// bit-identical and a crossing of two segments is the double nearest
    /// to the true crossing (#173). A point carrying a square root keeps
    /// [`XPoint::approx`].
    pub(crate) fn rounded(&self) -> Point2 {
        if matches!(self.exact.recipe, Recipe::Input(_)) || !self.is_rational() {
            // An input vertex is a double already.
            return self.approx;
        }
        let c = self.c();
        Point2::new(
            round_ratio(&c.xa, &c.w, self.approx.x),
            round_ratio(&c.ya, &c.w, self.approx.y),
        )
    }

    /// The identity of this point's shared coefficients: equal for copies
    /// of one point, and an equal identity proves equal points (#218).
    /// Points built separately may be equal with distinct identities.
    pub(crate) fn id(&self) -> usize {
        Arc::as_ptr(&self.exact) as usize
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
    Dot(&'a XPoint, &'a XPoint, &'a Direction),
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

/// A constant direction `(x, y)`, with sound boxes around its components.
#[derive(Debug)]
pub(crate) struct Direction {
    pub(crate) x: Dyadic,
    pub(crate) y: Dyadic,
    bx: Interval,
    by: Interval,
}

impl Direction {
    pub(crate) fn new(x: Dyadic, y: Dyadic) -> Self {
        Self {
            bx: x.enclosure(),
            by: y.enclosure(),
            x,
            y,
        }
    }

    /// The opposite direction.
    pub(crate) fn reversed(&self) -> Self {
        Self {
            x: self.x.neg(),
            y: self.y.neg(),
            bx: self.bx.neg(),
            by: self.by.neg(),
        }
    }
}

/// The travel direction of a curve, as a function of the point.
#[derive(Debug, Clone)]
pub(crate) enum Tangent {
    /// A constant direction (a segment traversed some way), shared by
    /// every piece of the segment (#218).
    Fixed(Arc<Direction>),
    /// `factor * perp(grad F)` for a circle `F`: `factor` positive is
    /// counter-clockwise travel.
    Circle(Arc<Circle>, Sign),
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
        let rational = p.is_rational();
        let p = p.c();
        let w = self.k(&p.w);
        if rational {
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
                let (dx, dy) = (e.k(&dir.x), e.k(&dir.y));
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
                // Only a circle's tangent depends on the point; two fixed
                // directions leave its coefficients unbuilt (#218).
                let on_circle = |tangent: &Tangent| matches!(tangent, Tangent::Circle(..));
                let [x, y, w] = if on_circle(u) || on_circle(v) {
                    e.point(at)?
                } else {
                    let one = e.k(&dy(1.0));
                    [one.clone(), one.clone(), one]
                };
                let vector = |tangent: &Tangent| -> [Nested<T>; 2] {
                    match tangent {
                        Tangent::Fixed(d) => [e.k(&d.x), e.k(&d.y)],
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
            let (dx, dy) = (dir.bx, dir.by);
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
        Pred::Tangents {
            u: Tangent::Fixed(u),
            v: Tangent::Fixed(v),
            cross,
            ..
        } => {
            if *cross {
                u.bx.mul(&v.by).sub(&u.by.mul(&v.bx))
            } else {
                u.bx.mul(&v.bx).add(&u.by.mul(&v.by))
            }
        }
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

/// Whether two points are one, when that is read off without arithmetic:
/// copies of one point, or two points whose boxes are the same single
/// point. `false` proves nothing.
pub(crate) fn plainly_same(a: &XPoint, b: &XPoint) -> bool {
    let exact = |p: &XPoint| p.bx.lo() == p.bx.hi() && p.by.lo() == p.by.hi();
    Arc::ptr_eq(&a.exact, &b.exact)
        || (exact(a) && exact(b) && a.bx.lo() == b.bx.lo() && a.by.lo() == b.by.lo())
}

/// Exact equality.
pub(crate) fn same_point(a: &XPoint, b: &XPoint) -> bool {
    if a.bx.disjoint(b.bx) || a.by.disjoint(b.by) {
        return false;
    }
    // Copies of one point share its coefficients (#218): a crossing is
    // built once and handed to both edges, so most pairs asked here end.
    if Arc::ptr_eq(&a.exact, &b.exact) {
        return true;
    }
    // Boxes that are single points hold their coordinate exactly: input
    // vertices, and every point equal to one.
    let exact = |p: &XPoint| p.bx.lo() == p.bx.hi() && p.by.lo() == p.by.hi();
    if exact(a) && exact(b) {
        return a.bx.lo() == b.bx.lo() && a.by.lo() == b.by.lo();
    }
    // Same radicand (or none) and proportional coefficients: equal, by
    // plain dyadic cross-multiplication. The common case (a crossing met
    // again from the other edge, a vertex shared by both operands) ends
    // here without building a tower. Not proportional does not prove
    // unequal (a square radicand, say), so that falls through.
    let rational = a.is_rational();
    let (a_, b_) = (a.c(), b.c());
    if a_.d == b_.d {
        let (a, b) = (a_, b_);
        let prop = |p: &Dyadic, q: &Dyadic| sgn(&p.mul(&b.w).sub(&q.mul(&a.w))) == Sign::Zero;
        if prop(&a.xa, &b.xa) && prop(&a.xb, &b.xb) && prop(&a.ya, &b.ya) && prop(&a.yb, &b.yb) {
            return true;
        }
        if rational {
            // Both rational: proportionality is also necessary.
            return false;
        }
    }
    sign(Pred::DiffX(a, b)) == Sign::Zero && sign(Pred::DiffY(a, b)) == Sign::Zero
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sign of `num / den - (a + b) / 2`, exactly.
    fn beyond(num: &Dyadic, den: &Dyadic, a: f64, b: f64) -> Sign {
        sgn(&num.sub(&dy(a).add(&dy(b)).mul(&dy(0.5)).mul(den)))
    }

    #[test]
    fn rational_points_round_to_the_nearest_double() {
        // Crossings of segments between random doubles: `rounded` must be
        // correctly rounded, and must differ from the quick `approx` for
        // some, or this test would not see a missing correction.
        let mut state = 0x173u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64 * 10.0 - 5.0
        };
        let mut corrected = 0;
        for _ in 0..20_000 {
            let (num, den) = (dy(next()).mul(&dy(next())), dy(next()).mul(&dy(next())));
            let (num, den) = if sgn(&den) == Sign::Negative {
                (num.neg(), den.neg())
            } else {
                (num, den)
            };
            if sgn(&den) == Sign::Zero {
                continue;
            }
            let point = XPoint::rational(num.clone(), num.clone(), den.clone());
            let value = point.rounded().x;
            assert_ne!(
                beyond(&num, &den, value.next_down(), value),
                Sign::Negative,
                "{value}"
            );
            assert_ne!(
                beyond(&num, &den, value, value.next_up()),
                Sign::Positive,
                "{value}"
            );
            if value != point.approx().x {
                corrected += 1;
            }
        }
        assert!(corrected > 0, "approx was always right: the test is blind");
    }
}
