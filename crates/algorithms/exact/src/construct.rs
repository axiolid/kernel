//! Exact constructions over `f64` inputs, decided without division.
//!
//! Each construction's result is kept symbolically -- as its inputs plus a
//! recipe -- and every question about it is a [`SignExpr`], answered by the
//! interval filter or, when that cannot decide, exactly. Nothing is rounded
//! until a caller asks for an approximate value for output.

use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Sign};

use crate::arith::{sign_product, Arith};
use crate::certify::{certify, filter, require_finite, ExactError, SignExpr};
use crate::root::{sign_root, Root2};

/// The line through two distinct points, parameterised `from + t*(to - from)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line {
    from: Point2,
    to: Point2,
}

impl Line {
    /// The line through `from` and `to`.
    pub fn new(from: Point2, to: Point2) -> Result<Self, ExactError> {
        require_finite(&[from.x, from.y, to.x, to.y])?;
        // Exact comparison on purpose: any two distinct floats define a line.
        if from.x == to.x && from.y == to.y {
            return Err(ExactError::DegenerateLine);
        }
        Ok(Self { from, to })
    }

    /// Parameter 0.
    #[must_use]
    pub const fn from(self) -> Point2 {
        self.from
    }

    /// Parameter 1.
    #[must_use]
    pub const fn to(self) -> Point2 {
        self.to
    }
}

/// A circle; a zero radius is allowed and behaves as a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle {
    centre: Point2,
    radius: f64,
}

impl Circle {
    /// The circle about `centre` with `radius`.
    pub fn new(centre: Point2, radius: f64) -> Result<Self, ExactError> {
        require_finite(&[centre.x, centre.y, radius])?;
        if radius < 0.0 {
            return Err(ExactError::NegativeRadius);
        }
        Ok(Self { centre, radius })
    }

    /// Centre.
    #[must_use]
    pub const fn centre(self) -> Point2 {
        self.centre
    }

    /// Radius.
    #[must_use]
    pub const fn radius(self) -> f64 {
        self.radius
    }
}

fn v<T: Arith>(value: f64) -> T {
    T::from_f64(value)
}

/// `orient2d(a, b, p)` as a polynomial.
fn orient<T: Arith>(a: Point2, b: Point2, px: &T, py: &T) -> T {
    let (ax, ay) = (v::<T>(a.x), v::<T>(a.y));
    let bax = v::<T>(b.x).sub(&ax);
    let bay = v::<T>(b.y).sub(&ay);
    bax.mul(&py.sub(&ay)).sub(&bay.mul(&px.sub(&ax)))
}

// ---------------------------------------------------------------- crossings

/// Where lines `first` and `second` cross, `(Nx/D, Ny/D)`, homogeneously.
fn crossing<T: Arith>(first: Line, second: Line) -> (T, T, T) {
    let (x1, y1) = (v::<T>(first.from.x), v::<T>(first.from.y));
    let (x2, y2) = (v::<T>(first.to.x), v::<T>(first.to.y));
    let (x3, y3) = (v::<T>(second.from.x), v::<T>(second.from.y));
    let (x4, y4) = (v::<T>(second.to.x), v::<T>(second.to.y));
    let d = x1
        .sub(&x2)
        .mul(&y3.sub(&y4))
        .sub(&y1.sub(&y2).mul(&x3.sub(&x4)));
    let c1 = x1.mul(&y2).sub(&y1.mul(&x2));
    let c2 = x3.mul(&y4).sub(&y3.mul(&x4));
    let nx = c1.mul(&x3.sub(&x4)).sub(&x1.sub(&x2).mul(&c2));
    let ny = c1.mul(&y3.sub(&y4)).sub(&y1.sub(&y2).mul(&c2));
    (nx, ny, d)
}

struct CrossingDenominator(Line, Line);

impl SignExpr for CrossingDenominator {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        crossing::<T>(self.0, self.1).2.sign()
    }
}

struct CrossingOrientation {
    first: Line,
    second: Line,
    a: Point2,
    b: Point2,
}

impl SignExpr for CrossingOrientation {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let (nx, ny, d) = crossing::<T>(self.first, self.second);
        // orient(a, b, N/D) * D, which has no division: substitute
        // p = N/D into orient and multiply through by D.
        let (ax, ay) = (v::<T>(self.a.x), v::<T>(self.a.y));
        let bax = v::<T>(self.b.x).sub(&ax);
        let bay = v::<T>(self.b.y).sub(&ay);
        let scaled = bax
            .mul(&ny.sub(&ay.mul(&d)))
            .sub(&bay.mul(&nx.sub(&ax.mul(&d))));
        // orient = scaled / D, so its sign is sign(scaled) * sign(D).
        Some(sign_product(scaled.sign()?, d.sign()?))
    }
}

/// Which side of the directed line `a -> b` the crossing of `first` and
/// `second` lies on: [`Sign::Positive`] for left.
///
/// Returns `Ok(None)` when the lines are parallel (no single crossing).
/// The crossing point is never rounded, so a crossing that lies exactly on
/// `a -> b` reports [`Sign::Zero`].
pub fn crossing_orientation(
    first: Line,
    second: Line,
    a: Point2,
    b: Point2,
) -> Result<Option<Sign>, ExactError> {
    require_finite(&[a.x, a.y, b.x, b.y])?;
    if certify(&CrossingDenominator(first, second))? == Sign::Zero {
        return Ok(None);
    }
    certify(&CrossingOrientation {
        first,
        second,
        a,
        b,
    })
    .map(Some)
}

/// The interval filter alone for [`crossing_orientation`], so the
/// escalation rate can be observed and tested. [`Certified::Uncertain`]
/// means the exact tier would run; parallel lines also report uncertain.
pub fn crossing_orientation_filter(first: Line, second: Line, a: Point2, b: Point2) -> Certified {
    if !a.x.is_finite() || !a.y.is_finite() || !b.x.is_finite() || !b.y.is_finite() {
        return filter(&NeverDecides);
    }
    filter(&CrossingOrientation {
        first,
        second,
        a,
        b,
    })
}

/// An expression no arithmetic decides: the filter's "uncertain" for input
/// the exact tier would refuse.
struct NeverDecides;

impl SignExpr for NeverDecides {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        None
    }
}

// ------------------------------------------------------- line meets circle

/// Which of the two solutions: `Minus` has the smaller parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Branch {
    /// `t = (-B - sqrt(disc)) / A`: the first hit along the line.
    Minus,
    /// `t = (-B + sqrt(disc)) / A`: the second hit along the line.
    Plus,
}

/// A point where a line meets a circle, held exactly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineHit {
    line: Line,
    circle: Circle,
    branch: Branch,
}

/// How a line meets a circle, hits in increasing parameter order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HitCount {
    /// The line misses the circle.
    Missed,
    /// The line touches the circle at exactly one point.
    Tangent(LineHit),
    /// The line crosses the circle at two points.
    Secant(LineHit, LineHit),
}

/// `A t^2 + 2 B t + C = 0` for the line's parameter, with
/// `disc = B^2 - A*C`, so `t = (-B +- sqrt(disc)) / A`.
fn quadratic<T: Arith>(line: Line, circle: Circle) -> (T, T, T) {
    let dx = v::<T>(line.to.x).sub(&v(line.from.x));
    let dy = v::<T>(line.to.y).sub(&v(line.from.y));
    let fx = v::<T>(line.from.x).sub(&v(circle.centre.x));
    let fy = v::<T>(line.from.y).sub(&v(circle.centre.y));
    let a = dx.square().add(&dy.square());
    let b = dx.mul(&fx).add(&dy.mul(&fy));
    let c = fx
        .square()
        .add(&fy.square())
        .sub(&v::<T>(circle.radius).square());
    let disc = b.square().sub(&a.mul(&c));
    (a, b, disc)
}

struct Discriminant(Line, Circle);

impl SignExpr for Discriminant {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        quadratic::<T>(self.0, self.1).2.sign()
    }
}

/// Where `line` meets `circle`, exactly: missed, tangent, or two hits.
///
/// Tangency is decided exactly, so a line that grazes the circle is never
/// reported as two nearly equal hits, or as a miss, by rounding.
pub fn line_circle_hits(line: Line, circle: Circle) -> Result<HitCount, ExactError> {
    let hit = |branch| LineHit {
        line,
        circle,
        branch,
    };
    Ok(match certify(&Discriminant(line, circle))? {
        Sign::Negative => HitCount::Missed,
        Sign::Zero => HitCount::Tangent(hit(Branch::Minus)),
        _ => HitCount::Secant(hit(Branch::Minus), hit(Branch::Plus)),
    })
}

impl LineHit {
    /// The line this hit lies on.
    #[must_use]
    pub const fn line(self) -> Line {
        self.line
    }

    /// The circle this hit lies on.
    #[must_use]
    pub const fn circle(self) -> Circle {
        self.circle
    }

    /// Which solution of the quadratic this is.
    #[must_use]
    pub const fn branch(self) -> Branch {
        self.branch
    }

    fn root_sign<T: Arith>(self) -> T {
        match self.branch {
            Branch::Minus => v::<T>(-1.0),
            Branch::Plus => v::<T>(1.0),
        }
    }

    /// The parameter as `(a + b*sqrt(c)) / d` in arithmetic `T`.
    fn parameter<T: Arith>(self) -> Root2<T> {
        let (a, b, disc) = quadratic::<T>(self.line, self.circle);
        Root2 {
            a: b.neg(),
            b: self.root_sign(),
            c: disc,
            d: a,
        }
    }

    /// Sign of `t - value`: whether this hit comes before (`Negative`),
    /// at, or after the line parameter `value`. With `0.0` and `1.0` this
    /// answers "is the hit inside the segment `from..to`".
    pub fn cmp_param(self, value: f64) -> Result<Sign, ExactError> {
        require_finite(&[value])?;
        certify(&HitVersusParam { hit: self, value })
    }

    /// Which side of the directed line `a -> b` this hit lies on.
    pub fn orientation(self, a: Point2, b: Point2) -> Result<Sign, ExactError> {
        require_finite(&[a.x, a.y, b.x, b.y])?;
        certify(&HitOrientation { hit: self, a, b })
    }

    /// The parameter rounded to `f64`, for output only.
    ///
    /// Never use this to make a decision; ask a sign question instead.
    #[must_use]
    pub fn approx_param(self) -> f64 {
        let (line, circle) = (self.line, self.circle);
        let (dx, dy) = (line.to.x - line.from.x, line.to.y - line.from.y);
        let (fx, fy) = (line.from.x - circle.centre.x, line.from.y - circle.centre.y);
        let a = dx * dx + dy * dy;
        let b = dx * fx + dy * fy;
        let c = fx * fx + fy * fy - circle.radius * circle.radius;
        let root = (b * b - a * c).max(0.0).sqrt();
        match self.branch {
            Branch::Minus => (-b - root) / a,
            Branch::Plus => (-b + root) / a,
        }
    }

    /// The point rounded to `f64`, for output only.
    #[must_use]
    pub fn approx_point(self) -> Point2 {
        let t = self.approx_param();
        let (from, to) = (self.line.from, self.line.to);
        Point2::new(from.x + t * (to.x - from.x), from.y + t * (to.y - from.y))
    }
}

struct HitVersusParam {
    hit: LineHit,
    value: f64,
}

impl SignExpr for HitVersusParam {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let mut root = self.hit.parameter::<T>();
        // t - value = (a - value*d + b*sqrt(c)) / d.
        root.a = root.a.sub(&v::<T>(self.value).mul(&root.d));
        root.sign()
    }
}

struct HitOrientation {
    hit: LineHit,
    a: Point2,
    b: Point2,
}

impl SignExpr for HitOrientation {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        // orient(a, b, from + t*dir) is linear in t: k0 + k1*t, with
        // k0 = orient(a, b, from) and k1 = cross(b - a, dir). Substituting
        // t = (-B + s*sqrt(disc)) / A and multiplying by A > 0:
        // (k0*A - k1*B) + s*k1*sqrt(disc).
        let line = self.hit.line;
        let (qa, qb, disc) = quadratic::<T>(line, self.hit.circle);
        let k0 = orient(self.a, self.b, &v::<T>(line.from.x), &v::<T>(line.from.y));
        let dx = v::<T>(line.to.x).sub(&v(line.from.x));
        let dy = v::<T>(line.to.y).sub(&v(line.from.y));
        let bax = v::<T>(self.b.x).sub(&v(self.a.x));
        let bay = v::<T>(self.b.y).sub(&v(self.a.y));
        let k1 = bax.mul(&dy).sub(&bay.mul(&dx));
        let rational = k0.mul(&qa).sub(&k1.mul(&qb));
        let irrational = self.hit.root_sign::<T>().mul(&k1);
        sign_root(&rational, &irrational, &disc)
    }
}

// ------------------------------------------------------ ordering along a line

struct HitOrder(LineHit, LineHit);

impl SignExpr for HitOrder {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        self.0.parameter::<T>().cmp_sign(&self.1.parameter::<T>())
    }
}

/// Sign of `t(first) - t(second)`: which hit comes first along their line.
///
/// Hits on different circles have different radicands; the comparison is
/// still exact ([`crate::sign_two_roots`]). Both hits must lie on the same
/// line, since parameters on different lines are not comparable.
pub fn compare_along(first: LineHit, second: LineHit) -> Result<Sign, ExactError> {
    if first.line != second.line {
        return Err(ExactError::DifferentLines);
    }
    if first.circle == second.circle {
        // Same quadratic: t = (-B -+ sqrt(disc)) / A with A = |to - from|^2
        // > 0, so Minus < Plus whenever both exist (disc > 0). Deciding this
        // structurally matters: evaluating it would subtract two equal
        // intervals, which never certify zero, and escalate every time.
        return Ok(match (first.branch, second.branch) {
            (Branch::Minus, Branch::Plus) => Sign::Negative,
            (Branch::Plus, Branch::Minus) => Sign::Positive,
            _ => Sign::Zero,
        });
    }
    certify(&HitOrder(first, second))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The same-circle shortcut in `compare_along` must agree with evaluating
    /// the comparison, which it exists to skip.
    #[test]
    fn same_circle_shortcut_agrees_with_evaluation() {
        let mut state = 0xC0FF_EE00_1234_5678u64;
        let mut f = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            ((state >> 11) as f64 / (1u64 << 53) as f64) * 200.0 - 100.0
        };
        let mut checked = 0;
        while checked < 500 {
            let line = Line::new(Point2::new(f(), f()), Point2::new(f(), f()));
            let circle = Circle::new(Point2::new(f(), f()), f().abs());
            let (Ok(line), Ok(circle)) = (line, circle) else {
                continue;
            };
            if let Ok(HitCount::Secant(a, b)) = line_circle_hits(line, circle) {
                for (x, y) in [(a, b), (b, a), (a, a), (b, b)] {
                    assert_eq!(compare_along(x, y), certify(&HitOrder(x, y)));
                }
                checked += 1;
            }
        }
    }
}
