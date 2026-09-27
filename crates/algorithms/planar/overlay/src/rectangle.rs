//! The minimum-area rectangle enclosing a point set (#182).
//!
//! # Exact choice, rounded output
//!
//! Some rectangle of least area has a side on an edge of the convex hull
//! (Freeman and Shapira), so rotating calipers try each hull edge `d` in
//! turn. Every decision is exact:
//!
//! - the hull comes from exact orientations;
//! - the calipers advance while the next vertex projects no less far along
//!   `d` (or along `d` turned a quarter) -- the sign of a dot product of
//!   differences of `f64`s, decided in intervals and else in dyadics;
//! - an orientation's area is `W H / |d|^2`, with `W` and `H` the spans of
//!   those projections, so two orientations compare exactly by
//!   `W_i H_i |d_j|^2` against `W_j H_j |d_i|^2`;
//! - ties (a square has two minimal orientations) are found exactly and
//!   broken by the least angle of the first axis, turned by quarter turns
//!   into `[0, 90)` degrees -- itself an exact comparison, so the answer
//!   does not depend on the order of the input.
//!
//! Only the output is rounded: the unit axes (a square root), the centre
//! and the half extents. [`RectangleEvidence::error`] bounds how far any
//! of them, and any corner, lies from the exact rectangle.

use axiolid_core::{Point2, Vec2};
use axiolid_exact::{certify, Arith, Dyadic, SignExpr};
use axiolid_guarantees::Sign;

/// A rectangle by its centre, two unit axes and the half extents along
/// them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientedRectangle {
    /// Centre.
    pub centre: Point2,
    /// Unit axes, counter-clockwise: the second is the first turned a
    /// quarter. The first points into `[0, 90)` degrees.
    pub axes: [Vec2; 2],
    /// Half the side lengths along the axes; one is zero for collinear
    /// input, both for a single point.
    pub half_extents: [f64; 2],
}

impl OrientedRectangle {
    /// Its area.
    #[must_use]
    pub fn area(&self) -> f64 {
        4.0 * self.half_extents[0] * self.half_extents[1]
    }

    /// The corners, counter-clockwise.
    #[must_use]
    pub fn corners(&self) -> [Point2; 4] {
        let u = self.axes[0] * self.half_extents[0];
        let v = self.axes[1] * self.half_extents[1];
        let c = self.centre;
        [c - u - v, c + u - v, c + u + v, c - u + v]
    }
}

/// What the construction did and how exact its output is.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct RectangleEvidence {
    /// Vertices of the convex hull.
    pub hull_vertices: usize,
    /// Distinct orientations of least area (a square has two, a generic
    /// set one). The one returned turns its first axis least from the
    /// x-axis.
    pub minimal_orientations: usize,
    /// A bound on the distance between any output point (the centre, a
    /// corner) or length (a half extent) and its exact value: a few ulps
    /// of the largest coordinate, from normalising the axes and projecting
    /// onto them. Zero when nothing needed rounding.
    pub error: f64,
}

/// The rectangle and its evidence.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct MinimumRectangle {
    /// The rectangle.
    pub rectangle: OrientedRectangle,
    /// What was rounded, and by how much at most.
    pub evidence: RectangleEvidence,
}

/// Why no rectangle was built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RectangleError {
    /// No points.
    Empty,
    /// A coordinate was not finite.
    NonFinite,
}

/// `(b - a) . d`, or `(b - a) . d⊥` with `d⊥` the quarter turn of `d`.
struct Along {
    a: Point2,
    b: Point2,
    d: Vec2,
    across: bool,
}

impl SignExpr for Along {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let (dx, dy) = if self.across {
            (f(self.d.y).neg(), f(self.d.x))
        } else {
            (f(self.d.x), f(self.d.y))
        };
        let ex = f(self.b.x).sub(&f(self.a.x));
        let ey = f(self.b.y).sub(&f(self.a.y));
        ex.mul(&dx).add(&ey.mul(&dy)).sign()
    }
}

/// Orientation of `c` against `a -> b`.
struct Orient {
    a: Point2,
    b: Point2,
    c: Point2,
}

impl SignExpr for Orient {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let (ux, uy) = (f(self.b.x).sub(&f(self.a.x)), f(self.b.y).sub(&f(self.a.y)));
        let (vx, vy) = (f(self.c.x).sub(&f(self.a.x)), f(self.c.y).sub(&f(self.a.y)));
        ux.mul(&vy).sub(&uy.mul(&vx)).sign()
    }
}

/// The inputs are finite, so the exact tier always decides.
fn sign<E: SignExpr>(e: &E) -> Sign {
    certify(e).unwrap_or(Sign::Zero)
}

/// The convex hull, counter-clockwise, without collinear vertices; the
/// two extremes for collinear input, one point for coincident input.
fn hull(points: &[Point2]) -> Vec<Point2> {
    let mut p = points.to_vec();
    p.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    p.dedup();
    if p.len() < 3 {
        return p;
    }
    let chain = |order: &mut dyn Iterator<Item = Point2>| {
        let mut out: Vec<Point2> = Vec::new();
        for c in order {
            while let [.., a, b] = out[..] {
                if sign(&Orient { a, b, c }) == Sign::Positive {
                    break;
                }
                out.pop();
            }
            out.push(c);
        }
        out.pop();
        out
    };
    let mut lower = chain(&mut p.iter().copied());
    lower.extend(chain(&mut p.iter().rev().copied()));
    lower
}

fn exact(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

/// A direction turned by quarter turns into `[0, 90)` degrees; exact.
fn canonical(mut d: Vec2) -> Vec2 {
    while !(d.x > 0.0 && d.y >= 0.0) {
        d = Vec2::new(d.y, -d.x);
    }
    d
}

/// Whether `a` turns strictly clockwise into `b`, exactly.
fn clockwise(a: Vec2, b: Vec2) -> bool {
    exact(a.x)
        .mul(&exact(b.y))
        .sub(&exact(a.y).mul(&exact(b.x)))
        .sign()
        == Some(Sign::Negative)
}

/// One orientation: a hull edge and the calipers' four vertices.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    d: Vec2,
    /// Hull vertices least and most along `d`, and most across it; the
    /// edge's own start is least across it.
    lo: usize,
    hi: usize,
    base: usize,
    top: usize,
}

impl Candidate {
    /// `W H` and `|d|^2`, exactly.
    fn area_parts(&self, h: &[Point2]) -> (Dyadic, Dyadic) {
        let (dx, dy) = (exact(self.d.x), exact(self.d.y));
        let span = |a: Point2, b: Point2, across: bool| {
            let (ex, ey) = (exact(b.x).sub(&exact(a.x)), exact(b.y).sub(&exact(a.y)));
            if across {
                ey.mul(&dx).sub(&ex.mul(&dy))
            } else {
                ex.mul(&dx).add(&ey.mul(&dy))
            }
        };
        let w = span(h[self.lo], h[self.hi], false);
        let t = span(h[self.base], h[self.top], true);
        (w.mul(&t), dx.mul(&dx).add(&dy.mul(&dy)))
    }
}

/// The minimum-area rectangle enclosing `points`.
///
/// # Errors
///
/// [`RectangleError::Empty`] for no points, [`RectangleError::NonFinite`]
/// for a coordinate that is not finite.
pub fn minimum_area_rectangle(points: &[Point2]) -> Result<MinimumRectangle, RectangleError> {
    if points.is_empty() {
        return Err(RectangleError::Empty);
    }
    if !points.iter().all(|p| p.is_finite()) {
        return Err(RectangleError::NonFinite);
    }
    let h = hull(points);
    let size = h
        .iter()
        .fold(0.0f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));
    let evidence = |minimal, error| RectangleEvidence {
        hull_vertices: h.len(),
        minimal_orientations: minimal,
        error,
    };
    let (rectangle, minimal) = match h.len() {
        1 => {
            let rectangle = OrientedRectangle {
                centre: h[0],
                axes: [Vec2::X, Vec2::Y],
                half_extents: [0.0, 0.0],
            };
            return Ok(MinimumRectangle {
                rectangle,
                evidence: evidence(1, 0.0),
            });
        }
        2 => {
            // Nothing lies across the segment. An even number of quarter
            // turns leaves the first axis along it, an odd number the second.
            let d = h[1] - h[0];
            let axis = canonical(d);
            let mut rectangle = fit(&h, axis);
            rectangle.half_extents[usize::from(axis == d || axis == -d)] = 0.0;
            rectangle.centre = Point2::new(0.5 * (h[0].x + h[1].x), 0.5 * (h[0].y + h[1].y));
            (rectangle, 1)
        }
        _ => {
            let (best, minimal) = calipers(&h);
            (fit(&h, canonical(best.d)), minimal)
        }
    };
    Ok(MinimumRectangle {
        rectangle,
        evidence: evidence(minimal, 64.0 * f64::EPSILON * size),
    })
}

/// The orientation of least area over a hull of three or more vertices,
/// and how many distinct orientations share it.
fn calipers(h: &[Point2]) -> (Candidate, usize) {
    let n = h.len();
    // Whether the vertex after `at` lies no less far (or, `larger` false,
    // no further) along `d` or across it.
    let step = |at: usize, d: Vec2, across: bool, larger: bool| {
        let s = sign(&Along {
            a: h[at],
            b: h[(at + 1) % n],
            d,
            across,
        });
        if larger {
            s != Sign::Negative
        } else {
            s != Sign::Positive
        }
    };
    let advance = |mut at: usize, d: Vec2, across: bool, larger: bool| {
        // A hull turns left, so a caliper only moves forward, and never
        // round the whole hull.
        for _ in 0..n {
            if !step(at, d, across, larger) {
                break;
            }
            at = (at + 1) % n;
        }
        at
    };
    let mut candidates = Vec::with_capacity(n);
    let (mut hi, mut top, mut lo) = (0, 0, 0);
    for base in 0..n {
        let d = h[(base + 1) % n] - h[base];
        if base == 0 {
            // First placement: each caliper starts where the previous one
            // stopped, walking round from the edge.
            hi = advance(0, d, false, true);
            top = advance(hi, d, true, true);
            lo = advance(top, d, false, false);
        } else {
            hi = advance(hi, d, false, true);
            top = advance(top, d, true, true);
            lo = advance(lo, d, false, false);
        }
        candidates.push(Candidate {
            d,
            lo,
            hi,
            base,
            top,
        });
    }
    let mut best = candidates[0];
    let mut best_parts = best.area_parts(h);
    let mut ties = vec![canonical(best.d)];
    for c in &candidates[1..] {
        let parts = c.area_parts(h);
        let order = parts
            .0
            .mul(&best_parts.1)
            .sub(&best_parts.0.mul(&parts.1))
            .sign();
        match order {
            Some(Sign::Negative) => {
                best = *c;
                best_parts = parts;
                ties = vec![canonical(c.d)];
            }
            Some(Sign::Zero) => {
                let axis = canonical(c.d);
                if clockwise(canonical(best.d), axis) {
                    best = *c;
                    best_parts = parts;
                }
                if !ties
                    .iter()
                    .any(|t| !clockwise(*t, axis) && !clockwise(axis, *t))
                {
                    ties.push(axis);
                }
            }
            _ => {}
        }
    }
    (best, ties.len())
}

/// The rectangle with first axis along `d` that encloses the hull,
/// rounded once per quantity.
fn fit(h: &[Point2], d: Vec2) -> OrientedRectangle {
    let l = d.x.hypot(d.y);
    let u = Vec2::new(d.x / l, d.y / l);
    let v = Vec2::new(-u.y, u.x);
    let span = |axis: Vec2| {
        h.iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                let t = p.x * axis.x + p.y * axis.y;
                (lo.min(t), hi.max(t))
            })
    };
    let ((u0, u1), (v0, v1)) = (span(u), span(v));
    let (mu, mv) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
    OrientedRectangle {
        centre: Point2::new(u.x * mu + v.x * mv, u.y * mu + v.y * mv),
        axes: [u, v],
        half_extents: [0.5 * (u1 - u0), 0.5 * (v1 - v0)],
    }
}
