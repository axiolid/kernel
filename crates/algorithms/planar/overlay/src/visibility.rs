//! The visibility polygon of a point inside a region with holes (#184).
//!
//! # A sweep with exact decisions
//!
//! Around the viewpoint `v` every boundary vertex has a direction. Between
//! two consecutive directions -- a wedge -- no vertex lies, and boundary
//! edges do not cross, so one edge is nearest throughout the wedge and
//! bounds what is seen there. The visibility polygon is the fan of those
//! pieces.
//!
//! Every choice is exact: the order of directions (a half-plane and an
//! orientation), which edges a wedge's middle ray meets, and which of them
//! is nearest along it. The middle ray points along `(w1 - v) + (w2 - v)`,
//! the sum of its two bounding directions, so it is exact too; a wedge is
//! narrower than a half-turn whenever the viewpoint is strictly inside, so
//! that sum points into it. Only the output points where a wedge's
//! boundary ray meets its edge -- the ends of shadows -- are rounded, once,
//! and the ring is presented like every other region operation.

use axiolid_core::{Point2, Tolerance};
use axiolid_exact::{certify, Arith, SignExpr};
use axiolid_guarantees::Sign;

use crate::arc::ArcRing;
use crate::arrangement::ArcArrangement;
use crate::region::Region;
use crate::OverlayError;

/// Why no visibility polygon was built.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum VisibilityError {
    /// The viewpoint is on the region's boundary or outside it.
    NotInside,
    /// The result failed the overlay's own checks.
    Overlay(OverlayError),
}

impl From<OverlayError> for VisibilityError {
    fn from(error: OverlayError) -> Self {
        Self::Overlay(error)
    }
}

impl Region {
    /// The part of the region in sight of `viewpoint`: every point the
    /// straight segment from the viewpoint reaches without leaving the
    /// region. Walls and holes cast shadows.
    ///
    /// # Errors
    ///
    /// [`VisibilityError::NotInside`] for a viewpoint on the boundary or
    /// outside.
    pub fn visibility_polygon(
        &self,
        viewpoint: Point2,
        tolerance: Tolerance,
    ) -> Result<Self, VisibilityError> {
        if !viewpoint.is_finite() || !strictly_inside(self, viewpoint) {
            return Err(VisibilityError::NotInside);
        }
        let edges: Vec<(Point2, Point2)> = self
            .boundary_rings()
            .iter()
            .flat_map(|r| {
                let n = r.points.len();
                (0..n).map(move |i| (r.points[i], r.points[(i + 1) % n]))
            })
            .filter(|(a, b)| a != b)
            .collect();
        let v = viewpoint;
        // One vertex per direction, counter-clockwise from +x.
        let mut around: Vec<Point2> = Vec::new();
        for &(a, _) in &edges {
            let at = around.partition_point(|&w| before(v, w, a));
            if at < around.len() && same_direction(v, around[at], a) {
                continue;
            }
            around.insert(at, a);
        }
        let k = around.len();
        let mut ring: Vec<Point2> = Vec::with_capacity(2 * k);
        for i in 0..k {
            let (w1, w2) = (around[i], around[(i + 1) % k]);
            let Some(edge) = nearest(v, w1, w2, &edges) else {
                // Every wedge meets the boundary when the viewpoint is inside.
                return Err(VisibilityError::NotInside);
            };
            for w in [w1, w2] {
                let p = on_ray(v, w, edge);
                if ring.last() != Some(&p) {
                    ring.push(p);
                }
            }
        }
        while ring.len() > 1 && ring.first() == ring.last() {
            ring.pop();
        }
        if ring.len() < 3 {
            return Ok(Self::empty());
        }
        let arrangement = ArcArrangement::new(&[ArcRing::from_points(&ring)], tolerance)?;
        Ok(crate::minkowski::region_of(
            &arrangement,
            |f| f[0],
            true,
            tolerance,
        )?)
    }
}

/// Where the ray from `v` towards `w` meets `edge`: an endpoint exactly
/// when it lies on the ray, else rounded once.
fn on_ray(v: Point2, w: Point2, (a, b): (Point2, Point2)) -> Point2 {
    for end in [a, b] {
        if sign(&Orient { a: v, b: w, c: end }) == Sign::Zero {
            return end;
        }
    }
    let (d, e) = (w - v, b - a);
    let t = (a - v).perp_dot(e) / d.perp_dot(e);
    v + d * t
}

/// The nearest edge along the wedge's middle ray from `v`, between the
/// directions to `w1` and `w2`.
fn nearest(
    v: Point2,
    w1: Point2,
    w2: Point2,
    edges: &[(Point2, Point2)],
) -> Option<(Point2, Point2)> {
    let mut best: Option<(Point2, Point2)> = None;
    for &edge in edges {
        if sign(&Hits { v, w1, w2, edge }) != Sign::Positive {
            continue;
        }
        if best.is_none_or(|b| {
            sign(&Nearer {
                v,
                w1,
                w2,
                near: edge,
                far: b,
            }) == Sign::Positive
        }) {
            best = Some(edge);
        }
    }
    best
}

/// Whether the direction from `v` to `a` comes strictly before that to
/// `b`, counter-clockwise from `+x`.
fn before(v: Point2, a: Point2, b: Point2) -> bool {
    let half = |p: Point2| u8::from(!(p.y > v.y || (p.y == v.y && p.x > v.x)));
    let (ha, hb) = (half(a), half(b));
    if ha != hb {
        return ha < hb;
    }
    sign(&Orient { a: v, b: a, c: b }) == Sign::Positive
}

fn same_direction(v: Point2, a: Point2, b: Point2) -> bool {
    !before(v, a, b) && !before(v, b, a)
}

/// Strictly inside the region: off every boundary ring and with an odd
/// number of rings round it, exactly.
fn strictly_inside(region: &Region, p: Point2) -> bool {
    let mut winding = 0usize;
    for ring in region.boundary_rings() {
        let n = ring.points.len();
        let mut crossings = 0usize;
        for i in 0..n {
            let (a, b) = (ring.points[i], ring.points[(i + 1) % n]);
            let s = sign(&Orient { a, b, c: p });
            let between = p.x >= a.x.min(b.x)
                && p.x <= a.x.max(b.x)
                && p.y >= a.y.min(b.y)
                && p.y <= a.y.max(b.y);
            if s == Sign::Zero && between {
                return false;
            }
            // Crossings of the rightward ray, half-open in y.
            let upward = a.y <= p.y && b.y > p.y;
            let downward = b.y <= p.y && a.y > p.y;
            if (upward && s == Sign::Positive) || (downward && s == Sign::Negative) {
                crossings += 1;
            }
        }
        winding += crossings % 2;
    }
    winding % 2 == 1
}

/// The inputs are finite, so the exact tier always decides.
fn sign<E: SignExpr>(e: &E) -> Sign {
    certify(e).unwrap_or(Sign::Zero)
}

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

/// `(p - q)` in `T`.
fn diff<T: Arith>(p: Point2, q: Point2) -> [T; 2] {
    [
        T::from_f64(p.x).sub(&T::from_f64(q.x)),
        T::from_f64(p.y).sub(&T::from_f64(q.y)),
    ]
}

fn cross<T: Arith>(a: &[T; 2], b: &[T; 2]) -> T {
    a[0].mul(&b[1]).sub(&a[1].mul(&b[0]))
}

/// The middle direction of a wedge, `(w1 - v) + (w2 - v)`.
fn middle<T: Arith>(v: Point2, w1: Point2, w2: Point2) -> [T; 2] {
    let (p, q) = (diff::<T>(w1, v), diff::<T>(w2, v));
    [p[0].add(&q[0]), p[1].add(&q[1])]
}

/// The ray from `v` along `u` meets the line of `(a, b)` at parameter
/// `N / D`: `N = (a - v) x (b - a)`, `D = u x (b - a)`.
fn parameter<T: Arith>(v: Point2, u: &[T; 2], (a, b): (Point2, Point2)) -> (T, T) {
    let e = diff::<T>(b, a);
    (cross(&diff::<T>(a, v), &e), cross(u, &e))
}

/// Positive when the wedge's middle ray crosses the edge's interior ahead
/// of the viewpoint.
struct Hits {
    v: Point2,
    w1: Point2,
    w2: Point2,
    edge: (Point2, Point2),
}

impl SignExpr for Hits {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let u = middle::<T>(self.v, self.w1, self.w2);
        let (a, b) = self.edge;
        let sa = cross(&u, &diff::<T>(a, self.v)).sign()?;
        let sb = cross(&u, &diff::<T>(b, self.v)).sign()?;
        let straddles = matches!(
            (sa, sb),
            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
        );
        if !straddles {
            return Some(Sign::Negative);
        }
        let (n, d) = parameter::<T>(self.v, &u, self.edge);
        let (sn, sd) = (n.sign()?, d.sign()?);
        Some(if sn != Sign::Zero && sn == sd {
            Sign::Positive
        } else {
            Sign::Negative
        })
    }
}

/// Positive when `near` is strictly nearer than `far` along the wedge's
/// middle ray; both are known to cross it ahead.
struct Nearer {
    v: Point2,
    w1: Point2,
    w2: Point2,
    near: (Point2, Point2),
    far: (Point2, Point2),
}

impl SignExpr for Nearer {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let u = middle::<T>(self.v, self.w1, self.w2);
        let (n1, d1) = parameter::<T>(self.v, &u, self.near);
        let (n2, d2) = parameter::<T>(self.v, &u, self.far);
        // t1 < t2  <=>  (n2 d1 - n1 d2) has the sign of d1 d2.
        let gap = n2.mul(&d1).sub(&n1.mul(&d2)).sign()?;
        let same = d1.sign()? == d2.sign()?;
        Some(match (gap, same) {
            (Sign::Zero, _) => Sign::Zero,
            (s, true) => s,
            (Sign::Positive, false) => Sign::Negative,
            (_, false) => Sign::Positive,
        })
    }
}
