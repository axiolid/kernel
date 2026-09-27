//! The minimum enclosing circle of a point set (#118).
//!
//! # Exact choice, enclosed output
//!
//! Welzl's algorithm, in its iterative form: points are visited in a fixed
//! pseudo-random order, and a point outside the current circle is put on
//! the boundary of the next one. The one geometric decision is whether a
//! point lies in the circle spanned by one, two or three support points,
//! and it is exact:
//!
//! - one support point: equality;
//! - two, the circle on them as diameter: the sign of `(p - a) . (p - b)`;
//! - three, their circumcircle: the incircle determinant against the
//!   triangle's orientation.
//!
//! Each is a polynomial in differences of the input `f64`s, decided in
//! intervals and else in dyadics, so the support set is the exact minimum
//! circle's. The visiting order is fixed, so the result does not depend on
//! chance; the circle does not depend on the input order either (the
//! minimum circle is unique), though which support points are reported can
//! when four or more lie on it.
//!
//! Only the output is rounded. The exact centre is enclosed from the
//! support points (a midpoint, or a circumcentre as a quotient of exact
//! dyadic polynomials), and the returned radius is rounded up so the
//! returned circle contains the exact one. [`CircleEvidence::error`] bounds
//! both the centre's distance from the exact centre and the radius's excess
//! over the exact radius.

use axiolid_core::Point2;
use axiolid_exact::{certify, Arith, Dyadic, Interval, SignExpr};
use axiolid_guarantees::Sign;

/// A circle by its centre and radius.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnclosingCircle {
    /// Centre.
    pub centre: Point2,
    /// Radius; zero for a single distinct point.
    pub radius: f64,
}

/// Which points determine the circle and how exact the output is.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct CircleEvidence {
    /// Indices into the input of the one, two or three points the exact
    /// minimum circle passes through and is determined by, ascending.
    pub support: Vec<usize>,
    /// A bound on the distance between the returned centre and the exact
    /// one, and on the returned radius's excess over the exact radius. The
    /// returned radius is never below the exact radius, so the returned
    /// circle contains every input point. Zero for a single point.
    pub error: f64,
}

/// The circle and its evidence.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct MinimumCircle {
    /// The circle.
    pub circle: EnclosingCircle,
    /// Its support and error bound.
    pub evidence: CircleEvidence,
}

/// Why no circle was built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CircleError {
    /// No points.
    Empty,
    /// A coordinate was not finite.
    NonFinite,
}

/// `(p - a) . (p - b)`: not positive exactly when `p` lies in the circle
/// on `a b` as diameter.
struct Diametral {
    a: Point2,
    b: Point2,
    p: Point2,
}

impl SignExpr for Diametral {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let (px, py) = (f(self.p.x), f(self.p.y));
        let ax = px.sub(&f(self.a.x));
        let ay = py.sub(&f(self.a.y));
        let bx = px.sub(&f(self.b.x));
        let by = py.sub(&f(self.b.y));
        ax.mul(&bx).add(&ay.mul(&by)).sign()
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

/// The incircle determinant: positive when `p` lies inside the circle
/// through `a b c` taken counter-clockwise.
struct InCircle {
    a: Point2,
    b: Point2,
    c: Point2,
    p: Point2,
}

impl SignExpr for InCircle {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let row = |q: Point2| {
            let x = f(q.x).sub(&f(self.p.x));
            let y = f(q.y).sub(&f(self.p.y));
            let l = x.mul(&x).add(&y.mul(&y));
            (x, y, l)
        };
        let (ax, ay, al) = row(self.a);
        let (bx, by, bl) = row(self.b);
        let (cx, cy, cl) = row(self.c);
        let minor = |x1: &T, y1: &T, x2: &T, y2: &T| x1.mul(y2).sub(&y1.mul(x2));
        al.mul(&minor(&bx, &by, &cx, &cy))
            .sub(&bl.mul(&minor(&ax, &ay, &cx, &cy)))
            .add(&cl.mul(&minor(&ax, &ay, &bx, &by)))
            .sign()
    }
}

/// The inputs are finite, so the exact tier always decides.
fn sign<E: SignExpr>(e: &E) -> Sign {
    certify(e).unwrap_or(Sign::Zero)
}

/// Whether `p` lies in (or on) the circle determined by `support`.
fn inside(points: &[Point2], support: &[usize], p: Point2) -> bool {
    match *support {
        [a] => points[a] == p,
        [a, b] => {
            sign(&Diametral {
                a: points[a],
                b: points[b],
                p,
            }) != Sign::Positive
        }
        [a, b, c] => {
            let (a, b, c) = (points[a], points[b], points[c]);
            let turn = sign(&Orient { a, b, c });
            let side = sign(&InCircle { a, b, c, p });
            side == Sign::Zero || side == turn
        }
        _ => unreachable!("a circle has one to three support points"),
    }
}

/// A fixed pseudo-random permutation of `0..n` (Fisher-Yates driven by
/// splitmix64 from a constant seed): the expected linear running time of
/// a random order, with a result that never changes between runs.
fn visiting_order(n: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut next = || {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    for i in (1..n).rev() {
        let j = (next() % (i as u64 + 1)) as usize;
        order.swap(i, j);
    }
    order
}

/// The support of the minimum circle, by Welzl's iterative algorithm.
fn welzl(points: &[Point2]) -> Vec<usize> {
    let order = visiting_order(points.len());
    let mut support = vec![order[0]];
    for i in 1..order.len() {
        let pi = order[i];
        if inside(points, &support, points[pi]) {
            continue;
        }
        // `pi` lies on the minimum circle of the points visited so far.
        support = vec![pi];
        for j in 0..i {
            let pj = order[j];
            if inside(points, &support, points[pj]) {
                continue;
            }
            // So do `pi` and `pj`.
            support = vec![pi, pj];
            for &pk in &order[..j] {
                if !inside(points, &support, points[pk]) {
                    support = vec![pi, pj, pk];
                }
            }
        }
    }
    support
}

fn exact(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

/// Enclosures of the exact centre of the circle on `support`.
fn centre_enclosure(points: &[Point2], support: &[usize]) -> [Interval; 2] {
    match *support {
        [a] => [Interval::point(points[a].x), Interval::point(points[a].y)],
        [a, b] => {
            let half = exact(0.5);
            let mid = |s: f64, t: f64| exact(s).add(&exact(t)).mul(&half).enclosure();
            [mid(points[a].x, points[b].x), mid(points[a].y, points[b].y)]
        }
        [a, b, c] => {
            let (a, b, c) = (points[a], points[b], points[c]);
            let (ux, uy) = (exact(b.x).sub(&exact(a.x)), exact(b.y).sub(&exact(a.y)));
            let (vx, vy) = (exact(c.x).sub(&exact(a.x)), exact(c.y).sub(&exact(a.y)));
            let uu = ux.mul(&ux).add(&uy.mul(&uy));
            let vv = vx.mul(&vx).add(&vy.mul(&vy));
            let d = ux.mul(&vy).sub(&uy.mul(&vx)).mul(&exact(2.0)).enclosure();
            let nx = uu.mul(&vy).sub(&vv.mul(&uy)).enclosure();
            let ny = vv.mul(&ux).sub(&uu.mul(&vx)).enclosure();
            [
                Interval::point(a.x).add(&nx.quotient(d)),
                Interval::point(a.y).add(&ny.quotient(d)),
            ]
        }
        _ => unreachable!("a circle has one to three support points"),
    }
}

/// An upper bound on the distance from `p` to any point of the box
/// `centre`, and a lower bound on the distance to the nearest.
fn distance_bounds(p: Point2, centre: &[Interval; 2]) -> (f64, f64) {
    let dx = Interval::point(p.x).sub(&centre[0]);
    let dy = Interval::point(p.y).sub(&centre[1]);
    let squared = dx.mul(&dx).add(&dy.mul(&dy));
    // `sqrt` is correctly rounded, so one step outward covers it.
    let low = squared.lo().max(0.0).sqrt().next_down().max(0.0);
    (low, squared.hi().sqrt().next_up())
}

/// The minimum circle enclosing `points`.
///
/// # Errors
///
/// [`CircleError::Empty`] for no points, [`CircleError::NonFinite`] for a
/// coordinate that is not finite.
pub fn minimum_enclosing_circle(points: &[Point2]) -> Result<MinimumCircle, CircleError> {
    if points.is_empty() {
        return Err(CircleError::Empty);
    }
    if !points.iter().all(|p| p.is_finite()) {
        return Err(CircleError::NonFinite);
    }
    let mut support = welzl(points);
    support.sort_unstable();
    if let [a] = *support {
        return Ok(MinimumCircle {
            circle: EnclosingCircle {
                centre: points[a],
                radius: 0.0,
            },
            evidence: CircleEvidence {
                support,
                error: 0.0,
            },
        });
    }
    let enclosure = centre_enclosure(points, &support);
    let mid = |i: Interval| i.lo() + 0.5 * (i.hi() - i.lo());
    let centre = Point2::new(mid(enclosure[0]), mid(enclosure[1]));
    // How far the returned centre can be from the exact one: the sum of
    // the per-axis gaps is at least their Euclidean length.
    let gap = |i: Interval, m: f64| (i.hi() - m).max(m - i.lo());
    let centre_error = (gap(enclosure[0], centre.x) + gap(enclosure[1], centre.y)).next_up();
    let (radius_low, radius_high) = distance_bounds(points[support[0]], &enclosure);
    let mut radius = (radius_high + centre_error).next_up();
    // The circle of that radius about the returned centre contains the
    // exact circle; confirm every point with enclosures all the same.
    let at = [Interval::point(centre.x), Interval::point(centre.y)];
    for &p in points {
        radius = radius.max(distance_bounds(p, &at).1);
    }
    Ok(MinimumCircle {
        circle: EnclosingCircle { centre, radius },
        evidence: CircleEvidence {
            support,
            error: (radius - radius_low).next_up(),
        },
    })
}
