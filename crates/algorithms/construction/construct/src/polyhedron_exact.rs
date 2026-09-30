//! Exact arithmetic for the planar-faced boolean in [`crate::polyhedron`]
//! (#199).
//!
//! Two jobs need more than f64 there, and both use dyadic numbers
//! (`axiolid-exact`), so every decision is a sign of an exact value:
//!
//! - **Split points.** A crossing is stored as f64 (ADR 0045), but as the
//!   double NEAREST to the exact point, so one exact point always gets the
//!   same bits and a crossing that is a double comes back exactly.
//! - **Classification points.** A fragment is classified by one point of its
//!   relative interior. A fragment one ULP wide has no double strictly
//!   inside it, so its f64 centroid rounds onto its own boundary -- often
//!   onto the other solid's boundary too -- and the answer is then about the
//!   wrong point. Such fragments are classified at an exact dyadic interior
//!   point instead.

use axiolid_core::Point3;
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;

/// A point with exact dyadic coordinates.
pub(crate) type ExactPoint = [Dyadic; 3];

/// Where a coplanar point sits relative to a face ring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RingPosition {
    Inside,
    OnEdge,
    Outside,
}

fn dy(v: f64) -> Option<Dyadic> {
    Dyadic::try_from_f64(v)
}

/// The exact value of a finite point.
pub(crate) fn exact(p: Point3) -> Option<ExactPoint> {
    Some([dy(p.x)?, dy(p.y)?, dy(p.z)?])
}

fn minus(p: &ExactPoint, q: &ExactPoint) -> ExactPoint {
    [p[0].sub(&q[0]), p[1].sub(&q[1]), p[2].sub(&q[2])]
}

fn cross(u: &ExactPoint, v: &ExactPoint) -> ExactPoint {
    [
        u[1].mul(&v[2]).sub(&u[2].mul(&v[1])),
        u[2].mul(&v[0]).sub(&u[0].mul(&v[2])),
        u[0].mul(&v[1]).sub(&u[1].mul(&v[0])),
    ]
}

fn dot(u: &ExactPoint, v: &ExactPoint) -> Dyadic {
    u[0].mul(&v[0]).add(&u[1].mul(&v[1])).add(&u[2].mul(&v[2]))
}

fn sign(value: &Dyadic) -> Sign {
    value
        .sign()
        .expect("a dyadic value always has a decided sign")
}

/// `a < b`, exactly.
fn less(a: &Dyadic, b: &Dyadic) -> bool {
    sign(&a.sub(b)) == Sign::Negative
}

/// The exact normal `(f1 - f0) x (f2 - f0)` of a face.
fn normal(face: &[Point3]) -> Option<ExactPoint> {
    let a = exact(face[0])?;
    Some(cross(
        &minus(&exact(face[1])?, &a),
        &minus(&exact(face[2])?, &a),
    ))
}

/// Orientation of `d` against the plane `a`, `b`, `c`, exactly.
///
/// Its sign convention is its own; the exact probe path compares signs
/// only with other signs from this function, never with `orient3d`.
fn orient3(a: &ExactPoint, b: &ExactPoint, c: &ExactPoint, d: &ExactPoint) -> Sign {
    sign(&dot(&minus(a, d), &cross(&minus(b, d), &minus(c, d))))
}

/// Which side of a face's plane an exact point lies on.
pub(crate) fn side(face: &[Point3], point: &ExactPoint) -> Option<Sign> {
    Some(orient3(
        &exact(face[0])?,
        &exact(face[1])?,
        &exact(face[2])?,
        point,
    ))
}

/// The turn of edge `a`-`b` seen along the segment `origin`-`far`.
pub(crate) fn turn(origin: &ExactPoint, far: Point3, a: Point3, b: Point3) -> Option<Sign> {
    Some(orient3(origin, &exact(far)?, &exact(a)?, &exact(b)?))
}

/// Where a coplanar exact point sits relative to a face ring.
///
/// The same parity rule as the f64 test in `polyhedron`, over exact values.
/// The projection drops an axis along which the EXACT normal is non-zero,
/// so the projected ring is never degenerate however thin the face is.
pub(crate) fn ring_position(face: &[Point3], point: &ExactPoint) -> Option<RingPosition> {
    let n = normal(face)?;
    // Prefer the dominant axis (the best-conditioned projection); any axis
    // with a non-zero exact component is correct.
    let approx = [0, 1, 2].map(|k| n[k].to_f64().abs());
    let mut axes = [0usize, 1, 2];
    axes.sort_by(|&i, &j| approx[j].total_cmp(&approx[i]));
    let drop = *axes.iter().find(|&&k| sign(&n[k]) != Sign::Zero)?;
    let (u, v) = ((drop + 1) % 3, (drop + 2) % 3);

    let ring: Vec<ExactPoint> = face.iter().map(|&p| exact(p)).collect::<Option<_>>()?;
    let (qu, qv) = (&point[u], &point[v]);
    let orient2 = |a: &ExactPoint, b: &ExactPoint| {
        sign(
            &b[u]
                .sub(&a[u])
                .mul(&qv.sub(&a[v]))
                .sub(&b[v].sub(&a[v]).mul(&qu.sub(&a[u]))),
        )
    };
    let within = |q: &Dyadic, a: &Dyadic, b: &Dyadic| {
        let (lo, hi) = if less(b, a) { (b, a) } else { (a, b) };
        !less(q, lo) && !less(hi, q)
    };

    for i in 0..ring.len() {
        let (a, b) = (&ring[i], &ring[(i + 1) % ring.len()]);
        if orient2(a, b) == Sign::Zero && within(qu, &a[u], &b[u]) && within(qv, &a[v], &b[v]) {
            return Some(RingPosition::OnEdge);
        }
    }
    let mut inside = false;
    for i in 0..ring.len() {
        let (a, b) = (&ring[i], &ring[(i + 1) % ring.len()]);
        if less(qv, &a[v]) != less(qv, &b[v]) {
            let s = orient2(a, b);
            let upward = less(&a[v], &b[v]);
            if (upward && s == Sign::Negative) || (!upward && s == Sign::Positive) {
                inside = !inside;
            }
        }
    }
    Some(if inside {
        RingPosition::Inside
    } else {
        RingPosition::Outside
    })
}

/// Whether two faces' outward normals point the same way, exactly.
pub(crate) fn normals_agree(first: &[Point3], second: &[Point3]) -> Option<bool> {
    Some(sign(&dot(&normal(first)?, &normal(second)?)) == Sign::Positive)
}

/// A point strictly inside a planar ring, exactly, or `None` when the ring
/// encloses no area.
///
/// Each fan triangle `(v0, vi, vi+1)` offers the dyadic point
/// `v0/4 + vi/4 + vi+1/2`, strictly inside the triangle when it has area.
/// The first such point certified strictly inside the ring is used, which
/// also serves a non-convex ring, whose first fan triangles may fall outside
/// it. The point lies in the ring's plane because it is an affine
/// combination of the ring's own coplanar vertices.
pub(crate) fn interior_point(ring: &[Point3]) -> Option<ExactPoint> {
    let quarter = dy(0.25)?;
    let half = dy(0.5)?;
    let first = exact(ring[0])?;
    for i in 1..ring.len() - 1 {
        let (b, c) = (exact(ring[i])?, exact(ring[i + 1])?);
        let area = cross(&minus(&b, &first), &minus(&c, &first));
        if area.iter().all(|a| sign(a) == Sign::Zero) {
            continue;
        }
        let point: ExactPoint =
            [0, 1, 2].map(|k| first[k].add(&b[k]).mul(&quarter).add(&c[k].mul(&half)));
        if ring_position(ring, &point)? == RingPosition::Inside {
            return Some(point);
        }
    }
    None
}

/// A nearby double for an exact point, for choosing a ray length only.
pub(crate) fn approx(point: &ExactPoint) -> Point3 {
    Point3::new(point[0].to_f64(), point[1].to_f64(), point[2].to_f64())
}

/// Where segment `a`-`b` meets the plane through `plane`'s first 3 points,
/// each coordinate the double nearest to the exact crossing.
///
/// The naive `a + (b - a) * t` rounds `t` and then rounds again, so the cut
/// through the plane `z = 1/3` could come back as `0.33333333333333326`.
/// Two splits that meet the same exact point then land ULPs apart, and a
/// later split through the pair emits a ring that encloses no area (#199).
/// Correct rounding is a function of the exact point alone, so every
/// construction of one exact point yields the same bits, and a crossing that
/// IS a double (every cut of an axis-aligned edge by an axis-aligned plane)
/// comes back exactly. This is not snapping: nothing moves towards a
/// neighbour; each coordinate is the best f64 for its own exact value.
///
/// `guess` is the f64 formula's answer; it only seeds the search.
pub(crate) fn plane_crossing(
    plane: &[Point3],
    a: Point3,
    b: Point3,
    guess: Point3,
) -> Option<Point3> {
    let normal = normal(plane)?;
    let (ea, eb) = (exact(a)?, exact(b)?);
    let along = minus(&eb, &ea);
    // p = a + (b - a) * reach / span, exactly.
    let mut span = dot(&normal, &along);
    let mut reach = dot(&normal, &minus(&exact(plane[0])?, &ea));
    match sign(&span) {
        Sign::Zero => return None,
        Sign::Negative => {
            span = span.neg();
            reach = reach.neg();
        }
        _ => {}
    }
    let mut out = [0.0; 3];
    for k in 0..3 {
        let (from, to) = ([a.x, a.y, a.z][k], [b.x, b.y, b.z][k]);
        out[k] = if sign(&along[k]) == Sign::Zero {
            // The segment does not move along this axis: the exact value is
            // the endpoint's own coordinate.
            from
        } else {
            let numerator = ea[k].mul(&span).add(&along[k].mul(&reach));
            // The crossing lies strictly between the endpoints, so the
            // answer is bracketed by them.
            let (lo, hi) = if from < to { (from, to) } else { (to, from) };
            nearest_ratio(&numerator, &span, [guess.x, guess.y, guess.z][k], lo, hi)?
        };
    }
    Some(Point3::new(out[0], out[1], out[2]))
}

/// A total-order key for finite doubles: `key(x) < key(y)` iff `x < y`.
fn key(x: f64) -> u64 {
    let bits = (x + 0.0).to_bits();
    if bits >> 63 == 1 {
        !bits
    } else {
        bits | (1 << 63)
    }
}

fn unkey(k: u64) -> f64 {
    f64::from_bits(if k >> 63 == 1 { k & !(1 << 63) } else { !k })
}

/// `numerator / denominator` (denominator positive) rounded to the nearest
/// double, ties to even, known to lie in `[lo, hi]`.
///
/// Each step compares the exact quotient with the midpoint between two
/// adjacent doubles. That midpoint is a dyadic, so the comparison is exact.
/// Starting from the f64 guess, a few steps settle the answer; a guess
/// further off than that is abandoned for bisection over the bracket, which
/// needs at most 64 comparisons, so a bad seed costs time, never
/// correctness.
fn nearest_ratio(
    numerator: &Dyadic,
    denominator: &Dyadic,
    guess: f64,
    lo: f64,
    hi: f64,
) -> Option<f64> {
    let half = dy(0.5)?;
    // Sign of `numerator / denominator - (x + y) / 2`.
    let beyond = |x: f64, y: f64| -> Option<Sign> {
        let mid = dy(x)?.add(&dy(y)?).mul(&half);
        Some(sign(&numerator.sub(&mid.mul(denominator))))
    };
    // Sign of `numerator / denominator - x`.
    let versus = |x: f64| -> Option<Sign> { Some(sign(&numerator.sub(&dy(x)?.mul(denominator)))) };
    let even = |r: f64| r.to_bits() & 1 == 0;

    const WALK: usize = 8;
    let mut r = if guess.is_finite() {
        guess.clamp(lo, hi)
    } else {
        lo
    };
    let mut settled = false;
    for _ in 0..WALK {
        let up = r.next_up();
        match beyond(r, up)? {
            Sign::Positive => r = up,
            Sign::Zero if !even(r) => r = up,
            _ => {
                settled = true;
                break;
            }
        }
    }
    if settled {
        for _ in 0..WALK {
            let down = r.next_down();
            match beyond(down, r)? {
                Sign::Negative => r = down,
                Sign::Zero if !even(r) => r = down,
                _ => return Some(r),
            }
        }
    }

    // Bisection: the largest double in the bracket not above the quotient,
    // then the nearer of it and its successor.
    let (mut low, mut high) = (key(lo), key(hi));
    if versus(lo)? == Sign::Negative || versus(hi)? == Sign::Positive {
        return None;
    }
    while high - low > 1 {
        let mid = low + (high - low) / 2;
        if versus(unkey(mid))? == Sign::Negative {
            high = mid;
        } else {
            low = mid;
        }
    }
    let floor = if versus(unkey(high))? == Sign::Negative {
        unkey(low)
    } else {
        unkey(high)
    };
    let next = floor.next_up();
    Some(match beyond(floor, next)? {
        Sign::Positive => next,
        Sign::Zero if !even(floor) => next,
        _ => floor,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quotient(n: f64, d: f64) -> (Dyadic, Dyadic) {
        (dy(n).unwrap(), dy(d).unwrap())
    }

    #[test]
    fn nearest_ratio_is_the_ieee_quotient_from_any_seed() {
        // IEEE division is correctly rounded, so it is the oracle; every
        // seed in the bracket -- including ones far beyond the walk -- must
        // reach it.
        for &(n, d) in &[(1.0, 3.0), (2.0, 3.0), (1.0, 7.0), (5.0, 9.0), (0.1, 0.3)] {
            let (num, den) = quotient(n, d);
            let want = n / d;
            for seed in [want, want.next_up(), 0.0, 1.0, f64::NAN, 0.5] {
                assert_eq!(
                    nearest_ratio(&num, &den, seed, 0.0, 1.0),
                    Some(want),
                    "{n}/{d} from seed {seed}"
                );
            }
        }
    }

    #[test]
    fn a_crossing_that_is_a_double_comes_back_exactly() {
        let t = 1.0 / 3.0;
        let plane = [
            Point3::new(0.0, 0.0, t),
            Point3::new(1.0, 0.0, t),
            Point3::new(0.0, 1.0, t),
        ];
        let (a, b) = (Point3::new(0.5, 0.5, 0.0), Point3::new(0.5, 0.5, 1.0));
        // A deliberately poor seed: the answer must not depend on it.
        let cut = plane_crossing(&plane, a, b, Point3::new(0.5, 0.5, 0.9)).unwrap();
        assert_eq!(cut, Point3::new(0.5, 0.5, t));
    }

    #[test]
    fn a_ring_enclosing_no_area_has_no_interior_point() {
        let z = 1.0 / 3.0;
        // The collapsed quad from #199: vertex pairs one ULP apart.
        let collapsed = [
            Point3::new(1.0, 1.0, z),
            Point3::new(2.0 / 3.0, 1.0, z),
            Point3::new(1.0, 1.0, z.next_down()),
            Point3::new(2.0 / 3.0, 1.0, z),
        ];
        assert_eq!(interior_point(&collapsed), None);
        // One ULP wide is thin, not collapsed: it has an exact interior.
        let sliver = [
            Point3::new(0.0, 1.0, z.next_down()),
            Point3::new(0.0, 1.0, z),
            Point3::new(1.0, 1.0, z),
            Point3::new(1.0, 1.0, z.next_down()),
        ];
        let inside = interior_point(&sliver).expect("a sliver has area");
        assert_eq!(ring_position(&sliver, &inside), Some(RingPosition::Inside));
    }

    #[test]
    fn ties_round_to_even_from_the_walk_and_from_bisection() {
        let one = dy(1.0).unwrap();
        let ulp = dy(f64::EPSILON).unwrap();
        let half_ulp = ulp.mul(&dy(0.5).unwrap());
        // 1 + ulp/2 lies midway between 1 (even) and 1 + ulp (odd).
        let low_tie = one.add(&half_ulp);
        // 1 + 3 ulp/2 lies midway between 1 + ulp (odd) and 1 + 2 ulp (even).
        let high_tie = one.add(&ulp).add(&half_ulp);
        let (even_low, even_high) = (1.0, 1.0 + 2.0 * f64::EPSILON);
        let seeds = |want: f64| [want, want.next_up(), want.next_down(), 1.0 + f64::EPSILON];
        for (tie, want) in [(&low_tie, even_low), (&high_tie, even_high)] {
            for seed in seeds(want) {
                assert_eq!(
                    nearest_ratio(tie, &one, seed, 0.0, 2.0),
                    Some(want),
                    "seed {seed}"
                );
            }
            // No seed: the walk from the bracket's end gives up and bisects.
            assert_eq!(nearest_ratio(tie, &one, f64::NAN, 0.0, 2.0), Some(want));
        }
    }

    #[test]
    fn an_interior_point_is_never_taken_on_the_ring() {
        // A non-convex ring whose first fan triangle offers (2, 1), which
        // lies on the edge (4, 1)-(1, 1); the next two offer points outside
        // the ring. Only the fourth, (0.75, 3.25), is strictly inside.
        let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
        let ring = [
            p(0.0, 4.0),
            p(0.0, 0.0),
            p(4.0, 0.0),
            p(4.0, 1.0),
            p(1.0, 1.0),
            p(1.0, 4.0),
        ];
        let point = interior_point(&ring).expect("the ring has area");
        assert_eq!(approx(&point), p(0.75, 3.25));
    }
}
