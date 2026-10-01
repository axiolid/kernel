//! `incircle` and `insphere`: is a point inside a circumscribed ball?
//!
//! These are the Delaunay predicates. `incircle(a, b, c, d)` asks whether `d`
//! lies inside the circle through `a`, `b`, `c`; `insphere` is the 3D analogue.
//! A zero means `d` lies exactly on the ball -- the cocircular/cospherical case
//! that makes a Delaunay triangulation non-unique and, if misjudged, produces
//! inverted or overlapping cells.
//!
//! Both are lifted determinants: adding a coordinate equal to the squared
//! distance from the origin turns "inside a ball" into "below a hyperplane".
//! That lift squares the operand magnitudes, so the error bound grows faster
//! than `orient*`'s and the filter fails sooner -- which is why the escalation
//! rate is measured rather than assumed.

use axiolid_core::{Point2, Point3};
use axiolid_guarantees::{Certified, Precision, Sign};

use crate::arithmetic::{expansion_product, expansion_sign, expansion_sum, negate_expansion};
use crate::expansion::two_diff;
use crate::orient3::{highest_bit_exponent, least_significant_bit_exponent};

/// Machine epsilon for binary64.
const EPSILON: f64 = f64::EPSILON / 2.0;

/// Relative error bound for the lifted 3x3 `incircle` determinant.
const INCIRCLE_ERROR_FACTOR: f64 = (10.0 + 96.0 * EPSILON) * EPSILON;

/// Relative error bound for the lifted 4x4 `insphere` determinant.
const INSPHERE_ERROR_FACTOR: f64 = (16.0 + 224.0 * EPSILON) * EPSILON;

/// Is `d` inside the circle through `a`, `b`, `c`?
///
/// [`Sign::Positive`] means inside when `a, b, c` are counter-clockwise.
/// Callers that cannot guarantee that orientation must normalise it first with
/// `orient2d`, because the sign of this determinant flips with it.
///
/// Always [`Certified::Certain`].
#[must_use]
pub fn incircle(a: Point2, b: Point2, c: Point2, d: Point2) -> Certified {
    match incircle_filter(a, b, c, d) {
        Certified::Certain { sign, .. } => Certified::exact_sign(sign),
        _ => Certified::exact_sign(incircle_exact(a, b, c, d)),
    }
}

/// The fast filter alone, exposed so escalation can be measured.
#[must_use]
pub fn incircle_filter(a: Point2, b: Point2, c: Point2, d: Point2) -> Certified {
    let (adx, ady) = (a.x - d.x, a.y - d.y);
    let (bdx, bdy) = (b.x - d.x, b.y - d.y);
    let (cdx, cdy) = (c.x - d.x, c.y - d.y);

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;

    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;

    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;

    let determinant =
        alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);

    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * alift
        + (cdxady.abs() + adxcdy.abs()) * blift
        + (adxbdy.abs() + bdxady.abs()) * clift;

    Certified::from_filter(
        determinant,
        INCIRCLE_ERROR_FACTOR * permanent,
        Precision::F64,
    )
}

/// Exact sign of the lifted `incircle` determinant.
///
/// Every coordinate difference is kept exactly, as a two-term expansion
/// (`two_diff`), and every product and sum after it is an expansion. The
/// differences used to be rounded first, which made this "exact" fallback
/// wrong exactly where the filter hands over -- nearly cocircular points
/// whose differences do not fit an `f64` -- and Delaunay flips driven by it
/// could cycle for ever (#190).
#[must_use]
fn incircle_exact(a: Point2, b: Point2, c: Point2, d: Point2) -> Sign {
    let (adx, ady) = (diff(a.x, d.x), diff(a.y, d.y));
    let (bdx, bdy) = (diff(b.x, d.x), diff(b.y, d.y));
    let (cdx, cdy) = (diff(c.x, d.x), diff(c.y, d.y));

    let bc = minor(&bdx, &cdy, &cdx, &bdy);
    let ca = minor(&cdx, &ady, &adx, &cdy);
    let ab = minor(&adx, &bdy, &bdx, &ady);

    let total = expansion_sum(
        &expansion_sum(&lift(&bc, &adx, &ady), &lift(&ca, &bdx, &bdy)),
        &lift(&ab, &cdx, &cdy),
    );
    expansion_sign(&total)
}

/// `p - q`, exactly, as an expansion.
#[must_use]
fn diff(p: f64, q: f64) -> Vec<f64> {
    let (d, err) = two_diff(p, q);
    let mut e = Vec::with_capacity(2);
    if err != 0.0 {
        e.push(err);
    }
    if d != 0.0 || e.is_empty() {
        e.push(d);
    }
    e
}

/// `p q - r s` over expansions, exactly.
#[must_use]
fn minor(p: &[f64], q: &[f64], r: &[f64], s: &[f64]) -> Vec<f64> {
    expansion_sum(
        &expansion_product(p, q),
        &negate_expansion(&expansion_product(r, s)),
    )
}

/// Multiply an expansion by `x*x + y*y`, exactly.
#[must_use]
fn lift(e: &[f64], x: &[f64], y: &[f64]) -> Vec<f64> {
    let square = expansion_sum(&expansion_product(x, x), &expansion_product(y, y));
    expansion_product(e, &square)
}

/// Is `e` inside the sphere through `a`, `b`, `c`, `d`?
///
/// [`Sign::Positive`] means inside when `a, b, c, d` are positively oriented
/// (`orient3d(a, b, c, d) > 0`). As with [`incircle`], the sign flips with the
/// base orientation, so a caller must normalise it.
///
/// Always [`Certified::Certain`].
#[must_use]
pub fn insphere(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> Certified {
    match insphere_filter(a, b, c, d, e) {
        Certified::Certain { sign, .. } => Certified::exact_sign(sign),
        _ => Certified::exact_sign(insphere_exact(a, b, c, d, e)),
    }
}

/// The fast filter alone, exposed so escalation can be measured.
#[must_use]
pub fn insphere_filter(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> Certified {
    let v = |p: Point3| (p.x - e.x, p.y - e.y, p.z - e.z);
    let (ax, ay, az) = v(a);
    let (bx, by, bz) = v(b);
    let (cx, cy, cz) = v(c);
    let (dx, dy, dz) = v(d);

    let (axby, bxay) = (ax * by, bx * ay);
    let (bxcy, cxby) = (bx * cy, cx * by);
    let (cxdy, dxcy) = (cx * dy, dx * cy);
    let (dxay, axdy) = (dx * ay, ax * dy);
    let (axcy, cxay) = (ax * cy, cx * ay);
    let (bxdy, dxby) = (bx * dy, dx * by);
    let ab = axby - bxay;
    let bc = bxcy - cxby;
    let cd = cxdy - dxcy;
    let da = dxay - axdy;
    let ac = axcy - cxay;
    let bd = bxdy - dxby;

    let abc = az * bc - bz * ac + cz * ab;
    let bcd = bz * cd - cz * bd + dz * bc;
    let cda = cz * da + dz * ac + az * cd;
    let dab = dz * ab + az * bd + bz * da;

    let alift = ax * ax + ay * ay + az * az;
    let blift = bx * bx + by * by + bz * bz;
    let clift = cx * cx + cy * cy + cz * cz;
    let dlift = dx * dx + dy * dy + dz * dz;

    let determinant = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    // Shewchuk's permanent: the same expression over the ABSOLUTE values of
    // every elementary product. Bounding by |abc| and friends instead -- the
    // rounded minors -- is not a bound at all: when a minor cancels, its
    // rounding error is not proportional to its value, and the filter
    // certified a non-zero sign for exactly cospherical points (#126).
    let (az, bz, cz, dz) = (az.abs(), bz.abs(), cz.abs(), dz.abs());
    let ab_plus = axby.abs() + bxay.abs();
    let bc_plus = bxcy.abs() + cxby.abs();
    let cd_plus = cxdy.abs() + dxcy.abs();
    let da_plus = dxay.abs() + axdy.abs();
    let ac_plus = axcy.abs() + cxay.abs();
    let bd_plus = bxdy.abs() + dxby.abs();
    let permanent = (cd_plus * bz + bd_plus * cz + bc_plus * dz) * alift
        + (da_plus * cz + ac_plus * dz + cd_plus * az) * blift
        + (ab_plus * dz + bd_plus * az + da_plus * bz) * clift
        + (bc_plus * az + ac_plus * bz + ab_plus * cz) * dlift;

    Certified::from_filter(
        determinant,
        INSPHERE_ERROR_FACTOR * permanent,
        Precision::F64,
    )
}

/// Exact sign of the lifted 4x4 `insphere` determinant.
///
/// Expands along the lifted column: each 3x3 minor is built from exact 2x2
/// cofactors, scaled by the remaining z difference, then by the squared
/// distance. Nothing is rounded between those steps.
#[must_use]
fn insphere_exact(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> Sign {
    if let Some(sign) = insphere_small_integer(a, b, c, d, e) {
        return sign;
    }
    insphere_expansion(a, b, c, d, e)
}

/// The expansion-arithmetic tier of [`insphere_exact`].
#[must_use]
fn insphere_expansion(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> Sign {
    // Differences exact, as for `incircle_exact` (#190).
    let v = |p: Point3| [diff(p.x, e.x), diff(p.y, e.y), diff(p.z, e.z)];
    let (a3, b3, c3, d3) = (v(a), v(b), v(c), v(d));

    let minor3 = |p: &[Vec<f64>; 3], q: &[Vec<f64>; 3], r: &[Vec<f64>; 3]| {
        let qr = minor(&q[0], &r[1], &r[0], &q[1]);
        let rp = minor(&r[0], &p[1], &p[0], &r[1]);
        let pq = minor(&p[0], &q[1], &q[0], &p[1]);
        expansion_sum(
            &expansion_sum(
                &expansion_product(&qr, &p[2]),
                &expansion_product(&rp, &q[2]),
            ),
            &expansion_product(&pq, &r[2]),
        )
    };

    let bcd = minor3(&b3, &c3, &d3);
    let cda = minor3(&c3, &d3, &a3);
    let dab = minor3(&d3, &a3, &b3);
    let abc = minor3(&a3, &b3, &c3);

    // Cofactor expansion signs alternate: +d -c +b -a.
    let total = expansion_sum(
        &expansion_sum(&lift3(&abc, &d3), &negate_expansion(&lift3(&dab, &c3))),
        &expansion_sum(&lift3(&cda, &b3), &negate_expansion(&lift3(&bcd, &a3))),
    );
    expansion_sign(&total)
}

/// Exact `insphere` in `i128`, for differences on a narrow dyadic grid.
///
/// Exactly degenerate inputs -- a cubic lattice, where every insphere test
/// is an exact zero -- reach the exact path on every call, and the
/// expansion arithmetic allocates throughout. When every coordinate
/// difference is itself exact in `f64` and all are integer multiples of one
/// power of two `2^k` with magnitude below `2^(k + 20)`, the determinant of
/// the scaled integers is at most `72 * 2^100` in magnitude and fits `i128`
/// exactly. Returns `None` (use the expansions) otherwise.
#[must_use]
fn insphere_small_integer(a: Point3, b: Point3, c: Point3, d: Point3, e: Point3) -> Option<Sign> {
    let mut differences = [[0.0f64; 3]; 4];
    for (row, p) in [a, b, c, d].into_iter().enumerate() {
        for (axis, (x, y)) in [(p.x, e.x), (p.y, e.y), (p.z, e.z)].into_iter().enumerate() {
            let (difference, error) = two_diff(x, y);
            if error != 0.0 || !difference.is_finite() {
                return None;
            }
            differences[row][axis] = difference;
        }
    }
    let nonzero = || differences.iter().flatten().filter(|v| **v != 0.0);
    let Some(lowest) = nonzero().map(|v| least_significant_bit_exponent(*v)).min() else {
        return Some(Sign::Zero);
    };
    let highest = nonzero().map(|v| highest_bit_exponent(*v)).max()?;
    if highest - lowest >= 20 {
        return None;
    }
    let scaled = |v: f64| -> i128 {
        if v == 0.0 {
            return 0;
        }
        // |v| = significand * 2^exponent exactly; shift onto the 2^lowest grid.
        let bits = v.abs().to_bits();
        let encoded = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1u64 << 52) - 1);
        let (significand, exponent) = if encoded == 0 {
            (fraction, -1074)
        } else {
            (fraction | (1u64 << 52), encoded - 1075)
        };
        let shift = exponent - lowest;
        let magnitude = if shift >= 0 {
            i128::from(significand) << shift
        } else {
            i128::from(significand >> -shift)
        };
        if v < 0.0 {
            -magnitude
        } else {
            magnitude
        }
    };
    let row = |r: usize| {
        let [x, y, z] = differences[r].map(scaled);
        [x, y, z, x * x + y * y + z * z]
    };
    let m = [row(0), row(1), row(2), row(3)];
    let det3 = |p: [i128; 4], q: [i128; 4], r: [i128; 4]| {
        p[0] * (q[1] * r[2] - r[1] * q[2]) - p[1] * (q[0] * r[2] - r[0] * q[2])
            + p[2] * (q[0] * r[1] - r[0] * q[1])
    };
    // The lifted determinant with the sign convention of `insphere`.
    let total = -m[0][3] * det3(m[1], m[2], m[3]) + m[1][3] * det3(m[0], m[2], m[3])
        - m[2][3] * det3(m[0], m[1], m[3])
        + m[3][3] * det3(m[0], m[1], m[2]);
    Some(match total.signum() {
        1 => Sign::Positive,
        -1 => Sign::Negative,
        _ => Sign::Zero,
    })
}

/// Multiply an expansion by `x*x + y*y + z*z`, exactly.
#[must_use]
fn lift3(e: &[f64], p: &[Vec<f64>; 3]) -> Vec<f64> {
    let square = expansion_sum(
        &expansion_sum(
            &expansion_product(&p[0], &p[0]),
            &expansion_product(&p[1], &p[1]),
        ),
        &expansion_product(&p[2], &p[2]),
    );
    expansion_product(e, &square)
}

/// Smallest exponent `k` with `2^k` a coordinate magnitude the exact path of
/// [`in_diametral_sphere`] accepts; zero is always accepted.
const DIAMETRAL_MIN_EXPONENT: i32 = -100;

/// Largest such exponent.
const DIAMETRAL_MAX_EXPONENT: i32 = 100;

/// Is `d` inside the diametral sphere of the triangle `a`, `b`, `c`?
///
/// The diametral sphere is the smallest sphere through `a`, `b`, `c`: its
/// centre is the triangle's circumcentre and its equator the circumcircle.
/// For `d` coplanar with the triangle this is therefore the coplanar
/// `incircle` test in 3D -- the question a 3D Delaunay triangulation asks of
/// a point lying exactly in the plane of a convex-hull face -- and, unlike
/// [`insphere`], it needs no orientation: [`Sign::Positive`] means strictly
/// inside, [`Sign::Negative`] strictly outside, [`Sign::Zero`] on the sphere.
///
/// The sign is that of
///
/// ```text
/// -( |n|^2 (w.w) - |u|^2 ((w x v).n) - |v|^2 ((u x w).n) )
/// ```
///
/// with `u = b - a`, `v = c - a`, `w = d - a` and `n = u x v`: `|n|^2` times
/// the power of `d` with respect to the sphere, a degree-6 polynomial in the
/// coordinates.
///
/// A degenerate triangle (collinear `a`, `b`, `c`) has no diametral sphere and
/// returns `Zero`. A forward running-error filter decides almost every case;
/// the exact expansion path decides the rest. The exact path needs every
/// non-zero coordinate magnitude in `[2^-100, 2^100]`, so that no product in
/// the degree-6 determinant leaves binary64's normal range; outside it, and
/// for non-finite input, the result is [`Certified::Uncertain`] rather than a
/// guess.
#[must_use]
pub fn in_diametral_sphere(a: Point3, b: Point3, c: Point3, d: Point3) -> Certified {
    match in_diametral_sphere_filter(a, b, c, d) {
        Certified::Certain { sign, .. } => Certified::exact_sign(sign),
        _ => {
            if [a, b, c, d].iter().all(|p| {
                [p.x, p.y, p.z].iter().all(|&x| {
                    x == 0.0
                        || (x.is_finite()
                            && x.abs() >= 2f64.powi(DIAMETRAL_MIN_EXPONENT)
                            && x.abs() <= 2f64.powi(DIAMETRAL_MAX_EXPONENT))
                })
            }) {
                Certified::exact_sign(in_diametral_sphere_exact(a, b, c, d))
            } else {
                Certified::Uncertain {
                    attempted: Precision::Exact,
                }
            }
        }
    }
}

/// A value and a bound on its absolute error, for a running-error filter.
#[derive(Clone, Copy)]
struct Bounded {
    value: f64,
    error: f64,
}

/// Unit roundoff, taken at twice its size so the `1 / (1 - eps)` factor of
/// the backward bound is absorbed.
const ROUNDOFF: f64 = f64::EPSILON;

impl Bounded {
    fn difference(p: f64, q: f64) -> Self {
        let value = p - q;
        Self {
            value,
            error: ROUNDOFF * value.abs(),
        }
    }

    fn add(self, other: Self) -> Self {
        let value = self.value + other.value;
        Self {
            value,
            error: self.error + other.error + ROUNDOFF * value.abs(),
        }
    }

    fn sub(self, other: Self) -> Self {
        let value = self.value - other.value;
        Self {
            value,
            error: self.error + other.error + ROUNDOFF * value.abs(),
        }
    }

    fn mul(self, other: Self) -> Self {
        let value = self.value * other.value;
        Self {
            value,
            error: self.value.abs() * other.error
                + other.value.abs() * self.error
                + self.error * other.error
                + ROUNDOFF * value.abs(),
        }
    }
}

/// The running-error filter alone, exposed so escalation can be measured.
///
/// Every operation carries a bound on its absolute error; the bound itself is
/// computed in rounded arithmetic, so the final bound is inflated by a
/// relative margin far larger than the rounding the bound's own evaluation
/// can incur. A non-zero coordinate difference outside `[2^-150, 2^150]`
/// (whose degree-6 products could leave the normal range) or a non-finite
/// value is left [`Certified::Uncertain`].
#[must_use]
pub fn in_diametral_sphere_filter(a: Point3, b: Point3, c: Point3, d: Point3) -> Certified {
    let vector = |p: Point3| {
        [
            Bounded::difference(p.x, a.x),
            Bounded::difference(p.y, a.y),
            Bounded::difference(p.z, a.z),
        ]
    };
    let (u, v, w) = (vector(b), vector(c), vector(d));
    let cross = |p: [Bounded; 3], q: [Bounded; 3]| {
        [
            p[1].mul(q[2]).sub(p[2].mul(q[1])),
            p[2].mul(q[0]).sub(p[0].mul(q[2])),
            p[0].mul(q[1]).sub(p[1].mul(q[0])),
        ]
    };
    let dot =
        |p: [Bounded; 3], q: [Bounded; 3]| p[0].mul(q[0]).add(p[1].mul(q[1])).add(p[2].mul(q[2]));
    let n = cross(u, v);
    let power = dot(n, n)
        .mul(dot(w, w))
        .sub(dot(u, u).mul(dot(cross(w, v), n)))
        .sub(dot(v, v).mul(dot(cross(u, w), n)));

    let tiny = 2f64.powi(-900);
    let operands = [u, v, w].into_iter().flatten().map(|x| x.value.abs());
    let representable = operands.clone().all(f64::is_finite)
        && operands
            .filter(|x| *x != 0.0)
            .all(|x| x >= 2f64.powi(-150) && x <= 2f64.powi(150))
        && power.value.is_finite()
        && power.error.is_finite();
    if !representable {
        return Certified::Uncertain {
            attempted: Precision::F64,
        };
    }
    // The bound's own evaluation rounds a few dozen times; 2^-40 relative
    // dwarfs that, and the absolute term covers underflow in the bound.
    let bound = power.error * (1.0 + 2f64.powi(-40)) + tiny;
    Certified::from_filter(-power.value, bound, Precision::F64)
}

/// Exact sign of the diametral-sphere determinant, over expansions.
#[must_use]
fn in_diametral_sphere_exact(a: Point3, b: Point3, c: Point3, d: Point3) -> Sign {
    let vector = |p: Point3| [diff(p.x, a.x), diff(p.y, a.y), diff(p.z, a.z)];
    let (u, v, w) = (vector(b), vector(c), vector(d));
    let cross = |p: &[Vec<f64>; 3], q: &[Vec<f64>; 3]| {
        [
            minor(&p[1], &q[2], &p[2], &q[1]),
            minor(&p[2], &q[0], &p[0], &q[2]),
            minor(&p[0], &q[1], &p[1], &q[0]),
        ]
    };
    let dot = |p: &[Vec<f64>; 3], q: &[Vec<f64>; 3]| {
        expansion_sum(
            &expansion_sum(
                &expansion_product(&p[0], &q[0]),
                &expansion_product(&p[1], &q[1]),
            ),
            &expansion_product(&p[2], &q[2]),
        )
    };
    let n = cross(&u, &v);
    let first = expansion_product(&dot(&n, &n), &dot(&w, &w));
    let second = expansion_product(&dot(&u, &u), &dot(&cross(&w, &v), &n));
    let third = expansion_product(&dot(&v, &v), &dot(&cross(&u, &w), &n));
    let power = expansion_sum(&first, &negate_expansion(&expansion_sum(&second, &third)));
    expansion_sign(&power).flip()
}

#[cfg(test)]
mod tests {
    use super::{insphere_expansion, insphere_small_integer};
    use axiolid_core::Point3;

    fn next(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    /// The integer tier agrees with the expansion tier wherever it answers,
    /// non-zero signs included (the filter usually settles those, so the
    /// public predicate alone would not exercise them), and it declines
    /// every input whose scaled integers could overflow `i128`.
    #[test]
    fn the_integer_tier_agrees_with_the_expansions() {
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let (mut answered, mut declined) = (0usize, 0usize);
        for round in 0..40_000u32 {
            // Magnitudes from a few units up to 2^30, on grids of 2^-8..2^8.
            let bits = 2 + round % 29;
            let unit = 2f64.powi((next(&mut state) % 17) as i32 - 8);
            let mut coordinate = || {
                let span = 1u64 << bits;
                ((next(&mut state) % (2 * span)) as f64 - span as f64) * unit
            };
            let mut point = || Point3::new(coordinate(), coordinate(), coordinate());
            let (a, b, c, d, e) = (point(), point(), point(), point(), point());
            let expected = insphere_expansion(a, b, c, d, e);
            match insphere_small_integer(a, b, c, d, e) {
                Some(sign) => {
                    assert_eq!(sign, expected, "{a:?} {b:?} {c:?} {d:?} {e:?}");
                    answered += 1;
                }
                None => declined += 1,
            }
        }
        assert!(
            answered > 10_000 && declined > 10_000,
            "{answered} / {declined}"
        );
    }
}
