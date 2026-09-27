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

    let ab = ax * by - bx * ay;
    let bc = bx * cy - cx * by;
    let cd = cx * dy - dx * cy;
    let da = dx * ay - ax * dy;
    let ac = ax * cy - cx * ay;
    let bd = bx * dy - dx * by;

    let abc = az * bc - bz * ac + cz * ab;
    let bcd = bz * cd - cz * bd + dz * bc;
    let cda = cz * da + dz * ac + az * cd;
    let dab = dz * ab + az * bd + bz * da;

    let alift = ax * ax + ay * ay + az * az;
    let blift = bx * bx + by * by + bz * bz;
    let clift = cx * cx + cy * cy + cz * cz;
    let dlift = dx * dx + dy * dy + dz * dz;

    let determinant = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    let permanent =
        (abc.abs() * dlift + dab.abs() * clift) + (cda.abs() * blift + bcd.abs() * alift);

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
