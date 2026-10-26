//! When two boundary edges meet: the validation contract's notion of
//! touching, which settling repairs by the same tests (#272).
//!
//! Non-adjacent edges meet when an end of one lies on the other within
//! the tolerance ([`on_segment`]), or when they cross properly, each with
//! one end strictly on each side of the other's line
//! ([`proper_crossing`]). Every decision is exact: the cross products are
//! signed by `axiolid_exact::certify`, and the box tests compare exact
//! differences, so a pair of edges decides the same wherever it lies.
//! Reading the side of a rounded cross product made an exact `0.0` a side:
//! two disjoint collinear edges of a settled ring then counted as
//! crossing, and `overlay` refused its own output (#272).
//!
//! `scripts/probe_edge_contact_mutants.py` lists the faults the tests
//! here and in `tests/edge_contact.rs` must catch.

use axiolid_core::Point2;
use axiolid_exact::{certify, Arith, SignExpr};
use axiolid_guarantees::Sign;

/// The sign of `(b - a) x (c - a) - offset`, exactly.
struct Cross {
    a: Point2,
    b: Point2,
    c: Point2,
    offset: f64,
}

impl SignExpr for Cross {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let (ux, uy) = (f(self.b.x).sub(&f(self.a.x)), f(self.b.y).sub(&f(self.a.y)));
        let (vx, vy) = (f(self.c.x).sub(&f(self.a.x)), f(self.c.y).sub(&f(self.a.y)));
        ux.mul(&vy).sub(&uy.mul(&vx)).sub(&f(self.offset)).sign()
    }
}

fn cross_sign(a: Point2, b: Point2, c: Point2, offset: f64) -> Sign {
    // Finite points never make the exact tier refuse.
    certify(&Cross { a, b, c, offset }).unwrap_or(Sign::Zero)
}

/// Which side of the directed line `a -> b` the point `c` lies on,
/// exactly (#272): `Positive` to the left.
pub(crate) fn side(a: Point2, b: Point2, c: Point2) -> Sign {
    cross_sign(a, b, c, 0.0)
}

/// Whether `p` lies on the segment `a`-`b` within `epsilon`: the
/// validation contract's notion of touching, which settling shares.
///
/// `|(b - a) x (p - a)| <= epsilon` and `p` inside the segment's box
/// widened by `epsilon`, both decided exactly (#272), so the answer does
/// not depend on where the segment lies. A point exactly on the line
/// through the segment touches it only inside its extent: two collinear
/// edges of a U-shape or a comb share a line without touching.
pub(crate) fn on_segment(a: Point2, b: Point2, p: Point2, epsilon: f64) -> bool {
    within_extent(a, b, p, epsilon)
        && cross_sign(a, b, p, epsilon) != Sign::Positive
        && cross_sign(a, b, p, -epsilon) != Sign::Negative
}

/// Whether `a`-`b` and `c`-`d` cross at one point inside both, exactly:
/// each segment has one end strictly on each side of the other's line
/// (#272). An end exactly on the other's line is no crossing here but
/// [`on_segment`]'s case; a rounded zero read as a side made two disjoint
/// collinear edges cross.
pub(crate) fn proper_crossing(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let opposite = |p: Sign, q: Sign| {
        matches!(
            (p, q),
            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
        )
    };
    opposite(side(a, b, c), side(a, b, d)) && opposite(side(c, d, a), side(c, d, b))
}

/// Whether non-adjacent boundary segments `a`-`b` and `c`-`d` meet within
/// `epsilon`: an end of one on the other ([`on_segment`]), or a proper
/// crossing ([`proper_crossing`]). Every decision is exact (#272), so a
/// pair decides the same wherever it lies. Collinear pairs fall to the
/// first case: overlapping or touching ones put an end on the other and
/// meet; disjoint ones have no end inside the other's extent and do not,
/// whatever the rounded cross products of their ends say.
pub(crate) fn segments_intersect(a: Point2, b: Point2, c: Point2, d: Point2, epsilon: f64) -> bool {
    // Boxes apart by more than `epsilon` cannot meet either way; this only
    // saves the exact signs.
    let apart = |p: f64, q: f64, r: f64, s: f64| {
        !at_least(r.max(s), p.min(q), -epsilon) || !at_least(p.max(q), r.min(s), -epsilon)
    };
    if apart(a.x, b.x, c.x, d.x) || apart(a.y, b.y, c.y, d.y) {
        return false;
    }
    on_segment(c, d, a, epsilon)
        || on_segment(c, d, b, epsilon)
        || on_segment(a, b, c, epsilon)
        || on_segment(a, b, d, epsilon)
        || proper_crossing(a, b, c, d)
}

/// Whether `x - y >= bound`, exactly: the rounded difference decides
/// unless it rounds onto `bound`, and then the sign of its rounding error
/// does. Inputs must be finite.
fn at_least(x: f64, y: f64, bound: f64) -> bool {
    let d = x - y;
    if d != bound {
        return d > bound;
    }
    // Two-sum: `x - y == d + error` exactly.
    let ny = -y;
    let v = d - x;
    let error = (x - (d - v)) + (ny - v);
    error >= 0.0
}

/// Whether `p` lies inside the axis-aligned box of segment `a`-`b`, widened
/// by `epsilon`, exactly. Combined with collinearity this puts `p` on the
/// segment. Written as `p.x >= min - epsilon` it would round `epsilon`
/// away far from the origin (#272).
fn within_extent(a: Point2, b: Point2, p: Point2, epsilon: f64) -> bool {
    at_least(p.x, a.x.min(b.x), -epsilon)
        && at_least(a.x.max(b.x), p.x, -epsilon)
        && at_least(p.y, a.y.min(b.y), -epsilon)
        && at_least(a.y.max(b.y), p.y, -epsilon)
}

#[cfg(test)]
mod tests {
    //! Every point lies on the grid of `2^-24` near the origin, so its
    //! coordinates times `2^24` are integers below `2^30`, and the cross
    //! products of their differences, of up to 61 bits, round in `f64`.
    //! Offsets of 1e5 and 5e6 translate such points exactly, and an
    //! `i128` oracle decides each pair independently of the code under
    //! test.
    use super::*;

    const SCALE: f64 = (1u64 << 24) as f64;
    const OFFSETS: [(f64, f64); 4] = [(0.0, 0.0), (1e5, -1e5), (5e6, 5e6), (-5e6, 1e5)];

    type I = (i64, i64);

    fn point((x, y): I, (dx, dy): (f64, f64)) -> Point2 {
        let p = Point2::new(x as f64 / SCALE + dx, y as f64 / SCALE + dy);
        assert_eq!(
            ((p.x - dx) * SCALE, (p.y - dy) * SCALE),
            (x as f64, y as f64),
            "the grid point must translate exactly"
        );
        p
    }

    fn cross_i(a: I, b: I, c: I) -> i128 {
        let (ux, uy) = (i128::from(b.0 - a.0), i128::from(b.1 - a.1));
        let (vx, vy) = (i128::from(c.0 - a.0), i128::from(c.1 - a.1));
        ux * vy - uy * vx
    }

    /// The oracle: whether `p` is within `cross_units` (a cross product, in
    /// grid units squared) of segment `a`-`b`'s line and inside its box
    /// widened by `box_units`.
    fn on_i(a: I, b: I, p: I, cross_units: i128, box_units: i64) -> bool {
        cross_i(a, b, p).abs() <= cross_units
            && p.0 >= a.0.min(b.0) - box_units
            && p.0 <= a.0.max(b.0) + box_units
            && p.1 >= a.1.min(b.1) - box_units
            && p.1 <= a.1.max(b.1) + box_units
    }

    fn proper_i(a: I, b: I, c: I, d: I) -> bool {
        cross_i(a, b, c).signum() * cross_i(a, b, d).signum() < 0
            && cross_i(c, d, a).signum() * cross_i(c, d, b).signum() < 0
    }

    fn meet_i(a: I, b: I, c: I, d: I, cross_units: i128, box_units: i64) -> bool {
        on_i(c, d, a, cross_units, box_units)
            || on_i(c, d, b, cross_units, box_units)
            || on_i(a, b, c, cross_units, box_units)
            || on_i(a, b, d, cross_units, box_units)
            || proper_i(a, b, c, d)
    }

    /// The strict branch as it was before #272: the side of a rounded
    /// cross product, `0.0` counted with the negatives.
    fn rounded_crossing(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
        let cross = |a: Point2, b: Point2, c: Point2| (b - a).perp_dot(c - a);
        (cross(a, b, c) > 0.0) != (cross(a, b, d) > 0.0)
            && (cross(c, d, a) > 0.0) != (cross(c, d, b) > 0.0)
    }

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 11
        }
        fn range(&mut self, n: i64) -> i64 {
            (self.next() % (2 * n as u64 + 1)) as i64 - n
        }
    }

    /// The point at parameter `t` (in 1/1024ths) along `a`-`b`, rounded to
    /// the grid, moved by `jiggle` units.
    fn along(a: I, b: I, t: i64, jiggle: I) -> I {
        (
            a.0 + (b.0 - a.0) * t / 1024 + jiggle.0,
            a.1 + (b.1 - a.1) * t / 1024 + jiggle.1,
        )
    }

    const COUNT: usize = 4000;

    /// Pairs near every degenerate case: the second segment's ends on, or
    /// a few grid steps off, the first one's line, before, at the ends of,
    /// inside and beyond its extent.
    fn cases() -> Vec<(I, I, I, I)> {
        let mut rng = Lcg(272);
        let big = 1i64 << 27;
        let params = [-700, -300, 0, 300, 512, 700, 1024, 1300, 1800];
        let mut out = Vec::new();
        for _ in 0..COUNT {
            let a = (rng.range(big), rng.range(big));
            let b = (rng.range(big), rng.range(big));
            let end = |rng: &mut Lcg| {
                let t = params[(rng.range(4) + 4) as usize];
                along(a, b, t, (rng.range(2), rng.range(2)))
            };
            let c = end(&mut rng);
            let d = if rng.next().is_multiple_of(3) {
                (rng.range(big), rng.range(big))
            } else {
                end(&mut rng)
            };
            out.push((a, b, c, d));
        }
        out
    }

    #[test]
    fn segments_meet_exactly_as_the_integer_oracle_says_at_every_offset() {
        let (mut rounded_wrong, mut meeting, mut proper) = (0, 0, 0);
        for (a, b, c, d) in cases() {
            for eps_units in [0i64, 1, 16] {
                let epsilon = eps_units as f64 / SCALE;
                // `|cross| <= epsilon` in grid units squared.
                let cross_units = i128::from(eps_units) << 24;
                let want = meet_i(a, b, c, d, cross_units, eps_units);
                meeting += usize::from(want);
                for offset in OFFSETS {
                    let (pa, pb) = (point(a, offset), point(b, offset));
                    let (pc, pd) = (point(c, offset), point(d, offset));
                    assert_eq!(
                        segments_intersect(pa, pb, pc, pd, epsilon),
                        want,
                        "{a:?}-{b:?} vs {c:?}-{d:?} at {offset:?}, epsilon {eps_units}"
                    );
                    if eps_units == 0 {
                        let exact = proper_i(a, b, c, d);
                        assert_eq!(proper_crossing(pa, pb, pc, pd), exact);
                        proper += usize::from(exact);
                        rounded_wrong += usize::from(rounded_crossing(pa, pb, pc, pd) != exact);
                    }
                }
            }
        }
        // Hard enough that rounded signs get some wrong, and varied enough
        // that every answer occurs often.
        assert!(rounded_wrong > 0, "no case defeats the rounded signs");
        assert!(proper > 400, "{proper}");
        assert!(meeting > 1000 && meeting < 3 * COUNT - 1000, "{meeting}");
    }

    #[test]
    fn collinear_edges_meet_only_when_their_extents_do() {
        // Exactly on one slanted line: the pieces for parameters 0..1 and
        // 2..3 are disjoint, 0..2 and 1..3 overlap, 0..1 and 1..2 share an
        // end.
        let (o, u) = ((3, -5), (40_000_001, 29_999_999));
        let at = |k: i64| (o.0 + k * u.0, o.1 + k * u.1);
        for offset in OFFSETS {
            let p = |k: i64| point(at(k), offset);
            assert!(!segments_intersect(p(0), p(1), p(2), p(3), 0.0));
            assert!(!segments_intersect(p(0), p(1), p(3), p(2), 1e-9));
            assert!(segments_intersect(p(0), p(2), p(1), p(3), 0.0));
            assert!(segments_intersect(p(0), p(1), p(1), p(2), 0.0));
            assert!(segments_intersect(p(1), p(0), p(2), p(1), 0.0));
            assert!(on_segment(p(0), p(2), p(1), 0.0));
            assert!(!on_segment(p(0), p(1), p(2), 0.0));
            assert!(!proper_crossing(p(0), p(2), p(1), p(3)));
        }
    }

    #[test]
    fn near_crossings_a_few_ulps_either_side_decide_exactly() {
        // `c`-`d` starts above the line of `a`-`b`, and `d` ends `k` ulps
        // of `y` above the point `m` of that line: the line runs to the
        // right, so `d` is left of it for `k > 0`, on it for `k == 0` and
        // right of it, crossing, for `k < 0`.
        let (a, b) = ((-3, 5), (100_000_005, 70_000_009));
        let m = along(a, b, 512, (0, 0));
        assert_eq!(cross_i(a, b, m), 0);
        let c = (m.0 - 9, m.1 + 13);
        assert!(cross_i(a, b, c) > 0);
        let nudge = |y: f64, k: i64| {
            let mut y = y;
            for _ in 0..k.abs() {
                y = if k > 0 { y.next_up() } else { y.next_down() };
            }
            y
        };
        for offset in OFFSETS {
            let (pa, pb, pc, pm) = (
                point(a, offset),
                point(b, offset),
                point(c, offset),
                point(m, offset),
            );
            for k in -3i64..=3 {
                let pd = Point2::new(pm.x, nudge(pm.y, k));
                let at = format!("k {k} at {offset:?}");
                let want = [Sign::Negative, Sign::Zero, Sign::Positive][(k.signum() + 1) as usize];
                assert_eq!(side(pa, pb, pd), want, "{at}");
                assert_eq!(proper_crossing(pa, pb, pc, pd), k < 0, "{at}");
                assert_eq!(segments_intersect(pa, pb, pc, pd, 0.0), k <= 0, "{at}");
            }
        }
    }

    #[test]
    fn differences_are_compared_exactly() {
        // `1 - 1e-17` rounds to `1`; the rounding error still decides.
        let y = 1e-17;
        assert!(at_least(1.0, y, 1.0 - f64::EPSILON / 2.0));
        assert!(!at_least(1.0, y, 1.0));
        assert!(at_least(1.0, -y, 1.0));
        // Far from the origin the tolerance is not rounded away.
        assert!(!at_least(5e6, 5e6 + 1e-9, -1e-15));
        assert!(at_least(5e6, 5e6, -1e-15));
    }

    #[test]
    fn the_widened_box_is_exact_far_from_the_origin() {
        // At 5e6 the coordinates step by 2^-30, and a tolerance of
        // `m * 2^-33` lies between those steps: `min - epsilon` rounds.
        let step = (2.0f64).powi(-30);
        let (a, b) = (Point2::new(5e6, 5e6), Point2::new(5e6 + 1.0, 5e6));
        for m in 1..40 {
            let epsilon = f64::from(m) * (2.0f64).powi(-33);
            for k in -6i32..=0 {
                let p = Point2::new(5e6 + f64::from(k) * step, 5e6);
                let inside = 8 * k >= -m;
                assert_eq!(within_extent(a, b, p, epsilon), inside, "k {k}, m {m}");
                assert_eq!(on_segment(a, b, p, epsilon), inside, "k {k}, m {m}");
                let q = Point2::new(p.x, p.y + 1.0);
                assert_eq!(
                    segments_intersect(p, q, a, b, epsilon),
                    inside,
                    "k {k}, m {m}"
                );
            }
        }
    }

    /// `(x, y)` with `p * y - q * x == 1` for coprime `p`, `q`, moved along
    /// `(p, q)` to near the foot of the origin.
    fn unit_step(p: i64, q: i64) -> I {
        let (mut r0, mut r1) = (p, q);
        let (mut s0, mut s1) = (1i64, 0i64);
        let (mut t0, mut t1) = (0i64, 1i64);
        while r1 != 0 {
            let k = r0 / r1;
            (r0, r1) = (r1, r0 - k * r1);
            (s0, s1) = (s1, s0 - k * s1);
            (t0, t1) = (t1, t0 - k * t1);
        }
        let (x, y) = (-t0 * r0, s0 * r0);
        let k = ((x * p + y * q) as f64 / (p * p + q * q) as f64).round() as i64;
        (x - k * p, y - k * q)
    }

    fn gcd(a: i64, b: i64) -> i64 {
        if b == 0 {
            a.abs()
        } else {
            gcd(b, a % b)
        }
    }

    #[test]
    fn the_tolerance_on_the_cross_product_is_exact() {
        // `h` is `r` grid units squared of cross product off the line of
        // a long segment, half way along it; the tolerance is one such
        // unit. The rounded cross product of such a point is mostly
        // rounding residue.
        let unit = 1.0 / (SCALE * SCALE);
        let mut rng = Lcg(0x272);
        let (mut tried, mut rounded_wrong) = (0, 0);
        while tried < 300 {
            let (p, q) = ((1 << 27) + rng.range(1 << 26), rng.range(1 << 27));
            if gcd(p, q) != 1 {
                continue;
            }
            tried += 1;
            let (x0, y0) = unit_step(p, q);
            let (a, b) = ((-p, -q), (p, q));
            for r in -3i64..=3 {
                let h = (r * x0, r * y0);
                assert_eq!(cross_i(a, b, h), 2 * i128::from(r));
                let near = cross_i(a, b, h).abs() <= 2;
                let (pa, pb, ph) = (
                    point(a, (0.0, 0.0)),
                    point(b, (0.0, 0.0)),
                    point(h, (0.0, 0.0)),
                );
                assert_eq!(
                    on_segment(pa, pb, ph, 2.0 * unit),
                    near,
                    "{a:?} {b:?} {h:?}"
                );
                let rounded = (pb - pa).perp_dot(ph - pa).abs() <= 2.0 * unit;
                rounded_wrong += usize::from(rounded != near);
            }
        }
        assert!(
            rounded_wrong > 0,
            "no point defeats the rounded cross product"
        );
    }
}
