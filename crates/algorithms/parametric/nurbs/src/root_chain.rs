//! Roots of a function along a parameter, certified with their multiplicity
//! by a Rolle chain (#397, #427, ADR 0077).
//!
//! The function is given by a [`Tower`]: enclosures of `H_0, H_1, ...` over
//! a parameter range, where each `H_{j+1}` is, up to a factor that never
//! vanishes, the derivative of `H_j`. A zero of `H_j` of multiplicity `m`
//! is then a zero of `H_{j+1}` of multiplicity `m - 1`, and if `H_k` has no
//! zero on a stretch, `H_j` has at most `k - j` there, counting
//! multiplicity.
//!
//! [`isolate`] splits the range until each leaf has an order `k` (at most
//! [`MAX_CONTACT_ORDER`]) whose `H_k` excludes zero, and [`chain`] climbs
//! down from `H_k`: between two zones of `H_{j+1}` the function `H_j` is
//! strictly monotone, so a sign change brackets exactly one simple zero,
//! and a zone whose enclosure straddles zero holds a zero of the
//! multiplicity the chain allows. Split points are signed at orders 0 to
//! [`GENERIC_ORDER`], so no zero sits on the end of a leaf: an end that
//! straddles zero refuses, except the curve's own ends, which report.
//! Everything that cannot be settled (a contact above the cap, one that
//! rounding hides, an exhausted budget) refuses by name and never guesses.
//!
//! This module is pure: it knows nothing of curves, and is tested against
//! polynomials.

use crate::interval::I;
use axiolid_core::Scalar;

/// The highest order of contact (multiplicity of a root) certified.
pub(crate) const MAX_CONTACT_ORDER: usize = 8;
/// The orders that must be signed at a split point.
pub(crate) const GENERIC_ORDER: usize = 2;
/// Nodes one isolation may visit before it refuses: the most any test of
/// the section, boolean and compile suites needed was 449, and this is the
/// first power of two (at least `2^14`) above sixteen times that.
pub(crate) const ROOT_BUDGET: usize = 1 << 14;
const MAX_DEPTH: u32 = 128;
/// The narrowest leaf, in parameter units.
const MIN_WIDTH: Scalar = 1.0 / (1u64 << 44) as Scalar;
/// Split points, as fractions of a stretch, middle first.
const SPLITS: [Scalar; 15] = [
    8.0 / 16.0,
    7.0 / 16.0,
    9.0 / 16.0,
    6.0 / 16.0,
    10.0 / 16.0,
    5.0 / 16.0,
    11.0 / 16.0,
    4.0 / 16.0,
    12.0 / 16.0,
    3.0 / 16.0,
    13.0 / 16.0,
    2.0 / 16.0,
    14.0 / 16.0,
    1.0 / 16.0,
    15.0 / 16.0,
];

/// Why the roots could not all be isolated with their multiplicity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RootRefusal {
    /// The work budget ran out.
    Budget,
    /// A contact above the cap, one that rounding hides, or an end of a
    /// stretch that sits on a root.
    Unresolved,
}

/// What an end of a stretch is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum End {
    /// The curve's own end: a root there is reported.
    Curve,
    /// A cut made by the isolation: no root may sit on it.
    Inner,
}

/// A root, as the range that holds it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Zone {
    pub(crate) lo: Scalar,
    pub(crate) hi: Scalar,
    /// The point reported for the root.
    pub(crate) at: Scalar,
    /// Its multiplicity (0 for the end of a stretch that is no root).
    pub(crate) m: usize,
    /// Whether the zone reaches an end of the stretch.
    pub(crate) end: Option<End>,
}

impl Zone {
    fn end_point(t: Scalar, kind: End) -> Self {
        Self {
            lo: t,
            hi: t,
            at: t,
            m: 0,
            end: Some(kind),
        }
    }

    /// One zone of the run of zones whose enclosures straddle zero.
    fn cluster(run: &[Zone]) -> Self {
        let mut at = &run[0];
        for q in run {
            if q.m > at.m {
                at = q;
            }
        }
        let curve = run.iter().any(|q| q.end == Some(End::Curve));
        let m = run.iter().map(|q| q.m).sum::<usize>() + 1;
        Self {
            lo: run[0].lo,
            hi: run[run.len() - 1].hi,
            at: at.at,
            m,
            end: curve.then_some(End::Curve),
        }
    }
}

/// A count of the nodes still allowed.
#[derive(Debug, Clone)]
pub(crate) struct Budget {
    left: usize,
}

impl Budget {
    pub(crate) fn new(nodes: usize) -> Self {
        Self { left: nodes }
    }

    fn spend(&mut self) -> Result<(), RootRefusal> {
        if self.left == 0 {
            return Err(RootRefusal::Budget);
        }
        self.left -= 1;
        Ok(())
    }
}

/// A function along a parameter, with its chain of derivatives.
pub(crate) trait Tower {
    /// Enclosures of `H_0..=H_order` along the curve, valid at every `t` in
    /// `[t0, t1]` (`t0 == t1`: one point). `None` where the function cannot
    /// be bounded.
    fn enclose(&self, t0: Scalar, t1: Scalar, order: usize) -> Option<Vec<I>>;

    /// Preferred split points strictly inside `(t0, t1)`, ascending.
    fn joins(&self, t0: Scalar, t1: Scalar) -> Vec<Scalar>;
}

/// Every root in `[t0, t1]` with its multiplicity, ascending.
pub(crate) fn isolate<T: Tower>(
    tower: &T,
    t0: Scalar,
    t1: Scalar,
    ends: [End; 2],
    budget: &mut Budget,
) -> Result<Vec<Zone>, RootRefusal> {
    let mut out = Vec::new();
    isolate_in(tower, t0, t1, ends, 0, budget, &mut out)?;
    Ok(out)
}

fn isolate_in<T: Tower>(
    tower: &T,
    t0: Scalar,
    t1: Scalar,
    ends: [End; 2],
    depth: u32,
    budget: &mut Budget,
    out: &mut Vec<Zone>,
) -> Result<(), RootRefusal> {
    budget.spend()?;
    let middle = 0.5 * (t0 + t1);
    let mut joins = tower.joins(t0, t1);
    joins.sort_by(|a, b| (a - middle).abs().total_cmp(&(b - middle).abs()));
    if let Some(&j) = joins.iter().find(|&&j| generic(tower, j)) {
        isolate_in(tower, t0, j, [ends[0], End::Inner], depth + 1, budget, out)?;
        return isolate_in(tower, j, t1, [End::Inner, ends[1]], depth + 1, budget, out);
    }
    let low = tower.enclose(t0, t1, 1).ok_or(RootRefusal::Unresolved)?;
    if low[0].excludes_zero() {
        return Ok(());
    }
    let k = if low[1].excludes_zero() {
        Some(1)
    } else {
        let mut e = tower.enclose(t0, t1, 2).ok_or(RootRefusal::Unresolved)?;
        if !e[2].excludes_zero() {
            e = tower
                .enclose(t0, t1, MAX_CONTACT_ORDER)
                .ok_or(RootRefusal::Unresolved)?;
        }
        (2..=MAX_CONTACT_ORDER).find(|&k| e[k].excludes_zero())
    };
    if let Some(k) = k {
        out.extend(chain(tower, t0, t1, k, ends)?);
        return Ok(());
    }
    if depth >= MAX_DEPTH || t1 - t0 <= MIN_WIDTH {
        return Err(RootRefusal::Unresolved);
    }
    for c in SPLITS {
        let m = t0 + (t1 - t0) * c;
        if generic(tower, m) {
            isolate_in(tower, t0, m, [ends[0], End::Inner], depth + 1, budget, out)?;
            return isolate_in(tower, m, t1, [End::Inner, ends[1]], depth + 1, budget, out);
        }
    }
    Err(RootRefusal::Unresolved)
}

/// Whether `t` is a point where orders 0 to [`GENERIC_ORDER`] are signed.
fn generic<T: Tower>(tower: &T, t: Scalar) -> bool {
    tower
        .enclose(t, t, GENERIC_ORDER)
        .is_some_and(|e| e.iter().all(|x| x.excludes_zero()))
}

fn positive(e: I) -> bool {
    e.lo > 0.0
}

/// The roots of `H_0` on `[a, b]`, given that `H_k` has none there.
fn chain<T: Tower>(
    tower: &T,
    a: Scalar,
    b: Scalar,
    k: usize,
    ends: [End; 2],
) -> Result<Vec<Zone>, RootRefusal> {
    let mut z: Vec<Zone> = Vec::new();
    for j in (0..k).rev() {
        let mut p: Vec<Zone> = Vec::with_capacity(z.len() + 2);
        if z.first().is_none_or(|q| q.lo > a) {
            p.push(Zone::end_point(a, ends[0]));
        }
        p.extend(z.iter().copied());
        if z.last().is_none_or(|q| q.hi < b) {
            p.push(Zone::end_point(b, ends[1]));
        }
        let mut e = Vec::with_capacity(p.len());
        for q in &p {
            let level = tower
                .enclose(q.lo, q.hi, j)
                .ok_or(RootRefusal::Unresolved)?;
            e.push(level[j]);
        }
        let straddles: Vec<bool> = e.iter().map(|x| !x.excludes_zero()).collect();
        let mut next: Vec<Zone> = Vec::new();
        let mut i = 0;
        while i < p.len() {
            if straddles[i] {
                let start = i;
                while i < p.len() && straddles[i] {
                    i += 1;
                }
                let run = &p[start..i];
                if run.iter().any(|q| q.end == Some(End::Inner)) {
                    return Err(RootRefusal::Unresolved);
                }
                let cluster = Zone::cluster(run);
                if start > 0 && i < p.len() {
                    let change = positive(e[start - 1]) != positive(e[i]);
                    if change != (cluster.m % 2 == 1) {
                        return Err(RootRefusal::Unresolved);
                    }
                }
                next.push(cluster);
            } else {
                if i + 1 < p.len() && !straddles[i + 1] && positive(e[i]) != positive(e[i + 1]) {
                    next.push(bisect(tower, p[i].hi, p[i + 1].lo, j, positive(e[i]))?);
                }
                i += 1;
            }
        }
        if next.iter().map(|q| q.m).sum::<usize>() > k - j {
            return Err(RootRefusal::Unresolved);
        }
        z = next;
    }
    Ok(z)
}

/// The simple zero of `H_j` between `lo` and `hi`, where it changes sign,
/// narrowed as far as its enclosures stay signed: at a point whose
/// enclosure straddles zero (the zero itself, or rounding) the next
/// candidate is tried, so the bracket never stops at a wide stretch.
fn bisect<T: Tower>(
    tower: &T,
    mut lo: Scalar,
    mut hi: Scalar,
    j: usize,
    positive_at_lo: bool,
) -> Result<Zone, RootRefusal> {
    for _ in 0..100 {
        let mut narrowed = false;
        for c in SPLITS {
            let m = lo + (hi - lo) * c;
            if m <= lo || m >= hi {
                continue;
            }
            let e = tower.enclose(m, m, j).ok_or(RootRefusal::Unresolved)?[j];
            if !e.excludes_zero() {
                continue;
            }
            if positive(e) == positive_at_lo {
                lo = m;
            } else {
                hi = m;
            }
            narrowed = true;
            break;
        }
        if !narrowed {
            break;
        }
    }
    Ok(Zone {
        lo,
        hi,
        at: 0.5 * (lo + hi),
        m: 1,
        end: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A polynomial, `H_j` its `j`-th derivative, evaluated by interval
    /// Horner.
    struct PolyTower {
        /// Coefficients, constant term first.
        coeffs: Vec<Scalar>,
        joins: Vec<Scalar>,
    }

    impl PolyTower {
        fn new(coeffs: Vec<Scalar>) -> Self {
            Self {
                coeffs,
                joins: Vec::new(),
            }
        }
    }

    fn mul(a: &[Scalar], b: &[Scalar]) -> Vec<Scalar> {
        let mut out = vec![0.0; a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        out
    }

    fn monomial(degree: usize) -> Vec<Scalar> {
        let mut c = vec![0.0; degree + 1];
        c[degree] = 1.0;
        c
    }

    impl Tower for PolyTower {
        fn enclose(&self, t0: Scalar, t1: Scalar, order: usize) -> Option<Vec<I>> {
            let x = I::new(t0, t1);
            let mut c = self.coeffs.clone();
            let mut out = Vec::new();
            for _ in 0..=order {
                let mut acc = I::point(0.0);
                for &k in c.iter().rev() {
                    acc = acc.mul(x).add(I::point(k));
                }
                out.push(acc);
                c = c
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(i, k)| k * i as Scalar)
                    .collect();
                if c.is_empty() {
                    c = vec![0.0];
                }
            }
            Some(out)
        }

        fn joins(&self, t0: Scalar, t1: Scalar) -> Vec<Scalar> {
            self.joins
                .iter()
                .copied()
                .filter(|&j| j > t0 && j < t1)
                .collect()
        }
    }

    fn run(tower: &PolyTower, t0: Scalar, t1: Scalar) -> Result<Vec<(Scalar, usize)>, RootRefusal> {
        let mut budget = Budget::new(ROOT_BUDGET);
        isolate(tower, t0, t1, [End::Inner, End::Inner], &mut budget)
            .map(|z| z.into_iter().map(|z| (z.at, z.m)).collect())
    }

    #[test]
    fn a_fourth_order_contact_has_multiplicity_four() {
        let roots = run(&PolyTower::new(monomial(4)), -1.0, 1.0).unwrap();
        assert_eq!(roots.len(), 1);
        assert_eq!(roots[0].1, 4);
        assert!(roots[0].0.abs() < 1e-12);
    }

    #[test]
    fn contacts_up_to_the_cap_are_certified_and_beyond_refused() {
        for m in [6, 8] {
            let roots = run(&PolyTower::new(monomial(m)), -0.7, 1.3).unwrap();
            assert_eq!(roots.len(), 1, "t^{m}");
            assert_eq!(roots[0].1, m);
            assert!(roots[0].0.abs() < 1e-9);
        }
        for m in [9, 10] {
            assert_eq!(
                run(&PolyTower::new(monomial(m)), -0.7, 1.3),
                Err(RootRefusal::Unresolved),
                "t^{m}"
            );
        }
    }

    #[test]
    fn roots_of_different_orders_keep_their_multiplicity() {
        // (t - 1/4)^2 (t + 1/2)^3
        let mut p = vec![1.0];
        for _ in 0..2 {
            p = mul(&p, &[-0.25, 1.0]);
        }
        for _ in 0..3 {
            p = mul(&p, &[0.5, 1.0]);
        }
        let roots = run(&PolyTower::new(p), -1.0, 1.0).unwrap();
        assert_eq!(roots.len(), 2);
        assert!((roots[0].0 + 0.5).abs() < 1e-3 && roots[0].1 == 3);
        assert!((roots[1].0 - 0.25).abs() < 1e-3 && roots[1].1 == 2);
    }

    #[test]
    fn a_contact_breaks_into_two_roots_or_none_by_its_sign() {
        let mut c = monomial(4);
        c[0] = 1e-12;
        assert_eq!(run(&PolyTower::new(c.clone()), -1.0, 1.0).unwrap(), vec![]);
        c[0] = -1e-12;
        let roots = run(&PolyTower::new(c), -1.0, 1.0).unwrap();
        assert_eq!(roots.len(), 2);
        for (t, m) in roots {
            assert_eq!(m, 1);
            assert!((t.abs() - 1e-3).abs() < 1e-9, "{t}");
        }
    }

    #[test]
    fn a_curve_end_reports_the_contact_and_an_inner_end_refuses() {
        let tower = PolyTower::new(monomial(4));
        let mut budget = Budget::new(ROOT_BUDGET);
        let zones = isolate(&tower, 0.0, 1.0, [End::Curve, End::Curve], &mut budget).unwrap();
        assert_eq!(zones.len(), 1);
        assert_eq!((zones[0].at, zones[0].m), (0.0, 4));
        let mut budget = Budget::new(ROOT_BUDGET);
        assert_eq!(
            isolate(&tower, 0.0, 1.0, [End::Inner, End::Curve], &mut budget).unwrap_err(),
            RootRefusal::Unresolved
        );
    }

    fn seven_simple() -> PolyTower {
        let mut p = monomial(1);
        for r in [0.3, 0.6, 0.9] {
            p = mul(&p, &[-r * r, 0.0, 1.0]);
        }
        for _ in 0..2 {
            p = mul(&p, &[1.0, 0.0, 1.0]);
        }
        PolyTower::new(p)
    }

    #[test]
    fn simple_roots_are_all_found_including_one_at_the_middle() {
        let roots = run(&seven_simple(), -1.0, 1.0).unwrap();
        let expected = [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9];
        assert_eq!(roots.len(), 7);
        for ((t, m), x) in roots.iter().zip(expected) {
            assert_eq!(*m, 1);
            assert!((t - x).abs() < 1e-9, "{t} vs {x}");
        }
    }

    #[test]
    fn the_budget_refuses_by_name() {
        let tower = seven_simple();
        let mut budget = Budget::new(2);
        assert_eq!(
            isolate(&tower, -1.0, 1.0, [End::Inner, End::Inner], &mut budget).unwrap_err(),
            RootRefusal::Budget
        );
    }

    #[test]
    fn the_chain_alone_reads_the_multiplicity_from_a_higher_order() {
        // t^4 (1 + t^2): H_6 never vanishes, and t = 0 is a root of order 4.
        let tower = PolyTower::new(mul(&monomial(4), &[1.0, 0.0, 1.0]));
        let zones = chain(&tower, -0.6, 0.9, 6, [End::Inner, End::Inner]).unwrap();
        assert_eq!(zones.len(), 1);
        assert_eq!(zones[0].m, 4);
    }

    #[test]
    fn a_join_at_a_root_is_not_used_to_split() {
        // t^2 (t - 3/4), with a preferred cut at the double root.
        let mut tower = PolyTower::new(mul(&monomial(2), &[-0.75, 1.0]));
        tower.joins = vec![0.0];
        let roots = run(&tower, -1.0, 1.0).unwrap();
        assert_eq!(roots.len(), 2);
        assert_eq!(roots[0].1, 2);
        assert!(roots[0].0.abs() < 1e-9);
        assert_eq!(roots[1].1, 1);
        assert!((roots[1].0 - 0.75).abs() < 1e-9);
    }

    /// A tower whose enclosures contradict each other: `H_1 = 2t` has one
    /// zero, at 0, yet `H_0` changes sign across it and straddles zero only
    /// in the middle. A zero of `H_1` makes it a double zero of `H_0` or
    /// none, so the chain cannot be right and must say so.
    struct Inconsistent;

    impl Tower for Inconsistent {
        fn enclose(&self, t0: Scalar, t1: Scalar, order: usize) -> Option<Vec<I>> {
            let h0 = if t1 < -0.25 {
                I::new(-2.0, -1.0)
            } else if t0 > 0.25 {
                I::new(1.0, 2.0)
            } else {
                I::new(-1.0, 1.0)
            };
            let all = [h0, I::new(2.0 * t0, 2.0 * t1), I::point(2.0)];
            Some(all[..=order.min(2)].to_vec())
        }

        fn joins(&self, _t0: Scalar, _t1: Scalar) -> Vec<Scalar> {
            Vec::new()
        }
    }

    #[test]
    fn a_zone_whose_parity_contradicts_its_neighbours_refuses() {
        let mut budget = Budget::new(ROOT_BUDGET);
        assert_eq!(
            isolate(
                &Inconsistent,
                -1.0,
                1.0,
                [End::Inner, End::Inner],
                &mut budget
            )
            .unwrap_err(),
            RootRefusal::Unresolved
        );
    }
}
