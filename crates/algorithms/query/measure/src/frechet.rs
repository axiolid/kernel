//! Fréchet distance between polylines (#147).
//!
//! The Fréchet distance is the least leash length that lets two walkers
//! traverse the two curves from start to end, each moving only forward. The
//! discrete version walks vertex to vertex; the continuous version walks
//! along the segments too, and is never larger.
//!
//! - [`discrete_frechet_distance`]: the Eiter-Mannila dynamic programme,
//!   `O(nm)` time and `O(m)` memory, iterative.
//! - [`frechet_at_most`]: the Alt-Godau decision, `O(nm)`: can the curves be
//!   walked with a leash of length `eps`?
//! - [`frechet_distance`]: the continuous distance. The answer is always one
//!   of the critical values of the free space -- a distance between two
//!   vertices, from a vertex to a segment, or from a point on a segment to
//!   two vertices of the other curve at once -- so it is the least critical
//!   value the decision accepts, found by binary search over the sorted
//!   critical values. When there are more than [`CANDIDATE_BUDGET`] of them
//!   (there are `O(n^2 m + n m^2)`), the range is first narrowed by bisection
//!   on the decision, so memory stays bounded.
//!
//! The `_2d` functions lift the points to `z = 0`, which changes no distance.
//!
//! # Rounding
//!
//! Everything is `f64`. Critical values are computed in floating point, and
//! the decision reads its free-space intervals with the same helpers, so at
//! a vertex-to-vertex or vertex-to-segment critical value the decision sees
//! exactly the distance the candidate was built from. A passage that opens
//! at a point equidistant to two vertices is decided by comparing interval
//! ends, which carry a few units of rounding; there the answer can be the
//! next larger critical value. Nothing here is certified.

use core::fmt;

use axiolid_core::{Point2, Point3};

/// Critical values held at once before the range is narrowed by bisection.
pub const CANDIDATE_BUDGET: usize = 1 << 20;

/// Why a Fréchet query was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrechetError {
    /// A polyline has no points.
    EmptyPolyline,
    /// A coordinate is NaN or infinite.
    NonFiniteInput,
    /// The leash length for [`frechet_at_most`] is negative or not finite.
    InvalidBound,
}

impl fmt::Display for FrechetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyPolyline => formatter.write_str("a polyline has no points"),
            Self::NonFiniteInput => formatter.write_str("polyline coordinates must be finite"),
            Self::InvalidBound => {
                formatter.write_str("the leash length must be finite and non-negative")
            }
        }
    }
}

impl std::error::Error for FrechetError {}

/// Discrete Fréchet distance between the vertex sequences of `a` and `b`.
///
/// # Errors
///
/// An empty polyline or a non-finite coordinate.
pub fn discrete_frechet_distance(a: &[Point3], b: &[Point3]) -> Result<f64, FrechetError> {
    check(a)?;
    check(b)?;
    Ok(discrete(a, b))
}

/// [`discrete_frechet_distance`] in the plane.
///
/// # Errors
///
/// As [`discrete_frechet_distance`].
pub fn discrete_frechet_distance_2d(a: &[Point2], b: &[Point2]) -> Result<f64, FrechetError> {
    discrete_frechet_distance(&lift(a), &lift(b))
}

/// Whether the continuous Fréchet distance between `a` and `b` is at most
/// `eps`.
///
/// # Errors
///
/// An empty polyline, a non-finite coordinate, or a negative or non-finite
/// `eps`.
pub fn frechet_at_most(a: &[Point3], b: &[Point3], eps: f64) -> Result<bool, FrechetError> {
    check(a)?;
    check(b)?;
    if !eps.is_finite() || eps < 0.0 {
        return Err(FrechetError::InvalidBound);
    }
    Ok(decide(a, b, eps))
}

/// [`frechet_at_most`] in the plane.
///
/// # Errors
///
/// As [`frechet_at_most`].
pub fn frechet_at_most_2d(a: &[Point2], b: &[Point2], eps: f64) -> Result<bool, FrechetError> {
    frechet_at_most(&lift(a), &lift(b), eps)
}

/// Continuous Fréchet distance between the polylines `a` and `b`.
///
/// # Errors
///
/// An empty polyline or a non-finite coordinate.
pub fn frechet_distance(a: &[Point3], b: &[Point3]) -> Result<f64, FrechetError> {
    check(a)?;
    check(b)?;
    Ok(continuous(a, b, CANDIDATE_BUDGET))
}

/// [`frechet_distance`] in the plane.
///
/// # Errors
///
/// As [`frechet_distance`].
pub fn frechet_distance_2d(a: &[Point2], b: &[Point2]) -> Result<f64, FrechetError> {
    frechet_distance(&lift(a), &lift(b))
}

pub(crate) fn check(points: &[Point3]) -> Result<(), FrechetError> {
    if points.is_empty() {
        return Err(FrechetError::EmptyPolyline);
    }
    if points.iter().any(|p| !p.is_finite()) {
        return Err(FrechetError::NonFiniteInput);
    }
    Ok(())
}

pub(crate) fn lift(points: &[Point2]) -> Vec<Point3> {
    points.iter().map(|p| Point3::new(p.x, p.y, 0.0)).collect()
}

pub(crate) fn distance(p: Point3, q: Point3) -> f64 {
    (p - q).length()
}

fn discrete(a: &[Point3], b: &[Point3]) -> f64 {
    // `row[j]` is the coupling cost of a[..=i] with b[..=j].
    let mut row = vec![0.0_f64; b.len()];
    for (i, &p) in a.iter().enumerate() {
        let mut diagonal = 0.0;
        for (j, &q) in b.iter().enumerate() {
            let here = distance(p, q);
            let above = row[j];
            let best = match (i, j) {
                (0, 0) => 0.0,
                (0, _) => row[j - 1],
                (_, 0) => above,
                _ => above.min(diagonal).min(row[j - 1]),
            };
            diagonal = above;
            row[j] = here.max(best);
        }
    }
    row[b.len() - 1]
}

/// The foot of `c` on the line through segment `(a, b)`: its parameter
/// (unclamped) and its distance to `c`. `None` for a zero-length segment.
pub(crate) fn foot(c: Point3, a: Point3, b: Point3) -> Option<(f64, f64)> {
    let d = b - a;
    let length_squared = d.length_squared();
    if length_squared == 0.0 {
        return None;
    }
    let t = (c - a).dot(d) / length_squared;
    Some((t, distance(a + d * t, c)))
}

/// A closed parameter interval on a segment, empty when `lo > hi`.
#[derive(Clone, Copy)]
struct Span {
    lo: f64,
    hi: f64,
}

impl Span {
    const EMPTY: Span = Span {
        lo: f64::INFINITY,
        hi: f64::NEG_INFINITY,
    };

    fn is_empty(self) -> bool {
        self.lo > self.hi
    }
}

/// The parameters of segment `(a, b)` within `eps` of `c`.
///
/// Endpoints are decided by the same vertex distance the vertex-to-vertex
/// critical values are built from, and the interior by the same foot
/// distance the vertex-to-segment ones are, so at a critical value the
/// interval opens exactly.
fn free(c: Point3, a: Point3, b: Point3, eps: f64) -> Span {
    let start = distance(a, c) <= eps;
    let end = distance(b, c) <= eps;
    let Some((t, h)) = foot(c, a, b) else {
        return if start {
            Span { lo: 0.0, hi: 1.0 }
        } else {
            Span::EMPTY
        };
    };
    if h > eps {
        return Span::EMPTY;
    }
    let w = (eps * eps - h * h).max(0.0).sqrt() / (b - a).length();
    let lo = if start { 0.0 } else { (t - w).max(0.0) };
    let hi = if end { 1.0 } else { (t + w).min(1.0) };
    if lo > hi {
        Span::EMPTY
    } else {
        Span { lo, hi }
    }
}

/// Alt-Godau: is the top-right corner of the free space reachable from the
/// bottom-left by a monotone path?
pub(crate) fn decide(a: &[Point3], b: &[Point3], eps: f64) -> bool {
    if distance(a[0], b[0]) > eps || distance(a[a.len() - 1], b[b.len() - 1]) > eps {
        return false;
    }
    let (p, q) = (a.len() - 1, b.len() - 1);
    if p == 0 || q == 0 {
        // One curve is a point: every point of the other must be within
        // `eps` of it, and a segment's farthest point is an endpoint.
        let (point, other) = if p == 0 { (a[0], b) } else { (b[0], a) };
        return other.iter().all(|&v| distance(point, v) <= eps);
    }
    // `left[j]`: reachable part of the left edge of cell (i, j), a span of
    // b's segment j at a's vertex i. `bottom`: reachable part of the bottom
    // edge of cell (i, j), a span of a's segment i at b's vertex j.
    let mut left: Vec<Span> = Vec::with_capacity(q);
    let mut open = true;
    for j in 0..q {
        let span = free(a[0], b[j], b[j + 1], eps);
        left.push(if open && span.lo == 0.0 {
            span
        } else {
            Span::EMPTY
        });
        open = open && !span.is_empty() && span.lo == 0.0 && span.hi == 1.0;
    }
    let mut bottom_open = true;
    for i in 0..p {
        // Bottom edge of cell (i, 0): along the start of b, as for `left`.
        let span = free(b[0], a[i], a[i + 1], eps);
        let mut bottom = if bottom_open && span.lo == 0.0 {
            span
        } else {
            Span::EMPTY
        };
        bottom_open = bottom_open && !span.is_empty() && span.lo == 0.0 && span.hi == 1.0;
        for j in 0..q {
            let l = left[j];
            let top = free(b[j + 1], a[i], a[i + 1], eps);
            let right = free(a[i + 1], b[j], b[j + 1], eps);
            let next_bottom = if !l.is_empty() {
                top
            } else if !bottom.is_empty() {
                Span {
                    lo: bottom.lo.max(top.lo),
                    hi: top.hi,
                }
            } else {
                Span::EMPTY
            };
            left[j] = if !bottom.is_empty() {
                right
            } else if !l.is_empty() {
                Span {
                    lo: l.lo.max(right.lo),
                    hi: right.hi,
                }
            } else {
                Span::EMPTY
            };
            bottom = next_bottom;
        }
        // `bottom` is now the top edge of cell (i, q - 1).
        if i == p - 1 {
            return (!bottom.is_empty() && bottom.hi == 1.0)
                || (!left[q - 1].is_empty() && left[q - 1].hi == 1.0);
        }
    }
    unreachable!("p >= 1 returns in the loop")
}

/// Every critical value: vertex-vertex, vertex-segment, and the points of a
/// segment equidistant to two vertices of the other curve.
fn for_each_candidate(a: &[Point3], b: &[Point3], mut visit: impl FnMut(f64)) {
    for &p in a {
        for &q in b {
            visit(distance(p, q));
        }
    }
    for (points, other) in [(a, b), (b, a)] {
        for &c in points {
            for s in other.windows(2) {
                if let Some((t, h)) = foot(c, s[0], s[1]) {
                    if (0.0..=1.0).contains(&t) {
                        visit(h);
                    }
                }
            }
        }
        for s in other.windows(2) {
            let (start, d) = (s[0], s[1] - s[0]);
            for k in 0..points.len() {
                for l in k + 1..points.len() {
                    let n = points[l] - points[k];
                    let along = n.dot(d);
                    if along == 0.0 {
                        continue;
                    }
                    let middle = (points[k] + points[l]) * 0.5;
                    let t = n.dot(middle - start) / along;
                    if (0.0..=1.0).contains(&t) {
                        let x = start + d * t;
                        visit(distance(x, points[k]).max(distance(x, points[l])));
                    }
                }
            }
        }
    }
}

fn continuous(a: &[Point3], b: &[Point3], budget: usize) -> f64 {
    // The endpoints must be coupled, so their distances bound from below.
    let mut lo = distance(a[0], b[0]).max(distance(a[a.len() - 1], b[b.len() - 1]));
    if decide(a, b, lo) {
        return lo;
    }
    // The discrete coupling, interpolated, is a continuous one.
    let mut hi = discrete(a, b);
    // Invariant: the answer lies in (lo, hi] and `decide(hi)` holds.
    loop {
        let mut count = 0usize;
        for_each_candidate(a, b, |c| {
            if c > lo && c <= hi {
                count += 1;
            }
        });
        let mid = lo + (hi - lo) * 0.5;
        if count <= budget || mid <= lo || mid >= hi {
            break;
        }
        if decide(a, b, mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let mut candidates = Vec::new();
    for_each_candidate(a, b, |c| {
        if c > lo && c <= hi {
            candidates.push(c);
        }
    });
    candidates.sort_by(f64::total_cmp);
    candidates.dedup();
    // Least accepted candidate; the decision is monotone in `eps`.
    let first = candidates.partition_point(|&c| !decide(a, b, c));
    candidates.get(first).copied().unwrap_or(hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic pseudo-random walk of `n` points in space.
    fn walk(seed: u64, n: usize) -> Vec<Point3> {
        let mut state = seed;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        };
        let mut p = Point3::ZERO;
        (0..n)
            .map(|_| {
                p += Point3::new(1.0 + next(), next(), next());
                p
            })
            .collect()
    }

    /// Narrowing the range by bisection before listing candidates finds
    /// the same critical value as listing them all at once.
    #[test]
    fn a_tiny_budget_bisects_to_the_same_answer() {
        for seed in 0..16 {
            let a = walk(seed, 12);
            let jitter = walk(seed + 50, 12);
            let b: Vec<Point3> = (0..12)
                .map(|i| {
                    if i == 0 || i == 11 {
                        a[i]
                    } else {
                        a[i] + (jitter[i] - jitter[i - 1]) * 0.4
                    }
                })
                .collect();
            let listed = continuous(&a, &b, usize::MAX);
            let mut bisections = 0;
            let mut lo = 0.0_f64;
            let mut hi = discrete(&a, &b);
            // The same narrowing as `continuous`, counted, to show it runs.
            loop {
                let mut count = 0usize;
                for_each_candidate(&a, &b, |c| {
                    if c > lo && c <= hi {
                        count += 1;
                    }
                });
                let mid = lo + (hi - lo) * 0.5;
                if count <= 2 || mid <= lo || mid >= hi {
                    break;
                }
                bisections += 1;
                if decide(&a, &b, mid) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            assert!(bisections > 0, "{seed}");
            assert_eq!(continuous(&a, &b, 2), listed, "{seed}");
            assert!(listed > 0.0);
        }
    }
}
