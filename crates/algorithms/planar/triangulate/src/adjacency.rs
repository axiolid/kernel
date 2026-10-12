// SPDX-License-Identifier: MPL-2.0

//! Twin pointers: the full computation, the local link, and the guard
//! that holds the two equal.
//!
//! Insertion and flipping keep the twin array current by re-linking only
//! the halfedges of the triangles they rewrite and of those triangles'
//! outside neighbours (#268). Recomputing it from the triangle list after
//! every insertion and flip, as the build did, cost O(n) per step and
//! O(n^2) or worse per triangulation.
//!
//! Local updates are easy to get subtly wrong: an operation that reuses a
//! triangle slot leaves the old neighbours of that slot pointing at
//! halfedges that now run between different vertices, and the corruption
//! only surfaces later as a point-location walk that leaves the mesh. So in
//! debug builds every operation ends with [`debug_check`], which compares
//! the maintained array with [`twins`] computed from scratch and panics on
//! the first difference. The comparison is O(n) per operation; outside
//! this crate's own unit tests it runs only up to [`FULL_CHECK_LIMIT`]
//! triangles, so downstream debug builds do not inherit the quadratic cost
//! the local updates removed.
//!
//! For a manifold triangulation, which every operation here preserves,
//! each halfedge has at most one reverse, so the full computation has a
//! single answer and "local equals full" means the build behaves exactly as
//! it did when it recomputed: same walks, same flips, same output.

use std::collections::HashMap;

use crate::mesh::NO_HALFEDGE;

/// Largest triangle count the debug check compares at in an ordinary debug
/// build; `None`, unlimited, in this crate's unit tests.
pub(crate) const FULL_CHECK_LIMIT: Option<usize> = if cfg!(test) { None } else { Some(1024) };

/// Every twin pointer, computed from the triangle list alone.
pub(crate) fn twins(triangles: &[u32]) -> Vec<u32> {
    pair(triangles).0
}

/// The twin pointers, and whether every directed edge occurs once.
///
/// When one occurs twice -- two overlapping triangles on the same edge, as
/// duplicate input points left before they were merged (#300) -- the
/// pairing is whichever the scan met last, and is not the only answer.
fn pair(triangles: &[u32]) -> (Vec<u32>, bool) {
    #[cfg(test)]
    counters::bump(&counters::FULL);
    let count = triangles.len();
    let mut halfedges = vec![NO_HALFEDGE; count];
    let mut seen: HashMap<(u32, u32), u32> = HashMap::with_capacity(count);
    let mut manifold = true;
    for e in 0..count {
        let t = e / 3;
        let i = e % 3;
        let from = triangles[3 * t + i];
        let to = triangles[3 * t + (i + 1) % 3];
        if let Some(&twin) = seen.get(&(to, from)) {
            halfedges[e] = twin;
            halfedges[twin as usize] = e as u32;
        } else if seen.insert((from, to), e as u32).is_some() {
            manifold = false;
        }
    }
    (halfedges, manifold)
}

/// Make `a` and `b` twins. `b` may be [`NO_HALFEDGE`], for a halfedge on
/// the boundary.
pub(crate) fn link(halfedges: &mut [u32], a: usize, b: u32) {
    halfedges[a] = b;
    if b != NO_HALFEDGE {
        halfedges[b as usize] = a as u32;
    }
}

/// Note one local adjacency update, for the tests that count them.
#[inline]
pub(crate) fn note_update() {
    #[cfg(test)]
    counters::bump(&counters::UPDATES);
}

/// How a maintained twin array compares with a full recomputation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Comparison {
    /// Every twin pointer agrees.
    Equal,
    /// The first halfedge whose twin differs.
    Differs(usize),
    /// A directed edge occurs twice, so the recomputation has no single
    /// answer to compare with.
    Ambiguous,
}

/// Compare the maintained twin array with a full recomputation.
pub(crate) fn compare(triangles: &[u32], halfedges: &[u32]) -> Comparison {
    let (full, manifold) = pair(triangles);
    if !manifold {
        return Comparison::Ambiguous;
    }
    if full.len() != halfedges.len() {
        return Comparison::Differs(full.len().min(halfedges.len()));
    }
    match full.iter().zip(halfedges).position(|(f, h)| f != h) {
        Some(e) => Comparison::Differs(e),
        None => Comparison::Equal,
    }
}

/// In debug builds, panic unless the maintained twin array equals a full
/// recomputation (up to [`FULL_CHECK_LIMIT`] triangles).
///
/// A triangle list with a repeated directed edge panics too. Coincident
/// input points are merged before they reach the mesh (#300), so no input
/// produces one, and a build that does is broken; this check used to skip
/// such lists, because duplicate points made them, and had nothing exact
/// to compare with.
#[inline]
pub(crate) fn debug_check(triangles: &[u32], halfedges: &[u32]) {
    let too_large = FULL_CHECK_LIMIT.is_some_and(|limit| triangles.len() / 3 > limit);
    if !cfg!(debug_assertions) || too_large {
        return;
    }
    #[cfg(test)]
    counters::bump(&counters::CHECKS);
    match compare(triangles, halfedges) {
        Comparison::Equal => {}
        Comparison::Ambiguous => panic!(
            "triangulation has a directed edge that occurs twice, \
             so its adjacency has no single rebuild"
        ),
        Comparison::Differs(e) => panic!(
            "triangulation adjacency diverged from a full rebuild at halfedge {e}: \
             maintained {:?}, rebuilt {:?}",
            halfedges.get(e),
            twins(triangles).get(e)
        ),
    }
}

/// Per-thread operation counters, so a test can show that adjacency is
/// maintained locally and checked after every update.
#[cfg(test)]
pub(crate) mod counters {
    use std::cell::Cell;
    use std::thread::LocalKey;

    thread_local! {
        /// Full twin computations, whether for maintenance or for checks.
        pub(crate) static FULL: Cell<usize> = const { Cell::new(0) };
        /// Local updates: splits and flips.
        pub(crate) static UPDATES: Cell<usize> = const { Cell::new(0) };
        /// Debug checks run.
        pub(crate) static CHECKS: Cell<usize> = const { Cell::new(0) };
    }

    pub(crate) fn bump(counter: &'static LocalKey<Cell<usize>>) {
        counter.with(|c| c.set(c.get() + 1));
    }

    /// Reset every counter on this thread.
    pub(crate) fn reset() {
        for counter in [&FULL, &UPDATES, &CHECKS] {
            counter.with(|c| c.set(0));
        }
    }

    /// `(full, updates, checks)` on this thread.
    pub(crate) fn read() -> (usize, usize, usize) {
        (
            FULL.with(Cell::get),
            UPDATES.with(Cell::get),
            CHECKS.with(Cell::get),
        )
    }
}

#[cfg(test)]
mod tests {
    use axiolid_core::Point2;

    use super::{compare, counters, debug_check, twins, Comparison};
    use crate::mesh::NO_HALFEDGE;
    use crate::{triangulate, Constraint};

    /// Seeded uniform points in `[0, side)^2`.
    fn random(count: usize, side: f64, seed: u64) -> Vec<Point2> {
        let mut s = seed;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            side * (s >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..count).map(|_| Point2::new(next(), next())).collect()
    }

    /// A closed polyline through `corners`, sampled every `spacing`, with
    /// its pieces as constraints numbered from `base`.
    fn outline(corners: &[(f64, f64)], spacing: f64, base: u32) -> (Vec<Point2>, Vec<Constraint>) {
        let mut points = Vec::new();
        for k in 0..corners.len() {
            let (a, b) = (corners[k], corners[(k + 1) % corners.len()]);
            let n = ((b.0 - a.0).hypot(b.1 - a.1) / spacing).ceil().max(1.0) as usize;
            for j in 0..n {
                let t = j as f64 / n as f64;
                points.push(Point2::new(a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t));
            }
        }
        let n = points.len() as u32;
        let constraints = (0..n)
            .map(|i| Constraint::new(base + i, base + (i + 1) % n))
            .collect();
        (points, constraints)
    }

    /// Triangulate, and return `(full, updates, checks)` for the call.
    fn counted(points: &[Point2], constraints: &[Constraint]) -> (usize, usize, usize) {
        counters::reset();
        triangulate(points, constraints).expect("triangulates");
        counters::read()
    }

    /// Every local update is followed by a check, and the only full
    /// computations are those checks and the one hand-over in `finish`:
    /// none maintains adjacency, so a release build, which skips the
    /// checks, computes it once per triangulation. The first triangle and
    /// its ghosts are one update for three points; every other point is at
    /// least one more.
    fn assert_local(name: &str, (full, updates, checks): (usize, usize, usize), points: usize) {
        assert!(
            updates + 2 >= points,
            "{name}: {updates} updates for {points} points"
        );
        if cfg!(debug_assertions) {
            assert_eq!(checks, updates, "{name}: an update went unchecked");
        }
        assert_eq!(
            full,
            checks + 1,
            "{name}: a full adjacency rebuild maintained the mesh"
        );
    }

    #[test]
    fn random_points_are_maintained_locally_and_checked_after_every_update() {
        for (count, seed) in [(3, 1), (50, 2), (400, 3), (1500, 4)] {
            let points = random(count, 10.0, 0x2680 + seed);
            assert_local("random", counted(&points, &[]), count);
        }
    }

    #[test]
    fn points_on_edges_are_split_locally() {
        // A grid puts most points exactly on an existing edge, so the
        // two-into-four split runs as often as the one-into-three.
        for side in [3usize, 7, 20] {
            let points: Vec<Point2> = (0..side * side)
                .map(|k| Point2::new((k % side) as f64, (k / side) as f64))
                .collect();
            assert_local("grid", counted(&points, &[]), points.len());
        }
    }

    #[test]
    fn recovery_and_restoring_flips_are_maintained_locally() {
        // An outer outline, a hole, and chords through a point cloud:
        // recovery flips its way through the cloud, and the restoring
        // pass flips the long triangles back.
        let (mut points, mut constraints) = outline(
            &[(-1.0, -1.0), (11.0, -1.0), (11.0, 11.0), (-1.0, 11.0)],
            0.5,
            0,
        );
        let base = points.len() as u32;
        let (hole, around) = outline(&[(4.0, 4.0), (4.0, 6.0), (6.5, 6.0), (6.0, 4.0)], 0.3, base);
        points.extend(hole);
        constraints.extend(around);
        for q in random(800, 10.0, 0x2681) {
            if !(q.x > 3.9 && q.x < 6.6 && q.y > 3.9 && q.y < 6.1) {
                points.push(q);
            }
        }
        for k in 0..4 {
            let a = points.len() as u32;
            points.push(Point2::new(-0.5, 0.7 + 0.8 * k as f64));
            points.push(Point2::new(3.5, 0.5 + 0.8 * k as f64));
            constraints.push(Constraint::new(a, a + 1));
        }
        let before = counted(&points, &[]).1;
        let with = counted(&points, &constraints);
        assert!(with.1 > before, "the constraints needed no flip");
        assert_local("constrained", with, points.len());
    }

    #[test]
    fn the_check_catches_a_stale_twin() {
        let points = random(30, 1.0, 0x2682);
        let tri = triangulate(&points, &[]).expect("triangulates");
        let mut halfedges = twins(&tri.triangles);
        assert_eq!(compare(&tri.triangles, &halfedges), Comparison::Equal);
        // A twin left pointing at the halfedge that used to hold its slot,
        // as a reused triangle slot leaves its old neighbours.
        let e = halfedges.iter().position(|&h| h != NO_HALFEDGE).unwrap();
        let h = halfedges[e] as usize;
        halfedges[e] = (3 * (h / 3) + (h % 3 + 1) % 3) as u32;
        assert_eq!(compare(&tri.triangles, &halfedges), Comparison::Differs(e));
    }

    #[test]
    fn duplicate_points_are_merged_and_checked_after_every_update() {
        // Before #300 a repeated point left overlapping triangles on one
        // directed edge, and the check skipped them. Merged, the mesh stays
        // manifold and every update is checked: a skip would now panic.
        let mut points: Vec<Point2> = (0..9)
            .map(|k| Point2::new(f64::from(k % 3), f64::from(k / 3)))
            .collect();
        points.extend_from_slice(&points.clone()[..3]);
        points.push(points[4]);
        points.extend(random(200, 2.0, 0x3000));
        points.extend_from_within(12..60);
        let (full, updates, checks) = counted(&points, &[]);
        if cfg!(debug_assertions) {
            assert_eq!(checks, updates);
        }
        assert_eq!(full, checks + 1);
        let tri = triangulate(&points, &[]).expect("triangulates");
        assert_eq!(compare(&tri.triangles, &tri.halfedges), Comparison::Equal);
        assert_eq!(tri.merged_points().len(), 3 + 1 + 48);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "occurs twice")]
    fn the_debug_check_panics_on_a_repeated_directed_edge() {
        assert_eq!(
            compare(&[0, 1, 2, 0, 1, 3], &[NO_HALFEDGE; 6]),
            Comparison::Ambiguous
        );
        debug_check(&[0, 1, 2, 0, 1, 3], &[NO_HALFEDGE; 6]);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "diverged from a full rebuild")]
    fn the_debug_check_panics_on_a_stale_twin() {
        let points = random(30, 1.0, 0x2683);
        let tri = triangulate(&points, &[]).expect("triangulates");
        let mut halfedges = tri.halfedges.clone();
        halfedges.swap(0, 1);
        debug_check(&tri.triangles, &halfedges);
    }
}
