// SPDX-License-Identifier: MPL-2.0

//! Helpers shared by the integration tests: a seeded point stream, the
//! exact convex hull, and the check that a triangulation covers it.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{incircle, orient2d};
use axiolid_triangulate::{Constraint, Triangulation};

/// A point's position as a hashable key; `+ 0.0` reads -0.0 as 0.0.
pub fn key(p: Point2) -> (u64, u64) {
    ((p.x + 0.0).to_bits(), (p.y + 0.0).to_bits())
}

pub fn sign(c: Certified) -> Sign {
    match c {
        Certified::Certain { sign, .. } => sign,
        other => panic!("undecided predicate: {other:?}"),
    }
}

/// A seeded xorshift stream of uniform values in `[0, 1)`.
pub struct Rng(pub u64);

impl Rng {
    pub fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

pub fn uniform(count: usize, seed: u64) -> Vec<Point2> {
    let mut r = Rng(seed);
    (0..count)
        .map(|_| Point2::new(r.unit(), r.unit()))
        .collect()
}

/// Every input point on the boundary of the exact convex hull, collinear
/// ones included.
pub fn hull_boundary(points: &[Point2]) -> HashSet<(u64, u64)> {
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    sorted.dedup();
    let mut corners: Vec<Point2> = Vec::new();
    for pass in 0..2 {
        let start = corners.len();
        let iter: Box<dyn Iterator<Item = &Point2>> = if pass == 0 {
            Box::new(sorted.iter())
        } else {
            Box::new(sorted.iter().rev())
        };
        for &p in iter {
            while corners.len() >= start + 2
                && sign(orient2d(
                    corners[corners.len() - 2],
                    corners[corners.len() - 1],
                    p,
                )) != Sign::Positive
            {
                corners.pop();
            }
            corners.push(p);
        }
        corners.pop();
    }
    let within = |a: f64, b: f64, x: f64| a.min(b) <= x && x <= a.max(b);
    let mut on = HashSet::new();
    for k in 0..corners.len() {
        let (a, b) = (corners[k], corners[(k + 1) % corners.len()]);
        for &p in &sorted {
            if sign(orient2d(a, b, p)) == Sign::Zero
                && within(a.x, b.x, p.x)
                && within(a.y, b.y, p.y)
            {
                on.insert(key(p));
            }
        }
    }
    on
}

/// Fail unless `tri` is a triangulation of the distinct points of
/// `points` whose boundary is their exact convex hull, Delaunay away from
/// its constraints.
pub fn assert_covers_hull(name: &str, points: &[Point2], tri: &Triangulation) {
    let at = |k: u32| points[k as usize];
    let mut edges: HashMap<(u32, u32), u32> = HashMap::new();
    let mut used = vec![false; points.len()];
    for t in tri.triangles().chunks_exact(3) {
        assert_eq!(
            sign(orient2d(at(t[0]), at(t[1]), at(t[2]))),
            Sign::Positive,
            "{name}: triangle {t:?} is not counter-clockwise"
        );
        for i in 0..3 {
            used[t[i] as usize] = true;
            let before = edges.insert((t[i], t[(i + 1) % 3]), t[(i + 2) % 3]);
            assert!(before.is_none(), "{name}: a directed edge occurs twice");
        }
    }
    let distinct: HashSet<(u64, u64)> = points.iter().map(|p| key(*p)).collect();
    let corners: HashSet<(u64, u64)> = (0..points.len())
        .filter(|&k| used[k])
        .map(|k| key(points[k]))
        .collect();
    assert_eq!(corners, distinct, "{name}: a point is not a corner");

    let mut boundary = HashSet::new();
    let mut out = HashMap::new();
    for (&(a, b), &apex) in &edges {
        match edges.get(&(b, a)) {
            Some(&other) => assert!(
                tri.constraints().contains(&Constraint::new(a, b))
                    || sign(incircle(at(a), at(b), at(apex), at(other))) != Sign::Positive,
                "{name}: edge ({a}, {b}) is not locally Delaunay"
            ),
            None => {
                assert!(out.insert(a, b).is_none(), "{name}: boundary pinches");
                boundary.insert(key(at(a)));
                // A hull edge has no point strictly outside it.
                assert!(
                    points
                        .iter()
                        .all(|&p| sign(orient2d(at(a), at(b), p)) != Sign::Negative),
                    "{name}: boundary edge ({a}, {b}) is not on the hull"
                );
            }
        }
    }
    // One closed walk through every boundary vertex.
    let first = *out.keys().next().expect("a boundary");
    let (mut v, mut steps) = (first, 0);
    loop {
        v = out[&v];
        steps += 1;
        if v == first {
            break;
        }
    }
    assert_eq!(steps, out.len(), "{name}: the boundary is not one loop");

    let hull = hull_boundary(points);
    assert_eq!(boundary, hull, "{name}: boundary is not the exact hull");
    let (n, h) = (distinct.len(), hull.len());
    assert_eq!(
        tri.triangle_count(),
        2 * n - 2 - h,
        "{name}: {n} points with {h} on the hull"
    );
}
