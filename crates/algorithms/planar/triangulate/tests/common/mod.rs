// SPDX-License-Identifier: MPL-2.0

//! Helpers shared by the integration tests: a seeded point stream, the
//! exact convex hull, and the check that a triangulation covers it.

#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{incircle, orient2d};
use axiolid_triangulate::{Constraint, Ring, RingKind, Triangulation};

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
    for t in tri.triangles().as_chunks::<3>().0 {
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

/// Every unconstrained edge between two triangles is locally Delaunay: the
/// vertex across it lies outside (or on) the circumcircle, exactly.
pub fn assert_constrained_delaunay(tri: &Triangulation) {
    let t = tri.triangles();
    let pts = tri.points();
    let constrained: std::collections::BTreeSet<Constraint> =
        tri.constraints().iter().copied().collect();
    let mut across: HashMap<(u32, u32), u32> = HashMap::new();
    for k in 0..tri.triangle_count() {
        for i in 0..3 {
            across.insert(
                (t[3 * k + i], t[3 * k + (i + 1) % 3]),
                t[3 * k + (i + 2) % 3],
            );
        }
    }
    for (&(u, v), &apex) in &across {
        if constrained.contains(&Constraint::new(u, v)) {
            continue;
        }
        if let Some(&other) = across.get(&(v, u)) {
            let s = incircle(
                pts[u as usize],
                pts[v as usize],
                pts[apex as usize],
                pts[other as usize],
            );
            let inside = matches!(
                s,
                Certified::Certain {
                    sign: Sign::Positive,
                    ..
                }
            );
            assert!(
                !inside,
                "edge {u}-{v} is not locally Delaunay: {:?} {:?} apex {:?} other {:?}",
                pts[u as usize], pts[v as usize], pts[apex as usize], pts[other as usize]
            );
        }
    }
}

type Pt = (i128, i128);

fn cross(o: Pt, a: Pt, b: Pt) -> i128 {
    (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
}

/// Whether `q` lies on the closed segment `a`-`b`.
fn on_segment(a: Pt, b: Pt, q: Pt) -> bool {
    cross(a, b, q) == 0
        && a.0.min(b.0) <= q.0
        && q.0 <= a.0.max(b.0)
        && a.1.min(b.1) <= q.1
        && q.1 <= a.1.max(b.1)
}

/// Winding number of the closed polygon `ring` about `q`, off its boundary.
fn winding(ring: &[Pt], q: Pt) -> i128 {
    let mut w = 0;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if a.1 <= q.1 && q.1 < b.1 && cross(a, b, q) > 0 {
            w += 1;
        } else if b.1 <= q.1 && q.1 < a.1 && cross(a, b, q) < 0 {
            w -= 1;
        }
    }
    w
}

/// An independent integer check of a domain triangulation.
///
/// Coordinates must be integers. Scaled by three, a triangle's centroid is
/// an integer point, so no rounding is involved. A triangle must have
/// membership 1 there (the sum of `+|winding|` over outer rings and
/// `-|winding|` over holes) and its centroid must be off every ring edge;
/// the doubled area of the triangles must equal the signed sum of the
/// doubled ring areas; and a point is a corner exactly when it is on a ring
/// or inside the domain, and in `outside_points` otherwise.
pub fn assert_domain(name: &str, points: &[Point2], rings: &[Ring], tri: &Triangulation) {
    let pts: Vec<Pt> = points.iter().map(|p| (p.x as i128, p.y as i128)).collect();
    let polys: Vec<(RingKind, Vec<Pt>)> = rings
        .iter()
        .map(|r| {
            (
                r.kind,
                r.vertices.iter().map(|&v| pts[v as usize]).collect(),
            )
        })
        .collect();
    let scaled = |q: Pt| (3 * q.0, 3 * q.1);
    let on_ring = |q: Pt| {
        polys.iter().any(|(_, ring)| {
            (0..ring.len())
                .any(|i| on_segment(scaled(ring[i]), scaled(ring[(i + 1) % ring.len()]), q))
        })
    };
    let membership = |q: Pt| -> i128 {
        polys
            .iter()
            .map(|(kind, ring)| {
                let big: Vec<Pt> = ring.iter().map(|&v| scaled(v)).collect();
                let w = winding(&big, q).abs();
                if *kind == RingKind::Outer {
                    w
                } else {
                    -w
                }
            })
            .sum()
    };

    let mut doubled = 0;
    let mut edges = HashSet::new();
    let mut used = vec![false; points.len()];
    for t in tri.triangles().as_chunks::<3>().0 {
        let c: Vec<Pt> = t.iter().map(|&v| pts[v as usize]).collect();
        let q = (c[0].0 + c[1].0 + c[2].0, c[0].1 + c[1].1 + c[2].1);
        assert!(
            cross(c[0], c[1], c[2]) > 0,
            "{name}: triangle {t:?} is not counter-clockwise"
        );
        assert!(
            !on_ring(q),
            "{name}: triangle {t:?} has its centroid on a ring"
        );
        assert_eq!(
            membership(q),
            1,
            "{name}: triangle {t:?} is outside the domain"
        );
        doubled += cross(c[0], c[1], c[2]);
        for i in 0..3 {
            used[t[i] as usize] = true;
            assert!(
                edges.insert((t[i], t[(i + 1) % 3])),
                "{name}: a directed edge occurs twice"
            );
        }
    }
    let area: i128 = polys
        .iter()
        .map(|(kind, ring)| {
            let a: i128 = (0..ring.len())
                .map(|i| {
                    let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
                    p.0 * q.1 - p.1 * q.0
                })
                .sum();
            if *kind == RingKind::Outer {
                a.abs()
            } else {
                -a.abs()
            }
        })
        .sum();
    assert_eq!(
        doubled, area,
        "{name}: the triangles do not cover the domain"
    );

    let merged: HashSet<u32> = tri.merged_points().iter().map(|m| m.input).collect();
    for (k, &q) in pts.iter().enumerate() {
        if merged.contains(&(k as u32)) {
            continue;
        }
        let q = scaled(q);
        let reported = tri.outside_points().contains(&(k as u32));
        assert_eq!(
            used[k], !reported,
            "{name}: point {k} is a corner xor outside"
        );
        // A point on a ring may touch no domain triangle (a corner that a
        // hole fills), so only an off-ring point has a stated fate.
        if !on_ring(q) {
            assert_eq!(
                used[k],
                membership(q) == 1,
                "{name}: point {k} is used iff inside"
            );
        }
    }
    for c in tri.constraints() {
        assert!(
            edges.contains(&(c.a, c.b)) || edges.contains(&(c.b, c.a)),
            "{name}: constraint {c:?} is not an edge"
        );
    }
}

/// The triangles as a sorted list, each rotated to start at its lowest
/// corner, so two triangulations compare as sets.
pub fn triangle_set(tri: &Triangulation) -> Vec<[u32; 3]> {
    let mut out: Vec<[u32; 3]> = tri
        .triangles()
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let s = (0..3).min_by_key(|&i| t[i]).unwrap();
            [t[s], t[(s + 1) % 3], t[(s + 2) % 3]]
        })
        .collect();
    out.sort_unstable();
    out
}
