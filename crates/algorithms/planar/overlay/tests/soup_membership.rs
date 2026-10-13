//! Unions of triangle soups against point membership (#218).
//!
//! The subdivision carries the rings containing a piece across a stop
//! where exactly one other segment crosses, and takes a segment no other
//! boundary runs along as carried by its own ring only; everywhere else it
//! counts windings from a sample. Soups snapped to a coarse grid put
//! vertices on other edges, edges along each other and several edges
//! through one point, where the shortcuts must give way; soups of random
//! doubles cross simply almost everywhere, where they apply. Either way a
//! point away from every edge lies in the union exactly when it lies in
//! some triangle. (Debug builds also check every shortcut against the
//! count it replaces.)

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{union_soup, Polygon, Ring};

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn cross(a: Point2, b: Point2, c: Point2) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Counter-clockwise triangles with vertices from `vertex`, none flat.
fn soup(rng: &mut Lcg, count: usize, mut vertex: impl FnMut(&mut Lcg) -> Point2) -> Vec<Ring> {
    let mut out = Vec::new();
    while out.len() < count {
        let mut p = [vertex(rng), vertex(rng), vertex(rng)];
        let area = cross(p[0], p[1], p[2]);
        if area.abs() < 1e-3 {
            continue;
        }
        if area < 0.0 {
            p.swap(1, 2);
        }
        out.push(Ring { points: p.to_vec() });
    }
    out
}

/// Distance from `p` to the segment `a b`.
fn distance(p: Point2, a: Point2, b: Point2) -> f64 {
    let d = b - a;
    let t = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
    (p - (a + d * t)).length()
}

fn edges(ring: &Ring) -> impl Iterator<Item = (Point2, Point2)> + '_ {
    let n = ring.points.len();
    (0..n).map(move |i| (ring.points[i], ring.points[(i + 1) % n]))
}

/// Crossing-number membership; `p` is far from every edge.
fn in_ring(ring: &Ring, p: Point2) -> bool {
    edges(ring).fold(false, |inside, (a, b)| {
        if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y) {
            !inside
        } else {
            inside
        }
    })
}

fn in_union(polygons: &[Polygon], p: Point2) -> bool {
    polygons
        .iter()
        .any(|poly| in_ring(&poly.outer, p) && !poly.holes.iter().any(|h| in_ring(h, p)))
}

/// Checks `points` sample points of the soup's box against the union.
fn check(name: &str, rings: &[Ring], rng: &mut Lcg, points: usize) {
    let union = union_soup(rings, Tolerance::new(1e-9, 1e-9).unwrap()).expect("valid rings");
    let (mut lo, mut hi) = (
        Point2::new(f64::MAX, f64::MAX),
        Point2::new(f64::MIN, f64::MIN),
    );
    for p in rings.iter().flat_map(|r| &r.points) {
        lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
        hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
    }
    let mut asked = 0;
    while asked < points {
        let p = Point2::new(
            lo.x + (hi.x - lo.x) * rng.next(),
            lo.y + (hi.y - lo.y) * rng.next(),
        );
        if rings
            .iter()
            .flat_map(edges)
            .any(|(a, b)| distance(p, a, b) < 1e-6)
        {
            continue;
        }
        let want = rings.iter().any(|r| in_ring(r, p));
        assert_eq!(in_union(&union, p), want, "{name}: point {p:?}");
        asked += 1;
    }
}

#[test]
fn grid_snapped_soups_match_membership() {
    // Vertices on a half-unit grid: shared vertices, vertices on other
    // edges, collinear overlaps and many edges through one point.
    let mut rng = Lcg(0x218);
    for case in 0..40 {
        let (count, extent) = ([6, 20, 60][case % 3], [3.0, 6.0][case % 2]);
        let rings = soup(&mut rng, count, |r| {
            let snap = |v: f64| (v * extent * 2.0).round() / 2.0;
            Point2::new(snap(r.next()), snap(r.next()))
        });
        check(&format!("grid case {case}"), &rings, &mut rng, 200);
    }
}

#[test]
fn random_soups_match_membership() {
    // Random doubles: nearly every stop a simple crossing of two segments.
    let mut rng = Lcg(7);
    for case in 0..12 {
        let count = [30, 120, 400][case % 3];
        let rings = soup(&mut rng, count, |r| {
            let (cx, cy) = (r.next() * 20.0, r.next() * 20.0);
            Point2::new(cx, cy)
        });
        // Keep triangles small enough to overlap a few others each.
        let rings: Vec<Ring> = rings
            .into_iter()
            .map(|ring| {
                let c = ring.points[0];
                Ring {
                    points: ring.points.iter().map(|p| c + (*p - c) * 0.25).collect(),
                }
            })
            .collect();
        check(&format!("random case {case}"), &rings, &mut rng, 300);
    }
}

#[test]
fn mixed_soups_match_membership() {
    // Random triangles over grid-snapped ones: simple crossings next to
    // stops where a random edge passes a grid vertex or overlaps nothing.
    let mut rng = Lcg(0x5eed);
    for case in 0..20 {
        let mut rings = soup(&mut rng, 15, |r| {
            Point2::new((r.next() * 8.0).round(), (r.next() * 8.0).round())
        });
        rings.extend(soup(&mut rng, 15, |r| {
            Point2::new(r.next() * 8.0, r.next() * 8.0)
        }));
        check(&format!("mixed case {case}"), &rings, &mut rng, 300);
    }
}

#[test]
fn a_vertex_within_the_tolerance_of_a_long_edge_snaps_onto_it() {
    // A triangle's tip 4 mm above a slab's 100 m top edge, among hundreds
    // of triangles in that edge's band of `x`: the snap finds the
    // tip through its grid (#218) and puts it on the edge (#222), so the
    // slab's outline runs through the tip.
    let mut rng = Lcg(0x222);
    let tip = Point2::new(50.0, 0.004);
    let mut rings = vec![
        Ring {
            points: vec![
                Point2::new(0.0, -10.0),
                Point2::new(100.0, -10.0),
                Point2::new(100.0, 0.0),
                Point2::new(0.0, 0.0),
            ],
        },
        Ring {
            points: vec![tip, Point2::new(51.0, 3.0), Point2::new(49.0, 3.0)],
        },
    ];
    rings.extend(soup(&mut rng, 300, |r| {
        Point2::new(r.next() * 100.0, 10.0 + r.next() * 40.0)
    }));
    let union = union_soup(&rings, Tolerance::new(0.01, 1e-9).unwrap()).expect("valid rings");
    let holding = union
        .iter()
        .filter(|p| p.outer.points.iter().any(|q| q.y < -5.0))
        .collect::<Vec<_>>();
    assert_eq!(holding.len(), 1);
    assert!(
        holding[0].outer.points.contains(&tip),
        "the tip lies on the slab's outline"
    );
}
