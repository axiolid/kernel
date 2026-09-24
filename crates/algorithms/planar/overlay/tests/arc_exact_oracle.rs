//! Exact arc boolean against independent oracles (#155).
//!
//! Scenes are snapped to a coarse grid on purpose: shared edges, identical
//! circles, tangent circles and vertices landing on the other boundary are
//! the cases a tolerance-based backend gets wrong, and they are frequent
//! here. Two oracles, neither calling the code under test:
//!
//! - Area identities: |A u B| + |A n B| = |A| + |B|, |A - B| = |A| - |A n B|,
//!   |A xor B| = |A u B| - |A n B|, with input areas from the closed-form
//!   `arc_ring_area`.
//! - Point membership: grid points classified against finely tessellated
//!   operands must match the result, away from every boundary.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{
    arc_overlay, arc_ring_area, validate_arc_ring, ArcOverlayResult, ArcRing, ArcVertex,
    OverlayOperation,
};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    /// A grid coordinate in [-2, 2], step 0.5.
    fn coord(&mut self) -> f64 {
        (self.below(9) as f64 - 4.0) * 0.5
    }
    fn radius(&mut self) -> f64 {
        (self.below(4) as f64 + 1.0) * 0.5
    }
}

fn v(x: f64, y: f64, b: f64) -> ArcVertex {
    ArcVertex::bulged(Point2::new(x, y), b)
}

fn shape(rng: &mut Rng) -> ArcRing {
    let (x, y) = (rng.coord(), rng.coord());
    let (w, h) = (rng.radius() * 2.0, rng.radius() * 2.0);
    match rng.below(6) {
        0 => ArcRing::circle(Point2::new(x, y), rng.radius()),
        1 => ArcRing::new(vec![
            v(x, y, 0.0),
            v(x + w, y, 0.0),
            v(x + w, y + h, 0.0),
            v(x, y + h, 0.0),
        ]),
        // Stadium: semicircular caps.
        2 => ArcRing::new(vec![
            v(x, y, 0.0),
            v(x + w, y, 1.0),
            v(x + w, y + h, 0.0),
            v(x, y + h, 1.0),
        ]),
        // A concave bite out of the bottom.
        3 => ArcRing::new(vec![
            v(x, y, -0.25),
            v(x + w, y, 0.0),
            v(x + w, y + h, 0.0),
            v(x, y + h, 0.0),
        ]),
        // A convex bulge on top.
        4 => ArcRing::new(vec![
            v(x, y, 0.0),
            v(x + w, y, 0.0),
            v(x + w, y + h, 0.5),
            v(x, y + h, 0.0),
        ]),
        // A circle drawn as four quarter arcs.
        _ => {
            let r = rng.radius();
            let b = (std::f64::consts::PI / 8.0).tan();
            ArcRing::new(vec![
                v(x + r, y, b),
                v(x, y + r, b),
                v(x - r, y, b),
                v(x, y - r, b),
            ])
        }
    }
}

/// Polyline approximation of a ring, fine enough that its error is far
/// below the membership margin.
fn tessellate(ring: &ArcRing) -> Vec<Point2> {
    let n = ring.vertices.len();
    let mut out = Vec::new();
    for i in 0..n {
        let a = ring.vertices[i];
        let b = ring.vertices[(i + 1) % n].point;
        out.push(a.point);
        if a.bulge == 0.0 {
            continue;
        }
        let theta = 4.0 * a.bulge.atan();
        let (dx, dy) = (b.x - a.point.x, b.y - a.point.y);
        let chord = (dx * dx + dy * dy).sqrt();
        let r = chord / (2.0 * (theta / 2.0).sin()).abs();
        let k = (1.0 - a.bulge * a.bulge) / (4.0 * a.bulge);
        let c = Point2::new(
            (a.point.x + b.x) / 2.0 - k * dy,
            (a.point.y + b.y) / 2.0 + k * dx,
        );
        let start = (a.point.y - c.y).atan2(a.point.x - c.x);
        let steps = ((theta.abs() / std::f64::consts::TAU) * 8192.0).ceil() as usize;
        for s in 1..steps {
            let t = start + theta * s as f64 / steps as f64;
            out.push(Point2::new(c.x + r * t.cos(), c.y + r * t.sin()));
        }
    }
    out
}

fn inside(poly: &[Point2], p: Point2) -> bool {
    let mut odd = false;
    let n = poly.len();
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        if (a.y > p.y) != (b.y > p.y) {
            let x = a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if p.x < x {
                odd = !odd;
            }
        }
    }
    odd
}

fn distance(poly: &[Point2], p: Point2) -> f64 {
    let n = poly.len();
    (0..n)
        .map(|i| {
            let (a, b) = (poly[i], poly[(i + 1) % n]);
            let (dx, dy) = (b.x - a.x, b.y - a.y);
            let len2 = dx * dx + dy * dy;
            let t = if len2 == 0.0 {
                0.0
            } else {
                (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0)
            };
            let (qx, qy) = (a.x + t * dx - p.x, a.y + t * dy - p.y);
            (qx * qx + qy * qy).sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

fn area(result: &ArcOverlayResult) -> f64 {
    result
        .regions
        .iter()
        .map(|r| {
            arc_ring_area(&r.outer).abs()
                - r.holes.iter().map(|h| arc_ring_area(h).abs()).sum::<f64>()
        })
        .sum()
}

fn rings(result: &ArcOverlayResult) -> Vec<Vec<Point2>> {
    result
        .regions
        .iter()
        .flat_map(|r| std::iter::once(&r.outer).chain(&r.holes))
        .map(tessellate)
        .collect()
}

const OPS: [OverlayOperation; 4] = [
    OverlayOperation::Union,
    OverlayOperation::Intersection,
    OverlayOperation::Difference,
    OverlayOperation::Xor,
];

fn expected(op: OverlayOperation, a: bool, b: bool) -> bool {
    match op {
        OverlayOperation::Union => a || b,
        OverlayOperation::Intersection => a && b,
        OverlayOperation::Difference => a && !b,
        OverlayOperation::Xor => a != b,
    }
}

fn scaled(ring: &ArcRing, s: f64) -> ArcRing {
    ArcRing::new(
        ring.vertices
            .iter()
            .map(|v| ArcVertex::bulged(Point2::new(v.point.x * s, v.point.y * s), v.bulge))
            .collect(),
    )
}

fn check_scene(a: &ArcRing, b: &ArcRing, label: &str) -> usize {
    check_scene_at(a, b, label, 1.0)
}

/// `a` and `b` are already at scale `s`; the membership grid and margins
/// follow it.
fn check_scene_at(a: &ArcRing, b: &ArcRing, label: &str, s: f64) -> usize {
    let tol = Tolerance::new(1e-9 * s, 1e-9).expect("valid tolerance");
    let results: Vec<ArcOverlayResult> = OPS
        .iter()
        .map(|&op| {
            arc_overlay(a, b, op, tol)
                .unwrap_or_else(|e| panic!("{label} {op:?}: refused with {e:?}"))
        })
        .collect();

    // Every output ring is a valid ring on its own.
    for (op, r) in OPS.iter().zip(&results) {
        for region in &r.regions {
            for ring in std::iter::once(&region.outer).chain(&region.holes) {
                validate_arc_ring(ring, tol)
                    .unwrap_or_else(|e| panic!("{label} {op:?}: invalid output ring {e:?}"));
            }
            assert!(
                arc_ring_area(&region.outer) > 0.0,
                "{label} {op:?}: outer not CCW"
            );
            for hole in &region.holes {
                assert!(arc_ring_area(hole) < 0.0, "{label} {op:?}: hole not CW");
            }
        }
    }

    let (area_a, area_b) = (arc_ring_area(a).abs(), arc_ring_area(b).abs());
    let [u, i, d, x] = [0, 1, 2, 3].map(|k| area(&results[k]));
    let eps = 1e-9 * (area_a + area_b).max(s * s);
    assert!(
        (u + i - area_a - area_b).abs() < eps,
        "{label}: |AuB|+|AnB| = {} vs {}",
        u + i,
        area_a + area_b
    );
    assert!(
        (d - (area_a - i)).abs() < eps,
        "{label}: |A-B| = {d} vs {}",
        area_a - i
    );
    assert!(
        (x - (u - i)).abs() < eps,
        "{label}: |AxB| = {x} vs {}",
        u - i
    );

    // Membership on a jittered grid, away from every boundary.
    let (ta, tb) = (tessellate(a), tessellate(b));
    let outs: Vec<Vec<Vec<Point2>>> = results.iter().map(rings).collect();
    let margin = 1e-5 * s;
    let mut checked = 0;
    for gx in 0..24 {
        for gy in 0..24 {
            let p = Point2::new((-4.3 + gx as f64 * 0.37) * s, (-4.3 + gy as f64 * 0.37) * s);
            if distance(&ta, p) < margin || distance(&tb, p) < margin {
                continue;
            }
            let (ia, ib) = (inside(&ta, p), inside(&tb, p));
            for (k, op) in OPS.iter().enumerate() {
                if outs[k].iter().any(|ring| distance(ring, p) < margin) {
                    continue;
                }
                let got = outs[k].iter().filter(|ring| inside(ring, p)).count() % 2 == 1;
                assert_eq!(
                    got,
                    expected(*op, ia, ib),
                    "{label} {op:?}: membership of ({}, {}) wrong",
                    p.x,
                    p.y
                );
                checked += 1;
            }
        }
    }
    checked
}

#[test]
fn random_grid_snapped_scenes_match_both_oracles() {
    let mut rng = Rng(0x05ee_d155);
    let mut checked = 0;
    for scene in 0..150 {
        let a = shape(&mut rng);
        let b = shape(&mut rng);
        checked += check_scene(&a, &b, &format!("scene {scene} a={a:?} b={b:?}"));
    }
    eprintln!("membership checks: {checked}");
    assert!(checked > 100_000, "only {checked} membership checks");
}

#[test]
fn an_operand_against_itself_and_its_redrawn_twin() {
    let disc = ArcRing::circle(Point2::new(0.5, -0.5), 1.5);
    check_scene(&disc, &disc, "identical discs");
    let square = ArcRing::from_points(&[
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(2.0, 2.0),
        Point2::new(0.0, 2.0),
    ]);
    // The same square, started at another vertex.
    let turned = ArcRing::from_points(&[
        Point2::new(2.0, 2.0),
        Point2::new(0.0, 2.0),
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
    ]);
    check_scene(&square, &turned, "same square, other start");
}

#[test]
fn tangent_and_concentric_circles() {
    let c = |x: f64, y: f64, r: f64| ArcRing::circle(Point2::new(x, y), r);
    check_scene(&c(0.0, 0.0, 1.0), &c(2.0, 0.0, 1.0), "externally tangent");
    check_scene(&c(0.0, 0.0, 2.0), &c(1.0, 0.0, 1.0), "internally tangent");
    check_scene(&c(0.0, 0.0, 2.0), &c(0.0, 0.0, 1.0), "concentric");
    check_scene(&c(0.0, 0.0, 1.0), &c(0.0, 0.0, 1.0), "same circle");
}

/// The same scenes on a decimal grid (step 0.05): coordinates binary cannot
/// hold, so constructed crossings land within rounding of vertices and
/// output rounding has to keep every ring valid.
#[test]
fn decimal_grid_scenes_match_both_oracles() {
    let mut rng = Rng(0x5eed_0155);
    let mut checked = 0;
    for scene in 0..120 {
        let a = scaled(&shape(&mut rng), 0.1);
        let b = scaled(&shape(&mut rng), 0.1);
        checked += check_scene_at(
            &a,
            &b,
            &format!("decimal scene {scene} a={a:?} b={b:?}"),
            0.1,
        );
    }
    assert!(checked > 80_000, "only {checked} membership checks");
}

/// A disc through the corners of a semicircular wall end, at 0.3 m wall
/// thickness: the crossings sit within rounding of the wall's vertices.
/// Rounding used to leave repeated vertices in the union and difference.
#[test]
fn a_disc_through_a_rounded_wall_end_gives_valid_rings() {
    let wall = ArcRing::new(vec![
        v(0.0, 0.0, 0.0),
        v(5.0, 0.0, 1.0),
        v(5.0, 0.3, 0.0),
        v(0.0, 0.3, 0.0),
    ]);
    let disc = ArcRing::circle(Point2::new(5.2, 0.15), 0.25);
    check_scene_at(&wall, &disc, "wall end", 1.0);
}

/// Regions and total vertex count of a union of shapes touching at one
/// vertex.
fn touching_union(a: &ArcRing, b: &ArcRing) -> (usize, Vec<usize>) {
    let result = arc_overlay(a, b, OverlayOperation::Union, Tolerance::METRE).expect("union");
    let mut sizes: Vec<usize> = result
        .regions
        .iter()
        .map(|r| r.outer.vertices.len())
        .collect();
    sizes.sort_unstable();
    (result.regions.len(), sizes)
}

/// Where two result boundaries meet at a vertex, linking must take the
/// sharpest left turn so each region comes out as its own simple ring. The
/// wrong turn gives one ring that touches itself (a figure eight), which
/// still has the right area and membership, so only the ring structure
/// shows it.
#[test]
fn shapes_touching_at_a_vertex_stay_separate_rings() {
    let square = ArcRing::from_points(&[
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(2.0, 2.0),
        Point2::new(0.0, 2.0),
    ]);
    // Opposite quadrant: one leaving piece turns left, the other right.
    let opposite = ArcRing::from_points(&[
        Point2::new(0.0, 0.0),
        Point2::new(0.0, -2.0),
        Point2::new(-2.0, -2.0),
        Point2::new(-2.0, 0.0),
    ]);
    assert_eq!(touching_union(&square, &opposite), (2, vec![4, 4]));
    // A wedge below the square's corner: both leaving pieces turn left, so
    // the choice is between two left turns.
    let wedge = ArcRing::from_points(&[
        Point2::new(0.0, 0.0),
        Point2::new(1.0, -2.0),
        Point2::new(2.0, -1.0),
    ]);
    assert_eq!(touching_union(&square, &wedge), (2, vec![3, 4]));
    // Union is symmetric, but linking starts from a different piece: with
    // the square as clip, the trace arrives at the shared corner down the
    // square's left edge and meets two left turns, so only this order
    // exercises the choice between them.
    assert_eq!(touching_union(&wedge, &square), (2, vec![3, 4]));
    assert_eq!(touching_union(&opposite, &square), (2, vec![4, 4]));
    // The same with a disc touching the square's corner from outside.
    let disc = ArcRing::new(vec![v(0.0, 0.0, 1.0), v(-1.0, -1.0, 1.0)]);
    let (regions, _) = touching_union(&square, &disc);
    assert_eq!(regions, 2);
}
