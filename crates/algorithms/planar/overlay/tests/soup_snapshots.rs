//! Exact outputs pinned bit for bit on plan-footprint soups (#315).
//!
//! The subdivision's speed-ups for soups of overlapping slender parts
//! (meetings merged by a sort, carriers searched among the edges sharing a
//! stretch, rings carried across any stop, repeated segments' crossings
//! shared, and an intersection's rings culled to the other operand's box)
//! change no output. Each case hashes every output coordinate's bits; the
//! hashes were taken from `axiolid-overlay` 0.3.14, before #315. A changed
//! hash is a changed output, which the exactness contract does not allow
//! without a stated reason.
//!
//! The soups are the ones axioval/engine measured: slender boxes
//! (1 m x 4 cm x 2 cm) at seeded turns and offsets, a torus and a lying
//! cylinder, each triangle projected to xy and wound counter-clockwise;
//! plus random triangles, whose segments no other boundary runs along, and
//! an arc subdivision with repeated and shared circles.

use axiolid_core::{Frame2, Point2, Tolerance};
use axiolid_overlay::{
    overlay, union_soup, ArcArrangement, ArcRing, ArcVertex, FillRule, OverlayInput,
    OverlayOperation, Polygon, Ring,
};

type P3 = [f64; 3];

/// FNV-1a over 64-bit words.
struct Hash(u64);

impl Hash {
    fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }

    fn word(&mut self, word: u64) {
        for byte in word.to_le_bytes() {
            self.0 ^= u64::from(byte);
            self.0 = self.0.wrapping_mul(0x0100_0000_01b3);
        }
    }

    fn polygons(&mut self, polygons: &[Polygon]) {
        self.word(polygons.len() as u64);
        for polygon in polygons {
            self.word(polygon.holes.len() as u64);
            for ring in std::iter::once(&polygon.outer).chain(&polygon.holes) {
                self.word(ring.points.len() as u64);
                for p in &ring.points {
                    self.word(p.x.to_bits());
                    self.word(p.y.to_bits());
                }
            }
        }
    }

    fn text(&mut self, text: &str) {
        for byte in text.bytes() {
            self.word(u64::from(byte));
        }
    }
}

/// `k` closed slender boxes as triangles, from a fixed xorshift seed.
fn assembly(k: usize) -> Vec<[P3; 3]> {
    let mut out = Vec::new();
    let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..k {
        let (cx, cy, cz, turn) = (next(), next(), next(), next() * std::f64::consts::PI);
        let (c, s) = (turn.cos(), turn.sin());
        let corner = |i: usize| {
            let x = if i & 1 == 0 { -0.5 } else { 0.5 };
            let y = if i & 2 == 0 { -0.02 } else { 0.02 };
            let z = if i & 4 == 0 { 0.0 } else { 0.02 };
            [cx + c * x - s * y, cy + s * x + c * y, cz + z]
        };
        let v: Vec<P3> = (0..8).map(corner).collect();
        for [a, b, c2, d] in [
            [0, 1, 3, 2],
            [4, 6, 7, 5],
            [0, 4, 5, 1],
            [2, 3, 7, 6],
            [0, 2, 6, 4],
            [1, 5, 7, 3],
        ] {
            out.push([v[a], v[b], v[c2]]);
            out.push([v[a], v[c2], v[d]]);
        }
    }
    out
}

/// A torus lying flat, `u` segments around, `v` around the tube.
fn torus(u: usize, v: usize) -> Vec<[P3; 3]> {
    let (big, small) = (1.0, 0.3);
    let at = |i: usize, j: usize| {
        let a = std::f64::consts::TAU * (i % u) as f64 / u as f64;
        let b = std::f64::consts::TAU * (j % v) as f64 / v as f64;
        let r = big + small * b.cos();
        [r * a.cos(), r * a.sin(), small * b.sin()]
    };
    let mut out = Vec::new();
    for i in 0..u {
        for j in 0..v {
            let (p, q, r, s) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            out.push([p, q, r]);
            out.push([p, r, s]);
        }
    }
    out
}

/// A cylinder lying along x, `n` segments around, `m` along, caps left
/// out.
fn lying_cylinder(n: usize, m: usize) -> Vec<[P3; 3]> {
    let at = |i: usize, k: usize| {
        let a = std::f64::consts::TAU * (i % n) as f64 / n as f64;
        [2.0 * k as f64 / m as f64, 0.2 * a.cos(), 0.2 * a.sin()]
    };
    let mut out = Vec::new();
    for i in 0..n {
        for k in 0..m {
            let (p, q, r, s) = (at(i, k), at(i + 1, k), at(i + 1, k + 1), at(i, k + 1));
            out.push([p, q, r]);
            out.push([p, r, s]);
        }
    }
    out
}

/// Random triangles near random centres, counter-clockwise.
fn random_triangles(count: usize) -> Vec<[P3; 3]> {
    let mut seed = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..count)
        .map(|_| {
            let (cx, cy) = (next(), next());
            let mut corner = || [cx + 0.2 * next() - 0.1, cy + 0.2 * next() - 0.1, 0.0];
            [corner(), corner(), corner()]
        })
        .collect()
}

/// Each triangle projected to xy and wound counter-clockwise, those
/// without area left out.
fn shadows(triangles: &[[P3; 3]]) -> Vec<Polygon> {
    triangles
        .iter()
        .filter_map(|t| {
            let mut points: Vec<Point2> = t.iter().map(|p| Point2::new(p[0], p[1])).collect();
            let (a, b, c) = (points[0], points[1], points[2]);
            let area = 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
            if area.abs() <= 1e-12 {
                return None;
            }
            if area < 0.0 {
                points.reverse();
            }
            Some(Polygon {
                outer: Ring { points },
                holes: Vec::new(),
            })
        })
        .collect()
}

fn input(polygons: Vec<Polygon>) -> OverlayInput {
    OverlayInput {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Point2::X,
            y: Point2::Y,
        },
        polygons,
    }
}

fn square(x: f64, y: f64, side: f64) -> Polygon {
    Polygon {
        outer: Ring {
            points: vec![
                Point2::new(x, y),
                Point2::new(x + side, y),
                Point2::new(x + side, y + side),
                Point2::new(x, y + side),
            ],
        },
        holes: Vec::new(),
    }
}

/// The hash of every overlay the consumer runs on one soup: the union of
/// the soup with itself, its union as a soup, and its intersection with,
/// and difference from, squares at `corners`.
fn soup_hash(soup: &[[P3; 3]], corners: &[(f64, f64, f64)]) -> u64 {
    let tolerance = Tolerance::new(1.0e-9, 1.0e-9).expect("a valid tolerance");
    let kept = shadows(soup);
    let rings: Vec<Ring> = kept.iter().map(|p| p.outer.clone()).collect();
    let body = input(kept);
    let mut hash = Hash::new();
    let union = overlay(
        &body,
        &body,
        OverlayOperation::Union,
        FillRule::NonZero,
        tolerance,
    )
    .expect("the union of valid shadows");
    hash.polygons(&union.polygons);
    hash.polygons(&union_soup(&rings, tolerance).expect("a valid soup"));
    for &(x, y, side) in corners {
        let small = input(vec![square(x, y, side)]);
        for (operation, fill) in [
            (OverlayOperation::Intersection, FillRule::NonZero),
            (OverlayOperation::Intersection, FillRule::EvenOdd),
            (OverlayOperation::Difference, FillRule::NonZero),
        ] {
            let result = overlay(&small, &body, operation, fill, tolerance)
                .expect("an overlay of valid operands");
            hash.polygons(&result.polygons);
        }
        // The soup as the subject: the cull drops subject rings too.
        let result = overlay(
            &body,
            &small,
            OverlayOperation::Intersection,
            FillRule::NonZero,
            tolerance,
        )
        .expect("an overlay of valid operands");
        hash.polygons(&result.polygons);
    }
    hash.0
}

/// Print every case's hash with `SNAPSHOT_PRINT=1`, to pin a deliberate
/// change.
fn check(name: &str, hash: u64, pinned: u64) {
    if std::env::var_os("SNAPSHOT_PRINT").is_some() {
        println!("{name}: {hash:#018x}");
        return;
    }
    assert_eq!(hash, pinned, "{name}: the output changed");
}

#[test]
fn slender_boxes_are_unchanged() {
    let corners = [(0.0, 0.0, 0.05), (0.45, 0.45, 0.05), (0.3, 0.2, 0.3)];
    check(
        "slender 25",
        soup_hash(&assembly(25), &corners),
        0x1c16_2366_3520_b19c,
    );
    check(
        "slender 60",
        soup_hash(&assembly(60), &corners),
        0x171f_5f87_36f1_a7d9,
    );
}

#[test]
fn meshes_are_unchanged() {
    check(
        "torus 24x12",
        soup_hash(&torus(24, 12), &[(1.0, 0.0, 0.05), (-0.2, -1.3, 0.4)]),
        0xdadc_2bf1_d5a1_4303,
    );
    check(
        "cylinder 24x8",
        soup_hash(
            &lying_cylinder(24, 8),
            &[(0.0, 0.0, 0.05), (1.0, -0.1, 0.5)],
        ),
        0x84c2_5fb6_e299_9904,
    );
}

#[test]
fn random_triangles_are_unchanged() {
    check(
        "random 300",
        soup_hash(&random_triangles(300), &[(0.2, 0.3, 0.05), (0.5, 0.5, 0.3)]),
        0x967f_0c52_eef4_a922,
    );
}

/// The arc subdivision shares the core: circles given twice, circles
/// crossing, half discs sharing an arc, and rectangles over them.
#[test]
fn arc_subdivision_is_unchanged() {
    let p = Point2::new;
    let half = |cx: f64, up: bool| {
        // A half disc of radius 1 about (cx, 0): the arc from one end of
        // the diameter to the other, then the diameter back.
        let (a, b) = if up {
            (p(cx + 1.0, 0.0), p(cx - 1.0, 0.0))
        } else {
            (p(cx - 1.0, 0.0), p(cx + 1.0, 0.0))
        };
        ArcRing::new(vec![ArcVertex::bulged(a, 1.0), ArcVertex::straight(b)])
    };
    let rect = |x0: f64, y0: f64, x1: f64, y1: f64| {
        ArcRing::new(vec![
            ArcVertex::straight(p(x0, y0)),
            ArcVertex::straight(p(x1, y0)),
            ArcVertex::straight(p(x1, y1)),
            ArcVertex::straight(p(x0, y1)),
        ])
    };
    let rings = vec![
        ArcRing::circle(p(0.0, 0.0), 1.0),
        ArcRing::circle(p(0.0, 0.0), 1.0),
        ArcRing::circle(p(0.7, 0.1), 1.0),
        half(0.0, true),
        half(0.0, false),
        half(0.5, true),
        rect(-0.5, -0.25, 1.5, 0.25),
        rect(-0.5, -0.25, 1.5, 0.25),
        rect(-1.0, 0.0, 0.0, 1.0),
    ];
    let tolerance = Tolerance::new(1.0e-9, 1.0e-9).expect("a valid tolerance");
    let arrangement = ArcArrangement::new(&rings, tolerance).expect("valid rings");
    let mut hash = Hash::new();
    for v in arrangement.vertices() {
        hash.word(v.x.to_bits());
        hash.word(v.y.to_bits());
    }
    for edge in arrangement.edges() {
        hash.word(edge.from as u64);
        hash.word(edge.to as u64);
        hash.word(edge.bulge.to_bits());
        hash.text(&format!("{:?}", edge.sources));
        for ring in 0..rings.len() {
            hash.word(u64::from(edge.inside_left(ring)));
            hash.word(u64::from(edge.inside_right(ring)));
        }
    }
    check("arcs", hash.0, 0x9c5f_b345_d813_8633);
}
