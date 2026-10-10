//! Collinear edges that share no point are not an intersection (#272).
//!
//! Operand validation asks whether non-adjacent edges meet. Its strict
//! test read the side of a rounded cross product, with `0.0` counted as a
//! side, so two disjoint edges on one line whose cross products rounded
//! to `0.0` and to a tiny positive residue counted as crossing, and
//! `overlay` refused its own settled output as `SelfIntersection`. Every
//! contact decision is now exact; the unit tests in `src/contact.rs`
//! check the predicate against an integer oracle, and these check rings.

use axiolid_core::{Frame2, Point2, Tolerance, Vec2};
use axiolid_overlay::{
    overlay, union_soup, FillRule, OverlayError, OverlayInput, OverlayOperation, Polygon, Region,
    Ring,
};

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

fn poly(points: &[(f64, f64)]) -> Polygon {
    Polygon {
        outer: Ring {
            points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        },
        holes: Vec::new(),
    }
}

fn union(polygons: Vec<Polygon>, tolerance: Tolerance) -> Result<Vec<Polygon>, OverlayError> {
    let input = OverlayInput {
        frame: frame(),
        polygons,
    };
    overlay(
        &input,
        &input,
        OverlayOperation::Union,
        FillRule::NonZero,
        tolerance,
    )
    .map(|result| result.polygons)
}

fn issue_tolerance() -> Tolerance {
    Tolerance::new(1e-15, 1e-15).unwrap()
}

/// The strict branch as it was before #272: the side of a rounded cross
/// product, `0.0` counted with the negatives.
fn rounded_crossing(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let cross = |a: Point2, b: Point2, c: Point2| (b - a).perp_dot(c - a);
    (cross(a, b, c) > 0.0) != (cross(a, b, d) > 0.0)
        && (cross(c, d, a) > 0.0) != (cross(c, d, b) > 0.0)
}

/// Three rectangles of a U in a frame turned by `angle` about `centre`:
/// posts `u` in [0, 0.3] and [1.2, 1.5], `v` in [0, 1]; a bar `u` in
/// [0, 1.5], `v` in [0.9, 1.2].
fn turned_u(angle: f64, centre: (f64, f64)) -> Vec<Polygon> {
    let (s, c) = angle.sin_cos();
    let at = |u: f64, v: f64| (centre.0 + u * c - v * s, centre.1 + u * s + v * c);
    let rect = |u0: f64, v0: f64, u1: f64, v1: f64| {
        poly(&[at(u0, v0), at(u1, v0), at(u1, v1), at(u0, v1)])
    };
    vec![
        rect(0.0, 0.0, 0.3, 1.0),
        rect(1.2, 0.0, 1.5, 1.0),
        rect(0.0, 0.9, 1.5, 1.2),
    ]
}

#[test]
fn the_issue_u_is_accepted_when_fed_back() {
    let pieces = vec![
        poly(&[
            (13.462, -2.583),
            (13.702424380556335, -2.4035672347810966),
            (13.10431516315999, -1.6021526329266438),
            (12.863890782603654, -1.7815853981455474),
        ]),
        poly(&[
            (14.423697522225343, -1.8652689391243857),
            (14.664121902781678, -1.685836173905482),
            (14.066012685385333, -0.8844215720510291),
            (13.825588304828997, -1.063854337269933),
        ]),
        poly(&[
            (12.923701704343289, -1.8617268583309925),
            (14.125823607124968, -0.9645630322364743),
            (13.946390841906064, -0.7241386516801386),
            (12.744268939124385, -1.6213024777746567),
        ]),
    ];
    let t = issue_tolerance();
    let output = union(pieces, t).unwrap();
    // The two collinear edges the issue names, still in the output, whose
    // rounded cross products read as a crossing.
    let ring = &output[0].outer.points;
    let edge = |p: (f64, f64), q: (f64, f64)| {
        let (p, q) = (Point2::new(p.0, p.1), Point2::new(q.0, q.1));
        assert!(ring.contains(&p) && ring.contains(&q));
        (p, q)
    };
    let (a, b) = edge(
        (14.664121902781678, -1.685836173905482),
        (14.125823607124968, -0.9645630322364743),
    );
    let (c, d) = edge(
        (14.066012685385333, -0.8844215720510291),
        (13.946390841906064, -0.7241386516801386),
    );
    assert!(rounded_crossing(a, b, c, d));
    let again = union(output.clone(), t).expect("settled output is a valid operand");
    assert_eq!(again, output);
    Region::new(output, t).expect("settled output is a valid region");
}

#[test]
fn turned_us_settle_into_valid_operands_near_and_far() {
    let t = issue_tolerance();
    let mut tried = 0;
    for centre in [
        (13.462, -2.583),
        (1e5 + 0.37, -1e5 + 0.11),
        (5e6 + 0.5, 5e6 - 0.25),
    ] {
        for step in 0..400 {
            let angle = (0.37 + 0.0173 * f64::from(step)).to_radians() + centre.0 * 1e-9;
            let output = union(turned_u(angle, centre), t)
                .unwrap_or_else(|e| panic!("U at {angle} about {centre:?}: {e:?}"));
            assert_eq!(output.len(), 1);
            union(output.clone(), t)
                .unwrap_or_else(|e| panic!("U at {angle} about {centre:?} fed back: {e:?}"));
            tried += 1;
        }
    }
    assert_eq!(tried, 1200);
}

// Rings on the grid of 2^-24 within 64 of the origin: offsets of 1e5 and
// 5e6 translate them exactly.
const SCALE: f64 = (1u64 << 24) as f64;
const OFFSETS: [(f64, f64); 4] = [(0.0, 0.0), (1e5, -1e5), (5e6, 5e6), (-5e6, 1e5)];

fn grid(points: &[(i64, i64)], (dx, dy): (f64, f64)) -> Vec<(f64, f64)> {
    points
        .iter()
        .map(|&(x, y)| {
            let p = (x as f64 / SCALE + dx, y as f64 / SCALE + dy);
            assert_eq!(
                ((p.0 - dx) * SCALE, (p.1 - dy) * SCALE),
                (x as f64, y as f64)
            );
            p
        })
        .collect()
}

fn back(polygons: &[Polygon], (dx, dy): (f64, f64)) -> Vec<Polygon> {
    let ring = |r: &Ring| Ring {
        points: r
            .points
            .iter()
            .map(|p| Point2::new(p.x - dx, p.y - dy))
            .collect(),
    };
    polygons
        .iter()
        .map(|p| Polygon {
            outer: ring(&p.outer),
            holes: p.holes.iter().map(ring).collect(),
        })
        .collect()
}

/// The ring's union with itself at every offset, translated back, is the
/// one near the origin.
fn same_everywhere(points: &[(i64, i64)], t: Tolerance) -> Result<Vec<Polygon>, OverlayError> {
    let near = union(vec![poly(&grid(points, OFFSETS[0]))], t);
    for offset in OFFSETS {
        let far = union(vec![poly(&grid(points, offset))], t);
        assert_eq!(far.map(|f| back(&f, offset)), near, "at {offset:?}");
    }
    near
}

/// `(x, y)` with `p * y - q * x == 1`, for coprime `p`, `q`.
fn unit_cross(p: i64, q: i64) -> (i64, i64) {
    // Extended Euclid on (p, q): s * p + r * q == 1, so x = -r, y = s.
    let (mut old_r, mut r) = (p, q);
    let (mut old_s, mut s) = (1i64, 0i64);
    let (mut old_t, mut t) = (0i64, 1i64);
    while r != 0 {
        let k = old_r / r;
        (old_r, r) = (r, old_r - k * r);
        (old_s, s) = (s, old_s - k * s);
        (old_t, t) = (t, old_t - k * t);
    }
    assert_eq!(old_r.abs(), 1, "p and q must be coprime");
    (-old_t * old_r, old_s * old_r)
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

/// U-shaped rings whose outer side is two edges `a`-`b` and `c`-`d`, a
/// gap apart, with `c` and `d` one or two grid units squared of cross
/// product off the line of `a`-`b`: the rounded cross products of such
/// ends are a few units of rounding residue or exactly `0.0`. Only rings
/// whose collinear edges the rounded strict test calls crossing are kept.
fn rounding_us() -> Vec<Vec<(i64, i64)>> {
    let mut state = 272u64;
    let mut next = |n: i64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) % (n as u64)) as i64
    };
    let mut out = Vec::new();
    for _ in 0..20_000 {
        if out.len() == 40 {
            break;
        }
        let (p, q) = ((1 << 26) + next(1 << 26), next(1 << 27) - (1 << 26));
        if gcd(p, q) != 1 {
            continue;
        }
        // The unit step off the line, moved along it to within half a
        // segment of the foot of `a`.
        let (x0, y0) = unit_cross(p, q);
        let k = ((x0 * p + y0 * q) as f64 / (p * p + q * q) as f64).round() as i64;
        let (x0, y0) = (x0 - k * p, y0 - k * q);
        let a = (next(1 << 20), next(1 << 20));
        let on = |t: i64, r: i64| (a.0 + t * p + r * x0, a.1 + t * q + r * y0);
        let b = on(1, 0);
        let (c, d) = (on(4, next(5) - 2), on(7, next(5) - 2));
        // Left of the line, `p` and `q` scaled down.
        let n = (-q / 4, p / 4);
        let plus = |v: (i64, i64), k: i64| (v.0 + k * n.0, v.1 + k * n.1);
        let ring = vec![a, b, plus(b, 1), plus(c, 1), c, d, plus(d, 3), plus(a, 3)];
        let f = grid(&ring, OFFSETS[0]);
        let pt = |i: usize| Point2::new(f[i].0, f[i].1);
        if rounded_crossing(pt(0), pt(1), pt(4), pt(5)) {
            out.push(ring);
        }
    }
    assert_eq!(out.len(), 40, "too few rings defeat the rounded signs");
    out
}

#[test]
fn collinear_disjoint_edges_with_rounded_zero_cross_products_are_valid() {
    let t = Tolerance::new(1e-12, 1e-12).unwrap();
    for ring in rounding_us() {
        let near = same_everywhere(&ring, t).unwrap_or_else(|e| panic!("{ring:?} refused: {e:?}"));
        // The union settles to one ring that is itself a valid operand.
        assert_eq!(near.len(), 1, "{ring:?}");
        union(near.clone(), t).unwrap_or_else(|e| panic!("{ring:?} fed back: {e:?}"));
    }
}

#[test]
fn collinear_overlapping_edges_are_refused_near_and_far() {
    let t = Tolerance::new(1e-12, 1e-12).unwrap();
    // On one slanted line: the ring runs out along it, up, back down a
    // post onto it, and back along the line over its own first edge.
    let u = (40_000_001, 29_999_999);
    let n = (-29_999_999 / 2, 40_000_001 / 2);
    let at = |k: i64, h: i64| (k * u.0 / 4 + h * n.0, k * u.1 / 4 + h * n.1);
    let overlapping = [
        at(0, 0),
        at(8, 0),
        at(8, 1),
        at(6, 1),
        at(6, 0),
        at(2, 0),
        at(2, 1),
        at(0, 1),
    ];
    assert_eq!(
        same_everywhere(&overlapping, t),
        Err(OverlayError::SelfIntersection)
    );
    // The same with a gap: the edge back along the line stops short of
    // the first one, a U.
    let disjoint = [
        at(0, 0),
        at(3, 0),
        at(3, 1),
        at(5, 1),
        at(5, 0),
        at(8, 0),
        at(8, 2),
        at(0, 2),
    ];
    assert!(same_everywhere(&disjoint, t).is_ok());
}

#[test]
fn endpoint_touches_and_crossings_are_refused_near_and_far() {
    let t = Tolerance::new(1e-12, 1e-12).unwrap();
    let u = (40_000_001, 29_999_999);
    let n = (-29_999_999, 40_000_001);
    let at = |k: i64, h: i64| (k * u.0 / 4 + h * n.0 / 4, k * u.1 / 4 + h * n.1 / 4);
    // A notch whose tip reaches down to the first edge: exactly onto it
    // (`at(4, 0)` is `u`, half way along `0`-`2u`), five grid units short
    // of it (`(-3, 4)` is nearly perpendicular to `u`, on the inside), and
    // as far past it.
    let notch = |tip: (i64, i64)| {
        [
            at(0, 0),
            at(8, 0),
            at(8, 4),
            at(5, 4),
            tip,
            at(3, 4),
            at(0, 4),
        ]
    };
    let tip = at(4, 0);
    assert_eq!(tip, u);
    assert_eq!(
        same_everywhere(&notch(tip), t),
        Err(OverlayError::SelfIntersection)
    );
    assert!(same_everywhere(&notch((tip.0 - 3, tip.1 + 4)), t).is_ok());
    assert_eq!(
        same_everywhere(&notch((tip.0 + 3, tip.1 - 4)), t),
        Err(OverlayError::SelfIntersection)
    );
    // A bow tie.
    let bow = [at(0, 0), at(8, 0), at(0, 4), at(8, 4)];
    assert_eq!(
        same_everywhere(&bow, t),
        Err(OverlayError::SelfIntersection)
    );
}

/// Random soups of turned rectangles and triangles, near and far from the
/// origin: each union is accepted, and so is the union fed back.
#[test]
fn settled_soups_round_trip() {
    let mut state = 0x272u64;
    let mut unit = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let t = Tolerance::new(1e-12, 1e-12).unwrap();
    for round in 0..60 {
        let base = [(0.0, 0.0), (1e5, -1e5), (5e6, 5e6)][round % 3];
        let mut rings = Vec::new();
        for _ in 0..6 {
            let centre = (base.0 + 4.0 * unit(), base.1 + 4.0 * unit());
            let (s, c) = (std::f64::consts::TAU * unit()).sin_cos();
            let (w, h) = (0.2 + 2.0 * unit(), 0.2 + 2.0 * unit());
            let at =
                |u: f64, v: f64| Point2::new(centre.0 + u * c - v * s, centre.1 + u * s + v * c);
            let points = if unit() < 0.5 {
                vec![at(0.0, 0.0), at(w, 0.0), at(w, h), at(0.0, h)]
            } else {
                vec![at(0.0, 0.0), at(w, 0.0), at(w * unit(), h)]
            };
            rings.push(Ring { points });
        }
        let soup = union_soup(&rings, t).unwrap_or_else(|e| panic!("round {round}: {e:?}"));
        let again =
            union(soup.clone(), t).unwrap_or_else(|e| panic!("round {round} fed back: {e:?}"));
        Region::new(soup, t).unwrap_or_else(|e| panic!("round {round} as a region: {e:?}"));
        union(again, t).unwrap_or_else(|e| panic!("round {round} fed back twice: {e:?}"));
    }
}

/// The ray-crossing test before #272: where the edge crosses the ray,
/// rounded.
fn rounded_contains(ring: &[(f64, f64)], p: (f64, f64)) -> bool {
    let mut inside = false;
    for (i, &a) in ring.iter().enumerate() {
        let b = ring[(i + 1) % ring.len()];
        if (a.1 > p.1) != (b.1 > p.1) && p.0 < (b.0 - a.0) * (p.1 - a.1) / (b.1 - a.1) + a.0 {
            inside = !inside;
        }
    }
    inside
}

#[test]
fn a_hole_corner_just_inside_a_slanted_outer_edge_is_inside() {
    // A turned square, and a hole whose corner `h` lies one grid unit
    // squared of cross product inside its first edge: a hair's breadth
    // from it, closer than the rounding of where that edge crosses a ray.
    let t = Tolerance::new(0.0, 0.0).unwrap();
    let mut state = 0x271u64;
    let mut next = |n: i64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) % (n as u64)) as i64
    };
    let mut rounded_wrong = 0;
    let mut tried = 0;
    while tried < 200 {
        let (p, q) = ((1 << 26) + next(1 << 26), next(1 << 26) + 1);
        if gcd(p, q) != 1 {
            continue;
        }
        tried += 1;
        let (x0, y0) = unit_cross(p, q);
        let k = ((x0 * p + y0 * q) as f64 / (p * p + q * q) as f64).round() as i64;
        let (x0, y0) = (x0 - k * p, y0 - k * q);
        let outer = [
            (0, 0),
            (2 * p, 2 * q),
            (2 * (p - q), 2 * (q + p)),
            (-2 * q, 2 * p),
        ];
        let inner = |s: i64, u: i64| ((s * p - u * q) / 4, (s * q + u * p) / 4);
        // `r == 1`: `h` is two units squared of cross product inside the
        // edge `(0, 0)`-`(2p, 2q)`, half way along it; `r == -1`: as far
        // outside.
        for r in [1, -1] {
            let h = (p + r * x0, q + r * y0);
            let hole = [h, inner(5, 3), inner(3, 3)];
            let polygon = |offset| Polygon {
                outer: poly(&grid(&outer, offset)).outer,
                holes: vec![poly(&grid(&hole, offset)).outer],
            };
            let near = union(vec![polygon(OFFSETS[0])], t);
            if r == 1 {
                assert!(near.is_ok(), "{outer:?} with {hole:?}: {near:?}");
            } else {
                assert_eq!(near, Err(OverlayError::HoleOutsideOuter));
            }
            for offset in OFFSETS {
                let far = union(vec![polygon(offset)], t);
                assert_eq!(far.map(|f| back(&f, offset)), near, "at {offset:?}");
                let g = grid(&outer, offset);
                let inside = rounded_contains(&g, grid(&[h], offset)[0]);
                rounded_wrong += usize::from(inside != (r == 1));
            }
        }
    }
    assert!(rounded_wrong > 0, "no hole defeats the rounded ray test");
}
