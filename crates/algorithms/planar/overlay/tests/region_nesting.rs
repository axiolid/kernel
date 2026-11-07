//! A region's rings are apart and nested (#448): `Region::new` refuses a
//! hole that crosses or leaves its outer ring, holes that overlap, and
//! polygons that overlap, by name, and accepts rings that only touch.
//!
//! `scripts/probe_region_nesting_mutants.py` lists the faults these tests
//! must catch.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{OverlayError, Polygon, Region, Ring};

fn t() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).unwrap()
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
    }
}

fn polygon(outer: &[(f64, f64)], holes: &[&[(f64, f64)]]) -> Polygon {
    Polygon {
        outer: ring(outer),
        holes: holes.iter().map(|h| ring(h)).collect(),
    }
}

fn verdict(polygons: Vec<Polygon>) -> Result<(), OverlayError> {
    Region::new(polygons, t()).map(|_| ())
}

const SQUARE: &[(f64, f64)] = &[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<(f64, f64)> {
    vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
}

const SEED_21351_HOLE: [(f64, f64); 54] = [
    (0.0, 1.8250074780188603),
    (0.07306154548026278, 2.0991951843234236),
    (0.11438570205247789, 1.8250074780188603),
    (0.1474844072045125, 2.1938799638241373),
    (0.1985538174806416, 1.8250074780188603),
    (0.22991781876370188, 2.223986074674239),
    (0.28528732254325995, 1.8250074780188603),
    (0.3207154829973691, 2.228547273896655),
    (0.3547453025024969, 1.8250074780188603),
    (0.4073068899993937, 2.1430084939613883),
    (0.43801704724942336, 1.8250074780188603),
    (0.49099392341724757, 2.221509428495067),
    (0.542718327899733, 1.8250074780188603),
    (0.5724573578308509, 2.236958347144714),
    (0.6216488767927268, 1.8250074780188603),
    (0.6646482460343657, 2.1277542426900213),
    (0.7021323305059841, 1.8250074780188603),
    (0.7388730972608316, 2.19578405838923),
    (0.7795191733522651, 1.8250074780188603),
    (0.8312433057541164, 2.1353581536026285),
    (0.8697436926877211, 1.8250074780188603),
    (0.9173382783376326, 2.251833728600097),
    (0.9428007705603633, 1.8250074780188603),
    (0.989578430304246, 2.295519651785309),
    (1.0497734043196263, 1.8250074780188603),
    (1.0917848174290699, 2.2884762747851033),
    (1.1138179878109766, 1.8250074780188603),
    (1.157215207039751, 2.2800108227061298),
    (1.1963548413723635, 1.8250074780188603),
    (1.2443524747444792, 2.1027975770481406),
    (1.2947553148090476, 1.8250074780188603),
    (1.3315885345421732, 2.1426502674755468),
    (1.3695600174830478, 1.8250074780188603),
    (1.4138763390738327, 2.101703744230068),
    (1.4694448582765696, 1.8250074780188603),
    (1.4949377691448262, 2.248565515905763),
    (1.5434544459382795, 1.8250074780188603),
    (1.5815483797729528, 2.197337998478368),
    (1.637109376240751, 1.8250074780188603),
    (1.6670534677799458, 2.197261132545407),
    (1.7166340558758815, 1.8250074780188603),
    (1.7166340558758815, 0.0),
    (1.7018529054072407, -0.202285916771616),
    (1.5747073371736333, 0.0),
    (1.374671011874045, -0.1886012581217287),
    (1.255652212436378, 0.0),
    (1.1172229364443222, -0.21598024533661572),
    (0.9600844436860753, 0.0),
    (0.7663577569770754, -0.25671560969723556),
    (0.5932937909584324, 0.0),
    (0.44315876263202025, -0.27893003272558586),
    (0.3334328650246222, 0.0),
    (0.15577189883950887, -0.15683685455927415),
    (0.0, 0.0),
];

/// Fuzz seed 21351 of #345's review: a toothed hole whose needles reach
/// y = 2.2955, past its outer ring's top edge at y = 2.1202. Fast
/// morphology's notch fill took the needles for the set's boundary and
/// `dilate_outer_with(0.9046, fast(1e-2))` left out (2.0172, 3.0256),
/// 0.8993 from the region, with a stated deviation of 4.7e-3.
#[test]
fn seed_21351_a_hole_crossing_its_outer_ring_is_refused() {
    let (lo, x1, y1) = (-0.2951669615187855, 2.011801017394667, 2.120174439537646);
    let outer = [(lo, lo), (x1, lo), (x1, y1), (lo, y1)];
    assert_eq!(
        verdict(vec![polygon(&outer, &[&SEED_21351_HOLE[..]])]),
        Err(OverlayError::HoleOutsideOuter)
    );
}

#[test]
fn holes_inside_their_outer_ring_are_accepted_however_they_touch_it() {
    for (name, hole) in [
        ("inside", rect(1.0, 1.0, 2.0, 2.0)),
        (
            "given counter-clockwise",
            vec![(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)],
        ),
        (
            "touching at a vertex",
            vec![(1.0, 4.0), (3.0, 3.0), (2.0, 2.0)],
        ),
        (
            "touching at a corner",
            vec![(4.0, 4.0), (2.0, 3.0), (3.0, 2.0)],
        ),
        (
            "touching every side",
            vec![(2.0, 0.0), (4.0, 2.0), (2.0, 4.0), (0.0, 2.0)],
        ),
        ("sharing a stretch from inside", rect(1.0, 0.0, 3.0, 1.0)),
        // On one line: encloses nothing, as a boolean reads it (#219).
        (
            "on one line, crossing",
            vec![(3.0, 2.0), (5.0, 2.0), (6.0, 2.0)],
        ),
    ] {
        assert_eq!(
            verdict(vec![polygon(SQUARE, &[&hole[..]])]),
            Ok(()),
            "{name}"
        );
    }
}

#[test]
fn a_hole_crossing_or_leaving_its_outer_ring_is_refused() {
    for (name, hole) in [
        ("crossing", rect(3.0, 1.0, 5.0, 2.0)),
        ("outside", rect(5.0, 1.0, 6.0, 2.0)),
        ("around the outer ring", rect(-1.0, -1.0, 5.0, 5.0)),
        // No two edges cross properly: it leaves through two vertices on
        // the top edge.
        (
            "through two vertices",
            vec![(2.0, 3.0), (3.0, 4.0), (2.0, 5.0), (1.0, 4.0)],
        ),
        ("sharing a stretch from outside", rect(1.0, -1.0, 3.0, 0.0)),
    ] {
        assert_eq!(
            verdict(vec![polygon(SQUARE, &[&hole[..]])]),
            Err(OverlayError::HoleOutsideOuter),
            "{name}"
        );
    }
    // An outer ring on one line encloses no hole.
    assert_eq!(
        verdict(vec![polygon(
            &[(0.0, 0.0), (4.0, 0.0), (2.0, 0.0)],
            &[&rect(1.0, 1.0, 2.0, 2.0)[..]]
        )]),
        Err(OverlayError::HoleOutsideOuter)
    );
}

#[test]
fn holes_of_one_polygon_may_touch_but_not_overlap() {
    let a = rect(1.0, 1.0, 2.0, 2.0);
    for (name, b, want) in [
        ("touching at a vertex", rect(2.0, 2.0, 3.0, 3.0), Ok(())),
        ("sharing an edge", rect(2.0, 1.0, 3.0, 2.0), Ok(())),
        (
            "crossing",
            rect(1.5, 0.5, 2.5, 3.0),
            Err(OverlayError::HolesOverlap),
        ),
        (
            "nested",
            rect(0.5, 0.5, 3.5, 3.5),
            Err(OverlayError::HolesOverlap),
        ),
        ("the same", a.clone(), Err(OverlayError::HolesOverlap)),
    ] {
        assert_eq!(
            verdict(vec![polygon(SQUARE, &[&a[..], &b[..]])]),
            want,
            "{name}"
        );
    }
}

#[test]
fn polygons_may_touch_but_not_overlap() {
    let island = |x0, y0, x1, y1| polygon(&rect(x0, y0, x1, y1), &[]);
    let holed = || polygon(SQUARE, &[&rect(1.0, 1.0, 3.0, 3.0)[..]]);
    for (name, polygons, want) in [
        (
            "sharing an edge",
            vec![island(0.0, 0.0, 4.0, 4.0), island(4.0, 0.0, 6.0, 4.0)],
            Ok(()),
        ),
        (
            "sharing part of an edge",
            vec![island(0.0, 0.0, 4.0, 4.0), island(4.0, 1.0, 6.0, 2.0)],
            Ok(()),
        ),
        (
            "touching at a corner",
            vec![island(0.0, 0.0, 4.0, 4.0), island(4.0, 4.0, 6.0, 6.0)],
            Ok(()),
        ),
        (
            "an island in a hole",
            vec![holed(), island(1.5, 1.5, 2.5, 2.5)],
            Ok(()),
        ),
        (
            "an island filling a hole",
            vec![holed(), island(1.0, 1.0, 3.0, 3.0)],
            Ok(()),
        ),
        (
            "an island touching its hole",
            vec![holed(), polygon(&[(1.0, 1.0), (2.0, 1.5), (1.5, 2.0)], &[])],
            Ok(()),
        ),
        (
            // Out through the stretch a hole shares with its outer ring:
            // half in the hole, half outside, in neither polygon twice.
            "through a hole's open side",
            vec![
                polygon(SQUARE, &[&rect(1.0, 0.0, 3.0, 2.0)[..]]),
                island(1.5, -1.0, 2.5, 1.0),
            ],
            Ok(()),
        ),
        (
            "inside another",
            vec![island(0.0, 0.0, 4.0, 4.0), island(1.0, 1.0, 2.0, 2.0)],
            Err(OverlayError::PolygonsOverlap),
        ),
        (
            "crossing",
            vec![island(0.0, 0.0, 4.0, 4.0), island(3.0, 1.0, 5.0, 2.0)],
            Err(OverlayError::PolygonsOverlap),
        ),
        (
            "the same",
            vec![island(0.0, 0.0, 4.0, 4.0), island(0.0, 0.0, 4.0, 4.0)],
            Err(OverlayError::PolygonsOverlap),
        ),
        (
            "inside, sharing an edge",
            vec![island(0.0, 0.0, 4.0, 4.0), island(0.0, 0.0, 2.0, 1.0)],
            Err(OverlayError::PolygonsOverlap),
        ),
    ] {
        assert_eq!(verdict(polygons), want, "{name}");
    }
}

/// What a polygon set of axis-aligned rectangles on the unit grid breaks,
/// read at the centre of every unit cell: every face of their arrangement
/// is a union of cells, so the centres see every face.
fn cell_oracle(rects: &[([i64; 4], Vec<[i64; 4]>)]) -> Vec<Result<(), OverlayError>> {
    let inside = |r: &[i64; 4], x: f64, y: f64| {
        x > r[0] as f64 && x < r[2] as f64 && y > r[1] as f64 && y < r[3] as f64
    };
    let (mut outside, mut crowded, mut overlap) = (false, false, false);
    for i in 0..9 {
        for j in 0..9 {
            let (x, y) = (f64::from(i) + 0.5, f64::from(j) + 0.5);
            let mut solid = 0;
            for (outer, holes) in rects {
                let held = holes.iter().filter(|h| inside(h, x, y)).count();
                crowded |= held >= 2;
                outside |= held >= 1 && !inside(outer, x, y);
                solid += usize::from(inside(outer, x, y) && held == 0);
            }
            overlap |= solid >= 2;
        }
    }
    let mut found = Vec::new();
    if outside {
        found.push(Err(OverlayError::HoleOutsideOuter));
    }
    if crowded {
        found.push(Err(OverlayError::HolesOverlap));
    }
    if overlap {
        found.push(Err(OverlayError::PolygonsOverlap));
    }
    if found.is_empty() {
        found.push(Ok(()));
    }
    found
}

/// Random rectangles on a 9 x 9 grid, as outer rings and holes (most holes
/// inside their outer ring), touching, sharing edges and overlapping in
/// every way the grid allows. Each case runs as given, or turned by the
/// scaled rotation `(3x - 4y, 4x + 3y)`, and moved to 5e6, with each ring
/// in either orientation: the verdict is one of the defects the oracle
/// finds, and `Ok` exactly when it finds none.
#[test]
fn the_verdict_agrees_with_a_cell_oracle_on_rectangles() {
    let mut s = 0x448u64;
    let mut next = |n: u64| {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s % n
    };
    let mut seen = std::collections::BTreeMap::new();
    for case in 0..4000u64 {
        let mut rects: Vec<([i64; 4], Vec<[i64; 4]>)> = Vec::new();
        for _ in 0..1 + next(3) {
            let x0 = next(7) as i64;
            let x1 = x0 + 1 + next((8 - x0) as u64) as i64;
            let y0 = next(7) as i64;
            let y1 = y0 + 1 + next((8 - y0) as u64) as i64;
            let outer = [x0, y0, x1, y1];
            let mut holes = Vec::new();
            for _ in 0..next(3) {
                let (lo, hi) = if next(5) == 0 {
                    ([0, 0], [8, 8])
                } else {
                    ([x0, y0], [x1, y1])
                };
                let hx0 = lo[0] + next((hi[0] - lo[0]) as u64) as i64;
                let hx1 = hx0 + 1 + next((hi[0] - hx0) as u64) as i64;
                let hy0 = lo[1] + next((hi[1] - lo[1]) as u64) as i64;
                let hy1 = hy0 + 1 + next((hi[1] - hy0) as u64) as i64;
                holes.push([hx0, hy0, hx1, hy1]);
            }
            rects.push((outer, holes));
        }
        let want = cell_oracle(&rects);
        let (turn, far) = (case % 2 == 1, if case % 3 == 0 { 5.0e6 } else { 0.0 });
        let at = |x: i64, y: i64| {
            let (u, v) = if turn {
                (3 * x - 4 * y, 4 * x + 3 * y)
            } else {
                (x, y)
            };
            Point2::new(u as f64 + far, v as f64 + far)
        };
        let mut flip = case;
        let mut to_ring = |r: &[i64; 4]| {
            let mut points = vec![
                at(r[0], r[1]),
                at(r[2], r[1]),
                at(r[2], r[3]),
                at(r[0], r[3]),
            ];
            flip = flip.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            if flip >> 63 == 1 {
                points.reverse();
            }
            Ring { points }
        };
        let polygons: Vec<Polygon> = rects
            .iter()
            .map(|(outer, holes)| Polygon {
                outer: to_ring(outer),
                holes: holes.iter().map(&mut to_ring).collect(),
            })
            .collect();
        let got = verdict(polygons);
        assert!(
            want.contains(&got),
            "case {case}: {got:?}, want one of {want:?}: {rects:?}"
        );
        *seen.entry(format!("{got:?}")).or_insert(0usize) += 1;
    }
    // Every verdict occurs often.
    assert_eq!(seen.len(), 4, "{seen:?}");
    assert!(seen.values().all(|&n| n >= 300), "{seen:?}");
}

/// The floor of `benches/ring_validation.rs`: 40 m by 30 m, teeth 0.5 mm
/// deep along both long walls.
fn floor(corners: usize) -> Vec<(f64, f64)> {
    let teeth = (corners - 4) / 4;
    let (w, h) = (40.0, 30.0);
    let step = w / (teeth as f64 + 1.0);
    let depth = 5.0e-4;
    let mut points = vec![(0.0, 0.0)];
    for k in 0..teeth {
        let x = step * (k as f64 + 1.0);
        points.push((x - 0.25 * step, 0.0));
        points.push((x, depth));
    }
    points.push((w, 0.0));
    points.push((w, h));
    for k in 0..teeth {
        let x = w - step * (k as f64 + 1.0);
        points.push((x + 0.25 * step, h));
        points.push((x, h - depth));
    }
    points.push((0.0, h));
    points
}

#[test]
fn a_floor_of_100000_corners_with_200_holes_is_checked() {
    let outer = floor(100_000);
    let mut holes: Vec<Vec<(f64, f64)>> = Vec::new();
    for i in 0..20 {
        for j in 0..10 {
            let (x, y) = (1.0 + 1.9 * f64::from(i), 1.0 + 2.8 * f64::from(j));
            holes.push(rect(x, y, x + 0.5, y + 0.5));
        }
    }
    let refs: Vec<&[(f64, f64)]> = holes.iter().map(Vec::as_slice).collect();
    assert_eq!(verdict(vec![polygon(&outer, &refs)]), Ok(()));
    // One hole moved across the toothed wall.
    holes[0] = rect(1.0, -0.25, 1.5, 0.25);
    let refs: Vec<&[(f64, f64)]> = holes.iter().map(Vec::as_slice).collect();
    assert_eq!(
        verdict(vec![polygon(&outer, &refs)]),
        Err(OverlayError::HoleOutsideOuter)
    );
}
