// SPDX-License-Identifier: MPL-2.0

//! Exactly coincident input points are merged and reported (#300).
//!
//! Before, a repeated point was split into the mesh at its twin's position:
//! zero-area and overlapping triangles, one directed edge used twice, and
//! nothing said. Now each coincident point is merged into the lowest input
//! index at its position, `merged_points` lists the merge, and a constraint
//! whose two endpoints coincide is refused by name. Every case here is held
//! to a valid triangulation of the distinct points that covers their exact
//! hull, with the merge report a brute force expects.

mod common;

use std::collections::HashMap;

use axiolid_core::Point2;
use axiolid_triangulate::{triangulate, Constraint, Triangulation, TriangulationError};
use common::{assert_covers_hull, uniform, Rng};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

/// The lowest input index at each point's position.
fn expected_kept(points: &[Point2]) -> Vec<u32> {
    let mut first: HashMap<(u64, u64), u32> = HashMap::new();
    points
        .iter()
        .enumerate()
        .map(|(k, q)| {
            // `+ 0.0` reads -0.0 as 0.0: the same position.
            let key = ((q.x + 0.0).to_bits(), (q.y + 0.0).to_bits());
            *first.entry(key).or_insert(k as u32)
        })
        .collect()
}

/// Fail unless `tri` merged exactly the coincident points of `points`,
/// uses only kept vertices, carries every constraint over to them, and
/// covers the exact hull.
fn assert_merged(name: &str, points: &[Point2], constraints: &[Constraint], tri: &Triangulation) {
    let kept = expected_kept(points);
    let report: Vec<(u32, u32)> = tri
        .merged_points()
        .iter()
        .map(|m| (m.input, m.kept))
        .collect();
    let expected: Vec<(u32, u32)> = (0..points.len() as u32)
        .filter(|&k| kept[k as usize] != k)
        .map(|k| (k, kept[k as usize]))
        .collect();
    assert_eq!(report, expected, "{name}: merge report");
    for &v in tri.triangles() {
        assert_eq!(kept[v as usize], v, "{name}: a merged index is a corner");
    }
    let mut carried: Vec<Constraint> = constraints
        .iter()
        .map(|c| Constraint::new(kept[c.a as usize], kept[c.b as usize]))
        .collect();
    carried.sort_unstable();
    carried.dedup();
    assert_eq!(tri.constraints(), &carried[..], "{name}: constraints");
    for c in &carried {
        assert!(
            tri.triangles()
                .chunks_exact(3)
                .any(|t| t.contains(&c.a) && t.contains(&c.b)),
            "{name}: constraint {c:?} is missing"
        );
    }
    assert_covers_hull(name, points, tri);
}

fn check(name: &str, points: &[Point2], constraints: &[Constraint]) -> Triangulation {
    let tri = triangulate(points, constraints).unwrap_or_else(|e| panic!("{name}: {e}"));
    assert_merged(name, points, constraints, &tri);
    tri
}

/// `points` with `copies` of its points, picked by `r`, appended and the
/// whole shuffled; returns the points and where each original went.
fn with_copies(points: &[Point2], copies: usize, r: &mut Rng) -> Vec<Point2> {
    let mut out = points.to_vec();
    for _ in 0..copies {
        let k = (r.unit() * points.len() as f64) as usize;
        out.push(points[k]);
    }
    for k in (1..out.len()).rev() {
        let j = (r.unit() * (k + 1) as f64) as usize;
        out.swap(k, j);
    }
    out
}

#[test]
fn grids_with_duplicates_triangulate_validly() {
    let mut r = Rng(0x3001);
    for side in [2usize, 3, 5, 12] {
        let grid: Vec<Point2> = (0..side * side)
            .map(|k| p((k % side) as f64, (k / side) as f64))
            .collect();
        // Every point twice, in input order: the 3x3 grid of #300.
        let mut doubled = grid.clone();
        doubled.extend_from_slice(&grid);
        let tri = check(&format!("doubled grid {side}"), &doubled, &[]);
        assert_eq!(tri.merged_points().len(), grid.len());
        // Some points many times, shuffled, so duplicates come first too.
        let copies = with_copies(&grid, 3 * grid.len(), &mut r);
        check(&format!("shuffled grid {side}"), &copies, &[]);
    }
    // The case of #300: a 3x3 grid with its first row repeated.
    let mut points: Vec<Point2> = (0..9)
        .map(|k| p(f64::from(k % 3), f64::from(k / 3)))
        .collect();
    points.extend_from_slice(&points.clone()[..3]);
    let tri = check("grid of #300", &points, &[]);
    assert_eq!(tri.triangle_count(), 8);
}

#[test]
fn random_sets_with_duplicates_triangulate_validly() {
    let mut r = Rng(0x3002);
    for (k, count) in [3usize, 4, 10, 100, 1000, 5000].into_iter().enumerate() {
        let points = uniform(count, 0x3100 + k as u64);
        let copies = with_copies(&points, count / 2 + 1, &mut r);
        check(&format!("{count} random"), &copies, &[]);
    }
    // Only three distinct points, many times over.
    let three = [p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0)];
    let many = with_copies(&three, 30, &mut r);
    let tri = check("three distinct", &many, &[]);
    assert_eq!(tri.triangle_count(), 1);
    // Negative zero is the same position as zero.
    let signed = [
        p(0.0, 0.0),
        p(1.0, 0.0),
        p(0.0, 1.0),
        p(-0.0, 0.0),
        p(1.0, -0.0),
    ];
    let tri = check("signed zeros", &signed, &[]);
    assert_eq!(tri.merged_points().len(), 2);
}

#[test]
fn every_point_at_one_position_is_refused_as_collinear() {
    let points = [p(2.0, 3.0); 5];
    assert_eq!(
        triangulate(&points, &[]).unwrap_err(),
        TriangulationError::AllPointsCollinear
    );
}

/// A closed polyline through `corners`, sampled every `spacing`, as
/// points numbered from `base` and the constraints between them.
fn ring(corners: &[(f64, f64)], spacing: f64, base: u32) -> (Vec<Point2>, Vec<Constraint>) {
    let mut points = Vec::new();
    for k in 0..corners.len() {
        let (a, b) = (corners[k], corners[(k + 1) % corners.len()]);
        let n = ((b.0 - a.0).hypot(b.1 - a.1) / spacing).ceil().max(1.0) as usize;
        for j in 0..n {
            let t = j as f64 / n as f64;
            points.push(p(a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t));
        }
    }
    let n = points.len() as u32;
    let constraints = (0..n)
        .map(|i| Constraint::new(base + i, base + (i + 1) % n))
        .collect();
    (points, constraints)
}

#[test]
fn duplicates_in_the_outline_holes_and_constraints_are_merged() {
    let mut r = Rng(0x3003);
    let (mut points, mut constraints) =
        ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 8.0), (0.0, 8.0)], 0.5, 0);
    let outline = points.len();

    // The outline's points again, as free points.
    let free: Vec<Point2> = points.iter().step_by(3).copied().collect();
    points.extend(free);

    // A hole whose ring is listed twice: once as free points, once with
    // its constraints.
    let base = points.len() as u32;
    let (hole, _) = ring(&[(3.0, 3.0), (3.0, 5.0), (6.0, 5.0), (6.0, 3.0)], 0.4, 0);
    points.extend(hole.iter().copied());
    let (copy, around) = ring(
        &[(3.0, 3.0), (3.0, 5.0), (6.0, 5.0), (6.0, 3.0)],
        0.4,
        base + hole.len() as u32,
    );
    points.extend(copy);
    constraints.extend(around);

    // A second region touching the outline at a corner, its ring through a
    // copy of that corner: the outline pinches there once merged.
    let (touch, around) = ring(
        &[(10.0, 8.0), (13.0, 8.5), (12.0, 11.0)],
        0.5,
        points.len() as u32,
    );
    points.extend(touch);
    constraints.extend(around);

    // Chords from the outline to the hole, through duplicate indices on
    // both ends, and the same chord again through the originals.
    let a = points.len() as u32;
    points.push(points[outline / 8]);
    points.push(hole[0]);
    constraints.push(Constraint::new(a, a + 1));
    constraints.push(Constraint::new((outline / 8) as u32, base));

    // Random interior points, some repeated.
    let cloud: Vec<Point2> = (0..150)
        .map(|_| p(0.2 + 9.6 * r.unit(), 0.2 + 7.6 * r.unit()))
        .filter(|q| !(q.x > 2.9 && q.x < 6.1 && q.y > 2.9 && q.y < 5.1))
        .collect();
    let cloud = with_copies(&cloud, 40, &mut r);
    points.extend(cloud);

    let tri = check("outline, hole, constraints", &points, &constraints);
    assert!(tri.merged_points().len() > 40 + hole.len());
}

#[test]
fn a_constraint_between_coincident_points_is_refused_by_name() {
    // A ring that repeats its first point at the end: the closing
    // constraint joins the point to its own copy.
    let (mut points, mut constraints) = ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0)], 1.0, 0);
    let last = points.len() as u32;
    points.push(points[0]);
    constraints.pop();
    constraints.push(Constraint::new(last - 1, last));
    constraints.push(Constraint::new(last, 0));
    let err = triangulate(&points, &constraints).unwrap_err();
    assert_eq!(
        err,
        TriangulationError::CoincidentConstraint { a: 0, b: last }
    );
    assert_eq!(
        err.to_string(),
        format!("constraint (0, {last}) joins two coincident points")
    );

    // Two free points at one position, joined by a constraint.
    let points = [
        p(0.0, 0.0),
        p(2.0, 0.0),
        p(1.0, 1.0),
        p(1.0, 0.5),
        p(1.0, 0.5),
    ];
    assert_eq!(
        triangulate(&points, &[Constraint::new(4, 3)]).unwrap_err(),
        TriangulationError::CoincidentConstraint { a: 3, b: 4 }
    );
    // Without that constraint, the same points triangulate.
    check("free pair", &points, &[Constraint::new(0, 4)]);
}
