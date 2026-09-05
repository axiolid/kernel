//! KNN and radius queries over point sets.
//!
//! The central test is agreement with brute force: an acceleration
//! structure is only correct if it returns exactly what an exhaustive scan
//! would, so every query test below has a brute-force oracle rather than a
//! hand-written expected list.

use axiolid_core::{Point3, Scalar};
use axiolid_spatial::{PointHit, PointIndex, PointQueryError};

/// Deterministic pseudo-random points, so a failure is reproducible.
fn scattered(count: usize) -> Vec<Point3> {
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as Scalar / (1u64 << 53) as Scalar
    };
    (0..count)
        .map(|_| Point3::new(next() * 10.0, next() * 4.0, next() * 7.0))
        .collect()
}

/// What an exhaustive scan would return within a radius.
fn brute_radius(points: &[Point3], query: Point3, radius: Scalar) -> Vec<usize> {
    let mut hits: Vec<usize> = points
        .iter()
        .enumerate()
        .filter(|(_, p)| p.is_finite() && (**p - query).length() <= radius)
        .map(|(i, _)| i)
        .collect();
    hits.sort_unstable();
    hits
}

/// What an exhaustive scan would return as the k nearest.
fn brute_nearest(points: &[Point3], query: Point3, k: usize) -> Vec<usize> {
    let mut all: Vec<(Scalar, usize)> = points
        .iter()
        .enumerate()
        .filter(|(_, p)| p.is_finite())
        .map(|(i, p)| ((*p - query).length(), i))
        .collect();
    all.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(a.1.cmp(&b.1))
    });
    all.into_iter().take(k).map(|(_, i)| i).collect()
}

fn indices(hits: &[PointHit]) -> Vec<usize> {
    hits.iter().map(|h| h.index).collect()
}

/// Radius search must agree with an exhaustive scan at every scale: too
/// small a search box silently loses far points, too large only costs time.
#[test]
fn radius_queries_agree_with_brute_force() {
    let points = scattered(400);
    let index = PointIndex::build(&points);
    let mut out = Vec::new();

    for query in [
        Point3::new(5.0, 2.0, 3.5),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(9.9, 3.9, 6.9),
        Point3::new(-5.0, -5.0, -5.0),
    ] {
        for radius in [0.0, 0.25, 1.0, 3.0, 25.0] {
            index.radius_into(query, radius, &mut out).expect("query");
            let mut got = indices(&out);
            got.sort_unstable();
            assert_eq!(
                got,
                brute_radius(&points, query, radius),
                "radius {radius} at {query:?} disagrees with brute force"
            );
        }
    }
}

/// KNN must agree with an exhaustive scan, including the ordering.
#[test]
fn nearest_queries_agree_with_brute_force() {
    let points = scattered(300);
    let index = PointIndex::build(&points);
    let mut out = Vec::new();

    for query in [
        Point3::new(5.0, 2.0, 3.5),
        Point3::new(-2.0, 1.0, 0.5),
        Point3::new(12.0, 6.0, 9.0),
    ] {
        for k in [1usize, 3, 12, 50] {
            index.nearest_into(query, k, &mut out).expect("query");
            assert_eq!(
                indices(&out),
                brute_nearest(&points, query, k),
                "k={k} at {query:?} disagrees with brute force"
            );
        }
    }
}

/// Distances are exact, not bounds. A caller may use them directly.
#[test]
fn reported_distances_are_exact() {
    let points = scattered(120);
    let index = PointIndex::build(&points);
    let query = Point3::new(3.0, 2.0, 1.0);
    let mut out = Vec::new();
    index.nearest_into(query, 10, &mut out).expect("query");

    for hit in &out {
        let actual = (points[hit.index] - query).length();
        assert!(
            (hit.distance - actual).abs() <= 1e-12,
            "reported {} but actual is {actual}",
            hit.distance
        );
    }
}

/// Degenerate inputs the issue calls out: empty, duplicate, coincident.
#[test]
fn an_empty_cloud_answers_rather_than_panicking() {
    let index = PointIndex::build(&[]);
    let mut out = Vec::new();
    index
        .radius_into(Point3::ZERO, 5.0, &mut out)
        .expect("empty radius query");
    assert!(out.is_empty());
    index
        .nearest_into(Point3::ZERO, 4, &mut out)
        .expect("empty knn query");
    assert!(out.is_empty());
    assert_eq!(index.nearest(Point3::ZERO).expect("nearest"), None);
}

/// Coincident points are all real answers: a duplicate is not a mistake to
/// be deduplicated behind the caller's back.
#[test]
fn duplicate_points_are_all_returned() {
    let points = vec![Point3::new(1.0, 1.0, 1.0); 6];
    let index = PointIndex::build(&points);
    let mut out = Vec::new();

    index
        .radius_into(Point3::new(1.0, 1.0, 1.0), 0.5, &mut out)
        .expect("query");
    assert_eq!(out.len(), 6, "every coincident point must be reported");
    assert!(out.iter().all(|h| h.distance == 0.0));

    index
        .nearest_into(Point3::new(1.0, 1.0, 1.0), 4, &mut out)
        .expect("query");
    assert_eq!(indices(&out), vec![0, 1, 2, 3], "ties break by index");
}

/// A cloud with no extent must not divide by zero when sizing cells.
#[test]
fn a_single_point_cloud_is_queryable() {
    let points = vec![Point3::new(2.0, 3.0, 4.0)];
    let index = PointIndex::build(&points);
    let hit = index.nearest(Point3::new(2.0, 3.0, 5.0)).expect("query");
    let hit = hit.expect("one point is nearest");
    assert_eq!(hit.index, 0);
    assert!((hit.distance - 1.0).abs() < 1e-12);
}

/// A planar sheet has zero extent on one axis; cell sizing must survive it.
#[test]
fn a_planar_cloud_is_queryable() {
    let points: Vec<Point3> = (0..40)
        .map(|i| Point3::new(i as Scalar * 0.5, (i % 7) as Scalar, 0.0))
        .collect();
    let index = PointIndex::build(&points);
    let mut out = Vec::new();
    index
        .radius_into(Point3::new(5.0, 3.0, 0.0), 2.0, &mut out)
        .expect("query");
    let mut got = indices(&out);
    got.sort_unstable();
    assert_eq!(got, brute_radius(&points, Point3::new(5.0, 3.0, 0.0), 2.0));
}

/// A non-finite position must not be indexed as if it sat at the origin,
/// which would corrupt every query near it. It is dropped and counted.
#[test]
fn non_finite_points_are_rejected_not_relocated() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(Scalar::NAN, 1.0, 1.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(Scalar::INFINITY, 0.0, 0.0),
    ];
    let index = PointIndex::build(&points);
    assert_eq!(index.rejected(), 2);

    let mut out = Vec::new();
    index
        .radius_into(Point3::ZERO, 0.5, &mut out)
        .expect("query");
    assert_eq!(
        indices(&out),
        vec![0],
        "only the finite point at the origin may be reported"
    );
}

/// Invalid queries are refused by name rather than returning an empty
/// result a caller would read as "nothing nearby".
#[test]
fn invalid_queries_are_refused_by_name() {
    let index = PointIndex::build(&scattered(20));
    let mut out = Vec::new();

    assert_eq!(
        index.radius_into(Point3::new(Scalar::NAN, 0.0, 0.0), 1.0, &mut out),
        Err(PointQueryError::NonFiniteQuery)
    );
    assert_eq!(
        index.nearest_into(Point3::new(0.0, Scalar::INFINITY, 0.0), 3, &mut out),
        Err(PointQueryError::NonFiniteQuery)
    );
    assert_eq!(
        index.radius_into(Point3::ZERO, -1.0, &mut out),
        Err(PointQueryError::InvalidRadius)
    );
    assert_eq!(
        index.radius_into(Point3::ZERO, Scalar::NAN, &mut out),
        Err(PointQueryError::InvalidRadius)
    );
}

/// Asking for more neighbours than exist returns everything, not an error:
/// that is a complete answer to the question asked.
#[test]
fn asking_for_more_neighbours_than_exist_returns_all_of_them() {
    let points = scattered(5);
    let index = PointIndex::build(&points);
    let mut out = Vec::new();
    index
        .nearest_into(Point3::ZERO, 50, &mut out)
        .expect("query");
    assert_eq!(out.len(), 5);
    assert_eq!(indices(&out), brute_nearest(&points, Point3::ZERO, 50));
}

/// The callback form allocates nothing per hit, so a caller can stream
/// results. Proven by counting through a closure with no vector at all.
#[test]
fn the_callback_form_needs_no_allocation() {
    let points = scattered(200);
    let index = PointIndex::build(&points);
    let query = Point3::new(5.0, 2.0, 3.0);

    let mut count = 0usize;
    let mut furthest: Scalar = 0.0;
    index
        .for_each_within(query, 2.0, |hit| {
            count += 1;
            furthest = furthest.max(hit.distance);
        })
        .expect("query");

    assert_eq!(count, brute_radius(&points, query, 2.0).len());
    assert!(furthest <= 2.0, "no hit may exceed the radius");
}

/// Repeating a query must give byte-identical results.
#[test]
fn queries_are_deterministic() {
    let points = scattered(250);
    let index = PointIndex::build(&points);
    let query = Point3::new(4.0, 1.5, 2.5);

    let mut first = Vec::new();
    let mut second = Vec::new();
    index.nearest_into(query, 15, &mut first).expect("query");
    index.nearest_into(query, 15, &mut second).expect("query");
    assert_eq!(first, second);

    let rebuilt = PointIndex::build(&points);
    let mut third = Vec::new();
    rebuilt.nearest_into(query, 15, &mut third).expect("query");
    assert_eq!(first, third, "a rebuilt index must answer identically");
}

/// A radius query is closed: a point exactly at distance r is inside.
/// Random points never land exactly on a boundary, so this needs points
/// placed there deliberately.
#[test]
fn a_point_exactly_on_the_radius_is_included() {
    let points = vec![
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
        Point3::new(0.0, 0.0, 2.0),
        Point3::new(2.000_000_1, 0.0, 0.0),
    ];
    let index = PointIndex::build(&points);
    let mut out = Vec::new();
    index
        .radius_into(Point3::ZERO, 2.0, &mut out)
        .expect("query");

    let mut got = indices(&out);
    got.sort_unstable();
    assert_eq!(
        got,
        vec![0, 1, 2],
        "points exactly at the radius are inside it; the one beyond is not"
    );
}
