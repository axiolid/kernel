//! Canonical form of output rings.

use axiolid_core::Point2;

/// Order of points by `x`, then `y` (total order on binary64).
pub fn point_order(left: &Point2, right: &Point2) -> std::cmp::Ordering {
    left.x
        .total_cmp(&right.x)
        .then_with(|| left.y.total_cmp(&right.y))
}

/// Lexicographic order of point sequences; a prefix sorts first.
pub fn sequence_order(left: &[Point2], right: &[Point2]) -> std::cmp::Ordering {
    left.iter()
        .zip(right)
        .map(|(a, b)| point_order(a, b))
        .find(|order| order.is_ne())
        .unwrap_or_else(|| left.len().cmp(&right.len()))
}

/// Index at which a ring starts to be canonical: at its lexicographically
/// smallest point; of several equal smallest points the rotation that orders
/// lowest wins.
pub fn canonical_start(points: &[Point2]) -> usize {
    let rotated = |start: usize| -> Vec<Point2> {
        points[start..]
            .iter()
            .chain(&points[..start])
            .copied()
            .collect()
    };
    let smallest = points.iter().min_by(|l, r| point_order(l, r)).copied();
    (0..points.len())
        .filter(|&start| Some(points[start]) == smallest)
        .min_by(|&left, &right| sequence_order(&rotated(left), &rotated(right)))
        .unwrap_or(0)
}

/// Rotate a ring to start at its lexicographically smallest point; of
/// several equal smallest points the rotation that orders lowest wins.
pub fn canonical_rotation(points: Vec<Point2>) -> Vec<Point2> {
    let start = canonical_start(&points);
    points[start..]
        .iter()
        .chain(&points[..start])
        .copied()
        .collect()
}
