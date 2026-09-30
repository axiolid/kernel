//! Output rings that the region's own validation accepts.
//!
//! The overlay backends and the offset may return rings with edges shorter
//! than the tolerance, rings that come back to a point they already passed
//! (two lobes touching, or a hole pinched off against the outer boundary),
//! and vertices within the tolerance of an edge they are not next to. Every
//! one of those is rejected by [`crate::validate_ring`], so feeding an
//! operation's result back into [`crate::Region::new`] or into another
//! operation failed (axioval, #188 follow-up).
//!
//! Settling repairs each ring by the same tests validation applies: short
//! edges are merged, a vertex touching another part of the ring is put on
//! it, and the ring is split wherever it passes a point twice. Each piece
//! keeps the side it bounds -- a piece wound like its ring bounds the same
//! kind of area, one wound the other way the opposite -- and holes are
//! given back to the outer ring that holds them. Every move is within the
//! tolerance.

use axiolid_core::{Point2, Tolerance};

use crate::{contains, cross, segments_intersect, signed, within_extent, Polygon, Ring};

/// Polygons whose every ring passes [`crate::validate_ring`].
pub(crate) fn settle(polygons: Vec<Polygon>, tolerance: Tolerance) -> Vec<Polygon> {
    let eps = tolerance.linear();
    let mut outers: Vec<Vec<Point2>> = Vec::new();
    let mut holes: Vec<Vec<Point2>> = Vec::new();
    for polygon in polygons {
        for (ring, hole) in std::iter::once((polygon.outer, false))
            .chain(polygon.holes.into_iter().map(|h| (h, true)))
        {
            let sense = signed(&ring) > 0.0;
            for piece in pieces(ring.points, eps) {
                let same = (area(&piece) > 0.0) == sense;
                if same != hole {
                    outers.push(piece);
                } else {
                    holes.push(piece);
                }
            }
        }
    }
    let mut out: Vec<Polygon> = outers
        .into_iter()
        .map(|points| Polygon {
            outer: Ring { points },
            holes: Vec::new(),
        })
        .collect();
    // Boxes of the outer rings: a point outside one, widened by the
    // tolerance, is not strictly inside that ring.
    let boxes: Vec<(Point2, Point2)> = out
        .iter()
        .map(|p| {
            p.outer.points.iter().fold(
                (
                    Point2::new(f64::INFINITY, f64::INFINITY),
                    Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
                ),
                |(lo, hi), q| (lo.min(*q), hi.max(*q)),
            )
        })
        .collect();
    for hole in holes {
        let ring = Ring { points: hole };
        // The smallest outer ring holding a vertex of the hole strictly.
        let owner = (0..out.len())
            .filter(|&i| {
                let (lo, hi) = boxes[i];
                ring.points.iter().any(|&p| {
                    p.x >= lo.x - eps
                        && p.x <= hi.x + eps
                        && p.y >= lo.y - eps
                        && p.y <= hi.y + eps
                        && strictly_inside(&out[i].outer, p, eps)
                })
            })
            .min_by(|&i, &j| {
                area(&out[i].outer.points)
                    .abs()
                    .total_cmp(&area(&out[j].outer.points).abs())
            });
        if let Some(i) = owner {
            out[i].holes.push(ring);
        }
    }
    out
}

fn area(points: &[Point2]) -> f64 {
    signed(&Ring {
        points: points.to_vec(),
    })
}

fn strictly_inside(ring: &Ring, p: Point2, eps: f64) -> bool {
    let n = ring.points.len();
    let near = (0..n).any(|i| {
        let (a, b) = (ring.points[i], ring.points[(i + 1) % n]);
        cross(a, b, p).abs() <= eps && within_extent(a, b, p, eps)
    });
    !near && contains(ring, p)
}

/// The ring repaired and cut into rings validation accepts.
fn pieces(points: Vec<Point2>, eps: f64) -> Vec<Vec<Point2>> {
    let mut done = Vec::new();
    let mut work = vec![points];
    let mut budget = 64 + 8 * work[0].len();
    while let Some(mut ring) = work.pop() {
        merge_short(&mut ring, eps);
        if ring.len() < 3 || area(&ring).abs() <= eps * eps {
            continue;
        }
        if let Some((i, j)) = split_point(&ring) {
            let (a, b) = cut(&ring, i, j);
            work.push(a);
            work.push(b);
            continue;
        }
        if budget == 0 {
            done.push(ring);
            continue;
        }
        budget -= 1;
        match touching(&ring, eps) {
            Some(fixed) => work.push(fixed),
            None => done.push(ring),
        }
    }
    done
}

/// Drop each vertex within `eps` of the one before it.
fn merge_short(ring: &mut Vec<Point2>, eps: f64) {
    let mut out: Vec<Point2> = Vec::with_capacity(ring.len());
    for &p in ring.iter() {
        if out.last().is_none_or(|q| (p - *q).length() > eps) {
            out.push(p);
        }
    }
    while out.len() > 1 && (out[0] - out[out.len() - 1]).length() <= eps {
        out.pop();
    }
    *ring = out;
}

/// Two positions of the ring at the same point.
fn split_point(ring: &[Point2]) -> Option<(usize, usize)> {
    // The first position of each point, by its bits (`-0.0` read as `0.0`,
    // as `==` does); the first pair in order is the least repeat's.
    let mut first: std::collections::HashMap<(u64, u64), usize> = std::collections::HashMap::new();
    let mut best: Option<(usize, usize)> = None;
    for (j, p) in ring.iter().enumerate() {
        let key = ((p.x + 0.0).to_bits(), (p.y + 0.0).to_bits());
        match first.get(&key) {
            Some(&i) => {
                if best.is_none_or(|b| (i, j) < b) {
                    best = Some((i, j));
                }
            }
            None => {
                first.insert(key, j);
            }
        }
    }
    best
}

fn cut(ring: &[Point2], i: usize, j: usize) -> (Vec<Point2>, Vec<Point2>) {
    let first = ring[i..j].to_vec();
    let mut second = ring[j..].to_vec();
    second.extend_from_slice(&ring[..i]);
    (first, second)
}

/// The first pair of edges validation would call intersecting, resolved:
/// a vertex near the other edge is put on it (onto its end if near that),
/// and edges crossing outright meet at their crossing. The ring then
/// passes that point twice and is split there next. `None` when no pair
/// touches.
/// The first pair of edges `(i, j)`, `i < j`, not neighbours, that meet
/// within `eps`. Only edges whose boxes, widened by `eps`, overlap can:
/// a sweep over the boxes finds them without asking every pair.
fn first_touch(ring: &[Point2], eps: f64) -> Option<(usize, usize)> {
    let n = ring.len();
    let edge = |i: usize| (ring[i], ring[(i + 1) % n]);
    let meets = |i: usize, j: usize| {
        let (first, second) = (i.min(j), i.max(j));
        if second == first + 1 || (first == 0 && second + 1 == n) {
            return false;
        }
        let ((a, b), (c, d)) = (edge(first), edge(second));
        segments_intersect(a, b, c, d, eps)
    };
    let low = |i: usize| edge(i).0.x.min(edge(i).1.x) - eps;
    let high = |i: usize| edge(i).0.x.max(edge(i).1.x) + eps;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| low(i).total_cmp(&low(j)));
    let mut active: Vec<usize> = Vec::new();
    let mut best: Option<(usize, usize)> = None;
    for &i in &order {
        active.retain(|&j| high(j) >= low(i));
        let (a, b) = edge(i);
        for &j in &active {
            let (c, d) = edge(j);
            if a.y.max(b.y) + eps < c.y.min(d.y) || c.y.max(d.y) + eps < a.y.min(b.y) {
                continue;
            }
            let pair = (i.min(j), i.max(j));
            if best.is_none_or(|b| pair < b) && meets(i, j) {
                best = Some(pair);
            }
        }
        active.push(i);
    }
    best
}

fn touching(ring: &[Point2], eps: f64) -> Option<Vec<Point2>> {
    let n = ring.len();
    {
        {
            let (i, j) = first_touch(ring, eps)?;
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            let (c, d) = (ring[j], ring[(j + 1) % n]);
            let on = |p: Point2, q: Point2, v: Point2| {
                cross(p, q, v).abs() <= eps && within_extent(p, q, v, eps)
            };
            let mut out = ring.to_vec();
            // A vertex of one edge on the other: the other edge passes it.
            // Snap to an end of that edge if the vertex is that close to it.
            let place = |out: &mut Vec<Point2>, v: Point2, at: usize, p: Point2, q: Point2| {
                if (v - p).length() <= eps {
                    out[at] = v;
                } else if (v - q).length() <= eps {
                    out[(at + 1) % n] = v;
                } else {
                    out.insert(at + 1, v);
                }
            };
            if on(c, d, a) {
                place(&mut out, a, j, c, d);
            } else if on(c, d, b) {
                place(&mut out, b, j, c, d);
            } else if on(a, b, c) {
                place(&mut out, c, i, a, b);
            } else if on(a, b, d) {
                place(&mut out, d, i, a, b);
            } else {
                // A proper crossing: both edges pass through it.
                let t = cross(c, d, a) / (cross(c, d, a) - cross(c, d, b));
                let x = a + (b - a) * t;
                out.insert(j + 1, x);
                out.insert(i + 1, x);
            }
            Some(out)
        }
    }
}
