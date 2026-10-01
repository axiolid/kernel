//! Convex pieces of a polygon with holes, for Minkowski sums with a
//! non-convex polygon (#145).
//!
//! The holes are bridged into the outer ring, the bridged ring is cut into
//! triangles by ear clipping, and triangles are merged across a diagonal
//! while the union stays convex (Hertel and Mehlhorn: at most four times the
//! fewest convex pieces). Every decision is an exact orientation, and every
//! piece's vertices are the polygon's own.
//!
//! The triangulation is then certified rather than trusted: the triangles
//! are counter-clockwise and the sum of their boundaries, as directed edges
//! with opposite pairs cancelling, is the polygon's boundary. The winding
//! number of that sum at a point counts the triangles holding it and equals
//! the polygon's own winding number, one inside and zero outside, so the
//! triangles tile the polygon exactly. Any step that cannot proceed, or a
//! certificate that fails, answers `None`, and the caller takes a route that
//! needs no decomposition.

use std::collections::HashMap;

use axiolid_core::Point2;
use axiolid_guarantees::Sign;

use crate::minkowski::orient;

/// Convex pieces, each counter-clockwise, whose union is the polygon with
/// outer ring `outer` and holes `holes` (either orientation); `None` when
/// no certified decomposition was found.
pub(crate) fn convex_parts(outer: &[Point2], holes: &[Vec<Point2>]) -> Option<Vec<Vec<Point2>>> {
    let outer = oriented(outer, true)?;
    let holes: Vec<Vec<Point2>> = holes
        .iter()
        .map(|h| oriented(h, false))
        .collect::<Option<_>>()?;
    let bridged = bridge(outer.clone(), &holes)?;
    // Where clipping starts decides the triangles, and so how many pieces
    // the merge leaves: try several starts and keep the fewest pieces.
    let n = bridged.len();
    let tries = if n <= 64 { n } else { 8 };
    let mut best: Option<Vec<Vec<Point2>>> = None;
    for t in 0..tries {
        let Some(triangles) = ear_clip(&bridged, t * n / tries) else {
            continue;
        };
        if !tiles(&triangles, &outer, &holes) {
            continue;
        }
        let pieces = merge(triangles);
        if best.as_ref().is_none_or(|b| pieces.len() < b.len()) {
            best = Some(pieces);
        }
    }
    best
}

/// The ring counter-clockwise (`ccw`) or clockwise, decided exactly at its
/// lowest, then leftmost, vertex.
fn oriented(ring: &[Point2], ccw: bool) -> Option<Vec<Point2>> {
    let n = ring.len();
    if n < 3 {
        return None;
    }
    let i = (0..n).min_by(|&a, &b| {
        ring[a]
            .y
            .total_cmp(&ring[b].y)
            .then(ring[a].x.total_cmp(&ring[b].x))
    })?;
    let turn = orient(ring[(i + n - 1) % n], ring[i], ring[(i + 1) % n]);
    if turn == Sign::Zero {
        return None;
    }
    let mut out = ring.to_vec();
    if (turn == Sign::Positive) != ccw {
        out.reverse();
    }
    Some(out)
}

/// Whether `p` lies on the closed segment `a b`, given that the three are
/// collinear.
fn within_box(a: Point2, b: Point2, p: Point2) -> bool {
    p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
}

/// Whether the closed segments `a b` and `c d` share a point other than an
/// endpoint common to both (by coordinates); a collinear overlap beyond a
/// common endpoint counts.
fn blocks(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let (o1, o2) = (orient(a, b, c), orient(a, b, d));
    let (o3, o4) = (orient(c, d, a), orient(c, d, b));
    let shared = |p: Point2| p == a || p == b;
    if shared(c) || shared(d) {
        // Incident: only a collinear overlap past the common point blocks.
        let (common, other) = if shared(c) { (c, d) } else { (d, c) };
        if shared(other) {
            return true;
        }
        let far = if common == a { b } else { a };
        return orient(common, far, other) == Sign::Zero
            && (within_box(common, far, other) || within_box(common, other, far));
    }
    if o1 != o2
        && o1 != Sign::Zero
        && o2 != Sign::Zero
        && o3 != o4
        && o3 != Sign::Zero
        && o4 != Sign::Zero
    {
        return true;
    }
    (o1 == Sign::Zero && within_box(a, b, c))
        || (o2 == Sign::Zero && within_box(a, b, d))
        || (o3 == Sign::Zero && within_box(c, d, a))
        || (o4 == Sign::Zero && within_box(c, d, b))
}

/// Whether the direction from `v` to `h` lies strictly inside the interior
/// wedge of a counter-clockwise ring at `v`, between `prev` and `next`.
fn in_cone(prev: Point2, v: Point2, next: Point2, h: Point2) -> bool {
    match orient(prev, v, next) {
        Sign::Positive => {
            orient(v, next, h) == Sign::Positive && orient(v, h, prev) == Sign::Positive
        }
        Sign::Negative => {
            !(orient(v, prev, h) != Sign::Negative && orient(v, h, next) != Sign::Negative)
        }
        _ => orient(prev, next, h) == Sign::Positive,
    }
}

/// The outer ring (counter-clockwise) with every hole (clockwise) spliced
/// in along a bridge from the hole's rightmost vertex to a vertex it sees
/// to its right: a weakly simple ring whose bridges are walked both ways.
fn bridge(mut ring: Vec<Point2>, holes: &[Vec<Point2>]) -> Option<Vec<Point2>> {
    // Rightmost hole first, so a later hole's bridge never needs a vertex
    // of a hole not yet spliced in.
    let right = |h: &Vec<Point2>| {
        (0..h.len())
            .max_by(|&a, &b| h[a].x.total_cmp(&h[b].x).then(h[a].y.total_cmp(&h[b].y)))
            .unwrap_or(0)
    };
    let mut order: Vec<usize> = (0..holes.len()).collect();
    order.sort_by(|&a, &b| {
        let (pa, pb) = (holes[a][right(&holes[a])], holes[b][right(&holes[b])]);
        pb.x.total_cmp(&pa.x).then(pb.y.total_cmp(&pa.y))
    });
    for (rank, &hi) in order.iter().enumerate() {
        let hole = &holes[hi];
        let at = right(hole);
        let h = hole[at];
        // Every edge a bridge must keep clear of: the ring so far, and
        // every hole not yet spliced in (this one included).
        let mut walls: Vec<(Point2, Point2)> = Vec::new();
        let ring_edges = |r: &[Point2]| {
            (0..r.len())
                .map(|i| (r[i], r[(i + 1) % r.len()]))
                .collect::<Vec<_>>()
        };
        walls.extend(ring_edges(&ring));
        for &other in &order[rank..] {
            walls.extend(ring_edges(&holes[other]));
        }
        let n = ring.len();
        let mut best: Option<(f64, usize)> = None;
        for j in 0..n {
            let v = ring[j];
            if v.x <= h.x {
                continue;
            }
            let d = (v - h).length_squared();
            if best.is_some_and(|(b, _)| b <= d) {
                continue;
            }
            if !in_cone(ring[(j + n - 1) % n], v, ring[(j + 1) % n], h) {
                continue;
            }
            if walls.iter().any(|&(p, q)| blocks(h, v, p, q)) {
                continue;
            }
            best = Some((d, j));
        }
        let (_, j) = best?;
        let v = ring[j];
        let mut spliced = Vec::with_capacity(n + hole.len() + 2);
        spliced.extend_from_slice(&ring[..=j]);
        spliced.extend((0..hole.len()).map(|k| hole[(at + k) % hole.len()]));
        spliced.push(h);
        spliced.push(v);
        spliced.extend_from_slice(&ring[j + 1..]);
        ring = spliced;
    }
    Some(ring)
}

/// Triangles of a (weakly) simple counter-clockwise ring by ear clipping:
/// a vertex is clipped when it turns left and no other vertex, by
/// coordinates distinct from the three corners, lies in the closed
/// triangle.
fn ear_clip(ring: &[Point2], start: usize) -> Option<Vec<[Point2; 3]>> {
    let mut live: Vec<usize> = (0..ring.len()).collect();
    let mut triangles = Vec::with_capacity(ring.len().saturating_sub(2));
    let mut at = start;
    let mut tried = 0;
    while live.len() > 3 {
        let n = live.len();
        let (ia, ib, ic) = ((at + n - 1) % n, at % n, (at + 1) % n);
        let (a, b, c) = (ring[live[ia]], ring[live[ib]], ring[live[ic]]);
        let ear = orient(a, b, c) == Sign::Positive
            && !live.iter().any(|&k| {
                let p = ring[k];
                p != a
                    && p != b
                    && p != c
                    && orient(a, b, p) != Sign::Negative
                    && orient(b, c, p) != Sign::Negative
                    && orient(c, a, p) != Sign::Negative
            });
        if ear {
            triangles.push([a, b, c]);
            live.remove(ib);
            at = if ib == 0 { 0 } else { ib - 1 };
            tried = 0;
        } else {
            at = (at + 1) % n;
            tried += 1;
            if tried > n {
                return None;
            }
        }
    }
    let (a, b, c) = (ring[live[0]], ring[live[1]], ring[live[2]]);
    if orient(a, b, c) != Sign::Positive {
        return None;
    }
    triangles.push([a, b, c]);
    Some(triangles)
}

/// A directed edge's key: its endpoints' bits (signed zeros merged), and
/// whether it runs against the key's order.
fn key(p: Point2, q: Point2) -> ((u64, u64, u64, u64), i32) {
    let bits = |r: Point2| ((r.x + 0.0).to_bits(), (r.y + 0.0).to_bits());
    let (pb, qb) = (bits(p), bits(q));
    if pb <= qb {
        ((pb.0, pb.1, qb.0, qb.1), 1)
    } else {
        ((qb.0, qb.1, pb.0, pb.1), -1)
    }
}

/// The certificate: the triangles' boundaries sum to the polygon's.
fn tiles(triangles: &[[Point2; 3]], outer: &[Point2], holes: &[Vec<Point2>]) -> bool {
    let mut chain: HashMap<(u64, u64, u64, u64), i32> = HashMap::new();
    let mut add = |p: Point2, q: Point2, s: i32| {
        let (k, d) = key(p, q);
        *chain.entry(k).or_insert(0) += s * d;
    };
    for t in triangles {
        if orient(t[0], t[1], t[2]) != Sign::Positive {
            return false;
        }
        for i in 0..3 {
            add(t[i], t[(i + 1) % 3], 1);
        }
    }
    for ring in std::iter::once(outer).chain(holes.iter().map(Vec::as_slice)) {
        for i in 0..ring.len() {
            add(ring[i], ring[(i + 1) % ring.len()], -1);
        }
    }
    chain.values().all(|&c| c == 0)
}

/// Merge triangles across shared diagonals while the union stays convex.
fn merge(triangles: Vec<[Point2; 3]>) -> Vec<Vec<Point2>> {
    let mut pieces: Vec<Option<Vec<Point2>>> = triangles.iter().map(|t| Some(t.to_vec())).collect();
    // Which piece holds each directed edge.
    let bits = |p: Point2, q: Point2| {
        let b = |r: Point2| ((r.x + 0.0).to_bits(), (r.y + 0.0).to_bits());
        (b(p), b(q))
    };
    let mut owner = HashMap::new();
    for (i, t) in triangles.iter().enumerate() {
        for e in 0..3 {
            owner.insert(bits(t[e], t[(e + 1) % 3]), i);
        }
    }
    for t in &triangles {
        for e in 0..3 {
            let (u, v) = (t[e], t[(e + 1) % 3]);
            let (Some(&i), Some(&j)) = (owner.get(&bits(u, v)), owner.get(&bits(v, u))) else {
                continue;
            };
            if i == j {
                continue;
            }
            let (Some(pi), Some(pj)) = (pieces[i].as_ref(), pieces[j].as_ref()) else {
                continue;
            };
            let Some(merged) = join(pi, pj, u, v) else {
                continue;
            };
            let (keep, gone) = (i.min(j), i.max(j));
            for k in 0..merged.len() {
                owner.insert(bits(merged[k], merged[(k + 1) % merged.len()]), keep);
            }
            owner.remove(&bits(u, v));
            owner.remove(&bits(v, u));
            pieces[keep] = Some(merged);
            pieces[gone] = None;
        }
    }
    pieces.into_iter().flatten().collect()
}

/// The union of convex `p` (holding the edge `u -> v`) and `q` (holding
/// `v -> u`), when it is convex.
fn join(p: &[Point2], q: &[Point2], u: Point2, v: Point2) -> Option<Vec<Point2>> {
    let iu = (0..p.len()).find(|&k| p[k] == u && p[(k + 1) % p.len()] == v)?;
    let jv = (0..q.len()).find(|&k| q[k] == v && q[(k + 1) % q.len()] == u)?;
    // p from v round to u, then q after u round to before v.
    let mut out: Vec<Point2> = (0..p.len()).map(|k| p[(iu + 1 + k) % p.len()]).collect();
    out.extend((0..q.len() - 2).map(|k| q[(jv + 2 + k) % q.len()]));
    let n = out.len();
    // Only the two joined corners changed: u at p.len() - 1, v at 0.
    let convex_at =
        |k: usize| orient(out[(k + n - 1) % n], out[k], out[(k + 1) % n]) != Sign::Negative;
    (convex_at(p.len() - 1) && convex_at(0)).then_some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pts(v: &[(f64, f64)]) -> Vec<Point2> {
        v.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    fn area(p: &[Point2]) -> f64 {
        (0..p.len())
            .map(|i| p[i].perp_dot(p[(i + 1) % p.len()]))
            .sum::<f64>()
            * 0.5
    }

    fn convex(p: &[Point2]) -> bool {
        let n = p.len();
        (0..n).all(|i| orient(p[i], p[(i + 1) % n], p[(i + 2) % n]) != Sign::Negative)
    }

    #[test]
    fn a_plus_is_three_pieces_from_any_start() {
        let plus = pts(&[
            (-1., -3.),
            (1., -3.),
            (1., -1.),
            (3., -1.),
            (3., 1.),
            (1., 1.),
            (1., 3.),
            (-1., 3.),
            (-1., 1.),
            (-3., 1.),
            (-3., -1.),
            (-1., -1.),
        ]);
        for r in 0..plus.len() {
            let mut p = plus.clone();
            p.rotate_left(r);
            let parts = convex_parts(&p, &[]).unwrap();
            assert_eq!(parts.len(), 3, "start {r}");
        }
    }

    #[test]
    fn an_l_shape_is_two_convex_pieces() {
        let l = pts(&[(0., 0.), (2., 0.), (2., 1.), (1., 1.), (1., 2.), (0., 2.)]);
        let parts = convex_parts(&l, &[]).unwrap();
        assert_eq!(parts.len(), 2);
        assert!(parts.iter().all(|p| convex(p)));
        assert_eq!(parts.iter().map(|p| area(p)).sum::<f64>(), 3.0);
    }

    #[test]
    fn a_square_with_a_square_hole_is_tiled() {
        let outer = pts(&[(0., 0.), (4., 0.), (4., 4.), (0., 4.)]);
        let hole = pts(&[(1., 1.), (3., 1.), (3., 3.), (1., 3.)]);
        let parts = convex_parts(&outer, &[hole]).unwrap();
        assert!(parts.iter().all(|p| convex(p)));
        assert_eq!(parts.iter().map(|p| area(p)).sum::<f64>(), 12.0);
        assert!(parts.len() <= 8, "{}", parts.len());
    }

    #[test]
    fn a_comb_with_two_holes_is_tiled() {
        let outer = pts(&[
            (0., 0.),
            (9., 0.),
            (9., 5.),
            (7., 5.),
            (7., 2.),
            (5., 2.),
            (5., 5.),
            (0., 5.),
        ]);
        let holes = vec![
            pts(&[(1., 1.), (2., 1.), (2., 2.), (1., 2.)]),
            pts(&[(1., 3.), (4., 3.), (2.5, 4.)]),
        ];
        let parts = convex_parts(&outer, &holes).unwrap();
        assert!(parts.iter().all(|p| convex(p)));
        let want = 9. * 5. - 2. * 3. - 1. - 1.5;
        assert_eq!(parts.iter().map(|p| area(p)).sum::<f64>(), want);
    }

    #[test]
    fn a_bridge_goes_round_a_hole_in_the_way() {
        // The nearest vertex right of the square hole is the tip of a
        // spike, behind the long thin triangular hole: the bridge must
        // pass that hole by.
        let outer = pts(&[
            (0., 0.),
            (10., 0.),
            (10., 5.4),
            (3., 5.5),
            (10., 5.6),
            (10., 10.),
            (0., 10.),
        ]);
        let holes = vec![
            pts(&[(1., 4.5), (2., 4.5), (2., 5.5), (1., 5.5)]),
            pts(&[(2.5, 1.), (2.7, 1.), (2.6, 9.)]),
        ];
        let parts = convex_parts(&outer, &holes).expect("a certified cut");
        assert!(parts.iter().all(|p| convex(p)));
        let total: f64 = parts.iter().map(|p| area(p)).sum();
        assert!((total - (100. - 0.7 - 1. - 0.8)).abs() < 1e-9, "{total}");
    }

    #[test]
    fn a_false_tiling_fails_its_certificate() {
        let square = pts(&[(0., 0.), (1., 0.), (1., 1.), (0., 1.)]);
        let t = [square[0], square[1], square[2]];
        assert!(!tiles(&[t], &square, &[]));
        assert!(!tiles(
            &[t, t, [square[0], square[2], square[3]]],
            &square,
            &[]
        ));
        assert!(tiles(&[t, [square[0], square[2], square[3]]], &square, &[]));
    }
}
