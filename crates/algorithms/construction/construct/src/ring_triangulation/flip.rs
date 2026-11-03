//! Flipping slivers out of a clipped triangulation (#415).
//!
//! Ear clipping ([`super::clip`]) cuts a vertex a few ulps off its
//! neighbours' line as a sliver ear, though the other diagonal of the
//! quadrilateral it leaves would not. So after clipping, each triangle that is a sliver (doubled area within `flat^2`) has its
//! longest edge flipped when the triangle across it and itself form a
//! strictly convex quadrilateral (both new triangles strictly
//! counter-clockwise by exact `orient2d`) and neither new triangle is a
//! sliver. Triangles are visited in index order and passes repeat until one
//! changes nothing, at most [`MAX_PASSES`]. Every flip removes a sliver and
//! makes none, so the sliver count falls and the pass ends; the bound only
//! guards the loop. Only edges shared by two triangles are flipped, so the
//! boundary and the bridges' two copies are never moved, and the caller
//! certifies the result again.

use std::collections::HashMap;

use axiolid_core::Point2;

use super::orient;

/// Passes over the triangles before the flipping gives up.
const MAX_PASSES: usize = 64;

/// Whether triangle `a p b` is a sliver: its doubled area is within
/// `flat^2`, the threshold a mesh audit reads as degenerate (#415). A `flat`
/// of `0` names no sliver: an exactly proved triangle is never refused for
/// the rounding of this estimate.
pub(crate) fn is_sliver(a: Point2, p: Point2, b: Point2, flat: f64) -> bool {
    flat > 0.0 && ((p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x)).abs() <= flat * flat
}

/// Whether triangle `a p b` is flat to rounding: its doubled area is within
/// `1e-9` of its longest side squared, the relative band `planar` reads a
/// run of corners as collinear in. Placement can round such a triangle to
/// zero area or invert it, whatever the tolerance (#415).
pub(crate) fn is_flat(a: Point2, p: Point2, b: Point2) -> bool {
    let double_area = ((p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x)).abs();
    let side = |u: Point2, v: Point2| (v.x - u.x).hypot(v.y - u.y);
    let longest = side(a, p).max(side(p, b)).max(side(b, a));
    double_area <= 1e-9 * longest * longest
}

/// Flip every sliver of `triangles` away that a convex neighbour allows.
pub(super) fn unsliver(points: &[Point2], triangles: &mut [[u32; 3]], flat: f64) {
    if flat <= 0.0 {
        return;
    }
    let at = |i: u32| points[i as usize];
    // The triangle holding each directed edge.
    let mut across: HashMap<(u32, u32), usize> = HashMap::with_capacity(3 * triangles.len());
    for (k, t) in triangles.iter().enumerate() {
        for e in 0..3 {
            across.insert((t[e], t[(e + 1) % 3]), k);
        }
    }
    for _ in 0..MAX_PASSES {
        let mut changed = false;
        for k in 0..triangles.len() {
            let t = triangles[k];
            if !is_sliver(at(t[0]), at(t[1]), at(t[2]), flat) {
                continue;
            }
            // The longest edge `u -> v`, the first of equals.
            let length = |e: usize| {
                let (u, v) = (at(t[e]), at(t[(e + 1) % 3]));
                (v.x - u.x).hypot(v.y - u.y)
            };
            let longest = (1..3).fold(0, |best, e| if length(e) > length(best) { e } else { best });
            let (u, v, w) = (t[longest], t[(longest + 1) % 3], t[(longest + 2) % 3]);
            let Some(&other) = across.get(&(v, u)) else {
                continue;
            };
            let s = triangles[other];
            let Some(x) = s.iter().copied().find(|&c| c != u && c != v) else {
                continue;
            };
            // The quadrilateral `u, x, v, w` is cut along `x w` instead.
            if orient(at(u), at(x), at(w)) <= 0 || orient(at(x), at(v), at(w)) <= 0 {
                continue;
            }
            if is_sliver(at(u), at(x), at(w), flat) || is_sliver(at(x), at(v), at(w), flat) {
                continue;
            }
            for t in [t, s] {
                for e in 0..3 {
                    across.remove(&(t[e], t[(e + 1) % 3]));
                }
            }
            triangles[k] = [u, x, w];
            triangles[other] = [x, v, w];
            for (index, t) in [(k, triangles[k]), (other, triangles[other])] {
                for e in 0..3 {
                    across.insert((t[e], t[(e + 1) % 3]), index);
                }
            }
            changed = true;
        }
        if !changed {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sliver `u v w` over the triangle `v u x` across its long edge,
    /// for a `flat` of 1 (doubled area within 1 is a sliver).
    fn pair(w: [f64; 2], x: [f64; 2]) -> (Vec<Point2>, Vec<[u32; 3]>) {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(10.0, 0.0),
            Point2::new(w[0], w[1]),
            Point2::new(x[0], x[1]),
        ];
        (points, vec![[0, 1, 2], [1, 0, 3]])
    }

    #[test]
    fn a_sliver_is_flipped_across_its_long_edge() {
        let (points, mut triangles) = pair([5.0, 0.08], [5.0, -3.0]);
        unsliver(&points, &mut triangles, 1.0);
        assert_eq!(triangles, vec![[0, 3, 2], [3, 1, 2]]);
    }

    #[test]
    fn a_flip_that_would_make_a_sliver_is_left_alone() {
        let (points, mut triangles) = pair([5.0, 0.08], [5.0, -0.02]);
        unsliver(&points, &mut triangles, 1.0);
        assert_eq!(triangles, vec![[0, 1, 2], [1, 0, 3]]);
    }

    #[test]
    fn a_flip_across_a_reflex_quadrilateral_is_left_alone() {
        // `x` lies far beyond `v`: the quadrilateral is reflex at `v`, and
        // `x w` would leave the polygon.
        let (points, mut triangles) = pair([5.0, 0.08], [100.0, -1.0]);
        unsliver(&points, &mut triangles, 1.0);
        assert_eq!(triangles, vec![[0, 1, 2], [1, 0, 3]]);
    }

    #[test]
    fn no_flat_length_flips_nothing() {
        let (points, mut triangles) = pair([5.0, 0.08], [5.0, -3.0]);
        unsliver(&points, &mut triangles, 0.0);
        assert_eq!(triangles, vec![[0, 1, 2], [1, 0, 3]]);
    }
}
