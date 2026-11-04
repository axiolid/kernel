//! Re-triangulating a triangle whose sides carry extra corners.
//!
//! Index-only: the corners are caller-chosen indices, and collinearity is
//! decided by which side a corner lies on, never from coordinates. The caller
//! guarantees that a corner listed on a side lies on that side, strictly
//! between its ends and in order.

/// Triangulate a proper triangle `t` whose side `i` (from `t[i]` to
/// `t[i+1]`) carries the corners `sides[i]` in order, into proper triangles
/// with the same winding, appended to `out` as index triples.
///
/// The triangle plus its side corners is a convex polygon whose points all
/// lie on three lines. Collinearity is decided from which sides a point is
/// on, not from coordinates: `t[i]` is on sides `i` and `i-1`, a side
/// corner only on its own side, and three points are collinear exactly when
/// they share a side. Ears are clipped while one exists whose three points
/// are not collinear and whose removal leaves a polygon that is not
/// collinear either; clipping such an ear from a convex polygon leaves a
/// convex polygon, so the loop always finds the next one and ends at a
/// proper triangle. `O(k^2)` for `k` points, and `k` is small. With no
/// corner on any side, `t` itself is appended.
///
/// Returns `false`, leaving `out` as it was, when no ear exists. The side
/// masks this function builds itself never leave the loop without one, so it
/// never returns `false` for them: properness of `t` (three distinct,
/// non-collinear corners) is the caller's precondition, and `false` reports
/// that it was violated. A caller then declines or keeps `t` unsplit; no
/// degenerate piece is ever emitted.
#[must_use]
pub fn split_triangle_at_side_corners(
    t: [usize; 3],
    sides: &[Vec<usize>; 3],
    out: &mut Vec<usize>,
) -> bool {
    let start = out.len();
    // (index, side-membership bitmask)
    let mut ring: Vec<(usize, u8)> =
        Vec::with_capacity(3 + sides.iter().map(Vec::len).sum::<usize>());
    for i in 0..3 {
        ring.push((t[i], (1 << i) | (1 << ((i + 2) % 3))));
        ring.extend(sides[i].iter().map(|&c| (c, 1u8 << i)));
    }
    let collinear =
        |masks: &mut dyn Iterator<Item = u8>| masks.fold(0b111u8, |acc, m| acc & m) != 0;
    while ring.len() > 3 {
        let n = ring.len();
        // Clip the first corner whose removal leaves a ring that is not
        // all on one line. That alone keeps every emitted ear proper: the
        // ring is a triangle with corners on its sides, so the first such
        // corner always joins two different lines (checked exhaustively
        // for up to four corners per side in the tests below).
        let ear = (0..n).find(|&k| {
            !collinear(
                &mut ring
                    .iter()
                    .enumerate()
                    .filter(|&(j, _)| j != k)
                    .map(|(_, p)| p.1),
            )
        });
        let Some(k) = ear else {
            out.truncate(start);
            return false;
        };
        let (prev, next) = ((k + n - 1) % n, (k + 1) % n);
        out.extend_from_slice(&[ring[prev].0, ring[k].0, ring[next].0]);
        ring.remove(k);
    }
    out.extend_from_slice(&[ring[0].0, ring[1].0, ring[2].0]);
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A right triangle with legs of 5 so integer corners fit on each side:
    /// the corner of side `i` at step `j` of `m` is `t[i] + (t[i+1]-t[i])*j/(m+1)`.
    const T: [[f64; 2]; 3] = [[0.0, 0.0], [5.0, 0.0], [0.0, 5.0]];

    fn area2(p: [f64; 2], q: [f64; 2], r: [f64; 2]) -> f64 {
        (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1])
    }

    /// Lay out the triangle and `counts[i]` corners on side `i`; returns the
    /// point table (triangle corners first) and the side lists.
    fn layout(counts: [usize; 3]) -> (Vec<[f64; 2]>, [Vec<usize>; 3]) {
        let mut points = T.to_vec();
        let mut sides: [Vec<usize>; 3] = Default::default();
        for i in 0..3 {
            let (a, b) = (T[i], T[(i + 1) % 3]);
            for j in 1..=counts[i] {
                let f = j as f64 / (counts[i] + 1) as f64;
                sides[i].push(points.len());
                points.push([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f]);
            }
        }
        (points, sides)
    }

    #[test]
    fn a_triangle_without_corners_is_kept() {
        let mut out = Vec::new();
        assert!(split_triangle_at_side_corners(
            [7, 8, 9],
            &Default::default(),
            &mut out
        ));
        assert_eq!(out, [7, 8, 9]);
    }

    /// Every mix of 0..=4 corners per side: only proper triangles with the
    /// parent's winding, covering exactly its area, `k + 1` of them, using
    /// every corner and no collinear piece, appended after existing output.
    #[test]
    fn every_mix_of_corners_gives_a_proper_cover_with_the_same_winding() {
        for counts in (0..125).map(|n| [n % 5, n / 5 % 5, n / 25]) {
            let (points, sides) = layout(counts);
            let mut out = vec![99];
            assert!(split_triangle_at_side_corners([0, 1, 2], &sides, &mut out));
            assert_eq!(out[0], 99, "{counts:?}: existing output kept");
            let pieces = &out[1..];
            let corners: usize = counts.iter().sum();
            assert_eq!(pieces.len(), 3 * (corners + 1), "{counts:?}");
            let parent = area2(T[0], T[1], T[2]);
            let mut total = 0.0;
            let mut seen = vec![false; points.len()];
            for piece in pieces.chunks(3) {
                let area = area2(points[piece[0]], points[piece[1]], points[piece[2]]);
                assert!(area * parent > 0.0, "{counts:?}: {piece:?} is not proper");
                total += area;
                piece.iter().for_each(|&v| seen[v] = true);
            }
            assert!((total - parent).abs() < 1e-9, "{counts:?}: area covered");
            assert!(seen.iter().all(|&s| s), "{counts:?}: every corner used");
        }
    }

    /// The output order is part of the contract.
    #[test]
    fn the_pieces_come_in_a_fixed_order() {
        let (_, sides) = layout([1, 0, 0]);
        let mut out = Vec::new();
        assert!(split_triangle_at_side_corners([0, 1, 2], &sides, &mut out));
        assert_eq!(out, [2, 0, 3, 3, 1, 2]);
    }

    /// A reversed triangle keeps its winding: the same call on mirrored
    /// indices emits pieces of the same (negative) orientation.
    #[test]
    fn a_clockwise_triangle_stays_clockwise() {
        let (points, sides) = layout([2, 1, 3]);
        let mut out = Vec::new();
        // Mirror: t = [0, 2, 1] with sides reversed accordingly.
        let mirrored: [Vec<usize>; 3] = [
            sides[2].iter().rev().copied().collect(),
            sides[1].iter().rev().copied().collect(),
            sides[0].iter().rev().copied().collect(),
        ];
        assert!(split_triangle_at_side_corners(
            [0, 2, 1],
            &mirrored,
            &mut out
        ));
        assert!(out
            .chunks(3)
            .all(|p| area2(points[p[0]], points[p[1]], points[p[2]]) < 0.0));
    }
}
