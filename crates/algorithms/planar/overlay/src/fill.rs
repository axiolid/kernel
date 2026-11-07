//! Notches filled for the fast disc morphology, with a stated excess
//! (#330).
//!
//! [`crate::minkowski`] builds a disc morphology from the sum `S + P` of a
//! set `S` (the region for a dilation, its complement for an erosion) with
//! a convex polygon `P`. Along a wall of many short teeth, after #323's
//! chord covers, what is left of each edge hull is a sliver at a tooth's
//! tip, and with the teeth closer than the slivers are wide each one
//! overlaps dozens of others: about 905,000 pieces of the subdivision for
//! 10,002 corners, and time and memory grew with the square of the
//! corners. The exact union of the slivers has cusps that are crossing
//! points, not vertex sums, so the exact route cannot do with fewer terms.
//! [`crate::MorphologyMode::Fast`] instead sums a set `S'` that holds `S`,
//! its notches between the teeth's tips filled by chords, and states how
//! much larger that makes the sum.
//!
//! # The notch and its excess
//!
//! A chord `a b` between two vertices of one ring of `S` is proved, with
//! #323's exact proof ([`crate::cover::proved`]) taken for the other side,
//! to lie in the closure of the complement of `S`: the vertices it spans
//! lie on the side of `S` or on it, it leaves both ends outside `S`'s
//! sector, and no other boundary crosses or touches it. It closes a
//! notch `N`, the polygon of the chord and the boundary it spans, and
//! `S' = S union N`. The ring keeps the chord as one edge instead of the
//! boundary it spans. A chord is used only when, in addition:
//!
//! - no other ring lies inside `N` (one vertex per ring near it is tested
//!   exactly: a ring meeting neither the chord nor the spanned boundary is
//!   wholly inside or wholly outside), so the rings of `S'` are nested as
//!   those of `S` are;
//! - every spanned vertex lies within `t` of the chord's line, `t` a lower
//!   bound on the radius of a disc about the origin holding `P` (the
//!   proposals reach no deeper, so a notch deeper than `t` stays: its slab
//!   bound below is a bound `excess` has, not one `fill` reaches);
//! - its excess `e`, below, is within the caller's budget.
//!
//! Two chords used never meet but at a vertex both end at. A chord of
//! another ring, or of another run of the same ring, ends outside `N` (a
//! ring other than the chord's lies wholly outside it, and the rest of the
//! chord's own ring leaves `a` outside `N`'s corner there, by the sector
//! test), and so would cross `N`'s boundary twice: the spanned boundary,
//! which its own proof refuses, or the chord, which a segment crosses at
//! most once. A vertex of either on the other its proof refuses too. So the
//! rings of `S'` stay simple, apart and nested as those of `S` are.
//!
//! Then `S' + B(t) <= S + B(t + e)` for the disc `B(t)` of radius `t`,
//! where for spanned vertices `u` at distance `d_u` from the line and at
//! `s_u` along it, `a` and `b` among them at distance 0,
//!
//! ```text
//! e = min( H,  min over subsets U:  max d_u (u in U) + g(U)^2 / (8 t) )
//! ```
//!
//! `H` the greatest `d_u` and `g(U)` the greatest gap between consecutive
//! `s_u` of `U` (with `a` and `b`). Proof: let `|x - n| <= t` with `n` in
//! `N`. If `x` lies in `N`, the perpendicular from `x` away from the chord
//! meets the spanned boundary, a part of `S`, within `H <= t`. Otherwise
//! the segment from `x` to `n` enters `N` through the spanned boundary
//! (a point of `S` within `t`) or through the chord at `q`, `|x - q| <= t`.
//! Let `f` be the foot of `x` on the line, `h <= t` its height. If `f`
//! lies beyond `b`, `|x - b| <= |x - q|`, and so beyond `a`. Otherwise
//! `f` lies between two
//! consecutive members of `U` and within `g / 2` along the line of one,
//! `u`, which lies on the far side of the line from `x`:
//! `|x - u|^2 <= (t + d_u)^2 + g^2 / 4`, so
//! `|x - u| <= t + d_u + g^2 / (8 t)`. And always `|x - s| <= t + H` for
//! the point `s` of `S` nearest `n`.
//!
//! # The slab bound
//!
//! A narrow notch has a second bound, one that does not depend on its
//! depth, for a disc of radius `t` cannot reach far into a needle. Let
//! every spanned vertex project onto the chord `a b`, of length `g`: then `N`, the polygon of
//! those vertices, lies in the slab between the perpendiculars to the line
//! at `a` and `b`, a convex set. Then
//!
//! ```text
//! e_W = sqrt(t^2 + g^2 / 4) - t = (g^2 / 4) / (sqrt(t^2 + g^2 / 4) + t)
//! ```
//!
//! is an excess whatever the depth (`fill` only proposes notches within
//! `t`, so for it the bound competes with the others). Proof: let
//! `|x - n| <= t` with `n` in `N`, and `L` the line through `n` parallel to the chord. The part of
//! `L` in `N` through `n` is a segment `p q` of length at most `g`, the
//! width of the slab. Its ends lie on the boundary of `N`; unless `L` is
//! the chord's own line (where they are `a` and `b`, in `S`) they cannot
//! lie on the chord, so they lie on the spanned boundary, a part of `S`.
//! Let `f` be the foot of `x` on `L` and `h <= t` the distance from `x` to
//! `L`. If `f` lies between `p` and `q`, the nearer end, say `p`, is at
//! most `g / 2` from `f` along `L`: `|x - p|^2 <= h^2 + g^2 / 4 <=
//! t^2 + g^2 / 4`. Otherwise `f` lies beyond one end, say `q`, and as `n`
//! lies between `p` and `q`, `|f - q| <= |f - n|`, so `|x - q|^2 = h^2 +
//! |f - q|^2 <= h^2 + |f - n|^2 = |x - n|^2 <= t^2`. Either way `x` lies
//! within `sqrt(t^2 + g^2 / 4)` of `S`, so `S' + B(t) <= S + B(t + e_W)`.
//! The slab is tested in `f64` with the same slack as the other bounds
//! (each vertex projects at least `slack` inside the chord's ends), and
//! `g` is the chord's length plus the slack. `e_W` joins the minimum above;
//! for a notch deeper than `t` the minimum is that of `e_W` and `H`, the
//! only bounds it has, and an overhanging vertex, beyond a chord's end,
//! with a deep notch is refused. Fast still keeps a needle deeper than `t`,
//! because nothing proposes a chord over it.
//!
//! Along a toothed wall the tips lie on the chord, so `e` is about the
//! square of the tips' spacing over `8 t`: 2 micrometres for teeth 1.8 mm
//! apart and a disc of 0.4 m (`t` 0.2 m), where `H`, the depth of the
//! teeth, is 0.5 mm. The values are computed in `f64` and widened by a
//! slack of 64 ulps of the largest coordinate, so `e` is an upper bound.
//!
//! A proposal that fails a condition is split in two, each part tried in
//! turn and proved again, so a refused chord gives way to shorter ones,
//! down to single edges, which the set keeps as they are. One too deep or
//! over the budget is split at its spanned vertex farthest from its line
//! (Douglas-Peucker); one not proved, or holding a ring, at its nearest,
//! which a chord along a wall's tips passes through. Of equally far
//! vertices the one nearest the middle is taken, so a long wall splits
//! into halves and costs about `n log n`.
//!
//! A chord that is not proved is tried without its start edge, without its
//! end edge, then without both, before it is split (#345). Along a wall
//! of teeth the chord wanted runs between two tips, but the split of a
//! proposal over a room's corner (too deep for the budget, split at the
//! corner) leaves a chord from the corner, a valley of the wall, along it.
//! That chord is not proved, for the tips between lie on the wrong side of
//! it, and neither is either half the nearest-vertex split makes: it halves
//! a chord between valleys at a valley again and again, so every tooth
//! between the corner and the far end stayed, about an eighth of a floor's
//! teeth, and the sum of that many teeth cost what #330's did. The parity
//! of the teeth decides whether the far end is a tip or a valley, hence
//! the three trims. The first one used replaces the chord: it is proved,
//! within the disc and the budget, with no ring in its notch, as any
//! chord is, so the excess and the side are those of the other chords.
//!
//! The proposals are [`cover::proposals`] for discs of radius `t`: chords
//! whose spanned vertices lie within `t` of the line (#340; #330 proposed
//! within `t / 16`, and so filled no needle deeper than that).
//!
//! # The result's side
//!
//! [`crate::minkowski`] sums `S'` with a polygon `P` inside `B(t')`,
//! `t' >= t`, so `S + P <= S' + P <= S + B(t' + e)`; `t` is half the
//! radius less the margin, and the budget is capped at `t`. A
//! circumscribed polygon (an inner erosion, an outer dilation) holds the
//! disc as before, so its side holds, and `e` is added to its deviation.
//! An inscribed one (an inner dilation, an outer erosion) is drawn in by
//! `e`, to `t' = r - margin - e >= t`, so that `S' + P` still lies within
//! the disc `B(r - margin)` that holds the exact route's polygon, and its
//! deviation `r - t' cos(pi / 64)` grows by about `e`. For an erosion the
//! summed set is the complement, so its notches are filled where the
//! region's teeth stand: the region loses them, and its own rings are
//! taken filled too, so the subdivision shares them with the sum's copy.
//! `scripts/probe_fast_morphology_mutants.py` lists the faults the tests
//! must catch.

use axiolid_core::Point2;
use axiolid_guarantees::Sign;

use crate::cover::{self, Chord, Index};
use crate::exact_arc::{orient_doubles as orient, Bounds};

/// The rings of `S'` and the greatest excess of a chord used.
#[derive(Debug, Clone)]
pub(crate) struct Filled {
    /// One ring per ring of `S`, in order: as given, or with notches
    /// filled.
    pub(crate) rings: Vec<Vec<Point2>>,
    /// The greatest excess `e` of a chord used; zero when none is.
    pub(crate) excess: f64,
    /// The number of chords used, for the tests.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) chords: usize,
}

/// Whether `q` lies inside the simple polygon `polygon`, by the parity of
/// crossings decided with exact orientations. `q` is not on its boundary.
fn inside(polygon: &[Point2], q: Point2) -> bool {
    let n = polygon.len();
    let mut odd = false;
    for i in 0..n {
        let (p, r) = (polygon[i], polygon[(i + 1) % n]);
        if (p.y > q.y) != (r.y > q.y) {
            let o = orient(p, r, q);
            if (r.y > p.y && o == Sign::Positive) || (r.y < p.y && o == Sign::Negative) {
                odd = !odd;
            }
        }
    }
    odd
}

/// The points of a chord's notch: its start, the vertices it spans, and
/// its end.
fn notch(points: &[Point2], chord: Chord) -> Vec<Point2> {
    let n = points.len();
    (0..=chord.span)
        .map(|k| points[(chord.from + k) % n])
        .collect()
}

/// Whether no ring other than `ring` lies inside the notch `polygon`.
fn encloses_none(ring: usize, polygon: &[Point2], index: &Index<'_>) -> bool {
    let query = Bounds::around(polygon, 0.0);
    let (mut found, mut scratch) = (Vec::new(), Vec::new());
    index
        .tree
        .query_into(|q| q.overlaps(&query), &mut found, &mut scratch);
    let mut tested: Vec<usize> = Vec::new();
    found.into_iter().all(|item| {
        let (r, e) = index.edges[item];
        if r == ring || tested.contains(&r) {
            return true;
        }
        tested.push(r);
        !inside(polygon, index.rings[r][e])
    })
}

/// The excess of the notch `polygon` (its chord from the first point to
/// the last) for discs of radius `t` or more, or `None` when a spanned
/// vertex lies farther than `t` from the chord's line and one projects
/// beyond the chord's ends (module docs).
fn excess(polygon: &[Point2], t: f64) -> Option<f64> {
    let (a, b) = (polygon[0], polygon[polygon.len() - 1]);
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len = dx.hypot(dy);
    if !(len > 0.0 && len.is_finite() && t > 0.0) {
        return None;
    }
    let scale = polygon
        .iter()
        .fold(0.0_f64, |m, p| m.max(p.x.abs()).max(p.y.abs()));
    let slack = 64.0 * f64::EPSILON * (scale + len);
    // Along the line and off it, every spanned vertex; the ends at 0.
    let mut at: Vec<(f64, f64)> = polygon[1..polygon.len() - 1]
        .iter()
        .map(|v| {
            let (ux, uy) = (v.x - a.x, v.y - a.y);
            ((ux * dx + uy * dy) / len, (ux * dy - uy * dx).abs() / len)
        })
        .collect();
    let deepest = at.iter().fold(0.0_f64, |m, &(_, d)| m.max(d)) + slack;
    // The slab: every spanned vertex projects strictly inside the chord.
    let slab = at.iter().all(|&(s, _)| s >= slack && s <= len - slack);
    if deepest.is_nan() || (deepest > t && !slab) {
        return None;
    }
    // The slab bound, which holds however deep the notch is.
    let g = len + slack;
    let wide = slab.then(|| {
        let h = 0.25 * g * g;
        h / ((t * t + h).sqrt() + t) * (1.0 + 1e-12)
    });
    if deepest > t {
        // The depth `H` bounds it too (module docs), the slab bound being
        // the smaller only for a narrow notch.
        return wide.map(|w| w.min(deepest)).filter(|e| e.is_finite());
    }
    at.push((0.0, 0.0));
    at.push((len, 0.0));
    at.sort_by(|p, q| p.0.total_cmp(&q.0));
    // The supports: vertices within `limit` of the line, for halving
    // limits; the best bound over them, or the notch's depth.
    let mut best = deepest;
    let mut limit = deepest;
    for _ in 0..48 {
        let (mut last, mut gap, mut top) = (None::<f64>, 0.0_f64, 0.0_f64);
        for &(s, d) in &at {
            if d <= limit {
                if let Some(l) = last {
                    gap = gap.max(s - l);
                }
                last = Some(s);
                top = top.max(d);
            }
        }
        let g = gap + 2.0 * slack;
        best = best.min(top + slack + g * g / (8.0 * t));
        if limit == 0.0 {
            break;
        }
        limit = if limit < slack { 0.0 } else { 0.5 * limit };
    }
    let mut bound = best * (1.0 + 1e-12);
    if let Some(w) = wide {
        bound = bound.min(w);
    }
    bound.is_finite().then_some(bound)
}

/// How a chord fares (`fill`): used with its excess, or refused.
enum Verdict {
    /// Proved, within the disc and the budget, and no ring in its notch.
    Used(f64),
    /// Not proved to lie outside the set.
    Unproved,
    /// Too deep for the disc, or over the budget.
    TooBig,
    /// Another ring lies in its notch.
    Holds,
}

/// Where to split a refused `chord`, as a step from its start: at the
/// spanned vertex farthest from its line when it was too deep or over the
/// budget (`deepest`), else (no proof, or a ring in its notch) at the
/// nearest, which a chord along a wall's tips or valleys passes through;
/// of equally far ones, the one nearest the middle, so that the parts are
/// about even.
fn split_step(points: &[Point2], chord: Chord, deepest: bool) -> usize {
    let polygon = notch(points, chord);
    let (a, b) = (polygon[0], polygon[polygon.len() - 1]);
    let off = |v: Point2| ((v.x - a.x) * (b.y - a.y) - (v.y - a.y) * (b.x - a.x)).abs();
    let from_middle = |i: usize| (2 * i).abs_diff(chord.span);
    (1..chord.span)
        .max_by(|&i, &j| {
            let (u, v) = (off(polygon[i]), off(polygon[j]));
            let far = if deepest {
                u.total_cmp(&v)
            } else {
                v.total_cmp(&u)
            };
            far.then(from_middle(j).cmp(&from_middle(i)))
                .then(j.cmp(&i))
        })
        .unwrap_or(chord.span / 2)
}

/// `S'` for the set bounded by `rings`, the set lying on side `sides[r]`
/// of ring `r`'s travel, for discs of radius `t` or more: each chord's
/// excess within `budget` (module docs). `linear` is the tolerance: a ring
/// filled to a zero area under it keeps its notches.
pub(crate) fn fill(
    rings: &[Vec<Point2>],
    sides: &[Sign],
    t: f64,
    budget: f64,
    linear: f64,
) -> Filled {
    let unchanged = || Filled {
        rings: rings.to_vec(),
        excess: 0.0,
        chords: 0,
    };
    if !(t > 0.0 && t.is_finite() && budget >= 0.0 && budget.is_finite()) {
        return unchanged();
    }
    let index = Index::new(rings);
    let mut out = rings.to_vec();
    let (mut most, mut count) = (0.0_f64, 0);
    for (r, ring) in rings.iter().enumerate() {
        if sides[r] == Sign::Zero {
            continue;
        }
        let boundary = cover::merged(ring);
        let points = &boundary.points;
        let n = points.len();
        // The chord lies outside the set: in the set on the other side.
        let outside = sides[r].flip();
        // The chord's verdict (module docs): used with its excess, or why not.
        let judge = |chord: Chord| {
            if !cover::proved(r, &boundary, chord, outside, &index) {
                return Verdict::Unproved;
            }
            let polygon = notch(points, chord);
            match excess(&polygon, t) {
                Some(e) if e <= budget => {
                    if encloses_none(r, &polygon, &index) {
                        Verdict::Used(e)
                    } else {
                        Verdict::Holds
                    }
                }
                _ => Verdict::TooBig,
            }
        };
        // The chords used and their excess.
        let mut used: Vec<(Chord, f64)> = Vec::new();
        let mut pending = cover::proposals(&boundary, outside, cover::Reach::within(t));
        pending.reverse();
        while let Some(chord) = pending.pop() {
            if chord.span < 2 {
                continue;
            }
            let verdict = judge(chord);
            if let Verdict::Used(e) = verdict {
                used.push((chord, e));
                continue;
            }
            // A chord not proved may have a valley of the wall for an end,
            // where a tip is wanted (module docs): it is tried without an
            // end edge, or both, in turn.
            if matches!(verdict, Verdict::Unproved) {
                let trimmed = [(1, 0), (0, 1), (1, 1)].into_iter().find_map(|(a, b)| {
                    let span = chord.span.checked_sub(a + b).filter(|&s| s >= 2)?;
                    let c = Chord {
                        from: (chord.from + a) % n,
                        span,
                    };
                    match judge(c) {
                        Verdict::Used(e) => Some((c, e)),
                        _ => None,
                    }
                });
                if let Some(used_chord) = trimmed {
                    used.push(used_chord);
                    continue;
                }
            }
            // Split, the first part tried first.
            let deepest = matches!(verdict, Verdict::TooBig);
            let k = split_step(points, chord, deepest);
            pending.push(Chord {
                from: (chord.from + k) % n,
                span: chord.span - k,
            });
            pending.push(Chord {
                from: chord.from,
                span: k,
            });
        }
        if used.is_empty() {
            continue;
        }
        let mut spanned = vec![false; n];
        for (chord, _) in &used {
            for k in 1..chord.span {
                spanned[(chord.from + k) % n] = true;
            }
        }
        let filled: Vec<Point2> = (0..n).filter(|&i| !spanned[i]).map(|i| points[i]).collect();
        if filled.len() < 3 || crate::zero_area(&filled, linear) {
            continue;
        }
        out[r] = filled;
        count += used.len();
        most = used.iter().fold(most, |m, &(_, e)| m.max(e));
    }
    Filled {
        rings: out,
        excess: most,
        chords: count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn points(p: &[(f64, f64)]) -> Vec<Point2> {
        p.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    /// A room's bottom wall with `teeth` teeth of `depth` into the room,
    /// from `(0, 0)` to `(4, 0)`, then two corners back: counter-clockwise.
    fn room(teeth: usize, depth: f64) -> Vec<Point2> {
        let mut ring: Vec<Point2> = (0..=2 * teeth)
            .map(|k| {
                let x = 4.0 * k as f64 / (2 * teeth) as f64;
                Point2::new(x, if k % 2 == 1 { depth } else { 0.0 })
            })
            .collect();
        ring.extend(points(&[(4.0, 3.0), (0.0, 3.0)]));
        ring
    }

    #[test]
    fn the_notches_between_the_tips_are_filled_with_a_small_excess() {
        // 200 teeth 0.01 deep, 0.02 apart: e = 0.02^2 / (8 * 0.5) = 1e-4,
        // a hundredth of the depth.
        let rings = vec![room(200, 0.01)];
        let close = |e: f64| (1e-4..1e-4 * (1.0 + 1e-6)).contains(&e);
        // The room itself (a dilation): the valleys are its outermost
        // corners, and one chord along them fills the notches.
        let room_side = fill(&rings, &[Sign::Positive], 0.5, 1e-3, 1e-9);
        assert_eq!(room_side.chords, 1);
        assert!(close(room_side.excess), "{}", room_side.excess);
        assert_eq!(
            room_side.rings[0],
            points(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)])
        );
        // Its complement (an erosion): the chord runs along the tips.
        let complement = fill(&rings, &[Sign::Negative], 0.5, 1e-3, 1e-9);
        assert_eq!(complement.chords, 1);
        assert!(close(complement.excess), "{}", complement.excess);
        assert_eq!(
            complement.rings[0],
            points(&[
                (0.0, 0.0),
                (0.01, 0.01),
                (3.99, 0.01),
                (4.0, 0.0),
                (4.0, 3.0),
                (0.0, 3.0)
            ])
        );
        // A budget below the excess: no chord, the ring as it was.
        for side in [Sign::Positive, Sign::Negative] {
            let none = fill(&rings, &[side], 0.5, 5e-5, 1e-9);
            assert_eq!((none.chords, none.excess), (0, 0.0));
            assert_eq!(none.rings, rings);
        }
        // A smaller disc sees farther into the notches.
        let small = fill(&rings, &[Sign::Positive], 0.25, 1e-3, 1e-9);
        assert!((2e-4..2.001e-4).contains(&small.excess), "{}", small.excess);
        // A zero-area ring or no side: nothing filled.
        assert_eq!(fill(&rings, &[Sign::Zero], 0.5, 1e-3, 1e-9).chords, 0);
    }

    #[test]
    fn the_excess_is_the_gap_between_supports_or_the_depth() {
        // Three supports on the line 0.1 apart, valleys 0.05 deep between.
        let notch = points(&[
            (0.0, 0.0),
            (0.05, 0.05),
            (0.1, 0.0),
            (0.15, 0.05),
            (0.2, 0.0),
        ]);
        let e = excess(&notch, 1.0).unwrap();
        let want = 0.1f64.powi(2) / 8.0;
        assert!((want..want * (1.0 + 1e-9)).contains(&e), "{e}");
        // Deeper than the disc, inside the slab: the depth, `H`, the
        // smaller of it and the slab bound here (0.068).
        let e = excess(&notch, 0.04).unwrap();
        assert!((0.05..0.05 * (1.0 + 1e-9)).contains(&e), "{e}");
        // With no support between, the depth bounds it.
        let deep = points(&[(0.0, 0.0), (0.1, 0.001), (0.2, 0.0)]);
        let e = excess(&deep, 1.0).unwrap();
        assert!((0.001..0.001 * (1.0 + 1e-9)).contains(&e), "{e}");
    }

    #[test]
    fn a_ring_inside_the_notch_is_found() {
        let wall = room(1, 0.5);
        let island = points(&[(1.9, 0.1), (2.1, 0.1), (2.0, 0.2)]);
        let outside = points(&[(1.9, 1.1), (2.1, 1.1), (2.0, 1.2)]);
        let index_rings = vec![wall.clone(), island.clone()];
        let index = Index::new(&index_rings);
        let polygon = points(&[(0.0, 0.0), (2.0, 0.5), (4.0, 0.0)]);
        assert!(!encloses_none(0, &polygon, &index));
        let index_rings = vec![wall, outside];
        let index = Index::new(&index_rings);
        assert!(encloses_none(0, &polygon, &index));
    }

    #[test]
    fn a_refused_chord_gives_way_to_shorter_ones() {
        // An island of the set inside the first tooth of the room's wall:
        // the chord along all the valleys would hold it. Split at the
        // valley in the middle, again and again towards the island, it
        // gives way to a few chords along the valleys, and the first tooth
        // keeps its notch.
        let island = points(&[(0.009, 0.002), (0.011, 0.002), (0.01, 0.004)]);
        let rings = vec![room(200, 0.01), island.clone()];
        let filled = fill(&rings, &[Sign::Positive; 2], 0.5, 1e-3, 1e-9);
        assert!((2..=10).contains(&filled.chords), "{}", filled.chords);
        let ring = &filled.rings[0];
        assert_eq!(ring[..3], points(&[(0.0, 0.0), (0.01, 0.01), (0.02, 0.0)]));
        // Every other corner of the wall lies on its line.
        assert!(ring.len() < 20 && ring[3..ring.len() - 2].iter().all(|p| p.y == 0.0));
        assert_eq!(
            ring[ring.len() - 3..],
            points(&[(4.0, 0.0), (4.0, 3.0), (0.0, 3.0)])
        );
        assert_eq!(filled.rings[1], island);
    }

    #[test]
    fn a_chord_another_ring_crosses_is_not_used() {
        // An island of the set straddling the line through the valleys
        // under the tooth at 0.21: no chord may run along that line past
        // it, and that tooth keeps its notch.
        let island = points(&[(0.21, -0.001), (0.211, 0.0), (0.21, 0.001), (0.209, 0.0)]);
        let rings = vec![room(200, 0.01), island.clone()];
        let filled = fill(&rings, &[Sign::Positive; 2], 0.5, 1e-3, 1e-9);
        assert!(filled.chords >= 2, "{}", filled.chords);
        assert!(filled.rings[0].contains(&Point2::new(0.21, 0.01)));
        assert!(filled.rings[0].len() < 20, "{}", filled.rings[0].len());
        assert_eq!(filled.rings[1], island);
    }

    /// A room whose bottom wall has `teeth` needles of `depth`, their tips
    /// `spacing` apart: counter-clockwise, `teeth * spacing` wide.
    fn needle_room(teeth: usize, depth: f64, spacing: f64) -> Vec<Point2> {
        let width = teeth as f64 * spacing;
        let mut ring: Vec<Point2> = (0..=2 * teeth)
            .map(|k| {
                let x = width * k as f64 / (2 * teeth) as f64;
                Point2::new(x, if k % 2 == 1 { depth } else { 0.0 })
            })
            .collect();
        ring.extend([Point2::new(width, 3.0), Point2::new(0.0, 3.0)]);
        ring
    }

    #[test]
    fn a_needle_wall_within_the_disc_is_filled_with_the_gap_excess() {
        // Teeth 0.05 deep (within `t`), tips 2e-4 apart, a disc of 0.2: the
        // chord along the valleys is filled, its excess the tips' gap
        // bound, `g^2 / (8 t)` = 2.5e-8.
        let rings = vec![needle_room(1000, 0.05, 2e-4)];
        let f = fill(&rings, &[Sign::Positive], 0.2, 1e-3, 1e-9);
        assert_eq!(f.chords, 1);
        assert!(f.excess < 3e-8, "{}", f.excess);
        assert_eq!(f.rings[0].len(), 4);
    }

    /// #345's floor: a square of 3 m with `corners` corners, 0.5 mm teeth
    /// on the bottom and top.
    pub(super) fn floor(corners: usize) -> Vec<Point2> {
        let n = (corners / 2 - 1) as u32;
        let at = |k: u32| 3.0 * f64::from(k) / f64::from(n);
        let tooth = |k: u32| if k % 2 == 1 { 5e-4 } else { 0.0 };
        let mut p: Vec<Point2> = (0..=n).map(|k| Point2::new(at(k), tooth(k))).collect();
        p.extend((0..=n).map(|k| Point2::new(3.0 - at(k), 3.0 - tooth(k))));
        p
    }

    #[test]
    fn a_floors_corner_chord_leaves_no_wall_of_teeth() {
        // #345: the floor's erosion fills the complement's notches. The
        // proposal over the room's corner (deeper than the budget) is split
        // at the corner, and the chord left, from the corner's valley to a
        // tip or a valley of the wall, is trimmed to the tips; whichever
        // the parity of the teeth, the wall is one chord, not one per tooth.
        for corners in [2_000, 2_002, 2_004, 20_000, 20_002] {
            for t in [0.7, 0.745, 0.8] {
                let rings = vec![floor(corners)];
                let f = fill(&rings, &[Sign::Negative], t, 1e-3, 1e-9);
                let what = format!("{corners} corners, t {t}");
                assert!(f.chords <= 8, "{what}: {} chords", f.chords);
                assert!(f.rings[0].len() <= 16, "{what}: {}", f.rings[0].len());
                assert!(f.excess > 0.0 && f.excess <= 1e-3, "{what}: {}", f.excess);
            }
        }
    }

    #[test]
    fn a_needle_wall_deeper_than_the_disc_is_kept() {
        // Teeth 0.25 deep, a disc of 0.2: no chord is proposed over them.
        let rings = vec![needle_room(100, 0.25, 2e-4)];
        let f = fill(&rings, &[Sign::Positive], 0.2, 1e-3, 1e-9);
        assert_eq!(f.chords, 0);
        assert_eq!(f.rings, rings);
    }

    #[test]
    fn the_slab_bound_holds_for_a_needle_deeper_than_the_disc() {
        // One valley 0.05 deep, g = 2e-4, t = 0.2 (inside the disc): the
        // slab bound, (g^2 / 4) / (sqrt(t^2 + g^2 / 4) + t).
        let t = 0.2_f64;
        let needle = points(&[(0.0, 0.0), (1e-4, 0.05), (2e-4, 0.0)]);
        let want = 1e-8 / ((t * t + 1e-8).sqrt() + t);
        let e = excess(&needle, t).unwrap();
        assert!((want..want * (1.0 + 1e-9)).contains(&e), "{e} {want}");
        // Deeper than the disc, inside the slab: bounded by the slab, not
        // refused (`H` = 0.5 is the larger).
        let deep = points(&[(0.0, 0.0), (1e-4, 0.5), (2e-4, 0.0)]);
        let e = excess(&deep, t).unwrap();
        assert!((want..want * (1.0 + 1e-9)).contains(&e), "{e} {want}");
        // An overhang past the chord's end, deeper than the disc: refused.
        let overhang = points(&[(0.0, 0.0), (-1e-4, 0.5), (2e-4, 0.0)]);
        assert!(excess(&overhang, t).is_none());
    }
}

#[cfg(test)]
mod diagnosis {
    //! #345: how the fast fill fares on a floor with sub-budget teeth, and
    //! what it leaves. Run with
    //! `cargo test -p axiolid-overlay --release --lib -- --ignored --nocapture diagnosis`.
    use super::tests::floor;
    use super::*;

    #[test]
    #[ignore = "diagnosis for #345, prints"]
    fn the_floor_with_sub_budget_teeth() {
        for corners in [2_000, 50_000, 100_000, 200_000] {
            for (what, side) in [("erosion", Sign::Negative), ("dilation", Sign::Positive)] {
                for r in [1.4_f64, 1.49, 1.6] {
                    // `t` is half the radius; the budget 1 mm.
                    let rings = vec![floor(corners)];
                    let f = fill(&rings, &[side], 0.5 * r, 1e-3, 1e-9);
                    println!(
                        "{corners} corners, r {r}, {what}: {} chords, ring {} corners, excess {:.1e}",
                        f.chords,
                        f.rings[0].len(),
                        f.excess
                    );
                }
            }
        }
    }
}
