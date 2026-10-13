//! Keyhole faces of a closed shell: no triangle edge on a seam (#283).
//!
//! [`seam`](super::seam) triangulates a keyhole's region and lets a seam
//! be a triangle edge, once each way, where the ear clipper puts one there.
//! Inside one face that is only a diagonal, but a seam is a segment the
//! ring runs along both ways, and other faces of a shell can run along it
//! too (a rim's two end faces at its cut): with the seam a triangle edge
//! as well, that segment carries four triangle uses and the shell is not a
//! two-manifold. So each seam a triangle edge lies on is taken out of the
//! triangulation, every decision an exact `orient2d` sign:
//!
//! 1. The two triangles on the seam `p -> q`, `(p, q, c)` and `(q, p, d)`,
//!    form the quadrilateral `p d q c`. Where it is strictly convex, its
//!    other diagonal `c d` replaces the seam: no point is added.
//! 2. Otherwise one corner of it is reflex, and no diagonal of the
//!    quadrilateral avoids the seam (a thin rim seamed at its corner
//!    forces this). A Steiner point `s` is placed strictly inside it, off
//!    the seam's line, where the four triangles `(p, d, s)`, `(d, q, s)`,
//!    `(q, c, s)` and `(c, p, s)` are all strictly counter-clockwise. The
//!    seam's midpoint `m` lies strictly inside all four half-planes, so
//!    points close enough to it do too: the candidates are `m` halved
//!    towards `c` or `d` one to [`MAX_HALVINGS`] times
//!    ([`Steiner::place`]), nearest the corners first, and the first that
//!    passes is halved once more ([`place`]) so that it lies inside by a
//!    margin, not on a line through two corners that only the rounding of
//!    a projection moved apart. If none passes, the seam is refused by
//!    name.
//!
//! A candidate is a chain of midpoints of the face's own points, so a
//! caller lifting the face back to space places it with the same chain
//! there: on a face in an axis plane, every midpoint of two points with
//! one coordinate in common keeps that coordinate exactly, and the point
//! lies exactly in the plane. Elsewhere it lies within the rounding of the
//! chain, which the caller reports.
//!
//! The result is certified against the cover it replaces
//! ([`certify_same_cover`]): every triangle strictly counter-clockwise, the
//! same boundary edges once each in the same direction, every other edge
//! once each way, two more triangles per Steiner point, and no triangle
//! edge on a seam. Positive triangles with the boundary of a certified
//! cover tile the same region exactly once, so a Steiner point adds no
//! area and moves no boundary.

use std::collections::HashMap;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;

use super::orient;

/// How often the seam's midpoint is halved towards a corner at most.
///
/// Fifty-three halvings bring a candidate within one unit in the last
/// place of the midpoint for any `f64` input; past that it is the
/// midpoint itself, on the seam.
pub(crate) const MAX_HALVINGS: u32 = 60;

/// A point a sealed triangulation adds: the midpoint `m` of the seam
/// `seam[0] -> seam[1]`, moved half way towards `toward`, `halvings`
/// times. Indices are into the triangulation's points, so `toward` may be
/// an earlier Steiner point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Steiner {
    pub(crate) seam: [u32; 2],
    pub(crate) toward: u32,
    pub(crate) halvings: u32,
}

impl Steiner {
    /// The point by its chain of midpoints, over any point type.
    pub(crate) fn place<T: Copy>(&self, at: impl Fn(u32) -> T, midpoint: impl Fn(T, T) -> T) -> T {
        let m = midpoint(at(self.seam[0]), at(self.seam[1]));
        let mut x = at(self.toward);
        for _ in 0..self.halvings {
            x = midpoint(m, x);
        }
        x
    }
}

/// The midpoint of two points, as [`Steiner::place`] chains it in 2D.
pub(crate) fn midpoint2(a: Point2, b: Point2) -> Point2 {
    Point2::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5)
}

/// Take every seam in `seams` (pairs of point indices, either way round)
/// out of `triangles`, a certified cover over `points`: flip it where its
/// quadrilateral is strictly convex, otherwise split that quadrilateral at
/// a Steiner point, appended to `points`. Returns the Steiner points in
/// the order they were appended, the new triangles certified.
pub(crate) fn seal(
    points: &mut Vec<Point2>,
    triangles: &mut Vec<[u32; 3]>,
    seams: &[[u32; 2]],
) -> GeomResult<Vec<Steiner>> {
    let before = triangles.clone();
    let mut added = Vec::new();
    // Which triangle runs each directed edge, kept current as triangles
    // are replaced.
    let mut owner: HashMap<(u32, u32), usize> = HashMap::with_capacity(3 * triangles.len());
    for (index, t) in triangles.iter().enumerate() {
        for k in 0..3 {
            owner.insert((t[k], t[(k + 1) % 3]), index);
        }
    }
    for &[p, q] in seams {
        let (Some(&one), Some(&other)) = (owner.get(&(p, q)), owner.get(&(q, p))) else {
            continue;
        };
        let c = third(triangles[one], p, q);
        let d = third(triangles[other], q, p);
        let at = |v: u32| points[v as usize];
        let replace: Vec<[u32; 3]> =
            if orient(at(c), at(p), at(d)) > 0 && orient(at(d), at(q), at(c)) > 0 {
                vec![[c, p, d], [d, q, c]]
            } else {
                let Some(steiner) = place(points, p, q, c, d) else {
                    return Err(GeomError::Degenerate(format!(
                        "profile seam from ({}, {}) to ({}, {}) is forced to be a triangle edge: \
                         no point inside its two triangles can take it out",
                        at(p).x,
                        at(p).y,
                        at(q).x,
                        at(q).y
                    )));
                };
                let s = u32::try_from(points.len()).map_err(|_| {
                    GeomError::InvalidInput(
                        "profile Steiner points exceed what a u32 index addresses".to_owned(),
                    )
                })?;
                let position = steiner.place(|v| points[v as usize], midpoint2);
                points.push(position);
                added.push(steiner);
                vec![[p, d, s], [d, q, s], [q, c, s], [c, p, s]]
            };
        for index in [one, other] {
            let t = triangles[index];
            for k in 0..3 {
                owner.remove(&(t[k], t[(k + 1) % 3]));
            }
        }
        let slots = [one, other];
        for (k, t) in replace.into_iter().enumerate() {
            let index = if k < 2 {
                triangles[slots[k]] = t;
                slots[k]
            } else {
                triangles.push(t);
                triangles.len() - 1
            };
            for j in 0..3 {
                owner.insert((t[j], t[(j + 1) % 3]), index);
            }
        }
    }
    certify_same_cover(points, &before, triangles, added.len(), seams)?;
    Ok(added)
}

/// The corner of `t` that is not `a` or `b`.
fn third(t: [u32; 3], a: u32, b: u32) -> u32 {
    t.into_iter().find(|&v| v != a && v != b).unwrap_or(t[0])
}

/// The candidate of the module notes for the quadrilateral `p d q c`:
/// one halving past the first that lies strictly inside it and off the
/// line `p q`, when that one passes too.
///
/// The plane a caller projects a face into is rounded: three corners on
/// one line in space can come out a hair apart in it, and a chain of
/// midpoints lands on such lines readily (half way between two corners,
/// or on a corner's line through the seam's midpoint). The first candidate
/// that passes can therefore sit on the quadrilateral's boundary in space
/// and only inside it by that rounding, a triangle of no area once lifted.
/// The midpoint `m` lies inside by the full clearance of the seam's
/// triangles, and the quadrilateral's open kernel is convex, so the
/// midpoint of `m` and a candidate in its closure lies strictly inside,
/// by at least half that clearance.
fn place(points: &[Point2], p: u32, q: u32, c: u32, d: u32) -> Option<Steiner> {
    let at = |v: u32| points[v as usize];
    let fits = |steiner: &Steiner| {
        let s = steiner.place(at, midpoint2);
        [(p, d), (d, q), (q, c), (c, p)]
            .into_iter()
            .all(|(a, b)| orient(at(a), at(b), s) > 0)
            && orient(at(p), at(q), s) != 0
    };
    for halvings in 1..=MAX_HALVINGS {
        for toward in [c, d] {
            let first = Steiner {
                seam: [p, q],
                toward,
                halvings,
            };
            if fits(&first) {
                let deeper = Steiner {
                    halvings: halvings + 1,
                    ..first
                };
                return Some(if fits(&deeper) { deeper } else { first });
            }
        }
    }
    None
}

/// Refuse `after` unless it covers what the certified cover `before` does:
/// every triangle strictly counter-clockwise, each directed edge once, the
/// same net boundary, `2 * steiner` more triangles, and no triangle edge
/// on a seam.
///
/// The net boundary is each directed edge's uses less its reverse's. Two
/// sets of positive triangles with the same net boundary cover every point
/// the same number of times, the winding number of that boundary; the
/// certified cover covers its region once, so `after` does too.
pub(crate) fn certify_same_cover(
    points: &[Point2],
    before: &[[u32; 3]],
    after: &[[u32; 3]],
    steiner: usize,
    seams: &[[u32; 2]],
) -> GeomResult<()> {
    let fail = |why: String| {
        Err(GeomError::Degenerate(format!(
            "keyhole triangulation failed its certificate: {why}"
        )))
    };
    if after.len() != before.len() + 2 * steiner {
        return fail(format!(
            "{} triangles, the cover it replaces has {} and {steiner} Steiner points add {}",
            after.len(),
            before.len(),
            2 * steiner
        ));
    }
    let net = |triangles: &[[u32; 3]]| {
        let mut net: HashMap<(u32, u32), i64> = HashMap::with_capacity(3 * triangles.len());
        for t in triangles {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let (key, sign) = if a < b { ((a, b), 1) } else { ((b, a), -1) };
                *net.entry(key).or_insert(0) += sign;
            }
        }
        net.retain(|_, n| *n != 0);
        net
    };
    let mut directed: HashMap<(u32, u32), u32> = HashMap::with_capacity(3 * after.len());
    for t in after {
        let corner = |k: usize| points[t[k] as usize];
        if orient(corner(0), corner(1), corner(2)) <= 0 {
            return fail(format!("triangle {t:?} is not strictly counter-clockwise"));
        }
        for k in 0..3 {
            let edge = (t[k], t[(k + 1) % 3]);
            let uses = directed.entry(edge).or_insert(0);
            *uses += 1;
            if *uses > 1 {
                return fail(format!("edge {}->{} is used twice", edge.0, edge.1));
            }
            if seams
                .iter()
                .any(|s| (s[0], s[1]) == edge || (s[1], s[0]) == edge)
            {
                return fail(format!("edge {}->{} lies on a seam", edge.0, edge.1));
            }
        }
    }
    if net(before) != net(after) {
        return fail("its boundary differs from the cover it replaces".to_owned());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y)
    }

    /// A thin rim's corner: the seam 0 -> 1 between triangles whose
    /// quadrilateral is reflex at 1, so it cannot be flipped.
    fn reflex() -> (Vec<Point2>, Vec<[u32; 3]>) {
        let points = vec![p(0.0, 0.0), p(1.0, 1.0), p(4.0, 1.0), p(1.0, 4.0)];
        // (0, 1, 3) and (1, 0, 2): the quadrilateral 0 2 1 3 turns right
        // at 1.
        (points, vec![[0, 1, 3], [1, 0, 2]])
    }

    #[test]
    fn a_convex_quadrilateral_flips_without_a_point() {
        let mut points = vec![p(0.0, 0.0), p(2.0, 2.0), p(2.0, 0.0), p(0.0, 2.0)];
        let mut triangles = vec![[0, 1, 3], [1, 0, 2]];
        let added = seal(&mut points, &mut triangles, &[[0, 1]]).expect("flips");
        assert!(added.is_empty());
        assert_eq!(points.len(), 4);
        assert_eq!(triangles, vec![[3, 0, 2], [2, 1, 3]]);
    }

    #[test]
    fn a_reflex_quadrilateral_takes_a_steiner_point_off_the_seam() {
        let (mut points, mut triangles) = reflex();
        let added = seal(&mut points, &mut triangles, &[[0, 1]]).expect("a point fits");
        assert_eq!(added.len(), 1);
        assert_eq!(triangles.len(), 4);
        let s = points[4];
        assert_ne!(orient(points[0], points[1], s), 0);
        assert_eq!(
            added[0].place(|v| points[v as usize], midpoint2),
            s,
            "the recipe places the point"
        );
    }

    #[test]
    fn the_certificate_names_each_defect() {
        let (mut points, before) = reflex();
        points.push(p(5.0, 5.0));
        let seams = [[0, 1]];
        let reason = |after: &[[u32; 3]], steiner: usize| match certify_same_cover(
            &points, &before, after, steiner, &seams,
        ) {
            Err(GeomError::Degenerate(message)) => message,
            other => panic!("expected a refusal, got {other:?}"),
        };
        assert!(reason(&before, 0).contains("lies on a seam"));
        assert!(reason(&before, 1).contains("Steiner points add 2"));
        assert!(reason(&[[0, 3, 1], [1, 2, 0]], 0).contains("not strictly counter-clockwise"));
        // The right count and winding, but over a different boundary.
        assert!(reason(&[[0, 2, 4], [0, 4, 3]], 0).contains("boundary differs"));
        assert!(reason(&[[0, 2, 4], [0, 2, 3]], 0).contains("used twice"));
    }
}
