//! Certified triangulation of a polygon with holes (#253).
//!
//! Three steps, every decision an exact `orient2d` sign:
//!
//! 1. [`validate`] refuses, by name, rings that do not bound a polygon with
//!    holes: a non-finite or repeated vertex, a ring that folds back on or
//!    crosses itself, holes that overlap or touch each other or the outer
//!    ring, a hole outside the outer ring or inside another hole. It also
//!    reads each ring's orientation, so rings may come either way round.
//! 2. [`bridge`] joins each hole to the outer boundary by a segment checked
//!    against every edge, giving one weakly simple ring whose bridge
//!    vertices appear twice.
//! 3. [`clip`] cuts ears from that ring. An ear is refused if any other
//!    vertex lies in its closed triangle, so a diagonal never runs through
//!    a vertex.
//!
//! The result is then certified, not trusted: every triangle strictly
//! counter-clockwise, every ring edge used exactly once in the direction
//! that keeps the polygon on its left, every other edge exactly once in
//! each direction. Positive triangles whose boundary is exactly the rings
//! tile the polygon once, so the certificate is the whole contract; a
//! triangulation that fails it is refused, never returned.
//!
//! Rings that touch at single points are refused under
//! [`PinchPolicy::Refuse`], which solids use: a pinch extrudes to a
//! non-manifold edge. Under [`PinchPolicy::Accept`], for 2D regions and
//! planar surface patches, [`pinch`] triangulates around them (#262).
//!
//! This replaced `earcut` (ADR 0083). earcut drops nodes where the bridged
//! ring runs straight on and lets only reflex nodes block an ear, so for two
//! holes in one horizontal band a cap triangle's edge ran along the band's
//! bottom line past both holes' inner corners: a T-junction, and an open
//! extrusion.

use std::collections::HashMap;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;
use axiolid_guarantees::Sign;
use axiolid_predicates::orient2d;

use crate::profile::PinchPolicy;

mod bridge;
mod clip;
mod pinch;
mod validate;

/// Triangulate `outer` with `holes`, over the vertex list `outer ++ holes`.
///
/// Rings may be given either way round; the triangles are always
/// counter-clockwise. Under [`PinchPolicy::Accept`] a point where rings
/// touch is referenced by its first index in `outer ++ holes`.
pub(crate) fn triangulate_rings(
    outer: &[Point2],
    holes: &[Vec<Point2>],
    policy: PinchPolicy,
) -> GeomResult<(Vec<Point2>, Vec<[u32; 3]>)> {
    let (points, rings) = gather(outer, holes)?;
    validate::check_rings(&points, &rings)?;
    let touches = validate::check_edges(&points, &rings, policy)?;
    if let Some(canon) = pinch::canonical(&points, &touches) {
        // Only `Accept` lets rings touch, so only it gets here.
        let triangles = pinch::triangulate(&points, &rings, &canon, &touches)?;
        return Ok((points, triangles));
    }
    let loops = validate::orient_and_place(&points, &rings)?;
    let mut polygon = bridge::bridge_holes(&points, &loops)?;
    let triangles = clip::clip_ears(&points, &mut polygon)?;
    certify(&points, &loops, &triangles, 1, loops.len() - 1)?;
    Ok((points, triangles))
}

/// The vertex list `outer ++ holes` with each ring's range in it, refused
/// when a `u32` cannot index it.
fn gather(
    outer: &[Point2],
    holes: &[Vec<Point2>],
) -> GeomResult<(Vec<Point2>, Vec<core::ops::Range<usize>>)> {
    let mut points = Vec::with_capacity(outer.len() + holes.iter().map(Vec::len).sum::<usize>());
    let mut rings = Vec::with_capacity(1 + holes.len());
    for ring in core::iter::once(outer).chain(holes.iter().map(Vec::as_slice)) {
        let start = points.len();
        points.extend_from_slice(ring);
        rings.push(start..points.len());
    }
    if u32::try_from(points.len()).is_err() {
        return Err(GeomError::InvalidInput(format!(
            "profile has {} vertices, more than a u32 index addresses",
            points.len()
        )));
    }
    Ok((points, rings))
}

/// The vertices lying inside another edge, as `(ring, edge, vertex)`, from
/// the validation [`triangulate_rings`] runs under [`PinchPolicy::Accept`]
/// (#265), refused as that validation refuses them.
pub(crate) fn ring_touches(
    outer: &[Point2],
    holes: &[Vec<Point2>],
) -> GeomResult<Vec<(usize, usize, usize)>> {
    let (points, rings) = gather(outer, holes)?;
    validate::check_rings(&points, &rings)?;
    let touches = validate::check_edges(&points, &rings, PinchPolicy::Accept)?;
    // A vertex inside an edge is found once per edge of its own that meets
    // it there; list it once, in the order validation first finds it.
    let mut listed = Vec::with_capacity(touches.len());
    for touch in touches {
        let entry = (touch.ring, touch.index, touch.vertex as usize);
        if !listed.contains(&entry) {
            listed.push(entry);
        }
    }
    Ok(listed)
}

/// One ring's vertex indices, in the order that keeps the polygon on the
/// left: the outer ring counter-clockwise, holes clockwise.
#[derive(Debug, Clone)]
struct Loop {
    vertices: Vec<u32>,
}

/// Exact orientation of `c` against the directed line `a -> b`: `1` left,
/// `-1` right, `0` on it.
///
/// `orient2d` escalates to exact arithmetic, so it always proves a sign.
/// A result it did not prove would read as `0`, which every caller treats
/// as the degenerate, refusing answer: an ear or a bridge is rejected, and
/// the certificate fails rather than passing a triangle it cannot prove.
fn orient(a: Point2, b: Point2, c: Point2) -> i8 {
    match orient2d(a, b, c).sign() {
        Some(Sign::Positive) => 1,
        Some(Sign::Negative) => -1,
        _ => 0,
    }
}

/// Whether `c`, collinear with `a` and `b`, lies on the closed segment.
fn within(a: Point2, b: Point2, c: Point2) -> bool {
    a.x.min(b.x) <= c.x && c.x <= a.x.max(b.x) && a.y.min(b.y) <= c.y && c.y <= a.y.max(b.y)
}

/// Whether the closed segments `p1 p2` and `q1 q2` share any point.
fn segments_touch(p1: Point2, p2: Point2, q1: Point2, q2: Point2) -> bool {
    let o1 = orient(p1, p2, q1);
    let o2 = orient(p1, p2, q2);
    let o3 = orient(q1, q2, p1);
    let o4 = orient(q1, q2, p2);
    if o1 * o2 < 0 && o3 * o4 < 0 {
        return true;
    }
    (o1 == 0 && within(p1, p2, q1))
        || (o2 == 0 && within(p1, p2, q2))
        || (o3 == 0 && within(q1, q2, p1))
        || (o4 == 0 && within(q1, q2, p2))
}

/// Whether direction `v -> p` leaves `v` into the polygon, for a vertex `v`
/// between `a` and `b` with the polygon on the left of `a -> v -> b`.
///
/// A direction along either edge is outside: it overlaps the boundary.
fn locally_inside(a: Point2, v: Point2, b: Point2, p: Point2) -> bool {
    let left_of_in = orient(a, v, p) > 0;
    let left_of_out = orient(v, b, p) > 0;
    if orient(a, v, b) > 0 {
        left_of_in && left_of_out
    } else {
        left_of_in || left_of_out
    }
}

/// Refuse a triangulation that does not tile the polygon exactly once.
///
/// `loops` are the boundary loops, the polygon on their left, bounding
/// `outers` connected parts with `holes` holes among them; a polygon with
/// `n` loop vertices then has `n + 2 holes - 2 outers` triangles.
fn certify(
    points: &[Point2],
    loops: &[Loop],
    triangles: &[[u32; 3]],
    outers: usize,
    holes: usize,
) -> GeomResult<()> {
    let fail = |why: String| {
        Err(GeomError::Degenerate(format!(
            "profile triangulation failed its certificate: {why}"
        )))
    };
    let vertex_count: usize = loops.iter().map(|ring| ring.vertices.len()).sum();
    let expected = (vertex_count + 2 * holes).saturating_sub(2 * outers);
    if triangles.len() != expected {
        let parts = if outers == 1 {
            String::new()
        } else {
            format!(" in {outers} parts")
        };
        return fail(format!(
            "{} triangles, a polygon with {vertex_count} vertices and {holes} holes{parts} has {expected}",
            triangles.len(),
        ));
    }
    let mut directed: HashMap<(u32, u32), u32> = HashMap::with_capacity(3 * triangles.len());
    for t in triangles {
        let corner = |k: usize| points[t[k] as usize];
        if orient(corner(0), corner(1), corner(2)) <= 0 {
            return fail(format!("triangle {t:?} is not strictly counter-clockwise"));
        }
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((a, b)).or_default() += 1;
        }
    }
    let mut boundary: HashMap<(u32, u32), ()> = HashMap::with_capacity(vertex_count);
    for ring in loops {
        let n = ring.vertices.len();
        for k in 0..n {
            let edge = (ring.vertices[k], ring.vertices[(k + 1) % n]);
            if directed.get(&edge) != Some(&1) || directed.contains_key(&(edge.1, edge.0)) {
                return fail(format!(
                    "ring edge {}->{} is not covered exactly once from inside",
                    edge.0, edge.1
                ));
            }
            boundary.insert(edge, ());
        }
    }
    for (&(a, b), &count) in &directed {
        if count != 1 {
            return fail(format!("edge {a}->{b} is used {count} times"));
        }
        if !boundary.contains_key(&(a, b)) && directed.get(&(b, a)) != Some(&1) {
            return fail(format!("interior edge {a}->{b} has no twin"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
        vec![
            Point2::new(x0, y0),
            Point2::new(x1, y0),
            Point2::new(x1, y1),
            Point2::new(x0, y1),
        ]
    }

    /// The issue's profile, laid out as `outer ++ holes` with the holes
    /// clockwise.
    fn issue_profile() -> (Vec<Point2>, Vec<Loop>) {
        let mut points = rect(-2.0, -2.0, 2.0, 2.0);
        for mut hole in [rect(-1.5, -0.5, -0.5, 0.5), rect(0.5, -0.5, 1.5, 0.5)] {
            hole.reverse();
            points.extend(hole);
        }
        let loops = vec![
            Loop {
                vertices: vec![0, 1, 2, 3],
            },
            Loop {
                vertices: vec![4, 5, 6, 7],
            },
            Loop {
                vertices: vec![8, 9, 10, 11],
            },
        ];
        (points, loops)
    }

    #[test]
    fn the_certificate_refuses_earcuts_t_junction() {
        // earcut's output for the issue's profile: right area, 12 triangles
        // instead of 14, and edge 10 -> 7 running past vertices 11 and 6.
        let (points, loops) = issue_profile();
        let earcut = [
            [6, 11, 8],
            [10, 7, 0],
            [10, 0, 1],
            [10, 1, 2],
            [3, 0, 7],
            [3, 7, 4],
            [3, 4, 5],
            [5, 6, 8],
            [9, 10, 2],
            [9, 2, 3],
            [3, 5, 8],
            [3, 8, 9],
        ];
        assert!(reason(&points, &loops, &earcut).contains("12 triangles"));
    }

    /// The certificate's refusal message, or a panic if it passed.
    fn reason(points: &[Point2], loops: &[Loop], triangles: &[[u32; 3]]) -> String {
        match certify(points, loops, triangles, 1, loops.len() - 1) {
            Err(GeomError::Degenerate(message)) => message,
            other => panic!("expected a certificate refusal, got {other:?}"),
        }
    }

    #[test]
    fn the_certificate_names_each_defect() {
        let (points, _) = issue_profile();
        let square = &points[..4];
        let loops = vec![Loop {
            vertices: vec![0, 1, 2, 3],
        }];
        assert!(certify(square, &loops, &[[0, 1, 2], [0, 2, 3]], 1, 0).is_ok());
        assert!(reason(square, &loops, &[[0, 2, 1], [0, 3, 2]])
            .contains("not strictly counter-clockwise"));
        // Right count, wrong cover: the second triangle overlaps the first.
        assert!(reason(square, &loops, &[[0, 1, 2], [1, 2, 3]])
            .contains("ring edge 1->2 is not covered exactly once"));

        // A convex hexagon: right count, every ring edge once, but the
        // middle triangle overlaps the three corner ones.
        let hexagon: Vec<Point2> = (0..6)
            .map(|k| {
                let t = f64::from(k) * core::f64::consts::FRAC_PI_3;
                Point2::new(t.cos(), t.sin())
            })
            .collect();
        let loops = vec![Loop {
            vertices: (0..6).collect(),
        }];
        assert!(certify(
            &hexagon,
            &loops,
            &[[0, 1, 2], [2, 3, 4], [4, 5, 0], [0, 2, 4]],
            1,
            0
        )
        .is_ok());
        assert!(reason(
            &hexagon,
            &loops,
            &[[0, 1, 2], [2, 3, 4], [4, 5, 0], [1, 3, 5]]
        )
        .contains("has no twin"));
    }

    #[test]
    fn the_issue_profile_triangulates_and_certifies() {
        let (points, _) = issue_profile();
        let (outer, holes) = (
            points[..4].to_vec(),
            vec![points[4..8].to_vec(), points[8..].to_vec()],
        );
        let (_, triangles) =
            triangulate_rings(&outer, &holes, PinchPolicy::Refuse).expect("triangulates");
        assert_eq!(triangles.len(), 14);
    }
}
