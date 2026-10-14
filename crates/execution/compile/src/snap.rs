//! A boolean's operands moved onto each other's faces by a rounding
//! residue before the mesh boolean cuts them (#276, narrowed in #291,
//! reaching the placement's rounding in #312).
//!
//! # Why
//!
//! A tool that stops a rounding error short of its host's face (an
//! opening exported `4.5e-15` short of the wall it crosses) leaves a skin
//! that thin, and a tool that reaches the same distance past a face leaves
//! a sliver. The mesh boolean decides on the `f64` numbers it is given, so
//! it keeps both: two faces `4.5e-15` apart. They measure correctly where
//! they are, but once the mesh is placed at georeferenced coordinates
//! (about `5.6e6`, where a coordinate rounds to about `1e-9`) the skin's two
//! faces cross, and the solid's volume is refused as self-intersecting.
//!
//! # Reach: rounding, not the tolerance
//!
//! Only a residue at the rounding scale of the coordinates is closed: of
//! the operands' own, or of those their result is placed at, whichever is
//! coarser. A coordinate moves by at most
//!
//! ```text
//! reach = min(eps / sqrt 3, REACH_EPSILONS * f64::EPSILON * S),   REACH_EPSILONS = 16
//! S     = max(S_operands, S_placed),   S_placed = P + S_operands  (P > 0)
//! ```
//!
//! where `eps` is the linear tolerance, so a vertex moves by at most
//! `sqrt 3 * reach <= eps`. `S_operands` is the largest coordinate
//! magnitude of either operand. `P` is the placed scale ([`placed_scale`],
//! `0` where nothing places the boolean): the world coordinate magnitude
//! the instances above it add, in the operands' own units, so `S_placed`
//! bounds the largest world coordinate magnitude the result reaches once
//! placed, in those units (up to the factor, at most `sqrt 3`, by which a
//! turn grows a box's largest coordinate, which the sixteen epsilons
//! absorb).
//!
//! A boolean computed in its host's own frame (`S_operands` about `3.67`
//! for a wall) and placed by an `Instance` `5.6e6` out is cut where a
//! coordinate rounds to `1e-15`, but measured where it rounds to `1e-9`:
//! a skin of `1e-12` is far above its operands' rounding and far below its
//! placement's, so it is kept, and its two faces cross once placed, and
//! the volume is refused as self-intersecting (#312). The reach is
//! therefore that of the placed coordinates. The compiler carries `P` down
//! through every `Instance` and `InstanceAtStation` with the tolerance and
//! the chord budget (`instance_local_budget`), and keys its cache on it,
//! so a source placed near and far is compiled once per scale. Sixteen relative epsilons are about 8 to 16
//! units in the last place of `S`: the rounding a coordinate accrues
//! through a few composed placements (each affine step at most about
//! `3 f64::EPSILON * S`) and a 15-significant-digit export (about
//! `4.5 f64::EPSILON` relative). The #276 skin, `4.5e-15` in a wall whose
//! coordinates reach `3.67`, is `5.5 f64::EPSILON * S`, inside the reach
//! (`1.3e-14`); at georeferenced coordinates (`S` about `6e6`, the
//! operands' own or their placement's) the reach is about `2e-8`, so a
//! residue of an ulp or a few there (`1e-9`) closes, and so does any
//! residue up to that of a wall cut in its own frame and placed there;
//! a skin of `1e-6` (fifty times the reach) is kept. A skin or sliver
//! above the reach is authored geometry: a `0.1`-`1 mm` skin under a
//! `1 mm` tolerance stays in the mesh (#291), however thin its tolerance
//! would allow it to be read.
//!
//! # What moves: an exact landing only
//!
//! A move by a rounding residue helps only if the vertex then lies exactly
//! on the face: one that lands within rounding of it is as likely to land
//! on its other side, and a face whose vertices land on both sides crosses
//! the face it should meet. The #276 snap did this to a real slab (#291):
//! its openings' bottoms were moved onto planes it computed through a
//! normalised normal (`0.9999999999999999`, not `1`) or that were tilted by a
//! rounded half turn, some vertices landing above the slab's bottom face
//! and some below, and the mesh boolean left a void tangent to that face,
//! refused. So the only faces anything is moved onto are those whose three
//! corners share one coordinate exactly (an axis-aligned plane `x_k = c`,
//! as every face of a host built in its own frame is), and the moved
//! coordinate becomes `c` itself, exactly. A face that is axis-aligned only
//! within rounding (a frame turned by a rounded half turn, a profile
//! edge's end evaluated as `a + (b - a)`) is cut as given.
//!
//! A coordinate moves by its value: where a vertex of the tool lies within
//! the reach of such a face of the subject (of its plane, and of its
//! triangle), every vertex of the tool with that coordinate value is moved
//! onto the plane, so a face of the tool on an axis-aligned plane moves
//! whole and stays planar, however little of it lies near the subject. A
//! vertex near two or three such planes (on an edge or a corner of the
//! subject) lands on all of them. Then the subject's coordinates are moved
//! onto the moved tool's faces the same way, so a large tool whose own
//! corners lie away from the subject (a slab cut off a whole face) is
//! settled from the subject's side. A value within the reach of two
//! parallel planes, or on one and near another (a skin of the other
//! operand thinner than the reach), stays.
//!
//! Each operand mesh is then the given one with every vertex moved by at
//! most `sqrt 3 * reach <= eps`, so every point of it (a convex
//! combination of its triangle's corners) by at most that: the mesh
//! boolean's result is the boolean of operands perturbed within the
//! tolerance, the reading ADR 0080 gives the exact boolean.
//!
//! # Never silently, never worse
//!
//! The largest move is reported: [`crate::deviation`] carries it as a
//! [`crate::DeviationBound::Certified`] contribution of the boolean path
//! under the detail [`crate::deviation::SNAPPED_OPERANDS`], through
//! instances (scaled) and enclosing booleans. At a linear tolerance of
//! zero nothing moves. A move that would turn or flatten a triangle of the
//! operand (a feature of its own thinner than the reach) is not made: that
//! operand is cut as given. And a snap never turns a working boolean into
//! a refusal: where the snapped operands' boolean is refused, or its result
//! touches itself, the compiler cuts the operands as given, and when that
//! succeeds it is the result, with no snap reported (#291,
//! `ReferenceMeshCompiler::build_boolean`).

use std::collections::HashMap;

use axiolid_core::{Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_mesh::TriMesh;

/// The reach in relative epsilons of the largest coordinate magnitude
/// (see the module notes).
pub(crate) const REACH_EPSILONS: Scalar = 16.0;

/// The operands after [`settle`]: each moved one, and the largest distance
/// any vertex moved.
#[derive(Debug, Clone)]
pub(crate) struct Settled {
    /// The subject, if any of its vertices moved.
    pub subject: Option<TriMesh>,
    /// The tool, if any of its vertices moved.
    pub tool: Option<TriMesh>,
    /// The largest move, `0` when nothing moved.
    pub moved: Scalar,
}

/// Significant bits [`placed_scale`] keeps (the leading one and two more):
/// a placed scale is rounded up to `1`, `1.25`, `1.5` or `1.75` times a
/// power of two, so a nested placement grows it by at most a quarter per
/// level, and instances at similar distances share a cache entry.
const PLACED_BITS: u32 = 3;

/// The placed scale of an instance's source: the world coordinate
/// magnitude `parent`, the placed scale of the instance itself, and
/// `transform`'s translation add, in the source's units (divided by the
/// transform's smallest column length, so a world rounding residue is never
/// read smaller than it is), rounded up to [`PLACED_BITS`] significant
/// bits. `None` when it is not finite. See the module notes (#312).
pub(crate) fn placed_scale(transform: Transform3, parent: Scalar) -> Option<Scalar> {
    let m = transform.matrix3;
    let columns = [m.x_axis, m.y_axis, m.z_axis].map(Vec3::length);
    let smallest = columns.into_iter().fold(Scalar::INFINITY, Scalar::min);
    let unit = if smallest > 0.0 {
        smallest
    } else {
        crate::deviation::stretch(transform)
    };
    let offset = parent + transform.translation.abs().max_element();
    let placed = offset / unit;
    if !placed.is_finite() || placed < 0.0 {
        return None;
    }
    let mask = (1_u64 << (52 - (PLACED_BITS - 1))) - 1;
    let bits = placed.to_bits();
    let rounded = if bits & mask == 0 {
        placed
    } else {
        Scalar::from_bits((bits | mask) + 1)
    };
    rounded.is_finite().then_some(rounded)
}

/// How far a coordinate of `subject` or `tool` may move: the rounding scale
/// of their coordinates, or of those their result is placed at (`placed`,
/// [`placed_scale`]), whichever is coarser, capped so that a vertex moves
/// by at most the linear tolerance, `0` at zero tolerance (see the module
/// notes).
pub(crate) fn reach(
    subject: &TriMesh,
    tool: &TriMesh,
    tolerance: Tolerance,
    placed: Scalar,
) -> Scalar {
    let eps = tolerance.linear();
    if eps <= 0.0 || !eps.is_finite() {
        return 0.0;
    }
    let operands = subject
        .positions
        .iter()
        .chain(&tool.positions)
        .map(|p| p.abs().max_element())
        .filter(|m| m.is_finite())
        .fold(0.0, Scalar::max);
    let placed = if placed > 0.0 && placed.is_finite() {
        placed + operands
    } else {
        0.0
    };
    let scale = operands.max(placed);
    (REACH_EPSILONS * Scalar::EPSILON * scale).min(eps / 3.0_f64.sqrt())
}

/// Move the tool's coordinates onto the subject's axis-aligned faces, then
/// the subject's onto the moved tool's, each by at most [`reach`] (see the
/// module notes); `placed` is the boolean's placed scale.
pub(crate) fn settle(
    subject: &TriMesh,
    tool: &TriMesh,
    tolerance: Tolerance,
    placed: Scalar,
) -> Settled {
    let reach = reach(subject, tool, tolerance, placed);
    let mut out = Settled {
        subject: None,
        tool: None,
        moved: 0.0,
    };
    if reach <= 0.0 {
        return out;
    }
    if let Some((moved_tool, by)) = onto(tool, subject, reach) {
        out.moved = out.moved.max(by);
        out.tool = Some(moved_tool);
    }
    let fixed = out.tool.as_ref().unwrap_or(tool);
    if let Some((moved_subject, by)) = onto(subject, fixed, reach) {
        out.moved = out.moved.max(by);
        out.subject = Some(moved_subject);
    }
    out
}

/// An axis-aligned triangle of the fixed operand: its plane `x_axis = at`,
/// and its box grown by the reach.
struct Face {
    corners: [Point3; 3],
    axis: usize,
    at: Scalar,
    lo: Point3,
    hi: Point3,
}

/// The axis whose coordinate all three corners share exactly, if the
/// triangle is not degenerate.
fn axis_of(corners: &[Point3; 3]) -> Option<usize> {
    let normal = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
    if normal == Vec3::ZERO || !normal.is_finite() {
        return None;
    }
    (0..3).find(|&k| corners[0][k] == corners[1][k] && corners[1][k] == corners[2][k])
}

/// `moving` with its coordinates within `reach` of the axis-aligned faces
/// of `fixed` set onto them, and the largest move; `None` when nothing
/// moved, or a move would turn or flatten one of `moving`'s triangles.
///
/// A coordinate moves by its value: where a vertex of `moving` lies within
/// `reach` of a face `x_k = c` of `fixed`, every vertex of `moving` whose
/// `k`-th coordinate has that same value gets `c`. So a face of `moving`
/// on an axis-aligned plane stays planar however little of it lies near
/// `fixed` (a door that reaches a micrometre below its wall's floor has
/// its end moved whole, not only the corners above the floor). A value
/// within `reach` of two planes `x_k = c` (a skin of `fixed` thinner than
/// the reach), or of a plane and already on another, stays.
fn onto(moving: &TriMesh, fixed: &TriMesh, reach: Scalar) -> Option<(TriMesh, Scalar)> {
    let (mlo, mhi) = bounds(&moving.positions)?;
    let grow = Vec3::splat(reach);
    let faces: Vec<Face> = fixed
        .indices
        .chunks_exact(3)
        .filter_map(|t| {
            let corners = [
                *fixed.positions.get(t[0] as usize)?,
                *fixed.positions.get(t[1] as usize)?,
                *fixed.positions.get(t[2] as usize)?,
            ];
            let lo = corners[0].min(corners[1]).min(corners[2]) - grow;
            let hi = corners[0].max(corners[1]).max(corners[2]) + grow;
            // Only faces the moving operand can reach.
            if lo.cmpgt(mhi).any() || hi.cmplt(mlo).any() {
                return None;
            }
            let axis = axis_of(&corners)?;
            Some(Face {
                corners,
                axis,
                at: corners[0][axis],
                lo,
                hi,
            })
        })
        .collect();
    if faces.is_empty() {
        return None;
    }
    // Per axis, the plane each coordinate value is held to: `None` once
    // two planes claim it.
    let mut held: [HashMap<u64, Option<Scalar>>; 3] = Default::default();
    let key = |v: Scalar| (v + 0.0).to_bits();
    for p in &moving.positions {
        for face in &faces {
            let k = face.axis;
            if p.cmplt(face.lo).any()
                || p.cmpgt(face.hi).any()
                || (face.at - p[k]).abs() > reach
                || distance_to_triangle(*p, &face.corners) > reach
            {
                continue;
            }
            held[k]
                .entry(key(p[k]))
                .and_modify(|at| {
                    if *at != Some(face.at) {
                        *at = None;
                    }
                })
                .or_insert(Some(face.at));
        }
    }
    let mut positions = moving.positions.clone();
    let mut moved: Scalar = 0.0;
    for p in &mut positions {
        let mut to = *p;
        for (k, held) in held.iter().enumerate() {
            if let Some(Some(at)) = held.get(&key(p[k])) {
                to[k] = *at;
            }
        }
        moved = moved.max((to - *p).length());
        *p = to;
    }
    if moved == 0.0 {
        return None;
    }
    // A move must keep every triangle facing the way it did, and keep its
    // area: an operand with a feature of its own thinner than the reach is
    // cut as given rather than folded.
    for t in moving.indices.chunks_exact(3) {
        let [a, b, c] = [t[0], t[1], t[2]].map(|i| i as usize);
        let before = (moving.positions[b] - moving.positions[a])
            .cross(moving.positions[c] - moving.positions[a]);
        let after = (positions[b] - positions[a]).cross(positions[c] - positions[a]);
        if before != after && after.dot(before) <= 0.0 {
            return None;
        }
    }
    let mut out = moving.clone();
    out.positions = positions;
    Some((out, moved))
}

/// The distance from `p` to the closed triangle `t`.
fn distance_to_triangle(p: Point3, t: &[Point3; 3]) -> Scalar {
    // Ericson, Real-Time Collision Detection, 5.1.5.
    let [a, b, c] = *t;
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return ap.length();
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return bp.length();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (p - (a + ab * v)).length();
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return cp.length();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (p - (a + ac * w)).length();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (b + (c - b) * w)).length();
    }
    let denom = va + vb + vc;
    if denom == 0.0 {
        // A degenerate triangle: its edges bound the distance.
        return [(a, b), (b, c), (c, a)]
            .iter()
            .map(|&(s, e)| segment_distance(p, s, e))
            .fold(Scalar::INFINITY, Scalar::min);
    }
    let (v, w) = (vb / denom, vc / denom);
    (p - (a + ab * v + ac * w)).length()
}

fn segment_distance(p: Point3, a: Point3, b: Point3) -> Scalar {
    let ab = b - a;
    let l = ab.length_squared();
    let t = if l > 0.0 {
        ((p - a).dot(ab) / l).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p - (a + ab * t)).length()
}

/// The box of `points`, `None` when empty or not finite.
fn bounds(points: &[Point3]) -> Option<(Point3, Point3)> {
    let first = *points.first()?;
    let (lo, hi) = points
        .iter()
        .fold((first, first), |(lo, hi), p| (lo.min(*p), hi.max(*p)));
    (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube(lo: Point3, hi: Point3) -> TriMesh {
        let p = |x: bool, y: bool, z: bool| {
            Point3::new(
                if x { hi.x } else { lo.x },
                if y { hi.y } else { lo.y },
                if z { hi.z } else { lo.z },
            )
        };
        let positions = vec![
            p(false, false, false),
            p(true, false, false),
            p(true, true, false),
            p(false, true, false),
            p(false, false, true),
            p(true, false, true),
            p(true, true, true),
            p(false, true, true),
        ];
        let indices = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ];
        TriMesh::new(positions, indices)
    }

    fn metre() -> Tolerance {
        Tolerance::METRE
    }

    #[test]
    fn the_reach_is_the_rounding_of_the_coordinates_capped_by_the_tolerance() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(Point3::new(1.0, 0.1, 0.5), Point3::new(2.0, 0.35, 2.5));
        let near = reach(&host, &tool, Tolerance::MILLIMETRE, 0.0);
        assert_eq!(near, REACH_EPSILONS * Scalar::EPSILON * 4.0);
        let far = Vec3::new(6.0e5, 5.6e6, 0.0);
        let placed = |m: &TriMesh| {
            let mut m = m.clone();
            m.positions.iter_mut().for_each(|p| *p += far);
            m
        };
        let georeferenced = reach(&placed(&host), &placed(&tool), Tolerance::MILLIMETRE, 0.0);
        // A few ulps there, about 2e-8: well under a micrometre.
        assert!((1e-9..1e-7).contains(&georeferenced), "{georeferenced}");
        // Never above the tolerance, nothing at zero.
        let tiny = reach(&host, &tool, Tolerance::new(1e-16, 1e-9).unwrap(), 0.0);
        assert_eq!(tiny, 1e-16 / 3.0_f64.sqrt());
        assert_eq!(reach(&host, &tool, Tolerance::ZERO, 0.0), 0.0);
    }

    #[test]
    fn the_reach_covers_the_coordinates_the_result_is_placed_at() {
        // A wall cut in its own frame and placed 5.6e6 out (#312): the
        // reach is that of the placed coordinates, about 2e-8.
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(Point3::new(1.0, 0.1, 0.5), Point3::new(2.0, 0.35, 2.5));
        let site = Transform3::from_translation(Vec3::new(6.0e5, 5.6e6, 204.3));
        let placed = placed_scale(site, 0.0).expect("finite");
        // Rounded up to 1.5 * 2^22, never below the translation.
        assert_eq!(placed, 1.5 * 4_194_304.0);
        let far = reach(&host, &tool, Tolerance::MILLIMETRE, placed);
        assert_eq!(far, REACH_EPSILONS * Scalar::EPSILON * (placed + 4.0));
        assert!((1e-8..3e-8).contains(&far), "{far}");
        // Unplaced, or placed near: the operands' own rounding.
        let near = reach(&host, &tool, Tolerance::MILLIMETRE, 0.0);
        assert_eq!(near, REACH_EPSILONS * Scalar::EPSILON * 4.0);
        // Still capped by the tolerance, and nothing at zero.
        let tiny = Tolerance::new(1e-10, 1e-9).unwrap();
        assert_eq!(reach(&host, &tool, tiny, placed), 1e-10 / 3.0_f64.sqrt());
        assert_eq!(reach(&host, &tool, Tolerance::ZERO, placed), 0.0);
    }

    #[test]
    fn the_placed_scale_composes_in_the_sources_units() {
        let at = |x: Scalar| Transform3::from_translation(Vec3::new(x, -2.0 * x, 0.5));
        // The largest translation component (here twice `x`), exact at a
        // value of three significant bits, rounded up otherwise.
        assert_eq!(placed_scale(at(3.0), 0.0), Some(6.0));
        assert_eq!(placed_scale(at(5.0), 0.0), Some(10.0));
        assert_eq!(placed_scale(at(4.5), 0.0), Some(10.0));
        assert_eq!(placed_scale(Transform3::IDENTITY, 0.0), Some(0.0));
        // A nested placement adds its translation to its parent's scale.
        assert_eq!(placed_scale(at(0.5), 6.0), Some(7.0));
        // A source scaled up is placed at a smaller scale in its own units,
        // read through its smallest column so a world residue is never read
        // smaller than it is: 6e6 / 2, rounded up to 1.5 * 2^21.
        let scaled = Transform3::from_cols(
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.0, 0.0, 8.0),
            Vec3::new(6.0e6, 0.0, 0.0),
        );
        assert_eq!(placed_scale(scaled, 0.0), Some(1.5 * 2_097_152.0));
        // Instances at similar distances share one scale.
        let a = placed_scale(at(5.61e6), 0.0);
        let b = placed_scale(at(5.89e6), 0.0);
        assert_eq!(a, b);
        assert_eq!(placed_scale(at(Scalar::INFINITY), 0.0), None);
    }

    #[test]
    fn a_tool_short_of_a_face_lands_on_it_exactly() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(Point3::new(1.0, 4.5e-15, 0.5), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, metre(), 0.0);
        let moved = s.tool.expect("the tool moves");
        assert!(moved.positions.iter().any(|p| p.y == 0.0));
        assert!(moved.positions.iter().all(|p| p.y == 0.0 || p.y == 0.35));
        assert!(s.moved > 0.0 && s.moved <= 1e-14);
        assert!(s.subject.is_none());
    }

    #[test]
    fn a_vertex_on_an_edge_lands_on_both_planes() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        // Bottom just above the host's floor face, end short of its front:
        // the corner near both lands on their common edge.
        let tool = cube(Point3::new(1.0, 4e-15, 4e-15), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, metre(), 0.0);
        let moved = s.tool.expect("the tool moves");
        assert!(moved
            .positions
            .iter()
            .filter(|p| p.y < 0.1)
            .all(|p| p.y == 0.0 && (p.z == 0.0 || p.z == 2.5)));
        assert!(s.moved > 4e-15 && s.moved < 6e-15, "{}", s.moved);
    }

    #[test]
    fn a_gap_above_the_rounding_stays_however_large_the_tolerance() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        for (gap, tolerance) in [
            (1e-12, metre()),
            (2e-6, metre()),
            (1e-6, Tolerance::MILLIMETRE),
            (5e-4, Tolerance::MILLIMETRE),
        ] {
            let tool = cube(Point3::new(1.0, gap, 0.5), Point3::new(2.0, 0.35, 2.5));
            let s = settle(&host, &tool, tolerance, 0.0);
            assert!(s.tool.is_none() && s.subject.is_none(), "{gap}");
            assert_eq!(s.moved, 0.0);
        }
    }

    #[test]
    fn nothing_moves_at_zero_tolerance() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(Point3::new(1.0, 4.5e-15, 0.5), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, Tolerance::ZERO, 0.0);
        assert!(s.tool.is_none() && s.subject.is_none());
    }

    #[test]
    fn a_large_tool_settles_the_subject_onto_its_face() {
        let host = cube(Point3::ZERO, Point3::ONE);
        // Cuts the top half off all but a skin a rounding error thick.
        let top = 1.0 - 4.0 * Scalar::EPSILON;
        let tool = cube(Point3::new(-1.0, -1.0, 0.5), Point3::new(2.0, 2.0, top));
        let s = settle(&host, &tool, metre(), 0.0);
        assert!(s.tool.is_none());
        let subject = s.subject.expect("the subject moves");
        assert!(subject.positions.iter().all(|p| p.z == 0.0 || p.z == top));
    }

    #[test]
    fn a_corner_near_three_faces_lands_on_their_corner() {
        let host = cube(Point3::ZERO, Point3::splat(2.0));
        let tool = cube(Point3::splat(4.5e-15), Point3::ONE);
        let s = settle(&host, &tool, metre(), 0.0);
        let moved = s.tool.expect("the tool moves");
        assert!(moved.positions.contains(&Point3::ZERO));
        assert!(s.moved <= 3.0_f64.sqrt() * reach(&host, &tool, metre(), 0.0));
    }

    #[test]
    fn a_face_mostly_beyond_the_face_it_meets_moves_whole() {
        // The tool stops a rounding error short of the host's front face
        // and reaches a micrometre below its floor: its bottom corners are
        // a micrometre from every face of the host, its top ones on the
        // front face. Its end moves whole, and stays planar (#291).
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(
            Point3::new(1.0, 4.5e-15, -1e-6),
            Point3::new(2.0, 0.35, 2.5),
        );
        let s = settle(&host, &tool, metre(), 0.0);
        let moved = s.tool.expect("the tool moves");
        assert!(moved.positions.iter().all(|p| p.y == 0.0 || p.y == 0.35));
        assert!(moved.positions.iter().any(|p| p.z == -1e-6));
    }

    #[test]
    fn a_tool_thinner_than_the_reach_is_not_flattened() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        // Both of its faces lie within the reach of the host's face:
        // moving both onto it would fold the tool flat.
        let tool = cube(Point3::new(1.0, -6e-15, 0.5), Point3::new(2.0, -2e-15, 2.5));
        let s = settle(&host, &tool, metre(), 0.0);
        assert!(s.tool.is_none());
    }

    #[test]
    fn a_vertex_between_parallel_planes_within_the_reach_stays() {
        // A subject slab thinner than the reach: no single plane to land on.
        let host = cube(Point3::ZERO, Point3::new(4.0, 1e-14, 3.0));
        let tool = cube(Point3::new(1.0, 5e-15, 0.5), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, metre(), 0.0);
        assert!(s.tool.is_none());
    }

    #[test]
    fn a_face_axis_aligned_only_within_rounding_is_not_landed_on() {
        // The host's front face runs from `y = 0` to `y = -1.4e-15` along
        // `x` (a profile edge's end evaluated as `a + (b - a)`): a tool
        // flush with one end of it, or a rounding error off, cannot land
        // exactly on its plane, and is cut as given (#291).
        let mut host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        for p in &mut host.positions {
            if p.x == 4.0 && p.y == 0.0 {
                p.y = -1.4e-15;
            }
        }
        for y in [0.0, 2e-15, -2e-15] {
            let tool = cube(Point3::new(1.0, y, 0.5), Point3::new(2.0, 0.35, 2.5));
            let s = settle(&host, &tool, metre(), 0.0);
            assert!(s.tool.is_none() && s.subject.is_none(), "{y}");
        }
    }
}
