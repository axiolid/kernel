//! A boolean's operands moved onto each other's faces within the linear
//! tolerance before the mesh boolean cuts them (#276).
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
//! # What moves, and by how much
//!
//! Every vertex of the tool within the linear tolerance `eps` of a
//! triangle of the subject is moved onto that triangle's plane; a vertex
//! near two or three such planes (on an edge or a corner of the subject)
//! onto their common line or point, by the smallest move that puts it on
//! all of them. Then every vertex of the subject within `eps` of the moved
//! tool is moved onto the tool's planes the same way, so a large tool
//! whose own corners lie away from the subject (a slab cut off a whole
//! face) is settled from the subject's side. A vertex is moved only when
//! the move is at most `eps`; one between two parallel planes within `eps`
//! (a skin of its own operand), or whose planes do not meet near it, stays.
//!
//! Each operand mesh is then the given one with every vertex moved by at
//! most `eps`, so every point of it (a convex combination of its
//! triangle's corners) by at most `eps`: the mesh boolean's result is the
//! boolean of operands perturbed within the tolerance, the reading ADR
//! 0080 gives the exact boolean. A skin or sliver thinner than `eps`
//! vanishes; one thicker than `eps` is out of reach and kept.
//!
//! # Exactly, where the faces allow it
//!
//! The move is computed so that an axis-aligned plane receives the vertex
//! exactly: the plane's unit normal is then exact, the residue `c - n . v`
//! is an exact `f64` difference of two close numbers, and the corrected
//! coordinate is the plane's own. A host built in its own frame, as a wall
//! with its openings subtracted before it is placed, is axis-aligned, so
//! the faces become exactly coplanar and the mesh boolean's coplanar rule
//! removes the skin. A plane under a general rotation has no `f64` points
//! exactly on it; the vertex then lands within rounding of it, which the
//! mesh boolean may still read as a skin that thin.
//!
//! # Never silently, never beyond the tolerance
//!
//! The largest move is reported: [`crate::deviation`] carries it as a
//! [`crate::DeviationBound::Certified`] contribution of the boolean path
//! under the detail [`crate::deviation::SNAPPED_OPERANDS`], through
//! instances (scaled) and enclosing booleans. At a linear tolerance of
//! zero nothing moves. A move that would turn or flatten a triangle of the
//! operand (a feature of its own thinner than `eps`) is not made: that
//! operand is cut as given.

use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;

/// Two unit normals are one direction, or opposite, within this sine.
/// Bookkeeping, not a decision: it only tells whether a second plane near
/// a vertex adds a constraint or repeats one.
const PARALLEL: Scalar = 1e-12;

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

/// Move the tool's vertices onto the subject's faces, then the subject's
/// onto the moved tool's, each by at most the linear tolerance (see the
/// module notes).
pub(crate) fn settle(subject: &TriMesh, tool: &TriMesh, tolerance: Tolerance) -> Settled {
    let eps = tolerance.linear();
    let mut out = Settled {
        subject: None,
        tool: None,
        moved: 0.0,
    };
    if eps <= 0.0 || !eps.is_finite() {
        return out;
    }
    if let Some((moved_tool, by)) = onto(tool, subject, eps) {
        out.moved = out.moved.max(by);
        out.tool = Some(moved_tool);
    }
    let fixed = out.tool.as_ref().unwrap_or(tool);
    if let Some((moved_subject, by)) = onto(subject, fixed, eps) {
        out.moved = out.moved.max(by);
        out.subject = Some(moved_subject);
    }
    out
}

/// A triangle of the fixed operand: its plane and its box grown by `eps`.
struct Face {
    corners: [Point3; 3],
    normal: Vec3,
    offset: Scalar,
    lo: Point3,
    hi: Point3,
}

/// `moving` with every vertex within `eps` of `fixed` moved onto it, and
/// the largest move; `None` when no vertex moved, or a move would turn or
/// flatten one of `moving`'s triangles.
fn onto(moving: &TriMesh, fixed: &TriMesh, eps: Scalar) -> Option<(TriMesh, Scalar)> {
    let (mlo, mhi) = bounds(&moving.positions)?;
    let grow = Vec3::splat(eps);
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
            let normal = (corners[1] - corners[0])
                .cross(corners[2] - corners[0])
                .try_normalize()?;
            Some(Face {
                corners,
                normal,
                offset: normal.dot(corners[0]),
                lo,
                hi,
            })
        })
        .collect();
    if faces.is_empty() {
        return None;
    }
    // Only vertices inside the reachable faces' (grown) box can move.
    let (flo, fhi) = faces
        .iter()
        .fold((faces[0].lo, faces[0].hi), |(lo, hi), f| {
            (lo.min(f.lo), hi.max(f.hi))
        });
    let mut positions = moving.positions.clone();
    let mut moved: Scalar = 0.0;
    for p in &mut positions {
        if p.cmplt(flo).any() || p.cmpgt(fhi).any() {
            continue;
        }
        if let Some(to) = target(*p, &faces, eps) {
            moved = moved.max((to - *p).length());
            *p = to;
        }
    }
    if moved == 0.0 {
        return None;
    }
    // A move must keep every triangle facing the way it did, and keep its
    // area: an operand with a feature of its own thinner than `eps` is cut
    // as given rather than folded.
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

/// Where `p` moves: onto every distinct plane of a face within `eps` of
/// it, by the smallest move, when that move is positive and at most `eps`.
fn target(p: Point3, faces: &[Face], eps: Scalar) -> Option<Point3> {
    // Up to three independent planes: `(normal, residue c - n . p)`.
    let mut planes: Vec<(Vec3, Scalar)> = Vec::with_capacity(3);
    for face in faces {
        if p.cmplt(face.lo).any() || p.cmpgt(face.hi).any() {
            continue;
        }
        let residue = face.offset - face.normal.dot(p);
        if residue.abs() > eps || distance_to_triangle(p, &face.corners) > eps {
            continue;
        }
        let mut repeated = false;
        for &(n, r) in &planes {
            if n.cross(face.normal).length() <= PARALLEL {
                // The same plane again (a neighbouring triangle), or a
                // parallel one: a vertex between two parallel planes
                // within `eps` stays where it is.
                let same = if n.dot(face.normal) > 0.0 {
                    (r - residue).abs()
                } else {
                    (r + residue).abs()
                };
                if same > PARALLEL * (1.0 + p.abs().max_element()) {
                    return None;
                }
                repeated = true;
                break;
            }
        }
        if repeated {
            continue;
        }
        if planes.len() == 2 {
            // A third plane through the line of the first two adds nothing
            // the line does not satisfy (checked below), or is not met.
            let (a, b) = (planes[0].0, planes[1].0);
            if a.cross(b).dot(face.normal).abs() <= PARALLEL {
                continue;
            }
        }
        if planes.len() == 3 {
            continue;
        }
        planes.push((face.normal, residue));
    }
    if planes.iter().all(|&(_, r)| r == 0.0) {
        return None;
    }
    let step = smallest_move(&planes)?;
    if step.length() > eps {
        return None;
    }
    let to = p + step;
    // Every plane met within `eps`, the dependent ones included.
    for face in faces {
        if p.cmplt(face.lo).any() || p.cmpgt(face.hi).any() {
            continue;
        }
        let residue = face.offset - face.normal.dot(p);
        if residue.abs() <= eps
            && distance_to_triangle(p, &face.corners) <= eps
            && (face.offset - face.normal.dot(to)).abs() > eps
        {
            return None;
        }
    }
    Some(to)
}

/// The smallest `d` with `n_i . d = r_i` for every plane: `d = sum l_i n_i`
/// with `G l = r`, `G` the normals' Gram matrix. For orthogonal unit
/// normals `G` is exactly the identity and `d` the sum of the exact
/// residues along their exact axes.
fn smallest_move(planes: &[(Vec3, Scalar)]) -> Option<Vec3> {
    let k = planes.len();
    let mut g = [[0.0; 4]; 3];
    for i in 0..k {
        for j in 0..k {
            g[i][j] = planes[i].0.dot(planes[j].0);
        }
        g[i][3] = planes[i].1;
    }
    // Gaussian elimination with partial pivoting on the small system.
    for col in 0..k {
        let pivot = (col..k).max_by(|&a, &b| g[a][col].abs().total_cmp(&g[b][col].abs()))?;
        if g[pivot][col].abs() <= PARALLEL {
            return None;
        }
        g.swap(col, pivot);
        for row in 0..k {
            if row != col && g[row][col] != 0.0 {
                let f = g[row][col] / g[col][col];
                for c in col..4 {
                    g[row][c] -= f * g[col][c];
                }
            }
        }
    }
    let mut d = Vec3::ZERO;
    for (i, &(n, _)) in planes.iter().enumerate() {
        d += n * (g[i][3] / g[i][i]);
    }
    d.is_finite().then_some(d)
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
    fn a_tool_short_of_a_face_lands_on_it_exactly() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(Point3::new(1.0, 4.5e-15, 0.5), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, metre());
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
        let tool = cube(Point3::new(1.0, 1e-12, 1e-12), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, metre());
        let moved = s.tool.expect("the tool moves");
        assert!(moved
            .positions
            .iter()
            .filter(|p| p.y < 0.1)
            .all(|p| p.y == 0.0 && (p.z == 0.0 || p.z == 2.5)));
        assert!(s.moved > 1e-12 && s.moved < 2e-12);
    }

    #[test]
    fn a_gap_beyond_the_tolerance_stays() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(Point3::new(1.0, 2e-6, 0.5), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, metre());
        assert!(s.tool.is_none() && s.subject.is_none());
        assert_eq!(s.moved, 0.0);
    }

    #[test]
    fn nothing_moves_at_zero_tolerance() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        let tool = cube(Point3::new(1.0, 4.5e-15, 0.5), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, Tolerance::ZERO);
        assert!(s.tool.is_none() && s.subject.is_none());
    }

    #[test]
    fn a_large_tool_settles_the_subject_onto_its_face() {
        let host = cube(Point3::ZERO, Point3::ONE);
        // Cuts the top half off all but a skin `1e-13` thick.
        let tool = cube(
            Point3::new(-1.0, -1.0, 0.5),
            Point3::new(2.0, 2.0, 1.0 - 1e-13),
        );
        let s = settle(&host, &tool, metre());
        assert!(s.tool.is_none());
        let subject = s.subject.expect("the subject moves");
        assert!(subject
            .positions
            .iter()
            .all(|p| p.z == 0.0 || p.z == 1.0 - 1e-13));
    }

    #[test]
    fn a_corner_further_than_the_tolerance_stays() {
        let host = cube(Point3::ZERO, Point3::splat(2.0));
        // Its corner is within the tolerance of three faces of the host,
        // but the host's corner is `0.9e-6 * sqrt 3` away.
        let tool = cube(Point3::splat(0.9e-6), Point3::ONE);
        let s = settle(&host, &tool, metre());
        let moved = s.tool.expect("the vertices near one face move");
        assert!(moved.positions.contains(&Point3::splat(0.9e-6)));
        assert!(s.moved <= 1e-6, "{}", s.moved);
    }

    #[test]
    fn a_tool_thinner_than_the_tolerance_is_not_flattened() {
        let host = cube(Point3::ZERO, Point3::new(4.0, 0.25, 3.0));
        // Both of its faces lie within the tolerance of the host's face:
        // moving both onto it would fold the tool flat.
        let tool = cube(Point3::new(1.0, -6e-7, 0.5), Point3::new(2.0, -2e-7, 2.5));
        let s = settle(&host, &tool, metre());
        assert!(s.tool.is_none());
    }

    #[test]
    fn a_vertex_between_parallel_planes_within_tolerance_stays() {
        // A subject slab thinner than the tolerance: no single plane to
        // land on.
        let host = cube(Point3::ZERO, Point3::new(4.0, 1e-7, 3.0));
        let tool = cube(Point3::new(1.0, 5e-8, 0.5), Point3::new(2.0, 0.35, 2.5));
        let s = settle(&host, &tool, metre());
        assert!(s.tool.is_none());
    }
}
