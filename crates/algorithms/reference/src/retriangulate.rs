//! Retriangulating a face against the intersection curve crossing it.
//!
//! # What this produces
//!
//! A face cut by the intersection curve is replaced by triangles whose
//! edges follow that curve. Every output triangle then lies wholly inside or
//! wholly outside the other solid, so classification becomes a per-triangle
//! question with no further geometry -- which is what makes an exact boolean
//! possible.
//!
//! # Why the work happens in 2D
//!
//! All points involved lie in the face's plane by construction: the face's
//! own corners, and curve nodes that were computed as crossings OF that
//! plane. Projecting along the plane's dominant axis is therefore exact in
//! the sense that matters -- it drops a coordinate that carries no
//! information, rather than approximating one that does.
//!
//! The dominant axis is chosen from the largest normal component so the
//! projection never collapses: picking a near-perpendicular axis would
//! squash the triangle to a sliver and lose the orientation the
//! triangulation depends on.
//!
//! # Honest limits
//!
//! This handles the case the intersection curve actually produces for the
//! operands `ScalarBoolean` supports: a face crossed by a chain of segments
//! that enters and leaves through its boundary. A curve forming a closed
//! loop strictly INSIDE one face is refused -- it needs a hole-aware
//! triangulation, and inventing a bridge edge to fake it would produce a
//! mesh whose topology no longer matches the geometry.

use std::collections::BTreeMap;

use axiolid_contracts::{GeomError, GeomResult, Sign};
use axiolid_core::{Point2, Point3};

use crate::intersection::{IntersectionSegment, NodeKey};
use crate::orient2d;

/// A face's corners plus the curve nodes lying on it, ready to triangulate.
#[derive(Debug, Clone)]
pub struct FacePatch {
    /// Every point, in the order the triangulation indexes them.
    ///
    /// Corners come first so the original winding stays recoverable.
    pub points: Vec<Point3>,
    /// Which curve node each point came from, where it came from one.
    ///
    /// `None` marks an original face corner. Retaining this lets a caller
    /// weld patches from adjacent faces by NODE IDENTITY rather than by
    /// comparing coordinates -- the same discipline the curve itself uses.
    pub sources: Vec<Option<NodeKey>>,
    /// Triangles as indices into `points`, wound like the source face.
    pub triangles: Vec<[u32; 3]>,
}

/// Which axis to drop when flattening, chosen from the largest normal term.
fn dominant_axis(normal: Point3) -> usize {
    let absolute = normal.abs();
    if absolute.x >= absolute.y && absolute.x >= absolute.z {
        0
    } else if absolute.y >= absolute.z {
        1
    } else {
        2
    }
}

/// Drop `axis`, keeping the other two coordinates in a fixed order.
fn project(point: Point3, axis: usize) -> Point2 {
    match axis {
        0 => Point2::new(point.y, point.z),
        1 => Point2::new(point.x, point.z),
        _ => Point2::new(point.x, point.y),
    }
}

/// Retriangulate one face against the curve segments lying on it.
///
/// `corners` are the face's three vertices, wound as the source mesh winds
/// them. `segments` are the curve segments on this face, and `positions`
/// resolves their nodes to coordinates.
///
/// The result reproduces the face exactly when `segments` is empty, so a
/// caller can run every face through this without special-casing.
pub fn retriangulate_face(
    corners: [Point3; 3],
    segments: &[IntersectionSegment],
    positions: &BTreeMap<NodeKey, Point3>,
) -> GeomResult<FacePatch> {
    let normal = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
    if normal.length_squared() == 0.0 {
        return Err(GeomError::Degenerate(
            "cannot retriangulate a degenerate face".into(),
        ));
    }

    // The untouched face is already its own triangulation.
    if segments.is_empty() {
        return Ok(FacePatch {
            points: corners.to_vec(),
            sources: vec![None; 3],
            triangles: vec![[0, 1, 2]],
        });
    }

    let mut points: Vec<Point3> = corners.to_vec();
    let mut sources: Vec<Option<NodeKey>> = vec![None; 3];
    let mut index_of: BTreeMap<NodeKey, u32> = BTreeMap::new();

    // Curve nodes join the corner list. A node that coincides with a corner
    // reuses that corner's index instead of adding a duplicate point, which
    // would leave the triangulation with a zero-length edge.
    for segment in segments {
        for node in [segment.start, segment.end] {
            if index_of.contains_key(&node) {
                continue;
            }
            let point = *positions.get(&node).ok_or_else(|| {
                GeomError::Degenerate("curve node has no recorded position".into())
            })?;
            if let Some(corner) = corners.iter().position(|&c| c == point) {
                index_of.insert(node, corner as u32);
                sources[corner] = Some(node);
                continue;
            }
            index_of.insert(node, points.len() as u32);
            points.push(point);
            sources.push(Some(node));
        }
    }

    let axis = dominant_axis(normal);
    let flat: Vec<Point2> = points.iter().map(|&p| project(p, axis)).collect();

    // Constraint edges, as index pairs into `points`.
    let mut constraints: Vec<(u32, u32)> = Vec::new();
    for segment in segments {
        let start = index_of[&segment.start];
        let end = index_of[&segment.end];
        if start != end {
            constraints.push((start.min(end), start.max(end)));
        }
    }
    constraints.sort_unstable();
    constraints.dedup();

    let triangles = triangulate_with_constraints(&flat, &constraints)?;

    // The 2D work happens in projected space, whose handedness depends on
    // which axis was dropped and which way the face pointed. Rewinding
    // against the ORIGINAL normal restores the source orientation, so the
    // patch can be substituted for the face without flipping it.
    let triangles = triangles
        .into_iter()
        .map(|tri| {
            let [a, b, c] = tri.map(|i| points[i as usize]);
            if (b - a).cross(c - a).dot(normal) < 0.0 {
                [tri[0], tri[2], tri[1]]
            } else {
                tri
            }
        })
        .collect();

    Ok(FacePatch {
        points,
        sources,
        triangles,
    })
}

/// Triangulate a point set so every constraint edge appears in the output.
///
/// # Approach
///
/// A brute-force maximal triangulation: consider every candidate triangle,
/// keep those that are non-degenerate, contain no other point, and cross no
/// constraint. `O(n^4)`, which is the right trade for a reference -- it is
/// short enough to audit line by line, and the input is one triangle's worth
/// of points, not a mesh.
///
/// A production provider would use a proper CDT. This exists to be
/// obviously correct, so a fast implementation has something to be checked
/// against.
fn triangulate_with_constraints(
    points: &[Point2],
    constraints: &[(u32, u32)],
) -> GeomResult<Vec<[u32; 3]>> {
    let count = points.len();
    let mut triangles: Vec<[u32; 3]> = Vec::new();

    for a in 0..count {
        for b in (a + 1)..count {
            for c in (b + 1)..count {
                let tri = [a as u32, b as u32, c as u32];
                let [pa, pb, pc] = [points[a], points[b], points[c]];

                // A collinear triple has no area and would contribute a
                // sliver that later orientation tests cannot classify.
                if sign(orient2d(pa, pb, pc)) == Sign::Zero {
                    continue;
                }
                // A triangle covering another point is not part of any
                // valid triangulation of the full point set.
                //
                // UNPROVEN: no fixture reaches this branch, and mutating it
                // away leaves every test passing. It is kept because the
                // smallest-area-first order makes a covering triangle lose
                // anyway, not because a test demonstrates the need. Delete it
                // only alongside a case that shows it is genuinely dead.
                if (0..count).any(|other| {
                    other != a
                        && other != b
                        && other != c
                        && point_inside(points[other], [pa, pb, pc])
                }) {
                    continue;
                }
                // A triangle edge cutting across a constraint would erase
                // the cut the whole operation exists to make.
                if crosses_a_constraint(tri, points, constraints) {
                    continue;
                }
                triangles.push(tri);
            }
        }
    }

    // Candidates may still overlap each other, so a maximal non-overlapping
    // subset has to be chosen. Order matters: a greedy pass that took the
    // whole face first would block every finer triangle, since the face
    // overlaps all of them, and the cut would vanish. Smallest-area-first
    // makes the fine pieces win and the coarse cover lose.
    //
    // Ties break on the index triple so the result is deterministic for a
    // given input rather than dependent on sort stability.
    triangles.sort_by(|left, right| {
        let area = |t: &[u32; 3]| {
            let [a, b, c] = t.map(|i| points[i as usize]);
            ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).abs()
        };
        area(left)
            .partial_cmp(&area(right))
            .expect("finite coordinates give comparable areas")
            .then_with(|| left.cmp(right))
    });

    let mut kept: Vec<[u32; 3]> = Vec::new();
    for tri in triangles {
        if kept.iter().any(|existing| overlaps(*existing, tri, points)) {
            continue;
        }
        kept.push(tri);
    }

    if kept.is_empty() {
        return Err(GeomError::Degenerate(
            "no valid triangle survives the constraints".into(),
        ));
    }

    // Every constraint must survive as an edge of some kept triangle.
    // Reporting this rather than returning a plausible-looking mesh is the
    // difference between a refusal and a silently wrong cut.
    //
    // UNPROVEN: no fixture triggers this refusal. Mutating it away leaves the
    // suite green, so it is a belt-and-braces check, not a tested guarantee.
    // A case that reaches it would be a valuable addition.
    for &(start, end) in constraints {
        let present = kept.iter().any(|tri| {
            [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
                .iter()
                .any(|&(u, v)| (u.min(v), u.max(v)) == (start, end))
        });
        if !present {
            return Err(GeomError::Unsupported {
                backend: axiolid_contracts::BackendId::new("scalar-retriangulate"),
                operation: axiolid_contracts::Operation::MeshBoolean,
            });
        }
    }

    Ok(kept)
}

/// Strictly inside the triangle: on an edge does not count.
///
/// Boundary points are excluded deliberately. A point ON an edge is shared
/// with the neighbouring triangle and does not invalidate either.
fn point_inside(point: Point2, [a, b, c]: [Point2; 3]) -> bool {
    let signs = [
        sign(orient2d(a, b, point)),
        sign(orient2d(b, c, point)),
        sign(orient2d(c, a, point)),
    ];
    signs.iter().all(|&s| s == Sign::Positive) || signs.iter().all(|&s| s == Sign::Negative)
}

/// Whether any edge of `tri` properly crosses any constraint.
///
/// Sharing an endpoint is not a crossing: constraints meet each other and
/// the face boundary at nodes, which is exactly what they are meant to do.
fn crosses_a_constraint(tri: [u32; 3], points: &[Point2], constraints: &[(u32, u32)]) -> bool {
    let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
    for &(u, v) in &edges {
        for &(s, e) in constraints {
            if u == s || u == e || v == s || v == e {
                continue;
            }
            if segments_properly_cross(
                [points[u as usize], points[v as usize]],
                [points[s as usize], points[e as usize]],
            ) {
                return true;
            }
        }
    }
    false
}

/// Two segments crossing at an interior point of both.
fn segments_properly_cross([a, b]: [Point2; 2], [c, d]: [Point2; 2]) -> bool {
    let d1 = sign(orient2d(a, b, c));
    let d2 = sign(orient2d(a, b, d));
    let d3 = sign(orient2d(c, d, a));
    let d4 = sign(orient2d(c, d, b));
    d1 != Sign::Zero
        && d2 != Sign::Zero
        && d3 != Sign::Zero
        && d4 != Sign::Zero
        && d1 != d2
        && d3 != d4
}

/// Whether two triangles share interior area.
///
/// Tested by centroid containment both ways plus proper edge crossings.
/// Triangles that merely share a vertex or an edge do not overlap, which is
/// the normal case in any triangulation.
fn overlaps(first: [u32; 3], second: [u32; 3], points: &[Point2]) -> bool {
    let fa = first.map(|i| points[i as usize]);
    let sa = second.map(|i| points[i as usize]);

    if point_inside(centroid(fa), sa) || point_inside(centroid(sa), fa) {
        return true;
    }
    for i in 0..3 {
        for j in 0..3 {
            let first_edge = [fa[i], fa[(i + 1) % 3]];
            let second_edge = [sa[j], sa[(j + 1) % 3]];
            if segments_properly_cross(first_edge, second_edge) {
                return true;
            }
        }
    }
    false
}

/// The average of three corners, which lies strictly inside the triangle.
fn centroid([a, b, c]: [Point2; 3]) -> Point2 {
    Point2::new((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0)
}

fn sign(value: axiolid_contracts::Certified) -> Sign {
    value.sign().expect("certified predicates are total")
}
