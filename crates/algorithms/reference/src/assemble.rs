//! Assembling an exact boolean result from retriangulated operands.
//!
//! # The three steps
//!
//! 1. [`intersection_segments`] finds WHERE the
//!    two surfaces cross.
//! 2. [`retriangulate_face`] rebuilds each cut face
//!    so the curve exists as mesh edges.
//! 3. This module decides which of the resulting pieces to keep.
//!
//! # Why step 2 makes step 3 easy
//!
//! After retriangulation no triangle straddles the other solid's surface:
//! each one lies wholly inside or wholly outside. So a single containment
//! test per triangle settles it, and the test can be taken at the centroid --
//! a point guaranteed to be in the triangle's interior, away from the
//! boundary where classification is ambiguous.
//!
//! Without step 2 this would be false: a triangle crossing the surface has no
//! single answer, and sampling it anywhere would be a guess.
//!
//! # Winding
//!
//! `Difference` keeps the subject's outside and the tool's inside, but the
//! tool's kept faces must be REVERSED: they become the cavity wall, and a
//! cavity's outward normal points into the removed volume. Getting this
//! wrong produces a mesh that looks right and has the wrong sign everywhere.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{BooleanOperator, Point3};
use axiolid_mesh::TriMesh;

use crate::boolean::contains_point_exact;
use crate::intersection::{intersection_segments, IntersectionSegment, NodeKey};
use crate::retriangulate::retriangulate_face;
use std::collections::BTreeMap;

/// Which side of the other solid a piece lies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    /// Strictly inside the other operand.
    Inside,
    /// Strictly outside it.
    Outside,
}

/// Rebuild one operand's faces against the curve, and classify each piece.
///
/// Returns the retriangulated triangles paired with the side they fall on,
/// so the caller can keep whichever the operation asks for.
fn split_and_classify(
    mesh: &TriMesh,
    other: &TriMesh,
    face_segments: &BTreeMap<u32, Vec<IntersectionSegment>>,
    positions: &BTreeMap<NodeKey, Point3>,
) -> GeomResult<Vec<([Point3; 3], Side)>> {
    let mut pieces = Vec::new();

    for face_index in 0..mesh.triangle_count() {
        let base = face_index * 3;
        let corners: [Point3; 3] = [
            mesh.positions[mesh.indices[base] as usize],
            mesh.positions[mesh.indices[base + 1] as usize],
            mesh.positions[mesh.indices[base + 2] as usize],
        ];

        let key = u32::try_from(face_index).unwrap_or(u32::MAX);
        let empty: Vec<IntersectionSegment> = Vec::new();
        let segments = face_segments.get(&key).unwrap_or(&empty);
        let patch = retriangulate_face(corners, segments, positions)?;

        for triangle in &patch.triangles {
            let [a, b, c] = triangle.map(|i| patch.points[i as usize]);
            // The centroid is interior to its OWN piece, but that says
            // nothing about where it falls relative to the other operand: it
            // can land exactly on the other's face, edge, or corner. Measured,
            // not assumed -- a corner-overlap subtraction puts a centroid on
            // exactly (2,4,2), the tool's corner edge.
            let centroid = Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            );
            // A centroid on the other operand's surface is unclassifiable, so
            // retry from points nudged toward each corner. These stay strictly
            // inside the piece -- so they classify the SAME piece -- while
            // moving off whatever feature the centroid landed on. Refusing
            // beats guessing: treating an on-surface centroid as "outside"
            // keeps a piece whose neighbours were dropped, leaving a hole that
            // still reports a plausible volume.
            let probes = [
                centroid,
                lerp(centroid, a, 0.25),
                lerp(centroid, b, 0.25),
                lerp(centroid, c, 0.25),
            ];
            let inside = probes
                .into_iter()
                .find_map(|probe| contains_point_exact(other, probe));
            let side = match inside {
                Some(true) => Side::Inside,
                Some(false) => Side::Outside,
                None => {
                    return Err(GeomError::Unsupported {
                        backend: BackendId::new("scalar-assemble"),
                        operation: Operation::MeshBoolean,
                    })
                }
            };
            pieces.push(([a, b, c], side));
        }
    }
    Ok(pieces)
}

/// Compute a boolean of two interpenetrating solids, exactly.
///
/// This is the case [`ScalarBoolean`](crate::ScalarBoolean) refuses: surfaces
/// that properly cross. Every decision is an exact predicate -- the curve from
/// `orient3d` signs, the retriangulation from `orient2d` signs, the
/// classification from ray parity -- so the result is not a tolerance
/// approximation of the answer, it is the answer.
///
/// # Errors
///
/// Propagates the curve's refusals: coplanar face overlap, degenerate faces,
/// and a curve that cannot be stitched. Refusing is deliberate; a boolean
/// that guesses in those cases is worse than one that declines.
pub fn exact_boolean(
    subject: &TriMesh,
    tool: &TriMesh,
    operation: BooleanOperator,
) -> GeomResult<TriMesh> {
    let curve = intersection_segments(subject, tool)?;

    let subject_pieces = split_and_classify(
        subject,
        tool,
        &curve.subject_face_segments,
        &curve.positions,
    )?;
    let tool_pieces =
        split_and_classify(tool, subject, &curve.tool_face_segments, &curve.positions)?;

    // Which side of each operand the operation keeps, and whether the tool's
    // kept faces have to be flipped.
    //
    // Difference keeps the subject's outside and the tool's inside, and the
    // tool's faces become the cavity wall: their outward normal must point
    // INTO the removed volume, so they are reversed. Union and Intersection
    // keep consistently-oriented faces from both, so they are not.
    let (keep_subject, keep_tool, flip_tool) = match operation {
        BooleanOperator::Union => (Side::Outside, Side::Outside, false),
        BooleanOperator::Intersection => (Side::Inside, Side::Inside, false),
        BooleanOperator::Difference => (Side::Outside, Side::Inside, true),
        _ => {
            return Err(axiolid_contracts::GeomError::Unsupported {
                backend: axiolid_contracts::BackendId::new("scalar-exact-boolean"),
                operation: Operation::MeshBoolean,
            })
        }
    };

    let mut positions: Vec<Point3> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    // Vertices are welded on exact coordinate bits, the same discipline the
    // curve uses for node identity. Two pieces that meet along the cut were
    // computed from the same arithmetic, so their shared corners are
    // bit-identical and join into one vertex -- leaving the result closed
    // rather than a shell of unconnected triangles.
    let mut welded: BTreeMap<[u64; 3], u32> = BTreeMap::new();
    let mut push = |point: Point3, positions: &mut Vec<Point3>| -> u32 {
        let bits = [point.x + 0.0, point.y + 0.0, point.z + 0.0].map(f64::to_bits);
        *welded.entry(bits).or_insert_with(|| {
            positions.push(point);
            (positions.len() - 1) as u32
        })
    };

    for (triangle, side) in &subject_pieces {
        if *side != keep_subject {
            continue;
        }
        for corner in triangle {
            let index = push(*corner, &mut positions);
            indices.push(index);
        }
    }
    for (triangle, side) in &tool_pieces {
        if *side != keep_tool {
            continue;
        }
        // Reversing swaps two corners, which flips the winding and so the
        // outward normal.
        let ordered = if flip_tool {
            [triangle[0], triangle[2], triangle[1]]
        } else {
            *triangle
        };
        for corner in &ordered {
            let index = push(*corner, &mut positions);
            indices.push(index);
        }
    }

    Ok(TriMesh::new(positions, indices))
}

/// A point a fraction of the way from `from` toward `to`.
///
/// Used to move a probe off a degenerate feature while keeping it strictly
/// inside the same piece, so it still classifies that piece.
fn lerp(from: Point3, to: Point3, t: f64) -> Point3 {
    Point3::new(
        from.x + (to.x - from.x) * t,
        from.y + (to.y - from.y) * t,
        from.z + (to.z - from.z) * t,
    )
}
