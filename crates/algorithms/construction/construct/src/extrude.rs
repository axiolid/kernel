//! Linear extrusion of a triangulated profile into a closed solid.
//!
//! The result must be watertight and outward-oriented, because that is exactly
//! what `axiolid-mesh-boolean-boolmesh` demands of its inputs. Getting the winding wrong here
//! produces a mesh that looks valid and computes wrong booleans -- the failure
//! mode that has already cost two debugging sessions.
//!
//! # A direction in the profile plane (#281)
//!
//! An offset with no component along the profile normal bounds no volume:
//! its "solid" is a sliver, or flat at exactly `z = 0`. The exact path
//! refuses an offset within tolerance of the plane, `|o.z| <= tolerance`,
//! as `"extrusion direction in the profile plane"`, and so does
//! [`extrude_profile`], by the same typed refusal, so the two paths agree
//! on what they build. [`extrude`] carries no tolerance, so it refuses
//! only an offset exactly in the plane, by the same name.

use axiolid_contracts::{GeomError, GeomResult, Operation, Sign};
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_reference::arithmetic::{
    expansion_sign, expansion_sum, grow_expansion, scale_expansion,
};
use axiolid_reference::expansion::two_product;

pub use crate::extrude_exact::extrude_profile_exact;

/// Extrude a triangulated 2D profile along `direction` by `depth`.
///
/// The profile lies in the local z = 0 plane. Caps use the triangulation;
/// sides are quads split into two triangles per boundary edge.
///
/// `boundary` lists the closed loops of the profile as index ranges into
/// `points`: each loop is a contiguous run, matching `profile::Rings` layout.
///
/// An offset exactly in the profile plane (`o.z == 0`) bounds no volume and
/// is refused as `"extrusion direction in the profile plane"`; with no
/// tolerance here, a sliver just off the plane is built.
/// [`extrude_profile`] refuses within its tolerance (see the module docs).
pub fn extrude(
    points: &[Point2],
    triangles: &[[u32; 3]],
    loops: &[core::ops::Range<usize>],
    direction: Vec3,
    depth: Scalar,
) -> GeomResult<TriMesh> {
    let offset = mesh_offset(direction, depth)?;
    if offset.z == 0.0 {
        return Err(in_profile_plane());
    }

    let n = points.len();
    let mut positions = Vec::with_capacity(n * 2);
    // Base ring first, then the offset ring: vertex i has its twin at i + n.
    positions.extend(points.iter().map(|p| Point3::new(p.x, p.y, 0.0)));
    positions.extend(points.iter().map(|p| Point3::new(p.x, p.y, 0.0) + offset));

    let mut indices: Vec<u32> = Vec::with_capacity(triangles.len() * 6 + n * 6);
    let top = n as u32;

    // Caps. The profile triangulation is counter-clockwise seen from +z, which
    // is outward for the TOP cap and inward for the bottom, so the bottom is
    // emitted reversed.
    for t in triangles {
        indices.extend_from_slice(&[t[0] + top, t[1] + top, t[2] + top]);
        indices.extend_from_slice(&[t[0], t[2], t[1]]);
    }

    // Sides. Each boundary edge (a -> b) becomes the quad a, b, b', a'.
    for range in loops {
        let len = range.len();
        if len < 3 {
            return Err(GeomError::InvalidInput(format!(
                "extrusion loop needs at least 3 vertices, got {len}"
            )));
        }
        for k in 0..len {
            let a = (range.start + k) as u32;
            let b = (range.start + (k + 1) % len) as u32;
            indices.extend_from_slice(&[a, b, b + top]);
            indices.extend_from_slice(&[a, b + top, a + top]);
        }
    }

    // Every winding above assumes the offset leaves the profile plane towards
    // +z. When it points below the plane the solid is the mirror image of that
    // case, so every triangle is inside-out (signed volume -area*depth) and a
    // boolean would refuse or silently invert it. Flipping each triangle once
    // restores outward orientation for caps and walls, outer and hole loops
    // alike. An offset IN the plane bounds no volume either way and was
    // refused above.
    if offset.z < 0.0 {
        for triangle in indices.chunks_exact_mut(3) {
            triangle.swap(1, 2);
        }
    }

    Ok(TriMesh::new(positions, indices))
}

/// Triangulate rings and extrude them in one step.
///
/// The loop layout must match `triangulate`'s vertex order exactly, so the
/// two are derived from the same `Rings` value here rather than by a caller
/// reconstructing the ranges.
///
/// A direction within `tolerance` of the profile plane (`|o.z| <=
/// tolerance` for the offset `o`) is refused as `"extrusion direction in
/// the profile plane"`, exactly as [`extrude_profile_exact`] refuses it.
pub fn extrude_profile(
    rings: &crate::profile::Rings,
    direction: Vec3,
    depth: Scalar,
    tolerance: axiolid_core::Tolerance,
) -> GeomResult<TriMesh> {
    // Refused before triangulating, by the exact path's rule (#281).
    if mesh_offset(direction, depth)?.z.abs() <= tolerance.linear() {
        return Err(in_profile_plane());
    }
    let (points, triangles) = crate::profile::triangulate(rings)?;
    let mut loops = Vec::with_capacity(1 + rings.holes.len());
    let mut start = 0usize;
    loops.push(start..rings.outer.len());
    start += rings.outer.len();
    for hole in &rings.holes {
        loops.push(start..start + hole.len());
        start += hole.len();
    }
    extrude(&points, &triangles, &loops, direction, depth)
}

/// The offset `direction / |direction| * depth`, validated.
fn mesh_offset(direction: Vec3, depth: Scalar) -> GeomResult<Vec3> {
    if !depth.is_finite() || depth <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "extrusion depth must be positive and finite, got {depth}"
        )));
    }
    if !direction.is_finite() || direction.length() <= 0.0 {
        return Err(GeomError::InvalidInput(
            "extrusion direction must be a finite non-zero vector".to_owned(),
        ));
    }
    let offset = direction.normalize() * depth;
    if !offset.is_finite() {
        return Err(GeomError::Degenerate(
            "extrusion direction could not be normalised".to_owned(),
        ));
    }
    Ok(offset)
}

/// The refusal both extrusion paths give an offset in the profile plane.
fn in_profile_plane() -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::Sweep,
        input: "extrusion direction in the profile plane",
    }
}

/// Whether a closed mesh is outward-oriented.
///
/// A face-counting majority does not work: a hollow section's inner wall
/// legitimately faces the opposite way from its outer wall, and for a thin
/// tube the two counts are comparable. Orientation is a property of the
/// enclosed volume, not of how faces point relative to a centre.
///
/// So the volume is summed exactly. Each tetrahedron about the reference point
/// contributes `a . (b x c)`, and those contributions are accumulated in
/// expansion arithmetic rather than f64, so the final sign is certified even
/// when the terms cancel catastrophically -- which is exactly what happens for
/// a thin plate or a large solid far from the origin.
///
/// Returns `None` when the mesh encloses exactly zero volume, which is not an
/// orientation and must not be reported as one.
#[must_use]
pub fn outward_orientation(mesh: &TriMesh) -> Option<bool> {
    if mesh.indices.len() < 12 {
        // Fewer than four triangles cannot bound a volume.
        return None;
    }
    let mut total: Vec<f64> = vec![0.0];
    for corner in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[corner[0] as usize];
        let b = mesh.positions[corner[1] as usize];
        let c = mesh.positions[corner[2] as usize];
        total = expansion_sum(&total, &triple_product(a, b, c));
    }
    match expansion_sign(&total) {
        Sign::Positive => Some(true),
        Sign::Negative => Some(false),
        Sign::Zero => None,
        // `Sign` is non-exhaustive; an unrecognised variant is not a verdict.
        _ => None,
    }
}

/// Exact `a . (b x c)` as an expansion: six times a tetrahedron's volume.
#[must_use]
fn triple_product(a: Point3, b: Point3, c: Point3) -> Vec<f64> {
    let term = |p: f64, q: f64, r: f64, s: f64, k: f64| {
        // k * (p*q - r*s), exactly.
        scale_expansion(&exact_difference_of_products(p, q, r, s), k)
    };
    let x = term(b.y, c.z, c.y, b.z, a.x);
    let y = term(b.z, c.x, c.z, b.x, a.y);
    let z = term(b.x, c.y, c.x, b.y, a.z);
    expansion_sum(&expansion_sum(&x, &y), &z)
}

/// Exact `p*q - r*s` as a four-term expansion.
#[must_use]
fn exact_difference_of_products(p: f64, q: f64, r: f64, s: f64) -> Vec<f64> {
    let (pq, pq_err) = two_product(p, q);
    let (rs, rs_err) = two_product(r, s);
    let e = grow_expansion(&[pq_err], -rs_err);
    let e = grow_expansion(&e, pq);
    grow_expansion(&e, -rs)
}
