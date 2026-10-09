//! Constant-distance offset and shelling of planar-faced solids (#78).
//!
//! # Which offset this is
//!
//! Offsetting by a sphere rounds every convex edge, which introduces curved
//! faces and leaves the planar family this crate refuses to leave. What is
//! implemented is the MITER offset: each face plane is pushed along its
//! normal by the distance, and the faces are extended until they meet again.
//! A cube offset outward by `d` becomes a cube of edge `a + 2d`, which is
//! the identity #78 verifies against.
//!
//! # Why vertices, not faces
//!
//! Pushing each face and re-deriving the solid works for convex operands and
//! silently fails for concave ones: at a reflex edge the pushed neighbours
//! overlap instead of meeting. Offsetting the VERTICES instead -- each moved
//! to the common point of its own incident pushed planes -- gives the miter
//! point in both cases, and preserves topology exactly, so the result has the
//! same face structure as its input and needs no re-assembly.

use crate::boolean_exact::unsupported;
use crate::polyhedron::Polyhedron;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Vec3};
use std::collections::BTreeMap;

/// Which way an offset moves the boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffsetDirection {
    /// Grow the solid.
    Outward,
    /// Shrink the solid. Collapses if the distance exceeds the local
    /// half-thickness, which is refused rather than emitted.
    Inward,
}

/// A face's supporting plane, as a unit normal and an offset along it.
#[derive(Debug, Clone, Copy)]
struct Plane {
    normal: Vec3,
    distance: f64,
}

/// Offset a planar-faced solid by a constant distance.
///
/// Every vertex moves to the intersection of its own incident face planes,
/// each pushed by `distance`. Topology is preserved, so the result has the
/// same faces as the input with new coordinates.
///
/// # Errors
///
/// Refuses a non-positive or non-finite distance, a vertex whose incident
/// planes do not meet in a single point, and any inward offset that collapses
/// the solid.
pub fn offset_solid(
    solid: &Polyhedron,
    distance: f64,
    direction: OffsetDirection,
) -> GeomResult<Polyhedron> {
    if !distance.is_finite() || distance <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "offset distance {distance} must be positive and finite"
        )));
    }
    let signed = match direction {
        OffsetDirection::Outward => distance,
        OffsetDirection::Inward => -distance,
    };

    let planes = face_planes(solid)?;
    let incidence = vertex_incidence(solid);

    let mut moved: BTreeMap<VertexKey, Point3> = BTreeMap::new();
    for (key, faces) in &incidence {
        let point = miter_point(&planes, faces, signed, key)?;
        moved.insert(*key, point);
    }

    let faces = solid
        .faces()
        .iter()
        .map(|face| {
            face.iter()
                .map(|&v| moved[&VertexKey::of(v)])
                .collect::<Vec<_>>()
        })
        .collect();

    let result = Polyhedron::new(faces)?;
    reject_collapse(solid, &result, &planes, signed)?;
    Ok(result)
}

/// Exact-coordinate key so vertices shared between faces are recognised.
///
/// Faces arrive as independent rings, so the same corner appears once per
/// incident face. Bit keying is correct because those copies come from the
/// same input literal or the same boolean split, never from separate
/// arithmetic, so no welding tolerance is invented here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct VertexKey([u64; 3]);

impl VertexKey {
    fn of(p: Point3) -> Self {
        // `+ 0.0` keys -0.0 as 0.0, the same corner (#269).
        Self([
            (p.x + 0.0).to_bits(),
            (p.y + 0.0).to_bits(),
            (p.z + 0.0).to_bits(),
        ])
    }
}

/// Unit normal and plane offset for each face, in face order.
fn face_planes(solid: &Polyhedron) -> GeomResult<Vec<Plane>> {
    solid
        .faces()
        .iter()
        .map(|face| {
            let raw = (face[1] - face[0]).cross(face[2] - face[0]);
            let length = raw.length();
            if length == 0.0 {
                return Err(unsupported("degenerate face has no offset direction"));
            }
            let normal = raw / length;
            Ok(Plane {
                normal,
                distance: normal.dot(face[0] - Point3::new(0.0, 0.0, 0.0)),
            })
        })
        .collect()
}

/// Which faces touch each vertex.
fn vertex_incidence(solid: &Polyhedron) -> BTreeMap<VertexKey, Vec<usize>> {
    let mut map: BTreeMap<VertexKey, Vec<usize>> = BTreeMap::new();
    for (index, face) in solid.faces().iter().enumerate() {
        for &v in face {
            let entry = map.entry(VertexKey::of(v)).or_default();
            if !entry.contains(&index) {
                entry.push(index);
            }
        }
    }
    map
}

/// Where a vertex lands when its incident planes are all pushed by `signed`.
///
/// Three independent planes meet in one point, which is the miter vertex.
/// This is what makes concave edges work: at a reflex edge the incident
/// planes still meet, and their common point is the correct inside corner,
/// whereas pushing faces and re-deriving the solid would leave them
/// overlapping.
///
/// A vertex with more than three incident faces is over-determined: the
/// planes only stay concurrent for special geometry, and if they do not, the
/// offset genuinely has no single answer there. Three are solved and the rest
/// verified, so a disagreement is refused rather than silently resolved by
/// picking the first three.
fn miter_point(
    planes: &[Plane],
    faces: &[usize],
    signed: f64,
    key: &VertexKey,
) -> GeomResult<Point3> {
    if faces.len() < 3 {
        return Err(unsupported(
            "a vertex needs 3 incident faces to determine an offset corner",
        ));
    }
    let (a, b, c) = (planes[faces[0]], planes[faces[1]], planes[faces[2]]);
    let point = intersect_three(a, b, c, signed)
        .ok_or_else(|| unsupported("incident face planes are parallel; offset corner undefined"))?;

    for &extra in &faces[3..] {
        let plane = planes[extra];
        let residual =
            plane.normal.dot(point - Point3::new(0.0, 0.0, 0.0)) - (plane.distance + signed);
        // Scale-relative: an absolute epsilon would reject large models and
        // accept tiny ones.
        let scale = point.x.abs().max(point.y.abs()).max(point.z.abs()).max(1.0);
        if residual.abs() > 1e-9 * scale {
            let _ = key;
            return Err(unsupported(
                "over-determined vertex: incident planes do not meet in one offset point",
            ));
        }
    }
    Ok(point)
}

/// Common point of three pushed planes, by Cramer's rule.
///
/// `None` when the planes are parallel or share a line, which has no unique
/// offset corner.
fn intersect_three(a: Plane, b: Plane, c: Plane, signed: f64) -> Option<Point3> {
    let bc = b.normal.cross(c.normal);
    let determinant = a.normal.dot(bc);
    // Scale-relative singularity test: the normals are unit vectors, so the
    // determinant is the parallelepiped volume they span and is dimensionless.
    if determinant.abs() < 1e-12 {
        return None;
    }
    let (da, db, dc) = (
        a.distance + signed,
        b.distance + signed,
        c.distance + signed,
    );
    let numerator = bc * da + c.normal.cross(a.normal) * db + a.normal.cross(b.normal) * dc;
    Some(Point3::new(0.0, 0.0, 0.0) + numerator / determinant)
}

/// Refuse an offset that turned the solid inside out.
///
/// An inward offset deeper than the local half-thickness makes opposing walls
/// cross. The signature is local and unmistakable: a face whose normal
/// reverses has been pushed past its opposite neighbour, so the boundary now
/// passes through itself.
///
/// Checking normals rather than volume is deliberate. A volume test needs the
/// result triangulated, and fanning a NON-CONVEX face emits triangles outside
/// the footprint, so the measurement would be wrong exactly on the operands
/// this offset exists to support. Normal preservation is a per-face invariant
/// that holds regardless of face convexity.
fn reject_collapse(
    input: &Polyhedron,
    result: &Polyhedron,
    planes: &[Plane],
    signed: f64,
) -> GeomResult<()> {
    for (before, after) in input.faces().iter().zip(result.faces()) {
        let n0 = (before[1] - before[0]).cross(before[2] - before[0]);
        let n1 = (after[1] - after[0]).cross(after[2] - after[0]);
        if n1.length() == 0.0 {
            return Err(unsupported(
                "offset collapsed a face to zero area; distance exceeds the local half-thickness",
            ));
        }
        if n0.dot(n1) <= 0.0 {
            return Err(unsupported(
                "offset reversed a face normal; the boundary passes through itself",
            ));
        }
    }

    // Normals alone are not enough, and neither is travel distance. Offset a
    // 2-cube inward by 1.5: the z=2 face lands at z=0.5 and the z=0 face at
    // z=1.5. Each kept its normal and each travelled exactly the requested
    // 1.5 -- yet they have SWAPPED SIDES and the solid is inside out.
    //
    // The collapse mode is opposing walls crossing, which is precisely the
    // "local half-thickness" the operation is bounded by. Two anti-parallel
    // faces bound a slab whose width is the sum of their plane offsets; an
    // inward offset narrows that slab by twice the distance, and a width that
    // reaches zero or goes negative means the walls have met or passed
    // through each other.
    //
    // Only anti-parallel pairs are tested. A containment or half-space test
    // would also catch this, but rejects non-convex operands: the L-prism's
    // reflex corner correctly moves OUTSIDE the offset planes of the faces
    // forming its notch, and that is the case this offset exists to support.
    for (i, a) in planes.iter().enumerate() {
        for b in &planes[i + 1..] {
            // Anti-parallel: normals opposed to within rounding of exactly
            // -1, which is the configuration that forms a bounded slab.
            if a.normal.dot(b.normal) > -1.0 + 1e-12 {
                continue;
            }
            let width_before = a.distance + b.distance;
            let width_after = width_before + 2.0 * signed;
            if width_before > 0.0 && width_after <= 1e-12 * width_before.max(1.0) {
                return Err(unsupported(
                    "offset closed the gap between opposing walls; \
                     distance exceeds the local half-thickness",
                ));
            }
        }
    }
    Ok(())
}

/// Hollow a solid to a stated wall thickness.
///
/// The shell is the difference between the solid and its inward offset, which
/// is why this needs #77's general boolean: the cavity of a non-convex solid
/// is not expressible any other way.
///
/// # Errors
///
/// Refuses a thickness that collapses the cavity, propagating the offset's
/// own refusal rather than returning a degenerate shell.
pub fn shell_solid(solid: &Polyhedron, thickness: f64) -> GeomResult<Polyhedron> {
    let cavity = offset_solid(solid, thickness, OffsetDirection::Inward)?;
    crate::polyhedron::boolean_polyhedra_exact(
        solid,
        &cavity,
        crate::polyhedron::BooleanOp::Difference,
    )
}
