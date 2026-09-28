#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Narrow-phase ray/triangle-mesh nearest-hit intersection.
//!
//! # Why this is its own package
//!
//! `axiolid-spatial` already owns the broad phase: [`SpatialIndex::visit_ray`]
//! walks a BVH and yields candidate keys. Without a narrow phase a caller gets
//! boxes and has to write ray/triangle themselves, which is how tolerance
//! policy fragments across consumers. This package closes that seam and owns
//! nothing else.
//!
//! It deliberately does not depend on `axiolid-spatial`: the narrow phase is
//! useful without an index, and an index is useful without this. Composition
//! happens at the call site by feeding candidate triangle indices into
//! [`nearest_hit_among`].
//!
//! # Boundary
//!
//! This package owns the intersection and the hit record. It does not own what
//! a ray *means*: sampling patterns, camera rigs, entity identity, or whether a
//! hit counts as an obstruction stay with the caller.
//!
//! # Fail closed, never silently miss
//!
//! A degenerate (zero-area) triangle has no well-defined ray intersection. This
//! package refuses with [`RayMeshError::DegenerateTriangle`] rather than
//! reporting a miss, because a silent miss is indistinguishable from real empty
//! space and quietly corrupts containment and visibility answers built on it.
//!
//! Ray direction is not required to be normalised, so the reported `t` is in
//! units of the supplied direction vector. That is stated rather than fixed up,
//! because normalising a caller's ray silently changes the meaning of every
//! distance they compare against.
//!
//! # What the tolerance decides
//!
//! `Tolerance::linear()` bounds only the parallel-ray rejection and the
//! barycentric slack that keeps edge and vertex hits. It never decides the
//! front/back branch: [`FaceSide`] comes from the certified `orient3d` sign.
//!
//! [`SpatialIndex::visit_ray`]: https://docs.rs/axiolid-spatial

use core::fmt;

use axiolid_core::{Point3, Ray3, Scalar, Tolerance};
use axiolid_guarantees::Sign;
use axiolid_mesh::TriangleMeshView;
use axiolid_predicates::orient3d;

/// Which side of a triangle the ray arrived from.
///
/// Determined by the certified orientation of the ray origin against the
/// triangle plane, not by the sign of a floating-point dot product, so a
/// grazing ray does not flip sides on rounding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceSide {
    /// The origin lies on the positive side of the triangle's winding normal.
    Front,
    /// The origin lies on the negative side of the triangle's winding normal.
    Back,
    /// The origin lies exactly in the triangle's plane.
    Coplanar,
}

/// One nearest-hit record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayHit3 {
    /// Parametric distance along the supplied (possibly unnormalised) direction.
    pub t: Scalar,
    /// Index of the hit triangle in the source mesh.
    pub triangle: usize,
    /// Barycentric coordinates `(u, v, w)` with `w = 1 - u - v`, ordered to
    /// match the triangle's stored corner order.
    pub barycentric: [Scalar; 3],
    /// Side the ray origin was on.
    pub side: FaceSide,
    /// Hit position reconstructed as `origin + direction * t`.
    pub point: Point3,
}

/// Fail-closed reasons a ray/mesh query cannot produce an answer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RayMeshError {
    /// A ray or mesh coordinate is NaN or infinite.
    NonFiniteInput,
    /// The ray direction is exactly zero, so no parametric distance exists.
    ZeroDirection,
    /// The tolerance policy is not usable for a parametric query.
    InvalidTolerance,
    /// A triangle references a position outside the mesh's position buffer.
    PositionIndexOutOfRange {
        /// Offending triangle.
        triangle: usize,
    },
    /// A triangle has zero area, so it has no defined ray intersection.
    ///
    /// Reported rather than skipped: a silent miss is indistinguishable from
    /// empty space.
    DegenerateTriangle {
        /// Offending triangle.
        triangle: usize,
    },
}

impl fmt::Display for RayMeshError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteInput => formatter.write_str("ray and mesh coordinates must be finite"),
            Self::ZeroDirection => formatter.write_str("ray direction must be non-zero"),
            Self::InvalidTolerance => {
                formatter.write_str("ray/mesh tolerance must be finite and non-negative")
            }
            Self::PositionIndexOutOfRange { triangle } => {
                write!(
                    formatter,
                    "triangle {triangle} references a missing position"
                )
            }
            Self::DegenerateTriangle { triangle } => {
                write!(formatter, "triangle {triangle} has zero area")
            }
        }
    }
}

impl std::error::Error for RayMeshError {}

/// Nearest hit over every triangle of `mesh`.
///
/// Prefer [`nearest_hit_among`] when a broad phase has already rejected most
/// triangles; this scans all of them.
pub fn nearest_hit(
    mesh: &impl TriangleMeshView,
    ray: &Ray3,
    tolerance: Tolerance,
) -> Result<Option<RayHit3>, RayMeshError> {
    nearest_hit_among(mesh, ray, tolerance, 0..mesh.triangle_count())
}

/// Nearest hit over caller-supplied candidate triangles.
///
/// This is the composition point with a broad phase: feed it the triangle
/// indices a BVH walk produced. Candidates may repeat and may arrive in any
/// order; the result does not depend on that order. A candidate index at or
/// beyond `mesh.triangle_count()` is skipped, not refused; a triangle that
/// references a missing position is refused with
/// [`RayMeshError::PositionIndexOutOfRange`].
///
/// # Determinism
///
/// Hits are ordered by `t`, then by triangle index. Two coplanar triangles
/// sharing an edge therefore resolve to the same triangle on every run and on
/// every platform, instead of depending on traversal order.
pub fn nearest_hit_among(
    mesh: &impl TriangleMeshView,
    ray: &Ray3,
    tolerance: Tolerance,
    candidates: impl IntoIterator<Item = usize>,
) -> Result<Option<RayHit3>, RayMeshError> {
    validate_ray(ray)?;
    validate_tolerance(tolerance)?;

    let mut best: Option<RayHit3> = None;
    for triangle in candidates {
        if triangle >= mesh.triangle_count() {
            continue;
        }
        let Some(hit) = triangle_hit(mesh, ray, tolerance, triangle)? else {
            continue;
        };
        if best.is_none_or(|current| is_closer(&hit, &current)) {
            best = Some(hit);
        }
    }
    Ok(best)
}

/// Intersect one triangle of `mesh`, reporting the hit or a certified miss.
pub fn triangle_hit(
    mesh: &impl TriangleMeshView,
    ray: &Ray3,
    tolerance: Tolerance,
    triangle: usize,
) -> Result<Option<RayHit3>, RayMeshError> {
    validate_ray(ray)?;
    validate_tolerance(tolerance)?;
    let corners = corners(mesh, triangle)?;
    intersect_triangle(ray, corners, tolerance, triangle)
}

/// Intersect a ray with a standalone triangle.
///
/// `triangle_index` only labels diagnostics; it is not used for geometry.
pub fn intersect_triangle(
    ray: &Ray3,
    corners: [Point3; 3],
    tolerance: Tolerance,
    triangle_index: usize,
) -> Result<Option<RayHit3>, RayMeshError> {
    validate_ray(ray)?;
    validate_tolerance(tolerance)?;
    if !corners.iter().all(|corner| corner.is_finite()) {
        return Err(RayMeshError::NonFiniteInput);
    }

    let [a, b, c] = corners;
    let edge1 = b - a;
    let edge2 = c - a;
    let normal = edge1.cross(edge2);
    // Exact zero area is a representation fact, not a tolerance question: a
    // degenerate triangle has no plane to intersect at any tolerance.
    if normal.length_squared() == 0.0 {
        return Err(RayMeshError::DegenerateTriangle {
            triangle: triangle_index,
        });
    }

    // Möller-Trumbore, double-sided. The determinant is compared against the
    // caller's linear tolerance scaled by the operand magnitudes, so a
    // parallel-in-plane ray is rejected in the model's units instead of against
    // a hidden epsilon.
    let pvec = ray.direction.cross(edge2);
    let determinant = edge1.dot(pvec);
    let parallel_bound = tolerance.linear() * edge1.length() * pvec.length();
    if determinant.abs() <= parallel_bound {
        return Ok(None);
    }

    let inverse = 1.0 / determinant;
    let tvec = ray.origin - a;
    let u = tvec.dot(pvec) * inverse;
    let qvec = tvec.cross(edge1);
    let v = ray.direction.dot(qvec) * inverse;
    let w = 1.0 - u - v;

    // Edge and vertex hits are kept: a ray grazing a shared edge must hit the
    // surface, not fall through it. The barycentric slack is the caller's
    // tolerance, not an invented constant.
    let slack = tolerance.linear();
    if u < -slack || v < -slack || w < -slack {
        return Ok(None);
    }

    let t = edge2.dot(qvec) * inverse;
    if t < 0.0 {
        return Ok(None);
    }
    if !t.is_finite() || !u.is_finite() || !v.is_finite() {
        return Err(RayMeshError::NonFiniteInput);
    }

    Ok(Some(RayHit3 {
        t,
        triangle: triangle_index,
        barycentric: [w, u, v],
        side: side_of(ray.origin, corners),
        point: ray.origin + ray.direction * t,
    }))
}

/// Certified side classification of the ray origin against a triangle plane.
///
/// `orient3d(a, b, c, d)` is positive when `d` lies opposite the side the
/// winding normal points to, so a positive origin sign means the ray reaches
/// the triangle from behind and strikes its back face.
fn side_of(origin: Point3, corners: [Point3; 3]) -> FaceSide {
    let [a, b, c] = corners;
    match orient3d(a, b, c, origin).sign() {
        Some(Sign::Positive) => FaceSide::Back,
        Some(Sign::Negative) => FaceSide::Front,
        Some(Sign::Zero) => FaceSide::Coplanar,
        // Non-finite coordinates are rejected before this point, and `Sign` is
        // `#[non_exhaustive]`, so anything unrecognised must not be guessed at.
        _ => FaceSide::Coplanar,
    }
}

fn is_closer(candidate: &RayHit3, current: &RayHit3) -> bool {
    match candidate.t.partial_cmp(&current.t) {
        Some(core::cmp::Ordering::Less) => true,
        Some(core::cmp::Ordering::Equal) => candidate.triangle < current.triangle,
        _ => false,
    }
}

fn corners(mesh: &impl TriangleMeshView, triangle: usize) -> Result<[Point3; 3], RayMeshError> {
    let indices = mesh.triangle(triangle);
    let mut corners = [Point3::ZERO; 3];
    for (slot, index) in corners.iter_mut().zip(indices) {
        let index = usize::try_from(index)
            .map_err(|_| RayMeshError::PositionIndexOutOfRange { triangle })?;
        if index >= mesh.position_count() {
            return Err(RayMeshError::PositionIndexOutOfRange { triangle });
        }
        *slot = mesh.position(index);
    }
    if !corners.iter().all(|corner| corner.is_finite()) {
        return Err(RayMeshError::NonFiniteInput);
    }
    Ok(corners)
}

fn validate_ray(ray: &Ray3) -> Result<(), RayMeshError> {
    if !ray.origin.is_finite() || !ray.direction.is_finite() {
        return Err(RayMeshError::NonFiniteInput);
    }
    if ray.direction.length_squared() == 0.0 {
        return Err(RayMeshError::ZeroDirection);
    }
    Ok(())
}

fn validate_tolerance(tolerance: Tolerance) -> Result<(), RayMeshError> {
    let linear = tolerance.linear();
    if !linear.is_finite() || linear < 0.0 {
        return Err(RayMeshError::InvalidTolerance);
    }
    Ok(())
}
