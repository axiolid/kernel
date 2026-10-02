//! Revolution of a profile about an axis.
//!
//! The profile is flattened to rings in its own plane, then each ring point
//! is rotated about the axis in steps chosen by the chord budget. A full
//! turn closes on itself; a partial turn is capped by the profile at each
//! end.
//!
//! Pappus gives the oracle: a full revolution has volume `2*pi*R*A`, where
//! `R` is the centroid's distance from the axis and `A` the profile area.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;

use crate::profile::Rings;

/// Rotation of `p` about the axis through `origin` along unit `dir`.
///
/// Rodrigues' formula. Written out rather than pulled from a matrix type
/// so the axis stays arbitrary: a revolution is not a Z-up operation.
pub(crate) fn rotate(p: Point3, origin: Point3, dir: Vec3, angle: Scalar) -> Point3 {
    let v = p - origin;
    let (s, c) = angle.sin_cos();
    origin + v * c + dir.cross(v) * s + dir * (dir.dot(v) * (1.0 - c))
}

/// Angular steps so the swept arc meets the chord budget.
///
/// The widest point of the profile governs: a ring point at distance `r`
/// from the axis traces a circle of that radius, and its sagitta is
/// `r(1 - cos(dtheta/2))`. Using the maximum radius means every other
/// point is sampled at least as finely.
pub(crate) fn steps(max_radius: Scalar, angle: Scalar, tol: Scalar) -> usize {
    if !(max_radius.is_finite() && max_radius > 0.0 && tol.is_finite() && tol > 0.0) {
        return 8;
    }
    let ratio = (1.0 - (tol / max_radius).min(1.0)).clamp(-1.0, 1.0);
    let per = 2.0 * ratio.acos().max(1e-9);
    ((angle.abs() / per).ceil() as usize).clamp(2, 4096)
}

/// Revolve a profile about an axis into a closed solid.
///
/// A full turn wraps its rings by index, so the seam shares vertices by
/// construction rather than by two samplings agreeing numerically. That is
/// the same structure `tessellate_primitive` uses for a cylinder, and it
/// avoids the trim-based path recorded as broken in issue #2.
pub fn revolve(
    rings: &Rings,
    axis_origin: Point3,
    axis_direction: Vec3,
    angle: Scalar,
    tolerance: Tolerance,
) -> GeomResult<TriMesh> {
    if !angle.is_finite() || angle == 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "revolution angle must be finite and non-zero, got {angle}"
        )));
    }
    let len = axis_direction.length();
    if !len.is_finite() || len <= 0.0 {
        return Err(GeomError::InvalidInput(
            "revolution axis must be a finite non-zero direction".to_owned(),
        ));
    }
    let dir = axis_direction / len;
    // Widest distance from the axis governs the angular step: a point at
    // radius r traces a circle of that radius.
    let mut max_r: Scalar = 0.0;
    for p in rings.outer.iter().chain(rings.holes.iter().flatten()) {
        let v = Point3::new(p.x, p.y, 0.0) - axis_origin;
        max_r = max_r.max((v - dir * dir.dot(v)).length());
    }
    let full = (angle.abs() - core::f64::consts::TAU).abs() <= 1e-9;
    let n = steps(max_r, angle, tolerance.linear());
    // A full turn emits n stations and wraps onto the first; a partial turn
    // emits n + 1 so both ends exist to be capped.
    let count = if full { n } else { n + 1 };
    let stations: Vec<crate::loft::Station> = (0..count)
        .rev()
        .map(|s| {
            let t = angle * (s as Scalar) / (n as Scalar);
            crate::loft::place(rings, |p| {
                rotate(Point3::new(p.x, p.y, 0.0), axis_origin, dir, t)
            })
        })
        .collect();
    let mut mesh = crate::loft::loft(rings, &stations, full)?;
    // The station order above (end of the sweep first) winds the walls
    // outward for a positive angle about `dir`. A negative angle sweeps the
    // other way round, so the same order winds them inward instead (#221):
    // `angle` and `dir` only matter through their product (the axis-angle
    // rotation vector), and `(dir, angle)` and `(-dir, -angle)` describe the
    // identical rotation, so the fix cannot key off `angle`'s stored sign
    // alone without also depending on `dir`'s. Settling it from the built
    // mesh's own orientation sidesteps that: whichever combination produced
    // it, a closed solid's sign is unambiguous.
    if matches!(crate::extrude::outward_orientation(&mesh), Some(false)) {
        for triangle in mesh.indices.chunks_exact_mut(3) {
            triangle.swap(1, 2);
        }
    }
    Ok(mesh)
}
