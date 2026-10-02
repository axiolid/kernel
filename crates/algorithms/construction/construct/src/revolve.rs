//! Revolution of a profile about an axis.
//!
//! The profile is flattened to rings in its own plane, then each ring point
//! is rotated about the axis in steps chosen by the chord budget. A full
//! turn closes on itself; a partial turn is capped by the profile at each
//! end.
//!
//! Pappus gives the oracle: a full revolution has volume `2*pi*R*A`, where
//! `R` is the centroid's distance from the axis and `A` the profile area.
//!
//! # How far the walls sit from the surface (#231)
//!
//! [`revolve`] bounds the distance from every point of the surface swept by
//! the GIVEN rings to its triangles by `tolerance.linear()`. Rings that
//! chord a curved profile add their own deviation on top: a caller that
//! wants the exact surface of revolution within a budget `c` flattens the
//! profile to `c_p` and revolves to `c - c_p`. The reference compiler takes
//! half each, which minimises the triangle count for a fixed sum: the steps
//! round the profile and round the axis each grow as the inverse square
//! root of their share, so their product is smallest at an even split. The
//! two shares add because a rotation is an isometry: a point of the exact
//! surface is the rotation of a profile point within `c_p` of the ring
//! polygon, so it is within `c_p` of the same rotation of the polygon.
//!
//! Each ring segment `q_a -> q_b` sweeps, between the stations at angles
//! `t` and `t + h`, the surface `R(s) q(x)` with `q(x) = (1-x) q_a + x q_b`.
//! For fixed `x` the point `R(s) q(x)` runs along an arc of radius
//! `rho(x)`, its distance from the axis, and the chord between the arc's
//! ends is the row at `x` of the bilinear patch over the quad
//! `R(t) q_a, R(t) q_b, R(t+h) q_a, R(t+h) q_b` (a rotation is linear in
//! the point). An arc of angle `h` stays within its sagitta
//! `rho (1 - cos(h/2))` of its chord, and `rho(x)` is the norm of an affine
//! function, so convex in `x` and largest at a ring vertex: every point of
//! the swept surface is within `rho_max (1 - cos(h/2))` of the patch, with
//! `rho_max` the largest vertex distance from the axis over every ring --
//! the outer equator, for a torus, never the tube centre. The patch is in
//! turn within `loft::quad_deviation` of the two triangles the
//! loft cuts it into. When the axis lies in the profile's plane every quad
//! is a planar trapezoid (its two chords are images of parallel radial
//! vectors under the same `R(t+h) - R(t)`), convex while the segment stays
//! on one side of the axis, and that term is zero: the bound is the sagitta
//! alone, whatever the sign of the Gaussian curvature (the inner side of a
//! torus is saddle-shaped, but its trapezoids are still flat). An axis off
//! the profile's plane sweeps hyperboloidal quads; their measured twist is
//! added, and the step shrinks until the sum fits.
//!
//! Every span of a revolution is the first span rotated, so the twist is
//! measured once.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_mesh::TriMesh;

use crate::profile::Rings;

/// The most steps a revolution or a circle takes before it refuses the
/// budget, rather than meet it with a coarser mesh that breaks the bound
/// the caller asked for.
pub(crate) const MAX_STEPS: usize = 4096;

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
/// point is sampled at least as finely. A budget that needs more than
/// [`MAX_STEPS`] steps is refused, never met by fewer.
pub(crate) fn steps(max_radius: Scalar, angle: Scalar, tol: Scalar) -> GeomResult<usize> {
    if !(tol.is_finite() && tol > 0.0) {
        return Err(GeomError::InvalidInput(format!(
            "chord budget must be positive and finite, got {tol}"
        )));
    }
    if !(max_radius.is_finite() && max_radius > 0.0) {
        return Ok(8);
    }
    let ratio = (1.0 - (tol / max_radius).min(1.0)).clamp(-1.0, 1.0);
    let per = 2.0 * ratio.acos();
    let n = (angle.abs() / per).ceil();
    // NaN (a zero step) counts as too many.
    if n.is_nan() || n > MAX_STEPS as Scalar {
        return Err(too_many_steps());
    }
    Ok((n as usize).max(2))
}

/// The refusal for a budget that needs more than [`MAX_STEPS`] steps.
pub(crate) fn too_many_steps() -> GeomError {
    GeomError::BudgetExceeded {
        resource: "revolution angular steps",
    }
}

/// The next step count to try when `n` misses the budget: a quarter more,
/// so a near miss costs little and a far one converges in a few rounds.
pub(crate) fn grow(n: usize) -> GeomResult<usize> {
    let next = n + n / 4 + 1;
    if next > MAX_STEPS {
        return Err(too_many_steps());
    }
    Ok(next)
}

/// Steps for a revolution of `rings` through `angle` whose walls stay
/// within `budget` of the surface the rings sweep (see the module notes).
fn bounded_steps(
    rings: &Rings,
    axis_origin: Point3,
    dir: Vec3,
    angle: Scalar,
    max_r: Scalar,
    budget: Scalar,
) -> GeomResult<usize> {
    let at = |t: Scalar| {
        crate::loft::place(rings, |p| {
            rotate(Point3::new(p.x, p.y, 0.0), axis_origin, dir, t)
        })
    };
    let mut n = steps(max_r, angle, budget)?;
    loop {
        let h = angle / n as Scalar;
        let sagitta = max_r * (1.0 - (0.5 * h).cos());
        // The loft walks the stations from the far end back, so each span
        // runs from angle `t + h` to `t`; every span is this one rotated.
        let twist = crate::loft::span_deviation(&at(h), &at(0.0));
        if sagitta + twist <= budget {
            return Ok(n);
        }
        n = grow(n)?;
    }
}

/// Revolve a profile about an axis into a closed solid.
///
/// A full turn wraps its rings by index, so the seam shares vertices by
/// construction rather than by two samplings agreeing numerically. That is
/// the same structure `tessellate_primitive` uses for a cylinder, and it
/// avoids the trim-based path recorded as broken in issue #2.
///
/// Every point of the surface the rings sweep lies within
/// `tolerance.linear()` of the mesh (see the module notes for the proof);
/// the rings' own distance from a curved profile is the caller's to add.
/// A budget that needs more than 4096 steps round the axis is refused with
/// `BudgetExceeded`.
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
    let n = bounded_steps(rings, axis_origin, dir, angle, max_r, tolerance.linear())?;
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
