//! Deterministic raw triangle-mesh measures.
use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_mesh::{audit_mesh, MeshHealth, TriangleMeshView};
use core::fmt;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceProperties {
    pub area: f64,
    pub centroid: Point3,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumeProperties {
    pub signed_volume: f64,
    pub centroid: Point3,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshMeasureError {
    MeshNotSurfaceUsable(MeshHealth),
    MeshNotVolumeUsable(MeshHealth),
    ZeroSurfaceArea,
    ZeroSignedVolume,
}
impl fmt::Display for MeshMeasureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("mesh is not suitable for requested raw measure")
    }
}
impl std::error::Error for MeshMeasureError {}

pub fn surface_properties<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    tolerance: Tolerance,
) -> Result<SurfaceProperties, MeshMeasureError> {
    let health = audit_mesh(mesh, tolerance);
    if !health.is_surface_usable() || health.degenerate_triangles != 0 {
        return Err(MeshMeasureError::MeshNotSurfaceUsable(health));
    }
    // Area uses edge differences, which are already local and well
    // conditioned. The centroid is not: it accumulates absolute positions,
    // so at large coordinates it loses the same precision `volume_properties`
    // did. Re-base it the same way.
    let base = mesh.position(0);
    let (mut area, mut weighted) = (0.0f64, Vec3::ZERO);
    for i in 0..mesh.triangle_count() {
        let [a, b, c] = mesh.triangle(i).map(|j| mesh.position(j as usize) - base);
        let weight = (b - a).cross(c - a).length() * 0.5;
        area += weight;
        weighted += (a + b + c) * (weight / 3.0);
    }
    if !area.is_finite() || area == 0.0 {
        return Err(MeshMeasureError::ZeroSurfaceArea);
    }
    Ok(SurfaceProperties {
        area,
        // Shift the local-frame centroid back to world coordinates.
        centroid: base + (weighted / area),
    })
}

pub fn volume_properties<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    tolerance: Tolerance,
) -> Result<VolumeProperties, MeshMeasureError> {
    let health = audit_mesh(mesh, tolerance);
    if !health.is_closed_two_manifold() {
        return Err(MeshMeasureError::MeshNotVolumeUsable(health));
    }
    // Sum about a local origin rather than the world origin.
    //
    // Volume and centroid are translation-invariant, so re-basing the operands
    // changes nothing mathematically -- but it changes the conditioning
    // completely. Evaluating `a . (b x c)` about the world origin scales every
    // term with the CUBE of the distance to it, so a 0.1 m box on a national
    // grid at 1e7 forms terms of order 1e21 whose cancellation has to produce
    // 1e-3. That measured 25% volume error and a centroid 8160 km off the
    // solid: a rule check reading either got a confidently wrong number.
    //
    // The first vertex is the local origin. Any point near the geometry works;
    // this one is free and always on the mesh.
    let base = mesh.position(0);
    let (mut volume, mut weighted) = (0.0f64, Vec3::ZERO);
    for i in 0..mesh.triangle_count() {
        let [a, b, c] = mesh.triangle(i).map(|j| mesh.position(j as usize) - base);
        let tetra = a.dot(b.cross(c)) / 6.0;
        volume += tetra;
        weighted += (a + b + c) * (tetra / 4.0);
    }
    if !volume.is_finite() || volume == 0.0 {
        return Err(MeshMeasureError::ZeroSignedVolume);
    }
    Ok(VolumeProperties {
        signed_volume: volume,
        // The sum ran in local coordinates, so the centroid comes back
        // relative to `base`. Shift it home.
        centroid: base + (weighted / volume),
    })
}

/// Second moments of the enclosed solid about the origin.
///
/// The `x`, `y`, `z` components are the integrals of `x^2`, `y^2` and
/// `z^2` over the enclosed volume, i.e. the diagonal of the second-moment
/// tensor for unit density. The classical inertia tensor diagonal is
/// `(Iyy + Izz, Ixx + Izz, Ixx + Iyy)` from these, so callers can derive
/// either convention without the provider guessing which one they meant.
///
/// # Method
///
/// The divergence theorem again, one order higher than the volume sum:
/// for a closed oriented triangulation the integral of `x^2` over the
/// enclosed region reduces to a sum over triangles of terms in the
/// vertices' coordinates. Each triangle contributes
/// `(nx / 60) * sum over the 10 symmetric monomials`, the standard
/// closed-form tetrahedral moment about the origin.
///
/// Requires the same closed two-manifold input as [`volume_properties`]:
/// an open shell has no enclosed region and the integral is meaningless,
/// so it is refused rather than summed into a plausible number.
pub fn second_moments<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    tolerance: Tolerance,
) -> Result<Vec3, MeshMeasureError> {
    let health = audit_mesh(mesh, tolerance);
    if !health.is_closed_two_manifold() {
        return Err(MeshMeasureError::MeshNotVolumeUsable(health));
    }
    let mut moments = Vec3::ZERO;
    for i in 0..mesh.triangle_count() {
        let [a, b, c] = mesh.triangle(i).map(|j| mesh.position(j as usize));
        // Signed volume of the tetrahedron on the origin, times 6.
        let six_v = a.dot(b.cross(c));
        for axis in 0..3 {
            let (pa, pb, pc) = (a[axis], b[axis], c[axis]);
            // Integral of t^2 over the tetrahedron, in barycentric closed
            // form: (a^2 + b^2 + c^2 + ab + ac + bc) / 60 per unit 6V.
            let quadratic = pa * pa + pb * pb + pc * pc + pa * pb + pa * pc + pb * pc;
            moments[axis] += six_v * quadratic / 60.0;
        }
    }
    if !moments.x.is_finite() || !moments.y.is_finite() || !moments.z.is_finite() {
        return Err(MeshMeasureError::ZeroSignedVolume);
    }
    Ok(moments)
}
