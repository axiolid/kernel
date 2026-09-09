//! Solid admissibility: what a mesh must satisfy to be a boolean operand.
//!
//! # Why this lives in the kernel
//!
//! Before this module, the only provider validated its own inputs inside its
//! adapter. That made the first backend the de-facto definition of "valid
//! solid": a second provider would have brought a second, silently different
//! definition, and callers would have seen admissibility change when dispatch
//! picked a different backend.
//!
//! Axiolid owns admissibility. The registry validates *before* dispatch, so a
//! provider never sees an operand the contract rejects. Providers must not
//! widen the set (accepting what Axiolid rejects) or narrow it (rejecting what
//! Axiolid accepts); the conformance suite checks both directions.

use axiolid_core::{Tolerance, Vec3};
use axiolid_mesh::TriMesh;

use axiolid_contracts::{GeomError, GeomResult};

/// How strictly an operand must be formed.
///
/// Levels are cumulative: each includes the checks below it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[non_exhaustive]
pub enum SolidRequirements {
    /// Indices in range, no NaN/infinite coordinates, at least one triangle.
    ///
    /// The floor. An operand failing this is malformed data, not geometry.
    Structural,
    /// Structural, plus a finite, non-zero enclosed signed volume.
    ///
    /// Volume accumulation that overflows or otherwise becomes non-finite is
    /// rejected rather than treated as evidence of a valid interior.
    ///
    /// A flat or self-cancelling shell has no interior, so no set operation on
    /// it has a defined meaning.
    Enclosing,
    /// Enclosing, plus outward orientation (positive signed volume).
    ///
    /// An inside-out operand would silently invert the operation, turning a
    /// difference into an intersection without any error.
    Oriented,
}

/// Why an operand was rejected, with the operand named.
///
/// `role` is `"subject"` or `"tool[3]"` so a caller learns *which* mesh was
/// wrong, not merely that some mesh was.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SolidRejection {
    /// Which operand failed.
    pub role: String,
    /// The level it failed at.
    pub level: SolidRequirements,
    /// Human-readable detail.
    pub detail: String,
}

impl SolidRequirements {
    /// Validate one operand at this level.
    pub fn validate(self, mesh: &TriMesh, role: &str) -> GeomResult<()> {
        mesh.validate_structure()
            .map_err(|error| GeomError::InvalidInput(format!("{role}: {error}")))?;
        if mesh.indices.is_empty() {
            return Err(GeomError::InvalidInput(format!(
                "{role}: mesh has no triangles"
            )));
        }
        if !mesh
            .positions
            .iter()
            .all(|point| point.x.is_finite() && point.y.is_finite() && point.z.is_finite())
        {
            return Err(GeomError::InvalidInput(format!(
                "{role}: mesh has a non-finite coordinate"
            )));
        }
        if self == Self::Structural {
            return Ok(());
        }

        let six_volume = six_signed_volume(mesh);
        if !six_volume.is_finite() {
            return Err(GeomError::InvalidInput(format!(
                "{role}: mesh has non-finite signed volume"
            )));
        }
        if six_volume == 0.0 {
            return Err(GeomError::Degenerate(format!(
                "{role}: mesh encloses zero signed volume, so it has no interior"
            )));
        }
        if self == Self::Enclosing {
            return Ok(());
        }

        if six_volume < 0.0 {
            return Err(GeomError::InvalidInput(format!(
                "{role}: mesh is inside-out (signed volume {:.6} < 0); \
                 boolean operations would silently invert",
                six_volume / 6.0
            )));
        }
        Ok(())
    }

    /// Validate a subject and its tools, naming each operand by index.
    pub fn validate_operands(self, subject: &TriMesh, tools: &[&TriMesh]) -> GeomResult<()> {
        self.validate(subject, "subject")?;
        for (index, tool) in tools.iter().enumerate() {
            self.validate(tool, &format!("tool[{index}]"))?;
        }
        Ok(())
    }
}

/// Six times the signed volume, computed about the mesh centroid.
///
/// Summing triple products of ABSOLUTE coordinates cancels catastrophically
/// when a small solid sits far from the origin: the terms scale with the
/// distance cubed while the answer stays the size of the solid, so the sign
/// becomes rounding noise and a well-formed mesh is misread as inside-out
/// (#99). Shifting to the centroid first makes every term the size of the
/// solid itself. Exact in real arithmetic -- signed volume is
/// translation-invariant -- and vastly better conditioned in f64.
///
/// The factor of six is left in deliberately: the sign and the zero test are
/// what matter, and dividing introduces rounding for no benefit.
pub(crate) fn six_signed_volume(mesh: &TriMesh) -> f64 {
    let n = mesh.positions.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for p in &mesh.positions {
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    let (cx, cy, cz) = (cx / n, cy / n, cz / n);
    let shifted = |i: u32| {
        let p = mesh.positions[i as usize];
        Vec3::new(p.x - cx, p.y - cy, p.z - cz)
    };
    mesh.indices
        .chunks_exact(3)
        .map(|triangle| {
            let a = shifted(triangle[0]);
            let b = shifted(triangle[1]);
            let c = shifted(triangle[2]);
            a.dot(b.cross(c))
        })
        .sum()
}

/// Enclosed volume of a validated operand, used by conformance invariants.
pub fn enclosed_volume(mesh: &TriMesh, _tolerance: Tolerance) -> f64 {
    six_signed_volume(mesh) / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::Point3;

    /// Unit cube at `origin`, outward-oriented.
    fn cube(origin: [f64; 3], size: f64) -> TriMesh {
        let (x, y, z) = (origin[0], origin[1], origin[2]);
        let s = size;
        let p = vec![
            Point3::new(x, y, z),
            Point3::new(x + s, y, z),
            Point3::new(x + s, y + s, z),
            Point3::new(x, y + s, z),
            Point3::new(x, y, z + s),
            Point3::new(x + s, y, z + s),
            Point3::new(x + s, y + s, z + s),
            Point3::new(x, y + s, z + s),
        ];
        let i = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ];
        TriMesh::new(p, i)
    }

    /// #99: a small solid far from the origin must not be called inside-out.
    ///
    /// About the origin the triple products are ~1e27 while the true answer
    /// is ~1e-18, so the sign was rounding noise and a valid cube was
    /// rejected. About the centroid every term is the size of the cube.
    #[test]
    fn a_small_solid_far_from_the_origin_is_not_inside_out() {
        for (offset, size) in [(1e9, 1e-6), (1e6, 1e-3), (1e12, 1.0)] {
            let far = cube([offset, offset, offset], size);
            let six = six_signed_volume(&far);
            assert!(
                six > 0.0,
                "cube of size {size} at {offset} must have positive signed volume, got {six}"
            );
            SolidRequirements::Oriented
                .validate(&far, "subject")
                .expect("a well-formed cube must validate wherever it sits");
        }
    }

    /// The guard must still catch a genuinely inverted solid.
    #[test]
    fn a_genuinely_inside_out_solid_is_still_refused() {
        let mut inverted = cube([1e9, 1e9, 1e9], 1e-6);
        for t in inverted.indices.chunks_exact_mut(3) {
            t.swap(1, 2);
        }
        assert!(six_signed_volume(&inverted) < 0.0);
        let err = SolidRequirements::Oriented
            .validate(&inverted, "subject")
            .expect_err("an inverted solid must still be refused");
        assert!(format!("{err}").contains("inside-out"));
    }

    /// Signed volume is translation-invariant, so placement must not
    /// change the measured magnitude beyond f64 noise.
    #[test]
    fn the_measured_volume_does_not_depend_on_placement() {
        let near = six_signed_volume(&cube([0.0, 0.0, 0.0], 1.0));
        let far = six_signed_volume(&cube([1e6, 1e6, 1e6], 1.0));
        assert!(
            (near - far).abs() <= 1e-6 * near.abs(),
            "near {near} vs far {far}"
        );
    }
}
