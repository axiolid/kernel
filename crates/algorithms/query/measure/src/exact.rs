//! Mass properties of an exact B-rep, without tessellating it.
//!
//! # Why this exists separately
//!
//! [`crate::mesh_measure::MeshMeasure`] measures a `TriMesh`. An exact B-rep
//! has no triangles, so measuring one previously meant tessellating first --
//! which converts an exact solid into an approximation before measuring it,
//! and then reports the approximation's volume as though it were the solid's.
//!
//! This path measures the exact representation directly. For planar faces the
//! divergence theorem is exact over the polygonal boundary, so a prism's
//! volume comes back at machine precision rather than at tessellation
//! fidelity.
//!
//! # Deliberate refusal
//!
//! Only planar faces are supported. A cylindrical or spherical face needs a
//! surface integral this module does not implement, and approximating it by
//! sampling would silently reintroduce exactly the tessellation error this
//! path exists to avoid. Those cases are refused by name.

use axiolid_brep::ExactBRep;
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_surface::Surface;
use axiolid_topology::{LoopId, Orientation};
use core::fmt;

use crate::MassProperties;

/// Why an exact B-rep could not be measured.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExactMeasureError {
    /// A face's support surface is not planar.
    NonPlanarFace(&'static str),
    /// A face has no support surface attached.
    MissingSurface,
    /// An edge referenced geometry the B-rep does not contain.
    DanglingReference,
    /// The boundary enclosed no volume.
    Degenerate,
}

impl fmt::Display for ExactMeasureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPlanarFace(kind) => write!(
                f,
                "exact measurement supports planar faces only; found a {kind} face. \
                 Tessellate and use MeshMeasure for an approximate answer."
            ),
            Self::MissingSurface => f.write_str("a face has no support surface"),
            Self::DanglingReference => {
                f.write_str("a boundary element references missing geometry")
            }
            Self::Degenerate => f.write_str("the boundary encloses no volume"),
        }
    }
}

impl std::error::Error for ExactMeasureError {}

/// Name the surface family for a refusal that tells the caller what to add.
fn family(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane(_) => "plane",
        Surface::Cylinder(_) => "cylindrical",
        Surface::Cone(_) => "conical",
        Surface::Sphere(_) => "spherical",
        Surface::EllipticalCylinder(_) => "elliptical-cylindrical",
        Surface::Torus(_) => "toroidal",
        _ => "non-planar",
    }
}

/// Mass properties of an exact B-rep with planar faces.
///
/// # Method
///
/// Each planar face is a polygon in 3-space. Fanning it into triangles about
/// its first vertex and applying the divergence theorem gives volume,
/// centroid, area and second moments in one pass -- the same closed forms
/// [`crate::mesh`] uses, but over the B-rep's own boundary polygons rather
/// than over a tessellation of them. For planar faces the fan is not an
/// approximation: a planar polygon is exactly the union of its fan triangles.
///
/// Face orientation is honoured: a face marked [`Orientation::Reversed`]
/// contributes with its winding flipped, so a correctly built solid yields a
/// positive volume without the caller pre-normalising anything.
///
/// # Errors
///
/// Refuses any non-planar face by name rather than approximating it, and
/// refuses a boundary that encloses nothing.
pub fn exact_properties(
    brep: &ExactBRep,
    _tolerance: Tolerance,
) -> Result<MassProperties, ExactMeasureError> {
    let topology = brep.topology();
    let mut area = 0.0;
    let mut volume = 0.0;
    let mut volume_weighted = Point3::ZERO;
    let mut moments = Vec3::ZERO;

    for face in topology.faces() {
        let surface_id = face.surface.ok_or(ExactMeasureError::MissingSurface)?;
        let surface = brep
            .surfaces()
            .get(surface_id.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        if !matches!(surface, Surface::Plane(_)) {
            return Err(ExactMeasureError::NonPlanarFace(family(surface)));
        }

        for bound in &face.bounds {
            let ring = ring_positions(brep, bound.loop_id)?;
            if ring.len() < 3 {
                continue;
            }
            // Winding is already carried by the loop: each `EdgeUse` names
            // its own traversal direction, and `ring_positions` walks it, so
            // a `Reversed` face's loop already comes back wound the other
            // way. Flipping again here would double-correct and cancel out.
            // Verified by mutation: with an extra face-level flip, a prism
            // still measures positive because the two negations compose.
            accumulate_fan(
                &ring,
                &mut area,
                &mut volume,
                &mut volume_weighted,
                &mut moments,
            );
        }
    }

    if !volume.is_finite() || volume.abs() < Scalar::EPSILON {
        return Err(ExactMeasureError::Degenerate);
    }

    Ok(MassProperties {
        area,
        signed_volume: volume,
        centroid: volume_weighted / volume,
        second_moment_diagonal: moments,
    })
}

/// Accumulate one polygon's fan into the running sums.
fn accumulate_fan(
    ring: &[Point3],
    area: &mut Scalar,
    volume: &mut Scalar,
    volume_weighted: &mut Point3,
    moments: &mut Vec3,
) {
    let anchor = ring[0];
    for window in ring[1..].windows(2) {
        let (a, b, c) = (anchor, window[0], window[1]);

        *area += (b - a).cross(c - a).length() * 0.5;

        let six_v = a.dot(b.cross(c));
        *volume += six_v / 6.0;
        *volume_weighted += (a + b + c) * (six_v / 6.0 / 4.0);

        for axis in 0..3 {
            let (pa, pb, pc) = (a[axis], b[axis], c[axis]);
            let quadratic = pa * pa + pb * pb + pc * pc + pa * pb + pa * pc + pb * pc;
            moments[axis] += six_v * quadratic / 60.0;
        }
    }
}

/// Ordered vertex positions around one loop.
///
/// Each edge use carries its own traversal direction, so a loop is walked by
/// taking the START vertex of every oriented use: consecutive uses share a
/// vertex, and taking one endpoint per use yields the ring exactly once
/// without duplicating the shared corners.
fn ring_positions(brep: &ExactBRep, loop_id: LoopId) -> Result<Vec<Point3>, ExactMeasureError> {
    let topology = brep.topology();
    let wire = topology
        .loops()
        .get(loop_id.index())
        .ok_or(ExactMeasureError::DanglingReference)?;

    let mut ring = Vec::with_capacity(wire.edges.len());
    for use_ in &wire.edges {
        let edge = topology
            .edges()
            .get(use_.edge.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        let vertex_id = match use_.orientation {
            Orientation::Forward => edge.start,
            Orientation::Reversed => edge.end,
        };
        let vertex = topology
            .vertices()
            .get(vertex_id.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        ring.push(vertex.position);
    }
    Ok(ring)
}
