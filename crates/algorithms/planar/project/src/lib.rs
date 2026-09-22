#![forbid(unsafe_code)]
//! Projection of triangle meshes onto a plane, and prism intersection.
//!
//! # The bridge between the 3D and 2D halves
//!
//! Triangle meshes live on one side of this kernel and planar booleans on the
//! other. Both halves exist; the fold that connects them did not, so a consumer
//! holding a mesh had to write the projection itself and pick its own winding
//! and degeneracy conventions.
//!
//! # Degenerate triangles are dropped explicitly
//!
//! A triangle seen edge-on projects to a zero-area sliver. It contributes
//! nothing to a union, but silently discarding it hides the difference between
//! "this mesh is edge-on" and "this mesh is empty".
//! The count is therefore reported in [`ProjectionEvidence`].
//!
//! # A projection is geometry, not a footprint
//!
//! [`project_mesh`] answers exactly one question: which points of the plane
//! does this mesh cover. It does not decide *which* mesh to project, *which*
//! plane counts as the reference, or which parts of a building should have
//! been included. Those are the questions that turn a projection into a
//! footprint, a shadow, a clash silhouette, or a formwork outline, and they
//! have different answers per consumer.
//!
//! Concretely, two downstream users asking for "the footprint" will disagree
//! about overhangs, balconies, and structure below grade. A kernel that picked
//! one of those answers would be wrong for the other user and would hide the
//! choice inside a function whose name implied there was nothing to choose.
//! The same split applies elsewhere in this kernel: a clearance query reports
//! a distance and refuses to call it adequate.
//!
//! So the mechanical half lives here -- project, union, preserve holes, report
//! evidence -- and the naming, plane selection, and inclusion rules live with
//! the consumer that has a reason to prefer one rule over another.

use axiolid_core::{PlaneFrame, Point2, Tolerance};
use axiolid_mesh::TriangleMeshView;
use axiolid_overlay::{
    overlay, union_soup, FillRule, OverlayError, OverlayInput, OverlayOperation, Polygon, Ring,
};

/// Why a projection could not be produced.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProjectionError {
    /// The projection plane basis was not orthonormal or not finite.
    ///
    /// No longer produced: `Plane` validates on construction, so an invalid
    /// basis cannot reach a projection. Retained because the enum is public
    /// and removing a variant is a breaking change.
    InvalidPlane,
    /// A mesh position was not finite.
    NonFinitePosition,
    /// A triangle referenced a position index the mesh does not have.
    IndexOutOfRange,
    /// The planar stage rejected the projected geometry.
    Planar(OverlayError),
}

impl From<OverlayError> for ProjectionError {
    fn from(error: OverlayError) -> Self {
        Self::Planar(error)
    }
}

/// What a projection actually did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct ProjectionEvidence {
    /// Triangles read from the source mesh.
    pub input_triangles: usize,
    /// Triangles whose projection had zero area and were dropped.
    ///
    /// Reported rather than silently skipped: a fully edge-on mesh projects to
    /// nothing, and that is a different fact from an empty mesh.
    pub degenerate_triangles: usize,
    /// Polygons in the resulting union.
    pub output_polygons: usize,
    /// Inner boundary components across all output polygons.
    pub output_holes: usize,
}

/// A projected footprint.
#[derive(Debug, Clone, PartialEq)]
pub struct Projection {
    /// The union of the projected triangles, holes preserved.
    pub polygons: Vec<Polygon>,
    /// What the projection did.
    pub evidence: ProjectionEvidence,
}

/// An orthonormal projection plane.
///
/// Alias for the core [`PlaneFrame`], which validates orthonormality on
/// construction using the ANGULAR tolerance. Orthonormality is a dimensionless
/// property; deciding it with the linear tolerance made the same skewed basis
/// valid in millimetres and invalid in metres.
pub type Plane = PlaneFrame;

/// Project a triangle mesh onto `plane` and union the projected triangles.
///
/// The result is a polygon set with holes, not an outline or a hull: a mesh
/// with a through-hole projects to a polygon that still has the hole.
///
/// Triangles are unioned pairwise into an accumulator rather than handed to the
/// backend as one soup, because the planar validator rejects self-intersecting
/// input and a raw triangle soup routinely overlaps itself.
pub fn project_mesh<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    plane: Plane,
    tolerance: Tolerance,
) -> Result<Projection, ProjectionError> {
    let triangles = collect_triangles(mesh, plane, tolerance)?;
    let polygons = union_all(triangles.rings, tolerance)?;
    let evidence = ProjectionEvidence {
        input_triangles: mesh.triangle_count(),
        degenerate_triangles: triangles.degenerate,
        output_polygons: polygons.len(),
        output_holes: polygons.iter().map(|p| p.holes.len()).sum(),
    };
    Ok(Projection { polygons, evidence })
}

struct Collected {
    rings: Vec<Ring>,
    degenerate: usize,
}

fn collect_triangles<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    plane: Plane,
    tolerance: Tolerance,
) -> Result<Collected, ProjectionError> {
    let mut rings = Vec::new();
    let mut degenerate = 0;
    let positions = mesh.position_count();

    for index in 0..mesh.triangle_count() {
        let corners = mesh.triangle(index);
        let mut projected = [Point2::new(0.0, 0.0); 3];
        for (slot, corner) in corners.iter().enumerate() {
            let corner = usize::try_from(*corner).map_err(|_| ProjectionError::IndexOutOfRange)?;
            if corner >= positions {
                return Err(ProjectionError::IndexOutOfRange);
            }
            let point = mesh.position(corner);
            if !point.is_finite() {
                return Err(ProjectionError::NonFinitePosition);
            }
            projected[slot] = plane.project(point);
        }

        // Twice the signed area. A triangle seen edge-on lands at zero here,
        // and is counted rather than quietly skipped.
        let cross = (projected[1] - projected[0]).perp_dot(projected[2] - projected[0]);
        // The ring validator rejects a zero-area ring, so the threshold here
        // matches the area threshold it applies rather than a looser one.
        if cross.abs() <= 2.0 * tolerance.linear().powi(2) {
            degenerate += 1;
            continue;
        }

        // Back-facing triangles project to clockwise rings. Orientation in 3D
        // is not a planar fact: a solid presents both facings to any plane, so
        // both must contribute positively to the footprint.
        let mut points = projected.to_vec();
        if cross < 0.0 {
            points.reverse();
        }
        rings.push(Ring { points });
    }

    Ok(Collected { rings, degenerate })
}

fn plane_frame() -> axiolid_core::Frame2 {
    axiolid_core::Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: axiolid_core::Vec2::new(1.0, 0.0),
        y: axiolid_core::Vec2::new(0.0, 1.0),
    }
}

/// Union already-CCW triangle rings into a normalised polygon set.
///
/// Delegates to `union_soup`, which validates each ring on its own but does
/// not require the set to be mutually disjoint. Folding pairwise through
/// `overlay` instead fails: the accumulator becomes self-touching as soon as
/// two triangles share an edge, and `overlay` rejects a self-intersecting
/// operand.
fn union_all(rings: Vec<Ring>, tolerance: Tolerance) -> Result<Vec<Polygon>, ProjectionError> {
    Ok(union_soup(&rings, tolerance)?)
}

/// Intersect a mesh footprint with a vertical prism over a 2D region.
///
/// The prism is infinite along the plane normal, so this is exactly the planar
/// intersection of the mesh footprint with `region`. Naming it separately keeps
/// the caller from having to know that identity, and leaves room for a bounded
/// prism later without changing the call site.
pub fn intersect_prism<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    plane: Plane,
    region: &[Polygon],
    tolerance: Tolerance,
) -> Result<Projection, ProjectionError> {
    let footprint = project_mesh(mesh, plane, tolerance)?;
    let subject = OverlayInput {
        frame: plane_frame(),
        polygons: footprint.polygons,
    };
    let clip = OverlayInput {
        frame: plane_frame(),
        polygons: region.to_vec(),
    };
    let clipped = overlay(
        &subject,
        &clip,
        OverlayOperation::Intersection,
        FillRule::NonZero,
        tolerance,
    )?;
    let evidence = ProjectionEvidence {
        output_polygons: clipped.polygons.len(),
        output_holes: clipped.polygons.iter().map(|p| p.holes.len()).sum(),
        ..footprint.evidence
    };
    Ok(Projection {
        polygons: clipped.polygons,
        evidence,
    })
}
