//! Mass properties of an exact B-rep, without tessellating it.
//!
//! # Why this exists separately
//!
//! [`crate::mesh_measure::MeshMeasure`] measures a `TriMesh`. An exact B-rep
//! has no triangles, so measuring one previously meant tessellating first --
//! which converts an exact solid into an approximation before measuring it,
//! and then reports the approximation's volume as though it were the solid's.
//!
//! This path measures the exact representation directly. For a planar face
//! bounded by straight edges the divergence theorem is exact over the
//! polygonal boundary, so a prism's volume comes back at machine precision
//! rather than at tessellation fidelity.
//!
//! A curved face, or a planar face with a curved edge, is integrated over its
//! own parameter domain by Green's theorem round the face's pcurves
//! (module `exact_face`): the exact surface and the exact trimming curves,
//! with adaptive Gauss-Kronrod quadrature held to a relative error near
//! machine precision. Nothing is faceted.
//!
//! # Deliberate refusal
//!
//! An unknown surface family, a face whose boundary does not enclose a domain
//! in its parameters, and an integral that does not converge are refused by
//! name. Approximating any of them would silently reintroduce exactly the
//! tessellation error this path exists to avoid.

use axiolid_brep::ExactBRep;
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::Curve3;
use axiolid_surface::Surface;
use axiolid_topology::{Face, LoopId, Orientation};
use core::fmt;

use crate::MassProperties;

/// Components integrated per face: area, volume, three first moments and
/// three second moments.
pub(crate) const COMPONENTS: usize = 8;

/// One value per component, in [`COMPONENTS`] order.
pub(crate) type Sums = [Scalar; COMPONENTS];

/// Why an exact B-rep could not be measured.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExactMeasureError {
    /// A face's support surface belongs to a family this module cannot
    /// integrate.
    NonPlanarFace(&'static str),
    /// A face has no support surface attached.
    MissingSurface,
    /// An edge referenced geometry the B-rep does not contain.
    DanglingReference,
    /// The boundary enclosed no volume.
    Degenerate,
    /// A face's pcurves do not bound a domain in its surface parameters.
    ParameterDomain(&'static str),
    /// A surface or pcurve could not be evaluated where the face needs it.
    Evaluation,
    /// A face integral did not reach its error bound.
    NotConverged,
}

impl fmt::Display for ExactMeasureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPlanarFace(kind) => write!(
                f,
                "exact measurement cannot integrate a {kind} face. \
                 Tessellate and use MeshMeasure for an approximate answer."
            ),
            Self::MissingSurface => f.write_str("a face has no support surface"),
            Self::DanglingReference => {
                f.write_str("a boundary element references missing geometry")
            }
            Self::Degenerate => f.write_str("the boundary encloses no volume"),
            Self::ParameterDomain(why) => write!(f, "a face cannot be integrated: {why}"),
            Self::Evaluation => f.write_str("a surface or pcurve could not be evaluated on a face"),
            Self::NotConverged => f.write_str(
                "a curved face integral did not reach its error bound; \
                 tessellate and use MeshMeasure for an approximate answer",
            ),
        }
    }
}

impl std::error::Error for ExactMeasureError {}

/// Name the surface family for a refusal that tells the caller what to add.
pub(crate) fn family(surface: &Surface) -> &'static str {
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

/// Mass properties of an exact B-rep.
///
/// # Method
///
/// A planar face bounded by straight edges is a polygon in 3-space. Fanning
/// it into triangles about its first vertex and applying the divergence
/// theorem gives volume, centroid and second moments in one pass -- the same
/// closed forms [`crate::mesh`] uses, but over the B-rep's own boundary
/// polygons rather than over a tessellation of them. For planar faces the fan
/// is not an approximation: a planar polygon is exactly the union of its fan
/// triangles.
///
/// Every other face -- cylinders, cones, spheres, tori, elliptical
/// cylinders, B-spline surfaces, and planar faces with an arc or ellipse on
/// their boundary -- contributes the same cone integrals over its own
/// parameter domain, by Green's theorem round its pcurves (see
/// `exact_face`). Both paths use the same fields, so a solid mixing
/// them sums consistently.
///
/// Orientation is honoured at every level: a face marked
/// [`Orientation::Reversed`], a reversed use of a face in its shell, and a
/// reversed bound each flip the loop's winding, exactly as `audit_brep`
/// reads them. A correctly built solid yields a positive volume wherever it
/// sits, without the caller pre-normalising anything. Every shell is summed,
/// so a void shell, used reversed, subtracts its volume.
///
/// # Errors
///
/// Refuses an unknown surface family, a face whose pcurves do not bound a
/// domain, an integral that does not converge, and a boundary that encloses
/// nothing -- never approximating any of them.
pub fn exact_properties(
    brep: &ExactBRep,
    tolerance: Tolerance,
) -> Result<MassProperties, ExactMeasureError> {
    let topology = brep.topology();
    let scale = characteristic_length(brep);
    // Consecutive pcurves must meet on the surface. The caller's tolerance
    // governs, with a floor for `Tolerance::ZERO` so rounding in the pcurve
    // evaluation itself is not read as a gap.
    let linear = tolerance.linear().max(1e-9 * scale);
    let mut area = 0.0;
    let mut volume = 0.0;
    let mut volume_weighted = Point3::ZERO;
    let mut moments = Vec3::ZERO;

    // How each face is used by the shell that holds it. A face outside
    // every shell (a bare face table) is taken as used forward.
    let mut shell_sense = vec![Orientation::Forward; topology.faces().len()];
    for shell in topology.shells() {
        for &(face_id, sense) in &shell.faces {
            if let Some(slot) = shell_sense.get_mut(face_id.index()) {
                *slot = sense;
            }
        }
    }

    for (face_index, face) in topology.faces().iter().enumerate() {
        let surface_id = face.surface.ok_or(ExactMeasureError::MissingSurface)?;
        let surface = brep
            .surfaces()
            .get(surface_id.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        let flip_face = (shell_sense[face_index] == Orientation::Reversed)
            ^ (face.orientation == Orientation::Reversed);

        if matches!(surface, Surface::Plane(_)) && straight_edged(brep, face)? {
            let mut vector_area = Vec3::ZERO;
            for bound in &face.bounds {
                let mut ring = ring_positions(brep, bound.loop_id)?;
                if ring.len() < 3 {
                    continue;
                }
                // Loops are wound in the support surface's own frame; the
                // face, its use in the shell, and the bound each flip that.
                // This is the convention `audit_brep` checks when it pairs
                // edge uses, so any audited closed solid measures correctly
                // under it.
                //
                // Ignoring the flips looked right for years because every
                // tested solid had its reversed faces in the plane z = 0,
                // where `int z n_z dA` is zero whichever way the face is
                // wound. A solid lifted off that plane exposes it
                // (`a_raised_solid_measures_the_same_as_one_on_the_ground`).
                if flip_face ^ (bound.orientation == Orientation::Reversed) {
                    ring.reverse();
                }
                accumulate_fan(
                    &ring,
                    &mut vector_area,
                    &mut volume,
                    &mut volume_weighted,
                    &mut moments,
                );
            }
            // Summed as vectors, a hole's opposite winding subtracts its
            // area, and a non-convex fan's overlapping triangles cancel.
            area += vector_area.length();
            continue;
        }

        let sums = crate::exact_face::face_sums(brep, face, surface, linear, scale)?;
        let sign = if flip_face { -1.0 } else { 1.0 };
        area += sums[0].abs();
        volume += sign * sums[1];
        volume_weighted += Vec3::new(sums[2], sums[3], sums[4]) * sign;
        moments += Vec3::new(sums[5], sums[6], sums[7]) * sign;
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

/// Whether every edge on the face's boundary is a straight line, so its
/// vertices alone state the boundary.
fn straight_edged(
    brep: &ExactBRep,
    face: &Face<axiolid_brep::SurfaceId>,
) -> Result<bool, ExactMeasureError> {
    let topology = brep.topology();
    for bound in &face.bounds {
        let wire = topology
            .loops()
            .get(bound.loop_id.index())
            .ok_or(ExactMeasureError::DanglingReference)?;
        for use_ in &wire.edges {
            let edge = topology
                .edges()
                .get(use_.edge.index())
                .ok_or(ExactMeasureError::DanglingReference)?;
            let curve = edge
                .curve
                .and_then(|id| brep.curves3().get(id.index()))
                .ok_or(ExactMeasureError::DanglingReference)?;
            if !matches!(curve, Curve3::Line(_)) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// A length on the solid's own scale, measured from the origin the cone
/// fields are taken about: the noise floor of every integral is set from it.
fn characteristic_length(brep: &ExactBRep) -> Scalar {
    let reach = brep
        .topology()
        .vertices()
        .iter()
        .map(|vertex| vertex.position.length())
        .fold(0.0, Scalar::max);
    let extent = brep
        .surfaces()
        .iter()
        .map(|surface| match surface {
            Surface::Cylinder(c) => c.frame.origin.length() + c.radius,
            Surface::EllipticalCylinder(c) => {
                c.frame.origin.length() + c.semi_axis_x.max(c.semi_axis_y)
            }
            Surface::Cone(c) => c.frame.origin.length() + c.radius.abs(),
            Surface::Sphere(s) => s.frame.origin.length() + s.radius,
            Surface::Torus(t) => t.frame.origin.length() + t.major_radius + t.minor_radius,
            _ => 0.0,
        })
        .filter(|value| value.is_finite())
        .fold(0.0, Scalar::max);
    let length = reach.max(extent);
    if length > 0.0 {
        length
    } else {
        1.0
    }
}

/// Accumulate one polygon's fan into the running sums.
fn accumulate_fan(
    ring: &[Point3],
    area: &mut Vec3,
    volume: &mut Scalar,
    volume_weighted: &mut Point3,
    moments: &mut Vec3,
) {
    let anchor = ring[0];
    for window in ring[1..].windows(2) {
        let (a, b, c) = (anchor, window[0], window[1]);

        *area += (b - a).cross(c - a) * 0.5;

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
