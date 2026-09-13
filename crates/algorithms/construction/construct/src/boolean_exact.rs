//! Exact boolean over axis-aligned prisms (#66).
//!
//! # Why this family, and why it is genuinely exact
//!
//! Exact boolean support previously reached only half-space-bounded
//! difference and intersection. Widening it to arbitrary solids needs a
//! general polyhedral boolean, which is a large algorithm in its own right.
//!
//! But there is a family where the general problem collapses to one the
//! kernel already solves EXACTLY: two prisms sharing an extrusion axis. Their
//! boolean is the 2D boolean of their cross-sections crossed with the boolean
//! of their height intervals. `axiolid-overlay` computes the planar part
//! exactly, so the result is exact -- not a tessellated approximation.
//!
//! That family is not a toy. A wall with a rectangular opening is exactly
//! this shape, and it is the dominant pattern in building models.
//!
//! # When the answer is NOT a prism
//!
//! The reduction only holds when the result is itself a prism:
//!
//! - Intersection: always. The heights intersect to one interval.
//! - Union: only when both operands span the same height. Otherwise the
//!   result is stepped -- two different cross-sections at two heights -- and
//!   a single prism cannot represent it.
//! - Difference: only when the tool spans at least the subject's full height.
//!   A tool ending mid-way leaves a stepped solid for the same reason.
//!
//! Those cases are refused rather than approximated by the nearest prism,
//! which would silently change the geometry.

use axiolid_brep::{ExactBRep, FaceName, Operand, SweptFace};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{BooleanOperator, Frame2, Point2, Scalar, Tolerance, Vec2, Vec3};
use axiolid_overlay::{
    arc_overlay, overlay, validate_arc_ring, ArcRing, FillRule, OverlayInput, OverlayOperation,
    Polygon, Ring,
};

use crate::boolean_provenance::{name_side_fragment, OperandRings};
use axiolid_brep_audit::geometric_audit;

use crate::extrude_arc::extrude_arc_ring;
use crate::extrude_exact::extrude_polygon_rings_named;
use crate::BACKEND_ID;

pub fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::MeshBoolean,
        input,
    }
}

/// A prism: a closed planar cross-section swept along +z.
///
/// Deliberately explicit rather than recovered from an `ExactBRep`. Reading a
/// prism back out of a general B-rep is its own inference problem, and
/// getting it wrong would silently mis-identify the operand.
#[derive(Debug, Clone, PartialEq)]
pub struct Prism {
    /// Outer boundary, counter-clockwise; further rings are holes.
    pub rings: Vec<Vec<Point2>>,
    /// Base height along z.
    pub bottom: Scalar,
    /// Top height along z.
    pub top: Scalar,
}

/// A prism whose cross-section may contain arc edges (ADR 0050).
///
/// Separate from [`Prism`] rather than replacing it: the polygon path is
/// exact through integer predicates, the arc path agrees with closed forms
/// to machine precision. Those are different guarantees, so a caller that
/// wants the stronger one must not be silently moved to the weaker one.
#[derive(Debug, Clone, PartialEq)]
pub struct ArcPrism {
    /// Cross-section boundary; arcs carry a non-zero bulge.
    pub section: ArcRing,
    /// Base height along z.
    pub bottom: Scalar,
    /// Top height along z.
    pub top: Scalar,
}

/// Exact boolean of two coaxial prisms.
///
/// Returns an exact B-rep, or a typed refusal naming why the result is not
/// itself a prism. Never falls back to a mesh.
pub fn boolean_prisms_exact(
    subject: &Prism,
    tool: &Prism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    validate(subject, "subject")?;
    validate(tool, "tool")?;

    // Height logic decides whether a prism can represent the answer at all,
    // so it is settled before any planar work.
    let (bottom, top) = resolve_span(
        (subject.bottom, subject.top),
        (tool.bottom, tool.top),
        operator,
        tolerance,
    )?;
    let operation = match operator {
        BooleanOperator::Intersection => OverlayOperation::Intersection,
        BooleanOperator::Union => OverlayOperation::Union,
        BooleanOperator::Difference => OverlayOperation::Difference,
        _ => return Err(unsupported("unknown exact prism boolean operator")),
    };

    let frame = Frame2 {
        origin: Vec2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    };
    let result = overlay(
        &OverlayInput {
            frame,
            polygons: to_polygons(subject),
        },
        &OverlayInput {
            frame,
            polygons: to_polygons(tool),
        },
        operation,
        FillRule::NonZero,
        tolerance,
    )
    .map_err(|error| GeomError::BackendContractViolation {
        backend: BACKEND_ID,
        detail: format!("exact prism cross-section overlay failed: {error:?}"),
    })?;

    if result.polygons.is_empty() {
        return Err(GeomError::Degenerate(
            "prism boolean produced an empty cross-section".to_owned(),
        ));
    }
    // A disconnected result is several solids, and one ExactBRep is one
    // solid. Returning just the first would silently discard material.
    if result.polygons.len() > 1 {
        return Err(unsupported(
            "exact prism boolean producing disconnected components",
        ));
    }

    let polygon = &result.polygons[0];
    let mut rings = Vec::with_capacity(1 + polygon.holes.len());
    rings.push(polygon.outer.points.clone());
    for hole in &polygon.holes {
        rings.push(hole.points.clone());
    }

    // `extrude_polygon_rings` builds from z = 0, so a band that does not start
    // there is not representable by it. Refusing is correct rather than
    // silently dropping the offset and returning a solid at the wrong height.
    if bottom.abs() > tolerance.linear() {
        return Err(unsupported(
            "exact prism boolean whose result does not start at z = 0",
        ));
    }
    // Every edge of an exact overlay lies on an edge of one input, so each
    // wall of the result is a fragment of an input wall and can say which.
    // Recovering it here -- from geometry, after the fact -- avoids
    // threading provenance through the overlay backend.
    let subject_rings = OperandRings {
        operand: Operand::Subject,
        rings: &subject.rings,
    };
    let tool_rings = OperandRings {
        operand: Operand::Tool,
        rings: &tool.rings,
    };
    let operands = [subject_rings, tool_rings];
    let mut solid =
        extrude_polygon_rings_named(&rings, Vec3::Z * (top - bottom), &mut |(start, end)| {
            name_side_fragment(start, end, &operands)
        })?;

    // The caps lie in the planes the height logic selected above, so they
    // are fragments of whichever operand supplied each bound.
    name_caps(&mut solid, subject, tool, bottom, top, tolerance);
    gate_geometry(solid, tolerance)
}

fn validate(prism: &Prism, role: &'static str) -> GeomResult<()> {
    if prism.rings.is_empty() {
        return Err(GeomError::InvalidInput(format!(
            "{role} prism has no cross-section rings"
        )));
    }
    for ring in &prism.rings {
        if ring.len() < 3 {
            return Err(GeomError::InvalidInput(format!(
                "{role} prism ring needs at least three points"
            )));
        }
        if !ring.iter().all(|p| p.x.is_finite() && p.y.is_finite()) {
            return Err(GeomError::InvalidInput(format!(
                "{role} prism ring has a non-finite point"
            )));
        }
    }
    if !prism.bottom.is_finite() || !prism.top.is_finite() {
        return Err(GeomError::InvalidInput(format!(
            "{role} prism heights must be finite"
        )));
    }
    if prism.top <= prism.bottom {
        return Err(GeomError::InvalidInput(format!(
            "{role} prism top must lie above its bottom"
        )));
    }
    Ok(())
}

fn to_polygons(prism: &Prism) -> Vec<Polygon> {
    let mut rings = prism.rings.iter();
    let outer = Ring {
        points: rings.next().cloned().unwrap_or_default(),
    };
    let holes = rings.map(|r| Ring { points: r.clone() }).collect();
    vec![Polygon { outer, holes }]
}

/// Name the result caps after the operand whose cap plane they lie in.
///
/// A coaxial boolean never tilts a cap, so each result cap is coplanar with
/// a cap of at least one operand. When both operands share the plane the
/// subject is named: the result is a fragment of both, and naming it after
/// the subject keeps the choice deterministic rather than order-dependent.
/// A cap matching neither operand is left unnamed rather than guessed.
fn name_caps(
    solid: &mut ExactBRep,
    subject: &Prism,
    tool: &Prism,
    bottom: Scalar,
    top: Scalar,
    tolerance: Tolerance,
) {
    let start = cap_operand(subject.bottom, tool.bottom, bottom, tolerance);
    let end = cap_operand(subject.top, tool.top, top, tolerance);
    solid.name_caps(
        start.map(|operand| FaceName::swept(SweptFace::StartCap).fragment(operand)),
        end.map(|operand| FaceName::swept(SweptFace::EndCap).fragment(operand)),
    );
}

/// Which operand a result cap height came from.
fn cap_operand(
    subject: Scalar,
    tool: Scalar,
    result: Scalar,
    tolerance: Tolerance,
) -> Option<Operand> {
    if tolerance.eq(subject, result) {
        Some(Operand::Subject)
    } else if tolerance.eq(tool, result) {
        Some(Operand::Tool)
    } else {
        None
    }
}

/// Exact boolean of two coaxial prisms whose sections may contain arcs.
///
/// # What is exact here
///
/// The height reduction is identical to [`boolean_prisms_exact`]: a
/// coaxial boolean is the planar boolean of the sections crossed with the
/// boolean of the height intervals. Arc edges survive as arcs, so a
/// cylindrical wall stays a `Cylinder` face rather than becoming a fan of
/// planar strips.
///
/// The planar part agrees with closed-form areas to machine precision
/// rather than being exact through integer predicates. That is a weaker
/// claim than the polygon path makes and is stated here so a caller can
/// choose deliberately.
pub fn boolean_arc_prisms_exact(
    subject: &ArcPrism,
    tool: &ArcPrism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    for (section, role) in [(&subject.section, "subject"), (&tool.section, "tool")] {
        validate_arc_ring(section, tolerance).map_err(|error| {
            GeomError::InvalidInput(format!("{role} arc prism section: {error:?}"))
        })?;
    }
    if !subject.bottom.is_finite()
        || !subject.top.is_finite()
        || !tool.bottom.is_finite()
        || !tool.top.is_finite()
    {
        return Err(GeomError::InvalidInput(
            "arc prism heights must be finite".to_owned(),
        ));
    }
    if subject.top <= subject.bottom || tool.top <= tool.bottom {
        return Err(GeomError::InvalidInput(
            "arc prism top must lie above its bottom".to_owned(),
        ));
    }

    let (bottom, top) = resolve_span(
        (subject.bottom, subject.top),
        (tool.bottom, tool.top),
        operator,
        tolerance,
    )?;

    let operation = match operator {
        BooleanOperator::Intersection => OverlayOperation::Intersection,
        BooleanOperator::Union => OverlayOperation::Union,
        BooleanOperator::Difference => OverlayOperation::Difference,
        _ => return Err(unsupported("unknown exact prism boolean operator")),
    };

    let result =
        arc_overlay(&subject.section, &tool.section, operation, tolerance).map_err(|error| {
            GeomError::BackendContractViolation {
                backend: BACKEND_ID,
                detail: format!("arc prism cross-section overlay failed: {error:?}"),
            }
        })?;

    if result.regions.is_empty() {
        return Err(GeomError::Degenerate(
            "arc prism boolean produced an empty cross-section".to_owned(),
        ));
    }
    // One ExactBRep is one solid, so several regions cannot be returned
    // without silently discarding material.
    if result.regions.len() > 1 {
        return Err(unsupported(
            "exact arc prism boolean producing disconnected components",
        ));
    }

    let region = &result.regions[0];
    // A hole needs a cap face carrying two bounds with arc loops, which
    // the arc extruder does not build. Refusing names the gap instead of
    // returning a solid with its opening filled in.
    if !region.holes.is_empty() {
        return Err(unsupported(
            "exact arc prism boolean whose result has an interior hole",
        ));
    }
    // The arc extruder builds from z = 0, so a band starting elsewhere
    // would come back at the wrong height.
    if bottom.abs() > tolerance.linear() {
        return Err(unsupported(
            "exact arc prism boolean whose result does not start at z = 0",
        ));
    }
    let solid = extrude_arc_ring(&region.outer, Vec3::Z * (top - bottom))?;
    gate_geometry(solid, tolerance)
}

/// The height span a coaxial boolean result occupies.
///
/// Shared by the polygon and arc paths: the height reduction does not
/// depend on what the cross-section looks like, so duplicating it would
/// invite the two paths to disagree about which spans are representable.
fn resolve_span(
    subject: (Scalar, Scalar),
    tool: (Scalar, Scalar),
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<(Scalar, Scalar)> {
    match operator {
        BooleanOperator::Intersection => {
            let bottom = subject.0.max(tool.0);
            let top = subject.1.min(tool.1);
            if top - bottom <= tolerance.linear() {
                return Err(GeomError::Degenerate(
                    "prism intersection is empty along the extrusion axis".to_owned(),
                ));
            }
            Ok((bottom, top))
        }
        BooleanOperator::Union => {
            // Differing spans give a stepped solid, which is not a prism.
            if !tolerance.eq(subject.0, tool.0) || !tolerance.eq(subject.1, tool.1) {
                return Err(unsupported(
                    "exact prism union with differing extrusion spans",
                ));
            }
            Ok(subject)
        }
        BooleanOperator::Difference => {
            // A tool that stops inside the subject leaves a step.
            if tool.0 > subject.0 + tolerance.linear() || tool.1 < subject.1 - tolerance.linear() {
                return Err(unsupported(
                    "exact prism difference with a tool shorter than the subject",
                ));
            }
            Ok(subject)
        }
        _ => Err(unsupported("unknown exact prism boolean operator")),
    }
}

/// Reject a boolean result whose geometry does not hold together.
///
/// # Why booleans specifically
///
/// A boolean is where pcurves get rebuilt against surfaces they did not
/// originally trim, so it is the operation most able to produce a solid that
/// is topologically perfect and geometrically wrong. Both properties held on
/// the cap loops of an earlier arc boolean, and nothing caught it.
///
/// The audit is cheap here because an exact analytic boolean returns a handful
/// of faces, not a mesh: a few evaluations per edge use.
fn gate_geometry(solid: ExactBRep, tolerance: Tolerance) -> GeomResult<ExactBRep> {
    let health = geometric_audit(&solid, tolerance);
    if health.is_consistent() {
        return Ok(solid);
    }
    // Report the measured deviation, not just the fact of failure: a caller
    // deciding whether their tolerance is wrong needs the number.
    let detail = match health.worst_error() {
        Some(error) => format!(
            "boolean result failed its geometric audit: {} defect(s), worst deviation {error:e}",
            health.defects().len()
        ),
        None => format!(
            "boolean result failed its geometric audit: {} defect(s)",
            health.defects().len()
        ),
    };
    Err(GeomError::BackendContractViolation {
        backend: BACKEND_ID,
        detail,
    })
}
