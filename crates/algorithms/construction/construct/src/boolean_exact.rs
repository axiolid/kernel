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
    arc_overlay, overlay, validate_arc_ring, ArcPolygon, ArcRing, FillRule, OverlayInput,
    OverlayOperation, Polygon, Ring,
};

use crate::boolean_provenance::{name_side_fragment, OperandRings};
use axiolid_brep_audit::geometric_audit;

use crate::extrude_arc::extrude_arc_rings_from;
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
///
/// A result that falls apart into several separate solids is refused here,
/// because one `ExactBRep` is one solid and returning just one piece would
/// silently discard material. [`boolean_prisms_exact_solids`] returns every
/// piece.
pub fn boolean_prisms_exact(
    subject: &Prism,
    tool: &Prism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let (bottom, top) = prism_span(subject, tool, operator, tolerance)?;
    let polygons = prism_sections(subject, tool, operator, tolerance)?;
    single_solid(
        polygons,
        "prism boolean produced an empty cross-section",
        "exact prism boolean producing disconnected components \
         (boolean_prisms_exact_solids returns every piece)",
        |polygon| prism_solid(polygon, subject, tool, (bottom, top), tolerance),
    )
}

/// Exact boolean of two coaxial prisms, one exact B-rep per separate piece.
///
/// Same reduction and refusals as [`boolean_prisms_exact`], except that a
/// result which falls apart into several solids (a difference cutting a
/// wall in two, an intersection of two separate islands) returns every
/// piece instead of refusing, and an empty result is an empty list rather
/// than an error.
///
/// Solids are ordered by the lowest vertex of their outer boundary, `x`
/// first and then `y`, so the order is stable across runs and does not
/// depend on the overlay's internal traversal.
pub fn boolean_prisms_exact_solids(
    subject: &Prism,
    tool: &Prism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<Vec<ExactBRep>> {
    let Some((bottom, top)) = empty_span_is_none(prism_span(subject, tool, operator, tolerance))?
    else {
        return Ok(Vec::new());
    };
    // `overlay` already orders polygons by their lowest outer vertex (its
    // documented contract), which is the order the arc path sorts into.
    let polygons = prism_sections(subject, tool, operator, tolerance)?;
    polygons
        .iter()
        .map(|polygon| prism_solid(polygon, subject, tool, (bottom, top), tolerance))
        .collect()
}

/// Validate both prisms and settle the height span of the result.
///
/// Height logic decides whether a prism can represent the answer at all,
/// so it is settled before any planar work.
fn prism_span(
    subject: &Prism,
    tool: &Prism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<(Scalar, Scalar)> {
    validate(subject, "subject")?;
    validate(tool, "tool")?;
    resolve_span(
        (subject.bottom, subject.top),
        (tool.bottom, tool.top),
        operator,
        tolerance,
    )
}

/// The planar boolean of the two cross-sections, one polygon per piece.
fn prism_sections(
    subject: &Prism,
    tool: &Prism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<Vec<Polygon>> {
    let operation = overlay_operation(operator)?;
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
    Ok(result.polygons)
}

/// Extrude one result piece and name its faces after the operands.
fn prism_solid(
    polygon: &Polygon,
    subject: &Prism,
    tool: &Prism,
    (bottom, top): (Scalar, Scalar),
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
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

/// The one solid a single-result boolean may return.
///
/// Empty and multi-piece results are refused with the messages callers
/// already match on: one `ExactBRep` is one solid, and returning just one
/// piece would silently discard material.
fn single_solid<T>(
    pieces: Vec<T>,
    empty: &'static str,
    disconnected: &'static str,
    build: impl FnOnce(&T) -> GeomResult<ExactBRep>,
) -> GeomResult<ExactBRep> {
    match pieces.as_slice() {
        [] => Err(GeomError::Degenerate(empty.to_owned())),
        [only] => build(only),
        _ => Err(unsupported(disconnected)),
    }
}

/// Map an empty height span to `None`; every other refusal stays an error.
///
/// For the multi-solid entry points an empty result is a valid answer (an
/// empty list), not a failure.
fn empty_span_is_none(span: GeomResult<(Scalar, Scalar)>) -> GeomResult<Option<(Scalar, Scalar)>> {
    match span {
        Ok(span) => Ok(Some(span)),
        Err(GeomError::Degenerate(_)) => Ok(None),
        Err(error) => Err(error),
    }
}

/// Order boundaries by their lowest vertex, `x` then `y`.
fn lowest_first(a: &[Point2], b: &[Point2]) -> std::cmp::Ordering {
    let key = |ring: &[Point2]| {
        ring.iter()
            .copied()
            .min_by(|p, q| p.x.total_cmp(&q.x).then(p.y.total_cmp(&q.y)))
    };
    match (key(a), key(b)) {
        (Some(p), Some(q)) => p.x.total_cmp(&q.x).then(p.y.total_cmp(&q.y)),
        (a, b) => a.is_some().cmp(&b.is_some()),
    }
}

fn overlay_operation(operator: BooleanOperator) -> GeomResult<OverlayOperation> {
    match operator {
        BooleanOperator::Intersection => Ok(OverlayOperation::Intersection),
        BooleanOperator::Union => Ok(OverlayOperation::Union),
        BooleanOperator::Difference => Ok(OverlayOperation::Difference),
        _ => Err(unsupported("unknown exact prism boolean operator")),
    }
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
/// The planar part is the exact arc overlay (ADR 0070): every
/// topological decision is exact for the given input, and crossing points
/// of two curves are rounded to `f64` once, in the output. Results may
/// carry holes (through-openings) and may start above `z = 0`.
///
/// # Refused
///
/// A result with several disconnected regions (one `ExactBRep` is one
/// solid; [`boolean_arc_prisms_exact_solids`] returns every piece) and the
/// stepped spans [`boolean_prisms_exact`] also refuses.
pub fn boolean_arc_prisms_exact(
    subject: &ArcPrism,
    tool: &ArcPrism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let span = arc_prism_span(subject, tool, operator, tolerance)?;
    let regions = arc_prism_sections(subject, tool, operator, tolerance)?;
    single_solid(
        regions,
        "arc prism boolean produced an empty cross-section",
        "exact arc prism boolean producing disconnected components \
         (boolean_arc_prisms_exact_solids returns every piece)",
        |region| arc_prism_solid(region, span, tolerance),
    )
}

/// Exact boolean of two coaxial arc prisms, one exact B-rep per piece.
///
/// Same reduction and refusals as [`boolean_arc_prisms_exact`], except that
/// a result which falls apart into several solids returns every piece, and
/// an empty result is an empty list rather than an error. Solids are
/// ordered as in [`boolean_prisms_exact_solids`]: by the lowest vertex of
/// their outer boundary, `x` first and then `y`.
pub fn boolean_arc_prisms_exact_solids(
    subject: &ArcPrism,
    tool: &ArcPrism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<Vec<ExactBRep>> {
    let Some(span) = empty_span_is_none(arc_prism_span(subject, tool, operator, tolerance))? else {
        return Ok(Vec::new());
    };
    let mut regions = arc_prism_sections(subject, tool, operator, tolerance)?;
    regions.sort_by(|a, b| lowest_first(&arc_points(&a.outer), &arc_points(&b.outer)));
    regions
        .iter()
        .map(|region| arc_prism_solid(region, span, tolerance))
        .collect()
}

fn arc_points(ring: &ArcRing) -> Vec<Point2> {
    ring.vertices.iter().map(|vertex| vertex.point).collect()
}

/// Validate both arc prisms and settle the height span of the result.
fn arc_prism_span(
    subject: &ArcPrism,
    tool: &ArcPrism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<(Scalar, Scalar)> {
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
    resolve_span(
        (subject.bottom, subject.top),
        (tool.bottom, tool.top),
        operator,
        tolerance,
    )
}

/// The exact planar boolean of the two arc sections, one region per piece.
fn arc_prism_sections(
    subject: &ArcPrism,
    tool: &ArcPrism,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> GeomResult<Vec<ArcPolygon>> {
    let operation = overlay_operation(operator)?;
    let result =
        arc_overlay(&subject.section, &tool.section, operation, tolerance).map_err(|error| {
            GeomError::BackendContractViolation {
                backend: BACKEND_ID,
                detail: format!("arc prism cross-section overlay failed: {error:?}"),
            }
        })?;
    Ok(result.regions)
}

/// Extrude one arc region between the result's heights.
fn arc_prism_solid(
    region: &ArcPolygon,
    (bottom, top): (Scalar, Scalar),
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    // Holes are through-openings: each becomes its own wall ring and a
    // second bound on both caps. The result may start above z = 0 (an
    // intersection with a raised tool), so the section is extruded from
    // its own base height.
    let mut rings = Vec::with_capacity(1 + region.holes.len());
    rings.push(region.outer.clone());
    rings.extend(region.holes.iter().cloned());
    let solid = extrude_arc_rings_from(&rings, bottom, Vec3::Z * (top - bottom))?;
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
