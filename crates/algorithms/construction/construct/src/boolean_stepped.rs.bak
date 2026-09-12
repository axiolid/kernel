//! Stepped union of two coaxial prisms (ADR 0050 follow-up).
//!
//! # Why a union with differing spans is not one prism
//!
//! [`boolean_prisms_exact`](crate::boolean_exact::boolean_prisms_exact)
//! refuses a union whose operands span different heights, because the
//! result is stepped: the cross-section CHANGES partway up, and a single
//! prism carries exactly one section.
//!
//! The refusal is honest but the shape is perfectly well defined. Cutting
//! the union at every height where an operand starts or stops leaves bands,
//! and WITHIN a band the active operand set is constant -- so each band is
//! a genuine prism whose section is the planar union of whatever is active
//! there. The stepped solid is that stack, and the decomposition is exact:
//! the planar work is the same overlay the single-prism path already uses.
//!
//! # What this module does and does not give you
//!
//! It returns the BANDS. Assembling them into one `ExactBRep` additionally
//! needs the ledge faces where the section changes, which is its own piece
//! of work; returning the exact decomposition is the honest half that a
//! caller can already use, and it is verifiable on its own terms.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Scalar, Tolerance, Vec2};
use axiolid_overlay::{overlay, FillRule, OverlayInput, OverlayOperation, Polygon, Ring};

use crate::boolean_exact::{unsupported, Prism};
use crate::BACKEND_ID;

/// One constant-section slab of a stepped result.
#[derive(Debug, Clone, PartialEq)]
pub struct Band {
    /// Cross-section rings: outer first, then holes.
    pub rings: Vec<Vec<axiolid_core::Point2>>,
    /// Base height of this slab.
    pub bottom: Scalar,
    /// Top height of this slab.
    pub top: Scalar,
}

/// Decompose a coaxial union into constant-section bands, bottom to top.
///
/// Returns one [`Band`] per height interval over which the set of active
/// operands does not change. A union whose operands span the same height
/// yields exactly one band, which is the case
/// [`boolean_prisms_exact`](crate::boolean_exact::boolean_prisms_exact)
/// already handles.
///
/// Operands that do not touch are refused: their union is two separate
/// solids, and a band stack describes one.
pub fn union_prisms_stepped(
    subject: &Prism,
    tool: &Prism,
    tolerance: Tolerance,
) -> GeomResult<Vec<Band>> {
    if tool.bottom > subject.top + tolerance.linear()
        || subject.bottom > tool.top + tolerance.linear()
    {
        return Err(unsupported(
            "stepped union of prisms that do not meet along the axis",
        ));
    }

    // Every height where an operand starts or stops is a potential
    // section change. Heights closer than tolerance are the SAME cut:
    // keeping both would emit a zero-thickness band that no solid can
    // represent.
    let mut cuts = vec![subject.bottom, subject.top, tool.bottom, tool.top];
    cuts.sort_by(|a, b| a.total_cmp(b));
    cuts.dedup_by(|a, b| tolerance.eq(*a, *b));

    let mut bands = Vec::with_capacity(cuts.len().saturating_sub(1));
    for pair in cuts.windows(2) {
        let (bottom, top) = (pair[0], pair[1]);
        // The midpoint decides membership: it is interior to the band, so
        // it cannot land on a boundary and give an ambiguous answer.
        let middle = 0.5 * (bottom + top);
        let in_subject = middle > subject.bottom && middle < subject.top;
        let in_tool = middle > tool.bottom && middle < tool.top;
        let rings = match (in_subject, in_tool) {
            (false, false) => continue,
            (true, false) => subject.rings.clone(),
            (false, true) => tool.rings.clone(),
            // Both active: the band's section is their planar union, which
            // is exactly the overlay the single-prism path performs.
            (true, true) => section_union(subject, tool, tolerance)?,
        };
        bands.push(Band { rings, bottom, top });
    }
    if bands.is_empty() {
        return Err(GeomError::Degenerate(
            "stepped union has no band of positive height".to_owned(),
        ));
    }
    Ok(bands)
}

/// Planar union of the two cross-sections.
fn section_union(
    subject: &Prism,
    tool: &Prism,
    tolerance: Tolerance,
) -> GeomResult<Vec<Vec<axiolid_core::Point2>>> {
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
        OverlayOperation::Union,
        FillRule::NonZero,
        tolerance,
    )
    .map_err(|error| GeomError::BackendContractViolation {
        backend: BACKEND_ID,
        detail: format!("stepped union cross-section overlay failed: {error:?}"),
    })?;
    if result.polygons.len() != 1 {
        return Err(unsupported(
            "stepped union band with a disconnected cross-section",
        ));
    }
    let polygon = &result.polygons[0];
    let mut rings = Vec::with_capacity(1 + polygon.holes.len());
    rings.push(polygon.outer.points.clone());
    for hole in &polygon.holes {
        rings.push(hole.points.clone());
    }
    Ok(rings)
}

fn to_polygons(prism: &Prism) -> Vec<Polygon> {
    let mut rings = prism.rings.iter();
    let outer = Ring {
        points: rings.next().cloned().unwrap_or_default(),
    };
    let holes = rings.map(|r| Ring { points: r.clone() }).collect();
    vec![Polygon { outer, holes }]
}
