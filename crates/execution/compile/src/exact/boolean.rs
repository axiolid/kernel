//! Exact booleans of placed operands (#228).
//!
//! A building model cuts openings from walls and slabs: a placed extrusion
//! minus one or more placed extrusions. A door or window opening is usually
//! extruded across the wall, perpendicular to the wall's own extrusion; a
//! shaft through a slab runs parallel to it. Both are compiled here, under
//! any rigid placement of either operand.
//!
//! # Dispatch
//!
//! - Two bare sharp rectangles extruded along `+z`, without placements,
//!   keep the prism path (`boolean_prisms_exact`): it decides the plan
//!   with integer predicates and takes every operator.
//! - A boolean of any operator but the symmetric difference (#376) runs the
//!   general exact B-rep boolean of ADR 0075
//!   (`axiolid_brep_boolean::boolean`) on the two operands' exact B-reps,
//!   each placed by [`ExactBRep::transformed`]. An operand, subject or
//!   tool, is whatever [`ExactCompilation::compile`] builds exactly (#378):
//!   an extrusion, a revolution (its cylinders, cones and tori), another
//!   boolean, or whatever else it builds exactly, under any rigid
//!   placement. Openings compose: `((wall - o1) - o2) - o3` is three such
//!   differences.
//!
//! The coaxial arc-prism path is not used for placed operands. It would
//! take only a single-ring subject in a shared frame, so it serves the first
//! opening of a slab and none after it, and the general boolean already
//! costs well under a millisecond per opening (see the crate's tests), so a
//! second route would add a second set of results to keep in agreement for
//! no measured gain.
//!
//! # Refused, by name
//!
//! An operand the exact compiler refuses (a primitive, a tapered
//! revolution, a mesh, a half-space as subject) is refused with its own
//! family named. So are a symmetric difference of placed operands, a scaled
//! or sheared placement, and every configuration the general boolean cannot
//! build exactly ([`BooleanError`]): `GeomError::UnsupportedInput` naming
//! the case. A boolean that leaves
//! nothing (a difference that removes the whole subject, disjoint
//! operands' intersection) is `GeomError::Degenerate`. Nothing is meshed.
//!
//! # Mesh operands, for a boolean's report only (#266)
//!
//! A compilation that flattens meshes (crate-internal, for
//! `crate::deviation`) reads an authored mesh, or a B-rep of planes and
//! warped faces without a surface, as a placed operand too, compiled to the
//! exact polyhedron of its triangles (`super::polyhedral`), as subject or
//! tool. The rest of the dispatch, and every refusal above, is
//! unchanged; the public compiler still refuses a mesh operand by name.

use axiolid_brep::ExactBRep;
use axiolid_brep_boolean::{boolean_with_report, BooleanError};
use axiolid_construct::boolean_exact::{boolean_prisms_exact, Prism};
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{BooleanOperator, Point2, Vec3};
use axiolid_model::{GeometryNode, NodeId, SolidOperation};
use axiolid_profile::Profile;

use super::{remap_construction_error, unsupported, ExactCompilation};

impl ExactCompilation<'_> {
    pub(super) fn compile_boolean(
        &mut self,
        left: NodeId,
        right: NodeId,
        operator: BooleanOperator,
    ) -> GeomResult<ExactBRep> {
        // A half-space tool is a clip (#234, `super::clip`).
        if let Some(tool) = self.clip_tool(right)? {
            return self.compile_clip(left, tool, operator);
        }
        if let (Some(subject), Some(tool)) = (self.bare_prism(left)?, self.bare_prism(right)?) {
            return boolean_prisms_exact(&subject, &tool, operator, self.options.tolerance())
                .map_err(remap_construction_error);
        }
        if operator == BooleanOperator::SymmetricDifference {
            return Err(unsupported("exact symmetric difference of placed operands"));
        }
        // Each operand is whatever the exact compiler compiles, read
        // through its placements (#378): a refusal is the operand's own,
        // named by its family.
        let subject = self.compile(left)?;
        let tool = self.compile(right)?;
        let (body, report) =
            boolean_with_report(&subject, &tool, operator, self.options.tolerance()).map_err(
                |error| match error {
                    BooleanError::EmptyResult => emptied(operator),
                    error => remap_boolean_error(error),
                },
            )?;
        self.boolean_report = Some(report);
        Ok(body)
    }

    /// The operand as a prism, when it is a bare sharp filled rectangle
    /// extruded along `+z`, and `None` for anything else.
    ///
    /// Recovering a prism from an arbitrary exact B-rep is its own inference
    /// problem, and guessing wrong would mis-identify the operand, so only
    /// this one unplaced shape is read as a prism.
    fn bare_prism(&self, id: NodeId) -> GeomResult<Option<Prism>> {
        let node = self.graph.get(id).ok_or_else(|| {
            GeomError::InvalidInput(format!("operand {id:?} does not belong to this graph"))
        })?;
        let GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction,
            depth,
        }) = *node
        else {
            return Ok(None);
        };
        if direction.normalize_or_zero().dot(Vec3::Z) < 1.0 - self.options.tolerance().linear() {
            return Ok(None);
        }
        let Some(GeometryNode::Profile(Profile::Rectangle(rectangle))) = self.graph.get(profile)
        else {
            return Ok(None);
        };
        if rectangle.thickness.is_some()
            || rectangle.outer_radius.is_some()
            || rectangle.inner_radius.is_some()
        {
            return Ok(None);
        }
        let (hx, hy) = (rectangle.x / 2.0, rectangle.y / 2.0);
        Ok(Some(Prism {
            rings: vec![vec![
                Point2::new(-hx, -hy),
                Point2::new(hx, -hy),
                Point2::new(hx, hy),
                Point2::new(-hx, hy),
            ]],
            bottom: 0.0,
            top: depth,
        }))
    }
}

/// A union, intersection or difference that leaves nothing: only the
/// difference removes "the subject".
fn emptied(operator: BooleanOperator) -> GeomError {
    GeomError::Degenerate(
        match operator {
            BooleanOperator::Difference => "exact boolean difference removes the whole subject",
            BooleanOperator::Intersection => "exact boolean intersection leaves nothing",
            BooleanOperator::Union => "exact boolean union leaves nothing",
            BooleanOperator::SymmetricDifference => {
                "exact boolean symmetric difference leaves nothing"
            }
            _ => "exact boolean leaves nothing",
        }
        .to_owned(),
    )
}

/// A general-boolean refusal, named for the compiler's caller.
pub(super) fn remap_boolean_error(error: BooleanError) -> GeomError {
    let input = match error {
        BooleanError::EmptyResult => {
            return GeomError::Degenerate(
                "exact boolean difference removes the whole subject".to_owned(),
            )
        }
        BooleanError::UnsupportedSection => {
            "exact boolean whose operands meet in a section curve the general boolean does not build"
        }
        BooleanError::UnsupportedTrim => {
            "exact boolean over a face trim the general boolean cannot intersect"
        }
        BooleanError::Undecided => {
            "exact boolean with a point too close to a face boundary to classify"
        }
        BooleanError::UnsupportedSplit => {
            "exact boolean over a face or section without an exact pcurve"
        }
        BooleanError::UnclosedSplit => "exact boolean whose split face pieces do not close",
        BooleanError::TangentSplit => "exact boolean with tangent pieces leaving one vertex",
        BooleanError::AmbiguousCavity => "exact boolean leaving a cavity inside no solid",
        BooleanError::Assembly => "exact boolean whose kept faces do not sew into a solid",
        BooleanError::NearCoincidence => {
            "exact boolean with features chained within tolerance over more than it"
        }
        BooleanError::UnsupportedContact => {
            "exact boolean whose plane touching a cylinder within tolerance is crossed where the contact cannot be placed"
        }
        BooleanError::Evaluation => "exact boolean over a curve or surface it cannot evaluate",
        BooleanError::DanglingReference => "exact boolean over a dangling operand handle",
        BooleanError::Measure(_) => "exact boolean over a face whose domain cannot be built",
        _ => "exact boolean configuration the general boolean refuses",
    };
    unsupported(input)
}
