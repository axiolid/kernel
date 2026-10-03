//! Plan distance, clearance and overlap between bodies of several exact
//! solids (#237): [`crate::plan_boundary_distance`] and its kin (#217) for
//! [`PlacedBody`] sides.
//!
//! # What is measured
//!
//! A body's shadow on the XY plane is the union of its items' shadows, and
//! each item's shadow is its boundary's (a vertical line through a point of
//! a solid meets its boundary). So the plan distance between two bodies is
//! the least horizontal distance between a boundary point of an item of one
//! and a boundary point of an item of the other: the minimum over item
//! pairs, found by ONE branch and bound over the elements of every item, as
//! [`super::body_boundary_distance`] in space. A pair of items far apart in
//! plan gets a large lower bound at the root and is never refined, and the
//! witnesses name the items they lie on.
//!
//! As for one pair, the distance is zero where the shadows overlap and
//! where one item's shadow lies inside another's: every boundary point of
//! the inner item then has a boundary point of the outer one directly above
//! or below it. Over several items this holds item pair by item pair, so a
//! body standing on another's footprint measures zero whichever of its
//! items does.
//!
//! # Items of one body
//!
//! How the items of one body lie against each other does not matter in
//! plan: their shadows may overlap freely (a column over its footing), and
//! the union of the shadows is the union of the boundaries' shadows however
//! the items touch, overlap or nest. So no layout is checked and nothing is
//! cut, unlike the Hausdorff queries in [`super`]; the items are only
//! merged and placed.
//!
//! # Placement
//!
//! The placement is applied before projecting: projection is along the
//! world's `z`, so a body turned about `z` measures as its plan turned, and
//! a tilted body casts the shadow of its tilted items. Witnesses are in the
//! world.

use axiolid_core::{Point2, Scalar, Tolerance};

use super::{distance, BodyDistance, BodyMeasureError, Items, PlacedBody};
use crate::exact_distance::{overlap_search, Clearance, Metric, PlanOverlap};

/// Plan distance between two bodies of exact solids -- between their
/// shadows on the XY plane -- to within `accuracy` (#237): the least over
/// item pairs, by one shared search.
///
/// Zero where an item's shadow overlaps or lies inside one of the other
/// body's, as [`crate::plan_boundary_distance`] for one pair; items of one
/// body may overlap in plan, or in space. When two planar faces show an
/// overlap of positive area, the interval is `[0, 0]` and the witnesses
/// project into it. `point_a` and `point_b` lie on the items `item_a` and
/// `item_b` name, their projections `upper` apart. A negative or NaN
/// accuracy is read as zero.
///
/// # Errors
///
/// [`BodyMeasureError::EmptyBody`] for an empty body,
/// [`BodyMeasureError::Placement`] for a placement that is not rigid, and
/// [`BodyMeasureError::Measure`] with what
/// [`crate::plan_boundary_distance`] refuses.
pub fn body_plan_boundary_distance<'a, 'b>(
    a: impl Into<PlacedBody<'a>>,
    b: impl Into<PlacedBody<'b>>,
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<BodyDistance, BodyMeasureError> {
    let accuracy = accuracy.max(0.0);
    distance(
        a.into(),
        b.into(),
        tolerance,
        Metric::Plan,
        &mut |lower, upper| upper - lower <= accuracy,
    )
}

/// [`body_plan_boundary_distance`] refined only until it clears `limit`.
///
/// # Errors
///
/// As [`body_plan_boundary_distance`].
pub fn body_plan_boundary_clearance<'a, 'b>(
    a: impl Into<PlacedBody<'a>>,
    b: impl Into<PlacedBody<'b>>,
    limit: Scalar,
    tolerance: Tolerance,
) -> Result<(BodyDistance, Clearance), BodyMeasureError> {
    let found = distance(
        a.into(),
        b.into(),
        tolerance,
        Metric::Plan,
        &mut |lower, upper| upper < limit || lower > limit,
    )?;
    let clearance = found.bounds.against(limit);
    Ok((found, clearance))
}

/// Whether the plan projections of two bodies overlap, and which items
/// show it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BodyPlanOverlap {
    /// Over a patch of positive area: `at` and a disc around it lie in the
    /// shadows of item `item_a` of the first body and item `item_b` of the
    /// second, shown by a planar face of each that is not vertical.
    Overlapping {
        /// A plan point inside both items' shadows.
        at: Point2,
        /// The item of the first body, by index in its slice.
        item_a: usize,
        /// The item of the second body, by index in its slice.
        item_b: usize,
    },
    /// Certainly apart: no item of one comes within `gap`, positive, of an
    /// item of the other in plan.
    Disjoint {
        /// A certified lower bound on the plan distance between the bodies.
        gap: Scalar,
    },
    /// Neither shown: the shadows may only touch, overlap where no planar
    /// face covers both, or lie closer than the search resolves.
    Undecided,
}

/// Whether the plan projections of two bodies overlap over a patch of
/// positive area -- an item of one over an item of the other -- are
/// certainly apart, or neither could be shown (#237).
///
/// As [`crate::plan_overlap`] for one pair, over every item pair by one
/// search: overlap is shown by two planar faces, one of an item of each
/// body, not vertical, whose shadows share an open patch. Items of one body
/// overlapping each other in plan are not an overlap between the bodies.
///
/// # Errors
///
/// As [`body_plan_boundary_distance`].
pub fn body_plan_overlap<'a, 'b>(
    a: impl Into<PlacedBody<'a>>,
    b: impl Into<PlacedBody<'b>>,
    tolerance: Tolerance,
) -> Result<BodyPlanOverlap, BodyMeasureError> {
    let (a, b) = Items::pair(a.into(), b.into())?;
    let (verdict, on) = overlap_search(&a.brep, &b.brep, tolerance)?;
    Ok(match (verdict, on) {
        (PlanOverlap::Overlapping { at }, Some((on_a, on_b))) => BodyPlanOverlap::Overlapping {
            at,
            item_a: a.item(on_a),
            item_b: b.item(on_b),
        },
        (PlanOverlap::Disjoint { gap }, _) => BodyPlanOverlap::Disjoint { gap },
        _ => BodyPlanOverlap::Undecided,
    })
}
