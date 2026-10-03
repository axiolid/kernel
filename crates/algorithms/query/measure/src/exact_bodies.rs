//! Certified boundary distance and Hausdorff distance between bodies made
//! of several exact solids (#229).
//!
//! A body here is a slice of [`ExactBRep`] items, one per solid: how an IFC
//! body with several items arrives. A slice, not one B-rep holding several
//! solids, because the items are produced and placed one by one, a
//! witness has to say which item it lies on, and the layout check below
//! works item against item. (One `ExactBRep` with several solids is
//! measured by the single-B-rep queries as the union of all its faces,
//! with none of the checks made here.)
//!
//! # Distance
//!
//! [`body_boundary_distance`] is the least distance between a point on the
//! boundary of an item of `a` and a point on the boundary of an item of
//! `b`: the minimum over item pairs, found by ONE branch and bound over the
//! elements of every item (as [`crate::boundary_distance`] for one pair),
//! so a pair of items far apart gets a large lower bound at the root and is
//! never refined. The witnesses name the items they lie on.
//!
//! Items may touch or overlap. The minimum over pairs is then still the
//! distance between the boundaries of the two unions, except where part of
//! one body lies wholly inside the other: a boundary point of an item that
//! is not on the union's boundary is interior to the union, and from a
//! point outside the union the segment to it crosses the union's boundary
//! at a point no farther away. As for one pair, containment is a separate
//! classification.
//!
//! # In plan
//!
//! [`body_plan_boundary_distance`], [`body_plan_boundary_clearance`] and
//! [`body_plan_overlap`] measure between the bodies' shadows on the XY
//! plane by the same shared search (#237, `plan.rs`). Items of one body may
//! overlap freely there, so nothing about their layout is checked.
//!
//! # Hausdorff
//!
//! [`body_boundary_hausdorff_distance`] and
//! [`one_sided_body_boundary_hausdorff`] measure the boundary of the UNION
//! of each body's items. Every body is first checked, item pair by item
//! pair, for a layout under which that boundary can be formed exactly:
//!
//! - **Apart**: shown by a plane with a certified positive gap between
//!   them -- a direction along which the exact projection ranges of every
//!   face patch and edge span of the two are disjoint, tried along the axes
//!   and every planar face's normal -- or, failing that, by a positive
//!   certified lower bound on the distance between their boundaries
//!   together with a boundary point of each outside an enclosing box of the
//!   other. The second test needs each item to be one solid (one outer
//!   shell, voids allowed): two such solids whose boundaries do not meet
//!   are disjoint unless the outer shell of one lies inside the other, and
//!   then its box lies in the other's box. However small the gap -- an
//!   exporter's 0.3 mm between parts -- both facing faces are on the
//!   union's boundary and are measured as such: a gap is never glued.
//! - **Touching without a shared patch**: a plane separates them to within
//!   rounding (`1e-9` relative to their size), and one of the two has no
//!   face that can hold a patch of that plane -- no planar face lying in it
//!   and no B-spline face reaching it (no other family holds an open patch
//!   of a plane). Interiors are then disjoint, and a boundary point of one
//!   that is not on the union's boundary would have a ball round it covered
//!   by the two, cut by the plane into halves, one in each: a disc of the
//!   plane on both boundaries, which one of them cannot hold. So contact
//!   along edges or at points leaves the union boundary the union of the
//!   boundaries. Planes are tried along the axes, every planar face's
//!   normal, and the sum of any two touching ones (two blocks meeting at an
//!   edge are separated by the diagonal plane, where neither has a face).
//! - **Exact face contact**: they touch on a plane where both have faces,
//!   every such face planar and lying, in the B-rep's own numbers, on one
//!   plane normal to a coordinate axis (its plane, vertices, lines and
//!   circles all at the very same coordinate along it; a placement turning
//!   about that axis keeps this), and bounded by lines and circles. Each
//!   face in contact is then cut down to its free region, what is left of
//!   it once the faces of other items lying against it are cut away, by an
//!   exact arrangement of their boundaries on the plane (`contact`).
//!
//! Why cutting is right: interiors are disjoint, so a point `p` of an
//! item's face `F` is off the union's boundary only if a ball round it is
//! covered by items. By the argument above, the half of the ball across
//! `F`'s plane is covered by other items whose faces then hold a disc of
//! that plane round `p` -- faces lying against `F`. So the union's boundary
//! is the closure of the free regions of all faces: those in no contact as
//! they are, and those in contact cut. Its edges are the edges of those
//! regions; an edge between two faces both cut away (the foot of a wall
//! shared by two blocks on a footing) is not on it and is dropped.
//!
//! Anything else is refused by name, never glued or snapped:
//! [`BodyMeasureError::ContactPlaneNotAxisNormal`] when two items touch, to
//! within rounding, on a plane no coordinate axis is normal to (a tilted
//! assembly, walls turned in plan: whether they meet, clear or overlap is
//! below rounding there); [`BodyMeasureError::ItemsNearlyShareFace`] with
//! the gap when they interpenetrate by no more than the caller's tolerance,
//! or their faces are on axis planes at coordinates that differ below
//! rounding; [`BodyMeasureError::ItemsShareFace`]
//! when such faces cannot be cut (a B-spline face, an elliptical edge);
//! and [`BodyMeasureError::ItemsOverlap`] when they were not shown apart or
//! touching: boundaries that cross or come closer than the search resolves
//! with no plane between them, one item possibly inside another, or curved
//! items touching (no planar face to separate them by).
//!
//! Gluing near contact and widening the result by the gap is not sound:
//! lifted 0.3 mm off its footing, a column's base disc and the disc of the
//! footing's top under it are both on the union's boundary, and their
//! centres are a whole radius from the boundary of the body glued shut. So
//! a gap is measured as a gap, and an interpenetration is refused.
//!
//! The one-sided distance is then measured by the search of
//! [`crate::one_sided_boundary_hausdorff`] on the union's boundary of each
//! body. Faces are matched across every pair of items, so a body moved by
//! a translation closes as fast as one solid does: every face patch is
//! held at `|t|` by its own translate, and the support point of each item
//! against `t` seeds the lower bound. A free region is bounded through the
//! face it was cut from, so an exact cut that falls out differently in the
//! last bit for a translated copy (a sliver along a wall turned in plan)
//! still closes at once (`exact_hausdorff/cut.rs`).
//!
//! # Placed bodies
//!
//! Every query takes each side as a [`PlacedBody`]: the items in the
//! body's own frame and one rigid placement. An IFC body's extruded items
//! stand on axis-normal planes in the body's frame however the body is
//! placed in the world, so contact is found and cut there, and the
//! placement then moves the cut boundary (and the items it was cut from)
//! before anything is measured: witnesses are in the world. A placement
//! cannot change what is measured, a rigid motion carrying the union's
//! boundary onto the placed union's. Items turned against each other in
//! the body's own frame still touch on planes no axis is normal to and are
//! refused as [`BodyMeasureError::ContactPlaneNotAxisNormal`]; a slice of
//! items converts to a placed body with the identity placement.

use std::borrow::Cow;
use std::collections::HashMap;

use core::fmt;

use axiolid_brep::{ExactBRep, ExactBRepBuilder, TransformError};
use axiolid_core::{Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_evaluate::evaluate3;
use axiolid_surface::Surface;

use crate::exact::ExactMeasureError;
use crate::exact_distance::{
    search, search_within, surface_of, Clearance, DistanceBounds, Metric, Shape, Side,
    OVERLAP_STEPS,
};
use crate::exact_hausdorff::{witnessed, Cut, MAX_SPLITS};
use crate::mesh_hausdorff::HausdorffBounds;

mod contact;
mod plan;

use contact::{Assembler, AxisPlane, Uncut};
pub use plan::{
    body_plan_boundary_clearance, body_plan_boundary_distance, body_plan_overlap, BodyPlanOverlap,
};

/// Which argument of a body query a refusal is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodySide {
    /// The first slice: `a`, or `from` of a one-sided query.
    First,
    /// The second slice: `b`, or `to` of a one-sided query.
    Second,
}

/// A body: its items in the body's own frame, and one rigid placement
/// that puts them in the world (#229).
///
/// Every `body_*` query takes one per side, and a slice of items (or a
/// vector or array of them) converts into one with the identity placement.
/// Face contact is found and cut in the body's own frame, where an IFC
/// body's extruded items usually stand on axis-normal planes, and the
/// placement then moves the cut boundary, so a body turned in any
/// direction still has its contact cut. Witnesses are in the world.
///
/// Items turned against each other within the body's own frame still touch
/// on planes no axis is normal to, and are refused by
/// [`BodyMeasureError::ContactPlaneNotAxisNormal`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlacedBody<'a> {
    /// The items, in the body's own frame.
    pub items: &'a [ExactBRep],
    /// The rigid placement of the body's frame in the world.
    pub placement: Transform3,
}

impl<'a> PlacedBody<'a> {
    /// `items` placed by `placement`.
    #[must_use]
    pub fn new(items: &'a [ExactBRep], placement: Transform3) -> Self {
        Self { items, placement }
    }
}

impl<'a> From<&'a [ExactBRep]> for PlacedBody<'a> {
    fn from(items: &'a [ExactBRep]) -> Self {
        Self::new(items, Transform3::IDENTITY)
    }
}

impl<'a> From<&'a Vec<ExactBRep>> for PlacedBody<'a> {
    fn from(items: &'a Vec<ExactBRep>) -> Self {
        Self::new(items, Transform3::IDENTITY)
    }
}

impl<'a, const N: usize> From<&'a [ExactBRep; N]> for PlacedBody<'a> {
    fn from(items: &'a [ExactBRep; N]) -> Self {
        Self::new(items, Transform3::IDENTITY)
    }
}

/// Why two bodies of exact solids could not be measured.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum BodyMeasureError {
    /// The body has no items.
    EmptyBody {
        /// Which body.
        body: BodySide,
    },
    /// An item could not be measured, as by the single-solid queries.
    Measure(ExactMeasureError),
    /// Two items of one body were not shown apart or touching (see the
    /// module docs): their boundaries may cross, one may lie inside the
    /// other, or they touch where no planar face separates them. The
    /// boundary of their union is not formed here.
    ItemsOverlap {
        /// Which body.
        body: BodySide,
        /// The first item, by index in the slice.
        first: usize,
        /// The second item, by index in the slice.
        second: usize,
    },
    /// Two items of one body touch where both have a face in the plane
    /// that separates them, and the shared patch cannot be cut exactly: a
    /// B-spline face reaches the plane, or a face in contact is bounded by
    /// an edge that is neither a line nor a circle.
    ItemsShareFace {
        /// Which body.
        body: BodySide,
        /// The first item, by index in the slice.
        first: usize,
        /// The second item, by index in the slice.
        second: usize,
    },
    /// Two items of one body are within the caller's tolerance of face
    /// contact without being in exact contact: they interpenetrate by up to
    /// `-gap`, or their faces are on planes normal to one axis but not at
    /// the same coordinate in their own numbers (a gap or an overlap below
    /// rounding). Gluing them would measure a boundary the items do not
    /// have, so the pair is refused rather than snapped (a certified
    /// positive gap is measured exactly instead: the two faces are on the
    /// union's boundary).
    ItemsNearlyShareFace {
        /// Which body.
        body: BodySide,
        /// The first item, by index in the slice.
        first: usize,
        /// The second item, by index in the slice.
        second: usize,
        /// The certified least separation along the plane between them:
        /// negative where they may interpenetrate, by at most its size.
        gap: Scalar,
    },
    /// Two items of one body touch, to within rounding, where both have
    /// faces on a plane that no coordinate axis is normal to: walls of
    /// items turned against each other, or an assembly placed before it was
    /// measured. Contact is cut only on axis-normal planes, where it can be
    /// shown exact; pass the items in the body's own frame with the
    /// placement in a [`PlacedBody`] instead, so the contact is cut before
    /// the placement turns it.
    ContactPlaneNotAxisNormal {
        /// Which body.
        body: BodySide,
        /// The first item, by index in the slice.
        first: usize,
        /// The second item, by index in the slice.
        second: usize,
    },
    /// The placement of a [`PlacedBody`] is not rigid, or its items cannot
    /// be moved by it.
    Placement {
        /// Which body.
        body: BodySide,
        /// Why the items could not be placed.
        error: TransformError,
    },
}

impl fmt::Display for BodyMeasureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBody { body } => write!(f, "the {body:?} body has no items"),
            Self::Measure(error) => write!(f, "an item could not be measured: {error}"),
            Self::ItemsOverlap {
                body,
                first,
                second,
            } => write!(
                f,
                "items {first} and {second} of the {body:?} body were not shown apart \
                 or touching; the boundary of their union is not formed"
            ),
            Self::ItemsShareFace {
                body,
                first,
                second,
            } => write!(
                f,
                "items {first} and {second} of the {body:?} body may share a patch of \
                 face that cannot be cut exactly"
            ),
            Self::ItemsNearlyShareFace {
                body,
                first,
                second,
                gap,
            } => write!(
                f,
                "items {first} and {second} of the {body:?} body are within the \
                 tolerance of face contact (gap {gap}) but not in exact contact"
            ),
            Self::ContactPlaneNotAxisNormal {
                body,
                first,
                second,
            } => write!(
                f,
                "items {first} and {second} of the {body:?} body touch on a plane no \
                 coordinate axis is normal to; pass them unplaced in a PlacedBody"
            ),
            Self::Placement { body, error } => {
                write!(f, "the {body:?} body could not be placed: {error}")
            }
        }
    }
}

impl std::error::Error for BodyMeasureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Measure(error) => Some(error),
            Self::Placement { error, .. } => Some(error),
            _ => None,
        }
    }
}

impl From<ExactMeasureError> for BodyMeasureError {
    fn from(error: ExactMeasureError) -> Self {
        Self::Measure(error)
    }
}

/// A certified distance between two bodies and the items its witnesses lie
/// on.
#[derive(Debug, Clone, PartialEq)]
pub struct BodyDistance {
    /// The interval and the witnesses, as for one pair.
    pub bounds: DistanceBounds,
    /// The item of the first body holding `bounds.point_a`.
    pub item_a: usize,
    /// The item of the second body holding `bounds.point_b`.
    pub item_b: usize,
}

/// A certified one-sided Hausdorff distance between two bodies and the
/// items its witnesses lie on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyHausdorffBounds {
    /// The interval and the witnesses, as for one pair.
    pub bounds: HausdorffBounds,
    /// The item of the measured body holding `bounds.point_from`.
    pub item_from: usize,
    /// The item of the other body holding `bounds.point_to`.
    pub item_to: usize,
}

/// The two-sided Hausdorff distance between the boundaries of two bodies
/// and both one-sided ones.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyHausdorff {
    /// `max(h(A, B), h(B, A))`, witnessed by the side with the larger lower
    /// bound.
    pub distance: BodyHausdorffBounds,
    /// `h(A, B)`, from the first body's items to the second's.
    pub forward: BodyHausdorffBounds,
    /// `h(B, A)`, from the second body's items to the first's.
    pub backward: BodyHausdorffBounds,
}

/// Distance between the boundaries of two bodies of exact solids, to within
/// `accuracy` (#229): the least over item pairs, by one shared search.
///
/// Items may touch or overlap (see the module docs). A negative or NaN
/// accuracy is read as zero.
///
/// # Errors
///
/// [`BodyMeasureError::EmptyBody`] for an empty slice, and
/// [`BodyMeasureError::Measure`] with what [`crate::boundary_distance`]
/// refuses.
pub fn body_boundary_distance<'a, 'b>(
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
        Metric::Space,
        &mut |lower, upper| upper - lower <= accuracy,
    )
}

/// [`body_boundary_distance`] refined only until it clears `limit`.
///
/// # Errors
///
/// As [`body_boundary_distance`].
pub fn body_boundary_clearance<'a, 'b>(
    a: impl Into<PlacedBody<'a>>,
    b: impl Into<PlacedBody<'b>>,
    limit: Scalar,
    tolerance: Tolerance,
) -> Result<(BodyDistance, Clearance), BodyMeasureError> {
    let found = distance(
        a.into(),
        b.into(),
        tolerance,
        Metric::Space,
        &mut |lower, upper| upper < limit || lower > limit,
    )?;
    let clearance = found.bounds.against(limit);
    Ok((found, clearance))
}

/// The least distance by `metric` over item pairs, by one search over the
/// placed items merged into one B-rep per body.
fn distance(
    a: PlacedBody<'_>,
    b: PlacedBody<'_>,
    tolerance: Tolerance,
    metric: Metric,
    done: &mut dyn FnMut(Scalar, Scalar) -> bool,
) -> Result<BodyDistance, BodyMeasureError> {
    let (a, b) = Items::pair(a, b)?;
    let found = search(&a.brep, &b.brep, tolerance, metric, done)?;
    Ok(BodyDistance {
        item_a: a.item(found.on.0),
        item_b: b.item(found.on.1),
        bounds: found.bounds,
    })
}

/// One-sided Hausdorff distance from the boundary of the union of `from`'s
/// items to that of `to`'s, to within `accuracy` (#229).
///
/// # Errors
///
/// [`BodyMeasureError::EmptyBody`] for an empty slice;
/// [`BodyMeasureError::ItemsOverlap`] or
/// [`BodyMeasureError::ItemsShareFace`] for two items of one body not shown
/// apart or touching without a shared patch (see the module docs); and
/// [`BodyMeasureError::Measure`] with what
/// [`crate::one_sided_boundary_hausdorff`] refuses.
pub fn one_sided_body_boundary_hausdorff<'a, 'b>(
    from: impl Into<PlacedBody<'a>>,
    to: impl Into<PlacedBody<'b>>,
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<BodyHausdorffBounds, BodyMeasureError> {
    one_sided_body_boundary_hausdorff_with_budget(from, to, accuracy, tolerance, MAX_SPLITS)
}

/// [`one_sided_body_boundary_hausdorff`] with the caller's cap on the
/// splits, as [`crate::one_sided_boundary_hausdorff_with_budget`].
///
/// # Errors
///
/// As [`one_sided_body_boundary_hausdorff`].
pub fn one_sided_body_boundary_hausdorff_with_budget<'a, 'b>(
    from: impl Into<PlacedBody<'a>>,
    to: impl Into<PlacedBody<'b>>,
    accuracy: Scalar,
    tolerance: Tolerance,
    max_splits: usize,
) -> Result<BodyHausdorffBounds, BodyMeasureError> {
    let source = Items::boundary(from.into(), BodySide::First, tolerance)?;
    let target = Items::boundary(to.into(), BodySide::Second, tolerance)?;
    one_sided(&source, &target, accuracy, tolerance, max_splits)
}

/// Two-sided Hausdorff distance between the boundaries of the unions of two
/// bodies' items, to within `accuracy` (#229).
///
/// # Errors
///
/// As [`one_sided_body_boundary_hausdorff`].
pub fn body_boundary_hausdorff_distance<'a, 'b>(
    a: impl Into<PlacedBody<'a>>,
    b: impl Into<PlacedBody<'b>>,
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<BodyHausdorff, BodyMeasureError> {
    let first = Items::boundary(a.into(), BodySide::First, tolerance)?;
    let second = Items::boundary(b.into(), BodySide::Second, tolerance)?;
    let forward = one_sided(&first, &second, accuracy, tolerance, MAX_SPLITS)?;
    let backward = one_sided(&second, &first, accuracy, tolerance, MAX_SPLITS)?;
    let mut distance = if backward.bounds.lower > forward.bounds.lower {
        backward
    } else {
        forward
    };
    distance.bounds.upper = forward.bounds.upper.max(backward.bounds.upper);
    Ok(BodyHausdorff {
        distance,
        forward,
        backward,
    })
}

fn one_sided(
    from: &Items<'_>,
    to: &Items<'_>,
    accuracy: Scalar,
    tolerance: Tolerance,
    max_splits: usize,
) -> Result<BodyHausdorffBounds, BodyMeasureError> {
    let cuts = match (&from.cut, &to.cut) {
        (None, None) => None,
        _ => Some((from.cut(), to.cut())),
    };
    let found = witnessed(
        &from.brep,
        &to.brep,
        accuracy,
        tolerance,
        max_splits,
        &from.edge_ranges(),
        cuts.as_ref().map(|(a, b)| (a, b)),
    )?;
    Ok(BodyHausdorffBounds {
        bounds: found.bounds,
        item_from: from.item(found.from),
        item_to: to.item(found.to),
    })
}

/// A body's items as one B-rep holding all their faces and edges, with the
/// first face and edge of each item to name it by.
struct Items<'a> {
    brep: Cow<'a, ExactBRep>,
    faces: Vec<usize>,
    edges: Vec<usize>,
    /// Where `brep` is the items' boundary cut down to free regions: the
    /// items uncut, and how each cut face and each face in contact relate.
    cut: Option<CutFrom<'a>>,
}

/// The items uncut, the face of them each cut face of the boundary came
/// from, and the faces lying against each face in contact, all by index
/// into the uncut items.
struct CutFrom<'a> {
    items: Cow<'a, ExactBRep>,
    origin: HashMap<usize, usize>,
    partners: HashMap<usize, Vec<usize>>,
}

impl<'a> Items<'a> {
    fn new(items: &'a [ExactBRep], body: BodySide) -> Result<Self, BodyMeasureError> {
        match items {
            [] => Err(BodyMeasureError::EmptyBody { body }),
            [one] => Ok(Self {
                brep: Cow::Borrowed(one),
                faces: vec![0],
                edges: vec![0],
                cut: None,
            }),
            _ => {
                let mut builder = ExactBRepBuilder::default();
                let (mut faces, mut edges) = (Vec::new(), Vec::new());
                for item in items {
                    faces.push(builder.topology_mut().faces().len());
                    edges.push(builder.topology_mut().edges().len());
                    builder.append(item, false);
                }
                // Every item is a valid B-rep, and so is their disjoint
                // union; a refusal here is a catalog the items lost.
                let brep = builder
                    .finish()
                    .map_err(|_| ExactMeasureError::DanglingReference)?;
                Ok(Self {
                    brep: Cow::Owned(brep),
                    faces,
                    edges,
                    cut: None,
                })
            }
        }
    }

    /// Both bodies' items, uncut, merged and placed.
    fn pair(a: PlacedBody<'a>, b: PlacedBody<'a>) -> Result<(Self, Self), BodyMeasureError> {
        let a = Self::new(a.items, BodySide::First)?.placed(a.placement, BodySide::First)?;
        let b = Self::new(b.items, BodySide::Second)?.placed(b.placement, BodySide::Second)?;
        Ok((a, b))
    }

    /// How the measured boundary was cut from the items: nothing cut when
    /// no faces were in contact.
    fn cut(&self) -> Cut<'_> {
        match &self.cut {
            Some(cut) => Cut {
                items: &cut.items,
                origin: cut.origin.clone(),
                partners: cut.partners.clone(),
            },
            None => Cut {
                items: &self.brep,
                origin: HashMap::new(),
                partners: HashMap::new(),
            },
        }
    }

    /// The boundary of the union of a body's items: their faces, with those
    /// in exact contact cut down to their free regions, once every pair of
    /// items is shown apart, touching or in exact contact.
    fn boundary(
        placed: PlacedBody<'a>,
        body: BodySide,
        tolerance: Tolerance,
    ) -> Result<Self, BodyMeasureError> {
        Self::cut_boundary(placed.items, body, tolerance)?.placed(placed.placement, body)
    }

    /// Move the measured boundary (and the items it was cut from) by a
    /// rigid placement; face and edge indices are kept.
    fn placed(self, placement: Transform3, body: BodySide) -> Result<Self, BodyMeasureError> {
        if placement == Transform3::IDENTITY {
            return Ok(self);
        }
        let place = |brep: &ExactBRep| {
            brep.transformed(&placement)
                .map_err(|error| BodyMeasureError::Placement { body, error })
        };
        let cut = match self.cut {
            Some(cut) => Some(CutFrom {
                items: Cow::Owned(place(&cut.items)?),
                origin: cut.origin,
                partners: cut.partners,
            }),
            None => None,
        };
        Ok(Self {
            brep: Cow::Owned(place(&self.brep)?),
            faces: self.faces,
            edges: self.edges,
            cut,
        })
    }

    /// The boundary of the union of the items, in their own frame.
    fn cut_boundary(
        items: &'a [ExactBRep],
        body: BodySide,
        tolerance: Tolerance,
    ) -> Result<Self, BodyMeasureError> {
        if items.is_empty() {
            return Err(BodyMeasureError::EmptyBody { body });
        }
        let partners = check_layout(items, body, tolerance)?;
        if partners.is_empty() {
            return Self::new(items, body);
        }
        let uncut = Self::new(items, body)?;
        let mut origin = HashMap::new();
        let mut assembler = Assembler::default();
        let (mut faces, mut edges) = (Vec::new(), Vec::new());
        for (index, item) in items.iter().enumerate() {
            let (first_face, first_edge) = assembler.counts();
            faces.push(first_face);
            edges.push(first_edge);
            assembler.next_item();
            for face in 0..item.topology().faces().len() {
                let cut = match partners.get(&(index, face)) {
                    Some((plane, against)) => {
                        let against: Vec<(&ExactBRep, usize)> =
                            against.iter().map(|&(j, g)| (&items[j], g)).collect();
                        contact::free_regions(item, face, &against, *plane)
                            .map_err(|()| BodyMeasureError::ItemsShareFace {
                                body,
                                first: index.min(partners[&(index, face)].1[0].0),
                                second: index.max(partners[&(index, face)].1[0].0),
                            })?
                            .map(|regions| (regions, *plane))
                    }
                    None => None,
                };
                match cut {
                    Some((regions, plane)) => {
                        for region in &regions {
                            origin.insert(assembler.counts().0, uncut.faces[index] + face);
                            assembler.add_region(region, plane);
                        }
                    }
                    None => assembler
                        .copy_face(item, face)
                        .ok_or(ExactMeasureError::DanglingReference)?,
                }
            }
        }
        let brep = assembler
            .finish()
            .ok_or(ExactMeasureError::DanglingReference)?;
        let partners = partners
            .iter()
            .map(|(&(item, face), (_, against))| {
                (
                    uncut.faces[item] + face,
                    against
                        .iter()
                        .map(|&(other, g)| uncut.faces[other] + g)
                        .collect(),
                )
            })
            .collect();
        Ok(Self {
            brep: Cow::Owned(brep),
            faces,
            edges,
            cut: Some(CutFrom {
                items: uncut.brep,
                origin,
                partners,
            }),
        })
    }

    /// Each item's edges, by index in the merged B-rep.
    fn edge_ranges(&self) -> Vec<core::ops::Range<usize>> {
        let total = self.brep.topology().edges().len();
        (0..self.edges.len())
            .map(|k| self.edges[k]..self.edges.get(k + 1).copied().unwrap_or(total))
            .collect()
    }

    /// The item an element of the merged B-rep belongs to.
    fn item(&self, shape: Shape) -> usize {
        let (starts, index) = match shape {
            Shape::Face { face, .. } => (&self.faces, face),
            Shape::Edge { edge, .. } => (&self.edges, edge),
        };
        starts
            .partition_point(|&start| start <= index)
            .saturating_sub(1)
    }
}

/// How two items of one body lie (see the module docs).
#[derive(Debug, Clone, PartialEq)]
enum Layout {
    /// Shown disjoint.
    Apart,
    /// Shown touching, with no patch of face shared.
    Touching,
    /// In exact face contact on `plane`: these faces of each lie on it.
    Contact {
        plane: AxisPlane,
        faces_a: Vec<usize>,
        faces_b: Vec<usize>,
    },
    /// Touching where both may hold a patch of the separating plane, on a
    /// face that is not planar.
    MayShareFace,
    /// Within `gap` of face contact (negative: interpenetrating), on axis
    /// planes that differ, or across a plane more than rounding deep.
    NearlyShareFace { gap: Scalar },
    /// Touching to within rounding where both have faces on a plane that
    /// no coordinate axis is normal to.
    NotAxisNormal,
    /// Neither shown.
    Undecided,
}

/// The faces of other items lying against each face in exact contact, by
/// item and face, with the plane they share.
type Partners = HashMap<(usize, usize), (AxisPlane, Vec<(usize, usize)>)>;

/// Refuse a body two of whose items are not shown apart, touching without
/// a shared patch, or in exact face contact; and list the faces in contact.
fn check_layout(
    items: &[ExactBRep],
    body: BodySide,
    tolerance: Tolerance,
) -> Result<Partners, BodyMeasureError> {
    let mut partners = Partners::new();
    if items.len() < 2 {
        return Ok(partners);
    }
    let linear = tolerance.linear().max(1e-12);
    let sides = items
        .iter()
        .map(|item| Side::new(item, linear, Metric::Space))
        .collect::<Result<Vec<_>, _>>()?;
    for first in 0..sides.len() {
        for second in first + 1..sides.len() {
            match layout(&sides[first], &sides[second], tolerance)? {
                Layout::Apart | Layout::Touching => {}
                Layout::Contact {
                    plane,
                    faces_a,
                    faces_b,
                } => {
                    for &fa in &faces_a {
                        for &fb in &faces_b {
                            for (from, to) in
                                [((first, fa), (second, fb)), ((second, fb), (first, fa))]
                            {
                                partners
                                    .entry(from)
                                    .or_insert_with(|| (plane, Vec::new()))
                                    .1
                                    .push(to);
                            }
                        }
                    }
                }
                Layout::NearlyShareFace { gap } => {
                    return Err(BodyMeasureError::ItemsNearlyShareFace {
                        body,
                        first,
                        second,
                        gap,
                    })
                }
                Layout::NotAxisNormal => {
                    return Err(BodyMeasureError::ContactPlaneNotAxisNormal {
                        body,
                        first,
                        second,
                    })
                }
                Layout::MayShareFace => {
                    return Err(BodyMeasureError::ItemsShareFace {
                        body,
                        first,
                        second,
                    })
                }
                Layout::Undecided => {
                    return Err(BodyMeasureError::ItemsOverlap {
                        body,
                        first,
                        second,
                    })
                }
            }
        }
    }
    Ok(partners)
}

/// The range of `d . x` over every element of a side: a sound enclosure of
/// the item's extent along `d`, which is its boundary's.
fn extent(side: &Side<'_>, d: Vec3) -> Result<(Scalar, Scalar), ExactMeasureError> {
    let mut range = (Scalar::INFINITY, Scalar::NEG_INFINITY);
    for element in &side.elements {
        let (lo, hi) = side.project(element, d)?;
        range = (range.0.min(lo), range.1.max(hi));
    }
    if range.0.is_finite() && range.1.is_finite() {
        Ok(range)
    } else {
        Err(crate::exact::EVALUATION)
    }
}

/// The unit normals of a side's planar faces.
fn plane_normals(side: &Side<'_>) -> Result<Vec<Vec3>, ExactMeasureError> {
    let mut normals = Vec::new();
    for face in side.brep.topology().faces() {
        if let Surface::Plane(plane) = surface_of(side.brep, face.surface)? {
            let n = plane.frame.x.cross(plane.frame.y);
            let length = n.length();
            if length.is_finite() && length > 0.0 {
                normals.push(n / length);
            }
        }
    }
    Ok(normals)
}

/// A plane `d . x = level` with `a` on its low side and `b` on its high
/// side, to within `slack`; `gap` is the certified least separation of
/// the two along `d`, negative where their ranges overlap.
#[derive(Debug, Clone, Copy)]
struct Plane {
    d: Vec3,
    level: Scalar,
    slack: Scalar,
    gap: Scalar,
}

/// How the items' ranges along the unit `d` lie: certainly apart (a
/// positive gap), within rounding of a plane between them, overlapping by
/// no more than `tolerance`, or overlapping more.
fn separate(
    a: &Side<'_>,
    b: &Side<'_>,
    d: Vec3,
    tolerance: Scalar,
) -> Result<Separation, ExactMeasureError> {
    let (a_lo, a_hi) = extent(a, d)?;
    let (b_lo, b_hi) = extent(b, d)?;
    let scale = a_lo.abs().max(a_hi.abs()).max(b_lo.abs()).max(b_hi.abs());
    let slack = 1e-9 * (1.0 + scale);
    let mut found = Separation::Overlapping;
    for (gap, d, level) in [
        (b_lo - a_hi, d, 0.5 * (a_hi + b_lo)),
        (a_lo - b_hi, -d, -0.5 * (b_hi + a_lo)),
    ] {
        if gap > 0.0 {
            return Ok(Separation::Apart);
        }
        let plane = Plane {
            d,
            level,
            slack,
            gap,
        };
        if gap >= -slack {
            found = Separation::Touching(plane);
        } else if gap >= -tolerance && !matches!(found, Separation::Touching(_)) {
            found = Separation::Near(Plane {
                slack: slack - gap,
                ..plane
            });
        }
    }
    Ok(found)
}

#[derive(Debug, Clone, Copy)]
enum Separation {
    Apart,
    Touching(Plane),
    Near(Plane),
    Overlapping,
}

/// The faces of the side that may hold an open patch of the plane: the
/// planar faces lying in it, and whether a B-spline face reaches it (no
/// other family holds an open patch of a plane).
fn patch_faces(side: &Side<'_>, plane: &Plane) -> Result<(Vec<usize>, bool), ExactMeasureError> {
    let band = (
        plane.level - 2.0 * plane.slack,
        plane.level + 2.0 * plane.slack,
    );
    let mut planar = Vec::new();
    let mut spline = false;
    for element in &side.elements {
        let Shape::Face { face, .. } = element.shape else {
            continue;
        };
        let surface = surface_of(side.brep, side.brep.topology().faces()[face].surface)?;
        let (lo, hi) = side.project(element, plane.d)?;
        match surface {
            Surface::Plane(_) if lo >= band.0 && hi <= band.1 => {
                if !planar.contains(&face) {
                    planar.push(face);
                }
            }
            Surface::BSpline(_) if hi >= band.0 && lo <= band.1 => spline = true,
            _ => {}
        }
    }
    Ok((planar, spline))
}

/// Whether some face of the side may hold an open patch of the plane.
fn may_hold_patch(side: &Side<'_>, plane: &Plane) -> Result<bool, ExactMeasureError> {
    let (planar, spline) = patch_faces(side, plane)?;
    Ok(spline || !planar.is_empty())
}

/// The layout of two items touching on `plane` where both may hold a patch
/// of it: exact face contact when every such face is planar and lies, in
/// its own numbers, on one axis plane; else near contact, or a face that
/// cannot be cut.
fn contact(a: &Side<'_>, b: &Side<'_>, plane: &Plane) -> Result<Layout, ExactMeasureError> {
    let (faces_a, spline_a) = patch_faces(a, plane)?;
    let (faces_b, spline_b) = patch_faces(b, plane)?;
    if spline_a || spline_b {
        return Ok(Layout::MayShareFace);
    }
    let mut common: Option<AxisPlane> = None;
    let mut off_level = false;
    let mut not_axis_normal = false;
    for (side, faces) in [(a, &faces_a), (b, &faces_b)] {
        for &face in faces {
            match (contact::axis_plane(side.brep, face), common) {
                (Err(Uncut::Edge), _) => return Ok(Layout::MayShareFace),
                (Err(Uncut::NotAxisNormal), _) => not_axis_normal = true,
                (Ok(found), None) => common = Some(found),
                (Ok(found), Some(seen)) if found == seen => {}
                _ => off_level = true,
            }
        }
    }
    // Touching to within rounding on a plane no axis is normal to: the
    // contact cannot be shown exact there, whatever the gap.
    if not_axis_normal {
        return Ok(Layout::NotAxisNormal);
    }
    // On axis planes, but not at one coordinate: a gap or an overlap below
    // rounding, or within the tolerance.
    if off_level {
        return Ok(Layout::NearlyShareFace { gap: plane.gap });
    }
    Ok(match common {
        Some(on) => Layout::Contact {
            plane: on,
            faces_a,
            faces_b,
        },
        None => Layout::NearlyShareFace { gap: plane.gap },
    })
}

fn layout(a: &Side<'_>, b: &Side<'_>, tolerance: Tolerance) -> Result<Layout, ExactMeasureError> {
    let mut candidates = vec![Vec3::X, Vec3::Y, Vec3::Z];
    for n in plane_normals(a)?.into_iter().chain(plane_normals(b)?) {
        // `n` and `-n` are one candidate: both senses are tried.
        if candidates
            .iter()
            .all(|seen| seen.dot(n).abs() < 1.0 - 1e-12)
        {
            candidates.push(n);
        }
    }
    let linear = tolerance.linear();
    let mut touching: Vec<Plane> = Vec::new();
    let mut near: Vec<Plane> = Vec::new();
    for d in candidates {
        match separate(a, b, d, linear)? {
            Separation::Apart => return Ok(Layout::Apart),
            Separation::Touching(plane) => touching.push(plane),
            Separation::Near(plane) => near.push(plane),
            Separation::Overlapping => {}
        }
    }
    // Two touching planes with `a` below both: so it is below their sum.
    let base = touching.len();
    for i in 0..base {
        for j in i + 1..base {
            let sum = touching[i].d + touching[j].d;
            let length = sum.length();
            if length > 1e-6 {
                match separate(a, b, sum / length, linear)? {
                    Separation::Apart => return Ok(Layout::Apart),
                    Separation::Touching(plane) => touching.push(plane),
                    Separation::Near(_) | Separation::Overlapping => {}
                }
            }
        }
    }
    for plane in &touching {
        if !may_hold_patch(a, plane)? || !may_hold_patch(b, plane)? {
            return Ok(Layout::Touching);
        }
    }
    // Touching where both hold a patch of the plane: face contact, cut
    // exactly or refused by name.
    if let Some(plane) = touching.first() {
        return contact(a, b, plane);
    }
    if boundaries_apart(a, b, tolerance)? && !nested(a, b)? {
        return Ok(Layout::Apart);
    }
    // Interpenetrating by no more than the tolerance across a plane where
    // both have a face: near contact, refused by name.
    for plane in &near {
        if may_hold_patch(a, plane)? && may_hold_patch(b, plane)? {
            return Ok(Layout::NearlyShareFace { gap: plane.gap });
        }
    }
    Ok(Layout::Undecided)
}

/// Whether the boundaries are certainly apart: a positive lower bound on
/// the distance between them.
fn boundaries_apart(
    a: &Side<'_>,
    b: &Side<'_>,
    tolerance: Tolerance,
) -> Result<bool, ExactMeasureError> {
    let found = search_within(
        a.brep,
        b.brep,
        tolerance,
        Metric::Space,
        OVERLAP_STEPS,
        &mut |lower, _| lower > 0.0,
    )?;
    Ok(found.bounds.lower > 0.0)
}

/// Whether one item may lie inside the other, for two items whose
/// boundaries are apart: unless each is one solid with a boundary point
/// outside an enclosing box of the other (see the module docs).
fn nested(a: &Side<'_>, b: &Side<'_>) -> Result<bool, ExactMeasureError> {
    let one_solid = |brep: &ExactBRep| {
        let topology = brep.topology();
        topology.solids().len() == 1
            || (topology.solids().is_empty() && topology.shells().len() == 1)
    };
    if !one_solid(a.brep) || !one_solid(b.brep) {
        return Ok(true);
    }
    Ok(!sticks_out(a, b)? || !sticks_out(b, a)?)
}

/// Whether a point certainly on `a`'s boundary lies outside an enclosing
/// box of `b`.
fn sticks_out(a: &Side<'_>, b: &Side<'_>) -> Result<bool, ExactMeasureError> {
    let mut boxes = [(0.0, 0.0); 3];
    for (axis, d) in [Vec3::X, Vec3::Y, Vec3::Z].into_iter().enumerate() {
        boxes[axis] = extent(b, d)?;
    }
    let outside = |p: Point3| {
        (0..3).any(|axis| {
            let (lo, hi) = boxes[axis];
            let margin = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
            p[axis] < lo - margin || p[axis] > hi + margin
        })
    };
    if a.elements
        .iter()
        .filter_map(|element| element.witness)
        .any(outside)
    {
        return Ok(true);
    }
    // Edge ends: a corner sticking out is often the only point that does.
    let topology = a.brep.topology();
    for (index, edge) in topology.edges().iter().enumerate() {
        let Some(curve) = edge.curve.and_then(|id| a.brep.curves3().get(id.index())) else {
            continue;
        };
        let Some(span) = topology
            .edge_id_at(index)
            .and_then(|id| a.brep.edge_interval(id))
        else {
            continue;
        };
        for t in [span.start, span.end] {
            if evaluate3(curve, t).is_ok_and(outside) {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    //! Placing a cut boundary moves the items it was cut from with it: the
    //! cut faces' bound (`exact_hausdorff/cut.rs`) matches those items
    //! between bodies, and left unplaced they would match bodies placed
    //! apart as though they coincided. A rigid motion makes some uncut face
    //! tie every cut face's farthest point, so no closed form sees it.

    use super::{contact, BodySide, Cow, CutFrom, HashMap, Items};
    use axiolid_core::{Point2, Transform3, Vec3};
    use axiolid_overlay::ArcRing;

    #[test]
    fn placing_a_cut_boundary_places_its_uncut_items() {
        let plane = contact::AxisPlane {
            axis: 2,
            level: 0.5,
        };
        let mut assembler = contact::Assembler::default();
        assembler.add_region(
            &[ArcRing::from_points(&[
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(1.0, 1.0),
                Point2::new(0.0, 1.0),
            ])],
            plane,
        );
        let sheet = assembler.finish().expect("a sheet");
        let items = Items {
            brep: Cow::Owned(sheet.clone()),
            faces: vec![0],
            edges: vec![0],
            cut: Some(CutFrom {
                items: Cow::Owned(sheet),
                origin: HashMap::new(),
                partners: HashMap::new(),
            }),
        };
        let shift = Vec3::new(2.0, -1.0, 3.0);
        let placed = items
            .placed(Transform3::from_translation(shift), BodySide::First)
            .expect("rigid");
        let first = |brep: &axiolid_brep::ExactBRep| brep.topology().vertices()[0].position;
        let cut = placed.cut.as_ref().expect("still cut");
        assert_eq!(first(&placed.brep), first(&cut.items));
        assert!((first(&cut.items).z - 3.5).abs() < 1e-12);
    }
}
