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
//! # Hausdorff
//!
//! [`body_boundary_hausdorff_distance`] and
//! [`one_sided_body_boundary_hausdorff`] measure the boundary of the UNION
//! of each body's items. They are computed on the union of the items'
//! boundaries, so every body is first checked, item pair by item pair, for
//! a layout under which the two are the same point set:
//!
//! - **Apart**: shown by a plane strictly between them -- a direction along
//!   which the exact projection ranges of every face patch and edge span of
//!   the two are disjoint, tried along the axes and every planar face's
//!   normal -- or, failing that, by a positive certified lower bound on the
//!   distance between their boundaries together with a boundary point of
//!   each outside an enclosing box of the other. The second test needs each
//!   item to be one solid (one outer shell, voids allowed): two such solids
//!   whose boundaries do not meet are disjoint unless the outer shell of
//!   one lies inside the other, and then its box lies in the other's box.
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
//!
//! Anything else is refused by name: [`BodyMeasureError::ItemsShareFace`]
//! when the items touch where both have a face in the separating plane --
//! a shared patch of face is interior to the union and NOT on its boundary,
//! and the union boundary is not formed here -- and
//! [`BodyMeasureError::ItemsOverlap`] when they were not shown apart or
//! touching: boundaries that cross or come closer than the search resolves
//! with no plane between them, one item possibly inside another, or curved
//! items touching (no planar face to separate them by). Nothing is assumed
//! disjoint.
//!
//! Once the layout is shown, the one-sided distance is
//! `max over p in dA of min over items B_j of d(p, dB_j)`, measured by the
//! search of [`crate::one_sided_boundary_hausdorff`] over all items at
//! once. Faces are matched across every pair of items, so a body moved by
//! a translation closes as fast as one solid does: every face patch is
//! held at `|t|` by its own translate, and the support point of the union
//! against `t` is `|t|` from the moved body.

use std::borrow::Cow;

use core::fmt;

use axiolid_brep::{ExactBRep, ExactBRepBuilder};
use axiolid_core::{Point3, Scalar, Tolerance, Vec3};
use axiolid_evaluate::evaluate3;
use axiolid_surface::Surface;

use crate::exact::ExactMeasureError;
use crate::exact_distance::{
    search, search_within, surface_of, Clearance, DistanceBounds, Metric, Shape, Side,
    OVERLAP_STEPS,
};
use crate::exact_hausdorff::{witnessed, MAX_SPLITS};
use crate::mesh_hausdorff::HausdorffBounds;

/// Which argument of a body query a refusal is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BodySide {
    /// The first slice: `a`, or `from` of a one-sided query.
    First,
    /// The second slice: `b`, or `to` of a one-sided query.
    Second,
}

/// Why two bodies of exact solids could not be measured.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
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
    /// that separates them, so they may share a patch of face: interior to
    /// their union and not on its boundary, which is not formed here.
    ItemsShareFace {
        /// Which body.
        body: BodySide,
        /// The first item, by index in the slice.
        first: usize,
        /// The second item, by index in the slice.
        second: usize,
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
                 face, which is not on the boundary of their union"
            ),
        }
    }
}

impl std::error::Error for BodyMeasureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Measure(error) => Some(error),
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
pub fn body_boundary_distance(
    a: &[ExactBRep],
    b: &[ExactBRep],
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<BodyDistance, BodyMeasureError> {
    let accuracy = accuracy.max(0.0);
    distance(a, b, tolerance, &mut |lower, upper| {
        upper - lower <= accuracy
    })
}

/// [`body_boundary_distance`] refined only until it clears `limit`.
///
/// # Errors
///
/// As [`body_boundary_distance`].
pub fn body_boundary_clearance(
    a: &[ExactBRep],
    b: &[ExactBRep],
    limit: Scalar,
    tolerance: Tolerance,
) -> Result<(BodyDistance, Clearance), BodyMeasureError> {
    let found = distance(a, b, tolerance, &mut |lower, upper| {
        upper < limit || lower > limit
    })?;
    let clearance = found.bounds.against(limit);
    Ok((found, clearance))
}

fn distance(
    a: &[ExactBRep],
    b: &[ExactBRep],
    tolerance: Tolerance,
    done: &mut dyn FnMut(Scalar, Scalar) -> bool,
) -> Result<BodyDistance, BodyMeasureError> {
    let a = Items::new(a, BodySide::First)?;
    let b = Items::new(b, BodySide::Second)?;
    let found = search(&a.brep, &b.brep, tolerance, Metric::Space, done)?;
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
pub fn one_sided_body_boundary_hausdorff(
    from: &[ExactBRep],
    to: &[ExactBRep],
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
pub fn one_sided_body_boundary_hausdorff_with_budget(
    from: &[ExactBRep],
    to: &[ExactBRep],
    accuracy: Scalar,
    tolerance: Tolerance,
    max_splits: usize,
) -> Result<BodyHausdorffBounds, BodyMeasureError> {
    let source = Items::new(from, BodySide::First)?;
    let target = Items::new(to, BodySide::Second)?;
    check_layout(from, BodySide::First, tolerance)?;
    check_layout(to, BodySide::Second, tolerance)?;
    one_sided(&source, &target, accuracy, tolerance, max_splits)
}

/// Two-sided Hausdorff distance between the boundaries of the unions of two
/// bodies' items, to within `accuracy` (#229).
///
/// # Errors
///
/// As [`one_sided_body_boundary_hausdorff`].
pub fn body_boundary_hausdorff_distance(
    a: &[ExactBRep],
    b: &[ExactBRep],
    accuracy: Scalar,
    tolerance: Tolerance,
) -> Result<BodyHausdorff, BodyMeasureError> {
    let first = Items::new(a, BodySide::First)?;
    let second = Items::new(b, BodySide::Second)?;
    check_layout(a, BodySide::First, tolerance)?;
    check_layout(b, BodySide::Second, tolerance)?;
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
    let found = witnessed(
        &from.brep,
        &to.brep,
        accuracy,
        tolerance,
        max_splits,
        &from.edge_ranges(),
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
}

impl<'a> Items<'a> {
    fn new(items: &'a [ExactBRep], body: BodySide) -> Result<Self, BodyMeasureError> {
        match items {
            [] => Err(BodyMeasureError::EmptyBody { body }),
            [one] => Ok(Self {
                brep: Cow::Borrowed(one),
                faces: vec![0],
                edges: vec![0],
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
                })
            }
        }
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Layout {
    /// Shown disjoint.
    Apart,
    /// Shown touching, with no patch of face shared.
    Touching,
    /// Touching where both may hold a patch of the separating plane.
    MayShareFace,
    /// Neither shown.
    Undecided,
}

/// Refuse a body two of whose items are not shown apart or touching without
/// a shared patch.
fn check_layout(
    items: &[ExactBRep],
    body: BodySide,
    tolerance: Tolerance,
) -> Result<(), BodyMeasureError> {
    if items.len() < 2 {
        return Ok(());
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
    Ok(())
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
/// side, to within `slack`.
#[derive(Debug, Clone, Copy)]
struct Plane {
    d: Vec3,
    level: Scalar,
    slack: Scalar,
}

/// How the items' ranges along the unit `d` lie: strictly apart either way
/// round, within rounding of a plane between them, or overlapping.
fn separate(a: &Side<'_>, b: &Side<'_>, d: Vec3) -> Result<Separation, ExactMeasureError> {
    let (a_lo, a_hi) = extent(a, d)?;
    let (b_lo, b_hi) = extent(b, d)?;
    let scale = a_lo.abs().max(a_hi.abs()).max(b_lo.abs()).max(b_hi.abs());
    let slack = 1e-9 * (1.0 + scale);
    let mut found = Separation::Overlapping;
    for (gap, d, level) in [
        (b_lo - a_hi, d, 0.5 * (a_hi + b_lo)),
        (a_lo - b_hi, -d, -0.5 * (b_hi + a_lo)),
    ] {
        if gap > slack {
            return Ok(Separation::Apart);
        }
        if gap >= -slack {
            found = Separation::Touching(Plane { d, level, slack });
        }
    }
    Ok(found)
}

#[derive(Debug, Clone, Copy)]
enum Separation {
    Apart,
    Touching(Plane),
    Overlapping,
}

/// Whether some face of the side may hold an open patch of the plane: a
/// planar face lying in it, or a B-spline face reaching it.
fn may_hold_patch(side: &Side<'_>, plane: &Plane) -> Result<bool, ExactMeasureError> {
    let band = (
        plane.level - 2.0 * plane.slack,
        plane.level + 2.0 * plane.slack,
    );
    for element in &side.elements {
        let Shape::Face { face, .. } = element.shape else {
            continue;
        };
        let surface = surface_of(side.brep, side.brep.topology().faces()[face].surface)?;
        let (lo, hi) = side.project(element, plane.d)?;
        let holds = match surface {
            Surface::Plane(_) => lo >= band.0 && hi <= band.1,
            Surface::BSpline(_) => hi >= band.0 && lo <= band.1,
            _ => false,
        };
        if holds {
            return Ok(true);
        }
    }
    Ok(false)
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
    let mut touching: Vec<Plane> = Vec::new();
    for d in candidates {
        match separate(a, b, d)? {
            Separation::Apart => return Ok(Layout::Apart),
            Separation::Touching(plane) => touching.push(plane),
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
                match separate(a, b, sum / length)? {
                    Separation::Apart => return Ok(Layout::Apart),
                    Separation::Touching(plane) => touching.push(plane),
                    Separation::Overlapping => {}
                }
            }
        }
    }
    for plane in &touching {
        if !may_hold_patch(a, plane)? || !may_hold_patch(b, plane)? {
            return Ok(Layout::Touching);
        }
    }
    if boundaries_apart(a, b, tolerance)? && !nested(a, b)? {
        return Ok(Layout::Apart);
    }
    Ok(if touching.is_empty() {
        Layout::Undecided
    } else {
        Layout::MayShareFace
    })
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
