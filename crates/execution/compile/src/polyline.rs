//! Curve nodes compiled to 3D polylines with a certified chord bound
//! (#311): the reference [`PolylineCompiler`].
//!
//! # What compiles
//!
//! - A `Curve3` over its domain, and a `Curve2` over its domain lifted to
//!   `z = 0` (the plane a 2D curve is placed from, as a curve placed at a
//!   station is).
//! - An `OpenProfile`: its path, which graph validation keeps a bounded
//!   open 2D curve.
//! - An `Instance` of any of these, and an `InstanceAtStation` whose frame
//!   is exact (a line basis, #264): the source's polyline moved by the
//!   placement.
//!
//! # The bound
//!
//! A polyline or a bounded line is its own exact curve: `0`. Any other
//! family is cut at the parameters where it may not be `C^1`
//! ([`axiolid_reference::bound::continuity_breaks2`]), and each smooth
//! piece is bisected until every span's certified chord bound
//! ([`axiolid_reference::bound::chord_bound2`], `chord_bound3`: the exact
//! sagitta of a circular arc, else `h^2 / 8 sup |c''|` from certified
//! derivative bounds) is within the budget. Every point of the curve on a
//! span lies within that span's bound of its chord, so the largest span
//! bound bounds the whole curve; that is the reported value. Vertices are
//! the curve evaluated at the span ends, so the rounding of evaluating them
//! is not part of the bound, as the rounding of computing a mesh's own
//! vertices is not (`crate::deviation`).
//!
//! A placement moves the source's polyline rigidly or scales it: the source
//! is compiled to the budget shrunk by an upper bound on the transform's
//! stretch ([`crate::deviation::shrunk_budget`]), and its bound grows by
//! that stretch, exactly as an instance's mesh deviation does (#258). It
//! also carries the rounding of applying the transform in `f64` (#295,
//! [`crate::deviation`]'s "Placement rounding"): the source leaves room for
//! it, the bound adds it, and a budget at or below it is refused with
//! [`crate::deviation::BUDGET_BELOW_PLACEMENT_ROUNDING`].
//!
//! # What is refused, by name
//!
//! An unbounded curve (a line, an elevated curve over a line plan), a
//! family whose flattening is not certified
//! ([`axiolid_reference::bound::certifies_flattening2`]: intrinsic, implicit,
//! banked and other numerically evaluated curves), a curve relation (trims,
//! composites, offsets), a curve placed at a numerically resolved station,
//! a collection, and any node that is not a curve.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Interval, Point3, Scalar, Transform3};
use axiolid_curve::{Curve2, Curve3};
use axiolid_mesh_compile_contract::{ChordDeviation, PolylineOutcome};
use axiolid_model::{GeometryGraph, GeometryNode, NodeId};

/// Nesting of placements and open profiles a root may hold.
const MAX_DEPTH: usize = 256;

/// Bisections of one span before its bound must fit.
const MAX_SPLITS: u32 = 48;

/// Vertices of one polyline.
const MAX_POINTS: usize = 1 << 20;

/// Compile `root` to a polyline within the options' chord budget.
pub(crate) fn compile(
    graph: &GeometryGraph,
    root: NodeId,
    options: &ExecutionOptions,
) -> GeomResult<PolylineOutcome> {
    let requested = crate::compiler::chord_error(options);
    let (points, bound) = resolve(graph, root, requested, options, 0)?;
    Ok(PolylineOutcome::new(
        points,
        ChordDeviation::new(requested, bound),
    ))
}

/// A node's polyline and its bound, within `chord` in its own units.
fn resolve(
    graph: &GeometryGraph,
    id: NodeId,
    chord: Scalar,
    options: &ExecutionOptions,
    depth: usize,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    if depth > MAX_DEPTH {
        return Err(GeomError::BudgetExceeded {
            resource: "polyline placement depth",
        });
    }
    match graph.get(id) {
        Some(GeometryNode::Curve3(curve)) => flatten3(curve, chord),
        Some(GeometryNode::Curve2(curve)) => flatten2(curve, chord),
        Some(GeometryNode::OpenProfile(profile)) => {
            resolve(graph, profile.path, chord, options, depth + 1)
        }
        Some(GeometryNode::Instance(instance)) => placed(
            graph,
            instance.source,
            instance.transform,
            chord,
            options,
            depth,
        ),
        Some(GeometryNode::InstanceAtStation(at)) => {
            let placement = crate::station::placement(graph, id)?;
            if !placement.exact {
                return Err(refused(
                    "a curve placed at a station whose frame is read at a numerically \
                     resolved distance: the placement has no certified bound",
                ));
            }
            placed(graph, at.source, placement.transform, chord, options, depth)
        }
        Some(GeometryNode::CurveRelation(_)) => Err(refused(
            "a curve relation (trim, composite, offset): only atomic curves, open profiles \
             of them and their placements compile to polylines",
        )),
        Some(GeometryNode::Collection(_)) => Err(refused(
            "a collection: a polyline is one curve, so compile each member on its own",
        )),
        Some(_) => Err(GeomError::InvalidInput(format!(
            "node {id:?} is not a curve: it has no polyline"
        ))),
        None => Err(GeomError::InvalidInput(format!(
            "node {id:?} does not belong to this graph"
        ))),
    }
}

/// The source's polyline moved by `transform`, compiled to the budget
/// shrunk by its stretch, less room for the rounding of applying the
/// transform in `f64` (#295), so the grown bound stays within `chord`.
/// As for a mesh (`crate::compiler`): the room starts at the floor at the
/// translation's scale plus a share for the source's own coordinates, and
/// the source is compiled once more with room for an eighth more than its
/// polyline rounded by where that carried the bound over `chord`; a
/// `chord` at or below the rounding is refused by name.
fn placed(
    graph: &GeometryGraph,
    source: NodeId,
    transform: Transform3,
    chord: Scalar,
    options: &ExecutionOptions,
    depth: usize,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    let (_, stretch) = crate::compiler::instance_local_tolerance(transform, options.tolerance())?;
    let mut room = crate::compiler::placement_room(transform, chord, 0.0)?;
    let mut retried = false;
    loop {
        let local = crate::deviation::shrunk_budget_within(chord, stretch, room);
        let (points, bound) = resolve(graph, source, local, options, depth + 1)?;
        let rounding = crate::deviation::placement_rounding(&points, transform);
        crate::compiler::refuse_below_rounding(chord, rounding)?;
        let grown = crate::deviation::grown(bound, stretch, rounding);
        if retried || rounding <= room || grown <= chord {
            let points = points
                .into_iter()
                .map(|point| transform.transform_point3(point))
                .collect();
            return Ok((points, grown));
        }
        room = room.max((rounding * 1.125).min(rounding + (chord - rounding) * 0.5));
        retried = true;
    }
}

fn refused(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: crate::BACKEND_ID,
        operation: Operation::CurveEvaluation,
        input,
    }
}

fn flatten2(curve: &Curve2, chord: Scalar) -> GeomResult<(Vec<Point3>, Scalar)> {
    let lift = |p: axiolid_core::Point2| Point3::new(p.x, p.y, 0.0);
    match curve {
        Curve2::Polyline(polyline) => Ok((
            own_vertices(polyline.points.iter().map(|&p| lift(p)), polyline.closed)?,
            0.0,
        )),
        Curve2::Line(_) => Err(refused("an unbounded 2D line: it has no finite polyline")),
        _ => {
            if !axiolid_reference::bound::certifies_flattening2(curve) {
                return Err(refused(uncertified2(curve)));
            }
            let domain = bounded_domain(axiolid_reference::curve::domain2(curve))?;
            let full_turn = matches!(curve, Curve2::Circle(_) | Curve2::Ellipse(_));
            bisected(
                |t| axiolid_reference::curve::evaluate2(curve, t).map(lift),
                |a, b| axiolid_reference::bound::chord_bound2(curve, a, b),
                domain,
                &axiolid_reference::bound::continuity_breaks2(curve, 1),
                chord,
                full_turn,
            )
        }
    }
}

fn flatten3(curve: &Curve3, chord: Scalar) -> GeomResult<(Vec<Point3>, Scalar)> {
    match curve {
        Curve3::Polyline(polyline) => Ok((
            own_vertices(polyline.points.iter().copied(), polyline.closed)?,
            0.0,
        )),
        Curve3::Line(_) => Err(refused("an unbounded 3D line: it has no finite polyline")),
        _ => {
            if !axiolid_reference::bound::certifies_flattening3(curve) {
                return Err(refused(uncertified3(curve)));
            }
            let domain = bounded_domain(axiolid_reference::curve::domain3(curve))?;
            let full_turn = matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_));
            bisected(
                |t| axiolid_reference::curve::evaluate3(curve, t),
                |a, b| axiolid_reference::bound::chord_bound3(curve, a, b),
                domain,
                &axiolid_reference::bound::continuity_breaks3(curve, 1),
                chord,
                full_turn,
            )
        }
    }
}

/// A polyline's own vertices, the first repeated last when it is closed.
fn own_vertices(points: impl Iterator<Item = Point3>, closed: bool) -> GeomResult<Vec<Point3>> {
    let mut points: Vec<Point3> = points.collect();
    if points.len() < 2 {
        return Err(GeomError::InvalidInput(
            "a polyline needs at least two points".into(),
        ));
    }
    if points.iter().any(|p| !p.is_finite()) {
        return Err(GeomError::InvalidInput(
            "a polyline point is not finite".into(),
        ));
    }
    if closed && points.first() != points.last() {
        points.push(points[0]);
    }
    Ok(points)
}

/// The curve's domain, refused unless finite and not empty.
fn bounded_domain(domain: Interval) -> GeomResult<Interval> {
    if !(domain.start.is_finite() && domain.end.is_finite()) {
        return Err(refused(
            "an unbounded curve (an elevated curve over a line plan): it has no finite polyline",
        ));
    }
    if domain.start == domain.end {
        return Err(GeomError::InvalidInput(
            "the curve's parameter domain is empty".into(),
        ));
    }
    Ok(domain)
}

/// Bisect each smooth piece of `domain` until every span's certified chord
/// bound fits `chord`; the polyline through the span ends and the largest
/// span bound. A full turn of a conic ends on its first vertex exactly.
fn bisected(
    eval: impl Fn(Scalar) -> GeomResult<Point3>,
    bound: impl Fn(Scalar, Scalar) -> Option<Scalar>,
    domain: Interval,
    breaks: &[Scalar],
    chord: Scalar,
    full_turn: bool,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    let (lo, hi) = if domain.start <= domain.end {
        (domain.start, domain.end)
    } else {
        (domain.end, domain.start)
    };
    let mut cuts: Vec<Scalar> = breaks
        .iter()
        .copied()
        .filter(|&t| t > lo && t < hi)
        .collect();
    cuts.sort_by(Scalar::total_cmp);
    cuts.dedup();
    let mut ends = Vec::with_capacity(cuts.len() + 2);
    ends.push(lo);
    ends.extend(cuts);
    ends.push(hi);

    let mut points = vec![eval(lo)?];
    let mut worst: Scalar = 0.0;
    for piece in ends.windows(2) {
        // Spans still to settle, the leftmost on top.
        let mut stack = vec![(piece[0], piece[1], 0_u32)];
        while let Some((a, b, splits)) = stack.pop() {
            let span = bound(a, b).ok_or_else(|| {
                refused("a curve span with no certified chord bound (a corner or invalid data)")
            })?;
            if span <= chord {
                points.push(eval(b)?);
                worst = worst.max(span);
                if points.len() > MAX_POINTS {
                    return Err(GeomError::BudgetExceeded {
                        resource: "polyline vertices",
                    });
                }
                continue;
            }
            let mid = 0.5 * (a + b);
            if splits >= MAX_SPLITS || !(mid > a && mid < b) {
                return Err(GeomError::BudgetExceeded {
                    resource: "polyline chord refinement",
                });
            }
            stack.push((mid, b, splits + 1));
            stack.push((a, mid, splits + 1));
        }
    }
    if full_turn && hi - lo >= core::f64::consts::TAU {
        let first = points[0];
        if let Some(last) = points.last_mut() {
            *last = first;
        }
    }
    Ok((points, worst))
}

/// The refusal of a 2D family whose flattening is not certified.
fn uncertified2(curve: &Curve2) -> &'static str {
    match curve {
        Curve2::Intrinsic(_) => {
            "an intrinsic 2D curve: it is integrated numerically, so its chords have no \
             certified bound"
        }
        Curve2::Implicit(_) => {
            "an implicit 2D curve outside a face: its chords have no certified bound"
        }
        Curve2::BSpline(_) => {
            "a rational 2D B-spline with a non-positive weight: its chords have no certified bound"
        }
        _ => "a 2D curve family whose chords have no certified bound",
    }
}

/// The refusal of a 3D family whose flattening is not certified.
fn uncertified3(curve: &Curve3) -> &'static str {
    match curve {
        Curve3::Intrinsic(_) => {
            "an intrinsic 3D curve: it is integrated numerically, so its chords have no \
             certified bound"
        }
        Curve3::Banked(_) => "a banked 3D curve: its chords have no certified bound",
        Curve3::Elevated(_) => {
            "an elevated 3D curve whose plan or profile has no certified chord bound"
        }
        Curve3::BSpline(_) => {
            "a rational 3D B-spline with a non-positive weight: its chords have no certified bound"
        }
        _ => "a 3D curve family whose chords have no certified bound",
    }
}
