//! Curve nodes compiled to 3D polylines with a certified chord bound
//! (#311): the reference [`PolylineCompiler`].
//!
//! # What compiles
//!
//! - A `Curve3` over its domain, and a `Curve2` over its domain lifted to
//!   `z = 0` (the plane a 2D curve is placed from, as a curve placed at a
//!   station is).
//! - A trim, a composite or a constant offset of such curves (#317),
//!   nested: see [Curve relations](#curve-relations).
//! - An `OpenProfile`: its path, which graph validation keeps a bounded
//!   open 2D curve, so a path that is a relation compiles too.
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
//! An intrinsic curve (#317) is bounded by its curvature law: it runs at
//! unit speed, so `|c''| = |k|`, and the law's closed-form supremum over a
//! span ([`axiolid_reference::bound::intrinsic_chord_bound2`], `3`) gives
//! `h^2 / 8 sup |k|`, across a seam of a piecewise law too. Its vertices
//! are read by the evaluator's quadrature, whose error is its stated
//! estimate, as an arc-length chain's intrinsic pieces' are. A 3D one
//! compiles only with its torsion law identically zero
//! (`Intrinsic3::is_planar`), as the certified station lengths take it
//! (#302): with torsion its frame is a truncated Magnus series.
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
//! # Curve relations
//!
//! - **A trim** of an atomic curve is that curve flattened between the trim
//!   parameters, with its own chord bounds. The parameters are read as a
//!   station basis reads them (`crate::station::trim_parameters`): a
//!   parameter, an arc length, or a point inverted exactly onto the basis;
//!   a closed conic wraps from the start selector by the trim's sense,
//!   across its seam if it must. The trim runs in its sense.
//! - **A composite** is each segment compiled to the budget, in its sense,
//!   end to end: the earlier segment's last vertex is the later one's
//!   first, one shared vertex, and so is the first vertex where the last
//!   segment closes onto it. A joint whose ends differ by `g`, at most the
//!   options' linear tolerance, moves that vertex by `g`, so the bound is
//!   the worst segment's plus the largest such `g`; when that misses the
//!   budget the segments are compiled once more within the budget less
//!   `g`. Ends further apart are refused by name: a polyline is one
//!   connected curve.
//! - **A constant offset** `B + d u` (`u` the unit left normal of a 2D
//!   basis, `normalise(V x B')` beside a 3D one with reference direction
//!   `V`) takes its basis's parameter, so a trim or a composite under it is
//!   offset span by span, the offset's left following the direction each
//!   piece is traversed in. Its chord bound is `h^2 / 8 sup |P''|` from the
//!   basis's derivative bounds through the unit-vector rule the certified
//!   station lengths use (#294, #307;
//!   [`axiolid_reference::bound::offset_chord_bound2`], `3`). Its pieces
//!   are cut where the basis may not be `C^2`, and a span whose tangent the
//!   bound cannot keep away from zero is halved.
//!
//! # What is refused, by name
//!
//! An unbounded curve (a line, an elevated curve over a line plan), a
//! family whose flattening is not certified
//! ([`axiolid_reference::bound::certifies_flattening2`]: implicit and banked
//! curves, a 3D intrinsic curve with torsion, other numerically evaluated
//! curves), a trim of anything but an atomic curve, a composite whose
//! segments do not meet, an offset of an offset or of a placed curve, an
//! offset across a corner of its basis or of a family it has no bound for,
//! an offset by stations, a surface or parameter curve, a curve placed at a
//! numerically resolved station, a collection, and any node that is not a
//! curve.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Interval, Point2, Point3, Scalar, Transform3, Vec3};
use axiolid_curve::{Curve2, Curve3};
use axiolid_mesh_compile_contract::{ChordDeviation, PolylineOutcome};
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraph, GeometryNode, NodeId, TrimSelector,
    TrimmingPreference,
};
use axiolid_reference::bound;
use axiolid_reference::curve::{domain2, domain3, evaluate2, evaluate3};
use axiolid_reference::station::StationCurve;

/// Nesting of placements, relations and open profiles a root may hold.
const MAX_DEPTH: usize = 256;

/// Bisections of one span before its bound must fit.
const MAX_SPLITS: u32 = 48;

/// Vertices of one polyline.
const MAX_POINTS: usize = 1 << 20;

/// Relative slack within which a trim parameter just past its basis's
/// domain (a point inverted onto an end) is read as the end.
const DOMAIN_SLACK: Scalar = 1e-9;

/// Why a composite whose segments do not meet is refused.
const COMPOSITE_GAP: &str = "a composite curve whose consecutive segments do not meet within the \
                             tolerance: a polyline is one connected curve";

/// Why a trim of a relation or a placement is refused.
const TRIM_OF_RELATION: &str = "a trim of a curve relation or of a placed curve: only an atomic \
                                curve's parameters are trimmed to a polyline";

/// Why an offset across a corner of its basis is refused.
const OFFSET_ACROSS_CORNER: &str = "an offset curve across a vertex or corner of its basis: the \
                                    offset is not continuous there";

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

/// A constant offset carried down to the curve it offsets (#317): its
/// distance, signed for the direction that curve is traversed in, and the
/// reference direction of a 3D offset.
#[derive(Debug, Clone, Copy)]
struct Offset {
    distance: Scalar,
    reference: Option<Vec3>,
}

impl Offset {
    /// The offset as the basis traversed backwards reads it: the left of
    /// the reversed direction (and `V x -T`) is the other side.
    fn reversed(self) -> Self {
        Self {
            distance: -self.distance,
            ..self
        }
    }
}

/// `offset` as a piece traversed `forward` or backwards reads it.
fn sensed(offset: Option<Offset>, forward: bool) -> Option<Offset> {
    if forward {
        offset
    } else {
        offset.map(Offset::reversed)
    }
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
        return Err(too_deep());
    }
    match graph.get(id) {
        Some(GeometryNode::Curve3(curve)) => atomic3(curve, None, chord, None),
        Some(GeometryNode::Curve2(curve)) => atomic2(curve, None, chord, None),
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
        Some(GeometryNode::CurveRelation(relation)) => {
            curve_relation(graph, relation, chord, options, depth, None)
        }
        Some(GeometryNode::Collection(_)) => Err(refused(
            "a collection: a polyline is one curve, so compile each member on its own",
        )),
        other => Err(not_a_curve(id, other.is_some())),
    }
}

/// A child of a relation: resolved as any root, or, under an offset, an
/// atomic curve or a relation of them offset piece by piece.
fn under(
    graph: &GeometryGraph,
    id: NodeId,
    chord: Scalar,
    options: &ExecutionOptions,
    depth: usize,
    offset: Option<Offset>,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    let Some(offset) = offset else {
        return resolve(graph, id, chord, options, depth);
    };
    if depth > MAX_DEPTH {
        return Err(too_deep());
    }
    match graph.get(id) {
        Some(GeometryNode::Curve2(curve)) => atomic2(curve, None, chord, Some(offset)),
        Some(GeometryNode::Curve3(curve)) => atomic3(curve, None, chord, Some(offset)),
        Some(GeometryNode::CurveRelation(relation)) => {
            curve_relation(graph, relation, chord, options, depth, Some(offset))
        }
        Some(GeometryNode::Instance(_) | GeometryNode::InstanceAtStation(_)) => Err(refused(
            "an offset of a placed curve: only atomic curves and their trims and composites \
             are offset to a polyline",
        )),
        other => Err(not_a_curve(id, other.is_some())),
    }
}

/// A trim, a composite or an offset (#317), under `offset` when one of its
/// ancestors is an offset.
fn curve_relation(
    graph: &GeometryGraph,
    relation: &CurveRelation,
    chord: Scalar,
    options: &ExecutionOptions,
    depth: usize,
    offset: Option<Offset>,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    match relation {
        CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement,
            preference,
        } => trimmed(
            graph,
            *basis,
            (start, end),
            *sense_agreement,
            *preference,
            chord,
            offset,
        ),
        CurveRelation::Composite { segments } => {
            composite(graph, segments, chord, options, depth, offset)
        }
        CurveRelation::Offset {
            basis,
            distance,
            reference_direction,
        } => {
            if offset.is_some() {
                return Err(refused(
                    "an offset of an offset curve: only atomic curves and their trims and \
                     composites are offset to a polyline",
                ));
            }
            if !distance.is_finite() {
                return Err(GeomError::InvalidInput(
                    "an offset curve's distance is not finite".into(),
                ));
            }
            let offset = Offset {
                distance: *distance,
                reference: *reference_direction,
            };
            under(graph, *basis, chord, options, depth + 1, Some(offset))
        }
        CurveRelation::OffsetByStations { .. } | CurveRelation::OffsetsAtStations { .. } => {
            Err(refused(
                "an offset by stations: its law is read along its basis's length, which no \
                 chord bound here reaches",
            ))
        }
        _ => Err(refused(
            "a surface curve or a parameter curve: only atomic curves and their trims, \
             composites and offsets compile to polylines",
        )),
    }
}

/// A trim of an atomic curve: the basis flattened between the trim
/// parameters, in the trim's sense.
fn trimmed(
    graph: &GeometryGraph,
    basis: NodeId,
    (start, end): (&[TrimSelector], &[TrimSelector]),
    sense: bool,
    preference: TrimmingPreference,
    chord: Scalar,
    offset: Option<Offset>,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    let offset = sensed(offset, sense);
    let (mut points, bound) = match graph.get(basis) {
        Some(GeometryNode::Curve2(curve)) => {
            let window = crate::station::trim_parameters(
                StationCurve::Two(curve),
                start,
                end,
                sense,
                preference,
            )?;
            atomic2(curve, Some(window), chord, offset)?
        }
        Some(GeometryNode::Curve3(curve)) => {
            let window = crate::station::trim_parameters(
                StationCurve::Three(curve),
                start,
                end,
                sense,
                preference,
            )?;
            atomic3(curve, Some(window), chord, offset)?
        }
        Some(_) => return Err(refused(TRIM_OF_RELATION)),
        None => return Err(not_a_curve(basis, false)),
    };
    if !sense {
        points.reverse();
    }
    Ok((points, bound))
}

/// A composite: its segments end to end, each joint one shared vertex.
fn composite(
    graph: &GeometryGraph,
    segments: &[CurveSegment],
    chord: Scalar,
    options: &ExecutionOptions,
    depth: usize,
    offset: Option<Offset>,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    if segments.is_empty() {
        return Err(GeomError::InvalidInput(
            "a composite curve has no segments".into(),
        ));
    }
    let tolerance = options.tolerance().linear();
    let mut budget = chord;
    loop {
        let mut points: Vec<Point3> = Vec::new();
        let (mut worst, mut gap): (Scalar, Scalar) = (0.0, 0.0);
        for segment in segments {
            let law = sensed(offset, segment.same_sense);
            let (mut child, bound) = under(graph, segment.curve, budget, options, depth + 1, law)?;
            if !segment.same_sense {
                child.reverse();
            }
            worst = worst.max(bound);
            if let (Some(&last), Some(&first)) = (points.last(), child.first()) {
                let joint = last.distance(first);
                if joint.is_nan() || joint > tolerance {
                    return Err(refused(COMPOSITE_GAP));
                }
                gap = gap.max(joint);
                child.remove(0);
            }
            points.append(&mut child);
            if points.len() > MAX_POINTS {
                return Err(GeomError::BudgetExceeded {
                    resource: "polyline vertices",
                });
            }
        }
        if segments.len() > 1 && points.len() > 2 {
            let (first, last) = (points[0], points[points.len() - 1]);
            let closing = first.distance(last);
            if closing > 0.0 && closing <= tolerance {
                gap = gap.max(closing);
                if let Some(end) = points.last_mut() {
                    *end = first;
                }
            }
        }
        let bound = if gap > 0.0 { worst + gap } else { worst };
        if bound <= chord || budget < chord || gap >= chord {
            return Ok((points, bound));
        }
        budget = chord - gap;
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

fn too_deep() -> GeomError {
    GeomError::BudgetExceeded {
        resource: "polyline placement depth",
    }
}

/// A node that is not a curve, or (`present` false) not in the graph.
fn not_a_curve(id: NodeId, present: bool) -> GeomError {
    GeomError::InvalidInput(if present {
        format!("node {id:?} is not a curve: it has no polyline")
    } else {
        format!("node {id:?} does not belong to this graph")
    })
}

fn lift(p: Point2) -> Point3 {
    Point3::new(p.x, p.y, 0.0)
}

/// `breaks` and the seams of an intrinsic curve's law, if any.
fn with_seams(mut breaks: Vec<Scalar>, seams: Option<Vec<Scalar>>) -> Vec<Scalar> {
    breaks.extend(seams.unwrap_or_default());
    breaks
}

/// A 2D curve over its domain or a trim's `window`, lifted to `z = 0`,
/// offset when `offset` is given.
fn atomic2(
    curve: &Curve2,
    window: Option<(Scalar, Scalar)>,
    chord: Scalar,
    offset: Option<Offset>,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    let periodic = matches!(curve, Curve2::Circle(_) | Curve2::Ellipse(_));
    match (curve, window, offset) {
        (Curve2::Polyline(polyline), None, None) => {
            return Ok((
                own_vertices(polyline.points.iter().map(|&p| lift(p)), polyline.closed)?,
                0.0,
            ))
        }
        (Curve2::Line(_), None, _) => {
            return Err(refused("an unbounded 2D line: it has no finite polyline"))
        }
        _ => {}
    }
    let unbounded = matches!(curve, Curve2::Line(_));
    let span = span(domain2(curve), window, periodic, unbounded)?;
    let full_turn = periodic && span.1 - span.0 >= core::f64::consts::TAU;
    let Some(offset) = offset else {
        if let Curve2::Intrinsic(intrinsic) = curve {
            if bound::intrinsic_chord_bound2(intrinsic, 0.0, intrinsic.length).is_none() {
                return Err(refused(uncertified2(curve)));
            }
            return bisected(
                |t| evaluate2(curve, t).map(lift),
                |a, b| bound::intrinsic_chord_bound2(intrinsic, a, b),
                span,
                &[],
                chord,
                false,
            );
        }
        if !bound::certifies_flattening2(curve) {
            return Err(refused(uncertified2(curve)));
        }
        return bisected(
            |t| evaluate2(curve, t).map(lift),
            |a, b| bound::chord_bound2(curve, a, b),
            span,
            &bound::continuity_breaks2(curve, 1),
            chord,
            full_turn,
        );
    };
    if offset.reference.is_some() {
        return Err(refused(
            "a 2D curve offset along a reference direction: a 2D offset runs along the left \
             normal",
        ));
    }
    let seams = match curve {
        Curve2::Intrinsic(intrinsic) => {
            if bound::intrinsic_chord_bound2(intrinsic, 0.0, intrinsic.length).is_none() {
                return Err(refused(uncertified2(curve)));
            }
            Some(intrinsic.curvature.seams_within(intrinsic.length))
        }
        Curve2::Chain(_) | Curve2::Implicit(_) => {
            return Err(refused(
                "an offset of an arc-length chain or an implicit 2D curve: its derivative \
                 bounds are not composed with an offset's",
            ))
        }
        _ if !bound::certifies_flattening2(curve) => return Err(refused(uncertified2(curve))),
        _ => None,
    };
    if inside(&bound::continuity_breaks2(curve, 1), span) {
        return Err(refused(OFFSET_ACROSS_CORNER));
    }
    let distance = offset.distance;
    bisected(
        |t| bound::offset_point2(curve, distance, t).map(lift),
        |a, b| Some(bound::offset_chord_bound2(curve, distance, a, b).unwrap_or(Scalar::INFINITY)),
        span,
        &with_seams(bound::continuity_breaks2(curve, 2), seams),
        chord,
        full_turn,
    )
}

/// [`atomic2`] for a 3D curve; an offset needs its reference direction.
fn atomic3(
    curve: &Curve3,
    window: Option<(Scalar, Scalar)>,
    chord: Scalar,
    offset: Option<Offset>,
) -> GeomResult<(Vec<Point3>, Scalar)> {
    let periodic = matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_));
    match (curve, window, offset) {
        (Curve3::Polyline(polyline), None, None) => {
            return Ok((
                own_vertices(polyline.points.iter().copied(), polyline.closed)?,
                0.0,
            ))
        }
        (Curve3::Line(_), None, _) => {
            return Err(refused("an unbounded 3D line: it has no finite polyline"))
        }
        _ => {}
    }
    let unbounded = matches!(curve, Curve3::Line(_));
    let span = span(domain3(curve), window, periodic, unbounded)?;
    let full_turn = periodic && span.1 - span.0 >= core::f64::consts::TAU;
    let Some(offset) = offset else {
        if let Curve3::Intrinsic(intrinsic) = curve {
            if !intrinsic.is_planar()
                || bound::intrinsic_chord_bound3(intrinsic, 0.0, intrinsic.length).is_none()
            {
                return Err(refused(uncertified3(curve)));
            }
            return bisected(
                |t| evaluate3(curve, t),
                |a, b| bound::intrinsic_chord_bound3(intrinsic, a, b),
                span,
                &[],
                chord,
                false,
            );
        }
        if !bound::certifies_flattening3(curve) {
            return Err(refused(uncertified3(curve)));
        }
        return bisected(
            |t| evaluate3(curve, t),
            |a, b| bound::chord_bound3(curve, a, b),
            span,
            &bound::continuity_breaks3(curve, 1),
            chord,
            full_turn,
        );
    };
    let Some(reference) = offset.reference else {
        return Err(refused(
            "a 3D curve offset without a reference direction: a 3D offset runs along V x T",
        ));
    };
    if matches!(curve, Curve3::Intrinsic(_)) {
        return Err(refused(
            "an offset of an intrinsic 3D curve: its derivative bounds are not composed with \
             an offset's",
        ));
    }
    if !bound::certifies_flattening3(curve) {
        return Err(refused(uncertified3(curve)));
    }
    if inside(&bound::continuity_breaks3(curve, 1), span) {
        return Err(refused(OFFSET_ACROSS_CORNER));
    }
    let distance = offset.distance;
    bisected(
        |t| bound::offset_point3(curve, distance, reference, t),
        |a, b| {
            Some(
                bound::offset_chord_bound3(curve, distance, reference, a, b)
                    .unwrap_or(Scalar::INFINITY),
            )
        },
        span,
        &bound::continuity_breaks3(curve, 2),
        chord,
        full_turn,
    )
}

/// Whether a break lies strictly inside `(lo, hi)`.
fn inside(breaks: &[Scalar], (lo, hi): (Scalar, Scalar)) -> bool {
    breaks.iter().any(|&t| t > lo && t < hi)
}

/// The parameters `lo < hi` a curve is flattened over: its whole domain,
/// refused unless finite and not empty, or a trim's `window`, which must
/// lie in the domain (to [`DOMAIN_SLACK`]) unless the curve is periodic
/// (a trim may cross its seam) or `unbounded`.
fn span(
    domain: Interval,
    window: Option<(Scalar, Scalar)>,
    periodic: bool,
    unbounded: bool,
) -> GeomResult<(Scalar, Scalar)> {
    let Some((lo, hi)) = window else {
        let domain = bounded_domain(domain)?;
        return Ok((domain.start.min(domain.end), domain.start.max(domain.end)));
    };
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err(GeomError::Degenerate(
            "a trimmed curve has an empty or non-finite parameter interval".into(),
        ));
    }
    if periodic || unbounded {
        return Ok((lo, hi));
    }
    let (start, end) = (domain.start.min(domain.end), domain.start.max(domain.end));
    let slack = |value: Scalar| DOMAIN_SLACK * value.abs().max(1.0);
    if (start.is_finite() && lo < start - slack(start))
        || (end.is_finite() && hi > end + slack(end))
    {
        return Err(GeomError::InvalidInput(
            "a trim's parameters lie outside its basis curve's domain".into(),
        ));
    }
    let (lo, hi) = (lo.max(start), hi.min(end));
    if lo >= hi {
        return Err(GeomError::Degenerate(
            "a trimmed curve has an empty parameter interval".into(),
        ));
    }
    Ok((lo, hi))
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

/// Bisect each smooth piece of `[lo, hi]` until every span's certified
/// chord bound fits `chord`; the polyline through the span ends and the
/// largest span bound. A full turn of a conic ends on its first vertex
/// exactly.
fn bisected(
    eval: impl Fn(Scalar) -> GeomResult<Point3>,
    bound: impl Fn(Scalar, Scalar) -> Option<Scalar>,
    (lo, hi): (Scalar, Scalar),
    breaks: &[Scalar],
    chord: Scalar,
    full_turn: bool,
) -> GeomResult<(Vec<Point3>, Scalar)> {
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
    if full_turn {
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
            "an intrinsic 2D curve whose curvature law has no closed-form bound: its chords \
             have no certified bound"
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
            "an intrinsic 3D curve with torsion, or whose curvature law has no closed-form \
             bound: its Frenet frame is integrated by a truncated series, so its chords have \
             no certified bound"
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
