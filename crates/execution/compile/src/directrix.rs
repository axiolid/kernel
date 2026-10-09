//! Safe resolution of graph-referenced 3D sweep directrices.

use axiolid_construct::sweep::SampledPath;
use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryNode, MasterRepresentation, NodeId, TrimSelector,
    TrimmingPreference,
};

mod elevated;
mod pieces;
pub(crate) use pieces::{exact_pieces, pieces};

const MAX_DEPTH: usize = 256;
const MAX_POINTS: usize = 1_000_000;
const MAX_FLATTEN_DEPTH: u32 = 16;

pub(crate) fn points(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Point3>> {
    sampled(graph, id, range, options).map(|path| path.points)
}

/// The directrix sampled to the options' chord budget, with its exact end
/// tangents when it is one smooth curve (#231, #232).
///
/// A circle or an ellipse, alone or trimmed, reports its tangents: those
/// are smooth everywhere, so a sweep may refine them until its walls fit
/// its budget. So does a B-spline with no corner knot (one of full
/// multiplicity) inside the swept span (#232). A line or polyline is its
/// own chords, and a composite or a cornered B-spline may carry corners
/// that no refinement removes, so those are swept as sampled.
pub(crate) fn sampled(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<SampledPath> {
    let (points, end_tangents) = resolve(graph, id, range, options, 0)?;
    if points.len() < 2 {
        return Err(GeomError::Degenerate(
            "a sweep directrix needs at least two points".into(),
        ));
    }
    Ok(SampledPath {
        points,
        end_tangents,
    })
}

/// Samples in the direction of travel, and the exact end tangents when
/// known (see [`sampled`]).
type Resolved = (Vec<Point3>, Option<[Vec3; 2]>);

/// Reverse a resolved path: the samples, and the tangents swapped and
/// negated so they still point along the travel.
fn reversed((mut points, ends): Resolved) -> Resolved {
    points.reverse();
    (points, ends.map(|[start, end]| [-end, -start]))
}

/// Unit tangents of a smooth curve at the ends of `[start, end]`, in the
/// direction of increasing parameter: a conic, or a B-spline with no
/// corner knot strictly inside (#232); `None` for any other family.
fn smooth_ends(curve: &axiolid_curve::Curve3, start: Scalar, end: Scalar) -> Option<[Vec3; 2]> {
    let smooth = match curve {
        axiolid_curve::Curve3::Circle(_) | axiolid_curve::Curve3::Ellipse(_) => true,
        axiolid_curve::Curve3::BSpline(_) => {
            let (lo, hi) = (start.min(end), start.max(end));
            !axiolid_reference::bound::continuity_breaks3(curve, 1)
                .into_iter()
                .any(|knot| knot > lo && knot < hi)
        }
        // Smooth unless the grade breaks inside the span (#252).
        curve if elevated::is_elevated(curve) => return elevated::smooth_ends(curve, start, end),
        _ => false,
    };
    if !smooth {
        return None;
    }
    let unit = |t| {
        axiolid_reference::curve::derivative3(curve, t)
            .ok()
            .map(Vec3::normalize_or_zero)
            .filter(|v| *v != Vec3::ZERO && v.is_finite())
    };
    Some([unit(start)?, unit(end)?])
}

fn resolve(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
    depth: usize,
) -> GeomResult<Resolved> {
    if depth > MAX_DEPTH {
        return Err(GeomError::BudgetExceeded {
            resource: "directrix relation depth",
        });
    }
    match graph.get(id) {
        Some(GeometryNode::Curve3(curve)) => sample_curve(curve, range, options),
        Some(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve {
            curve_3d,
            master: MasterRepresentation::Curve3d,
            ..
        })) => resolve(graph, *curve_3d, range, options, depth + 1),
        // Selecting curve_3d for a pcurve master can move a seam to the wrong
        // side of a periodic surface. Until pcurve evaluation and Both-agreement
        // validation exist, refusing these representations is the safe result.
        Some(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve { .. })) => {
            Err(unsupported_curve_evaluation())
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement,
            preference,
        })) => {
            // A point selector can only be inverted against a curve, so the
            // basis is resolved first. A relation basis (trimmed-of-trimmed)
            // has no single analytic curve to invert against; a parameter
            // selector still works there, a point selector is refused by name.
            let basis_curve = match graph.get(*basis) {
                Some(GeometryNode::Curve3(curve)) => Some(curve),
                _ => None,
            };
            let (a, b) = match basis_curve {
                Some(curve) => (
                    parameter(
                        start,
                        *preference,
                        "start",
                        curve,
                        *sense_agreement,
                        options.tolerance(),
                    )?,
                    parameter(
                        end,
                        *preference,
                        "end",
                        curve,
                        *sense_agreement,
                        options.tolerance(),
                    )?,
                ),
                None => (
                    parameter_only(start, *preference, "start")?,
                    parameter_only(end, *preference, "end")?,
                ),
            };
            if a == b {
                return Err(GeomError::Degenerate(
                    "trimmed directrix has an empty interval".into(),
                ));
            }
            // A periodic basis may be trimmed ACROSS its seam: the curve runs
            // from `a` the way `sense_agreement` says, wrapping if it must.
            // Sorting `a` and `b` would pick the complementary arc (#168).
            if let Some((curve, period)) = basis_curve.and_then(|c| period_of(c).map(|p| (c, p))) {
                let out =
                    sample_periodic_trim(curve, period, a, b, *sense_agreement, range, options)?;
                return Ok(if *sense_agreement { out } else { reversed(out) });
            }
            let selected = range.unwrap_or((a, b));
            let lo = a.min(b);
            let hi = a.max(b);
            let slack = options.tolerance().linear();
            if selected.0.min(selected.1) < lo - slack || selected.0.max(selected.1) > hi + slack {
                return Err(GeomError::InvalidInput(
                    "sweep range exceeds trimmed directrix".into(),
                ));
            }
            let out = resolve(graph, *basis, Some(selected), options, depth + 1)?;
            Ok(if *sense_agreement { out } else { reversed(out) })
        }
        // A composite may turn a corner at a joint, which no refinement
        // removes, so it reports no end tangents.
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { segments })) => {
            let mut out = Vec::new();
            for segment in segments {
                let (mut child, _) = resolve(graph, segment.curve, None, options, depth + 1)?;
                if !segment.same_sense {
                    child.reverse();
                }
                stitch(&mut out, child, options.tolerance().linear())?;
                if out.len() > MAX_POINTS {
                    return Err(GeomError::BudgetExceeded {
                        resource: "directrix points",
                    });
                }
            }
            match range {
                Some((start, end)) => {
                    trim_by_length(&out, start, end, options.tolerance().linear())
                        .map(|points| (points, None))
                }
                None => Ok((out, None)),
            }
        }
        // Offsets at stations along a basis curve, sampled (#241).
        Some(GeometryNode::CurveRelation(relation @ CurveRelation::OffsetByStations { .. })) => {
            crate::station::offset_curve_points(graph, relation, range, options)
                .map(|points| (points, None))
        }
        // A curve placed at a station (#264): its source sampled in its own
        // coordinates (a 2D one in `z = 0`), then moved rigidly.
        Some(GeometryNode::InstanceAtStation(placed)) => {
            let transform = crate::station::placement(graph, id)?.transform;
            let (points, ends) = match graph.get(placed.source) {
                Some(GeometryNode::Curve2(curve)) => sample_curve2(curve, range, options)?,
                _ => resolve(graph, placed.source, range, options, depth + 1).map_err(|error| {
                    match error {
                        GeomError::InvalidInput(detail) if detail.contains("not a 3D curve") => {
                            GeomError::UnsupportedInput {
                                backend: crate::BACKEND_ID,
                                operation: axiolid_contracts::Operation::CurveEvaluation,
                                input: "a 2D curve relation placed at a station, as a sweep \
                                        directrix: only an atomic 2D curve is lifted",
                            }
                        }
                        other => other,
                    }
                })?,
            };
            Ok((
                points
                    .into_iter()
                    .map(|p| transform.transform_point3(p))
                    .collect(),
                ends.map(|ends| ends.map(|v| transform.transform_vector3(v))),
            ))
        }
        Some(GeometryNode::CurveRelation(_)) => Err(unsupported_curve_evaluation()),
        Some(_) => Err(GeomError::InvalidInput(format!(
            "sweep directrix {id:?} is not a 3D curve"
        ))),
        None => Err(GeomError::InvalidInput(format!(
            "directrix {id:?} is outside the graph"
        ))),
    }
}

fn unsupported_curve_evaluation() -> GeomError {
    GeomError::Unsupported {
        backend: axiolid_contracts::BackendId::new("scalar-compile"),
        operation: axiolid_contracts::Operation::CurveEvaluation,
    }
}

fn parameter(
    selectors: &[TrimSelector],
    preference: TrimmingPreference,
    label: &str,
    basis: &axiolid_curve::Curve3,
    sense: bool,
    tolerance: axiolid_core::Tolerance,
) -> GeomResult<Scalar> {
    // A Cartesian selector names a POINT. Some formats can only state a trim
    // that way -- a three-point arc knows its endpoints, not their parameters
    // -- so the point is inverted against the basis rather than refused.
    // Inversion is exact or it refuses; it never projects an off-curve point.
    // An arc length is a parameter-kind selector, resolved against the basis
    // by quadrature (`axiolid_reference::arc_parameter`); its failure is
    // reported, never replaced by another selector.
    let first_parameter = || -> GeomResult<Option<Scalar>> {
        selectors
            .iter()
            .find(|selector| is_parameter_kind(selector))
            .map_or(Ok(None), |selector| {
                parameter_kind(selector, Some(basis), sense, label)
            })
    };
    let selected = match preference {
        TrimmingPreference::Parameter => first_parameter()?,
        TrimmingPreference::Unspecified => match selectors.first() {
            Some(selector) if is_parameter_kind(selector) => {
                parameter_kind(selector, Some(basis), sense, label)?
            }
            _ => invert_first_point(selectors, basis, tolerance),
        },
        TrimmingPreference::Cartesian => match invert_first_point(selectors, basis, tolerance) {
            Some(value) => Some(value),
            None => first_parameter()?,
        },
    };
    selected.filter(|value| value.is_finite()).ok_or_else(|| {
        GeomError::InvalidInput(format!(
            "trimmed directrix {label} needs a parameter selector, or a point \
             selector that lies on the basis curve"
        ))
    })
}

/// Whether a selector names a position by a measure along the curve
/// (a parameter or an arc length) rather than by a point.
fn is_parameter_kind(selector: &TrimSelector) -> bool {
    matches!(
        selector,
        TrimSelector::Parameter(_) | TrimSelector::ArcLength(_)
    )
}

/// A parameter-kind selector read as a basis parameter.
///
/// An arc length is measured from the basis parameter `0` in the trim's
/// sense, to the tolerance `axiolid_reference::arc_parameter` states. It
/// needs an analytic basis to measure along: on a relation basis it is
/// refused by name rather than read as a parameter.
fn parameter_kind(
    selector: &TrimSelector,
    basis: Option<&axiolid_curve::Curve3>,
    sense: bool,
    label: &str,
) -> GeomResult<Option<Scalar>> {
    match selector {
        TrimSelector::Parameter(value) => Ok(Some(*value)),
        TrimSelector::ArcLength(length) => {
            let Some(curve) = basis else {
                return Err(GeomError::InvalidInput(format!(
                    "trimmed directrix {label} is an arc-length selector, but its basis \
                     is a curve relation with no single curve to measure along"
                )));
            };
            if !length.is_finite() {
                return Err(GeomError::InvalidInput(format!(
                    "trimmed directrix {label} arc length must be finite"
                )));
            }
            let signed = if sense { *length } else { -*length };
            axiolid_reference::arc_parameter::parameter_at_arc_length3(curve, 0.0, signed)
                .map(Some)
                .map_err(|error| {
                    GeomError::InvalidInput(format!(
                        "trimmed directrix {label} arc length {length} does not resolve on \
                         its basis: {error}"
                    ))
                })
        }
        _ => Ok(None),
    }
}

/// Parameter selectors only, for a basis with no invertible analytic curve.
fn parameter_only(
    selectors: &[TrimSelector],
    preference: TrimmingPreference,
    label: &str,
) -> GeomResult<Scalar> {
    let selected = match preference {
        TrimmingPreference::Parameter | TrimmingPreference::Unspecified => selectors
            .iter()
            .find(|selector| is_parameter_kind(selector))
            .map_or(Ok(None), |selector| {
                parameter_kind(selector, None, true, label)
            })?,
        // The basis is a relation, so there is no analytic curve to invert a
        // point against. Refusing names that, rather than silently using a
        // parameter the file did not designate as authoritative.
        TrimmingPreference::Cartesian => None,
    };
    selected.filter(|value| value.is_finite()).ok_or_else(|| {
        GeomError::InvalidInput(format!(
            "trimmed directrix {label} needs a finite parameter selector; its \
             basis is a curve relation, so a point selector cannot be inverted"
        ))
    })
}

/// First point selector inverted against the basis, if one resolves.
fn invert_first_point(
    selectors: &[TrimSelector],
    basis: &axiolid_curve::Curve3,
    tolerance: axiolid_core::Tolerance,
) -> Option<Scalar> {
    selectors.iter().find_map(|selector| match selector {
        TrimSelector::Point3(point) => {
            axiolid_reference::curve::invert3(basis, *point, tolerance).ok()
        }
        _ => None,
    })
}

/// Parameter period of a closed conic basis; `None` for anything else.
fn period_of(curve: &axiolid_curve::Curve3) -> Option<Scalar> {
    match curve {
        axiolid_curve::Curve3::Circle(_) | axiolid_curve::Curve3::Ellipse(_) => {
            Some(core::f64::consts::TAU)
        }
        _ => None,
    }
}

/// Relative slack for "at most one period": authoring tools write a quarter
/// bend as `(3pi/2, 2pi + 1e-15)`, and that must stay a quarter, not wrap.
const PERIOD_SLACK: Scalar = 1e-9;

/// Sample a trim of a periodic basis in its UNWRAPPED basis interval.
///
/// The trimmed curve runs from `a` forward (sense) or backward (against) to
/// the next occurrence of `b`, so its basis interval is `[a, a + span]` or
/// `[a - span, a]`, which may extend past the natural `[0, period]`; a conic
/// evaluates the same there. Returned in increasing basis parameter; the
/// caller reverses for `!sense`, as for every other basis.
///
/// A sweep `range` is read in the trimmed curve's own parameter, which is
/// the basis parameter: each end is taken modulo the period into the
/// unwrapped interval, so `(330, 30)`, `(330, 390)` and `(-30, 30)` degrees
/// all name the same sub-arc of a `315 -> 45` trim. An end that lands on no
/// point of the arc is refused, as a range outside any trim is.
fn sample_periodic_trim(
    curve: &axiolid_curve::Curve3,
    period: Scalar,
    a: Scalar,
    b: Scalar,
    sense: bool,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Resolved> {
    let (lo, hi) = periodic_trim_interval(period, a, b, sense, range, options)?;
    let points = axiolid_reference::curve::flatten3(
        curve,
        axiolid_core::Interval { start: lo, end: hi },
        crate::compiler::chord_error(options),
        MAX_FLATTEN_DEPTH,
    )?;
    Ok((points, smooth_ends(curve, lo, hi)))
}

/// The unwrapped basis interval `[lo, hi]` of a periodic trim, narrowed to
/// a sweep `range` when one is given (see [`sample_periodic_trim`]).
fn periodic_trim_interval(
    period: Scalar,
    a: Scalar,
    b: Scalar,
    sense: bool,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<(Scalar, Scalar)> {
    let travelled = if sense { b - a } else { a - b };
    let span = if travelled > 0.0 && travelled <= period * (1.0 + PERIOD_SLACK) {
        travelled
    } else {
        match travelled.rem_euclid(period) {
            // `a` and `b` differ by whole turns: the trim is a full turn.
            0.0 => period,
            wrapped => wrapped,
        }
    };
    let (lo, hi) = if sense { (a, a + span) } else { (a - span, a) };
    match range {
        None => Ok((lo, hi)),
        Some((start, end)) => {
            if !(start.is_finite() && end.is_finite()) {
                return Err(GeomError::InvalidInput(
                    "sweep parameter range must be finite".into(),
                ));
            }
            if start == end {
                return Err(GeomError::Degenerate(
                    "sweep parameter range is empty".into(),
                ));
            }
            let slack = options.tolerance().linear();
            let into_arc = |value: Scalar| -> GeomResult<Scalar> {
                let shifted = lo + (value - lo).rem_euclid(period);
                // Just below `lo`, a value wraps to just below `lo + period`:
                // read it as `lo` rather than refuse a rounding residue.
                let shifted = if shifted > lo + period - slack {
                    lo
                } else {
                    shifted
                };
                if shifted > hi + slack {
                    return Err(GeomError::InvalidInput(
                        "sweep range exceeds trimmed directrix".into(),
                    ));
                }
                Ok(shifted.min(hi))
            };
            let (s, e) = (into_arc(start)?, into_arc(end)?);
            // An end ON the trim's far end maps back to `lo` when the arc is
            // a full turn; a range is never empty by that accident.
            let (s, e) = if s == e { (lo, hi) } else { (s, e) };
            Ok((s.min(e), s.max(e)))
        }
    }
}

fn sample_curve(
    curve: &axiolid_curve::Curve3,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Resolved> {
    let domain = sweep_domain(
        axiolid_reference::curve::domain3(curve),
        matches!(curve, axiolid_curve::Curve3::Line(_)),
        range,
        options,
    )?;
    if elevated::is_elevated(curve) {
        elevated::check_bounded(domain)?;
    }
    let points = axiolid_reference::curve::flatten3(
        curve,
        domain,
        crate::compiler::chord_error(options),
        MAX_FLATTEN_DEPTH,
    )?;
    Ok((points, smooth_ends(curve, domain.start, domain.end)))
}

/// A 2D curve sampled as [`sample_curve`] samples a 3D one, lying in
/// `z = 0` (#264: the source of a curve placed at a station). It reports
/// no end tangents, so a sweep along it is swept as sampled.
fn sample_curve2(
    curve: &axiolid_curve::Curve2,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Resolved> {
    let domain = sweep_domain(
        axiolid_reference::curve::domain2(curve),
        matches!(curve, axiolid_curve::Curve2::Line(_)),
        range,
        options,
    )?;
    let points = axiolid_reference::curve::flatten2(
        curve,
        domain,
        crate::compiler::chord_error(options),
        MAX_FLATTEN_DEPTH,
    )?;
    Ok((
        points
            .into_iter()
            .map(|p| Point3::new(p.x, p.y, 0.0))
            .collect(),
        None,
    ))
}

/// The parameter interval a sweep `range` selects of a curve whose
/// natural domain is `natural`; an `unbounded` curve (a line) takes the
/// range as given.
fn sweep_domain(
    natural: axiolid_core::Interval,
    unbounded: bool,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<axiolid_core::Interval> {
    Ok(match range {
        None => natural,
        Some((start, end)) => {
            if !(start.is_finite() && end.is_finite()) {
                return Err(GeomError::InvalidInput(
                    "sweep parameter range must be finite".into(),
                ));
            }
            if start == end {
                return Err(GeomError::Degenerate(
                    "sweep parameter range is empty".into(),
                ));
            }
            let mut lo = natural.start.min(natural.end);
            let mut hi = natural.start.max(natural.end);
            let rlo = start.min(end);
            let rhi = start.max(end);
            if unbounded {
                lo = rlo;
                hi = rhi;
            }
            let slack = options.tolerance().linear();
            if rlo < lo - slack || rhi > hi + slack {
                return Err(GeomError::InvalidInput(
                    "sweep parameter range falls outside curve domain".into(),
                ));
            }
            axiolid_core::Interval {
                start: rlo.max(lo),
                end: rhi.min(hi),
            }
        }
    })
}

fn stitch(target: &mut Vec<Point3>, mut child: Vec<Point3>, tolerance: Scalar) -> GeomResult<()> {
    if target.is_empty() {
        target.append(&mut child);
        return Ok(());
    }
    let gap = target
        .last()
        .unwrap()
        .distance(*child.first().ok_or_else(|| {
            GeomError::Degenerate("composite directrix segment has no points".into())
        })?);
    if gap > tolerance {
        return Err(GeomError::InvalidInput(format!(
            "composite directrix has a {gap} unit gap"
        )));
    }
    child.remove(0);
    target.append(&mut child);
    Ok(())
}

fn trim_by_length(
    points: &[Point3],
    start: Scalar,
    end: Scalar,
    tolerance: Scalar,
) -> GeomResult<Vec<Point3>> {
    if !(start.is_finite() && end.is_finite()) {
        return Err(GeomError::InvalidInput(
            "composite range must be finite".into(),
        ));
    }
    if start == end {
        return Err(GeomError::Degenerate("composite range is empty".into()));
    }
    let mut cumulative = Vec::with_capacity(points.len());
    cumulative.push(0.0);
    for edge in points.windows(2) {
        let next = cumulative.last().copied().unwrap() + edge[0].distance(edge[1]);
        cumulative.push(next);
    }
    let total = cumulative.last().copied().unwrap_or(0.0);
    let lo = start.min(end);
    let hi = start.max(end);
    if lo < -tolerance || hi > total + tolerance {
        return Err(GeomError::InvalidInput(format!(
            "composite range ({start}, {end}) exceeds length {total}"
        )));
    }
    let lo = lo.max(0.0);
    let hi = hi.min(total);
    let mut out = vec![point_at_length(points, &cumulative, lo)?];
    for (&distance, &point) in cumulative
        .iter()
        .zip(points)
        .skip(1)
        .take(points.len().saturating_sub(2))
    {
        if distance > lo && distance < hi {
            out.push(point);
        }
    }
    out.push(point_at_length(points, &cumulative, hi)?);
    if start > end {
        out.reverse();
    }
    Ok(out)
}

fn point_at_length(points: &[Point3], cumulative: &[Scalar], target: Scalar) -> GeomResult<Point3> {
    for (index, limits) in cumulative.windows(2).enumerate() {
        if target <= limits[1] {
            let span = limits[1] - limits[0];
            if span == 0.0 {
                return Ok(points[index]);
            }
            return Ok(points[index].lerp(points[index + 1], (target - limits[0]) / span));
        }
    }
    points
        .last()
        .copied()
        .ok_or_else(|| GeomError::Degenerate("directrix has no points".into()))
}

/// A directrix an exact swept disk can follow: one straight segment or one
/// circular arc (#223). Returned by [`crate::exact_directrix`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ExactDirectrix {
    /// From the first point to the second.
    Segment(Point3, Point3),
    /// Over an angle span of the circle, start to end.
    Arc(axiolid_curve::Circle3, axiolid_core::Interval),
}

/// Resolve a directrix to one segment or one arc, read with the same
/// parameter conventions as [`points`]; anything with corners, or curved
/// other than a circle, is refused by name through `unsupported`.
pub(crate) fn exact(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
    unsupported: fn(&'static str) -> GeomError,
) -> GeomResult<ExactDirectrix> {
    exact_at(graph, id, range, options, unsupported, 0)
}

fn finite_range((start, end): (Scalar, Scalar)) -> GeomResult<(Scalar, Scalar)> {
    if !(start.is_finite() && end.is_finite()) {
        return Err(GeomError::InvalidInput(
            "sweep parameter range must be finite".into(),
        ));
    }
    if start == end {
        return Err(GeomError::Degenerate(
            "sweep parameter range is empty".into(),
        ));
    }
    Ok((start, end))
}

fn exact_at(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
    unsupported: fn(&'static str) -> GeomError,
    depth: usize,
) -> GeomResult<ExactDirectrix> {
    use axiolid_curve::Curve3;
    use core::f64::consts::TAU;
    if depth > MAX_DEPTH {
        return Err(GeomError::BudgetExceeded {
            resource: "directrix relation depth",
        });
    }
    let slack = options.tolerance().linear();
    let other = "exact swept disk along a directrix other than one segment or one arc";
    match graph.get(id) {
        Some(GeometryNode::Curve3(Curve3::Line(line))) => {
            let range =
                range.ok_or_else(|| unsupported("exact swept disk along an unbounded line"))?;
            let (start, end) = finite_range(range)?;
            Ok(ExactDirectrix::Segment(
                line.origin + line.direction * start,
                line.origin + line.direction * end,
            ))
        }
        Some(GeometryNode::Curve3(Curve3::Polyline(polyline)))
            if polyline.points.len() == 2 && !polyline.closed && range.is_none() =>
        {
            Ok(ExactDirectrix::Segment(
                polyline.points[0],
                polyline.points[1],
            ))
        }
        Some(GeometryNode::Curve3(Curve3::Circle(circle))) => {
            let span = match range {
                None => axiolid_core::Interval::new(0.0, TAU),
                Some(range) => {
                    // As the mesh path reads it: an increasing span inside
                    // the circle's natural domain `[0, 2 pi]`.
                    let (start, end) = finite_range(range)?;
                    let (lo, hi) = (start.min(end), start.max(end));
                    if lo < -slack || hi > TAU + slack {
                        return Err(GeomError::InvalidInput(
                            "sweep parameter range falls outside curve domain".into(),
                        ));
                    }
                    axiolid_core::Interval::new(lo.max(0.0), hi.min(TAU))
                }
            };
            Ok(ExactDirectrix::Arc(*circle, span))
        }
        // A straight plan under a constant grade is a segment (#252).
        Some(GeometryNode::Curve3(curve)) => {
            elevated::exact_segment(curve, range, options, unsupported)
                .unwrap_or_else(|| Err(unsupported(other)))
        }
        Some(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve {
            curve_3d,
            master: MasterRepresentation::Curve3d,
            ..
        })) => exact_at(graph, *curve_3d, range, options, unsupported, depth + 1),
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { segments }))
            if segments.len() == 1 && range.is_none() =>
        {
            // One segment has no corners, and its sense does not change
            // the swept solid.
            exact_at(
                graph,
                segments[0].curve,
                None,
                options,
                unsupported,
                depth + 1,
            )
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement,
            preference,
        })) => {
            let Some(GeometryNode::Curve3(curve)) = graph.get(*basis) else {
                return Err(unsupported(
                    "exact swept disk along a trim of a curve relation",
                ));
            };
            let a = parameter(
                start,
                *preference,
                "start",
                curve,
                *sense_agreement,
                options.tolerance(),
            )?;
            let b = parameter(
                end,
                *preference,
                "end",
                curve,
                *sense_agreement,
                options.tolerance(),
            )?;
            if a == b {
                return Err(GeomError::Degenerate(
                    "trimmed directrix has an empty interval".into(),
                ));
            }
            match curve {
                Curve3::Line(line) => {
                    let (lo, hi) = (a.min(b), a.max(b));
                    let (s, e) = match range {
                        None => (a, b),
                        Some(range) => {
                            let (s, e) = finite_range(range)?;
                            if s.min(e) < lo - slack || s.max(e) > hi + slack {
                                return Err(GeomError::InvalidInput(
                                    "sweep range exceeds trimmed directrix".into(),
                                ));
                            }
                            (s.clamp(lo, hi), e.clamp(lo, hi))
                        }
                    };
                    Ok(ExactDirectrix::Segment(
                        line.origin + line.direction * s,
                        line.origin + line.direction * e,
                    ))
                }
                Curve3::Circle(circle) => {
                    let (lo, hi) =
                        periodic_trim_interval(TAU, a, b, *sense_agreement, range, options)?;
                    Ok(ExactDirectrix::Arc(
                        *circle,
                        axiolid_core::Interval::new(lo, hi),
                    ))
                }
                _ => Err(unsupported(
                    "exact swept disk along a trim of a curve other than a line or circle",
                )),
            }
        }
        Some(GeometryNode::CurveRelation(_)) => Err(unsupported(other)),
        Some(GeometryNode::InstanceAtStation(_)) => Err(unsupported(
            "exact swept disk along a curve placed at a station",
        )),
        Some(_) => Err(GeomError::InvalidInput(format!(
            "sweep directrix {id:?} is not a 3D curve"
        ))),
        None => Err(GeomError::InvalidInput(format!(
            "directrix {id:?} is outside the graph"
        ))),
    }
}

/// What a sweep's directrix is, for its deviation report (#232).
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum DirectrixKind {
    /// One straight segment.
    Segment,
    /// One circular arc.
    Arc(axiolid_curve::Circle3),
    /// One smooth curve over a parameter span, as the mesh path samples it:
    /// an ellipse or a B-spline.
    Smooth(axiolid_curve::Curve3, axiolid_core::Interval),
    /// Anything else, named.
    Other(&'static str),
}

/// Classify a directrix the way [`sampled`] reads it.
///
/// A segment or an arc is what [`exact`] resolves. An ellipse or a
/// B-spline, alone, trimmed or read through a curve-3D surface curve, is
/// [`DirectrixKind::Smooth`] over the span [`sampled`] flattens; a trim of
/// a periodic ellipse uses the same unwrapped interval. Everything else
/// (polylines, composites, relations of relations) is named.
pub(crate) fn kind(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<DirectrixKind> {
    fn marker(_: &'static str) -> GeomError {
        GeomError::Degenerate("not one segment or one arc".to_owned())
    }
    if let Ok(exact) = exact(graph, id, range, options, marker) {
        return Ok(match exact {
            ExactDirectrix::Segment(..) => DirectrixKind::Segment,
            ExactDirectrix::Arc(circle, _) => DirectrixKind::Arc(circle),
        });
    }
    smooth_kind(graph, id, range, options, 0)
}

fn smooth_kind(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
    depth: usize,
) -> GeomResult<DirectrixKind> {
    use axiolid_curve::Curve3;
    if depth > MAX_DEPTH {
        return Err(GeomError::BudgetExceeded {
            resource: "directrix relation depth",
        });
    }
    let smooth = |curve: &Curve3| {
        matches!(curve, Curve3::Ellipse(_) | Curve3::BSpline(_)) || elevated::is_elevated(curve)
    };
    // An elevated curve is one smooth curve only where its grade does not
    // break (#252).
    let classified = |curve: &Curve3, span: axiolid_core::Interval| {
        if elevated::is_elevated(curve) && elevated::has_corner(curve, span.start, span.end) {
            DirectrixKind::Other("elevated directrix with a grade break")
        } else {
            DirectrixKind::Smooth(curve.clone(), span)
        }
    };
    Ok(match graph.get(id) {
        Some(GeometryNode::Curve3(curve)) if smooth(curve) => {
            classified(curve, clamped_span(curve, range)?)
        }
        Some(GeometryNode::Curve3(Curve3::Polyline(_))) => {
            DirectrixKind::Other("polyline directrix")
        }
        Some(GeometryNode::Curve3(_)) => DirectrixKind::Other("directrix curve family"),
        Some(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve {
            curve_3d,
            master: MasterRepresentation::Curve3d,
            ..
        })) => smooth_kind(graph, *curve_3d, range, options, depth + 1)?,
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement,
            preference,
        })) => {
            let Some(GeometryNode::Curve3(curve)) = graph.get(*basis) else {
                return Ok(DirectrixKind::Other("trim of a curve relation"));
            };
            if !smooth(curve) {
                return Ok(DirectrixKind::Other("trimmed directrix curve family"));
            }
            let a = parameter(
                start,
                *preference,
                "start",
                curve,
                *sense_agreement,
                options.tolerance(),
            )?;
            let b = parameter(
                end,
                *preference,
                "end",
                curve,
                *sense_agreement,
                options.tolerance(),
            )?;
            if let Some(period) = period_of(curve) {
                let (lo, hi) =
                    periodic_trim_interval(period, a, b, *sense_agreement, range, options)?;
                DirectrixKind::Smooth(curve.clone(), axiolid_core::Interval::new(lo, hi))
            } else {
                let (s, e) = range.unwrap_or((a, b));
                classified(curve, clamped_span(curve, Some((s, e)))?)
            }
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { .. })) => {
            DirectrixKind::Other("composite directrix")
        }
        Some(GeometryNode::CurveRelation(_)) => DirectrixKind::Other("directrix relation"),
        Some(GeometryNode::InstanceAtStation(_)) => {
            DirectrixKind::Other("directrix placed at a station")
        }
        _ => DirectrixKind::Other("not a 3D curve"),
    })
}

/// The span [`sample_curve`] flattens: the natural domain, narrowed to a
/// range when one is given.
fn clamped_span(
    curve: &axiolid_curve::Curve3,
    range: Option<(Scalar, Scalar)>,
) -> GeomResult<axiolid_core::Interval> {
    let natural = axiolid_reference::curve::domain3(curve);
    Ok(match range {
        None => natural,
        Some((start, end)) => {
            let (start, end) = finite_range((start, end))?;
            let lo = natural.start.min(natural.end);
            let hi = natural.start.max(natural.end);
            axiolid_core::Interval::new(start.min(end).max(lo), start.max(end).min(hi))
        }
    })
}
