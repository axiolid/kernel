//! Safe resolution of graph-referenced 3D sweep directrices.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point3, Scalar};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryNode, MasterRepresentation, NodeId, TrimSelector,
    TrimmingPreference,
};

const MAX_DEPTH: usize = 256;
const MAX_POINTS: usize = 1_000_000;
const MAX_FLATTEN_DEPTH: u32 = 16;

pub(crate) fn points(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Point3>> {
    let points = resolve(graph, id, range, options, 0)?;
    if points.len() < 2 {
        return Err(GeomError::Degenerate(
            "a sweep directrix needs at least two points".into(),
        ));
    }
    Ok(points)
}

fn resolve(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
    depth: usize,
) -> GeomResult<Vec<Point3>> {
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
                    parameter(start, *preference, "start", curve, options.tolerance())?,
                    parameter(end, *preference, "end", curve, options.tolerance())?,
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
                let mut out =
                    sample_periodic_trim(curve, period, a, b, *sense_agreement, range, options)?;
                if !sense_agreement {
                    out.reverse();
                }
                return Ok(out);
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
            let mut out = resolve(graph, *basis, Some(selected), options, depth + 1)?;
            if !sense_agreement {
                out.reverse();
            }
            Ok(out)
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { segments })) => {
            let mut out = Vec::new();
            for segment in segments {
                let mut child = resolve(graph, segment.curve, None, options, depth + 1)?;
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
                }
                None => Ok(out),
            }
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
    tolerance: axiolid_core::Tolerance,
) -> GeomResult<Scalar> {
    // A Cartesian selector names a POINT. Some formats can only state a trim
    // that way -- a three-point arc knows its endpoints, not their parameters
    // -- so the point is inverted against the basis rather than refused.
    // Inversion is exact or it refuses; it never projects an off-curve point.
    let selected = match preference {
        TrimmingPreference::Parameter => selectors.iter().find_map(as_parameter),
        TrimmingPreference::Unspecified => selectors
            .first()
            .and_then(as_parameter)
            .or_else(|| invert_first_point(selectors, basis, tolerance)),
        TrimmingPreference::Cartesian => invert_first_point(selectors, basis, tolerance)
            .or_else(|| selectors.iter().find_map(as_parameter)),
    };
    selected.filter(|value| value.is_finite()).ok_or_else(|| {
        GeomError::InvalidInput(format!(
            "trimmed directrix {label} needs a parameter selector, or a point \
             selector that lies on the basis curve"
        ))
    })
}

/// Parameter selectors only, for a basis with no invertible analytic curve.
fn parameter_only(
    selectors: &[TrimSelector],
    preference: TrimmingPreference,
    label: &str,
) -> GeomResult<Scalar> {
    let selected = match preference {
        TrimmingPreference::Parameter | TrimmingPreference::Unspecified => {
            selectors.iter().find_map(as_parameter)
        }
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
) -> GeomResult<Vec<Point3>> {
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
    let (lo, hi) = match range {
        None => (lo, hi),
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
            (s.min(e), s.max(e))
        }
    };
    axiolid_reference::curve::flatten3(
        curve,
        axiolid_core::Interval { start: lo, end: hi },
        crate::compiler::chord_error(options),
        MAX_FLATTEN_DEPTH,
    )
}

fn as_parameter(selector: &TrimSelector) -> Option<Scalar> {
    match selector {
        TrimSelector::Parameter(value) => Some(*value),
        _ => None,
    }
}

fn sample_curve(
    curve: &axiolid_curve::Curve3,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Point3>> {
    let natural = axiolid_reference::curve::domain3(curve);
    let domain = match range {
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
            if matches!(curve, axiolid_curve::Curve3::Line(_)) {
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
    };
    axiolid_reference::curve::flatten3(
        curve,
        domain,
        crate::compiler::chord_error(options),
        MAX_FLATTEN_DEPTH,
    )
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
