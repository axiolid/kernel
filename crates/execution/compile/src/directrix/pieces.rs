//! A swept disk's directrix read as exact segments and arcs (#232).
//!
//! Lines, polylines, circles, their trims and composites of them resolve
//! to [`PathPiece`]s with the same parameter conventions as the sampled
//! reading in the parent module, so `axiolid_construct::pipe` can bound
//! the tube against its exact centreline and mitre its corners. Any
//! other family (an ellipse, a B-spline) yields `None`, and the caller
//! sweeps the sampled path as before.

use axiolid_construct::pipe::PathPiece;
use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Interval, Point3, Scalar};
use axiolid_curve::{Circle3, Curve3};
use axiolid_model::{CurveRelation, GeometryGraph, GeometryNode, MasterRepresentation, NodeId};

use super::{
    parameter, parameter_only, period_of, periodic_trim_interval, unsupported_curve_evaluation,
    MAX_DEPTH,
};

/// The directrix as exact pieces in the direction of travel, or `None`
/// when it holds a curve other than a line, a polyline or a circle, or is
/// one circle or arc alone (that keeps the conic sweep of #231).
pub(crate) fn pieces(
    graph: &GeometryGraph,
    id: NodeId,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Option<Vec<PathPiece>>> {
    let pieces = resolve(graph, id, range, options, 0)?;
    Ok(pieces.filter(|p| !matches!(p.as_slice(), [PathPiece::Arc { .. }])))
}

/// The directrix as exact pieces in the direction of travel, a lone arc
/// included, or `None` when it holds a curve other than a line, a
/// polyline or a circle (#263: its joins are seams).
pub(crate) fn exact_pieces(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
) -> GeomResult<Option<Vec<PathPiece>>> {
    resolve(graph, id, None, options, 0)
}

type Resolved = Option<Vec<PathPiece>>;

fn reversed(pieces: Resolved) -> Resolved {
    pieces.map(|p| p.iter().rev().map(PathPiece::reversed).collect())
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
        Some(GeometryNode::Curve3(curve)) => curve_pieces(curve, range, options),
        Some(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve {
            curve_3d,
            master: MasterRepresentation::Curve3d,
            ..
        })) => resolve(graph, *curve_3d, range, options, depth + 1),
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
            if let Some((curve, period)) = basis_curve.and_then(|c| period_of(c).map(|p| (c, p))) {
                let Curve3::Circle(circle) = curve else {
                    return Ok(None);
                };
                let (lo, hi) =
                    periodic_trim_interval(period, a, b, *sense_agreement, range, options)?;
                let out = arc(circle, lo, hi).map(|p| vec![p]);
                return Ok(if *sense_agreement { out } else { reversed(out) });
            }
            let selected = range.unwrap_or((a, b));
            let slack = options.tolerance().linear();
            if selected.0.min(selected.1) < a.min(b) - slack
                || selected.0.max(selected.1) > a.max(b) + slack
            {
                return Err(GeomError::InvalidInput(
                    "sweep range exceeds trimmed directrix".into(),
                ));
            }
            let out = resolve(graph, *basis, Some(selected), options, depth + 1)?;
            Ok(if *sense_agreement { out } else { reversed(out) })
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { segments })) => {
            let mut out: Vec<PathPiece> = Vec::new();
            let slack = options.tolerance().linear();
            for segment in segments {
                let Some(child) = resolve(graph, segment.curve, None, options, depth + 1)? else {
                    return Ok(None);
                };
                let child = if segment.same_sense {
                    Some(child)
                } else {
                    reversed(Some(child))
                }
                .unwrap_or_default();
                if let (Some(last), Some(first)) = (out.last(), child.first()) {
                    let gap = last.end_point().distance(first.start_point());
                    if gap > slack {
                        return Err(GeomError::InvalidInput(format!(
                            "composite directrix has a {gap} unit gap"
                        )));
                    }
                }
                out.extend(child);
            }
            match range {
                None => Ok(Some(out)),
                Some((start, end)) => trim_by_length(&out, start, end, slack).map(Some),
            }
        }
        // No exact pieces: swept as sampled (#241).
        Some(GeometryNode::CurveRelation(CurveRelation::OffsetByStations { .. })) => Ok(None),
        Some(GeometryNode::CurveRelation(_)) => Err(unsupported_curve_evaluation()),
        Some(_) => Err(GeomError::InvalidInput(format!(
            "sweep directrix {id:?} is not a 3D curve"
        ))),
        None => Err(GeomError::InvalidInput(format!(
            "directrix {id:?} is outside the graph"
        ))),
    }
}

/// The interval of `curve` a sweep `range` selects, read as the sampled
/// path reads it: sorted, inside the natural domain up to the linear
/// tolerance (a line has none), clamped to it.
fn domain(
    curve: &Curve3,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Interval> {
    let natural = axiolid_reference::curve::domain3(curve);
    let Some((start, end)) = range else {
        return Ok(natural);
    };
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
    let (rlo, rhi) = (start.min(end), start.max(end));
    if matches!(curve, Curve3::Line(_)) {
        return Ok(Interval::new(rlo, rhi));
    }
    let (lo, hi) = (
        natural.start.min(natural.end),
        natural.start.max(natural.end),
    );
    let slack = options.tolerance().linear();
    if rlo < lo - slack || rhi > hi + slack {
        return Err(GeomError::InvalidInput(
            "sweep parameter range falls outside curve domain".into(),
        ));
    }
    Ok(Interval::new(rlo.max(lo), rhi.min(hi)))
}

fn curve_pieces(
    curve: &Curve3,
    range: Option<(Scalar, Scalar)>,
    options: &ExecutionOptions,
) -> GeomResult<Resolved> {
    match curve {
        Curve3::Line(line) => {
            let d = domain(curve, range, options)?;
            Ok(Some(vec![PathPiece::Segment {
                start: line.origin + line.direction * d.start,
                end: line.origin + line.direction * d.end,
            }]))
        }
        Curve3::Polyline(polyline) => {
            // Its last section would have to be mitred back onto its first,
            // and the frames carried round a loop that leaves its plane do
            // not return to themselves; not built yet (#245).
            if polyline.closed {
                return Err(GeomError::UnsupportedInput {
                    backend: crate::BACKEND_ID,
                    operation: axiolid_contracts::Operation::Sweep,
                    input: "a swept disk along a closed polyline (the mitre where it closes)",
                });
            }
            let d = domain(curve, range, options)?;
            let points = &polyline.points;
            let requested = d.end - d.start;
            let segments = points.len().saturating_sub(1);
            // As `flatten3`: a normalised range on a longer polyline would
            // silently drop its vertices.
            if segments > 1 && requested <= 1.0 {
                return Err(GeomError::InvalidInput(format!(
                    "polyline domain {d:?} spans {requested} of {segments} segments; a \
                     polyline parameter is one unit per segment"
                )));
            }
            let mut out = Vec::with_capacity(segments);
            for k in 0..segments {
                let (lo, hi) = (d.start.max(k as Scalar), d.end.min((k + 1) as Scalar));
                if hi <= lo {
                    continue;
                }
                let (p, q) = (points[k], points[k + 1]);
                let (start, end) = (
                    p + (q - p) * (lo - k as Scalar),
                    p + (q - p) * (hi - k as Scalar),
                );
                // A repeated vertex has no direction to sweep along.
                if start != end {
                    out.push(PathPiece::Segment { start, end });
                }
            }
            if out.is_empty() {
                return Err(GeomError::Degenerate(
                    "a sweep directrix needs at least two points".into(),
                ));
            }
            Ok(Some(out))
        }
        Curve3::Circle(circle) => {
            let d = domain(curve, range, options)?;
            Ok(arc(circle, d.start, d.end).map(|p| vec![p]))
        }
        _ => Ok(None),
    }
}

/// The arc of `circle` from parameter `lo` to `hi`, or `None` when its
/// frame is not orthonormal (then it is an ellipse in all but name).
fn arc(circle: &Circle3, lo: Scalar, hi: Scalar) -> Option<PathPiece> {
    let (x, y) = (circle.frame.x, circle.frame.y);
    let unit = |v: axiolid_core::Vec3| (v.length() - 1.0).abs() <= 1e-9;
    if !(unit(x) && unit(y) && x.dot(y).abs() <= 1e-9 && circle.radius > 0.0) {
        return None;
    }
    let point: Point3 = circle.frame.origin + (x * lo.cos() + y * lo.sin()) * circle.radius;
    Some(PathPiece::Arc {
        centre: circle.frame.origin,
        axis: x.cross(y).normalize(),
        start: point,
        angle: hi - lo,
    })
}

/// The pieces between arc lengths `start` and `end`, walked backwards
/// when `start > end`: the sampled reading's length trim, measured on the
/// exact pieces.
fn trim_by_length(
    pieces: &[PathPiece],
    start: Scalar,
    end: Scalar,
    tolerance: Scalar,
) -> GeomResult<Vec<PathPiece>> {
    if !(start.is_finite() && end.is_finite()) {
        return Err(GeomError::InvalidInput(
            "composite range must be finite".into(),
        ));
    }
    if start == end {
        return Err(GeomError::Degenerate("composite range is empty".into()));
    }
    let total: Scalar = pieces.iter().map(PathPiece::length).sum();
    let (lo, hi) = (start.min(end), start.max(end));
    if lo < -tolerance || hi > total + tolerance {
        return Err(GeomError::InvalidInput(format!(
            "composite range ({start}, {end}) exceeds length {total}"
        )));
    }
    let (lo, hi) = (lo.max(0.0), hi.min(total));
    let mut out = Vec::new();
    let mut at = 0.0;
    for piece in pieces {
        let length = piece.length();
        let (a, b) = (lo.max(at), hi.min(at + length));
        if b > a && length > 0.0 {
            out.push(piece.between((a - at) / length, (b - at) / length));
        }
        at += length;
    }
    if start > end {
        out = out.iter().rev().map(PathPiece::reversed).collect();
    }
    Ok(out)
}
