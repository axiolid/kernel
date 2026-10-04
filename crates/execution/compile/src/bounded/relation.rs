//! Curve-relation boundaries of a curve-bounded plane (#255).
//!
//! A boundary may be a composite of segments or a trim of a basis curve,
//! as IFC `IfcCompositeCurve` and `IfcTrimmedCurve` boundaries are. These
//! are read by the sweep directrix reader ([`crate::directrix::points`]),
//! which already resolves composites (each segment's sense, the joint gap
//! rule) and trims (every [`TrimSelector`] kind, periodic seams, arc
//! lengths through `axiolid_reference::arc_parameter`); nothing about that
//! reading is repeated here.
//!
//! The directrix reader takes 3D curves, while a boundary lives in the
//! plane's 2D parameters. So the boundary's relation tree is first copied
//! into a small graph of its own with every 2D curve lifted to `z = 0` and
//! every 2D point selector with it: a line, circle, ellipse, polyline or
//! B-spline evaluates to the same `(u, v)` lifted as it does in 2D, at the
//! same parameter, so the lift changes no point and no parameter. A 2D
//! family with no 3D twin, a relation kind other than a composite or a
//! trim, and a reference to anything but a curve are refused by name while
//! copying.
//!
//! The deviation of such a boundary is that of its leaves, read on the
//! lifted copy exactly as the atomic boundaries are read: straight leaves
//! are exact, certified families are within the chord budget, any other
//! family is unbounded by name. A composite adds its widest joint gap
//! (at most the linear tolerance), which the reader stitches shut.

use std::collections::HashMap;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{BSplineCurve3, Circle3, Curve2, Curve3, Ellipse3, Line3, Polyline3};
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId,
    TrimSelector,
};

use crate::deviation::DeviationBound;

/// Deepest relation nesting copied, as for directrices.
const MAX_DEPTH: usize = 256;

/// The points of a curve-relation boundary, in travel order, in the
/// plane's parameters lifted to `z = 0` (a 3D leaf keeps its own `z`, which
/// the caller checks).
pub(super) fn points(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
) -> GeomResult<Vec<Point3>> {
    let (lifted, root) = lift(graph, id)?;
    crate::directrix::points(&lifted, root, None, options).map_err(|error| match error {
        GeomError::InvalidInput(message) => GeomError::InvalidInput(format!(
            "curve-bounded boundary {id:?}, read as a curve path: {message}"
        )),
        GeomError::Degenerate(message) => GeomError::Degenerate(format!(
            "curve-bounded boundary {id:?}, read as a curve path: {message}"
        )),
        other => other,
    })
}

/// The deviation of a curve-relation boundary, before the plane frame's
/// stretch (see the module notes).
pub(crate) fn bound(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
) -> DeviationBound {
    match lift(graph, id) {
        Ok((lifted, root)) => leaf_bound(&lifted, root, options, 0),
        // The compiler refuses the boundary for the same reason.
        Err(_) => DeviationBound::Unbounded("curve-bounded plane boundary relation"),
    }
}

fn leaf_bound(
    graph: &GeometryGraph,
    id: NodeId,
    options: &ExecutionOptions,
    depth: usize,
) -> DeviationBound {
    if depth > MAX_DEPTH {
        return DeviationBound::Unbounded("curve-bounded plane boundary relation depth");
    }
    match graph.get(id) {
        Some(GeometryNode::Curve3(Curve3::Line(_) | Curve3::Polyline(_))) => {
            DeviationBound::Proven(0.0)
        }
        Some(GeometryNode::Curve3(curve))
            if axiolid_reference::bound::certifies_flattening3(curve) =>
        {
            DeviationBound::Proven(crate::compiler::chord_error(options))
        }
        Some(GeometryNode::Curve3(_)) => {
            DeviationBound::Unbounded("curve-bounded plane boundary family")
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed { basis, .. })) => {
            leaf_bound(graph, *basis, options, depth + 1)
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { segments })) => {
            let worst = segments
                .iter()
                .map(|segment| leaf_bound(graph, segment.curve, options, depth + 1))
                .fold(DeviationBound::Proven(0.0), DeviationBound::worst);
            match (worst, widest_joint(graph, segments, options)) {
                (DeviationBound::Proven(d), Some(gap)) => DeviationBound::Proven(d + gap),
                (DeviationBound::Proven(_), None) => {
                    DeviationBound::Unbounded("curve-bounded plane boundary joints")
                }
                (other, _) => other,
            }
        }
        // `lift` copies nothing else.
        _ => DeviationBound::Unbounded("curve-bounded plane boundary relation"),
    }
}

/// The widest gap between one composite segment's end and the next one's
/// start, each read in its travel direction: the reader stitches such a
/// gap shut, so the mesh's edge there lies up to that far from the curve.
fn widest_joint(
    graph: &GeometryGraph,
    segments: &[CurveSegment],
    options: &ExecutionOptions,
) -> Option<Scalar> {
    let mut ends = Vec::with_capacity(segments.len());
    for segment in segments {
        let points = crate::directrix::points(graph, segment.curve, None, options).ok()?;
        let (first, last) = (*points.first()?, *points.last()?);
        ends.push(if segment.same_sense {
            (first, last)
        } else {
            (last, first)
        });
    }
    Some(
        ends.windows(2)
            .map(|pair| pair[0].1.distance(pair[1].0))
            .fold(0.0, Scalar::max),
    )
}

/// Copy the boundary's relation tree into a graph of its own, every 2D
/// curve and point lifted to `z = 0`; returns that graph and its root.
fn lift(graph: &GeometryGraph, id: NodeId) -> GeomResult<(GeometryGraph, NodeId)> {
    let mut builder = GeometryGraphBuilder::new();
    let mut copied = HashMap::new();
    let root = lift_node(graph, id, id, &mut builder, &mut copied, 0)?;
    let lifted = builder.finish(vec![root]).map_err(|error| {
        GeomError::InvalidInput(format!(
            "curve-bounded boundary {id:?} cannot be read as a curve path: {error}"
        ))
    })?;
    Ok((lifted, root))
}

fn lift_node(
    graph: &GeometryGraph,
    boundary: NodeId,
    id: NodeId,
    builder: &mut GeometryGraphBuilder,
    copied: &mut HashMap<NodeId, NodeId>,
    depth: usize,
) -> GeomResult<NodeId> {
    if depth > MAX_DEPTH {
        return Err(GeomError::BudgetExceeded {
            resource: "curve-bounded boundary relation depth",
        });
    }
    if let Some(&done) = copied.get(&id) {
        return Ok(done);
    }
    let node = match graph.get(id) {
        Some(GeometryNode::Curve2(curve)) => GeometryNode::Curve3(lift_curve(curve)?),
        Some(GeometryNode::Curve3(curve)) => GeometryNode::Curve3(curve.clone()),
        Some(GeometryNode::CurveRelation(CurveRelation::Composite { segments })) => {
            let mut lifted = Vec::with_capacity(segments.len());
            for segment in segments {
                lifted.push(CurveSegment {
                    curve: lift_node(graph, boundary, segment.curve, builder, copied, depth + 1)?,
                    ..*segment
                });
            }
            GeometryNode::CurveRelation(CurveRelation::Composite { segments: lifted })
        }
        Some(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start,
            end,
            sense_agreement,
            preference,
        })) => GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis: lift_node(graph, boundary, *basis, builder, copied, depth + 1)?,
            start: start.iter().map(lift_selector).collect(),
            end: end.iter().map(lift_selector).collect(),
            sense_agreement: *sense_agreement,
            preference: *preference,
        }),
        Some(GeometryNode::CurveRelation(relation)) => {
            return Err(GeomError::UnsupportedInput {
                backend: crate::BACKEND_ID,
                operation: Operation::Tessellation,
                input: relation_name(relation),
            })
        }
        Some(_) => {
            return Err(GeomError::InvalidInput(format!(
                "curve-bounded boundary {boundary:?} refers to {id:?}, which is not a curve node"
            )))
        }
        None => {
            return Err(GeomError::InvalidInput(format!(
                "curve-bounded boundary {boundary:?} refers to {id:?}, which does not belong \
                 to this graph"
            )))
        }
    };
    let pushed = builder.push(node).map_err(|error| {
        GeomError::InvalidInput(format!(
            "curve-bounded boundary {boundary:?} cannot be read as a curve path: {error}"
        ))
    })?;
    copied.insert(id, pushed);
    Ok(pushed)
}

/// The refusal's name for a relation kind a boundary cannot be.
fn relation_name(relation: &CurveRelation) -> &'static str {
    match relation {
        CurveRelation::Offset { .. } => "curve-bounded plane boundary: an offset curve relation",
        CurveRelation::SurfaceCurve { .. } => {
            "curve-bounded plane boundary: a surface curve relation"
        }
        CurveRelation::ParameterCurve { .. } => {
            "curve-bounded plane boundary: a parameter curve relation"
        }
        CurveRelation::OffsetByStations { .. } => {
            "curve-bounded plane boundary: a curve offset at stations"
        }
        _ => "curve-bounded plane boundary: this curve relation kind",
    }
}

fn lift_point(p: Point2) -> Point3 {
    Point3::new(p.x, p.y, 0.0)
}

fn lift_vector(v: Vec2) -> Vec3 {
    Vec3::new(v.x, v.y, 0.0)
}

fn lift_frame(frame: &Frame2) -> Frame3 {
    let (x, y) = (lift_vector(frame.x), lift_vector(frame.y));
    Frame3 {
        origin: lift_point(frame.origin),
        x,
        y,
        z: x.cross(y),
    }
}

fn lift_selector(selector: &TrimSelector) -> TrimSelector {
    match selector {
        TrimSelector::Point2(p) => TrimSelector::Point3(lift_point(*p)),
        other => *other,
    }
}

/// A 2D curve as the 3D curve in `z = 0` that evaluates to the same points
/// at the same parameters.
fn lift_curve(curve: &Curve2) -> GeomResult<Curve3> {
    Ok(match curve {
        Curve2::Line(line) => Curve3::Line(Line3 {
            origin: lift_point(line.origin),
            direction: lift_vector(line.direction),
        }),
        Curve2::Circle(circle) => Curve3::Circle(Circle3 {
            frame: lift_frame(&circle.frame),
            radius: circle.radius,
        }),
        Curve2::Ellipse(ellipse) => Curve3::Ellipse(Ellipse3 {
            frame: lift_frame(&ellipse.frame),
            semi_axis_x: ellipse.semi_axis_x,
            semi_axis_y: ellipse.semi_axis_y,
        }),
        Curve2::Polyline(polyline) => Curve3::Polyline(Polyline3 {
            points: polyline.points.iter().copied().map(lift_point).collect(),
            closed: polyline.closed,
        }),
        Curve2::BSpline(spline) => Curve3::BSpline(BSplineCurve3 {
            degree: spline.degree,
            control_points: spline
                .control_points
                .iter()
                .copied()
                .map(lift_point)
                .collect(),
            knots: spline.knots.clone(),
            multiplicities: spline.multiplicities.clone(),
            weights: spline.weights.clone(),
            closed: spline.closed,
            self_intersect: spline.self_intersect,
            knot_spec: spline.knot_spec,
        }),
        _ => {
            return Err(GeomError::UnsupportedInput {
                backend: crate::BACKEND_ID,
                operation: Operation::Tessellation,
                input: "curve-bounded plane boundary: a curve relation over a 2D curve family \
                        other than a line, circle, ellipse, polyline or B-spline",
            })
        }
    })
}
