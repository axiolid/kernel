//! The edge pass: every edge's samples are fixed before any face (#336).
//!
//! Two faces meet along an edge only where both run along it through the
//! same vertices. Before #336 the first face to reach an edge sampled it
//! and later faces took its samples, but a planar face with straight edges
//! took only the edge's endpoints: where the curved face across had split a
//! straight ruling into pieces (its spans are sized by the surface's
//! curvature and turning, #310), the shell had T-junctions and stayed open,
//! and which face sampled an edge depended on the order the shell listed
//! them in.
//!
//! So the samples are fixed here, in one pass over the shells, before a
//! single face is meshed:
//!
//! 1. Every face that takes the curved path states what it needs of each
//!    edge it uses: through a pcurve, the count [`super::edge_sample_count`]
//!    gives (the trim's 3D midpoint test and its certified chord bound);
//!    without one, the count [`super::inversion::sample_edge`] settles on
//!    for this face's surface (the midpoint test, the quarter-turn guard and
//!    the certified span terms). Each count is a power of two over uniform
//!    spans of one parameterisation, so the largest is the strictest: its
//!    samples refine every other face's.
//! 2. Each edge is sampled once at that count. With no pcurve on any use it
//!    is sampled along its 3D trace, in its own direction, and keeps its
//!    span bounds on that curve. When a use has a pcurve, the samples are
//!    that pcurve's image on its face, the use with the smallest pcurve node
//!    taken; every other pcurve use must agree with them within the chord
//!    budget, as before. A planar face states no requirement: one segment
//!    between the vertices is all it needs.
//! 3. The samples are interned as mesh vertices, endpoints welded to the
//!    topology vertices, and stored as the shared-edge cache. Every face
//!    then takes them: a curved face inverts or re-evaluates them, a planar
//!    face (straight edges included) inserts the interior samples into its
//!    rings as corners, so each piece of the edge is used once from each
//!    side.
//!
//! The outcome depends on the faces using an edge, not on their order. An
//! edge only planar faces use is not sampled here; such faces keep their
//! vertex polygons, exactly as before.

use std::collections::BTreeMap;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Vec3;
use axiolid_model::NodeId;
use axiolid_topology::Orientation;

use super::inversion::{self, EdgeTrace, SurfaceInverter};
use super::{
    edge_sample_count, edge_samples, face_surface, takes_curved_path, trim_samples,
    with_curved_face_transaction, EdgeSamples, FaceBound, FaceContext, SpanBounds, WeldedVertices,
};
use axiolid_mesh::TriMesh;

/// What the faces using one edge need of it.
struct EdgePlan<'g> {
    /// The strictest number of uniform spans any face asked for.
    segments: usize,
    /// The pcurve use the samples are taken through, when a use has one.
    pcurve: Option<PcurveSource<'g>>,
}

/// A pcurve use of an edge: its trim, its face's surface and its direction.
struct PcurveSource<'g> {
    node: NodeId,
    trim: &'g axiolid_curve::Curve2,
    surface: &'g axiolid_surface::Surface,
    sense: Orientation,
}

/// Fix and intern the samples of every edge a curved face of `shells`
/// uses, into `cache` (the module notes give the rules).
pub(super) fn sample_shell_edges(
    mesh: &mut TriMesh,
    ctx: &FaceContext<'_>,
    shells: &[(&axiolid_topology::Shell, bool)],
    welded: &mut WeldedVertices,
    cache: &mut EdgeSamples,
    total_curved_records: &mut usize,
) -> GeomResult<()> {
    let plans = plan_edges(ctx, shells)?;
    with_curved_face_transaction(
        mesh,
        welded,
        cache,
        total_curved_records,
        ctx.chord_error,
        |state| {
            for (index, plan) in &plans {
                let edge = ctx
                    .brep
                    .edges()
                    .get(*index)
                    .ok_or_else(|| GeomError::InvalidInput("edge missing".to_string()))?;
                let id = ctx
                    .brep
                    .edge_id_at(*index)
                    .ok_or_else(|| GeomError::InvalidInput("edge missing".to_string()))?;
                let (points, sense, spans) = planned_samples(ctx, edge, plan)?;
                let (start, end) = if sense == Orientation::Forward {
                    (edge.start, edge.end)
                } else {
                    (edge.end, edge.start)
                };
                edge_samples(state, id, sense, start, end, &points, spans)?;
            }
            Ok(())
        },
    )
}

/// The strictest requirement of every curved face on each edge it uses,
/// keyed by edge index so the edges are interned in a fixed order.
fn plan_edges<'g>(
    ctx: &FaceContext<'g>,
    shells: &[(&axiolid_topology::Shell, bool)],
) -> GeomResult<BTreeMap<usize, EdgePlan<'g>>> {
    let (brep, graph) = (ctx.brep, ctx.graph);
    let mut plans: BTreeMap<usize, EdgePlan<'g>> = BTreeMap::new();
    for &(face_id, _) in shells.iter().flat_map(|(shell, _)| shell.faces.iter()) {
        let face = brep
            .faces()
            .get(face_id.index())
            .ok_or_else(|| GeomError::InvalidInput("face missing".to_string()))?;
        let Some(surface) = face_surface(graph, face)? else {
            continue;
        };
        if !takes_curved_path(brep, graph, face, surface) {
            continue;
        }
        let bound = FaceBound::new(surface);
        let mut inverter: Option<SurfaceInverter<'_>> = None;
        for face_bound in &face.bounds {
            let wire = brep
                .loops()
                .get(face_bound.loop_id.index())
                .ok_or_else(|| GeomError::InvalidInput("loop missing".to_string()))?;
            for use_ in &wire.edges {
                let edge = brep
                    .edges()
                    .get(use_.edge.index())
                    .ok_or_else(|| GeomError::InvalidInput("edge missing".to_string()))?;
                let (segments, source) = match use_.pcurve {
                    Some(node) => {
                        let Some(axiolid_model::GeometryNode::Curve2(trim)) = graph.get(node)
                        else {
                            return Err(GeomError::InvalidInput(
                                "edge pcurve must reference a Curve2 node".to_string(),
                            ));
                        };
                        let segments = edge_sample_count(trim, surface, ctx.chord_error)?;
                        let source = PcurveSource {
                            node,
                            trim,
                            surface,
                            sense: use_.orientation,
                        };
                        (segments, Some(source))
                    }
                    None => {
                        let inverter = match inverter.as_mut() {
                            Some(inverter) => inverter,
                            None => inverter.insert(SurfaceInverter::new(surface, ctx.tolerance)?),
                        };
                        let trace = inversion::edge_trace(ctx, edge)?;
                        let (points, _, _) = inversion::sample_edge(&trace, inverter, &bound, ctx)?;
                        (points.len().saturating_sub(1).max(1), None)
                    }
                };
                let plan = plans.entry(use_.edge.index()).or_insert(EdgePlan {
                    segments,
                    pcurve: None,
                });
                plan.segments = plan.segments.max(segments);
                if let Some(source) = source {
                    if plan
                        .pcurve
                        .as_ref()
                        .is_none_or(|kept| source.node < kept.node)
                    {
                        plan.pcurve = Some(source);
                    }
                }
            }
        }
    }
    Ok(plans)
}

/// One planned edge's samples, the direction they run in and their span
/// bounds on the edge's 3D curve.
fn planned_samples(
    ctx: &FaceContext<'_>,
    edge: &axiolid_topology::Edge<NodeId>,
    plan: &EdgePlan<'_>,
) -> GeomResult<(Vec<Vec3>, Orientation, SpanBounds)> {
    let n = plan.segments;
    if let Some(source) = &plan.pcurve {
        let points = trim_samples(source.trim, n)?
            .iter()
            .map(|p| axiolid_reference::surface::evaluate(source.surface, p.x, p.y))
            .collect::<GeomResult<Vec<_>>>()?;
        return Ok((points, source.sense, Err(inversion::THROUGH_PCURVE)));
    }
    let trace: EdgeTrace<'_> = inversion::edge_trace(ctx, edge)?;
    let points = (0..=n)
        .map(|i| trace.at(i as axiolid_core::Scalar / n as axiolid_core::Scalar))
        .collect::<GeomResult<Vec<_>>>()?;
    let spans = inversion::spans_of(&trace, &points).ok_or(inversion::NO_CHORD_BOUND);
    Ok((points, Orientation::Forward, spans))
}
