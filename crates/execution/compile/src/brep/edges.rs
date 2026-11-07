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
//!    samples refine every other face's. An edge with no pcurve on any use
//!    is cut at its 3D curve's derivative breaks (#342), a polyline's
//!    corners for one, and its requirement is one count per piece between
//!    them, the largest per piece, their sum within the segment cap; an edge
//!    some use runs through a pcurve is one uniform count over the whole
//!    edge, its breaks ignored.
//! 2. Each edge is sampled once at that count. With no pcurve on any use it
//!    is sampled along its 3D trace, in its own direction, and keeps its
//!    span bounds on that curve. When a use has a pcurve, the samples are
//!    that pcurve's image on its face, the use with the smallest pcurve node
//!    taken; every other pcurve use has them located on its own pcurve
//!    ([`super::locate`], #413), within the chord budget, and the samples are
//!    halved until the located spans of every such use meet the budget too,
//!    up to the segment cap. Along its 3D trace each break is a sample, bit for
//!    bit, so a straight piece is one exact span. A planar face states no requirement: one segment
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

use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Vec3;
use axiolid_topology::Orientation;

use super::inversion::{self, EdgeTrace, SurfaceInverter};
use super::locate::{locate_samples, Limits};
use super::view::Pcurve;
use super::{
    edge_sample_count, edge_samples, span_fits, takes_curved_path, trim_samples,
    with_curved_face_transaction, EdgeSamples, FaceBound, FaceContext, SpanBounds, WeldedVertices,
    MAX_CURVED_EDGE_SEGMENTS,
};
use axiolid_mesh::TriMesh;
use axiolid_reference::bound::SurfaceBoundOracle;

/// What the faces using one edge need of it.
struct EdgePlan<'g> {
    /// The strictest sampling any face asked for.
    segments: Requirement,
    /// The pcurve use the samples are taken through, when a use has one.
    pcurve: Option<PcurveSource<'g>>,
    /// Every pcurve use of the edge, the one sampled included: the samples
    /// are located on each (#413).
    uses: Vec<PcurveSource<'g>>,
}

/// How finely an edge is sampled.
enum Requirement {
    /// Uniform spans of the whole edge, the largest count asked for: an edge
    /// some use runs through a pcurve, whose samples are that pcurve's.
    Uniform(usize),
    /// Per piece between the 3D curve's derivative breaks (#342), uniform
    /// spans inside it, the largest count asked for: an edge with no pcurve
    /// on any use, sampled along its 3D trace.
    Pieces(Vec<usize>),
}

/// A pcurve use of an edge: its trim, its face's surface and its direction.
#[derive(Clone, Copy)]
struct PcurveSource<'g> {
    pcurve: Pcurve<'g>,
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
        (ctx.chord_error, ctx.endpoint_limit()),
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
    let brep = ctx.brep;
    let mut faces: Vec<(
        &axiolid_surface::Surface,
        Vec<&axiolid_topology::EdgeUse<Pcurve<'g>>>,
    )> = Vec::new();
    for &(face_id, _) in shells.iter().flat_map(|(shell, _)| shell.faces.iter()) {
        let face = brep
            .faces()
            .get(face_id.index())
            .ok_or_else(|| GeomError::InvalidInput("face missing".to_string()))?;
        let Some(surface) = face.surface else {
            continue;
        };
        if !takes_curved_path(brep, face, surface) {
            continue;
        }
        let mut uses: Vec<&axiolid_topology::EdgeUse<Pcurve<'g>>> = Vec::new();
        for face_bound in &face.bounds {
            let wire = brep
                .loops()
                .get(face_bound.loop_id.index())
                .ok_or_else(|| GeomError::InvalidInput("loop missing".to_string()))?;
            uses.extend(wire.edges.iter());
        }
        faces.push((surface, uses));
    }
    // An edge some use runs through a pcurve is sampled through it, whole.
    let through_pcurve: BTreeSet<usize> = faces
        .iter()
        .flat_map(|(_, uses)| uses.iter())
        .filter(|use_| use_.pcurve.is_some())
        .map(|use_| use_.edge.index())
        .collect();
    let mut plans: BTreeMap<usize, EdgePlan<'g>> = BTreeMap::new();
    for (surface, uses) in faces {
        let bound = FaceBound::new(surface);
        let mut inverter: Option<SurfaceInverter<'_>> = None;
        for use_ in uses {
            let edge = brep
                .edges()
                .get(use_.edge.index())
                .ok_or_else(|| GeomError::InvalidInput("edge missing".to_string()))?;
            let (segments, source) = match use_.pcurve {
                Some(pcurve) => {
                    let segments = edge_sample_count(&pcurve, surface, ctx.chord_error)?;
                    let source = PcurveSource {
                        pcurve,
                        surface,
                        sense: use_.orientation,
                    };
                    (Requirement::Uniform(segments), Some(source))
                }
                None => {
                    let inverter = match inverter.as_mut() {
                        Some(inverter) => inverter,
                        None => inverter.insert(SurfaceInverter::new(surface, ctx.tolerance)?),
                    };
                    let trace = inversion::edge_trace(ctx, edge)?;
                    let uniform = through_pcurve.contains(&use_.edge.index());
                    let breaks = if uniform { Vec::new() } else { trace.breaks() };
                    let (_, _, _, counts) =
                        inversion::sample_edge(&trace, inverter, &bound, ctx, &breaks)?;
                    if uniform {
                        (Requirement::Uniform(counts[0]), None)
                    } else {
                        (Requirement::Pieces(counts), None)
                    }
                }
            };
            match plans.entry(use_.edge.index()) {
                Entry::Vacant(slot) => {
                    slot.insert(EdgePlan {
                        segments,
                        pcurve: source,
                        uses: source.into_iter().collect(),
                    });
                }
                Entry::Occupied(slot) => {
                    let plan = slot.into_mut();
                    plan.segments.merge(segments)?;
                    if let Some(source) = source {
                        plan.uses.push(source);
                        if plan
                            .pcurve
                            .as_ref()
                            .is_none_or(|kept| source.pcurve.key < kept.pcurve.key)
                        {
                            plan.pcurve = Some(source);
                        }
                    }
                }
            }
        }
    }
    Ok(plans)
}

impl Requirement {
    /// The stricter of two requirements on one edge: the larger count, per
    /// piece. Counts are powers of two, so the larger refines the smaller.
    fn merge(&mut self, other: Self) -> GeomResult<()> {
        match (self, other) {
            (Self::Uniform(kept), Self::Uniform(n)) => *kept = (*kept).max(n),
            (Self::Pieces(kept), Self::Pieces(counts)) if kept.len() == counts.len() => {
                for (kept, n) in kept.iter_mut().zip(counts) {
                    *kept = (*kept).max(n);
                }
                // The largest per piece may exceed the cap that each use met.
                if kept.iter().sum::<usize>() > MAX_CURVED_EDGE_SEGMENTS {
                    return Err(GeomError::BudgetExceeded {
                        resource: "curved edge samples",
                    });
                }
            }
            _ => {
                return Err(GeomError::Degenerate(
                    "the faces using an edge disagree on how it is cut into pieces".to_owned(),
                ))
            }
        }
        Ok(())
    }
}

/// One planned edge's samples, the direction they run in and their span
/// bounds on the edge's 3D curve.
fn planned_samples(
    ctx: &FaceContext<'_>,
    edge: &axiolid_topology::Edge<&axiolid_curve::Curve3>,
    plan: &EdgePlan<'_>,
) -> GeomResult<(Vec<Vec3>, Orientation, SpanBounds)> {
    let (Some(source), Requirement::Uniform(n)) = (&plan.pcurve, &plan.segments) else {
        let (None, Requirement::Pieces(counts)) = (&plan.pcurve, &plan.segments) else {
            return Err(GeomError::Degenerate(
                "an edge's sampling plan and its pcurve use disagree".to_owned(),
            ));
        };
        let trace: EdgeTrace<'_> = inversion::edge_trace(ctx, edge)?;
        let fractions = inversion::piece_fractions(&trace.breaks(), counts);
        let points = fractions
            .iter()
            .map(|&s| trace.at(s))
            .collect::<GeomResult<Vec<_>>>()?;
        let spans =
            inversion::spans_of(&trace, &fractions, &points).ok_or(inversion::NO_CHORD_BOUND);
        return Ok((points, Orientation::Forward, spans));
    };
    // Halve the sampling use's span until the other uses' located spans
    // meet the budget too, up to the cap (#413): a located span can be
    // longer than the sampling use's own, and a boundary edge is never
    // split by the face's interior refinement.
    let mut n = *n;
    loop {
        let points = trim_samples(&source.pcurve, n)?
            .iter()
            .map(|p| axiolid_reference::surface::evaluate(source.surface, p.x, p.y))
            .collect::<GeomResult<Vec<_>>>()?;
        if n >= MAX_CURVED_EDGE_SEGMENTS || located_within_budget(ctx, plan, source, &points)? {
            return Ok((points, source.sense, Err(inversion::THROUGH_PCURVE)));
        }
        n *= 2;
    }
}

/// Whether every other pcurve use of the edge, with `points` (in the
/// direction of `source`) located on it, has each span's certified chord
/// bound within the chord budget. A use whose bound cannot be had counts as
/// within: it is reported as it stands.
fn located_within_budget(
    ctx: &FaceContext<'_>,
    plan: &EdgePlan<'_>,
    source: &PcurveSource<'_>,
    points: &[Vec3],
) -> GeomResult<bool> {
    let limits = Limits {
        interior: ctx.chord_error,
        end: ctx.endpoint_limit(),
    };
    for other in &plan.uses {
        if other.pcurve.key == source.pcurve.key
            && other.sense == source.sense
            && std::ptr::eq(other.surface, source.surface)
        {
            continue;
        }
        let Some(oracle) = SurfaceBoundOracle::new(other.surface) else {
            continue;
        };
        let mut targets = points.to_vec();
        if other.sense != source.sense {
            targets.reverse();
        }
        let located = locate_samples(&other.pcurve, other.surface, &targets, limits)?;
        let breaks = axiolid_reference::bound::continuity_breaks2(other.pcurve.curve, 1);
        for span in located.ts.windows(2) {
            let fits = span_fits(
                &oracle,
                (other.pcurve.curve, &breaks),
                span[0],
                span[1],
                ctx.chord_error,
            );
            if fits == Some(false) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
