//! Closed regions of a solid's cut (ADR 0092, #372).
//!
//! The cut is bounded by the cut chains and by the region edges the cut lies
//! across (`request.rs`). Both are directed with the material on the left, so
//! together they form a directed graph in which every node has as many
//! outgoing as incoming edges. Which edge links to which is decided on source
//! topology and exact source signs, never on coordinates:
//!
//! 1. Every cut segment is an edge in its recorded direction, tagged `Cut`;
//!    every bounding region edge is a vertex-to-vertex edge tagged `Region`.
//! 2. A node with unequal degrees breaks the topology. At a vertex with two
//!    or more outgoing edges the boundary meets itself (a pinch, #373). Its
//!    edges are ordered counter-clockwise in the frame by exact signs (the
//!    half-plane on the projected axes, then `orient`), outgoing and incoming
//!    edges alternate, and each incoming edge continues with the next
//!    outgoing edge clockwise, the one that opens its material wedge. Edges
//!    that do not alternate are refused as crossing rings, two edges along
//!    one ray as a pinch.
//! 3. Each closed walk is split wherever it comes back to a vertex, so no
//!    ring repeats a vertex: corners that touch give rings that touch there,
//!    and a hole that touches its outer ring stays a hole that touches it.
//!    Each region is one connected piece of the cut's inside, written as
//!    `axiolid-overlay` writes regions (ADR 0104). A chain that passes a pinch
//!    (ADR 0091) must stay on one ring, or the pinch is refused by name.
//! 4. A node of a ring is dropped when both its edges are cut segments and
//!    they are exactly collinear (`chains::is_droppable`), or when both are
//!    region edges and the region ring itself dropped the vertex. A junction
//!    of a cut segment and a region edge is a vertex of both outputs and
//!    stays.
//! 5. Every node gets its exact point (`plane_points`): a source vertex in the
//!    plane or the exact crossing of a source edge with it. `ring_nesting`
//!    decides on those points whether a ring is an outer ring or a hole and
//!    groups them; the returned coordinates only label the output.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation, Sign};
use axiolid_mesh_section_contract::{
    CutRegion, CutRegionLimits, REFUSED_CROSSING_CUT_RINGS, REFUSED_PINCHED_CUT_REGION,
};

use super::chains::{chain_points, has_adjacent_duplicate, is_droppable, ChainContext};
use super::plane_points::PlanePoints;
use super::ring_nesting::{nest_rings, TracedRing};
use super::{internal_topology_error, NodeKey, ScalarSection, Segment};

/// What an edge of the cut graph is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Link {
    /// A cut segment.
    Cut,
    /// An edge of an in-plane region that the cut lies across.
    Region,
}

/// What ring tracing needs besides the chain context.
pub(super) struct CutInputs<'a> {
    pub(super) chains: &'a ChainContext<'a>,
    pub(super) segments: &'a BTreeSet<Segment>,
    /// Region edges that bound the cut, directed with the material on the left.
    pub(super) region_edges: &'a BTreeSet<(u32, u32)>,
    /// Source vertices the in-plane regions keep on their rings.
    pub(super) region_ring_vertices: &'a BTreeSet<u32>,
}

/// One directed edge of the cut graph, with the material on its left.
#[derive(Debug, Clone, Copy)]
struct CutEdge {
    from: NodeKey,
    to: NodeKey,
    link: Link,
}

fn refused(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: ScalarSection::ID,
        operation: Operation::MeshPlaneSection,
        input,
    }
}

/// The edges of the cut graph: cut segments in their recorded direction, then
/// bounding region edges, both in their set order.
fn cut_edges(inputs: &CutInputs<'_>) -> GeomResult<Vec<CutEdge>> {
    let mut edges = Vec::with_capacity(inputs.segments.len() + inputs.region_edges.len());
    for &segment in inputs.segments {
        let Segment(left, right) = segment;
        let forward = *inputs
            .chains
            .directions
            .get(&segment)
            .ok_or_else(internal_topology_error)?;
        let (from, to) = if forward {
            (left, right)
        } else {
            (right, left)
        };
        edges.push(CutEdge {
            from,
            to,
            link: Link::Cut,
        });
    }
    for &(from, to) in inputs.region_edges {
        edges.push(CutEdge {
            from: NodeKey::Vertex(from),
            to: NodeKey::Vertex(to),
            link: Link::Region,
        });
    }
    Ok(edges)
}

/// Incoming and outgoing edges of one node.
#[derive(Default)]
struct Incident {
    incoming: Vec<usize>,
    outgoing: Vec<usize>,
}

/// The edge that follows each edge, before rings are split at repeated
/// vertices, and the incidences of every node.
fn successors(
    ctx: &ChainContext<'_>,
    edges: &[CutEdge],
    exact: &mut PlanePoints,
    options: &ExecutionOptions,
) -> GeomResult<(Vec<usize>, BTreeMap<NodeKey, Incident>)> {
    let mut nodes = BTreeMap::<NodeKey, Incident>::new();
    for (index, edge) in edges.iter().enumerate() {
        nodes.entry(edge.from).or_default().outgoing.push(index);
        nodes.entry(edge.to).or_default().incoming.push(index);
    }
    let mut next = vec![usize::MAX; edges.len()];
    for (&node, incident) in &nodes {
        options.check_cancelled()?;
        if incident.incoming.len() != incident.outgoing.len() {
            return Err(internal_topology_error());
        }
        match (node, incident.outgoing.as_slice()) {
            (_, &[only]) => next[incident.incoming[0]] = only,
            (NodeKey::Edge(_), _) => return Err(internal_topology_error()),
            (NodeKey::Vertex(_), _) => {
                for (into, onward) in pair_at_pinch(ctx, edges, exact, node, incident)? {
                    next[into] = onward;
                }
            }
        }
    }
    Ok((next, nodes))
}

/// An edge at a pinch vertex, seen from the vertex.
struct Spoke {
    edge: usize,
    /// Whether the edge ends at the vertex.
    incoming: bool,
    /// The exact point of its other end.
    towards: usize,
    /// Whether that end lies in the lower half-plane around the vertex:
    /// below it on the second projected axis, or level and not after it on
    /// the first.
    lower: bool,
}

/// Pair the edges at a vertex where the boundary meets itself: around the
/// vertex, counter-clockwise in the frame, outgoing and incoming edges
/// alternate, and the material lies in the wedge from an outgoing edge to the
/// next incoming one. Each incoming edge continues with the outgoing edge that
/// opens its wedge, the next outgoing edge clockwise.
fn pair_at_pinch(
    ctx: &ChainContext<'_>,
    edges: &[CutEdge],
    exact: &mut PlanePoints,
    node: NodeKey,
    incident: &Incident,
) -> GeomResult<Vec<(usize, usize)>> {
    let here = exact.node(ctx.mesh, ctx.classifications, node)?;
    let mut spokes = Vec::with_capacity(2 * incident.outgoing.len());
    for (list, incoming) in [(&incident.incoming, true), (&incident.outgoing, false)] {
        for &edge in list {
            let far = if incoming {
                edges[edge].from
            } else {
                edges[edge].to
            };
            let towards = exact.node(ctx.mesh, ctx.classifications, far)?;
            if exact.same(towards, here)? {
                return Err(GeomError::Degenerate(
                    "cut region collapsed under projection".into(),
                ));
            }
            let upper = match exact.compare(towards, here, 1)? {
                Ordering::Greater => true,
                Ordering::Equal => exact.compare(towards, here, 0)? == Ordering::Greater,
                Ordering::Less => false,
            };
            spokes.push(Spoke {
                edge,
                incoming,
                towards,
                lower: !upper,
            });
        }
    }
    // Insertion sort: the comparison is fallible and a pinch has few edges.
    for index in 1..spokes.len() {
        let mut at = index;
        while at > 0 && spoke_before(exact, here, &spokes[at], &spokes[at - 1])? {
            spokes.swap(at, at - 1);
            at -= 1;
        }
    }
    let count = spokes.len();
    for index in 0..count {
        if spokes[index].incoming == spokes[(index + 1) % count].incoming {
            return Err(refused(REFUSED_CROSSING_CUT_RINGS));
        }
    }
    Ok(spokes
        .iter()
        .enumerate()
        .filter(|(_, spoke)| spoke.incoming)
        .map(|(index, spoke)| (spoke.edge, spokes[(index + count - 1) % count].edge))
        .collect())
}

/// Whether spoke `a` comes before spoke `b` counter-clockwise from the
/// direction of the frame's positive first projected axis.
fn spoke_before(exact: &PlanePoints, here: usize, a: &Spoke, b: &Spoke) -> GeomResult<bool> {
    if a.lower != b.lower {
        return Ok(!a.lower);
    }
    match exact.orient(here, a.towards, b.towards)? {
        Sign::Positive => Ok(true),
        Sign::Negative => Ok(false),
        // Two edges leave along one ray: their pairing is not decided.
        _ => Err(refused(REFUSED_PINCHED_CUT_REGION)),
    }
}

/// The closed walks of the successor map, from the smallest unused edge.
fn walks(next: &[usize], options: &ExecutionOptions) -> GeomResult<Vec<Vec<usize>>> {
    let mut used = vec![false; next.len()];
    let mut walks = Vec::new();
    for start in 0..next.len() {
        if used[start] {
            continue;
        }
        options.check_cancelled()?;
        let mut walk = Vec::new();
        let mut edge = start;
        loop {
            if used[edge] {
                return Err(internal_topology_error());
            }
            used[edge] = true;
            walk.push(edge);
            edge = *next.get(edge).ok_or_else(internal_topology_error)?;
            if edge == start {
                break;
            }
        }
        walks.push(walk);
    }
    Ok(walks)
}

/// Split a closed walk wherever it comes back to a vertex it left before, so
/// that no ring visits a vertex twice.
fn split_walk(edges: &[CutEdge], walk: &[usize], rings: &mut Vec<Vec<usize>>) {
    let mut stack: Vec<usize> = Vec::with_capacity(walk.len());
    let mut position = BTreeMap::<NodeKey, usize>::new();
    for &edge in walk {
        let node = edges[edge].from;
        if let Some(&at) = position.get(&node) {
            let ring = stack.split_off(at);
            for &done in &ring {
                position.remove(&edges[done].from);
            }
            rings.push(ring);
        }
        position.insert(node, stack.len());
        stack.push(edge);
    }
    rings.push(stack);
}

/// A chain that passes through a pinch vertex is one chain (ADR 0091), so
/// its two segments must stay consecutive on one ring.
fn check_chains_continue(
    edges: &[CutEdge],
    nodes: &BTreeMap<NodeKey, Incident>,
    rings: &[Vec<usize>],
) -> GeomResult<()> {
    let mut after = vec![usize::MAX; edges.len()];
    for ring in rings {
        for index in 0..ring.len() {
            after[ring[index]] = ring[(index + 1) % ring.len()];
        }
    }
    for incident in nodes.values() {
        if incident.outgoing.len() < 2 {
            continue;
        }
        let cut = |list: &[usize]| -> Vec<usize> {
            list.iter()
                .copied()
                .filter(|&edge| edges[edge].link == Link::Cut)
                .collect()
        };
        match (
            cut(&incident.incoming).as_slice(),
            cut(&incident.outgoing).as_slice(),
        ) {
            ([], []) | ([_], []) | ([], [_]) => {}
            ([into], [onward]) if after[*into] == *onward => {}
            _ => return Err(refused(REFUSED_PINCHED_CUT_REGION)),
        }
    }
    Ok(())
}

/// The cut regions of a solid and the candidate pairs compared to nest them.
pub(super) fn build_cut_regions(
    inputs: &CutInputs<'_>,
    exact: &mut PlanePoints,
    extra: Vec<TracedRing>,
    limits: CutRegionLimits,
    options: &ExecutionOptions,
) -> GeomResult<(Vec<CutRegion>, usize)> {
    let ctx = inputs.chains;
    let edges = cut_edges(inputs)?;
    let (next, nodes) = successors(ctx, &edges, exact, options)?;
    let mut traced = Vec::new();
    for walk in walks(&next, options)? {
        split_walk(&edges, &walk, &mut traced);
    }
    check_chains_continue(&edges, &nodes, &traced)?;

    let mut rings = Vec::new();
    for ring_edges in &traced {
        options.check_cancelled()?;
        let count = ring_edges.len();
        let path: Vec<(NodeKey, Link)> = (0..count)
            .map(|index| {
                (
                    edges[ring_edges[index]].from,
                    edges[ring_edges[(index + count - 1) % count]].link,
                )
            })
            .collect();
        let ring = simplified_ring(inputs, &path)?;
        if ring.len() < 3 {
            return Err(GeomError::Degenerate(
                "cut region collapsed below three vertices".into(),
            ));
        }
        let points = chain_points(ctx.mesh, ctx.frame, ctx.classifications, &ring)?;
        // A fail-closed check of the returned output, which the registry
        // requires; not a geometric decision (those are on the exact points).
        if has_adjacent_duplicate(&points, true) {
            return Err(GeomError::Degenerate(
                "cut region collapsed under projection".into(),
            ));
        }
        let mut nodes = Vec::with_capacity(ring.len());
        for &node in &ring {
            nodes.push(exact.node(ctx.mesh, ctx.classifications, node)?);
        }
        for index in 0..nodes.len() {
            if exact.same(nodes[index], nodes[(index + 1) % nodes.len()])? {
                return Err(GeomError::Degenerate(
                    "cut region collapsed under projection".into(),
                ));
            }
        }
        // A ring whose vertices are all collinear has no area; the nesting
        // relies on a turn at its smallest vertex.
        let mut turns = false;
        for &node in &nodes[2..] {
            if exact.orient(nodes[0], nodes[1], node)? != Sign::Zero {
                turns = true;
                break;
            }
        }
        if !turns {
            return Err(GeomError::Degenerate(
                "cut region ring has no signed area".into(),
            ));
        }
        rings.push(TracedRing { nodes, points });
    }
    rings.extend(extra);
    nest_rings(exact, rings, limits, options)
}

/// The nodes of a traced ring that remain after dropping the redundant ones.
fn simplified_ring(inputs: &CutInputs<'_>, path: &[(NodeKey, Link)]) -> GeomResult<Vec<NodeKey>> {
    let count = path.len();
    let mut kept = Vec::with_capacity(count);
    for (index, &(node, entered)) in path.iter().enumerate() {
        let (next, leaves) = path[(index + 1) % count];
        let previous = path[(index + count - 1) % count].0;
        let droppable = match (entered, leaves) {
            (Link::Cut, Link::Cut) => is_droppable(inputs.chains, previous, node, next)?,
            (Link::Region, Link::Region) => match node {
                NodeKey::Vertex(vertex) => !inputs.region_ring_vertices.contains(&vertex),
                NodeKey::Edge(_) => return Err(internal_topology_error()),
            },
            // A junction is a vertex of both the chain and the region.
            _ => false,
        };
        if !droppable {
            kept.push(node);
        }
    }
    Ok(kept)
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec3};
    use axiolid_mesh::TriMesh;
    use axiolid_mesh_section_contract::SectionRole;

    use super::super::regions::Projection;
    use super::super::{classify_vertices, ExactSectionPlane};
    use super::*;

    /// Source vertices on the plane z = 0, all cut-segment nodes:
    /// 0 = (0,0), 1 = (1,0), 2 = (2,0), 3 = (2,1), 4 = (0,1), 5 = (-1,0).
    /// Vertex 1 lies exactly between 0 and 2, on the line y = 0.
    fn positions() -> Vec<Point3> {
        [
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (0.0, 1.0),
            (-1.0, 0.0),
        ]
        .map(|(x, y)| Point3::new(x, y, 0.0))
        .to_vec()
    }

    fn vertex(index: u32) -> NodeKey {
        NodeKey::Vertex(index)
    }

    /// Cut paths (closed or open) and region runs over `positions`, traced
    /// with the given ring vertices and region vertices.
    fn trace(
        positions: Vec<Point3>,
        cut_paths: &[(&[u32], bool)],
        region_runs: &[&[u32]],
        ring_vertices: &[u32],
        region_vertices: &[u32],
    ) -> GeomResult<(Vec<CutRegion>, usize)> {
        let mesh = TriMesh::new(positions, vec![]);
        let frame = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let plane = ExactSectionPlane::new(frame).expect("plane");
        let mut segments = BTreeSet::new();
        let mut directions = BTreeMap::new();
        for &(path, closed) in cut_paths {
            let count = if closed { path.len() } else { path.len() - 1 };
            for index in 0..count {
                let (from, to) = (vertex(path[index]), vertex(path[(index + 1) % path.len()]));
                let segment = Segment::new(from, to).expect("segment");
                directions.insert(segment, from == segment.0);
                segments.insert(segment);
            }
        }
        let region_edges: BTreeSet<(u32, u32)> = region_runs
            .iter()
            .flat_map(|run| run.windows(2).map(|pair| (pair[0], pair[1])))
            .collect();
        let region_vertices: BTreeSet<u32> = region_vertices.iter().copied().collect();
        let ring_vertices: BTreeSet<u32> = ring_vertices.iter().copied().collect();
        let classifications = classify_vertices(&mesh, plane).expect("classes");
        let chains = ChainContext {
            mesh: &mesh,
            frame,
            classifications: &classifications,
            projection: Projection::new(&plane).expect("projection"),
            directions: &directions,
            coplanar_edges: &BTreeSet::new(),
            region_vertices: &region_vertices,
            role: SectionRole::Solid,
            boundary_edges: &BTreeSet::new(),
            boundary_vertices: &BTreeSet::new(),
        };
        build_cut_regions(
            &CutInputs {
                chains: &chains,
                segments: &segments,
                region_edges: &region_edges,
                region_ring_vertices: &ring_vertices,
            },
            &mut PlanePoints::new(&plane, chains.projection).expect("points"),
            Vec::new(),
            CutRegionLimits::new(10, 100, 1_000),
            &ExecutionOptions::new(Tolerance::METRE),
        )
    }

    fn ring(points: &[(f64, f64)]) -> Vec<Point2> {
        points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    #[test]
    fn a_chain_closes_along_a_region_edge() {
        // The chain runs 2 -> 3 -> 4 -> 0; the region edge 0 -> 2 closes it.
        let (regions, _) = trace(
            positions(),
            &[(&[2, 3, 4, 0], false)],
            &[&[0, 2]],
            &[0, 2],
            &[0, 2],
        )
        .expect("cut");
        assert_eq!(regions.len(), 1);
        assert_eq!(
            regions[0].outer,
            ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
        );
        assert!(regions[0].holes.is_empty());
    }

    #[test]
    fn a_collinear_region_run_vertex_is_kept_only_if_the_region_ring_kept_it() {
        let run: &[&[u32]] = &[&[0, 1, 2]];
        let chain: &[(&[u32], bool)] = &[(&[2, 3, 4, 0], false)];
        let (dropped, _) = trace(positions(), chain, run, &[0, 2], &[0, 1, 2]).expect("cut");
        assert_eq!(
            dropped[0].outer,
            ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
        );
        let (kept, _) = trace(positions(), chain, run, &[0, 1, 2], &[0, 1, 2]).expect("cut");
        assert_eq!(
            kept[0].outer,
            ring(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
        );
    }

    #[test]
    fn a_collinear_junction_of_a_chain_and_a_region_edge_stays() {
        // Region edge 0 -> 1, then the chain 1 -> 2 -> 3 -> 4 -> 0. Vertex 1
        // lies between 0 and 2 but is where the region ends.
        let (regions, _) = trace(
            positions(),
            &[(&[1, 2, 3, 4, 0], false)],
            &[&[0, 1]],
            &[0, 1],
            &[0, 1],
        )
        .expect("cut");
        assert_eq!(
            regions[0].outer,
            ring(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])
        );
    }

    #[test]
    fn a_ring_with_no_area_is_degenerate() {
        // Vertex 1 is protected, so the three collinear vertices stay.
        let flat = trace(positions(), &[(&[0, 1, 2], true)], &[], &[], &[1]);
        assert!(matches!(flat, Err(GeomError::Degenerate(_))), "{flat:?}");
    }

    #[test]
    fn a_pinch_of_cut_segments_alone_is_refused_by_name() {
        // Two triangles share vertex 1, where four cut segments meet: a chain
        // would pass it twice (chain assembly refuses it as a saddle first).
        let pinch = trace(
            positions(),
            &[(&[0, 1, 4], true), (&[1, 2, 3], true)],
            &[],
            &[],
            &[],
        );
        assert!(matches!(
            pinch,
            Err(GeomError::UnsupportedInput {
                input: REFUSED_PINCHED_CUT_REGION,
                ..
            })
        ));
    }

    #[test]
    fn a_chain_end_without_a_region_edge_breaks_the_topology() {
        let open = trace(
            positions(),
            &[(&[2, 3, 4, 0], false)],
            &[],
            &[0, 2],
            &[0, 2],
        );
        assert!(matches!(
            open,
            Err(GeomError::BackendContractViolation { .. })
        ));
    }

    #[test]
    fn a_ring_far_from_the_origin_is_an_outer_ring_by_its_exact_area() {
        // Counter-clockwise, 1e8 from the origin, where the float area has
        // the wrong sign; a float test would call it a hole and refuse it.
        let triangle = vec![
            Point3::new(100000002.228125, 100000001.5671875, 0.0),
            Point3::new(100000001.9375, 100000001.275, 0.0),
            Point3::new(100000002.615625, 100000000.5578125, 0.0),
        ];
        let (regions, _) = trace(triangle, &[(&[0, 1, 2], true)], &[], &[], &[]).expect("cut");
        assert_eq!(regions.len(), 1);
        let reversed = trace(
            vec![
                Point3::new(100000002.228125, 100000001.5671875, 0.0),
                Point3::new(100000001.9375, 100000001.275, 0.0),
                Point3::new(100000002.615625, 100000000.5578125, 0.0),
            ],
            &[(&[2, 1, 0], true)],
            &[],
            &[],
            &[],
        );
        assert!(
            reversed.is_err(),
            "a clockwise ring alone is a hole without an outer ring"
        );
    }

    fn at(points: &[(f64, f64)]) -> Vec<Point3> {
        points
            .iter()
            .map(|&(x, y)| Point3::new(x, y, 0.0))
            .collect()
    }

    fn refused_as(result: &GeomResult<(Vec<CutRegion>, usize)>, name: &str) -> bool {
        matches!(result, Err(GeomError::UnsupportedInput { input, .. }) if *input == name)
    }

    #[test]
    fn two_rings_pinched_at_a_region_corner_split_there() {
        // 2 = (1,1) is a corner of both: A runs 0 1 2 3, B runs 2 4 5 6. The
        // region edges 1->2, 2->3, 6->2 and 2->4 meet at 2.
        let points = at(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (2.0, 1.0),
            (2.0, 2.0),
            (1.0, 2.0),
        ]);
        let corners = [1, 2, 3, 4, 6];
        let (regions, _) = trace(
            points,
            &[(&[3, 0, 1], false), (&[4, 5, 6], false)],
            &[&[1, 2, 3], &[6, 2, 4]],
            &corners,
            &corners,
        )
        .expect("cut");
        assert_eq!(regions.len(), 2);
        assert_eq!(
            regions[0].outer,
            ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
        );
        assert_eq!(
            regions[1].outer,
            ring(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)])
        );
    }

    #[test]
    fn three_rings_pinched_at_one_vertex_split_there() {
        let points = at(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (-2.0, 0.0),
            (-2.0, -2.0),
            (0.0, -2.0),
        ]);
        let corners = [0, 1, 2, 3, 4, 5, 6];
        let (regions, _) = trace(
            points,
            &[(&[1, 2], false), (&[3, 4], false), (&[5, 6], false)],
            &[&[2, 0, 1], &[4, 0, 3], &[6, 0, 5]],
            &corners,
            &corners,
        )
        .expect("cut");
        assert_eq!(
            regions
                .iter()
                .map(|region| region.outer.clone())
                .collect::<Vec<_>>(),
            [
                ring(&[(-2.0, -2.0), (0.0, -2.0), (0.0, 0.0)]),
                ring(&[(-2.0, 0.0), (0.0, 0.0), (0.0, 2.0)]),
                ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0)]),
            ]
        );
    }

    #[test]
    fn a_hole_pinched_against_its_outer_ring_is_a_hole_that_touches_it() {
        // The closed chain passes 1 = (1.5,0); the region edges 1->5, 5->6,
        // 6->1 bound a triangular hole whose corner is 1.
        let points = at(&[
            (0.0, 0.0),
            (1.5, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (0.0, 3.0),
            (0.5, 2.0),
            (2.5, 2.0),
        ]);
        let (regions, _) = trace(
            points,
            &[(&[0, 1, 2, 3, 4], true)],
            &[&[1, 5, 6, 1]],
            &[1, 5, 6],
            &[1, 5, 6],
        )
        .expect("cut");
        assert_eq!(regions.len(), 1);
        assert_eq!(
            regions[0].outer,
            ring(&[(0.0, 0.0), (1.5, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)])
        );
        assert_eq!(
            regions[0].holes,
            [ring(&[(0.5, 2.0), (2.5, 2.0), (1.5, 0.0)])]
        );
    }

    #[test]
    fn a_chain_parted_at_a_pinch_is_refused_by_name() {
        // The chain 2 -> 0 -> 3 passes 0, where P (0 1 2) and Q (0 3 4) touch:
        // its two segments end up on different rings.
        let points = at(&[
            (0.0, 0.0),
            (-1.0, -1.0),
            (0.0, -1.0),
            (1.0, 0.0),
            (1.0, 1.0),
        ]);
        let all = [0, 1, 2, 3, 4];
        let parted = trace(
            points,
            &[(&[2, 0, 3], false)],
            &[&[0, 1, 2], &[3, 4, 0]],
            &all,
            &all,
        );
        assert!(
            refused_as(&parted, REFUSED_PINCHED_CUT_REGION),
            "{parted:?}"
        );
    }

    #[test]
    fn edges_leaving_a_pinch_along_one_ray_are_refused_by_name() {
        // A leaves 0 towards 1 = (2,0); B enters 0 from 3 = (1,0).
        let points = at(&[(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 0.0), (0.0, -1.0)]);
        let all = [0, 1, 2, 3, 4];
        let ray = trace(points, &[], &[&[0, 1, 2, 0], &[3, 0, 4, 3]], &all, &all);
        assert!(refused_as(&ray, REFUSED_PINCHED_CUT_REGION), "{ray:?}");
    }

    #[test]
    fn two_outgoing_edges_side_by_side_at_a_pinch_are_crossing_rings() {
        // A leaves 0 at 0 degrees and returns at 90; B leaves at 45 and
        // returns at 135: their material overlaps between 45 and 90.
        let points = at(&[(0.0, 0.0), (2.0, 0.0), (0.0, 2.0), (1.0, 1.0), (-1.0, 1.0)]);
        let all = [0, 1, 2, 3, 4];
        let overlap = trace(points, &[], &[&[0, 1, 2, 0], &[0, 3, 4, 0]], &all, &all);
        assert!(
            refused_as(&overlap, REFUSED_CROSSING_CUT_RINGS),
            "{overlap:?}"
        );
    }

    /// A square whose bottom and top pass 1 = (1.5,0) and 4 = (1.5,3), and a
    /// diamond 1 6 4 7 that touches it at both.
    fn split_square() -> Vec<Point3> {
        at(&[
            (0.0, 0.0),
            (1.5, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (1.5, 3.0),
            (0.0, 3.0),
            (0.5, 1.5),
            (2.5, 1.5),
        ])
    }

    #[test]
    fn a_hole_touching_its_outer_ring_twice_leaves_two_regions() {
        // All region edges: the diamond cuts the square's material in two, and
        // each piece is one region with a connected inside.
        let all = [0, 1, 2, 3, 4, 5, 6, 7];
        let (regions, _) = trace(
            split_square(),
            &[],
            &[&[0, 1, 2, 3, 4, 5, 0], &[1, 6, 4, 7, 1]],
            &all,
            &all,
        )
        .expect("cut");
        assert_eq!(
            regions
                .iter()
                .map(|region| (region.outer.clone(), region.holes.len()))
                .collect::<Vec<_>>(),
            [
                (
                    ring(&[(0.0, 0.0), (1.5, 0.0), (0.5, 1.5), (1.5, 3.0), (0.0, 3.0)]),
                    0
                ),
                (
                    ring(&[(1.5, 0.0), (3.0, 0.0), (3.0, 3.0), (1.5, 3.0), (2.5, 1.5)]),
                    0
                ),
            ]
        );
    }

    #[test]
    fn a_closed_chain_through_both_contacts_of_such_a_hole_is_refused_by_name() {
        // The same cut with the square as one closed chain: the chain would
        // run along two rings.
        let all = [0, 1, 2, 3, 4, 5, 6, 7];
        let parted = trace(
            split_square(),
            &[(&[0, 1, 2, 3, 4, 5], true)],
            &[&[1, 6, 4, 7, 1]],
            &all,
            &all,
        );
        assert!(
            refused_as(&parted, REFUSED_PINCHED_CUT_REGION),
            "{parted:?}"
        );
    }

    #[test]
    fn a_zero_length_edge_at_a_pinch_is_degenerate() {
        // B leaves 0 for 3, a second vertex at the same position.
        let points = at(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (0.0, 2.0),
            (0.0, 0.0),
            (-2.0, 0.0),
            (0.0, -2.0),
        ]);
        let all = [0, 1, 2, 3, 4, 5];
        let zero = trace(points, &[], &[&[0, 1, 2, 0], &[0, 3, 4, 5, 0]], &all, &all);
        assert!(matches!(zero, Err(GeomError::Degenerate(_))), "{zero:?}");
    }
}
