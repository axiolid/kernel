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
//! 2. A node with unequal degrees breaks the topology. A vertex with two
//!    outgoing edges is a pinch, where the boundary touches itself, and is
//!    refused by name.
//! 3. Each ring is traced from its smallest node. A node is dropped when both
//!    its edges are cut segments and they are exactly collinear
//!    (`chains::is_droppable`), or when both are region edges and the region
//!    ring itself dropped the vertex. A junction of a cut segment and a region
//!    edge is a vertex of both outputs and stays.
//! 4. Every node gets its exact point (`plane_points`): a source vertex in the
//!    plane or the exact crossing of a source edge with it. `ring_nesting`
//!    decides on those points whether a ring is an outer ring or a hole and
//!    groups them; the returned coordinates only label the output.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation, Sign};
use axiolid_mesh_section_contract::{CutRegion, CutRegionLimits, REFUSED_PINCHED_CUT_REGION};

use super::chains::{chain_points, has_adjacent_duplicate, is_droppable, ChainContext};
use super::plane_points::PlanePoints;
use super::ring_nesting::{nest_rings, TracedRing};
use super::{internal_topology_error, ExactSectionPlane, NodeKey, ScalarSection, Segment};

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
    pub(super) plane: ExactSectionPlane,
    pub(super) segments: &'a BTreeSet<Segment>,
    /// Region edges that bound the cut, directed with the material on the left.
    pub(super) region_edges: &'a BTreeSet<(u32, u32)>,
    /// Source vertices the in-plane regions keep on their rings.
    pub(super) region_ring_vertices: &'a BTreeSet<u32>,
}

type Adjacency = BTreeMap<NodeKey, Vec<(NodeKey, Link)>>;

fn add_edge(
    outgoing: &mut Adjacency,
    incoming: &mut Adjacency,
    from: NodeKey,
    to: NodeKey,
    link: Link,
) {
    outgoing.entry(from).or_default().push((to, link));
    incoming.entry(to).or_default().push((from, link));
}

/// The cut regions of a solid and the candidate pairs compared to nest them.
pub(super) fn build_cut_regions(
    inputs: &CutInputs<'_>,
    limits: CutRegionLimits,
    options: &ExecutionOptions,
) -> GeomResult<(Vec<CutRegion>, usize)> {
    let ctx = inputs.chains;
    let mut outgoing = Adjacency::new();
    let mut incoming = Adjacency::new();
    for &segment in inputs.segments {
        let Segment(left, right) = segment;
        let forward = *ctx
            .directions
            .get(&segment)
            .ok_or_else(internal_topology_error)?;
        let (from, to) = if forward {
            (left, right)
        } else {
            (right, left)
        };
        add_edge(&mut outgoing, &mut incoming, from, to, Link::Cut);
    }
    for &(from, to) in inputs.region_edges {
        add_edge(
            &mut outgoing,
            &mut incoming,
            NodeKey::Vertex(from),
            NodeKey::Vertex(to),
            Link::Region,
        );
    }

    for node in outgoing.keys().chain(incoming.keys()) {
        let out = outgoing.get(node).map_or(0, Vec::len);
        let inn = incoming.get(node).map_or(0, Vec::len);
        if out != inn {
            return Err(internal_topology_error());
        }
        if out > 1 {
            return match node {
                NodeKey::Vertex(_) => Err(GeomError::UnsupportedInput {
                    backend: ScalarSection::ID,
                    operation: Operation::MeshPlaneSection,
                    input: REFUSED_PINCHED_CUT_REGION,
                }),
                NodeKey::Edge(_) => Err(internal_topology_error()),
            };
        }
    }

    let mut exact = PlanePoints::new(&inputs.plane, ctx.projection)?;
    let mut visited = BTreeSet::<NodeKey>::new();
    let mut rings = Vec::new();
    for (&start, onward) in &outgoing {
        if visited.contains(&start) {
            continue;
        }
        options.check_cancelled()?;
        // Each node with the link that enters it.
        let mut path = Vec::new();
        let mut node = start;
        let mut entered = incoming[&start][0].1;
        let mut link = onward[0].1;
        let mut next = onward[0].0;
        loop {
            visited.insert(node);
            path.push((node, entered));
            entered = link;
            node = next;
            if node == start {
                break;
            }
            if path.len() > outgoing.len() {
                return Err(internal_topology_error());
            }
            (next, link) = outgoing[&node][0];
        }

        let ring = simplified_ring(inputs, &path)?;
        if ring.len() < 3 {
            return Err(GeomError::Degenerate(
                "cut region collapsed below three vertices".into(),
            ));
        }
        let points = chain_points(ctx.mesh, ctx.frame, ctx.classifications, &ring)?;
        if has_adjacent_duplicate(&points, true) {
            return Err(GeomError::Degenerate(
                "cut region collapsed under projection".into(),
            ));
        }
        let mut nodes = Vec::with_capacity(ring.len());
        for &node in &ring {
            nodes.push(exact.node(ctx.mesh, node)?);
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
    nest_rings(&exact, rings, limits, options)
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
    use super::super::ExactSectionPlane;
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
        let chains = ChainContext {
            mesh: &mesh,
            frame,
            classifications: &[],
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
                plane,
                segments: &segments,
                region_edges: &region_edges,
                region_ring_vertices: &ring_vertices,
            },
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
    fn a_pinch_is_refused_by_name() {
        // Two triangles share vertex 1: it has two outgoing edges.
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
}
