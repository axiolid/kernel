//! Assembly of cut segments into plane-local polylines.
//!
//! `assemble_contours` is the legacy path: every node has degree two and every
//! polyline is a closed contour. `assemble_chains` serves requests (ADR 0091):
//! it also yields open chains that end on an in-plane region, and refuses a
//! vertex where more than two cut segments meet.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation, Sign};
use axiolid_core::{Frame3, Point2};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    ChainEnd, SectionChain, SectionContour, SectionLimits, REFUSED_SADDLE_VERTEX,
};

use super::regions::{canonical_rotation, corners, strictly_between, Projection};
use super::{
    internal_topology_error, node_point, point_order, project_local, sequence_order, signed_area,
    simplify_collinear, Classification, EdgeKey, NodeKey, ScalarSection, Segment,
};

pub(super) fn assemble_contours(
    mesh: &TriMesh,
    frame: Frame3,
    classifications: &[Classification],
    segments: &BTreeSet<Segment>,
    limits: SectionLimits,
    options: &ExecutionOptions,
) -> GeomResult<Vec<SectionContour>> {
    let mut adjacency = BTreeMap::<NodeKey, Vec<NodeKey>>::new();
    for &Segment(left, right) in segments {
        adjacency.entry(left).or_default().push(right);
        adjacency.entry(right).or_default().push(left);
    }
    for neighbours in adjacency.values_mut() {
        neighbours.sort_unstable();
        if neighbours.len() != 2 {
            return Err(GeomError::NotManifold(format!(
                "section graph has degree {}, expected 2",
                neighbours.len()
            )));
        }
    }

    let mut visited = BTreeSet::<Segment>::new();
    let mut contours = Vec::new();
    for &start in adjacency.keys() {
        let start_is_complete = adjacency[&start].iter().try_fold(true, |complete, &next| {
            let segment = Segment::new(start, next)?;
            Ok::<bool, GeomError>(complete && visited.contains(&segment))
        })?;
        if start_is_complete {
            continue;
        }
        options.check_cancelled()?;
        if contours.len() >= limits.max_contours {
            return Err(GeomError::BudgetExceeded {
                resource: "section contours",
            });
        }
        let mut nodes = Vec::new();
        let mut previous = None;
        let mut current = start;
        loop {
            if nodes.len() >= limits.max_output_vertices {
                return Err(GeomError::BudgetExceeded {
                    resource: "section output vertices",
                });
            }
            nodes.push(current);
            let neighbours = adjacency
                .get(&current)
                .ok_or_else(internal_topology_error)?;
            let next = match previous {
                None => neighbours[0],
                Some(previous) if neighbours[0] == previous => neighbours[1],
                Some(_) => neighbours[0],
            };
            let edge = Segment::new(current, next)?;
            if !visited.insert(edge) && next != start {
                return Err(GeomError::NotManifold(
                    "section graph revisited an edge before closing a contour".into(),
                ));
            }
            previous = Some(current);
            current = next;
            if current == start {
                break;
            }
            if nodes.len() > segments.len() {
                return Err(internal_topology_error());
            }
        }
        if nodes.len() < 3 {
            return Err(GeomError::Degenerate(
                "section contour has fewer than three source-topology nodes".into(),
            ));
        }
        let mut points = Vec::new();
        points
            .try_reserve_exact(nodes.len())
            .map_err(|_| GeomError::BudgetExceeded { resource: "memory" })?;
        for node in nodes {
            let world = node_point(mesh, classifications, node)?;
            let local = world - frame.origin;
            let point = Point2::new(local.dot(frame.x), local.dot(frame.y));
            if !point.is_finite() {
                return Err(GeomError::Degenerate(
                    "section projection exceeded finite arithmetic".into(),
                ));
            }
            points.push(point);
        }
        simplify_collinear(&mut points, options.tolerance().linear());
        if points.len() < 3 {
            return Err(GeomError::Degenerate(
                "section contour collapsed below three vertices".into(),
            ));
        }
        let area = signed_area(&points);
        if !area.is_finite() || area == 0.0 {
            return Err(GeomError::Degenerate(
                "section contour has no finite signed area".into(),
            ));
        }
        if area < 0.0 {
            points[1..].reverse();
        }
        contours.push(SectionContour::new(points));
    }
    contours.sort_by(|left, right| point_order(&left.points[0], &right.points[0]));
    Ok(contours)
}

/// What chain assembly needs besides the segments.
pub(super) struct ChainContext<'a> {
    pub(super) mesh: &'a TriMesh,
    pub(super) frame: Frame3,
    pub(super) classifications: &'a [Classification],
    pub(super) projection: Projection,
    /// Per segment: whether the cut runs from its first node to its second,
    /// with the solid on the left (counter-clockwise in the frame). Decided
    /// from the exact vertex signs of the segment's triangle.
    pub(super) directions: &'a BTreeMap<Segment, bool>,
    /// Crossing edges whose two triangles are exactly coplanar: the cut points
    /// on either side of such an edge's node are collinear with it.
    pub(super) coplanar_edges: &'a BTreeSet<EdgeKey>,
    /// Vertices of in-plane triangles.
    pub(super) region_vertices: &'a BTreeSet<u32>,
}

/// Open and closed chains of the cut segments of a request (ADR 0091).
///
/// A vertex node of degree one ends an open chain; it must be a vertex of an
/// in-plane triangle, where the chain meets a region boundary. A vertex node
/// of degree above two is a saddle and is refused by name.
///
/// Simplification and orientation are exact, decided on source data. A closed
/// chain is made counter-clockwise by the segment directions. A node is
/// dropped only when its two cut segments are exactly collinear: a vertex node
/// between two vertex nodes with a zero `orient2d`, or an edge node whose two
/// triangles are exactly coplanar (`orient3d`). Vertices of in-plane triangles
/// are never dropped, so a chain end stays bit-identical to the region vertex
/// it meets, and the ends of an open chain are never dropped.
pub(super) fn assemble_chains(
    ctx: &ChainContext<'_>,
    segments: &BTreeSet<Segment>,
    limits: SectionLimits,
    options: &ExecutionOptions,
) -> GeomResult<Vec<SectionChain>> {
    let (mesh, frame, classifications) = (ctx.mesh, ctx.frame, ctx.classifications);
    let region_vertices = ctx.region_vertices;
    let mut adjacency = BTreeMap::<NodeKey, Vec<NodeKey>>::new();
    for &Segment(left, right) in segments {
        adjacency.entry(left).or_default().push(right);
        adjacency.entry(right).or_default().push(left);
    }
    for (node, neighbours) in &mut adjacency {
        neighbours.sort_unstable();
        match (node, neighbours.len()) {
            (_, 1) | (_, 2) => {}
            (NodeKey::Vertex(_), _) => {
                return Err(GeomError::UnsupportedInput {
                    backend: ScalarSection::ID,
                    operation: Operation::MeshPlaneSection,
                    input: REFUSED_SADDLE_VERTEX,
                })
            }
            (NodeKey::Edge(_), _) => return Err(internal_topology_error()),
        }
    }

    let is_region_vertex = |node: NodeKey| matches!(node, NodeKey::Vertex(vertex) if region_vertices.contains(&vertex));
    let mut visited = BTreeSet::<Segment>::new();
    let mut chains = Vec::new();
    let mut chain_vertices = 0usize;

    // Open chains first, from their degree-one ends in sorted order.
    for (&start, neighbours) in &adjacency {
        if neighbours.len() != 1 || visited.contains(&Segment::new(start, neighbours[0])?) {
            continue;
        }
        options.check_cancelled()?;
        let mut nodes = vec![start];
        let mut previous = start;
        let mut current = neighbours[0];
        loop {
            visited.insert(Segment::new(previous, current)?);
            nodes.push(current);
            let onward = &adjacency[&current];
            if onward.len() == 1 {
                break;
            }
            let next = if onward[0] == previous {
                onward[1]
            } else {
                onward[0]
            };
            previous = current;
            current = next;
            if nodes.len() > segments.len() + 1 {
                return Err(internal_topology_error());
            }
        }
        // Landing 1: an open chain can only end on an in-plane region.
        if !(is_region_vertex(start) && is_region_vertex(current)) {
            return Err(internal_topology_error());
        }
        let nodes = simplify_nodes(ctx, nodes, false)?;
        let mut points = chain_points(mesh, frame, classifications, &nodes)?;
        if points.len() < 2 || has_adjacent_duplicate(&points, false) {
            return Err(GeomError::Degenerate(
                "section chain collapsed under projection".into(),
            ));
        }
        match point_order(&points[0], &points[points.len() - 1]) {
            std::cmp::Ordering::Equal => {
                return Err(GeomError::Degenerate(
                    "open section chain starts and ends at one point".into(),
                ))
            }
            std::cmp::Ordering::Greater => points.reverse(),
            std::cmp::Ordering::Less => {}
        }
        push_chain(
            &mut chains,
            &mut chain_vertices,
            SectionChain::open(points, ChainEnd::InPlaneRegion, ChainEnd::InPlaneRegion),
            limits,
        )?;
    }

    // What remains is closed: every node has degree two.
    for &start in adjacency.keys() {
        let complete = adjacency[&start].iter().try_fold(true, |complete, &next| {
            Ok::<bool, GeomError>(complete && visited.contains(&Segment::new(start, next)?))
        })?;
        if complete {
            continue;
        }
        options.check_cancelled()?;
        let mut nodes = Vec::new();
        let mut previous = None;
        let mut current = start;
        loop {
            nodes.push(current);
            let neighbours = &adjacency[&current];
            let next = match previous {
                None => neighbours[0],
                Some(previous) if neighbours[0] == previous => neighbours[1],
                Some(_) => neighbours[0],
            };
            visited.insert(Segment::new(current, next)?);
            previous = Some(current);
            current = next;
            if current == start {
                break;
            }
            if nodes.len() > segments.len() {
                return Err(internal_topology_error());
            }
        }
        if nodes.len() < 3 {
            return Err(GeomError::Degenerate(
                "section contour has fewer than three source-topology nodes".into(),
            ));
        }
        orient_nodes(ctx, &mut nodes)?;
        let nodes = simplify_nodes(ctx, nodes, true)?;
        if nodes.len() < 3 {
            return Err(GeomError::Degenerate(
                "section contour collapsed below three vertices".into(),
            ));
        }
        let points = chain_points(mesh, frame, classifications, &nodes)?;
        if has_adjacent_duplicate(&points, true) {
            return Err(GeomError::Degenerate(
                "section chain collapsed under projection".into(),
            ));
        }
        push_chain(
            &mut chains,
            &mut chain_vertices,
            SectionChain::closed(canonical_rotation(points)),
            limits,
        )?;
    }

    chains.sort_by(|left, right| {
        sequence_order(&left.points, &right.points)
            .then_with(|| right.is_closed().cmp(&left.is_closed()))
    });
    Ok(chains)
}

fn push_chain(
    chains: &mut Vec<SectionChain>,
    chain_vertices: &mut usize,
    chain: SectionChain,
    limits: SectionLimits,
) -> GeomResult<()> {
    if chains.len() >= limits.max_contours {
        return Err(GeomError::BudgetExceeded {
            resource: "section contours",
        });
    }
    *chain_vertices += chain.points.len();
    if *chain_vertices > limits.max_output_vertices {
        return Err(GeomError::BudgetExceeded {
            resource: "section output vertices",
        });
    }
    chains.push(chain);
    Ok(())
}

fn chain_points(
    mesh: &TriMesh,
    frame: Frame3,
    classifications: &[Classification],
    nodes: &[NodeKey],
) -> GeomResult<Vec<Point2>> {
    nodes
        .iter()
        .map(|&node| project_local(frame, node_point(mesh, classifications, node)?))
        .collect()
}

fn has_adjacent_duplicate(points: &[Point2], ring: bool) -> bool {
    let pairs = if ring {
        points.len()
    } else {
        points.len().saturating_sub(1)
    };
    (0..pairs).any(|index| points[index] == points[(index + 1) % points.len()])
}

/// Make a closed chain run counter-clockwise by its segment directions; the
/// directions of all its segments must agree.
fn orient_nodes(ctx: &ChainContext<'_>, nodes: &mut [NodeKey]) -> GeomResult<()> {
    let mut forward = None;
    for index in 0..nodes.len() {
        let (from, to) = (nodes[index], nodes[(index + 1) % nodes.len()]);
        let segment = Segment::new(from, to)?;
        let runs_forward = ctx
            .directions
            .get(&segment)
            .map(|&first_to_second| first_to_second == (from == segment.0))
            .ok_or_else(internal_topology_error)?;
        match forward {
            None => forward = Some(runs_forward),
            Some(known) if known != runs_forward => {
                return Err(GeomError::NotManifold(
                    "section chain segments disagree on their direction".into(),
                ))
            }
            Some(_) => {}
        }
    }
    if forward == Some(false) {
        nodes.reverse();
    }
    Ok(())
}

/// Whether the cut segments on both sides of a node are exactly collinear,
/// so the node can go without moving the polyline.
fn is_droppable(
    ctx: &ChainContext<'_>,
    previous: NodeKey,
    node: NodeKey,
    next: NodeKey,
) -> GeomResult<bool> {
    match node {
        NodeKey::Edge(edge) => Ok(ctx.coplanar_edges.contains(&edge)),
        NodeKey::Vertex(vertex) => {
            if ctx.region_vertices.contains(&vertex) {
                return Ok(false);
            }
            let (NodeKey::Vertex(before), NodeKey::Vertex(after)) = (previous, next) else {
                return Ok(false);
            };
            let [a, b, c] = corners(ctx.mesh, [before, vertex, after])?;
            Ok(ctx.projection.orient(a, b, c)? == Sign::Zero && strictly_between(a, b, c))
        }
    }
}

/// Drop every node of a chain whose two cut segments are exactly collinear.
/// Collinearity is judged on the original neighbours: a dropped neighbour lies
/// on the same line, so one pass suffices. An open chain keeps its ends.
fn simplify_nodes(
    ctx: &ChainContext<'_>,
    nodes: Vec<NodeKey>,
    ring: bool,
) -> GeomResult<Vec<NodeKey>> {
    let last = nodes.len() - 1;
    let mut kept = Vec::with_capacity(nodes.len());
    for (index, &node) in nodes.iter().enumerate() {
        if !ring && (index == 0 || index == last) {
            kept.push(node);
            continue;
        }
        let previous = nodes[(index + nodes.len() - 1) % nodes.len()];
        let next = nodes[(index + 1) % nodes.len()];
        if !is_droppable(ctx, previous, node, next)? {
            kept.push(node);
        }
    }
    Ok(kept)
}

#[cfg(test)]
mod tests {
    use axiolid_core::{Point3, Tolerance, Vec3};

    use super::super::ExactSectionPlane;
    use super::*;

    /// Source vertices on the plane z = 0, all cut-segment nodes:
    /// 0 = (0,0), 1 = (1,0), 2 = (2,0), 3 = (2,1), 4 = (0,1), 5 = (-1,0).
    /// Vertex 1 lies exactly between 0 and 2, on the line y = 0.
    struct Fixture {
        mesh: TriMesh,
        frame: Frame3,
        projection: Projection,
    }

    fn fixture() -> Fixture {
        let positions = [
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (0.0, 1.0),
            (-1.0, 0.0),
        ]
        .map(|(x, y)| Point3::new(x, y, 0.0))
        .to_vec();
        let frame = Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        };
        let plane = ExactSectionPlane::new(frame).expect("plane");
        Fixture {
            mesh: TriMesh::new(positions, vec![]),
            frame,
            projection: Projection::new(&plane).expect("projection"),
        }
    }

    fn vertex(index: u32) -> NodeKey {
        NodeKey::Vertex(index)
    }

    /// Chains of the given vertex path, run in the given direction.
    fn chains(path: &[u32], closed: bool, region_vertices: &[u32]) -> Vec<SectionChain> {
        let fixture = fixture();
        let region_vertices: BTreeSet<u32> = region_vertices.iter().copied().collect();
        let mut segments = BTreeSet::new();
        let mut directions = BTreeMap::new();
        let count = if closed { path.len() } else { path.len() - 1 };
        for index in 0..count {
            let (from, to) = (vertex(path[index]), vertex(path[(index + 1) % path.len()]));
            let segment = Segment::new(from, to).expect("segment");
            directions.insert(segment, from == segment.0);
            segments.insert(segment);
        }
        let coplanar_edges = BTreeSet::new();
        let ctx = ChainContext {
            mesh: &fixture.mesh,
            frame: fixture.frame,
            classifications: &[],
            projection: fixture.projection,
            directions: &directions,
            coplanar_edges: &coplanar_edges,
            region_vertices: &region_vertices,
        };
        let options = ExecutionOptions::new(Tolerance::METRE);
        assemble_chains(
            &ctx,
            &segments,
            SectionLimits::new(100, 100, 100, 10),
            &options,
        )
        .expect("chains")
    }

    fn points(chain: &SectionChain) -> Vec<(f64, f64)> {
        chain.points.iter().map(|p| (p.x, p.y)).collect()
    }

    #[test]
    fn a_closed_chain_drops_an_exactly_collinear_vertex_unless_it_is_protected() {
        // 0 -> 1 -> 2 -> 3 -> 4 is counter-clockwise; vertex 1 is collinear.
        let path = [0, 1, 2, 3, 4];
        let plain = chains(&path, true, &[]);
        assert_eq!(
            points(&plain[0]),
            [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)]
        );
        let kept = chains(&path, true, &[1]);
        assert_eq!(
            points(&kept[0]),
            [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)]
        );
    }

    #[test]
    fn an_open_chain_keeps_its_ends_and_a_protected_vertex_it_passes_through() {
        // The ends are region vertices; the middle one is collinear.
        let plain = chains(&[0, 1, 2], false, &[0, 2]);
        assert_eq!(points(&plain[0]), [(0.0, 0.0), (2.0, 0.0)]);
        let through = chains(&[0, 1, 2], false, &[0, 1, 2]);
        assert_eq!(points(&through[0]), [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
    }

    #[test]
    fn the_ends_of_an_open_chain_are_never_candidates() {
        // 0 lies between 5 and 1, so a ring would drop it; an open chain's
        // end must stay.
        let fixture = fixture();
        let region_vertices = BTreeSet::new();
        let coplanar_edges = BTreeSet::new();
        let directions = BTreeMap::new();
        let ctx = ChainContext {
            mesh: &fixture.mesh,
            frame: fixture.frame,
            classifications: &[],
            projection: fixture.projection,
            directions: &directions,
            coplanar_edges: &coplanar_edges,
            region_vertices: &region_vertices,
        };
        let nodes = vec![vertex(0), vertex(1), vertex(5)];
        assert_eq!(
            simplify_nodes(&ctx, nodes.clone(), false).expect("nodes"),
            nodes
        );
    }
}
