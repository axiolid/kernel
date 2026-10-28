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

use super::regions::canonical_rotation;
use super::{
    internal_topology_error, node_point, point_order, project_local, sequence_order, signed_area,
    simplify_collinear, Classification, NodeKey, ScalarSection, Segment,
};
use crate::orientation::orient2d;

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

/// Open and closed chains of the cut segments of a request (ADR 0091).
///
/// A vertex node of degree one ends an open chain; it must be a vertex of an
/// in-plane triangle, where the chain meets a region boundary. A vertex node
/// of degree above two is a saddle and is refused by name. Points of vertices
/// listed in `region_vertices` are never simplified away, so a chain end
/// stays bit-identical to the region vertex it meets.
pub(super) fn assemble_chains(
    mesh: &TriMesh,
    frame: Frame3,
    classifications: &[Classification],
    segments: &BTreeSet<Segment>,
    region_vertices: &BTreeSet<u32>,
    limits: SectionLimits,
    options: &ExecutionOptions,
) -> GeomResult<Vec<SectionChain>> {
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
        let mut points = chain_points(mesh, frame, classifications, &nodes)?;
        let mut protected: Vec<bool> = nodes.iter().map(|&node| is_region_vertex(node)).collect();
        simplify_chain(
            &mut points,
            &mut protected,
            options.tolerance().linear(),
            false,
        );
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
        let mut points = chain_points(mesh, frame, classifications, &nodes)?;
        let mut protected: Vec<bool> = nodes.iter().map(|&node| is_region_vertex(node)).collect();
        simplify_chain(
            &mut points,
            &mut protected,
            options.tolerance().linear(),
            true,
        );
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

/// The legacy collinear simplification, with `protected` points kept and, for
/// an open chain, its two ends kept as well.
fn simplify_chain(points: &mut Vec<Point2>, protected: &mut Vec<bool>, tolerance: f64, ring: bool) {
    loop {
        let minimum = if ring { 3 } else { 2 };
        if points.len() <= minimum {
            return;
        }
        let range = if ring {
            0..points.len()
        } else {
            1..points.len() - 1
        };
        let mut removed = false;
        for index in range {
            if protected[index] {
                continue;
            }
            let previous = points[(index + points.len() - 1) % points.len()];
            let current = points[index];
            let next = points[(index + 1) % points.len()];
            let chord = next - previous;
            let scale = chord.length();
            let cross = (current - previous).perp_dot(chord).abs();
            let within_tolerance = scale > 0.0 && cross <= tolerance * scale;
            if orient2d(previous, current, next).sign() == Some(Sign::Zero) || within_tolerance {
                points.remove(index);
                protected.remove(index);
                removed = true;
                break;
            }
        }
        if !removed {
            return;
        }
    }
}
