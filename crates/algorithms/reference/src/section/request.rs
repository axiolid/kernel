//! `ScalarSection::section_request` (ADR 0091): faces lying in the plane
//! become regions, and cut curves are transversal crossings only.
//!
//! The vertex classification and segment building are the legacy path's.
//! What differs is the treatment of an edge lying in the plane, by its
//! incident triangles:
//!
//! - two off-plane triangles whose third vertices lie on opposite sides: a
//!   cut segment between the edge's ends;
//! - two off-plane triangles on one side: a touching edge, counted and not
//!   emitted;
//! - an in-plane triangle: a region boundary, never a cut.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{PlaneSection, PlaneSectionEvidence, SectionRequest};

use super::chains::assemble_chains;
use super::regions::{build_regions, check_overlaps, InPlaneSet, InPlaneTriangle, Projection};
use super::{
    check_source_limits, classify_vertices, insert_segment, internal_topology_error, mesh_index,
    opposite, EdgeKey, ExactSectionPlane, NodeKey, Segment,
};

/// A triangle incident to a mesh edge that lies in the plane.
#[derive(Debug, Clone, Copy)]
enum Incidence {
    /// An off-plane triangle with its third vertex on this side.
    Off(Sign),
    /// An in-plane triangle.
    In,
}

pub(super) fn plane_section(
    mesh: &TriMesh,
    request: &SectionRequest,
    options: &ExecutionOptions,
) -> GeomResult<PlaneSection> {
    options.check_cancelled()?;
    let limits = request.limits;
    check_source_limits(mesh, limits)?;

    let frame = request.frame;
    let plane = ExactSectionPlane::new(frame)?;
    let projection = Projection::new(&plane)?;
    let classifications = classify_vertices(mesh, plane)?;

    let mut segments = BTreeSet::<Segment>::new();
    let mut on_plane_edges = BTreeMap::<EdgeKey, Vec<Incidence>>::new();
    let mut in_plane = Vec::new();
    for (triangle_index, triangle) in mesh.triangles().enumerate() {
        options.check_cancelled()?;
        let vertices = [
            mesh_index(triangle[0])?,
            mesh_index(triangle[1])?,
            mesh_index(triangle[2])?,
        ];
        let signs = vertices.map(|index| classifications[index].sign);
        let zero_count = signs.iter().filter(|&&sign| sign == Sign::Zero).count();
        match zero_count {
            3 => {
                in_plane.push(InPlaneTriangle::new(
                    mesh,
                    projection,
                    triangle_index,
                    triangle,
                )?);
                for (from, to) in [(0, 1), (1, 2), (2, 0)] {
                    on_plane_edges
                        .entry(EdgeKey::new(triangle[from], triangle[to]))
                        .or_default()
                        .push(Incidence::In);
                }
            }
            2 => {
                let mut zeros = triangle
                    .into_iter()
                    .zip(signs)
                    .filter_map(|(index, sign)| (sign == Sign::Zero).then_some(index));
                let left = zeros.next().ok_or_else(internal_topology_error)?;
                let right = zeros.next().ok_or_else(internal_topology_error)?;
                let third = signs
                    .into_iter()
                    .find(|&sign| sign != Sign::Zero)
                    .ok_or_else(internal_topology_error)?;
                on_plane_edges
                    .entry(EdgeKey::new(left, right))
                    .or_default()
                    .push(Incidence::Off(third));
            }
            1 => {
                let zero_corner = (0..3)
                    .find(|&corner| signs[corner] == Sign::Zero)
                    .ok_or_else(internal_topology_error)?;
                let (first, second) = match zero_corner {
                    0 => (1, 2),
                    1 => (0, 2),
                    2 => (0, 1),
                    _ => return Err(internal_topology_error()),
                };
                if opposite(signs[first], signs[second]) {
                    let segment = Segment::new(
                        NodeKey::Vertex(triangle[zero_corner]),
                        NodeKey::Edge(EdgeKey::new(triangle[first], triangle[second])),
                    )?;
                    insert_segment(&mut segments, segment, limits)?;
                }
            }
            0 => {
                let mut crossing = [None, None];
                let mut crossing_count = 0usize;
                for (left, right) in [(0, 1), (1, 2), (2, 0)] {
                    if opposite(signs[left], signs[right]) {
                        if crossing_count >= crossing.len() {
                            return Err(internal_topology_error());
                        }
                        crossing[crossing_count] =
                            Some(NodeKey::Edge(EdgeKey::new(triangle[left], triangle[right])));
                        crossing_count += 1;
                    }
                }
                match (crossing[0], crossing[1]) {
                    (Some(first), Some(second)) => {
                        insert_segment(&mut segments, Segment::new(first, second)?, limits)?;
                    }
                    (None, None) => {}
                    _ => return Err(internal_topology_error()),
                }
            }
            _ => return Err(internal_topology_error()),
        }
    }

    let mut touching_edges = 0usize;
    for (edge, incidences) in on_plane_edges {
        options.check_cancelled()?;
        match incidences.as_slice() {
            [Incidence::Off(first), Incidence::Off(second)] => {
                if opposite(*first, *second) {
                    insert_segment(
                        &mut segments,
                        Segment::new(NodeKey::Vertex(edge.0), NodeKey::Vertex(edge.1))?,
                        limits,
                    )?;
                } else {
                    touching_edges += 1;
                }
            }
            [Incidence::In, Incidence::Off(_)] | [Incidence::Off(_), Incidence::In] => {}
            [Incidence::In, Incidence::In] => {}
            _ => {
                return Err(GeomError::NotManifold(format!(
                    "on-plane mesh edge {:?} has {} incident triangles",
                    edge,
                    incidences.len()
                )))
            }
        }
    }

    let in_plane_triangles = in_plane.len();
    let region_vertices: BTreeSet<u32> = in_plane
        .iter()
        .flat_map(|triangle| triangle.vertices)
        .collect();
    let set = InPlaneSet::new(in_plane)?;
    let in_plane_pairs = check_overlaps(
        mesh,
        projection,
        &set,
        request.region_limits.max_in_plane_pairs,
        options,
    )?;
    let segment_vertices: BTreeSet<u32> = segments
        .iter()
        .flat_map(|&Segment(left, right)| [left, right])
        .filter_map(|node| match node {
            NodeKey::Vertex(vertex) => Some(vertex),
            NodeKey::Edge(_) => None,
        })
        .collect();
    let regions = build_regions(
        mesh,
        frame,
        projection,
        &set,
        &segment_vertices,
        request.region_limits,
        options,
    )?;
    let chains = assemble_chains(
        mesh,
        frame,
        &classifications,
        &segments,
        &region_vertices,
        limits,
        options,
    )?;

    let evidence = PlaneSectionEvidence::input_mesh(mesh.triangle_count())
        .with_source_counts(in_plane_triangles, touching_edges, in_plane_pairs)
        .with_output(&chains, &regions);
    Ok(PlaneSection::new(frame, chains, regions, evidence))
}
