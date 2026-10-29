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
//!
//! A `SectionRole::Sheet` may be open. A mesh edge lying in the plane with one
//! incident triangle is then either a region boundary (the triangle is in the
//! plane) or a touching edge (it is off the plane), and a chain can end on a
//! boundary edge or boundary vertex of the mesh. A solid refuses an edge with
//! one incident triangle.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation, Sign};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    PlaneSection, PlaneSectionEvidence, SectionRequest, SectionRole,
};

use super::chains::{assemble_chains, ChainContext};
use super::regions::{build_regions, check_overlaps, InPlaneSet, InPlaneTriangle, Projection};
use super::{
    check_source_limits, classify_vertices, insert_segment, internal_topology_error, mesh_index,
    opposite, EdgeKey, ExactSectionPlane, NodeKey, ScalarSection, Segment,
};
use crate::orient3::orient3d;

/// A triangle incident to a mesh edge that lies in the plane.
#[derive(Debug, Clone, Copy)]
enum Incidence {
    /// An off-plane triangle with its third vertex on this side; the cut along
    /// the edge runs from the first vertex to the second.
    Off(Sign, u32, u32),
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
    let sheet = match request.role {
        SectionRole::Solid => false,
        SectionRole::Sheet => true,
        _ => {
            return Err(GeomError::Unsupported {
                backend: ScalarSection::ID,
                operation: Operation::MeshPlaneSection,
            })
        }
    };
    let boundary = if sheet {
        Boundary::of(mesh, options)?
    } else {
        Boundary::default()
    };

    let mut segments = BTreeSet::<Segment>::new();
    let mut on_plane_edges = BTreeMap::<EdgeKey, Vec<Incidence>>::new();
    let mut in_plane = Vec::new();
    // Per segment, whether the cut runs from its first node to its second with
    // the solid on its left. The triangles are outward oriented, so the cut
    // leaves a triangle's boundary where the winding goes from above the plane
    // (negative `orient3d`) to below it, and enters where it goes back up.
    let mut directions = BTreeMap::<Segment, bool>::new();
    // Third vertices of the triangles on each crossing edge.
    let mut crossing_thirds = BTreeMap::<EdgeKey, Vec<u32>>::new();
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
                let corner = (0..3)
                    .find(|&corner| signs[corner] != Sign::Zero)
                    .ok_or_else(internal_topology_error)?;
                let third = signs[corner];
                // The winding order of the two on-plane vertices. Above the
                // plane the cut follows it; below, it runs against it.
                let (first, second) = (triangle[(corner + 1) % 3], triangle[(corner + 2) % 3]);
                let (from, to) = if third == Sign::Negative {
                    (first, second)
                } else {
                    (second, first)
                };
                on_plane_edges
                    .entry(EdgeKey::new(left, right))
                    .or_default()
                    .push(Incidence::Off(third, from, to));
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
                    let vertex = NodeKey::Vertex(triangle[zero_corner]);
                    let edge = EdgeKey::new(triangle[first], triangle[second]);
                    let segment = Segment::new(vertex, NodeKey::Edge(edge))?;
                    insert_segment(&mut segments, segment, limits)?;
                    // The winding runs zero corner, next, previous; the edge
                    // is crossed from the next to the previous vertex.
                    let next = (zero_corner + 1) % 3;
                    if signs[next] == Sign::Negative {
                        record_direction(&mut directions, NodeKey::Edge(edge), vertex)?;
                    } else {
                        record_direction(&mut directions, vertex, NodeKey::Edge(edge))?;
                    }
                    crossing_thirds
                        .entry(edge)
                        .or_default()
                        .push(triangle[zero_corner]);
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
                        let edge = EdgeKey::new(triangle[left], triangle[right]);
                        // The cut starts on the edge that goes from above to below.
                        crossing[crossing_count] =
                            Some((NodeKey::Edge(edge), signs[left] == Sign::Negative));
                        crossing_thirds
                            .entry(edge)
                            .or_default()
                            .push(triangle[3 - left - right]);
                        crossing_count += 1;
                    }
                }
                match (crossing[0], crossing[1]) {
                    (Some((first, first_starts)), Some((second, _))) => {
                        insert_segment(&mut segments, Segment::new(first, second)?, limits)?;
                        if first_starts {
                            record_direction(&mut directions, first, second)?;
                        } else {
                            record_direction(&mut directions, second, first)?;
                        }
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
            [Incidence::Off(first, from, to), Incidence::Off(second, _, _)] => {
                if opposite(*first, *second) {
                    insert_segment(
                        &mut segments,
                        Segment::new(NodeKey::Vertex(edge.0), NodeKey::Vertex(edge.1))?,
                        limits,
                    )?;
                    record_direction(
                        &mut directions,
                        NodeKey::Vertex(*from),
                        NodeKey::Vertex(*to),
                    )?;
                } else {
                    touching_edges += 1;
                }
            }
            [Incidence::In, Incidence::Off(..)] | [Incidence::Off(..), Incidence::In] => {}
            [Incidence::In, Incidence::In] => {}
            // Boundary edges of a sheet: a region boundary, or an edge the
            // sheet touches the plane in.
            [Incidence::In] if sheet => {}
            [Incidence::Off(..)] if sheet => touching_edges += 1,
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
    let coplanar_edges = coplanar_crossing_edges(mesh, &crossing_thirds)?;
    let chains = assemble_chains(
        &ChainContext {
            mesh,
            frame,
            classifications: &classifications,
            projection,
            directions: &directions,
            coplanar_edges: &coplanar_edges,
            region_vertices: &region_vertices,
            role: request.role,
            boundary_edges: &boundary.edges,
            boundary_vertices: &boundary.vertices,
        },
        &segments,
        limits,
        options,
    )?;

    let evidence = PlaneSectionEvidence::input_mesh(mesh.triangle_count())
        .with_source_counts(in_plane_triangles, touching_edges, in_plane_pairs)
        .with_role(request.role, boundary.edges.len())
        .with_output(&chains, &regions);
    Ok(PlaneSection::new(frame, chains, regions, evidence))
}

/// The boundary of a sheet: edges used by one triangle, and their vertices.
#[derive(Debug, Default)]
struct Boundary {
    edges: BTreeSet<EdgeKey>,
    vertices: BTreeSet<u32>,
}

impl Boundary {
    fn of(mesh: &TriMesh, options: &ExecutionOptions) -> GeomResult<Self> {
        let mut uses = BTreeMap::<EdgeKey, usize>::new();
        for triangle in mesh.triangles() {
            options.check_cancelled()?;
            for (from, to) in [(0, 1), (1, 2), (2, 0)] {
                *uses
                    .entry(EdgeKey::new(triangle[from], triangle[to]))
                    .or_default() += 1;
            }
        }
        let mut boundary = Self::default();
        for (edge, count) in uses {
            match count {
                1 => {
                    boundary.edges.insert(edge);
                    boundary.vertices.extend([edge.0, edge.1]);
                }
                2 => {}
                _ => {
                    return Err(GeomError::NotManifold(format!(
                        "mesh edge {edge:?} has {count} incident triangles"
                    )))
                }
            }
        }
        Ok(boundary)
    }
}

/// Record that the cut runs from `from` to `to`.
fn record_direction(
    directions: &mut BTreeMap<Segment, bool>,
    from: NodeKey,
    to: NodeKey,
) -> GeomResult<()> {
    let segment = Segment::new(from, to)?;
    directions.insert(segment, from == segment.0);
    Ok(())
}

/// Crossing edges whose two triangles are exactly coplanar: the four vertices
/// have a zero `orient3d`.
fn coplanar_crossing_edges(
    mesh: &TriMesh,
    thirds: &BTreeMap<EdgeKey, Vec<u32>>,
) -> GeomResult<BTreeSet<EdgeKey>> {
    let mut coplanar = BTreeSet::new();
    for (edge, thirds) in thirds {
        let [first, second] = thirds.as_slice() else {
            continue;
        };
        let point = |index: u32| -> GeomResult<_> { Ok(mesh.positions[mesh_index(index)?]) };
        let certified = orient3d(
            point(edge.0)?,
            point(edge.1)?,
            point(*first)?,
            point(*second)?,
        );
        if certified.sign() == Some(Sign::Zero) {
            coplanar.insert(*edge);
        }
    }
    Ok(coplanar)
}
