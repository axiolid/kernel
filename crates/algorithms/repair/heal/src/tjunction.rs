//! Splitting triangle edges at vertices that lie on them (T-junctions).
//!
//! A cut opening leaves a long edge on one face matched by two short edges
//! on the other: the mesh looks closed, but no edge is shared, so the
//! oriented-cycle audit counts three unbalanced edges and no volume is
//! measured. Splitting the long side at the vertex restores the pairing.
//!
//! Only edges the audit finds unbalanced are split, and only at vertices
//! that are endpoints of such an edge, so a closed mesh is never touched.
//! A vertex lies on an edge when it is strictly inside it (`0 < t < 1`),
//! at most `tolerance.linear()` from its line and further than
//! `tolerance.linear()` from both endpoints. No vertex is created or moved,
//! so per-vertex data is untouched; corner data follows the split corners.
//!
//! The split is declined, and the triangle kept as it was, when it would
//! give a wrong answer: the triangle's own third corner lies on its edge, a
//! piece would be degenerate or face against the parent, or the two ends of
//! a split side hold different values in a corner-indexed buffer, so no
//! single value could be given to the new corner, or a per-vertex value at
//! the inserted vertex is not what interpolating its side's endpoints gives.
//! Declines are counted.
//!
//! Recommended pipeline before `axiolid_measure::cycle_volume_properties`:
//! weld vertices, drop degenerate elements, split T-junctions, with the same
//! tolerance.

use std::collections::{BTreeMap, BTreeSet};
use std::ops::ControlFlow;

use axiolid_core::{Aabb, Point3, Tolerance};
use axiolid_mesh::{split_triangle_at_side_corners, TriMesh};
use axiolid_spatial::{Bvh, SpatialIndex, SpatialItem};

use crate::carry;
use crate::mesh::MeshHealError;

/// What [`split_t_junctions`] did.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct TJunctionSplit {
    /// Vertex-on-edge incidences found on unbalanced edges, counted once per
    /// undirected edge however many triangles use it.
    pub t_junctions: usize,
    /// Triangles that were split.
    pub split_triangles: usize,
    /// Triangles added: each split triangle's pieces minus one.
    pub triangles_added: usize,
    /// Triangles with a vertex on an edge that were kept unsplit.
    pub declined_triangles: usize,
    /// Largest distance from a used vertex to the edge it was placed on.
    /// The split triangles stay within this (at most the tolerance) of the
    /// originals; it is reported, never refused.
    pub max_deviation: f64,
}

/// Split triangle edges at vertices lying on them within `tolerance`.
///
/// The first piece of a split triangle keeps its index; the other pieces are
/// appended in ascending parent order, so triangle indices stay stable.
/// Returns the mesh unchanged when nothing is split.
///
/// # Errors
/// [`MeshHealError::RaggedIndices`] when the index buffer is not a whole
/// number of triangles.
pub fn split_t_junctions(
    mesh: &TriMesh,
    tolerance: Tolerance,
) -> Result<(TriMesh, TJunctionSplit), MeshHealError> {
    if !mesh.indices.len().is_multiple_of(3) {
        return Err(MeshHealError::RaggedIndices {
            indices: mesh.indices.len(),
        });
    }
    let linear = tolerance.linear();
    let limit = linear.powi(4);
    let count = mesh.indices.len() / 3;
    let mut report = TJunctionSplit {
        t_junctions: 0,
        split_triangles: 0,
        triangles_added: 0,
        declined_triangles: 0,
        max_deviation: 0.0,
    };

    // The audit's usability: valid indices, finite positions, and a doubled
    // area above the tolerance (compared squared).
    let usable: Vec<Option<[u32; 3]>> = (0..count)
        .map(|t| {
            let ids = [
                mesh.indices[t * 3],
                mesh.indices[t * 3 + 1],
                mesh.indices[t * 3 + 2],
            ];
            let [a, b, c] = ids.map(|i| mesh.positions.get(i as usize).copied());
            let (a, b, c) = (a?, b?, c?);
            if !(a.is_finite() && b.is_finite() && c.is_finite()) {
                return None;
            }
            ((b - a).cross(c - a).length_squared() > limit).then_some(ids)
        })
        .collect();

    // Net direction per undirected edge over usable triangles.
    let mut net: BTreeMap<(u32, u32), i64> = BTreeMap::new();
    for ids in usable.iter().flatten() {
        for s in 0..3 {
            let (u, v) = (ids[s], ids[(s + 1) % 3]);
            *net.entry((u.min(v), u.max(v))).or_default() += if u < v { 1 } else { -1 };
        }
    }
    let unbalanced: Vec<(u32, u32)> = net
        .iter()
        .filter(|&(_, &winding)| winding != 0)
        .map(|(&edge, _)| edge)
        .collect();
    let candidates: BTreeSet<u32> = unbalanced.iter().flat_map(|&(a, b)| [a, b]).collect();
    if candidates.is_empty() {
        return Ok((mesh.clone(), report));
    }

    let index = Bvh::build(
        candidates
            .iter()
            .map(|&v| SpatialItem::new(v, Aabb::from_point(mesh.positions[v as usize]))),
    );

    // Vertices on each unbalanced edge, with their distance from it, in
    // ascending (parameter, vertex index) order from the low endpoint.
    let mut on_edge: BTreeMap<(u32, u32), Vec<(u32, f64)>> = BTreeMap::new();
    // Parameter of each vertex along its edge, from the low endpoint.
    let mut parameter: BTreeMap<(u32, u32, u32), f64> = BTreeMap::new();
    for &(a, b) in &unbalanced {
        let (pa, pb) = (mesh.positions[a as usize], mesh.positions[b as usize]);
        let mut probe = Aabb::from_point(pa);
        probe.extend(pb);
        probe.min = Point3::new(
            probe.min.x - linear,
            probe.min.y - linear,
            probe.min.z - linear,
        );
        probe.max = Point3::new(
            probe.max.x + linear,
            probe.max.y + linear,
            probe.max.z + linear,
        );
        let ab = pb - pa;
        let length_squared = ab.length_squared();
        let mut found: Vec<(f64, u32, f64)> = Vec::new();
        index.visit_aabb(&probe, &mut |&m: &u32| {
            if m != a && m != b {
                let pm = mesh.positions[m as usize];
                let t = (pm - pa).dot(ab) / length_squared;
                let distance = ab.cross(pm - pa).length() / length_squared.sqrt();
                if t > 0.0
                    && t < 1.0
                    && distance <= linear
                    && (pm - pa).length() > linear
                    && (pm - pb).length() > linear
                {
                    found.push((t, m, distance));
                }
            }
            ControlFlow::Continue(())
        });
        if found.is_empty() {
            continue;
        }
        found.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        report.t_junctions += found.len();
        for &(t, m, _) in &found {
            parameter.insert((a, b, m), t);
        }
        on_edge.insert((a, b), found.into_iter().map(|(_, m, d)| (m, d)).collect());
    }
    if on_edge.is_empty() {
        return Ok((mesh.clone(), report));
    }

    let aligned = carry::corner_buffers_aligned(mesh);
    let mut indices = mesh.indices.clone();
    let mut source: Vec<usize> = (0..indices.len()).collect();
    for (t, ids) in usable.iter().enumerate() {
        let Some(ids) = *ids else {
            continue;
        };
        // Points on side `s`, from `ids[s]` to `ids[s + 1]`.
        let sides: [Vec<(u32, f64)>; 3] = std::array::from_fn(|s| {
            let (u, v) = (ids[s], ids[(s + 1) % 3]);
            let mut points = on_edge
                .get(&(u.min(v), u.max(v)))
                .cloned()
                .unwrap_or_default();
            if u > v {
                points.reverse();
            }
            points
        });
        if sides.iter().all(Vec::is_empty) {
            continue;
        }
        let own_corner_on_edge =
            (0..3).any(|s| sides[s].iter().any(|&(m, _)| m == ids[(s + 2) % 3]));
        let corners_split_cleanly = aligned
            && (0..3).all(|s| {
                sides[s].is_empty() || carry::corners_agree(mesh, t * 3 + s, t * 3 + (s + 1) % 3)
            });
        // New corners reuse a vertex that already carries its own per-vertex
        // data, so that data must be what interpolation would have given.
        let vertex_data_consistent = (0..3).all(|s| {
            let (u, v) = (ids[s], ids[(s + 1) % 3]);
            let (low, high) = (u.min(v), u.max(v));
            sides[s].iter().all(|&(m, _)| {
                parameter
                    .get(&(low, high, m))
                    .is_some_and(|&at| carry::vertex_data_agrees(mesh, low, high, m, at))
            })
        });
        let pieces = (!own_corner_on_edge && corners_split_cleanly && vertex_data_consistent)
            .then(|| triangulate(t, ids, &sides, mesh, limit))
            .flatten();
        let Some(Pieces {
            vertices,
            corners,
            tokens,
        }) = pieces
        else {
            report.declined_triangles += 1;
            continue;
        };
        let pieces = tokens.len() / 3;
        for (slot, &token) in tokens.iter().enumerate() {
            let (vertex, corner) = (vertices[token], corners[token]);
            if slot < 3 {
                indices[t * 3 + slot] = vertex;
                source[t * 3 + slot] = corner;
            } else {
                indices.push(vertex);
                source.push(corner);
            }
        }
        report.split_triangles += 1;
        report.triangles_added += pieces - 1;
        for &(_, distance) in sides.iter().flatten() {
            report.max_deviation = report.max_deviation.max(distance);
        }
    }
    if report.split_triangles == 0 {
        return Ok((mesh.clone(), report));
    }
    let mut out = mesh.clone();
    out.indices = indices;
    carry::remap_corners(&mut out, &source);
    Ok((out, report))
}

/// A split triangle as tokens: `0..3` are its corners, the rest its side
/// points in side order.
struct Pieces {
    /// Vertex of every token.
    vertices: Vec<u32>,
    /// Source corner of every token; a side point takes the corner its side
    /// starts at.
    corners: Vec<usize>,
    /// The pieces as tokens, three per piece.
    tokens: Vec<usize>,
}

/// Triangulate triangle `t` with its side points; `None` when a piece would
/// be degenerate or face against the parent.
fn triangulate(
    t: usize,
    ids: [u32; 3],
    sides: &[Vec<(u32, f64)>; 3],
    mesh: &TriMesh,
    limit: f64,
) -> Option<Pieces> {
    let mut vertices = ids.to_vec();
    let mut corners: Vec<usize> = (0..3).map(|i| t * 3 + i).collect();
    let mut side_tokens: [Vec<usize>; 3] = Default::default();
    for (s, points) in sides.iter().enumerate() {
        for &(m, _) in points {
            side_tokens[s].push(vertices.len());
            vertices.push(m);
            corners.push(t * 3 + s);
        }
    }
    let mut tokens = Vec::new();
    if !split_triangle_at_side_corners([0, 1, 2], &side_tokens, &mut tokens) || tokens.len() <= 3 {
        return None;
    }
    let point = |token: usize| mesh.positions[vertices[token] as usize];
    let parent = (point(1) - point(0)).cross(point(2) - point(0));
    let sound = tokens.chunks(3).all(|piece| {
        let normal = (point(piece[1]) - point(piece[0])).cross(point(piece[2]) - point(piece[0]));
        normal.length_squared() > limit && normal.dot(parent) > 0.0
    });
    sound.then_some(Pieces {
        vertices,
        corners,
        tokens,
    })
}
