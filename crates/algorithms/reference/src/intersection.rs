//! Intersection segments between two triangle meshes, and the polylines they
//! assemble into.
//!
//! # The capability this adds
//!
//! [`ScalarBoolean`] is exact and total for disjoint,
//! nested, and identical operands, and refuses everything else. The single
//! missing capability is resolving surfaces that properly cross, which needs
//! the intersection curve first. This module computes it.
//!
//! # Why nodes are symbolic, not coordinates
//!
//! An intersection point is named by the *source topology that produced it*
//! -- a vertex index, or the pair of faces and the edge that crossed -- never
//! by its computed position. Two faces sharing an edge then produce byte-
//! identical node names, so stitching segments into a polyline is exact
//! integer matching with no tolerance anywhere.
//!
//! Matching on coordinates instead would need an epsilon, and an epsilon in
//! the stitching step is precisely how a boolean develops cracks: two
//! segments that should share an endpoint fail to join, and the curve opens.
//! `ScalarSection` already uses this technique for plane cuts; this reuses it
//! for the mesh-mesh case.
//!
//! # Honest limits
//!
//! Coplanar face pairs are refused, not approximated. Their intersection is
//! an area rather than a curve, and resolving it needs a 2D overlap policy
//! the caller must choose. Refusing keeps this module's output meaning
//! exactly one thing.

use std::collections::{BTreeMap, BTreeSet};

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation, Sign};
use axiolid_core::Point3;
use axiolid_mesh::TriMesh;

use crate::boolean::ScalarBoolean;
use crate::{orient3d, triangle_triangle_relation, TriangleTriangleRelation};

/// Which operand a piece of source topology belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Operand {
    /// The first mesh passed to [`intersection_segments`].
    Subject,
    /// The second mesh passed to [`intersection_segments`].
    Tool,
}

/// An undirected mesh edge, named by its two vertex indices in sorted order.
///
/// Sorting is what makes the name canonical: the two faces sharing an edge
/// visit it in opposite directions, and both must produce the same key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeKey {
    operand: Operand,
    low: u32,
    high: u32,
}

impl EdgeKey {
    fn new(operand: Operand, first: u32, second: u32) -> Self {
        Self {
            operand,
            low: first.min(second),
            high: first.max(second),
        }
    }
}

/// An intersection point, named by the source topology that produced it.
///
/// Never by coordinates: see the module docs for why that distinction
/// decides whether the assembled curve can crack.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NodeKey {
    /// An original mesh vertex lying exactly on the other surface.
    Vertex {
        /// Which operand owns the vertex.
        operand: Operand,
        /// Its index in that operand's position array.
        index: u32,
    },
    /// An edge of one operand crossing the other operand's surface.
    ///
    /// Identified by the edge ALONE, deliberately. The crossed triangle is
    /// not part of the name: a closed surface is triangulated arbitrarily,
    /// so one puncture point can sit on a shared triangle edge and be
    /// reported once per incident triangle. Including the face index would
    /// give that single point two names, and the curve would fragment into
    /// disconnected two-node pieces instead of closing into a loop.
    ///
    /// An edge crossing a plane it is NOT parallel to punctures it exactly
    /// once, so the pierced plane completes the name. The plane is identified
    /// by its own geometry rather than by a triangle index: a flat side of a
    /// solid is triangulated arbitrarily, and naming the triangle would give
    /// one puncture two names whenever it lands on a shared triangle edge.
    ///
    /// This matters for a through-cut. An edge that enters one side of a slab
    /// and leaves the other punctures the surface TWICE; identifying the node
    /// by the edge alone would collapse both punctures into a single name and
    /// corrupt the curve into degree-3 nodes.
    EdgeSurface {
        /// The crossing edge.
        edge: EdgeKey,
        /// Where along that edge the puncture lies, as exact coordinate bits.
        at: PointKey,
    },
}

/// A point named by its exact coordinate bits.
///
/// Used only to disambiguate several punctures along ONE edge. It is not a
/// coordinate-tolerance match: two punctures are the same node only when
/// their coordinates are bit-identical, which they are when the same
/// arithmetic produced them from the same inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PointKey {
    bits: [u64; 3],
}

impl PointKey {
    fn new(point: Point3) -> Self {
        Self {
            // `+ 0.0` folds `-0.0` into `0.0` so equal points share bits.
            bits: [point.x + 0.0, point.y + 0.0, point.z + 0.0].map(f64::to_bits),
        }
    }
}

/// One segment of the intersection curve, joining two nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct IntersectionSegment {
    /// One endpoint.
    pub start: NodeKey,
    /// The other endpoint.
    pub end: NodeKey,
}

impl IntersectionSegment {
    /// Build a segment, rejecting one that collapses to a single node.
    fn new(start: NodeKey, end: NodeKey) -> GeomResult<Self> {
        if start == end {
            return Err(GeomError::Degenerate(
                "an intersection segment collapsed to one source-topology node".into(),
            ));
        }
        // Canonical order, so a segment found twice compares equal.
        Ok(Self {
            start: start.min(end),
            end: start.max(end),
        })
    }
}

/// Exact orientation sign of `point` against the plane of `triangle`.
fn plane_sign(triangle: [Point3; 3], point: Point3) -> Sign {
    orient3d(triangle[0], triangle[1], triangle[2], point)
        .sign()
        .expect("certified predicates are total")
}

/// Nodes where one face's edges meet the other face's plane.
///
/// Returns at most two: a triangle is convex, so its boundary crosses a
/// plane in at most two places. A vertex exactly on the plane contributes a
/// `Vertex` node; an edge straddling it contributes an `EdgeFace` node.
///
/// This finds where the edges meet the *plane*, which is a superset of where
/// they meet the *triangle*. The caller filters to the triangle.
fn crossing_nodes(
    face_vertices: [u32; 3],
    face_operand: Operand,
    face_points: [Point3; 3],
    other: [Point3; 3],
) -> Vec<(NodeKey, Point3)> {
    let signs = face_points.map(|point| plane_sign(other, point));
    let mut nodes = Vec::new();

    for corner in 0..3 {
        let next = (corner + 1) % 3;

        // A vertex ON the plane is itself an intersection point, and is
        // named by its own index so both operands agree on it.
        if signs[corner] == Sign::Zero {
            nodes.push((
                NodeKey::Vertex {
                    operand: face_operand,
                    index: face_vertices[corner],
                },
                face_points[corner],
            ));
            continue;
        }

        // A straddling edge crosses once, strictly between its endpoints.
        // `Zero` at `next` is handled when that corner is visited, so only
        // a genuine sign flip counts here -- otherwise the point is emitted
        // twice under two different names.
        if signs[next] != Sign::Zero && signs[corner] != signs[next] {
            // Which puncture along this edge: the parameter of the crossing,
            // as exact bits. An edge that pierces the other surface several
            // times (a through-cut) gets a distinct name per puncture, while
            // the SAME puncture reported from several incident triangles of a
            // flat face gets one name, because the parameter is identical.
            let crossing = plane_crossing(face_points[corner], face_points[next], other);
            let key = NodeKey::EdgeSurface {
                edge: EdgeKey::new(face_operand, face_vertices[corner], face_vertices[next]),
                at: PointKey::new(crossing),
            };
            nodes.push((key, crossing));
        }
    }
    nodes
}

/// Where segment `start`->`end` meets the plane of `triangle`.
///
/// The caller has already PROVEN a crossing exists using exact predicates;
/// this only computes the position. That separation is deliberate: the
/// decision is exact, and the coordinate is the best binary64 approximation
/// of a point whose existence is certain. Rounding moves the point slightly,
/// never invents or removes it, because the node's identity comes from its
/// symbolic name rather than this value.
fn plane_crossing(start: Point3, end: Point3, triangle: [Point3; 3]) -> Point3 {
    let normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
    let start_height = (start - triangle[0]).dot(normal);
    let end_height = (end - triangle[0]).dot(normal);
    let span = start_height - end_height;
    if span == 0.0 {
        // Unreachable for a proven crossing: opposite exact signs cannot
        // produce equal heights. Returning the midpoint keeps the function
        // total rather than panicking inside a geometry kernel.
        return start.midpoint(end);
    }
    start + (end - start) * (start_height / span)
}

/// Whether `point` lies within `triangle`, given it is already on its plane.
///
/// Uses the same exact sign test as the rest of the module: the point is
/// inside when it is on the same side of all three edges, with `Zero`
/// accepted so boundary contact counts as inside.
fn point_in_triangle(point: Point3, triangle: [Point3; 3]) -> bool {
    let normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
    let mut positive = false;
    let mut negative = false;
    for corner in 0..3 {
        let next = (corner + 1) % 3;
        // Build a tetrahedron from the edge, the point, and the face normal;
        // its orientation says which side of the edge the point is on.
        let apex = triangle[corner] + normal;
        match plane_sign([triangle[corner], triangle[next], apex], point) {
            Sign::Positive => positive = true,
            Sign::Negative => negative = true,
            // `Zero` means the point is exactly on this edge's plane, which
            // is boundary contact and counts as inside. `Sign` is
            // non-exhaustive, so an unknown future variant is treated the
            // same rather than silently changing the answer.
            _ => {}
        }
    }
    !(positive && negative)
}

/// The intersection curve between two meshes, as segments plus positions.
#[derive(Debug, Clone, Default)]
pub struct IntersectionCurve {
    /// Segments, deduplicated and in canonical order.
    pub segments: Vec<IntersectionSegment>,
    /// Position of every node named by a segment.
    pub positions: BTreeMap<NodeKey, Point3>,
}

/// Compute the intersection curve between two triangle meshes.
///
/// `O(n*m)`: every face pair is tested. This is the reference implementation,
/// so it is written to be obviously right rather than fast -- a BVH here
/// would be a second thing to get wrong. A production provider adds one.
///
/// # Errors
///
/// Refuses coplanar face pairs, whose intersection is an area rather than a
/// curve and needs a 2D overlap policy the caller must choose.
pub fn intersection_segments(subject: &TriMesh, tool: &TriMesh) -> GeomResult<IntersectionCurve> {
    let mut segments = BTreeSet::new();
    let mut positions = BTreeMap::new();

    for subject_face in 0..subject.triangle_count() {
        let (subject_indices, subject_points) = face(subject, subject_face)?;
        for tool_face in 0..tool.triangle_count() {
            let (tool_indices, tool_points) = face(tool, tool_face)?;

            match triangle_triangle_relation(subject_points, tool_points) {
                // No shared point: nothing to record.
                TriangleTriangleRelation::Disjoint => continue,
                // An area, not a curve. Refuse rather than pick a policy.
                TriangleTriangleRelation::Coplanar => {
                    // An area, not a curve: resolving it needs a 2D overlap
                    // policy the caller must choose, so refuse by contract.
                    return Err(GeomError::Unsupported {
                        backend: BackendId::new("scalar-intersection"),
                        operation: Operation::MeshBoolean,
                    });
                }
                // A degenerate source face has no well-defined plane.
                TriangleTriangleRelation::DegenerateTriangle => {
                    return Err(GeomError::Degenerate(
                        "a source face is degenerate; heal the mesh before intersecting".into(),
                    ))
                }
                // Both contribute segments: `Touching` includes edge-on-face
                // contact, which is a real part of the curve.
                TriangleTriangleRelation::Proper | TriangleTriangleRelation::Touching => {}
            }

            let subject_normal = (subject_points[1] - subject_points[0])
                .cross(subject_points[2] - subject_points[0]);
            let tool_normal =
                (tool_points[1] - tool_points[0]).cross(tool_points[2] - tool_points[0]);

            // Nodes contributed by each face crossing the other's plane,
            // filtered to those actually inside the other triangle.
            let mut nodes = Vec::new();
            for (key, point) in crossing_nodes(
                subject_indices,
                Operand::Subject,
                subject_points,
                tool_points,
            ) {
                if point_in_triangle(point, tool_points) {
                    nodes.push((key, point));
                }
            }
            for (key, point) in
                crossing_nodes(tool_indices, Operand::Tool, tool_points, subject_points)
            {
                if point_in_triangle(point, subject_points) {
                    nodes.push((key, point));
                }
            }

            // The two triangles' planes meet in a line; each triangle clips
            // that line to an interval, and the curve here is the OVERLAP of
            // those two intervals.
            //
            // Up to four nodes arrive: each operand's edges can puncture the
            // other's triangle. Four is the normal transverse case, not an
            // error -- discarding it was what cracked the curve into
            // disconnected two-node pieces.
            //
            // Ordering the nodes ALONG the intersection line and taking the
            // middle two yields exactly the shared interval: the outer two
            // are each outside the other triangle.
            nodes.sort_by(|left, right| left.0.cmp(&right.0));
            nodes.dedup_by(|left, right| left.0 == right.0);
            if nodes.len() < 2 {
                // A single point of contact contributes no segment.
                for (key, point) in nodes {
                    positions.insert(key, point);
                }
                continue;
            }

            // Direction of the plane-plane intersection line.
            let axis = subject_normal.cross(tool_normal);
            if axis.length_squared() == 0.0 {
                // Parallel planes that are not coplanar cannot cross; the
                // coplanar case was refused above.
                for (key, point) in nodes {
                    positions.insert(key, point);
                }
                continue;
            }
            nodes.sort_by(|left, right| {
                let left_t = left.1.dot(axis);
                let right_t = right.1.dot(axis);
                left_t
                    .partial_cmp(&right_t)
                    .expect("finite coordinates give an orderable projection")
            });
            let interval = if nodes.len() == 2 {
                [nodes[0], nodes[1]]
            } else {
                [nodes[nodes.len() / 2 - 1], nodes[nodes.len() / 2]]
            };
            for (key, point) in nodes.iter().copied() {
                positions.insert(key, point);
            }
            if interval[0].0 == interval[1].0 {
                continue;
            }

            let segment = IntersectionSegment::new(interval[0].0, interval[1].0)?;
            segments.insert(segment);
        }
    }

    Ok(IntersectionCurve {
        segments: segments.into_iter().collect(),
        positions,
    })
}

/// A connected run of the intersection curve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polyline {
    /// Nodes in traversal order.
    pub nodes: Vec<NodeKey>,
    /// Whether the run returns to its first node.
    ///
    /// A closed loop is the normal result for two closed solids. An open
    /// run means the curve reached a boundary, which a closed operand
    /// should not have.
    pub closed: bool,
}

/// Stitch segments into maximal connected polylines.
///
/// Pure integer graph traversal: nodes are symbolic names, so joining two
/// segments is an equality test rather than a distance comparison. No
/// tolerance is involved at any point.
///
/// # Errors
///
/// Refuses a node of degree three or more. On a clean pair of closed
/// surfaces the intersection curve is a 1-manifold, so every node has one
/// or two neighbours; a branch means the input is self-intersecting or
/// non-manifold, and continuing would silently pick one arbitrary path.
pub fn assemble_polylines(segments: &[IntersectionSegment]) -> GeomResult<Vec<Polyline>> {
    let mut adjacency: BTreeMap<NodeKey, Vec<NodeKey>> = BTreeMap::new();
    for segment in segments {
        adjacency
            .entry(segment.start)
            .or_default()
            .push(segment.end);
        adjacency
            .entry(segment.end)
            .or_default()
            .push(segment.start);
    }
    for (node, neighbours) in &mut adjacency {
        neighbours.sort_unstable();
        neighbours.dedup();
        if neighbours.len() > 2 {
            return Err(GeomError::NotManifold(format!(
                "intersection curve branches at {node:?} with degree {}",
                neighbours.len()
            )));
        }
        // The intersection of two CLOSED surfaces is a set of closed rings, so
        // every node must have exactly two neighbours. A degree-1 node means
        // the curve was cut short -- in practice because an edge of one
        // operand crosses an edge of the other at the same point, and each
        // operand named that puncture from its own side. The two names do not
        // join, so the ring opens.
        //
        // Refuse rather than return the broken run. Merging coincident names
        // by coordinate was tried and rejected: it welded genuinely distinct
        // nodes in the corner-overlap case, replacing a visible failure with a
        // quiet wrong answer.
        if neighbours.len() < 2 {
            return Err(GeomError::Unsupported {
                backend: ScalarBoolean::ID,
                operation: Operation::MeshBoolean,
            });
        }
    }

    let mut visited = BTreeSet::new();
    let mut polylines = Vec::new();

    // Open runs first: starting from a degree-1 node walks the whole run in
    // one pass. Starting mid-run would produce two half-runs instead.
    let endpoints: Vec<NodeKey> = adjacency
        .iter()
        .filter(|(_, neighbours)| neighbours.len() == 1)
        .map(|(node, _)| *node)
        .collect();
    for start in endpoints {
        if visited.contains(&start) {
            continue;
        }
        polylines.push(walk(start, &adjacency, &mut visited, false));
    }

    // Whatever remains is a cycle: every node has degree two.
    let cycle_starts: Vec<NodeKey> = adjacency.keys().copied().collect();
    for start in cycle_starts {
        if visited.contains(&start) {
            continue;
        }
        polylines.push(walk(start, &adjacency, &mut visited, true));
    }

    Ok(polylines)
}

/// Walk one connected run from `start`, marking nodes visited.
fn walk(
    start: NodeKey,
    adjacency: &BTreeMap<NodeKey, Vec<NodeKey>>,
    visited: &mut BTreeSet<NodeKey>,
    closed: bool,
) -> Polyline {
    let mut nodes = vec![start];
    visited.insert(start);
    let mut current = start;
    let mut previous = None;

    loop {
        let Some(neighbours) = adjacency.get(&current) else {
            break;
        };
        // Step to the neighbour we did not arrive from. On a cycle both are
        // unvisited at the first step, so the choice of direction is
        // arbitrary but consistent -- `adjacency` is sorted.
        let next = neighbours
            .iter()
            .copied()
            .find(|candidate| Some(*candidate) != previous && !visited.contains(candidate));
        let Some(next) = next else {
            break;
        };
        nodes.push(next);
        visited.insert(next);
        previous = Some(current);
        current = next;
    }

    Polyline { nodes, closed }
}

/// Vertex indices and positions of one face.
fn face(mesh: &TriMesh, index: usize) -> GeomResult<([u32; 3], [Point3; 3])> {
    let base = index * 3;
    let indices: [u32; 3] = mesh
        .indices
        .get(base..base + 3)
        .ok_or_else(|| GeomError::Degenerate("face index out of range".into()))?
        .try_into()
        .map_err(|_| GeomError::Degenerate("face index slice is not three wide".into()))?;
    let mut points = [Point3::ZERO; 3];
    for (slot, vertex) in indices.iter().enumerate() {
        points[slot] = *mesh
            .positions
            .get(*vertex as usize)
            .ok_or_else(|| GeomError::Degenerate("face references a missing vertex".into()))?;
    }
    Ok((indices, points))
}
