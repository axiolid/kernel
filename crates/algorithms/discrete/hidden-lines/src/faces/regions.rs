//! Regions from the labelled arrangement.
//!
//! A piece is a boundary piece when its two labels differ. For each source
//! `K` its half-edges are the boundary pieces with `K` on exactly one side,
//! directed so that `K` is on the left. `trace` closes them into rings and
//! nests holes. A ring vertex is kept unless exactly two boundary pieces, of
//! any source, meet there in a straight line, so neighbouring regions share
//! bit-identical vertices.
//!
//! The triangles of a region are those frontmost on its sides. The sides of
//! `K` are joined into components: the two sides of a wedge at a vertex (one
//! face of the arrangement) and the two sides of a piece that carries `K` on
//! both. A component without a boundary half-edge is a cluster that
//! dissolved inside a region, and is placed by exact point location.

use std::collections::BTreeMap;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};
use axiolid_core::Point2;
use axiolid_exact_arrangement::{
    canonical_rotation, sequence_order, trace, Arrangement, ExactPoints, HalfEdge, Location,
    PieceId, PointId, TracedRing, VertexId,
};
use axiolid_mesh_hidden_lines_contract::{FaceSource, VisibleRegion};

use super::label::{Key, Label};
use super::{Error, Faces, Work};

/// Limits on the output.
#[derive(Debug, Clone, Copy)]
pub(crate) struct OutputLimits {
    pub max_regions: usize,
    pub max_region_vertices: usize,
}

/// The side `left` (`false`: right) of a piece, as an index.
fn side(piece: usize, left: bool) -> usize {
    2 * piece + usize::from(!left)
}

/// Union-find over a sorted list of side indices.
struct Components {
    nodes: Vec<usize>,
    parent: Vec<usize>,
}

impl Components {
    fn new(mut nodes: Vec<usize>) -> Self {
        nodes.sort_unstable();
        nodes.dedup();
        let parent = (0..nodes.len()).collect();
        Self { nodes, parent }
    }

    fn index(&self, node: usize) -> usize {
        self.nodes
            .binary_search(&node)
            .expect("every side of a source is a node")
    }

    fn find(&mut self, mut at: usize) -> usize {
        while self.parent[at] != at {
            self.parent[at] = self.parent[self.parent[at]];
            at = self.parent[at];
        }
        at
    }

    fn join(&mut self, a: usize, b: usize) {
        let (a, b) = (self.index(a), self.index(b));
        let (a, b) = (self.find(a), self.find(b));
        if a != b {
            self.parent[a.max(b)] = a.min(b);
        }
    }
}

/// The faces' points of one traced ring.
fn ring_points(arrangement: &Arrangement, ring: &TracedRing) -> Vec<PointId> {
    ring.vertices
        .iter()
        .map(|vertex| arrangement.vertices[vertex.0].point)
        .collect()
}

/// Whether the two boundary pieces at `at` run in a straight line through it.
fn straight(
    points: &ExactPoints,
    before: PointId,
    at: PointId,
    after: PointId,
) -> GeomResult<bool> {
    Ok(points.orient(before, at, after)? == Sign::Zero
        && points.strictly_between(before, at, after)?)
}

/// The vertices a ring must keep: where the number of boundary pieces is not
/// two, or two that do not run in a straight line.
fn kept_vertices(
    points: &ExactPoints,
    arrangement: &Arrangement,
    boundary: &[bool],
) -> GeomResult<Vec<bool>> {
    let mut around: Vec<Vec<VertexId>> = vec![Vec::new(); arrangement.vertices.len()];
    for (index, piece) in arrangement.pieces.iter().enumerate() {
        if boundary[index] {
            around[piece.ends[0].0].push(piece.ends[1]);
            around[piece.ends[1].0].push(piece.ends[0]);
        }
    }
    let mut keep = vec![true; around.len()];
    for (vertex, others) in around.iter().enumerate() {
        if let [a, b] = others[..] {
            let point = |v: VertexId| arrangement.vertices[v.0].point;
            keep[vertex] = !straight(
                points,
                point(a),
                arrangement.vertices[vertex].point,
                point(b),
            )?;
        }
    }
    Ok(keep)
}

/// The inconsistency of labels inside one face of the arrangement.
fn disagreement(faces: &Faces<'_>, a: &Label, b: &Label) -> Error {
    let message = GeomError::Degenerate("visible-face labels disagree inside a face".into());
    match a.keys.first().or(b.keys.first()) {
        Some(&(mesh, _)) => Error::mesh(faces.ids[mesh], message),
        None => Error::from(message),
    }
}

/// The regions of every source.
pub(crate) fn build(
    faces: &Faces<'_>,
    points: &ExactPoints,
    arrangement: &Arrangement,
    labels: &[[Label; 2]],
    limits: OutputLimits,
    work: &mut Work,
    options: &ExecutionOptions,
) -> Result<Vec<VisibleRegion>, Error> {
    let boundary: Vec<bool> = labels.iter().map(|[l, r]| l.keys != r.keys).collect();
    let keep = kept_vertices(points, arrangement, &boundary)?;

    // Half-edges per source, with the source on their left.
    let mut half_edges: BTreeMap<Key, Vec<HalfEdge>> = BTreeMap::new();
    for (index, [left, right]) in labels.iter().enumerate() {
        if !boundary[index] {
            continue;
        }
        let piece = PieceId(index);
        for key in &left.keys {
            if right.keys.binary_search(key).is_err() {
                half_edges.entry(*key).or_default().push(HalfEdge {
                    piece,
                    forward: true,
                });
            }
        }
        for key in &right.keys {
            if left.keys.binary_search(key).is_err() {
                half_edges.entry(*key).or_default().push(HalfEdge {
                    piece,
                    forward: false,
                });
            }
        }
    }

    // The joins of sides, per source.
    let label_of = |at: usize| &labels[at / 2][at % 2];
    let mut joins: BTreeMap<Key, Vec<(usize, usize)>> = BTreeMap::new();
    for vertex in 0..arrangement.vertices.len() {
        let fan = arrangement.fan(points, VertexId(vertex))?;
        for k in 0..fan.len() {
            let ((p, p_out), (q, q_out)) = (fan[k], fan[(k + 1) % fan.len()]);
            // The counter-clockwise side of `p` and the clockwise side of `q`
            // face one wedge.
            let a = side(p.0, p_out);
            let b = side(q.0, !q_out);
            let (la, lb) = (label_of(a), label_of(b));
            if la.keys != lb.keys {
                return Err(disagreement(faces, la, lb));
            }
            for key in &la.keys {
                joins.entry(*key).or_default().push((a, b));
            }
        }
    }
    for (index, [left, right]) in labels.iter().enumerate() {
        for key in &left.keys {
            if right.keys.binary_search(key).is_ok() {
                joins
                    .entry(*key)
                    .or_default()
                    .push((side(index, true), side(index, false)));
            }
        }
    }

    let mut regions: Vec<VisibleRegion> = Vec::new();
    let mut vertex_total = 0usize;
    for (key, half_edges) in &half_edges {
        options.check_cancelled()?;
        let (mesh, source) = *key;
        let id = faces.ids[mesh];
        let blame = |error: GeomError| Error::mesh(id, error);
        let (rings, traced) = trace(
            points,
            arrangement,
            half_edges,
            &|vertex| keep[vertex.0],
            &mut work.budget,
            options,
        )
        .map_err(blame)?;

        // The region of each boundary side, then of every side by component.
        let mut region_of_ring = vec![usize::MAX; rings.len()];
        for (region, found) in traced.iter().enumerate() {
            for ring in std::iter::once(&found.outer).chain(&found.holes) {
                region_of_ring[*ring] = region;
            }
        }
        let mut boundary_sides: Vec<(usize, usize)> = Vec::new();
        for (ring, traced_ring) in rings.iter().enumerate() {
            for half_edge in &traced_ring.half_edges {
                boundary_sides.push((
                    side(half_edge.piece.0, half_edge.forward),
                    region_of_ring[ring],
                ));
            }
        }
        let own_joins = joins.get(key).map_or(&[][..], Vec::as_slice);
        let mut nodes: Vec<usize> = own_joins.iter().flat_map(|&(a, b)| [a, b]).collect();
        nodes.extend(boundary_sides.iter().map(|&(s, _)| s));
        let mut components = Components::new(nodes);
        for &(a, b) in own_joins {
            components.join(a, b);
        }
        let mut region_of_component: BTreeMap<usize, usize> = BTreeMap::new();
        for &(s, region) in &boundary_sides {
            let root = components.index(s);
            let root = components.find(root);
            if region_of_component
                .insert(root, region)
                .is_some_and(|r| r != region)
            {
                return Err(blame(GeomError::Degenerate(
                    "a visible region's sides reach two rings".into(),
                )));
            }
        }

        // Ring points, for placing floating components.
        let outlines: Vec<(Vec<PointId>, Vec<Vec<PointId>>)> = traced
            .iter()
            .map(|found| {
                (
                    ring_points(arrangement, &rings[found.outer]),
                    found
                        .holes
                        .iter()
                        .map(|&hole| ring_points(arrangement, &rings[hole]))
                        .collect(),
                )
            })
            .collect();

        let mut triangles: Vec<Vec<u32>> = vec![Vec::new(); traced.len()];
        for index in 0..components.nodes.len() {
            let root = components.find(index);
            let node = components.nodes[index];
            let region = match region_of_component.get(&root) {
                Some(&region) => region,
                None => {
                    let piece = &arrangement.pieces[node / 2];
                    let at = arrangement.vertices[piece.ends[0].0].point;
                    let placed = place(points, &outlines, at, work).map_err(blame)?;
                    region_of_component.insert(root, placed);
                    placed
                }
            };
            for &t in label_of(node).fronts_of(key) {
                triangles[region].push((t - faces.firsts[mesh]) as u32);
            }
        }

        for (region, found) in traced.iter().enumerate() {
            if regions.len() >= limits.max_regions {
                return Err(Error::from(GeomError::BudgetExceeded {
                    resource: "visible-face regions",
                }));
            }
            let round = |ring: &TracedRing| round_ring(points, arrangement, ring).map_err(blame);
            let outer = round(&rings[found.outer])?;
            let mut holes = found
                .holes
                .iter()
                .map(|&hole| round(&rings[hole]))
                .collect::<Result<Vec<_>, _>>()?;
            holes.sort_by(|a, b| sequence_order(a, b));
            vertex_total += outer.len() + holes.iter().map(Vec::len).sum::<usize>();
            if vertex_total > limits.max_region_vertices {
                return Err(Error::from(GeomError::BudgetExceeded {
                    resource: "visible-face region vertices",
                }));
            }
            let mut listed = std::mem::take(&mut triangles[region]);
            listed.sort_unstable();
            listed.dedup();
            if source == FaceSource::Cap {
                listed.clear();
            }
            regions.push(VisibleRegion::new(id, source, outer, holes, listed));
        }
    }

    let position = |id| faces.ids.iter().position(|&m| m == id);
    regions.sort_by(|a, b| {
        position(a.mesh)
            .cmp(&position(b.mesh))
            .then(a.source.cmp(&b.source))
            .then_with(|| sequence_order(&a.outer, &b.outer))
    });
    Ok(regions)
}

/// The region of `key` that holds a floating component: the one whose outer
/// ring contains the point and none of whose holes does.
fn place(
    points: &ExactPoints,
    outlines: &[(Vec<PointId>, Vec<Vec<PointId>>)],
    at: PointId,
    work: &mut Work,
) -> GeomResult<usize> {
    for (region, (outer, holes)) in outlines.iter().enumerate() {
        work.budget.count()?;
        if points.locate(at, outer)? != Location::Inside {
            continue;
        }
        let mut inside_hole = false;
        for hole in holes {
            work.budget.count()?;
            if points.locate(at, hole)? != Location::Outside {
                inside_hole = true;
                break;
            }
        }
        if !inside_hole {
            return Ok(region);
        }
    }
    Err(GeomError::Degenerate(
        "a visible region's cluster lies in no region".into(),
    ))
}

/// A ring rounded once per vertex, in canonical rotation. A ring that
/// collapses under rounding is refused.
fn round_ring(
    points: &ExactPoints,
    arrangement: &Arrangement,
    ring: &TracedRing,
) -> GeomResult<Vec<Point2>> {
    let mut rounded = Vec::with_capacity(ring.vertices.len());
    for vertex in &ring.vertices {
        rounded.push(points.round(arrangement.vertices[vertex.0].point)?);
    }
    check_collapse(&rounded)?;
    Ok(canonical_rotation(rounded))
}

/// Refuse a ring with fewer than three points or two equal consecutive ones.
pub(crate) fn check_collapse(ring: &[Point2]) -> GeomResult<()> {
    let repeats = ring
        .iter()
        .zip(ring.iter().cycle().skip(1))
        .any(|(a, b)| a == b);
    if ring.len() < 3 || repeats {
        return Err(GeomError::Degenerate(
            "visible region collapsed under rounding".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_ring_that_collapses_under_rounding_is_refused() {
        let at = |x: f64, y: f64| Point2::new(x, y);
        let ok = [at(0.0, 0.0), at(1.0, 0.0), at(0.0, 1.0)];
        assert!(check_collapse(&ok).is_ok());
        // Two consecutive vertices one ulp apart round to the same point.
        let next = f64::from_bits(1.0f64.to_bits() + 1);
        let apart = [at(0.0, 0.0), at(1.0, 0.0), at(next, 0.0), at(0.0, 1.0)];
        assert!(check_collapse(&apart).is_ok(), "distinct after rounding");
        let rounded = [at(0.0, 0.0), at(1.0, 0.0), at(1.0, 0.0), at(0.0, 1.0)];
        assert!(matches!(
            check_collapse(&rounded),
            Err(GeomError::Degenerate(_))
        ));
        // The wrap-around pair counts, and fewer than three points do not
        // make a ring.
        let wrapped = [at(0.0, 0.0), at(1.0, 0.0), at(0.0, 0.0)];
        assert!(check_collapse(&wrapped).is_err());
        assert!(check_collapse(&ok[..2]).is_err());
    }
}
