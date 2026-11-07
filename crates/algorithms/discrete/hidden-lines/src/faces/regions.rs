//! Regions from the labelled arrangement.
//!
//! A piece is a boundary piece when its two labels differ. For each source
//! `K` its half-edges are the boundary pieces with `K` on exactly one side,
//! directed so that `K` is on the left. `trace` closes them into rings and
//! nests holes. A ring vertex is kept unless exactly two boundary pieces, of
//! any source, meet there in a straight line, so neighbouring regions share
//! bit-identical vertices.
//!
//! Those exact rings are then snap-rounded (ADR 0100). Every ring edge, the
//! straight run between two kept vertices, goes once into
//! `axiolid_exact_arrangement::snap_round`, which routes it through the
//! rounded position of every vertex whose rounding cell it passes through. A
//! source's boundary is netted on the snapped fragments (opposite traversals
//! cancel, so a part of zero width vanishes), traced again into rings by the
//! same rule and certified. A snapped region lists the triangles of every
//! exact region one of its pieces was cut from. Exact rings and regions with
//! no surviving boundary are counted as collapsed.
//!
//! The triangles of an exact region are those frontmost on its sides. The sides of
//! `K` are joined into components: the two sides of a wedge at a vertex (one
//! face of the arrangement) and the two sides of a piece that carries `K` on
//! both. A component without a boundary half-edge is a cluster that
//! dissolved inside a region, and is placed by exact point location.

use std::collections::BTreeMap;

use axiolid_contracts::{DisplacementBound, ExecutionOptions, GeomError, GeomResult, Sign};
use axiolid_core::Point2;
use axiolid_exact_arrangement::{
    canonical_rotation, check_winding, sequence_order, snap_round, trace, Arrangement, ExactPoints,
    HalfEdge, Location, PieceId, PointId, Segment, Snapped, TracedRegion, TracedRing, VertexId,
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

/// The regions and what snap rounding did to them.
pub(crate) struct Built {
    pub regions: Vec<VisibleRegion>,
    pub snap_displacement: DisplacementBound,
    pub rerouted_edges: usize,
    pub collapsed_regions: usize,
    pub collapsed_holes: usize,
}

/// The regions of every source: the exact regions, snap-rounded to the
/// binary64 grid (ADR 0100).
///
/// Every edge of an exact ring (the straight run between two kept vertices)
/// goes into one snap rounding. Each source's chain of edges is netted on the
/// snapped fragments, so a part that became zero-width cancels, and its
/// boundary is traced again into rings and certified.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build(
    faces: &Faces<'_>,
    points: &ExactPoints,
    arrangement: &Arrangement,
    labels: &[[Label; 2]],
    limits: OutputLimits,
    declared_bytes: usize,
    work: &mut Work,
    options: &ExecutionOptions,
) -> Result<Built, Error> {
    let boundary: Vec<bool> = labels.iter().map(|[l, r]| l.keys != r.keys).collect();
    let keep = kept_vertices(points, arrangement, &boundary)?;
    let sources = exact_sources(
        faces,
        points,
        arrangement,
        labels,
        &boundary,
        &keep,
        work,
        options,
    )?;

    // The edges of all exact rings, once each, and each source's chain over
    // them with the exact ring of every element.
    let vertex_point = |vertex: VertexId| arrangement.vertices[vertex.0].point;
    let mut edge_of: BTreeMap<(VertexId, VertexId), usize> = BTreeMap::new();
    let mut edges: Vec<Segment> = Vec::new();
    let mut chains: Vec<Vec<(usize, bool)>> = Vec::with_capacity(sources.len());
    let mut chain_rings: Vec<Vec<usize>> = Vec::with_capacity(sources.len());
    for source in &sources {
        let (mut chain, mut rings) = (Vec::new(), Vec::new());
        for (index, ring) in source.rings.iter().enumerate() {
            let count = ring.vertices.len();
            for at in 0..count {
                let (from, to) = (ring.vertices[at], ring.vertices[(at + 1) % count]);
                let ends = (from.min(to), from.max(to));
                let edge = *edge_of.entry(ends).or_insert_with(|| {
                    edges.push(Segment {
                        ends: [vertex_point(ends.0), vertex_point(ends.1)],
                    });
                    edges.len() - 1
                });
                chain.push((edge, from < to));
                rings.push(index);
            }
        }
        chains.push(chain);
        chain_rings.push(rings);
    }
    let snapped = snap_round(points, &edges, &mut work.budget, declared_bytes, options)?;

    // The net boundary of each source, and the vertices to keep.
    let mut boundary_after = vec![false; snapped.arrangement.pieces.len()];
    let mut nets = Vec::with_capacity(sources.len());
    for (source, chain) in sources.iter().zip(&chains) {
        let id = faces.ids[source.key.0];
        let net = snapped
            .class_boundary(chain)
            .map_err(|error| Error::mesh(id, error))?;
        for half_edge in &net.half_edges {
            boundary_after[half_edge.piece.0] = true;
        }
        nets.push(net);
    }
    let keep_after = kept_vertices(&snapped.points, &snapped.arrangement, &boundary_after)?;

    let mut built = Built {
        regions: Vec::new(),
        snap_displacement: snapped.displacement,
        rerouted_edges: snapped.rerouted,
        collapsed_regions: 0,
        collapsed_holes: 0,
    };
    let mut vertex_total = 0usize;
    for (((source, net), element_rings), chain) in
        sources.iter().zip(&nets).zip(&chain_rings).zip(&chains)
    {
        options.check_cancelled()?;
        let (mesh, kind) = source.key;
        let id = faces.ids[mesh];
        let blame = |error: GeomError| Error::mesh(id, error);

        // Which exact rings and regions survive, and which exact regions
        // each snapped piece was cut from. Each ring is netted alone: a ring
        // that collapses onto another surviving edge of the class leaves
        // nothing of its own, and its triangles are not attributed to it.
        let ring_pieces =
            ring_boundaries(&snapped, chain, element_rings, source.rings.len()).map_err(blame)?;
        let ring_survives: Vec<bool> = ring_pieces
            .iter()
            .map(|pieces| !pieces.is_empty())
            .collect();
        let mut region_of_ring = vec![usize::MAX; source.rings.len()];
        for (region, found) in source.traced.iter().enumerate() {
            for ring in std::iter::once(&found.outer).chain(&found.holes) {
                region_of_ring[*ring] = region;
            }
            if ring_survives[found.outer] {
                built.collapsed_holes += found
                    .holes
                    .iter()
                    .filter(|&&hole| !ring_survives[hole])
                    .count();
            } else {
                built.collapsed_regions += 1;
            }
        }
        let mut regions_of_piece: BTreeMap<PieceId, Vec<usize>> = BTreeMap::new();
        for (ring, pieces) in ring_pieces.iter().enumerate() {
            for half_edge in pieces {
                regions_of_piece
                    .entry(half_edge.piece)
                    .or_default()
                    .push(region_of_ring[ring]);
            }
        }

        let (rings, traced) = trace(
            &snapped.points,
            &snapped.arrangement,
            &net.half_edges,
            &|vertex| keep_after[vertex.0],
            &mut work.budget,
            options,
        )
        .map_err(blame)?;
        check_winding(
            &snapped.points,
            &snapped.arrangement,
            &rings,
            &traced,
            &mut work.budget,
            options,
        )
        .map_err(blame)?;

        for found in &traced {
            if built.regions.len() >= limits.max_regions {
                return Err(Error::from(GeomError::BudgetExceeded {
                    resource: "visible-face regions",
                }));
            }
            let positions = |ring: &TracedRing| -> Result<Vec<Point2>, Error> {
                let at = ring
                    .vertices
                    .iter()
                    .map(|&vertex| snapped.point(vertex))
                    .collect::<GeomResult<Vec<_>>>()
                    .map_err(blame)?;
                Ok(canonical_rotation(at))
            };
            let outer = positions(&rings[found.outer])?;
            let holes = found
                .holes
                .iter()
                .map(|&hole| positions(&rings[hole]))
                .collect::<Result<Vec<_>, _>>()?;
            let holes = ordered_holes(holes);
            vertex_total += outer.len() + holes.iter().map(Vec::len).sum::<usize>();
            if vertex_total > limits.max_region_vertices {
                return Err(Error::from(GeomError::BudgetExceeded {
                    resource: "visible-face region vertices",
                }));
            }
            // The triangles of every exact region a piece of it was cut from.
            let mut listed: Vec<u32> = Vec::new();
            for &ring in std::iter::once(&found.outer).chain(&found.holes) {
                for half_edge in &rings[ring].half_edges {
                    for &region in regions_of_piece
                        .get(&half_edge.piece)
                        .map_or(&[][..], Vec::as_slice)
                    {
                        listed.extend(&source.triangles[region]);
                    }
                }
            }
            listed.sort_unstable();
            listed.dedup();
            if kind == FaceSource::Cap {
                listed.clear();
            }
            built
                .regions
                .push(VisibleRegion::new(id, kind, outer, holes, listed));
        }
    }

    let position = |id| faces.ids.iter().position(|&m| m == id);
    built.regions.sort_by(|a, b| {
        position(a.mesh)
            .cmp(&position(b.mesh))
            .then(a.source.cmp(&b.source))
            .then_with(|| sequence_order(&a.outer, &b.outer))
    });
    Ok(built)
}

/// A source's exact regions: its rings, their nesting and, per region, the
/// triangles frontmost on its sides.
struct ExactSource {
    key: Key,
    rings: Vec<TracedRing>,
    traced: Vec<TracedRegion>,
    triangles: Vec<Vec<u32>>,
}

/// The exact regions of every source, with the vertices `keep` asks for.
#[allow(clippy::too_many_arguments)]
fn exact_sources(
    faces: &Faces<'_>,
    points: &ExactPoints,
    arrangement: &Arrangement,
    labels: &[[Label; 2]],
    boundary: &[bool],
    keep: &[bool],
    work: &mut Work,
    options: &ExecutionOptions,
) -> Result<Vec<ExactSource>, Error> {
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

    let mut sources = Vec::with_capacity(half_edges.len());
    for (key, half_edges) in &half_edges {
        options.check_cancelled()?;
        let (mesh, _) = *key;
        let blame = |error: GeomError| Error::mesh(faces.ids[mesh], error);
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

        for listed in &mut triangles {
            listed.sort_unstable();
            listed.dedup();
        }
        sources.push(ExactSource {
            key: *key,
            rings,
            traced,
            triangles,
        });
    }
    Ok(sources)
}

/// The net boundary of each exact ring on its own: `chain` holds the elements
/// of every ring of a source and `element_rings` the ring of each element.
fn ring_boundaries(
    snapped: &Snapped,
    chain: &[(usize, bool)],
    element_rings: &[usize],
    rings: usize,
) -> GeomResult<Vec<Vec<HalfEdge>>> {
    (0..rings)
        .map(|ring| {
            let own: Vec<(usize, bool)> = chain
                .iter()
                .zip(element_rings)
                .filter(|&(_, &of)| of == ring)
                .map(|(&element, _)| element)
                .collect();
            Ok(snapped.class_boundary(&own)?.half_edges)
        })
        .collect()
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

/// The holes of a region sorted lexicographically by their points.
fn ordered_holes(mut holes: Vec<Vec<Point2>>) -> Vec<Vec<Point2>> {
    holes.sort_by(|a, b| sequence_order(a, b));
    holes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn holes_are_sorted_by_their_points() {
        let ring = |x: f64| vec![Point2::new(x, 0.0), Point2::new(x + 1.0, 0.0)];
        let sorted = ordered_holes(vec![ring(5.0), ring(1.0), ring(3.0)]);
        assert_eq!(sorted, vec![ring(1.0), ring(3.0), ring(5.0)]);
    }

    #[test]
    fn a_ring_that_collapses_onto_a_surviving_edge_of_its_class_leaves_nothing() {
        use axiolid_contracts::ExecutionOptions;
        use axiolid_core::Tolerance;
        use axiolid_exact::{Arith, Dyadic};
        use axiolid_exact_arrangement::PairBudget;

        let d = |v: f64| Dyadic::try_from_f64(v).expect("finite");
        // 2^-1080, below the smallest binary64 number: it rounds to zero.
        let tiny = d(2f64.powi(-540)).mul(&d(2f64.powi(-540)));
        let mut points = ExactPoints::new();
        let mut at = |x: f64, y: f64| points.push_flat(Point2::new(x, y)).unwrap();
        let square = [at(0.0, 0.0), at(4.0, 0.0), at(4.0, 4.0), at(0.0, 4.0)];
        // A counter-clockwise sliver below the square's bottom edge.
        let below = |x: f64, depth: f64, points: &mut ExactPoints| {
            points
                .push([d(x), d(-depth).mul(&tiny), d(1.0)])
                .expect("point")
        };
        let sliver = [
            below(1.0, 3.0, &mut points),
            below(3.0, 3.0, &mut points),
            below(2.0, 1.0, &mut points),
        ];
        let mut edges = Vec::new();
        let mut chain = Vec::new();
        let mut element_rings = Vec::new();
        for (ring, corners) in [&square[..], &sliver[..]].into_iter().enumerate() {
            for i in 0..corners.len() {
                edges.push(Segment {
                    ends: [corners[i], corners[(i + 1) % corners.len()]],
                });
                chain.push((edges.len() - 1, true));
                element_rings.push(ring);
            }
        }
        let options = ExecutionOptions::new(Tolerance::METRE);
        let snapped = snap_round(
            &points,
            &edges,
            &mut PairBudget::new(100_000, "test pairs"),
            0,
            &options,
        )
        .expect("snapped");
        // Netted together, the sliver's elements still carry the square's
        // bottom edge, so they would be reported as surviving.
        let together = snapped.class_boundary(&chain).expect("class");
        assert!(together.survives[4..].iter().any(|&kept| kept));
        let alone = ring_boundaries(&snapped, &chain, &element_rings, 2).expect("rings");
        assert_eq!(
            alone[0].len(),
            7,
            "the square keeps its boundary, split at the sliver"
        );
        assert!(alone[1].is_empty(), "the sliver has none of its own");
    }
}
