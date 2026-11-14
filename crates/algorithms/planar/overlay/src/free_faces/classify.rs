//! Which faces are free, and which free faces are closed gaps.
//!
//! A face of the arrangement of inputs and bridges is read as obstacle or
//! free by ring membership, as ADR 0104 reads a region: crossing a piece
//! flips, once per carrier, the ring an obstacle edge belongs to; a barrier
//! or a bridge flips nothing. The walk starts in the unbounded face, which
//! is in no ring, and goes over shared pieces depth first, undoing each flip
//! on its way back. A face is in a polygon when it is in the outer ring and
//! in none of the holes; it is an obstacle when it is in any polygon.
//!
//! A bounded free face is a closed gap when it fits in a strip narrower than
//! the closing width: for some edge of the convex hull of its exact ring,
//! every hull vertex is strictly nearer to that edge's line than the width.
//! The scan for an edge starts at the vertex opposite it and stops at the
//! first that is not nearer, so a room is refused after one test per edge.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};
use axiolid_exact_arrangement::{
    trace, Arrangement, ExactPoints, Faces, HalfEdge, PairBudget, PieceId, PointId, TracedRing,
    VertexId,
};

use super::bridges::piece_sources;
use super::{Model, PieceSource};

/// A bounded free face, traced exactly.
pub(super) struct Class {
    /// Whether it is a closed gap rather than a face.
    pub gap: bool,
    /// Its rings: one outer ring and its holes.
    pub rings: Vec<TracedRing>,
    /// The index in `rings` of the outer ring.
    pub outer: usize,
}

pub(super) struct Classified {
    /// The bounded free faces, in the order of their numbers.
    pub classes: Vec<Class>,
    /// Per bridge found, whether it lies in free space.
    pub bridge_alive: Vec<bool>,
    /// Per piece, the sources of its carriers.
    pub piece_sources: Vec<Vec<PieceSource>>,
}

/// Whether the collinear `at` lies strictly between its two neighbours.
pub(super) fn straight(
    points: &ExactPoints,
    before: PointId,
    at: PointId,
    after: PointId,
) -> GeomResult<bool> {
    Ok(points.orient(before, at, after)? == Sign::Zero
        && points.strictly_between(before, at, after)?)
}

/// The vertices a ring must keep: unless exactly two boundary pieces meet
/// there in a straight line and carry equal `labels`.
pub(super) fn kept_vertices<L: PartialEq>(
    points: &ExactPoints,
    arrangement: &Arrangement,
    boundary: &[bool],
    labels: &[L],
) -> GeomResult<Vec<bool>> {
    let mut around: Vec<Vec<(VertexId, usize)>> = vec![Vec::new(); arrangement.vertices.len()];
    for (index, piece) in arrangement.pieces.iter().enumerate() {
        if boundary[index] {
            around[piece.ends[0].0].push((piece.ends[1], index));
            around[piece.ends[1].0].push((piece.ends[0], index));
        }
    }
    let point = |vertex: VertexId| arrangement.vertices[vertex.0].point;
    let mut keep = vec![true; around.len()];
    for (vertex, others) in around.iter().enumerate() {
        if let [(a, first), (b, second)] = others[..] {
            keep[vertex] = labels[first] != labels[second]
                || !straight(
                    points,
                    point(a),
                    arrangement.vertices[vertex].point,
                    point(b),
                )?;
        }
    }
    Ok(keep)
}

/// The rings of the polygons, as flags per ring, and the polygon each is of.
struct Membership {
    /// Per segment, its ring when it is an obstacle edge.
    ring_of: Vec<Option<usize>>,
    polygon_of: Vec<usize>,
    hole: Vec<bool>,
    outer_of: Vec<usize>,
    inside: Vec<bool>,
    holes_in: Vec<u32>,
    covering: usize,
}

impl Membership {
    fn new(sources: &[PieceSource]) -> Self {
        let mut index = BTreeMap::<(usize, usize), usize>::new();
        let mut member = Self {
            ring_of: Vec::with_capacity(sources.len()),
            polygon_of: Vec::new(),
            hole: Vec::new(),
            outer_of: Vec::new(),
            inside: Vec::new(),
            holes_in: Vec::new(),
            covering: 0,
        };
        for source in sources {
            let PieceSource::Obstacle { polygon, ring, .. } = *source else {
                member.ring_of.push(None);
                continue;
            };
            let next = member.polygon_of.len();
            let at = *index.entry((polygon, ring)).or_insert(next);
            if at == next {
                member.polygon_of.push(polygon);
                member.hole.push(ring > 0);
                member.inside.push(false);
                if member.outer_of.len() <= polygon {
                    member.outer_of.resize(polygon + 1, usize::MAX);
                    member.holes_in.resize(polygon + 1, 0);
                }
                if ring == 0 {
                    member.outer_of[polygon] = at;
                }
            }
            member.ring_of.push(Some(at));
        }
        member
    }

    fn covered(&self, polygon: usize) -> bool {
        self.outer_of[polygon] != usize::MAX
            && self.inside[self.outer_of[polygon]]
            && self.holes_in[polygon] == 0
    }

    /// Cross the rings carried by the segments of a piece.
    fn cross(&mut self, segments: &[usize]) {
        for &segment in segments {
            let Some(ring) = self.ring_of[segment] else {
                continue;
            };
            let polygon = self.polygon_of[ring];
            let before = self.covered(polygon);
            self.inside[ring] = !self.inside[ring];
            if self.hole[ring] {
                if self.inside[ring] {
                    self.holes_in[polygon] += 1;
                } else {
                    self.holes_in[polygon] -= 1;
                }
            }
            match (before, self.covered(polygon)) {
                (false, true) => self.covering += 1,
                (true, false) => self.covering -= 1,
                _ => {}
            }
        }
    }
}

fn inconsistent(message: &str) -> GeomError {
    GeomError::Degenerate(message.into())
}

/// Per face, whether an obstacle covers it.
fn obstacle_faces(
    model: &Model,
    arrangement: &Arrangement,
    faces: &Faces,
    options: &ExecutionOptions,
) -> GeomResult<Vec<bool>> {
    let mut member = Membership::new(&model.sources);
    let mut adjacent: Vec<Vec<usize>> = vec![Vec::new(); faces.count()];
    for piece in 0..arrangement.pieces.len() {
        let (left, right) = (faces.left(PieceId(piece)), faces.right(PieceId(piece)));
        if left != right {
            adjacent[left.0].push(piece);
            adjacent[right.0].push(piece);
        }
    }
    let mut obstacle = vec![false; faces.count()];
    let mut visited = vec![false; faces.count()];
    visited[Faces::UNBOUNDED.0] = true;
    // Per face on the path: its next neighbour to try, and the piece it was
    // entered over.
    let mut stack: Vec<(usize, usize, Option<usize>)> = vec![(Faces::UNBOUNDED.0, 0, None)];
    while let Some(&(face, at, via)) = stack.last() {
        options.check_cancelled()?;
        let Some(&piece) = adjacent[face].get(at) else {
            stack.pop();
            if let Some(piece) = via {
                member.cross(&arrangement.pieces[piece].segments);
            }
            continue;
        };
        if let Some(top) = stack.last_mut() {
            top.1 += 1;
        }
        let id = PieceId(piece);
        let other = if faces.left(id).0 == face {
            faces.right(id).0
        } else {
            faces.left(id).0
        };
        if visited[other] {
            continue;
        }
        visited[other] = true;
        member.cross(&arrangement.pieces[piece].segments);
        obstacle[other] = member.covering > 0;
        stack.push((other, 0, Some(piece)));
    }
    if visited.contains(&false) {
        return Err(inconsistent(
            "a face cannot be reached from the unbounded one",
        ));
    }
    Ok(obstacle)
}

/// The convex hull of exact points, counter-clockwise, without repeats.
fn hull(points: &ExactPoints, ring: &[PointId]) -> GeomResult<Vec<PointId>> {
    let mut sorted = ring.to_vec();
    let mut failure = None;
    sorted.sort_by(|&a, &b| {
        points.lexicographic(a, b).unwrap_or_else(|error| {
            failure.get_or_insert(error);
            Ordering::Equal
        })
    });
    if let Some(error) = failure {
        return Err(error);
    }
    sorted.dedup();
    if sorted.len() < 3 {
        return Ok(sorted);
    }
    let mut chain: Vec<PointId> = Vec::new();
    let run = |chain: &mut Vec<PointId>, point: PointId, floor: usize| -> GeomResult<()> {
        while chain.len() > floor
            && points.orient(chain[chain.len() - 2], chain[chain.len() - 1], point)?
                != Sign::Positive
        {
            chain.pop();
        }
        chain.push(point);
        Ok(())
    };
    for &point in &sorted {
        run(&mut chain, point, 1)?;
    }
    let floor = chain.len();
    for &point in sorted.iter().rev().skip(1) {
        run(&mut chain, point, floor)?;
    }
    chain.pop();
    Ok(chain)
}

/// Whether the face with exact outer ring `ring` fits in a strip narrower
/// than `gap`.
fn fits_in_strip(points: &ExactPoints, ring: &[PointId], gap: f64) -> GeomResult<bool> {
    let hull = hull(points, ring)?;
    let count = hull.len();
    for edge in 0..count {
        let (start, end) = (hull[edge], hull[(edge + 1) % count]);
        let mut all = true;
        for step in 0..count {
            let vertex = hull[(edge + count / 2 + step) % count];
            if points.compare_line_distance(vertex, start, end, gap)? != Ordering::Less {
                all = false;
                break;
            }
        }
        if all {
            return Ok(true);
        }
    }
    Ok(false)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn classify(
    model: &Model,
    arrangement: &Arrangement,
    faces: &Faces,
    inputs: usize,
    bridges: usize,
    gap: f64,
    budget: &mut PairBudget,
    options: &ExecutionOptions,
) -> GeomResult<Classified> {
    let obstacle = obstacle_faces(model, arrangement, faces, options)?;
    let mut bridge_alive = Vec::with_capacity(bridges);
    for bridge in 0..bridges {
        let (piece, _) = *arrangement.paths[inputs + bridge]
            .first()
            .ok_or_else(|| inconsistent("a bridge has no piece"))?;
        bridge_alive.push(!obstacle[faces.left(piece).0]);
    }

    // A bounded free face is a class.
    let mut class_of: Vec<Option<usize>> = vec![None; faces.count()];
    let mut count = 0;
    for face in 1..faces.count() {
        if !obstacle[face] {
            class_of[face] = Some(count);
            count += 1;
        }
    }
    let mut boundary = vec![false; arrangement.pieces.len()];
    let mut half_edges: Vec<Vec<HalfEdge>> = vec![Vec::new(); count];
    for (index, flag) in boundary.iter_mut().enumerate() {
        let piece = PieceId(index);
        let (left, right) = (
            class_of[faces.left(piece).0],
            class_of[faces.right(piece).0],
        );
        if left == right {
            continue;
        }
        *flag = true;
        if let Some(class) = left {
            half_edges[class].push(HalfEdge {
                piece,
                forward: true,
            });
        }
        if let Some(class) = right {
            half_edges[class].push(HalfEdge {
                piece,
                forward: false,
            });
        }
    }
    let sources = piece_sources(arrangement, &model.sources);
    let keep = kept_vertices(&model.points, arrangement, &boundary, &sources)?;

    let mut classes = Vec::with_capacity(count);
    for half_edges in &half_edges {
        options.check_cancelled()?;
        let (rings, regions) = trace(
            &model.points,
            arrangement,
            half_edges,
            &|vertex| keep[vertex.0],
            budget,
            options,
        )?;
        let [region] = &regions[..] else {
            return Err(inconsistent(
                "a free face has not exactly one outer boundary",
            ));
        };
        let outline: Vec<PointId> = rings[region.outer]
            .vertices
            .iter()
            .map(|vertex| arrangement.vertices[vertex.0].point)
            .collect();
        let gap = gap > 0.0 && fits_in_strip(&model.points, &outline, gap)?;
        classes.push(Class {
            gap,
            outer: region.outer,
            rings,
        });
    }
    Ok(Classified {
        classes,
        bridge_alive,
        piece_sources: sources,
    })
}
