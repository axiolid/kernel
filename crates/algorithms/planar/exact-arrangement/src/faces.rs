//! The faces of the complement of an arrangement.
//!
//! The pieces of an [`Arrangement`] cut the plane into faces. Each piece has
//! a face on either side, and [`faces`] numbers them:
//!
//! 1. Contours: a half-edge is a piece with a direction and its face on the
//!    left. At a vertex the half-edge that follows an incoming one is the
//!    first met turning clockwise from the way it came (the rule of
//!    [`trace`](crate::trace)), so the half-edges fall into orbits, one per
//!    boundary walk. A dangling piece is walked on both its sides.
//! 2. Components: the pieces join into connected components. The extreme
//!    vertex of a component (smallest by `x`, then `y`) is a convex corner,
//!    so exactly one wedge between consecutive pieces there is wider than a
//!    half turn. The contour that leaves along the piece that starts that
//!    wedge has the unbounded side of the component on its left: it is the
//!    component's outer contour.
//! 3. Faces: every other contour bounds one bounded face of its component,
//!    numbered from 1 by the smallest half-edge of the contour.
//! 4. Nesting: a component that lies inside a bounded face of another one
//!    has that face on the sides of its outer contour. Of the faces that
//!    hold its extreme vertex the innermost wins, which is the one whose
//!    component has the greatest extreme vertex, because a nested component
//!    lies strictly inside the walk around the face.
//!
//! A piece with one face on both sides is dangling, or part of a tree.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Sign};

use crate::arrangement::{Arrangement, PairBudget, PieceId, VertexId};
use crate::points::{inconsistent, ExactPoints, Location, PointId, RingBox};

/// A face of the complement of an [`Arrangement`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FaceId(pub usize);

/// The face on each side of every piece.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Faces {
    /// Per piece, the faces left and right of `ends[0]` to `ends[1]`.
    sides: Vec<[FaceId; 2]>,
    count: usize,
}

impl Faces {
    /// The unbounded face.
    pub const UNBOUNDED: FaceId = FaceId(0);

    /// The number of faces, the unbounded one included.
    pub fn count(&self) -> usize {
        self.count
    }

    /// The face left of `piece`, running from `ends[0]` to `ends[1]`.
    pub fn left(&self, piece: PieceId) -> FaceId {
        self.sides[piece.0][0]
    }

    /// The face right of `piece`, running from `ends[0]` to `ends[1]`.
    pub fn right(&self, piece: PieceId) -> FaceId {
        self.sides[piece.0][1]
    }
}

/// The half-edge of `piece` that runs from `ends[0]` when `forward`:
/// `2 * piece`, plus one for the other direction.
fn half_edge(piece: PieceId, forward: bool) -> usize {
    2 * piece.0 + usize::from(!forward)
}

/// Find the root of `item` in a union-find forest.
fn find(parent: &mut [usize], mut item: usize) -> usize {
    while parent[item] != item {
        parent[item] = parent[parent[item]];
        item = parent[item];
    }
    item
}

/// Number the faces of the complement of `arrangement`.
///
/// Face 0 is the unbounded face. Bounded faces are numbered from 1,
/// ascending by the smallest half-edge (`2 * piece`, plus one for the
/// backward direction) of their boundary walk. With `V` vertices, `E`
/// pieces and `C` connected components, the count `F` satisfies
/// `V - E + F = 1 + C`.
///
/// Each containment decided between components counts one pair against
/// `budget`; cancellation is checked per vertex. A component whose extreme
/// vertex does not decide its outer contour is `Degenerate`.
pub fn faces(
    points: &ExactPoints,
    arrangement: &Arrangement,
    budget: &mut PairBudget,
    options: &ExecutionOptions,
) -> GeomResult<Faces> {
    let piece_count = arrangement.pieces.len();
    let mut next = vec![usize::MAX; 2 * piece_count];
    let mut fans = Vec::with_capacity(arrangement.vertices.len());
    for vertex in 0..arrangement.vertices.len() {
        options.check_cancelled()?;
        let fan = arrangement.fan(points, VertexId(vertex))?;
        for (index, &(piece, forward)) in fan.iter().enumerate() {
            // The piece arrives here along the twin of the one that leaves,
            // and the boundary goes on along the spoke clockwise of it.
            let (turn, turn_forward) = fan[(index + fan.len() - 1) % fan.len()];
            next[half_edge(piece, forward) ^ 1] = half_edge(turn, turn_forward);
        }
        fans.push(fan);
    }
    let tail = |half: usize| {
        let ends = arrangement.pieces[half / 2].ends;
        ends[usize::from(half % 2 == 1)].0
    };

    // Contours, numbered by their smallest half-edge.
    let mut contour_of = vec![usize::MAX; 2 * piece_count];
    let mut contours = Vec::<Vec<usize>>::new();
    for start in 0..2 * piece_count {
        if contour_of[start] != usize::MAX {
            continue;
        }
        let mut walk = Vec::new();
        let mut half = start;
        loop {
            contour_of[half] = contours.len();
            walk.push(half);
            half = next[half];
            if half == start {
                break;
            }
            if half == usize::MAX || contour_of[half] != usize::MAX {
                return Err(inconsistent());
            }
        }
        contours.push(walk);
    }

    // Components, each with its extreme vertex.
    let mut parent: Vec<usize> = (0..arrangement.vertices.len()).collect();
    for piece in &arrangement.pieces {
        let (low, high) = (
            find(&mut parent, piece.ends[0].0),
            find(&mut parent, piece.ends[1].0),
        );
        parent[low.max(high)] = low.min(high);
    }
    let mut extreme = vec![usize::MAX; arrangement.vertices.len()];
    for (vertex, fan) in fans.iter().enumerate() {
        if fan.is_empty() {
            continue;
        }
        let root = find(&mut parent, vertex);
        let best = extreme[root];
        if best == usize::MAX
            || points
                .lexicographic(
                    arrangement.vertices[vertex].point,
                    arrangement.vertices[best].point,
                )?
                .is_lt()
        {
            extreme[root] = vertex;
        }
    }
    let roots: Vec<usize> = (0..extreme.len())
        .filter(|&root| extreme[root] != usize::MAX)
        .collect();

    // The outer contour of each component.
    let mut is_outer = vec![false; contours.len()];
    for &root in &roots {
        let vertex = extreme[root];
        let fan = &fans[vertex];
        let here = arrangement.vertices[vertex].point;
        let far = |&(piece, forward): &(PieceId, bool)| {
            arrangement.vertices[arrangement.pieces[piece.0].ends[usize::from(forward)].0].point
        };
        let mut wide = None;
        for (index, spoke) in fan.iter().enumerate() {
            let after = &fan[(index + 1) % fan.len()];
            let turns_clockwise =
                fan.len() == 1 || points.orient(here, far(spoke), far(after))? == Sign::Negative;
            if turns_clockwise {
                if wide.is_some() {
                    return Err(GeomError::Degenerate(
                        "a component has two wedges wider than a half turn at its extreme vertex"
                            .into(),
                    ));
                }
                wide = Some(index);
            }
        }
        let index = wide.ok_or_else(|| {
            GeomError::Degenerate("a component has no wedge wider than a half turn".into())
        })?;
        let (piece, forward) = fan[index];
        is_outer[contour_of[half_edge(piece, forward)]] = true;
    }

    // Bounded faces, then the faces around the nested components.
    let mut face_of = vec![0usize; contours.len()];
    let mut count = 1;
    for (contour, outer) in is_outer.iter().enumerate() {
        if !outer {
            face_of[contour] = count;
            count += 1;
        }
    }
    if roots.len() > 1 {
        let walks: Vec<Vec<PointId>> = contours
            .iter()
            .map(|walk| {
                walk.iter()
                    .map(|&half| arrangement.vertices[tail(half)].point)
                    .collect()
            })
            .collect();
        let boxes: Vec<RingBox> = walks.iter().map(|walk| points.ring_box(walk)).collect();
        let component_of = |contour: usize| find_root(&parent, tail(contours[contour][0]));
        let outer_contours: Vec<usize> = (0..contours.len()).filter(|&c| is_outer[c]).collect();
        for &outer in &outer_contours {
            options.check_cancelled()?;
            let root = component_of(outer);
            let vertex = extreme[root];
            let here = arrangement.vertices[vertex].point;
            let here_box = points.ring_box(&[here]);
            let mut innermost: Option<(usize, usize)> = None;
            for (contour, walk) in walks.iter().enumerate() {
                if is_outer[contour]
                    || component_of(contour) == root
                    || !here_box.possibly_within(&boxes[contour])
                {
                    continue;
                }
                options.check_cancelled()?;
                budget.count()?;
                match points.locate(here, walk)? {
                    Location::Inside => {}
                    Location::Outside => continue,
                    Location::Boundary => return Err(inconsistent()),
                }
                let holder = extreme[component_of(contour)];
                let deeper = match innermost {
                    None => true,
                    Some((_, known)) => points
                        .lexicographic(
                            arrangement.vertices[holder].point,
                            arrangement.vertices[known].point,
                        )?
                        .is_gt(),
                };
                if deeper {
                    innermost = Some((contour, holder));
                }
            }
            if let Some((contour, _)) = innermost {
                face_of[outer] = face_of[contour];
            }
        }
    }

    let sides = (0..piece_count)
        .map(|piece| {
            [
                FaceId(face_of[contour_of[2 * piece]]),
                FaceId(face_of[contour_of[2 * piece + 1]]),
            ]
        })
        .collect();
    Ok(Faces { sides, count })
}

/// The root of `item` without compressing the paths.
fn find_root(parent: &[usize], mut item: usize) -> usize {
    while parent[item] != item {
        item = parent[item];
    }
    item
}
