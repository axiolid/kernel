//! Rings that touch at single points, under [`PinchPolicy::Accept`] (#262).
//!
//! A plan footprint that is the union of two shadows touching at a corner,
//! a hole whose corner sits on the outer ring or on another hole, an outer
//! ring pinched at a vertex: each bounds a valid region, but not a polygon
//! with holes in the sense of [`validate`](super::validate), and no single
//! ring through the pinch can be ear clipped (its two visits of the pinch
//! point bound a neck of zero width). The region is rebuilt instead, every
//! decision exact:
//!
//! 1. Coincident vertices become one ([`canonical`]): every triangle uses
//!    the first index of its point in `outer ++ holes`. A vertex inside
//!    another ring's edge (a [`Touch`]) is inserted into that edge, so every
//!    contact is a shared vertex.
//! 2. Each ring is split at the vertices it visits twice into simple
//!    loops. A loop's depth is the number of loops containing it; even
//!    depths bound the region from outside and are turned
//!    counter-clockwise, odd ones are holes, turned clockwise. A hole ring
//!    must still lie inside the outer ring and outside every other hole.
//! 3. At a shared vertex the loops' edges, sorted by exact angle, must
//!    alternate leaving and arriving: then the region near it is a set of
//!    wedges, each from a leaving edge counter-clockwise to the next
//!    arriving one, and anything else means the rings cross there. Each
//!    wedge becomes one node, joining that arriving edge to that leaving
//!    edge. The nodes then form the boundary cycles of the region's
//!    connected parts: two shadows touching at a corner become two cycles,
//!    a hole touching the outer ring merges with it into one.
//! 4. Each cycle that bounds from outside is bridged to the hole cycles it
//!    contains and clipped as before ([`bridge`](super::bridge),
//!    [`clip`](super::clip)).
//!
//! The certificate is the same, over the split loops: positive triangles
//! whose boundary is exactly the loops' edges, every other edge twinned,
//! and `n + 2h - 2c` triangles for `n` loop vertices, `h` hole cycles and
//! `c` outer cycles (each outer cycle with `k` vertices has interior
//! angles summing to `(k - 2) pi`, each hole cycle `(k + 2) pi`, and a
//! triangle takes `pi`). A pinch vertex is counted once per visit, which
//! is the "pinch vertex counted twice" of the issue.

use std::collections::{HashMap, HashSet};

use core::cmp::Ordering;
use core::ops::Range;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::Point2;

use super::validate::{name, refuse, ring_turns_left, strictly_inside, Touch};
use super::{bridge, by_value, certify, clip, lexicographic, locally_inside, orient, Loop};

/// The first index of each point's coordinates, if any two vertices
/// coincide or `touches` is not empty; `None` when the rings touch nowhere
/// and the ordinary path applies.
pub(super) fn canonical(points: &[Point2], touches: &[Touch]) -> Option<Vec<u32>> {
    let mut first: HashMap<(u64, u64), u32> = HashMap::with_capacity(points.len());
    let mut shared = false;
    let canon: Vec<u32> = points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            // `+ 0.0` turns -0.0 into 0.0, which compare equal.
            let key = ((p.x + 0.0).to_bits(), (p.y + 0.0).to_bits());
            let index = *first.entry(key).or_insert(i as u32);
            shared |= index != i as u32;
            index
        })
        .collect();
    (shared || !touches.is_empty()).then_some(canon)
}

/// One simple loop of a ring split at its repeated vertices.
#[derive(Debug, Clone)]
struct Piece {
    /// The ring it came from: 0 the outer ring, `h + 1` hole `h`.
    ring: usize,
    vertices: Vec<u32>,
}

/// A loop vertex visit: a node of one boundary cycle.
#[derive(Debug, Clone, Copy)]
struct Node {
    vertex: u32,
    prev: usize,
    next: usize,
}

/// Triangulate rings that touch at single points, given the canonical
/// indices from [`canonical`] and the touches from validation. `source`
/// names each ring: the input ring it is part of (#270).
pub(super) fn triangulate(
    points: &[Point2],
    rings: &[Range<usize>],
    source: &[usize],
    canon: &[u32],
    touches: &[Touch],
) -> GeomResult<Vec<[u32; 3]>> {
    let at = |v: u32| points[v as usize];

    // 1. Canonical rings with every touching vertex inserted in its edge.
    let mut inserted: HashMap<(usize, usize), Vec<u32>> = HashMap::new();
    for touch in touches {
        let list = inserted.entry((touch.ring, touch.index)).or_default();
        let vertex = canon[touch.vertex as usize];
        if !list.contains(&vertex) {
            list.push(vertex);
        }
    }
    let mut pieces: Vec<Piece> = Vec::new();
    for (r, ring) in rings.iter().enumerate() {
        let len = ring.len();
        let mut walk: Vec<u32> = Vec::with_capacity(len);
        for k in 0..len {
            let (a, b) = (ring.start + k, ring.start + (k + 1) % len);
            walk.push(canon[a]);
            if let Some(list) = inserted.get_mut(&(r, k)) {
                let (pa, pb) = (points[a], points[b]);
                // Collinear and strictly inside `a b`: order along the
                // axis on which the edge moves.
                list.sort_by(|&u, &v| {
                    let (pu, pv) = (at(u), at(v));
                    let along = if pa.x != pb.x {
                        by_value(pu.x, pv.x)
                    } else {
                        by_value(pu.y, pv.y)
                    };
                    if (pa.x != pb.x && pb.x < pa.x) || (pa.x == pb.x && pb.y < pa.y) {
                        along.reverse()
                    } else {
                        along
                    }
                });
                walk.extend_from_slice(list);
            }
        }
        // 2. Split at repeated vertices: each return to a vertex on the
        // stack closes the loop above it.
        let mut stack: Vec<u32> = Vec::with_capacity(walk.len());
        let mut position: HashMap<u32, usize> = HashMap::new();
        for &v in &walk {
            if let Some(&i) = position.get(&v) {
                for w in &stack[i + 1..] {
                    position.remove(w);
                }
                let vertices = stack.split_off(i);
                stack.push(v);
                pieces.push(Piece {
                    ring: source[r],
                    vertices,
                });
            } else {
                position.insert(v, stack.len());
                stack.push(v);
            }
        }
        pieces.push(Piece {
            ring: source[r],
            vertices: stack,
        });
    }
    if let Some(piece) = pieces.iter().find(|p| p.vertices.len() < 3) {
        return Err(GeomError::Degenerate(format!(
            "profile {} splits into a loop of {} vertices at a pinch",
            name(piece.ring),
            piece.vertices.len()
        )));
    }

    // Each loop is turned on its own: a ring that crosses itself at a
    // vertex (a figure eight through a corner, as even-odd fill reports a
    // pinch) bounds the same two lobes as one that touches there.
    let depth = nest(points, &pieces)?;
    for (piece, &d) in pieces.iter_mut().zip(&depth) {
        let ring: Vec<Point2> = piece.vertices.iter().map(|&v| at(v)).collect();
        if ring_turns_left(&ring) != (d % 2 == 0) {
            piece.vertices.reverse();
        }
    }

    // 3. One node per visit, then one node per wedge at shared vertices.
    let mut nodes: Vec<Node> = Vec::new();
    for piece in &pieces {
        let first = nodes.len();
        let n = piece.vertices.len();
        for (k, &vertex) in piece.vertices.iter().enumerate() {
            nodes.push(Node {
                vertex,
                prev: first + (k + n - 1) % n,
                next: first + (k + 1) % n,
            });
        }
    }
    let mut visits: HashMap<u32, Vec<usize>> = HashMap::new();
    for (index, node) in nodes.iter().enumerate() {
        visits.entry(node.vertex).or_default().push(index);
    }
    let mut shared: Vec<(u32, Vec<usize>)> =
        visits.into_iter().filter(|(_, v)| v.len() > 1).collect();
    shared.sort_unstable_by_key(|(vertex, _)| *vertex);
    for (vertex, at_vertex) in &shared {
        join_wedges(points, &mut nodes, *vertex, at_vertex)?;
    }

    // The boundary cycles.
    let mut cycle_of = vec![usize::MAX; nodes.len()];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for start in 0..nodes.len() {
        if cycle_of[start] != usize::MAX {
            continue;
        }
        let mut cycle = Vec::new();
        let mut node = start;
        while cycle_of[node] == usize::MAX {
            cycle_of[node] = cycles.len();
            cycle.push(node);
            node = nodes[node].next;
        }
        cycles.push(cycle);
    }
    // A cycle bounds its part from outside when every visit of its
    // lexicographically smallest vertex turns strictly left: every other
    // point lies to the right of or above it, so the part's wedges there
    // are convex, while the outside of a hole cycle wraps round it. A cycle
    // of hole edges can bound from outside too, where holes touching in a
    // ring close off a part of the region.
    let outer: Vec<bool> = cycles
        .iter()
        .map(|cycle| {
            let lowest = cycle
                .iter()
                .map(|&n| at(nodes[n].vertex))
                .min_by(|&a, &b| lexicographic(a, b))
                .unwrap_or(Point2::new(0.0, 0.0));
            cycle
                .iter()
                .filter(|&&n| at(nodes[n].vertex) == lowest)
                .all(|&n| {
                    let node = nodes[n];
                    orient(
                        at(nodes[node.prev].vertex),
                        lowest,
                        at(nodes[node.next].vertex),
                    ) > 0
                })
        })
        .collect();

    // 4. Each hole cycle goes to the innermost outer cycle around it,
    // read at one of its vertices off that cycle: rings do not cross, so
    // such a vertex is off its boundary and places the whole hole cycle.
    // (A container never shares a vertex with its hole: they would have
    // joined into one cycle.)
    let rings_of: Vec<Vec<Point2>> = cycles
        .iter()
        .map(|cycle| cycle.iter().map(|&n| at(nodes[n].vertex)).collect())
        .collect();
    let vertices_of: Vec<HashSet<u32>> = cycles
        .iter()
        .map(|cycle| cycle.iter().map(|&n| nodes[n].vertex).collect())
        .collect();
    // Whether cycle `o` strictly contains cycle `c`, read at a vertex of
    // `c` that `o` does not visit.
    let around = |o: usize, c: usize| {
        cycles[c]
            .iter()
            .map(|&n| nodes[n].vertex)
            .find(|v| !vertices_of[o].contains(v))
            .is_some_and(|v| strictly_inside(&rings_of[o], at(v)))
    };
    let mut holes_of: Vec<Vec<usize>> = vec![Vec::new(); cycles.len()];
    for h in (0..cycles.len()).filter(|&h| !outer[h]) {
        let candidates: Vec<usize> = (0..cycles.len())
            .filter(|&o| outer[o] && around(o, h))
            .collect();
        // Candidates nest; the innermost contains none of the others.
        let container = candidates
            .iter()
            .copied()
            .find(|&o| candidates.iter().all(|&d| d == o || !around(o, d)));
        let Some(o) = container else {
            return Err(GeomError::Degenerate(
                "profile pinch leaves a hole cycle outside every outer cycle".to_owned(),
            ));
        };
        holes_of[o].push(h);
    }

    let mut triangles = Vec::new();
    let mut hole_count = 0;
    let mut outer_count = 0;
    for (o, cycle) in cycles.iter().enumerate() {
        if !outer[o] {
            continue;
        }
        outer_count += 1;
        hole_count += holes_of[o].len();
        let loops: Vec<Loop> = core::iter::once(cycle)
            .chain(holes_of[o].iter().map(|&h| &cycles[h]))
            .map(|c| Loop {
                vertices: c.iter().map(|&n| nodes[n].vertex).collect(),
            })
            .collect();
        let mut polygon = bridge::bridge_holes(points, &loops)?;
        triangles.extend(clip::clip_ears(points, &mut polygon)?);
    }
    let loops: Vec<Loop> = pieces
        .into_iter()
        .map(|piece| Loop {
            vertices: piece.vertices,
        })
        .collect();
    certify(points, &loops, &triangles, outer_count, hole_count)?;
    Ok(triangles)
}

/// Each piece's depth: how many other pieces contain it. Also refuses a
/// hole outside the outer ring or inside another hole, as validation does
/// for rings that touch nowhere.
///
/// Pieces do not cross, so one point of a piece off another's boundary
/// places it; when every vertex of piece `j` is a vertex of piece `i`,
/// `j`'s first edge decides by the side of `i`'s corner it leaves into.
fn nest(points: &[Point2], pieces: &[Piece]) -> GeomResult<Vec<usize>> {
    let at = |v: u32| points[v as usize];
    let rings: Vec<Vec<Point2>> = pieces
        .iter()
        .map(|p| p.vertices.iter().map(|&v| at(v)).collect())
        .collect();
    let boxes: Vec<(Point2, Point2)> = rings
        .iter()
        .map(|ring| {
            ring.iter().fold((ring[0], ring[0]), |(lo, hi), p| {
                (
                    Point2::new(lo.x.min(p.x), lo.y.min(p.y)),
                    Point2::new(hi.x.max(p.x), hi.y.max(p.y)),
                )
            })
        })
        .collect();
    let counter_clockwise: Vec<bool> = rings.iter().map(|r| ring_turns_left(r)).collect();
    let contains = |i: usize, j: usize| -> bool {
        let ((lo_i, hi_i), (lo_j, hi_j)) = (boxes[i], boxes[j]);
        if lo_j.x < lo_i.x || lo_j.y < lo_i.y || hi_j.x > hi_i.x || hi_j.y > hi_i.y {
            return false;
        }
        let own: HashMap<u32, usize> = pieces[i]
            .vertices
            .iter()
            .enumerate()
            .map(|(k, &v)| (v, k))
            .collect();
        if let Some(&v) = pieces[j].vertices.iter().find(|v| !own.contains_key(v)) {
            return strictly_inside(&rings[i], at(v));
        }
        // Every vertex shared: read the sector of `i` at `j`'s first vertex,
        // with `i`'s interior on its left.
        let (v, w) = (pieces[j].vertices[0], pieces[j].vertices[1]);
        let n = pieces[i].vertices.len();
        let k = own[&v];
        let (mut before, mut after) = (
            pieces[i].vertices[(k + n - 1) % n],
            pieces[i].vertices[(k + 1) % n],
        );
        if !counter_clockwise[i] {
            core::mem::swap(&mut before, &mut after);
        }
        locally_inside(at(before), at(v), at(after), at(w))
    };
    let mut depth = vec![0usize; pieces.len()];
    let mut inside: Vec<Vec<usize>> = vec![Vec::new(); pieces.len()];
    for j in 0..pieces.len() {
        for i in 0..pieces.len() {
            if i != j && contains(i, j) {
                depth[j] += 1;
                inside[j].push(i);
            }
        }
    }
    for (j, piece) in pieces.iter().enumerate() {
        if piece.ring == 0 {
            continue;
        }
        let h = piece.ring - 1;
        if !inside[j].iter().any(|&i| pieces[i].ring == 0) {
            return refuse(format!("profile hole {h} lies outside the outer ring"));
        }
        if let Some(&i) = inside[j]
            .iter()
            .find(|&&i| pieces[i].ring != 0 && pieces[i].ring != piece.ring)
        {
            return refuse(format!(
                "profile hole {h} lies inside hole {}",
                pieces[i].ring - 1
            ));
        }
    }
    Ok(depth)
}

/// Relink the nodes visiting `vertex` so that each one is a wedge of the
/// region: arriving along one edge and leaving along the next one
/// clockwise from it. Refuses the vertex when leaving and arriving edges do
/// not alternate around it, which is where rings cross.
fn join_wedges(
    points: &[Point2],
    nodes: &mut [Node],
    vertex: u32,
    at_vertex: &[usize],
) -> GeomResult<()> {
    let at = |v: u32| points[v as usize];
    let centre = at(vertex);
    // (far end, leaves here, node)
    let mut spokes: Vec<(Point2, bool, usize)> = Vec::with_capacity(2 * at_vertex.len());
    for &n in at_vertex {
        spokes.push((at(nodes[nodes[n].next].vertex), true, n));
        spokes.push((at(nodes[nodes[n].prev].vertex), false, n));
    }
    // Counter-clockwise from the positive x axis: the upper half-plane
    // first (the axis itself included), then the lower one.
    let half = |p: Point2| u8::from(!(p.y > centre.y || (p.y == centre.y && p.x > centre.x)));
    let mut tie = false;
    spokes.sort_by(|a, b| {
        half(a.0)
            .cmp(&half(b.0))
            .then_with(|| match orient(centre, a.0, b.0) {
                1 => Ordering::Less,
                -1 => Ordering::Greater,
                _ => {
                    // Two edges leaving one way overlap, which validation
                    // refuses; kept as a refusal rather than an arbitrary order.
                    tie = true;
                    Ordering::Equal
                }
            })
    });
    let count = spokes.len();
    let alternates = (0..count).all(|i| spokes[i].1 != spokes[(i + 1) % count].1);
    if tie || !alternates {
        return refuse(format!(
            "profile rings cross at their shared vertex ({}, {})",
            centre.x, centre.y
        ));
    }
    let next: Vec<usize> = at_vertex.iter().map(|&n| nodes[n].next).collect();
    let original_next = |n: usize| next[at_vertex.iter().position(|&m| m == n).unwrap_or(0)];
    for i in 0..count {
        let (leaving, arriving) = (spokes[i], spokes[(i + 1) % count]);
        if !leaving.1 {
            continue;
        }
        // The wedge from `leaving` counter-clockwise to `arriving`: the
        // arriving node now leaves along `leaving`'s edge.
        let (into, out) = (arriving.2, original_next(leaving.2));
        nodes[into].next = out;
        nodes[out].prev = into;
    }
    Ok(())
}
