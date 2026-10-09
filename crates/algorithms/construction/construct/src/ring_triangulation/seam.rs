//! Rings that join a hole by a seam traversed both ways, under
//! [`PinchPolicy::Accept`](super::PinchPolicy) (#270).
//!
//! A planar face can carry its hole in one ring: the outer boundary, a
//! seam to the hole, the hole, and the same seam back -- a keyhole. The
//! seam is one segment used once in each direction, so as a boundary it
//! cancels: the ring bounds exactly what its other edges bound, the outer
//! polygon minus the hole. It is the weakly simple ring that
//! [`bridge`](super::bridge) builds before ear clipping, arriving with the
//! bridge already in place. Validation alone refuses it, because the seam
//! and its reverse overlap.
//!
//! The ring is rebuilt without its seams instead, every decision exact:
//!
//! 1. A seam ([`find`]) is a pair of edges of one ring with the same two
//!    endpoints, coordinates compared by value, in opposite directions,
//!    when no other edge of any ring has those endpoints. A ring's seams
//!    must nest like brackets along it; two that interleave are not
//!    seams, and their edges stay to be refused as an overlap.
//! 2. No other edge, and no other seam, may meet a seam except at its
//!    endpoints ([`check_clear`]): an edge that crosses a seam, runs along
//!    it or ends inside it is refused by name.
//! 3. Removing a seam's two edges splits its ring in two ([`split`]): the
//!    edges strictly between them, and the rest. Each part closes on its
//!    own, since the first copy ends where the second starts, and seams
//!    nested inside split it further. A part of one or two vertices, or
//!    one folding back where a seam left it, has no area and is refused.
//!    So is a part of no edges, seams in a row (`a -> b -> c` out and
//!    back): their shared end `b` would be in no part, and a triangle edge
//!    could run past it, a T-junction for the faces that share it.
//! 4. The parts are triangulated as rings of their own on the
//!    [`pinch`](super::pinch) path, named after the ring they came from.
//!    That path orients each loop by its nesting depth, so a hole drawn
//!    either way round is a hole, and a seam joining two loops side by
//!    side leaves two parts, each bounded from outside. Parts that touch
//!    at points are accepted there as any pinch is; parts that cross or
//!    overlap are refused.
//!
//! The certificate is the pinch path's over the parts: positive triangles
//! whose boundary is exactly the parts' edges, every other edge twinned,
//! and `n + 2h - 2c` triangles for the parts' `n` vertices. A seam's two
//! copies cancel as boundary chains, so the parts' edges bound what the
//! ring bounds and the certificate is the ring's. A seam's endpoints are
//! repeated vertices: each part keeps its own visit, and a triangle uses
//! the point's first index in `outer ++ holes`, as at a pinch. A seam may
//! be a triangle edge, used once in each direction, or not one at all;
//! some regions force it (a thin rim seamed at its corner), and faces
//! sharing the seam then meet four triangle uses of it.
//!
//! [`PinchPolicy::Refuse`](super::PinchPolicy) refuses a ring with a seam
//! by name ([`refuse_for_solid`]): an extruded keyhole puts two
//! coincident wall faces on the seam.

use std::collections::HashMap;

use core::ops::Range;

use axiolid_contracts::GeomResult;
use axiolid_core::Point2;

use super::validate::{check_edges, folds_back, name, refuse, Touch};
use super::{by_value, orient, pinch, segments_touch, within, PinchPolicy};

/// A seam of ring `ring`: its edge `open` and its later edge `close` run
/// between the same two points in opposite directions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Seam {
    pub(super) ring: usize,
    pub(super) open: usize,
    pub(super) close: usize,
}

/// A point's coordinates as a hash key.
type Key = (u64, u64);

/// A point's key, `-0.0` the same as `0.0`.
fn key(p: Point2) -> Key {
    ((p.x + 0.0).to_bits(), (p.y + 0.0).to_bits())
}

/// Every seam of `rings`, in ring order and then by `open`.
pub(super) fn find(points: &[Point2], rings: &[Range<usize>]) -> Vec<Seam> {
    // Each undirected segment's uses: (ring, edge, from its lower key).
    type Use = (usize, usize, bool);
    let mut uses: HashMap<(Key, Key), Vec<Use>> = HashMap::new();
    for (r, ring) in rings.iter().enumerate() {
        let n = ring.len();
        for k in 0..n {
            let (a, b) = (
                key(points[ring.start + k]),
                key(points[ring.start + (k + 1) % n]),
            );
            let forward = a < b;
            let segment = if forward { (a, b) } else { (b, a) };
            uses.entry(segment).or_default().push((r, k, forward));
        }
    }
    let mut pairs: Vec<Seam> = uses
        .into_values()
        .filter_map(|list| match list[..] {
            [(r, i, fi), (s, j, fj)] if r == s && fi != fj => Some(Seam {
                ring: r,
                open: i.min(j),
                close: i.max(j),
            }),
            _ => None,
        })
        .collect();
    pairs.sort_unstable_by_key(|s| (s.ring, s.open));

    // Seams must nest like brackets along their ring: walking its edges,
    // a seam closing while a later one is still open interleaves with it,
    // and with every other seam opened after it and still open.
    let mut good = vec![true; pairs.len()];
    let mut events: Vec<(usize, usize, usize)> = Vec::with_capacity(2 * pairs.len());
    for (id, seam) in pairs.iter().enumerate() {
        events.push((seam.ring, seam.open, id));
        events.push((seam.ring, seam.close, id));
    }
    events.sort_unstable();
    let mut stack: Vec<usize> = Vec::new();
    let mut ring = usize::MAX;
    for (r, k, id) in events {
        if r != ring {
            ring = r;
            stack.clear();
        }
        if pairs[id].open == k {
            stack.push(id);
            continue;
        }
        let Some(at) = stack.iter().rposition(|&other| other == id) else {
            continue;
        };
        if at + 1 != stack.len() {
            for &other in &stack[at..] {
                good[other] = false;
            }
        }
        stack.remove(at);
    }
    pairs
        .into_iter()
        .zip(good)
        .filter_map(|(seam, good)| good.then_some(seam))
        .collect()
}

/// The refusal of a ring with a seam for a solid cap.
pub(super) fn refuse_for_solid<T>(seam: &Seam) -> GeomResult<T> {
    refuse(format!(
        "profile {} runs along a seam both ways (edges {} and {}), which a solid cannot extrude",
        name(seam.ring),
        seam.open,
        seam.close
    ))
}

/// Triangulate rings with seams: [`check_clear`], [`split`], then the
/// pinch path over the parts. Triangles use each point's first index in
/// `points`.
pub(super) fn triangulate(
    points: &[Point2],
    rings: &[Range<usize>],
    seams: &[Seam],
) -> GeomResult<Vec<[u32; 3]>> {
    let (parts, touches) = validated_parts(points, rings, seams)?;
    let identity = |n: usize| (0..n as u32).collect::<Vec<u32>>();
    let canon =
        pinch::canonical(&parts.points, &touches).unwrap_or_else(|| identity(parts.points.len()));
    let triangles =
        pinch::triangulate(&parts.points, &parts.rings, &parts.source, &canon, &touches)?;
    // Back to `points`: a part's point is a point of `points`, and each
    // point to its first index there. Both maps keep distinct points
    // distinct, so the certified triangulation is only relabelled.
    let first = pinch::canonical(points, &[]).unwrap_or_else(|| identity(points.len()));
    Ok(triangles
        .into_iter()
        .map(|t| t.map(|v| first[parts.back[v as usize] as usize]))
        .collect())
}

/// [`check_clear`], [`split`], then the parts' own validation: the parts
/// and the vertices lying inside their edges.
fn validated_parts(
    points: &[Point2],
    rings: &[Range<usize>],
    seams: &[Seam],
) -> GeomResult<(Parts, Vec<Touch>)> {
    check_clear(points, rings, seams)?;
    let parts = split(points, rings, seams)?;
    let touches = check_edges(
        &parts.points,
        &parts.rings,
        &parts.source,
        PinchPolicy::Accept,
    )?;
    Ok((parts, touches))
}

/// The vertices lying inside another edge of rings with seams, as
/// `(ring, edge, vertex)` of the input, which [`triangulate`] inserts into
/// those edges (#265). Every edge of a part is an edge of its source ring
/// that is not a seam, and no vertex lies inside a seam.
pub(super) fn touches(
    points: &[Point2],
    rings: &[Range<usize>],
    seams: &[Seam],
) -> GeomResult<Vec<(usize, usize, usize)>> {
    let (parts, touches) = validated_parts(points, rings, seams)?;
    Ok(touches
        .into_iter()
        .map(|touch| {
            let source = parts.source[touch.ring];
            let start = parts.back[parts.rings[touch.ring].start + touch.index] as usize;
            (
                source,
                start - rings[source].start,
                parts.back[touch.vertex as usize] as usize,
            )
        })
        .collect())
}

/// Each seam edge by `(ring, edge)`: `true` for its first copy, `false`
/// for its second.
fn roles(seams: &[Seam]) -> HashMap<(usize, usize), bool> {
    let mut role = HashMap::with_capacity(2 * seams.len());
    for seam in seams {
        role.insert((seam.ring, seam.open), true);
        role.insert((seam.ring, seam.close), false);
    }
    role
}

/// Refuse an edge or another seam meeting a seam anywhere but at the
/// seam's endpoints.
///
/// Edges are swept in order of their smallest x, as in
/// [`check_edges`], and only pairs with a seam among them are tested.
pub(super) fn check_clear(
    points: &[Point2],
    rings: &[Range<usize>],
    seams: &[Seam],
) -> GeomResult<()> {
    // The first copy stands for the seam; the second is left out.
    let role = roles(seams);
    // (a, b, ring, is a seam)
    let mut segments: Vec<(Point2, Point2, usize, bool)> = Vec::with_capacity(points.len());
    for (r, ring) in rings.iter().enumerate() {
        let n = ring.len();
        for k in 0..n {
            let is_seam = match role.get(&(r, k)) {
                Some(false) => continue,
                Some(true) => true,
                None => false,
            };
            segments.push((
                points[ring.start + k],
                points[ring.start + (k + 1) % n],
                r,
                is_seam,
            ));
        }
    }
    segments.sort_by(|e, f| by_value(e.0.x.min(e.1.x), f.0.x.min(f.1.x)));
    for (i, e) in segments.iter().enumerate() {
        let right = e.0.x.max(e.1.x);
        let (low, high) = (e.0.y.min(e.1.y), e.0.y.max(e.1.y));
        for f in &segments[i + 1..] {
            if f.0.x.min(f.1.x) > right {
                break;
            }
            if !(e.3 || f.3) || f.0.y.max(f.1.y) < low || f.0.y.min(f.1.y) > high {
                continue;
            }
            let (seam, other) = if e.3 { (e, f) } else { (f, e) };
            let clear = meets_only_at_ends(seam.0, seam.1, other.0, other.1)
                && (!other.3 || meets_only_at_ends(other.0, other.1, seam.0, seam.1));
            if !clear {
                let what = if seam.2 == other.2 {
                    "itself".to_owned()
                } else {
                    name(other.2)
                };
                return refuse(format!(
                    "profile {} has a seam from ({}, {}) to ({}, {}) that crosses or overlaps {what}",
                    name(seam.2),
                    seam.0.x,
                    seam.0.y,
                    seam.1.x,
                    seam.1.y
                ));
            }
        }
    }
    Ok(())
}

/// Whether segment `c d` meets the seam `a b` nowhere but at `a` or `b`.
fn meets_only_at_ends(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    if !segments_touch(a, b, c, d) {
        return true;
    }
    if orient(a, b, c) == 0 && orient(a, b, d) == 0 {
        // On one line: only an end shared, with the far ends on opposite
        // sides of it, touches in a single point.
        for (p, far_seam) in [(a, b), (b, a)] {
            for (q, far_other) in [(c, d), (d, c)] {
                if p == q {
                    return far_seam != far_other && !folds_back(far_seam, p, far_other);
                }
            }
        }
        return false;
    }
    // Not on one line: they meet in one point, which must be `a` or `b`.
    let on = |p: Point2| orient(c, d, p) == 0 && within(c, d, p);
    on(a) || on(b)
}

/// Rings split at their seams: the parts' points concatenated, each
/// part's range and source ring, and each point's index in the input.
#[derive(Debug, Clone)]
pub(super) struct Parts {
    pub(super) points: Vec<Point2>,
    pub(super) rings: Vec<Range<usize>>,
    pub(super) source: Vec<usize>,
    pub(super) back: Vec<u32>,
}

/// Split every ring at its seams into the parts the seams separate.
pub(super) fn split(
    points: &[Point2],
    rings: &[Range<usize>],
    seams: &[Seam],
) -> GeomResult<Parts> {
    let role = roles(seams);
    let mut parts = Parts {
        points: Vec::with_capacity(points.len()),
        rings: Vec::new(),
        source: Vec::new(),
        back: Vec::with_capacity(points.len()),
    };
    for (r, ring) in rings.iter().enumerate() {
        // Edge lists: the ring's remainder at the bottom, one list per
        // seam still open above it, with the seam's first edge. Each edge
        // is named by its start.
        let n = ring.len();
        let mut open: Vec<(usize, Vec<usize>)> = vec![(n, Vec::new())];
        let mut closed: Vec<Vec<usize>> = Vec::new();
        for k in 0..n {
            match role.get(&(r, k)) {
                Some(true) => open.push((k, Vec::new())),
                Some(false) => {
                    let Some((first, edges)) = open.pop() else {
                        continue;
                    };
                    if edges.is_empty() {
                        // Only seams between the two copies: seams in a
                        // row, whose shared end is in no part, so a
                        // triangle edge could run past it.
                        return refuse(format!(
                            "profile {} has seams in a row through vertex {}",
                            name(r),
                            (first + 1) % n
                        ));
                    }
                    closed.push(edges);
                }
                None => {
                    if let Some((_, top)) = open.last_mut() {
                        top.push(k);
                    }
                }
            }
        }
        closed.extend(open.into_iter().map(|(_, edges)| edges));
        for edges in closed {
            if edges.len() < 3 {
                return refuse(format!(
                    "profile {} splits at a seam into a loop of {} vertices",
                    name(r),
                    edges.len()
                ));
            }
            let m = edges.len();
            for (i, &k) in edges.iter().enumerate() {
                let at = |j: usize| points[ring.start + edges[(i + j) % m]];
                if folds_back(at(m - 1), at(0), at(1)) {
                    return refuse(format!(
                        "profile {} folds back on itself at vertex {k} once its seam is removed",
                        name(r)
                    ));
                }
            }
            let start = parts.points.len();
            for k in edges {
                parts.points.push(points[ring.start + k]);
                parts.back.push((ring.start + k) as u32);
            }
            parts.rings.push(start..parts.points.len());
            parts.source.push(r);
        }
    }
    Ok(parts)
}
