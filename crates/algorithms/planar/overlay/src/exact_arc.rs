//! Exact arc-aware boolean (#155, ADR 0070).
//!
//! Every topological decision — where boundaries cross, the order of
//! crossings along an edge, whether a piece of one boundary lies inside,
//! outside or on the other region, how pieces link into rings, and which
//! ring is a hole of which — is an exact sign computed by `axiolid-exact`
//! (interval filter first, exact arithmetic only when the filter cannot
//! decide). No tolerance enters any of them, so the answer cannot depend
//! on drawing units.
//!
//! Constructed crossing points carry a square root; they are rounded to
//! `f64` exactly once, when the output rings are written.
//!
//! # Steps
//!
//! 1. Split every edge at every point where it meets an edge of the other
//!    operand (overlaps contribute their end points, so overlapping pieces
//!    end up identical).
//! 2. Pick an exact rational point strictly inside each piece.
//! 3. Classify that point against the other operand: on its boundary
//!    (shared piece, same or opposite direction), inside, or outside. The
//!    inside test is a crossing-number count over y-monotone pieces with
//!    the half-open rule, so vertices and tangencies need no special case.
//! 4. Keep pieces by the operation's rule, reversing the ones that bound
//!    the result from the other side.
//! 5. Link kept pieces into rings, taking the rightmost turn at vertices
//!    where several leave, which yields minimal rings.
//! 6. Counter-clockwise rings are outer boundaries; each clockwise ring is
//!    a hole of the smallest outer that contains it.

mod edge;
mod point;

use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;

use crate::arc::{arc_ring_area, ArcRing, ArcVertex};
use crate::{OverlayError, OverlayOperation};

use edge::{crossings, Carrier, Edge};
use point::{cmp_y, dy, orient, same_point, sign, Circle, Pred, Tangent, XPoint};

/// Which operand a piece came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Subject,
    Clip,
}

/// Where a piece lies relative to the other operand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Status {
    Inside,
    Outside,
    SharedSame,
    SharedOpposite,
}

/// A piece of an input edge between consecutive split points.
#[derive(Debug, Clone)]
struct Piece {
    from: XPoint,
    to: XPoint,
    /// An exact point strictly inside the piece.
    sample: XPoint,
    tangent: Tangent,
    /// `Some` for arcs: the circle and the travel turn (after reversal).
    arc: Option<(Circle, Sign)>,
    /// The original bulge when the piece is a whole input edge.
    whole_bulge: Option<f64>,
}

impl Piece {
    fn reversed(self) -> Self {
        let tangent = match self.tangent {
            Tangent::Fixed(x, y) => Tangent::Fixed(x.neg(), y.neg()),
            Tangent::Circle(c, turn) => Tangent::Circle(c, turn.flip()),
        };
        Self {
            from: self.to,
            to: self.from,
            sample: self.sample,
            tangent,
            arc: self.arc.map(|(c, turn)| (c, turn.flip())),
            whole_bulge: self.whole_bulge.map(|b| -b),
        }
    }

    fn bulge(&self) -> f64 {
        let Some((circle, turn)) = &self.arc else {
            return 0.0;
        };
        if let Some(bulge) = self.whole_bulge {
            return bulge;
        }
        let c = circle.approx_centre();
        let angle = |p: Point2| (p.y - c.y).atan2(p.x - c.x);
        let t = if *turn == Sign::Negative { -1.0 } else { 1.0 };
        let tau = std::f64::consts::TAU;
        let a0 = angle(self.from.approx());
        let sweep = (t * (angle(self.to.approx()) - a0)).rem_euclid(tau);
        let to_sample = (t * (angle(self.sample.approx()) - a0)).rem_euclid(tau);
        // A sub-arc too short to resolve in f64 can wrap to nearly a full
        // turn; the sample, which lies inside the piece, tells them apart.
        let sweep = if to_sample > sweep { 0.0 } else { sweep };
        t * (sweep / 4.0).tan()
    }
}

/// One y-monotone part of an input edge, for crossing counts.
struct Mono {
    a: XPoint,
    b: XPoint,
    arc: Option<(Circle, Sign)>,
}

fn edges_of(ring: &ArcRing) -> Vec<Edge> {
    let n = ring.vertices.len();
    (0..n)
        .map(|i| {
            let from = ring.vertices[i];
            let to = ring.vertices[(i + 1) % n];
            Edge::new(from.point, to.point, from.bulge)
        })
        .collect()
}

/// Split an edge into y-monotone parts at the circle's top and bottom.
fn monotone(edge: &Edge) -> Vec<Mono> {
    let Carrier::Arc { circle, turn, .. } = &edge.carrier else {
        return vec![Mono {
            a: edge.p0.clone(),
            b: edge.p1.clone(),
            arc: None,
        }];
    };
    let mut cuts: Vec<XPoint> = [Sign::Positive, Sign::Negative]
        .into_iter()
        .map(|which| circle.extreme(which))
        .filter(|x| edge.holds(x) && !same_point(x, &edge.p0) && !same_point(x, &edge.p1))
        .collect();
    if cuts.len() == 2 && edge.order(&cuts[0], &cuts[1]) == Sign::Positive {
        cuts.swap(0, 1);
    }
    let mut stops = vec![edge.p0.clone()];
    stops.extend(cuts);
    stops.push(edge.p1.clone());
    stops
        .windows(2)
        .map(|w| Mono {
            a: w[0].clone(),
            b: w[1].clone(),
            arc: Some((circle.clone(), *turn)),
        })
        .collect()
}

/// Which side of a monotone part `s` lies on, relative to travel `a -> b`:
/// positive is left. Only called for `s` strictly within the part's height
/// range (half-open) and not on it.
fn side_of(part: &Mono, s: &XPoint) -> Sign {
    let chord = orient(&part.a, &part.b, s);
    let Some((circle, turn)) = &part.arc else {
        return chord;
    };
    // The arc lies on side -turn of its chord, and the region between
    // chord and arc is the disc restricted to that side: a monotone part is
    // at most a semicircle.
    let bulge = turn.flip();
    if chord != bulge && chord != Sign::Zero {
        return chord;
    }
    match sign(Pred::OnCircle(s, circle)) {
        Sign::Negative => bulge.flip(),
        Sign::Positive => bulge,
        _ => Sign::Zero,
    }
}

/// Winding number of `other`'s boundary around `s`, which is not on it.
fn winding(s: &XPoint, parts: &[Mono]) -> i64 {
    let mut total = 0;
    for part in parts {
        let ya = cmp_y(&part.a, s);
        let yb = cmp_y(&part.b, s);
        let up = ya != Sign::Positive && yb == Sign::Positive;
        let down = yb != Sign::Positive && ya == Sign::Positive;
        if !(up || down) {
            continue;
        }
        let side = side_of(part, s);
        if up && side == Sign::Positive {
            total += 1;
        } else if down && side == Sign::Negative {
            total -= 1;
        }
    }
    total
}

fn tangent_of(edge: &Edge) -> Tangent {
    match &edge.carrier {
        Carrier::Segment => {
            let (x, y) = edge_direction(edge);
            Tangent::Fixed(x, y)
        }
        Carrier::Arc { circle, turn, .. } => Tangent::Circle(circle.clone(), *turn),
    }
}

fn edge_direction(edge: &Edge) -> (Dyadic, Dyadic) {
    let a = edge.p0.approx();
    let b = edge.p1.approx();
    // Input vertices are exact doubles, so this difference is exact in
    // dyadic arithmetic.
    (dy(b.x).sub(&dy(a.x)), dy(b.y).sub(&dy(a.y)))
}

/// Approximate edge parameter of a point on the edge, by `f64` bisection
/// on the approximate positions. Only a seed: the sample is verified
/// exactly.
fn approx_param(edge: &Edge, p: Point2) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    let p0 = edge.p0.approx();
    let turn = match &edge.carrier {
        Carrier::Segment => None,
        Carrier::Arc { turn, .. } => Some(if *turn == Sign::Negative { -1.0 } else { 1.0 }),
    };
    let before = |q: Point2| -> bool {
        // Is q (at some parameter) before p along the edge?
        match turn {
            None => {
                let d = edge.p1.approx();
                let (dx, dy) = (d.x - p0.x, d.y - p0.y);
                (q.x - p0.x) * dx + (q.y - p0.y) * dy < (p.x - p0.x) * dx + (p.y - p0.y) * dy
            }
            Some(t) => {
                let o = (q.x - p0.x) * (p.y - p0.y) - (q.y - p0.y) * (p.x - p0.x);
                t * o > 0.0
            }
        }
    };
    for _ in 0..44 {
        let mid = 0.5 * (lo + hi);
        if before(edge.approx_at(mid)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// An exact rational point strictly between `from` and `to` on `edge`.
///
/// Seeded at the midpoint of the pieces' approximate parameters, rounded
/// to a short dyadic so the sample stays cheap in every later predicate;
/// verified exactly, with exact bisection over the whole edge as the
/// fallback when the seed is not strictly inside.
fn sample(edge: &Edge, from: &XPoint, to: &XPoint) -> XPoint {
    let inside =
        |p: &XPoint| edge.order(p, from) == Sign::Positive && edge.order(p, to) == Sign::Negative;
    let (a, b) = (
        approx_param(edge, from.approx()),
        approx_param(edge, to.approx()),
    );
    let mid = 0.5 * (a + b);
    for bits in [12, 24, 40] {
        let scale = f64::from(1u32 << (bits / 2)) * f64::from(1u32 << (bits - bits / 2));
        let u = (mid * scale).round() / scale;
        if u > 0.0 && u < 1.0 {
            let p = edge.point_at(&dy(u));
            if inside(&p) {
                return p;
            }
        }
    }
    // The piece is too short for the f64 seed (constructed crossings can lie
    // within rounding of a vertex). Bracket it near the seed, verified
    // exactly, so the exact bisection below runs a handful of steps instead
    // of resolving the whole edge from [0, 1].
    let bracket = |target: f64, want_before: bool| -> Dyadic {
        let fallback = if want_before { dy(0.0) } else { dy(1.0) };
        let mut step = f64::EPSILON;
        while step < 1.0 {
            let u = if want_before {
                target - step
            } else {
                target + step
            };
            if !(0.0..=1.0).contains(&u) {
                return fallback;
            }
            let p = edge.point_at(&dy(u));
            let ok = if want_before {
                edge.order(&p, from) != Sign::Positive
            } else {
                edge.order(&p, to) != Sign::Negative
            };
            if ok {
                return dy(u);
            }
            step *= 16.0;
        }
        fallback
    };
    let (mut lo, mut hi) = (bracket(a.min(b), true), bracket(a.max(b), false));
    loop {
        let mid = lo.add(&hi).mul(&dy(0.5));
        let p = edge.point_at(&mid);
        if edge.order(&p, from) != Sign::Positive {
            lo = mid;
        } else if edge.order(&p, to) != Sign::Negative {
            hi = mid;
        } else {
            return p;
        }
    }
}

/// Split both operands' edges and return their pieces.
fn pieces(own: &[Edge], other: &[Edge], ring: &ArcRing) -> Vec<Piece> {
    let mut out = Vec::new();
    for (index, edge) in own.iter().enumerate() {
        let mut stops = vec![edge.p0.clone(), edge.p1.clone()];
        for theirs in other {
            for x in crossings(edge, theirs) {
                if !stops.iter().any(|y| same_point(y, &x)) {
                    stops.push(x);
                }
            }
        }
        stops.sort_by(|a, b| match edge.order(a, b) {
            Sign::Negative => std::cmp::Ordering::Less,
            Sign::Positive => std::cmp::Ordering::Greater,
            _ => std::cmp::Ordering::Equal,
        });
        let split = stops.len() > 2;
        let arc = match &edge.carrier {
            Carrier::Segment => None,
            Carrier::Arc { circle, turn, .. } => Some((circle.clone(), *turn)),
        };
        for w in stops.windows(2) {
            out.push(Piece {
                sample: sample(edge, &w[0], &w[1]),
                from: w[0].clone(),
                to: w[1].clone(),
                tangent: tangent_of(edge),
                arc: arc.clone(),
                whole_bulge: (!split && arc.is_some()).then(|| ring.vertices[index].bulge),
            });
        }
    }
    out
}

fn classify(piece: &Piece, other: &[Edge], parts: &[Mono]) -> Status {
    for edge in other {
        if edge.contains(&piece.sample) {
            let theirs = tangent_of(edge);
            let dot = sign(Pred::Tangents {
                at: &piece.sample,
                u: &piece.tangent,
                v: &theirs,
                cross: false,
            });
            return if dot == Sign::Positive {
                Status::SharedSame
            } else {
                Status::SharedOpposite
            };
        }
    }
    if winding(&piece.sample, parts) != 0 {
        Status::Inside
    } else {
        Status::Outside
    }
}

/// Whether a piece is kept, and whether reversed.
fn keep(operation: OverlayOperation, side: Side, status: Status) -> Option<bool> {
    use OverlayOperation as Op;
    use Status as S;
    let subject = side == Side::Subject;
    match (operation, status) {
        (Op::Union, S::Outside) => Some(false),
        (Op::Intersection, S::Inside) => Some(false),
        (Op::Union | Op::Intersection, S::SharedSame) if subject => Some(false),
        (Op::Difference, S::Outside) if subject => Some(false),
        (Op::Difference, S::Inside) if !subject => Some(true),
        (Op::Difference, S::SharedOpposite) if subject => Some(false),
        (Op::Xor, S::Outside) => Some(false),
        (Op::Xor, S::Inside) => Some(true),
        _ => None,
    }
}

/// Preference rank of leaving along `d` after arriving along `din`.
///
/// Result boundaries keep their region on the left, so the minimal ring is
/// traced by taking, at every vertex, the first leaving direction met when
/// rotating clockwise from "straight back": left turns, then straight on,
/// then right turns, then back. Lower rank wins.
fn turn_rank(at: &XPoint, din: &Tangent, d: &Tangent) -> u8 {
    let cross = sign(Pred::Tangents {
        at,
        u: din,
        v: d,
        cross: true,
    });
    match cross {
        Sign::Positive => 0,
        Sign::Negative => 2,
        _ => {
            let dot = sign(Pred::Tangents {
                at,
                u: din,
                v: d,
                cross: false,
            });
            if dot == Sign::Negative {
                3
            } else {
                1
            }
        }
    }
}

/// Link pieces into closed rings.
///
/// Two pieces leaving one vertex along the same tangent (curves touching
/// tangentially at a vertex of the result) are ranked equal; the first is
/// taken. Rings stay closed either way; only how touching rings are
/// grouped could differ, which is a known limit (ADR 0070).
fn link(mut pool: Vec<Piece>) -> Result<Vec<Vec<Piece>>, OverlayError> {
    let mut rings = Vec::new();
    while let Some(first) = pool.pop() {
        let start = first.from.clone();
        let mut ring = vec![first];
        loop {
            let last = ring.last().expect("ring has a first piece");
            if same_point(&last.to, &start) {
                break;
            }
            let at = last.to.clone();
            let din = last.tangent.clone();
            let mut best: Option<(usize, u8)> = None;
            for (index, cand) in pool.iter().enumerate() {
                if !same_point(&cand.from, &at) {
                    continue;
                }
                let rank = turn_rank(&at, &din, &cand.tangent);
                let better = match best {
                    None => true,
                    Some((_, best_rank)) if rank < best_rank => true,
                    // Within the left or the right half, the more
                    // counter-clockwise direction comes first.
                    Some((b, best_rank)) if rank == best_rank && rank != 1 && rank != 3 => {
                        sign(Pred::Tangents {
                            at: &at,
                            u: &pool[b].tangent,
                            v: &cand.tangent,
                            cross: true,
                        }) == Sign::Positive
                    }
                    _ => false,
                };
                if better {
                    best = Some((index, rank));
                }
            }
            // With simple operands every vertex has as many kept pieces
            // leaving as arriving; a dead end means an operand crosses
            // itself.
            let (index, _) = best.ok_or(OverlayError::SelfIntersection)?;
            ring.push(pool.swap_remove(index));
        }
        rings.push(ring);
    }
    Ok(rings)
}

fn to_ring(pieces: &[Piece]) -> ArcRing {
    ArcRing {
        vertices: pieces
            .iter()
            .map(|p| ArcVertex::bulged(p.from.approx(), p.bulge()))
            .collect(),
    }
}

/// Monotone parts of a ring made of pieces, for containment tests.
fn ring_parts(pieces: &[Piece]) -> Vec<Mono> {
    let mut out = Vec::new();
    for p in pieces {
        let Some((circle, turn)) = &p.arc else {
            out.push(Mono {
                a: p.from.clone(),
                b: p.to.clone(),
                arc: None,
            });
            continue;
        };
        // A piece is part of one arc; cut it at the circle's extremes that
        // fall strictly inside it. On a sub-arc from P to Q with travel
        // turn t, X lies strictly inside iff it is on the circle and on
        // side -t of the chord P -> Q.
        let mut cuts: Vec<XPoint> = [Sign::Positive, Sign::Negative]
            .into_iter()
            .map(|which| circle.extreme(which))
            .filter(|x| orient(&p.from, &p.to, x) == turn.flip())
            .collect();
        if cuts.len() == 2 && orient(&p.from, &cuts[0], &cuts[1]) != *turn {
            cuts.swap(0, 1);
        }
        let mut stops = vec![p.from.clone()];
        stops.extend(cuts);
        stops.push(p.to.clone());
        for w in stops.windows(2) {
            out.push(Mono {
                a: w[0].clone(),
                b: w[1].clone(),
                arc: Some((circle.clone(), *turn)),
            });
        }
    }
    out
}

/// The exact boolean of two validated rings, both counter-clockwise.
pub(crate) fn boolean(
    subject: &ArcRing,
    clip: &ArcRing,
    operation: OverlayOperation,
) -> Result<Vec<(ArcRing, Vec<ArcRing>)>, OverlayError> {
    let a_edges = edges_of(subject);
    let b_edges = edges_of(clip);
    let a_parts: Vec<Mono> = a_edges.iter().flat_map(monotone).collect();
    let b_parts: Vec<Mono> = b_edges.iter().flat_map(monotone).collect();

    let t0 = std::time::Instant::now();
    let mut kept = Vec::new();
    for (side, own, other, parts, ring) in [
        (Side::Subject, &a_edges, &b_edges, &b_parts, subject),
        (Side::Clip, &b_edges, &a_edges, &a_parts, clip),
    ] {
        for piece in pieces(own, other, ring) {
            let status = classify(&piece, other, parts);
            if let Some(reverse) = keep(operation, side, status) {
                kept.push(if reverse { piece.reversed() } else { piece });
            }
        }
    }

    let t1 = t0.elapsed();
    let c1: Vec<usize> = point::COUNTS
        .iter()
        .map(|c| c.swap(0, std::sync::atomic::Ordering::Relaxed))
        .collect();
    let rings = link(kept)?;
    let t2 = t0.elapsed();
    let c2: Vec<usize> = point::COUNTS
        .iter()
        .map(|c| c.swap(0, std::sync::atomic::Ordering::Relaxed))
        .collect();
    if std::env::var("ZZ").is_ok() {
        eprintln!("split+classify {t1:?} {c1:?} | link {:?} {c2:?}", t2 - t1);
    }
    let mut outers: Vec<(ArcRing, f64, Vec<Mono>)> = Vec::new();
    let mut holes: Vec<(ArcRing, XPoint)> = Vec::new();
    for pieces in &rings {
        let ring = to_ring(pieces);
        let area = arc_ring_area(&ring);
        if area > 0.0 {
            outers.push((ring, area, ring_parts(pieces)));
        } else {
            holes.push((ring, pieces[0].sample.clone()));
        }
    }
    let mut regions: Vec<(ArcRing, Vec<ArcRing>)> = outers
        .iter()
        .map(|(ring, _, _)| (ring.clone(), Vec::new()))
        .collect();
    for (hole, probe) in holes {
        let owner = outers
            .iter()
            .enumerate()
            .filter(|(_, (_, _, parts))| winding(&probe, parts) != 0)
            .min_by(|x, y| x.1 .1.total_cmp(&y.1 .1))
            .map(|(index, _)| index);
        // A hole with no containing outer cannot come out of a boolean of
        // simple rings; keep it rather than lose area silently.
        match owner {
            Some(index) => regions[index].1.push(hole),
            None => {
                if let Some(region) = regions.first_mut() {
                    region.1.push(hole);
                }
            }
        }
    }
    if std::env::var("ZZ").is_ok() {
        let c3: Vec<usize> = point::COUNTS
            .iter()
            .map(|c| c.swap(0, std::sync::atomic::Ordering::Relaxed))
            .collect();
        eprintln!("  nest {:?} {c3:?}", t0.elapsed());
    }
    Ok(regions)
}
