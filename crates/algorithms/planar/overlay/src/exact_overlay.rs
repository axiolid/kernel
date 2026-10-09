//! Exact straight-edge boolean (#173).
//!
//! The polygon path used to hand its rings to an integer backend that maps
//! the operands' bounding box onto an `i32` grid and back, so every output
//! coordinate, even an input vertex the operation does not move, came back
//! snapped to a step of about 1.5e-8 of the extent. It now runs on the
//! exact subdivision the arc path uses (`exact_arc::arrangement`; a
//! straight ring is an arc ring without arcs):
//!
//! 1. Cut the plane by every ring of both operands at once. Where
//!    boundaries cross, which pieces coincide and which rings contain the
//!    region on each side of a piece are exact signs.
//! 2. Count, on each side of every piece, each operand's winding number:
//!    the rings containing that side, each counted `+1` when given
//!    counter-clockwise and `-1` when clockwise. The fill rule reads the
//!    count, the operation combines the two operands.
//! 3. Keep the pieces with the result on exactly one side, facing it, and
//!    link them into regions. A vertex where the boundary goes straight on
//!    (it only marked where another ring touched) is dropped, exactly.
//! 4. Round once: an input vertex is a double already and comes back
//!    bit-identical, and a crossing is the double nearest to it.
//!
//! Before step 1, each operand's rings are reduced as a chain: a winding
//! number is a sum over directed edges, so an edge given once each way
//! between the same two exact points cancels without changing any winding
//! number off it. A mesh given as a soup of its triangles shrinks to its
//! outline this way. What is left is split into cycles through distinct
//! vertices; if every one is simple, decided exactly, they replace the
//! operand's rings, else the rings are kept as given.
//!
//! The tolerance enters once, before step 1: features closer than it are
//! taken as touching and snapped together (see [`snap_within`], #222).
//! Past that it enters no decision. It validates operands, and
//! the result is settled like every other output ([`crate::settle`]), which
//! merges edges that rounding left shorter than the tolerance.

use std::collections::HashMap;

use axiolid_core::Point2;
use axiolid_guarantees::Sign;

use crate::arc::{ArcRing, ArcVertex};
use crate::exact_arc::{arrangement, orient_doubles};
use crate::{orientation, FillRule, OverlayError, OverlayOperation, Polygon, Ring};

/// Which operand a ring belongs to, and how it counts.
struct Counted {
    clip: bool,
    /// `+1` for a counter-clockwise ring as given, `-1` for clockwise.
    weight: i64,
}

fn filled(fill: FillRule, winding: i64) -> bool {
    match fill {
        FillRule::EvenOdd => winding % 2 != 0,
        FillRule::NonZero => winding != 0,
        FillRule::Positive => winding > 0,
        FillRule::Negative => winding < 0,
    }
}

fn combined(operation: OverlayOperation, subject: bool, clip: bool) -> bool {
    match operation {
        OverlayOperation::Intersection => subject && clip,
        OverlayOperation::Union => subject || clip,
        OverlayOperation::Difference => subject && !clip,
        OverlayOperation::Xor => subject != clip,
    }
}

/// The boolean of two sets of polygons, each ring counted by its own
/// orientation. Rings must be simple (validated); they may overlap each
/// other, within an operand too.
///
/// Returns the result's rings as `(outer, holes)`, before canonical
/// orientation and settling.
pub(crate) fn boolean(
    subject: &[Polygon],
    clip: &[Polygon],
    operation: OverlayOperation,
    fill: FillRule,
    snap: f64,
) -> Result<Vec<(Ring, Vec<Ring>)>, OverlayError> {
    boolean_reduced(subject, clip, operation, fill, snap, true)
}

/// [`boolean`], with the chain reduction on or off.
fn boolean_reduced(
    subject: &[Polygon],
    clip: &[Polygon],
    operation: OverlayOperation,
    fill: FillRule,
    snap: f64,
    reduced: bool,
) -> Result<Vec<(Ring, Vec<Ring>)>, OverlayError> {
    let mut all: Vec<(bool, Vec<Point2>)> = Vec::new();
    for (is_clip, polygons) in [(false, subject), (true, clip)] {
        for ring in polygons
            .iter()
            .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
        {
            all.push((is_clip, ring.points.clone()));
        }
    }
    snap_within(&mut all, snap);
    let mut rings: Vec<ArcRing> = Vec::new();
    let mut counted: Vec<Counted> = Vec::new();
    for is_clip in [false, true] {
        let given: Vec<Vec<Point2>> = all
            .iter()
            .filter(|(clip, _)| *clip == is_clip)
            .map(|(_, points)| points.clone())
            // A ring on one line -- the shadow of a vertical face, say --
            // encloses nothing and changes no winding number (#219).
            .filter(|points| points.len() >= 3 && !on_one_line(points))
            .collect();
        let rings_of = if reduced { reduce(given) } else { given };
        for mut points in rings_of {
            // Exact, so a sliver far from the origin keeps its winding
            // (#274).
            let positive = orientation(&points) == Sign::Positive;
            if !positive {
                points.reverse();
            }
            rings.push(ArcRing {
                vertices: points.into_iter().map(ArcVertex::straight).collect(),
            });
            counted.push(Counted {
                clip: is_clip,
                weight: if positive { 1 } else { -1 },
            });
        }
    }
    if rings.is_empty() {
        return Ok(Vec::new());
    }

    let raw = arrangement::build(&rings);
    // Whether the result holds the region on one side of a piece: the
    // rings containing that side are the ones containing the piece, and
    // each carrier whose inside faces that side.
    let result = |edge: &arrangement::RawEdge, left: bool| {
        let (mut subject, mut clip) = (0, 0);
        let rings = edge.inside.iter().copied().chain(
            edge.sources
                .iter()
                .filter(|source| source.2 == left)
                .map(|source| source.0),
        );
        for ring in rings {
            let count = &counted[ring];
            if count.clip {
                clip += count.weight;
            } else {
                subject += count.weight;
            }
        }
        combined(operation, filled(fill, subject), filled(fill, clip))
    };
    let keep: Vec<Option<bool>> = raw
        .edges
        .iter()
        .map(|edge| {
            let (left, right) = (result(edge, true), result(edge, false));
            // Kept pieces face the result: reversed when it lies right.
            (left != right).then_some(right)
        })
        .collect();

    // Only the vertices of kept pieces are rounded: rounding a crossing
    // exactly is not free, and most vertices of a soup are discarded.
    let positions: std::cell::RefCell<HashMap<usize, Point2>> = Default::default();
    let to_ring = |uses: &arrangement::RingUses| {
        let n = uses.len();
        let points = (0..n)
            .filter(|&i| !raw.straight_on(uses[(i + n - 1) % n], uses[i]))
            .map(|i| {
                let (piece, reversed) = uses[i];
                let edge = &raw.edges[piece];
                let vertex = if reversed { edge.to } else { edge.from };
                *positions
                    .borrow_mut()
                    .entry(vertex)
                    .or_insert_with(|| raw.vertex_position(vertex))
            })
            .collect();
        Ring { points }
    };
    Ok(raw
        .regions(&keep)?
        .iter()
        .map(|(outer, holes)| (to_ring(outer), holes.iter().map(to_ring).collect()))
        .collect())
}

/// Take features closer than `tolerance` as touching (#222): a vertex that
/// near an earlier one moves onto it, and a vertex that near another ring's
/// edge is inserted into it. Coordinates computed two ways -- a door's jamb
/// at `4.199999999999999` beside a wall face at `4.2` -- then meet exactly,
/// as the caller meant, where the exact boolean would keep a gap one ulp
/// wide. Inputs with nothing that near another feature come back
/// unchanged, bit for bit. A ring snapping would pinch onto itself keeps
/// its own vertices.
fn snap_within(rings: &mut [(bool, Vec<Point2>)], tolerance: f64) {
    if !(tolerance.is_finite() && tolerance > 0.0) {
        return;
    }
    let cell = |p: Point2| -> Option<(i64, i64)> {
        let (x, y) = ((p.x / tolerance).floor(), (p.y / tolerance).floor());
        (x.abs() < 1e17 && y.abs() < 1e17).then_some((x as i64, y as i64))
    };
    let original: Vec<Vec<Point2>> = rings.iter().map(|(_, r)| r.clone()).collect();
    // Vertices onto earlier ones, through a grid of `tolerance` cells.
    let mut grid: HashMap<(i64, i64), Vec<Point2>> = HashMap::new();
    let mut moved = false;
    for (_, ring) in rings.iter_mut() {
        for v in ring.iter_mut() {
            let Some((cx, cy)) = cell(*v) else {
                return;
            };
            let mut found = None;
            'search: for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(reps) = grid.get(&(cx + dx, cy + dy)) {
                        for r in reps {
                            if (*r - *v).length() <= tolerance {
                                found = Some(*r);
                                break 'search;
                            }
                        }
                    }
                }
            }
            match found {
                Some(r) => {
                    if identity(r) != identity(*v) {
                        moved = true;
                    }
                    *v = r;
                }
                None => grid.entry((cx, cy)).or_default().push(*v),
            }
        }
    }
    // Vertices onto the edges of other rings, found by x range.
    let mut reps: Vec<(Point2, usize)> = Vec::new();
    for (index, (_, ring)) in rings.iter().enumerate() {
        reps.extend(ring.iter().map(|p| (*p, index)));
    }
    reps.sort_by(|a, b| a.0.x.total_cmp(&b.0.x).then(a.0.y.total_cmp(&b.0.y)));
    reps.dedup_by(|a, b| identity(a.0) == identity(b.0) && a.1 == b.1);
    for (index, entry) in rings.iter_mut().enumerate() {
        let ring = entry.1.clone();
        let n = ring.len();
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let (p, q) = (ring[i], ring[(i + 1) % n]);
            out.push(p);
            let d = q - p;
            let length2 = d.dot(d);
            if length2.is_nan() || length2 <= 0.0 {
                continue;
            }
            let (lo, hi) = (p.x.min(q.x) - tolerance, p.x.max(q.x) + tolerance);
            let first = reps.partition_point(|r| r.0.x < lo);
            let mut inside: Vec<(f64, Point2)> = Vec::new();
            for &(r, owner) in reps[first..].iter().take_while(|r| r.0.x <= hi) {
                if owner == index || identity(r) == identity(p) || identity(r) == identity(q) {
                    continue;
                }
                let t = (r - p).dot(d) / length2;
                if !(t > 0.0 && t < 1.0) {
                    continue;
                }
                let foot = p + d * t;
                if (r - foot).length() <= tolerance
                    && !inside.iter().any(|(_, s)| identity(*s) == identity(r))
                {
                    inside.push((t, r));
                }
            }
            if !inside.is_empty() {
                moved = true;
                inside.sort_by(|a, b| a.0.total_cmp(&b.0));
                out.extend(inside.into_iter().map(|(_, r)| r));
            }
        }
        entry.1 = out;
    }
    if !moved {
        return;
    }
    for (index, (_, ring)) in rings.iter_mut().enumerate() {
        // Consecutive duplicates go; a spike back along an edge goes.
        ring.dedup_by(|a, b| identity(*a) == identity(*b));
        while ring.len() > 1 && identity(ring[0]) == identity(ring[ring.len() - 1]) {
            ring.pop();
        }
        let mut changed = true;
        while changed && ring.len() >= 3 {
            changed = false;
            let n = ring.len();
            for i in 0..n {
                if identity(ring[i]) == identity(ring[(i + 2) % n]) {
                    let drop = [(i + 1) % n, (i + 2) % n];
                    let (a, b) = (drop[0].max(drop[1]), drop[0].min(drop[1]));
                    ring.remove(a);
                    ring.remove(b);
                    changed = true;
                    break;
                }
            }
        }
        // A pinched ring keeps its own vertices.
        let mut seen = std::collections::HashSet::new();
        if !ring.iter().all(|p| seen.insert(identity(*p))) {
            *ring = original[index].clone();
        }
    }
}

/// Whether every vertex of a ring lies on one line, exactly: then the
/// ring encloses no area, whatever its rounded area says, and a boolean
/// leaves it out. The exact subdivision assumes simple rings, and a ring
/// running back over itself along a line is not one (#219).
pub(crate) fn on_one_line(points: &[Point2]) -> bool {
    let Some(&a) = points.first() else {
        return true;
    };
    let Some(&b) = points.iter().find(|p| **p != a) else {
        return true;
    };
    points
        .iter()
        .all(|&p| orient_doubles(a, b, p) == Sign::Zero)
}

/// An exact point's identity: its coordinates' bits, `-0.0` read as `0.0`.
fn identity(p: Point2) -> (u64, u64) {
    ((p.x + 0.0).to_bits(), (p.y + 0.0).to_bits())
}

/// One operand's rings, as a chain with every edge given both ways
/// cancelled, split into cycles through distinct vertices -- or the rings
/// as given, when nothing cancels or a cycle is not simple.
fn reduce(rings: Vec<Vec<Point2>>) -> Vec<Vec<Point2>> {
    let mut ids: HashMap<(u64, u64), usize> = HashMap::new();
    let mut points: Vec<Point2> = Vec::new();
    // Net multiplicity of each edge, `+` from the lower vertex id.
    let mut net: HashMap<(usize, usize), i64> = HashMap::new();
    let mut total = 0usize;
    for ring in &rings {
        let mut id = |p: Point2| {
            *ids.entry(identity(p)).or_insert_with(|| {
                points.push(p);
                points.len() - 1
            })
        };
        let n = ring.len();
        let vertices: Vec<usize> = ring.iter().map(|p| id(*p)).collect();
        for i in 0..n {
            let (a, b) = (vertices[i], vertices[(i + 1) % n]);
            if a == b {
                continue;
            }
            total += 1;
            let (edge, step) = if a < b { ((a, b), 1) } else { ((b, a), -1) };
            *net.entry(edge).or_default() += step;
        }
    }
    let left: i64 = net.values().map(|m| m.abs()).sum();
    if usize::try_from(left).is_ok_and(|left| left == total) {
        return rings;
    }
    let mut edges: Vec<((usize, usize), i64)> = net.into_iter().filter(|(_, m)| *m != 0).collect();
    edges.sort_unstable();
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); points.len()];
    for ((lo, hi), m) in edges {
        let (from, to) = if m > 0 { (lo, hi) } else { (hi, lo) };
        for _ in 0..m.unsigned_abs() {
            out[from].push(to);
        }
    }
    // Every vertex has as many edges in as out, so a walk only stops where
    // it started. A walk returning to a vertex it passed closes a cycle
    // there, which is cut off: each cycle meets each vertex once.
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    let mut position = vec![usize::MAX; points.len()];
    for start in 0..points.len() {
        while !out[start].is_empty() {
            let mut path = vec![start];
            position[start] = 0;
            let mut at = start;
            loop {
                let Some(next) = out[at].pop() else {
                    return rings;
                };
                if position[next] == usize::MAX {
                    position[next] = path.len();
                    path.push(next);
                    at = next;
                    continue;
                }
                let from = position[next];
                let cycle: Vec<usize> = path.drain(from..).collect();
                for &v in &cycle {
                    position[v] = usize::MAX;
                }
                cycles.push(cycle);
                if path.is_empty() {
                    break;
                }
                position[next] = path.len();
                path.push(next);
                at = next;
            }
        }
    }
    let cycles: Vec<Vec<Point2>> = cycles
        .into_iter()
        .map(|cycle| cycle.into_iter().map(|v| points[v]).collect())
        .collect();
    if cycles.iter().all(|cycle| simple(cycle)) {
        cycles
    } else {
        rings
    }
}

/// Whether the closed segments `a-b` and `c-d` share a point, exactly.
fn segments_meet(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let within = |p: Point2, q: Point2, r: Point2| {
        r.x >= p.x.min(q.x) && r.x <= p.x.max(q.x) && r.y >= p.y.min(q.y) && r.y <= p.y.max(q.y)
    };
    let (o1, o2) = (orient_doubles(a, b, c), orient_doubles(a, b, d));
    let (o3, o4) = (orient_doubles(c, d, a), orient_doubles(c, d, b));
    let opposite = |p: Sign, q: Sign| {
        matches!(
            (p, q),
            (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
        )
    };
    if opposite(o1, o2) && opposite(o3, o4) {
        return true;
    }
    (o1 == Sign::Zero && within(a, b, c))
        || (o2 == Sign::Zero && within(a, b, d))
        || (o3 == Sign::Zero && within(c, d, a))
        || (o4 == Sign::Zero && within(c, d, b))
}

/// Whether a ring through distinct vertices is simple, exactly: edges
/// apart from their neighbours, and each pair of neighbours sharing only
/// their common vertex.
fn simple(ring: &[Point2]) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let edge = |i: usize| (ring[i], ring[(i + 1) % n]);
    let mut order: Vec<usize> = (0..n).collect();
    let low = |i: usize| edge(i).0.x.min(edge(i).1.x);
    order.sort_by(|&i, &j| low(i).total_cmp(&low(j)));
    let mut active: Vec<usize> = Vec::new();
    for &i in &order {
        let (a, b) = edge(i);
        active.retain(|&j| {
            let (c, d) = edge(j);
            c.x.max(d.x) >= low(i)
        });
        for &j in &active {
            let (c, d) = edge(j);
            if a.y.max(b.y) < c.y.min(d.y) || c.y.max(d.y) < a.y.min(b.y) {
                continue;
            }
            let (first, second) = (i.min(j), i.max(j));
            let neighbours = second == first + 1 || (first == 0 && second == n - 1);
            if neighbours {
                // The shared vertex, and the far ends: on one line and on
                // one side of it, the edges overlap.
                let (v, u, w) = if second == first + 1 {
                    (ring[second], ring[first], ring[(second + 1) % n])
                } else {
                    (ring[0], ring[1], ring[n - 1])
                };
                if orient_doubles(u, v, w) == Sign::Zero {
                    let ahead = |p: Point2| ((p.x - v.x).signum(), (p.y - v.y).signum());
                    let (du, dw) = (ahead(u), ahead(w));
                    if (du.0 == dw.0 && u.x != v.x) || (du.1 == dw.1 && u.y != v.y) {
                        return false;
                    }
                }
            } else if segments_meet(a, b, c, d) {
                return false;
            }
        }
        active.push(i);
    }
    true
}

#[cfg(test)]
mod tests {
    //! The exact boolean against the integer backend it replaced: the same
    //! regions, up to that backend's grid.
    use super::*;
    use axiolid_core::Point2;
    use i_overlay::core::{fill_rule::FillRule as Fill, overlay_rule::OverlayRule};
    use i_overlay::float::single::SingleFloatOverlay;

    /// A small deterministic generator, so failures reproduce.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// A star-shaped ring around `(cx, cy)`: simple by construction.
    fn star(rng: &mut Lcg, cx: f64, cy: f64, n: usize, clockwise: bool) -> Ring {
        let mut points: Vec<Point2> = (0..n)
            .map(|i| {
                let angle = std::f64::consts::TAU * (i as f64 + 0.8 * rng.next()) / n as f64;
                let radius = 0.4 + rng.next();
                Point2::new(cx + radius * angle.cos(), cy + radius * angle.sin())
            })
            .collect();
        if clockwise {
            points.reverse();
        }
        Ring { points }
    }

    /// Even-odd membership over all rings of a result.
    fn inside(polygons: &[(Ring, Vec<Ring>)], p: Point2) -> bool {
        let crosses = |ring: &Ring| {
            let n = ring.points.len();
            (0..n)
                .filter(|&i| {
                    let (a, b) = (ring.points[i], ring.points[(i + 1) % n]);
                    (a.y > p.y) != (b.y > p.y)
                        && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x
                })
                .count()
                % 2
                == 1
        };
        polygons
            .iter()
            .flat_map(|(outer, holes)| std::iter::once(outer).chain(holes))
            .filter(|ring| crosses(ring))
            .count()
            % 2
            == 1
    }

    fn backend(polygons: &[Polygon]) -> Vec<Vec<Vec<[f64; 2]>>> {
        polygons
            .iter()
            .map(|p| {
                std::iter::once(&p.outer)
                    .chain(&p.holes)
                    .map(|r| r.points.iter().map(|q| [q.x, q.y]).collect())
                    .collect()
            })
            .collect()
    }

    fn distance_to_boundary(polygons: &[Polygon], p: Point2) -> f64 {
        let mut best = f64::INFINITY;
        for ring in polygons
            .iter()
            .flat_map(|q| std::iter::once(&q.outer).chain(&q.holes))
        {
            let n = ring.points.len();
            for i in 0..n {
                let (a, b) = (ring.points[i], ring.points[(i + 1) % n]);
                let t = ((p - a).dot(b - a) / (b - a).length_squared()).clamp(0.0, 1.0);
                best = best.min((a + (b - a) * t - p).length());
            }
        }
        best
    }

    /// A jittered grid of `n x n` cells from `(x0, y0)`, two triangles
    /// each, clockwise when asked: a mesh as a soup of its triangles.
    fn mesh(rng: &mut Lcg, n: usize, x0: f64, y0: f64, clockwise: bool) -> Vec<Polygon> {
        let at: Vec<Vec<Point2>> = (0..=n)
            .map(|i| {
                (0..=n)
                    .map(|j| {
                        let jitter = if i % n == 0 || j % n == 0 { 0.0 } else { 0.04 };
                        Point2::new(
                            x0 + i as f64 * 0.25 + jitter * (rng.next() - 0.5),
                            y0 + j as f64 * 0.25 + jitter * (rng.next() - 0.5),
                        )
                    })
                    .collect()
            })
            .collect();
        let mut out = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let (a, b, c, d) = (at[i][j], at[i + 1][j], at[i + 1][j + 1], at[i][j + 1]);
                for mut points in [vec![a, b, c], vec![a, c, d]] {
                    if clockwise {
                        points.reverse();
                    }
                    out.push(Polygon {
                        outer: Ring { points },
                        holes: Vec::new(),
                    });
                }
            }
        }
        out
    }

    /// Cancelling a soup's shared edges changes no result: every operation
    /// and fill rule, on meshes (overlapping, touching at a corner, with a
    /// cell missing, clockwise), triangles, and polygons with holes, gives
    /// the same polygons, bit for bit, with the reduction as without.
    #[test]
    fn the_reduction_changes_nothing() {
        let mut rng = Lcg(198);
        let operations = [
            OverlayOperation::Intersection,
            OverlayOperation::Union,
            OverlayOperation::Difference,
            OverlayOperation::Xor,
        ];
        let fills = [
            FillRule::EvenOdd,
            FillRule::NonZero,
            FillRule::Positive,
            FillRule::Negative,
        ];
        let mut reduced_somewhere = false;
        for case in 0..64 {
            let mut subject = mesh(&mut rng, 3 + case % 4, 0.0, 0.0, case % 5 == 4);
            if case % 3 == 0 {
                // A cell missing: the outline has a hole.
                subject.remove(2 * (1 + case % 3) + 1);
                subject.remove(2 * (1 + case % 3));
            }
            if case % 4 == 1 {
                // A second mesh touching the first at a corner only.
                subject.extend(mesh(&mut rng, 2, -0.5, -0.5, false));
            }
            let mut clip = if case % 2 == 0 {
                let x0 = 0.3 + 0.1 * rng.next();
                mesh(&mut rng, 3, x0, 0.2, case % 7 == 6)
            } else {
                (0..4)
                    .map(|_| {
                        let (cx, cy, cw) = (rng.next(), rng.next(), rng.next() < 0.3);
                        Polygon {
                            outer: star(&mut rng, cx, cy, 3, cw),
                            holes: Vec::new(),
                        }
                    })
                    .collect()
            };
            if case % 6 == 5 {
                clip.push(Polygon {
                    outer: Ring {
                        points: vec![
                            Point2::new(-1.0, -1.0),
                            Point2::new(2.0, -1.0),
                            Point2::new(2.0, 2.0),
                            Point2::new(-1.0, 2.0),
                        ],
                    },
                    holes: vec![Ring {
                        points: vec![
                            Point2::new(0.25, 0.25),
                            Point2::new(0.25, 0.5),
                            Point2::new(0.5, 0.5),
                            Point2::new(0.5, 0.25),
                        ],
                    }],
                });
            }
            let given = |polygons: &[Polygon]| -> Vec<Vec<Point2>> {
                polygons
                    .iter()
                    .flat_map(|p| std::iter::once(&p.outer).chain(&p.holes))
                    .map(|r| r.points.clone())
                    .collect()
            };
            reduced_somewhere |= reduce(given(&subject)).len() < given(&subject).len();
            let operation = operations[case % 4];
            let fill = fills[(case / 4) % 4];
            let with = boolean_reduced(&subject, &clip, operation, fill, 0.0, true).unwrap();
            let without = boolean_reduced(&subject, &clip, operation, fill, 0.0, false).unwrap();
            let canonical = |rings: Vec<(Ring, Vec<Ring>)>| {
                crate::canonical_polygons(rings)
                    .into_iter()
                    .map(|p| {
                        std::iter::once(p.outer)
                            .chain(p.holes)
                            .map(|r| {
                                r.points
                                    .iter()
                                    .map(|q| (q.x.to_bits(), q.y.to_bits()))
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            };
            assert_eq!(
                canonical(with),
                canonical(without),
                "case {case}: {operation:?} {fill:?}"
            );
        }
        assert!(reduced_somewhere);
    }

    /// A mesh reduces to its outline: one ring, and one more per hole.
    #[test]
    fn a_mesh_reduces_to_its_outline() {
        let mut rng = Lcg(1);
        let mut soup = mesh(&mut rng, 6, 0.0, 0.0, false);
        let rings = |polygons: &[Polygon]| -> Vec<Vec<Point2>> {
            polygons.iter().map(|p| p.outer.points.clone()).collect()
        };
        let outline = reduce(rings(&soup));
        assert_eq!(outline.len(), 1);
        assert_eq!(outline[0].len(), 24);
        // A cell missing inside: its two triangles out, a hole in.
        let k = 2 * (6 * 2 + 2);
        soup.drain(k..k + 2);
        assert_eq!(reduce(rings(&soup)).len(), 2);
        // Crossing triangles cancel nothing: kept as given.
        let crossing = vec![
            vec![
                Point2::new(0.0, 0.0),
                Point2::new(2.0, 0.0),
                Point2::new(1.0, 2.0),
            ],
            vec![
                Point2::new(0.0, 1.0),
                Point2::new(2.0, 1.0),
                Point2::new(1.0, -1.0),
            ],
        ];
        assert_eq!(reduce(crossing.clone()), crossing);
    }

    /// A ring whose vertex touches another of its edges is not simple,
    /// however it is started or run; one that only comes near is.
    #[test]
    fn a_ring_touching_itself_is_not_simple() {
        let p = Point2::new;
        let touching = [
            p(0.0, 0.0),
            p(4.0, 0.0),
            p(4.0, 2.0),
            p(2.0, 0.0),
            p(0.0, 2.0),
        ];
        let near = [
            p(0.0, 0.0),
            p(4.0, 0.0),
            p(4.0, 2.0),
            p(2.0, 0.5),
            p(0.0, 2.0),
        ];
        for shift in 0..5 {
            for reversed in [false, true] {
                let turn = |ring: &[Point2]| {
                    let mut out: Vec<Point2> = (0..5).map(|i| ring[(i + shift) % 5]).collect();
                    if reversed {
                        out.reverse();
                    }
                    out
                };
                assert!(!simple(&turn(&touching)), "{shift} {reversed}");
                assert!(simple(&turn(&near)), "{shift} {reversed}");
            }
        }
        // Neighbours folding back along one line: only the neighbour
        // check sees it, every edge of a triangle being a neighbour.
        assert!(!simple(&[p(0.0, 0.0), p(2.0, 0.0), p(1.0, 0.0)]));
        assert!(!simple(&[p(0.0, 0.0), p(0.0, 2.0), p(0.0, 1.0)]));
    }

    /// Two simple rings sharing an edge each way, the second reaching back
    /// over the first: the cycle left crosses itself, so the reduction
    /// keeps the rings as given -- and the result is the same either way.
    #[test]
    fn a_crossing_cycle_keeps_the_rings() {
        let p = Point2::new;
        let square = vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
        let hook = vec![
            p(1.0, 1.0),
            p(1.0, 0.0),
            p(2.0, 0.0),
            p(2.0, 2.0),
            p(0.5, 2.0),
            p(0.5, 0.5),
            p(0.7, 0.5),
            p(0.7, 1.5),
        ];
        let polygon = |points: &Vec<Point2>| Polygon {
            outer: Ring {
                points: points.clone(),
            },
            holes: Vec::new(),
        };
        let given = vec![square.clone(), hook.clone()];
        assert_eq!(reduce(given.clone()), given);
        let soup = [polygon(&square), polygon(&hook)];
        let with = boolean_reduced(
            &soup,
            &[],
            OverlayOperation::Union,
            FillRule::NonZero,
            0.0,
            true,
        );
        let without = boolean_reduced(
            &soup,
            &[],
            OverlayOperation::Union,
            FillRule::NonZero,
            0.0,
            false,
        );
        assert_eq!(format!("{with:?}"), format!("{without:?}"));
    }

    #[test]
    fn agrees_with_the_integer_backend() {
        let mut rng = Lcg(173);
        let operations = [
            (OverlayOperation::Intersection, OverlayRule::Intersect),
            (OverlayOperation::Union, OverlayRule::Union),
            (OverlayOperation::Difference, OverlayRule::Difference),
            (OverlayOperation::Xor, OverlayRule::Xor),
        ];
        let fills = [
            (FillRule::EvenOdd, Fill::EvenOdd),
            (FillRule::NonZero, Fill::NonZero),
            (FillRule::Positive, Fill::Positive),
            (FillRule::Negative, Fill::Negative),
        ];
        let mut compared = 0;
        for case in 0..256 {
            let mut operand = |count: usize| -> Vec<Polygon> {
                (0..count)
                    .map(|_| {
                        let (cx, cy) = (2.0 * rng.next(), 2.0 * rng.next());
                        // Every eighth case has long rings, which the
                        // subdivision indexes rather than scans.
                        let n = if case % 8 == 7 {
                            40 + (rng.next() * 40.0) as usize
                        } else {
                            3 + (rng.next() * 9.0) as usize
                        };
                        let clockwise = rng.next() < 0.3;
                        Polygon {
                            outer: star(&mut rng, cx, cy, n, clockwise),
                            holes: Vec::new(),
                        }
                    })
                    .collect()
            };
            let subject = operand(1 + case % 3);
            let clip = operand(1 + case % 2);
            let (operation, rule) = operations[case % 4];
            let (fill, backend_fill) = fills[(case / 4) % 4];
            let exact = boolean(&subject, &clip, operation, fill, 0.0).expect("simple rings link");
            let grid: Vec<(Ring, Vec<Ring>)> = backend(&subject)
                .overlay(&backend(&clip), rule, backend_fill)
                .into_iter()
                .map(|shape| {
                    let mut rings = shape.into_iter().map(|r| Ring {
                        points: r.into_iter().map(|p| Point2::new(p[0], p[1])).collect(),
                    });
                    let outer = rings.next().expect("a shape has an outer ring");
                    (outer, rings.collect())
                })
                .collect();
            let all: Vec<Polygon> = subject.iter().chain(&clip).cloned().collect();
            for _ in 0..200 {
                let p = Point2::new(4.0 * rng.next() - 1.0, 4.0 * rng.next() - 1.0);
                // Away from every input boundary, the grid cannot matter.
                if distance_to_boundary(&all, p) < 1e-6 {
                    continue;
                }
                compared += 1;
                assert_eq!(
                    inside(&exact, p),
                    inside(&grid, p),
                    "case {case}: {operation:?} {fill:?} at {p:?}"
                );
            }
        }
        assert!(compared > 40_000, "{compared}");
    }
}
