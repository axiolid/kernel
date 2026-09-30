//! Flat convex patches of the mesh measured TO.
//!
//! The piece bound `min over T of max over corners of d(v, T)` uses one
//! triangle at a time, so a piece over the seam between two triangles of a
//! flat face stays loose until it is split below the requested accuracy;
//! along every seam of a finely triangulated wall that is millions of
//! pieces. A patch is a set of triangles whose union is, up to a certified
//! excess, a convex set `K`: then `d(., union) <= d(., K) + excess`, the
//! distance to `K` is convex, and its maximum over a piece is at a corner,
//! where `d(v, K) <= min over members of d(v, T)`.
//!
//! # Why the hull is within `2 H` of the members
//!
//! Take a normal `n` and project along it onto the plane through a point of
//! the patch. Moving a point along `n` leaves every triple product
//! `n . (u x w)` unchanged, so the orientation tests below decide the
//! projected configuration exactly as they are evaluated on the unprojected
//! points. When they show the projected members' union is convex, a point
//! `y` of the hull of the members projects into the hull of the projected
//! corners, which is that union, so into some projected member; the same
//! barycentric point of the unprojected member is within `H` of the
//! projection, which is within `H` of `y`, where `H` bounds every corner's
//! height above the plane. Signs are only taken when they clear the triple
//! product's rounding (`16 eps |n| |u| |w|`, several times the standard-model
//! bound), and heights are rounded up.
//!
//! Two kinds are recognised, both common on the flat faces of building
//! models: two triangles sharing an edge that form a convex quadrilateral,
//! and the closed fan of triangles around an interior vertex whose link is a
//! convex polygon winding once.

use std::collections::HashMap;

use axiolid_core::Point3;

const EPS: f64 = f64::EPSILON;

/// Triangles whose union lies within `excess` of its convex hull's reach:
/// every point of the hull is within `excess` of a member.
pub(super) struct Patch {
    pub(super) members: Vec<usize>,
    pub(super) excess: f64,
}

type Bits = [u64; 3];

fn bits(p: Point3) -> Bits {
    p.to_array().map(f64::to_bits)
}

/// A normal with its length, for certified orientation signs and heights.
struct Frame {
    n: Point3,
    length: f64,
}

impl Frame {
    fn new(n: Point3) -> Option<Self> {
        let length = n.length();
        (length > 0.0 && length.is_finite()).then_some(Self { n, length })
    }

    /// The sign of `n . (u x w)`, or `None` within its rounding.
    fn sign(&self, u: Point3, w: Point3) -> Option<i8> {
        let value = self.n.dot(u.cross(w));
        let error = 16.0 * EPS * self.length * u.length() * w.length();
        if value > error {
            Some(1)
        } else if value < -error {
            Some(-1)
        } else {
            None
        }
    }

    /// An upper bound on the height of `p` above the plane through `origin`.
    fn height(&self, origin: Point3, p: Point3) -> f64 {
        let offset = p - origin;
        (self.n.dot(offset).abs() + 8.0 * EPS * self.length * offset.length())
            / (self.length * (1.0 - 4.0 * EPS))
    }

    /// The certified hull excess `2 H` for corners above `origin`.
    fn excess(&self, origin: Point3, corners: impl Iterator<Item = Point3>) -> f64 {
        let height = corners.map(|p| self.height(origin, p)).fold(0.0, f64::max);
        2.0 * height * (1.0 + 4.0 * EPS)
    }
}

/// Flat convex patches of `triangles` whose excess is at most `flatness`,
/// and for each triangle the patches it belongs to. Thin triangles are left
/// out. The order is deterministic.
pub(super) fn flat_patches(
    triangles: &[[Point3; 3]],
    thin: &[bool],
    flatness: f64,
) -> (Vec<Patch>, Vec<Vec<usize>>) {
    let mut edges: HashMap<(Bits, Bits), Vec<usize>> = HashMap::new();
    let mut stars: HashMap<Bits, Vec<usize>> = HashMap::new();
    for (index, corners) in triangles.iter().enumerate() {
        if thin[index] {
            continue;
        }
        for slot in 0..3 {
            let (p, q) = (bits(corners[slot]), bits(corners[(slot + 1) % 3]));
            edges.entry((p.min(q), p.max(q))).or_default().push(index);
            stars.entry(p).or_default().push(index);
        }
    }
    let mut edges: Vec<_> = edges.into_iter().collect();
    edges.sort_unstable_by_key(|(edge, _)| *edge);
    let mut stars: Vec<_> = stars.into_iter().collect();
    stars.sort_unstable_by_key(|(vertex, _)| *vertex);

    let mut patches = Vec::new();
    for ((p, q), members) in edges {
        for (at, &first) in members.iter().enumerate() {
            for &second in &members[at + 1..] {
                if let Some(excess) = pair_excess(triangles[first], triangles[second], p, q) {
                    patches.push(Patch {
                        members: vec![first, second],
                        excess,
                    });
                }
            }
        }
    }
    for (vertex, members) in stars {
        if members.len() < 3 {
            continue;
        }
        if let Some(excess) = fan_excess(triangles, &members, vertex) {
            patches.push(Patch { members, excess });
        }
    }
    patches.retain(|patch| patch.excess <= flatness);

    let mut of = vec![Vec::new(); triangles.len()];
    for (id, patch) in patches.iter().enumerate() {
        for &member in &patch.members {
            of[member].push(id);
        }
    }
    (patches, of)
}

/// The corner of `corners` that is neither `p` nor `q`.
fn apex(corners: [Point3; 3], p: Bits, q: Bits) -> Option<Point3> {
    let mut rest = corners
        .into_iter()
        .filter(|&c| bits(c) != p && bits(c) != q);
    let apex = rest.next()?;
    rest.next().is_none().then_some(apex)
}

/// `abc` and `abd` sharing `ab`: the excess when the projected quadrilateral
/// `a c b d` is convex, which holds exactly when `c`, `d` lie on opposite
/// sides of `ab` and `a`, `b` on opposite sides of `cd`.
fn pair_excess(first: [Point3; 3], second: [Point3; 3], p: Bits, q: Bits) -> Option<f64> {
    let a = first.into_iter().find(|&c| bits(c) == p)?;
    let b = first.into_iter().find(|&c| bits(c) == q)?;
    let c = apex(first, p, q)?;
    let d = apex(second, p, q)?;
    let frame = Frame::new((b - a).cross(c - a))?;
    let across_ab = frame.sign(b - a, c - a)? == 1 && frame.sign(b - a, d - a)? == -1;
    let across_cd = frame.sign(d - c, a - c)? * frame.sign(d - c, b - c)? == -1;
    (across_ab && across_cd).then(|| frame.excess(a, [b, c, d].into_iter()))
}

/// The closed fan of `members` around `vertex`: the excess when its link is
/// a convex polygon that winds once around the vertex, each fan triangle
/// turning the same way. The union is then the link polygon, convex.
fn fan_excess(triangles: &[[Point3; 3]], members: &[usize], vertex: Bits) -> Option<f64> {
    let v = triangles[members[0]]
        .into_iter()
        .find(|&c| bits(c) == vertex)?;
    // The link edge opposite the vertex in each member.
    let mut link: Vec<(Point3, Point3)> = Vec::with_capacity(members.len());
    for &member in members {
        let mut rest = triangles[member].into_iter().filter(|&c| bits(c) != vertex);
        let (l, r) = (rest.next()?, rest.next()?);
        if rest.next().is_some() || bits(l) == bits(r) {
            return None;
        }
        link.push((l, r));
    }
    // Chain the link edges into one cycle; every link vertex must be shared
    // by exactly two of them.
    let mut cycle = vec![link[0].0, link[0].1];
    let mut used = vec![false; link.len()];
    used[0] = true;
    for _ in 1..link.len() {
        let tail = bits(*cycle.last()?);
        let mut next = None;
        for (at, &(l, r)) in link.iter().enumerate() {
            if used[at] {
                continue;
            }
            let step = if bits(l) == tail {
                Some(r)
            } else if bits(r) == tail {
                Some(l)
            } else {
                None
            };
            if let Some(step) = step {
                if next.is_some() {
                    return None;
                }
                next = Some((at, step));
            }
        }
        let (at, step) = next?;
        used[at] = true;
        cycle.push(step);
    }
    // Closed: back at the start, and no link vertex visited twice.
    if bits(cycle.pop()?) != bits(cycle[0]) {
        return None;
    }
    let mut seen: Vec<Bits> = cycle.iter().map(|&p| bits(p)).collect();
    seen.sort_unstable();
    seen.dedup();
    if seen.len() != cycle.len() {
        return None;
    }

    let count = cycle.len();
    let normal = (0..count).fold(Point3::ZERO, |sum, i| {
        sum + (cycle[i] - v).cross(cycle[(i + 1) % count] - v)
    });
    let frame = Frame::new(normal)?;
    let unit = frame.n / frame.length;
    let mut turning = 0.0;
    for i in 0..count {
        let (prev, here, next) = (
            cycle[(i + count - 1) % count],
            cycle[i],
            cycle[(i + 1) % count],
        );
        // Each fan triangle turns the same way about the vertex, and the
        // link turns that way at every corner.
        if frame.sign(here - v, next - v)? != 1 || frame.sign(here - prev, next - here)? != 1 {
            return None;
        }
        let (a, b) = (here - v, next - v);
        turning += unit
            .dot(a.cross(b))
            .atan2(a.dot(b) - unit.dot(a) * unit.dot(b));
    }
    // The fan angles sum to 2 pi times the winding number; rounding moves
    // the sum by far less than the one radian allowed here.
    if (turning - core::f64::consts::TAU).abs() >= 1.0 {
        return None;
    }
    Some(frame.excess(v, cycle.into_iter()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z)
    }

    #[test]
    fn a_convex_quadrilateral_pairs_and_a_dart_does_not() {
        let (a, b, c, d) = (
            p(0.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
        );
        let flat = pair_excess([a, b, c], [a, b, d], bits(a), bits(b));
        assert!(flat.is_some_and(|excess| excess < 1e-13));
        // `d` on the same side of `ab` as `c`: the union is not convex.
        let same_side = pair_excess([a, b, c], [a, b, p(0.9, 0.2, 0.0)], bits(a), bits(b));
        assert!(same_side.is_none());
        // A dart: `a` reflex, so `a` and `b` lie on one side of `cd`.
        let dart = pair_excess([a, b, c], [a, b, p(-1.0, -0.2, 0.0)], bits(a), bits(b));
        assert!(dart.is_none());
        // Folded out of plane: accepted, with the fold's height as excess.
        let fold = pair_excess([a, b, c], [a, b, p(0.0, 1.0, 0.25)], bits(a), bits(b));
        assert!(fold.is_some_and(|excess| excess >= 0.5));
    }

    #[test]
    fn a_doubly_wound_fan_is_refused() {
        // A heptagram {7/2} link: every turn is convex and every fan
        // triangle positive, but it winds twice and leaves the caps out.
        let v = p(0.0, 0.0, 0.0);
        let ring: Vec<Point3> = (0..7)
            .map(|i| {
                let angle = core::f64::consts::TAU * f64::from(i) / 7.0;
                p(angle.cos(), angle.sin(), 0.0)
            })
            .collect();
        let triangles: Vec<[Point3; 3]> = (0..7)
            .map(|i| [v, ring[(2 * i) % 7], ring[(2 * i + 2) % 7]])
            .collect();
        let members: Vec<usize> = (0..7).collect();
        assert!(fan_excess(&triangles, &members, bits(v)).is_none());
        // Winding once, the same ring is a convex fan.
        let once: Vec<[Point3; 3]> = (0..7).map(|i| [v, ring[i], ring[(i + 1) % 7]]).collect();
        assert!(fan_excess(&once, &members, bits(v)).is_some());
    }

    #[test]
    fn a_fan_with_a_reflex_link_corner_is_refused() {
        // Six ring points, one pulled in towards the centre: the fan still
        // winds once, but its union is a dented hexagon.
        let v = p(0.0, 0.0, 0.0);
        let ring: Vec<Point3> = (0..6)
            .map(|i| {
                let angle = core::f64::consts::TAU * f64::from(i) / 6.0;
                let radius = if i == 0 { 0.3 } else { 1.0 };
                p(radius * angle.cos(), radius * angle.sin(), 0.0)
            })
            .collect();
        let triangles: Vec<[Point3; 3]> = (0..6).map(|i| [v, ring[i], ring[(i + 1) % 6]]).collect();
        let members: Vec<usize> = (0..6).collect();
        assert!(fan_excess(&triangles, &members, bits(v)).is_none());
    }
}
