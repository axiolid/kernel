//! Lofting a sequence of station rings into a solid.
//!
//! Every sweep family reduces to the same shape: place the profile at a
//! series of stations along a path, stitch consecutive stations into
//! walls, and cap the ends unless the path closes on itself. Extrusion is
//! two stations; revolution is an arc of them; a sectioned spine supplies
//! its own.
//!
//! Writing that once means winding, hole orientation and cap pairing are
//! fixed in one place rather than re-derived per family.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_mesh::TriMesh;

use crate::profile::Rings;

/// One station: the profile's rings already placed in 3D.
pub struct Station {
    /// Outer ring followed by each hole, in `Rings` order.
    pub loops: Vec<Vec<Point3>>,
}

/// Place `rings` into 3D with a per-point mapping.
///
/// The mapping takes the profile's 2D point and returns its 3D image, so a
/// caller expresses only what its family does: extrusion translates,
/// revolution rotates, a spine sweep applies a frame.
pub fn place(rings: &Rings, mut f: impl FnMut(Point2) -> Point3) -> Station {
    let mut loops = Vec::with_capacity(1 + rings.holes.len());
    loops.push(rings.outer.iter().map(|p| f(*p)).collect());
    for hole in &rings.holes {
        loops.push(hole.iter().map(|p| f(*p)).collect());
    }
    Station { loops }
}

/// Stitch stations into a closed solid.
///
/// `closed` wraps the last station onto the first, which welds a periodic
/// seam by index rather than by two samplings agreeing numerically. Open
/// lofts are capped with the profile triangulation at each end.
///
/// All stations must share the profile's ring structure: the walls pair
/// vertices by position, so a differing ring length has no meaningful
/// pairing and is refused rather than silently truncated.
pub fn loft(rings: &Rings, stations: &[Station], closed: bool) -> GeomResult<TriMesh> {
    if stations.len() < 2 {
        return Err(GeomError::InvalidInput(format!(
            "a loft needs at least two stations, got {}",
            stations.len()
        )));
    }
    let shape: Vec<usize> = stations[0].loops.iter().map(|r| r.len()).collect();
    for s in stations {
        let this: Vec<usize> = s.loops.iter().map(|r| r.len()).collect();
        if this != shape {
            return Err(GeomError::InvalidInput(
                "loft stations must share the profile's ring structure".to_owned(),
            ));
        }
    }
    let per_station: usize = shape.iter().sum();
    let mut positions: Vec<Point3> = Vec::with_capacity(stations.len() * per_station);
    for s in stations {
        for ring in &s.loops {
            positions.extend(ring.iter().copied());
        }
    }
    // A closed loft wraps onto station 0; an open one stops one short.
    let spans = if closed {
        stations.len()
    } else {
        stations.len() - 1
    };
    let mut indices: Vec<u32> = Vec::new();
    for s in 0..spans {
        let a0 = s * per_station;
        let b0 = ((s + 1) % stations.len()) * per_station;
        let mut base = 0usize;
        for m in shape.iter() {
            for k in 0..*m {
                let kn = (k + 1) % *m;
                let (a, b) = ((a0 + base + k) as u32, (a0 + base + kn) as u32);
                let (c, d) = ((b0 + base + k) as u32, (b0 + base + kn) as u32);
                // Both rings wind the same way here. A hole ring is already
                // stored clockwise (profile.rs keeps holes reversed), which
                // is what turns its wall normal inward; flipping the index
                // order as well would invert it a second time.
                indices.extend([a, b, d, a, d, c]);
            }
            base += *m;
        }
    }
    // Caps. A closed loft needs none: the wall meets itself. An open one is
    // bounded by the profile at each end, wound opposite so both face out.
    if !closed {
        let (points, tris) = crate::profile::triangulate(rings)?;
        if points.len() != per_station {
            return Err(GeomError::Degenerate(format!(
                "cap triangulation produced {} points for {per_station} ring points",
                points.len()
            )));
        }
        let last = ((stations.len() - 1) * per_station) as u32;
        for t in &tris {
            indices.extend([t[0], t[2], t[1]]);
            indices.extend([last + t[0], last + t[1], last + t[2]]);
        }
    }
    Ok(TriMesh::new(positions, indices))
}

/// How far the bilinear patch over one wall quad can sit from the two
/// triangles [`loft`] splits it into.
///
/// The quad joins ring points `a -> b` on one station to `c -> d` on the
/// next, and [`loft`] cuts it along the diagonal `a-d` into `(a, b, d)` and
/// `(a, d, c)`. The patch `B(x, y) = (1-y)((1-x)a + xb) + y((1-x)c + xd)`
/// is the surface the chords between the two stations span: its row at
/// fixed `x` is the chord between the two placements of one profile point,
/// which is why the sweeps bound the exact surface against it. Every bound
/// below holds pointwise, so the smallest is returned:
///
/// - `|T| / 4` with the twist `T = a - b - c + d`. On the half `x >= y`,
///   `B - L = -y(1-x) T` where `L` is the point of `(a, b, d)` with the
///   same barycentric weights, and `y(1-x) <= 1/4` there; the other half is
///   symmetric. Always valid.
/// - `|b - a|`: `B` is within `(1-y)x|b - a|` of the point
///   `(1-y)a + y((1-x)c + xd)` of triangle `(a, d, c)`. Likewise `|c - a|`
///   against `(a, b, d)`. These make a quad that collapses to a triangle on
///   an axis (a ring point on a revolution axis) exact.
/// - `2|d - c|`: `B` is within `xy|d - c|` of the point
///   `(1-y)((1-x)a + xb) + yc` of the plane triangle `(a, b, c)`, and moving
///   that triangle's `c` to `d` moves the point by at most `|d - c|`.
///   Likewise `2|d - b|`.
/// - The flattening bound. Take a plane with unit normal `n` and heights
///   `h(p) = n.p - k`, project the corners onto it, and suppose the
///   projected quad `a', b', d', c'` is convex. The projected patch is a
///   convex combination of the projected corners, so it lies in that quad,
///   which the two plane triangles `(a', b', d')`, `(a', d', c')` cover:
///   each `B(x, y)` has a point `Z` there with `Z = B(x, y)` in the plane.
///   Lifting `Z` back with the same weights gives a point of `(a, b, d)` or
///   `(a, d, c)`, and the two heights are weighted means of corner heights,
///   so the distance is at most the spread of the corner heights. Two
///   planes are used. The plane through `a, b, c` leaves only `d` off it,
///   at the height `g`: bound `|g|`. The plane parallel to both diagonals,
///   halfway between them, puts `a, d` at `e` and `b, c` at `-e`, with
///   `2|e|` the distance between the diagonals: bound `2|e|`, half of
///   `|g|` for a twisted quad (where `|g| = 4|e|`). A planar convex quad,
///   the wall a revolution builds about an axis in the profile's plane,
///   gets 0.
///
/// The convexity test runs in floating point with a relative slack of
/// `1e-12`: a quad that fails it by less leaves an uncovered sliver some
/// `1e-12` of its size wide, far below any chord budget.
pub(crate) fn quad_deviation(a: Point3, b: Point3, c: Point3, d: Point3) -> Scalar {
    let twist = (a - b - c + d).length() / 4.0;
    let mut best = twist
        .min((b - a).length())
        .min((c - a).length())
        .min(2.0 * (d - c).length())
        .min(2.0 * (d - b).length());
    // The plane through a, b, c: only d leaves it.
    if let Some(n) = unit((b - a).cross(c - a)) {
        let g = (d - a).dot(n);
        if convex_projection([a, b, d - n * g, c], n) {
            best = best.min(g.abs());
        }
    }
    // The plane halfway between the diagonals a-d and b-c.
    if let Some(n) = unit((d - a).cross(c - b)) {
        let e = 0.5 * (a - b).dot(n);
        let (a2, d2) = (a - n * e, d - n * e);
        let (b2, c2) = (b + n * e, c + n * e);
        if convex_projection([a2, b2, d2, c2], n) {
            best = best.min(2.0 * e.abs());
        }
    }
    best
}

/// `v` scaled to unit length, or `None` when it has none to scale.
fn unit(v: Vec3) -> Option<Vec3> {
    let length = v.length();
    (length > 0.0 && length.is_finite()).then(|| v / length)
}

/// Whether the polygon `p0 -> p1 -> p2 -> p3`, lying in a plane with unit
/// normal `n`, is convex: it turns the same way, about `n`, at every
/// corner (either way round), with a relative slack of `1e-12`.
fn convex_projection(corners: [Point3; 4], n: Vec3) -> bool {
    let scale = (0..4)
        .map(|k| (corners[(k + 1) % 4] - corners[k]).length())
        .fold(0.0, Scalar::max);
    let slack = 1e-12 * scale * scale;
    let turns: Vec<Scalar> = (0..4)
        .map(|k| {
            let (p, q, r) = (corners[k], corners[(k + 1) % 4], corners[(k + 2) % 4]);
            (q - p).cross(r - q).dot(n)
        })
        .collect();
    turns.iter().all(|t| *t >= -slack) || turns.iter().all(|t| *t <= slack)
}

/// The largest [`quad_deviation`] over the walls between two stations.
///
/// Pairs ring points exactly as [`loft`] does: point `k` with `k + 1`,
/// wrapping, on every ring.
pub(crate) fn span_deviation(from: &Station, to: &Station) -> Scalar {
    let mut worst: Scalar = 0.0;
    for (ring, next) in from.loops.iter().zip(&to.loops) {
        let m = ring.len().min(next.len());
        for k in 0..m {
            let kn = (k + 1) % m;
            worst = worst.max(quad_deviation(ring[k], ring[kn], next[k], next[kn]));
        }
    }
    worst
}

/// An orthonormal frame at a point on a path.
///
/// A sweep needs a full frame, not just a tangent: the profile's x and y
/// axes have to be carried along the path, and how they are carried is
/// exactly what distinguishes the sweep families from each other.
pub struct Frame {
    /// Frame origin on the directrix.
    pub origin: Point3,
    /// Image of the profile's +x.
    pub x: Vec3,
    /// Image of the profile's +y.
    pub y: Vec3,
}

impl Frame {
    /// Frame whose z is `tangent` and whose x is `reference` made
    /// perpendicular to it.
    ///
    /// Refuses a reference parallel to the tangent instead of silently
    /// picking a fallback axis: the caller supplied a direction that cannot
    /// orient the profile, and quietly substituting one rotates the section
    /// by an arbitrary angle.
    pub fn from_reference(origin: Point3, tangent: Vec3, reference: Vec3) -> GeomResult<Self> {
        let t = tangent.normalize_or_zero();
        if t == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "sweep tangent must be a non-zero direction".to_owned(),
            ));
        }
        let x = (reference - t * t.dot(reference)).normalize_or_zero();
        if x == Vec3::ZERO {
            return Err(GeomError::InvalidInput(
                "sweep reference direction must not be parallel to the directrix".to_owned(),
            ));
        }
        Ok(Self {
            origin,
            x,
            y: t.cross(x),
        })
    }
}

/// Place a profile point using a frame.
pub fn at(frame: &Frame, p: Point2) -> Point3 {
    frame.origin + frame.x * p.x + frame.y * p.y
}

/// Linear blend of two profile rings.
///
/// Tapered families interpolate between a start and an end profile, so the
/// blend belongs here rather than in each of them.
pub fn blend(a: Point2, b: Point2, t: Scalar) -> Point2 {
    Point2::new(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t)
}

/// Loft where the two ends carry different profiles.
///
/// [`loft`] caps both ends from one ring set, which is wrong the moment the
/// ends differ: a tapered solid needs each cap triangulated from its own
/// profile. The walls are identical, so only the caps are rebuilt here.
pub fn loft_tapered(start: &Rings, end: &Rings, stations: &[Station]) -> GeomResult<TriMesh> {
    let mut mesh = loft(start, stations, false)?;
    let per_station: usize = stations[0].loops.iter().map(|r| r.len()).sum();
    // Drop the caps `loft` built from `start` and rebuild the far one from
    // `end`. Wall triangles come first, so truncating to the wall count is
    // exact rather than a search.
    let spans = stations.len() - 1;
    let ring_edges: usize = stations[0].loops.iter().map(|r| r.len()).sum();
    let wall_indices = spans * ring_edges * 6;
    mesh.indices.truncate(wall_indices);
    let (near_pts, near_tris) = crate::profile::triangulate(start)?;
    let (far_pts, far_tris) = crate::profile::triangulate(end)?;
    if near_pts.len() != per_station || far_pts.len() != per_station {
        return Err(GeomError::Degenerate(
            "tapered cap triangulation disagrees with the station rings".to_owned(),
        ));
    }
    let last = (spans * per_station) as u32;
    for t in &near_tris {
        mesh.indices.extend([t[0], t[2], t[1]]);
    }
    for t in &far_tris {
        mesh.indices.extend([last + t[0], last + t[1], last + t[2]]);
    }
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_measure::proximity::closest_point_on_triangle;

    /// The largest distance from a dense sampling of the bilinear patch to
    /// the loft's two triangles, measured by brute force.
    fn measured(a: Point3, b: Point3, c: Point3, d: Point3) -> Scalar {
        let n = 64;
        let mut worst: Scalar = 0.0;
        for i in 0..=n {
            for j in 0..=n {
                let (x, y) = (i as Scalar / n as Scalar, j as Scalar / n as Scalar);
                let p = (a * (1.0 - x) + b * x) * (1.0 - y) + (c * (1.0 - x) + d * x) * y;
                let near = [[a, b, d], [a, d, c]]
                    .iter()
                    .filter_map(|t| closest_point_on_triangle(p, *t).ok())
                    .map(|q| p.distance(q))
                    .fold(Scalar::INFINITY, Scalar::min);
                worst = worst.max(near);
            }
        }
        worst
    }

    /// A deterministic spread of quads: planar convex, planar concave
    /// (the loft's diagonal outside the quad), twisted, and collapsed.
    fn quads() -> Vec<[Point3; 4]> {
        let mut seed: u64 = 0x2310_5eed;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((seed >> 11) as Scalar / (1u64 << 53) as Scalar) * 2.0 - 1.0
        };
        let mut out = vec![
            // Planar trapezoid: a revolution's wall about an in-plane axis.
            [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.2, 1.0, 0.0),
            ],
            // Planar, but `d` folded back across the diagonal a-d: concave.
            [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.2, 0.2, 0.0),
            ],
            // A saddle: the hyperbolic paraboloid z = xy.
            [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 1.0),
            ],
            // Collapsed onto an axis: a and b coincide.
            [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.3),
            ],
        ];
        for _ in 0..200 {
            let mut q = [Point3::ZERO; 4];
            for p in &mut q {
                *p = Point3::new(next(), next(), 0.3 * next());
            }
            out.push(q);
        }
        out
    }

    #[test]
    fn the_quad_bound_holds_for_every_quad() {
        for [a, b, c, d] in quads() {
            let bound = quad_deviation(a, b, c, d);
            let got = measured(a, b, c, d);
            assert!(
                got <= bound * (1.0 + 1e-9) + 1e-12,
                "quad {a:?} {b:?} {c:?} {d:?}: patch {got} from its triangles, bound {bound}"
            );
        }
    }

    #[test]
    fn a_planar_convex_quad_costs_nothing() {
        let [a, b, c, d] = quads()[0];
        assert!(quad_deviation(a, b, c, d) <= 1e-15);
    }

    #[test]
    fn a_twisted_quad_is_bounded_tightly() {
        // z = xy: the patch's centre is |T| / 4 = 1/4 above the diagonal's
        // midpoint, and 1/(4 sqrt 2) from the triangles' planes.
        let [a, b, c, d] = quads()[2];
        let bound = quad_deviation(a, b, c, d);
        assert!((bound - 0.25).abs() <= 1e-12, "{bound}");
        let got = measured(a, b, c, d);
        assert!((got - 0.25 / 2f64.sqrt()).abs() <= 1e-3, "{got}");
    }
}
