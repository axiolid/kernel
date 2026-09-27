//! Planar regions of a triangle mesh (#131).
//!
//! # Grown, then certified
//!
//! Regions grow from the largest unassigned triangle across shared edges,
//! taking a neighbour when its normal is within the angle of the region's
//! and its corners within the distance of the region's plane; the plane is
//! refitted (area-weighted normal and centroid) as the region grows.
//!
//! Growing decides in floating point; what is reported is then proven. For
//! the plane given -- a point and a normal, both `f64`, so an exact plane --
//! [`DetectedPlane::deviation`] bounds the distance of every corner of every
//! member triangle from it, computed in outward-rounded intervals. A
//! triangle that pushes the bound past the requested distance (plus the
//! few ulps any fitted plane is off, so exactly flat regions hold together
//! at distance zero) is peeled off the region until the bound holds. So
//! every region's corners, and so its whole surface, lie within
//! `deviation` of its plane.
//!
//! [`DetectedPlane::coplanar`] says more when it holds: every corner lies
//! on one plane exactly, decided by exact `orient3d`.

use axiolid_core::{Point3, Vec3};
use axiolid_exact::{certify, Arith, Dyadic, SignExpr};
use axiolid_guarantees::Sign;
use axiolid_mesh::TriangleMeshView;
use std::collections::HashMap;

/// How far a region may depart from a plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneTolerance {
    /// Greatest distance of a corner from the region's plane.
    pub distance: f64,
    /// Greatest angle, in radians, between a member triangle's normal and
    /// the plane's, while growing.
    pub angle: f64,
}

/// One planar region.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct DetectedPlane {
    /// Member triangle indices, ascending.
    pub triangles: Vec<u32>,
    /// A point of the plane: the members' area-weighted centroid.
    pub point: Point3,
    /// Unit normal (to rounding), the members' area-weighted mean, on the
    /// side the triangles face.
    pub normal: Vec3,
    /// A certified upper bound on the distance of any member corner from
    /// the plane through `point` with normal `normal`: at most the
    /// requested distance plus a few ulps of the region's coordinates (a
    /// lone triangle's is its own rounding).
    pub deviation: f64,
    /// Whether every member corner lies on one plane exactly.
    pub coplanar: bool,
    /// Summed area of the members (rounded).
    pub area: f64,
}

/// Why no segmentation was made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PlaneError {
    /// A coordinate is not finite.
    NonFinite,
    /// The distance or angle is negative or not finite.
    InvalidTolerance,
}

/// The mesh cut into planar regions, largest first; triangles of no area
/// belong to none.
///
/// # Errors
///
/// [`PlaneError`] for non-finite input or tolerance.
pub fn detect_planes<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    tolerance: PlaneTolerance,
) -> Result<Vec<DetectedPlane>, PlaneError> {
    let PlaneTolerance { distance, angle } = tolerance;
    if !(distance.is_finite() && distance >= 0.0 && angle.is_finite() && angle >= 0.0) {
        return Err(PlaneError::InvalidTolerance);
    }
    if (0..mesh.position_count()).any(|i| !mesh.position(i).is_finite()) {
        return Err(PlaneError::NonFinite);
    }
    let n = mesh.triangle_count();
    let corners: Vec<[u64; 3]> = (0..n).map(|t| mesh.triangle(t)).collect();
    let at = |i: u64| mesh.position(i as usize);
    let tri: Vec<[Point3; 3]> = corners.iter().map(|c| c.map(at)).collect();
    let normals: Vec<Vec3> = tri
        .iter()
        .map(|t| (t[1] - t[0]).cross(t[2] - t[0]))
        .collect();
    let areas: Vec<f64> = normals.iter().map(|v| 0.5 * v.length()).collect();
    // Edge neighbours by vertex index.
    let mut edges: HashMap<(u64, u64), Vec<usize>> = HashMap::new();
    for (t, c) in corners.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (c[k], c[(k + 1) % 3]);
            edges.entry((a.min(b), a.max(b))).or_default().push(t);
        }
    }
    let neighbours = |t: usize| -> Vec<usize> {
        let c = corners[t];
        let mut out = Vec::new();
        for k in 0..3 {
            let (a, b) = (c[k], c[(k + 1) % 3]);
            for &o in &edges[&(a.min(b), a.max(b))] {
                if o != t && !out.contains(&o) {
                    out.push(o);
                }
            }
        }
        out
    };
    let cos_limit = angle.min(std::f64::consts::PI).cos();
    let mut order: Vec<usize> = (0..n).filter(|&t| areas[t] > 0.0).collect();
    order.sort_by(|&a, &b| areas[b].total_cmp(&areas[a]).then(a.cmp(&b)));
    let mut owner = vec![usize::MAX; n];
    let mut out: Vec<DetectedPlane> = Vec::new();
    for &seed in &order {
        if owner[seed] != usize::MAX {
            continue;
        }
        let id = out.len();
        let mut members = vec![seed];
        owner[seed] = id;
        let mut fit = Fit::of(&members, &tri, &normals, &areas);
        let mut frontier = vec![seed];
        let mut since_refit = 0;
        while let Some(t) = frontier.pop() {
            for o in neighbours(t) {
                if owner[o] != usize::MAX || areas[o] == 0.0 {
                    continue;
                }
                let unit = normals[o] / (2.0 * areas[o]);
                if unit.dot(fit.normal) < cos_limit {
                    continue;
                }
                if tri[o].iter().any(|&v| fit.distance(v) > distance) {
                    continue;
                }
                owner[o] = id;
                members.push(o);
                frontier.push(o);
                since_refit += 1;
                if since_refit >= 16 {
                    fit = Fit::of(&members, &tri, &normals, &areas);
                    since_refit = 0;
                }
            }
        }
        // Certify; peel the worst member off while the bound fails.
        loop {
            fit = Fit::of(&members, &tri, &normals, &areas);
            let bounds: Vec<f64> = members.iter().map(|&t| fit.bound(&tri[t])).collect();
            let (worst, bound) = bounds.iter().enumerate().fold(
                (0, 0.0f64),
                |(wi, wb), (i, &b)| if b > wb { (i, b) } else { (wi, wb) },
            );
            if bound <= distance + fit.slack || members.len() == 1 {
                // A single triangle's own plane holds it; its bound is
                // rounding only, and is reported as it is.
                let mut sorted: Vec<u32> = members.iter().map(|&t| t as u32).collect();
                sorted.sort_unstable();
                let coplanar = coplanar(members.iter().flat_map(|&t| tri[t]));
                out.push(DetectedPlane {
                    triangles: sorted,
                    point: fit.point,
                    normal: fit.normal,
                    deviation: bound,
                    coplanar,
                    area: members.iter().map(|&t| areas[t]).sum(),
                });
                break;
            }
            let t = members.swap_remove(worst);
            owner[t] = usize::MAX;
        }
        // Peeled triangles are seeded again in their turn; any peeled
        // before their seed's turn came are picked up by the loop below.
    }
    // Triangles peeled after their turn in `order` passed.
    let mut left: Vec<usize> = order
        .iter()
        .copied()
        .filter(|&t| owner[t] == usize::MAX)
        .collect();
    while let Some(t) = left.pop() {
        if owner[t] != usize::MAX {
            continue;
        }
        let fit = Fit::of(&[t], &tri, &normals, &areas);
        owner[t] = out.len();
        out.push(DetectedPlane {
            triangles: vec![t as u32],
            point: fit.point,
            normal: fit.normal,
            deviation: fit.bound(&tri[t]),
            coplanar: true,
            area: areas[t],
        });
    }
    out.sort_by(|a, b| {
        b.area
            .total_cmp(&a.area)
            .then(a.triangles[0].cmp(&b.triangles[0]))
    });
    Ok(out)
}

/// A region's plane: area-weighted centroid and unit mean normal, and
/// the rounding a plane fitted to it cannot avoid.
struct Fit {
    point: Point3,
    normal: Vec3,
    slack: f64,
}

impl Fit {
    fn of(members: &[usize], tri: &[[Point3; 3]], normals: &[Vec3], areas: &[f64]) -> Self {
        let mut total = 0.0;
        let mut weighted = Vec3::ZERO;
        let mut direction = Vec3::ZERO;
        // About the first corner, for conditioning far from the origin.
        let base = tri[members[0]][0];
        for &t in members {
            let [a, b, c] = tri[t];
            weighted += ((a - base) + (b - base) + (c - base)) * (areas[t] / 3.0);
            total += areas[t];
            // `normals[t]` has length twice the area: an area weighting.
            direction += normals[t];
        }
        let reach = members
            .iter()
            .flat_map(|&t| tri[t])
            .map(|v| (v - base).abs().max_element())
            .fold(0.0, f64::max);
        let far = base.abs().max_element() + reach;
        Self {
            point: base + weighted / total,
            normal: direction.normalize(),
            slack: 32.0 * f64::EPSILON * far,
        }
    }

    fn distance(&self, v: Point3) -> f64 {
        self.normal.dot(v - self.point).abs()
    }

    /// A certified upper bound on the distance of each corner from the
    /// plane through `point` with normal `normal`, both taken exactly:
    /// `|n . (v - p)| / |n|`, in outward-rounded intervals.
    fn bound(&self, t: &[Point3; 3]) -> f64 {
        let n = [self.normal.x, self.normal.y, self.normal.z].map(Iv::point);
        let norm2 = n[0].mul(n[0]).add(n[1].mul(n[1])).add(n[2].mul(n[2]));
        t.iter()
            .map(|v| {
                let d = [v.x, v.y, v.z];
                let p = [self.point.x, self.point.y, self.point.z];
                let dot = (0..3)
                    .map(|k| n[k].mul(Iv::point(d[k]).sub(Iv::point(p[k]))))
                    .fold(Iv::point(0.0), Iv::add);
                let top = dot.lo.abs().max(dot.hi.abs());
                // |dot| / sqrt(norm2), rounded up.
                let q = (top * top).next_up() / norm2.lo;
                q.next_up().sqrt().next_up()
            })
            .fold(0.0, f64::max)
    }
}

/// Whether the points all lie on one plane, exactly.
fn coplanar(points: impl Iterator<Item = Point3>) -> bool {
    let points: Vec<Point3> = points.collect();
    let Some((a, b, c)) = spanning(&points) else {
        return true;
    };
    points
        .iter()
        .all(|&d| certify(&Orient3 { p: [a, b, c, d] }).ok() == Some(Sign::Zero))
}

/// Three points not on one line, if any.
fn spanning(points: &[Point3]) -> Option<(Point3, Point3, Point3)> {
    let a = points[0];
    let b = *points.iter().find(|&&p| p != a)?;
    let c = *points.iter().find(|&&p| {
        let (u, v) = (b - a, p - a);
        let cross = [
            exact(u.y)
                .mul(&exact(v.z))
                .sub(&exact(u.z).mul(&exact(v.y))),
            exact(u.z)
                .mul(&exact(v.x))
                .sub(&exact(u.x).mul(&exact(v.z))),
            exact(u.x)
                .mul(&exact(v.y))
                .sub(&exact(u.y).mul(&exact(v.x))),
        ];
        cross.iter().any(|x| x.sign() != Some(Sign::Zero))
    })?;
    Some((a, b, c))
}

fn exact(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

struct Orient3 {
    p: [Point3; 4],
}

impl SignExpr for Orient3 {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let q = |i: usize| {
            let p = self.p[i];
            [T::from_f64(p.x), T::from_f64(p.y), T::from_f64(p.z)]
        };
        let (a, b, c, d) = (q(0), q(1), q(2), q(3));
        let sub = |x: &[T; 3], y: &[T; 3]| [x[0].sub(&y[0]), x[1].sub(&y[1]), x[2].sub(&y[2])];
        let (u, v, w) = (sub(&b, &a), sub(&c, &a), sub(&d, &a));
        let n = [
            u[1].mul(&v[2]).sub(&u[2].mul(&v[1])),
            u[2].mul(&v[0]).sub(&u[0].mul(&v[2])),
            u[0].mul(&v[1]).sub(&u[1].mul(&v[0])),
        ];
        n[0].mul(&w[0])
            .add(&n[1].mul(&w[1]))
            .add(&n[2].mul(&w[2]))
            .sign()
    }
}

/// An outward-rounded interval.
#[derive(Debug, Clone, Copy)]
struct Iv {
    lo: f64,
    hi: f64,
}

impl Iv {
    fn point(v: f64) -> Self {
        Self { lo: v, hi: v }
    }

    fn outward(lo: f64, hi: f64) -> Self {
        Self {
            lo: lo.next_down(),
            hi: hi.next_up(),
        }
    }

    fn add(self, o: Self) -> Self {
        Self::outward(self.lo + o.lo, self.hi + o.hi)
    }

    fn sub(self, o: Self) -> Self {
        Self::outward(self.lo - o.hi, self.hi - o.lo)
    }

    fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self::outward(
            p.iter().copied().fold(f64::INFINITY, f64::min),
            p.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    }
}
