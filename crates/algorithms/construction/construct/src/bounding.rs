//! Bounding volumes of 3D point sets (#118): the minimum enclosing sphere
//! and a containing oriented box.
//!
//! # Minimum enclosing sphere: exact choice, enclosed output
//!
//! Welzl's algorithm, iterative, over a fixed pseudo-random visiting order.
//! Whether a point lies in the sphere spanned by one to four support
//! points is decided exactly (intervals, then dyadics):
//!
//! - one support point: equality;
//! - two: the sign of `(p - a) . (p - b)`;
//! - three, the least sphere through them (centred in their plane): with
//!   `u = b - a`, `v = c - a`, `w = u x v` and
//!   `N = |u|^2 (v x w) + |v|^2 (w x u)` (so the centre is
//!   `a + N / 2|w|^2`), the sign of `|p - a|^2 |w|^2 - (p - a) . N`;
//! - four, their circumsphere: with `D = 2 u . (v x w)` and
//!   `M = |u|^2 (v x w) + |v|^2 (w x u) + |w|^2 (u x v)`, the sign of
//!   `|p - a|^2 D - 2 (p - a) . M` against the sign of `D`.
//!
//! So the support set is the exact minimum sphere's. The centre is enclosed
//! from the same exact numerators and denominators, and the radius is
//! rounded up, so the returned sphere contains the exact one and every
//! input point. [`SphereEvidence::error`] bounds the centre's distance from
//! the exact centre and the radius's excess over the exact radius.
//!
//! # Oriented box: certified containment, not certified optimality
//!
//! [`oriented_bounding_box`] tries these orientations and keeps the least
//! volume:
//!
//! - the axis-aligned box;
//! - the principal axes of the points' covariance;
//! - for each of the world axes, the principal axes, every face normal of
//!   the exact convex hull ([`crate::hull::convex_hull`]) and, for flat
//!   input, the plane's normal: that normal as one axis and the exact
//!   minimum-area rectangle ([`axiolid_overlay::minimum_area_rectangle`])
//!   of the points projected across it for the other two.
//!
//! What is certified is containment: for every input point `p` and every
//! axis, `|(p - centre) . axes[i]| <= half_extents[i]` holds exactly for the
//! returned `f64` values -- the extents are measured in outward-rounded
//! intervals. The volume is no more than the axis-aligned box's. What is
//! **not** claimed is the global minimum volume: an optimal box need not
//! have a face flush with a hull face (O'Rourke's exact algorithm, cubic in
//! the hull size, is not implemented), so the result is a good box, not
//! the best one.

use axiolid_core::{Point2, Point3, Vec3};
use axiolid_exact::{certify, Arith, Dyadic, Interval, SignExpr};
use axiolid_guarantees::Sign;

/// Why no bounding volume was built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum BoundingError {
    /// No points.
    Empty,
    /// A coordinate was not finite.
    NonFinite,
}

fn validate(points: &[Point3]) -> Result<(), BoundingError> {
    if points.is_empty() {
        return Err(BoundingError::Empty);
    }
    if !points.iter().all(|p| p.is_finite()) {
        return Err(BoundingError::NonFinite);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Minimum enclosing sphere
// ---------------------------------------------------------------------------

/// A sphere by its centre and radius.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnclosingSphere {
    /// Centre.
    pub centre: Point3,
    /// Radius; zero for a single distinct point.
    pub radius: f64,
}

/// Which points determine the sphere and how exact the output is.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct SphereEvidence {
    /// Indices into the input of the one to four points the exact minimum
    /// sphere passes through and is determined by, ascending.
    pub support: Vec<usize>,
    /// A bound on the distance between the returned centre and the exact
    /// one, and on the returned radius's excess over the exact radius. The
    /// returned radius is never below the exact radius, so the returned
    /// sphere contains every input point. Zero for a single point.
    pub error: f64,
}

/// The sphere and its evidence.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct MinimumSphere {
    /// The sphere.
    pub sphere: EnclosingSphere,
    /// Its support and error bound.
    pub evidence: SphereEvidence,
}

type V<T> = [T; 3];

fn diff<T: Arith>(a: Point3, b: Point3) -> V<T> {
    let f = T::from_f64;
    [
        f(a.x).sub(&f(b.x)),
        f(a.y).sub(&f(b.y)),
        f(a.z).sub(&f(b.z)),
    ]
}

fn dot<T: Arith>(a: &V<T>, b: &V<T>) -> T {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn cross<T: Arith>(a: &V<T>, b: &V<T>) -> V<T> {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

fn scale<T: Arith>(s: &T, a: &V<T>) -> V<T> {
    [s.mul(&a[0]), s.mul(&a[1]), s.mul(&a[2])]
}

fn plus<T: Arith>(a: &V<T>, b: &V<T>) -> V<T> {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}

/// The centre of the least sphere through `a b c` is `a + N / den`.
fn triangle_centre<T: Arith>(a: Point3, b: Point3, c: Point3) -> (V<T>, T) {
    let u: V<T> = diff(b, a);
    let v: V<T> = diff(c, a);
    let w = cross(&u, &v);
    let n = plus(
        &scale(&dot(&u, &u), &cross(&v, &w)),
        &scale(&dot(&v, &v), &cross(&w, &u)),
    );
    let ww = dot(&w, &w);
    (n, ww.add(&ww))
}

/// The circumcentre of `a b c d` is `a + M / D`.
fn tetrahedron_centre<T: Arith>(a: Point3, b: Point3, c: Point3, d: Point3) -> (V<T>, T) {
    let u: V<T> = diff(b, a);
    let v: V<T> = diff(c, a);
    let w: V<T> = diff(d, a);
    let m = plus(
        &plus(
            &scale(&dot(&u, &u), &cross(&v, &w)),
            &scale(&dot(&v, &v), &cross(&w, &u)),
        ),
        &scale(&dot(&w, &w), &cross(&u, &v)),
    );
    let det = dot(&u, &cross(&v, &w));
    (m, det.add(&det))
}

/// `(p - a) . (p - b)`.
struct Diametral {
    a: Point3,
    b: Point3,
    p: Point3,
}

impl SignExpr for Diametral {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        dot::<T>(&diff(self.p, self.a), &diff(self.p, self.b)).sign()
    }
}

/// `|p - a|^2 den - 2 (p - a) . N` for the centre `a + N / den`: its sign
/// against `den`'s says whether `p` is outside.
struct Beyond {
    support: [Point3; 4],
    count: usize,
    p: Point3,
}

impl Beyond {
    fn parts<T: Arith>(&self) -> (T, T) {
        let [a, b, c, d] = self.support;
        let (n, den) = if self.count == 3 {
            triangle_centre::<T>(a, b, c)
        } else {
            tetrahedron_centre::<T>(a, b, c, d)
        };
        let q: V<T> = diff(self.p, a);
        let lhs = dot(&q, &q).mul(&den);
        let rhs = dot(&q, &n);
        (lhs.sub(&rhs.add(&rhs)), den)
    }
}

impl SignExpr for Beyond {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        self.parts::<T>().0.sign()
    }
}

struct Denominator(Beyond);

impl SignExpr for Denominator {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        self.0.parts::<T>().1.sign()
    }
}

/// The inputs are finite, so the exact tier always decides.
fn sign<E: SignExpr>(e: &E) -> Sign {
    certify(e).unwrap_or(Sign::Zero)
}

/// Whether `p` lies in (or on) the sphere determined by `support`.
fn inside(points: &[Point3], support: &[usize], p: Point3) -> bool {
    match *support {
        [a] => points[a] == p,
        [a, b] => {
            sign(&Diametral {
                a: points[a],
                b: points[b],
                p,
            }) != Sign::Positive
        }
        _ => {
            let mut s = [points[support[0]]; 4];
            for (slot, &i) in s.iter_mut().zip(support) {
                *slot = points[i];
            }
            let beyond = Beyond {
                support: s,
                count: support.len(),
                p,
            };
            let side = sign(&beyond);
            let den = sign(&Denominator(beyond));
            // `|p - C|^2 - |a - C|^2` has the sign of `side * den`.
            side == Sign::Zero || den == Sign::Zero || side != den
        }
    }
}

/// A fixed pseudo-random permutation of `0..n` (Fisher-Yates driven by
/// splitmix64 from a constant seed).
fn visiting_order(n: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    let mut state: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut next = || {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    for i in (1..n).rev() {
        let j = (next() % (i as u64 + 1)) as usize;
        order.swap(i, j);
    }
    order
}

/// The support of the minimum sphere, by Welzl's iterative algorithm.
fn welzl(points: &[Point3]) -> Vec<usize> {
    let order = visiting_order(points.len());
    let mut support = vec![order[0]];
    for i in 1..order.len() {
        let pi = order[i];
        if inside(points, &support, points[pi]) {
            continue;
        }
        support = vec![pi];
        for j in 0..i {
            let pj = order[j];
            if inside(points, &support, points[pj]) {
                continue;
            }
            support = vec![pi, pj];
            for k in 0..j {
                let pk = order[k];
                if inside(points, &support, points[pk]) {
                    continue;
                }
                support = vec![pi, pj, pk];
                for &pl in &order[..k] {
                    if !inside(points, &support, points[pl]) {
                        support = vec![pi, pj, pk, pl];
                    }
                }
            }
        }
    }
    support
}

/// Enclosures of the exact centre of the sphere on `support`.
fn centre_enclosure(points: &[Point3], support: &[usize]) -> [Interval; 3] {
    let p = |i: usize| points[support[i]];
    let a = p(0);
    let at = [a.x, a.y, a.z];
    let (n, den): (V<Dyadic>, Dyadic) = match support.len() {
        1 => return at.map(Interval::point),
        2 => {
            let b = p(1);
            let half = Dyadic::from_f64(0.5);
            let mid = |s: f64, t: f64| {
                Dyadic::from_f64(s)
                    .add(&Dyadic::from_f64(t))
                    .mul(&half)
                    .enclosure()
            };
            return [mid(a.x, b.x), mid(a.y, b.y), mid(a.z, b.z)];
        }
        3 => triangle_centre(a, p(1), p(2)),
        _ => tetrahedron_centre(a, p(1), p(2), p(3)),
    };
    let den = den.enclosure();
    let mut out = [Interval::point(0.0); 3];
    for k in 0..3 {
        out[k] = Interval::point(at[k]).add(&n[k].enclosure().quotient(den));
    }
    out
}

/// Bounds on the distance from `p` to the nearest and farthest points of
/// the box `centre`.
fn distance_bounds(p: Point3, centre: &[Interval; 3]) -> (f64, f64) {
    let at = [p.x, p.y, p.z];
    let mut squared = Interval::point(0.0);
    for k in 0..3 {
        let d = Interval::point(at[k]).sub(&centre[k]);
        squared = squared.add(&d.mul(&d));
    }
    // `sqrt` is correctly rounded, so one step outward covers it.
    let low = squared.lo().max(0.0).sqrt().next_down().max(0.0);
    (low, squared.hi().sqrt().next_up())
}

/// The minimum sphere enclosing `points`.
///
/// # Errors
///
/// [`BoundingError::Empty`] for no points, [`BoundingError::NonFinite`] for
/// a coordinate that is not finite.
pub fn minimum_enclosing_sphere(points: &[Point3]) -> Result<MinimumSphere, BoundingError> {
    validate(points)?;
    let mut support = welzl(points);
    support.sort_unstable();
    if let [a] = *support {
        return Ok(MinimumSphere {
            sphere: EnclosingSphere {
                centre: points[a],
                radius: 0.0,
            },
            evidence: SphereEvidence {
                support,
                error: 0.0,
            },
        });
    }
    let enclosure = centre_enclosure(points, &support);
    let mid = |i: Interval| i.lo() + 0.5 * (i.hi() - i.lo());
    let centre = Point3::new(mid(enclosure[0]), mid(enclosure[1]), mid(enclosure[2]));
    // The sum of the per-axis gaps is at least their Euclidean length.
    let gap = |i: Interval, m: f64| (i.hi() - m).max(m - i.lo());
    let centre_error =
        (gap(enclosure[0], centre.x) + gap(enclosure[1], centre.y) + gap(enclosure[2], centre.z))
            .next_up()
            .next_up();
    let (radius_low, radius_high) = distance_bounds(points[support[0]], &enclosure);
    let mut radius = (radius_high + centre_error).next_up();
    // That sphere contains the exact one; confirm every point with
    // enclosures all the same.
    let at = [centre.x, centre.y, centre.z].map(Interval::point);
    for &p in points {
        radius = radius.max(distance_bounds(p, &at).1);
    }
    Ok(MinimumSphere {
        sphere: EnclosingSphere { centre, radius },
        evidence: SphereEvidence {
            support,
            error: (radius - radius_low).next_up(),
        },
    })
}

// ---------------------------------------------------------------------------
// Oriented bounding box
// ---------------------------------------------------------------------------

/// A box by its centre, three unit axes and the half extents along them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrientedBox {
    /// Centre.
    pub centre: Point3,
    /// Axes, right-handed and orthonormal up to rounding (see
    /// [`BoxEvidence::orthogonality`]).
    pub axes: [Vec3; 3],
    /// Half the side lengths along the axes.
    pub half_extents: [f64; 3],
}

impl OrientedBox {
    /// Its volume.
    #[must_use]
    pub fn volume(&self) -> f64 {
        8.0 * self.half_extents[0] * self.half_extents[1] * self.half_extents[2]
    }

    /// The eight corners: bit `k` of the index picks the sign along axis
    /// `k`.
    #[must_use]
    pub fn corners(&self) -> [Point3; 8] {
        std::array::from_fn(|i| {
            let mut p = self.centre;
            for k in 0..3 {
                let s = if i & (1 << k) == 0 { -1.0 } else { 1.0 };
                p += self.axes[k] * (s * self.half_extents[k]);
            }
            p
        })
    }
}

/// How the box was chosen and how far its axes are from orthonormal.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct BoxEvidence {
    /// Orientations tried.
    pub candidates: usize,
    /// The volume of the axis-aligned box, computed the same way; the
    /// returned box's [`OrientedBox::volume`] is never above it.
    pub axis_aligned_volume: f64,
    /// An upper bound on `|axes[i] . axes[j] - [i == j]|` over all pairs,
    /// measured exactly: how far the rounded axes are from orthonormal.
    pub orthogonality: f64,
}

/// The box and its evidence.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct OrientedBoundingBox {
    /// The box. Every input point lies in it, exactly.
    pub bounding_box: OrientedBox,
    /// How it was chosen.
    pub evidence: BoxEvidence,
}

fn unit(v: Vec3) -> Option<Vec3> {
    let length = v.length();
    (length > 0.0 && length.is_finite()).then(|| v / length)
}

/// A right-handed frame whose third axis is `normal` and whose first lies
/// along `first` projected off it.
fn frame(normal: Vec3, first: Vec3) -> Option<[Vec3; 3]> {
    let n = unit(normal)?;
    let a = unit(first - n * first.dot(n))?;
    let b = unit(n.cross(a))?;
    Some([a, b, n])
}

/// Any unit vector perpendicular to `n`.
fn perpendicular(n: Vec3) -> Vec3 {
    let helper = if n.x.abs() <= n.y.abs() && n.x.abs() <= n.z.abs() {
        Vec3::X
    } else if n.y.abs() <= n.z.abs() {
        Vec3::Y
    } else {
        Vec3::Z
    };
    unit(n.cross(helper)).unwrap_or(Vec3::X)
}

/// The frame with `normal` as its third axis and the minimum-area
/// rectangle of `points` projected across it for the other two.
fn flush_frame(points: &[Point3], normal: Vec3) -> Option<[Vec3; 3]> {
    let n = unit(normal)?;
    let e1 = perpendicular(n);
    let e2 = n.cross(e1);
    let projected: Vec<Point2> = points
        .iter()
        .map(|p| Point2::new(p.dot(e1), p.dot(e2)))
        .collect();
    let rectangle = axiolid_overlay::minimum_area_rectangle(&projected).ok()?;
    let [r, _] = rectangle.rectangle.axes;
    frame(n, e1 * r.x + e2 * r.y)
}

/// The eigenvectors of the points' covariance, by cyclic Jacobi rotations.
fn principal_axes(points: &[Point3]) -> [Vec3; 3] {
    let n = points.len() as f64;
    let mean = points.iter().fold(Vec3::ZERO, |s, p| s + *p) / n;
    let mut c = [[0.0f64; 3]; 3];
    for p in points {
        let d = *p - mean;
        let d = [d.x, d.y, d.z];
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] += d[i] * d[j];
            }
        }
    }
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..32 {
        let off = c[0][1].abs() + c[0][2].abs() + c[1][2].abs();
        if off == 0.0 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if c[p][q] == 0.0 {
                continue;
            }
            let theta = (c[q][q] - c[p][p]) / (2.0 * c[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let cs = 1.0 / (t * t + 1.0).sqrt();
            let sn = t * cs;
            for row in &mut c {
                let (a, b) = (row[p], row[q]);
                row[p] = cs * a - sn * b;
                row[q] = sn * a + cs * b;
            }
            for k in 0..3 {
                let (a, b) = (c[p][k], c[q][k]);
                c[p][k] = cs * a - sn * b;
                c[q][k] = sn * a + cs * b;
            }
            for row in &mut v {
                let (a, b) = (row[p], row[q]);
                row[p] = cs * a - sn * b;
                row[q] = sn * a + cs * b;
            }
        }
    }
    std::array::from_fn(|k| Vec3::new(v[0][k], v[1][k], v[2][k]))
}

/// The volume of the box along `axes` holding `points`, in plain floating
/// point, with `pad` added to every side: only for ranking candidates,
/// never for the returned box. The pad, a tiny fraction of the points'
/// size, makes flat input rank by area and collinear input by length
/// instead of by rounding noise in a vanishing volume.
fn ranking_volume(points: &[Point3], axes: &[Vec3; 3], pad: f64) -> f64 {
    let mut volume = 1.0;
    for a in axes {
        let (low, high) = points
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), p| {
                let s = p.dot(*a);
                (l.min(s), h.max(s))
            });
        volume *= high - low + pad;
    }
    volume
}

/// The box along `axes` holding every point. Projections are exact
/// dyadics, and the half extents are rounded up from exact values, so
/// `|(p - centre) . axes[i]| <= half_extents[i]` holds exactly; an extent
/// that is representable (any axis-aligned box on representable
/// coordinates) is returned exactly.
fn fit(points: &[Point3], axes: [Vec3; 3]) -> OrientedBox {
    let e = Dyadic::from_f64;
    let project = |p: Point3, a: Vec3| {
        e(p.x)
            .mul(&e(a.x))
            .add(&e(p.y).mul(&e(a.y)))
            .add(&e(p.z).mul(&e(a.z)))
    };
    let is_less = |x: &Dyadic, y: &Dyadic| x.sub(y).sign() == Some(Sign::Negative);
    let mut low: Vec<Dyadic> = axes.iter().map(|a| project(points[0], *a)).collect();
    let mut high = low.clone();
    for &p in &points[1..] {
        for k in 0..3 {
            let s = project(p, axes[k]);
            if is_less(&s, &low[k]) {
                low[k] = s;
            } else if is_less(&high[k], &s) {
                high[k] = s;
            }
        }
    }
    let half = e(0.5);
    let mut centre = Point3::ZERO;
    for k in 0..3 {
        centre += axes[k] * low[k].add(&high[k]).mul(&half).to_f64();
    }
    let half_extents = std::array::from_fn(|k| {
        let c = project(centre, axes[k]);
        let (above, below) = (high[k].sub(&c), c.sub(&low[k]));
        let widest = if is_less(&above, &below) {
            below
        } else {
            above
        };
        round_up(&widest).max(0.0)
    });
    OrientedBox {
        centre,
        axes,
        half_extents,
    }
}

/// The least double no smaller than `d`.
fn round_up(d: &Dyadic) -> f64 {
    let mut x = d.to_f64();
    let e = Dyadic::from_f64;
    while e(x).sub(d).sign() == Some(Sign::Negative) {
        x = x.next_up();
    }
    while e(x.next_down()).sub(d).sign() != Some(Sign::Negative) && x.next_down() >= 0.0 {
        x = x.next_down();
    }
    x
}

/// `max |a_i . a_j - [i == j]|`, measured exactly and rounded up.
fn orthogonality(axes: &[Vec3; 3]) -> f64 {
    let e = Dyadic::from_f64;
    let mut worst = 0.0f64;
    for i in 0..3 {
        for j in i..3 {
            let (a, b) = (axes[i], axes[j]);
            let mut d = e(a.x)
                .mul(&e(b.x))
                .add(&e(a.y).mul(&e(b.y)))
                .add(&e(a.z).mul(&e(b.z)));
            if i == j {
                d = d.sub(&e(1.0));
            }
            if d.sign() == Some(Sign::Negative) {
                d = d.neg();
            }
            worst = worst.max(d.enclosure().hi());
        }
    }
    worst
}

/// A box holding every point, of volume no more than the axis-aligned box.
///
/// See the module documentation for the orientations tried and for what is
/// and is not certified: containment is, minimum volume is not.
///
/// # Errors
///
/// [`BoundingError::Empty`] for no points, [`BoundingError::NonFinite`] for
/// a coordinate that is not finite.
pub fn oriented_bounding_box(points: &[Point3]) -> Result<OrientedBoundingBox, BoundingError> {
    validate(points)?;
    // Candidate volumes are compared on the hull's vertices, which hold the
    // same extremes; the chosen box is then fitted to every point.
    let hull = crate::hull::convex_hull(points).ok();
    let extreme: &[Point3] = hull.as_ref().map_or(points, |h| &h.positions);

    let world = [Vec3::X, Vec3::Y, Vec3::Z];
    let principal = principal_axes(points);
    let mut frames: Vec<[Vec3; 3]> = vec![world];
    if let Some(f) = frame(principal[0].cross(principal[1]), principal[0]) {
        frames.push(f);
    }
    let mut normals: Vec<Vec3> = world.to_vec();
    normals.extend(principal);
    match &hull {
        Some(h) => {
            for t in h.indices.chunks_exact(3) {
                let [a, b, c] = [0, 1, 2].map(|k| h.positions[t[k] as usize]);
                normals.push((b - a).cross(c - a));
            }
        }
        None => normals.push(flat_normal(points)),
    }
    let mut seen: Vec<Vec3> = Vec::new();
    for normal in normals {
        let Some(n) = unit(normal) else { continue };
        // A normal and its opposite give the same box.
        let n = if (n.x, n.y, n.z) < (0.0, 0.0, 0.0) {
            -n
        } else {
            n
        };
        if seen.contains(&n) {
            continue;
        }
        seen.push(n);
        if let Some(f) = flush_frame(extreme, n) {
            frames.push(f);
        }
    }

    let axis_aligned_volume = fit(points, world).volume();
    let mut best = world;
    // The widest side of the axis-aligned box.
    let size = world
        .iter()
        .map(|a| {
            let along = extreme.iter().map(|p| p.dot(*a));
            along.clone().fold(f64::NEG_INFINITY, f64::max) - along.fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max);
    let pad = 1e-9 * size;
    let mut best_volume = ranking_volume(extreme, &world, pad);
    for f in &frames[1..] {
        let volume = ranking_volume(extreme, f, pad);
        if volume < best_volume {
            best = *f;
            best_volume = volume;
        }
    }
    let mut chosen = fit(points, best);
    // Fitting every point can only widen a box by rounding; never let that
    // push it past the axis-aligned one.
    if chosen.volume() > axis_aligned_volume {
        best = world;
        chosen = fit(points, world);
    }
    Ok(OrientedBoundingBox {
        bounding_box: chosen,
        evidence: BoxEvidence {
            candidates: frames.len(),
            axis_aligned_volume,
            orthogonality: orthogonality(&best),
        },
    })
}

/// A normal of flat (coplanar or collinear) input: the cross product of
/// the direction to the farthest point and to the point farthest from
/// that line; for collinear input, any perpendicular of the line.
fn flat_normal(points: &[Point3]) -> Vec3 {
    let a = points[0];
    let farthest = |from: &dyn Fn(Point3) -> f64| {
        points
            .iter()
            .copied()
            .max_by(|p, q| from(*p).total_cmp(&from(*q)))
            .unwrap_or(a)
    };
    let b = farthest(&|p| (p - a).length_squared());
    let d = b - a;
    let c = farthest(&|p| (p - a).cross(d).length_squared());
    let n = d.cross(c - a);
    if unit(n).is_some() {
        n
    } else if let Some(d) = unit(d) {
        perpendicular(d)
    } else {
        Vec3::Z
    }
}
