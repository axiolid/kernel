//! Certified volume shared by two closed meshes, and the volume of their
//! difference (#183).
//!
//! # Prisms, not a boolean
//!
//! No mesh is built. For a closed, consistently oriented, non-self-
//! intersecting triangle mesh `A` lying above a plane `z = z0`, a point is
//! inside `A` exactly when the faces above it, counted `+1` for faces
//! facing up and `-1` for faces facing down, sum to one. So, almost
//! everywhere,
//!
//! ```text
//! 1_A = sum over faces f of  s_f * 1{ (x, y) in f', z0 < z < h_f(x, y) }
//! ```
//!
//! with `f'` the face's shadow on the plane, `h_f` the height of its plane
//! and `s_f` the sign of the shadow's orientation (vertical faces cast no
//! shadow and drop out). Multiplying two such sums and integrating,
//!
//! ```text
//! vol(A and B) = sum over f in A, g in B of  s_f s_g * integral over f' and g' of (min(h_f, h_g) - z0)
//! ```
//!
//! Each term is a convex polygon -- the overlap of two triangles --
//! split by the line where the two planes cross, and a linear function
//! integrated over each piece.
//!
//! # Exact decisions, enclosed arithmetic
//!
//! Every vertex of every piece before the split is an input vertex or the
//! crossing of two input edge lines, so which side of a line or of the
//! plane crossing it lies on is the sign of a polynomial in the input
//! coordinates: decided by interval arithmetic and else exactly in
//! dyadics. The numbers -- crossing points, heights, areas -- are then
//! computed in outward-rounded intervals, so the sum is an interval that
//! contains the true volume. Touching bodies and shared faces cancel
//! exactly in the decisions; their volume comes out as a small interval
//! around zero, and is clamped to be no less than zero.

use axiolid_core::{Aabb, Point2, Point3, Tolerance};
use axiolid_exact::{certify, Arith, SignExpr};
use axiolid_guarantees::Sign;
use axiolid_heal::self_intersections;
use axiolid_mesh::{audit_mesh, TriangleMeshView};
use axiolid_spatial::{Bvh, SpatialIndex, SpatialItem};
use core::ops::ControlFlow;

/// A closed interval of volumes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumeInterval {
    /// No greater than the true volume, and never negative.
    pub lower: f64,
    /// No less than the true volume.
    pub upper: f64,
}

impl VolumeInterval {
    /// Whether `volume` lies in the interval.
    #[must_use]
    pub fn contains(&self, volume: f64) -> bool {
        self.lower <= volume && volume <= self.upper
    }

    /// `upper - lower`.
    #[must_use]
    pub fn width(&self) -> f64 {
        self.upper - self.lower
    }
}

/// Which operand an error is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Operand {
    /// The first mesh given.
    First,
    /// The second mesh given.
    Second,
}

/// Why no volume was bracketed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum OverlapError {
    /// A coordinate is not finite.
    NonFinite {
        /// Which mesh.
        operand: Operand,
    },
    /// The mesh is not a closed, consistently wound two-manifold without
    /// degenerate triangles, so it encloses no definite solid.
    NotClosed {
        /// Which mesh.
        operand: Operand,
    },
    /// Two triangles of the mesh intersect other than along shared
    /// vertices, so points inside it are not well defined.
    SelfIntersecting {
        /// Which mesh.
        operand: Operand,
        /// The first pair found, lower index first.
        triangles: [u32; 2],
    },
    /// The mesh's enclosed volume is not certified nonzero.
    NoVolume {
        /// Which mesh.
        operand: Operand,
    },
}

/// The volume enclosed by a closed mesh, certified. Either winding is
/// accepted.
///
/// # Errors
///
/// [`OverlapError`] for a non-finite, open or self-intersecting mesh, or
/// one whose volume is not certified nonzero.
pub fn enclosed_volume<M: TriangleMeshView + ?Sized>(
    mesh: &M,
) -> Result<VolumeInterval, OverlapError> {
    let solid = Solid::new(mesh, Operand::First, None)?;
    Ok(clamp(solid.volume, 0.0, f64::INFINITY))
}

/// The volume two closed meshes share, certified.
///
/// # Errors
///
/// [`OverlapError`] naming the operand that is non-finite, open, self-
/// intersecting or without volume.
pub fn intersection_volume<A, B>(first: &A, second: &B) -> Result<VolumeInterval, OverlapError>
where
    A: TriangleMeshView + ?Sized,
    B: TriangleMeshView + ?Sized,
{
    let (a, b) = solids(first, second)?;
    Ok(shared(&a, &b))
}

/// The volume of the first mesh outside the second, certified.
///
/// # Errors
///
/// As [`intersection_volume`].
pub fn difference_volume<A, B>(first: &A, second: &B) -> Result<VolumeInterval, OverlapError>
where
    A: TriangleMeshView + ?Sized,
    B: TriangleMeshView + ?Sized,
{
    let (a, b) = solids(first, second)?;
    let both = shared(&a, &b);
    let rest = a.volume.sub(Iv::new(both.lower, both.upper));
    Ok(clamp(rest, 0.0, a.volume.hi))
}

fn solids<A, B>(first: &A, second: &B) -> Result<(Solid, Solid), OverlapError>
where
    A: TriangleMeshView + ?Sized,
    B: TriangleMeshView + ?Sized,
{
    let base = lowest(first).min(lowest(second));
    let a = Solid::new(first, Operand::First, Some(base))?;
    let b = Solid::new(second, Operand::Second, Some(base))?;
    Ok((a, b))
}

/// The shared volume, clamped to what either solid holds.
fn shared(a: &Solid, b: &Solid) -> VolumeInterval {
    let items: Vec<SpatialItem<u32>> = b
        .faces
        .iter()
        .enumerate()
        .map(|(i, f)| SpatialItem::new(i as u32, f.shadow()))
        .collect();
    let bvh = Bvh::build(items);
    let mut total = Iv::point(0.0);
    for f in &a.faces {
        bvh.visit_aabb(&f.shadow(), &mut |j: &u32| {
            let g = &b.faces[*j as usize];
            let term = overlap(f, g, a.base);
            total = if f.up == g.up {
                total.add(term)
            } else {
                total.sub(term)
            };
            ControlFlow::Continue(())
        });
    }
    clamp(total, 0.0, a.volume.hi.min(b.volume.hi))
}

fn clamp(v: Iv, low: f64, high: f64) -> VolumeInterval {
    VolumeInterval {
        lower: v.lo.max(low).min(high),
        upper: v.hi.min(high).max(low),
    }
}

fn lowest<M: TriangleMeshView + ?Sized>(mesh: &M) -> f64 {
    (0..mesh.position_count())
        .map(|i| mesh.position(i).z)
        .fold(f64::INFINITY, f64::min)
}

/// A face that casts a shadow, its corners counter-clockwise in the
/// plane, and whether its outward side (after orienting the solid) faces
/// up.
#[derive(Debug, Clone, Copy)]
struct Face {
    p: [Point3; 3],
    up: bool,
}

impl Face {
    fn shadow(&self) -> Aabb {
        let mut b = Aabb::from_point(Point3::new(self.p[0].x, self.p[0].y, 0.0));
        for q in &self.p[1..] {
            b.extend(Point3::new(q.x, q.y, 0.0));
        }
        b
    }

    fn flat(&self, i: usize) -> Point2 {
        Point2::new(self.p[i].x, self.p[i].y)
    }

    /// The edge lines, each from a corner to the next.
    fn lines(&self) -> [Line; 3] {
        [0, 1, 2].map(|i| Line {
            a: self.flat(i),
            b: self.flat((i + 1) % 3),
        })
    }

    /// The plane's normal, enclosed: `(p1 - p0) x (p2 - p0)`.
    fn normal(&self) -> [Iv; 3] {
        let d = |i: usize, k: usize| Iv::point(self.p[i][k]).sub(Iv::point(self.p[0][k]));
        let (u, v) = ([d(1, 0), d(1, 1), d(1, 2)], [d(2, 0), d(2, 1), d(2, 2)]);
        [
            u[1].mul(v[2]).sub(u[2].mul(v[1])),
            u[2].mul(v[0]).sub(u[0].mul(v[2])),
            u[0].mul(v[1]).sub(u[1].mul(v[0])),
        ]
    }

    /// Height of the plane above `(x, y)`, less `base`.
    fn height(&self, at: [Iv; 2], base: f64) -> Iv {
        let n = self.normal();
        let dx = at[0].sub(Iv::point(self.p[0].x));
        let dy = at[1].sub(Iv::point(self.p[0].y));
        let lift = n[0].mul(dx).add(n[1].mul(dy)).div(n[2]);
        Iv::point(self.p[0].z).sub(Iv::point(base)).sub(lift)
    }
}

struct Solid {
    faces: Vec<Face>,
    volume: Iv,
    base: f64,
}

impl Solid {
    fn new<M: TriangleMeshView + ?Sized>(
        mesh: &M,
        operand: Operand,
        base: Option<f64>,
    ) -> Result<Self, OverlapError> {
        if (0..mesh.position_count()).any(|i| !mesh.position(i).is_finite()) {
            return Err(OverlapError::NonFinite { operand });
        }
        if !audit_mesh(mesh, Tolerance::ZERO).is_closed_two_manifold() {
            return Err(OverlapError::NotClosed { operand });
        }
        if let Some(pair) = self_intersections(mesh).first() {
            return Err(OverlapError::SelfIntersecting {
                operand,
                triangles: [pair.first, pair.second],
            });
        }
        let base = base.unwrap_or_else(|| lowest(mesh));
        let mut faces = Vec::new();
        for t in 0..mesh.triangle_count() {
            let [a, b, c] = mesh.triangle(t).map(|i| mesh.position(i as usize));
            let flat = |p: Point3| Point2::new(p.x, p.y);
            match sign(&Orient {
                a: flat(a),
                b: flat(b),
                c: flat(c),
            }) {
                Sign::Positive => faces.push(Face {
                    p: [a, b, c],
                    up: true,
                }),
                Sign::Negative => faces.push(Face {
                    p: [a, c, b],
                    up: false,
                }),
                _ => {}
            }
        }
        // Volume: each shadow's area times the mean height above the base.
        let mut volume = Iv::point(0.0);
        for f in &faces {
            let area = triangle_area(f.flat(0), f.flat(1), f.flat(2));
            let rise = [0, 1, 2]
                .map(|i| Iv::point(f.p[i].z).sub(Iv::point(base)))
                .into_iter()
                .fold(Iv::point(0.0), Iv::add)
                .div(Iv::point(3.0));
            let term = area.mul(rise);
            volume = if f.up {
                volume.add(term)
            } else {
                volume.sub(term)
            };
        }
        // Wound inward: turn every face over.
        if volume.hi < 0.0 {
            volume = Iv::point(0.0).sub(volume);
            for f in &mut faces {
                f.up = !f.up;
            }
        } else if volume.lo <= 0.0 {
            return Err(OverlapError::NoVolume { operand });
        }
        Ok(Self {
            faces,
            volume,
            base,
        })
    }
}

fn triangle_area(a: Point2, b: Point2, c: Point2) -> Iv {
    let (ux, uy) = (
        Iv::point(b.x).sub(Iv::point(a.x)),
        Iv::point(b.y).sub(Iv::point(a.y)),
    );
    let (vx, vy) = (
        Iv::point(c.x).sub(Iv::point(a.x)),
        Iv::point(c.y).sub(Iv::point(a.y)),
    );
    ux.mul(vy).sub(uy.mul(vx)).mul(Iv::point(0.5))
}

/// The integral of `min(h_f, h_g) - base` over the overlap of the two
/// shadows.
fn overlap(f: &Face, g: &Face, base: f64) -> Iv {
    // The overlap: f's shadow clipped by each of g's edge lines, exactly.
    let mut poly: Vec<(Vertex, Line)> = (0..3)
        .map(|i| (Vertex::Input(f.flat(i)), f.lines()[i]))
        .collect();
    for line in g.lines() {
        poly = clip(&poly, line);
        if poly.len() < 3 {
            return Iv::point(0.0);
        }
    }
    // Which plane is lower at each vertex, exactly.
    let signs: Vec<Sign> = poly
        .iter()
        .map(|(v, _)| {
            sign(&Lower {
                f: *f,
                g: *g,
                v: *v,
            })
        })
        .collect();
    let points: Vec<[Iv; 2]> = poly.iter().map(|(v, _)| v.enclose()).collect();
    if signs.iter().all(|s| *s == Sign::Zero) {
        // Coplanar here: one plane over the whole overlap.
        return integral(&points, f, base);
    }
    let rise = |p: [Iv; 2]| f.height(p, base).sub(g.height(p, base));
    let below = split(&points, &signs, Sign::Negative, &rise);
    let above = split(&points, &signs, Sign::Positive, &rise);
    integral(&below, f, base).add(integral(&above, g, base))
}

/// The part of a convex polygon where the sign is `keep` or zero, cut
/// where the plane difference changes sign.
fn split(
    points: &[[Iv; 2]],
    signs: &[Sign],
    keep: Sign,
    rise: &dyn Fn([Iv; 2]) -> Iv,
) -> Vec<[Iv; 2]> {
    let n = points.len();
    let mut out = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let (si, sj) = (signs[i], signs[j]);
        if si == keep || si == Sign::Zero {
            out.push(points[i]);
        }
        let opposite = matches!(
            (si, sj),
            (Sign::Negative, Sign::Positive) | (Sign::Positive, Sign::Negative)
        );
        if opposite {
            // Where along the edge the difference is zero: in (0, 1).
            let (ri, rj) = (rise(points[i]), rise(points[j]));
            let t = ri.div(ri.sub(rj)).within(0.0, 1.0);
            let at = |k: usize| points[i][k].add(t.mul(points[j][k].sub(points[i][k])));
            out.push([at(0), at(1)]);
        }
    }
    out
}

/// The integral of the face's height above the base over a convex
/// polygon: a fan of triangles, each its area times its mean height.
fn integral(points: &[[Iv; 2]], face: &Face, base: f64) -> Iv {
    if points.len() < 3 {
        return Iv::point(0.0);
    }
    let heights: Vec<Iv> = points.iter().map(|p| face.height(*p, base)).collect();
    let mut total = Iv::point(0.0);
    for i in 1..points.len() - 1 {
        let (a, b, c) = (points[0], points[i], points[i + 1]);
        let area = b[0]
            .sub(a[0])
            .mul(c[1].sub(a[1]))
            .sub(b[1].sub(a[1]).mul(c[0].sub(a[0])))
            .mul(Iv::point(0.5));
        let mean = heights[0]
            .add(heights[i])
            .add(heights[i + 1])
            .div(Iv::point(3.0));
        total = total.add(area.mul(mean));
    }
    total
}

/// A line through two input points, directed from `a` to `b`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Line {
    a: Point2,
    b: Point2,
}

/// A vertex of an overlap: an input corner, or where two input lines
/// cross.
#[derive(Debug, Clone, Copy)]
enum Vertex {
    Input(Point2),
    Cross(Line, Line),
}

impl Vertex {
    /// Homogeneous coordinates `(X, Y, W)` in the arithmetic `T`: the
    /// point is `(X / W, Y / W)`, with `W` the lines' cross product for a
    /// crossing.
    fn homogeneous<T: Arith>(&self) -> [T; 3] {
        let f = T::from_f64;
        match *self {
            Vertex::Input(p) => [f(p.x), f(p.y), f(1.0)],
            Vertex::Cross(l1, l2) => {
                let d1 = [f(l1.b.x).sub(&f(l1.a.x)), f(l1.b.y).sub(&f(l1.a.y))];
                let d2 = [f(l2.b.x).sub(&f(l2.a.x)), f(l2.b.y).sub(&f(l2.a.y))];
                let w = d1[0].mul(&d2[1]).sub(&d1[1].mul(&d2[0]));
                let e = [f(l2.a.x).sub(&f(l1.a.x)), f(l2.a.y).sub(&f(l1.a.y))];
                let n = e[0].mul(&d2[1]).sub(&e[1].mul(&d2[0]));
                [
                    f(l1.a.x).mul(&w).add(&d1[0].mul(&n)),
                    f(l1.a.y).mul(&w).add(&d1[1].mul(&n)),
                    w,
                ]
            }
        }
    }

    fn enclose(&self) -> [Iv; 2] {
        match *self {
            Vertex::Input(p) => [Iv::point(p.x), Iv::point(p.y)],
            Vertex::Cross(l1, l2) => {
                let d1 = [
                    Iv::point(l1.b.x).sub(Iv::point(l1.a.x)),
                    Iv::point(l1.b.y).sub(Iv::point(l1.a.y)),
                ];
                let d2 = [
                    Iv::point(l2.b.x).sub(Iv::point(l2.a.x)),
                    Iv::point(l2.b.y).sub(Iv::point(l2.a.y)),
                ];
                let w = d1[0].mul(d2[1]).sub(d1[1].mul(d2[0]));
                let e = [
                    Iv::point(l2.a.x).sub(Iv::point(l1.a.x)),
                    Iv::point(l2.a.y).sub(Iv::point(l1.a.y)),
                ];
                let t = e[0].mul(d2[1]).sub(e[1].mul(d2[0])).div(w);
                [
                    Iv::point(l1.a.x).add(d1[0].mul(t)),
                    Iv::point(l1.a.y).add(d1[1].mul(t)),
                ]
            }
        }
    }
}

/// Sutherland-Hodgman against the closed left side of `line`, keeping
/// the line under each edge so every new vertex is a crossing of two
/// input lines.
fn clip(poly: &[(Vertex, Line)], line: Line) -> Vec<(Vertex, Line)> {
    let n = poly.len();
    let signs: Vec<Sign> = poly
        .iter()
        .map(|(v, _)| sign(&SideOf { line, v: *v }))
        .collect();
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..n {
        let (v, edge) = poly[i];
        let (si, sj) = (signs[i], signs[(i + 1) % n]);
        if si != Sign::Negative {
            // From a vertex on the line to the next one out, the kept
            // boundary runs along the line to where it comes back in.
            let along = si == Sign::Zero && sj == Sign::Negative;
            out.push((v, if along { line } else { edge }));
        }
        match (si, sj) {
            (Sign::Positive, Sign::Negative) => out.push((Vertex::Cross(edge, line), line)),
            (Sign::Negative, Sign::Positive) => out.push((Vertex::Cross(edge, line), edge)),
            _ => {}
        }
    }
    out
}

/// The inputs are finite, so the exact tier always decides.
fn sign<E: SignExpr>(e: &E) -> Sign {
    certify(e).unwrap_or(Sign::Zero)
}

struct Orient {
    a: Point2,
    b: Point2,
    c: Point2,
}

impl SignExpr for Orient {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        SideOf {
            line: Line {
                a: self.a,
                b: self.b,
            },
            v: Vertex::Input(self.c),
        }
        .sign_in::<T>()
    }
}

/// Which side of `line` a vertex lies on: positive to the left.
struct SideOf {
    line: Line,
    v: Vertex,
}

impl SignExpr for SideOf {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let [x, y, w] = self.v.homogeneous::<T>();
        let (a, b) = (self.line.a, self.line.b);
        let (dx, dy) = (f(b.x).sub(&f(a.x)), f(b.y).sub(&f(a.y)));
        let rx = x.sub(&f(a.x).mul(&w));
        let ry = y.sub(&f(a.y).mul(&w));
        let s = dx.mul(&ry).sub(&dy.mul(&rx)).sign()?;
        Some(times(s, w.sign()?))
    }
}

/// The sign of `h_f - h_g` at a vertex: negative where `f`'s plane is the
/// lower. Both shadows are counter-clockwise, so both planes' normals
/// point up and multiplying through by them keeps the sign.
struct Lower {
    f: Face,
    g: Face,
    v: Vertex,
}

impl SignExpr for Lower {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        let f = T::from_f64;
        let [x, y, w] = self.v.homogeneous::<T>();
        let normal = |face: &Face| {
            let d = |i: usize, k: usize| f(face.p[i][k]).sub(&f(face.p[0][k]));
            let (u, v) = ([d(1, 0), d(1, 1), d(1, 2)], [d(2, 0), d(2, 1), d(2, 2)]);
            [
                u[1].mul(&v[2]).sub(&u[2].mul(&v[1])),
                u[2].mul(&v[0]).sub(&u[0].mul(&v[2])),
                u[0].mul(&v[1]).sub(&u[1].mul(&v[0])),
            ]
        };
        let (nf, ng) = (normal(&self.f), normal(&self.g));
        // W n_z h(x, y) = W n_z p0.z - n_x (X - W p0.x) - n_y (Y - W p0.y).
        let lift = |face: &Face, n: &[T; 3]| {
            let p = face.p[0];
            w.mul(&n[2])
                .mul(&f(p.z))
                .sub(&n[0].mul(&x.sub(&w.mul(&f(p.x)))))
                .sub(&n[1].mul(&y.sub(&w.mul(&f(p.y)))))
        };
        let difference = lift(&self.f, &nf)
            .mul(&ng[2])
            .sub(&lift(&self.g, &ng).mul(&nf[2]));
        Some(times(difference.sign()?, w.sign()?))
    }
}

fn times(a: Sign, b: Sign) -> Sign {
    match (a, b) {
        (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
        (x, y) if x == y => Sign::Positive,
        _ => Sign::Negative,
    }
}

/// An outward-rounded interval: every operation rounds to nearest and
/// steps each bound one float outward, which covers the rounding error.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Iv {
    lo: f64,
    hi: f64,
}

impl Iv {
    fn point(v: f64) -> Self {
        Self { lo: v, hi: v }
    }

    fn new(lo: f64, hi: f64) -> Self {
        Self { lo, hi }
    }

    fn outward(lo: f64, hi: f64) -> Self {
        if lo.is_nan() || hi.is_nan() {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        Self::new(lo.next_down(), hi.next_up())
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
        let lo = p.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = p.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::outward(lo, hi)
    }

    /// Quotient; the whole line if the divisor may be zero.
    fn div(self, o: Self) -> Self {
        if o.lo <= 0.0 && o.hi >= 0.0 {
            return Self::new(f64::NEG_INFINITY, f64::INFINITY);
        }
        let q = [
            self.lo / o.lo,
            self.lo / o.hi,
            self.hi / o.lo,
            self.hi / o.hi,
        ];
        let lo = q.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = q.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::outward(lo, hi)
    }

    /// Intersected with `[low, high]`, where the value is known to lie.
    fn within(self, low: f64, high: f64) -> Self {
        Self::new(self.lo.max(low), self.hi.min(high))
    }
}
