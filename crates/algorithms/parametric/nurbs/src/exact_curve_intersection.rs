//! Exact intersections of analytic curves with elementary surfaces and
//! with each other (#119).
//!
//! # Method
//!
//! A line `o + t d` or a conic `o + x a cos(theta) + y b sin(theta)` is
//! written as a polynomial vector over one parameter: `t` for a line, the
//! half-angle `w = tan(theta / 2)` for a conic, using
//! `cos = (1 - w^2) / (1 + w^2)` and `sin = 2 w / (1 + w^2)`. Substituting
//! that into the implicit equation of the other operand and clearing the
//! positive denominator gives an integer polynomial in one variable, whose
//! real roots are isolated exactly by Sturm sequences (`axiolid-exact`).
//!
//! Every decision is an exact sign on the given `f64` input:
//!
//! - **containment**: the polynomial is identically zero, so the curve
//!   lies in the other operand;
//! - **count and order**: exact root isolation;
//! - **tangency**: a root's multiplicity, decided by exact signs of
//!   derivatives at the root;
//! - **the conic point at `theta = pi`**, where `w` is infinite: the
//!   polynomial loses degree, and the lost degree is that point's
//!   multiplicity.
//!
//! # Which surface
//!
//! Implicit equations are derived from the same parametrisations
//! `axiolid-evaluate` uses (`o + x r cos u + y r sin u + z v` for a
//! cylinder, and so on), with frame axes taken exactly as given, not
//! assumed orthonormal. Local coordinates come from Cramer's rule, so the
//! equations stay polynomial in the input. A cone's slope is the `f64`
//! value `tan(semi_angle)`, the value the evaluator also uses; the answer
//! is exact for that cone. A cone point counts only on the nappe the
//! evaluator covers (`radius + v * slope >= 0`).
//!
//! Only returned parameters are exact. Points are rounded, for output.

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_curve::{Curve2, Curve3};
use axiolid_evaluate::evaluate3;
use axiolid_exact::{Arith, Dyadic, IntPoly, RealRoot};
use axiolid_guarantees::Sign;
use axiolid_surface::Surface;

/// Why an exact curve intersection was not computed.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ExactCurveRefusal {
    /// The curve is not a line, circle or ellipse. B-spline curves go
    /// through the certified numeric path instead.
    UnsupportedCurve,
    /// The surface is a B-spline surface or an unknown kind.
    UnsupportedSurface,
    /// An input coordinate, radius or angle is NaN or infinite.
    NonFinite,
    /// A frame is singular, a direction is zero, or a radius is not
    /// positive: there is no curve or surface to intersect.
    Degenerate,
    /// The curve lies in the cone's quadric but crosses its apex, so only
    /// part of it is on the modelled nappe: the overlap is a ray, which
    /// is neither a finite point set nor the whole curve.
    PartialOverlap,
}

/// Where on the first curve an intersection lies, exactly.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ExactCurveParameter {
    /// A line's own parameter `t` (point `origin + t * direction`).
    Line(RealRoot),
    /// A conic's half-angle parameter `w = tan(theta / 2)`.
    HalfAngle(RealRoot),
    /// A conic's point at `theta = pi`, where `w` is infinite.
    Antipode,
}

impl ExactCurveParameter {
    /// The curve parameter as a double, for output: `t` for a line, the
    /// angle `theta` in `[0, 2 pi)` for a conic.
    #[must_use]
    pub fn approx(&self) -> f64 {
        match self {
            Self::Line(t) => t.approx(),
            Self::HalfAngle(w) => {
                let theta = 2.0 * w.approx().atan();
                if theta < 0.0 {
                    theta + std::f64::consts::TAU
                } else {
                    theta
                }
            }
            Self::Antipode => std::f64::consts::PI,
        }
    }
}

/// One intersection point.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ExactCurveHit {
    /// Its parameter on the first curve, exactly.
    pub parameter: ExactCurveParameter,
    /// Its multiplicity: 1 for a transverse crossing, 2 or more where the
    /// curve touches the other operand.
    pub multiplicity: usize,
    /// The point, rounded from the curve at the parameter (output only).
    /// Planar results have `z = 0`.
    pub point: Point3,
}

impl ExactCurveHit {
    /// Whether the curve touches rather than crosses here.
    #[must_use]
    pub fn is_tangent(&self) -> bool {
        self.multiplicity >= 2
    }
}

/// The exact intersection of a curve with another operand.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ExactCurveIntersection {
    /// The whole curve lies in the other operand.
    Contained,
    /// Finitely many points, ordered along the first curve: by `t` for a
    /// line, by `theta` in `[0, 2 pi)` for a conic. Empty when they miss.
    Points(Vec<ExactCurveHit>),
}

/// The exact intersection of an analytic curve with an elementary surface.
///
/// Curves: line, circle, ellipse. Surfaces: plane, cylinder, elliptical
/// cylinder, cone, sphere, torus.
///
/// # Errors
///
/// [`ExactCurveRefusal`] names the unsupported or malformed operand.
pub fn exact_curve_surface_intersection(
    curve: &Curve3,
    surface: &Surface,
) -> Result<ExactCurveIntersection, ExactCurveRefusal> {
    let param = Param::of(curve)?;
    let locus = surface_locus(surface, &param)?;
    let result = solve(curve, &param, &locus);
    if matches!(result, ExactCurveIntersection::Contained) {
        if let Some(condition) = &locus.condition {
            if !nonnegative_everywhere(condition) {
                return Err(ExactCurveRefusal::PartialOverlap);
            }
        }
    }
    Ok(result)
}

/// Whether an integer polynomial is `>= 0` on the whole real line (and so
/// also at a conic's antipode, its formal leading term).
///
/// A sign change needs a root of odd multiplicity; without one the sign
/// away from roots is the leading coefficient's.
fn nonnegative_everywhere(c: &DPoly) -> bool {
    let int = c.to_int();
    if int.is_zero() {
        return true;
    }
    let odd = int
        .real_roots()
        .iter()
        .any(|root| multiplicity(root, &int) % 2 == 1);
    !odd && c.lead_sign() != Sign::Negative
}

/// The exact intersection of two analytic space curves.
///
/// Parameters are those of `first`; swap the arguments for `second`'s.
///
/// # Errors
///
/// [`ExactCurveRefusal`] names the unsupported or malformed operand.
pub fn exact_curve_curve_intersection3(
    first: &Curve3,
    second: &Curve3,
) -> Result<ExactCurveIntersection, ExactCurveRefusal> {
    let param = Param::of(first)?;
    let locus = curve_locus(second, &param)?;
    Ok(solve(first, &param, &locus))
}

/// The exact intersection of two analytic plane curves.
///
/// Parameters are those of `first`; points carry `z = 0`.
///
/// # Errors
///
/// [`ExactCurveRefusal`] names the unsupported or malformed operand.
pub fn exact_curve_curve_intersection2(
    first: &Curve2,
    second: &Curve2,
) -> Result<ExactCurveIntersection, ExactCurveRefusal> {
    exact_curve_curve_intersection3(&lift(first)?, &lift(second)?)
}

// --- polynomials over dyadics ----------------------------------------------

/// A polynomial with dyadic coefficients, lowest degree first.
#[derive(Debug, Clone)]
struct DPoly(Vec<Dyadic>);

impl DPoly {
    fn constant(c: Dyadic) -> Self {
        Self(vec![c])
    }

    fn add(&self, other: &Self) -> Self {
        let n = self.0.len().max(other.0.len());
        Self(
            (0..n)
                .map(|i| match (self.0.get(i), other.0.get(i)) {
                    (Some(a), Some(b)) => a.add(b),
                    (Some(a), None) | (None, Some(a)) => a.clone(),
                    (None, None) => Dyadic::zero(),
                })
                .collect(),
        )
    }

    fn scale(&self, c: &Dyadic) -> Self {
        Self(self.0.iter().map(|a| a.mul(c)).collect())
    }

    fn sub(&self, other: &Self) -> Self {
        self.add(&other.scale(&dy(-1.0)))
    }

    fn mul(&self, other: &Self) -> Self {
        if self.0.is_empty() || other.0.is_empty() {
            return Self(Vec::new());
        }
        let mut out = vec![Dyadic::zero(); self.0.len() + other.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in other.0.iter().enumerate() {
                out[i + j] = out[i + j].add(&a.mul(b));
            }
        }
        Self(out)
    }

    fn square(&self) -> Self {
        self.mul(self)
    }

    fn to_int(&self) -> IntPoly {
        IntPoly::from_dyadic(&self.0)
    }

    /// Sign of the highest non-zero coefficient (zero for the zero
    /// polynomial).
    fn lead_sign(&self) -> Sign {
        self.0
            .iter()
            .rev()
            .filter_map(Arith::sign)
            .find(|s| *s != Sign::Zero)
            .unwrap_or(Sign::Zero)
    }
}

fn dy(value: f64) -> Dyadic {
    Dyadic::from_f64(value)
}

type V3 = [Dyadic; 3];

fn v3(v: Vec3) -> V3 {
    [dy(v.x), dy(v.y), dy(v.z)]
}

fn dot(a: &V3, b: &V3) -> Dyadic {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn cross(a: &V3, b: &V3) -> V3 {
    [
        a[1].mul(&b[2]).sub(&a[2].mul(&b[1])),
        a[2].mul(&b[0]).sub(&a[0].mul(&b[2])),
        a[0].mul(&b[1]).sub(&a[1].mul(&b[0])),
    ]
}

fn is_zero(v: &Dyadic) -> bool {
    v.sign() == Some(Sign::Zero)
}

fn finite(values: &[f64]) -> Result<(), ExactCurveRefusal> {
    if values.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(ExactCurveRefusal::NonFinite)
    }
}

fn positive(values: &[f64]) -> Result<(), ExactCurveRefusal> {
    finite(values)?;
    if values.iter().all(|v| *v > 0.0) {
        Ok(())
    } else {
        Err(ExactCurveRefusal::Degenerate)
    }
}

fn frame_finite(f: &Frame3) -> Result<(), ExactCurveRefusal> {
    finite(&[
        f.origin.x, f.origin.y, f.origin.z, f.x.x, f.x.y, f.x.z, f.y.x, f.y.y, f.y.z, f.z.x, f.z.y,
        f.z.z,
    ])
}

// --- the first curve as a polynomial vector --------------------------------

/// `X(s) = n(s) / w(s)` with `w > 0` for every real `s`.
struct Param {
    n: [DPoly; 3],
    w: DPoly,
    /// Degree of `w`: 0 for a line, 2 for a conic (half-angle form).
    w_degree: usize,
}

impl Param {
    fn of(curve: &Curve3) -> Result<Self, ExactCurveRefusal> {
        match curve {
            Curve3::Line(l) => {
                finite(&[
                    l.origin.x,
                    l.origin.y,
                    l.origin.z,
                    l.direction.x,
                    l.direction.y,
                    l.direction.z,
                ])?;
                let (o, d) = (v3(l.origin), v3(l.direction));
                if d.iter().all(is_zero) {
                    return Err(ExactCurveRefusal::Degenerate);
                }
                Ok(Self {
                    n: [0, 1, 2].map(|i| DPoly(vec![o[i].clone(), d[i].clone()])),
                    w: DPoly::constant(dy(1.0)),
                    w_degree: 0,
                })
            }
            Curve3::Circle(c) => Self::conic(&c.frame, c.radius, c.radius),
            Curve3::Ellipse(e) => Self::conic(&e.frame, e.semi_axis_x, e.semi_axis_y),
            _ => Err(ExactCurveRefusal::UnsupportedCurve),
        }
    }

    /// `o + x a cos + y b sin` with the half-angle substitution, times
    /// `1 + w^2`: `n = o (1 + w^2) + x a (1 - w^2) + y b (2 w)`.
    fn conic(frame: &Frame3, a: f64, b: f64) -> Result<Self, ExactCurveRefusal> {
        frame_finite(frame)?;
        positive(&[a, b])?;
        let (o, x, y) = (v3(frame.origin), v3(frame.x), v3(frame.y));
        if cross(&x, &y).iter().all(is_zero) {
            return Err(ExactCurveRefusal::Degenerate);
        }
        let (a, b) = (dy(a), dy(b));
        let n = [0, 1, 2].map(|i| {
            let xa = x[i].mul(&a);
            DPoly(vec![
                o[i].add(&xa),
                dy(2.0).mul(&y[i]).mul(&b),
                o[i].sub(&xa),
            ])
        });
        Ok(Self {
            n,
            w: DPoly(vec![dy(1.0), Dyadic::zero(), dy(1.0)]),
            w_degree: 2,
        })
    }

    /// The homogenised linear form `u . (X - o)`, times `w`.
    fn linear(&self, u: &V3, o: &V3) -> DPoly {
        let mut out = self.w.scale(&dot(u, o).neg());
        for (n, c) in self.n.iter().zip(u) {
            out = out.add(&n.scale(c));
        }
        out
    }
}

// --- the other operand as equations ----------------------------------------

/// Polynomial equations `p = 0` with their degree `k` in `X`, plus an
/// optional side condition `c >= 0` (the cone's nappe).
struct Locus {
    equations: Vec<(DPoly, usize)>,
    condition: Option<DPoly>,
}

/// Local coordinates by Cramer's rule, times the frame determinant `det`:
/// `alpha = det(P, y, z)`, `beta = det(x, P, z)`, `gamma = det(x, y, P)`.
struct Local {
    alpha: DPoly,
    beta: DPoly,
    gamma: DPoly,
    det: Dyadic,
}

fn local(frame: (&V3, &V3, &V3, &V3), param: &Param) -> Result<Local, ExactCurveRefusal> {
    let (o, x, y, z) = frame;
    let det = dot(x, &cross(y, z));
    if is_zero(&det) {
        return Err(ExactCurveRefusal::Degenerate);
    }
    Ok(Local {
        alpha: param.linear(&cross(y, z), o),
        beta: param.linear(&cross(z, x), o),
        gamma: param.linear(&cross(x, y), o),
        det,
    })
}

fn frame_of(f: &Frame3) -> Result<(V3, V3, V3, V3), ExactCurveRefusal> {
    frame_finite(f)?;
    Ok((v3(f.origin), v3(f.x), v3(f.y), v3(f.z)))
}

fn surface_locus(surface: &Surface, param: &Param) -> Result<Locus, ExactCurveRefusal> {
    let one = |p: DPoly, k: usize| Locus {
        equations: vec![(p, k)],
        condition: None,
    };
    // (det * w)^2, the scale every radius term carries.
    let dw2 = |det: &Dyadic| param.w.scale(det).square();
    match surface {
        Surface::Plane(p) => {
            let (o, x, y, _) = frame_of(&p.frame)?;
            let normal = cross(&x, &y);
            if normal.iter().all(is_zero) {
                return Err(ExactCurveRefusal::Degenerate);
            }
            Ok(one(param.linear(&normal, &o), 1))
        }
        Surface::Cylinder(c) => {
            positive(&[c.radius])?;
            let f = frame_of(&c.frame)?;
            let l = local((&f.0, &f.1, &f.2, &f.3), param)?;
            let r2 = dy(c.radius).square();
            let p = l
                .alpha
                .square()
                .add(&l.beta.square())
                .sub(&dw2(&l.det).scale(&r2));
            Ok(one(p, 2))
        }
        Surface::EllipticalCylinder(c) => {
            positive(&[c.semi_axis_x, c.semi_axis_y])?;
            let f = frame_of(&c.frame)?;
            let l = local((&f.0, &f.1, &f.2, &f.3), param)?;
            let (a2, b2) = (dy(c.semi_axis_x).square(), dy(c.semi_axis_y).square());
            let p = l
                .alpha
                .square()
                .scale(&b2)
                .add(&l.beta.square().scale(&a2))
                .sub(&dw2(&l.det).scale(&a2.mul(&b2)));
            Ok(one(p, 2))
        }
        Surface::Sphere(s) => {
            positive(&[s.radius])?;
            let f = frame_of(&s.frame)?;
            let l = local((&f.0, &f.1, &f.2, &f.3), param)?;
            let r2 = dy(s.radius).square();
            let p = l
                .alpha
                .square()
                .add(&l.beta.square())
                .add(&l.gamma.square())
                .sub(&dw2(&l.det).scale(&r2));
            Ok(one(p, 2))
        }
        Surface::Cone(c) => {
            let slope = c.semi_angle.tan();
            finite(&[c.radius, c.semi_angle, slope])?;
            let f = frame_of(&c.frame)?;
            let l = local((&f.0, &f.1, &f.2, &f.3), param)?;
            // Local radius times det: R det w + slope gamma.
            let reach = param
                .w
                .scale(&dy(c.radius).mul(&l.det))
                .add(&l.gamma.scale(&dy(slope)));
            let p = l.alpha.square().add(&l.beta.square()).sub(&reach.square());
            // R + v slope >= 0 with v = gamma / (det w): the sign of
            // `reach` times the sign of det (w is positive).
            let condition = if l.det.sign() == Some(Sign::Negative) {
                reach.scale(&dy(-1.0))
            } else {
                reach
            };
            Ok(Locus {
                equations: vec![(p, 2)],
                condition: Some(condition),
            })
        }
        Surface::Torus(t) => {
            positive(&[t.major_radius, t.minor_radius])?;
            let f = frame_of(&t.frame)?;
            let l = local((&f.0, &f.1, &f.2, &f.3), param)?;
            let (big2, small2) = (dy(t.major_radius).square(), dy(t.minor_radius).square());
            let planar = l.alpha.square().add(&l.beta.square());
            let dw2 = dw2(&l.det);
            // (|L|^2 + R^2 - r^2)^2 = 4 R^2 (Lx^2 + Ly^2), times det^4 w^4.
            let s = planar
                .add(&l.gamma.square())
                .add(&dw2.scale(&big2.sub(&small2)));
            let p = s.square().sub(&dw2.mul(&planar).scale(&dy(4.0).mul(&big2)));
            Ok(one(p, 4))
        }
        _ => Err(ExactCurveRefusal::UnsupportedSurface),
    }
}

fn curve_locus(curve: &Curve3, param: &Param) -> Result<Locus, ExactCurveRefusal> {
    match curve {
        Curve3::Line(l) => {
            finite(&[
                l.origin.x,
                l.origin.y,
                l.origin.z,
                l.direction.x,
                l.direction.y,
                l.direction.z,
            ])?;
            let (o, d) = (v3(l.origin), v3(l.direction));
            if d.iter().all(is_zero) {
                return Err(ExactCurveRefusal::Degenerate);
            }
            // (X - o) x d = 0, one linear equation per component.
            let z = Dyadic::zero;
            let rows = [
                [z(), d[2].clone(), d[1].neg()],
                [d[2].neg(), z(), d[0].clone()],
                [d[1].clone(), d[0].neg(), z()],
            ];
            Ok(Locus {
                equations: rows.iter().map(|u| (param.linear(u, &o), 1)).collect(),
                condition: None,
            })
        }
        Curve3::Circle(c) => conic_locus(&c.frame, c.radius, c.radius, param),
        Curve3::Ellipse(e) => conic_locus(&e.frame, e.semi_axis_x, e.semi_axis_y, param),
        _ => Err(ExactCurveRefusal::UnsupportedCurve),
    }
}

/// A conic is its plane (`gamma = 0` with `z = x cross y`) and, within it,
/// `b^2 alpha^2 + a^2 beta^2 = a^2 b^2 det^2`.
fn conic_locus(frame: &Frame3, a: f64, b: f64, param: &Param) -> Result<Locus, ExactCurveRefusal> {
    positive(&[a, b])?;
    let (o, x, y, _) = frame_of(frame)?;
    let z = cross(&x, &y);
    let l = local((&o, &x, &y, &z), param)?;
    let (a2, b2) = (dy(a).square(), dy(b).square());
    let ring = l
        .alpha
        .square()
        .scale(&b2)
        .add(&l.beta.square().scale(&a2))
        .sub(&param.w.scale(&l.det).square().scale(&a2.mul(&b2)));
    Ok(Locus {
        equations: vec![(l.gamma, 1), (ring, 2)],
        condition: None,
    })
}

fn lift(curve: &Curve2) -> Result<Curve3, ExactCurveRefusal> {
    use axiolid_curve::{Circle3, Ellipse3, Line3};
    let frame = |f: &axiolid_core::Frame2| Frame3 {
        origin: Point3::new(f.origin.x, f.origin.y, 0.0),
        x: Vec3::new(f.x.x, f.x.y, 0.0),
        y: Vec3::new(f.y.x, f.y.y, 0.0),
        z: Vec3::Z,
    };
    Ok(match curve {
        Curve2::Line(l) => Curve3::Line(Line3 {
            origin: Point3::new(l.origin.x, l.origin.y, 0.0),
            direction: Vec3::new(l.direction.x, l.direction.y, 0.0),
        }),
        Curve2::Circle(c) => Curve3::Circle(Circle3 {
            frame: frame(&c.frame),
            radius: c.radius,
        }),
        Curve2::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            frame: frame(&e.frame),
            semi_axis_x: e.semi_axis_x,
            semi_axis_y: e.semi_axis_y,
        }),
        _ => return Err(ExactCurveRefusal::UnsupportedCurve),
    })
}

// --- solving -----------------------------------------------------------------

/// Multiplicity of `root` in `g`: the order of the first derivative that
/// does not vanish there.
fn multiplicity(root: &RealRoot, g: &IntPoly) -> usize {
    let mut order = 1;
    let mut q = g.derivative();
    while !q.is_zero() && root.sign_of(&q) == Sign::Zero {
        order += 1;
        q = q.derivative();
    }
    order
}

fn solve(curve: &Curve3, param: &Param, locus: &Locus) -> ExactCurveIntersection {
    // Each equation becomes an integer polynomial whose FORMAL degree is
    // k * deg(w): its coefficient of s^(k deg w) is the equation at the
    // conic's point at infinity (theta = pi).
    let polys: Vec<(IntPoly, usize)> = locus
        .equations
        .iter()
        .map(|(p, k)| (p.to_int(), k * param.w_degree))
        .filter(|(p, _)| !p.is_zero())
        .collect();
    if polys.is_empty() {
        return ExactCurveIntersection::Contained;
    }
    // Common roots of all equations are the roots of their gcd.
    let g = polys
        .iter()
        .skip(1)
        .fold(polys[0].0.clone(), |acc, (p, _)| acc.gcd(p));
    let condition_int = locus.condition.as_ref().map(DPoly::to_int);
    let allowed = |root: &RealRoot| {
        condition_int
            .as_ref()
            .is_none_or(|c| c.is_zero() || root.sign_of(c) != Sign::Negative)
    };

    let is_line = param.w_degree == 0;
    let mut hits: Vec<(i8, ExactCurveHit)> = Vec::new();
    for root in g.real_roots() {
        if !allowed(&root) {
            continue;
        }
        let multiplicity = multiplicity(&root, &g);
        // Order along the conic by theta in [0, 2 pi): w >= 0 first, then
        // the antipode, then w < 0.
        let band = if is_line || root.cmp_dyadic(&Dyadic::zero()) != Sign::Negative {
            0
        } else {
            2
        };
        let parameter = if is_line {
            ExactCurveParameter::Line(root)
        } else {
            ExactCurveParameter::HalfAngle(root)
        };
        hits.push((band, hit(curve, parameter, multiplicity)));
    }
    if !is_line {
        // A common zero at infinity: every equation lost degree. Its
        // multiplicity is the smallest loss.
        let lost = polys
            .iter()
            .map(|(p, formal)| formal - p.degree().unwrap_or(0))
            .min()
            .unwrap_or(0);
        let on_nappe = locus.condition.as_ref().is_none_or(|c| {
            // The condition's formal degree is deg(w) = 2; its leading
            // coefficient is its value at the antipode.
            c.0.get(2)
                .is_none_or(|lead| lead.sign() != Some(Sign::Negative))
        });
        if lost >= 1 && on_nappe {
            hits.push((1, hit(curve, ExactCurveParameter::Antipode, lost)));
        }
    }
    // Roots come in increasing order; a stable sort by band keeps it.
    hits.sort_by_key(|(band, _)| *band);
    ExactCurveIntersection::Points(hits.into_iter().map(|(_, h)| h).collect())
}

fn hit(curve: &Curve3, parameter: ExactCurveParameter, multiplicity: usize) -> ExactCurveHit {
    let point =
        evaluate3(curve, parameter.approx()).unwrap_or(Point3::new(f64::NAN, f64::NAN, f64::NAN));
    ExactCurveHit {
        parameter,
        multiplicity,
        point,
    }
}
