//! Exact sections of a cylinder or cone by a quadric (ADR 0076).
//!
//! # Identity
//!
//! A ruled carrier is linear in `v` at every angle `u`:
//! `P(u, v) = A0(u) + v A1(u)`, each `A` of the form
//! `e0 + e_c cos(u) + e_s sin(u)`. A quadric is `Q(p) = p^T M p + 2 q^T p + k`.
//! Substituting,
//!
//! ```text
//! Q(P(u, v)) = a(u) v^2 + b(u) v + c(u),
//! a = A1^T M A1,  b = 2 (A1^T M A0 + q^T A1),  c = A0^T M A0 + 2 q^T A0 + k,
//! ```
//!
//! trigonometric polynomials of degree at most two. The section is where one
//! of the two roots in `v` exists: where `D(u) = b^2 - 4ac >= 0`.
//!
//! # What is exact
//!
//! Every coefficient is built in dyadic arithmetic from the operands' own
//! doubles, so `a`, `b`, `c` and `D` are the exact trigonometric polynomials
//! of the given surfaces. Under `t = tan(u/2)`, `D (1 + t^2)^4` is an integer
//! polynomial of degree 8 whose real roots are isolated exactly
//! (`axiolid-exact`, Sturm): which spans of angle carry the curve, whether it
//! is one loop, two loops or none, and whether the surfaces only touch are
//! decided without rounding. The stored coefficients and span ends are the
//! only rounded values.
//!
//! # Scope
//!
//! Carriers: cylinder, elliptical cylinder, and a cone cut by a plane.
//! Quadrics: plane, sphere, cylinder, elliptical cylinder. A cone cut by a
//! curved quadric needs its nappe decided along a square-root branch and is
//! refused by name, as are tori and B-spline surfaces.

use axiolid_core::{Interval, Scalar};
use axiolid_curve::{Branch, Curve3, QuadraticGraph2, RuledCarrier, RuledSection3, Trig2};
use axiolid_exact::{Arith, Dyadic, IntPoly};
use axiolid_guarantees::Sign;
use axiolid_surface::Surface;

use crate::exact_surface_intersection::{
    Derivation, ExactIntersectionCurve, ExactIntersectionRefusal,
};

pub(crate) type D3 = [Dyadic; 3];

/// `e0 + ec cos(u) + es sin(u)`, exactly.
#[derive(Clone)]
struct TrigVec {
    e0: D3,
    ec: D3,
    es: D3,
}

/// A degree-2 trigonometric polynomial, exactly: `[1, cos, sin, cos2, sin2]`.
pub(crate) type ETrig = [Dyadic; 5];

/// `p^T M p + 2 q^T p + k`, exactly.
struct Quadric {
    m: [D3; 3],
    q: D3,
    k: Dyadic,
}

pub(crate) fn refuse() -> ExactIntersectionRefusal {
    ExactIntersectionRefusal::DegenerateFrame
}

pub(crate) fn d(value: Scalar) -> Result<Dyadic, ExactIntersectionRefusal> {
    Dyadic::try_from_f64(value).ok_or_else(refuse)
}

pub(crate) fn d3(v: axiolid_core::Vec3) -> Result<D3, ExactIntersectionRefusal> {
    Ok([d(v.x)?, d(v.y)?, d(v.z)?])
}

pub(crate) fn zero() -> Dyadic {
    Dyadic::zero()
}

/// A small integer, exactly (every such value is a double).
pub(crate) fn int(value: i64) -> Dyadic {
    Dyadic::try_from_f64(value as f64).expect("small integers are finite")
}

fn half() -> Dyadic {
    Dyadic::try_from_f64(0.5).expect("finite")
}

pub(crate) fn dot(a: &D3, b: &D3) -> Dyadic {
    a[0].mul(&b[0]).add(&a[1].mul(&b[1])).add(&a[2].mul(&b[2]))
}

fn scale(a: &D3, s: &Dyadic) -> D3 {
    [a[0].mul(s), a[1].mul(s), a[2].mul(s)]
}

fn add3(a: &D3, b: &D3) -> D3 {
    [a[0].add(&b[0]), a[1].add(&b[1]), a[2].add(&b[2])]
}

fn zero3() -> D3 {
    [zero(), zero(), zero()]
}

fn apply(m: &[D3; 3], v: &D3) -> D3 {
    [dot(&m[0], v), dot(&m[1], v), dot(&m[2], v)]
}

/// `x^T M y` for trigonometric-linear `x` and `y`: degree 2.
fn bilinear(m: &[D3; 3], x: &TrigVec, y: &TrigVec) -> ETrig {
    let f = |a: &D3, b: &D3| dot(a, &apply(m, b));
    let cc = f(&x.ec, &y.ec);
    let ss = f(&x.es, &y.es);
    let cs = f(&x.ec, &y.es).add(&f(&x.es, &y.ec));
    // cos^2 = (1 + cos2)/2, sin^2 = (1 - cos2)/2, cos sin = sin2/2.
    [
        f(&x.e0, &y.e0).add(&cc.add(&ss).mul(&half())),
        f(&x.e0, &y.ec).add(&f(&x.ec, &y.e0)),
        f(&x.e0, &y.es).add(&f(&x.es, &y.e0)),
        cc.sub(&ss).mul(&half()),
        cs.mul(&half()),
    ]
}

/// `q^T x` for trigonometric-linear `x`: degree 1.
fn linear(q: &D3, x: &TrigVec) -> ETrig {
    [dot(q, &x.e0), dot(q, &x.ec), dot(q, &x.es), zero(), zero()]
}

fn trig_add(a: &ETrig, b: &ETrig) -> ETrig {
    [
        a[0].add(&b[0]),
        a[1].add(&b[1]),
        a[2].add(&b[2]),
        a[3].add(&b[3]),
        a[4].add(&b[4]),
    ]
}

fn trig_scale(a: &ETrig, s: &Dyadic) -> ETrig {
    [
        a[0].mul(s),
        a[1].mul(s),
        a[2].mul(s),
        a[3].mul(s),
        a[4].mul(s),
    ]
}

pub(crate) fn rounded(a: &ETrig) -> Trig2 {
    Trig2 {
        constant: a[0].to_f64(),
        cos: a[1].to_f64(),
        sin: a[2].to_f64(),
        cos2: a[3].to_f64(),
        sin2: a[4].to_f64(),
    }
}

pub(crate) fn is_zero(a: &ETrig) -> bool {
    a.iter().all(|c| c.sign() == Some(Sign::Zero))
}

/// Polynomial product, lowest degree first.
pub(crate) fn pmul(a: &[Dyadic], b: &[Dyadic]) -> Vec<Dyadic> {
    let mut out = vec![zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = out[i + j].add(&x.mul(y));
        }
    }
    out
}

pub(crate) fn padd(a: &[Dyadic], b: &[Dyadic]) -> Vec<Dyadic> {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| {
            let x = a.get(i).cloned().unwrap_or_else(zero);
            let y = b.get(i).cloned().unwrap_or_else(zero);
            x.add(&y)
        })
        .collect()
}

pub(crate) fn psub(a: &[Dyadic], b: &[Dyadic]) -> Vec<Dyadic> {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| {
            let x = a.get(i).cloned().unwrap_or_else(zero);
            let y = b.get(i).cloned().unwrap_or_else(zero);
            x.sub(&y)
        })
        .collect()
}

/// `T(u) (1 + t^2)^2` as a polynomial in `t = tan(u/2)`, degree 4.
pub(crate) fn in_half_angle(a: &ETrig) -> Vec<Dyadic> {
    // 1 -> 1 + 2t^2 + t^4; cos -> 1 - t^4; sin -> 2t + 2t^3;
    // cos2 -> 1 - 6t^2 + t^4; sin2 -> 4t - 4t^3.
    let basis: [[i64; 5]; 5] = [
        [1, 0, 2, 0, 1],
        [1, 0, 0, 0, -1],
        [0, 2, 0, 2, 0],
        [1, 0, -6, 0, 1],
        [0, 4, 0, -4, 0],
    ];
    let mut out = vec![zero(); 5];
    for (coefficient, row) in a.iter().zip(basis) {
        for (slot, weight) in out.iter_mut().zip(row) {
            *slot = slot.add(&coefficient.mul(&int(weight)));
        }
    }
    out
}

/// Value of an exact trigonometric polynomial at `u = pi`, exactly.
pub(crate) fn at_pi(a: &ETrig) -> Dyadic {
    a[0].sub(&a[1]).add(&a[3])
}

/// Whether a frame spans space: every axis finite and the three independent.
/// A degenerate frame describes a degenerate parameterisation, and an
/// implicit equation that ignores the frame would describe a different
/// surface from the one evaluation draws.
pub(crate) fn frame_spans(frame: &axiolid_core::Frame3) -> Result<(), ExactIntersectionRefusal> {
    let (x, y, z) = (d3(frame.x)?, d3(frame.y)?, d3(frame.z)?);
    let cross = [
        x[1].mul(&y[2]).sub(&x[2].mul(&y[1])),
        x[2].mul(&y[0]).sub(&x[0].mul(&y[2])),
        x[0].mul(&y[1]).sub(&x[1].mul(&y[0])),
    ];
    d3(frame.origin)?;
    if dot(&cross, &z).sign() == Some(Sign::Zero) {
        return Err(refuse());
    }
    Ok(())
}

pub(crate) fn surface_frame(surface: &Surface) -> Option<&axiolid_core::Frame3> {
    match surface {
        Surface::Plane(p) => Some(&p.frame),
        Surface::Cylinder(c) => Some(&c.frame),
        Surface::EllipticalCylinder(c) => Some(&c.frame),
        Surface::Cone(c) => Some(&c.frame),
        Surface::Sphere(s) => Some(&s.frame),
        Surface::Torus(t) => Some(&t.frame),
        _ => None,
    }
}

/// The carrier as `A0 + v A1`, with the matching curve value.
fn carrier(
    surface: &Surface,
) -> Result<Option<(TrigVec, TrigVec, RuledCarrier)>, ExactIntersectionRefusal> {
    let (frame, rx, ry, slope) = match surface {
        Surface::Cylinder(c) => (c.frame, c.radius, c.radius, 0.0),
        Surface::EllipticalCylinder(c) => (c.frame, c.semi_axis_x, c.semi_axis_y, 0.0),
        Surface::Cone(c) => (c.frame, c.radius, c.radius, c.semi_angle.tan()),
        _ => return Ok(None),
    };
    let (o, x, y, z) = (d3(frame.origin)?, d3(frame.x)?, d3(frame.y)?, d3(frame.z)?);
    let (rx_e, ry_e, s) = (d(rx)?, d(ry)?, d(slope)?);
    let a0 = TrigVec {
        e0: o,
        ec: scale(&x, &rx_e),
        es: scale(&y, &ry_e),
    };
    let a1 = TrigVec {
        e0: z,
        ec: scale(&x, &s),
        es: scale(&y, &s),
    };
    Ok(Some((
        a0,
        a1,
        RuledCarrier {
            frame,
            x_radius: rx,
            y_radius: ry,
            slope,
        },
    )))
}

fn outer(a: &D3, b: &D3) -> [D3; 3] {
    [
        [a[0].mul(&b[0]), a[0].mul(&b[1]), a[0].mul(&b[2])],
        [a[1].mul(&b[0]), a[1].mul(&b[1]), a[1].mul(&b[2])],
        [a[2].mul(&b[0]), a[2].mul(&b[1]), a[2].mul(&b[2])],
    ]
}

fn madd(a: &[D3; 3], b: &[D3; 3]) -> [D3; 3] {
    [add3(&a[0], &b[0]), add3(&a[1], &b[1]), add3(&a[2], &b[2])]
}

fn mscale(a: &[D3; 3], s: &Dyadic) -> [D3; 3] {
    [scale(&a[0], s), scale(&a[1], s), scale(&a[2], s)]
}

fn identity(s: &Dyadic) -> [D3; 3] {
    [
        [s.clone(), zero(), zero()],
        [zero(), s.clone(), zero()],
        [zero(), zero(), s.clone()],
    ]
}

/// Centre a quadric `(p - o)^T M (p - o) + lin . (p - o) + k0` at `o`.
fn centred(m: [D3; 3], o: &D3, k0: Dyadic) -> Quadric {
    let mo = apply(&m, o);
    Quadric {
        q: scale(&mo, &int(-1)),
        k: dot(o, &mo).add(&k0),
        m,
    }
}

/// The implicit equation of the quadric operand, in the operand frame's
/// axes as given (orthonormal up to their own rounding).
fn quadric(surface: &Surface) -> Result<Option<Quadric>, ExactIntersectionRefusal> {
    Ok(Some(match surface {
        Surface::Plane(p) => {
            let (o, n) = (d3(p.frame.origin)?, d3(p.frame.z)?);
            // Q = n . (p - o): linear, so M = 0 and 2q = n.
            Quadric {
                m: [zero3(), zero3(), zero3()],
                q: scale(&n, &half()),
                k: dot(&n, &o).neg(),
            }
        }
        Surface::Sphere(s) => {
            let (o, r) = (d3(s.frame.origin)?, d(s.radius)?);
            centred(identity(&int(1)), &o, r.square().neg())
        }
        Surface::Cylinder(c) => {
            // |w|^2 |p - o|^2 - ((p - o) . w)^2 - r^2 |w|^2.
            let (o, w, r) = (d3(c.frame.origin)?, d3(c.frame.z)?, d(c.radius)?);
            let ww = dot(&w, &w);
            let m = madd(&identity(&ww), &mscale(&outer(&w, &w), &int(-1)));
            centred(m, &o, r.square().mul(&ww).neg())
        }
        Surface::EllipticalCylinder(c) => {
            // b^2 ((p - o) . X)^2 + a^2 ((p - o) . Y)^2 - a^2 b^2.
            let (o, x, y) = (d3(c.frame.origin)?, d3(c.frame.x)?, d3(c.frame.y)?);
            let (a2, b2) = (d(c.semi_axis_x)?.square(), d(c.semi_axis_y)?.square());
            let m = madd(&mscale(&outer(&x, &x), &b2), &mscale(&outer(&y, &y), &a2));
            centred(m, &o, a2.mul(&b2).neg())
        }
        Surface::Cone(c) => {
            // With d = p - o and h = d . Z: radial^2 = (r + s h)^2, i.e.
            // d^T (I - (1 + s^2) Z Z^T) d - 2 r s Z . d - r^2 = 0. Both
            // nappes; the modelled one is selected by `nappes`.
            let (o, z) = (d3(c.frame.origin)?, d3(c.frame.z)?);
            let (r, sl) = (d(c.radius)?, d(c.semi_angle.tan())?);
            let m = madd(
                &identity(&int(1)),
                &mscale(&outer(&z, &z), &int(1).add(&sl.square()).neg()),
            );
            let mut quadric = centred(m, &o, r.square().neg());
            let rs = r.mul(&sl);
            quadric.q = add3(&quadric.q, &scale(&z, &rs.neg()));
            quadric.k = quadric.k.add(&int(2).mul(&rs).mul(&dot(&z, &o)));
            quadric
        }
        _ => return Ok(None),
    }))
}

/// Conditions `h0(u) + h1(u) v >= 0` a point must meet to lie on the
/// modelled nappe of each cone taking part: the carrier's own radius
/// `r + s v`, and the other cone's `r' + s' Z' . (P - O')`.
fn nappes(
    carrier_surface: &Surface,
    other: &Surface,
    a0: &TrigVec,
    a1: &TrigVec,
) -> Result<Vec<(ETrig, ETrig)>, ExactIntersectionRefusal> {
    let constant = |value: Dyadic| [value, zero(), zero(), zero(), zero()];
    let mut out = Vec::new();
    if let Surface::Cone(c) = carrier_surface {
        out.push((constant(d(c.radius)?), constant(d(c.semi_angle.tan())?)));
    }
    if let Surface::Cone(c) = other {
        let (o, z) = (d3(c.frame.origin)?, d3(c.frame.z)?);
        let (r, sl) = (d(c.radius)?, d(c.semi_angle.tan())?);
        let sz = scale(&z, &sl);
        let mut h0 = linear(&sz, a0);
        h0[0] = h0[0].add(&r).sub(&dot(&sz, &o));
        out.push((h0, linear(&sz, a1)));
    }
    Ok(out)
}

/// Angles in `(-pi, pi)` where an exact half-angle polynomial vanishes, in
/// increasing order, as `(t_lo, t_hi, u)`: the isolating interval in `t` and
/// the root's angle rounded once.
pub(crate) fn angle_roots(poly: &[Dyadic]) -> Vec<(Dyadic, Dyadic, Scalar)> {
    let poly = IntPoly::from_dyadic(poly);
    if poly.is_zero() {
        return Vec::new();
    }
    let width = Dyadic::try_from_f64(2.0_f64.powi(-60)).expect("finite");
    poly.real_roots()
        .into_iter()
        .map(|mut root| {
            root.refine_to_width(&width);
            let (lo, hi) = root.bounds();
            let u = 2.0 * root.approx().atan();
            (lo.clone(), hi.clone(), u)
        })
        .collect()
}

/// A dyadic `t` strictly between two isolated roots.
fn between(left: &(Dyadic, Dyadic, Scalar), right: &(Dyadic, Dyadic, Scalar)) -> Dyadic {
    left.1.add(&right.0).mul(&half())
}

pub(crate) fn sign_at(poly: &IntPoly, t: &Dyadic) -> Sign {
    poly.sign_at(t)
}

/// Spans of angle where `positive` holds, from the roots of the
/// polynomials that can change it and a test at one dyadic point per span.
/// Each span is `[start, end]` with `end > start`, possibly running past
/// `pi`.
///
/// `u = pi` is `t = infinity` under the half-angle map, so a root there is
/// never isolated as a polynomial root: it shows as the polynomial's degree
/// dropping. `vanishes_at_pi` says whether any break is zero there exactly,
/// and then `pi` is a break too. Spans beside it are sampled beyond the
/// outermost finite root.
pub(crate) fn spans(
    breaks: &[Vec<Dyadic>],
    positive: &dyn Fn(&Dyadic) -> Option<bool>,
    positive_at_pi: bool,
    vanishes_at_pi: bool,
) -> Vec<(Scalar, Scalar)> {
    let mut roots: Vec<(Dyadic, Dyadic, Scalar)> =
        breaks.iter().flat_map(|p| angle_roots(p)).collect();
    roots.sort_by(|a, b| a.2.total_cmp(&b.2));
    roots.dedup_by(|a, b| (a.2 - b.2).abs() < 1e-15);
    let one = int(1);
    let mut out = Vec::new();
    if roots.is_empty() {
        let whole = if vanishes_at_pi {
            positive(&zero()) == Some(true)
        } else {
            positive_at_pi
        };
        if whole {
            out.push((-core::f64::consts::PI, core::f64::consts::PI));
        }
        return out;
    }
    for pair in roots.windows(2) {
        if positive(&between(&pair[0], &pair[1])) == Some(true) {
            out.push((pair[0].2, pair[1].2));
        }
    }
    let (first, last) = (&roots[0], &roots[roots.len() - 1]);
    if vanishes_at_pi {
        // (last, pi) and (-pi, first) are separate spans.
        if positive(&last.1.add(&one)) == Some(true) {
            out.push((last.2, core::f64::consts::PI));
        }
        if positive(&first.0.sub(&one)) == Some(true) {
            out.push((-core::f64::consts::PI, first.2));
        }
    } else if positive_at_pi {
        out.push((last.2, first.2 + core::f64::consts::TAU));
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

/// The exact section of a ruled carrier by a quadric, when the pair is in
/// scope; `Ok(None)` when it is not.
pub(crate) fn ruled_section(
    first: &Surface,
    second: &Surface,
) -> Result<Option<ExactIntersectionCurve>, ExactIntersectionRefusal> {
    // A cylinder carries in preference to a cone, so a cone's nappe is
    // decided as the other operand where it can be.
    let pick = |carrier_surface: &Surface, other: &Surface| {
        matches!(
            (carrier_surface, other),
            (
                Surface::Cylinder(_) | Surface::EllipticalCylinder(_),
                Surface::Plane(_)
                    | Surface::Sphere(_)
                    | Surface::Cylinder(_)
                    | Surface::EllipticalCylinder(_)
                    | Surface::Cone(_)
            ) | (
                Surface::Cone(_),
                Surface::Plane(_) | Surface::Sphere(_) | Surface::Cone(_)
            )
        )
    };
    let (carrier_surface, other) = if pick(first, second) {
        (first, second)
    } else if pick(second, first) {
        (second, first)
    } else {
        return Ok(None);
    };
    for surface in [carrier_surface, other] {
        if let Some(frame) = surface_frame(surface) {
            frame_spans(frame)?;
        }
    }
    let Some((a0, a1, curve_carrier)) = carrier(carrier_surface)? else {
        return Ok(None);
    };
    let Some(q) = quadric(other)? else {
        return Ok(None);
    };

    let a = bilinear(&q.m, &a1, &a1);
    let b = trig_scale(
        &trig_add(&bilinear(&q.m, &a1, &a0), &linear(&q.q, &a1)),
        &int(2),
    );
    let c = {
        let quad = bilinear(&q.m, &a0, &a0);
        let lin = trig_scale(&linear(&q.q, &a0), &int(2));
        let mut sum = trig_add(&quad, &lin);
        sum[0] = sum[0].add(&q.k);
        sum
    };
    let graph = |branch: Branch| QuadraticGraph2 {
        a: rounded(&a),
        b: rounded(&b),
        c: rounded(&c),
        branch,
    };
    let piece = |branch: Branch| {
        Curve3::RuledSection(RuledSection3 {
            carrier: curve_carrier,
            graph: graph(branch),
        })
    };

    let (ta, tb, tc) = (in_half_angle(&a), in_half_angle(&b), in_half_angle(&c));
    let mut branches = Vec::new();
    let mut spans_out = Vec::new();

    if is_zero(&a) && matches!(other, Surface::Cone(_)) {
        // A cone whose generators all run parallel to the carrier's rulings
        // in the quadratic term: not a graph this path splits by nappe.
        return Ok(None);
    }
    if is_zero(&a) {
        // Linear in v (a plane): v = -c / b wherever b != 0, and on a cone
        // only on the modelled nappe, r + s v >= 0, i.e. (r b - s c) b >= 0.
        if is_zero(&b) {
            // The plane contains every ruling direction: parallel to a
            // cylinder's axis, handled by the dedicated rulings path.
            return Ok(None);
        }
        let slope = d(curve_carrier.slope)?;
        let radius = d(curve_carrier.x_radius)?;
        let nappe_trig = {
            let mut t = trig_scale(&b, &radius);
            let sc = trig_scale(&c, &slope);
            for (slot, value) in t.iter_mut().zip(sc.iter()) {
                *slot = slot.sub(value);
            }
            t
        };
        // A plane through a cone's apex cuts rulings (or only the apex):
        // lines of constant angle, not a graph over the angle. Left to the
        // closed-form path, which names it.
        if curve_carrier.slope != 0.0 && is_zero(&nappe_trig) {
            return Ok(None);
        }
        let tb_poly = IntPoly::from_dyadic(&tb);
        let tn = in_half_angle(&nappe_trig);
        let tn_poly = IntPoly::from_dyadic(&tn);
        let cone = curve_carrier.slope != 0.0;
        let good = |t: &Dyadic| -> Option<bool> {
            let sb = sign_at(&tb_poly, t);
            if sb == Sign::Zero {
                return None;
            }
            if !cone {
                return Some(true);
            }
            let sn = sign_at(&tn_poly, t);
            Some(sn == Sign::Zero || (sn == sb))
        };
        let at_pi_ok = {
            let sb = at_pi(&b).sign().unwrap_or(Sign::Zero);
            sb != Sign::Zero
                && (!cone || {
                    let sn = at_pi(&nappe_trig).sign().unwrap_or(Sign::Zero);
                    sn == Sign::Zero || sn == sb
                })
        };
        let mut breaks = vec![tb.clone()];
        if cone {
            breaks.push(tn.clone());
        }
        let vanishes = at_pi(&b).sign() == Some(Sign::Zero)
            || (cone && at_pi(&nappe_trig).sign() == Some(Sign::Zero));
        for (start, end) in spans(&breaks, &good, at_pi_ok, vanishes) {
            // Branch whose rationalised root stays finite: sign of b.
            let mid = 0.5 * (start + end);
            let branch = if rounded(&b).value(mid) > 0.0 {
                Branch::Plus
            } else {
                Branch::Minus
            };
            branches.push(piece(branch));
            spans_out.push(Some(Interval::new(start, end)));
        }
    } else {
        // Every polynomial below is its trigonometric value times
        // (1 + t^2)^k with the same k where they are compared, so signs and
        // roots are those of the values themselves.
        // D (1 + t^2)^4 = B^2 - 4 A C.
        let discriminant = psub(&pmul(&tb, &tb), &pmul(&trig_scale_poly(&ta, 4), &tc));
        let d_poly = IntPoly::from_dyadic(&discriminant);
        if d_poly.is_zero() {
            return Err(ExactIntersectionRefusal::NotRegularCurve);
        }
        let a_poly = IntPoly::from_dyadic(&ta);
        // A nappe condition h0 + h1 v >= 0 along v = (-b + sign sqrt D) / 2a
        // has the sign of a * (P + sign Q sqrt D) with
        // P = 2 a h0 - b h1 and Q = h1, all scaled to (1 + t^2)^4.
        let conditions: Vec<(Vec<Dyadic>, Vec<Dyadic>, ETrig, ETrig)> =
            nappes(carrier_surface, other, &a0, &a1)?
                .into_iter()
                .map(|(h0, h1)| {
                    let (h0t, h1t) = (in_half_angle(&h0), in_half_angle(&h1));
                    let p = psub(&trig_scale_poly(&pmul(&ta, &h0t), 2), &pmul(&tb, &h1t));
                    (p, h1t, h0, h1)
                })
                .collect();
        let mut breaks = vec![discriminant.clone()];
        if !is_const(&a) {
            breaks.push(ta.clone());
        }
        for (p, q, _, _) in &conditions {
            breaks.push(psub(&pmul(p, p), &pmul(&pmul(q, q), &discriminant)));
        }
        // Values at u = pi, exactly, for the span through it.
        let (ap, bp, cp) = (at_pi(&a), at_pi(&b), at_pi(&c));
        let d_pi = bp.square().sub(&int(4).mul(&ap).mul(&cp));
        let conditions_pi: Vec<(Dyadic, Dyadic)> = conditions
            .iter()
            .map(|(_, _, h0, h1)| {
                let (h0p, h1p) = (at_pi(h0), at_pi(h1));
                (int(2).mul(&ap).mul(&h0p).sub(&bp.mul(&h1p)), h1p)
            })
            .collect();
        let vanishes = d_pi.sign() == Some(Sign::Zero)
            || ap.sign() == Some(Sign::Zero)
            || conditions_pi
                .iter()
                .any(|(p, q)| p.square().sub(&q.square().mul(&d_pi)).sign() == Some(Sign::Zero));
        // Whether a branch holds at one point, from exact values there.
        let holds =
            |branch: Branch, a_s: Sign, d_v: &Dyadic, cond: &[(Dyadic, Dyadic)]| -> Option<bool> {
                if d_v.sign()? != Sign::Positive || a_s == Sign::Zero {
                    return Some(false);
                }
                for (p, q) in cond {
                    let q = match branch {
                        Branch::Plus => q.clone(),
                        Branch::Minus => q.neg(),
                    };
                    let s = axiolid_exact::sign_root(p, &q, d_v)?;
                    // The condition's sign is sign(a) * s: negative when the
                    // two signs are strictly opposite.
                    let opposite = matches!(
                        (a_s, s),
                        (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
                    );
                    if opposite {
                        return Some(false);
                    }
                }
                Some(true)
            };
        let mut any = false;
        for branch in [Branch::Plus, Branch::Minus] {
            let positive = |t: &Dyadic| {
                let d_v = poly_value(&discriminant, t);
                let a_s = sign_at(&a_poly, t);
                let cond: Vec<(Dyadic, Dyadic)> = conditions
                    .iter()
                    .map(|(p, q, _, _)| (poly_value(p, t), poly_value(q, t)))
                    .collect();
                holds(branch, a_s, &d_v, &cond)
            };
            let at_pi_ok = holds(
                branch,
                ap.sign().unwrap_or(Sign::Zero),
                &d_pi,
                &conditions_pi,
            ) == Some(true);
            for (start, end) in spans(&breaks, &positive, at_pi_ok, vanishes) {
                any = true;
                branches.push(piece(branch));
                spans_out.push(Some(Interval::new(start, end)));
            }
        }
        if !any {
            // No span: apart, touching at isolated points, or meeting only
            // on the nappes the operands do not model.
            return Err(
                if angle_roots(&discriminant).is_empty() && d_pi.sign() != Some(Sign::Positive) {
                    ExactIntersectionRefusal::Disjoint
                } else if conditions.is_empty() {
                    ExactIntersectionRefusal::NotRegularCurve
                } else {
                    ExactIntersectionRefusal::Disjoint
                },
            );
        }
    }

    if branches.is_empty() {
        return Err(ExactIntersectionRefusal::Disjoint);
    }
    Ok(Some(ExactIntersectionCurve::with_spans(
        branches,
        spans_out,
        Derivation::RuledQuadricSection,
    )))
}

/// Exact value of a dyadic-coefficient polynomial at `t`.
fn poly_value(poly: &[Dyadic], t: &Dyadic) -> Dyadic {
    poly.iter()
        .rev()
        .fold(zero(), |acc, coefficient| acc.mul(t).add(coefficient))
}

/// Whether a trigonometric polynomial is constant.
fn is_const(a: &ETrig) -> bool {
    a[1..].iter().all(|c| c.sign() == Some(Sign::Zero))
}

fn trig_scale_poly(p: &[Dyadic], s: i64) -> Vec<Dyadic> {
    p.iter().map(|c| c.mul(&int(s))).collect()
}
