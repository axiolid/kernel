//! Exact sections of a torus by a plane or sphere (ADR 0076).
//!
//! # Identity
//!
//! With `P(u, v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z`:
//!
//! - a plane `n . (p - o) = 0` gives
//!   `(R + r cos v)(n.X cos u + n.Y sin u) = -(n . (O - o) + r n.Z sin v)`;
//! - a sphere `|p - c|^2 = rho^2`, with `d = O - c`, gives
//!   `2 (R + r cos v)(d.X cos u + d.Y sin u)
//!    = rho^2 - |d|^2 - R^2 - r^2 - 2 R r cos v - 2 r d.Z sin v`.
//!
//! Both are `A(v) cos u + B(v) sin u = C(v)` with `A`, `B`, `C` of degree one
//! in `v`, solved for `u` wherever `E = A^2 + B^2 - C^2 >= 0`.
//!
//! # What is exact
//!
//! Coefficients are dyadic in the operands' doubles; under `t = tan(v/2)`,
//! `E` and the polynomial `B^2 C^2 - A^2 E` (zero where a solution crosses
//! `u = pi`, where the returned angle wraps) are integer polynomials whose
//! roots are isolated exactly. Spans are split at both, so each piece is a
//! continuous pcurve and its existence is never a rounding question.
//!
//! # Scope
//!
//! Ring tori (`R > r`) against planes and spheres off the torus axis. On
//! the axis the closed-form coaxial path applies; cylinders, cones and
//! other tori give quartics in `u` and are refused by name.

use axiolid_core::Interval;
use axiolid_curve::{AngleGraph2, Branch, Curve3, TorusCarrier, TorusSection3};
use axiolid_exact::{Arith, Dyadic, IntPoly};
use axiolid_guarantees::Sign;
use axiolid_surface::{Surface, Torus};

use crate::exact_surface_intersection::{
    Derivation, ExactIntersectionCurve, ExactIntersectionRefusal,
};
use crate::ruled_section::{
    angle_roots, at_pi, d, d3, dot, frame_spans, in_half_angle, int, is_zero, pmul, psub, rounded,
    sign_at, spans, surface_frame, zero, ETrig,
};

fn trig1(constant: Dyadic, cos: Dyadic, sin: Dyadic) -> ETrig {
    [constant, cos, sin, zero(), zero()]
}

/// The exact section of a torus by a plane or sphere, when the pair is in
/// scope; `Ok(None)` when it is not.
pub(crate) fn torus_section(
    first: &Surface,
    second: &Surface,
) -> Result<Option<ExactIntersectionCurve>, ExactIntersectionRefusal> {
    let (torus, other) = match (first, second) {
        (Surface::Torus(t), other @ (Surface::Plane(_) | Surface::Sphere(_)))
        | (other @ (Surface::Plane(_) | Surface::Sphere(_)), Surface::Torus(t)) => (t, other),
        _ => return Ok(None),
    };
    for surface in [first, second] {
        if let Some(frame) = surface_frame(surface) {
            frame_spans(frame)?;
        }
    }
    let Torus {
        frame,
        major_radius,
        minor_radius,
    } = *torus;
    // A spindle or horn torus crosses its own axis: not this graph.
    if major_radius.is_nan() || major_radius <= minor_radius.abs() || minor_radius <= 0.0 {
        return Ok(None);
    }
    let (o, x, y, z) = (d3(frame.origin)?, d3(frame.x)?, d3(frame.y)?, d3(frame.z)?);
    let (big, small) = (d(major_radius)?, d(minor_radius)?);

    let (a, b, c) = match other {
        Surface::Plane(p) => {
            let (po, n) = (d3(p.frame.origin)?, d3(p.frame.z)?);
            let (nx, ny, nz) = (dot(&n, &x), dot(&n, &y), dot(&n, &z));
            let k0 = dot(&n, &o).sub(&dot(&n, &po));
            (
                trig1(big.mul(&nx), small.mul(&nx), zero()),
                trig1(big.mul(&ny), small.mul(&ny), zero()),
                trig1(k0.neg(), zero(), small.mul(&nz).neg()),
            )
        }
        Surface::Sphere(sphere) => {
            let (centre, rho) = (d3(sphere.frame.origin)?, d(sphere.radius)?);
            let dd = [
                o[0].sub(&centre[0]),
                o[1].sub(&centre[1]),
                o[2].sub(&centre[2]),
            ];
            let (dx, dy, dz) = (dot(&dd, &x), dot(&dd, &y), dot(&dd, &z));
            let two = int(2);
            (
                trig1(two.mul(&big).mul(&dx), two.mul(&small).mul(&dx), zero()),
                trig1(two.mul(&big).mul(&dy), two.mul(&small).mul(&dy), zero()),
                trig1(
                    rho.square()
                        .sub(&dot(&dd, &dd))
                        .sub(&big.square())
                        .sub(&small.square()),
                    two.mul(&big).mul(&small).neg(),
                    two.mul(&small).mul(&dz).neg(),
                ),
            )
        }
        _ => return Ok(None),
    };
    if is_zero(&a) && is_zero(&b) {
        // On the torus axis: circles, from the coaxial closed form.
        return Ok(None);
    }

    let (ta, tb, tc) = (in_half_angle(&a), in_half_angle(&b), in_half_angle(&c));
    // E (1 + t^2)^4 and (B^2 C^2 - A^2 E)(1 + t^2)^8.
    let rr = crate::ruled_section::padd(&pmul(&ta, &ta), &pmul(&tb, &tb));
    let e = psub(&rr, &pmul(&tc, &tc));
    let cut = psub(
        &pmul(&pmul(&tb, &tb), &pmul(&tc, &tc)),
        &pmul(&pmul(&ta, &ta), &e),
    );
    let e_poly = IntPoly::from_dyadic(&e);
    if e_poly.is_zero() {
        return Err(ExactIntersectionRefusal::NotRegularCurve);
    }
    let (ap, bp, cp) = (at_pi(&a), at_pi(&b), at_pi(&c));
    let e_pi = ap.square().add(&bp.square()).sub(&cp.square());
    let cut_pi = bp.square().mul(&cp.square()).sub(&ap.square().mul(&e_pi));
    let vanishes = e_pi.sign() == Some(Sign::Zero) || cut_pi.sign() == Some(Sign::Zero);
    let positive = |t: &Dyadic| Some(sign_at(&e_poly, t) == Sign::Positive);
    let found = spans(
        &[e.clone(), cut],
        &positive,
        e_pi.sign() == Some(Sign::Positive),
        vanishes,
    );
    if found.is_empty() {
        return Err(if angle_roots(&e).is_empty() {
            ExactIntersectionRefusal::Disjoint
        } else {
            ExactIntersectionRefusal::NotRegularCurve
        });
    }

    let carrier = TorusCarrier {
        frame,
        major_radius,
        minor_radius,
    };
    let mut branches = Vec::new();
    let mut spans_out = Vec::new();
    for (start, end) in found {
        for branch in [Branch::Plus, Branch::Minus] {
            branches.push(Curve3::TorusSection(TorusSection3 {
                torus: carrier,
                graph: AngleGraph2 {
                    a: rounded(&a),
                    b: rounded(&b),
                    c: rounded(&c),
                    branch,
                },
            }));
            spans_out.push(Some(Interval::new(start, end)));
        }
    }
    Ok(Some(ExactIntersectionCurve::with_spans(
        branches,
        spans_out,
        Derivation::TorusAngleSection,
    )))
}
