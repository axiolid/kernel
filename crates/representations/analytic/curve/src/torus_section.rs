//! Where a plane or sphere cuts a torus, read in the torus's own parameters
//! (ADR 0076).
//!
//! A torus is not ruled, but at a fixed tube angle `v` its points form a
//! circle about the axis, and a plane or sphere meets that circle where
//!
//! ```text
//! A(v) cos(u) + B(v) sin(u) = C(v),
//! ```
//!
//! with `A`, `B`, `C` trigonometric in `v`. So the section is `u` as a
//! function of `v`, in closed form, wherever `E = A^2 + B^2 - C^2 >= 0`:
//!
//! ```text
//! cos u = (A C - s B sqrt E) / (A^2 + B^2),
//! sin u = (B C + s A sqrt E) / (A^2 + B^2),     s = +1 or -1.
//! ```
//!
//! [`AngleGraph2`] holds that as a pcurve; [`TorusSection3`] the same curve
//! in space. The angle is returned in `(-pi, pi]`; a piece is built so it
//! never crosses `u = pi`, where that range wraps.

use axiolid_core::{Frame3, Point3, Scalar, Vec3};

use crate::quadric_section::{Branch, Trig2};

/// The solution `u(t)` of `a(t) cos u + b(t) sin u = c(t)` chosen by
/// `branch`, as a graph over its second coordinate: the point at `t` is
/// `(u(t), t)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AngleGraph2 {
    /// Coefficient of `cos u`.
    pub a: Trig2,
    /// Coefficient of `sin u`.
    pub b: Trig2,
    /// Right-hand side.
    pub c: Trig2,
    /// Which of the two solutions.
    pub branch: Branch,
}

impl AngleGraph2 {
    /// The angle at `t`, in `(-pi, pi]`, or `None` where no solution exists.
    #[must_use]
    pub fn angle(&self, t: Scalar) -> Option<Scalar> {
        let (a, b, c) = (self.a.value(t), self.b.value(t), self.c.value(t));
        let rr = a * a + b * b;
        if rr == 0.0 {
            return None;
        }
        let e = rr - c * c;
        // An exact zero at a branch end rounds either way; the span says
        // the point is on the curve.
        let size =
            |t: &Trig2| t.constant.abs() + t.cos.abs() + t.sin.abs() + t.cos2.abs() + t.sin2.abs();
        let scale = size(&self.a).powi(2) + size(&self.b).powi(2) + size(&self.c).powi(2);
        if e < -1e-12 * scale {
            return None;
        }
        let root = e.max(0.0).sqrt();
        let s = self.branch.sign();
        let u = (b * c + s * a * root).atan2(a * c - s * b * root);
        u.is_finite().then_some(u)
    }

    /// `du/dt`, by differentiating `a cos u + b sin u - c = 0` implicitly;
    /// `None` at a branch end, where the tangent is parallel to `u`.
    #[must_use]
    pub fn slope(&self, t: Scalar) -> Option<Scalar> {
        let u = self.angle(t)?;
        let (s, c) = u.sin_cos();
        let f_t = self.a.derivative(t) * c + self.b.derivative(t) * s - self.c.derivative(t);
        let f_u = -self.a.value(t) * s + self.b.value(t) * c;
        let slope = -f_t / f_u;
        slope.is_finite().then_some(slope)
    }

    /// `d2u/dt2`.
    #[must_use]
    pub fn bend(&self, t: Scalar) -> Option<Scalar> {
        let u = self.angle(t)?;
        let p = self.slope(t)?;
        let (s, c) = u.sin_cos();
        let (a, b) = (self.a.value(t), self.b.value(t));
        let (da, db) = (self.a.derivative(t), self.b.derivative(t));
        let f_u = -a * s + b * c;
        let f_tt = self.a.second(t) * c + self.b.second(t) * s - self.c.second(t);
        let f_tu = -da * s + db * c;
        let f_uu = -a * c - b * s;
        let bend = -(f_tt + 2.0 * f_tu * p + f_uu * p * p) / f_u;
        bend.is_finite().then_some(bend)
    }

    /// Whether every coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.a.is_finite() && self.b.is_finite() && self.c.is_finite()
    }
}

/// A torus, as the curve needs it, in `axiolid_surface::Torus`'s
/// parameterisation: `O + (R + r cos v)(cos u X + sin u Y) + r sin v Z`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorusCarrier {
    /// Local frame; `z` is the axis.
    pub frame: Frame3,
    /// Distance from the axis to the tube centre.
    pub major_radius: Scalar,
    /// Tube radius.
    pub minor_radius: Scalar,
}

impl TorusCarrier {
    /// The torus point at `(u, v)`.
    #[must_use]
    pub fn point(&self, u: Scalar, v: Scalar) -> Point3 {
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let ring = self.major_radius + self.minor_radius * cv;
        self.frame.origin
            + self.frame.x * (ring * cu)
            + self.frame.y * (ring * su)
            + self.frame.z * (self.minor_radius * sv)
    }

    /// `(dP/du, dP/dv)`.
    #[must_use]
    pub fn partials(&self, u: Scalar, v: Scalar) -> (Vec3, Vec3) {
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let ring = self.major_radius + self.minor_radius * cv;
        let r = self.minor_radius;
        (
            self.frame.x * (-ring * su) + self.frame.y * (ring * cu),
            self.frame.x * (-r * sv * cu) + self.frame.y * (-r * sv * su) + self.frame.z * (r * cv),
        )
    }
}

/// A plane's or sphere's section of a torus, in space: the torus point at
/// `(graph.angle(t), t)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorusSection3 {
    /// The torus the curve lies on.
    pub torus: TorusCarrier,
    /// The curve in the torus's parameters.
    pub graph: AngleGraph2,
}

impl TorusSection3 {
    /// Point at `t`.
    #[must_use]
    pub fn point(&self, t: Scalar) -> Option<Point3> {
        Some(self.torus.point(self.graph.angle(t)?, t))
    }

    /// Tangent at `t`: `P_u du/dt + P_v`.
    #[must_use]
    pub fn tangent(&self, t: Scalar) -> Option<Vec3> {
        let u = self.graph.angle(t)?;
        let (pu, pv) = self.torus.partials(u, t);
        Some(pu * self.graph.slope(t)? + pv)
    }

    /// Second derivative at `t`: `P_uu u'^2 + 2 P_uv u' + P_vv + P_u u''`.
    #[must_use]
    pub fn bend(&self, t: Scalar) -> Option<Vec3> {
        let u = self.graph.angle(t)?;
        let p = self.graph.slope(t)?;
        let q = self.graph.bend(t)?;
        let (su, cu) = u.sin_cos();
        let (sv, cv) = t.sin_cos();
        let (x, y, z) = (self.torus.frame.x, self.torus.frame.y, self.torus.frame.z);
        let r = self.torus.minor_radius;
        let ring = self.torus.major_radius + r * cv;
        let p_u = x * (-ring * su) + y * (ring * cu);
        let p_uu = x * (-ring * cu) + y * (-ring * su);
        let p_uv = x * (r * sv * su) + y * (-r * sv * cu);
        let p_vv = x * (-r * cv * cu) + y * (-r * cv * su) + z * (-r * sv);
        Some(p_uu * (p * p) + p_uv * (2.0 * p) + p_vv + p_u * q)
    }
}
