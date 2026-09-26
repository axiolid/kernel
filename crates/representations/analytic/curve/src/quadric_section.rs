//! Where a quadric cuts a cylinder or cone, read in the ruled surface's own
//! parameters (ADR 0076).
//!
//! A cylinder, elliptical cylinder or cone is RULED: at a fixed angle `u`
//! its points run along a straight line in `v`,
//!
//! ```text
//! P(u, v) = O + (rx + s v) cos(u) X + (ry + s v) sin(u) Y + v Z,
//! ```
//!
//! so substituting it into any quadric's equation gives, for each `u`, a
//! quadratic in `v`:
//!
//! ```text
//! a(u) v^2 + b(u) v + c(u) = 0,
//! ```
//!
//! with `a`, `b`, `c` trigonometric polynomials of degree at most two. The
//! section curve is the graph of one root of that quadratic over the angles
//! where its discriminant is not negative. [`QuadraticGraph2`] holds that
//! graph as a pcurve and [`RuledSection3`] the same curve in space. Both are
//! exact: nothing is sampled or fitted, and a plane's cut ([`Sinusoid2`]'s
//! case, `a = 0`) is the degenerate member of the family.
//!
//! [`Sinusoid2`]: crate::Sinusoid2

use axiolid_core::{Frame3, Point3, Scalar, Vec3};

/// `constant + cos cos(t) + sin sin(t) + cos2 cos(2t) + sin2 sin(2t)`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Trig2 {
    /// Constant term.
    pub constant: Scalar,
    /// Coefficient of `cos(t)`.
    pub cos: Scalar,
    /// Coefficient of `sin(t)`.
    pub sin: Scalar,
    /// Coefficient of `cos(2t)`.
    pub cos2: Scalar,
    /// Coefficient of `sin(2t)`.
    pub sin2: Scalar,
}

impl Trig2 {
    /// Value at `t`.
    #[must_use]
    pub fn value(&self, t: Scalar) -> Scalar {
        let (s, c) = t.sin_cos();
        let (s2, c2) = (2.0 * t).sin_cos();
        self.constant + self.cos * c + self.sin * s + self.cos2 * c2 + self.sin2 * s2
    }

    /// First derivative at `t`.
    #[must_use]
    pub fn derivative(&self, t: Scalar) -> Scalar {
        let (s, c) = t.sin_cos();
        let (s2, c2) = (2.0 * t).sin_cos();
        -self.cos * s + self.sin * c - 2.0 * self.cos2 * s2 + 2.0 * self.sin2 * c2
    }

    /// Second derivative at `t`.
    #[must_use]
    pub fn second(&self, t: Scalar) -> Scalar {
        let (s, c) = t.sin_cos();
        let (s2, c2) = (2.0 * t).sin_cos();
        -self.cos * c - self.sin * s - 4.0 * self.cos2 * c2 - 4.0 * self.sin2 * s2
    }

    /// Whether every coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.constant.is_finite()
            && self.cos.is_finite()
            && self.sin.is_finite()
            && self.cos2.is_finite()
            && self.sin2.is_finite()
    }
}

/// Which root of the quadratic a graph follows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Branch {
    /// `v = (-b + sqrt(b^2 - 4ac)) / 2a`.
    Plus,
    /// `v = (-b - sqrt(b^2 - 4ac)) / 2a`.
    Minus,
}

impl Branch {
    /// `+1` or `-1`.
    #[must_use]
    pub fn sign(self) -> Scalar {
        match self {
            Self::Plus => 1.0,
            Self::Minus => -1.0,
        }
    }
}

/// The graph of one root of `a(t) v^2 + b(t) v + c(t) = 0`, parameterised by
/// its first coordinate: the point at `t` is `(t, v(t))`.
///
/// Where `a` vanishes the quadratic is linear and the branch that stays
/// finite is `v = -c / b`; the evaluation below is the rationalised form
/// `v = 2c / (-b - sign sqrt(D))`, which is finite there and free of the
/// cancellation `-b + sqrt(D)` suffers when `4ac` is small.
///
/// The graph is only defined where `D = b^2 - 4ac >= 0`. Which spans those
/// are is decided when the graph is constructed (exactly, from the
/// operands' own numbers) and carried by the edge interval that uses it;
/// evaluating outside them is refused, not extrapolated.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticGraph2 {
    /// Coefficient of `v^2`.
    pub a: Trig2,
    /// Coefficient of `v`.
    pub b: Trig2,
    /// Constant coefficient.
    pub c: Trig2,
    /// Which root.
    pub branch: Branch,
}

impl QuadraticGraph2 {
    /// The discriminant `b^2 - 4ac` at `t`.
    #[must_use]
    pub fn discriminant(&self, t: Scalar) -> Scalar {
        let b = self.b.value(t);
        b * b - 4.0 * self.a.value(t) * self.c.value(t)
    }

    /// Height at `t`, or `None` where the root does not exist or diverges.
    #[must_use]
    pub fn height(&self, t: Scalar) -> Option<Scalar> {
        let (a, b, c) = (self.a.value(t), self.b.value(t), self.c.value(t));
        let d = b * b - 4.0 * a * c;
        // A discriminant a few ulps below zero at a branch end is rounding
        // of an exact zero; the caller's span says the point is on the
        // curve, so it is read as zero rather than refused. At the end both
        // terms are themselves near zero, so the bound is taken from the
        // coefficients' sizes, not from the values there.
        let size =
            |t: &Trig2| t.constant.abs() + t.cos.abs() + t.sin.abs() + t.cos2.abs() + t.sin2.abs();
        let scale = size(&self.b).powi(2) + 4.0 * size(&self.a) * size(&self.c);
        if d < -1e-12 * scale {
            return None;
        }
        let root = d.max(0.0).sqrt();
        let denominator = -b - self.branch.sign() * root;
        if denominator != 0.0 {
            let v = 2.0 * c / denominator;
            if v.is_finite() {
                return Some(v);
            }
        }
        // `-b - sign sqrt(D)` vanishes only where this branch's root is the
        // other formula's: `(-b + sign sqrt(D)) / 2a`.
        if a != 0.0 {
            let v = (-b + self.branch.sign() * root) / (2.0 * a);
            return v.is_finite().then_some(v);
        }
        None
    }

    /// `dv/dt` at `t`, by differentiating the quadratic implicitly; `None`
    /// where the height is undefined or the tangent is vertical (a branch
    /// end, where `2 a v + b = 0`).
    #[must_use]
    pub fn slope(&self, t: Scalar) -> Option<Scalar> {
        let v = self.height(t)?;
        let f_t = self.a.derivative(t) * v * v + self.b.derivative(t) * v + self.c.derivative(t);
        let f_v = 2.0 * self.a.value(t) * v + self.b.value(t);
        let slope = -f_t / f_v;
        slope.is_finite().then_some(slope)
    }

    /// `d2v/dt2` at `t`.
    #[must_use]
    pub fn bend(&self, t: Scalar) -> Option<Scalar> {
        let v = self.height(t)?;
        let p = self.slope(t)?;
        let f_v = 2.0 * self.a.value(t) * v + self.b.value(t);
        let f_tt = self.a.second(t) * v * v + self.b.second(t) * v + self.c.second(t);
        let f_tv = 2.0 * self.a.derivative(t) * v + self.b.derivative(t);
        let f_vv = 2.0 * self.a.value(t);
        let bend = -(f_tt + 2.0 * f_tv * p + f_vv * p * p) / f_v;
        bend.is_finite().then_some(bend)
    }

    /// Whether every coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.a.is_finite() && self.b.is_finite() && self.c.is_finite()
    }
}

/// A ruled carrier surface, as the curve needs it: a cylinder
/// (`rx = ry`, `slope = 0`), an elliptical cylinder (`slope = 0`) or a cone
/// (`rx = ry`, `slope = tan(semi-angle)`), in the same parameterisation as
/// the matching `axiolid_surface` family.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RuledCarrier {
    /// Local frame; `z` is the axis.
    pub frame: Frame3,
    /// Radius along local `x` at `v = 0`.
    pub x_radius: Scalar,
    /// Radius along local `y` at `v = 0`.
    pub y_radius: Scalar,
    /// How fast both radii grow with `v`.
    pub slope: Scalar,
}

impl RuledCarrier {
    /// The carrier's point at `(u, v)`.
    #[must_use]
    pub fn point(&self, u: Scalar, v: Scalar) -> Point3 {
        let (s, c) = u.sin_cos();
        self.frame.origin
            + self.frame.x * ((self.x_radius + self.slope * v) * c)
            + self.frame.y * ((self.y_radius + self.slope * v) * s)
            + self.frame.z * v
    }

    /// `(dP/du, dP/dv)` at `(u, v)`.
    #[must_use]
    pub fn partials(&self, u: Scalar, v: Scalar) -> (Vec3, Vec3) {
        let (s, c) = u.sin_cos();
        (
            self.frame.x * (-(self.x_radius + self.slope * v) * s)
                + self.frame.y * ((self.y_radius + self.slope * v) * c),
            self.frame.x * (self.slope * c) + self.frame.y * (self.slope * s) + self.frame.z,
        )
    }

    /// Whether every number is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.frame.origin.is_finite()
            && self.frame.x.is_finite()
            && self.frame.y.is_finite()
            && self.frame.z.is_finite()
            && self.x_radius.is_finite()
            && self.y_radius.is_finite()
            && self.slope.is_finite()
    }
}

/// A quadric's section of a ruled carrier, in space: the point at `t` is the
/// carrier's point at `(t, graph.height(t))`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RuledSection3 {
    /// The ruled surface the curve lies on.
    pub carrier: RuledCarrier,
    /// The curve in the carrier's parameters.
    pub graph: QuadraticGraph2,
}

impl RuledSection3 {
    /// Point at `t`, or `None` outside the graph's spans.
    #[must_use]
    pub fn point(&self, t: Scalar) -> Option<Point3> {
        Some(self.carrier.point(t, self.graph.height(t)?))
    }

    /// Tangent at `t`: `P_u + P_v dv/dt`.
    #[must_use]
    pub fn tangent(&self, t: Scalar) -> Option<Vec3> {
        let v = self.graph.height(t)?;
        let (pu, pv) = self.carrier.partials(t, v);
        Some(pu + pv * self.graph.slope(t)?)
    }

    /// Second derivative at `t`: `P_uu + 2 P_uv v' + P_v v''` (the carrier
    /// is linear in `v`, so `P_vv = 0`).
    #[must_use]
    pub fn bend(&self, t: Scalar) -> Option<Vec3> {
        let v = self.graph.height(t)?;
        let slope = self.graph.slope(t)?;
        let bend = self.graph.bend(t)?;
        let (s, c) = t.sin_cos();
        let k = &self.carrier;
        let p_uu = k.frame.x * (-(k.x_radius + k.slope * v) * c)
            + k.frame.y * (-(k.y_radius + k.slope * v) * s);
        let p_uv = k.frame.x * (-k.slope * s) + k.frame.y * (k.slope * c);
        let (_, p_v) = k.partials(t, v);
        Some(p_uu + p_uv * (2.0 * slope) + p_v * bend)
    }
}
