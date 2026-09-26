//! Extrema between points, curves and surfaces (#119, B7): the smallest
//! distance between two pieces of geometry, certified.
//!
//! A piece is a point, a curve over a span, or a surface over a box of its
//! parameters. The search is branch and bound over the pieces' parameter
//! ranges. Each sub-range has a box certain to hold its image:
//! - lines and patches of a plane map affinely;
//! - circles, ellipses, cylinders, cones, spheres and tori have exact
//!   ranges of their harmonics;
//! - B-splines are held by the control hull of the restricted Bezier
//!   pieces (rational weights by interval division);
//! - traced section curves by their certified cells and their carrier's
//!   image of them.
//!
//! The distance between two boxes bounds the pair from below. The distance
//! between two evaluated points bounds it from above, and those points
//! witness it. The search stops when the bounds are within the accuracy
//! asked for. Nothing is sampled on trust: a pair is discarded only when its
//! lower bound exceeds the best distance found.

use axiolid_core::{Interval, Point2, Point3, Scalar, Vec3};
use axiolid_curve::implicit::{bound_simple, partial, Cell};
use axiolid_curve::{Axis, Carrier, Curve3, ImplicitCurve2};
use axiolid_evaluate::evaluate3;
use axiolid_surface::Surface;
use core::f64::consts::{PI, TAU};

use crate::field::carrier_of;

/// One piece of geometry.
#[derive(Debug, Clone, Copy)]
pub enum Piece<'a> {
    /// A point.
    Point(Point3),
    /// A curve over a span of its parameter.
    Curve {
        /// The curve.
        curve: &'a Curve3,
        /// The span.
        span: Interval,
    },
    /// A surface over a box of its parameters.
    Patch {
        /// The surface.
        surface: &'a Surface,
        /// Lower corner of the parameter box.
        lo: Point2,
        /// Upper corner.
        hi: Point2,
    },
}

/// The smallest distance between two pieces, bracketed, with points on
/// each at the upper bound.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extremum {
    /// Certain lower bound on the distance.
    pub lower: Scalar,
    /// Distance between the witnesses: an upper bound.
    pub upper: Scalar,
    /// Point on the first piece.
    pub on_first: Point3,
    /// Point on the second piece.
    pub on_second: Point3,
    /// Its parameters on the first piece (`t`, or `(u, v)`; `NaN` for a
    /// point).
    pub first: Point2,
    /// Its parameters on the second piece.
    pub second: Point2,
}

/// Why an extremum could not be found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ExtremaRefusal {
    /// A curve or surface family this module has no image bounds for.
    Unsupported,
    /// A piece could not be evaluated.
    Evaluation,
    /// The search exceeded its work budget before reaching the accuracy.
    Budget,
}

/// An axis-aligned box.
#[derive(Debug, Clone, Copy)]
struct Aabb {
    lo: Vec3,
    hi: Vec3,
}

impl Aabb {
    fn point(p: Point3) -> Self {
        Self { lo: p, hi: p }
    }

    fn join(self, o: Self) -> Self {
        Self {
            lo: self.lo.min(o.lo),
            hi: self.hi.max(o.hi),
        }
    }

    fn distance(&self, o: &Self) -> Scalar {
        let gap = (self.lo - o.hi).max(o.lo - self.hi).max(Vec3::ZERO);
        gap.length()
    }

    fn diameter(&self) -> Scalar {
        (self.hi - self.lo).length()
    }

    fn widen(self, by: Scalar) -> Self {
        Self {
            lo: self.lo - Vec3::splat(by),
            hi: self.hi + Vec3::splat(by),
        }
    }
}

/// An interval of reals, with the arithmetic the bounds need.
#[derive(Debug, Clone, Copy)]
struct Iv(Scalar, Scalar);

impl Iv {
    fn add(self, o: Iv) -> Iv {
        Iv(self.0 + o.0, self.1 + o.1)
    }

    fn mul(self, o: Iv) -> Iv {
        let p = [self.0 * o.0, self.0 * o.1, self.1 * o.0, self.1 * o.1];
        Iv(
            p.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        )
    }

    fn scale(self, c: Scalar) -> Iv {
        let (a, b) = (self.0 * c, self.1 * c);
        Iv(a.min(b), a.max(b))
    }
}

fn cos_iv(a: Scalar, b: Scalar) -> Iv {
    if b - a >= TAU {
        return Iv(-1.0, 1.0);
    }
    let (ca, cb) = (a.cos(), b.cos());
    let (mut lo, mut hi) = (ca.min(cb), ca.max(cb));
    if ((a / TAU).ceil()) * TAU <= b {
        hi = 1.0;
    }
    if (((a - PI) / TAU).ceil()) * TAU + PI <= b {
        lo = -1.0;
    }
    Iv(lo, hi)
}

fn sin_iv(a: Scalar, b: Scalar) -> Iv {
    cos_iv(a - 0.5 * PI, b - 0.5 * PI)
}

/// `origin + sum_k axis_k * factor_k` over intervals of the factors.
fn combine(origin: Point3, terms: &[(Vec3, Iv)]) -> Aabb {
    let mut lo = origin;
    let mut hi = origin;
    for (axis, factor) in terms {
        for k in 0..3 {
            let c = factor.scale(axis[k]);
            lo[k] += c.0;
            hi[k] += c.1;
        }
    }
    Aabb { lo, hi }
}

/// A box holding the carrier's points over a parameter box.
fn carrier_box(carrier: &Carrier, lo: Point2, hi: Point2) -> Option<Aabb> {
    let (cu, su) = (cos_iv(lo.x, hi.x), sin_iv(lo.x, hi.x));
    let bx = match carrier {
        Carrier::Plane(f) => {
            let mut b = Aabb::point(f.origin + f.x * lo.x + f.y * lo.y);
            for (u, v) in [(hi.x, lo.y), (lo.x, hi.y), (hi.x, hi.y)] {
                b = b.join(Aabb::point(f.origin + f.x * u + f.y * v));
            }
            b
        }
        Carrier::Ruled(k) => {
            let v = Iv(lo.y, hi.y);
            let rx = v.scale(k.slope).add(Iv(k.x_radius, k.x_radius));
            let ry = v.scale(k.slope).add(Iv(k.y_radius, k.y_radius));
            combine(
                k.frame.origin,
                &[
                    (k.frame.x, rx.mul(cu)),
                    (k.frame.y, ry.mul(su)),
                    (k.frame.z, v),
                ],
            )
        }
        Carrier::Sphere { frame, radius } => {
            let (cv, sv) = (cos_iv(lo.y, hi.y), sin_iv(lo.y, hi.y));
            let r = Iv(*radius, *radius);
            combine(
                frame.origin,
                &[
                    (frame.x, r.mul(cv).mul(cu)),
                    (frame.y, r.mul(cv).mul(su)),
                    (frame.z, r.mul(sv)),
                ],
            )
        }
        Carrier::Torus(t) => {
            let (cv, sv) = (cos_iv(lo.y, hi.y), sin_iv(lo.y, hi.y));
            let ring = cv
                .scale(t.minor_radius)
                .add(Iv(t.major_radius, t.major_radius));
            combine(
                t.frame.origin,
                &[
                    (t.frame.x, ring.mul(cu)),
                    (t.frame.y, ring.mul(su)),
                    (t.frame.z, sv.scale(t.minor_radius)),
                ],
            )
        }
        Carrier::Spline(b) => return spline_box(b, lo, hi),
    };
    // Rounding in the sums.
    let size = bx.lo.abs().max(bx.hi.abs()).max_element();
    Some(bx.widen(32.0 * Scalar::EPSILON * (1.0 + size)))
}

/// A box holding a B-spline surface over a parameter box: each coordinate's
/// homogeneous field bounded by its restricted Bernstein coefficients,
/// divided by the weight's bounds.
#[allow(clippy::needless_range_loop)]
fn spline_box(b: &axiolid_curve::BSplineSurface, lo: Point2, hi: Point2) -> Option<Aabb> {
    let fields = crate::spline_field::homogeneous_fields(b)?;
    let cell = Cell { lo, hi };
    let w = bound_simple(&fields[3], &cell);
    if w.lo <= 0.0 {
        return None;
    }
    let mut out = Aabb {
        lo: Vec3::ZERO,
        hi: Vec3::ZERO,
    };
    for k in 0..3 {
        let x = bound_simple(&fields[k], &cell);
        let q = [x.lo / w.lo, x.lo / w.hi, x.hi / w.lo, x.hi / w.hi];
        out.lo[k] = q.iter().copied().fold(Scalar::INFINITY, Scalar::min);
        out.hi[k] = q.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max);
    }
    Some(out)
}

/// The box `[free0, free1] x [...]` of carrier parameters certain to hold
/// an implicit curve's cell for `s` in `[s0, s1]`.
fn cell_hull(
    curve: &ImplicitCurve2,
    index: usize,
    s0: Scalar,
    s1: Scalar,
) -> Option<(Point2, Point2)> {
    let cell = &curve.cells[index];
    let free = |s: Scalar| cell.from + (cell.to - cell.from) * s;
    let (f0, f1) = (free(s0), free(s1));
    let (f_lo, f_hi) = (f0.min(f1), f0.max(f1));
    let place = |fr: Scalar, w: Scalar| match cell.axis {
        Axis::U => Point2::new(fr, w),
        Axis::V => Point2::new(w, fr),
    };
    let w0 = match cell.axis {
        Axis::U => curve.point(index as Scalar + s0)?.y,
        Axis::V => curve.point(index as Scalar + s0)?.x,
    };
    // A bridge is straight: the box of its ends.
    // A bridge lies in its bracket about the straight line.
    if cell.bridge.is_some() {
        let (w_lo, w_hi) = cell.solved_range();
        let (a, b) = (place(f_lo, w_lo), place(f_hi, w_hi));
        return Some((a.min(b), a.max(b)));
    }
    let (d_free, d_solved) = match cell.axis {
        Axis::U => (partial(&curve.field, true), partial(&curve.field, false)),
        Axis::V => (partial(&curve.field, false), partial(&curve.field, true)),
    };
    let whole = {
        let (a, b) = (place(f_lo, cell.low), place(f_hi, cell.high));
        Cell {
            lo: a.min(b),
            hi: a.max(b),
        }
    };
    let fr = bound_simple(&d_free, &whole);
    let so = bound_simple(&d_solved, &whole);
    let (w_lo, w_hi) = if so.straddles_zero() {
        (cell.low, cell.high)
    } else {
        let floor = so.lo.abs().min(so.hi.abs());
        let reach = fr.lo.abs().max(fr.hi.abs()) / floor * (f_hi - f_lo);
        ((w0 - reach).max(cell.low), (w0 + reach).min(cell.high))
    };
    let (a, b) = (place(f_lo, w_lo), place(f_hi, w_hi));
    Some((a.min(b), a.max(b)))
}

/// The range of `a cos t + b sin t` over `[t0, t1]`, exactly.
fn harmonic(a: Scalar, b: Scalar, t0: Scalar, t1: Scalar) -> Iv {
    let r = a.hypot(b);
    if r == 0.0 {
        return Iv(0.0, 0.0);
    }
    let phi = b.atan2(a);
    cos_iv(t0 - phi, t1 - phi).scale(r)
}

/// The range of `n . P` over the carrier's points on a parameter box: the
/// direct bound intersected with the mean-value form
/// `f(c) + [f_u] [-h_u, h_u] + [f_v] [-h_v, h_v]`, which is second order
/// where `n` is normal to the surface -- exactly where a closest point is,
/// so the search there needs boxes only as small as the accuracy's square
/// root.
fn carrier_projection(carrier: &Carrier, n: Vec3, lo: Point2, hi: Point2) -> Option<Iv> {
    let direct = carrier_projection_direct(carrier, n, lo, hi)?;
    let Some((fu, fv)) = carrier_gradient(carrier, n, lo, hi) else {
        return Some(direct);
    };
    let c = (lo + hi) * 0.5;
    let (hu, hv) = (0.5 * (hi.x - lo.x), 0.5 * (hi.y - lo.y));
    let jet = carrier.jet(c.x, c.y);
    let f = n.dot(jet.point);
    let m = fu.scale(1.0).mul(Iv(-hu, hu)).add(fv.mul(Iv(-hv, hv)));
    let slack = 32.0 * Scalar::EPSILON * (1.0 + jet.point.length());
    let mean = Iv(f + m.0 - slack, f + m.1 + slack);
    Some(Iv(direct.0.max(mean.0), direct.1.min(mean.1)))
}

/// Ranges of `n . P_u` and `n . P_v` over a parameter box.
fn carrier_gradient(carrier: &Carrier, n: Vec3, lo: Point2, hi: Point2) -> Option<(Iv, Iv)> {
    let (cu, su) = (cos_iv(lo.x, hi.x), sin_iv(lo.x, hi.x));
    Some(match carrier {
        Carrier::Plane(f) => {
            let (a, b) = (n.dot(f.x), n.dot(f.y));
            (Iv(a, a), Iv(b, b))
        }
        Carrier::Ruled(k) => {
            let f = &k.frame;
            let (a, b, c) = (n.dot(f.x), n.dot(f.y), n.dot(f.z));
            let v = Iv(lo.y, hi.y);
            let rx = v.scale(k.slope).add(Iv(k.x_radius, k.x_radius));
            let ry = v.scale(k.slope).add(Iv(k.y_radius, k.y_radius));
            // P_u = -rx sin u X + ry cos u Y; P_v = s cos u X + s sin u Y + Z.
            let pu = rx.mul(su.scale(-a)).add(ry.mul(cu.scale(b)));
            let pv = harmonic(a * k.slope, b * k.slope, lo.x, hi.x).add(Iv(c, c));
            (pu, pv)
        }
        Carrier::Sphere { frame, radius } => {
            let (a, b, c) = (n.dot(frame.x), n.dot(frame.y), n.dot(frame.z));
            let (cv, sv) = (cos_iv(lo.y, hi.y), sin_iv(lo.y, hi.y));
            // P_u = r cos v (-sin u X + cos u Y); P_v = r (-sin v (cos u X +
            // sin u Y) + cos v Z).
            let pu = cv.mul(harmonic(b, -a, lo.x, hi.x)).scale(*radius);
            let pv = sv
                .mul(harmonic(a, b, lo.x, hi.x))
                .scale(-1.0)
                .add(cv.scale(c))
                .scale(*radius);
            (pu, pv)
        }
        Carrier::Torus(t) => {
            let f = &t.frame;
            let (a, b, c) = (n.dot(f.x), n.dot(f.y), n.dot(f.z));
            let (cv, sv) = (cos_iv(lo.y, hi.y), sin_iv(lo.y, hi.y));
            let r = t.minor_radius;
            let ring = cv.scale(r).add(Iv(t.major_radius, t.major_radius));
            let pu = ring.mul(harmonic(b, -a, lo.x, hi.x));
            let pv = sv
                .mul(harmonic(a, b, lo.x, hi.x))
                .scale(-r)
                .add(cv.scale(c * r));
            (pu, pv)
        }
        Carrier::Spline(_) => return None,
    })
}

/// The direct bound of [`carrier_projection`].
fn carrier_projection_direct(carrier: &Carrier, n: Vec3, lo: Point2, hi: Point2) -> Option<Iv> {
    let pad = |iv: Iv, scale: Scalar| {
        let m = 32.0 * Scalar::EPSILON * (1.0 + scale);
        Iv(iv.0 - m, iv.1 + m)
    };
    Some(match carrier {
        Carrier::Plane(f) => {
            let o = n.dot(f.origin);
            let (a, b) = (n.dot(f.x), n.dot(f.y));
            let (u, v) = (Iv(lo.x, hi.x).scale(a), Iv(lo.y, hi.y).scale(b));
            pad(Iv(o, o).add(u).add(v), o.abs() + a.abs() + b.abs())
        }
        Carrier::Ruled(k) => {
            let f = &k.frame;
            let (a, b, c) = (n.dot(f.x), n.dot(f.y), n.dot(f.z));
            let o = n.dot(f.origin);
            let v = Iv(lo.y, hi.y);
            let iv = if k.slope == 0.0 && k.x_radius == k.y_radius {
                // A cylinder: separable, exact.
                harmonic(a * k.x_radius, b * k.y_radius, lo.x, hi.x).add(v.scale(c))
            } else {
                let rx = v.scale(k.slope).add(Iv(k.x_radius, k.x_radius));
                let ry = v.scale(k.slope).add(Iv(k.y_radius, k.y_radius));
                rx.mul(cos_iv(lo.x, hi.x).scale(a))
                    .add(ry.mul(sin_iv(lo.x, hi.x).scale(b)))
                    .add(v.scale(c))
            };
            pad(
                Iv(o, o).add(iv),
                o.abs() + k.x_radius + k.y_radius + v.1.abs(),
            )
        }
        Carrier::Sphere { frame, radius } => {
            let (a, b, c) = (n.dot(frame.x), n.dot(frame.y), n.dot(frame.z));
            let o = n.dot(frame.origin);
            let h = harmonic(a, b, lo.x, hi.x);
            let iv = cos_iv(lo.y, hi.y)
                .mul(h)
                .add(sin_iv(lo.y, hi.y).scale(c))
                .scale(*radius);
            pad(Iv(o, o).add(iv), o.abs() + radius)
        }
        Carrier::Torus(t) => {
            let f = &t.frame;
            let (a, b, c) = (n.dot(f.x), n.dot(f.y), n.dot(f.z));
            let o = n.dot(f.origin);
            let h = harmonic(a, b, lo.x, hi.x);
            let ring = cos_iv(lo.y, hi.y)
                .scale(t.minor_radius)
                .add(Iv(t.major_radius, t.major_radius));
            let iv = ring
                .mul(h)
                .add(sin_iv(lo.y, hi.y).scale(c * t.minor_radius));
            pad(Iv(o, o).add(iv), o.abs() + t.major_radius + t.minor_radius)
        }
        Carrier::Spline(b) => {
            let bx = spline_box(b, lo, hi)?;
            project_box(&bx, n)
        }
    })
}

/// The range of `n . P` over a box.
fn project_box(b: &Aabb, n: Vec3) -> Iv {
    let mut lo = 0.0;
    let mut hi = 0.0;
    for k in 0..3 {
        let (x, y) = (n[k] * b.lo[k], n[k] * b.hi[k]);
        lo += x.min(y);
        hi += x.max(y);
    }
    Iv(lo, hi)
}

/// The range of `n . C(t)` over `[t0, t1]` of a curve.
fn curve_projection(curve: &Curve3, n: Vec3, t0: Scalar, t1: Scalar) -> Option<Iv> {
    let (a, b) = (t0.min(t1), t0.max(t1));
    match curve {
        Curve3::Line(l) => {
            let (p, q) = (
                n.dot(l.origin + l.direction * a),
                n.dot(l.origin + l.direction * b),
            );
            Some(Iv(p.min(q), p.max(q)))
        }
        // Exact: one harmonic.
        Curve3::Circle(c) => {
            let o = n.dot(c.frame.origin);
            let h = harmonic(
                n.dot(c.frame.x) * c.radius,
                n.dot(c.frame.y) * c.radius,
                a,
                b,
            );
            Some(Iv(o, o).add(h))
        }
        Curve3::Ellipse(e) => {
            let o = n.dot(e.frame.origin);
            let h = harmonic(
                n.dot(e.frame.x) * e.semi_axis_x,
                n.dot(e.frame.y) * e.semi_axis_y,
                a,
                b,
            );
            Some(Iv(o, o).add(h))
        }
        Curve3::ImplicitSection(s) => {
            let mut out: Option<Iv> = None;
            for index in 0..s.curve.cells.len() {
                let (c0, c1) = (index as Scalar, index as Scalar + 1.0);
                let (lo, hi) = (a.max(c0), b.min(c1));
                if hi < lo {
                    continue;
                }
                let (p, q) = cell_hull(&s.curve, index, lo - c0, hi - c0)?;
                let iv = carrier_projection(&s.carrier, n, p, q)?;
                out = Some(out.map_or(iv, |o: Iv| Iv(o.0.min(iv.0), o.1.max(iv.1))));
            }
            out
        }
        _ => curve_box(curve, a, b).map(|bx| project_box(&bx, n)),
    }
}

/// A box holding a curve's image over `[t0, t1]`.
fn curve_box(curve: &Curve3, t0: Scalar, t1: Scalar) -> Option<Aabb> {
    let (a, b) = (t0.min(t1), t0.max(t1));
    match curve {
        Curve3::Line(l) => Some(
            Aabb::point(l.origin + l.direction * a).join(Aabb::point(l.origin + l.direction * b)),
        ),
        Curve3::Circle(c) => Some(combine(
            c.frame.origin,
            &[
                (c.frame.x, cos_iv(a, b).scale(c.radius)),
                (c.frame.y, sin_iv(a, b).scale(c.radius)),
            ],
        )),
        Curve3::Ellipse(e) => Some(combine(
            e.frame.origin,
            &[
                (e.frame.x, cos_iv(a, b).scale(e.semi_axis_x)),
                (e.frame.y, sin_iv(a, b).scale(e.semi_axis_y)),
            ],
        )),
        Curve3::ImplicitSection(s) => {
            let mut out: Option<Aabb> = None;
            for index in 0..s.curve.cells.len() {
                let (c0, c1) = (index as Scalar, index as Scalar + 1.0);
                let (lo, hi) = (a.max(c0), b.min(c1));
                if hi < lo {
                    continue;
                }
                let (p, q) = cell_hull(&s.curve, index, lo - c0, hi - c0)?;
                let bx = carrier_box(&s.carrier, p, q)?;
                out = Some(out.map_or(bx, |o| o.join(bx)));
            }
            out
        }
        Curve3::BSpline(_) => {
            // Its control hull over the span, by the curve's own Bezier
            // pieces: bounded conservatively by the whole control polygon
            // restricted to the span's knot range.
            crate::spline_field::curve_hull(curve, a, b).map(|(lo, hi)| Aabb { lo, hi })
        }
        _ => None,
    }
}

/// A piece reduced to what the search needs: a parameter range, how to
/// bound its image and where its points are.
#[derive(Debug, Clone)]
enum Part {
    Point(Point3),
    Curve {
        curve: Curve3,
        lo: Scalar,
        hi: Scalar,
    },
    Patch {
        surface: Surface,
        carrier: Carrier,
        lo: Point2,
        hi: Point2,
    },
}

impl Part {
    fn of(piece: &Piece<'_>) -> Result<Part, ExtremaRefusal> {
        Ok(match piece {
            Piece::Point(p) => Part::Point(*p),
            Piece::Curve { curve, span } => {
                // The ADR 0076 graphs are bounded through their traced form.
                let curve = match curve {
                    Curve3::RuledSection(_) | Curve3::TorusSection(_) => Curve3::ImplicitSection(
                        crate::implicit_ops::implicit_view(curve, *span)
                            .ok_or(ExtremaRefusal::Unsupported)?,
                    ),
                    other => (*other).clone(),
                };
                let (lo, hi) = match (&curve, piece) {
                    (
                        Curve3::ImplicitSection(s),
                        Piece::Curve {
                            curve: original, ..
                        },
                    ) if !matches!(original, Curve3::ImplicitSection(_)) => (0.0, s.curve.end()),
                    _ => (span.start.min(span.end), span.start.max(span.end)),
                };
                Part::Curve { curve, lo, hi }
            }
            Piece::Patch { surface, lo, hi } => Part::Patch {
                surface: (*surface).clone(),
                carrier: carrier_of(surface).ok_or(ExtremaRefusal::Unsupported)?,
                lo: *lo,
                hi: *hi,
            },
        })
    }
}

/// A sub-range of a part.
#[derive(Debug, Clone, Copy)]
enum Range2 {
    Point,
    Curve(Scalar, Scalar),
    Patch(Point2, Point2),
}

impl Range2 {
    fn whole(part: &Part) -> Range2 {
        match part {
            Part::Point(_) => Range2::Point,
            Part::Curve { lo, hi, .. } => Range2::Curve(*lo, *hi),
            Part::Patch { lo, hi, .. } => Range2::Patch(*lo, *hi),
        }
    }

    fn bounds(&self, part: &Part) -> Option<Aabb> {
        match (self, part) {
            (Range2::Point, Part::Point(p)) => Some(Aabb::point(*p)),
            (Range2::Curve(a, b), Part::Curve { curve, .. }) => curve_box(curve, *a, *b),
            (Range2::Patch(lo, hi), Part::Patch { carrier, .. }) => carrier_box(carrier, *lo, *hi),
            _ => None,
        }
    }

    /// The range of `n . P` over the part's points in this range.
    fn projection(&self, part: &Part, n: Vec3) -> Option<Iv> {
        match (self, part) {
            (Range2::Point, Part::Point(p)) => {
                let x = n.dot(*p);
                Some(Iv(x, x))
            }
            (Range2::Curve(a, b), Part::Curve { curve, .. }) => curve_projection(curve, n, *a, *b),
            (Range2::Patch(lo, hi), Part::Patch { carrier, .. }) => {
                carrier_projection(carrier, n, *lo, *hi)
            }
            _ => None,
        }
    }

    /// The part's point in the middle of the range, and its parameters.
    fn sample(&self, part: &Part) -> Option<(Point3, Point2)> {
        match (self, part) {
            (Range2::Point, Part::Point(p)) => Some((*p, Point2::splat(Scalar::NAN))),
            (Range2::Curve(a, b), Part::Curve { curve, .. }) => {
                let t = 0.5 * (a + b);
                Some((evaluate3(curve, t).ok()?, Point2::new(t, Scalar::NAN)))
            }
            (Range2::Patch(lo, hi), Part::Patch { surface, .. }) => {
                let c = (*lo + *hi) * 0.5;
                Some((
                    axiolid_evaluate::surface::evaluate(surface, c.x, c.y).ok()?,
                    c,
                ))
            }
            _ => None,
        }
    }

    fn split(&self) -> Vec<Range2> {
        match self {
            Range2::Point => vec![Range2::Point],
            Range2::Curve(a, b) => {
                let m = 0.5 * (a + b);
                vec![Range2::Curve(*a, m), Range2::Curve(m, *b)]
            }
            Range2::Patch(lo, hi) => {
                let d = *hi - *lo;
                if d.x >= d.y {
                    let m = 0.5 * (lo.x + hi.x);
                    vec![
                        Range2::Patch(*lo, Point2::new(m, hi.y)),
                        Range2::Patch(Point2::new(m, lo.y), *hi),
                    ]
                } else {
                    let m = 0.5 * (lo.y + hi.y);
                    vec![
                        Range2::Patch(*lo, Point2::new(hi.x, m)),
                        Range2::Patch(Point2::new(lo.x, m), *hi),
                    ]
                }
            }
        }
    }

    fn is_point(&self) -> bool {
        matches!(self, Range2::Point)
    }
}

/// The smallest distance between two pieces, to within `accuracy`.
///
/// # Errors
///
/// A family without image bounds, an evaluation failure, or a search that
/// needs more than its budget of pairs.
pub fn minimum_distance(
    first: &Piece<'_>,
    second: &Piece<'_>,
    accuracy: Scalar,
) -> Result<Extremum, ExtremaRefusal> {
    let (pa, pb) = (Part::of(first)?, Part::of(second)?);
    let (ra, rb) = (Range2::whole(&pa), Range2::whole(&pb));
    let sample = |r: &Range2, p: &Part| r.sample(p).ok_or(ExtremaRefusal::Evaluation);
    let (sa, sb) = (sample(&ra, &pa)?, sample(&rb, &pb)?);
    let mut best = Extremum {
        lower: 0.0,
        upper: (sa.0 - sb.0).length(),
        on_first: sa.0,
        on_second: sb.0,
        first: sa.1,
        second: sb.1,
    };
    let bounds = |r: &Range2, p: &Part| r.bounds(p).ok_or(ExtremaRefusal::Unsupported);
    // A lower bound for a pair: the gap between their boxes, or better the
    // gap between their projections on the line joining their boxes'
    // centres, which closes on the true distance as the pieces shrink.
    let lower_bound = |x: &Range2, y: &Range2, bx: &Aabb, by: &Aabb| -> Scalar {
        let boxes = bx.distance(by);
        let d = (bx.lo + bx.hi) * 0.5 - (by.lo + by.hi) * 0.5;
        let length = d.length();
        if length == 0.0 {
            return boxes;
        }
        let n = d / length;
        match (x.projection(&pa, n), y.projection(&pb, n)) {
            (Some(px), Some(py)) => boxes.max(px.0 - py.1),
            _ => boxes,
        }
    };
    // Pending pairs with their lower bounds, the smallest searched first.
    let (ba, bb) = (bounds(&ra, &pa)?, bounds(&rb, &pb)?);
    let mut pending: Vec<(Scalar, Range2, Range2)> =
        vec![(lower_bound(&ra, &rb, &ba, &bb), ra, rb)];
    let mut work = 0usize;
    loop {
        // The smallest lower bound still open.
        let Some(index) = pending
            .iter()
            .enumerate()
            .min_by(|x, y| x.1 .0.total_cmp(&y.1 .0))
            .map(|(i, _)| i)
        else {
            best.lower = best.upper;
            return Ok(best);
        };
        let (lower, a, b) = pending.swap_remove(index);
        best.lower = lower.min(best.upper);
        if best.upper - lower <= accuracy {
            return Ok(best);
        }
        work += 1;
        if work > 200_000 {
            return Err(ExtremaRefusal::Budget);
        }
        // Split the side whose image is larger.
        let (da, db) = (bounds(&a, &pa)?.diameter(), bounds(&b, &pb)?.diameter());
        let (split_first, parts) = if (da >= db && !a.is_point()) || b.is_point() {
            (true, a.split())
        } else {
            (false, b.split())
        };
        for part in parts {
            let (x, y) = if split_first { (part, b) } else { (a, part) };
            let (bx, by) = (bounds(&x, &pa)?, bounds(&y, &pb)?);
            let low = lower_bound(&x, &y, &bx, &by);
            if low > best.upper {
                continue;
            }
            let (p, q) = (sample(&x, &pa)?, sample(&y, &pb)?);
            let d = (p.0 - q.0).length();
            if d < best.upper {
                best.upper = d;
                best.on_first = p.0;
                best.on_second = q.0;
                best.first = p.1;
                best.second = q.1;
            }
            if low <= best.upper {
                pending.push((low, x, y));
            }
        }
        pending.retain(|(low, _, _)| *low <= best.upper);
    }
}
