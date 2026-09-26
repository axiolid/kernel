//! Curve/curve intersection for the section families (#119, B8;
//! ADR 0077).
//!
//! **In the plane.** Every plane curve of the section families is a
//! field's zero set: a line or conic by its implicit equation, a
//! [`Sinusoid2`] by `v - w(u)`, a [`QuadraticGraph2`] by its quadratic, an
//! [`AngleGraph2`] by `A cos u + B sin u - C`, an [`ImplicitCurve2`] by its
//! own field. One piece is taken as traced cells, over its span (traced
//! with certified topology when it is not already an implicit curve), and
//! the other's field has its roots isolated along those cells with interval
//! bounds. Each root is kept when it lies on both pieces within their
//! spans: a graph's other branch shares its field, and is filtered here.
//!
//! **In space.** Two space curves meet only where one crosses a surface the
//! other lies on: the first is intersected with the second's carrier (a
//! traced section's surface, a conic's plane, a line's two planes) by the
//! certified curve/surface routines, and each point is kept when it lies on
//! the second curve within `tolerance`. In space the question is itself
//! only meaningful up to a tolerance: two curves built from rounded
//! numbers rarely meet exactly.
//!
//! [`Sinusoid2`]: axiolid_curve::Sinusoid2
//! [`QuadraticGraph2`]: axiolid_curve::QuadraticGraph2
//! [`AngleGraph2`]: axiolid_curve::AngleGraph2
//! [`ImplicitCurve2`]: axiolid_curve::ImplicitCurve2

use axiolid_core::{Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Basis, Carrier, Curve2, Curve3, Field2, ImplicitCurve2, SeriesField2, Trig2};
use axiolid_evaluate::curve::{evaluate2, evaluate3, locate2, locate3};
use axiolid_surface::{Plane, Surface};

use crate::exact_curve_intersection::{ExactCurveIntersection, ExactCurveRefusal};

/// One point where two curves meet, with its parameter on each.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveHit {
    /// Parameter on the first curve.
    pub first: Scalar,
    /// Parameter on the second curve.
    pub second: Scalar,
    /// The point (on the first curve at `first`).
    pub point: Point3,
    /// 1 where the curves cross, 2 or more where they touch.
    pub multiplicity: usize,
}

fn series(u: Basis, v: Basis, coefficients: Vec<Vec<Scalar>>) -> Field2 {
    Field2::Series(SeriesField2 { u, v, coefficients })
}

fn fourier(t: &Trig2) -> Vec<Scalar> {
    vec![t.constant, t.cos, t.sin, t.cos2, t.sin2]
}

/// A plane curve's defining field, or `None` for a family without one.
fn plane_field(curve: &Curve2) -> Option<Field2> {
    let p = Basis::Power;
    let f = Basis::Fourier;
    Some(match curve {
        Curve2::Line(l) => {
            let n = Vec2::new(-l.direction.y, l.direction.x);
            if n.length() == 0.0 {
                return None;
            }
            let n = n.normalize();
            series(p, p, vec![vec![-n.dot(l.origin), n.y], vec![n.x, 0.0]])
        }
        Curve2::Circle(c) => conic_field(c.frame, c.radius, c.radius)?,
        Curve2::Ellipse(e) => conic_field(e.frame, e.semi_axis_x, e.semi_axis_y)?,
        Curve2::Sinusoid(w) => series(
            f,
            p,
            vec![vec![-w.mean, 1.0], vec![-w.cosine, 0.0], vec![-w.sine, 0.0]],
        ),
        Curve2::QuadraticGraph(g) => {
            let (a, b, c) = (fourier(&g.a), fourier(&g.b), fourier(&g.c));
            series(f, p, (0..5).map(|i| vec![c[i], b[i], a[i]]).collect())
        }
        Curve2::AngleGraph(g) => {
            let (a, b, c) = (fourier(&g.a), fourier(&g.b), fourier(&g.c));
            series(f, f, vec![c.iter().map(|x| -x).collect(), a, b])
        }
        Curve2::Implicit(c) => c.field.clone(),
        _ => return None,
    })
}

/// `b^2 x^2 + a^2 y^2 - a^2 b^2` in an orthonormal frame's coordinates.
fn conic_field(frame: axiolid_core::Frame2, a: Scalar, b: Scalar) -> Option<Field2> {
    let (x, y) = (frame.x, frame.y);
    if (x.length() - 1.0).abs() > 1e-12
        || (y.length() - 1.0).abs() > 1e-12
        || x.dot(y).abs() > 1e-12
    {
        return None;
    }
    // x_l = X . (p - o), y_l = Y . (p - o): affine in (u, v).
    let lin = |axis: Vec2| [-axis.dot(frame.origin), axis.x, axis.y];
    let (lx, ly) = (lin(x), lin(y));
    // Square of c0 + c1 u + c2 v, as coefficients [i][j] of u^i v^j.
    let square = |l: [Scalar; 3]| {
        let mut c = vec![vec![0.0; 3]; 3];
        c[0][0] = l[0] * l[0];
        c[1][0] = 2.0 * l[0] * l[1];
        c[0][1] = 2.0 * l[0] * l[2];
        c[2][0] = l[1] * l[1];
        c[1][1] = 2.0 * l[1] * l[2];
        c[0][2] = l[2] * l[2];
        c
    };
    let (sx, sy) = (square(lx), square(ly));
    let mut c = vec![vec![0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = b * b * sx[i][j] + a * a * sy[i][j];
        }
    }
    c[0][0] -= a * a * b * b;
    Some(series(Basis::Power, Basis::Power, c))
}

/// A plane curve's piece as traced cells of its own field, over a window
/// holding it.
fn traced_piece(curve: &Curve2, span: Interval) -> Option<ImplicitCurve2> {
    if let Curve2::Implicit(c) = curve {
        let (a, b) = (span.start.min(span.end), span.start.max(span.end));
        return c.sub(a, b, None);
    }
    let field = plane_field(curve)?;
    let n = 128;
    let at = |i: usize| {
        evaluate2(
            curve,
            span.start + (span.end - span.start) * i as Scalar / n as Scalar,
        )
        .ok()
    };
    let samples: Vec<Point2> = (0..=n).map(at).collect::<Option<_>>()?;
    let (mut lo, mut hi) = (samples[0], samples[0]);
    for p in &samples {
        lo = lo.min(*p);
        hi = hi.max(*p);
    }
    let pad = (hi - lo) * 0.1 + Vec2::splat(1e-3 * (1.0 + (hi - lo).length()));
    let curves = crate::implicit_trace::trace(
        &field,
        axiolid_curve::implicit::Cell {
            lo: lo - pad,
            hi: hi + pad,
        },
        crate::implicit_trace::Periodic { u: false, v: false },
    )
    .ok()?;
    let closed = (samples[0] - samples[n]).length() <= 1e-12 * (1.0 + samples[0].length());
    crate::implicit_ops::extract_stretch(
        &curves,
        (false, false),
        samples[0],
        [samples[n / 3], samples[2 * n / 3]],
        samples[n],
        closed,
    )
}

/// Whether `t` lies in `span`, either way round.
fn within(t: Scalar, span: Interval) -> bool {
    let (lo, hi) = (span.start.min(span.end), span.start.max(span.end));
    let slack = 1e-9 * (1.0 + lo.abs().max(hi.abs()));
    t >= lo - slack && t <= hi + slack
}

/// Where two plane curves meet, each over its span.
///
/// Lines, circles, ellipses, sinusoids, quadratic and angle graphs and
/// implicit curves, in any pairing. Hits are ordered along the first curve.
///
/// # Errors
///
/// `UnsupportedCurve` when neither curve has a field (a B-spline or a
/// lifted pcurve: B-spline pairs take the certified NURBS tier), or when a
/// piece's trace is refused.
pub fn section_curve_curve_intersection2(
    first: &Curve2,
    first_span: Interval,
    second: &Curve2,
    second_span: Interval,
    tolerance: Tolerance,
) -> Result<Vec<CurveCurveHit>, ExactCurveRefusal> {
    // Trace one piece, and find the other's field's roots along it.
    let (traced, other, swapped) = match (plane_field(first), plane_field(second)) {
        (_, Some(f)) => (traced_piece(first, first_span), f, false),
        (Some(f), None) => (traced_piece(second, second_span), f, true),
        (None, None) => return Err(ExactCurveRefusal::UnsupportedCurve),
    };
    let traced = match traced {
        Some(t) => t,
        None if plane_field(if swapped { second } else { first }).is_none() => {
            return Err(ExactCurveRefusal::UnsupportedCurve)
        }
        None => return Err(ExactCurveRefusal::UnsupportedCurve),
    };
    let roots = crate::implicit_ops::roots_along(&traced, &other);
    let mut hits = Vec::new();
    for (t, multiplicity) in roots {
        let Some(p) = traced.point(t) else { continue };
        let (a, b) = match (locate2(first, p, tolerance), locate2(second, p, tolerance)) {
            (Ok(a), Ok(b)) => (a, b),
            _ => continue,
        };
        if !within(a, first_span) || !within(b, second_span) {
            continue;
        }
        hits.push(CurveCurveHit {
            first: a,
            second: b,
            point: Point3::new(p.x, p.y, 0.0),
            multiplicity,
        });
    }
    hits.sort_by(|x, y| x.first.total_cmp(&y.first));
    hits.dedup_by(|x, y| (x.first - y.first).abs() <= 1e-9 * (1.0 + x.first.abs()));
    Ok(hits)
}

/// The surface a space curve lies on that best separates it from others:
/// a traced section's carrier, a conic's plane; a line needs two planes.
fn carriers_of(curve: &Curve3) -> Option<Vec<Surface>> {
    let plane = |origin: Point3, z: Vec3| -> Surface {
        let z = z.normalize();
        let helper = if z.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
        let x = helper.cross(z).normalize();
        Surface::Plane(Plane {
            frame: Frame3 {
                origin,
                x,
                y: z.cross(x),
                z,
            },
        })
    };
    Some(match curve {
        Curve3::Line(l) => {
            let d = l.direction.normalize();
            let helper = if d.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
            let n = d.cross(helper).normalize();
            vec![plane(l.origin, n), plane(l.origin, d.cross(n))]
        }
        Curve3::Circle(c) => vec![plane(c.frame.origin, c.frame.x.cross(c.frame.y))],
        Curve3::Ellipse(e) => vec![plane(e.frame.origin, e.frame.x.cross(e.frame.y))],
        Curve3::ImplicitSection(s) => vec![surface_of(&s.carrier)?],
        Curve3::RuledSection(r) => vec![surface_of(&Carrier::Ruled(r.carrier))?],
        Curve3::TorusSection(t) => vec![surface_of(&Carrier::Torus(t.torus))?],
        _ => return None,
    })
}

/// A carrier as the surface it copies.
pub(crate) fn surface_of(carrier: &Carrier) -> Option<Surface> {
    use axiolid_surface::{Cone, Cylinder, EllipticalCylinder, Sphere, Torus};
    Some(match carrier {
        Carrier::Plane(f) => Surface::Plane(Plane { frame: *f }),
        Carrier::Ruled(k) => {
            if k.slope != 0.0 {
                if k.x_radius != k.y_radius {
                    return None;
                }
                Surface::Cone(Cone {
                    frame: k.frame,
                    radius: k.x_radius,
                    semi_angle: k.slope.atan(),
                })
            } else if k.x_radius == k.y_radius {
                Surface::Cylinder(Cylinder {
                    frame: k.frame,
                    radius: k.x_radius,
                })
            } else {
                Surface::EllipticalCylinder(EllipticalCylinder {
                    frame: k.frame,
                    semi_axis_x: k.x_radius,
                    semi_axis_y: k.y_radius,
                })
            }
        }
        Carrier::Sphere { frame, radius } => Surface::Sphere(Sphere {
            frame: *frame,
            radius: *radius,
        }),
        Carrier::Torus(t) => Surface::Torus(Torus {
            frame: t.frame,
            major_radius: t.major_radius,
            minor_radius: t.minor_radius,
        }),
        Carrier::Spline(b) => Surface::BSpline((**b).clone()),
    })
}

/// Where two space curves meet, each over its span, within `tolerance`.
///
/// The first is intersected with each surface the second lies on
/// (certified curve/surface intersection), and a point is kept where it
/// lies on the second curve within `tolerance` and inside both spans.
///
/// # Errors
///
/// `UnsupportedCurve` for a family neither routine handles.
pub fn section_curve_curve_intersection3(
    first: &Curve3,
    first_span: Interval,
    second: &Curve3,
    second_span: Interval,
    tolerance: Tolerance,
) -> Result<Vec<CurveCurveHit>, ExactCurveRefusal> {
    let carriers = carriers_of(second).ok_or(ExactCurveRefusal::UnsupportedCurve)?;
    let mut hits = Vec::new();
    for carrier in &carriers {
        let found = match first {
            Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_) => {
                crate::exact_curve_intersection::exact_curve_surface_intersection(first, carrier)
            }
            _ => {
                crate::implicit_ops::section_curve_surface_intersection(first, first_span, carrier)
            }
        };
        let points = match found {
            Ok(ExactCurveIntersection::Points(points)) => points,
            // The first curve lies on this carrier: the other carriers (a
            // line's second plane) or the membership test decide.
            Ok(_) => continue,
            Err(e) => return Err(e),
        };
        for hit in points {
            let Ok(b) = locate3(second, hit.point, tolerance) else {
                continue;
            };
            let Ok(on) = evaluate3(second, b) else {
                continue;
            };
            if (on - hit.point).length() > tolerance.linear() {
                continue;
            }
            let a = hit.parameter.approx();
            if !within(a, first_span) || !within(b, second_span) {
                continue;
            }
            hits.push(CurveCurveHit {
                first: a,
                second: b,
                point: hit.point,
                multiplicity: hit.multiplicity,
            });
        }
    }
    hits.sort_by(|x, y| x.first.total_cmp(&y.first));
    hits.dedup_by(|x, y| (x.first - y.first).abs() <= 1e-9 * (1.0 + x.first.abs()));
    Ok(hits)
}
