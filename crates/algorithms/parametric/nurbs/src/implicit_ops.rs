//! Working with section curves through their implicit form (ADR 0077):
//! pcurves of a section on any analytic face, and a curve of the section
//! families against a surface.
//!
//! Every section family is a field's zero set on a carrier: an
//! `ImplicitSection3` is one by definition. A ruled section (ADR 0076) is
//! the root branch of `a v^2 + b v + c = 0` on its ruled carrier, and a
//! torus section is where `A cos u + B sin u = C` on its torus. Each is
//! traced once in its own carrier's parameters. The stretch that matches
//! the curve is cut out by its end points, so the operations below have one
//! implementation for all of them.

use axiolid_core::{Interval, Point2, Scalar, Vec2};
use axiolid_curve::implicit::Cell;
use axiolid_curve::{
    Axis, Basis, Carrier, Curve3, Field2, ImplicitCurve2, ImplicitSection3, SeriesField2, Trig2,
};
use axiolid_surface::Surface;
use core::f64::consts::{PI, TAU};

use crate::curve_tower::CurveTower;
use crate::exact_curve_intersection::{
    ExactCurveHit, ExactCurveIntersection, ExactCurveParameter, ExactCurveRefusal, Isolated,
};
use crate::exact_surface_intersection::ExactIntersectionRefusal;
use crate::field::{carrier_of, section_field};
use crate::implicit_trace::Periodic;
use crate::root_chain::{isolate, Budget, End, RootRefusal, Zone, ROOT_BUDGET};

/// Every component of the section of `surface` by `other` in a window of
/// `surface`'s parameters, as curves in those parameters.
///
/// A window a whole turn wide in a periodic parameter wraps.
///
/// # Errors
///
/// As [`crate::implicit_surface_intersection`]; `Disjoint` is not an
/// error here (the list is empty).
pub fn trace_section_pcurves(
    surface: &Surface,
    other: &Surface,
    window: (Point2, Point2),
) -> Result<Vec<ImplicitCurve2>, ExactIntersectionRefusal> {
    let carrier = carrier_of(surface).ok_or(ExactIntersectionRefusal::UnsupportedPair)?;
    let field = section_field(&carrier, other).ok_or(ExactIntersectionRefusal::UnsupportedPair)?;
    // A quarter of a section's budget: the section itself is traced
    // already, and a face without these pcurves reads it back instead.
    trace_field(&carrier, &field, window, crate::implicit_trace::BUDGET / 4)
}

fn trace_field(
    carrier: &Carrier,
    field: &Field2,
    window: (Point2, Point2),
    budget: usize,
) -> Result<Vec<ImplicitCurve2>, ExactIntersectionRefusal> {
    let (pu, pv) = carrier.periodic();
    let (lo, hi) = window;
    let wraps = |periodic: bool, a: Scalar, b: Scalar| {
        periodic && ((b - a) - TAU).abs() <= 1e-12 * (1.0 + a.abs())
    };
    let periodic = Periodic {
        u: wraps(pu, lo.x, hi.x),
        v: wraps(pv, lo.y, hi.y),
    };
    crate::implicit_trace::trace_within(field, Cell { lo, hi }, periodic, budget)
        .map_err(crate::implicit_section::refusal_of)
}

/// The stretch of one of `curves` that runs from `start` through `first`
/// and then `second` to `end` (all in parameters, each possibly a whole
/// period off), as a curve of its own running from `start` to `end`.
/// `closed` asks for a whole loop from `start` back to itself, in the
/// direction that meets `first` before `second`.
#[must_use]
pub fn extract_stretch(
    curves: &[ImplicitCurve2],
    periodic: (bool, bool),
    start: Point2,
    [first, second]: [Point2; 2],
    end: Point2,
    closed: bool,
) -> Option<ImplicitCurve2> {
    let shifts = |p: Point2| -> Vec<Point2> {
        let ku: &[Scalar] = if periodic.0 {
            &[0.0, -1.0, 1.0, -2.0, 2.0]
        } else {
            &[0.0]
        };
        let kv: &[Scalar] = if periodic.1 {
            &[0.0, -1.0, 1.0, -2.0, 2.0]
        } else {
            &[0.0]
        };
        let mut out = Vec::new();
        for a in ku {
            for b in kv {
                out.push(Point2::new(p.x + a * TAU, p.y + b * TAU));
            }
        }
        out
    };
    let find = |curve: &ImplicitCurve2, p: Point2| -> Option<Scalar> {
        let eps = 1e-7 * (1.0 + p.x.abs().max(p.y.abs()));
        let mut best: Option<(Scalar, Scalar)> = None;
        for q in shifts(p) {
            let Some(t) = curve.parameter_of(q) else {
                continue;
            };
            let Some(on) = curve.point(t) else {
                continue;
            };
            // Compare up to whole periods.
            let mut d = on - q;
            if periodic.0 {
                d.x -= (d.x / TAU).round() * TAU;
            }
            if periodic.1 {
                d.y -= (d.y / TAU).round() * TAU;
            }
            let miss = d.length();
            if miss <= eps && best.is_none_or(|(m, _)| miss < m) {
                best = Some((miss, t));
            }
        }
        best.map(|(_, t)| t)
    };
    for curve in curves {
        let (Some(ts), Some(t1), Some(t2), Some(te)) = (
            find(curve, start),
            find(curve, first),
            find(curve, second),
            find(curve, end),
        ) else {
            continue;
        };
        let n = curve.end();
        let closure = curve.closure(periodic.0, periodic.1);
        let ahead = |t: Scalar| (t - ts).rem_euclid(n);
        let stretch = if closed {
            // The whole loop from `start` round to itself.
            let whole = curve.rotated(ts, closure?)?;
            if ahead(t1) < ahead(t2) {
                whole
            } else {
                whole.reversed()
            }
        } else if let Some(closure) = closure {
            // On a loop: forward if the inner points come in order before
            // `end` going forward from `start`.
            if ahead(t1) < ahead(t2) && ahead(t2) < ahead(te) {
                curve.sub(ts, te, Some(closure))?
            } else {
                curve.sub(te, ts, Some(closure))?.reversed()
            }
        } else if ts <= t1 && t1 <= t2 && t2 <= te {
            curve.sub(ts, te, None)?
        } else if te <= t2 && t2 <= t1 && t1 <= ts {
            curve.sub(te, ts, None)?.reversed()
        } else {
            continue;
        };
        return Some(stretch);
    }
    None
}

/// The carrier and defining field of a section family, or `None` for
/// other curves.
fn defining_field(curve: &Curve3) -> Option<(Carrier, Field2)> {
    let fourier = |t: &Trig2| vec![t.constant, t.cos, t.sin, t.cos2, t.sin2];
    match curve {
        Curve3::ImplicitSection(s) => Some((s.carrier.clone(), s.curve.field.clone())),
        Curve3::RuledSection(r) => {
            // c(u) + b(u) v + a(u) v^2, u harmonic, v a power.
            let (a, b, c) = (
                fourier(&r.graph.a),
                fourier(&r.graph.b),
                fourier(&r.graph.c),
            );
            let coefficients = (0..5).map(|i| vec![c[i], b[i], a[i]]).collect();
            Some((
                Carrier::Ruled(r.carrier),
                Field2::Series(SeriesField2 {
                    u: Basis::Fourier,
                    v: Basis::Power,
                    coefficients,
                }),
            ))
        }
        Curve3::TorusSection(t) => {
            // A(v) cos u + B(v) sin u - C(v), both harmonic.
            let (a, b, c) = (
                fourier(&t.graph.a),
                fourier(&t.graph.b),
                fourier(&t.graph.c),
            );
            let coefficients = vec![c.iter().map(|x| -x).collect(), a, b];
            Some((
                Carrier::Torus(t.torus),
                Field2::Series(SeriesField2 {
                    u: Basis::Fourier,
                    v: Basis::Fourier,
                    coefficients,
                }),
            ))
        }
        _ => None,
    }
}

/// A section-family curve over `span` as an [`ImplicitSection3`], with the
/// same point set and direction; `None` for other families or where the
/// trace is refused.
#[must_use]
pub fn implicit_view(curve: &Curve3, span: Interval) -> Option<ImplicitSection3> {
    if let Curve3::ImplicitSection(s) = curve {
        let closed = (span.end - span.start).abs() >= s.curve.end() - 1e-12;
        let sub = if closed {
            s.curve.clone()
        } else {
            let (a, b) = (span.start.min(span.end), span.start.max(span.end));
            let piece = s.curve.sub(a, b, None)?;
            if span.end < span.start {
                piece.reversed()
            } else {
                piece
            }
        };
        return Some(ImplicitSection3 {
            carrier: s.carrier.clone(),
            curve: sub,
        });
    }
    let (carrier, field) = defining_field(curve)?;
    let at = |t: Scalar| -> Option<Point2> {
        let p = axiolid_evaluate::evaluate3(curve, t).ok()?;
        let (u, v) = carrier.parameters(p);
        Some(Point2::new(u, v))
    };
    // Carrier parameters along the span, unwrapped, to size the window.
    let n = 64;
    let mut samples = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let t = span.start + (span.end - span.start) * i as Scalar / n as Scalar;
        let mut p = at(t)?;
        if let Some(q) = samples.last().copied() {
            let q: Point2 = q;
            let (pu, pv) = carrier.periodic();
            if pu {
                p.x += ((q.x - p.x) / TAU).round() * TAU;
            }
            if pv {
                p.y += ((q.y - p.y) / TAU).round() * TAU;
            }
        }
        samples.push(p);
    }
    let (mut lo, mut hi) = (samples[0], samples[0]);
    for p in &samples {
        lo = lo.min(*p);
        hi = hi.max(*p);
    }
    let pad = (hi - lo) * 0.25 + Vec2::splat(0.05);
    let (mut lo, mut hi) = (lo - pad, hi + pad);
    // Never more than a turn in a periodic parameter.
    let (pu, pv) = carrier.periodic();
    if pu && hi.x - lo.x >= TAU {
        let c = 0.5 * (lo.x + hi.x);
        (lo.x, hi.x) = (c - PI, c + PI);
    }
    if pv && hi.y - lo.y >= TAU {
        let c = 0.5 * (lo.y + hi.y);
        (lo.y, hi.y) = (c - PI, c + PI);
    }
    let curves = trace_field(&carrier, &field, (lo, hi), crate::implicit_trace::BUDGET).ok()?;
    // A loop closes in space, whatever whole turns its parameters make.
    let (first, last) = (
        axiolid_evaluate::evaluate3(curve, span.start).ok()?,
        axiolid_evaluate::evaluate3(curve, span.end).ok()?,
    );
    let closed = (first - last).length() <= 1e-9 * (1.0 + first.length());
    let stretch = extract_stretch(
        &curves,
        (pu, pv),
        samples[0],
        [samples[n / 3], samples[2 * n / 3]],
        samples[n],
        closed,
    )?;
    Some(ImplicitSection3 {
        carrier,
        curve: stretch,
    })
}

/// Where a section-family curve meets a surface, with parameters on the
/// curve itself.
///
/// The curve is read as an implicit curve on its carrier (over `span`, or
/// its whole traced extent for an [`ImplicitSection3`]); the surface's
/// equation on that carrier is a second field, and its roots along the
/// curve are isolated with interval bounds by a Rolle chain on the
/// derivatives along the curve (ADR 0077). A root where the second field
/// changes sign is a crossing (multiplicity 1); a root of higher
/// multiplicity, up to eight, is a touch or a flat crossing and reports its
/// order; a curve on which the second field vanishes identically is
/// `Contained`.
///
/// # Errors
///
/// `UnsupportedCurve` for a curve that is not a section family or whose
/// trace was refused, `UnsupportedSurface` for a B-spline surface,
/// `Undecided` where the roots could not all be isolated with their
/// multiplicity.
pub fn section_curve_surface_intersection(
    curve: &Curve3,
    span: Interval,
    surface: &Surface,
) -> Result<ExactCurveIntersection, ExactCurveRefusal> {
    let view = implicit_view(curve, span).ok_or(ExactCurveRefusal::UnsupportedCurve)?;
    let other =
        section_field(&view.carrier, surface).ok_or(ExactCurveRefusal::UnsupportedSurface)?;
    let own = &view.curve;
    // Contained: the second field vanishes along the whole curve, to the
    // rounding of its own size.
    let scale = 1e-9 * other.magnitude().max(1.0);
    let samples: Vec<Scalar> = (0..=32)
        .filter_map(|i| {
            let p = own.point(own.end() * i as Scalar / 32.0)?;
            Some(other.value(p))
        })
        .collect();
    if samples.len() == 33 && samples.iter().all(|h| h.abs() <= scale) {
        return Ok(ExactCurveIntersection::Contained);
    }
    let (pu, pv) = view.carrier.periodic();
    let closed = own.closure(pu, pv).is_some();
    let roots = roots_along(own, &other, closed).map_err(|_| ExactCurveRefusal::Undecided)?;
    let mut hits = Vec::new();
    for (t, multiplicity) in roots {
        let Some(p) = own.point(t) else { continue };
        let point = view.carrier.jet(p.x, p.y).point;
        // Back to the curve's own parameter.
        let parameter = match curve {
            Curve3::ImplicitSection(_) => {
                if span.end < span.start {
                    span.start - t
                } else {
                    span.start + t
                }
            }
            _ => {
                match axiolid_evaluate::curve::invert3(curve, point, axiolid_core::Tolerance::METRE)
                {
                    Ok(x) => nearest_turn(x, span),
                    Err(_) => continue,
                }
            }
        };
        hits.push(ExactCurveHit {
            parameter: ExactCurveParameter::Certified(Isolated::new(parameter)),
            multiplicity,
            point,
        });
    }
    Ok(ExactCurveIntersection::Points(hits))
}

/// The roots of `other` along an implicit curve, by the curve's parameter,
/// each with its multiplicity, ascending: 1 where the field changes sign,
/// more where it touches zero or crosses it with a flat contact.
///
/// The stretches of consecutive cells are isolated by a Rolle chain on the
/// intrinsic derivatives along the curve ([`crate::root_chain`],
/// [`crate::curve_tower`]); a bridge cell is scanned for sign changes (its
/// roots are not certified). A `closed` curve's start and end are one
/// point, so a root there is reported once.
///
/// # Errors
///
/// The refusal of the isolation: a contact above order eight, one hidden by
/// rounding, or an exhausted work budget.
pub(crate) fn roots_along(
    own: &ImplicitCurve2,
    other: &Field2,
    closed: bool,
) -> Result<Vec<(Scalar, usize)>, RootRefusal> {
    roots_along_within(own, other, closed, ROOT_BUDGET)
}

/// [`roots_along`] within a work budget of `budget` enclosures.
pub(crate) fn roots_along_within(
    own: &ImplicitCurve2,
    other: &Field2,
    closed: bool,
    budget: usize,
) -> Result<Vec<(Scalar, usize)>, RootRefusal> {
    let tower = CurveTower::new(own, other);
    let mut budget = Budget::new(budget);
    let cells = own.cells.len();
    let mut zones: Vec<Zone> = Vec::new();
    let mut index = 0;
    while index < cells {
        if own.cells[index].bridge.is_some() {
            for (t, m) in scan_bridge(own, other, index) {
                zones.push(Zone {
                    lo: t,
                    hi: t,
                    at: t,
                    m,
                    end: None,
                });
            }
            index += 1;
            continue;
        }
        let start = index;
        while index < cells && own.cells[index].bridge.is_none() {
            index += 1;
        }
        let ends = [
            if start == 0 { End::Curve } else { End::Inner },
            if index == cells {
                End::Curve
            } else {
                End::Inner
            },
        ];
        zones.extend(isolate(
            &tower,
            start as Scalar,
            index as Scalar,
            ends,
            &mut budget,
        )?);
    }
    zones.sort_by(|a, b| a.at.total_cmp(&b.at));
    // A loop's seam is one point, met at both ends of the parameter.
    if closed && zones.len() > 1 {
        let (first, last) = (zones[0], zones[zones.len() - 1]);
        if first.lo == 0.0 && last.hi == cells as Scalar {
            zones.pop();
            zones[0].m = first.m.max(last.m);
        }
    }
    let roots: Vec<(Scalar, usize)> = zones.iter().map(|z| (z.at, z.m)).collect();
    debug_assert!(roots.windows(2).all(|w| w[0].0 < w[1].0));
    Ok(roots)
}

/// The sign changes of `other` along a bridge cell, by 64 samples. A
/// bridge into a crossing is too short and too close to the singular point
/// for bounds, so its roots are scanned and not certified.
// TODO(#441): certify the roots along bridge cells.
fn scan_bridge(own: &ImplicitCurve2, other: &Field2, index: usize) -> Vec<(Scalar, usize)> {
    let cell = &own.cells[index];
    let h = |s: Scalar| -> Option<Scalar> {
        let free = cell.from + (cell.to - cell.from) * s;
        let solved = own.solve_cell(cell, free)?;
        let p = match cell.axis {
            Axis::U => Point2::new(free, solved),
            Axis::V => Point2::new(solved, free),
        };
        Some(other.value(p))
    };
    let mut out = Vec::new();
    let n = 64;
    let at = |k: usize| k as Scalar / n as Scalar;
    let mut last = h(at(0));
    for k in 1..=n {
        let now = h(at(k));
        if let (Some(a), Some(b)) = (last, now) {
            if a == 0.0 && k == 1 {
                out.push((at(0), 1));
            } else if (a < 0.0) != (b < 0.0) && b != 0.0 {
                out.push((bisect_bridge(&h, at(k - 1), at(k), a), 1));
            } else if b == 0.0 {
                out.push((at(k), 1));
            }
        }
        last = now;
    }
    let mut roots: Vec<(Scalar, usize)> = out
        .into_iter()
        .map(|(s, m)| (index as Scalar + s, m))
        .collect();
    roots.dedup_by(|a, b| (a.0 - b.0).abs() <= 1e-9);
    roots
}

fn bisect_bridge(
    h: &dyn Fn(Scalar) -> Option<Scalar>,
    mut s0: Scalar,
    mut s1: Scalar,
    a: Scalar,
) -> Scalar {
    let negative = a < 0.0;
    for _ in 0..100 {
        let m = 0.5 * (s0 + s1);
        if m <= s0 || m >= s1 {
            break;
        }
        match h(m) {
            Some(0.0) => return m,
            Some(v) if (v < 0.0) == negative => s0 = m,
            Some(_) => s1 = m,
            None => break,
        }
    }
    0.5 * (s0 + s1)
}

/// `x` moved by whole turns into `span` when that is possible.
fn nearest_turn(x: Scalar, span: Interval) -> Scalar {
    let (lo, hi) = (span.start.min(span.end), span.start.max(span.end));
    for k in [0.0, 1.0, -1.0, 2.0, -2.0] {
        let y = x + k * TAU;
        if y >= lo - 1e-9 && y <= hi + 1e-9 {
            return y;
        }
    }
    x
}

/// A line, circle or ellipse against a B-spline surface: the curve is the
/// meeting of two analytic surfaces (a line of two planes, a circle of its
/// plane and a sphere, an ellipse of its plane and an elliptical cylinder).
/// The first is traced on the spline (certified, over its Bernstein
/// patches); the second's roots along each traced piece are the hits.
///
/// # Errors
///
/// `UnsupportedCurve` for other families, `UnsupportedSurface` for a
/// surface that is not a B-spline or a spline whose trace is refused.
pub fn conic_spline_intersection(
    curve: &Curve3,
    surface: &Surface,
) -> Result<ExactCurveIntersection, ExactCurveRefusal> {
    use axiolid_core::Frame3;
    use axiolid_surface::{EllipticalCylinder, Plane, Sphere};
    if !matches!(surface, Surface::BSpline(_)) {
        return Err(ExactCurveRefusal::UnsupportedSurface);
    }
    let plane = |origin: axiolid_core::Point3, z: axiolid_core::Vec3| -> Surface {
        let z = z.normalize();
        let helper = if z.x.abs() < 0.9 {
            axiolid_core::Vec3::X
        } else {
            axiolid_core::Vec3::Y
        };
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
    let (first, second) = match curve {
        Curve3::Line(l) => {
            let d = l.direction.normalize();
            let helper = if d.x.abs() < 0.9 {
                axiolid_core::Vec3::X
            } else {
                axiolid_core::Vec3::Y
            };
            let n1 = d.cross(helper).normalize();
            (plane(l.origin, n1), plane(l.origin, d.cross(n1)))
        }
        Curve3::Circle(c) => (
            plane(c.frame.origin, c.frame.x.cross(c.frame.y)),
            Surface::Sphere(Sphere {
                frame: c.frame,
                radius: c.radius,
            }),
        ),
        Curve3::Ellipse(e) => (
            plane(e.frame.origin, e.frame.x.cross(e.frame.y)),
            Surface::EllipticalCylinder(EllipticalCylinder {
                frame: Frame3 {
                    z: e.frame.x.cross(e.frame.y).normalize(),
                    ..e.frame
                },
                semi_axis_x: e.semi_axis_x,
                semi_axis_y: e.semi_axis_y,
            }),
        ),
        _ => return Err(ExactCurveRefusal::UnsupportedCurve),
    };
    let pieces = match crate::implicit_section::implicit_surface_intersection(surface, &first, None)
    {
        Ok(pieces) => pieces,
        Err(ExactIntersectionRefusal::Disjoint) => {
            return Ok(ExactCurveIntersection::Points(Vec::new()))
        }
        Err(_) => return Err(ExactCurveRefusal::UnsupportedSurface),
    };
    let mut hits = Vec::new();
    for piece in pieces {
        let end = piece.curve.end();
        let traced = Curve3::ImplicitSection(piece);
        match section_curve_surface_intersection(&traced, Interval::new(0.0, end), &second)? {
            ExactCurveIntersection::Points(found) => {
                for hit in found {
                    let t = match curve {
                        Curve3::Line(l) => {
                            (hit.point - l.origin).dot(l.direction) / l.direction.length_squared()
                        }
                        _ => match axiolid_evaluate::curve::invert3(
                            curve,
                            hit.point,
                            axiolid_core::Tolerance::METRE,
                        ) {
                            Ok(t) => t,
                            Err(_) => continue,
                        },
                    };
                    hits.push(ExactCurveHit {
                        parameter: ExactCurveParameter::Certified(Isolated::new(t)),
                        ..hit
                    });
                }
            }
            // The traced piece lies on the second surface too: the curve
            // runs in the spline.
            _ => return Ok(ExactCurveIntersection::Contained),
        }
    }
    hits.sort_by(|a, b| a.parameter.approx().total_cmp(&b.parameter.approx()));
    Ok(ExactCurveIntersection::Points(hits))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::implicit_trace::trace_within;
    use axiolid_curve::implicit::Basis;

    fn power(terms: &[(usize, usize, Scalar)]) -> Field2 {
        let mut c = vec![vec![0.0; 10]; 10];
        for &(i, j, x) in terms {
            c[i][j] += x;
        }
        Field2::Series(SeriesField2 {
            u: Basis::Power,
            v: Basis::Power,
            coefficients: c,
        })
    }

    fn traced(field: &Field2, lo: Point2, hi: Point2) -> ImplicitCurve2 {
        let mut curves = trace_within(
            field,
            Cell { lo, hi },
            Periodic { u: false, v: false },
            1 << 20,
        )
        .expect("a traced curve");
        assert_eq!(curves.len(), 1);
        curves.remove(0)
    }

    fn window() -> (Point2, Point2) {
        (Point2::new(-1.0, -1.0), Point2::new(1.0, 1.0))
    }

    /// `v - 0.3 u^2` and `F + u^k`: on the curve the second field is `u^k`.
    fn graph(k: usize) -> (ImplicitCurve2, Field2) {
        let f = [(0, 1, 1.0), (2, 0, -0.3)];
        let mut h = f.to_vec();
        h.push((k, 0, 1.0));
        let (lo, hi) = window();
        (traced(&power(&f), lo, hi), power(&h))
    }

    #[test]
    fn a_contact_of_any_order_up_to_eight_reports_its_multiplicity() {
        for k in 2..=8 {
            let (curve, h) = graph(k);
            let roots = roots_along(&curve, &h, false).expect("decided");
            assert_eq!(roots.len(), 1, "u^{k}: {roots:?}");
            assert_eq!(roots[0].1, k, "u^{k}");
            let p = curve.point(roots[0].0).unwrap();
            assert!(p.x.abs() < 1e-6 && p.y.abs() < 1e-6, "u^{k}: {p:?}");
        }
    }

    #[test]
    fn a_contact_above_the_cap_is_undecided() {
        let (curve, h) = graph(9);
        assert!(roots_along(&curve, &h, false).is_err());
    }

    #[test]
    fn a_contact_at_the_end_of_the_curve_is_reported() {
        let f = power(&[(0, 1, 1.0), (2, 0, -0.3)]);
        let h = power(&[(0, 1, 1.0), (2, 0, -0.3), (4, 0, 1.0)]);
        let curve = traced(&f, Point2::new(-1.0, -1.0), Point2::new(0.0, 1.0));
        let roots = roots_along(&curve, &h, false).expect("decided");
        assert_eq!(roots.len(), 1, "{roots:?}");
        assert_eq!(roots[0].1, 4);
        assert!(curve.point(roots[0].0).unwrap().x.abs() < 1e-6);
    }

    #[test]
    fn the_tangent_at_a_loops_seam_touches_it_once() {
        let f = power(&[(2, 0, 1.0), (0, 2, 1.0), (0, 0, -0.25)]);
        let (lo, hi) = window();
        let curve = traced(&f, lo, hi);
        let p = curve.point(0.0).unwrap();
        let h = power(&[(1, 0, p.x), (0, 1, p.y), (0, 0, -(p.x * p.x + p.y * p.y))]);
        let roots = roots_along(&curve, &h, true).expect("decided");
        assert_eq!(roots.len(), 1, "{roots:?}");
        assert_eq!(roots[0].1, 2);
    }

    #[test]
    fn the_budget_refuses_by_name() {
        // v = 0 against u times four factors u^2 - r^2: nine simple roots,
        // more than the cap of the chain, so the leaves must split.
        let f = power(&[(0, 1, 1.0)]);
        let mut p = vec![0.0, 1.0];
        for r in [0.2f64, 0.5, 0.8, 0.9] {
            let mut q = vec![0.0; p.len() + 2];
            for (i, c) in p.iter().enumerate() {
                q[i + 2] += c;
                q[i] -= c * r * r;
            }
            p = q;
        }
        let terms: Vec<(usize, usize, Scalar)> =
            p.iter().enumerate().map(|(i, c)| (i, 0, *c)).collect();
        let h = power(&terms);
        let (lo, hi) = window();
        let curve = traced(&f, lo, hi);
        assert_eq!(roots_along(&curve, &h, false).expect("decided").len(), 9);
        assert_eq!(
            roots_along_within(&curve, &h, false, 2),
            Err(RootRefusal::Budget)
        );
    }

    /// Prints the most enclosures a decided case of this module's tests
    /// needed, for the calibration of `ROOT_BUDGET` (sixteen times the peak
    /// of every suite must stay below it). Run with
    /// `cargo test -p axiolid-nurbs --lib peak_enclosures -- --ignored --nocapture`.
    #[test]
    #[ignore = "calibration, not a check of behaviour"]
    fn peak_enclosures_of_decided_cases() {
        for k in 2..=8 {
            let (curve, h) = graph(k);
            roots_along(&curve, &h, false).expect("decided");
        }
        let peak = crate::root_chain::PEAK.load(std::sync::atomic::Ordering::Relaxed);
        println!("peak enclosures of the section cases: {peak}");
        assert!(16 * peak <= ROOT_BUDGET, "{peak}");
    }
}
