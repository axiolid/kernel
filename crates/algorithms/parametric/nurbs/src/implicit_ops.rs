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
use axiolid_curve::implicit::{bound_simple, partial, Cell, Range};
use axiolid_curve::{
    Axis, Basis, Carrier, Curve3, Field2, ImplicitCell, ImplicitCurve2, ImplicitSection3, Trig2,
};
use axiolid_surface::Surface;
use core::f64::consts::{PI, TAU};

use crate::exact_curve_intersection::{
    ExactCurveHit, ExactCurveIntersection, ExactCurveParameter, ExactCurveRefusal, Isolated,
};
use crate::exact_surface_intersection::ExactIntersectionRefusal;
use crate::field::{carrier_of, section_field};
use crate::implicit_trace::{trace, Periodic};

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
    trace_field(&carrier, &field, window)
}

fn trace_field(
    carrier: &Carrier,
    field: &Field2,
    window: (Point2, Point2),
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
    trace(field, Cell { lo, hi }, periodic).map_err(|_| ExactIntersectionRefusal::NotRegularCurve)
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
        Curve3::ImplicitSection(s) => Some((s.carrier, s.curve.field.clone())),
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
                Field2 {
                    u: Basis::Fourier,
                    v: Basis::Power,
                    coefficients,
                },
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
                Field2 {
                    u: Basis::Fourier,
                    v: Basis::Fourier,
                    coefficients,
                },
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
            carrier: s.carrier,
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
    let curves = trace_field(&carrier, &field, (lo, hi)).ok()?;
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
/// equation on that carrier is a second field, and its roots along each
/// cell are isolated with interval bounds over boxes certain to hold the
/// curve. A root where the second field changes sign is a crossing; one
/// where it only touches zero is reported with multiplicity 2; a curve on
/// which the second field vanishes identically is `Contained`.
///
/// # Errors
///
/// `UnsupportedCurve` for a curve that is not a section family or whose
/// trace was refused, `UnsupportedSurface` for a B-spline surface.
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
    let d_u = partial(&own.field, true);
    let d_v = partial(&own.field, false);
    let h_u = partial(&other, true);
    let h_v = partial(&other, false);
    let mut roots: Vec<(Scalar, usize)> = Vec::new();
    for (index, cell) in own.cells.iter().enumerate() {
        let job = CellRoots {
            curve: own,
            cell,
            other: &other,
            d_u: &d_u,
            d_v: &d_v,
            h_u: &h_u,
            h_v: &h_v,
        };
        let mut found = Vec::new();
        job.roots(0.0, 1.0, 0, &mut found);
        for (s, m) in found {
            roots.push((index as Scalar + s, m));
        }
    }
    roots.sort_by(|a, b| a.0.total_cmp(&b.0));
    roots.dedup_by(|a, b| (a.0 - b.0).abs() <= 1e-9);
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

/// The roots of a second field along one cell of an implicit curve.
struct CellRoots<'a> {
    curve: &'a ImplicitCurve2,
    cell: &'a ImplicitCell,
    other: &'a Field2,
    d_u: &'a Field2,
    d_v: &'a Field2,
    h_u: &'a Field2,
    h_v: &'a Field2,
}

impl CellRoots<'_> {
    fn free(&self, s: Scalar) -> Scalar {
        self.cell.from + (self.cell.to - self.cell.from) * s
    }

    fn place(&self, free: Scalar, solved: Scalar) -> Point2 {
        match self.cell.axis {
            Axis::U => Point2::new(free, solved),
            Axis::V => Point2::new(solved, free),
        }
    }

    fn solved(&self, free: Scalar) -> Option<Scalar> {
        let one = ImplicitCurve2 {
            field: self.curve.field.clone(),
            cells: vec![ImplicitCell {
                from: free,
                to: free,
                ..*self.cell
            }],
        };
        let p = one.point(0.0)?;
        Some(match self.cell.axis {
            Axis::U => p.y,
            Axis::V => p.x,
        })
    }

    fn h(&self, s: Scalar) -> Option<Scalar> {
        let free = self.free(s);
        Some(self.other.value(self.place(free, self.solved(free)?)))
    }

    /// A box certain to hold the curve for `s` in `[s0, s1]`, and bounds of
    /// the solved parameter's slope there.
    fn hull(&self, s0: Scalar, s1: Scalar) -> Option<(Cell, Range)> {
        let (f0, f1) = (self.free(s0), self.free(s1));
        let (f_lo, f_hi) = (f0.min(f1), f0.max(f1));
        let w0 = self.solved(f0)?;
        let make = |w_lo: Scalar, w_hi: Scalar| {
            let (a, b) = (self.place(f_lo, w_lo), self.place(f_hi, w_hi));
            Cell {
                lo: a.min(b),
                hi: a.max(b),
            }
        };
        let (d_free, d_solved) = match self.cell.axis {
            Axis::U => (self.d_u, self.d_v),
            Axis::V => (self.d_v, self.d_u),
        };
        let whole = make(self.cell.low, self.cell.high);
        let free = bound_simple(d_free, &whole);
        let solved = bound_simple(d_solved, &whole);
        if solved.straddles_zero() {
            return None;
        }
        let floor = solved.lo.abs().min(solved.hi.abs());
        let big = free.lo.abs().max(free.hi.abs());
        let reach = big / floor * (f_hi - f_lo);
        let hull = make(
            (w0 - reach).max(self.cell.low),
            (w0 + reach).min(self.cell.high),
        );
        // Slope bound of the solved parameter over the hull: -F_free / F_solved.
        let free = bound_simple(d_free, &hull);
        let solved = bound_simple(d_solved, &hull);
        let q = [
            -free.lo / solved.lo,
            -free.lo / solved.hi,
            -free.hi / solved.lo,
            -free.hi / solved.hi,
        ];
        let slope = Range {
            lo: q.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            hi: q.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        };
        Some((hull, slope))
    }

    fn roots(&self, s0: Scalar, s1: Scalar, depth: u32, out: &mut Vec<(Scalar, usize)>) {
        let Some((hull, slope)) = self.hull(s0, s1) else {
            return;
        };
        if !bound_simple(self.other, &hull).straddles_zero() {
            return;
        }
        let (Some(a), Some(b)) = (self.h(s0), self.h(s1)) else {
            return;
        };
        // d h / d free = H_free + H_solved * slope.
        let (h_free, h_solved) = match self.cell.axis {
            Axis::U => (self.h_u, self.h_v),
            Axis::V => (self.h_v, self.h_u),
        };
        let hf = bound_simple(h_free, &hull);
        let hs = bound_simple(h_solved, &hull);
        let p = [
            hs.lo * slope.lo,
            hs.lo * slope.hi,
            hs.hi * slope.lo,
            hs.hi * slope.hi,
        ];
        let lo = hf.lo + p.iter().copied().fold(Scalar::INFINITY, Scalar::min);
        let hi = hf.hi + p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max);
        let monotone = lo > 0.0 || hi < 0.0;
        if monotone {
            if a == 0.0 {
                out.push((s0, 1));
            } else if (a < 0.0) != (b < 0.0) {
                out.push((self.bisect(s0, s1, a), 1));
            }
            return;
        }
        if depth >= 44 || s1 - s0 <= 1e-13 {
            if (a < 0.0) != (b < 0.0) {
                out.push((self.bisect(s0, s1, a), 1));
            } else {
                // Touching: the second field reaches zero without crossing.
                let scale = 1e-9 * self.other.magnitude().max(1.0);
                if a.abs().min(b.abs()) <= scale {
                    out.push((0.5 * (s0 + s1), 2));
                }
            }
            return;
        }
        let m = 0.5 * (s0 + s1);
        self.roots(s0, m, depth + 1, out);
        self.roots(m, s1, depth + 1, out);
    }

    fn bisect(&self, mut s0: Scalar, mut s1: Scalar, a: Scalar) -> Scalar {
        let negative = a < 0.0;
        for _ in 0..100 {
            let m = 0.5 * (s0 + s1);
            if m <= s0 || m >= s1 {
                break;
            }
            match self.h(m) {
                Some(0.0) => return m,
                Some(h) if (h < 0.0) == negative => s0 = m,
                Some(_) => s1 = m,
                None => break,
            }
        }
        0.5 * (s0 + s1)
    }
}
