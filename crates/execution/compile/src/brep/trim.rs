//! Edge curves trimmed at their vertices (#319, #333).
//!
//! An edge runs from its start vertex to its end vertex in its curve's own
//! sense, and its vertices bound it, not its curve's range: IFC trims an
//! unbounded `IfcLine`, a whole `IfcCircle` or an `IfcBSplineCurveWithKnots`
//! at `EdgeStart` and `EdgeEnd`. A curve whose whole domain already runs from
//! the start vertex to the end vertex is the edge as it stands. Otherwise
//! each vertex is located on the curve and the edge is the span between
//! them ([`vertex_span`]):
//!
//! - A line, a circle or an ellipse locates its vertices by closed-form
//!   inversion, which checks the residual against the linear tolerance.
//! - A B-spline locates them by the certified global projection
//!   ([`axiolid_nurbs::project_curve3_certified`]), the curve analogue of
//!   [`axiolid_nurbs::invert_surface_certified`]: a vertex whose certified
//!   distance lower bound exceeds the linear tolerance is off the curve, and
//!   one whose surviving minimiser cover is not one connected interval
//!   (across the seam of a closed curve) has no unique parameter. Both are
//!   refused by name. The search runs at a thousandth of the tolerance, and
//!   since it stops as soon as its distance gap is met, its cover is then
//!   refined by certified speed bounds to a fixed parameter resolution
//!   ([`minimiser_cover`]) before uniqueness is judged, as the surface
//!   inversion judges it on cells within its parameter tolerance. Newton
//!   then sharpens the parameter inside the connected cover.
//! - A polyline locates them exactly: per segment, the closest point and
//!   the parameters within the linear tolerance of the vertex, in closed
//!   form, under the same uniqueness rule.
//! - A section curve, or any other family, is refused by name: its
//!   inversion takes the first of several turns naming the point, so the
//!   span between two vertices is not determined by it.
//!
//! On an open curve the span runs between the two parameters, whichever
//! comes first (the line's rule). On a closed curve (a circle, an ellipse, a
//! closed polyline, a B-spline whose `closed` seam
//! [`axiolid_nurbs::PeriodicCurve3`] verifies) it runs forward in the curve's
//! direction from the start vertex to the end vertex, across the seam where
//! it must, and is the whole curve from it for an edge with one vertex. Two
//! vertices at one point of a closed curve name either no span or the whole
//! curve, and are refused as ambiguous.
//!
//! # Span bounds across breaks
//!
//! A span's certified chord bound comes from the curve's own bounds
//! ([`axiolid_reference::bound::chord_bound3`]), which need the span free of
//! derivative breaks: a polyline's corners, a B-spline's full-multiplicity
//! knots, a closed curve's seam. A span over such points is split there.
//! Each piece lies within its own chord bound `d_i` of its own chord, and
//! that chord's points lie within `h`, the largest distance of a split point
//! from the span's chord, of the span's chord (the distance to a segment is
//! convex along a segment), so the span is within `max d_i + h` of its chord,
//! and its length is at most the sum of its pieces'.
//!
//! The same breaks tell the edge pass where to sample (#342):
//! [`CurveSpan::break_fractions`] gives them as fractions of the span, each
//! the float whose parameter is the break's own, so a sample there is the
//! break exactly and no span between samples crosses one.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Interval, Scalar, Tolerance, Vec3};
use axiolid_curve::{BSplineCurve3, Curve3, Polyline3};
use axiolid_nurbs::CertifiedProjectionOptions;
use axiolid_topology::VertexId;

use super::inversion::SpanBound;
use super::segment_distance;

/// Cells the certified B-spline projection may generate per vertex.
const CURVE_PROJECTION_NODES: u32 = 100_000;

/// Subdivision depth of one Bezier segment in that projection.
const CURVE_PROJECTION_DEPTH: u16 = 64;

/// The fraction of a B-spline's domain its minimiser cover is refined to.
const CURVE_PARAMETER_RESOLUTION: Scalar = 1e-6;

/// Cells that refinement may visit per vertex.
const CURVE_COVER_WORK: usize = 100_000;

/// Newton steps sharpening a certified B-spline parameter.
const POLISH_STEPS: usize = 8;

/// The part of an edge curve its vertices bound.
#[derive(Debug, Clone)]
pub(super) struct CurveSpan {
    /// From the start vertex's parameter to the end vertex's, unwrapped: on
    /// a closed curve the end may lie up to one period past the domain.
    span: Interval,
    /// A closed curve's native domain, into which parameters past it wrap.
    period: Option<Interval>,
    /// The curve's own derivative breaks, ascending: found once here, as
    /// the edge pass asks for a point and a bound many times.
    breaks: Vec<Scalar>,
}

impl CurveSpan {
    /// `span` of `curve`, wrapping into `period` on a closed curve.
    pub(super) fn new(curve: &Curve3, span: Interval, period: Option<Interval>) -> Self {
        let mut breaks = axiolid_reference::bound::continuity_breaks3(curve, 1);
        breaks.sort_by(Scalar::total_cmp);
        breaks.dedup();
        Self {
            span,
            period,
            breaks,
        }
    }

    /// The unwrapped parameter at fraction `s` of the span.
    fn parameter(&self, s: Scalar) -> Scalar {
        self.span.start + (self.span.end - self.span.start) * s
    }

    /// The parameter whose copy one period on is `t`, where that is exact:
    /// a break `b` with `b + period == t`, or the domain start for the
    /// seam `domain.end`.
    fn exact_native(&self, t: Scalar) -> Option<Scalar> {
        let domain = self.period?;
        if t == domain.end {
            return Some(domain.start);
        }
        let period = domain.end - domain.start;
        let i = self.breaks.partition_point(|b| b + period < t);
        self.breaks.get(i).copied().filter(|b| b + period == t)
    }

    /// `t` moved into the native domain, as the end of a piece whose
    /// midpoint is `mid`.
    ///
    /// Shifting by `t - period` is inexact for a non-dyadic parameter, and a
    /// piece that then ends past an interior C0 knot hides the knot inside
    /// it (#426). The cuts are the only parameters whose native value is
    /// known: a copy `b + period` of a break `b` maps back to `b`, and the
    /// seam `domain.end` to `domain.start`. Each is recognised by float
    /// equality with the very value `cuts` produced, so the decision is
    /// exact; any other parameter is shifted.
    fn native(&self, t: Scalar, mid: Scalar) -> Scalar {
        match self.period {
            Some(domain) if mid > domain.end => self
                .exact_native(t)
                .unwrap_or(t - (domain.end - domain.start)),
            _ => t,
        }
    }

    /// The native parameter of the sample at unwrapped `t`, inside the
    /// native interval of the piece [`CurveSpan::bound`] bounds around it.
    ///
    /// A sample between two cuts shifts by `t - period`, which may land a
    /// ulp outside that piece, past the knot it follows or precedes, where
    /// the piece's bound does not reach. The piece's ends are exact (as
    /// [`CurveSpan::native`] gives them) and rounding is monotone, so
    /// clamping the shift into them is exact and the bound covers the sample.
    fn sample(&self, t: Scalar) -> Scalar {
        let Some(domain) = self.period.filter(|domain| t > domain.end) else {
            return t;
        };
        if let Some(exact) = self.exact_native(t) {
            return exact;
        }
        let period = domain.end - domain.start;
        let hi = self.span.start.max(self.span.end);
        let i = self.breaks.partition_point(|b| b + period < t);
        let lower = match i.checked_sub(1).map(|j| self.breaks[j]) {
            Some(b) if b + period > domain.end => b,
            _ => domain.start,
        };
        let upper = match self.breaks.get(i) {
            Some(&b) if b + period < hi => b,
            _ => self.native(hi, hi),
        };
        (t - period).max(lower).min(upper)
    }

    /// The point at fraction `s` of the span.
    pub(super) fn evaluate(&self, curve: &Curve3, s: Scalar) -> GeomResult<Vec3> {
        let t = self.parameter(s);
        axiolid_reference::curve::evaluate3(curve, self.sample(t))
    }

    /// The unwrapped parameters strictly between `lo` and `hi` where the
    /// curve has a derivative break: its own breaks, their copies one
    /// period on, and a closed curve's seam.
    fn cuts(&self, lo: Scalar, hi: Scalar) -> Vec<Scalar> {
        let breaks = &self.breaks;
        let mut cuts: Vec<Scalar> = breaks.clone();
        if let Some(domain) = self.period {
            let period = domain.end - domain.start;
            cuts.extend(breaks.iter().map(|t| t + period));
            cuts.push(domain.end);
        }
        cuts.retain(|&t| lo < t && t < hi);
        cuts.sort_by(Scalar::total_cmp);
        cuts.dedup();
        cuts
    }

    /// The fractions of the span at which the curve has a derivative break
    /// (#342): strictly between 0 and 1, ascending, without repeats. Each
    /// is the fraction whose parameter is the break's own, so that a sample
    /// there is the break exactly and no piece between them spans one.
    pub(super) fn break_fractions(&self) -> Vec<Scalar> {
        let length = self.span.end - self.span.start;
        let (lo, hi) = if length >= 0.0 {
            (self.span.start, self.span.end)
        } else {
            (self.span.end, self.span.start)
        };
        let mut fractions: Vec<Scalar> = self
            .cuts(lo, hi)
            .into_iter()
            .map(|t| self.exact_fraction(t, (t - self.span.start) / length))
            .filter(|&s| 0.0 < s && s < 1.0)
            .collect();
        fractions.sort_by(Scalar::total_cmp);
        fractions.dedup();
        fractions
    }

    /// The float next to `guess`, within a few steps, whose parameter is
    /// `t` exactly; `guess` itself where none is.
    fn exact_fraction(&self, t: Scalar, guess: Scalar) -> Scalar {
        let (mut down, mut up) = (guess, guess);
        for _ in 0..4 {
            if self.parameter(down) == t {
                return down;
            }
            if self.parameter(up) == t {
                return up;
            }
            down = down.next_down();
            up = up.next_up();
        }
        guess
    }

    /// The certified bounds of the piece between fractions `a` and `b`,
    /// whose ends are `pa` and `pb`, split at the curve's breaks (the module
    /// notes give the argument); `None` where the family has no bound.
    pub(super) fn bound(
        &self,
        curve: &Curve3,
        a: Scalar,
        b: Scalar,
        pa: Vec3,
        pb: Vec3,
    ) -> Option<SpanBound> {
        let (ta, tb) = (self.parameter(a), self.parameter(b));
        let (lo, hi) = if ta <= tb { (ta, tb) } else { (tb, ta) };
        let cuts = self.cuts(lo, hi);
        let ends: Vec<Scalar> = core::iter::once(lo)
            .chain(cuts.iter().copied())
            .chain(core::iter::once(hi))
            .collect();
        let (mut chord, mut arc): (Scalar, Scalar) = (0.0, 0.0);
        for piece in ends.windows(2) {
            let mid = 0.5 * (piece[0] + piece[1]);
            let (n0, n1) = (self.native(piece[0], mid), self.native(piece[1], mid));
            let d = axiolid_reference::bound::chord_bound3(curve, n0, n1)?;
            arc += if d == 0.0 {
                let p0 = axiolid_reference::curve::evaluate3(curve, n0).ok()?;
                let p1 = axiolid_reference::curve::evaluate3(curve, n1).ok()?;
                (p1 - p0).length()
            } else {
                local_arc(curve, n0, n1)?
            };
            chord = chord.max(d);
        }
        let corner = cuts
            .iter()
            .map(|&cut| {
                axiolid_reference::curve::evaluate3(curve, self.native(cut, cut))
                    .map(|point| segment_distance(point, pa, pb))
            })
            .try_fold(0.0, |worst: Scalar, off| {
                off.map(|off| worst.max(if off.is_nan() { Scalar::INFINITY } else { off }))
            })
            .ok()?;
        let chord = chord + corner;
        if !(chord.is_finite() && arc.is_finite()) {
            return None;
        }
        if chord == 0.0 {
            return Some(SpanBound::straight(pa, pb));
        }
        Some(SpanBound { chord, arc })
    }
}

/// A certified bound on the arc length of `curve` between `n0` and `n1`:
/// the smaller of the speed bound over the span times its length, and the
/// speed at its middle plus the second-derivative bound over half of it,
/// `|C'(m)| h + sup|C''| h^2 / 4` with `h = |n1 - n0|` (`|C'|` is Lipschitz
/// with constant `sup|C''|`). A conic's speed bound is its largest semi-axis
/// over any span, while its second-derivative bound follows the span, so
/// the second form is the tighter one along an ellipse's flat side.
fn local_arc(curve: &Curve3, n0: Scalar, n1: Scalar) -> Option<Scalar> {
    let bounds = axiolid_reference::bound::curve_derivative_bounds3(curve, n0, n1)?;
    let h = (n1 - n0).abs();
    let global = bounds.first * h;
    let local = axiolid_reference::curve::derivative3(curve, 0.5 * (n0 + n1))
        .ok()
        .map(|speed| {
            (speed.length() * h + 0.25 * bounds.second * h * h)
                * (1.0 + axiolid_reference::bound::ROUNDING)
        })
        .filter(|local| local.is_finite());
    Some(local.map_or(global, |local| local.min(global)))
}

/// The span of `curve` an edge from `start` to `end` covers (the module
/// notes give the rules).
pub(super) fn vertex_span(
    curve: &Curve3,
    (start_id, start): (VertexId, Vec3),
    (end_id, end): (VertexId, Vec3),
    tolerance: Tolerance,
) -> GeomResult<CurveSpan> {
    let linear = tolerance.linear();
    let (t0, t1, period) = match curve {
        Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_) => {
            let locate = |vertex: Vec3, which: &str| {
                axiolid_reference::curve::invert3(curve, vertex, tolerance).map_err(|error| {
                    GeomError::InvalidInput(format!(
                        "the {which} vertex of an edge without a p-curve is not on its curve: \
                         {error}"
                    ))
                })
            };
            let period = (!matches!(curve, Curve3::Line(_))).then_some(Interval {
                start: 0.0,
                end: core::f64::consts::TAU,
            });
            (locate(start, "start")?, locate(end, "end")?, period)
        }
        _ => {
            let domain = axiolid_reference::curve::domain3(curve);
            let gap = |parameter: Scalar, vertex: Vec3| -> GeomResult<Scalar> {
                Ok((axiolid_reference::curve::evaluate3(curve, parameter)? - vertex).length())
            };
            let bounded =
                domain.start.is_finite() && domain.end.is_finite() && domain.end > domain.start;
            if bounded && gap(domain.start, start)? <= linear && gap(domain.end, end)? <= linear {
                return Ok(CurveSpan::new(curve, domain, None));
            }
            match curve {
                Curve3::BSpline(spline) => {
                    let period = axiolid_nurbs::PeriodicCurve3::try_new(spline, tolerance)
                        .ok()
                        .map(|view| {
                            let (start, end) = view.domain();
                            Interval { start, end }
                        });
                    let locate = |vertex: Vec3, which: &str| {
                        locate_on_spline(curve, spline, vertex, which, period, tolerance)
                    };
                    (locate(start, "start")?, locate(end, "end")?, period)
                }
                Curve3::Polyline(polyline) => {
                    let period = polyline.closed.then_some(domain);
                    let locate = |vertex: Vec3, which: &str| {
                        locate_on_polyline(polyline, vertex, which, period, linear)
                    };
                    (locate(start, "start")?, locate(end, "end")?, period)
                }
                _ => {
                    return Err(GeomError::InvalidInput(
                        "an edge curve without a p-curve neither runs its whole domain from the \
                         edge's start vertex to its end vertex nor is a line, conic, B-spline or \
                         polyline its vertices can be located on unambiguously"
                            .to_owned(),
                    ))
                }
            }
        }
    };
    let span = match period {
        None => Interval { start: t0, end: t1 },
        Some(domain) => {
            let turn = domain.end - domain.start;
            let sweep = if start_id == end_id {
                turn
            } else if (end - start).length() <= linear {
                return Err(GeomError::InvalidInput(
                    "an edge without a p-curve joins two vertices at one point of its closed \
                     curve: it is either no arc or the whole turn, so its span is ambiguous"
                        .to_owned(),
                ));
            } else {
                (t1 - t0).rem_euclid(turn)
            };
            Interval {
                start: t0,
                end: t0 + sweep,
            }
        }
    };
    if !(span.start.is_finite() && span.end.is_finite() && span.start != span.end) {
        return Err(GeomError::Degenerate(
            "an edge without a p-curve has no length between its vertices on its curve".to_owned(),
        ));
    }
    // A circle or ellipse evaluates past its turn by itself.
    let period = period.filter(|_| !matches!(curve, Curve3::Circle(_) | Curve3::Ellipse(_)));
    Ok(CurveSpan::new(curve, span, period))
}

fn off_curve(which: &str, distance: Scalar) -> GeomError {
    GeomError::InvalidInput(format!(
        "the {which} vertex of an edge without a p-curve is not on its curve: it lies at least \
         {distance:e} from it"
    ))
}

fn ambiguous(which: &str) -> GeomError {
    GeomError::InvalidInput(format!(
        "the {which} vertex of an edge without a p-curve has no unique parameter on its curve: \
         the curve passes through it more than once, so its span is ambiguous"
    ))
}

/// The connected interval of `intervals` holding `at`, when they form one
/// connected cover: on a closed curve's `period` the cover may run across
/// the seam, as two linear pieces touching the domain's two ends. `None`
/// when the cover has two separated parts.
fn unique_cover(
    intervals: &[(Scalar, Scalar)],
    period: Option<Interval>,
    at: Scalar,
) -> Option<(Scalar, Scalar)> {
    let mut sorted = intervals.to_vec();
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut parts: Vec<(Scalar, Scalar)> = Vec::new();
    for (lo, hi) in sorted {
        match parts.last_mut() {
            Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
            _ => parts.push((lo, hi)),
        }
    }
    let connected = match (parts.as_slice(), period) {
        ([_], _) => true,
        ([first, last], Some(domain)) => first.0 <= domain.start && last.1 >= domain.end,
        _ => false,
    };
    if !connected {
        return None;
    }
    parts
        .iter()
        .copied()
        .find(|&(lo, hi)| lo <= at && at <= hi)
        .or_else(|| parts.first().copied())
}

/// The parameter naming `vertex` on a B-spline edge curve, certified unique
/// (the module notes give the rules).
fn locate_on_spline(
    curve: &Curve3,
    spline: &BSplineCurve3,
    vertex: Vec3,
    which: &str,
    period: Option<Interval>,
    tolerance: Tolerance,
) -> GeomResult<Scalar> {
    let linear = tolerance.linear();
    let run = |gap: Scalar| -> GeomResult<axiolid_nurbs::CurveProjectionCertificate3> {
        let gap = Tolerance::new(gap, tolerance.angular()).map_err(|error| {
            GeomError::InvalidInput(format!("curve projection tolerance: {error:?}"))
        })?;
        let options =
            CertifiedProjectionOptions::new(gap, CURVE_PROJECTION_NODES, CURVE_PROJECTION_DEPTH)?;
        axiolid_nurbs::project_curve3_certified(spline, vertex, options)
    };
    let certificate = match run(linear * 1e-3) {
        Err(GeomError::BudgetExceeded { .. }) => run(linear)?,
        outcome => outcome?,
    };
    if certificate.distance_lower_bound > linear {
        return Err(off_curve(which, certificate.distance_lower_bound));
    }
    let cells: Vec<(Scalar, Scalar)> = certificate
        .possible_minimizer_intervals
        .iter()
        .map(|interval| (interval.start, interval.end))
        .collect();
    let threshold = certificate.distance_upper_bound + linear * 1e-3;
    let intervals = minimiser_cover(curve, vertex, cells, threshold)?;
    let (lo, hi) =
        unique_cover(&intervals, period, certificate.parameter).ok_or_else(|| ambiguous(which))?;
    let t = polish(curve, vertex, certificate.parameter.clamp(lo, hi), (lo, hi))?;
    let residual = (axiolid_reference::curve::evaluate3(curve, t)? - vertex).length();
    if residual.is_nan() || residual > linear {
        return Err(off_curve(which, residual));
    }
    Ok(t)
}

/// The certified search's cover refined to the parameter resolution.
///
/// The search stops once its distance gap is met, which for a point on
/// the curve is at once: its cells are then whole Bezier segments whose
/// boxes hold the point, and a segment the curve never brings near the
/// point survives beside the one that does. Each cell is therefore split
/// until it is narrower than [`CURVE_PARAMETER_RESOLUTION`] of the domain,
/// and dropped where `|C(m) - p| - L w / 2`, with `L` the certified speed
/// bound over the cell, `m` its midpoint and `w` its width, exceeds
/// `threshold`: every point of it is then farther than any minimiser.
fn minimiser_cover(
    curve: &Curve3,
    point: Vec3,
    cells: Vec<(Scalar, Scalar)>,
    threshold: Scalar,
) -> GeomResult<Vec<(Scalar, Scalar)>> {
    let domain = axiolid_reference::curve::domain3(curve);
    let resolution = (domain.end - domain.start).abs() * CURVE_PARAMETER_RESOLUTION;
    let mut pending = cells;
    let mut kept = Vec::new();
    let mut work = 0_usize;
    while let Some((a, b)) = pending.pop() {
        work += 1;
        if work > CURVE_COVER_WORK {
            return Err(GeomError::BudgetExceeded {
                resource: "certified B-spline location of an edge vertex",
            });
        }
        let mid = 0.5 * (a + b);
        let near = (axiolid_reference::curve::evaluate3(curve, mid)? - point).length();
        let lower = axiolid_reference::bound::curve_derivative_bounds3(curve, a, b)
            .map(|speed| near - speed.first * 0.5 * (b - a));
        if lower.is_some_and(|lower| lower > threshold) {
            continue;
        }
        if b - a <= resolution {
            kept.push((a, b));
        } else {
            pending.push((a, mid));
            pending.push((mid, b));
        }
    }
    Ok(kept)
}

/// Newton on `|C(t) - p|^2` inside `[lo, hi]`, keeping the best parameter.
fn polish(
    curve: &Curve3,
    point: Vec3,
    start: Scalar,
    (lo, hi): (Scalar, Scalar),
) -> GeomResult<Scalar> {
    let residual = |t: Scalar| -> GeomResult<Scalar> {
        Ok((axiolid_reference::curve::evaluate3(curve, t)? - point).length())
    };
    let (mut best, mut best_residual) = (start, residual(start)?);
    let mut current = start;
    for _ in 0..POLISH_STEPS {
        if best_residual == 0.0 {
            break;
        }
        let r = axiolid_reference::curve::evaluate3(curve, current)? - point;
        let d = axiolid_reference::curve::derivative3(curve, current)?;
        let dd = axiolid_reference::curve::second_derivative3(curve, current)?;
        let slope = d.dot(d) + r.dot(dd);
        if !(slope.is_finite() && slope > 0.0) {
            break;
        }
        current = (current - r.dot(d) / slope).clamp(lo, hi);
        let next = residual(current)?;
        if next < best_residual {
            best = current;
            best_residual = next;
        } else {
            break;
        }
    }
    Ok(best)
}

/// The parameter naming `vertex` on a polyline edge curve, located exactly:
/// the closest point over its segments, unique when the parameters within
/// `linear` of the vertex form one connected cover.
fn locate_on_polyline(
    polyline: &Polyline3,
    vertex: Vec3,
    which: &str,
    period: Option<Interval>,
    linear: Scalar,
) -> GeomResult<Scalar> {
    let points = &polyline.points;
    let count = points.len();
    let segments = if polyline.closed {
        count
    } else {
        count.saturating_sub(1)
    };
    if count < 2 || segments == 0 {
        return Err(GeomError::Degenerate(
            "an edge's polyline has no segment".to_owned(),
        ));
    }
    let mut nearest: Option<(Scalar, Scalar)> = None;
    let mut within: Vec<(Scalar, Scalar)> = Vec::new();
    for i in 0..segments {
        let (a, b) = (points[i], points[(i + 1) % count]);
        let along = b - a;
        let length = along.length_squared();
        let offset = vertex - a;
        let base = i as Scalar;
        let f = if length > 0.0 {
            (offset.dot(along) / length).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let distance = (a + along * f - vertex).length();
        if nearest.is_none_or(|(best, _)| distance < best) {
            nearest = Some((distance, base + f));
        }
        if length == 0.0 {
            if offset.length() <= linear {
                within.push((base, base + 1.0));
            }
            continue;
        }
        // |offset - f along|^2 <= linear^2, a quadratic in f.
        let centre = offset.dot(along) / length;
        let across = offset.length_squared() - offset.dot(along) * centre;
        let reach = (linear * linear - across) / length;
        if reach >= 0.0 {
            let half = reach.sqrt();
            let (lo, hi) = ((centre - half).max(0.0), (centre + half).min(1.0));
            if lo <= hi {
                within.push((base + lo, base + hi));
            }
        }
    }
    let (distance, t) = nearest.ok_or_else(|| off_curve(which, Scalar::INFINITY))?;
    if within.is_empty() || distance.is_nan() || distance > linear {
        return Err(off_curve(which, distance));
    }
    unique_cover(&within, period, t).ok_or_else(|| ambiguous(which))?;
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::Point3;

    fn vertices(count: usize) -> Vec<VertexId> {
        let mut brep = axiolid_topology::BRep::<axiolid_model::NodeId>::default();
        (0..count)
            .map(|_| {
                brep.add_vertex(axiolid_topology::Vertex {
                    position: Point3::ZERO,
                })
            })
            .collect()
    }

    fn error_text(result: GeomResult<impl core::fmt::Debug>) -> String {
        format!("{:?}", result.expect_err("refused"))
    }

    /// Over a flat ellipse's vertex the arc bound follows the short
    /// semi-axis, not the long one, and still covers the arc.
    #[test]
    fn a_span_s_arc_bound_follows_the_curve_s_local_speed() {
        let ellipse = Curve3::Ellipse(axiolid_curve::Ellipse3 {
            frame: axiolid_core::Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            semi_axis_x: 1.2,
            semi_axis_y: 0.1,
        });
        for (n0, n1) in [(-0.05, 0.05), (0.0, 0.2), (1.0, 1.6)] {
            let bound = local_arc(&ellipse, n0, n1).unwrap();
            let point = |t: Scalar| axiolid_reference::curve::evaluate3(&ellipse, t).unwrap();
            let arc: Scalar = (0..1000)
                .map(|i| {
                    let t = |k: Scalar| n0 + (n1 - n0) * k / 1000.0;
                    (point(t(i as Scalar + 1.0)) - point(t(i as Scalar))).length()
                })
                .sum();
            assert!(arc <= bound, "[{n0}, {n1}]: arc {arc} above {bound}");
            assert!(
                bound <= 1.5 * arc,
                "[{n0}, {n1}]: {bound} for an arc of {arc}"
            );
        }
    }

    /// A closed square from one vertex a third of the way along its first
    /// side breaks at its corners, seam included, each fraction naming its
    /// corner's parameter exactly (#342).
    #[test]
    fn a_closed_polyline_breaks_at_its_corners_and_its_seam() {
        let p = |x: Scalar, y: Scalar| Point3::new(x, y, 0.0);
        let square = Curve3::Polyline(Polyline3 {
            points: vec![p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)],
            closed: true,
        });
        let [a] = vertices(1)[..] else { unreachable!() };
        let start = (a, p(0.3, 0.0));
        let span = vertex_span(&square, start, start, Tolerance::MILLIMETRE).unwrap();
        let fractions = span.break_fractions();
        assert_eq!(fractions.len(), 4, "{fractions:?}");
        for (s, corner) in fractions.iter().zip([1.0, 2.0, 3.0, 4.0]) {
            assert_eq!(span.parameter(*s), corner);
            assert!(0.0 < *s && *s < 1.0);
        }
        assert!(fractions.windows(2).all(|w| w[0] < w[1]));
    }

    /// The fraction naming a break is nudged to the float whose parameter is
    /// the break's own where plain division misses it by a rounding.
    #[test]
    fn a_break_fraction_names_its_parameter_exactly() {
        let p = |x: Scalar, y: Scalar| Point3::new(x, y, 0.0);
        let bent = Curve3::Polyline(Polyline3 {
            points: vec![p(0.0, 0.0), p(1.0, 1.0), p(2.0, 0.0)],
            closed: false,
        });
        let span = CurveSpan::new(
            &bent,
            Interval {
                start: 0.05808157514973589,
                end: 1.501464907695806,
            },
            None,
        );
        let plain = (1.0 - span.span.start) / (span.span.end - span.span.start);
        assert_ne!(
            span.parameter(plain),
            1.0,
            "the fixture needs a rounding miss"
        );
        let fractions = span.break_fractions();
        assert_eq!(fractions.len(), 1);
        assert_eq!(span.parameter(fractions[0]), 1.0);
    }

    /// An open polyline run backwards between two interior vertices has its
    /// fractions counted along the edge, ascending.
    #[test]
    fn a_reversed_polyline_span_breaks_in_edge_order() {
        let p = |x: Scalar, y: Scalar| Point3::new(x, y, 0.0);
        let zigzag = Curve3::Polyline(Polyline3 {
            points: vec![p(0.0, 0.0), p(1.0, 1.0), p(2.0, 0.0), p(3.0, 1.0)],
            closed: false,
        });
        let [a, b] = vertices(2)[..] else {
            unreachable!()
        };
        let span = vertex_span(
            &zigzag,
            (a, p(2.5, 0.5)),
            (b, p(0.5, 0.5)),
            Tolerance::MILLIMETRE,
        )
        .unwrap();
        let fractions = span.break_fractions();
        assert_eq!(fractions.len(), 2, "{fractions:?}");
        assert_eq!(span.parameter(fractions[0]), 2.0);
        assert_eq!(span.parameter(fractions[1]), 1.0);
        assert!(fractions[0] < fractions[1]);
        // A line has none.
        let line = Curve3::Polyline(Polyline3 {
            points: vec![p(0.0, 0.0), p(1.0, 0.0)],
            closed: false,
        });
        let span = vertex_span(
            &line,
            (a, p(0.0, 0.0)),
            (b, p(1.0, 0.0)),
            Tolerance::MILLIMETRE,
        )
        .unwrap();
        assert!(span.break_fractions().is_empty());
    }

    /// Two segments crossing at the vertex: the parameter is not unique.
    #[test]
    fn a_vertex_where_a_curve_crosses_itself_is_ambiguous() {
        let corners = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let crossing = Point3::new(0.5, 0.5, 0.0);
        let polyline = Polyline3 {
            points: corners.clone(),
            closed: false,
        };
        let text = error_text(locate_on_polyline(&polyline, crossing, "start", None, 1e-6));
        assert!(text.contains("no unique parameter"), "{text}");
        let spline = BSplineCurve3 {
            degree: 1,
            control_points: corners,
            knots: vec![0.0, 1.0, 2.0, 3.0],
            multiplicities: vec![2, 1, 1, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: axiolid_curve::KnotSpec::Unspecified,
        };
        let curve = Curve3::BSpline(spline.clone());
        let text = error_text(locate_on_spline(
            &curve,
            &spline,
            crossing,
            "start",
            None,
            Tolerance::MILLIMETRE,
        ));
        assert!(text.contains("no unique parameter"), "{text}");
        // Away from the crossing the same curve locates exactly.
        let t = locate_on_spline(
            &curve,
            &spline,
            Point3::new(0.25, 0.25, 0.0),
            "start",
            None,
            Tolerance::MILLIMETRE,
        )
        .unwrap();
        assert!((t - 0.25).abs() < 1e-12, "{t}");
    }

    /// The certified search names a B-spline vertex to its cell; the polish
    /// brings it to rounding.
    #[test]
    fn a_bspline_vertex_parameter_is_polished() {
        let w = core::f64::consts::FRAC_1_SQRT_2;
        let spline = BSplineCurve3 {
            degree: 2,
            control_points: vec![
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            knots: vec![0.0, 1.0],
            multiplicities: vec![3, 3],
            weights: Some(vec![1.0, w, 1.0]),
            closed: false,
            self_intersect: None,
            knot_spec: axiolid_curve::KnotSpec::Unspecified,
        };
        let curve = Curve3::BSpline(spline.clone());
        let angle: Scalar = 0.7;
        let vertex = Point3::new(angle.cos(), angle.sin(), 0.0);
        let t =
            locate_on_spline(&curve, &spline, vertex, "end", None, Tolerance::MILLIMETRE).unwrap();
        let on = axiolid_reference::curve::evaluate3(&curve, t).unwrap();
        assert!((on - vertex).length() < 1e-12, "{t}: {on:?}");
    }

    /// A polyline vertex is located on its segment in closed form.
    #[test]
    fn a_polyline_vertex_is_located_exactly() {
        let polyline = Polyline3 {
            points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
            ],
            closed: false,
        };
        let t =
            locate_on_polyline(&polyline, Point3::new(2.0, 0.5, 0.0), "end", None, 1e-6).unwrap();
        assert_eq!(t, 1.25);
        // At a corner both segments name it, as one connected cover.
        let t =
            locate_on_polyline(&polyline, Point3::new(2.0, 0.0, 0.0), "end", None, 1e-6).unwrap();
        assert_eq!(t, 1.0);
        let text = error_text(locate_on_polyline(
            &polyline,
            Point3::new(1.0, 0.1, 0.0),
            "end",
            None,
            1e-6,
        ));
        assert!(text.contains("is not on its curve"), "{text}");
    }

    /// A sample past the seam lies on its own side of every knot copy, and
    /// a knot copy itself is the knot exactly (#426).
    #[test]
    fn a_sample_past_the_seam_stays_in_its_piece() {
        let knots = [0.1, 0.2, 1.0 / 3.0];
        let span = CurveSpan {
            span: Interval {
                start: 0.9,
                end: 1.9,
            },
            period: Some(Interval {
                start: 0.0,
                end: 1.0,
            }),
            breaks: knots.to_vec(),
        };
        for knot in knots {
            let copy = knot + 1.0;
            assert_ne!(copy - 1.0, knot, "the fixture needs an inexact shift");
            assert_eq!(span.sample(copy), knot, "the copy of {knot}");
            let (mut below, mut above) = (copy, copy);
            for _ in 0..4 {
                below = below.next_down();
                above = above.next_up();
                assert!(
                    span.sample(below) <= knot,
                    "{below} before the copy of {knot}"
                );
                assert!(
                    span.sample(above) >= knot,
                    "{above} after the copy of {knot}"
                );
            }
        }
        let mut t = 1.0;
        assert_eq!(span.sample(t), 1.0, "the seam is the domain end itself");
        for _ in 0..4 {
            t = t.next_up();
            assert!(span.sample(t) >= 0.0, "{t} after the seam");
        }
        assert_eq!(span.sample(0.9), 0.9, "before the seam");
    }

    /// A span across a closed polyline's seam is split there: the corner
    /// at the seam bounds how far the span strays from its chord.
    #[test]
    fn a_span_across_the_seam_is_bounded_by_the_seam_corner() {
        let square = Curve3::Polyline(Polyline3 {
            points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            closed: true,
        });
        let span = CurveSpan::new(
            &square,
            Interval {
                start: 3.5,
                end: 4.5,
            },
            Some(Interval {
                start: 0.0,
                end: 4.0,
            }),
        );
        let (a, b) = (
            span.evaluate(&square, 0.0).unwrap(),
            span.evaluate(&square, 1.0).unwrap(),
        );
        assert!((a - Point3::new(0.0, 0.5, 0.0)).length() < 1e-15, "{a:?}");
        assert!((b - Point3::new(0.5, 0.0, 0.0)).length() < 1e-15, "{b:?}");
        let bound = span.bound(&square, 0.0, 1.0, a, b).unwrap();
        let corner = 0.5 * core::f64::consts::FRAC_1_SQRT_2;
        assert!((bound.chord - corner).abs() < 1e-12, "{bound:?}");
        assert!((bound.arc - 1.0).abs() < 1e-12, "{bound:?}");
    }

    /// A cover touching both ends of a closed curve's domain is connected
    /// across its seam; two separate parts are not.
    #[test]
    fn a_cover_is_connected_across_a_closed_curve_seam() {
        let domain = Interval {
            start: 0.0,
            end: 4.0,
        };
        let seam = [(3.9, 4.0), (0.0, 0.1)];
        assert_eq!(unique_cover(&seam, Some(domain), 0.05), Some((0.0, 0.1)));
        assert_eq!(unique_cover(&seam, None, 0.05), None);
        assert_eq!(
            unique_cover(&[(1.0, 1.1), (2.0, 2.1)], Some(domain), 1.0),
            None
        );
    }

    /// Two vertices at one point of a closed polyline name no span or the
    /// whole loop; one vertex names the whole loop from it.
    #[test]
    fn two_vertices_at_one_point_of_a_closed_polyline_are_ambiguous() {
        let square = Curve3::Polyline(Polyline3 {
            points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            closed: true,
        });
        let ids = vertices(2);
        let point = Point3::new(0.5, 0.0, 0.0);
        let text = error_text(vertex_span(
            &square,
            (ids[0], point),
            (ids[1], point),
            Tolerance::MILLIMETRE,
        ));
        assert!(text.contains("ambiguous"), "{text}");
        let whole = vertex_span(
            &square,
            (ids[0], point),
            (ids[0], point),
            Tolerance::MILLIMETRE,
        )
        .unwrap();
        assert_eq!((whole.span.start, whole.span.end), (0.5, 4.5));
    }
}
