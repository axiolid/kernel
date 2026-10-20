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
#[derive(Debug, Clone, Copy)]
pub(super) struct CurveSpan {
    /// From the start vertex's parameter to the end vertex's, unwrapped: on
    /// a closed curve the end may lie up to one period past the domain.
    span: Interval,
    /// A closed curve's native domain, into which parameters past it wrap.
    period: Option<Interval>,
}

impl CurveSpan {
    /// The unwrapped parameter at fraction `s` of the span.
    fn parameter(&self, s: Scalar) -> Scalar {
        self.span.start + (self.span.end - self.span.start) * s
    }

    /// `t` moved into the native domain of a piece whose midpoint is `mid`.
    fn native(&self, t: Scalar, mid: Scalar) -> Scalar {
        match self.period {
            Some(domain) if mid > domain.end => t - (domain.end - domain.start),
            _ => t,
        }
    }

    /// The point at fraction `s` of the span.
    pub(super) fn evaluate(&self, curve: &Curve3, s: Scalar) -> GeomResult<Vec3> {
        let t = self.parameter(s);
        axiolid_reference::curve::evaluate3(curve, self.native(t, t))
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
        let breaks = axiolid_reference::bound::continuity_breaks3(curve, 1);
        let mut cuts: Vec<Scalar> = breaks.clone();
        if let Some(domain) = self.period {
            let period = domain.end - domain.start;
            cuts.extend(breaks.iter().map(|t| t + period));
            cuts.push(domain.end);
        }
        cuts.retain(|&t| lo < t && t < hi);
        cuts.sort_by(Scalar::total_cmp);
        cuts.dedup();
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
                axiolid_reference::bound::curve_derivative_bounds3(curve, n0, n1)?.first
                    * (n1 - n0).abs()
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
                return Ok(CurveSpan {
                    span: domain,
                    period: None,
                });
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
    Ok(CurveSpan { span, period })
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
        let span = CurveSpan {
            span: Interval {
                start: 3.5,
                end: 4.5,
            },
            period: Some(Interval {
                start: 0.0,
                end: 4.0,
            }),
        };
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
