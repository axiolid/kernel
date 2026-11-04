//! Lower an arbitrary contour profile onto the arc-ring extruder (ADR 0053).
//!
//! # Why this is an adapter, not new geometry
//!
//! [`ContourProfile`](axiolid_profile::ContourProfile) carries segments of
//! any [`Curve2`] kind. The arc-ring
//! extruder already builds exact walls for the two kinds that matter here --
//! `Line2` becomes a planar wall, `Circle2` becomes a cylindrical one -- so
//! lowering a contour to an [`ArcRing`] is a translation, not a construction.
//!
//! A `Polyline2` segment is exact too: one straight vertex per polyline
//! vertex it passes (see `lower_polyline`).
//!
//! Anything else is REFUSED rather than sampled. Tessellating an ellipse or a
//! spline into short chords would produce a solid that looks right, passes
//! every closure check, and silently is not the requested shape.
//!
//! # Closure at `Tolerance::ZERO` (#250)
//!
//! A ring is a cycle of vertices, each carrying the bulge of the edge that
//! leaves it, so it is closed by construction: the vertex at a joint is ONE
//! value, the end of the edge entering it and the start of the edge leaving
//! it. It is the leaving segment's start -- a line's stored origin bit for
//! bit, an arc's circle at its first parameter -- and nothing is averaged
//! or moved onto it.
//!
//! The segments themselves cannot meet exactly in floating point, however
//! carefully the producer shares its points (the section router computes
//! every corner and tangent point once and hands that one value to both
//! segments). A line stored as `origin + t direction` ends at
//! `origin + direction`, which rounds, and its stored direction is already
//! the rounded difference of the two points; an arc ends at the `cos` and
//! `sin` of a transcendental sweep. No choice of `Line2` or `Circle2` makes
//! both ends of every segment bit-exact. So the joint check compares the
//! entering segment's evaluated end with the leaving segment's start
//! against the larger of the caller's tolerance and the rounding of those
//! two evaluations: eight machine epsilons of the magnitudes they are
//! computed from, a few ulps of the coordinates. That is what lets a
//! section with non-dyadic sizes (IPE 300 left 2.6e-18 at the flange tip)
//! lower at `Tolerance::ZERO`, while a contour open by more than its own
//! rounding is still refused there. The check only decides WHETHER a ring
//! is built, never what it is: a ring that lowered before is bit-identical.

use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Point2, Scalar, Tolerance};
use axiolid_curve::{Circle2, Curve2, Polyline2};
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_profile::{Contour, ProfileSegment};

/// Convert one closed contour into an arc ring.
///
/// Returns the ring in contour order, one vertex per joint. Segment
/// endpoints must already meet, within `tolerance` or within the rounding of
/// their own evaluation, whichever is larger (see the module notes on
/// `Tolerance::ZERO`); gaps are reported rather than closed, because a
/// silently bridged gap changes the profile the caller asked for.
pub fn contour_to_arc_ring(contour: &Contour, tolerance: Tolerance) -> GeomResult<ArcRing> {
    lower_contour(contour, tolerance)
        .map(|lowered| lowered.ring)
        .map_err(ContourFault::into_geom_error)
}

/// A contour lowered to its arc ring, with what a validator needs to map the
/// ring back to the contour (#352).
#[derive(Debug, Clone)]
pub(crate) struct LoweredContour {
    pub(crate) ring: ArcRing,
    /// The contour segment each ring edge came from: a segment sweeping half
    /// a turn or more contributes several consecutive edges.
    pub(crate) segment_of_edge: Vec<usize>,
    /// The largest joint gap the lowering closed, closing joint included.
    pub(crate) max_gap: Scalar,
}

/// Why a contour does not lower, with the segment at fault.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ContourFault {
    /// A ring needs two vertices; `count` segments lowered to `vertices`.
    TooFewSegments {
        count: usize,
        vertices: usize,
    },
    UnsupportedSegment {
        segment: usize,
        kind: &'static str,
    },
    DegenerateSegment {
        segment: usize,
        reason: &'static str,
    },
    /// The end of `segment` and the start of the next are `gap` apart, beyond
    /// the slack. `closing` is the joint from the last segment to the first.
    Gap {
        segment: usize,
        gap: Scalar,
        at: Point2,
        closing: bool,
    },
}

impl ContourFault {
    /// The error `contour_to_arc_ring` has always returned for this fault.
    fn into_geom_error(self) -> GeomError {
        match self {
            Self::TooFewSegments { count, vertices } => GeomError::InvalidInput(format!(
                "a closed contour needs at least two ring vertices, its {count} segment(s) give {vertices}"
            )),
            Self::UnsupportedSegment { kind, .. } => GeomError::UnsupportedInput {
                backend: crate::BACKEND_ID,
                operation: Operation::Sweep,
                input: kind,
            },
            Self::DegenerateSegment { reason, .. } => GeomError::Degenerate(reason.to_owned()),
            Self::Gap {
                gap,
                at,
                closing: false,
                ..
            } => {
                GeomError::InvalidInput(format!("contour segments leave a gap of {gap} at {at:?}"))
            }
            Self::Gap {
                gap, closing: true, ..
            } => GeomError::InvalidInput(format!(
                "contour does not close: {gap} from last segment end to start"
            )),
        }
    }
}

/// Lower one closed contour, stopping at the first fault.
pub(crate) fn lower_contour(
    contour: &Contour,
    tolerance: Tolerance,
) -> Result<LoweredContour, ContourFault> {
    if contour.segments.is_empty() {
        return Err(ContourFault::TooFewSegments {
            count: 0,
            vertices: 0,
        });
    }

    let mut vertices = Vec::with_capacity(contour.segments.len());
    let mut segment_of_edge = Vec::with_capacity(contour.segments.len());
    let mut previous_end: Option<(Point2, Scalar)> = None;
    let mut first_scale = 0.0;
    let mut max_gap: Scalar = 0.0;

    for (index, segment) in contour.segments.iter().enumerate() {
        let lowered = lower_segment(segment, index)?;
        let start = lowered.vertices[0].point;
        if let Some((previous, previous_scale)) = previous_end {
            let gap = (start - previous).length();
            if gap > joint_slack(tolerance, previous_scale, lowered.start_scale) {
                return Err(ContourFault::Gap {
                    segment: index - 1,
                    gap,
                    at: start,
                    closing: false,
                });
            }
            max_gap = max_gap.max(gap);
        } else {
            first_scale = lowered.start_scale;
        }
        segment_of_edge.extend(std::iter::repeat_n(index, lowered.vertices.len()));
        vertices.extend(lowered.vertices);
        previous_end = Some((lowered.end, lowered.end_scale));
    }

    // A ring is a cycle of at least two vertices. Segments are not counted:
    // one closed polyline or one full circle is a whole ring by itself.
    if vertices.len() < 2 {
        return Err(ContourFault::TooFewSegments {
            count: contour.segments.len(),
            vertices: vertices.len(),
        });
    }

    // The contour must close back onto its own first vertex.
    if let (Some((last, last_scale)), Some(first)) = (previous_end, vertices.first()) {
        let gap = (first.point - last).length();
        if gap > joint_slack(tolerance, last_scale, first_scale) {
            return Err(ContourFault::Gap {
                segment: contour.segments.len() - 1,
                gap,
                at: first.point,
                closing: true,
            });
        }
        max_gap = max_gap.max(gap);
    }

    Ok(LoweredContour {
        ring: ArcRing::new(vertices),
        segment_of_edge,
        max_gap,
    })
}

/// Rounding of one evaluated segment endpoint, per unit of the magnitudes
/// it is computed from.
///
/// A line's end `origin + direction t` rounds twice (product and sum), and
/// its stored direction is already the rounded difference of the two points
/// the producer shared. An arc's end `centre + x r cos + y r sin` adds the
/// `sin`/`cos` error and frame axes that are unit only to an ulp or two
/// after normalisation. Each step is at most half an ulp of what it
/// touches, so eight machine epsilons (~1.8e-15) of the operand magnitudes
/// bound one side of a joint with room; a contour placed by an affine map
/// adds the map's own few roundings and stays inside it. A joint closer
/// than that is the representation's rounding, not a gap anyone drew.
const JOINT_ROUNDING: Scalar = 8.0 * Scalar::EPSILON;

/// The largest gap a joint may show: the caller's tolerance, or the
/// rounding of the two endpoint evaluations that meet there when that is
/// larger (at `Tolerance::ZERO`, or a tolerance below an ulp of the
/// coordinates).
pub(crate) fn joint_slack(tolerance: Tolerance, entering: Scalar, leaving: Scalar) -> Scalar {
    tolerance
        .linear()
        .max(JOINT_ROUNDING * (entering + leaving))
}

/// The magnitudes one segment's oriented end (`at_end`) or start is evaluated
/// from, which bound its rounding (see `JOINT_ROUNDING`); `point` is the
/// flattened endpoint, the fallback for a family evaluated otherwise.
pub(crate) fn endpoint_scale(segment: &ProfileSegment, at_end: bool, point: Point2) -> Scalar {
    let t = if at_end == segment.same_sense {
        segment.domain.end
    } else {
        segment.domain.start
    };
    match &segment.curve {
        Curve2::Line(line) => line.origin.length() + line.direction.length() * t.abs(),
        Curve2::Circle(circle) => circle.frame.origin.length() + circle.radius.abs(),
        _ => point.x.abs() + point.y.abs(),
    }
}

/// One segment lowered: the ring vertices it contributes (each with the
/// bulge of the edge leaving it), its evaluated end point, and the
/// magnitudes its start and end are evaluated from, which bound their
/// rounding (see `JOINT_ROUNDING`).
struct LoweredSegment {
    vertices: Vec<ArcVertex>,
    end: Point2,
    start_scale: Scalar,
    end_scale: Scalar,
}

/// One segment as the ring vertices it contributes and its end point.
fn lower_segment(segment: &ProfileSegment, index: usize) -> Result<LoweredSegment, ContourFault> {
    let (from, to) = if segment.same_sense {
        (segment.domain.start, segment.domain.end)
    } else {
        (segment.domain.end, segment.domain.start)
    };

    match &segment.curve {
        Curve2::Line(line) => {
            let start = line.origin + line.direction * from;
            let end = line.origin + line.direction * to;
            let origin = line.origin.length();
            let direction = line.direction.length();
            Ok(LoweredSegment {
                vertices: vec![ArcVertex {
                    point: start,
                    bulge: 0.0,
                }],
                end,
                start_scale: origin + direction * from.abs(),
                end_scale: origin + direction * to.abs(),
            })
        }
        Curve2::Circle(circle) => {
            let (vertices, end) = lower_arc(circle, from, to, index)?;
            let scale = circle.frame.origin.length() + circle.radius.abs();
            Ok(LoweredSegment {
                vertices,
                end,
                start_scale: scale,
                end_scale: scale,
            })
        }
        Curve2::Polyline(polyline) => lower_polyline(polyline, from, to, index),
        other => Err(ContourFault::UnsupportedSegment {
            segment: index,
            kind: unsupported_curve_name(other),
        }),
    }
}

/// A polyline segment as one straight ring vertex per polyline vertex it
/// passes, and its end point.
///
/// The parameter is one unit per polyline segment, as in the evaluator: the
/// domain is `[0, segments]` with `segments` the point count (closed) or one
/// less (open). At an integer parameter the stored point is taken bit for
/// bit, so a ring built from a polyline carries exactly its points; between
/// two integers the point is interpolated. The segment emits its start, then
/// every integer parameter strictly between `from` and `to` in traversal
/// order, all with bulge zero. A domain outside the points is refused, not
/// clamped.
fn lower_polyline(
    polyline: &Polyline2,
    from: Scalar,
    to: Scalar,
    index: usize,
) -> Result<LoweredSegment, ContourFault> {
    let degenerate = |reason: &'static str| ContourFault::DegenerateSegment {
        segment: index,
        reason,
    };
    let points = &polyline.points;
    let count = points.len();
    if count < 2 || points.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
        return Err(degenerate(
            "polyline profile segment has fewer than two finite points",
        ));
    }
    if from == to {
        return Err(degenerate(
            "polyline profile segment has an empty parameter range",
        ));
    }
    let segments = if polyline.closed { count } else { count - 1 } as Scalar;
    let inside = |t: Scalar| (0.0..=segments).contains(&t);
    if !inside(from) || !inside(to) {
        return Err(degenerate(
            "polyline profile segment's domain lies outside its points",
        ));
    }

    let evaluate = |t: Scalar| {
        let i = t.floor();
        if i == t {
            return points[i as usize % count];
        }
        let (a, b) = (points[i as usize], points[(i as usize + 1) % count]);
        a + (b - a) * (t - i)
    };
    let start = evaluate(from);
    let end = evaluate(to);

    let mut vertices = vec![ArcVertex {
        point: start,
        bulge: 0.0,
    }];
    let mut push = |t: Scalar| {
        let point = evaluate(t);
        let last = vertices[vertices.len() - 1].point;
        if point.x.to_bits() != last.x.to_bits() || point.y.to_bits() != last.y.to_bits() {
            vertices.push(ArcVertex { point, bulge: 0.0 });
        }
    };
    if from < to {
        let mut k = from.floor() + 1.0;
        while k < to {
            push(k);
            k += 1.0;
        }
    } else {
        let mut k = from.ceil() - 1.0;
        while k > to {
            push(k);
            k -= 1.0;
        }
    }
    // A closed polyline may repeat its first point; the segment's end then
    // equals its last interior vertex, and the next segment supplies it.
    if vertices.len() > 1 {
        let last = vertices[vertices.len() - 1].point;
        if last.x.to_bits() == end.x.to_bits() && last.y.to_bits() == end.y.to_bits() {
            vertices.pop();
        }
    }

    Ok(LoweredSegment {
        vertices,
        end,
        start_scale: start.x.abs() + start.y.abs(),
        end_scale: end.x.abs() + end.y.abs(),
    })
}

/// A circular segment as ring vertices with the bulges the extruder
/// expects, and its end point.
///
/// `bulge` is `tan(sweep / 4)` for the SIGNED sweep measured in world
/// orientation. The frame's handedness matters: a left-handed frame runs the
/// parameter backwards relative to the plane, so the world sweep is the
/// negation of the parameter sweep. Ignoring that flips the arc onto the
/// wrong side of its chord.
///
/// A bulge determines its arc only below half a turn, so a segment sweeping
/// half a turn or more is split into `k` equal sub-arcs, each below half a
/// turn: an IFC arch is commonly one semicircle (ADR 0053, amended by
/// #228). Each split vertex is the circle evaluated at its parameter,
/// exactly as the segment's own ends are, so every sub-arc carries the same
/// circle to the same rounding an unsplit arc does (a few ulps of the
/// circle's coordinates); nothing is fitted. A segment sweeping more than a
/// whole turn overlaps itself and is refused.
fn lower_arc(
    circle: &Circle2,
    from: Scalar,
    to: Scalar,
    index: usize,
) -> Result<(Vec<ArcVertex>, Point2), ContourFault> {
    let evaluate = |t: Scalar| {
        let (sin, cos) = t.sin_cos();
        circle.frame.origin
            + circle.frame.x * (circle.radius * cos)
            + circle.frame.y * (circle.radius * sin)
    };
    let start = evaluate(from);
    let end = evaluate(to);

    let handedness = circle.frame.x.perp_dot(circle.frame.y);
    if handedness == 0.0 {
        return Err(ContourFault::DegenerateSegment {
            segment: index,
            reason: "circular profile segment has a degenerate frame",
        });
    }
    let sweep = (to - from) * handedness.signum();
    if sweep == 0.0 {
        return Err(ContourFault::DegenerateSegment {
            segment: index,
            reason: "circular profile segment has an empty parameter range",
        });
    }
    if sweep.abs() > core::f64::consts::TAU {
        return Err(ContourFault::UnsupportedSegment {
            segment: index,
            kind: "circular profile segment sweeping more than a whole turn",
        });
    }
    // Equal sub-arcs, each below half a turn.
    let pieces = (sweep.abs() / core::f64::consts::PI).floor() as usize + 1;
    let bulge = (sweep / (4.0 * pieces as Scalar)).tan();
    let vertices = (0..pieces)
        .map(|k| ArcVertex {
            point: if k == 0 {
                start
            } else {
                evaluate(from + (to - from) * k as Scalar / pieces as Scalar)
            },
            bulge,
        })
        .collect();
    Ok((vertices, end))
}

fn unsupported_curve_name(curve: &Curve2) -> &'static str {
    match curve {
        Curve2::Ellipse(_) => "elliptical contour segment",
        Curve2::BSpline(_) => "B-spline contour segment",
        Curve2::Intrinsic(_) => "intrinsic contour segment",
        Curve2::Chain(_) => "arc-length chain contour segment",
        _ => "contour segment of an unsupported curve kind",
    }
}

/// Signed area of a bulge-encoded ring; positive when counter-clockwise.
///
/// The polygon shoelace over the vertices plus each arc's circular-segment
/// term, so a ring whose shape is decided by its arcs is classified by its
/// actual geometry rather than by its chords. A crescent thin enough that its
/// chord polygon winds the other way would otherwise be misread.
pub fn arc_ring_signed_area(ring: &ArcRing) -> Scalar {
    let count = ring.vertices.len();
    let mut total = 0.0;
    for index in 0..count {
        let here = ring.vertices[index];
        let next = ring.vertices[(index + 1) % count];
        total += here.point.perp_dot(next.point);
        if here.bulge != 0.0 {
            let sweep = 4.0 * here.bulge.atan();
            let chord = (next.point - here.point).length();
            let half = (sweep.abs() / 2.0).sin();
            if half > 0.0 {
                let radius = chord / (2.0 * half);
                total += radius * radius * (sweep - sweep.sin());
            }
        }
    }
    total / 2.0
}

/// Reverse an arc ring in place, preserving its arcs.
///
/// Reversing the vertex order alone is NOT enough: a bulge belongs to the
/// edge LEAVING its vertex, so after reversal each vertex must take the
/// negated bulge of what was previously its predecessor. Getting this wrong
/// flips every arc to the wrong side while the ring still closes.
pub fn reverse_arc_ring(ring: &ArcRing) -> ArcRing {
    let count = ring.vertices.len();
    let mut vertices = Vec::with_capacity(count);
    for index in (0..count).rev() {
        let previous = (index + count - 1) % count;
        vertices.push(ArcVertex {
            point: ring.vertices[index].point,
            bulge: -ring.vertices[previous].bulge,
        });
    }
    ArcRing { vertices }
}

/// Force a ring to the winding a boundary role requires.
///
/// Outer boundaries run counter-clockwise and holes run clockwise, so that a
/// ring's wall normals point out of the material in both cases.
pub fn orient_arc_ring(ring: &ArcRing, counter_clockwise: bool) -> GeomResult<ArcRing> {
    let area = arc_ring_signed_area(ring);
    if area == 0.0 {
        return Err(GeomError::Degenerate(
            "contour ring encloses no area".to_owned(),
        ));
    }
    if (area > 0.0) == counter_clockwise {
        Ok(ring.clone())
    } else {
        Ok(reverse_arc_ring(ring))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::{Frame2, Interval, Vec2};
    use axiolid_curve::Line2;

    fn arc(from: Scalar, to: Scalar) -> ProfileSegment {
        ProfileSegment {
            curve: Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::new(0.0, 0.0),
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                radius: 1.0,
            }),
            domain: Interval::new(from, to),
            same_sense: true,
        }
    }

    #[test]
    fn the_sub_edges_of_a_split_arc_map_to_their_segment() {
        let chord = ProfileSegment {
            curve: Curve2::Line(Line2 {
                origin: Point2::new(0.0, -1.0),
                direction: Point2::new(1.0, 1.0) - Point2::new(0.0, 0.0),
            }),
            domain: Interval::UNIT,
            same_sense: true,
        };
        // Three quarters of a turn is two sub-arcs; the line closes it.
        let contour = Contour::new(vec![arc(0.0, 1.5 * core::f64::consts::PI), chord]);
        let lowered = lower_contour(&contour, Tolerance::METRE).expect("lowers");
        assert_eq!(lowered.segment_of_edge, vec![0, 0, 1]);
        assert_eq!(lowered.ring.vertices.len(), 3);
    }
}
