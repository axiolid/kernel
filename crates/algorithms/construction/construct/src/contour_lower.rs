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
//! Anything else is REFUSED rather than sampled. Tessellating an ellipse or a
//! spline into short chords would produce a solid that looks right, passes
//! every closure check, and silently is not the requested shape.

use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Point2, Scalar, Tolerance};
use axiolid_curve::{Circle2, Curve2};
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_profile::{Contour, ProfileSegment};

/// Convert one closed contour into an arc ring.
///
/// Returns the ring in contour order. Segment endpoints must already meet;
/// gaps are reported rather than closed, because a silently bridged gap
/// changes the profile the caller asked for.
pub fn contour_to_arc_ring(contour: &Contour, tolerance: Tolerance) -> GeomResult<ArcRing> {
    if contour.segments.len() < 2 {
        return Err(GeomError::InvalidInput(format!(
            "a closed contour needs at least two segments, got {}",
            contour.segments.len()
        )));
    }

    let mut vertices = Vec::with_capacity(contour.segments.len());
    let mut previous_end: Option<Point2> = None;

    for segment in &contour.segments {
        let (pieces, end) = lower_segment(segment)?;
        let start = pieces[0].point;
        if let Some(previous) = previous_end {
            let gap = (start - previous).length();
            if gap > tolerance.linear() {
                return Err(GeomError::InvalidInput(format!(
                    "contour segments leave a gap of {gap} at {start:?}"
                )));
            }
        }
        vertices.extend(pieces);
        previous_end = Some(end);
    }

    // The contour must close back onto its own first vertex.
    if let (Some(last), Some(first)) = (previous_end, vertices.first()) {
        let gap = (first.point - last).length();
        if gap > tolerance.linear() {
            return Err(GeomError::InvalidInput(format!(
                "contour does not close: {gap} from last segment end to start"
            )));
        }
    }

    Ok(ArcRing::new(vertices))
}

/// One segment as the ring vertices it contributes (each with the bulge of
/// the edge leaving it) and its end point.
fn lower_segment(segment: &ProfileSegment) -> GeomResult<(Vec<ArcVertex>, Point2)> {
    let (from, to) = if segment.same_sense {
        (segment.domain.start, segment.domain.end)
    } else {
        (segment.domain.end, segment.domain.start)
    };

    match &segment.curve {
        Curve2::Line(line) => {
            let start = line.origin + line.direction * from;
            let end = line.origin + line.direction * to;
            Ok((
                vec![ArcVertex {
                    point: start,
                    bulge: 0.0,
                }],
                end,
            ))
        }
        Curve2::Circle(circle) => lower_arc(circle, from, to),
        other => Err(GeomError::UnsupportedInput {
            backend: crate::BACKEND_ID,
            operation: Operation::Sweep,
            input: unsupported_curve_name(other),
        }),
    }
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
fn lower_arc(circle: &Circle2, from: Scalar, to: Scalar) -> GeomResult<(Vec<ArcVertex>, Point2)> {
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
        return Err(GeomError::Degenerate(
            "circular profile segment has a degenerate frame".to_owned(),
        ));
    }
    let sweep = (to - from) * handedness.signum();
    if sweep == 0.0 {
        return Err(GeomError::Degenerate(
            "circular profile segment has an empty parameter range".to_owned(),
        ));
    }
    if sweep.abs() > core::f64::consts::TAU {
        return Err(GeomError::UnsupportedInput {
            backend: crate::BACKEND_ID,
            operation: Operation::Sweep,
            input: "circular profile segment sweeping more than a whole turn",
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
        Curve2::Polyline(_) => "polyline contour segment",
        Curve2::BSpline(_) => "B-spline contour segment",
        Curve2::Intrinsic(_) => "intrinsic contour segment",
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
