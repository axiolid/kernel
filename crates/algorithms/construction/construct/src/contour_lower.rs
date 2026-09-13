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
        let (start, end, bulge) = lower_segment(segment)?;
        if let Some(previous) = previous_end {
            let gap = (start - previous).length();
            if gap > tolerance.linear() {
                return Err(GeomError::InvalidInput(format!(
                    "contour segments leave a gap of {gap} at {start:?}"
                )));
            }
        }
        vertices.push(ArcVertex {
            point: start,
            bulge,
        });
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

/// One segment as (start, end, bulge of the edge leaving start).
fn lower_segment(segment: &ProfileSegment) -> GeomResult<(Point2, Point2, Scalar)> {
    let (from, to) = if segment.same_sense {
        (segment.domain.start, segment.domain.end)
    } else {
        (segment.domain.end, segment.domain.start)
    };

    match &segment.curve {
        Curve2::Line(line) => {
            let start = line.origin + line.direction * from;
            let end = line.origin + line.direction * to;
            Ok((start, end, 0.0))
        }
        Curve2::Circle(circle) => lower_arc(circle, from, to),
        other => Err(GeomError::UnsupportedInput {
            backend: crate::BACKEND_ID,
            operation: Operation::Sweep,
            input: unsupported_curve_name(other),
        }),
    }
}

/// A circular segment as start, end, and the bulge the extruder expects.
///
/// `bulge` is `tan(sweep / 4)` for the SIGNED sweep measured in world
/// orientation. The frame's handedness matters: a left-handed frame runs the
/// parameter backwards relative to the plane, so the world sweep is the
/// negation of the parameter sweep. Ignoring that flips the arc onto the
/// wrong side of its chord.
fn lower_arc(circle: &Circle2, from: Scalar, to: Scalar) -> GeomResult<(Point2, Point2, Scalar)> {
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
    // A single bulge cannot express a half turn or more: the tangent blows up
    // at pi and the chord no longer determines the arc.
    if sweep.abs() >= core::f64::consts::PI {
        return Err(GeomError::UnsupportedInput {
            backend: crate::BACKEND_ID,
            operation: Operation::Sweep,
            input: "circular profile segment sweeping half a turn or more",
        });
    }
    Ok((start, end, (sweep / 4.0).tan()))
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
