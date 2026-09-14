//! Exact centre-line offsetting (ADR 0056).
//!
//! `center_line.rs` offsets the FLATTENED path, which is correct for the
//! tessellating pipeline but bakes a chord budget into the result. The exact
//! extruder cannot accept that: flattening here would produce a solid that
//! closes, audits clean, and is not the requested shape.
//!
//! # What can be offset exactly
//!
//! Only curve kinds whose offset is the SAME kind, measured rather than
//! assumed:
//!
//! - a line offsets to a parallel line;
//! - a circular arc offsets to a CONCENTRIC arc of radius `r -/+ d`, keeping
//!   its centre, frame and parameter domain, so endpoints correspond at equal
//!   parameters;
//! - an ellipse does NOT offset to an ellipse (best-fit residual 0.074 for a
//!   3:1 ellipse at distance 0.3), and neither do splines in general.
//!
//! The last group is refused. Sampling it into chords is exactly the silent
//! approximation this path exists to avoid.

use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Interval, Point2, Scalar, Tolerance, Vec2};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_profile::{CenterLineProfile, Contour, ContourProfile, ProfileSegment};

use crate::BACKEND_ID;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// One path segment reduced to what offsetting needs.
#[derive(Debug, Clone, Copy)]
enum Piece {
    Line {
        from: Point2,
        to: Point2,
    },
    Arc {
        centre: Point2,
        x: Vec2,
        y: Vec2,
        radius: Scalar,
        start: Scalar,
        end: Scalar,
        /// Signed sweep in WORLD orientation, not parameter orientation.
        sweep: Scalar,
    },
}

impl Piece {
    fn start_point(&self) -> Point2 {
        match self {
            Piece::Line { from, .. } => *from,
            Piece::Arc {
                centre,
                x,
                y,
                radius,
                start,
                ..
            } => conic(*centre, *x, *y, *radius, *start),
        }
    }

    fn end_point(&self) -> Point2 {
        match self {
            Piece::Line { to, .. } => *to,
            Piece::Arc {
                centre,
                x,
                y,
                radius,
                end,
                ..
            } => conic(*centre, *x, *y, *radius, *end),
        }
    }

    /// Unit tangent at the start, in traversal direction.
    fn start_tangent(&self) -> Vec2 {
        match self {
            Piece::Line { from, to } => (*to - *from).normalize_or_zero(),
            Piece::Arc {
                x, y, start, end, ..
            } => arc_tangent(*x, *y, *start, *end > *start),
        }
    }

    /// Unit tangent at the end, in traversal direction.
    fn end_tangent(&self) -> Vec2 {
        match self {
            Piece::Line { from, to } => (*to - *from).normalize_or_zero(),
            Piece::Arc {
                x, y, start, end, ..
            } => arc_tangent(*x, *y, *end, *end > *start),
        }
    }
}

fn conic(centre: Point2, x: Vec2, y: Vec2, radius: Scalar, t: Scalar) -> Point2 {
    let (s, c) = t.sin_cos();
    centre + x * (radius * c) + y * (radius * s)
}

fn arc_tangent(x: Vec2, y: Vec2, t: Scalar, increasing: bool) -> Vec2 {
    let (s, c) = t.sin_cos();
    let d = x * -s + y * c;
    let d = if increasing { d } else { -d };
    d.normalize_or_zero()
}
/// Read one path segment, refusing kinds with no exact offset.
fn read_piece(segment: &ProfileSegment) -> GeomResult<Piece> {
    let (from, to) = if segment.same_sense {
        (segment.domain.start, segment.domain.end)
    } else {
        (segment.domain.end, segment.domain.start)
    };
    match &segment.curve {
        Curve2::Line(line) => Ok(Piece::Line {
            from: line.origin + line.direction * from,
            to: line.origin + line.direction * to,
        }),
        Curve2::Circle(circle) => {
            let handedness = circle.frame.x.perp_dot(circle.frame.y);
            if handedness == 0.0 {
                return Err(GeomError::Degenerate(
                    "centre-line arc has a degenerate frame".to_owned(),
                ));
            }
            Ok(Piece::Arc {
                centre: circle.frame.origin,
                x: circle.frame.x,
                y: circle.frame.y,
                radius: circle.radius,
                start: from,
                end: to,
                // World sweep, not parameter sweep: a left-handed frame runs
                // the parameter backwards, and the offset SIDE depends on the
                // world turn direction.
                sweep: (to - from) * handedness.signum(),
            })
        }
        Curve2::Ellipse(_) => Err(unsupported(
            "centre-line elliptical segment has no elliptical offset",
        )),
        _ => Err(unsupported(
            "centre-line segment of a kind with no exact offset",
        )),
    }
}

/// Offset one piece sideways by a signed distance.
///
/// `distance` is positive to the LEFT of the traversal direction.
fn offset_piece(piece: &Piece, distance: Scalar) -> GeomResult<ProfileSegment> {
    match piece {
        Piece::Line { from, to } => {
            let along = (*to - *from).normalize_or_zero();
            if along == Vec2::ZERO {
                return Err(GeomError::Degenerate(
                    "centre line has a zero-length segment".to_owned(),
                ));
            }
            let normal = Vec2::new(-along.y, along.x);
            let start = *from + normal * distance;
            let end = *to + normal * distance;
            Ok(ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: start,
                    direction: end - start,
                }),
                domain: Interval::UNIT,
                same_sense: true,
            })
        }
        Piece::Arc {
            centre,
            x,
            y,
            radius,
            start,
            end,
            sweep,
        } => {
            // Measured: the left normal of a counter-clockwise arc points at
            // its centre, so the left offset SHRINKS the radius; clockwise
            // grows it. Centre, frame and domain are unchanged, which is what
            // makes the offset endpoints correspond at equal parameters.
            let offset_radius = radius - distance * sweep.signum();
            if offset_radius <= 0.0 {
                return Err(GeomError::Degenerate(format!(
                    "centre-line half-width {} collapses an arc of radius {radius}",
                    distance.abs()
                )));
            }
            Ok(ProfileSegment {
                curve: Curve2::Circle(Circle2 {
                    frame: axiolid_core::Frame2 {
                        origin: *centre,
                        x: *x,
                        y: *y,
                    },
                    radius: offset_radius,
                }),
                domain: Interval::new(*start, *end),
                same_sense: true,
            })
        }
    }
}
/// Resolve a centre-line profile into an exact closed contour.
///
/// The boundary is the left offset walked forward, a butt end cap, the right
/// offset walked back, and a butt start cap. Butt caps are used because the
/// source states a width and an extent, not an end treatment: a round or
/// square cap would add material the author never declared.
pub fn center_line_contour(
    profile: &CenterLineProfile,
    tolerance: Tolerance,
) -> GeomResult<ContourProfile> {
    // Explicit non-positive test so a NaN half-width is refused too.
    if profile.half_width <= 0.0 || profile.half_width.is_nan() {
        return Err(GeomError::Degenerate(format!(
            "centre line half-width must be positive, got {}",
            profile.half_width
        )));
    }
    if profile.path.segments.is_empty() {
        return Err(GeomError::Degenerate(
            "centre line path has no segments".to_owned(),
        ));
    }

    let pieces = profile
        .path
        .segments
        .iter()
        .map(read_piece)
        .collect::<GeomResult<Vec<_>>>()?;

    // The path must be connected and must NOT be closed: a closed path
    // denotes an annulus, which needs a hole rather than a single ring, and
    // silently treating it as open would weld the two ends together.
    let eps = tolerance.linear();
    for pair in pieces.windows(2) {
        let gap = (pair[1].start_point() - pair[0].end_point()).length();
        if gap > eps {
            return Err(GeomError::InvalidInput(format!(
                "centre line path is disconnected by {gap}"
            )));
        }
    }
    let first = pieces.first().expect("at least one segment");
    let last = pieces.last().expect("at least one segment");
    if (last.end_point() - first.start_point()).length() <= eps {
        return Err(unsupported(
            "closed centre-line path denotes an annulus, not a single ring",
        ));
    }

    // A tangent reversal makes the two offsets cross, producing a ring whose
    // enclosed area depends on where it self-intersects. Refused rather than
    // emitting a spike.
    for pair in pieces.windows(2) {
        let incoming = pair[0].end_tangent();
        let outgoing = pair[1].start_tangent();
        if incoming.dot(outgoing) <= -1.0 + 1e-12 {
            return Err(GeomError::Degenerate(
                "centre line reverses on itself".to_owned(),
            ));
        }
    }

    let half = profile.half_width;
    let mut segments = Vec::with_capacity(2 * pieces.len() + 2);

    // RIGHT side forward, then LEFT side back. Verified by shoelace: walking
    // the left side first traces the boundary clockwise, which builds the
    // solid inside-out -- the volume comes out correct in magnitude and
    // negative in sign, so only a SIGNED check catches it.
    for piece in &pieces {
        segments.push(offset_piece(piece, -half)?);
    }
    // End cap: straight across the width, right to left.
    let end_centre = last.end_point();
    let end_normal = left_normal(last.end_tangent());
    segments.push(straight(
        end_centre - end_normal * half,
        end_centre + end_normal * half,
    ));
    // Left side, backward: reverse the order AND each segment's sense.
    for piece in pieces.iter().rev() {
        let mut offset = offset_piece(piece, half)?;
        offset.same_sense = !offset.same_sense;
        segments.push(offset);
    }
    // Start cap closes the ring, left to right.
    let start_centre = first.start_point();
    let start_normal = left_normal(first.start_tangent());
    segments.push(straight(
        start_centre + start_normal * half,
        start_centre - start_normal * half,
    ));

    Ok(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

fn left_normal(tangent: Vec2) -> Vec2 {
    Vec2::new(-tangent.y, tangent.x)
}

fn straight(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}
