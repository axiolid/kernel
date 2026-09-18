//! Exact distance-to-parameter conversion for polylines.
//!
//! A polyline is the one refused family whose arc length is a finite
//! sum rather than an integral: locating a distance is a running sum
//! until the interval is found, then one linear interpolation. The
//! only transcendental is the `sqrt` in each segment length, which
//! `Line` already relies on and reports as `ArcLength3d`, so refusing
//! the sequence while accepting each element was not defensible on
//! exactness grounds (kernel#107).
//!
//! Three edges have a plausible-looking wrong answer, so each is
//! pinned deliberately rather than left to fall out of the code:
//!
//! - A SEAM is two-valued. At an interior vertex the tangent jumps,
//!   so this returns the parameter of the OUTGOING segment: a
//!   distance that lands exactly on a vertex reads the heading the
//!   curve is about to take, not the one it arrived with.
//! - A ZERO-LENGTH segment (repeated identical points) has no
//!   direction. Skipping it would silently change the
//!   parameterisation, so it is refused instead.
//! - A CLOSED polyline's wrap segment is real length and is
//!   included, but distance still may not exceed the total: it is
//!   clamped by refusal, never wrapped around.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point3, Scalar};
use axiolid_curve::Polyline3;

/// Index pairs of the segments of `points`, honouring `closed`.
fn spans(count: usize, closed: bool) -> Vec<(usize, usize)> {
    if closed {
        (0..count).map(|i| (i, (i + 1) % count)).collect()
    } else {
        (0..count.saturating_sub(1)).map(|i| (i, i + 1)).collect()
    }
}

/// Total 3D length of `polyline`, or an error if any segment is
/// degenerate or non-finite.
pub fn polyline_length(polyline: &Polyline3) -> GeomResult<Scalar> {
    let mut total = 0.0;
    for (i, j) in spans(polyline.points.len(), polyline.closed) {
        total += segment_length(polyline.points[i], polyline.points[j])?;
    }
    Ok(total)
}

/// Length of one segment, refusing the degenerate cases by name.
fn segment_length(a: Point3, b: Point3) -> GeomResult<Scalar> {
    let length = (b - a).length();
    if !length.is_finite() {
        return Err(GeomError::Degenerate(
            "polyline segment has a non-finite length".to_string(),
        ));
    }
    if length <= 0.0 {
        return Err(GeomError::Degenerate(
            "polyline has a zero-length segment, so distance along it is ambiguous".to_string(),
        ));
    }
    Ok(length)
}

/// Convert a distance along `polyline` to its native parameter.
///
/// The parameter runs `[0, segment_count]` with the integer part
/// selecting the segment, matching `curve::evaluate3`.
pub fn polyline_parameter(polyline: &Polyline3, distance: Scalar) -> GeomResult<Scalar> {
    let spans = spans(polyline.points.len(), polyline.closed);
    if spans.is_empty() {
        return Err(GeomError::Degenerate(format!(
            "polyline with {} points has no evaluable segment",
            polyline.points.len()
        )));
    }

    // Every segment is measured up front so a degenerate one is
    // refused even when the requested distance stops short of it:
    // the curve is ill-defined as a whole, not just past that point.
    let mut lengths = Vec::with_capacity(spans.len());
    for (i, j) in &spans {
        lengths.push(segment_length(polyline.points[*i], polyline.points[*j])?);
    }
    let total: Scalar = lengths.iter().sum();

    if !distance.is_finite() || distance < 0.0 || distance > total {
        return Err(GeomError::InvalidInput(format!(
            "distance {distance} is outside the polyline's length [0, {total}]"
        )));
    }

    // Running sum until the interval is found. Both `<` and `<=` yield the
    // same parameter at a vertex (`t` lands exactly on the integer either
    // way); what actually pins the seam is `polyline_span`'s `floor`, which
    // sends an integral `t` to the segment STARTING there. `<=` is kept so
    // the loop resolves a vertex distance without relying on the trailing
    // fallback below.
    let mut run = 0.0;
    for (index, length) in lengths.iter().enumerate() {
        if distance <= run + length {
            let local = (distance - run) / length;
            return Ok(index as Scalar + local);
        }
        run += length;
    }

    // Only reachable when rounding leaves `distance` a hair above
    // the accumulated sum; the domain end is the exact answer.
    Ok(spans.len() as Scalar)
}
