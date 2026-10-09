//! Seams of a station's basis curve: where they lie, which piece a station
//! on one reads, and the mitre between the two (#263, ADR 0082 amendment).
//!
//! # Where the seams are
//!
//! [`station_seams2`] and [`station_seams3`] list the places strictly
//! inside a curve where two of its pieces meet, each as a
//! [`StationSeam`]: its distance in the curve's station measure (the
//! [parent module](super)'s convention), its native parameter, whether
//! the stored data guarantees the same section frame on both sides
//! ([`StationSeam::smooth`]), and whether the distance is read from stored
//! data alone ([`StationSeam::exact`]). They are read without evaluating
//! the curve:
//!
//! - a polyline's interior vertices, at the running sum of its segment
//!   lengths (a repeated vertex is one seam), not smooth;
//! - a B-spline's interior knots of multiplicity at least its degree
//!   (where it may turn a corner), not smooth, and NOT exact: the arc
//!   length to a knot is a quadrature ([`crate::arc_parameter`]);
//! - an intrinsic curve's curvature (and, in 3D, torsion) law seams, and
//!   an arc-length chain's joins and its intrinsic pieces' seams, at their
//!   stored arc lengths, smooth: the heading is continuous by construction;
//! - an elevated curve's plan seams (smooth, as above) and its profile's
//!   seams, where the grade may jump (not smooth); a banked curve's, and
//!   the seams of its cant law (the roll may jump) and its pivot law (the
//!   grade of its point path may jump), all at their stored plan
//!   distances.
//!
//! A line, a circle and an ellipse have none. [`exact_station_seams2`] and
//! [`exact_station_seams3`] refuse, by a typed
//! [`GeomError::UnsupportedInput`], a curve with a seam that is not exact.
//!
//! # Which piece a station on a seam reads
//!
//! A station whose distance lies within [`ARC_LENGTH_TOLERANCE`]` * max(1,
//! s)` of a seam that is not smooth is ON it, and is read at the seam's
//! own distance from the piece [`SeamSide`] names: the piece that starts
//! there ([`SeamSide::Outgoing`], what every evaluator reads by default) or
//! the one that ends there ([`SeamSide::Incoming`]). The incoming piece is
//! read as the curve TRUNCATED at the seam, at its end: a polyline's
//! segment ending at the vertex, a B-spline's span ending at the knot
//! (the reversed spline's outgoing span), an elevated curve whose profile,
//! and a banked curve whose cant and pivot laws, stop at the seam.
//!
//! # Mitre
//!
//! Where a run of sections or offsets crosses a seam whose two tangents
//! differ, its section there stands in the [`Mitre`] plane: through the
//! seam's point, normal to the bisector `n = (t_in + t_out) / |t_in +
//! t_out|`. A section point is placed twice, once in each side's frame,
//! each projected along its own side's tangent onto the plane, which is
//! where that side's extrusion of the section meets it; the two agree
//! whenever the outgoing frame is the incoming one turned about `t_in x
//! t_out` (a plan turn, a grade break on a straight plan), and their
//! midpoint is taken otherwise. A longitudinal offset moves the section
//! along `n`. Tangents within [`SEAM_TANGENT_TOLERANCE`] of each other
//! need no mitre, and a near reversal, `|t_in + t_out| / 2` (the cosine of
//! half the turn) at most [`MITRE_TOLERANCE`], has no usable plane and is
//! refused by name.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Point2, Point3, Scalar, Vec3};
use axiolid_curve::{
    BSplineCurve, Banked3, CantLaw, ChainPiece2, Curve2, Curve3, Elevated3, ElevationLaw, SeamSide,
};

use super::{
    banked_frame, invalid, planar_frame, reference_up_frame, station_length2, station_length3,
    SectionFrame,
};
use crate::arc_parameter::{arc_length2, arc_length3, ARC_LENGTH_TOLERANCE};
use crate::curve::{derivative2, derivative3, domain2, domain3, evaluate2, evaluate3};

/// Largest angle, in radians, between the incoming and outgoing tangents
/// of a seam that still counts as one direction: no mitre is cut there.
pub const SEAM_TANGENT_TOLERANCE: Scalar = 1e-9;

/// Smallest cosine of half the turn at a seam, `|t_in + t_out| / 2`, that
/// a [`Mitre`] admits; below it the tangents nearly reverse and the mitre
/// plane would stretch a section without bound.
pub const MITRE_TOLERANCE: Scalar = 1e-6;

/// A place inside a curve where two of its pieces meet; see the
/// [module documentation](self).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StationSeam {
    /// Distance from the curve's start in its station measure: plan
    /// distance on an elevated or banked curve, arc length otherwise.
    pub distance: Scalar,
    /// The curve's own parameter there: the vertex index of a polyline,
    /// the knot of a B-spline, the distance itself on every family whose
    /// parameter is its measure.
    pub parameter: Scalar,
    /// Whether the stored data guarantees one section frame on both sides
    /// (a chain join, an intrinsic law's seam). `false` means the frame
    /// MAY jump: evaluate both sides to see whether it does.
    pub smooth: bool,
    /// Whether `distance` is read from stored data alone (sums of stored
    /// lengths and of segment lengths); `false` for a quadrature.
    pub exact: bool,
}

impl StationSeam {
    /// A seam at `distance` (the station measure) and `parameter` (the
    /// curve's own), `smooth` when the stored data guarantees one frame on
    /// both sides and `exact` when `distance` needs no quadrature.
    #[must_use]
    pub const fn new(distance: Scalar, parameter: Scalar, smooth: bool, exact: bool) -> Self {
        Self {
            distance,
            parameter,
            smooth,
            exact,
        }
    }

    /// A seam whose parameter is its distance.
    const fn measured(distance: Scalar, smooth: bool) -> Self {
        Self::new(distance, distance, smooth, true)
    }
}

/// The arc-length tolerance around `distance`: what counts as ON a seam.
fn slack(distance: Scalar) -> Scalar {
    ARC_LENGTH_TOLERANCE * distance.abs().max(1.0)
}

/// Seams strictly inside `(0, length)`, sorted, those within tolerance of
/// each other merged (smooth and exact only if every merged one is).
pub(super) fn tidy(mut seams: Vec<StationSeam>, length: Option<Scalar>) -> Vec<StationSeam> {
    seams.retain(|seam| {
        seam.distance.is_finite()
            && seam.distance > slack(0.0)
            && length.is_none_or(|l| seam.distance < l - slack(l))
    });
    seams.sort_by(|a, b| a.distance.total_cmp(&b.distance));
    let mut out: Vec<StationSeam> = Vec::with_capacity(seams.len());
    for seam in seams {
        match out.last_mut() {
            Some(last) if seam.distance - last.distance <= slack(seam.distance) => {
                last.smooth &= seam.smooth;
                last.exact &= seam.exact;
            }
            _ => out.push(seam),
        }
    }
    out
}

/// The non-degenerate segments of a polyline's points: `(start vertex,
/// end vertex, length)`, a closed polyline's wrap included.
fn polyline_segments<P>(
    points: &[P],
    closed: bool,
    length: impl Fn(P, P) -> Scalar,
) -> Vec<(usize, usize, Scalar)>
where
    P: Copy,
{
    let count = points.len();
    let spans: Vec<(usize, usize)> = if closed {
        (0..count).map(|i| (i, (i + 1) % count)).collect()
    } else {
        (0..count.saturating_sub(1)).map(|i| (i, i + 1)).collect()
    };
    spans
        .into_iter()
        .map(|(i, j)| (i, j, length(points[i], points[j])))
        .collect()
}

/// A polyline's seams: the vertex where each non-degenerate segment
/// starts after an earlier non-degenerate one, at the running sum of the
/// segment lengths. The parameter is that segment's start vertex.
fn polyline_seams(segments: &[(usize, usize, Scalar)]) -> Vec<StationSeam> {
    let mut out = Vec::new();
    let mut run = 0.0;
    let mut started = false;
    for &(start, _, length) in segments {
        if length > 0.0 {
            if started {
                out.push(StationSeam::new(run, start as Scalar, false, true));
            }
            started = true;
        }
        run += length;
    }
    out
}

/// The segment a polyline vertex seam reads on `side`: the first
/// non-degenerate one starting at `vertex` or later, or the last ending at
/// it or earlier. Its index is its start vertex.
fn polyline_side(
    segments: &[(usize, usize, Scalar)],
    vertex: usize,
    side: SeamSide,
) -> GeomResult<(usize, usize)> {
    let at = segments
        .iter()
        .position(|&(start, _, _)| start == vertex)
        .ok_or_else(|| invalid("station: the seam is not a vertex of the polyline".into()))?;
    let found = match side {
        SeamSide::Incoming => segments[..at].iter().rev().find(|s| s.2 > 0.0),
        _ => segments[at..].iter().find(|s| s.2 > 0.0),
    };
    found.map(|&(i, j, _)| (i, j)).ok_or_else(|| {
        invalid("station: the polyline has no segment on that side of the seam".into())
    })
}

/// Interior knots of a B-spline where it may turn a corner (multiplicity
/// at least its degree), inside its domain.
fn corner_knots(breaks: Vec<Scalar>, lo: Scalar, hi: Scalar) -> Vec<Scalar> {
    breaks.into_iter().filter(|&t| t > lo && t < hi).collect()
}

/// The same B-spline traversed backwards: knots mirrored about the sum of
/// the first and last, every list reversed. Its outgoing span at the
/// mirrored knot `first + last - t` is this spline's incoming span at `t`.
fn reversed_spline<P: Clone>(spline: &BSplineCurve<P>) -> (BSplineCurve<P>, Scalar) {
    let first = spline.knots.first().copied().unwrap_or(0.0);
    let last = spline.knots.last().copied().unwrap_or(0.0);
    let pivot = first + last;
    let mut reversed = spline.clone();
    reversed.knots = spline.knots.iter().rev().map(|t| pivot - t).collect();
    reversed.multiplicities.reverse();
    reversed.control_points.reverse();
    if let Some(weights) = reversed.weights.as_mut() {
        weights.reverse();
    }
    (reversed, pivot)
}

/// The seams of a 2D curve, every one located (see the
/// [module documentation](self)); a B-spline's corner knots by quadrature,
/// flagged not exact.
///
/// # Errors
///
/// A family a station cannot be measured on, a malformed length, and the
/// quadrature's refusals, by name.
pub fn station_seams2(curve: &Curve2) -> GeomResult<Vec<StationSeam>> {
    let length = station_length2(curve)?;
    let seams = match curve {
        Curve2::Polyline(polyline) => polyline_seams(&polyline_segments(
            &polyline.points,
            polyline.closed,
            |a: Point2, b: Point2| (b - a).length(),
        )),
        Curve2::BSpline(_) => {
            let domain = domain2(curve);
            let (lo, hi) = (domain.start.min(domain.end), domain.start.max(domain.end));
            corner_knots(crate::bound::continuity_breaks2(curve, 1), lo, hi)
                .into_iter()
                .map(|t| {
                    Ok(StationSeam::new(
                        arc_length2(curve, lo, t)?,
                        t,
                        false,
                        false,
                    ))
                })
                .collect::<GeomResult<_>>()?
        }
        Curve2::Intrinsic(intrinsic) => intrinsic
            .curvature
            .seams_within(intrinsic.length)
            .into_iter()
            .map(|s| StationSeam::measured(s, true))
            .collect(),
        Curve2::Chain(chain) => chain
            .joins()
            .unwrap_or_default()
            .into_iter()
            .chain(crate::chain::chain_curvature_seams(chain))
            .map(|s| StationSeam::measured(s, true))
            .collect(),
        _ => Vec::new(),
    };
    Ok(tidy(seams, length))
}

/// The seams of a 3D curve, every one located; see [`station_seams2`].
///
/// # Errors
///
/// As [`station_seams2`].
pub fn station_seams3(curve: &Curve3) -> GeomResult<Vec<StationSeam>> {
    let length = station_length3(curve)?;
    let seams = match curve {
        Curve3::Polyline(polyline) => polyline_seams(&polyline_segments(
            &polyline.points,
            polyline.closed,
            |a: Point3, b: Point3| (b - a).length(),
        )),
        Curve3::BSpline(_) => {
            let domain = domain3(curve);
            let (lo, hi) = (domain.start.min(domain.end), domain.start.max(domain.end));
            corner_knots(crate::bound::continuity_breaks3(curve, 1), lo, hi)
                .into_iter()
                .map(|t| {
                    Ok(StationSeam::new(
                        arc_length3(curve, lo, t)?,
                        t,
                        false,
                        false,
                    ))
                })
                .collect::<GeomResult<_>>()?
        }
        Curve3::Intrinsic(intrinsic) => intrinsic
            .curvature
            .seams_within(intrinsic.length)
            .into_iter()
            .chain(intrinsic.torsion.seams_within(intrinsic.length))
            .map(|s| StationSeam::measured(s, true))
            .collect(),
        Curve3::Elevated(elevated) => elevated_seams(elevated)?,
        Curve3::Banked(banked) => {
            let mut seams = elevated_seams(&banked.base)?;
            seams.extend(
                banked
                    .cant
                    .seams()
                    .into_iter()
                    .chain(banked.pivot.seams())
                    .map(|s| StationSeam::measured(s, false)),
            );
            seams
        }
        _ => Vec::new(),
    };
    Ok(tidy(seams, length))
}

/// An elevated curve's plan seams and its profile's.
fn elevated_seams(curve: &Elevated3) -> GeomResult<Vec<StationSeam>> {
    let mut seams = match curve.plan.as_ref() {
        // The plan's own seams are in its arc length, which is the plan
        // distance; a line or a circle has none.
        plan @ (Curve2::Intrinsic(_) | Curve2::Chain(_)) => station_seams2(plan)?,
        _ => Vec::new(),
    };
    seams.extend(
        crate::elevation::elevation_seams(&curve.elevation)
            .into_iter()
            .map(|s| StationSeam::measured(s, false)),
    );
    Ok(seams)
}

fn not_exact(seams: Vec<StationSeam>) -> GeomResult<Vec<StationSeam>> {
    if seams.iter().any(|seam| !seam.exact) {
        return Err(GeomError::UnsupportedInput {
            backend: BackendId::new("axiolid-evaluate"),
            operation: Operation::CurveEvaluation,
            input: "exact seam distances of a B-spline with a corner knot: its arc length to \
                    the knot is a quadrature",
        });
    }
    Ok(seams)
}

/// A chain whose parametric piece turns a corner inside itself: that
/// corner is a seam no stored length locates.
fn chain_hides_corners(curve: &Curve2) -> bool {
    match curve {
        Curve2::Chain(chain) => chain.pieces.iter().any(|piece| match piece {
            ChainPiece2::Parametric { curve, .. } => {
                !crate::bound::continuity_breaks2(curve, 1).is_empty()
            }
            _ => false,
        }),
        _ => false,
    }
}

fn hidden_corner() -> GeomError {
    GeomError::UnsupportedInput {
        backend: BackendId::new("axiolid-evaluate"),
        operation: Operation::CurveEvaluation,
        input: "exact seam distances of a chain whose parametric piece turns a corner inside \
                itself: no stored length locates it",
    }
}

/// [`station_seams2`], refusing a curve with a seam whose distance is not
/// read from stored data alone.
///
/// # Errors
///
/// As [`station_seams2`], and [`GeomError::UnsupportedInput`] for a
/// B-spline with a corner knot or a chain whose parametric piece has one.
pub fn exact_station_seams2(curve: &Curve2) -> GeomResult<Vec<StationSeam>> {
    if chain_hides_corners(curve) {
        return Err(hidden_corner());
    }
    not_exact(station_seams2(curve)?)
}

/// [`station_seams3`], refusing a curve with a seam whose distance is not
/// read from stored data alone.
///
/// # Errors
///
/// As [`exact_station_seams2`]; an elevated or banked curve's chain plan
/// is checked the same way.
pub fn exact_station_seams3(curve: &Curve3) -> GeomResult<Vec<StationSeam>> {
    let plan = match curve {
        Curve3::Elevated(elevated) => Some(elevated.plan.as_ref()),
        Curve3::Banked(banked) => Some(banked.base.plan.as_ref()),
        _ => None,
    };
    if plan.is_some_and(chain_hides_corners) {
        return Err(hidden_corner());
    }
    not_exact(station_seams3(curve)?)
}

// --- a station on a seam ------------------------------------------------------

/// The seam a station at the admitted `distance` lies on, if any: a seam
/// that is not smooth within the arc-length tolerance.
fn seam_near(seams: Vec<StationSeam>, distance: Scalar) -> Option<StationSeam> {
    seams
        .into_iter()
        .filter(|seam| !seam.smooth && (seam.distance - distance).abs() <= slack(distance))
        .min_by(|a, b| {
            (a.distance - distance)
                .abs()
                .total_cmp(&(b.distance - distance).abs())
        })
}

/// Whether a family can have a seam a station reads two ways; the others
/// skip locating seams altogether.
fn sided2(curve: &Curve2) -> bool {
    matches!(curve, Curve2::Polyline(_) | Curve2::BSpline(_))
}

pub(super) fn sided3(curve: &Curve3) -> bool {
    matches!(
        curve,
        Curve3::Polyline(_) | Curve3::BSpline(_) | Curve3::Elevated(_) | Curve3::Banked(_)
    )
}

/// The section of a 2D curve ON a seam near `distance`, read from `side`;
/// `None` when no seam is near.
pub(super) fn seam_section2(
    curve: &Curve2,
    distance: Scalar,
    side: SeamSide,
) -> GeomResult<Option<SectionFrame>> {
    if !sided2(curve) {
        return Ok(None);
    }
    let Some(seam) = seam_near(station_seams2(curve)?, distance) else {
        return Ok(None);
    };
    let frame = match curve {
        Curve2::Polyline(polyline) => {
            let segments = polyline_segments(&polyline.points, polyline.closed, |a: Point2, b| {
                (b - a).length()
            });
            let vertex = seam.parameter as usize;
            let (i, j) = polyline_side(&segments, vertex, side)?;
            planar_frame(
                polyline.points[vertex],
                polyline.points[j] - polyline.points[i],
            )?
        }
        Curve2::BSpline(spline) => {
            let point = evaluate2(curve, seam.parameter)?;
            let tangent = match side {
                SeamSide::Incoming => {
                    let (reversed, pivot) = reversed_spline(spline);
                    -derivative2(&Curve2::BSpline(reversed), pivot - seam.parameter)?
                }
                _ => derivative2(curve, seam.parameter)?,
            };
            planar_frame(point, tangent)?
        }
        _ => return Ok(None),
    };
    Ok(Some(frame))
}

/// The seam a station at the admitted `distance` along a 3D curve lies ON
/// and reads two ways, if any: the rule [`seam_section3`] applies, and
/// the curve-evaluation provider's sided queries with it (#286).
pub(crate) fn seam_at3(curve: &Curve3, distance: Scalar) -> GeomResult<Option<StationSeam>> {
    if !sided3(curve) {
        return Ok(None);
    }
    Ok(seam_near(station_seams3(curve)?, distance))
}

/// The section of a 3D curve ON a seam near `distance`, read from `side`;
/// `None` when no seam is near.
pub(super) fn seam_section3(
    curve: &Curve3,
    distance: Scalar,
    side: SeamSide,
) -> GeomResult<Option<SectionFrame>> {
    let Some(seam) = seam_at3(curve, distance)? else {
        return Ok(None);
    };
    let incoming = side == SeamSide::Incoming;
    let at = seam.distance;
    let frame = match curve {
        Curve3::Polyline(polyline) => {
            let segments = polyline_segments(&polyline.points, polyline.closed, |a: Point3, b| {
                (b - a).length()
            });
            let vertex = seam.parameter as usize;
            let (i, j) = polyline_side(&segments, vertex, side)?;
            reference_up_frame(
                polyline.points[vertex],
                polyline.points[j] - polyline.points[i],
            )?
        }
        Curve3::BSpline(spline) => {
            let point = evaluate3(curve, seam.parameter)?;
            let tangent = if incoming {
                let (reversed, pivot) = reversed_spline(spline);
                -derivative3(&Curve3::BSpline(reversed), pivot - seam.parameter)?
            } else {
                derivative3(curve, seam.parameter)?
            };
            reference_up_frame(point, tangent)?
        }
        Curve3::Elevated(elevated) => {
            let truncated;
            let read = if incoming {
                truncated = Elevated3 {
                    plan: elevated.plan.clone(),
                    elevation: elevation_ending_at(&elevated.elevation, at),
                };
                &truncated
            } else {
                elevated
            };
            reference_up_frame(
                crate::arc_length::elevated_point(read, at)?,
                crate::arc_length::elevated_tangent(read, at)?,
            )?
        }
        Curve3::Banked(banked) => {
            let truncated;
            let read = if incoming {
                truncated = Banked3 {
                    base: Elevated3 {
                        plan: banked.base.plan.clone(),
                        elevation: elevation_ending_at(&banked.base.elevation, at),
                    },
                    cant: law_ending_at(&banked.cant, at),
                    pivot: law_ending_at(&banked.pivot, at),
                    ..banked.clone()
                };
                &truncated
            } else {
                banked
            };
            banked_frame(read, at)?
        }
        _ => return Ok(None),
    };
    Ok(Some(frame))
}

/// A profile cut at `distance`: the pieces starting before it (outside
/// the tolerance), nested laws cut the same way, so the piece read at
/// `distance` is the one ENDING there.
fn elevation_ending_at(law: &ElevationLaw, distance: Scalar) -> ElevationLaw {
    match law {
        ElevationLaw::Piecewise { breaks, laws } if laws.len() == breaks.len() + 1 => {
            let keep = breaks.partition_point(|b| *b < distance - slack(distance));
            let start = if keep == 0 { 0.0 } else { breaks[keep - 1] };
            let mut kept = laws[..keep].to_vec();
            kept.push(elevation_ending_at(&laws[keep], distance - start));
            ElevationLaw::Piecewise {
                breaks: breaks[..keep].to_vec(),
                laws: kept,
            }
        }
        other => other.clone(),
    }
}

/// A cant or pivot law cut at `distance`: the pieces starting before it,
/// so its last piece, closed at its far end, is the one ending there.
fn law_ending_at(law: &CantLaw, distance: Scalar) -> CantLaw {
    let mut start = 0.0;
    let mut pieces = Vec::with_capacity(law.pieces.len());
    for piece in &law.pieces {
        if !pieces.is_empty() && start >= distance - slack(distance) {
            break;
        }
        pieces.push(piece.clone());
        start += piece.length;
    }
    CantLaw::new(pieces)
}

// --- mitre --------------------------------------------------------------------

/// The plane a section stands in where a run crosses a seam whose
/// tangents differ; see the [module documentation](self#mitre).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mitre {
    /// The section frame of the piece ending at the seam.
    pub incoming: SectionFrame,
    /// The section frame of the piece starting there.
    pub outgoing: SectionFrame,
    /// Unit normal of the mitre plane, `(t_in + t_out) / |t_in + t_out|`.
    pub normal: Vec3,
}

impl Mitre {
    /// The mitre between the two frames of a seam, or `None` when their
    /// tangents agree within [`SEAM_TANGENT_TOLERANCE`].
    ///
    /// # Errors
    ///
    /// Tangents that nearly reverse (the cosine of half the turn at most
    /// [`MITRE_TOLERANCE`]), and two frames that do not share their point
    /// (a seam that is not position continuous), by name.
    pub fn between(incoming: &SectionFrame, outgoing: &SectionFrame) -> GeomResult<Option<Self>> {
        let (t_in, t_out) = (incoming.tangent, outgoing.tangent);
        let turn = t_in.cross(t_out).length().atan2(t_in.dot(t_out));
        if !turn.is_finite() {
            return Err(invalid("station: a seam's tangents are not finite".into()));
        }
        if turn <= SEAM_TANGENT_TOLERANCE {
            return Ok(None);
        }
        let sum = t_in + t_out;
        let half = 0.5 * sum.length();
        if half <= MITRE_TOLERANCE {
            return Err(GeomError::Degenerate(format!(
                "station: the curve turns back on itself at a seam (by {turn} rad), so a run \
                 across it has no mitre plane"
            )));
        }
        let gap = (outgoing.point - incoming.point).length();
        let scale = incoming.point.abs().max_element().max(1.0);
        if gap > ARC_LENGTH_TOLERANCE * scale {
            return Err(invalid(format!(
                "station: the curve's pieces are {gap} apart at a seam, so a run across it has \
                 no mitre"
            )));
        }
        Ok(Some(Self {
            incoming: *incoming,
            outgoing: *outgoing,
            normal: sum / (2.0 * half),
        }))
    }

    /// A section point in the mitre plane: `incoming` and `outgoing` are
    /// the same section point placed in each side's frame (without a
    /// longitudinal offset); each is projected along its side's tangent
    /// onto the plane, and the midpoint is moved `longitudinal` along the
    /// normal.
    #[must_use]
    pub fn place(&self, incoming: Point3, outgoing: Point3, longitudinal: Scalar) -> Point3 {
        let project = |frame: &SectionFrame, point: Point3| {
            let height = (point - frame.point).dot(self.normal);
            point - frame.tangent * (height / frame.tangent.dot(self.normal))
        };
        let a = project(&self.incoming, incoming);
        let b = project(&self.outgoing, outgoing);
        (a + b) * 0.5 + self.normal * longitudinal
    }
}
