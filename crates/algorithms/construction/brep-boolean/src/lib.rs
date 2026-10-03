#![forbid(unsafe_code)]

//! General exact booleans between exact B-reps with analytic faces.
//!
//! The general-fuse pipeline of ADR 0075, built in stages:
//!
//! - [`section_edges`]: where two operands' faces cross, as exact curves
//!   trimmed exactly to both faces;
//! - [`split_face`]: one face cut along its section edges into regions, in
//!   the face's own parameters, with exact pcurves;
//! - [`boolean`]: the regions classified against the other solid, selected
//!   by operator and sewn into the result.
//!
//! Operands may touch as well as cross: faces sharing a patch of one
//! surface, sections running along existing edges, tangent contact, and
//! solids meeting along an edge or at a point are all handled, not
//! refused.
//!
//! Nothing is meshed or fitted. A configuration a step cannot build is
//! refused by name ([`BooleanError`]).
//!
//! # What a result guarantees
//!
//! Every surface of the result is an operand's surface, unchanged. Every
//! new edge lies on an exact section curve (a closed form or a certified
//! trace); vertices are those curves evaluated in `f64`.
//!
//! **Exact decisions.** Which section two supports have (closed forms,
//! ruled, torus and traced sections), where a section crosses the surface
//! next to a face's edge (`exact_curve_surface_intersection`), whether a
//! point lies in a face (certified membership; too close to decide is
//! refused, not guessed) and whether it lies in a solid (ray parity over
//! exact ray/surface intersections).
//!
//! **Decisions within tolerance.** Operands built or placed separately meet
//! faces that agree only up to rounding (#228). These decisions use the
//! caller's [`Tolerance`] and nothing else: `eps = tolerance.linear()` for
//! distances, `alpha = tolerance.angular()` for directions.
//!
//! - two supports are one surface when their normals or axes agree within
//!   `alpha` and their offsets within `eps`;
//! - two supports touch rather than cross where their normals agree within
//!   `alpha`;
//! - a point lies on an edge, or on a pole, within `eps`;
//! - a section runs along an edge when every point of the edge is proved
//!   within `eps` of it (two lines, or two circles, in closed form);
//! - a section that touches the surface next to an edge, whose double root
//!   rounding split or lost, is cut where it meets the edge itself, within
//!   `eps`;
//! - cuts on a section within `eps` of each other are one cut, a boundary
//!   split within `eps` of the boundary piece's end is that end, and two
//!   evaluations of a vertex within `eps` are one vertex;
//! - a section is an iso-parameter curve of a face within `eps`;
//! - a plane is parallel or perpendicular to a cylinder's axis, or touches
//!   the cylinder, when that moves the plane by at most `eps` over where
//!   the section can matter.
//!
//! **Exact first (#236).** A reading about the operands' own surfaces --
//! two supports are one, a plane is parallel or perpendicular to a
//! cylinder's axis or touches it -- is first put to an exact predicate on
//! the operands' `f64` numbers (dyadic arithmetic, `axiolid-exact`). Faces
//! that are exactly coplanar, parallel or perpendicular, as operands placed
//! by matrices with entries `0` and `+-1` meet, are decided exactly and use
//! no tolerance. A reading about constructed points (a cut, a vertex, a
//! point on an edge) compares `f64` evaluations of exact curves, so a
//! residue up to `2^-40` of the operands' extent (about four thousand units
//! in the last place) is the rounding of one exact point, read without the
//! tolerance; only a larger residue is a decision within tolerance.
//!
//! **Guarantee.** [`boolean_with_report`] returns a [`BooleanReport`] of the
//! within-tolerance decisions that fired, by kind, with the furthest each
//! moved or turned the operands. When it is empty
//! ([`BooleanReport::is_exact`]), the result is the exact boolean of the
//! operands as given: every decision was exact, and new vertices are exact
//! points evaluated in `f64`. At [`Tolerance::ZERO`] nothing can be read
//! within tolerance, so every result is exact; operands that only miss each
//! other by rounding are then refused, not guessed. When the report is not
//! empty, the result is the exact boolean of operands whose faces were
//! moved by at most `eps` (and, for a coincidence or contact of
//! directions, turned by at most `alpha`), with every surface and curve
//! exact for those; the report bounds how far. A reading that no single
//! such perturbation explains is refused instead: cuts chained within `eps`
//! of each other over more than `eps` are [`BooleanError::NearCoincidence`].
//! Features further apart than the tolerance always go through the exact
//! predicates: a gap of ten tolerances stays a gap.
//!
//! Bookkeeping in a face's parameters (ordering pieces leaving a vertex,
//! naming a point's parameters on the surface it was evaluated from) uses
//! fixed relative slacks far below any tolerance. It never decides
//! geometry: an order it cannot settle is refused
//! ([`BooleanError::TangentSplit`], [`BooleanError::UnclosedSplit`]).

mod assemble;
mod bounds;
mod classify;
mod predicate;
mod report;
mod seams;
mod section;
mod split;
mod support;

pub use report::{BooleanReport, ToleranceDecision, ToleranceDecisionKind};
pub use section::{section_edges, SectionEdge};
pub use split::{split_face, Piece, PieceSource, Region};

use axiolid_brep::ExactBRep;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_evaluate::surface::{locate, normal};
use axiolid_surface::Surface;
use axiolid_topology::Orientation;

/// The exact boolean of two exact B-rep solids.
///
/// ADR 0075 stages 1 and 2: faces on planes, cylinders, elliptical
/// cylinders, cones, spheres and tori, meeting in any section
/// `exact_surface_intersection` builds (lines, conics, ruled, torus and
/// traced sections).
///
/// # Errors
///
/// Anything a step refuses (see [`BooleanError`]), or a result with nothing
/// left.
pub fn boolean(
    a: &ExactBRep,
    b: &ExactBRep,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> Result<ExactBRep, BooleanError> {
    boolean_with_report(a, b, operator, tolerance).map(|(result, _)| result)
}

/// [`boolean`], with the within-tolerance decisions it took (#236).
///
/// An exact report ([`BooleanReport::is_exact`]) means no decision used the
/// tolerance: the result is the exact boolean of `a` and `b` as given (see
/// "What a result guarantees"). Otherwise the report names each kind of
/// reading that fired and the furthest it moved or turned the operands,
/// within `tolerance`. At [`Tolerance::ZERO`] nothing can be read within
/// tolerance, so a result is always exact: operands placed with exact axis
/// matrices (entries `0` and `+-1`) meet in exactly coincident, parallel and
/// perpendicular faces, which the exact predicates decide.
///
/// # Errors
///
/// As [`boolean`].
pub fn boolean_with_report(
    a: &ExactBRep,
    b: &ExactBRep,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> Result<(ExactBRep, BooleanReport), BooleanError> {
    let session = report::open(&[a, b]);
    let result = run(a, b, operator, tolerance)?;
    Ok((result, session.finish()))
}

fn run(
    a: &ExactBRep,
    b: &ExactBRep,
    operator: BooleanOperator,
    tolerance: Tolerance,
) -> Result<ExactBRep, BooleanError> {
    // Faces that wind round their surface without a seam get one first.
    let seamed_a = seams::with_seams(a, tolerance)?;
    let seamed_b = seams::with_seams(b, tolerance)?;
    let a = seamed_a.as_ref().unwrap_or(a);
    let b = seamed_b.as_ref().unwrap_or(b);
    let edges = section_edges(a, b, tolerance)?;
    let solid_a = classify::Solid::new(a, tolerance)?;
    let solid_b = classify::Solid::new(b, tolerance)?;
    let mut kept = Vec::new();
    let cuts: Vec<axiolid_core::Point3> = edges.iter().flat_map(|e| [e.start, e.end]).collect();
    for (operand, other, other_solid, first) in [(a, b, &solid_b, true), (b, a, &solid_a, false)] {
        let topology = operand.topology();
        for index in 0..topology.faces().len() {
            let face = topology
                .face_id_at(index)
                .ok_or(BooleanError::DanglingReference)?;
            let record = &topology.faces()[index];
            let surface: Surface = record
                .surface
                .and_then(|id| operand.surfaces().get(id.index()))
                .ok_or(BooleanError::DanglingReference)?
                .clone();
            let mine: Vec<SectionEdge> = edges
                .iter()
                .filter(|e| {
                    if first {
                        e.face_a == face
                    } else {
                        e.face_b == face
                    }
                })
                .cloned()
                .collect();
            // The other operand's faces on this face's surface: a region
            // may lie on one of them rather than inside or outside.
            let coincident: Vec<usize> = (0..other.topology().faces().len())
                .filter(|&f| {
                    other.topology().faces()[f]
                        .surface
                        .and_then(|s| other.surfaces().get(s.index()))
                        .is_some_and(|s| support::same_support(&surface, s, tolerance))
                })
                .collect();
            for region in split_face(operand, face, &mine, first, &cuts, tolerance)? {
                let sign = match record.orientation {
                    Orientation::Forward if !region.against => 1.0,
                    Orientation::Reversed if region.against => 1.0,
                    _ => -1.0,
                };
                let (keep, flip) = decide(
                    &region,
                    &surface,
                    sign,
                    other_solid,
                    &coincident,
                    operator,
                    first,
                    tolerance,
                )?;
                if keep {
                    kept.push(assemble::Kept {
                        surface: surface.clone(),
                        orientation: record.orientation,
                        region,
                        flip,
                    });
                }
            }
        }
    }
    assemble::assemble(&kept, tolerance)
}

/// Whether a region is kept, and whether it bounds the result from its
/// other side: classified at the first interior point that decides.
#[allow(clippy::too_many_arguments)]
fn decide(
    region: &Region,
    surface: &Surface,
    sign: f64,
    other_solid: &classify::Solid<'_>,
    coincident: &[usize],
    operator: BooleanOperator,
    first: bool,
    tolerance: Tolerance,
) -> Result<(bool, bool), BooleanError> {
    let mut last = BooleanError::Undecided;
    for point in classify::interior_points(region, surface)? {
        match classify_point(
            point,
            surface,
            sign,
            other_solid,
            coincident,
            operator,
            first,
            tolerance,
        ) {
            Err(BooleanError::Undecided) => last = BooleanError::Undecided,
            other => return other,
        }
    }
    Err(last)
}

#[allow(clippy::too_many_arguments)]
fn classify_point(
    point: axiolid_core::Point3,
    surface: &Surface,
    sign: f64,
    other_solid: &classify::Solid<'_>,
    coincident: &[usize],
    operator: BooleanOperator,
    first: bool,
    tolerance: Tolerance,
) -> Result<(bool, bool), BooleanError> {
    if let Some(theirs) = other_solid.on_face(point, coincident, tolerance)? {
        // On the other solid's boundary: kept once, from the first
        // operand, where the operator leaves a boundary.
        let (u, v) = locate(surface, point, report::floored(tolerance))
            .map_err(|_| BooleanError::Evaluation)?;
        let ours = normal(surface, u, v).map_err(|_| BooleanError::Evaluation)? * sign;
        let same = ours.dot(theirs) > 0.0;
        let keep = first
            && match operator {
                BooleanOperator::Union | BooleanOperator::Intersection => same,
                BooleanOperator::Difference => !same,
                _ => return Err(BooleanError::UnsupportedSection),
            };
        return Ok((keep, false));
    }
    let inside = other_solid.contains(point, tolerance)?;
    Ok(match (operator, first) {
        (BooleanOperator::Union, _) => (!inside, false),
        (BooleanOperator::Intersection, _) => (inside, false),
        (BooleanOperator::Difference, true) => (!inside, false),
        (BooleanOperator::Difference, false) => (inside, true),
        _ => return Err(BooleanError::UnsupportedSection),
    })
}

use axiolid_measure::ExactMeasureError;
use core::fmt;

/// Why a boolean step could not be carried out exactly.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BooleanError {
    /// Two faces' supports intersect in a curve this stage does not take
    /// (anything but a line, circle or ellipse).
    UnsupportedSection,
    /// A face is trimmed by an edge or pcurve family this stage cannot
    /// intersect or classify against.
    UnsupportedTrim,
    /// A point was too close to a face boundary to classify.
    Undecided,
    /// A curve or surface could not be evaluated or inverted.
    Evaluation,
    /// A handle referenced missing geometry.
    DanglingReference,
    /// A face domain could not be built.
    Measure(ExactMeasureError),
    /// A face whose surface, or a section on it, has no exact pcurve in
    /// this stage.
    UnsupportedSplit,
    /// The pieces of a split face do not close into loops.
    UnclosedSplit,
    /// Two pieces leave a vertex of a split face in the same direction and
    /// bend alike there, so no order between them can be read off.
    TangentSplit,
    /// The operation leaves nothing: an intersection of solids that only
    /// touch, or a difference that removes everything. An exact B-rep
    /// cannot be empty, so this is the answer, not a refusal.
    EmptyResult,
    /// A cavity of the result lies inside none of its solids.
    AmbiguousCavity,
    /// The kept faces did not sew into a valid exact B-rep.
    Assembly,
    /// Features lie within tolerance of each other in a chain longer than
    /// the tolerance, so no single perturbation within tolerance makes them
    /// one: reading them as one would move a feature further than the
    /// caller allowed (#228).
    NearCoincidence,
}

impl fmt::Display for BooleanError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSection => {
                f.write_str("a face pair meets in a curve this stage does not build")
            }
            Self::UnsupportedTrim => {
                f.write_str("a face is trimmed by a curve family this stage cannot intersect")
            }
            Self::Undecided => f.write_str("a point lies too close to a face boundary to classify"),
            Self::Evaluation => f.write_str("a curve or surface could not be evaluated"),
            Self::DanglingReference => f.write_str("a handle references missing geometry"),
            Self::Measure(error) => write!(f, "a face domain could not be built: {error}"),
            Self::UnsupportedSplit => {
                f.write_str("a face or section has no exact pcurve in this stage")
            }
            Self::UnclosedSplit => f.write_str("the pieces of a split face do not close"),
            Self::TangentSplit => {
                f.write_str("two pieces leave a vertex of a split face in one direction")
            }
            Self::EmptyResult => f.write_str("the operation leaves nothing"),
            Self::AmbiguousCavity => {
                f.write_str("a cavity of the result lies inside none of its solids")
            }
            Self::Assembly => f.write_str("the kept faces did not sew into a valid exact B-rep"),
            Self::NearCoincidence => f.write_str(
                "features lie within tolerance of each other in a chain longer than the tolerance",
            ),
        }
    }
}

impl std::error::Error for BooleanError {}
