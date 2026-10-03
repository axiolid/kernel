//! Which within-tolerance decisions a boolean took (#236).
//!
//! Every decision the crate docs list under "decisions within tolerance"
//! goes through one of the helpers here, which either settles it without
//! the tolerance or records it in the report of the boolean that is
//! running:
//!
//! - **Support decisions** ([`support`]) are about the operands' own
//!   surfaces: two supports are one, a plane is parallel or perpendicular to
//!   a cylinder's axis or touches it. An exact predicate on the operands'
//!   numbers (`axiolid-exact`, see `crate::predicate`) settles the exactly
//!   coincident, parallel and perpendicular cases; only a reading the exact
//!   predicate does not confirm is recorded.
//! - **Point decisions** ([`near`], [`near_angle`]) are about constructed
//!   points and curves: a cut on an edge, cuts or vertices read as one, a
//!   section running along an edge. Those values are rounded `f64`
//!   evaluations of exact curves, so a residue up to the rounding floor
//!   ([`rounding`]: `2^-40` of the operands' extent, about `4000` units in
//!   the last place) is the rounding of one exact point, not a decision; a
//!   residue above it and within the caller's tolerance is recorded.
//!
//! The report lives in a per-thread session opened by the public entry
//! points ([`crate::boolean_with_report`], [`crate::boolean`],
//! [`crate::section_edges`], [`crate::split_face`]) with the caller's
//! tolerance; an entry point called inside another's session joins it.
//! Nothing is shared across threads and nothing outlives the call.
//!
//! **Nothing is read at a zero tolerance (#251).** A support reading the
//! exact predicate rejects moves or turns the operands by a positive
//! amount, however small its `f64` measure: one whose measure rounds to
//! exactly `0` is not exact. So a tolerance with a zero part takes no such
//! reading ([`support`]): at [`Tolerance::ZERO`] the exact answer stands
//! and the general closed form decides. The session enforces the contract
//! too: a decision recorded at [`Tolerance::ZERO`], or beyond the caller's
//! tolerance, is a bug -- a debug assertion -- and in a release build the
//! entry point refuses it ([`BooleanError::ToleranceExceeded`]) rather than
//! return a result its report misdescribes.

use std::cell::RefCell;

use axiolid_brep::ExactBRep;
use axiolid_core::{Scalar, Tolerance};
use axiolid_surface::Surface;

use crate::BooleanError;

/// What a within-tolerance decision read (#236).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ToleranceDecisionKind {
    /// Two faces' supports read as one surface: normals or axes within the
    /// angular tolerance, offsets within the linear one.
    CoincidentSupports,
    /// Two surfaces read as touching, not crossing, where their normals
    /// agree within the angular tolerance.
    Contact,
    /// A point read as lying on an edge, a curve, a pole or a surface.
    IncidentPoint,
    /// A section read as running along a boundary edge.
    SectionAlongEdge,
    /// A section touching the surface next to an edge, cut where it meets
    /// the edge itself.
    TangentCrossing,
    /// Cuts on a curve, a boundary split and the piece's end, or two
    /// evaluations of a vertex, read as one point.
    MergedPoints,
    /// A section read as an iso-parameter curve of a face (a meridian, a
    /// latitude, a ruling or a circle about an axis).
    IsoCurve,
    /// A plane read as parallel to a cylinder's axis.
    PlaneParallelToAxis,
    /// A plane read as perpendicular to a cylinder's axis.
    PlanePerpendicularToAxis,
    /// A plane read as touching a cylinder along a ruling.
    PlaneTouchesCylinder,
}

/// The worst reading of one kind a boolean made within tolerance.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToleranceDecision {
    /// What was read.
    pub kind: ToleranceDecisionKind,
    /// The largest distance a feature was moved by a reading of this kind,
    /// in the operands' length unit (zero for a purely angular reading).
    pub linear: Scalar,
    /// The largest angle a direction was turned by a reading of this kind,
    /// in radians (zero for a purely linear reading).
    pub angular: Scalar,
}

/// The within-tolerance decisions one boolean took (#236).
///
/// Empty means no such decision fired: the result is the exact boolean of
/// the operands as given, its new vertices and curves evaluated in `f64`
/// (see the crate docs). Otherwise it is the exact boolean of operands
/// moved by at most [`Self::linear`] and turned by at most
/// [`Self::angular`], both within the caller's tolerance.
///
/// Either way, constructed points closer than [`Self::rounding_floor`]
/// were read as one point without a decision, and are not reported (#244).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BooleanReport {
    decisions: Vec<ToleranceDecision>,
    /// The operands' extent the rounding floor was scaled from (#244).
    extent: Scalar,
}

impl BooleanReport {
    /// One entry per kind of reading that fired, in [`ToleranceDecisionKind`]
    /// order, each with its worst magnitude.
    #[must_use]
    pub fn decisions(&self) -> &[ToleranceDecision] {
        &self.decisions
    }

    /// Whether no within-tolerance decision fired: the result is the exact
    /// boolean of the operands as given.
    #[must_use]
    pub fn is_exact(&self) -> bool {
        self.decisions.is_empty()
    }

    /// The largest distance any reading moved a feature (zero when exact).
    #[must_use]
    pub fn linear(&self) -> Scalar {
        self.decisions
            .iter()
            .map(|d| d.linear)
            .fold(0.0, Scalar::max)
    }

    /// The largest angle any reading turned a direction (zero when exact).
    #[must_use]
    pub fn angular(&self) -> Scalar {
        self.decisions
            .iter()
            .map(|d| d.angular)
            .fold(0.0, Scalar::max)
    }

    /// Whether a reading of `kind` fired.
    #[must_use]
    pub fn contains(&self, kind: ToleranceDecisionKind) -> bool {
        self.decisions.iter().any(|d| d.kind == kind)
    }

    /// The absolute rounding floor that applied to this result (#244), in
    /// the operands' length unit: [`ROUNDING_FACTOR`] times
    /// [`Self::extent`].
    ///
    /// A residue between constructed points (cuts, vertices, a point on an
    /// edge) up to this floor was read as one exact point evaluated twice,
    /// without the tolerance, and is in no [`ToleranceDecision`]; only a
    /// larger one is a decision within tolerance. So a consumer widening
    /// distances on the result widens them by this floor too, even when
    /// [`Self::is_exact`] holds. Zero when no general boolean produced the
    /// report (the default report).
    #[must_use]
    pub fn rounding_floor(&self) -> Scalar {
        ROUNDING_FACTOR * self.extent
    }

    /// The extent [`Self::rounding_floor`] was scaled from (#244): the
    /// largest coordinate magnitude of the operands as given to the
    /// boolean, the larger of the two operands' readings, each the largest
    /// of
    ///
    /// - `max(|x|, |y|, |z|)` over every topology vertex, and
    /// - per plane, cylinder, elliptical cylinder, cone, sphere or torus
    ///   surface, `max(|x|, |y|, |z|)` of its frame origin plus its size:
    ///   zero for a plane, `|radius|` for a cylinder, a cone (its radius at
    ///   the frame origin) or a sphere, the larger semi-axis for an
    ///   elliptical cylinder, `|major + minor|` for a torus; other surface
    ///   families add nothing beyond their vertices.
    ///
    /// Coordinates are read about the operands' own origin, so the extent
    /// grows with their distance from it. An operand whose reading is not
    /// finite counts as zero. After [`Self::merged`], the largest extent of the merged
    /// reports.
    #[must_use]
    pub fn extent(&self) -> Scalar {
        self.extent
    }

    /// This report and `other` together: what a chain of booleans read
    /// (each kind once, with the worst magnitude of either).
    ///
    /// The merged rounding floor is the larger of the two
    /// ([`Self::extent`] is the larger extent), sound for every boolean of
    /// the chain: each one's floor is at most it. A rigid placement of a
    /// result leaves distances between its points unchanged, so a report
    /// carried through one keeps its floor as is (#244).
    #[must_use]
    pub fn merged(mut self, other: &BooleanReport) -> BooleanReport {
        for d in &other.decisions {
            self.add(d.kind, d.linear, d.angular);
        }
        self.extent = self.extent.max(other.extent);
        self
    }

    fn add(&mut self, kind: ToleranceDecisionKind, linear: Scalar, angular: Scalar) {
        match self.decisions.binary_search_by(|d| d.kind.cmp(&kind)) {
            Ok(i) => {
                let d = &mut self.decisions[i];
                d.linear = d.linear.max(linear);
                d.angular = d.angular.max(angular);
            }
            Err(i) => self.decisions.insert(
                i,
                ToleranceDecision {
                    kind,
                    linear,
                    angular,
                },
            ),
        }
    }
}

/// The rounding floor relative to the operands' extent: `2^-40` (#236,
/// #244).
///
/// The absolute linear floor of a boolean is this factor times the
/// operands' extent; [`BooleanReport::rounding_floor`] returns it and
/// [`BooleanReport::extent`] defines the extent. Angular residues have a
/// floor of this factor itself, in radians, unscaled.
pub const ROUNDING_FACTOR: Scalar = 1.0 / 1_099_511_627_776.0;

/// The relative rounding floor ([`ROUNDING_FACTOR`]).
pub(crate) const ROUNDING: Scalar = ROUNDING_FACTOR;

struct Session {
    /// The caller's tolerance: every recorded decision lies within it, and
    /// at [`Tolerance::ZERO`] none may be recorded (#251).
    tolerance: Tolerance,
    /// Whether a decision the tolerance does not admit was asked for: the
    /// entry point refuses ([`BooleanError::ToleranceExceeded`]).
    exceeded: bool,
    /// The operands' extent ([`BooleanReport::extent`]).
    extent: Scalar,
    /// The linear rounding floor: [`ROUNDING`] times the operands' extent.
    rounding: Scalar,
    report: BooleanReport,
}

thread_local! {
    static SESSION: RefCell<Option<Session>> = const { RefCell::new(None) };
}

/// The open session of one entry point; dropping it closes the session.
pub(crate) struct Guard {
    owner: bool,
}

impl Guard {
    /// The report of the session, closing it. A guard that joined an outer
    /// session returns an empty report: the outer one holds the decisions.
    ///
    /// # Errors
    ///
    /// [`BooleanError::ToleranceExceeded`] when a decision the caller's
    /// tolerance does not admit was asked for (see [`record`]).
    pub(crate) fn finish(self) -> Result<BooleanReport, BooleanError> {
        if !self.owner {
            return Ok(BooleanReport::default());
        }
        SESSION
            .with(|s| {
                s.borrow_mut().as_mut().map(|s| {
                    if s.exceeded {
                        return Err(BooleanError::ToleranceExceeded);
                    }
                    Ok(BooleanReport {
                        extent: s.extent,
                        ..std::mem::take(&mut s.report)
                    })
                })
            })
            .unwrap_or_else(|| Ok(BooleanReport::default()))
    }
}

impl Drop for Guard {
    fn drop(&mut self) {
        if self.owner {
            SESSION.with(|s| *s.borrow_mut() = None);
        }
    }
}

/// Open a session for the given operands at the caller's tolerance, or
/// join the one already open.
pub(crate) fn open(operands: &[&ExactBRep], tolerance: Tolerance) -> Guard {
    let owner = SESSION.with(|s| s.borrow().is_none());
    if owner {
        let extent = operands.iter().map(|b| extent(b)).fold(0.0, Scalar::max);
        SESSION.with(|s| {
            *s.borrow_mut() = Some(Session {
                tolerance,
                exceeded: false,
                extent,
                rounding: ROUNDING * extent,
                report: BooleanReport::default(),
            });
        });
    }
    Guard { owner }
}

/// The largest coordinate magnitude of a B-rep: its vertices, and each
/// surface's origin widened by its size.
fn extent(brep: &ExactBRep) -> Scalar {
    let mut out: Scalar = 0.0;
    for v in brep.topology().vertices() {
        out = out.max(v.position.abs().max_element());
    }
    for s in brep.surfaces() {
        let (origin, size) = match s {
            Surface::Plane(p) => (p.frame.origin, 0.0),
            Surface::Cylinder(c) => (c.frame.origin, c.radius),
            Surface::EllipticalCylinder(c) => (c.frame.origin, c.semi_axis_x.max(c.semi_axis_y)),
            Surface::Cone(c) => (c.frame.origin, c.radius.abs()),
            Surface::Sphere(c) => (c.frame.origin, c.radius),
            Surface::Torus(t) => (t.frame.origin, t.major_radius + t.minor_radius),
            _ => continue,
        };
        out = out.max(origin.abs().max_element() + size.abs());
    }
    if out.is_finite() {
        out
    } else {
        0.0
    }
}

/// The linear rounding floor of the open session (zero without one).
pub(crate) fn rounding() -> Scalar {
    SESSION.with(|s| s.borrow().as_ref().map_or(0.0, |s| s.rounding))
}

/// The tolerance with its parts raised to the rounding floor: for naming a
/// rounded point on the curve or surface it was evaluated from, which is
/// bookkeeping, not a decision.
pub(crate) fn floored(tolerance: Tolerance) -> Tolerance {
    Tolerance::new(
        tolerance.linear().max(rounding()),
        tolerance.angular().max(ROUNDING),
    )
    .unwrap_or(tolerance)
}

/// Record a within-tolerance reading in the open session.
///
/// A reading must lie within the session's tolerance, and at
/// [`Tolerance::ZERO`] there is none (#236, #251): the helpers below never
/// ask otherwise. A call that does is a bug -- a debug assertion -- and in
/// a release build marks the session, whose entry point then refuses
/// ([`Guard::finish`]) instead of returning a result the report
/// misdescribes.
pub(crate) fn record(kind: ToleranceDecisionKind, linear: Scalar, angular: Scalar) {
    SESSION.with(|s| {
        if let Some(s) = s.borrow_mut().as_mut() {
            let admitted = s.tolerance != Tolerance::ZERO
                && linear <= s.tolerance.linear()
                && angular <= s.tolerance.angular();
            // The crate's own unit tests take the release build's refusal.
            #[cfg(not(test))]
            debug_assert!(
                admitted,
                "{kind:?} ({linear}, {angular}) recorded at tolerance {:?} (#251)",
                s.tolerance
            );
            if admitted {
                s.report.add(kind, linear, angular);
            } else {
                s.exceeded = true;
            }
        }
    });
}

/// A point decision on a linear residue: within the rounding floor it is
/// one point, unrecorded; within the caller's tolerance it is recorded.
pub(crate) fn near(kind: ToleranceDecisionKind, residue: Scalar, tolerance: Tolerance) -> bool {
    let within = residue <= tolerance.linear().max(rounding());
    if !within {
        return false;
    }
    if residue > rounding() {
        record(kind, residue, 0.0);
    }
    true
}

/// A point decision on an angular residue, as [`near`].
pub(crate) fn near_angle(
    kind: ToleranceDecisionKind,
    residue: Scalar,
    tolerance: Tolerance,
) -> bool {
    let within = residue <= tolerance.angular().max(ROUNDING);
    if !within {
        return false;
    }
    if residue > ROUNDING {
        record(kind, 0.0, residue);
    }
    true
}

/// How a support decision came out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Reading {
    /// The exact predicate holds: no tolerance was used.
    Exact,
    /// Not exactly, but within the caller's tolerance: recorded.
    Within,
    /// Neither.
    Apart,
}

impl Reading {
    /// Exactly, or within tolerance.
    pub(crate) fn holds(self) -> bool {
        self != Reading::Apart
    }
}

/// A support decision: `linear` and `angular` are the `f64` measures of how
/// far the reading moves and turns the operands. Within the caller's
/// tolerance (or within rounding of zero), the exact predicate is asked:
/// confirmed, the reading is exact and unrecorded; otherwise it is taken
/// within tolerance and recorded, or refused beyond it.
///
/// A reading the exact predicate rejects moves or turns the operands by a
/// positive amount even when its `f64` measure is exactly `0` (a rotation's
/// rounding cancels in `f64`, not in the numbers given). A tolerance part of
/// zero admits no positive amount, so such a reading is refused there,
/// whatever its measure: at [`Tolerance::ZERO`] the exact answer stands
/// (#251).
pub(crate) fn support(
    kind: ToleranceDecisionKind,
    linear: Scalar,
    angular: Scalar,
    tolerance: Tolerance,
    exact: impl FnOnce() -> bool,
) -> Reading {
    let (eps, alpha) = (tolerance.linear(), tolerance.angular());
    if !(linear <= eps.max(rounding()) && angular <= alpha.max(ROUNDING)) {
        return Reading::Apart;
    }
    if exact() {
        return Reading::Exact;
    }
    if eps > 0.0 && alpha > 0.0 && linear <= eps && angular <= alpha {
        record(kind, linear, angular);
        return Reading::Within;
    }
    Reading::Apart
}

#[cfg(test)]
mod tests {
    use super::*;
    use ToleranceDecisionKind::*;

    #[test]
    fn a_report_keeps_each_kind_once_with_its_worst_magnitude() {
        let mut r = BooleanReport::default();
        assert!(r.is_exact());
        r.add(MergedPoints, 1e-9, 0.0);
        r.add(CoincidentSupports, 2e-9, 1e-12);
        r.add(MergedPoints, 3e-9, 0.0);
        assert_eq!(r.decisions().len(), 2);
        assert_eq!(r.decisions()[0].kind, CoincidentSupports);
        assert_eq!(r.linear(), 3e-9);
        assert_eq!(r.angular(), 1e-12);
        let merged = BooleanReport::default().merged(&r);
        assert_eq!(merged, r);
    }

    #[test]
    fn merged_reports_keep_the_larger_rounding_floor() {
        let small = BooleanReport {
            extent: 3.0,
            ..BooleanReport::default()
        };
        let large = BooleanReport {
            extent: 10.0,
            ..BooleanReport::default()
        };
        for merged in [small.clone().merged(&large), large.clone().merged(&small)] {
            assert_eq!(merged.extent(), 10.0);
            assert_eq!(merged.rounding_floor(), ROUNDING_FACTOR * 10.0);
            assert!(merged.is_exact());
        }
        assert_eq!(BooleanReport::default().rounding_floor(), 0.0);
    }

    #[test]
    fn a_session_reports_the_floor_it_applied() {
        let guard = open(&[], Tolerance::ZERO);
        assert_eq!(rounding(), 0.0);
        let report = guard.finish().expect("nothing recorded");
        assert_eq!((report.extent(), report.rounding_floor()), (0.0, 0.0));
        assert_eq!(ROUNDING_FACTOR, (2.0 as Scalar).powi(-40));
    }

    #[test]
    fn readings_within_rounding_are_unrecorded_and_beyond_tolerance_refused() {
        let tol = Tolerance::METRE;
        let guard = open(&[], tol);
        // No operands: the floor is zero, so any residue is a decision.
        assert!(near(MergedPoints, 0.0, tol));
        assert!(near(MergedPoints, 5e-7, tol));
        assert!(!near(MergedPoints, 2e-6, tol));
        assert!(!near(MergedPoints, Scalar::NAN, tol));
        // Exact support readings are not recorded; others are.
        let exact = support(CoincidentSupports, 1e-7, 0.0, tol, || true);
        assert_eq!(exact, Reading::Exact);
        let within = support(PlaneParallelToAxis, 1e-7, 0.0, tol, || false);
        assert_eq!(within, Reading::Within);
        let apart = support(PlaneParallelToAxis, 1e-5, 0.0, tol, || true);
        assert_eq!(apart, Reading::Apart);
        let report = guard.finish().expect("within tolerance");
        assert_eq!(report.decisions().len(), 2);
        assert!(report.contains(MergedPoints));
        assert!(report.contains(PlaneParallelToAxis));
        assert!(!report.contains(CoincidentSupports));
    }

    #[test]
    fn a_reading_the_exact_predicate_rejects_is_never_taken_at_a_zero_part() {
        // #251: an `f64` measure of exactly zero is not exact. At
        // `Tolerance::ZERO`, or with a zero part, only the exact predicate
        // can confirm a reading; nothing is recorded.
        let partial = [
            Tolerance::ZERO,
            Tolerance::new(1e-6, 0.0).expect("valid"),
            Tolerance::new(0.0, 1e-9).expect("valid"),
        ];
        for tol in partial {
            let guard = open(&[], tol);
            for kind in [
                CoincidentSupports,
                PlaneParallelToAxis,
                PlanePerpendicularToAxis,
                PlaneTouchesCylinder,
            ] {
                assert_eq!(support(kind, 0.0, 0.0, tol, || false), Reading::Apart);
                assert_eq!(support(kind, 0.0, 0.0, tol, || true), Reading::Exact);
            }
            let report = guard.finish().expect("nothing recorded");
            assert!(report.is_exact(), "{tol:?}: {report:?}");
        }
        // Both parts positive: a zero measure the predicate rejects is
        // still taken, and reported (unchanged).
        let guard = open(&[], Tolerance::METRE);
        let within = support(PlanePerpendicularToAxis, 0.0, 0.0, Tolerance::METRE, || {
            false
        });
        assert_eq!(within, Reading::Within);
        let report = guard.finish().expect("within tolerance");
        assert!(report.contains(PlanePerpendicularToAxis));
    }

    #[test]
    fn points_within_rounding_are_never_recorded_at_zero_tolerance() {
        let guard = open(&[], Tolerance::ZERO);
        assert!(near(MergedPoints, 0.0, Tolerance::ZERO));
        assert!(!near(MergedPoints, 1e-300, Tolerance::ZERO));
        assert!(near_angle(Contact, 0.0, Tolerance::ZERO));
        assert!(!near_angle(Contact, 2.0 * ROUNDING, Tolerance::ZERO));
        assert!(guard.finish().expect("nothing recorded").is_exact());
    }

    #[test]
    fn a_decision_recorded_at_zero_tolerance_is_refused() {
        // The helpers never record at `Tolerance::ZERO`; a direct record is
        // a bug: a debug assertion outside these tests, and a named refusal.
        let guard = open(&[], Tolerance::ZERO);
        record(PlanePerpendicularToAxis, 0.0, 0.0);
        assert_eq!(guard.finish(), Err(BooleanError::ToleranceExceeded));
    }

    #[test]
    fn a_decision_beyond_the_tolerance_is_refused() {
        let guard = open(&[], Tolerance::METRE);
        record(MergedPoints, 2e-6, 0.0);
        assert_eq!(guard.finish(), Err(BooleanError::ToleranceExceeded));
    }

    #[test]
    fn a_joined_session_keeps_the_outer_tolerance() {
        let outer = open(&[], Tolerance::METRE);
        let inner = open(&[], Tolerance::ZERO);
        // The outer session's tolerance applies to a joined entry point.
        assert!(near(MergedPoints, 5e-7, Tolerance::METRE));
        assert!(inner.finish().expect("joined").is_exact());
        assert!(outer.finish().expect("within").contains(MergedPoints));
    }
}
