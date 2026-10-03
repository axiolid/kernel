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
//! [`crate::section_edges`], [`crate::split_face`]); an entry point called
//! inside another's session joins it. Nothing is shared across threads and
//! nothing outlives the call.

use std::cell::RefCell;

use axiolid_brep::ExactBRep;
use axiolid_core::{Scalar, Tolerance};
use axiolid_surface::Surface;

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
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BooleanReport {
    decisions: Vec<ToleranceDecision>,
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

    /// This report and `other` together: what a chain of booleans read
    /// (each kind once, with the worst magnitude of either).
    #[must_use]
    pub fn merged(mut self, other: &BooleanReport) -> BooleanReport {
        for d in &other.decisions {
            self.add(d.kind, d.linear, d.angular);
        }
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

/// The rounding floor, relative to the operands' extent: `2^-40`.
pub(crate) const ROUNDING: Scalar = 1.0 / 1_099_511_627_776.0;

struct Session {
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
    pub(crate) fn finish(self) -> BooleanReport {
        if !self.owner {
            return BooleanReport::default();
        }
        SESSION
            .with(|s| {
                s.borrow_mut()
                    .as_mut()
                    .map(|s| std::mem::take(&mut s.report))
            })
            .unwrap_or_default()
    }
}

impl Drop for Guard {
    fn drop(&mut self) {
        if self.owner {
            SESSION.with(|s| *s.borrow_mut() = None);
        }
    }
}

/// Open a session for the given operands, or join the one already open.
pub(crate) fn open(operands: &[&ExactBRep]) -> Guard {
    let owner = SESSION.with(|s| s.borrow().is_none());
    if owner {
        let extent = operands.iter().map(|b| extent(b)).fold(0.0, Scalar::max);
        SESSION.with(|s| {
            *s.borrow_mut() = Some(Session {
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
pub(crate) fn record(kind: ToleranceDecisionKind, linear: Scalar, angular: Scalar) {
    SESSION.with(|s| {
        if let Some(s) = s.borrow_mut().as_mut() {
            s.report.add(kind, linear, angular);
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
    if linear <= eps && angular <= alpha {
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
    fn readings_within_rounding_are_unrecorded_and_beyond_tolerance_refused() {
        let guard = open(&[]);
        let tol = Tolerance::METRE;
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
        let report = guard.finish();
        assert_eq!(report.decisions().len(), 2);
        assert!(report.contains(MergedPoints));
        assert!(report.contains(PlaneParallelToAxis));
        assert!(!report.contains(CoincidentSupports));
    }
}
