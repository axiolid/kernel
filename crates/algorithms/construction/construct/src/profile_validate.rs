//! Exact validation of profiles on the contour the exact path builds (#352).
//!
//! [`ExactProfileValidator`] implements [`ValidateProfile`] for
//! [`ContourProfile`] and [`Profile`]. It decides closure, self-contact and
//! hole nesting, and refuses with a typed [`ProfileFinding`] instead of
//! guessing.
//!
//! # Why the lowered ring, not the authored curves
//!
//! A `ProfileSegment` stores a `Circle2` (an `f64` frame and radius) and an
//! `f64` angle domain, so an arc end is `centre + r (cos t, sin t)`. By the
//! Lindemann-Weierstrass theorem `cos t` and `sin t` are transcendental for
//! every nonzero dyadic `t`: an arc-to-arc joint almost never closes exactly,
//! and an arc end almost never lies exactly on a line. An exact predicate on
//! the literal data would answer "open contour" and "off the axis", both
//! correctly and both uselessly. Intent can only be recovered through a
//! declared interpretation, and the exact path already has one:
//! [`contour_to_arc_ring`](crate::contour_lower::contour_to_arc_ring).
//!
//! - A joint is ONE `f64` point, the start of the segment leaving it.
//! - A joint closes when the entering end lies within
//!   `max(tolerance.linear(), 8 eps * operand magnitudes)` of it.
//! - Each arc becomes a bulge, split into equal sub-arcs of less than half a
//!   turn.
//!
//! The extruder, the revolver and the half-space path build from that ring.
//! So the validator judges the lowered ring and "valid" means valid as the
//! exact path will build it. Each ring edge maps back to the segment it came
//! from, and findings name segments, not ring edges.
//!
//! # Exact, then tolerance-free
//!
//! Self-contact, shared boundary and hole nesting are decided on the lowered
//! ring by `axiolid_overlay::check_arc_polygon`: exact signs over the ring's
//! `f64` data, no tolerance in any contact or nesting decision (ADR 0070).
//! The tolerance enters only where it does in the exact path: joint closure,
//! and the per-ring degeneracies of `validate_arc_ring` (repeated vertex,
//! zero-radius arc, zero area).
//!
//! What the lowering moves is reported, not hidden:
//! [`ProfileReport::max_joint_gap`] is the largest joint gap it closed. The
//! lowered boundary lies within `2 * max_joint_gap` plus rounding of the
//! authored one:
//!
//! - Each lowered arc is a similarity image of the authored arc, its two
//!   ends moved by at most the gap at their joints.
//! - For a sub-arc of less than half a turn, a point of the image moves by
//!   at most `|d0| + |d1 - d0|`, where `d0` and `d1` are the movements of
//!   its ends; both are at most the gap, which gives the bound.
//!
//! Two parts of the authored boundary that come closer than that, without
//! touching, are decided for the lowered ring.
//!
//! # Near-contacts (opt-in, #388)
//!
//! Parts that miss each other by less than a distance are valid exactly.
//! A consumer that treats them as touching asks
//! [`ExactProfileValidator::validate_with`] for a clearance (usually
//! `tolerance.linear()`) and gets a [`ProfileFinding::Clearance`] naming the
//! two segments, with a certified lower bound on their distance. A pair is
//! reported when its distance has a local minimum below the clearance away
//! from where they meet, or when it forms an exact cusp (they meet and leave
//! in the same direction); the neighbourhood of a meeting point at a nonzero
//! angle, a narrow wedge, is not. The check is
//! `axiolid_overlay::check_arc_polygon_with`'s: exact, and a pair it cannot
//! decide is reported with `undecided` set, never cleared. It is off in
//! [`ValidateProfile::validate`], so the extruder and every existing caller
//! see the same report as before. The lowered ring is judged, as everywhere
//! here: the authored parts may be `2 * max_joint_gap` closer or farther.
//!
//! # Refused and reported
//!
//! Refused with a finding, checking nothing more of that contour: too few
//! segments, a segment family the lowering refuses (ellipse, B-spline,
//! chains, ...), a degenerate frame, an empty parameter range, more than a
//! whole turn, and an open joint. A contour that does not lower makes the
//! profile report only those findings, because nesting against a missing
//! ring would be a wrong answer. TODO(#389): polyline segments are refused
//! like the other families although they lower exactly.
//!
//! Reported, with the segments involved: a self-contact (crossing, pinch,
//! tangency or fold-back), a shared stretch of boundary between contours, a
//! hole outside the outer contour, and holes whose interiors meet. An
//! isolated point where two contours meet is a [`ProfileContact`], not a
//! finding: [`ProfileReport::is_valid`] accepts it and
//! [`ProfileReport::is_strict`] refuses it. Contacts that together cut the
//! interior apart (a hole touching the outer contour at several points, a
//! closed chain of holes) are a `DisconnectedInterior` finding.

use std::collections::BTreeMap;

use axiolid_contracts::GeomError;
use axiolid_core::{Point2, Scalar, Tolerance};
use axiolid_overlay::{
    check_arc_polygon_with, ArcCheckOptions, ArcPolygon, ArcPolygonDefect, ArcRingRef, OverlayError,
};
use axiolid_profile::{Contour, ContourProfile, Profile, ValidateProfile};

use crate::contour_lower::{lower_contour, ContourFault, LoweredContour};
use crate::extrude_exact::profile_to_contour;

/// Which contour of a [`ContourProfile`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ContourRef {
    Outer,
    Hole(usize),
}

impl From<ArcRingRef> for ContourRef {
    fn from(ring: ArcRingRef) -> Self {
        match ring {
            ArcRingRef::Outer => Self::Outer,
            ArcRingRef::Hole(k) => Self::Hole(k),
        }
    }
}

/// Why a profile is not valid as the exact path would build it.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum ProfileFinding {
    /// The profile does not lower to a contour (family or dimensions);
    /// nothing else was checked.
    NotLowered(GeomError),
    TooFewSegments {
        contour: ContourRef,
        count: usize,
    },
    /// A segment family the exact lowering refuses (ellipse, polyline,
    /// B-spline, ...), or an arc of more than a whole turn: `kind` names it.
    UnsupportedSegment {
        contour: ContourRef,
        segment: usize,
        kind: &'static str,
    },
    /// Degenerate frame or empty parameter range.
    DegenerateSegment {
        contour: ContourRef,
        segment: usize,
        reason: &'static str,
    },
    /// The end of `segment` and the start of the next (`(segment + 1) % n`,
    /// so `n - 1` is the closing joint) are `gap` apart, beyond the slack.
    Gap {
        contour: ContourRef,
        segment: usize,
        gap: Scalar,
    },
    /// A `validate_arc_ring` refusal of the lowered ring.
    DegenerateRing {
        contour: ContourRef,
        error: OverlayError,
    },
    /// Two segments of one contour share a point other than the joint
    /// between consecutive ones. `segments` are segment indices, ascending;
    /// both name one segment when a single arc touches itself. `at` is the
    /// first such point, rounded once.
    SelfContact {
        contour: ContourRef,
        segments: [usize; 2],
        at: Point2,
    },
    /// Two contours share a stretch of boundary. `contours` ascending.
    BoundaryOverlap {
        contours: [ContourRef; 2],
    },
    /// Part of a hole lies outside the outer contour.
    HoleOutsideOuter {
        hole: usize,
    },
    /// Contacts between contours cut the interior into separate pieces.
    /// `contours` ascending, every contour of the cut component.
    DisconnectedInterior {
        contours: Vec<ContourRef>,
    },
    /// A defect of the arc polygon check this version does not map to a
    /// finding of its own, such as `NestingUndecided`: nesting could not be
    /// decided, so the profile is not valid.
    Polygon(ArcPolygonDefect),
    /// Two holes' interiors meet. `holes` ascending.
    HolesOverlap {
        holes: [usize; 2],
    },
    /// Two segments are closer than the requested clearance
    /// ([`ExactProfileValidator::validate_with`]): their distance has a
    /// local minimum strictly between zero and it away from their meeting
    /// points, or they form an exact cusp (they meet and leave in the same
    /// direction). The neighbourhood of a meeting point at a nonzero angle (a
    /// narrow wedge) is not reported. `segments` are
    /// `(contour, segment)` pairs, ascending. `lower_bound` is a certified
    /// lower bound on their distance; `undecided` is false when the exact
    /// test proved the distance is below the clearance and true when it
    /// could not decide and the pair is reported rather than cleared.
    Clearance {
        segments: [(ContourRef, usize); 2],
        lower_bound: Scalar,
        undecided: bool,
    },
}

/// A point where two different contours meet.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileContact {
    /// The two contours, ascending.
    pub contours: [ContourRef; 2],
    /// The meeting point, rounded once.
    pub at: Point2,
}

/// The result of [`ExactProfileValidator`].
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ProfileReport {
    pub findings: Vec<ProfileFinding>,
    pub contacts: Vec<ProfileContact>,
    /// Largest joint gap closed by lowering; the validated boundary lies
    /// within `2 * max_joint_gap` plus rounding of the authored one.
    pub max_joint_gap: Scalar,
}

impl ProfileReport {
    /// No finding. Isolated contacts between contours are allowed.
    pub fn is_valid(&self) -> bool {
        self.findings.is_empty()
    }

    /// No finding and no contact between contours.
    pub fn is_strict(&self) -> bool {
        self.findings.is_empty() && self.contacts.is_empty()
    }
}

/// Validates a profile on the arc ring the exact path lowers it to.
///
/// See the module documentation for what is decided and what is refused.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ExactProfileValidator;

impl ExactProfileValidator {
    /// [`ValidateProfile::validate`] with options: the default options are
    /// exactly `validate`; `ArcCheckOptions::with_clearance` adds
    /// [`ProfileFinding::Clearance`] for segments closer than the distance.
    pub fn validate_with(
        &self,
        profile: &Profile,
        tolerance: Tolerance,
        options: ArcCheckOptions,
    ) -> ProfileReport {
        match profile_to_contour(profile, tolerance) {
            Ok(contours) => validate_contours_with(&contours, tolerance, options),
            Err(error) => ProfileReport {
                findings: vec![ProfileFinding::NotLowered(error)],
                ..ProfileReport::default()
            },
        }
    }

    /// [`ExactProfileValidator::validate_with`] for a [`ContourProfile`].
    pub fn validate_contours_with(
        &self,
        profile: &ContourProfile,
        tolerance: Tolerance,
        options: ArcCheckOptions,
    ) -> ProfileReport {
        validate_contours_with(profile, tolerance, options)
    }
}

impl ValidateProfile<ContourProfile> for ExactProfileValidator {
    type Report = ProfileReport;

    fn validate(&self, profile: &ContourProfile, tolerance: Tolerance) -> ProfileReport {
        validate_contours(profile, tolerance)
    }
}

impl ValidateProfile<Profile> for ExactProfileValidator {
    type Report = ProfileReport;

    fn validate(&self, profile: &Profile, tolerance: Tolerance) -> ProfileReport {
        self.validate_with(profile, tolerance, ArcCheckOptions::new())
    }
}

pub(crate) fn validate_contours(profile: &ContourProfile, tolerance: Tolerance) -> ProfileReport {
    validate_contours_with(profile, tolerance, ArcCheckOptions::new())
}

fn validate_contours_with(
    profile: &ContourProfile,
    tolerance: Tolerance,
    options: ArcCheckOptions,
) -> ProfileReport {
    let mut report = ProfileReport::default();
    let mut lowered: Vec<LoweredContour> = Vec::with_capacity(1 + profile.holes.len());
    let mut broken = false;
    let contours: Vec<(ContourRef, &Contour)> =
        std::iter::once((ContourRef::Outer, &profile.outer))
            .chain(
                profile
                    .holes
                    .iter()
                    .enumerate()
                    .map(|(k, hole)| (ContourRef::Hole(k), hole)),
            )
            .collect();
    for (name, contour) in contours {
        match lower_contour(contour, tolerance) {
            Ok(done) => {
                report.max_joint_gap = report.max_joint_gap.max(done.max_gap);
                lowered.push(done);
            }
            Err(fault) => {
                report.findings.push(finding_of(name, fault));
                broken = true;
            }
        }
    }
    if broken {
        return report;
    }

    let mut rings = lowered.iter().map(|l| l.ring.clone());
    let polygon = ArcPolygon {
        outer: rings.next().expect("the outer contour lowered"),
        holes: rings.collect(),
    };
    let check = check_arc_polygon_with(&polygon, tolerance, options);
    let segment = |ring: ArcRingRef, edge: usize| match ring {
        ArcRingRef::Outer => lowered[0].segment_of_edge[edge],
        ArcRingRef::Hole(k) => lowered[k + 1].segment_of_edge[edge],
    };
    // Clearance defects name edges; several edges of one arc make one
    // finding, so they are merged per segment pair and sorted.
    let mut near: BTreeMap<[(ContourRef, usize); 2], (Scalar, bool)> = BTreeMap::new();
    for defect in check.defects {
        report.findings.push(match defect {
            ArcPolygonDefect::Clearance {
                edges,
                lower_bound,
                undecided,
            } => {
                let mut pair = edges.map(|(ring, edge)| (ring.into(), segment(ring, edge)));
                pair.sort_unstable();
                let entry = near.entry(pair).or_insert((lower_bound, undecided));
                entry.0 = entry.0.min(lower_bound);
                entry.1 &= undecided;
                continue;
            }
            ArcPolygonDefect::Ring { ring, error } => ProfileFinding::DegenerateRing {
                contour: ring.into(),
                error,
            },
            ArcPolygonDefect::SelfContact { ring, edges, at } => {
                let mut segments = edges.map(|edge| segment(ring, edge));
                segments.sort_unstable();
                ProfileFinding::SelfContact {
                    contour: ring.into(),
                    segments,
                    at,
                }
            }
            ArcPolygonDefect::BoundaryOverlap { rings } => ProfileFinding::BoundaryOverlap {
                contours: rings.map(ContourRef::from),
            },
            ArcPolygonDefect::HoleOutsideOuter { hole } => {
                ProfileFinding::HoleOutsideOuter { hole }
            }
            ArcPolygonDefect::HolesOverlap { holes } => ProfileFinding::HolesOverlap { holes },
            ArcPolygonDefect::DisconnectedInterior { rings } => {
                ProfileFinding::DisconnectedInterior {
                    contours: rings.into_iter().map(ContourRef::from).collect(),
                }
            }
            // `ArcPolygonDefect` is non-exhaustive: a defect this version
            // does not know is carried whole, never dropped.
            other => ProfileFinding::Polygon(other),
        });
    }
    report.findings.extend(
        near.into_iter().map(
            |(segments, (lower_bound, undecided))| ProfileFinding::Clearance {
                segments,
                lower_bound,
                undecided,
            },
        ),
    );
    report.contacts = check
        .contacts
        .into_iter()
        .map(|contact| ProfileContact {
            contours: contact.rings.map(ContourRef::from),
            at: contact.at,
        })
        .collect();
    report
}

fn finding_of(contour: ContourRef, fault: ContourFault) -> ProfileFinding {
    match fault {
        ContourFault::TooFewSegments { count } => ProfileFinding::TooFewSegments { contour, count },
        ContourFault::UnsupportedSegment { segment, kind } => ProfileFinding::UnsupportedSegment {
            contour,
            segment,
            kind,
        },
        ContourFault::DegenerateSegment { segment, reason } => ProfileFinding::DegenerateSegment {
            contour,
            segment,
            reason,
        },
        ContourFault::Gap { segment, gap, .. } => ProfileFinding::Gap {
            contour,
            segment,
            gap,
        },
    }
}
