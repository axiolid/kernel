//! Conformance checks for providers of
//! [`MeshPlaneSection::section_request`](crate::MeshPlaneSection::section_request).
//!
//! The suite runs analytic fixtures with hand-derived expectations
//! (`fixtures.rs`): faces in the plane, a post rising out of a parapet whose
//! top lies in the plane, edge and point contacts, and overlapping in-plane
//! faces. A provider whose first case returns `Unsupported` does not claim
//! the method: the report passes with `claimed == false`. Every case runs
//! twice and the results must be identical.

use core::fmt;

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::{Point2, Tolerance};

use crate::conformance::{ConformanceFailure, ConformanceSuite};
use crate::fixtures::{cases, Case, Expected, Output};
use crate::request::{
    InPlaneRegion, PlaneSection, PlaneSectionEvidence, RegionLimits, SectionChain, SectionRequest,
};
use crate::{MeshPlaneSection, SectionLimits};

/// Absolute tolerance for comparing output points with the fixtures.
const POINT_TOLERANCE: f64 = 1e-9;

/// One conformance failure of the region suite.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum RegionConformanceFailure {
    /// The provider's legacy `section` fails the base suite.
    Base(ConformanceFailure),
    /// A case returned an error that is not its expected refusal.
    ReturnedError {
        /// Fixture name.
        case: &'static str,
        /// The error.
        error: String,
    },
    /// The chains differ from the fixture's.
    WrongChains {
        /// Fixture name.
        case: &'static str,
    },
    /// The regions differ from the fixture's.
    WrongRegions {
        /// Fixture name.
        case: &'static str,
    },
    /// The evidence does not match the fixture and the output.
    WrongEvidence {
        /// Fixture name.
        case: &'static str,
    },
    /// A case must be refused and was not, or was refused otherwise.
    WrongRefusal {
        /// Fixture name.
        case: &'static str,
        /// The expected input family.
        expected: &'static str,
        /// What happened instead.
        actual: String,
    },
    /// Two runs of one case differ.
    NonDeterministic {
        /// Fixture name.
        case: &'static str,
    },
}

/// Result of [`RegionConformanceSuite::run`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct RegionConformanceReport {
    /// Every failure found.
    pub failures: Vec<RegionConformanceFailure>,
    /// Whether the provider implements `section_request`. When `false` the
    /// suite ran nothing else and passes.
    pub claimed: bool,
}

impl RegionConformanceReport {
    /// Whether no case failed.
    pub fn is_success(&self) -> bool {
        self.failures.is_empty()
    }
}

impl fmt::Display for RegionConformanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if !self.claimed {
            write!(f, "not claimed")
        } else if self.is_success() {
            write!(f, "conformant")
        } else {
            write!(f, "{:?}", self.failures)
        }
    }
}

/// The conformance suite for `section_request`.
pub struct RegionConformanceSuite;

impl RegionConformanceSuite {
    /// Run every fixture against `provider`.
    pub fn run(provider: &dyn MeshPlaneSection) -> RegionConformanceReport {
        let mut report = RegionConformanceReport {
            failures: Vec::new(),
            claimed: true,
        };
        let limits = SectionLimits::new(4096, 8192, 16384, 64);
        let region_limits = RegionLimits::new(64, 4096, 1 << 20);
        let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"));
        for (index, case) in cases().iter().enumerate() {
            let request = SectionRequest::new(case.frame, limits, region_limits);
            let first = provider.section_request(&case.mesh, &request, &options);
            if index == 0 && matches!(first, Err(GeomError::Unsupported { .. })) {
                report.claimed = false;
                return report;
            }
            let second = provider.section_request(&case.mesh, &request, &options);
            check_case(case, first, second, &mut report);
        }
        for failure in ConformanceSuite::run(provider).failures {
            report
                .failures
                .push(RegionConformanceFailure::Base(failure));
        }
        report
    }
}

fn check_case(
    case: &Case,
    first: Result<PlaneSection, GeomError>,
    second: Result<PlaneSection, GeomError>,
    report: &mut RegionConformanceReport,
) {
    let name = case.name;
    let stable = match (&first, &second) {
        (Ok(a), Ok(b)) => a == b,
        (Err(a), Err(b)) => a.to_string() == b.to_string(),
        _ => false,
    };
    if !stable {
        report
            .failures
            .push(RegionConformanceFailure::NonDeterministic { case: name });
    }
    match (&case.expected, first) {
        (Expected::Output(expected), Ok(got))
        | (Expected::OutputOrRefused(_, expected), Ok(got)) => {
            check_output(case, expected, &got, report);
        }
        (Expected::OutputOrRefused(family, _), Err(error)) if is_refusal(&error, family) => {}
        (Expected::Output(_), Err(error)) => {
            report
                .failures
                .push(RegionConformanceFailure::ReturnedError {
                    case: name,
                    error: error.to_string(),
                });
        }
        (Expected::OutputOrRefused(family, _), Err(error)) => {
            report
                .failures
                .push(RegionConformanceFailure::WrongRefusal {
                    case: name,
                    expected: family,
                    actual: error.to_string(),
                });
        }
    }
}

fn is_refusal(error: &GeomError, family: &str) -> bool {
    matches!(
        error,
        GeomError::UnsupportedInput { operation: Operation::MeshPlaneSection, input, .. }
            if *input == family
    )
}

fn check_output(
    case: &Case,
    expected: &Output,
    got: &PlaneSection,
    report: &mut RegionConformanceReport,
) {
    let name = case.name;
    if got.frame != case.frame
        || !chains_match(&expected.chains, &got.chains)
        || !open_ends_are_region_vertices(got)
    {
        report
            .failures
            .push(RegionConformanceFailure::WrongChains { case: name });
    }
    if !regions_match(&expected.regions, &got.regions) {
        report
            .failures
            .push(RegionConformanceFailure::WrongRegions { case: name });
    }
    let evidence = got.evidence;
    let derived = PlaneSectionEvidence::input_mesh(case.mesh.triangle_count())
        .with_source_counts(
            expected.in_plane_triangles,
            expected.touching_edges,
            evidence.in_plane_pairs,
        )
        .with_output(&got.chains, &got.regions);
    if evidence != derived || !evidence.is_derived_from_input_mesh() {
        report
            .failures
            .push(RegionConformanceFailure::WrongEvidence { case: name });
    }
}

/// Both end points of every open chain equal a region vertex bit for bit:
/// chains and regions are projected by one function.
fn open_ends_are_region_vertices(section: &PlaneSection) -> bool {
    let on_region = |point: &Point2| {
        section.regions.iter().any(|region| {
            std::iter::once(&region.outer)
                .chain(&region.holes)
                .any(|ring| ring.iter().any(|vertex| vertex == point))
        })
    };
    section
        .chains
        .iter()
        .filter(|chain| !chain.is_closed())
        .all(|chain| {
            chain.points.first().is_some_and(on_region)
                && chain.points.last().is_some_and(on_region)
        })
}

fn chains_match(expected: &[SectionChain], got: &[SectionChain]) -> bool {
    expected.len() == got.len()
        && expected
            .iter()
            .zip(got)
            .all(|(e, g)| e.ends == g.ends && points_match(&e.points, &g.points))
}

fn regions_match(expected: &[InPlaneRegion], got: &[InPlaneRegion]) -> bool {
    expected.len() == got.len()
        && expected.iter().zip(got).all(|(e, g)| {
            e.side == g.side
                && points_match(&e.outer, &g.outer)
                && e.holes.len() == g.holes.len()
                && e.holes
                    .iter()
                    .zip(&g.holes)
                    .all(|(eh, gh)| points_match(eh, gh))
        })
}

fn points_match(expected: &[Point2], got: &[Point2]) -> bool {
    expected.len() == got.len()
        && expected.iter().zip(got).all(|(e, g)| {
            (e.x - g.x).abs() <= POINT_TOLERANCE && (e.y - g.y).abs() <= POINT_TOLERANCE
        })
}
