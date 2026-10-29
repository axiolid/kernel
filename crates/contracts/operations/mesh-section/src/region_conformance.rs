//! Conformance checks for providers of
//! [`MeshPlaneSection::section_request`](crate::MeshPlaneSection::section_request).
//!
//! The suite runs analytic fixtures with hand-derived expectations
//! (`fixtures.rs`): faces in the plane, a post rising out of a parapet whose
//! top lies in the plane, edge and point contacts, overlapping in-plane
//! faces, and sheets (terrain, a sheet in the plane, an annulus, a saddle). A provider whose first case returns `Unsupported` does not claim
//! the method: the report passes with `claimed == false`. Every case runs
//! twice and the results must be identical.
//!
//! Every case runs a second time with cut regions requested (ADR 0092). A
//! provider whose first such request returns `Unsupported`, or an answer
//! without cut regions (it predates them), does not claim them
//! (`cut_regions_claimed == false`) and the rest is skipped. Otherwise
//! the answer must equal the plain one, plus the fixture's cut regions and
//! the cut evidence; a plain refusal must be repeated; a sheet is invalid
//! input; and the cut rings must cover the chains bit for bit.

use core::fmt;

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::{Point2, Tolerance};

use crate::conformance::{ConformanceFailure, ConformanceSuite};
use crate::cut::{CutRegion, CutRegionLimits};
use crate::fixtures::{cases, Case, Cut, Expected, Output};
use crate::request::{
    ChainEnd, InPlaneRegion, PlaneSection, PlaneSectionEvidence, RegionLimits, SectionChain,
    SectionRequest, SectionRole,
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
    /// The cut regions, or the answer around them, differ from the fixture's.
    WrongCutRegions {
        /// Fixture name.
        case: &'static str,
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
    /// Whether the provider also sections `SectionRole::Sheet`. When
    /// `false`, the first sheet case returned `Unsupported` and the sheet
    /// cases were skipped; the report can still pass.
    pub sheets_claimed: bool,
    /// Whether the provider also traces cut regions. When `false`, the first
    /// cut-region request returned `Unsupported` and the cut checks were
    /// skipped; the report can still pass.
    pub cut_regions_claimed: bool,
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
            let sheets = if self.sheets_claimed {
                ""
            } else {
                " (sheets not claimed)"
            };
            let cuts = if self.cut_regions_claimed {
                ""
            } else {
                " (cut regions not claimed)"
            };
            write!(f, "conformant{sheets}{cuts}")
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
            sheets_claimed: true,
            cut_regions_claimed: true,
        };
        let mut seen_sheet = false;
        let limits = SectionLimits::new(4096, 8192, 16384, 64);
        let region_limits = RegionLimits::new(64, 4096, 1 << 20);
        let cut_limits = CutRegionLimits::new(64, 4096, 1 << 20);
        let mut seen_solid = false;
        let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"));
        for (index, case) in cases().iter().enumerate() {
            if case.role == SectionRole::Sheet && !report.sheets_claimed {
                continue;
            }
            let request =
                SectionRequest::new(case.frame, limits, region_limits).with_role(case.role);
            let first = provider.section_request(&case.mesh, &request, &options);
            if case.role == SectionRole::Sheet
                && !seen_sheet
                && matches!(first, Err(GeomError::Unsupported { .. }))
            {
                report.sheets_claimed = false;
                continue;
            }
            seen_sheet |= case.role == SectionRole::Sheet;
            if index == 0 && matches!(first, Err(GeomError::Unsupported { .. })) {
                report.claimed = false;
                return report;
            }
            let second = provider.section_request(&case.mesh, &request, &options);
            let plain = first.clone();
            check_case(case, first, second, &mut report);
            if !report.cut_regions_claimed {
                continue;
            }
            let cut_request = request.with_cut_regions(cut_limits);
            let first = provider.section_request(&case.mesh, &cut_request, &options);
            if case.role == SectionRole::Solid
                && !seen_solid
                && match &first {
                    Err(error) => matches!(error, GeomError::Unsupported { .. }),
                    // A provider that predates cut regions ignores the request.
                    Ok(section) => section.cut_regions.is_none(),
                }
            {
                report.cut_regions_claimed = false;
                continue;
            }
            seen_solid |= case.role == SectionRole::Solid;
            let second = provider.section_request(&case.mesh, &cut_request, &options);
            check_cut_case(case, &plain, first, second, &mut report);
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
        (Expected::OutputOrRefused(family, _) | Expected::Refused(family), Err(error))
            if is_refusal(&error, family) => {}
        (Expected::Refused(family), Ok(_)) => {
            report
                .failures
                .push(RegionConformanceFailure::WrongRefusal {
                    case: name,
                    expected: family,
                    actual: "a section was returned".into(),
                });
        }
        (Expected::Output(_), Err(error)) => {
            report
                .failures
                .push(RegionConformanceFailure::ReturnedError {
                    case: name,
                    error: error.to_string(),
                });
        }
        (Expected::OutputOrRefused(family, _) | Expected::Refused(family), Err(error)) => {
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

/// Check the answer to a cut-region request against the plain answer.
fn check_cut_case(
    case: &Case,
    plain: &Result<PlaneSection, GeomError>,
    first: Result<PlaneSection, GeomError>,
    second: Result<PlaneSection, GeomError>,
    report: &mut RegionConformanceReport,
) {
    let name = case.name;
    let wrong = RegionConformanceFailure::WrongCutRegions { case: name };
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
    match (&case.cut, first) {
        (Cut::Invalid, Err(GeomError::InvalidInput(_))) => {}
        (Cut::Refused(family), Err(error)) if is_refusal(&error, family) => {}
        (Cut::Regions(_), Err(error)) if plain.as_ref().is_err_and(|p| same_error(p, &error)) => {}
        (Cut::Refused(family), Err(error)) => {
            report
                .failures
                .push(RegionConformanceFailure::WrongRefusal {
                    case: name,
                    expected: family,
                    actual: error.to_string(),
                });
        }
        (_, Err(error)) => {
            report
                .failures
                .push(RegionConformanceFailure::ReturnedError {
                    case: name,
                    error: error.to_string(),
                });
        }
        (Cut::Regions(expected), Ok(got)) => {
            if !plain
                .as_ref()
                .is_ok_and(|plain| cut_answer_matches(case, plain, expected, &got))
            {
                report.failures.push(wrong);
            }
        }
        (Cut::Refused(_) | Cut::Invalid, Ok(_)) => report.failures.push(wrong),
    }
}

fn same_error(left: &GeomError, right: &GeomError) -> bool {
    left.to_string() == right.to_string()
}

/// Whether a cut answer is the plain one plus the expected cut regions, with
/// matching evidence, and its rings cover the chains.
fn cut_answer_matches(
    case: &Case,
    plain: &PlaneSection,
    expected: &[CutRegion],
    got: &PlaneSection,
) -> bool {
    let Some(cut) = got.cut_regions.as_deref() else {
        return false;
    };
    let evidence = got.evidence;
    let derived = PlaneSectionEvidence::input_mesh(case.mesh.triangle_count())
        .with_source_counts(
            plain.evidence.in_plane_triangles,
            plain.evidence.touching_edges,
            plain.evidence.in_plane_pairs,
        )
        .with_role(case.role, case.boundary_edges)
        .with_output(&got.chains, &got.regions)
        .with_cut_pairs(evidence.cut_pairs)
        .with_cut_output(cut);
    got.frame == plain.frame
        && got.chains == plain.chains
        && got.regions == plain.regions
        && cut_regions_match(expected, cut)
        && evidence == derived
        && cut_rings_cover_the_chains(got, cut)
}

fn cut_regions_match(expected: &[CutRegion], got: &[CutRegion]) -> bool {
    expected.len() == got.len()
        && expected.iter().zip(got).all(|(e, g)| {
            points_match(&e.outer, &g.outer)
                && e.holes.len() == g.holes.len()
                && e.holes
                    .iter()
                    .zip(&g.holes)
                    .all(|(eh, gh)| points_match(eh, gh))
        })
}

/// The bit-for-bit relations between chains, regions and cut rings: every
/// closed chain is a ring; every open chain with both ends on a region is a
/// contiguous run of exactly one ring, forward or reversed; every ring
/// vertex is a chain point or a region vertex.
fn cut_rings_cover_the_chains(section: &PlaneSection, cut: &[CutRegion]) -> bool {
    let rings: Vec<&Vec<Point2>> = cut
        .iter()
        .flat_map(|region| std::iter::once(&region.outer).chain(&region.holes))
        .collect();
    let closed_chains_are_rings = section
        .chains
        .iter()
        .filter(|chain| chain.is_closed())
        .all(|chain| rings.iter().any(|ring| **ring == chain.points));
    let open_chains_are_runs = section
        .chains
        .iter()
        .filter(|chain| chain.ends == Some([ChainEnd::InPlaneRegion; 2]))
        .all(|chain| {
            let mut backwards = chain.points.clone();
            backwards.reverse();
            rings
                .iter()
                .filter(|ring| is_run(ring, &chain.points) || is_run(ring, &backwards))
                .count()
                == 1
        });
    let known = |point: &Point2| {
        section
            .chains
            .iter()
            .any(|chain| chain.points.contains(point))
            || section.regions.iter().any(|region| {
                std::iter::once(&region.outer)
                    .chain(&region.holes)
                    .any(|ring| ring.contains(point))
            })
    };
    closed_chains_are_rings
        && open_chains_are_runs
        && rings.iter().all(|ring| ring.iter().all(known))
}

/// Whether `run` is a contiguous cyclic run of `ring`, in ring order.
fn is_run(ring: &[Point2], run: &[Point2]) -> bool {
    run.len() <= ring.len()
        && (0..ring.len()).any(|start| {
            run.iter()
                .enumerate()
                .all(|(offset, point)| ring[(start + offset) % ring.len()] == *point)
        })
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
        .with_role(case.role, case.boundary_edges)
        .with_output(&got.chains, &got.regions);
    if evidence != derived || !evidence.is_derived_from_input_mesh() {
        report
            .failures
            .push(RegionConformanceFailure::WrongEvidence { case: name });
    }
}

/// Every end of an open chain declared to end on a region equals a region
/// vertex bit for bit: chains and regions are projected by one function.
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
        .filter_map(|chain| Some((chain, chain.ends?)))
        .all(|(chain, ends)| {
            let on_region_if_declared = |end: ChainEnd, point: Option<&Point2>| {
                end != ChainEnd::InPlaneRegion || point.is_some_and(on_region)
            };
            on_region_if_declared(ends[0], chain.points.first())
                && on_region_if_declared(ends[1], chain.points.last())
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
