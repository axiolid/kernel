//! Mutation probe: does the section conformance suite actually bite?
use axiolid_mesh_section_contract::conformance::{ConformanceSuite, RegionConformanceSuite};
use axiolid_mesh_section_contract::MeshPlaneSection;
use axiolid_reference::section::ScalarSection;

#[test]
fn scalar_section_satisfies_the_geometry_conformance_suite() {
    let report = ConformanceSuite::run(&ScalarSection);
    assert!(report.is_success(), "conformance failures: {report}");
}

#[test]
fn scalar_section_satisfies_the_region_conformance_suite() {
    let report = RegionConformanceSuite::run(&ScalarSection);
    assert!(
        report.claimed,
        "the reference provider implements section_request"
    );
    assert!(report.is_success(), "region conformance failures: {report}");
}

/// The reference provider, but refusing sheets as a solids-only provider would.
#[derive(Debug)]
struct SolidsOnly;

impl axiolid_contracts::Backend for SolidsOnly {
    fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
        ScalarSection.descriptor()
    }
}

impl MeshPlaneSection for SolidsOnly {
    fn section(
        &self,
        mesh: &axiolid_mesh::TriMesh,
        frame: axiolid_core::Frame3,
        limits: axiolid_mesh_section_contract::SectionLimits,
        options: &axiolid_contracts::ExecutionOptions,
    ) -> axiolid_contracts::GeomResult<axiolid_mesh_section_contract::SectionOutcome> {
        ScalarSection.section(mesh, frame, limits, options)
    }

    fn section_request(
        &self,
        mesh: &axiolid_mesh::TriMesh,
        request: &axiolid_mesh_section_contract::SectionRequest,
        options: &axiolid_contracts::ExecutionOptions,
    ) -> axiolid_contracts::GeomResult<axiolid_mesh_section_contract::PlaneSection> {
        if request.role == axiolid_mesh_section_contract::SectionRole::Sheet {
            return Err(axiolid_contracts::GeomError::Unsupported {
                backend: ScalarSection::ID,
                operation: axiolid_contracts::Operation::MeshPlaneSection,
            });
        }
        ScalarSection.section_request(mesh, request, options)
    }
}

#[test]
fn a_provider_without_sheets_passes_and_reports_it() {
    let report = RegionConformanceSuite::run(&SolidsOnly);
    assert!(report.claimed && report.is_success(), "{report}");
    assert!(!report.sheets_claimed);
    let full = RegionConformanceSuite::run(&ScalarSection);
    assert!(full.sheets_claimed);
}
