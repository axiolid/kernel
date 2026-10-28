//! Mutation probe: does the section conformance suite actually bite?
use axiolid_mesh_section_contract::conformance::{ConformanceSuite, RegionConformanceSuite};
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
