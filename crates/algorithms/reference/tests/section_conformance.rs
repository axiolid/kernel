//! Mutation probe: does the section conformance suite actually bite?
use axiolid_mesh_section_contract::conformance::ConformanceSuite;
use axiolid_reference::section::ScalarSection;

#[test]
fn scalar_section_satisfies_the_geometry_conformance_suite() {
    let report = ConformanceSuite::run(&ScalarSection);
    assert!(report.is_success(), "conformance failures: {report}");
}
