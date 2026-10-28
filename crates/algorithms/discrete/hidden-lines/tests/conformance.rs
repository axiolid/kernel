//! The shared conformance suite, run against the exact scalar provider.

use axiolid_hidden_lines::ScalarHiddenLines;
use axiolid_mesh_hidden_lines_contract::conformance::ConformanceSuite;

#[test]
fn scalar_hidden_lines_passes_the_conformance_suite() {
    let report = ConformanceSuite::run(&ScalarHiddenLines::new());
    assert!(report.is_success(), "{report}");
}
