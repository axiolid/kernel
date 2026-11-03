//! The shared visible-face conformance suite, run against the exact scalar
//! provider and its test variants.

use axiolid_hidden_lines::{ScalarVisibleFaces, ScanVisibleFaces, UnprunedVisibleFaces};
use axiolid_mesh_hidden_lines_contract::conformance::faces::FaceConformanceSuite;

#[test]
fn scalar_visible_faces_passes_the_conformance_suite() {
    let report = FaceConformanceSuite::run(&ScalarVisibleFaces::new());
    assert!(report.is_success(), "{report}");
}

#[test]
fn the_scan_and_unpruned_variants_pass_the_conformance_suite() {
    for report in [
        FaceConformanceSuite::run(&ScanVisibleFaces),
        FaceConformanceSuite::run(&UnprunedVisibleFaces),
    ] {
        assert!(report.is_success(), "{report}");
    }
}
