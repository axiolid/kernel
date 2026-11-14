//! See-through triangles (ADR 0111). Every provider claims them and passes
//! the see-through suite; flags that are all opaque change nothing, evidence
//! and exact work included. A provider that claims them and lets them hide
//! is refused at registration.

mod common;

use axiolid_contracts::{Backend, BackendDescriptor, ExecutionOptions, GeomError};
use axiolid_core::{Point3, Tolerance};
use axiolid_dispatch::{MeshHiddenLinesRegistry, MeshVisibleFacesRegistry};
use axiolid_hidden_lines::{
    ScalarHiddenLines, ScalarVisibleFaces, ScanHiddenLines, ScanVisibleFaces, UnprunedVisibleFaces,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::conformance::{
    faces::FaceConformanceFailure, see_through::SeeThroughConformanceSuite, ConformanceFailure,
};
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineError, HiddenLineMesh, HiddenLineOutcome, HiddenLineRequest,
    MeshHiddenLines, MeshId, MeshRole, MeshVisibleFaces, Occlusion, TriangleOcclusion, View,
    VisibleFaceOutcome, VisibleFaceRequest,
};
use common::{
    box_mesh, plan_frame, scene, solids, try_draw_with, try_faces, views, FACE_LIMITS, LIMITS,
};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

/// A provider that claims see-through triangles and lets them hide: it draws
/// the request with every flag dropped.
#[derive(Debug)]
struct Ignores;

fn opaque<'a>(meshes: &[HiddenLineMesh<'a>]) -> Vec<HiddenLineMesh<'a>> {
    meshes
        .iter()
        .map(|entry| {
            let mut entry = *entry;
            entry.occlusion = TriangleOcclusion::AllOpaque;
            entry
        })
        .collect()
}

impl Backend for Ignores {
    fn descriptor(&self) -> BackendDescriptor {
        ScalarHiddenLines.descriptor()
    }
}

impl MeshHiddenLines for Ignores {
    fn supports_occlusion(&self, _occlusion: Occlusion) -> bool {
        true
    }

    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        let meshes = opaque(request.meshes);
        let mut request = request.clone();
        request.meshes = &meshes;
        ScalarHiddenLines.hidden_lines(&request, options)
    }
}

impl MeshVisibleFaces for Ignores {
    fn supports_occlusion(&self, _occlusion: Occlusion) -> bool {
        true
    }

    fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError> {
        let meshes = opaque(request.meshes);
        let mut request = request.clone();
        request.meshes = &meshes;
        ScalarVisibleFaces.visible_faces(&request, options)
    }
}

#[test]
fn register_conformant_refuses_an_unverified_see_through_claim() {
    let mut lines = MeshHiddenLinesRegistry::new();
    let report = lines
        .register_conformant(0, Ignores)
        .expect_err("the claim is false");
    assert!(report.failures.iter().any(|f| matches!(
        f,
        ConformanceFailure::ReturnedError {
            case: "see-through claim",
            ..
        }
    )));
    assert_eq!(lines.providers().count(), 0);

    let mut faces = MeshVisibleFacesRegistry::new();
    let report = faces
        .register_conformant(0, Ignores)
        .expect_err("the claim is not verified");
    assert!(report.failures.iter().any(|f| matches!(
        f,
        FaceConformanceFailure::ReturnedError {
            case: "see-through claim",
            ..
        }
    )));
    assert_eq!(faces.providers().count(), 0);
}

#[test]
fn the_hidden_line_providers_pass_the_see_through_suite() {
    let lines: [&dyn MeshHiddenLines; 2] = [&ScalarHiddenLines, &ScanHiddenLines];
    for provider in lines {
        let report = SeeThroughConformanceSuite::run(provider);
        assert!(report.claimed);
        assert!(report.is_success(), "{report}");
    }
    let mut registry = MeshHiddenLinesRegistry::new();
    registry
        .register_conformant(0, ScalarHiddenLines::new())
        .unwrap_or_else(|report| panic!("{report}"));
}

/// The seeded scenes of the oracle in every view, with and without a near
/// plane through the middle of one box: flags that are all opaque give the
/// same outcome, evidence and escalations included, as no flags.
#[test]
fn opaque_flags_change_no_hidden_line_drawing() {
    for seed in 1..=40 {
        let meshes = scene(seed);
        let plain = solids(&meshes);
        let opaque: Vec<Vec<Occlusion>> = meshes
            .iter()
            .map(|m| vec![Occlusion::Opaque; m.triangle_count()])
            .collect();
        let flagged: Vec<HiddenLineMesh<'_>> = plain
            .iter()
            .zip(&opaque)
            .map(|(entry, flags)| entry.with_occlusion(flags))
            .collect();
        for (name, frame) in views() {
            for near in [None, Some(0.5)] {
                let view = near.map_or(View::new(frame), |near| View::new(frame).with_near(near));
                let edges = EdgeSelection::new(Some(0.5));
                let a = try_draw_with(&ScalarHiddenLines, &plain, view, edges);
                let b = try_draw_with(&ScalarHiddenLines, &flagged, view, edges);
                assert_eq!(a, b, "seed {seed}, {name}, near {near:?}");
                if let Ok(outcome) = &b {
                    assert_eq!(outcome.evidence.see_through_triangles, 0);
                    assert_eq!(outcome.evidence.seen_through_segments, 0);
                    assert!(outcome.segments.iter().all(|s| s.seen_through == 0));
                }
            }
        }
    }
}

/// A piece that collapses under rounding joins its neighbours only when they
/// agree in visibility and count.
#[test]
fn a_collapsed_piece_joins_neighbours_of_one_count_only() {
    // The floor's edge y = 4 from x = 4 to 0 is under pane B from x = 4 to
    // 3 + u/4 (u the ulp of 3: B's left edge runs from (3, 0) to (3 + u, 16)),
    // under nothing to x = 3, a part whose ends both round to 3, and under
    // two stacked panes from 3 to 1.
    let ulp = f64::from_bits(3.0_f64.to_bits() + 1) - 3.0;
    let floor = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(4.0, 4.0, 0.0),
            Point3::new(0.0, 4.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let b = TriMesh::new(
        vec![
            Point3::new(3.0, 0.0, 2.0),
            Point3::new(6.0, 8.0, 2.0),
            Point3::new(3.0 + ulp, 16.0, 2.0),
        ],
        vec![0, 1, 2],
    );
    let pane = |z: f64| {
        TriMesh::new(
            vec![
                Point3::new(1.0, 1.0, z),
                Point3::new(3.0, 1.0, z),
                Point3::new(3.0, 5.0, z),
                Point3::new(1.0, 5.0, z),
            ],
            vec![0, 1, 2, 0, 2, 3],
        )
    };
    let (low, high) = (pane(1.0), pane(1.5));
    let (one, two) = ([Occlusion::SeeThrough], [Occlusion::SeeThrough; 2]);
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &floor),
        HiddenLineMesh::new(MeshId(2), MeshRole::Sheet, &b).with_occlusion(&one),
        HiddenLineMesh::new(MeshId(3), MeshRole::Sheet, &low).with_occlusion(&two),
        HiddenLineMesh::new(MeshId(4), MeshRole::Sheet, &high).with_occlusion(&two),
    ];
    let request = HiddenLineRequest::new(
        &meshes,
        View::new(plan_frame()),
        EdgeSelection::new(None),
        LIMITS,
    );
    let outcome = ScalarHiddenLines
        .hidden_lines(&request, &options())
        .expect("draws");
    let pieces: Vec<_> = outcome
        .segments
        .iter()
        .filter(|s| s.source.mesh == MeshId(1) && s.source.vertices == [2, 3])
        .map(|s| (s.seen_through, s.start.x, s.end.x))
        .collect();
    assert_eq!(pieces, [(1, 4.0, 3.0), (2, 3.0, 1.0), (0, 1.0, 0.0)]);
    assert!(outcome.evidence.collapsed_segments > 0);
}

#[test]
fn a_short_occlusion_slice_is_refused_by_the_provider() {
    let mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let flags = [Occlusion::Opaque; 11];
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &mesh),
        HiddenLineMesh::new(MeshId(2), MeshRole::Solid, &mesh).with_occlusion(&flags),
    ];
    let request = HiddenLineRequest::new(
        &meshes,
        View::new(plan_frame()),
        EdgeSelection::new(None),
        LIMITS,
    );
    let error = ScalarHiddenLines
        .hidden_lines(&request, &options())
        .expect_err("refused");
    assert_eq!(error.mesh, Some(MeshId(2)));
    assert!(
        matches!(error.error, GeomError::InvalidInput(_)),
        "{error:?}"
    );
}

#[test]
fn the_visible_face_providers_pass_the_see_through_suite() {
    let faces: [&dyn MeshVisibleFaces; 3] = [
        &ScalarVisibleFaces,
        &ScanVisibleFaces,
        &UnprunedVisibleFaces,
    ];
    for provider in faces {
        let report = SeeThroughConformanceSuite::run_faces(provider);
        assert!(report.claimed);
        assert!(report.is_success(), "{report}");
    }
    let mut registry = MeshVisibleFacesRegistry::new();
    registry
        .register_conformant(0, ScalarVisibleFaces::new())
        .unwrap_or_else(|report| panic!("{report}"));
}

/// The seeded scenes of the face oracle in every view, with and without a
/// near plane: flags that are all opaque give the same outcome, evidence and
/// escalations included, as no flags.
#[test]
fn opaque_flags_change_no_visible_faces() {
    for seed in 0..12 {
        let meshes = scene(seed);
        let plain = solids(&meshes);
        let flags: Vec<Vec<Occlusion>> = meshes
            .iter()
            .map(|m| vec![Occlusion::Opaque; m.triangle_count()])
            .collect();
        let flagged: Vec<HiddenLineMesh<'_>> = plain
            .iter()
            .zip(&flags)
            .map(|(entry, flags)| entry.with_occlusion(flags))
            .collect();
        for (name, frame) in views() {
            for near in [None, Some(0.5)] {
                let view = near.map_or(View::new(frame), |near| View::new(frame).with_near(near));
                let a = try_faces(&ScalarVisibleFaces, &plain, view, None);
                let b = try_faces(&ScalarVisibleFaces, &flagged, view, None);
                assert_eq!(a, b, "seed {seed}, {name}, near {near:?}");
                if let Ok(outcome) = &b {
                    assert_eq!(outcome.evidence.see_through_triangles, 0);
                    assert_eq!(outcome.evidence.behind_surfaces, 0);
                    assert!(outcome.regions.iter().all(|r| r.behind.is_empty()));
                }
            }
        }
    }
}

#[test]
fn a_short_occlusion_slice_or_a_mixed_group_is_refused_by_the_face_provider() {
    let mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let short = [Occlusion::SeeThrough; 11];
    let mut mixed = [Occlusion::Opaque; 12];
    mixed[3] = Occlusion::SeeThrough;
    for flags in [&short[..], &mixed[..]] {
        let meshes = [
            HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &mesh),
            HiddenLineMesh::new(MeshId(2), MeshRole::Solid, &mesh).with_occlusion(flags),
        ];
        let request = VisibleFaceRequest::new(
            &meshes,
            View::new(plan_frame()),
            FACE_LIMITS.with_max_behind_surfaces(16),
        );
        let error = ScalarVisibleFaces
            .visible_faces(&request, &options())
            .expect_err("refused");
        assert_eq!(error.mesh, Some(MeshId(2)));
        assert!(
            matches!(error.error, GeomError::InvalidInput(_)),
            "{error:?}"
        );
    }
}
