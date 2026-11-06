//! The two oracle scenes that binary64 rounding used to refuse (#410, ADR
//! 0100) are answered, and the evidence says what snap rounding did.
//!
//! The scenes are the seeded rotated boxes of `oracle.rs` in the oblique view
//! (seed 1) and the isometric view (seed 9), each with the near plane at the
//! middle of the depth range. Their evidence counts are pinned here; the
//! conformance suite pins frozen copies that do not depend on `sin_cos`.

mod common;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame3, Point3, Tolerance};
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    HiddenLineMesh, MeshId, MeshRole, MeshVisibleFaces, OrthographicView, VisibleFaceEvidence,
    VisibleFaceRequest,
};
use common::topology::assert_valid_topology;
use common::{scene, solids, try_faces, views, FACE_LIMITS};

/// The evidence and the region count of one scene; its topology is checked.
fn snapped(seed: u64, view: &str) -> (VisibleFaceEvidence, usize) {
    let (_, frame): (_, Frame3) = views()
        .into_iter()
        .find(|(name, _)| *name == view)
        .expect("the view");
    let boxes = scene(seed);
    let depth = |p: &axiolid_core::Point3| (*p - frame.origin).dot(frame.z);
    let depths = boxes.iter().flat_map(|m| m.positions.iter().map(depth));
    let (lo, hi) = depths.fold((f64::MAX, f64::MIN), |(l, h), w| (l.min(w), h.max(w)));
    let view = OrthographicView::new(frame).with_near((lo + hi) / 2.0);
    let meshes = solids(&boxes);
    let outcome = try_faces(&ScalarVisibleFaces::new(), &meshes, view, None)
        .unwrap_or_else(|error| panic!("seed {seed}, {view:?}: {error}"));
    assert_valid_topology(&outcome);
    (outcome.evidence, outcome.regions.len())
}

#[test]
fn seed_1_oblique_with_a_near_plane_is_snapped() {
    let (evidence, regions) = snapped(1, "oblique");
    // A sliver region of three points that round to one binary64 point.
    assert_eq!(evidence.snap_displacement.get(), 2f64.powi(-52));
    assert_eq!(
        (
            evidence.rerouted_edges,
            evidence.collapsed_regions,
            evidence.collapsed_holes,
            regions
        ),
        (0, 1, 0, 8)
    );
}

#[test]
fn seed_9_isometric_with_a_near_plane_is_snapped() {
    let (evidence, regions) = snapped(9, "isometric");
    // Two edges run through the rounding cell of another vertex, and a
    // sliver region is thinner than one binary64 step.
    assert_eq!(evidence.snap_displacement.get(), 2f64.powi(-53));
    assert_eq!(
        (
            evidence.rerouted_edges,
            evidence.collapsed_regions,
            evidence.collapsed_holes,
            regions
        ),
        (2, 1, 0, 23)
    );
}

/// A sheet a few binary64 steps wide in front of another, in an oblique view,
/// under a tolerance so small that the audit lets it be a triangle: its
/// corners are exact, non-binary64 image points that round to one cell. It
/// leaves a hole of no width in the sheet behind it and a region of no width
/// of its own. Both collapse and are counted; the sheet behind keeps its
/// outer ring.
#[test]
fn a_hole_and_its_occluder_thinner_than_a_cell_collapse() {
    let (_, frame) = views()
        .into_iter()
        .find(|(name, _)| *name == "oblique")
        .expect("the oblique view");
    let behind = TriMesh::new(
        vec![
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(-1.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let step = |x: f64| f64::from_bits(x.to_bits() + 1);
    let (x, y) = (0.3, 0.2);
    let tiny = TriMesh::new(
        vec![
            Point3::new(x, y, 1.0),
            Point3::new(step(x), y, 1.0),
            Point3::new(x, step(y), 1.0),
        ],
        vec![0, 1, 2],
    );
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &behind),
        HiddenLineMesh::new(MeshId(2), MeshRole::Sheet, &tiny),
    ];
    let request = VisibleFaceRequest::new(&meshes, OrthographicView::new(frame), FACE_LIMITS);
    let options = ExecutionOptions::new(Tolerance::new(1e-40, 1e-40).expect("tolerance"));
    let outcome = ScalarVisibleFaces::new()
        .visible_faces(&request, &options)
        .unwrap_or_else(|error| panic!("{error}"));
    assert_valid_topology(&outcome);
    let evidence = outcome.evidence;
    assert_eq!(
        (evidence.collapsed_regions, evidence.collapsed_holes),
        (1, 1)
    );
    assert!(evidence.snap_displacement.get() > 0.0);
    assert_eq!(outcome.regions.len(), 1);
    assert_eq!(outcome.regions[0].mesh, MeshId(1));
    assert!(outcome.regions[0].holes.is_empty());
}
