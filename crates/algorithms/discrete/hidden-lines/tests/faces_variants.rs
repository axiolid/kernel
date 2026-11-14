//! The linear-scan broad phase and the unpruned arrangement give the regions
//! of the provider, on seeded box scenes in several views.
//!
//! A generic rotated scene can put two exact vertices of a ring within one
//! ulp of each other; the contract refuses that as `Degenerate`, and then
//! every variant must refuse it alike.

mod common;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::Tolerance;
use axiolid_hidden_lines::{ScalarVisibleFaces, ScanVisibleFaces, UnprunedVisibleFaces};
use axiolid_mesh_hidden_lines_contract::{
    HiddenLineMesh, MeshVisibleFaces, Occlusion, View, VisibleFaceRequest,
};
use common::{scene, solids, try_faces, views, FACE_LIMITS};

#[test]
fn the_variants_give_equal_outcomes() {
    let (mut answered, mut refused) = (0, 0);
    for seed in 0..20 {
        let boxes = scene(seed);
        let meshes = solids(&boxes);
        // Plan, an elevation and the oblique view.
        for (name, frame) in views()
            .into_iter()
            .filter(|(n, _)| ["plan", "-Y", "oblique"].contains(n))
        {
            let view = View::new(frame);
            let base = try_faces(&ScalarVisibleFaces::new(), &meshes, view, None);
            for (variant, other) in [
                ("scan", try_faces(&ScanVisibleFaces, &meshes, view, None)),
                (
                    "unpruned",
                    try_faces(&UnprunedVisibleFaces, &meshes, view, None),
                ),
            ] {
                let context = format!("seed {seed}, {name}, {variant}");
                let (Ok(a), Ok(b)) = (&base, &other) else {
                    assert_eq!(base.as_ref().err(), other.as_ref().err(), "{context}");
                    continue;
                };
                assert_eq!(a.regions, b.regions, "{context}");
                // Evidence that legitimately differs: how the arrangement was
                // reached, not what it says. `tied_sides` counts sides of
                // pieces, so only whether there is a tie must agree.
                let (a, b) = (a.evidence, b.evidence);
                assert_eq!(
                    (
                        a.source_triangles,
                        a.edge_on_triangles,
                        a.removed_triangles,
                        a.tied_sides > 0,
                        a.regions,
                        a.cap_regions
                    ),
                    (
                        b.source_triangles,
                        b.edge_on_triangles,
                        b.removed_triangles,
                        b.tied_sides > 0,
                        b.regions,
                        b.cap_regions
                    ),
                    "{context}"
                );
            }
            if base.is_ok() {
                answered += 1;
            } else {
                refused += 1;
            }
        }
    }
    assert!(answered >= 50, "answered {answered}, refused {refused}");
}

/// The same with every third box wholly see-through (ADR 0111), with and
/// without a near plane: the prune must keep every edge across which a stack
/// changes, so the regions and their `behind` lists agree.
#[test]
fn the_variants_agree_on_see_through_stacks() {
    let (mut answered, mut stacked) = (0, 0);
    for seed in 0..20 {
        let boxes = scene(seed);
        let flags: Vec<Vec<Occlusion>> = boxes
            .iter()
            .enumerate()
            .map(|(k, mesh)| {
                let flag = if k % 3 == 0 {
                    Occlusion::SeeThrough
                } else {
                    Occlusion::Opaque
                };
                vec![flag; mesh.triangle_count()]
            })
            .collect();
        let meshes: Vec<HiddenLineMesh<'_>> = solids(&boxes)
            .into_iter()
            .zip(&flags)
            .map(|(entry, flags)| entry.with_occlusion(flags))
            .collect();
        for (name, frame) in views()
            .into_iter()
            .filter(|(n, _)| ["plan", "-Y", "oblique"].contains(n))
        {
            for near in [None, Some(0.5)] {
                let view = near.map_or(View::new(frame), |near| View::new(frame).with_near(near));
                let ask = |provider: &dyn MeshVisibleFaces| {
                    let request = VisibleFaceRequest::new(
                        &meshes,
                        view,
                        FACE_LIMITS.with_max_behind_surfaces(1_000_000),
                    );
                    provider.visible_faces(&request, &options())
                };
                let base = ask(&ScalarVisibleFaces::new());
                for (variant, other) in [
                    ("scan", ask(&ScanVisibleFaces)),
                    ("unpruned", ask(&UnprunedVisibleFaces)),
                ] {
                    let context = format!("seed {seed}, {name}, near {near:?}, {variant}");
                    let (Ok(a), Ok(b)) = (&base, &other) else {
                        assert_eq!(base.as_ref().err(), other.as_ref().err(), "{context}");
                        continue;
                    };
                    assert_eq!(a.regions, b.regions, "{context}");
                    assert_eq!(
                        (
                            a.evidence.behind_surfaces,
                            a.evidence.see_through_triangles,
                            a.evidence.tied_sides > 0
                        ),
                        (
                            b.evidence.behind_surfaces,
                            b.evidence.see_through_triangles,
                            b.evidence.tied_sides > 0
                        ),
                        "{context}"
                    );
                }
                if let Ok(outcome) = &base {
                    answered += 1;
                    stacked += usize::from(outcome.evidence.behind_surfaces > 0);
                }
            }
        }
    }
    assert!(answered >= 100, "answered {answered}");
    assert!(stacked >= 50, "stacked {stacked} of {answered}");
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}
