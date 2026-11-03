//! The linear-scan broad phase and the unpruned arrangement give the regions
//! of the provider, on seeded box scenes in several views.
//!
//! A generic rotated scene can put two exact vertices of a ring within one
//! ulp of each other; the contract refuses that as `Degenerate`, and then
//! every variant must refuse it alike.

mod common;

use axiolid_hidden_lines::{ScalarVisibleFaces, ScanVisibleFaces, UnprunedVisibleFaces};
use axiolid_mesh_hidden_lines_contract::OrthographicView;
use common::{scene, solids, try_faces, views};

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
            let view = OrthographicView::new(frame);
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
