//! Refinement reports each channel's fate truthfully.
//!
//! A report that says `Preserved` must be backed by the channel actually
//! being on the returned mesh; one that says `Dropped` must mean it is not.

use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_mesh::{AttributeChannel, AttributeFate, Blend, DropReason, NormalAttribute, TriMesh};
use axiolid_refine::{refine, RefineTarget};

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

fn triangle_with(blend: Blend) -> TriMesh {
    let p = |x: f64, y: f64| Point3::new(x, y, 0.0);
    let mut mesh = TriMesh::new(vec![p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0)], vec![0, 1, 2]);
    mesh.attributes
        .push(AttributeChannel::new("c", vec![1.0, 2.0, 3.0], 1, blend));
    mesh.normals = Some(NormalAttribute {
        values: vec![Vec3::Z; 3],
        indices: None,
    });
    mesh
}

#[test]
fn zero_levels_returns_every_channel_it_reports_preserved() {
    for blend in [Blend::Linear, Blend::Nearest, Blend::None] {
        let mesh = triangle_with(blend);
        let (out, report) =
            refine(&mesh, RefineTarget::Uniform { levels: 0 }, None, tol()).expect("refine");
        assert_eq!(report.vertices_added, 0);
        assert_eq!(
            report.attribute_fates,
            vec![("c".to_owned(), AttributeFate::Preserved)]
        );
        assert_eq!(
            out.attributes, mesh.attributes,
            "{blend:?}: reported Preserved, must be present"
        );
        assert_eq!(
            out.normals, mesh.normals,
            "{blend:?}: normals are exact too"
        );
        out.validate_structure().expect("valid");
    }
}

#[test]
fn a_real_refinement_drops_what_it_reports_dropped() {
    for (blend, reason) in [
        (Blend::Linear, DropReason::ProviderLimitation),
        (Blend::None, DropReason::NotBlendable),
    ] {
        let (out, report) = refine(
            &triangle_with(blend),
            RefineTarget::Uniform { levels: 1 },
            None,
            tol(),
        )
        .expect("refine");
        assert_eq!(
            report.attribute_fates,
            vec![("c".to_owned(), AttributeFate::Dropped(reason))]
        );
        assert!(out.attributes.is_empty());
        out.validate_structure().expect("valid");
    }
}
