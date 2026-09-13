use std::f64::consts::TAU;

use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_contracts::{GeomError, Operation};
use axiolid_core::{Tolerance, Vec2, Vec3};
use axiolid_curve::Curve2;
use axiolid_evaluate::{evaluate2, evaluate3};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_surface::Surface;
use axiolid_topology::{audit_brep, Orientation};

fn rectangle(thickness: Option<f64>) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 4.0,
        y: 2.0,
        thickness,
        outer_radius: None,
        inner_radius: None,
    })
}

fn assert_complete_and_closed(value: &axiolid_construct::ExactBRep) {
    let topology = value.topology();
    assert!(audit_brep(topology).is_closed_manifold());
    assert_eq!(topology.shells().len(), 1);
    assert!(topology.shells()[0].closed);
    assert_eq!(topology.solids().len(), 1);

    for index in 0..topology.edges().len() {
        let edge = topology.edge_id_at(index).expect("dense edge");
        assert!(value.edge_interval(edge).is_some(), "edge {index} interval");
    }
    for (loop_index, wire) in topology.loops().iter().enumerate() {
        let loop_id = topology.loop_id_at(loop_index).expect("dense loop");
        for use_index in 0..wire.edges.len() {
            assert!(
                value.pcurve_interval(loop_id, use_index).is_some(),
                "loop {loop_index} use {use_index} pcurve interval"
            );
        }
    }
}

fn assert_pcurve_surface_agreement(brep: &axiolid_construct::ExactBRep) {
    let topology = brep.topology();
    for face in topology.faces() {
        let surface_id = face.surface.expect("exact face support");
        let surface = &brep.surfaces()[surface_id.index()];
        for bound in &face.bounds {
            let wire = &topology.loops()[bound.loop_id.index()];
            for (use_index, edge_use) in wire.edges.iter().enumerate() {
                let pcurve_id = edge_use.pcurve.expect("exact pcurve");
                let pcurve = &brep.curves2()[pcurve_id.index()];
                let pspan = brep
                    .pcurve_interval(bound.loop_id, use_index)
                    .expect("pcurve span");
                let edge = &topology.edges()[edge_use.edge.index()];
                let curve_id = edge.curve.expect("exact edge support");
                let curve = &brep.curves3()[curve_id.index()];
                let espan = brep.edge_interval(edge_use.edge).expect("edge span");

                for alpha in [0.0, 0.5, 1.0] {
                    let p_t = pspan.start + (pspan.end - pspan.start) * alpha;
                    let e_t = match edge_use.orientation {
                        Orientation::Forward => espan.start + (espan.end - espan.start) * alpha,
                        Orientation::Reversed => espan.end + (espan.start - espan.end) * alpha,
                    };
                    let uv = evaluate2(pcurve, p_t).expect("pcurve evaluation");
                    let from_surface = axiolid_evaluate::surface::evaluate(surface, uv.x, uv.y)
                        .expect("surface evaluation");
                    let from_curve = evaluate3(curve, e_t).expect("edge evaluation");
                    assert!(
                        from_surface.distance(from_curve) <= 1.0e-10,
                        "pcurve/surface mismatch: {from_surface:?} vs {from_curve:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn filled_rectangle_extrusion_is_an_exact_planar_solid() {
    let value = extrude_profile_exact(
        &rectangle(None),
        Vec3::new(0.25, 0.5, 1.0),
        3.0,
        Tolerance::METRE,
    )
    .expect("exact rectangle extrusion");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
    assert_eq!(value.topology().vertices().len(), 8);
    assert_eq!(value.topology().edges().len(), 12);
    assert_eq!(value.topology().faces().len(), 6);
    assert_eq!(value.surfaces().len(), 6);
    assert!(value
        .surfaces()
        .iter()
        .all(|surface| matches!(surface, Surface::Plane(_))));
}

#[test]
fn hollow_rectangle_uses_cap_inner_bounds_and_through_passage_faces() {
    let value = extrude_profile_exact(&rectangle(Some(0.25)), Vec3::Z, 2.0, Tolerance::METRE)
        .expect("exact hollow rectangle extrusion");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
    assert_eq!(value.topology().vertices().len(), 16);
    assert_eq!(value.topology().edges().len(), 24);
    assert_eq!(value.topology().faces().len(), 10);
    assert_eq!(value.topology().faces()[0].bounds.len(), 2);
    assert_eq!(value.topology().faces()[1].bounds.len(), 2);
    assert_eq!(value.topology().solids()[0].voids.len(), 0);
}

#[test]
fn axial_circle_extrusion_preserves_cylinder_and_distinct_seam_charts() {
    let value = extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius: 2.0,
            thickness: None,
        }),
        Vec3::Z,
        5.0,
        Tolerance::METRE,
    )
    .expect("exact cylinder");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
    assert_eq!(value.topology().vertices().len(), 2);
    assert_eq!(value.topology().edges().len(), 3);
    assert_eq!(value.topology().faces().len(), 3);
    assert_eq!(
        value
            .surfaces()
            .iter()
            .filter(|surface| matches!(surface, Surface::Cylinder(_)))
            .count(),
        1
    );
    assert_eq!(
        value
            .surfaces()
            .iter()
            .filter(|surface| matches!(surface, Surface::Plane(_)))
            .count(),
        2
    );

    let side = &value.topology().loops()[2];
    let seam_u = side
        .edges
        .iter()
        .filter_map(|edge_use| edge_use.pcurve)
        .filter_map(|id| match &value.curves2()[id.index()] {
            Curve2::Line(line) if line.direction.x == 0.0 => Some(line.origin.x),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(seam_u.len(), 2);
    assert!(seam_u.contains(&0.0));
    assert!(seam_u.contains(&TAU));
}

#[test]
fn exact_extrusion_refuses_families_whose_supports_are_not_populated() {
    let rounded = Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 1.0,
        thickness: None,
        outer_radius: Some(0.1),
        inner_radius: None,
    });
    let annulus = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: Some(0.2),
    });
    for (profile, family) in [
        (rounded, "rounded rectangle extrusion"),
        (annulus, "annular circle extrusion"),
    ] {
        let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
            .expect_err("unsupported exact family must refuse");
        assert!(matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::Sweep,
                input,
                ..
            } if input == family
        ));
    }
}

#[test]
fn oblique_circle_extrusion_refuses_instead_of_mislabeling_a_cylinder() {
    let profile = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let error = extrude_profile_exact(&profile, Vec3::new(1.0, 0.0, 1.0), 1.0, Tolerance::METRE)
        .expect_err("oblique circle needs a non-cylindrical swept support");

    assert!(matches!(
        error,
        GeomError::UnsupportedInput {
            operation: Operation::Sweep,
            input: "oblique circle extrusion",
            ..
        }
    ));
}

#[test]
fn sub_tolerance_oblique_circle_refuses_instead_of_snapping_geometry() {
    let profile = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });
    let direction = Vec3::new(Tolerance::METRE.linear() * 0.5, 0.0, 1.0);
    let error = extrude_profile_exact(&profile, direction, 1.0, Tolerance::METRE)
        .expect_err("an exact constructor must not tolerance-snap an oblique direction");

    assert!(matches!(
        error,
        GeomError::UnsupportedInput {
            operation: Operation::Sweep,
            input: "oblique circle extrusion",
            ..
        }
    ));
}

#[test]
fn nonzero_direction_magnitude_is_not_a_model_length() {
    let direction = Vec3::Z * (Tolerance::METRE.linear() * 0.5);
    let value = extrude_profile_exact(
        &Profile::Circle(CircleProfile {
            radius: 1.0,
            thickness: None,
        }),
        direction,
        1.0,
        Tolerance::METRE,
    )
    .expect("direction magnitude is unitless and normalization preserves the axis");

    assert_complete_and_closed(&value);
    assert_pcurve_surface_agreement(&value);
}

#[test]
fn finite_nonzero_extreme_direction_magnitudes_normalize_safely() {
    let profile = Profile::Circle(CircleProfile {
        radius: 1.0,
        thickness: None,
    });

    for magnitude in [f64::MIN_POSITIVE, f64::MAX] {
        let value = extrude_profile_exact(
            &profile,
            Vec3::new(0.0, 0.0, magnitude),
            1.0,
            Tolerance::METRE,
        )
        .expect("finite nonzero unitless directions must normalize without range loss");

        assert_complete_and_closed(&value);
        assert_pcurve_surface_agreement(&value);
    }
}

#[test]
fn mixed_extreme_direction_refuses_nonzero_component_loss() {
    let error = extrude_profile_exact(
        &rectangle(None),
        Vec3::new(f64::MIN_POSITIVE, 0.0, f64::MAX),
        1.0,
        Tolerance::METRE,
    )
    .expect_err("normalization must not silently erase a nonzero direction component");

    assert!(matches!(error, GeomError::Degenerate(_)));
}

#[test]
fn invalid_exact_extrusion_inputs_are_not_reclassified_as_unsupported() {
    let error = extrude_profile_exact(
        &rectangle(None),
        Vec2::ZERO.extend(0.0),
        1.0,
        Tolerance::METRE,
    )
    .expect_err("zero direction");
    assert!(matches!(error, GeomError::InvalidInput(_)));
}
