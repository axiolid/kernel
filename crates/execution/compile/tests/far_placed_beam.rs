//! A designed tangent contact placed far from the origin (#380).
//!
//! The beam of #243, whose round web hole touches the top flange's inner
//! face, under a placement that turns about `z` and moves it. The rounding
//! of the placement makes the tangent contact read as a graze of a few
//! thousandths of a radian that is no section: it keeps its chord, as it
//! does without #380, so the mesh cut refuses it no more often than it did.
//! The cut's own refusals ("meet tangentially") at the rounding of the
//! coordinates were counted per scale before #380: one of twelve turns at
//! 12.5 m, none at a few hundred metres, one at a kilometre, five at UTM.

use std::f64::consts::FRAC_PI_2;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance, Transform3, Vec3};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{CircleProfile, Profile, SectionProfile};

fn graph(building: Transform3) -> (GeometryGraph, NodeId) {
    let mut g = GeometryGraphBuilder::new();
    let ext = |g: &mut GeometryGraphBuilder, profile: Profile, depth: f64| {
        let profile = g.push(GeometryNode::Profile(profile)).unwrap();
        g.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
        .unwrap()
    };
    let (depth, web, flange, radius) = (0.3, 0.0071, 0.0107, 0.05);
    let section = Profile::Section(SectionProfile::I {
        depth,
        width: 0.15,
        web_thickness: web,
        flange_thickness: flange,
        fillet_radius: None,
        flange_edge_radius: None,
        flange_slope: None,
    });
    let beam = ext(&mut g, section, 4.0);
    let beam = g
        .push(GeometryNode::Instance(Instance {
            source: beam,
            transform: building,
        }))
        .unwrap();
    let hole = ext(
        &mut g,
        Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        }),
        web,
    );
    let across =
        Transform3::from_translation(Vec3::new(-web / 2.0, depth / 2.0 - flange - radius, 2.0))
            * Transform3::from_rotation_y(FRAC_PI_2);
    let hole = g
        .push(GeometryNode::Instance(Instance {
            source: hole,
            transform: building * across,
        }))
        .unwrap();
    let root = g
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: beam,
            right: hole,
            operator: BooleanOperator::Difference,
        }))
        .unwrap();
    (g.finish(vec![root]).unwrap(), root)
}

/// The meshes refused of twelve turns of the beam about `z` at `at`.
fn refused(at: Vec3) -> usize {
    let options = ExecutionOptions::new(Tolerance::new(1e-3, 1e-9).unwrap())
        .with_chord_error(1e-3)
        .unwrap();
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    (0..12)
        .filter(|&k| {
            let turn = Transform3::from_rotation_z(0.6 + 0.05 * k as f64);
            let (graph, root) = graph(Transform3::from_translation(at) * turn);
            compiler.compile_mesh(&graph, root, &options).is_err()
        })
        .count()
}

#[test]
fn a_tangent_hole_placed_far_away_is_refused_no_more_often_than_before() {
    for (at, before) in [
        (Vec3::new(12.5, -4.0, 3.2), 1),
        (Vec3::new(250.0, -400.0, 3.2), 0),
        (Vec3::new(1000.0, 2000.0, 10.0), 1),
        (Vec3::new(3.6e5, 5.8e6, 120.0), 5),
    ] {
        let now = refused(at);
        assert!(
            now <= before,
            "{at:?}: {now} of 12 refused, {before} before"
        );
    }
}
