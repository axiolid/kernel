//! An exact placed difference on a rolled-section beam at `Tolerance::ZERO`
//! (#250).
//!
//! IPE sizes are decimal, so the I profile's contour segments meet only to
//! the rounding of their own evaluation. Lowering refused that at
//! `Tolerance::ZERO` ("contour segments leave a gap of 2.6e-18") before any
//! boolean ran. Now the beam lowers, and a round web hole clear of the root
//! fillets -- it crosses only the web's two planar faces, transversally --
//! cuts it exactly: an empty report and the closed-form volume.

use std::f64::consts::PI;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Mat3, Tolerance, Transform3, Vec3};
use axiolid_measure::exact_properties;
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{CircleProfile, Profile, SectionProfile};

/// IPE 300: h, b, tw, tf, root radius r.
const H: f64 = 0.3;
const B: f64 = 0.15;
const TW: f64 = 0.0071;
const TF: f64 = 0.0107;
const R: f64 = 0.015;
const LEN: f64 = 2.0;
/// The web hole: radius 50 mm on the section's centre line, far below the
/// fillets, which start `H / 2 - TF - R = 0.1243` from it.
const HOLE: f64 = 0.05;

fn ipe(fillet: Option<f64>) -> Profile {
    Profile::Section(SectionProfile::I {
        depth: H,
        width: B,
        web_thickness: TW,
        flange_thickness: TF,
        fillet_radius: fillet,
        flange_edge_radius: None,
        flange_slope: None,
    })
}

fn area(fillet: Option<f64>) -> f64 {
    let fillets = fillet.map_or(0.0, |r| 4.0 * r * r * (1.0 - PI / 4.0));
    2.0 * B * TF + (H - 2.0 * TF) * TW + fillets
}

/// The beam `[0, LEN]` along `z` under `placement`, minus a round hole along
/// `x` through the web at mid-length, from `start` for `depth`.
fn beam_with_hole(
    fillet: Option<f64>,
    (start, depth): (f64, f64),
    placement: Transform3,
) -> (GeometryGraph, NodeId) {
    let mut g = GeometryGraphBuilder::new();
    let mut push = |node| g.push(node).expect("a valid node");
    let section = push(GeometryNode::Profile(ipe(fillet)));
    let beam = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile: section,
        direction: Vec3::Z,
        depth: LEN,
    }));
    let beam = push(GeometryNode::Instance(Instance {
        source: beam,
        transform: placement,
    }));
    let disk = push(GeometryNode::Profile(Profile::Circle(CircleProfile {
        radius: HOLE,
        thickness: None,
    })));
    let rod = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile: disk,
        direction: Vec3::Z,
        depth,
    }));
    // The rod's `+z` turned onto `+x` with entries 0 and +-1.
    let onto_x = Transform3::from_mat3(Mat3::from_cols(-Vec3::Z, Vec3::Y, Vec3::X));
    let rod = push(GeometryNode::Instance(Instance {
        source: rod,
        transform: placement
            * Transform3::from_translation(Vec3::new(start, 0.0, LEN / 2.0))
            * onto_x,
    }));
    let body = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: beam,
        right: rod,
        operator: BooleanOperator::Difference,
    }));
    (g.finish(vec![body]).expect("a valid graph"), body)
}

/// Placements with exact entries: none, and a quarter or half turn with a
/// dyadic offset.
fn exact_placements() -> [Transform3; 3] {
    let quarter = Mat3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z);
    let half = Mat3::from_cols(-Vec3::X, -Vec3::Y, Vec3::Z);
    [
        Transform3::IDENTITY,
        Transform3::from_mat3_translation(quarter, Vec3::new(12.5, -4.0, 3.25)),
        Transform3::from_mat3_translation(half, Vec3::new(-7.75, 2.5, 0.0)),
    ]
}

fn audited(brep: &ExactBRep) {
    let health = geometric_audit(brep, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(brep.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
}

#[test]
fn an_ipe_beam_extrudes_at_zero_tolerance() {
    for fillet in [None, Some(R)] {
        let mut g = GeometryGraphBuilder::new();
        let section = g.push(GeometryNode::Profile(ipe(fillet))).expect("node");
        let beam = g
            .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
                profile: section,
                direction: Vec3::Z,
                depth: LEN,
            }))
            .expect("node");
        let graph = g.finish(vec![beam]).expect("graph");
        let (brep, report) = ReferenceExactCompiler::new()
            .compile_exact_with_report(&graph, beam, &ExecutionOptions::new(Tolerance::ZERO))
            .unwrap_or_else(|e| panic!("fillet {fillet:?}: {e}"));
        audited(&brep);
        assert!(report.is_exact(), "{report:?}");
        let volume = exact_properties(&brep, Tolerance::METRE)
            .expect("measurable")
            .signed_volume;
        let expected = area(fillet) * LEN;
        assert!(
            (volume - expected).abs() <= 1e-12 * expected,
            "fillet {fillet:?}: volume {volume}, closed form {expected}"
        );
    }
}

#[test]
fn a_web_hole_clear_of_the_fillets_cuts_an_ipe_beam_exactly_at_zero_tolerance() {
    // Through the web only, and right across the section past both flange
    // tips: either way the rod meets the web's two faces and nothing else.
    let runs = [(-TW, 2.0 * TW), (-0.75 * B, 1.5 * B)];
    for fillet in [None, Some(R)] {
        for run in runs {
            for placement in exact_placements() {
                let (graph, body) = beam_with_hole(fillet, run, placement);
                let (brep, report) = ReferenceExactCompiler::new()
                    .compile_exact_with_report(
                        &graph,
                        body,
                        &ExecutionOptions::new(Tolerance::ZERO),
                    )
                    .unwrap_or_else(|e| panic!("fillet {fillet:?}, run {run:?}: {e}"));
                audited(&brep);
                assert!(
                    report.is_exact(),
                    "fillet {fillet:?}, run {run:?}: {report:?}"
                );
                let volume = exact_properties(&brep, Tolerance::METRE)
                    .expect("measurable")
                    .signed_volume;
                // The hole's cylindrical walls are integrated numerically
                // by the measure, to ~1e-12 relative, not to the bit.
                let expected = area(fillet) * LEN - PI * HOLE * HOLE * TW;
                assert!(
                    (volume - expected).abs() <= 1e-11 * expected,
                    "fillet {fillet:?}, run {run:?}: volume {volume}, closed form {expected}"
                );
            }
        }
    }
}
