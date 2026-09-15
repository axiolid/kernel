//! A torsion curve works as a graph directrix (ADR 0062).
//!
//! The end-to-end claim: `CurveRelation` and the sweep machinery are
//! GENERIC over curve family, so dispatching `Curve3::Intrinsic` in
//! `evaluate3` makes them work on a torsion curve with no change to the
//! relation code. A unit test on the evaluator cannot show that; this
//! walks the real graph and compiles a mesh.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame3, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{CurvatureLaw, Curve3, Intrinsic3};
use axiolid_measure::volume_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraphBuilder, GeometryNode, SolidOperation, Transition,
    TrimSelector, TrimmingPreference,
};
use axiolid_profile::{Profile, RectangleProfile};

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn rect() -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 0.4,
        y: 0.6,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// Helix of radius `a`, pitch `b`, as an arc-length natural equation.
fn helix(length: Scalar) -> Intrinsic3 {
    let (a, b): (Scalar, Scalar) = (3.0, 1.5);
    let c = a.hypot(b);
    Intrinsic3::new(
        Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::new(1.0, 0.0, 0.0),
            y: Vec3::new(0.0, 1.0, 0.0),
            z: Vec3::new(0.0, 0.0, 1.0),
        },
        CurvatureLaw::circular(a / (c * c)),
        CurvatureLaw::circular(b / (c * c)),
        length,
    )
}

/// Sweep `rect` along a directrix node built by `build`.
fn sweep(
    build: impl FnOnce(&mut GeometryGraphBuilder) -> axiolid_model::NodeId,
) -> axiolid_contracts::GeomResult<axiolid_mesh::TriMesh> {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect())).unwrap();
    let directrix = build(&mut b);
    let swept = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::FixedReferenceSweep {
                profile,
                directrix,
                reference_direction: Vec3::Z,
                parameter_range: None,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![swept]).unwrap();
    compiler().compile_mesh(&graph, swept, &options())
}

/// Volume of a swept solid, or a panic naming the failure.
fn volume_of(mesh: &axiolid_mesh::TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .expect("a swept solid must be closed and two-manifold")
        .signed_volume
        .abs()
}

#[test]
fn a_torsion_curve_is_a_valid_sweep_directrix() {
    // Before Intrinsic was dispatched in evaluate3, this refused outright.
    // Pappus gives the volume: section area times the centroid's travel,
    // and for an arc-length curve that travel IS the length.
    let length = 20.0;
    let mesh = sweep(|b| {
        b.push(GeometryNode::Curve3(Curve3::Intrinsic(helix(length))))
            .unwrap()
    })
    .expect("a torsion curve is a valid directrix");

    let exact = 0.4 * 0.6 * length;
    let ratio = volume_of(&mesh) / exact;
    // A swept solid on a curved path is not exactly the prism volume --
    // the section sweeps slightly less on the inside of the bend. Bounded
    // rather than pinned, so the test states a real property.
    assert!(
        (0.90..=1.02).contains(&ratio),
        "swept volume ratio {ratio} against the prism estimate"
    );
}

#[test]
fn a_trimmed_torsion_curve_sweeps_only_the_kept_span() {
    // CurveRelation::Trimmed is generic: it trims by PARAMETER, and for an
    // arc-length curve the parameter is arc length. Half the span must
    // sweep about half the volume -- which also proves the relation is
    // reading the curve's real domain, not a unit interval.
    let length = 20.0;
    let mesh = sweep(|b| {
        let basis = b
            .push(GeometryNode::Curve3(Curve3::Intrinsic(helix(length))))
            .unwrap();
        b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start: vec![TrimSelector::Parameter(4.0)],
            end: vec![TrimSelector::Parameter(14.0)],
            sense_agreement: true,
            preference: TrimmingPreference::Parameter,
        }))
        .unwrap()
    })
    .expect("a trimmed torsion curve is a valid directrix");

    let exact = 0.4 * 0.6 * 10.0;
    let ratio = volume_of(&mesh) / exact;
    assert!(
        (0.90..=1.02).contains(&ratio),
        "trimmed sweep ratio {ratio}: the trim must cut the span, not the whole curve"
    );
}

#[test]
fn a_composite_of_torsion_curves_stitches_end_to_end() {
    // CurveRelation::Composite stitches children by proximity of endpoints.
    // Two helices of the same law, the second starting where the first
    // ends, must sweep about the sum of their lengths.
    let first = helix(8.0);
    let mut second = helix(8.0);
    // Anchor the second at the first's endpoint, with the frame it has
    // there, so the two genuinely meet.
    second.start = axiolid_reference::frenet_frame(&first, 8.0).expect("end frame");

    let mesh = sweep(|b| {
        let a = b
            .push(GeometryNode::Curve3(Curve3::Intrinsic(first)))
            .unwrap();
        let c = b
            .push(GeometryNode::Curve3(Curve3::Intrinsic(second)))
            .unwrap();
        b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
            segments: vec![
                CurveSegment {
                    curve: a,
                    same_sense: true,
                    transition: Transition::ContinuousSameGradient,
                },
                CurveSegment {
                    curve: c,
                    same_sense: true,
                    transition: Transition::ContinuousSameGradient,
                },
            ],
        }))
        .unwrap()
    })
    .expect("a composite of torsion curves is a valid directrix");

    let exact = 0.4 * 0.6 * 16.0;
    let ratio = volume_of(&mesh) / exact;
    assert!(
        (0.90..=1.02).contains(&ratio),
        "composite sweep ratio {ratio}: both halves must contribute"
    );
}
