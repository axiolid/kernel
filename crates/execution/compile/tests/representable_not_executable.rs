//! Representable is not executable (#20).
//!
//! If the neutral model could only hold what a provider can execute today,
//! importers would have to discard data at the door -- and the loss would be
//! silent. The model must therefore carry geometry no provider can yet build,
//! and executing it must produce a typed refusal that NAMES the gap.
//!
//! Two properties, and the pair is the point:
//!
//! - the data survives a round trip through the graph unchanged, so nothing is
//!   quietly dropped on the way in;
//! - asking to execute it refuses by input family, so a caller learns which
//!   capability is missing rather than that "something" failed.
//!
//! A refusal that says only "unsupported operation" fails the second half: it
//! cannot tell a caller whether to register a revolution provider or a sweep
//! provider, which is the actionable part.

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::{Point3, Tolerance, Vec3};
use axiolid_curve::{Curve3, Polyline3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, SolidOperation};
use axiolid_profile::{Profile, RectangleProfile};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn rect() -> Profile {
    Profile::Rectangle(RectangleProfile {
        x: 1.0,
        y: 1.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// A revolution round-trips through the graph with every field intact.
///
/// No exact provider can build one, but the importer must still be able to
/// put it in the model and read it back byte-identical. Losing the axis or
/// the angle here would be an unrecoverable import loss.
#[test]
fn a_revolution_survives_a_round_trip_through_the_model() {
    let mut builder = GeometryGraphBuilder::new();
    let profile = builder.push(GeometryNode::Profile(rect())).unwrap();
    let axis_origin = Point3::new(2.0, 0.0, 0.0);
    let axis_direction = Vec3::Z;
    let angle = 1.234_567_89;

    let solid = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Revolution {
            profile,
            axis_origin,
            axis_direction,
            angle,
        }))
        .unwrap();
    let graph = builder.finish(vec![solid]).unwrap();

    let Some(GeometryNode::SolidOperation(SolidOperation::Revolution {
        axis_origin: stored_origin,
        axis_direction: stored_direction,
        angle: stored_angle,
        ..
    })) = graph.get(solid)
    else {
        panic!("the model must still hold the revolution it was given");
    };

    assert_eq!(*stored_origin, axis_origin, "axis origin was not preserved");
    assert_eq!(
        *stored_direction, axis_direction,
        "axis direction was not preserved"
    );
    assert_eq!(
        *stored_angle, angle,
        "angle must survive exactly, not within a tolerance"
    );
}

/// Executing a still-unsupported family refuses, naming that family.
///
/// The named family is the actionable word: it tells a caller which provider
/// to register. A bare "unsupported operation" would not.
///
/// This deliberately uses a family with no exact provider. Exact revolution
/// landed for the offset full-turn case, so revolution is no longer a valid
/// example -- a test asserting it is refused would now be asserting the
/// absence of a capability that exists.
#[test]
fn executing_an_unsupported_family_refuses_by_name() {
    let mut builder = GeometryGraphBuilder::new();
    let profile = builder.push(GeometryNode::Profile(rect())).unwrap();
    let directrix = builder
        .push_value(Curve3::Polyline(Polyline3 {
            points: vec![Point3::ZERO, Point3::new(0.0, 0.0, 3.0)],
            closed: false,
        }))
        .unwrap();
    let solid = builder
        .push(GeometryNode::SolidOperation(
            SolidOperation::FixedReferenceSweep {
                profile,
                directrix,
                reference_direction: Vec3::X,
                parameter_range: None,
            },
        ))
        .unwrap();
    let graph = builder.finish(vec![solid]).unwrap();

    let error = ReferenceExactCompiler::new()
        .compile_exact(&graph, solid, &options())
        .expect_err("the exact compiler has no sweep provider wired yet");

    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::GraphCompilation,
                input: "fixed-reference sweep",
                ..
            }
        ),
        "the refusal must name the missing family, got {error:?}"
    );
}

/// Distinct unsupported families produce DISTINCT refusals.
///
/// This is the test that would fail against a compiler collapsing everything
/// into one generic "unsupported" answer. If two different missing
/// capabilities are indistinguishable to a caller, the diagnostic carries no
/// information and the caller cannot act on it.
#[test]
fn different_unsupported_families_are_named_differently() {
    type BuildFamily = Box<dyn Fn(&mut GeometryGraphBuilder) -> axiolid_model::NodeId>;
    let families: Vec<(&str, BuildFamily)> = vec![
        (
            "tapered revolution",
            Box::new(|b: &mut GeometryGraphBuilder| {
                let start_profile = b.push(GeometryNode::Profile(rect())).unwrap();
                let end_profile = b.push(GeometryNode::Profile(rect())).unwrap();
                b.push(GeometryNode::SolidOperation(
                    SolidOperation::TaperedRevolution {
                        start_profile,
                        end_profile,
                        axis_origin: Point3::new(-3.0, 0.0, 0.0),
                        axis_direction: Vec3::Y,
                        angle: 1.0,
                    },
                ))
                .unwrap()
            }),
        ),
        (
            "fixed-reference sweep",
            Box::new(|b: &mut GeometryGraphBuilder| {
                let profile = b.push(GeometryNode::Profile(rect())).unwrap();
                let directrix = b
                    .push_value(Curve3::Polyline(Polyline3 {
                        points: vec![Point3::ZERO, Point3::new(0.0, 0.0, 3.0)],
                        closed: false,
                    }))
                    .unwrap();
                b.push(GeometryNode::SolidOperation(
                    SolidOperation::FixedReferenceSweep {
                        profile,
                        directrix,
                        reference_direction: Vec3::X,
                        parameter_range: None,
                    },
                ))
                .unwrap()
            }),
        ),
    ];

    let mut seen: Vec<&str> = Vec::new();
    for (expected, build) in families {
        let mut builder = GeometryGraphBuilder::new();
        let root = build(&mut builder);
        let graph = builder.finish(vec![root]).unwrap();

        let error = ReferenceExactCompiler::new()
            .compile_exact(&graph, root, &options())
            .expect_err("no exact provider exists for this family yet");

        let GeomError::UnsupportedInput { input, .. } = error else {
            panic!("expected a typed input-family refusal, got {error:?}");
        };
        assert_eq!(input, expected, "refusal named the wrong family");
        assert!(
            !seen.contains(&input),
            "two families shared the refusal name {input:?}, so the diagnostic is not actionable"
        );
        seen.push(input);
    }
    assert_eq!(seen.len(), 2);
}

/// Every declared solid family has a distinct, non-empty diagnostic name.
///
/// `SolidOperation` is `#[non_exhaustive]`, so the mesh compiler keeps a
/// catch-all arm for families added later. That arm used to answer with a bare
/// `Unsupported { operation: Sweep }`, which tells a caller that "a sweep"
/// failed but not WHICH capability is missing -- revolution, swept disk and
/// fixed-reference sweep were indistinguishable. It now reports the family,
/// matching what the exact path already did.
///
/// The naming table is what both paths depend on, so this pins it directly:
/// every family distinct, none blank, none a generic placeholder.
#[test]
fn every_solid_family_has_a_distinct_diagnostic_name() {
    let names = axiolid_mesh_compile::SOLID_FAMILY_NAMES;
    assert!(
        names.len() >= 10,
        "expected every declared family, got {names:?}"
    );

    let mut sorted: Vec<&str> = names.to_vec();
    sorted.sort_unstable();
    let before = sorted.len();
    sorted.dedup();
    assert_eq!(
        before,
        sorted.len(),
        "two families share a diagnostic name, so a caller cannot tell them apart"
    );

    for name in names {
        assert!(!name.is_empty(), "a family has no diagnostic name");
        assert!(
            !name.contains("unknown") && !name.contains("unsupported"),
            "{name:?} is a placeholder, not an actionable family name"
        );
    }
}
