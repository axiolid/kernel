//! An elevated alignment centreline survives a round trip through the graph
//! (ADR 0060, issue #105).
//!
//! The issue's second acceptance criterion. Storing the composition as an
//! ordinary `Curve3` is the whole point: if the graph had to learn a special
//! case for it, the value would not really be a curve.

use axiolid_core::{Frame2, Point2, Vec2};
use axiolid_curve::{CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2};
use axiolid_model::{GeometryGraphBuilder, GeometryNode};

/// A clothoid plan with a parabolic vertical profile: a real alignment shape.
fn centreline() -> Curve3 {
    let (radius, length) = (300.0, 120.0);
    let plan = Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::new(1.0, 0.0),
            y: Vec2::new(0.0, 1.0),
        },
        CurvatureLaw::clothoid(0.0, 1.0 / radius, length),
        length,
    ));
    Curve3::Elevated(Elevated3::new(
        plan,
        ElevationLaw::parabolic(100.0, 0.02, -0.02, length),
    ))
}

#[test]
fn an_elevated_centreline_round_trips_through_the_graph() {
    let original = centreline();
    let mut builder = GeometryGraphBuilder::default();
    let node = builder
        .push(GeometryNode::Curve3(original.clone()))
        .expect("an elevated curve is a valid geometry node");
    // Finishing validates the graph, so this also proves the new variant is
    // accepted as a well-formed 3D curve rather than merely stored.
    let graph = builder.finish(vec![node]).expect("the graph validates");

    let stored = graph.get(node).expect("the node is present");
    let GeometryNode::Curve3(recovered) = stored else {
        panic!("a curve3 node must come back as a curve3 node");
    };

    // Structural equality: both laws recovered unchanged, not refitted.
    assert_eq!(
        &original, recovered,
        "the composition must survive storage byte for byte"
    );

    let Curve3::Elevated(elevated) = recovered else {
        panic!("the elevated variant must survive");
    };
    let Curve2::Intrinsic(plan) = elevated.plan.as_ref() else {
        panic!("the plan must still be an exact intrinsic curve, not a fitted spline");
    };
    // The spiral is still a spiral: curvature ramps from 0 to 1/300.
    assert!(
        !plan.is_straight(),
        "the transition spiral must not have been flattened"
    );
    assert_eq!(plan.length, 120.0);
}
