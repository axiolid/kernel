//! Persistent names must survive the rebuilds that invalidate arena indices.
//!
//! The point of a name is that it outlives the assembly it was minted in.
//! An arena index does not: rebuild with a different face count and slot 7
//! is a different face. These tests pin the difference.

use axiolid_brep::{EdgeName, FaceName, SweptFace};
use axiolid_construct::feature::{fillet_extruded_profile, EdgeSelector, FeatureSize};
use axiolid_core::{Point2, Scalar, Tolerance, Vec3};
use axiolid_profile::{Profile, RectangleProfile};

fn rectangle(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The blend face is named after the edge it replaced, not its own slot.
#[test]
fn a_blend_is_named_after_the_edge_it_replaced() {
    let solid = fillet_extruded_profile(
        &rectangle(4.0, 6.0),
        Vec3::Z,
        2.0,
        EdgeSelector::NearestCorner(Point2::new(2.0, 3.0)),
        FeatureSize::ConstantRadius(0.5),
        Tolerance::METRE,
    )
    .expect("a rectangle fillet is supported");

    // The +x/+y corner sits between walls 1 and 2 of the four-sided ring.
    let expected = FaceName::blend(EdgeName::between(
        FaceName::swept(SweptFace::Side(1)),
        FaceName::swept(SweptFace::Side(2)),
    ));
    assert!(
        solid.face_by_name(&expected).is_some(),
        "the blend must carry the name of the edge it replaced"
    );
}

/// A name minted on one solid selects the same edge after an upstream edit.
///
/// This is the property an arena index cannot provide, and the reason the
/// fillet previously had to re-locate its target by nearest-point search.
#[test]
fn a_name_survives_a_profile_resize() {
    let corner_edge = EdgeName::between(
        FaceName::swept(SweptFace::Side(1)),
        FaceName::swept(SweptFace::Side(2)),
    );

    // Same NAME, two different profiles. The corner is at a different
    // coordinate in each, so a positional selector would need updating.
    for (x, y) in [(4.0, 6.0), (10.0, 3.0)] {
        let solid = fillet_extruded_profile(
            &rectangle(x, y),
            Vec3::Z,
            2.0,
            EdgeSelector::Named(corner_edge.clone()),
            FeatureSize::ConstantRadius(0.5),
            Tolerance::METRE,
        )
        .expect("a named edge resolves on any rectangle");

        let blend = FaceName::blend(corner_edge.clone());
        assert!(
            solid.face_by_name(&blend).is_some(),
            "the {x} x {y} solid must carry the blend named by that edge"
        );
    }
}

/// Naming an edge from either side gives the same name.
///
/// Without this an edge has two names, and a feature applied from the other
/// side silently misses.
#[test]
fn an_edge_name_does_not_depend_on_which_face_is_mentioned_first() {
    let a = FaceName::swept(SweptFace::Side(1));
    let b = FaceName::swept(SweptFace::Side(2));
    assert_eq!(
        EdgeName::between(a.clone(), b.clone()),
        EdgeName::between(b, a),
        "the edge where two faces meet is one edge, so it is one name"
    );
}

/// Both selector spellings of the same corner produce the same solid.
#[test]
fn a_named_selector_agrees_with_the_positional_one() {
    let named = fillet_extruded_profile(
        &rectangle(4.0, 6.0),
        Vec3::Z,
        2.0,
        EdgeSelector::Named(EdgeName::between(
            FaceName::swept(SweptFace::Side(1)),
            FaceName::swept(SweptFace::Side(2)),
        )),
        FeatureSize::ConstantRadius(0.5),
        Tolerance::METRE,
    )
    .expect("named selection is supported");
    let positional = fillet_extruded_profile(
        &rectangle(4.0, 6.0),
        Vec3::Z,
        2.0,
        EdgeSelector::NearestCorner(Point2::new(2.0, 3.0)),
        FeatureSize::ConstantRadius(0.5),
        Tolerance::METRE,
    )
    .expect("positional selection is supported");

    assert_eq!(
        named.topology().faces().len(),
        positional.topology().faces().len(),
        "the two spellings must select the same corner"
    );
    let blend = FaceName::blend(EdgeName::between(
        FaceName::swept(SweptFace::Side(1)),
        FaceName::swept(SweptFace::Side(2)),
    ));
    assert!(named.face_by_name(&blend).is_some());
    assert!(positional.face_by_name(&blend).is_some());
}

/// An unresolvable name is refused, not silently snapped to a near corner.
#[test]
fn an_unknown_name_is_refused_rather_than_guessed() {
    let result = fillet_extruded_profile(
        &rectangle(4.0, 6.0),
        Vec3::Z,
        2.0,
        EdgeSelector::Named(EdgeName::between(
            FaceName::swept(SweptFace::Side(7)),
            FaceName::swept(SweptFace::Side(9)),
        )),
        FeatureSize::ConstantRadius(0.5),
        Tolerance::METRE,
    );
    assert!(
        result.is_err(),
        "a name with no matching edge must refuse, not pick the nearest"
    );
}
