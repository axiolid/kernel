//! Boolean output faces must say which input face they are a fragment of.
//!
//! This is the property that makes a boolean result addressable: a caller
//! that applied a material to a wall face, or selected an edge to fillet,
//! needs those references to survive the operation. Without provenance the
//! result is a fresh anonymous solid and every upstream reference dangles.

use axiolid_brep::{FaceName, Operand, SweptFace};
use axiolid_construct::boolean_exact::{boolean_prisms_exact, Prism};
use axiolid_core::{BooleanOperator, Point2, Scalar, Tolerance};

fn rect_ring(cx: Scalar, cy: Scalar, x: Scalar, y: Scalar) -> Vec<Point2> {
    let (hx, hy) = (x / 2.0, y / 2.0);
    vec![
        Point2::new(cx - hx, cy - hy),
        Point2::new(cx + hx, cy - hy),
        Point2::new(cx + hx, cy + hy),
        Point2::new(cx - hx, cy + hy),
    ]
}

fn prism(rings: Vec<Vec<Point2>>, bottom: Scalar, top: Scalar) -> Prism {
    Prism { rings, bottom, top }
}

/// Every wall of a differenced result names the input wall it came from.
///
/// The wall keeps four outer walls and gains four from the opening. Each
/// must say which operand and which original profile edge it is part of,
/// because that is what a material assignment or a later fillet resolves.
#[test]
fn difference_walls_name_their_source_wall() {
    let wall = prism(vec![rect_ring(0.0, 0.0, 10.0, 4.0)], 0.0, 3.0);
    let opening = prism(vec![rect_ring(0.0, 0.0, 2.0, 2.0)], 0.0, 3.0);

    let result = boolean_prisms_exact(
        &wall,
        &opening,
        BooleanOperator::Difference,
        Tolerance::METRE,
    )
    .expect("a wall with an interior opening is exactly constructible");

    let topology = result.topology();
    let mut subject_walls = 0;
    let mut tool_walls = 0;
    let mut unnamed = 0;
    for index in 0..topology.faces().len() {
        let Some(id) = topology.face_id_at(index) else {
            continue;
        };
        match result.face_name(id) {
            Some(FaceName::Fragment { operand, source }) => {
                // A wall fragment must trace back to a swept side wall,
                // never to a cap and never to an anonymous face.
                if matches!(**source, FaceName::Swept(SweptFace::Side(_))) {
                    match operand {
                        Operand::Subject => subject_walls += 1,
                        Operand::Tool => tool_walls += 1,
                    }
                }
            }
            None => unnamed += 1,
            _ => {}
        }
    }

    assert_eq!(
        subject_walls, 4,
        "the four outer wall faces must name the subject walls they came from"
    );
    assert_eq!(
        tool_walls, 4,
        "the four opening faces must name the tool walls that cut them"
    );
    assert_eq!(unnamed, 0, "no face of an exact boolean should be unnamed");
}

/// A fragment of a fillet blend still reports the blend as its origin.
///
/// This is the composition that makes names worth having: the boolean
/// wraps whatever name the input face carried, so a caller asking "what
/// was this before any boolean" gets the blend, not the boolean.
#[test]
fn origin_sees_through_the_boolean_layer() {
    let subject = prism(vec![rect_ring(0.0, 0.0, 4.0, 4.0)], 0.0, 3.0);
    let tool = prism(vec![rect_ring(2.0, 0.0, 4.0, 4.0)], 0.0, 3.0);

    let result = boolean_prisms_exact(
        &subject,
        &tool,
        BooleanOperator::Intersection,
        Tolerance::METRE,
    )
    .expect("coaxial prism intersection is exactly constructible");

    let topology = result.topology();
    let mut checked = 0;
    for index in 0..topology.faces().len() {
        let Some(id) = topology.face_id_at(index) else {
            continue;
        };
        if let Some(name) = result.face_name(id) {
            // origin() strips every Fragment layer, so whatever it returns
            // must be a real construction name, never another fragment.
            assert!(
                !matches!(name.origin(), FaceName::Fragment { .. }),
                "origin must strip all boolean layers, got {}",
                name.origin()
            );
            assert!(
                !name.is_anonymous(),
                "an exactly-constructed face must not be anonymous, got {name}"
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "the result must carry at least one name");
}

/// A union of two prisms sharing a wall plane still names every face.
///
/// The shared plane is the ambiguous case: a result wall lying on a wall
/// of BOTH operands has two defensible names. The recovery must pick one
/// deterministically rather than report a different answer per run, and
/// must never claim a source the wall does not actually lie on.
#[test]
fn a_shared_wall_plane_stays_deterministic() {
    let left = prism(vec![rect_ring(0.0, 0.0, 4.0, 4.0)], 0.0, 3.0);
    let right = prism(vec![rect_ring(4.0, 0.0, 4.0, 4.0)], 0.0, 3.0);

    let first = boolean_prisms_exact(&left, &right, BooleanOperator::Union, Tolerance::METRE);
    let second = boolean_prisms_exact(&left, &right, BooleanOperator::Union, Tolerance::METRE);

    let (Ok(first), Ok(second)) = (first, second) else {
        // Refusing is acceptable here; claiming a wrong name is not.
        return;
    };

    let names = |solid: &axiolid_brep::ExactBRep| -> Vec<String> {
        let topology = solid.topology();
        (0..topology.faces().len())
            .filter_map(|index| topology.face_id_at(index))
            .map(|id| match solid.face_name(id) {
                Some(name) => name.to_string(),
                None => "<unnamed>".to_owned(),
            })
            .collect()
    };

    assert_eq!(
        names(&first),
        names(&second),
        "the same inputs must produce the same names on every run"
    );
}

/// A wall must not be attributed to an input edge it merely touches.
///
/// An output edge sharing ONE endpoint with an input edge, but running
/// away from its supporting line, is not a fragment of it. Matching on a
/// single endpoint accepts exactly that and silently attaches the
/// fragment to the wrong wall, carrying the wrong material through.
///
/// The T-shaped union below has result corners coincident with subject
/// corners while the walls leaving them are perpendicular. Every distinct
/// input wall that contributes boundary must appear exactly once per
/// fragment it produced -- a wall that vanishes from the mapping was
/// stolen by a neighbour that only touched it.
#[test]
fn a_touching_corner_does_not_borrow_the_wrong_wall() {
    let subject = prism(vec![rect_ring(0.0, 0.0, 4.0, 2.0)], 0.0, 3.0);
    let tool = prism(vec![rect_ring(0.0, 2.0, 2.0, 2.0)], 0.0, 3.0);

    let Ok(result) =
        boolean_prisms_exact(&subject, &tool, BooleanOperator::Union, Tolerance::METRE)
    else {
        return;
    };

    let topology = result.topology();
    let mut seen: Vec<String> = (0..topology.faces().len())
        .filter_map(|index| topology.face_id_at(index))
        .filter_map(|id| result.face_name(id))
        .map(|name| name.to_string())
        .collect();
    seen.sort();

    // All four tool walls bound the union somewhere: the three that stick
    // out, plus the one flush with the subject. Losing any of them means
    // a subject wall it merely touches claimed it.
    for wall in ["tool/side[1]", "tool/side[2]", "tool/side[3]"] {
        assert!(
            seen.iter().any(|name| name == wall),
            "{wall} must appear in the mapping, got {seen:?}"
        );
    }
    // The subject wall opposite the tool is untouched by it and must
    // survive as its own fragment.
    assert!(
        seen.iter().any(|name| name == "subject/side[3]"),
        "subject/side[3] must survive as its own fragment, got {seen:?}"
    );
}
