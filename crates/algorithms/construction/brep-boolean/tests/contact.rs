//! Operands that touch rather than cross (#167, ADR 0075 stage 2): shared
//! face patches, sections along existing edges, and tangent contact.
//!
//! Same oracles as `boolean.rs`: audit, closed manifold, exact volume from
//! closed forms, and `|A u B| + |A n B| = |A| + |B|`.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{boolean, BooleanError};
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Point2, Tolerance};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn solid(section: ArcRing, bottom: f64, top: f64) -> ExactBRep {
    let p = ArcPrism {
        section,
        bottom,
        top,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a solid")
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> ArcRing {
    ArcRing::from_points(&[
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ])
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

fn attempt(a: &ExactBRep, b: &ExactBRep, op: BooleanOperator) -> Result<f64, BooleanError> {
    let result = boolean(a, b, op, tol())?;
    let health = geometric_audit(&result, tol());
    assert!(health.is_consistent(), "{op:?}: {:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{op:?}: {topology:?}");
    Ok(volume(&result))
}

fn run(a: &ExactBRep, b: &ExactBRep, op: BooleanOperator) -> f64 {
    attempt(a, b, op).unwrap_or_else(|e| panic!("{op:?}: {e}"))
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

/// All three operators; `None` for an empty intersection or difference.
fn check(
    a: &ExactBRep,
    b: &ExactBRep,
    union: f64,
    intersection: Option<f64>,
    difference: Option<f64>,
) {
    close("union", run(a, b, BooleanOperator::Union), union);
    for (op, want) in [
        (BooleanOperator::Intersection, intersection),
        (BooleanOperator::Difference, difference),
    ] {
        match (attempt(a, b, op), want) {
            (Ok(got), Some(want)) => close(&format!("{op:?}"), got, want),
            (Err(BooleanError::EmptyResult), None) => {}
            (got, want) => panic!("{op:?}: expected {want:?}, got {got:?}"),
        }
    }
    let i = intersection.unwrap_or(0.0);
    close("identity", union + i, volume(a) + volume(b));
}

#[test]
fn boxes_stacked_with_an_offset_share_part_of_a_face() {
    // B sits on A's roof; they share the square [1, 2]^2 at z = 1.
    let a = solid(square(0.0, 0.0, 2.0, 2.0), 0.0, 1.0);
    let b = solid(square(1.0, 1.0, 3.0, 3.0), 1.0, 2.0);
    check(&a, &b, 8.0, None, Some(4.0));
    let union = boolean(&a, &b, BooleanOperator::Union, tol()).expect("union");
    assert_eq!(
        union.topology().solids().len(),
        1,
        "joined through the patch"
    );
}

#[test]
fn boxes_side_by_side_on_one_floor_overlap_in_coplanar_faces() {
    // Floors and roofs coplanar and facing the same way.
    let a = solid(square(0.0, 0.0, 2.0, 2.0), 0.0, 1.0);
    let b = solid(square(1.0, 1.0, 3.0, 3.0), 0.0, 1.0);
    check(&a, &b, 7.0, Some(1.0), Some(3.0));
}

#[test]
fn a_box_against_itself() {
    let a = solid(square(0.0, 0.0, 2.0, 3.0), 0.0, 1.0);
    check(&a, &a.clone(), 6.0, Some(6.0), None);
}

#[test]
fn boxes_touching_along_a_whole_face() {
    let a = solid(square(0.0, 0.0, 2.0, 2.0), 0.0, 1.0);
    let b = solid(square(2.0, 0.0, 3.0, 2.0), 0.0, 1.0);
    check(&a, &b, 6.0, None, Some(4.0));
}

#[test]
fn a_box_plane_through_a_diamond_prism_edge() {
    // The box's face x = 0 contains two of the diamond's vertical edges,
    // and the diamond lies on both sides of it there.
    let diamond = solid(
        ArcRing::from_points(&[
            Point2::new(0.0, -1.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(-1.0, 0.0),
        ]),
        0.0,
        1.0,
    );
    let block = solid(square(0.0, -2.0, 2.0, 2.0), -1.0, 2.0);
    check(&diamond, &block, 25.0, Some(1.0), Some(1.0));
}

#[test]
fn a_pipe_touching_a_box_side_from_outside() {
    let r = 0.5;
    let a = solid(square(0.0, 0.0, 2.0, 2.0), 0.0, 1.0);
    let pipe = solid(ArcRing::circle(Point2::new(2.0 + r, 1.0), r), -1.0, 2.0);
    let vp = PI * r * r * 3.0;
    check(&a, &pipe, 4.0 + vp, None, Some(4.0));
}

#[test]
fn a_pipe_touching_a_box_side_from_inside() {
    // Its wall is tangent to the face x = 2, and its circles on the roof and
    // floor touch those faces' edges.
    let r = 0.5;
    let a = solid(square(0.0, 0.0, 2.0, 2.0), 0.0, 1.0);
    let pipe = solid(ArcRing::circle(Point2::new(2.0 - r, 1.0), r), -1.0, 2.0);
    let disc = PI * r * r;
    check(&a, &pipe, 4.0 + 2.0 * disc, Some(disc), Some(4.0 - disc));
}

#[test]
fn boxes_touching_along_an_edge_stay_two_manifold_solids() {
    // They share only the vertical edge x = y = 1: four faces meet there.
    let a = solid(square(0.0, 0.0, 1.0, 1.0), 0.0, 1.0);
    let b = solid(square(1.0, 1.0, 2.0, 2.0), 0.0, 1.0);
    check(&a, &b, 2.0, None, Some(1.0));
    let union = boolean(&a, &b, BooleanOperator::Union, tol()).expect("union");
    assert_eq!(union.topology().solids().len(), 2);
}

#[test]
fn boxes_touching_at_a_corner() {
    let a = solid(square(0.0, 0.0, 1.0, 1.0), 0.0, 1.0);
    let b = solid(square(1.0, 1.0, 2.0, 2.0), 1.0, 2.0);
    check(&a, &b, 2.0, None, Some(1.0));
}

#[test]
fn an_l_shaped_union_that_touches_itself_along_an_edge() {
    // A and B overlap, and B also touches A along an edge elsewhere: one
    // solid, with four faces meeting along that edge. Here: A an L of two
    // boxes [0,2]x[0,1] and [0,1]x[0,2]; B = [1,2]x[1,2] touches the L's
    // inner corner along the edge x = y = 1 and along two faces.
    let l = boolean(
        &solid(square(0.0, 0.0, 2.0, 1.0), 0.0, 1.0),
        &solid(square(0.0, 0.0, 1.0, 2.0), 0.0, 1.0),
        BooleanOperator::Union,
        tol(),
    )
    .expect("an L");
    close("L", volume(&l), 3.0);
    let b = solid(square(1.0, 1.0, 2.0, 2.0), 0.0, 1.0);
    check(&l, &b, 4.0, None, Some(3.0));
}
