//! Booleans whose result encloses a cavity (#120).
//!
//! A tool buried inside the subject -- inside its plan and between its caps
//! -- leaves a closed cavity. The result is one solid whose void shell faces
//! into the cavity. It used to be refused because the mesh compiler dropped
//! void shells; it now tessellates them.
//!
//! Oracles are closed forms from the inputs, measured signed with
//! `exact_properties`, which sums every shell.

use axiolid_brep::ExactBRep;
use axiolid_construct::boolean_exact::{
    boolean_arc_prisms_exact, boolean_prisms_exact, boolean_prisms_exact_solids, ArcPrism, Prism,
};
use axiolid_core::{BooleanOperator, Point2, Tolerance};
use axiolid_overlay::ArcRing;
use axiolid_surface::Surface;

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Point2> {
    vec![
        Point2::new(x0, y0),
        Point2::new(x1, y0),
        Point2::new(x1, y1),
        Point2::new(x0, y1),
    ]
}

fn prism(ring: Vec<Point2>, bottom: f64, top: f64) -> Prism {
    Prism {
        rings: vec![ring],
        bottom,
        top,
    }
}

fn audit(solid: &ExactBRep) {
    let health = axiolid_brep_audit::geometric_audit(solid, tol());
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topo = axiolid_topology::audit_brep(solid.topology());
    assert!(topo.is_closed_manifold(), "{topo:?}");
}

fn volume(solid: &ExactBRep) -> f64 {
    axiolid_measure::exact_properties(solid, tol())
        .expect("measurable")
        .signed_volume
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-11 * want.abs(),
        "{what}: expected {want}, got {got}"
    );
}

#[test]
fn a_buried_box_leaves_a_void_shell() {
    let block = prism(rect(0.0, 0.0, 4.0, 3.0), 1.0, 3.0);
    let buried = prism(rect(1.0, 1.0, 2.0, 2.0), 1.5, 2.25);
    let solid = boolean_prisms_exact(&block, &buried, BooleanOperator::Difference, tol())
        .expect("a cavity is representable");
    audit(&solid);

    let solids = solid.topology().solids();
    assert_eq!(solids.len(), 1);
    assert_eq!(solids[0].voids.len(), 1, "one cavity, one void shell");
    close("volume", volume(&solid), 4.0 * 3.0 * 2.0 - 0.75);

    // The multi-solid entry point returns the same single piece.
    let pieces = boolean_prisms_exact_solids(&block, &buried, BooleanOperator::Difference, tol())
        .expect("one piece");
    assert_eq!(pieces.len(), 1);
    close("piece volume", volume(&pieces[0]), 24.0 - 0.75);
}

#[test]
fn a_round_cavity_keeps_its_cylinder_wall() {
    let (cx, cy, r) = (2.0, 1.5, 0.75);
    let block = ArcPrism {
        section: ArcRing::from_points(&rect(0.0, 0.0, 4.0, 3.0)),
        bottom: 0.0,
        top: 2.0,
    };
    let bore = ArcPrism {
        section: ArcRing::circle(Point2::new(cx, cy), r),
        bottom: 0.5,
        top: 1.25,
    };
    let solid = boolean_arc_prisms_exact(&block, &bore, BooleanOperator::Difference, tol())
        .expect("a round cavity is representable");
    audit(&solid);

    let solids = solid.topology().solids();
    assert_eq!(solids[0].voids.len(), 1);
    let cylinders = solid
        .surfaces()
        .iter()
        .filter(|surface| matches!(surface, Surface::Cylinder(_)))
        .count();
    assert!(cylinders >= 1, "the cavity wall stays a cylinder");
    close("volume", volume(&solid), 24.0 - PI * r * r * 0.75);
}
