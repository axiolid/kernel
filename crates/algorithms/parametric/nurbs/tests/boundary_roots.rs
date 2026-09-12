//! Roots lying exactly ON a patch domain edge.
//!
//! Ordinary Krawczyk certification proves a root exists by mapping a
//! parameter box strictly inside itself. A root sitting exactly on a face of
//! that box can never satisfy the strict inequality, so subdivision halves
//! the box forever and the query returns `Unresolved` at any tolerance.
//!
//! The boundary-restricted certificate fixes the parameter that is pinned to
//! the domain edge and certifies the reduced system in the remaining two
//! unknowns, where the root IS interior.

use axiolid_core::Point3;
use axiolid_curve::{BSplineCurve3, KnotSpec};
use axiolid_nurbs::{
    intersect_curve_surface_certified, CertifiedCurveSurfaceIntersection3,
    CertifiedCurveSurfaceIntersectionOptions,
};
use axiolid_surface::BSplineSurface;

/// A unit plane patch spanning x,y in [-1,1] at z=0.
///
/// Its own parameter domain is u in [10,12], v in [-4,0], deliberately
/// unrelated to the model coordinates so a test cannot accidentally pass by
/// confusing parameter space with model space.
fn plane() -> BSplineSurface {
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(-1.0, -1.0, 0.0), Point3::new(-1.0, 1.0, 0.0)],
            vec![Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
        u_knots: vec![10.0, 12.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![-4.0, 0.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

/// A vertical line piercing z=0 at (x, y).
fn vertical_line(x: f64, y: f64) -> BSplineCurve3 {
    BSplineCurve3 {
        degree: 1,
        control_points: vec![Point3::new(x, y, -1.0), Point3::new(x, y, 2.0)],
        knots: vec![-2.0, 2.0],
        multiplicities: vec![2, 2],
        weights: None,
        knot_spec: KnotSpec::Unspecified,
        closed: false,
        self_intersect: None,
    }
}

fn options(tolerance: f64) -> CertifiedCurveSurfaceIntersectionOptions {
    CertifiedCurveSurfaceIntersectionOptions::new(tolerance, 100_000, 64)
        .expect("valid certified policy")
}

/// A root strictly inside the patch still certifies, unchanged.
///
/// Guards the ordinary path: the boundary-restricted attempt must not
/// disturb, duplicate, or replace an interior certificate.
#[test]
fn an_interior_root_still_certifies_exactly_once() {
    let result =
        intersect_curve_surface_certified(&vertical_line(0.25, 0.5), &plane(), options(1.0e-6))
            .expect("certified query terminates");

    match result {
        CertifiedCurveSurfaceIntersection3::Complete { intersections, .. } => {
            assert_eq!(
                intersections.len(),
                1,
                "one piercing point means exactly one certified root"
            );
        }
        other => panic!("interior root must certify, got {other:?}"),
    }
}

/// A root exactly ON the patch's v_start edge now certifies.
///
/// `y = -1.0` is the model image of the patch's own `v = -4` domain edge, so
/// the piercing point lands precisely on the boundary. Before the
/// boundary-restricted certificate this returned `Unresolved` at every
/// tolerance, because the Krawczyk image could not be proven strictly inside
/// a box whose face the root sits on.
#[test]
fn a_root_on_the_v_start_domain_edge_certifies() {
    let result =
        intersect_curve_surface_certified(&vertical_line(0.25, -1.0), &plane(), options(1.0e-6))
            .expect("certified query terminates");

    match result {
        CertifiedCurveSurfaceIntersection3::Complete { intersections, .. } => {
            assert_eq!(
                intersections.len(),
                1,
                "the boundary piercing point must be certified, not abandoned"
            );
            let root = &intersections[0];
            assert!(
                (root.point.y - (-1.0)).abs() < 1.0e-9,
                "certified point must lie on the edge, got y={}",
                root.point.y
            );
        }
        other => panic!("a boundary root must now certify, got {other:?}"),
    }
}

/// A root exactly on the u_start edge certifies too.
///
/// Exercises the other pinned axis: `x = -1.0` is the image of `u = 10`.
/// Without this, a bug that only handled the v family would pass unnoticed.
#[test]
fn a_root_on_the_u_start_domain_edge_certifies() {
    let result =
        intersect_curve_surface_certified(&vertical_line(-1.0, 0.25), &plane(), options(1.0e-6))
            .expect("certified query terminates");

    match result {
        CertifiedCurveSurfaceIntersection3::Complete { intersections, .. } => {
            assert_eq!(intersections.len(), 1, "u-edge root must be certified");
            let root = &intersections[0];
            assert!(
                (root.point.x - (-1.0)).abs() < 1.0e-9,
                "certified point must lie on the u_start edge, got x={}",
                root.point.x
            );
        }
        other => panic!("a u-edge boundary root must certify, got {other:?}"),
    }
}

/// A line missing the patch entirely is still proven disjoint.
///
/// The boundary machinery must not manufacture a root from a near miss: the
/// line passes outside the patch, so the answer is an empty Complete, not a
/// certificate.
#[test]
fn a_line_outside_the_patch_certifies_no_root() {
    let result =
        intersect_curve_surface_certified(&vertical_line(3.0, 3.0), &plane(), options(1.0e-6))
            .expect("certified query terminates");

    match result {
        CertifiedCurveSurfaceIntersection3::Complete { intersections, .. } => {
            assert!(
                intersections.is_empty(),
                "a line outside the patch has no roots, got {}",
                intersections.len()
            );
        }
        other => panic!("a clear miss must be proven empty, got {other:?}"),
    }
}

/// A line meeting the edge's supporting line but missing the patch is refused.
///
/// The reduced 2x2 solve only enforces two of the three coordinate
/// equations. This line crosses z=0 at y=-1 -- on the v_start edge's
/// supporting line -- but at x=5.0, far outside the patch's x extent. If the
/// discarded third row were assumed rather than checked, the solve would
/// report a root that does not exist on the patch.
#[test]
fn a_root_on_the_edge_line_but_off_the_patch_is_refused() {
    let result =
        intersect_curve_surface_certified(&vertical_line(5.0, -1.0), &plane(), options(1.0e-6))
            .expect("certified query terminates");

    match result {
        CertifiedCurveSurfaceIntersection3::Complete { intersections, .. } => {
            assert!(
                intersections.is_empty(),
                "a point off the patch must not be certified, got {} root(s) at {:?}",
                intersections.len(),
                intersections.first().map(|root| root.point)
            );
        }
        // A conservative Unresolved is also acceptable here: the contract is
        // that no FALSE root is ever certified, not that every miss is proven.
        CertifiedCurveSurfaceIntersection3::Unresolved { intersections, .. } => {
            assert!(
                intersections.is_empty(),
                "an unresolved outcome must still not certify a false root"
            );
        }
        other => panic!("unexpected outcome {other:?}"),
    }
}
