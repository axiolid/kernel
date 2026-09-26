//! The general exact boolean, stage 1 (#167, ADR 0075).
//!
//! Every result must audit clean (faces agree with their pcurves and share
//! edges consistently) and measure exactly (`exact_properties`, ADR 0073)
//! to the closed form from the inputs. Vertical-column cases are also
//! checked against the column builder (ADR 0072), an independent exact
//! path, and every pair against `|A u B| + |A n B| = |A| + |B|`.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, clip_arc_prism_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Plane3, Point2, Point3, Tolerance, Vec3};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use axiolid_primitive::HalfSpace;

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn prism(section: ArcRing, bottom: f64, top: f64) -> ArcPrism {
    ArcPrism {
        section,
        bottom,
        top,
    }
}

fn solid(section: ArcRing, bottom: f64, top: f64) -> ExactBRep {
    let p = prism(section, bottom, top);
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

fn run(a: &ExactBRep, b: &ExactBRep, op: BooleanOperator) -> ExactBRep {
    let result = boolean(a, b, op, tol()).unwrap_or_else(|e| panic!("{op:?}: {e}"));
    let health = geometric_audit(&result, tol());
    assert!(health.is_consistent(), "{op:?}: {:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{op:?}: {topology:?}");
    result
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

/// All three operators, the identity, and closed forms where given.
fn check(a: &ExactBRep, b: &ExactBRep, union: f64, intersection: f64, difference: f64) {
    let u = volume(&run(a, b, BooleanOperator::Union));
    let i = volume(&run(a, b, BooleanOperator::Intersection));
    let d = volume(&run(a, b, BooleanOperator::Difference));
    close("union", u, union);
    close("intersection", i, intersection);
    close("difference", d, difference);
    close("identity", u + i, volume(a) + volume(b));
    close("difference = A - (A n B)", d, volume(a) - i);
}

#[test]
fn a_pipe_through_a_box() {
    let block = solid(square(-1.0, -1.0, 1.0, 1.0), 0.0, 2.0);
    let r = 0.5;
    let pipe = solid(ArcRing::circle(Point2::new(0.2, -0.1), r), -1.0, 3.0);
    let disc = PI * r * r;
    check(
        &block,
        &pipe,
        8.0 + 4.0 * disc - 2.0 * disc,
        2.0 * disc,
        8.0 - 2.0 * disc,
    );
}

#[test]
fn a_pipe_at_a_box_corner_agrees_with_the_column_builder() {
    let (r, cx, cy) = (0.5, 1.0, 0.9);
    let box_prism = prism(square(-1.0, -1.0, 1.0, 1.0), 0.0, 2.0);
    let pipe_prism = prism(ArcRing::circle(Point2::new(cx, cy), r), -1.0, 3.0);
    let columns = |op| {
        volume(&boolean_arc_prisms_exact(&box_prism, &pipe_prism, op, tol()).expect("columns"))
    };
    let block = solid(square(-1.0, -1.0, 1.0, 1.0), 0.0, 2.0);
    let pipe = solid(ArcRing::circle(Point2::new(cx, cy), r), -1.0, 3.0);
    check(
        &block,
        &pipe,
        columns(BooleanOperator::Union),
        columns(BooleanOperator::Intersection),
        columns(BooleanOperator::Difference),
    );
}

#[test]
fn a_pipe_through_a_sloped_roof() {
    let column = prism(square(-2.0, -2.0, 2.0, 2.0), 0.0, 10.0);
    let roof = HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, 3.0),
            normal: Vec3::new(-0.4, 0.2, 1.0),
        },
        agreement: false,
    };
    let block = clip_arc_prism_exact(&column, &roof, tol()).expect("a sloped block");
    let (cx, cy, r) = (0.3, 0.2, 0.6);
    let pipe = solid(ArcRing::circle(Point2::new(cx, cy), r), -1.0, 12.0);
    // The block over the pipe's circle averages the roof's height at the
    // circle's centre.
    let mean = 3.0 + 0.4 * cx - 0.2 * cy;
    let under = PI * r * r * mean;
    let (va, vb) = (volume(&block), volume(&pipe));
    check(&block, &pipe, va + vb - under, under, va - under);
}

/// `int_{-h}^{h} sqrt(a^2 - z^2) dz`.
fn chord_integral(a: f64, h: f64) -> f64 {
    let f = |z: f64| 0.5 * z * (a * a - z * z).sqrt() + 0.5 * a * a * (z / a).asin();
    f(h) - f(-h)
}

#[test]
fn a_box_crossed_by_a_ring_on_another_axis() {
    // Not a column: the ring's cylinders run along world y while the box's
    // walls are vertical. The ring is a rectangle x in [4, 6], y in
    // [-1.5, 1.5] about the axis through (-5, 0, 0) along y: distance
    // rho = sqrt((x + 5)^2 + z^2) from 4 to 6, 3 wide. The box spans
    // x in [-2, 2] (x + 5 in [3, 7]), y in [-1, 1], z in [-1, 1].
    use axiolid_construct::revolve_exact::revolve_profile_exact;
    use axiolid_profile::{Profile, RectangleProfile};
    let ring = revolve_profile_exact(
        &Profile::Rectangle(RectangleProfile {
            x: 2.0,
            y: 3.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Point3::new(-5.0, 0.0, 0.0),
        Vec3::Y,
        std::f64::consts::TAU,
        tol(),
    )
    .expect("a ring");
    let block = solid(square(-2.0, -1.0, 2.0, 1.0), -1.0, 1.0);
    // Over z in [-1, 1] the annulus 4 <= rho <= 6 stays within the box's
    // x range, so the intersection is that band of it, 2 deep in y.
    let band = chord_integral(6.0, 1.0) - chord_integral(4.0, 1.0);
    let intersection = 2.0 * band;
    let (va, vb) = (volume(&block), volume(&ring));
    close("ring volume", vb, std::f64::consts::TAU * 5.0 * 6.0);
    check(
        &block,
        &ring,
        va + vb - intersection,
        intersection,
        va - intersection,
    );
}
