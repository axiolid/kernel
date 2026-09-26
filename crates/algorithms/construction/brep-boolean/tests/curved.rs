//! The general boolean on cone and torus faces, and on sections with no
//! closed form (#167 stage 2, ADR 0075 and 0077).
//!
//! Oracles: closed-form volumes where a plane through the axis or across
//! it halves or slices a solid of revolution (Pappus), the identity
//! `|A u B| + |A n B| = |A| + |B|`, and `|A - B| = |A| - |A n B|`, all
//! measured exactly (ADR 0073) on results that audit clean.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_core::{BooleanOperator, Frame2, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

fn tol() -> Tolerance {
    Tolerance::METRE
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

fn close(what: &str, got: f64, want: f64, rel: f64) {
    assert!(
        (got - want).abs() <= rel * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

/// All three operators against the identity, and closed forms where given.
fn check(a: &ExactBRep, b: &ExactBRep, intersection: Option<f64>, rel: f64) -> f64 {
    let u = volume(&run(a, b, BooleanOperator::Union));
    let i = volume(&run(a, b, BooleanOperator::Intersection));
    let d = volume(&run(a, b, BooleanOperator::Difference));
    if let Some(want) = intersection {
        close("intersection", i, want, rel);
    }
    let (va, vb) = (volume(a), volume(b));
    close("identity", u + i, va + vb, rel);
    close("difference", d, va - i, rel);
    i
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn arc(centre: Point2, radius: f64, from: f64, to: f64) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }),
        domain: Interval::new(from, to),
        same_sense: true,
    }
}

fn revolve(segments: Vec<ProfileSegment>) -> ExactBRep {
    let profile = Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    });
    let solid =
        revolve_profile_exact(&profile, Point3::ZERO, Vec3::Y, TAU, tol()).expect("revolves");
    let health = geometric_audit(&solid, tol());
    assert!(health.is_consistent(), "{:?}", health.defects());
    solid
}

/// A torus about the world y axis: tube radius `r` about a circle of
/// radius `big`.
fn ring(big: f64, r: f64) -> ExactBRep {
    let c = Point2::new(big, 0.0);
    let q = 0.5 * PI;
    revolve(vec![
        arc(c, r, 0.0, q),
        arc(c, r, q, PI),
        arc(c, r, PI, 3.0 * q),
        arc(c, r, 3.0 * q, TAU),
    ])
}

/// A box from a column prism: `x`, `y` (world) and `z` ranges.
fn block(x: (f64, f64), y: (f64, f64), z: (f64, f64)) -> ExactBRep {
    let p = ArcPrism {
        section: ArcRing::from_points(&[
            Point2::new(x.0, y.0),
            Point2::new(x.1, y.0),
            Point2::new(x.1, y.1),
            Point2::new(x.0, y.1),
        ]),
        bottom: z.0,
        top: z.1,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a box")
}

#[test]
fn a_box_holding_half_a_torus() {
    // The box's face x = 0 is a plane through the torus's axis (world y):
    // it cuts two tube circles, and the box holds exactly half the torus.
    let (big, r) = (4.0, 1.0);
    let torus = ring(big, r);
    close("torus", volume(&torus), 2.0 * PI * PI * big * r * r, 1e-10);
    let half = block((0.0, 6.0), (-6.0, 6.0), (-6.0, 6.0));
    check(&torus, &half, Some(PI * PI * big * r * r), 1e-9);
}

/// Area of a circle of radius `r` beyond a chord at distance `d`.
fn segment(r: f64, d: f64) -> f64 {
    r * r * (d / r).acos() - d * (r * r - d * d).sqrt()
}

#[test]
fn a_slab_across_the_tube_of_a_torus() {
    // The box above the plane y = 0.3 (world), across the whole ring: by
    // Pappus the torus above it is the tube segment beyond the chord at
    // 0.3, swept round at its centroid radius, which is the ring's own.
    let (big, r) = (4.0, 1.0);
    let torus = ring(big, r);
    let slab = block((-6.0, 6.0), (0.3, 6.0), (-6.0, 6.0));
    check(&torus, &slab, Some(TAU * big * segment(r, 0.3)), 1e-9);
}

#[test]
fn half_a_revolved_cone_ring() {
    // A trapezoid about the y axis: a cone wall outside, a cylinder inside.
    let solid = revolve(vec![
        line(Point2::new(2.0, 0.0), Point2::new(4.0, 0.0)),
        line(Point2::new(4.0, 0.0), Point2::new(3.0, 2.0)),
        line(Point2::new(3.0, 2.0), Point2::new(2.0, 2.0)),
        line(Point2::new(2.0, 2.0), Point2::new(2.0, 0.0)),
    ]);
    let v = volume(&solid);
    let half = block((0.0, 6.0), (-6.0, 6.0), (-6.0, 6.0));
    check(&solid, &half, Some(0.5 * v), 1e-9);
    // And a slab across the cone: the ring between heights 0.5 and 2.
    let slab = block((-6.0, 6.0), (0.5, 6.0), (-6.0, 6.0));
    // Frustum above y = 0.5: outer radius 3.75 -> 3, inner 2, height 1.5.
    let frustum = PI * 1.5 / 3.0 * (3.75 * 3.75 + 3.75 * 3.0 + 3.0 * 3.0);
    check(&solid, &slab, Some(frustum - PI * 4.0 * 1.5), 1e-9);
}

#[test]
fn a_pipe_through_the_tube_of_a_torus() {
    // A vertical pipe (world z) through the ring where the ring itself runs
    // along z: the torus/cylinder section has no closed form and is traced
    // (ADR 0077).
    let (big, r) = (4.0, 1.0);
    let torus = ring(big, r);
    let p = ArcPrism {
        section: ArcRing::circle(Point2::new(big, 0.1), 0.5),
        bottom: -8.0,
        top: 8.0,
    };
    let pipe = boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).unwrap();
    let i = check(&torus, &pipe, None, 1e-8);
    assert!(i > 0.0 && i < volume(&pipe));
}
