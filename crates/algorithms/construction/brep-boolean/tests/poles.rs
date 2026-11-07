//! Booleans on faces that reach a pole or apex (#379).
//!
//! A cone's apex and a sphere's pole are single vertices whose wall loop
//! runs `seam, [rim], seam reversed`: in parameters it leaves a gap along the
//! pole line, and the loop's u jumps by an arbitrary amount there. That jump
//! is not a turn round the surface, so the boolean must not give the face a
//! second seam for it. Each volume is a closed form from the inputs.

mod common;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_core::{BooleanOperator, Interval, Point2, Point3, Tolerance, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

fn audited(brep: &ExactBRep, what: &str) {
    let health = geometric_audit(brep, tol());
    assert!(health.is_consistent(), "{what}: {:?}", health.defects());
    let topology = axiolid_topology::audit_brep(brep.topology());
    assert!(topology.is_closed_manifold(), "{what}: {topology:?}");
}

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

fn polygon(points: &[(f64, f64)]) -> Profile {
    let points: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    let segments = (0..points.len())
        .map(|i| {
            let (from, to) = (points[i], points[(i + 1) % points.len()]);
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: from,
                    direction: to - from,
                }),
                domain: Interval::UNIT,
                same_sense: true,
            }
        })
        .collect();
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

fn revolved(profile: &Profile, angle: f64) -> ExactBRep {
    let solid = revolve_profile_exact(profile, Point3::ZERO, Vec3::Y, angle, tol())
        .expect("a revolution touching the axis");
    audited(&solid, "operand");
    solid
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-10 * want.abs(),
        "{what}: expected {want}, got {got}"
    );
}

fn cut(a: &ExactBRep, b: &ExactBRep, want: f64, what: &str) {
    let result = boolean(a, b, BooleanOperator::Difference, tol())
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    audited(&result, what);
    close(what, volume(&result), want);
}

/// The slab `y >= 1` and an off-axis block, cut from a three-quarter cone.
#[test]
fn a_three_quarter_cone_loses_a_slab_and_a_block() {
    let cone = revolved(&polygon(&[(0.0, 0.0), (1.5, 0.0), (0.0, 2.0)]), 1.5 * PI);
    close("operand", volume(&cone), 0.75 * PI * 1.5 * 1.5 * 2.0 / 3.0);
    cut(
        &cone,
        &block((-6.0, 6.0), (1.0, 6.0), (-6.0, 6.0)),
        3.0925052683774528,
        "slab",
    );
    cut(
        &cone,
        &block((0.2, 0.7), (0.1, 0.6), (-0.25, 0.3)),
        3.4717917352885173,
        "block",
    );
}
