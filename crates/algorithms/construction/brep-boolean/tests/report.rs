//! Exact booleans at `Tolerance::ZERO`, and the report of what a boolean
//! read within tolerance (#236).
//!
//! A wall minus an opening across it: through, flush (its caps in the
//! wall's faces) and blind. Placed by matrices whose entries are exactly
//! `0` and `+-1`, the operands' faces are exactly coplanar or perpendicular,
//! the exact predicates decide every contact, and the difference succeeds
//! at `Tolerance::ZERO` with nothing read within tolerance. Built from
//! sines and cosines under a general placement, the flush caps miss the
//! wall's faces by rounding: read as coplanar at a positive tolerance, and
//! reported. At zero tolerance nothing is read under any placement (#251).

use std::f64::consts::FRAC_PI_2;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::{
    boolean, boolean_with_report, BooleanError, BooleanReport, ToleranceDecisionKind,
    ROUNDING_FACTOR,
};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{BooleanOperator, Mat3, Tolerance, Transform3, Vec3};
use axiolid_measure::exact_properties;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_surface::Surface;

/// Wall length, thickness and height.
const L: f64 = 6.0;
const T: f64 = 0.25;
const H: f64 = 3.0;
/// Opening width and height.
const W: f64 = 1.0;
const OH: f64 = 1.5;

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// A quarter turn about `x`, every entry exactly `0` or `+-1`.
fn exact_across() -> Transform3 {
    Transform3::from_mat3(Mat3::from_cols(Vec3::X, Vec3::Z, -Vec3::Y))
}

/// The same turn from sines and cosines: `cos(pi / 2)` leaves `6e-17`.
fn rounded_across() -> Transform3 {
    Transform3::from_rotation_x(FRAC_PI_2)
}

/// Placements with exact entries: none, a quarter turn about `z`, and a
/// half turn about `x` (the wall upside down), each with an offset.
fn exact_placements() -> [Transform3; 3] {
    [
        Transform3::IDENTITY,
        Transform3::from_mat3_translation(
            Mat3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z),
            Vec3::new(12.5, -4.0, 3.25),
        ),
        Transform3::from_mat3_translation(
            Mat3::from_cols(Vec3::X, -Vec3::Y, -Vec3::Z),
            Vec3::new(-1.5, 0.75, 9.0),
        ),
    ]
}

/// A rigid placement with no axis left aligned.
fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

fn wall(p: Transform3) -> ExactBRep {
    extrude_profile_exact(&rect(L, T), Vec3::Z, H, Tolerance::ZERO)
        .expect("a wall")
        .transformed(&p)
        .expect("rigid")
}

/// An opening `W x OH` centred at wall `(0.5, 1.25)` in `x, z`, from wall
/// `y = start` along `-y` for `depth`.
fn opening(depth: f64, start: f64, across: Transform3, p: Transform3) -> ExactBRep {
    extrude_profile_exact(&rect(W, OH), Vec3::Z, depth, Tolerance::ZERO)
        .expect("an opening")
        .transformed(&(p * Transform3::from_translation(Vec3::new(0.5, start, 1.25)) * across))
        .expect("rigid")
}

/// Through, flush and blind (to the wall's middle plane): depth, start and
/// the volume removed.
fn cases() -> [(&'static str, f64, f64, f64); 3] {
    [
        ("through", T + 0.5, T / 2.0 + 0.25, W * OH * T),
        ("flush", T, T / 2.0, W * OH * T),
        ("blind", 0.25 + T / 2.0, T / 2.0 + 0.25, W * OH * T / 2.0),
    ]
}

/// The difference with its report, audited and measured against `removed`
/// within `slack`.
fn cut(
    a: &ExactBRep,
    b: &ExactBRep,
    tolerance: Tolerance,
    removed: f64,
    slack: f64,
) -> BooleanReport {
    let (result, report) =
        boolean_with_report(a, b, BooleanOperator::Difference, tolerance).expect("difference");
    let health = geometric_audit(&result, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    let expected = L * T * H - removed;
    let measured = exact_properties(&result, Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    assert!(
        (measured - expected).abs() <= slack + 1e-12 * expected,
        "volume {measured}, expected {expected}"
    );
    report
}

#[test]
fn exact_axis_openings_are_exact_at_zero_tolerance() {
    for p in exact_placements() {
        for (name, depth, start, removed) in cases() {
            let report = cut(
                &wall(p),
                &opening(depth, start, exact_across(), p),
                Tolerance::ZERO,
                removed,
                0.0,
            );
            assert!(report.is_exact(), "{name}: {report:?}");
            assert_eq!((report.linear(), report.angular()), (0.0, 0.0));
        }
    }
}

#[test]
fn exact_axis_openings_read_nothing_at_a_positive_tolerance() {
    // The exact predicates come first: exactly coplanar caps are exact
    // whatever the tolerance.
    for p in exact_placements() {
        for (name, depth, start, removed) in cases() {
            let report = cut(
                &wall(p),
                &opening(depth, start, exact_across(), p),
                Tolerance::METRE,
                removed,
                0.0,
            );
            assert!(report.is_exact(), "{name}: {report:?}");
        }
    }
}

#[test]
fn a_rounded_flush_opening_reads_its_caps_within_tolerance() {
    let eps = Tolerance::METRE.linear();
    for p in [Transform3::IDENTITY, general()] {
        let (_, depth, start, removed) = cases()[1];
        let report = cut(
            &wall(p),
            &opening(depth, start, rounded_across(), p),
            Tolerance::METRE,
            removed,
            W * OH * eps,
        );
        assert!(report.contains(ToleranceDecisionKind::CoincidentSupports));
        assert!(
            report.linear() > 0.0 && report.linear() <= eps,
            "{report:?}"
        );
        assert!(report.angular() <= Tolerance::METRE.angular());
        // At zero tolerance nothing may be read: refused, never guessed.
        boolean(
            &wall(p),
            &opening(depth, start, rounded_across(), p),
            BooleanOperator::Difference,
            Tolerance::ZERO,
        )
        .expect_err("rounded caps cannot be read as coplanar at zero tolerance");
    }
}

#[test]
fn a_rounded_through_opening_crosses_the_wall_exactly() {
    // Crossing faces need no reading: under a general rotation a through
    // opening meets the wall's faces transversally, and the result is the
    // exact difference of the operands as given, at either tolerance.
    let p = general();
    let (_, depth, start, removed) = cases()[0];
    for tolerance in [Tolerance::ZERO, Tolerance::METRE] {
        let report = cut(
            &wall(p),
            &opening(depth, start, rounded_across(), p),
            tolerance,
            removed,
            1e-12,
        );
        assert!(report.is_exact(), "{report:?}");
    }
}

#[test]
fn exactly_parallel_caps_a_fraction_of_the_tolerance_apart_are_read_and_reported() {
    // Exact axes, but the opening stops a tenth of the tolerance short of
    // the far face: exactly parallel, not coplanar. At a positive tolerance
    // the cap is read as flush, and reported with its offset; at zero the
    // skin is kept, exactly.
    let eps = Tolerance::METRE.linear();
    let skin = 0.1 * eps;
    let (depth, start) = (T - skin, T / 2.0);
    let removed = W * OH * T;
    let report = cut(
        &wall(Transform3::IDENTITY),
        &opening(depth, start, exact_across(), Transform3::IDENTITY),
        Tolerance::METRE,
        removed,
        W * OH * eps,
    );
    assert!(report.contains(ToleranceDecisionKind::CoincidentSupports));
    assert!((report.linear() - skin).abs() <= 1e-3 * skin, "{report:?}");
    let report = cut(
        &wall(Transform3::IDENTITY),
        &opening(depth, start, exact_across(), Transform3::IDENTITY),
        Tolerance::ZERO,
        W * OH * (T - skin),
        1e-15,
    );
    assert!(report.is_exact(), "{report:?}");
}

#[test]
fn boolean_is_the_result_of_boolean_with_report() {
    let p = general();
    let (_, depth, start, _) = cases()[1];
    let (a, b) = (wall(p), opening(depth, start, rounded_across(), p));
    let plain = boolean(&a, &b, BooleanOperator::Difference, Tolerance::METRE).expect("cut");
    let (reported, _) =
        boolean_with_report(&a, &b, BooleanOperator::Difference, Tolerance::METRE).expect("cut");
    assert_eq!(plain, reported);
}

#[test]
fn the_report_carries_the_rounding_floor_of_its_operands() {
    // #244. The wall `[-3, 3] x [-1/8, 1/8] x [0, 3]` and an opening inside
    // it: the extent is the wall's largest coordinate, `L / 2`. Placed by
    // the quarter turn offset by `(12.5, -4, 3.25)`, the wall's `y = -1/8`
    // face lands at `x = 12.625`. An exact result at zero tolerance still
    // exposes its floor.
    let placements = exact_placements();
    for (p, extent) in [(placements[0], L / 2.0), (placements[1], 12.5 + T / 2.0)] {
        for tolerance in [Tolerance::ZERO, Tolerance::METRE] {
            let (_, depth, start, removed) = cases()[1];
            let report = cut(
                &wall(p),
                &opening(depth, start, exact_across(), p),
                tolerance,
                removed,
                0.0,
            );
            assert!(report.is_exact(), "{report:?}");
            assert_eq!(report.extent(), extent, "{report:?}");
            assert_eq!(report.rounding_floor(), ROUNDING_FACTOR * extent);
        }
    }
    // The larger operand sets the extent: an opening reaching `y = -10`.
    let (_, _, _, removed) = cases()[0];
    let report = cut(
        &wall(Transform3::IDENTITY),
        &opening(20.0, 10.0, exact_across(), Transform3::IDENTITY),
        Tolerance::ZERO,
        removed,
        0.0,
    );
    assert_eq!(report.extent(), 10.0, "{report:?}");
    assert_eq!(report.rounding_floor(), 10.0 * ROUNDING_FACTOR);
}

const R_HOLE: f64 = 0.5;

/// A round hole of radius `R_HOLE` across the wall at `(0.5, 1.25)` in
/// `x, z`, from wall `y = start` along `-y` for `depth`.
fn round_hole(depth: f64, start: f64, p: Transform3) -> ExactBRep {
    let profile = Profile::Circle(CircleProfile {
        radius: R_HOLE,
        thickness: None,
    });
    extrude_profile_exact(&profile, Vec3::Z, depth, Tolerance::ZERO)
        .expect("a hole")
        .transformed(
            &(p * Transform3::from_translation(Vec3::new(0.5, start, 1.25)) * exact_across()),
        )
        .expect("rigid")
}

#[test]
fn a_round_hole_under_a_general_placement_is_exact_at_zero_tolerance() {
    // #251. Placed by one general rotation, the hole's axis and the wall's
    // faces' normals are the same `f64` vector: their `f64` cross product
    // is exactly zero, while the exact numbers (a frame orthonormal only to
    // rounding) are not perpendicular. That reading used to be taken at
    // zero tolerance and reported with `linear` 0. Now the exact answer
    // stands: the general closed form cuts the hole, the report is empty.
    let p = general();
    let (w, h) = (wall(p), round_hole(T + 0.5, T / 2.0 + 0.25, p));
    let axis = h
        .surfaces()
        .iter()
        .find_map(|s| match s {
            Surface::Cylinder(c) => Some(c.frame.z.normalize()),
            _ => None,
        })
        .expect("a cylinder");
    let zero_cross = w.surfaces().iter().any(|s| match s {
        Surface::Plane(q) => q.frame.z.normalize().cross(axis).length() == 0.0,
        _ => false,
    });
    assert!(zero_cross, "the f64 reading must be perpendicular");
    let disc = std::f64::consts::PI * R_HOLE * R_HOLE;
    for (depth, start, removed) in [
        (T + 0.5, T / 2.0 + 0.25, disc * T),
        (0.25 + T / 2.0, T / 2.0 + 0.25, disc * T / 2.0),
    ] {
        let hole = round_hole(depth, start, p);
        let report = cut(&w, &hole, Tolerance::ZERO, removed, 1e-9);
        assert!(report.is_exact(), "depth {depth}: {report:?}");
        assert_eq!((report.linear(), report.angular()), (0.0, 0.0));
        // At a positive tolerance the reading is still taken, and reported.
        let report = cut(&w, &hole, Tolerance::METRE, removed, 1e-9);
        assert!(
            report.contains(ToleranceDecisionKind::PlanePerpendicularToAxis),
            "{report:?}"
        );
    }
}

/// A deterministic stream of numbers in `[0, 1)` (SplitMix64).
struct Stream(u64);

impl Stream {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A random rigid placement: a turn about a random axis, offset up to
    /// twenty units.
    fn placement(&mut self) -> Transform3 {
        let mut c = || 2.0 * self.next() - 1.0;
        let axis = Vec3::new(c(), c(), c()).normalize_or(Vec3::Z);
        let angle = std::f64::consts::PI * c();
        let offset = Vec3::new(20.0 * c(), 20.0 * c(), 20.0 * c());
        Transform3::from_translation(offset) * Transform3::from_axis_angle(axis, angle)
    }
}

#[test]
fn reports_at_zero_tolerance_are_empty_under_random_general_placements() {
    // #236, #251: at zero tolerance a boolean is exact or refused by name,
    // whatever the placement. Rectangular openings (exact and rounded axis
    // matrices) and round holes, through, flush and blind, under random
    // rigid placements, by every operator.
    //
    // Through and blind tools cross the wall's faces transversally (the
    // round hole's cylinder meets them perpendicularly in `f64`, not in the
    // numbers given) and always succeed; a flush tool's caps miss the
    // wall's faces by the placement's rounding and may only be refused.
    let mut stream = Stream(0x251);
    let mut exact = 0;
    for _ in 0..24 {
        let p = stream.placement();
        let w = wall(p);
        for (name, depth, start, _) in cases() {
            for tool in [
                opening(depth, start, exact_across(), p),
                opening(depth, start, rounded_across(), p),
                round_hole(depth, start, p),
            ] {
                for operator in [
                    BooleanOperator::Difference,
                    BooleanOperator::Union,
                    BooleanOperator::Intersection,
                ] {
                    match boolean_with_report(&w, &tool, operator, Tolerance::ZERO) {
                        Ok((_, report)) => {
                            assert!(report.is_exact(), "{name} {p:?} {operator:?}: {report:?}");
                            exact += 1;
                        }
                        Err(e) => {
                            assert_ne!(e, BooleanError::ToleranceExceeded);
                            assert_eq!(name, "flush", "{p:?} {operator:?}: {e}");
                        }
                    }
                }
            }
        }
    }
    assert_eq!(exact, 24 * 2 * 3 * 3);
}
