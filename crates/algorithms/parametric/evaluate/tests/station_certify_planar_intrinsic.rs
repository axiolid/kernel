//! Certified joint distances on and beside a planar 3D intrinsic curve
//! (#302): a clothoid stored as a `Curve3::Intrinsic` with no torsion in a
//! level plane.
//!
//! The reference here is independent of the kernel's closed-form speed:
//! the curve's heading and position are integrated by classical RK4 from
//! its curvature law, the offset point is built from them in components
//! (`Z x T` for the lateral, `+Z` for up), and the offset's length is the
//! sum of its chords at two step sizes, Richardson-extrapolated. That
//! agrees with the closed form `L - d theta` of a constant offset to well
//! inside `1e-12`, which the first test checks.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{
    CurvatureLaw, Curve3, CurvePath, Intrinsic3, OffsetFrame, OffsetLaw, PathCurve, PathOffset,
    PathOffsets, PathPiece,
};
use axiolid_curve_evaluate_contract::{CurveEvaluator, DistanceInterval};
use axiolid_evaluate::station::{
    offset_pieces, CompositeBasis, StationCurve, StationPiece, UNCERTIFIED_LENGTH,
};
use axiolid_evaluate::ReferenceCurveEvaluator;

// --- the curve ------------------------------------------------------------------

/// The curvature at `s` of [`clothoid_law`]: straight to `1 / 300` over
/// 100 m, then the arc it reaches.
fn kappa(s: Scalar) -> Scalar {
    if s <= 100.0 {
        s / 30_000.0
    } else {
        1.0 / 300.0
    }
}

fn clothoid_law() -> CurvatureLaw {
    CurvatureLaw::piecewise(
        vec![100.0],
        vec![
            CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 100.0),
            CurvatureLaw::Constant {
                curvature: 1.0 / 300.0,
            },
        ],
    )
}

/// A level start frame at height 2: turning left (`y = +Y`) or right (`y =
/// -Y`, the binormal `-Z`).
fn level(left: bool) -> Frame3 {
    let side = if left { 1.0 } else { -1.0 };
    Frame3 {
        origin: Point3::new(5.0, -3.0, 2.0),
        x: Vec3::X,
        y: Vec3::Y * side,
        z: Vec3::Z * side,
    }
}

fn spiral(start: Frame3, torsion: CurvatureLaw) -> Curve3 {
    Curve3::Intrinsic(Intrinsic3::new(start, clothoid_law(), torsion, 150.0))
}

fn planar(left: bool) -> Curve3 {
    spiral(level(left), CurvatureLaw::Constant { curvature: 0.0 })
}

// --- an independent reference ---------------------------------------------------

/// The displacement of the offset point from the curve's, given its unit
/// tangent `t` and the law `(lateral, vertical)` there.
fn displacement(t: Vec3, lateral: Scalar, vertical: Scalar) -> Vec3 {
    let left = Vec3::Z.cross(t).normalize();
    left * lateral + Vec3::Z * vertical
}

/// The length over `[a, b]` of the offset by `law(s) = (lateral,
/// vertical)` beside the curve from `start`: RK4 on `theta' = k`, `P' =
/// cos(theta) x + sin(theta) y`, its chords summed at `n` and `2 n` steps
/// and Richardson-extrapolated. Each chord is the curve's RK4 increment
/// plus the change of the displacement, so no rounding of the position's
/// magnitude enters it, and the chords are summed compensated.
fn reference(
    start: Frame3,
    a: Scalar,
    b: Scalar,
    law: &dyn Fn(Scalar) -> (Scalar, Scalar),
) -> Scalar {
    let chords = |n: usize| -> Scalar {
        let h = 150.0 / n as Scalar;
        let tangent = |theta: Scalar| start.x * theta.cos() + start.y * theta.sin();
        let derivative = |s: Scalar, theta: Scalar| (kappa(s), tangent(theta));
        let inside = |s: Scalar| s >= a - 1e-9 && s <= b + 1e-9;
        let mut theta: Scalar = 0.0;
        let (mut total, mut carry) = (0.0, 0.0);
        for i in 0..n {
            let s = h * i as Scalar;
            let (k1, v1) = derivative(s, theta);
            let (k2, v2) = derivative(s + 0.5 * h, theta + 0.5 * h * k1);
            let (k3, v3) = derivative(s + 0.5 * h, theta + 0.5 * h * k2);
            let (k4, v4) = derivative(s + h, theta + h * k3);
            let next = theta + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            let step = (v1 + v2 * 2.0 + v3 * 2.0 + v4) * (h / 6.0);
            let end = h * (i + 1) as Scalar;
            if inside(s) && inside(end) {
                let (l0, w0) = law(s);
                let (l1, w1) = law(end);
                let chord = step + displacement(tangent(next), l1, w1)
                    - displacement(tangent(theta), l0, w0);
                let term = chord.length() - carry;
                let sum = total + term;
                carry = (sum - total) - term;
                total = sum;
            }
            theta = next;
        }
        total
    };
    // Steps that land on 100 and 150 exactly.
    let (coarse, fine) = (chords(3 * (1 << 14)), chords(3 * (1 << 15)));
    (4.0 * fine - coarse) / 3.0
}

fn certifies(interval: DistanceInterval, truth: Scalar, tolerance: Scalar, what: &str) {
    assert!(
        interval.lower <= interval.upper && interval.width() <= tolerance,
        "{what}: {interval:?} is wider than {tolerance:e}"
    );
    assert!(
        interval.contains(truth),
        "{what}: {interval:?} does not contain {truth} (off by {:e})",
        if truth < interval.lower {
            interval.lower - truth
        } else {
            truth - interval.upper
        }
    );
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, input: &str, what: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input: found, .. }) => assert_eq!(found, input, "{what}"),
        other => panic!("{what}: expected a refusal naming {input:?}, got {other:?}"),
    }
}

/// An offset by `distance` to the left in plan: a 3D offset needs its
/// reference direction, here `+Z`.
fn upright(distance: Scalar) -> OffsetLaw {
    OffsetLaw::Directed {
        distance,
        reference_direction: Vec3::Z,
    }
}

/// The whole of an offset of `base` over `[a, b]`.
fn offset_path_piece(base: &Curve3, a: Scalar, b: Scalar, law: OffsetLaw) -> PathPiece {
    PathPiece::whole(PathCurve::Offset(Box::new(PathOffset::new(
        PathPiece::new(PathCurve::Three(base.clone()), a, b),
        law,
    ))))
}

// --- the tests --------------------------------------------------------------------

#[test]
fn the_reference_agrees_with_the_closed_form() {
    // Beside a clothoid turning by 100 / 600 rad at `d` the offset runs `1 -
    // d k` times as fast: its length is `100 - d / 6`.
    for (left, d) in [(true, 5.0), (true, -5.0), (false, 5.0)] {
        let sign = if left { 1.0 } else { -1.0 };
        let found = reference(level(left), 0.0, 100.0, &|_| (d, 0.0));
        let truth = 100.0 - sign * d / 6.0;
        assert!((found - truth).abs() <= 1e-12, "{found} vs {truth}");
    }
}

#[test]
fn a_planar_3d_clothoid_measures_its_joints_exactly() {
    let curve = planar(true);
    let pieces = vec![
        StationPiece::between(StationCurve::Three(&curve), 0.0, 100.0).unwrap(),
        StationPiece::between(StationCurve::Three(&curve), 100.0, 150.0).unwrap(),
    ];
    let basis = CompositeBasis::new(pieces).unwrap();
    let joints = basis.joint_intervals(1e-9).unwrap();
    assert_eq!(joints, vec![DistanceInterval::point(100.0)]);
    let path = CurvePath::new(vec![
        PathPiece::new(PathCurve::Three(curve.clone()), 0.0, 100.0),
        PathPiece::new(PathCurve::Three(curve.clone()), 100.0, 150.0),
    ]);
    let joints = ReferenceCurveEvaluator::new()
        .path_joint_distances(&path, 1e-9)
        .unwrap();
    certifies(joints[0], 100.0, 1e-9, "path");
}

#[test]
fn a_constant_offset_of_a_planar_3d_clothoid_is_certified() {
    for left in [true, false] {
        let curve = planar(left);
        let whole = StationPiece::whole(StationCurve::Three(&curve)).unwrap();
        for distance in [5.0, -5.0] {
            let law = upright(distance);
            let truth = reference(level(left), 0.0, 100.0, &|_| (distance, 0.0));
            let basis = CompositeBasis::new(offset_pieces(&[whole], law).unwrap()).unwrap();
            let path = CurvePath::new(vec![
                offset_path_piece(&curve, 0.0, 100.0, law),
                offset_path_piece(&curve, 100.0, 150.0, law),
            ]);
            let what = format!("left {left}, d {distance}");
            for tolerance in [1e-6, 1e-9] {
                certifies(
                    basis.joint_intervals(tolerance).unwrap()[0],
                    truth,
                    tolerance,
                    &what,
                );
                let joints = ReferenceCurveEvaluator::new()
                    .path_joint_distances(&path, tolerance)
                    .unwrap();
                certifies(joints[0], truth, tolerance, &what);
            }
        }
    }
}

#[test]
fn an_offset_by_distances_of_a_planar_3d_clothoid_is_certified() {
    // Widening from 1 to 3 and rising from 0 to 2 over the clothoid, then
    // held beside the arc: measured in 3D, so the rise lengthens it.
    let widening = OffsetLaw::Linear {
        start: PathOffsets::new(1.0, 0.0, 0.0),
        end: PathOffsets::new(3.0, 2.0, 0.0),
        frame: OffsetFrame::Section,
    };
    let held = OffsetLaw::Linear {
        start: PathOffsets::new(3.0, 2.0, 0.0),
        end: PathOffsets::new(3.0, 2.0, 0.0),
        frame: OffsetFrame::Section,
    };
    for left in [true, false] {
        let curve = planar(left);
        let first = StationPiece::between(StationCurve::Three(&curve), 0.0, 100.0).unwrap();
        let second = StationPiece::between(StationCurve::Three(&curve), 100.0, 150.0).unwrap();
        let mut pieces = offset_pieces(&[first], widening).unwrap();
        pieces.extend(offset_pieces(&[second], held).unwrap());
        let basis = CompositeBasis::new(pieces).unwrap();
        let truth = reference(level(left), 0.0, 100.0, &|s| {
            (1.0 + 2.0 * s / 100.0, 2.0 * s / 100.0)
        });
        // The rise lengthens it: about `100 (sqrt(1 + 0.02^2) - 1)`.
        let flat = reference(level(left), 0.0, 100.0, &|s| (1.0 + 2.0 * s / 100.0, 0.0));
        assert!(truth - flat > 0.015);
        let path = CurvePath::new(vec![
            offset_path_piece(&curve, 0.0, 100.0, widening),
            offset_path_piece(&curve, 100.0, 150.0, held),
        ]);
        for tolerance in [1e-6, 1e-9] {
            let what = format!("left {left}");
            certifies(
                basis.joint_intervals(tolerance).unwrap()[0],
                truth,
                tolerance,
                &what,
            );
            let joints = ReferenceCurveEvaluator::new()
                .path_joint_distances(&path, tolerance)
                .unwrap();
            certifies(joints[0], truth, tolerance, &what);
        }
    }
}

#[test]
fn a_twisted_or_tilted_intrinsic_is_refused_by_name() {
    let law = upright(5.0);
    // A torsion law not identically zero: a space curve.
    let twisted = spiral(level(true), CurvatureLaw::Constant { curvature: 1e-3 });
    // No torsion, in a plane tilted about the start tangent (rolled) or
    // about the start normal (pitched): the lateral of its reference-up
    // frame leaves that plane. By a large angle, and by `1e-7`, a tilt
    // the frame's own orthonormality in plan does not show (its axes'
    // plan lengths are `1 - 5e-15`).
    let flat = || CurvatureLaw::Constant { curvature: 0.0 };
    let tilted = |x: Vec3, y: Vec3| {
        spiral(
            Frame3 {
                origin: Point3::new(5.0, -3.0, 2.0),
                x,
                y,
                z: x.cross(y),
            },
            flat(),
        )
    };
    let mut planes = Vec::new();
    for (angle, roll, pitch) in [
        (0.9, "rolled", "pitched"),
        (1e-7, "rolled slightly", "pitched slightly"),
    ] {
        let (c, s): (Scalar, Scalar) = (Scalar::cos(angle), Scalar::sin(angle));
        planes.push((roll, tilted(Vec3::X, Vec3::new(0.0, c, s))));
        planes.push((pitch, tilted(Vec3::new(c, 0.0, s), Vec3::Y)));
    }
    // A torsion law that is zero only once evaluated (a frozen sine that
    // cancels its mean) is not read as planar either.
    let frozen = spiral(
        level(true),
        CurvatureLaw::Sinusoid {
            mean: -1.0,
            amplitude: 1.0,
            angular_frequency: 0.0,
            phase: core::f64::consts::FRAC_PI_2,
        },
    );
    planes.push(("twisted", twisted));
    planes.push(("frozen", frozen));
    for (what, curve) in planes {
        let whole = StationPiece::whole(StationCurve::Three(&curve)).unwrap();
        let basis = CompositeBasis::new(offset_pieces(&[whole], law).unwrap()).unwrap();
        refused(basis.joint_intervals(1e-6), UNCERTIFIED_LENGTH, what);
        let path = CurvePath::new(vec![
            offset_path_piece(&curve, 0.0, 100.0, law),
            offset_path_piece(&curve, 100.0, 150.0, law),
        ]);
        refused(
            ReferenceCurveEvaluator::new().path_joint_distances(&path, 1e-6),
            UNCERTIFIED_LENGTH,
            what,
        );
    }
}
