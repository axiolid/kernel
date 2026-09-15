//! Space curves from curvature and torsion (ADR 0061).
//!
//! The helix is checked against its CLOSED FORM, and the planar case against
//! the 2D path that was itself pinned to Fresnel. Both are independent of the
//! integrator: a scheme checked only against another run of itself proves
//! nothing.

use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{CurvatureLaw, Intrinsic2, Intrinsic3};
use axiolid_evaluate::{frenet_frame, frenet_point, frenet_tangent, intrinsic_point};

fn start_frame() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

/// A helix of radius `a` and pitch parameter `b` has constant curvature and
/// torsion, so it is the exactly-solvable case.
fn helix_laws(a: Scalar, b: Scalar) -> (Scalar, Scalar) {
    let c = a.hypot(b);
    (a / (c * c), b / (c * c))
}

#[test]
fn a_helix_matches_its_closed_form() {
    // The helix axis is the DARBOUX vector (tau, 0, k)/|omega|, not +z: the
    // binormal starts along +z but precesses. With constant omega the motion
    // is a screw about that axis, so
    //   p(s) = u (u.e1) s + (e1 - u (u.e1)) sin(w s)/w + (u x e1)(1 - cos(w s))/w
    // is the exact closed form, independent of the integrator.
    let (a, b): (Scalar, Scalar) = (3.0, 1.5);
    let c = a.hypot(b);
    let (k, tau) = helix_laws(a, b);
    let length = 2.0 * core::f64::consts::TAU * c; // two full turns
    let curve = Intrinsic3::new(
        start_frame(),
        CurvatureLaw::circular(k),
        CurvatureLaw::circular(tau),
        length,
    );

    let omega = Vec3::new(tau, 0.0, k);
    let rate = omega.length();
    let axis = omega / rate;
    let along = axis.dot(Vec3::X);

    for fraction in [0.0, 0.125, 0.25, 0.5, 0.75, 1.0] {
        let s = fraction * length;
        let turn = rate * s;
        let want = axis * (along * s)
            + (Vec3::X - axis * along) * (turn.sin() / rate)
            + axis.cross(Vec3::X) * ((1.0 - turn.cos()) / rate);
        let want = Point3::new(want.x, want.y, want.z);
        let got = frenet_point(&curve, s).expect("a helix evaluates");
        let error = (got - want).length();
        assert!(
            error < 1e-9,
            "helix at s={s}: got {got:?} want {want:?} error {error:e}"
        );
    }
}

#[test]
fn a_helix_tangent_matches_its_closed_form() {
    let (a, b): (Scalar, Scalar) = (2.0, 0.75);
    let c = a.hypot(b);
    let (k, tau) = helix_laws(a, b);
    let length = core::f64::consts::TAU * c;
    let curve = Intrinsic3::new(
        start_frame(),
        CurvatureLaw::circular(k),
        CurvatureLaw::circular(tau),
        length,
    );

    let s = 0.4 * length;
    // Rodrigues: the tangent is e1 rotated about the Darboux axis by |omega| s.
    let omega = Vec3::new(tau, 0.0, k);
    let rate = omega.length();
    let axis = omega / rate;
    let turn = rate * s;
    let want = Vec3::X * turn.cos()
        + axis.cross(Vec3::X) * turn.sin()
        + axis * (axis.dot(Vec3::X)) * (1.0 - turn.cos());
    let got = frenet_tangent(&curve, s).expect("a helix has a tangent");
    assert!(
        (got - want).length() < 1e-9,
        "tangent: got {got:?} want {want:?}"
    );
}

#[test]
fn zero_torsion_reproduces_the_planar_curve_exactly() {
    // The whole point of the 3D path is that it must not disagree with the
    // 2D one where both apply. The 2D clothoid was pinned against Fresnel.
    let (radius, length) = (300.0, 120.0);
    let law = CurvatureLaw::clothoid(0.0, 1.0 / radius, length);
    let spatial = Intrinsic3::new(start_frame(), law.clone(), CurvatureLaw::straight(), length);
    let planar = Intrinsic2::new(
        Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        law,
        length,
    );

    for fraction in [0.25, 0.5, 1.0] {
        let s = fraction * length;
        let flat = intrinsic_point(&planar, s).expect("planar evaluates");
        let space = frenet_point(&spatial, s).expect("spatial evaluates");
        assert!(
            (space.x - flat.x).abs() < 1e-9 && (space.y - flat.y).abs() < 1e-9,
            "at s={s}: 3D ({}, {}) vs 2D ({}, {})",
            space.x,
            space.y,
            flat.x,
            flat.y
        );
        assert!(
            space.z.abs() < 1e-12,
            "zero torsion must stay in the start plane, got z={}",
            space.z
        );
    }
}

#[test]
fn the_frame_stays_orthonormal_on_a_non_commuting_law() {
    // tau/k is NOT constant here, so the generators at different arc lengths
    // do not commute and exp(int Omega) is not the solution. This is the case
    // that separates a Lie-group integrator from a naive one.
    let curve = Intrinsic3::new(
        start_frame(),
        CurvatureLaw::clothoid(0.01, 0.09, 40.0),
        CurvatureLaw::Sinusoid {
            mean: 0.005,
            amplitude: 0.02,
            angular_frequency: 0.3,
            phase: 0.0,
        },
        40.0,
    );

    let frame = frenet_frame(&curve, 40.0).expect("a varying law evaluates");
    let dots = [
        frame.x.dot(frame.y),
        frame.y.dot(frame.z),
        frame.z.dot(frame.x),
    ];
    for dot in dots {
        assert!(
            dot.abs() < 1e-12,
            "frame axes must stay perpendicular: {dot:e}"
        );
    }
    for axis in [frame.x, frame.y, frame.z] {
        assert!(
            (axis.length() - 1.0).abs() < 1e-12,
            "frame axes must stay unit: {}",
            axis.length()
        );
    }
    // Right-handed, and not merely by sign: T x N must BE B, which also
    // pins the frame ordering rather than just its chirality.
    let cross = frame.x.cross(frame.y);
    assert!(
        (cross - frame.z).length() < 1e-12,
        "T x N must equal B, got {cross:?} vs {:?}",
        frame.z
    );
}

/// A helix whose Darboux axis and rate match the LOCAL generator at `s`.
///
/// Frenet with constant `omega` is a screw motion, and this is its exact
/// closed form. Used to build a reference for a varying law by stepping in
/// many small pieces, each locally a screw -- an integrator of a different
/// kind from the one under test, so agreement is evidence rather than
/// tautology.
fn screw_step(
    omega: Vec3,
    h: Scalar,
    tangent: Vec3,
    normal: Vec3,
    binormal: Vec3,
) -> (Vec3, [Vec3; 3]) {
    let rate = omega.length();
    if rate < 1e-14 {
        return (tangent * h, [tangent, normal, binormal]);
    }
    let axis = omega / rate;
    let turn = rate * h;
    // Body-frame displacement of a screw about `axis` over `h`.
    let along = axis.dot(Vec3::X);
    let local = axis * (along * h)
        + (Vec3::X - axis * along) * (turn.sin() / rate)
        + axis.cross(Vec3::X) * ((1.0 - turn.cos()) / rate);
    let world = tangent * local.x + normal * local.y + binormal * local.z;
    // Rodrigues on each body axis.
    let rotate = |v: Vec3| {
        v * turn.cos() + axis.cross(v) * turn.sin() + axis * (axis.dot(v)) * (1.0 - turn.cos())
    };
    let (rx, ry, rz) = (rotate(Vec3::X), rotate(Vec3::Y), rotate(Vec3::Z));
    let to_world = |v: Vec3| tangent * v.x + normal * v.y + binormal * v.z;
    (world, [to_world(rx), to_world(ry), to_world(rz)])
}

#[test]
fn the_commutator_term_is_load_bearing_on_a_non_commuting_law() {
    // tau/k varies, so generators at different arc lengths do not commute and
    // the Magnus commutator term is what makes the scheme fourth order.
    //
    // Reference: 20,000 locally-exact screw steps. Each step uses the closed
    // form for CONSTANT omega, so the reference never uses the Magnus
    // expansion at all -- dropping the commutator cannot move it.
    let length = 40.0;
    let k = CurvatureLaw::clothoid(0.01, 0.09, length);
    let tau = CurvatureLaw::circular(0.05);
    let curve = Intrinsic3::new(start_frame(), k, tau, length);

    let steps = 20_000;
    let h = length / steps as Scalar;
    let mut position = Vec3::ZERO;
    let mut axes = [Vec3::X, Vec3::Y, Vec3::Z];
    for step in 0..steps {
        let s = (step as Scalar + 0.5) * h; // midpoint of this step
        let curvature = 0.01 + (0.09 - 0.01) * s / length;
        let (delta, next) = screw_step(
            Vec3::new(0.05, 0.0, curvature),
            h,
            axes[0],
            axes[1],
            axes[2],
        );
        position += delta;
        axes = next;
    }

    let got = frenet_point(&curve, length).expect("evaluates");
    let error = (Vec3::new(got.x, got.y, got.z) - position).length();
    // With the commutator the error is ~5e-6; dropping it gives ~2.7e-2, a
    // factor of 5000. This bound sits between the two, so the test states a
    // real accuracy claim rather than merely 'it runs'.
    assert!(
        error < 1e-4,
        "Magnus must track the locally-exact screw reference, error {error:e}"
    );
}

#[test]
fn a_non_commuting_law_differs_from_the_frozen_generator_answer() {
    // Guards the commutator term. If someone replaces the Magnus step with
    // exp(int Omega), this curve is where it shows: a helix would not.
    let length = 40.0;
    let curve = Intrinsic3::new(
        start_frame(),
        CurvatureLaw::clothoid(0.01, 0.09, length),
        CurvatureLaw::circular(0.05),
        length,
    );
    // exp(int Omega) over the whole span in one go: total turning about each
    // axis, exponentiated once. That is the wrong answer here.
    let bend = Intrinsic2::new(
        Frame2 {
            origin: Point2::new(0.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        CurvatureLaw::clothoid(0.01, 0.09, length),
        length,
    )
    .total_turning()
    .expect("turning integrates");
    let twist = 0.05 * length;
    let naive_angle = (bend * bend + twist * twist).sqrt();

    let frame = frenet_frame(&curve, length).expect("evaluates");
    // The true tangent, versus the one a single exponential would give.
    let axis = Vec3::new(twist, 0.0, bend) / naive_angle;
    let cos = naive_angle.cos();
    let naive_tangent = Vec3::X * cos
        + axis.cross(Vec3::X) * naive_angle.sin()
        + axis * (axis.dot(Vec3::X)) * (1.0 - cos);
    let gap = (frame.x - naive_tangent).length();
    assert!(
        gap > 1e-3,
        "a non-commuting law must differ from the frozen-generator answer, gap={gap:e}"
    );
}

#[test]
fn a_helix_is_recognised_and_a_general_curve_is_not() {
    let helix = Intrinsic3::new(
        start_frame(),
        CurvatureLaw::circular(0.1),
        CurvatureLaw::circular(0.02),
        10.0,
    );
    assert!(helix.is_helical(), "constant k and tau is a helix");
    assert!(!helix.is_planar(), "non-zero torsion is not planar");

    let flat = Intrinsic3::new(
        start_frame(),
        CurvatureLaw::clothoid(0.0, 0.1, 10.0),
        CurvatureLaw::straight(),
        10.0,
    );
    assert!(flat.is_planar(), "zero torsion is planar");
    assert!(!flat.is_helical(), "a varying curvature is not a helix");
}

#[test]
fn a_non_orthonormal_start_frame_is_refused() {
    let skewed = Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: Vec3::X,
        y: Vec3::new(0.5, 0.5, 0.0), // neither unit nor perpendicular
        z: Vec3::Z,
    };
    let curve = Intrinsic3::new(
        skewed,
        CurvatureLaw::circular(0.1),
        CurvatureLaw::circular(0.01),
        10.0,
    );
    assert!(
        frenet_point(&curve, 5.0).is_err(),
        "a start frame that is not a rotation must be refused, not fixed up"
    );

    // Orthonormal but LEFT-handed: every length and dot product is right, only
    // the handedness is wrong, so this is the case a norms-only check misses.
    let mirrored = Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: -Vec3::Z,
    };
    let flipped = Intrinsic3::new(
        mirrored,
        CurvatureLaw::circular(0.1),
        CurvatureLaw::circular(0.01),
        10.0,
    );
    assert!(
        frenet_point(&flipped, 5.0).is_err(),
        "a left-handed start frame must be refused"
    );
}

#[test]
fn arc_length_outside_the_curve_is_refused() {
    let curve = Intrinsic3::new(
        start_frame(),
        CurvatureLaw::circular(0.1),
        CurvatureLaw::circular(0.01),
        10.0,
    );
    assert!(frenet_point(&curve, -1.0).is_err(), "negative arc length");
    assert!(frenet_point(&curve, 10.5).is_err(), "past the end");
    assert!(
        frenet_point(&curve, 10.0).is_ok(),
        "the end itself is inside"
    );
}
