//! Certified positions of elevated and banked curves over an intrinsic
//! plan (#393).
//!
//! The points of a gradient curve and of a banked one over a clothoid are
//! held to a reference written here and independent of the rule under test
//! (ADR 0060): the Fresnel integrals' power series and the closed-form
//! profile and pivot. The chord bound of a curve over a fast harmonic plan,
//! whose evaluated points are centimetres off, is held to a compensated
//! Simpson reference, and the families whose plan is a line or a circle
//! keep the exact-chord bound bit for bit.

use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, Circle2, CurvatureLaw,
    Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2, Line2,
};
use axiolid_evaluate::bound::{certifies_flattening3, chord_bound3, ElevatedPositionTable};
use axiolid_evaluate::curve::evaluate3;
use axiolid_evaluate::elevated::{banked_chord_bound, elevated_chord_bound};

fn world() -> Frame2 {
    Frame2 {
        origin: Point2::new(20.0, 10.0),
        x: Vec2::new(0.8, 0.6),
        y: Vec2::new(-0.6, 0.8),
    }
}

/// `int_0^s (cos(a u^2), sin(a u^2)) du` by its power series.
fn fresnel(a: Scalar, s: Scalar) -> Vec2 {
    let (mut x, mut y) = (0.0, 0.0);
    let mut factorial = 1.0;
    for n in 0..40 {
        if n > 0 {
            factorial *= n as Scalar;
        }
        let value = a.powi(n) * s.powi(2 * n + 1) / (factorial * (2 * n + 1) as Scalar);
        match n % 4 {
            0 => x += value,
            1 => y += value,
            2 => x -= value,
            _ => y -= value,
        }
    }
    Vec2::new(x, y)
}

fn place(local: Vec2, z: Scalar) -> Point3 {
    let w = world();
    Point3::new(
        w.origin.x + w.x.x * local.x + w.y.x * local.y,
        w.origin.y + w.x.y * local.x + w.y.y * local.y,
        z,
    )
}

fn clothoid_plan() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        world(),
        CurvatureLaw::clothoid(0.0, 1.0 / 120.0, 60.0),
        60.0,
    ))
}

/// The clothoid's heading is `d^2 / 14400`.
fn clothoid_point(d: Scalar, z: Scalar) -> Point3 {
    place(fresnel(1.0 / 14400.0, d), z)
}

fn gradient_curve() -> Curve3 {
    Curve3::Elevated(Elevated3::new(
        clothoid_plan(),
        ElevationLaw::parabolic(50.0, 0.05, -0.03, 60.0),
    ))
}

fn canted_curve() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(
            clothoid_plan(),
            ElevationLaw::parabolic(50.0, 0.05, -0.03, 60.0),
        ),
        CantLaw::new(vec![
            CantPiece::sine(30.0, 0.0, 0.12),
            CantPiece::constant(30.0, 0.12),
        ]),
        CantLaw::new(vec![
            CantPiece::sine(30.0, 0.0, 0.06),
            CantPiece::constant(30.0, 0.06),
        ]),
        1.5,
        BankConvention::TangentRotation,
    ))
}

/// `z = 50 + 0.05 d - (0.08 / 120) d^2`, the parabolic profile.
fn parabola(d: Scalar) -> Scalar {
    50.0 + 0.05 * d - 0.08 / 120.0 * d * d
}

/// The pivot of [`canted_curve`]: `0.06 (xi - sin(2 pi xi) / 2 pi)`, then
/// `0.06`.
fn pivot(d: Scalar) -> Scalar {
    if d >= 30.0 {
        return 0.06;
    }
    let xi = d / 30.0;
    let tau = core::f64::consts::TAU;
    0.06 * (xi - (tau * xi).sin() / tau)
}

fn assert_certified(name: &str, curve: &Curve3, height: impl Fn(Scalar) -> Scalar, target: Scalar) {
    let table = ElevatedPositionTable::new(curve, target, 60.0).expect("a table");
    let max = table.max_error();
    assert!(max <= target, "{name}: max_error {max:e} above {target:e}");
    let mut worst: Scalar = 0.0;
    for k in 0..=120 {
        let d = 0.5 * k as Scalar;
        let (point, error) = table.point3(d).expect("a point");
        let exact = clothoid_point(d, height(d));
        let gap = point.distance(exact);
        worst = worst.max(gap);
        assert!(
            gap <= error + 1e-13,
            "{name}: d = {d}: {gap:e} from the exact point, above its error {error:e}"
        );
        assert!(
            error <= max,
            "{name}: d = {d}: error {error:e} above {max:e}"
        );
    }
    eprintln!("{name}: worst {worst:e}, max_error {max:e}");
}

#[test]
fn points_of_a_gradient_curve_over_a_clothoid_are_certified() {
    assert_certified("gradient", &gradient_curve(), parabola, 1e-9);
}

#[test]
fn points_of_a_canted_curve_over_a_clothoid_are_certified() {
    assert_certified("canted", &canted_curve(), |d| parabola(d) + pivot(d), 1e-9);
}

#[test]
fn the_table_refuses_what_it_cannot_certify() {
    let line = Curve2::Line(Line2 {
        origin: Point2::new(0.0, 0.0),
        direction: Vec2::new(1.0, 0.0),
    });
    let on_line = Curve3::Elevated(Elevated3::new(
        line,
        ElevationLaw::constant_grade(0.0, 0.01),
    ));
    assert!(ElevatedPositionTable::new(&on_line, 1e-6, 10.0).is_none());
    let curve = gradient_curve();
    assert!(ElevatedPositionTable::new(&curve, 1e-6, 61.0).is_none());
    assert!(ElevatedPositionTable::new(&curve, 0.0, 60.0).is_none());
    assert!(ElevatedPositionTable::new(&curve, 1e-6, 0.0).is_none());
    let table = ElevatedPositionTable::new(&curve, 1e-6, 30.0).expect("a table");
    assert!(table.point3(30.5).is_none());
    assert!(table.point3(-0.5).is_none());
    // A plain space line is not an elevated curve.
    let space = Curve3::Line(axiolid_curve::Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    });
    assert!(ElevatedPositionTable::new(&space, 1e-6, 1.0).is_none());
}

// --- the fast harmonic ---------------------------------------------------------

const AMPLITUDE: Scalar = 0.01;
const FREQUENCY: Scalar = 50.0;

fn harmonic_plan(amplitude: Scalar) -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        world(),
        CurvatureLaw::Sinusoid {
            mean: 0.0,
            amplitude,
            angular_frequency: FREQUENCY,
            phase: 0.0,
        },
        100.0,
    ))
}

fn harmonic_elevated(amplitude: Scalar) -> Curve3 {
    Curve3::Elevated(Elevated3::new(
        harmonic_plan(amplitude),
        ElevationLaw::constant_grade(50.0, 0.04),
    ))
}

fn harmonic_banked(amplitude: Scalar) -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(
            harmonic_plan(amplitude),
            ElevationLaw::constant_grade(50.0, 0.04),
        ),
        CantLaw::new(vec![CantPiece::constant(100.0, 0.05)]),
        CantLaw::zero(100.0),
        1.5,
        BankConvention::TangentRotation,
    ))
}

/// A Neumaier-compensated running sum.
#[derive(Clone, Copy, Default)]
struct Sum {
    total: Scalar,
    carry: Scalar,
}

impl Sum {
    fn add(&mut self, x: Scalar) {
        let t = self.total + x;
        self.carry += if self.total.abs() >= x.abs() {
            (self.total - t) + x
        } else {
            (x - t) + self.total
        };
        self.total = t;
    }

    fn value(self) -> Scalar {
        self.total + self.carry
    }
}

/// The exact points of the harmonic plan (`theta = a (1 - cos(w s))`,
/// `a = amplitude / w`) at the multiples of `1 / 4096` in `[from, to]`, by
/// one pass of compensated Simpson with step `1 / 8192`, placed by the
/// world frame at `z = 50 + 0.04 d`.
fn harmonic_points(amplitude: Scalar, from: Scalar, to: Scalar) -> Vec<(Scalar, Point3)> {
    let a = amplitude / FREQUENCY;
    let f = |s: Scalar| {
        let theta = a * (1.0 - (FREQUENCY * s).cos());
        (theta.cos(), theta.sin())
    };
    let step = 1.0 / 8192.0;
    let (mut x, mut y) = (Sum::default(), Sum::default());
    let mut out = Vec::new();
    let panels = (to / (2.0 * step)).round() as usize;
    for m in 0..=panels {
        let s = 2.0 * step * m as Scalar;
        if s >= from {
            out.push((s, place(Vec2::new(x.value(), y.value()), 50.0 + 0.04 * s)));
        }
        let (f0, f1, f2) = (f(s), f(s + step), f(s + 2.0 * step));
        x.add(step / 3.0 * (f0.0 + 4.0 * f1.0 + f2.0));
        y.add(step / 3.0 * (f0.1 + 4.0 * f1.1 + f2.1));
    }
    out
}

/// The distance of `p` from the segment `a b`.
fn segment_distance(p: Point3, a: Point3, b: Point3) -> Scalar {
    let ab: Vec3 = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
    p.distance(a + ab * t)
}

/// The worst distance of the exact curve over `[lo, hi]` from the chord of
/// the evaluated ends.
fn worst_chord_gap(curve: &Curve3, amplitude: Scalar, lo: Scalar, hi: Scalar) -> Scalar {
    let a = evaluate3(curve, lo).unwrap();
    let b = evaluate3(curve, hi).unwrap();
    harmonic_points(amplitude, lo, hi)
        .into_iter()
        .filter(|&(s, _)| s <= hi)
        .map(|(_, p)| segment_distance(p, a, b))
        .fold(0.0, Scalar::max)
}

#[test]
fn a_chord_bound_covers_the_evaluated_ends_of_a_fast_harmonic_plan() {
    let (lo, hi) = (96.4375, 96.5625);
    for (name, curve) in [
        ("elevated", harmonic_elevated(AMPLITUDE)),
        ("banked", harmonic_banked(AMPLITUDE)),
    ] {
        let exact = worst_chord_gap(&curve, AMPLITUDE, lo, hi);
        let bound = chord_bound3(&curve, lo, hi).expect("a bound");
        eprintln!("{name}: exact {exact:e}, bound {bound:e}");
        // The exact curve lies 2.8348e-3 from the chord of the evaluated
        // ends, 145 times the exact-chord bound `h^2 / 8 k`, 1.95e-5.
        assert!(exact >= 2.8348e-3, "{name}: the reference moved: {exact:e}");
        assert!(
            bound >= exact,
            "{name}: {bound:e} below the exact {exact:e}"
        );
        assert!(bound <= 2.9e-3, "{name}: {bound:e} is not tight");
        // At a point the bound is the evaluated point's distance from the
        // exact one, 1.4e-2 at 96.5.
        let at = chord_bound3(&curve, 96.5, 96.5).expect("a bound at a point");
        assert!(at >= 1.39e-2, "{name}: {at:e} at a point");
    }
}

#[test]
fn a_chord_bound_over_a_smaller_amplitude_still_covers_the_exact_chord() {
    let (lo, hi) = (96.4375, 96.5625);
    let curve = harmonic_elevated(0.002);
    let exact = worst_chord_gap(&curve, 0.002, lo, hi);
    let bound = chord_bound3(&curve, lo, hi).expect("a bound");
    eprintln!("exact {exact:e}, bound {bound:e}");
    // On main: 2.5e-6 against 5.68e-4.
    assert!(exact >= 5.6e-4, "the reference moved: {exact:e}");
    assert!(bound >= exact, "{bound:e} below the exact {exact:e}");
}

#[test]
fn a_chord_bound_over_a_line_or_a_circle_plan_is_the_exact_chord_bound() {
    let line = Curve2::Line(Line2 {
        origin: Point2::new(20.0, 10.0),
        direction: Vec2::new(0.8, 0.6),
    });
    let arc = Curve2::Circle(Circle2 {
        frame: world(),
        radius: 150.0,
    });
    for plan in [line, arc] {
        let elevated = Elevated3::new(plan, ElevationLaw::parabolic(50.0, 0.05, -0.03, 60.0));
        let banked = Banked3::new(
            elevated.clone(),
            CantLaw::new(vec![CantPiece::constant(60.0, 0.1)]),
            CantLaw::zero(60.0),
            1.5,
            BankConvention::TangentRotation,
        );
        for (a, b) in [(1.0, 1.5), (10.0, 12.0), (20.0, 20.25)] {
            let got = chord_bound3(&Curve3::Elevated(elevated.clone()), a, b);
            assert_eq!(
                got.map(Scalar::to_bits),
                elevated_chord_bound(&elevated, a, b).map(Scalar::to_bits)
            );
            let got = chord_bound3(&Curve3::Banked(banked.clone()), a, b);
            assert_eq!(
                got.map(Scalar::to_bits),
                banked_chord_bound(&banked, a, b).map(Scalar::to_bits)
            );
        }
    }
}

/// A chain plan's evaluated ends carry no certified gap (#438): the exact
/// chord bound under-reports (1.95e-5 against 2.8e-3 over this span, and
/// zero at a point against 1.4e-2), so the curve is not claimed certified.
#[test]
fn a_chain_plan_is_not_certified_for_flattening() {
    let chain = Curve2::Chain(Chain2::new(
        world(),
        vec![
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::Constant { curvature: 0.0 },
                length: 20.0,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::Sinusoid {
                    mean: 0.0,
                    amplitude: AMPLITUDE,
                    angular_frequency: FREQUENCY,
                    phase: 0.0,
                },
                length: 100.0,
            },
        ],
    ));
    let curve = Curve3::Elevated(Elevated3::new(
        chain,
        ElevationLaw::constant_grade(50.0, 0.04),
    ));
    assert!(!certifies_flattening3(&curve));
    let bound = chord_bound3(&curve, 116.4375, 116.5625).expect("the exact chord");
    assert!(bound < 1e-4, "{bound:e}");
}
