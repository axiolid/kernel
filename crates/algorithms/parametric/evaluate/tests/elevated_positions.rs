//! Certified positions of elevated and banked curves over an intrinsic
//! plan (#393) or a chain plan (#438).
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
use axiolid_evaluate::bound::{
    certifies_flattening3, chord_bound3, ElevatedPositionTable, SeamSide,
};
use axiolid_evaluate::curve::{derivative3, evaluate3, second_derivative3};
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

// --- chain plans (#438) ----------------------------------------------------------

/// A straight 20 m, then the fast harmonic over 100 m: the evaluator's point
/// is 1.4e-2 off at 116.5.
fn harmonic_chain() -> Curve2 {
    Curve2::Chain(Chain2::new(
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
    ))
}

fn chain_elevated() -> Curve3 {
    Curve3::Elevated(Elevated3::new(
        harmonic_chain(),
        ElevationLaw::constant_grade(50.0, 0.04),
    ))
}

fn chain_banked() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(harmonic_chain(), ElevationLaw::constant_grade(50.0, 0.04)),
        CantLaw::new(vec![CantPiece::constant(120.0, 0.05)]),
        CantLaw::zero(120.0),
        1.5,
        BankConvention::TangentRotation,
    ))
}

/// The exact points of the chain at the multiples of `1 / 4096` in
/// `[from, to]`: the straight's 20 m along the start tangent, then the
/// harmonic's.
fn chain_points(from: Scalar, to: Scalar) -> Vec<(Scalar, Point3)> {
    let shift = Vec3::new(20.0 * world().x.x, 20.0 * world().x.y, 0.04 * 20.0);
    harmonic_points(AMPLITUDE, from - 20.0, to - 20.0)
        .into_iter()
        .map(|(s, p)| (s + 20.0, p + shift))
        .collect()
}

#[test]
fn a_chord_bound_over_a_chain_plan_covers_the_evaluated_ends() {
    // The exact-chord bound is 1.95e-5 here, the exact arc 2.8348e-3 from
    // the chord of the evaluated ends: the chain's evaluated ends are as
    // far off as the harmonic's alone.
    let (lo, hi) = (116.4375, 116.5625);
    for (name, curve) in [("elevated", chain_elevated()), ("banked", chain_banked())] {
        // A banked curve is never claimed whole; the elevated one is.
        assert_eq!(certifies_flattening3(&curve), name == "elevated", "{name}");
        let a = evaluate3(&curve, lo).unwrap();
        let b = evaluate3(&curve, hi).unwrap();
        let exact = chain_points(lo, hi)
            .into_iter()
            .filter(|&(s, _)| s <= hi)
            .map(|(_, p)| {
                // The banked curve rides the pivot, zero here, and the
                // elevated height is the same closed form.
                segment_distance(p, a, b)
            })
            .fold(0.0, Scalar::max);
        let bound = chord_bound3(&curve, lo, hi).expect("a bound");
        assert!(exact >= 2.8348e-3, "{name}: the reference moved: {exact:e}");
        assert!(
            bound >= exact,
            "{name}: {bound:e} below the exact {exact:e}"
        );
        assert!(bound <= 2.9e-3, "{name}: {bound:e} is not tight");
        let at = chord_bound3(&curve, 116.5, 116.5).expect("a bound at a point");
        assert!(
            (1.39e-2..=1.5e-2).contains(&at),
            "{name}: {at:e} at a point"
        );
    }
}

#[test]
fn the_chain_plan_points_of_a_gradient_curve_are_certified() {
    // As the intrinsic plan's: every point of the table within its own
    // error of the exact one, and the error under the table's maximum.
    let curve = chain_elevated();
    let table = ElevatedPositionTable::new(&curve, 1e-6, 119.0).expect("a table");
    assert!(table.max_error() > 0.0);
    for (s, exact) in chain_points(0.0, 119.0).into_iter().step_by(509) {
        let (point, error) = table.point3(s).expect("a point");
        assert!(error <= table.max_error());
        assert!(
            point.distance(exact) <= error + 1e-12,
            "s = {s}: {:e} from the exact point, above its error {error:e}",
            point.distance(exact)
        );
    }
}

/// The cubic `y = x^3 / (6 R X)` over `X` as a Bezier, with its arc length
/// by the binomial series (as `tests/arc_length_chain.rs`).
fn cubic(a: Scalar, reach: Scalar) -> Curve2 {
    Curve2::BSpline(axiolid_curve::BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(reach / 3.0, 0.0),
            Point2::new(2.0 * reach / 3.0, 0.0),
            Point2::new(reach, a * reach.powi(3)),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: axiolid_curve::KnotSpec::PiecewiseBezier,
    })
}

fn cubic_length(a: Scalar, x: Scalar) -> Scalar {
    let z = 9.0 * a * a * x.powi(4);
    let (mut binomial, mut power, mut sum) = (1.0, 1.0, 0.0);
    for n in 0..200 {
        let term = binomial * power / (4.0 * n as Scalar + 1.0);
        sum += term;
        if term.abs() < 1e-30 {
            break;
        }
        binomial *= (0.5 - n as Scalar) / (n as Scalar + 1.0);
        power *= z;
    }
    x * sum
}

fn cubic_abscissa(a: Scalar, length: Scalar) -> Scalar {
    let mut x = length;
    for _ in 0..50 {
        x -= (cubic_length(a, x) - length) / (1.0 + 9.0 * a * a * x.powi(4)).sqrt();
    }
    x
}

const STRAIGHT: Scalar = 50.0;
const TRANSITION: Scalar = 60.0;
const RADIUS: Scalar = 300.0;

fn layout() -> Curve2 {
    Curve2::Chain(Chain2::new(
        world(),
        vec![
            ChainPiece2::Parametric {
                curve: Curve2::Line(Line2 {
                    origin: Vec2::ZERO,
                    direction: Vec2::X,
                }),
                start: 0.0,
                length: STRAIGHT,
            },
            ChainPiece2::Parametric {
                curve: cubic(1.0 / (6.0 * RADIUS * TRANSITION), TRANSITION),
                start: 0.0,
                length: TRANSITION,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::circular(1.0 / RADIUS),
                length: 40.0,
            },
        ],
    ))
}

/// The layout's exact local point, heading and curvature at `d`.
fn layout_exact(d: Scalar) -> (Vec2, Scalar, Scalar) {
    let a = 1.0 / (6.0 * RADIUS * TRANSITION);
    if d < STRAIGHT {
        return (Vec2::new(d, 0.0), 0.0, 0.0);
    }
    if d < STRAIGHT + TRANSITION {
        let x = cubic_abscissa(a, d - STRAIGHT);
        let slope = 3.0 * a * x * x;
        return (
            Vec2::new(STRAIGHT + x, a * x.powi(3)),
            slope.atan(),
            6.0 * a * x / (1.0 + slope * slope).powf(1.5),
        );
    }
    let x = cubic_abscissa(a, TRANSITION);
    let heading = (3.0 * a * x * x).atan();
    let join = Vec2::new(STRAIGHT + x, a * x.powi(3));
    let phi = (d - STRAIGHT - TRANSITION) / RADIUS;
    let arc = Vec2::new(RADIUS * phi.sin(), RADIUS * (1.0 - phi.cos()));
    let (s, c) = heading.sin_cos();
    (
        join + Vec2::new(c * arc.x - s * arc.y, s * arc.x + c * arc.y),
        heading + phi,
        1.0 / RADIUS,
    )
}

#[test]
fn an_elevated_curve_over_a_cubic_chain_reads_its_series_plan() {
    // Gradient 0.01 into 0.0148 over the chain's 150 m: z = 10 + 0.01 d +
    // 0.0048 d^2 / 300.
    let curve = Curve3::Elevated(Elevated3::new(
        layout(),
        ElevationLaw::parabolic(10.0, 0.01, 0.0148, 150.0),
    ));
    let table = ElevatedPositionTable::new(&curve, 1e-9, 150.0).expect("a table");
    let [e0, e1, e2] = table.jet_errors();
    assert_eq!(e0, table.max_error());
    assert!(e0 > 0.0 && e1 > 0.0 && e2 > 0.0);
    let w = world();
    for k in 0..=150 {
        let d = k as Scalar;
        let (local, heading, kappa) = layout_exact(d);
        let z = 10.0 + 0.01 * d + 0.0048 * d * d / 300.0;
        let exact = place(local, z);
        let (point, derivative, second) = {
            let (p, c1, c2) = table.jet3_on(d, SeamSide::Outgoing).expect("a jet");
            (p, c1, c2)
        };
        let (read, error) = table.point3(d).expect("a point");
        assert_eq!(read, point);
        assert!(error <= e0);
        assert!(point.distance(exact) <= error, "d = {d}");
        // c' = (tangent, z'), c'' = (k n, z'').
        let tangent = w.x * heading.cos() + w.y * heading.sin();
        let normal = w.x * -heading.sin() + w.y * heading.cos();
        let grade = 0.01 + 0.0048 * d / 150.0;
        let exact1 = Vec3::new(tangent.x, tangent.y, grade);
        let exact2 = Vec3::new(normal.x * kappa, normal.y * kappa, 0.0048 / 150.0);
        assert!((derivative - exact1).length() <= e1, "c' at {d}");
        assert!((second - exact2).length() <= e2, "c'' at {d}");
    }
}

#[test]
fn the_jet_over_an_intrinsic_plan_is_the_curves_own() {
    let curve = gradient_curve();
    let table = ElevatedPositionTable::new(&curve, 1e-6, 60.0).expect("a table");
    let [_, e1, e2] = table.jet_errors();
    assert_eq!((e1, e2), (0.0, 0.0));
    for d in [0.0, 12.5, 30.0, 59.0] {
        let (point, c1, c2) = table.jet3_on(d, SeamSide::Outgoing).expect("a jet");
        assert_eq!(point, table.point3(d).unwrap().0);
        assert_eq!(c1, derivative3(&curve, d).unwrap());
        assert_eq!(c2, second_derivative3(&curve, d).unwrap());
    }
}
