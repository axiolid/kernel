//! Certified frame-law bounds along intrinsic curves (#361): a 2D
//! clothoid, a level planar space curve, and an elevated curve over a
//! piecewise plan of a line, a clothoid, an arc, a clothoid and a line,
//! with a parabolic sag and a circular-arc crest.
//!
//! Every reference is independent of the kernel's quadratures: the plan's
//! points from the Fresnel integrals' power series and closed forms, its
//! headings and the profile's heights and grades in closed form, and the
//! frames built here as the section frames are defined. The derivative
//! bounds are checked by finite differences of that frame; each frame the
//! kernel reads must lie within the reported errors of it.

use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, CurvatureLaw, Curve2, Curve3,
    Elevated3, ElevationLaw, Intrinsic2, Intrinsic3, RailSide,
};
use axiolid_evaluate::station::{
    piece_frame_bounds, station_section2, station_section3, FrameLaw, FrameLawBounds, FrameReading,
    SectionFrame, StationCurve, StationPiece,
};

const TOLERANCE: Scalar = 1e-6;

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
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

fn rotate(angle: Scalar, v: Vec2) -> Vec2 {
    let (s, c) = angle.sin_cos();
    Vec2::new(c * v.x - s * v.y, s * v.x + c * v.y)
}

#[derive(Clone, Copy)]
enum Piece {
    Line(Scalar),
    Rising(Scalar, Scalar),
    Arc(Scalar, Scalar),
    Falling(Scalar, Scalar),
}

impl Piece {
    fn length(self) -> Scalar {
        match self {
            Piece::Line(l) | Piece::Rising(l, _) | Piece::Arc(l, _) | Piece::Falling(l, _) => l,
        }
    }

    fn law(self) -> CurvatureLaw {
        match self {
            Piece::Line(_) => CurvatureLaw::Constant { curvature: 0.0 },
            Piece::Rising(l, r) => CurvatureLaw::clothoid(0.0, 1.0 / r, l),
            Piece::Arc(_, r) => CurvatureLaw::Constant { curvature: 1.0 / r },
            Piece::Falling(l, r) => CurvatureLaw::clothoid(1.0 / r, 0.0, l),
        }
    }

    fn turn(self) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Rising(l, r) | Piece::Falling(l, r) => 0.5 * l / r,
            Piece::Arc(l, r) => l / r,
        }
    }

    /// The heading turned over the first `s`.
    fn heading(self, s: Scalar) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Rising(l, r) => 0.5 * s * s / (r * l),
            Piece::Arc(_, r) => s / r,
            Piece::Falling(l, r) => self.turn() - 0.5 * (l - s) * (l - s) / (r * l),
        }
    }

    fn displacement(self, s: Scalar) -> Vec2 {
        match self {
            Piece::Line(_) => Vec2::new(s, 0.0),
            Piece::Rising(l, r) => fresnel(0.5 / (r * l), s),
            Piece::Arc(_, r) => Vec2::new(r * (s / r).sin(), r * (1.0 - (s / r).cos())),
            Piece::Falling(l, r) => {
                let a = 0.5 / (r * l);
                let tail = fresnel(a, l) - fresnel(a, l - s);
                rotate(self.turn(), Vec2::new(tail.x, -tail.y))
            }
        }
    }
}

struct Alignment {
    curve: Intrinsic2,
    pieces: Vec<Piece>,
}

impl Alignment {
    fn new(pieces: Vec<Piece>) -> Self {
        let mut breaks = Vec::new();
        let mut length = 0.0;
        for piece in &pieces[..pieces.len() - 1] {
            length += piece.length();
            breaks.push(length);
        }
        length += pieces.last().unwrap().length();
        let law = CurvatureLaw::piecewise(breaks, pieces.iter().map(|p| p.law()).collect());
        Self {
            curve: Intrinsic2::new(frame(), law, length),
            pieces,
        }
    }

    /// The exact point and heading at `s`, in the start frame's plane.
    fn exact(&self, s: Scalar) -> (Vec2, Scalar) {
        let (mut at, mut heading, mut begin) = (Vec2::ZERO, 0.0, 0.0);
        for piece in &self.pieces {
            let inside = (s - begin).clamp(0.0, piece.length());
            at += rotate(heading, piece.displacement(inside));
            if s <= begin + piece.length() {
                return (at, heading + piece.heading(inside));
            }
            heading += piece.turn();
            begin += piece.length();
        }
        unreachable!("past the end")
    }

    /// The point and unit tangent in the world.
    fn world(&self, s: Scalar) -> (Vec2, Vec2) {
        let (local, heading) = self.exact(s);
        let f = frame();
        (
            f.origin + f.x * local.x + f.y * local.y,
            f.x * heading.cos() + f.y * heading.sin(),
        )
    }
}

fn alignment() -> Alignment {
    Alignment::new(vec![
        Piece::Line(40.0),
        Piece::Rising(60.0, 300.0),
        Piece::Arc(80.0, 300.0),
        Piece::Falling(60.0, 300.0),
        Piece::Line(40.0),
    ])
}

/// A parabolic sag over `[0, 120]` into a circular crest of radius 4000.
fn profile() -> ElevationLaw {
    ElevationLaw::Piecewise {
        breaks: vec![120.0],
        laws: vec![
            ElevationLaw::parabolic(10.0, 0.01, 0.0148, 120.0),
            ElevationLaw::circular_arc(11.488, 0.0148, -4000.0),
        ],
    }
}

/// The profile's height and grade at `d`, in closed form.
fn profile_at(d: Scalar) -> (Scalar, Scalar) {
    if d <= 120.0 {
        return (10.0 + 0.01 * d + 2e-5 * d * d, 0.01 + 4e-5 * d);
    }
    let (r, t0) = (-4000.0_f64, 0.0148_f64.atan());
    let t = (t0.sin() + (d - 120.0) / r).asin();
    (11.488 + r * (t0.cos() - t.cos()), t.tan())
}

/// A reference frame, built as the section frames are defined.
#[derive(Clone, Copy)]
struct Reference {
    point: Vec3,
    tangent: Vec3,
    lateral: Vec3,
    up: Vec3,
}

fn reference_up(point: Vec3, tangent: Vec3) -> Reference {
    let right = tangent.cross(Vec3::Z);
    let right = right / right.length();
    Reference {
        point,
        tangent,
        lateral: -right,
        up: right.cross(tangent),
    }
}

fn planar(point: Vec3, tangent: Vec3) -> Reference {
    Reference {
        point,
        tangent,
        lateral: Vec3::new(-tangent.y, tangent.x, 0.0),
        up: Vec3::Z,
    }
}

fn plan_reading(r: Reference) -> Reference {
    let t = Vec3::new(r.tangent.x, r.tangent.y, 0.0);
    planar(r.point, t / t.length())
}

fn elevated_at(a: &Alignment, d: Scalar) -> Reference {
    let (p, t) = a.world(d);
    let (z, g) = profile_at(d);
    let size = (1.0 + g * g).sqrt();
    reference_up(
        Vec3::new(p.x, p.y, z),
        Vec3::new(t.x / size, t.y / size, g / size),
    )
}

fn plan_at(a: &Alignment, d: Scalar) -> Reference {
    let (p, t) = a.world(d);
    planar(Vec3::new(p.x, p.y, 0.0), Vec3::new(t.x, t.y, 0.0))
}

type Axis = fn(&Reference) -> Vec3;

/// The bounds against first and second differences of the reference frame
/// over a lattice of `[lo, hi]`.
fn assert_bounded(
    reference: impl Fn(Scalar) -> Reference,
    b: &FrameLawBounds,
    lo: Scalar,
    hi: Scalar,
) {
    let h = 1e-3 * (hi - lo);
    for i in 1..20 {
        let s = lo + (hi - lo) * Scalar::from(i) / 20.0;
        let (a, m, c) = (reference(s - h), reference(s), reference(s + h));
        let first = |f: Axis| ((f(&c) - f(&a)) / (2.0 * h)).length();
        let second = |f: Axis| ((f(&c) - 2.0 * f(&m) + f(&a)) / (h * h)).length();
        let slack = |bound: Scalar| bound * (1.0 + 1e-3) + 1e-4;
        let axes: [(Axis, [Scalar; 2], &str); 4] = [
            (|f| f.point, b.point, "point"),
            (|f| f.tangent, b.tangent, "tangent"),
            (|f| f.lateral, b.lateral, "lateral"),
            (|f| f.up, b.up, "up"),
        ];
        for (f, [b1, b2], name) in axes {
            assert!(first(f) <= slack(b1), "{name}' {} > {b1} at {s}", first(f));
            assert!(
                second(f) <= slack(b2),
                "{name}'' {} > {b2} at {s}",
                second(f)
            );
        }
    }
}

/// Everything the module promises of the bounds over `[lo, hi]`: the
/// derivatives, an exact measure, and the kernel's reads within the errors
/// of the reference.
fn check(
    law: &FrameLaw<'_, '_>,
    reference: impl Fn(Scalar) -> Reference + Copy,
    read: impl Fn(Scalar) -> SectionFrame,
    lo: Scalar,
    hi: Scalar,
) {
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let bounds = law.bounds(lo, hi, reading).expect("bounded");
        let wanted = move |s: Scalar| match reading {
            FrameReading::Section => reference(s),
            _ => plan_reading(reference(s)),
        };
        assert_bounded(wanted, &bounds, lo, hi);
        assert_eq!(bounds.measure, 0.0);
        assert!(
            bounds.point_error <= 2.0 * TOLERANCE,
            "{}",
            bounds.point_error
        );
        assert!(bounds.axis_error <= 1e-11, "{}", bounds.axis_error);
        for s in [lo, hi] {
            let frame = match reading {
                FrameReading::Section => read(s),
                _ => read(s).plan().unwrap(),
            };
            let exact = wanted(s);
            let off = (frame.point - Point3::ZERO - exact.point).length();
            assert!(
                off <= bounds.point_error,
                "{off} > {} at {s}",
                bounds.point_error
            );
            for (got, want) in [
                (frame.tangent, exact.tangent),
                (frame.lateral, exact.lateral),
                (frame.up, exact.up),
            ] {
                assert!((got - want).length() <= bounds.axis_error + 1e-15);
            }
        }
    }
}

// Spans inside the plan's pieces (seams at 40, 100, 180, 240) and the
// profile's (120).
const SPANS: [(Scalar, Scalar); 7] = [
    (5.0, 35.0),
    (45.0, 95.0),
    (105.0, 118.0),
    (125.0, 175.0),
    (185.0, 235.0),
    (245.0, 275.0),
    (60.0, 61.0),
];

#[test]
fn a_clothoid_basis_is_bounded() {
    let c = 1.0 / 1800.0;
    let curve = Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::clothoid(0.0, 60.0 * c, 60.0),
        60.0,
    ));
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    let reference = |s: Scalar| {
        let f = frame();
        let local = fresnel(0.5 * c, s);
        let p = f.origin + f.x * local.x + f.y * local.y;
        let t = f.x * (0.5 * c * s * s).cos() + f.y * (0.5 * c * s * s).sin();
        planar(Vec3::new(p.x, p.y, 0.0), Vec3::new(t.x, t.y, 0.0))
    };
    for (lo, hi) in [(2.0, 20.0), (30.0, 58.0), (10.0, 10.5)] {
        check(
            &law,
            reference,
            |s| station_section2(&curve, s).unwrap(),
            lo,
            hi,
        );
    }
}

/// Composite Simpson, 2^16 panels, of `(cos theta, sin theta)` over `[0, s]`.
fn simpson(theta: &dyn Fn(Scalar) -> Scalar, s: Scalar) -> Vec2 {
    let n = 1usize << 16;
    let h = s / n as Scalar;
    let (mut x, mut y) = (0.0, 0.0);
    for i in 0..=n {
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let u = h * i as Scalar;
        x += w * theta(u).cos();
        y += w * theta(u).sin();
    }
    Vec2::new(x * h / 3.0, y * h / 3.0)
}

/// A cosine ramp and a sine-corrected transition (harmonic laws) are
/// bounded, their points from a Simpson reference.
#[test]
fn harmonic_transitions_are_bounded() {
    let (l, kappa) = (60.0, 1.0 / 100.0);
    let pi = core::f64::consts::PI;
    let tau = core::f64::consts::TAU;
    let cosine = CurvatureLaw::Sinusoid {
        mean: 0.5 * kappa,
        amplitude: -0.5 * kappa,
        angular_frequency: pi / l,
        phase: core::f64::consts::FRAC_PI_2,
    };
    let corrected = CurvatureLaw::sine_corrected_transition(0.0, kappa, l);
    let cosine_theta = move |u: Scalar| 0.5 * kappa * (u - l / pi * (pi * u / l).sin());
    let corrected_theta = move |u: Scalar| {
        kappa * u * u / (2.0 * l) + kappa * l / (tau * tau) * ((tau * u / l).cos() - 1.0)
    };
    let cases: [(CurvatureLaw, &dyn Fn(Scalar) -> Scalar); 2] =
        [(cosine, &cosine_theta), (corrected, &corrected_theta)];
    for (law, theta) in cases {
        let curve = Curve2::Intrinsic(Intrinsic2::new(frame(), law, l));
        let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
        let reference = |s: Scalar| {
            let f = frame();
            let local = simpson(theta, s);
            let p = f.origin + f.x * local.x + f.y * local.y;
            let t = f.x * theta(s).cos() + f.y * theta(s).sin();
            planar(Vec3::new(p.x, p.y, 0.0), Vec3::new(t.x, t.y, 0.0))
        };
        for (lo, hi) in [(2.0, 20.0), (25.0, 58.0)] {
            check(
                &law,
                reference,
                |s| station_section2(&curve, s).unwrap(),
                lo,
                hi,
            );
        }
    }
}

#[test]
fn a_level_space_curve_is_bounded() {
    let c = 1.0 / 1800.0;
    let start = Frame3 {
        origin: Point3::new(3.0, -1.0, 4.0),
        x: Vec3::new(0.6, 0.8, 0.0),
        y: Vec3::new(-0.8, 0.6, 0.0),
        z: Vec3::Z,
    };
    let curve = Curve3::Intrinsic(Intrinsic3::new(
        start,
        CurvatureLaw::clothoid(0.0, 60.0 * c, 60.0),
        CurvatureLaw::Constant { curvature: 0.0 },
        60.0,
    ));
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let reference = |s: Scalar| {
        let local = fresnel(0.5 * c, s);
        let p = start.origin + start.x * local.x + start.y * local.y;
        let t = start.x * (0.5 * c * s * s).cos() + start.y * (0.5 * c * s * s).sin();
        planar(p, t)
    };
    for (lo, hi) in [(2.0, 20.0), (30.0, 58.0)] {
        check(
            &law,
            reference,
            |s| station_section3(&curve, s).unwrap(),
            lo,
            hi,
        );
    }
}

#[test]
fn a_gradient_curve_over_the_alignment_is_bounded() {
    let a = alignment();
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Intrinsic(a.curve.clone()),
        profile(),
    ));
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    for (lo, hi) in SPANS {
        check(
            &law,
            |d| elevated_at(&a, d),
            |d| station_section3(&curve, d).unwrap(),
            lo,
            hi,
        );
    }
    // A plan seam inside, and the profile's seam inside: split by the caller.
    for (lo, hi) in [(30.0, 50.0), (110.0, 130.0), (230.0, 250.0)] {
        assert!(law.bounds(lo, hi, FrameReading::Section).is_none());
    }
    // The plan alone, the same alignment.
    let plan = Curve2::Intrinsic(a.curve.clone());
    let law = FrameLaw::atomic(StationCurve::Two(&plan), None, TOLERANCE).unwrap();
    for (lo, hi) in [(5.0, 35.0), (45.0, 95.0), (105.0, 175.0), (245.0, 275.0)] {
        check(
            &law,
            |d| plan_at(&a, d),
            |d| station_section2(&plan, d).unwrap(),
            lo,
            hi,
        );
    }
}

/// The cant over the alignment, in closed form: a ramp on the first line
/// (so a roll forgotten there is seen), a second ramp across the first
/// clothoid, a plateau on the arc, and the way back.
const CANT_B: Scalar = 1.5;

fn cant_pieces(scale: Scalar) -> Vec<CantPiece> {
    vec![
        CantPiece::linear(40.0, 0.0, 0.05 * scale),
        CantPiece::linear(60.0, 0.05 * scale, 0.1 * scale),
        CantPiece::constant(80.0, 0.1 * scale),
        CantPiece::linear(60.0, 0.1 * scale, 0.0),
        CantPiece::constant(40.0, 0.0),
    ]
}

/// The cant's value and rate at `d`, from the pieces above.
fn cant_at(d: Scalar) -> (Scalar, Scalar) {
    match d {
        d if d < 40.0 => (0.05 * d / 40.0, 0.05 / 40.0),
        d if d < 100.0 => (0.05 + 0.05 * (d - 40.0) / 60.0, 0.05 / 60.0),
        d if d < 180.0 => (0.1, 0.0),
        d if d < 240.0 => (0.1 * (240.0 - d) / 60.0, -0.1 / 60.0),
        _ => (0.0, 0.0),
    }
}

fn canted(convention: BankConvention, pivot: bool) -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(Curve2::Intrinsic(alignment().curve), profile()),
        CantLaw::new(cant_pieces(1.0)),
        if pivot {
            // Rotation about the low rail: the point rises with half the cant.
            CantLaw::new(cant_pieces(0.5))
        } else {
            CantLaw::zero(280.0)
        },
        CANT_B,
        convention,
    ))
}

/// The rolled section frame at `d`, built as the section frames are
/// defined: the point path's tangent `(p, G)`, the horizontal left normal
/// `n`, `u = t x n`, rolled by `rho` from the cant.
fn canted_at(a: &Alignment, d: Scalar, convention: BankConvention, pivot: bool) -> Reference {
    let (p, t) = a.world(d);
    let (z, g) = profile_at(d);
    let (cant, rate) = cant_at(d);
    let (e, e_rate) = if pivot {
        (0.5 * cant, 0.5 * rate)
    } else {
        (0.0, 0.0)
    };
    let slope = g + e_rate;
    let derivative = Vec3::new(t.x, t.y, slope);
    let tangent = derivative / derivative.length();
    let normal = Vec3::new(-t.y, t.x, 0.0);
    let square_up = tangent.cross(normal);
    let sine = match convention {
        BankConvention::TangentRotation => cant / CANT_B,
        _ => cant * (1.0 + slope * slope).sqrt() / CANT_B,
    };
    let roll = sine.asin();
    let (sin, cos) = roll.sin_cos();
    Reference {
        point: Vec3::new(p.x, p.y, z + e),
        tangent,
        lateral: cos * normal + sin * square_up,
        up: -sin * normal + cos * square_up,
    }
}

#[test]
fn a_canted_curve_over_the_alignment_is_bounded() {
    let a = alignment();
    for (convention, pivot) in [
        (BankConvention::TangentRotation, false),
        (BankConvention::VerticalRise, false),
        (BankConvention::TangentRotation, true),
    ] {
        let curve = canted(convention, pivot);
        let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
        for (lo, hi) in SPANS {
            check(
                &law,
                |d| canted_at(&a, d, convention, pivot),
                |d| station_section3(&curve, d).unwrap(),
                lo,
                hi,
            );
        }
        // A plan seam, the profile's and the cant's inside: split by the
        // caller.
        for (lo, hi) in [(30.0, 50.0), (110.0, 130.0), (230.0, 250.0)] {
            assert!(law.bounds(lo, hi, FrameReading::Section).is_none());
        }
    }
}

/// A Viennese bend of the bank angle (0 to 0.1 rad over 40..140) with the
/// pivot held about the right rail, in closed form: `psi`, `D = b sin(psi)`
/// and the pivot `e = D / 2`, `e' = D' / 2`.
fn viennese_cant(d: Scalar) -> (Scalar, Scalar) {
    let (psi, rate) = if d <= 40.0 {
        (0.0, 0.0)
    } else if d >= 140.0 {
        (0.1, 0.0)
    } else {
        let xi: Scalar = (d - 40.0) / 100.0;
        let shape = xi.powi(4) * (35.0 - 84.0 * xi + 70.0 * xi * xi - 20.0 * xi.powi(3));
        (0.1 * shape, 0.1 * 140.0 * (xi * (1.0 - xi)).powi(3) / 100.0)
    };
    (CANT_B * psi.sin(), CANT_B * psi.cos() * rate)
}

fn viennese_curve() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(Curve2::Intrinsic(alignment().curve), profile()),
        CantLaw::new(vec![
            // Angle pieces throughout: a constant piece is a height.
            CantPiece::viennese_bend(40.0, 0.0, 0.0),
            CantPiece::viennese_bend(100.0, 0.0, 0.1),
            CantPiece::viennese_bend(140.0, 0.1, 0.1),
        ]),
        CantLaw::new(vec![CantPiece::about_rail(280.0, RailSide::Right, 0.0)]),
        CANT_B,
        BankConvention::VerticalRise,
    ))
}

fn viennese_at(a: &Alignment, d: Scalar) -> Reference {
    let (p, t) = a.world(d);
    let (z, g) = profile_at(d);
    let (cant, rate) = viennese_cant(d);
    let slope = g + 0.5 * rate;
    let tangent = Vec3::new(t.x, t.y, slope) / (1.0 + slope * slope).sqrt();
    let normal = Vec3::new(-t.y, t.x, 0.0);
    let square_up = tangent.cross(normal);
    let sine = cant * (1.0 + slope * slope).sqrt() / CANT_B;
    let (sin, cos) = sine.asin().sin_cos();
    Reference {
        point: Vec3::new(p.x, p.y, z + 0.5 * cant),
        tangent,
        lateral: cos * normal + sin * square_up,
        up: -sin * normal + cos * square_up,
    }
}

/// A Viennese bend (an angle piece) and a held-rail pivot derived from it,
/// under the vertical-rise convention, over the alignment (#396).
#[test]
fn a_viennese_bend_with_a_held_rail_pivot_over_the_alignment_is_bounded() {
    let a = alignment();
    let curve = viennese_curve();
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    // Inside the plan's, the profile's and the cant's pieces (seams 40, 100,
    // 120, 140, 180, 240).
    for (lo, hi) in [
        (5.0, 35.0),
        (45.0, 95.0),
        (105.0, 118.0),
        (122.0, 138.0),
        (145.0, 175.0),
        (185.0, 235.0),
        (245.0, 275.0),
    ] {
        check(
            &law,
            |d| viennese_at(&a, d),
            |d| station_section3(&curve, d).unwrap(),
            lo,
            hi,
        );
    }
}

fn clothoid_plan() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
        60.0,
    ))
}

/// Every family left unbounded by name stays so, and `piece_frame_bounds`
/// bounds no intrinsic curve at all: only the reader, which compares the
/// frames it reads, may certify them.
#[test]
fn what_is_not_certified_stays_unbounded() {
    let plan = clothoid_plan();
    let banked = Curve3::Banked(Banked3::new(
        Elevated3::new(plan.clone(), ElevationLaw::level(0.0)),
        CantLaw::new(vec![CantPiece::constant(60.0, 0.15)]),
        CantLaw::zero(60.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    let chain_plan = Curve2::Chain(Chain2::new(
        frame(),
        vec![ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
            length: 60.0,
        }],
    ));
    let banked_chain = Curve3::Banked(Banked3::new(
        Elevated3::new(chain_plan, ElevationLaw::level(0.0)),
        CantLaw::new(vec![CantPiece::constant(60.0, 0.15)]),
        CantLaw::zero(60.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    let chain = Curve2::Chain(Chain2::new(
        frame(),
        vec![ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
            length: 60.0,
        }],
    ));
    let stretched = Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::new(2.0, 0.0),
            y: Vec2::new(0.0, 2.0),
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
        60.0,
    ));
    let stretched_plan =
        Curve3::Elevated(Elevated3::new(stretched.clone(), ElevationLaw::level(0.0)));
    let law = CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0);
    let flat = CurvatureLaw::Constant { curvature: 0.0 };
    let level = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::new(0.6, 0.8, 0.0),
        y: Vec3::new(-0.8, 0.6, 0.0),
        z: Vec3::Z,
    };
    let tilted_frame = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::new(0.6, 0.64, 0.48),
        y: Vec3::new(-0.8, 0.48, 0.36),
        z: Vec3::new(0.0, -0.6, 0.8),
    };
    let tilted = Curve3::Intrinsic(Intrinsic3::new(
        tilted_frame,
        law.clone(),
        flat.clone(),
        60.0,
    ));
    let twisted = Curve3::Intrinsic(Intrinsic3::new(
        level,
        law.clone(),
        CurvatureLaw::Constant { curvature: 0.01 },
        60.0,
    ));
    // Tilted by a hair: its plane projection is orthonormal to 1e-12, but
    // the curve is not level.
    let hair = Frame3 {
        origin: Point3::ZERO,
        x: Vec3::new(0.6, 0.8, 1e-9),
        y: Vec3::new(-0.8, 0.6, 0.0),
        z: Vec3::new(-6e-10, -8e-10, 1.0),
    };
    let barely = Curve3::Intrinsic(Intrinsic3::new(hair, law.clone(), flat.clone(), 60.0));
    let cases: [(&str, StationCurve<'_>); 7] = [
        ("banked over a chain", StationCurve::Three(&banked_chain)),
        ("chain", StationCurve::Two(&chain)),
        ("stretched", StationCurve::Two(&stretched)),
        ("stretched, elevated", StationCurve::Three(&stretched_plan)),
        ("tilted", StationCurve::Three(&tilted)),
        ("tilted by a hair", StationCurve::Three(&barely)),
        ("twisted", StationCurve::Three(&twisted)),
    ];
    for (name, curve) in cases {
        let law = FrameLaw::atomic(curve, None, TOLERANCE).unwrap();
        for reading in [FrameReading::Section, FrameReading::Plan] {
            assert!(law.bounds(10.0, 20.0, reading).is_none(), "{name}");
        }
    }
    // The public closed-form bounds never certify an intrinsic curve.
    let elevated = Curve3::Elevated(Elevated3::new(plan.clone(), profile()));
    let level_curve = Curve3::Intrinsic(Intrinsic3::new(level, law, flat, 60.0));
    let families: [(&str, StationCurve<'_>); 5] = [
        ("plan", StationCurve::Two(&plan)),
        ("level", StationCurve::Three(&level_curve)),
        ("elevated", StationCurve::Three(&elevated)),
        ("banked", StationCurve::Three(&banked)),
        ("tilted", StationCurve::Three(&tilted)),
    ];
    for (name, curve) in families {
        let piece = StationPiece::between(curve, 10.0, 20.0).unwrap();
        for reading in [FrameReading::Section, FrameReading::Plan] {
            assert!(
                piece_frame_bounds(&piece, 0.0, 10.0, reading).is_none(),
                "{name}"
            );
        }
    }
}
