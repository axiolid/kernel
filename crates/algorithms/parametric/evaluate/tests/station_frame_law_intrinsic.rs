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
    Elevated3, ElevationLaw, Intrinsic2, Intrinsic3,
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
            assert!(off <= bounds.point_error, "{off} > {}", bounds.point_error);
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
    let cases: [(&str, StationCurve<'_>); 6] = [
        ("banked over a clothoid", StationCurve::Three(&banked)),
        ("chain", StationCurve::Two(&chain)),
        ("stretched", StationCurve::Two(&stretched)),
        ("stretched, elevated", StationCurve::Three(&stretched_plan)),
        ("tilted", StationCurve::Three(&tilted)),
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
