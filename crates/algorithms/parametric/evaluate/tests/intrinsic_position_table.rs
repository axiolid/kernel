//! The certified position table of intrinsic curves (#361).
//!
//! Every expectation is independent of the Gauss-Legendre rule under test
//! (ADR 0060): an arc in closed form, a clothoid from the Fresnel
//! integrals' power series, a line, clothoid, arc, clothoid, line composed
//! from those, and a Bloss transition against a composite Simpson rule.
//! Each returned point must lie within its reported error of the exact one,
//! at loose and tight targets, and the error stay within the target.

use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{CurvatureLaw, Intrinsic2};
use axiolid_evaluate::bound::IntrinsicPositionTable;

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    }
}

fn place(frame: &Frame2, local: Vec2) -> Point2 {
    frame.origin + frame.x * local.x + frame.y * local.y
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

/// One piece of an alignment: a line, a clothoid rising from zero to
/// `1/radius`, an arc of `radius`, or a clothoid falling back to zero.
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

    /// The heading turned over the whole piece.
    fn turn(self) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Rising(l, r) | Piece::Falling(l, r) => 0.5 * l / r,
            Piece::Arc(l, r) => l / r,
        }
    }

    /// The displacement over its first `s`, starting along `+x`.
    fn displacement(self, s: Scalar) -> Vec2 {
        match self {
            Piece::Line(_) => Vec2::new(s, 0.0),
            Piece::Rising(l, r) => fresnel(0.5 / (r * l), s),
            Piece::Arc(_, r) => Vec2::new(r * (s / r).sin(), r * (1.0 - (s / r).cos())),
            Piece::Falling(l, r) => {
                // The heading is the rising clothoid's reversed, from the
                // other end: `theta_f(s) = turn - theta_r(l - s)`.
                let a = 0.5 / (r * l);
                let tail = fresnel(a, l) - fresnel(a, l - s);
                rotate(self.turn(), Vec2::new(tail.x, -tail.y))
            }
        }
    }
}

/// An alignment of pieces in closed form: the curve and its exact point.
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

    fn exact(&self, s: Scalar) -> Point2 {
        let (mut at, mut heading, mut begin) = (Vec2::ZERO, 0.0, 0.0);
        for piece in &self.pieces {
            let inside = (s - begin).clamp(0.0, piece.length());
            at += rotate(heading, piece.displacement(inside));
            if s <= begin + piece.length() {
                break;
            }
            heading += piece.turn();
            begin += piece.length();
        }
        place(&frame(), Vec2::new(at.x, at.y))
    }
}

fn gentle() -> Alignment {
    Alignment::new(vec![
        Piece::Line(40.0),
        Piece::Rising(60.0, 300.0),
        Piece::Arc(80.0, 300.0),
        Piece::Falling(60.0, 300.0),
        Piece::Line(40.0),
    ])
}

fn samples(length: Scalar, count: usize) -> Vec<Scalar> {
    (0..=count)
        .map(|i| length * i as Scalar / count as Scalar)
        .collect()
}

#[test]
fn an_arc_is_within_its_certified_error() {
    let radius = 12.0;
    let arc = Intrinsic2::new(
        frame(),
        CurvatureLaw::Constant {
            curvature: 1.0 / radius,
        },
        40.0,
    );
    for target in [1e-3, 1e-6, 1e-9] {
        let table = IntrinsicPositionTable::new2(&arc, target, 40.0).unwrap();
        for s in samples(40.0, 97) {
            let (p, error) = table.point2(s).unwrap();
            let a = s / radius;
            let exact = place(
                &frame(),
                Vec2::new(radius * a.sin(), radius * (1.0 - a.cos())),
            );
            assert!((p - exact).length() <= error, "{target}: {s}");
            assert!(error <= target, "{error} > {target}");
        }
    }
}

#[test]
fn a_clothoid_is_within_its_certified_error() {
    let c = 1.0 / 1800.0;
    let spiral = Intrinsic2::new(frame(), CurvatureLaw::clothoid(0.0, 60.0 * c, 60.0), 60.0);
    for target in [1e-3, 1e-6, 1e-9] {
        let table = IntrinsicPositionTable::new2(&spiral, target, 60.0).unwrap();
        for s in samples(60.0, 61) {
            let (p, error) = table.point2(s).unwrap();
            let exact = place(&frame(), fresnel(0.5 * c, s));
            assert!((p - exact).length() <= error, "{target}: {s}");
            assert!(error <= target);
        }
    }
}

#[test]
fn a_composed_alignment_is_within_its_certified_error() {
    let alignment = gentle();
    for target in [1e-3, 1e-6, 1e-9] {
        let table = IntrinsicPositionTable::new2(&alignment.curve, target, 280.0).unwrap();
        let mut seen = samples(280.0, 233);
        // The seams, and a hair either side.
        for seam in [40.0, 100.0, 180.0, 240.0] {
            seen.extend([seam - 1e-9, seam, seam + 1e-9]);
        }
        for s in seen {
            let (p, error) = table.point2(s).unwrap();
            let off = (p - alignment.exact(s)).length();
            assert!(off <= error, "{target}: {s}: {off} > {error}");
            assert!(error <= target, "{error} > {target}");
        }
    }
}

#[test]
fn a_table_over_part_of_a_curve_reads_only_that_part() {
    let alignment = gentle();
    let table = IntrinsicPositionTable::new2(&alignment.curve, 1e-6, 150.0).unwrap();
    let (p, error) = table.point2(150.0).unwrap();
    assert!((p - alignment.exact(150.0)).length() <= error);
    assert!(table.point2(150.5).is_none());
    assert!(table.point2(-1e-3).is_none());
    assert!(table.point2(Scalar::NAN).is_none());
}

#[test]
fn a_bloss_transition_is_within_its_certified_error() {
    // k = K (3 (s/L)^2 - 2 (s/L)^3), theta = K (s^3 / L^2 - s^4 / (2 L^3)).
    let (big_k, l) = (1.0 / 250.0, 90.0);
    let law = CurvatureLaw::Polynomial {
        coefficients: vec![0.0, 0.0, 3.0 * big_k / (l * l), -2.0 * big_k / (l * l * l)],
    };
    let curve = Intrinsic2::new(frame(), law, l);
    let theta = |u: Scalar| big_k * (u.powi(3) / (l * l) - u.powi(4) / (2.0 * l * l * l));
    // Composite Simpson, 2^18 panels.
    let simpson = |s: Scalar| {
        let n = 1usize << 18;
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
    };
    for target in [1e-3, 1e-6, 1e-9] {
        let table = IntrinsicPositionTable::new2(&curve, target, l).unwrap();
        for s in samples(l, 12) {
            let (p, error) = table.point2(s).unwrap();
            let exact = place(&frame(), simpson(s));
            assert!((p - exact).length() <= error, "{target}: {s}");
            assert!(error <= target);
        }
    }
}

#[test]
fn reads_are_bitwise_equal_in_any_order() {
    let alignment = gentle();
    let table = IntrinsicPositionTable::new2(&alignment.curve, 1e-9, 280.0).unwrap();
    let at = samples(280.0, 300);
    let read = |s: Scalar| {
        let (p, e) = table.point2(s).unwrap();
        [p.x.to_bits(), p.y.to_bits(), e.to_bits()]
    };
    let forward: Vec<_> = at.iter().map(|&s| read(s)).collect();
    let backward: Vec<_> = at.iter().rev().map(|&s| read(s)).collect();
    for (i, bits) in forward.iter().enumerate() {
        assert_eq!(bits, &backward[at.len() - 1 - i]);
    }
    // A fixed pseudo-random order.
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    for _ in 0..600 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let i = (state >> 33) as usize % at.len();
        assert_eq!(read(at[i]), forward[i]);
    }
}

#[test]
fn a_planar_space_curve_is_its_plane_curve() {
    use axiolid_core::{Frame3, Point3, Vec3};
    use axiolid_curve::Intrinsic3;
    let c = 1.0 / 1800.0;
    let z = Vec3::new(1.0, -2.0, 2.0).normalize();
    let x = Vec3::new(2.0, 1.0, 0.0).normalize();
    let start = Frame3 {
        origin: Point3::new(5.0, -3.0, 1.5),
        x,
        y: z.cross(x),
        z,
    };
    let law = CurvatureLaw::clothoid(0.0, 60.0 * c, 60.0);
    let flat = CurvatureLaw::Constant { curvature: 0.0 };
    let space = Intrinsic3::new(start, law.clone(), flat, 60.0);
    let table = IntrinsicPositionTable::new3(&space, 1e-6, 60.0).unwrap();
    for s in samples(60.0, 40) {
        let (p, error) = table.point3(s).unwrap();
        let local = fresnel(0.5 * c, s);
        let exact = start.origin + start.x * local.x + start.y * local.y;
        assert!((p - exact).length() <= error);
        assert!(error <= 1e-6);
    }
    assert!(table.point2(1.0).is_none());
    let twisted = Intrinsic3::new(start, law, CurvatureLaw::Constant { curvature: 0.01 }, 60.0);
    assert!(IntrinsicPositionTable::new3(&twisted, 1e-6, 60.0).is_none());
}

#[test]
fn an_over_cap_table_and_a_bad_end_are_refused() {
    let tight = Intrinsic2::new(frame(), CurvatureLaw::Constant { curvature: 1.0 }, 1.0e6);
    assert!(IntrinsicPositionTable::new2(&tight, 1e-12, 1.0e6).is_none());
    let alignment = gentle();
    let length = alignment.curve.length;
    assert!(IntrinsicPositionTable::new2(&alignment.curve, 1e-6, length + 1.0).is_none());
    assert!(IntrinsicPositionTable::new2(&alignment.curve, 1e-6, 0.0).is_none());
    assert!(IntrinsicPositionTable::new2(&alignment.curve, 0.0, 10.0).is_none());
    assert!(IntrinsicPositionTable::new2(&alignment.curve, 1e-6, Scalar::NAN).is_none());
}

#[test]
fn a_partial_panel_is_charged_its_own_remainder() {
    // One panel covers the whole arc at this target, so a point inside it
    // is all partial panel: the rule's error there, about `w^5 k^4 / 4320`,
    // is far above the rounding and must be in the reported error.
    let radius = 100.0;
    let arc = Intrinsic2::new(
        frame(),
        CurvatureLaw::Constant {
            curvature: 1.0 / radius,
        },
        10.0,
    );
    let table = IntrinsicPositionTable::new2(&arc, 1e-3, 10.0).unwrap();
    assert_eq!(table.panels(), 1);
    let mut charged: Scalar = 0.0;
    for s in samples(10.0, 20) {
        let (p, error) = table.point2(s).unwrap();
        let a = s / radius;
        let exact = place(
            &frame(),
            Vec2::new(radius * a.sin(), radius * (1.0 - a.cos())),
        );
        let off = (p - exact).length();
        assert!(off <= error, "{s}: {off} > {error}");
        charged = charged.max(off);
    }
    // The rule really errs here, well past rounding.
    assert!(charged > 1e-10, "{charged}");
}

#[test]
fn the_reported_error_never_falls_along_the_curve() {
    // Checkpoints carry the remainder so far: an error that fell at one
    // would be a remainder forgotten.
    let alignment = gentle();
    let table = IntrinsicPositionTable::new2(&alignment.curve, 1e-9, 280.0).unwrap();
    assert!(table.panels() > 128, "{}", table.panels());
    let mut last = 0.0;
    for s in samples(280.0, 4000) {
        let (_, error) = table.point2(s).unwrap();
        assert!(error >= last, "{s}: {error} < {last}");
        last = error;
    }
}

#[test]
fn a_ten_kilometre_alignment_needs_few_panels() {
    let alignment = Alignment::new(vec![
        Piece::Line(2000.0),
        Piece::Rising(300.0, 1000.0),
        Piece::Arc(3000.0, 1000.0),
        Piece::Falling(300.0, 1000.0),
        Piece::Line(4400.0),
    ]);
    let table = IntrinsicPositionTable::new2(&alignment.curve, 2.5e-7, 10_000.0).unwrap();
    assert!(table.panels() < 200_000, "{}", table.panels());
    for s in samples(10_000.0, 211) {
        let (p, error) = table.point2(s).unwrap();
        let off = (p - alignment.exact(s)).length();
        assert!(off <= error, "{s}: {off} > {error}");
        assert!(error <= 2.5e-7);
    }
}
