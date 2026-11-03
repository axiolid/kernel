//! The error of one station's frame against the exact frame law (#405):
//! `FrameLaw::error_at` on circles, polylines, intrinsic plans alone and
//! under a profile, and composites of them.
//!
//! Every reference is independent of the kernel's quadratures: circles and
//! polylines in closed form, the alignment's points from the Fresnel
//! integrals' power series and closed forms (as in
//! `station_frame_law_intrinsic.rs`), the ellipse's perimeter by a fine
//! Simpson rule. The kernel's read frame must lie within the error of the
//! exact frame at the distance asked for and at the distances its
//! resolution may have put it at.

use std::f64::consts::{PI, TAU};

use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, Circle2,
    Circle3, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Ellipse2, Intrinsic2, KnotSpec,
    Line2, OffsetLaw, Polyline2, SeamSide,
};
use axiolid_evaluate::arc_parameter::ARC_LENGTH_TOLERANCE;
use axiolid_evaluate::station::{
    offset_pieces, station_section2_on, station_section3_on, CompositeBasis, FrameLaw,
    FrameReading, SectionFrame, StationCurve, StationFrameError, StationPiece,
};

const TOLERANCE: Scalar = 1e-6;
const SIDES: [SeamSide; 2] = [SeamSide::Outgoing, SeamSide::Incoming];
const READINGS: [FrameReading; 2] = [FrameReading::Section, FrameReading::Plan];

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    }
}

/// An exact frame, built as the section frames are defined.
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

/// The error covers the distance of the read frame from the exact one.
fn covers(error: StationFrameError, read: &SectionFrame, exact: Reference, what: &str) {
    let point = (read.point - Point3::ZERO - exact.point).length();
    assert!(point <= error.point, "{what}: point {point} > {error:?}");
    for (name, got, want) in [
        ("tangent", read.tangent, exact.tangent),
        ("lateral", read.lateral, exact.lateral),
        ("up", read.up, exact.up),
    ] {
        let off = (got - want).length();
        assert!(off <= error.axis, "{what}: {name} {off} > {error:?}");
    }
}

fn read_as(frame: SectionFrame, reading: FrameReading) -> SectionFrame {
    match reading {
        FrameReading::Section => frame,
        _ => frame.plan().unwrap(),
    }
}

fn exact_as(exact: Reference, reading: FrameReading) -> Reference {
    match reading {
        FrameReading::Section => exact,
        _ => plan_reading(exact),
    }
}

fn arc_slack(distance: Scalar) -> Scalar {
    ARC_LENGTH_TOLERANCE * distance.abs().max(1.0)
}

#[test]
fn a_circle_station_is_charged_its_distance_resolution() {
    let (r, s) = (30.0, 150.0);
    let f = frame();
    let circle = Curve2::Circle(Circle2 {
        frame: f,
        radius: r,
    });
    let exact = |s: Scalar| {
        let (sin, cos) = (s / r).sin_cos();
        let p = f.origin + r * (f.x * cos + f.y * sin);
        let t = f.y * cos - f.x * sin;
        planar(Vec3::new(p.x, p.y, 0.0), Vec3::new(t.x, t.y, 0.0))
    };
    let law = FrameLaw::atomic(StationCurve::Two(&circle), None, TOLERANCE).unwrap();
    let delta = 2.0 * arc_slack(s);
    for reading in READINGS {
        let error = law.error_at(s, SeamSide::Outgoing, reading).unwrap();
        let read = read_as(
            station_section2_on(&circle, s, SeamSide::Outgoing).unwrap(),
            reading,
        );
        for at in [s - delta, s, s + delta] {
            covers(error, &read, exact_as(exact(at), reading), "circle");
        }
        // Not looser than the resolution and the rounding.
        assert!(
            error.point >= delta && error.point <= 4.0 * delta,
            "{error:?}"
        );
        assert!(
            error.axis >= delta / r && error.axis <= 4.0 * delta / r,
            "{error:?}"
        );
    }

    // A tilted circle in space, the reference-up frame.
    let (r, s) = (5.0, 7.0);
    let x = Vec3::new(1.0, 0.0, 0.0);
    let y = Vec3::new(0.0, 0.8, 0.6);
    let tilted = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(1.0, 2.0, 3.0),
            x,
            y,
            z: x.cross(y),
        },
        radius: r,
    });
    let exact = |s: Scalar| {
        let (sin, cos) = (s / r).sin_cos();
        let p = Vec3::new(1.0, 2.0, 3.0) + r * (x * cos + y * sin);
        reference_up(p, y * cos - x * sin)
    };
    let law = FrameLaw::atomic(StationCurve::Three(&tilted), None, TOLERANCE).unwrap();
    let delta = 2.0 * arc_slack(s);
    for reading in READINGS {
        let error = law.error_at(s, SeamSide::Outgoing, reading).unwrap();
        let read = read_as(
            station_section3_on(&tilted, s, SeamSide::Outgoing).unwrap(),
            reading,
        );
        for at in [s - delta, s, s + delta] {
            covers(error, &read, exact_as(exact(at), reading), "tilted circle");
        }
    }
}

/// `int_0^s (cos(a u^2), sin(a u^2)) du` by its power series.
fn fresnel(a: Scalar, s: Scalar) -> Vec2 {
    let (mut x, mut y) = (0.0, 0.0);
    // `(a s^2)^n / n!`, built term by term: the powers alone overflow.
    let mut term = 1.0;
    for n in 0..120 {
        if n > 0 {
            term *= a * s * s / Scalar::from(n);
        }
        let value = s * term / Scalar::from(2 * n + 1);
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
    /// A clothoid with curvature `c s`.
    Ramp(Scalar, Scalar),
    Rising(Scalar, Scalar),
    Arc(Scalar, Scalar),
    Falling(Scalar, Scalar),
}

impl Piece {
    fn length(self) -> Scalar {
        match self {
            Piece::Line(l)
            | Piece::Ramp(l, _)
            | Piece::Rising(l, _)
            | Piece::Arc(l, _)
            | Piece::Falling(l, _) => l,
        }
    }

    fn law(self) -> CurvatureLaw {
        match self {
            Piece::Line(_) => CurvatureLaw::Constant { curvature: 0.0 },
            Piece::Ramp(_, c) => CurvatureLaw::Polynomial {
                coefficients: vec![0.0, c],
            },
            Piece::Rising(l, r) => CurvatureLaw::clothoid(0.0, 1.0 / r, l),
            Piece::Arc(_, r) => CurvatureLaw::Constant { curvature: 1.0 / r },
            Piece::Falling(l, r) => CurvatureLaw::clothoid(1.0 / r, 0.0, l),
        }
    }

    fn turn(self) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Ramp(l, c) => 0.5 * c * l * l,
            Piece::Rising(l, r) | Piece::Falling(l, r) => 0.5 * l / r,
            Piece::Arc(l, r) => l / r,
        }
    }

    /// The heading turned over the first `s`.
    fn heading(self, s: Scalar) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Ramp(_, c) => 0.5 * c * s * s,
            Piece::Rising(l, r) => 0.5 * s * s / (r * l),
            Piece::Arc(_, r) => s / r,
            Piece::Falling(l, r) => self.turn() - 0.5 * (l - s) * (l - s) / (r * l),
        }
    }

    fn displacement(self, s: Scalar) -> Vec2 {
        match self {
            Piece::Line(_) => Vec2::new(s, 0.0),
            Piece::Ramp(_, c) => fresnel(0.5 * c, s),
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

    fn plan_at(&self, d: Scalar) -> Reference {
        let (p, t) = self.world(d);
        planar(Vec3::new(p.x, p.y, 0.0), Vec3::new(t.x, t.y, 0.0))
    }

    fn elevated_at(&self, d: Scalar) -> Reference {
        let (p, t) = self.world(d);
        let (z, g) = profile_at(d);
        let size = (1.0 + g * g).sqrt();
        reference_up(
            Vec3::new(p.x, p.y, z),
            Vec3::new(t.x / size, t.y / size, g / size),
        )
    }
}

/// Seams of the plan at 40, 100, 180 and 240.
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

#[test]
fn a_clothoid_station_is_charged_its_table() {
    let a = Alignment::new(vec![Piece::Ramp(120.0, 1.0 / 300.0)]);
    let plan = Curve2::Intrinsic(a.curve.clone());
    let law = FrameLaw::atomic(StationCurve::Two(&plan), None, TOLERANCE).unwrap();
    for reading in READINGS {
        let error = law.error_at(90.0, SeamSide::Outgoing, reading).unwrap();
        let read = read_as(
            station_section2_on(&plan, 90.0, SeamSide::Outgoing).unwrap(),
            reading,
        );
        covers(error, &read, exact_as(a.plan_at(90.0), reading), "clothoid");
        assert!(
            error.point <= 2.0 * TOLERANCE && error.axis <= 1e-11,
            "{error:?}"
        );
    }
}

#[test]
fn an_alignment_station_is_charged_on_and_off_its_seams() {
    let a = alignment();
    let plan = Curve2::Intrinsic(a.curve.clone());
    let elevated = Curve3::Elevated(Elevated3::new(plan.clone(), profile()));
    let under_plan = FrameLaw::atomic(StationCurve::Two(&plan), None, TOLERANCE).unwrap();
    let under_profile = FrameLaw::atomic(StationCurve::Three(&elevated), None, TOLERANCE).unwrap();
    // Inside a piece, on the plan's seams (40, 100, 180, 240) and on the
    // profile's (120).
    for d in [
        10.0, 70.0, 40.0, 100.0, 120.0, 150.0, 180.0, 200.0, 240.0, 270.0,
    ] {
        for side in SIDES {
            for reading in READINGS {
                let what = format!("alignment at {d} {side:?} {reading:?}");
                let error = under_plan.error_at(d, side, reading).expect(&what);
                let read = read_as(station_section2_on(&plan, d, side).unwrap(), reading);
                covers(error, &read, exact_as(a.plan_at(d), reading), &what);
                assert!(
                    error.point <= 2.0 * TOLERANCE && error.axis <= 1e-11,
                    "{what}"
                );

                let error = under_profile.error_at(d, side, reading).expect(&what);
                let read = read_as(station_section3_on(&elevated, d, side).unwrap(), reading);
                covers(error, &read, exact_as(a.elevated_at(d), reading), &what);
                assert!(
                    error.point <= 2.0 * TOLERANCE && error.axis <= 1e-11,
                    "{what}"
                );
            }
        }
    }

    // The same alignment as the one piece of a composite.
    let piece = StationPiece::whole(StationCurve::Three(&elevated)).unwrap();
    let basis = CompositeBasis::new(vec![piece]).unwrap();
    let law = basis.frame_law(TOLERANCE).unwrap();
    for d in [70.0, 120.0, 200.0] {
        let error = law
            .error_at(d, SeamSide::Outgoing, FrameReading::Section)
            .unwrap();
        let read = basis.section_on(d, SeamSide::Outgoing).unwrap();
        covers(error, &read, a.elevated_at(d), "composite alignment");
    }
}

#[test]
fn a_gradient_over_a_line_is_charged_its_rounding_alone() {
    // The plan line has no length: the frame is read on a span about the
    // station, and a plan distance is the parameter.
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::new(5.0, -2.0),
            direction: Vec2::X,
        }),
        ElevationLaw::constant_grade(100.0, 0.02),
    ));
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let size = (1.0_f64 + 0.02 * 0.02).sqrt();
    for d in [0.0, 0.5, 40.0, 1.0e5] {
        for reading in READINGS {
            let error = law.error_at(d, SeamSide::Outgoing, reading).unwrap();
            let read = read_as(
                station_section3_on(&curve, d, SeamSide::Outgoing).unwrap(),
                reading,
            );
            let exact = reference_up(
                Vec3::new(5.0 + d, -2.0, 100.0 + 0.02 * d),
                Vec3::new(1.0, 0.0, 0.02) / size,
            );
            covers(error, &read, exact_as(exact, reading), "gradient");
            assert!(error.axis <= 1e-12 && error.point <= 1e-8, "{error:?}");
        }
    }
    assert!(law
        .error_at(-1.0, SeamSide::Outgoing, FrameReading::Section)
        .is_none());
}

fn polyline() -> Curve2 {
    Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::ZERO,
            Point2::new(30.0, 0.0),
            Point2::new(30.0, 40.0),
            Point2::new(0.0, 40.0),
        ],
        closed: false,
    })
}

#[test]
fn a_polyline_station_is_exact_on_a_vertex_and_refused_beside_it() {
    let curve = polyline();
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    let (before, after) = (Vec3::X, Vec3::Y);
    for (side, tangent) in [(SeamSide::Incoming, before), (SeamSide::Outgoing, after)] {
        let error = law.error_at(30.0, side, FrameReading::Section).unwrap();
        let read = station_section2_on(&curve, 30.0, side).unwrap();
        let exact = planar(Vec3::new(30.0, 0.0, 0.0), tangent);
        covers(error, &read, exact, "vertex");
    }
    // Within the resolution of the vertex but not on it: which side the
    // exact frame stands on is not known.
    let slack = arc_slack(30.0);
    for d in [30.0 + 1.5 * slack, 30.0 - 1.5 * slack] {
        assert!(law
            .error_at(d, SeamSide::Outgoing, FrameReading::Section)
            .is_none());
    }
    // Clear of it by more than the span the resolution reads.
    let delta = 2.0 * slack;
    for (d, tangent, y) in [
        (30.0 + 3.0 * delta, after, 3.0 * delta),
        (30.0 - 3.0 * delta, before, 0.0),
    ] {
        let error = law
            .error_at(d, SeamSide::Outgoing, FrameReading::Section)
            .unwrap();
        let read = station_section2_on(&curve, d, SeamSide::Outgoing).unwrap();
        let exact = if tangent == after {
            planar(Vec3::new(30.0, y, 0.0), tangent)
        } else {
            planar(Vec3::new(d, 0.0, 0.0), tangent)
        };
        covers(error, &read, exact, "beside a vertex");
    }
}

/// The perimeter of the ellipse about `(30, 0)` with semi-axes 10 and 6,
/// by a fine Simpson rule.
fn perimeter() -> Scalar {
    let speed = |t: Scalar| (100.0 * t.sin().powi(2) + 36.0 * t.cos().powi(2)).sqrt();
    let n = 20_000;
    let h = TAU / n as Scalar;
    let mut sum = speed(0.0) + speed(TAU);
    for i in 1..n {
        sum += speed(i as Scalar * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    sum * h / 3.0
}

fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(30.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    })
}

#[test]
fn a_circle_after_an_ellipse_carries_its_joint() {
    let ellipse = ellipse();
    // Radius 5 about `(35, 0)`: it starts where the ellipse ends, heading
    // `+y`.
    let circle = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(35.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 5.0,
    });
    let basis = CompositeBasis::new(vec![
        StationPiece::whole(StationCurve::Two(&ellipse)).unwrap(),
        StationPiece::whole(StationCurve::Two(&circle)).unwrap(),
    ])
    .unwrap();
    let law = basis.frame_law(TOLERANCE).unwrap();
    let joint = basis.joint_intervals(TOLERANCE).unwrap()[0];
    let truth = perimeter();
    assert!(joint.lower <= truth && truth <= joint.upper);
    let d = truth + 4.0;
    let exact = |joint: Scalar| {
        let (sin, cos) = ((d - joint) / 5.0).sin_cos();
        let p = Vec3::new(35.0 + 5.0 * cos, 5.0 * sin, 0.0);
        planar(p, Vec3::new(-sin, cos, 0.0))
    };
    let error = law
        .error_at(d, SeamSide::Outgoing, FrameReading::Section)
        .unwrap();
    let read = basis.section_on(d, SeamSide::Outgoing).unwrap();
    // The joint truly where the ellipse ends, and anywhere certified.
    for at in [truth, joint.lower, joint.upper] {
        covers(error, &read, exact(at), "after an ellipse");
    }
    // The joint's uncertainty is charged, and no more than the tolerance
    // allows.
    let width = joint.upper - joint.lower;
    assert!(
        error.point >= 0.5 * width && error.point <= 4.0 * TOLERANCE,
        "{error:?}"
    );
    // A frame on the ellipse itself is read by quadrature: not charged.
    assert!(law
        .error_at(5.0, SeamSide::Outgoing, FrameReading::Section)
        .is_none());

    // On the joint the side named is read; beside it, within the joint's
    // uncertainty but not on it, the exact frame might stand on the
    // ellipse's side: refused. Clear of that, charged again.
    let start = basis.pieces()[0].length();
    let delta = law
        .error_at(start, SeamSide::Outgoing, FrameReading::Section)
        .expect("on the joint")
        .point;
    assert!(delta > 1e-9, "{delta}");
    for at in [0.5 * delta, 1.5 * delta] {
        assert!(
            law.error_at(start + at, SeamSide::Outgoing, FrameReading::Section)
                .is_none(),
            "{at}"
        );
    }
    assert!(law
        .error_at(
            start + 3.0 * delta,
            SeamSide::Outgoing,
            FrameReading::Section
        )
        .is_some());
}

#[test]
fn a_closed_form_read_is_charged_the_rounding_of_its_data() {
    // A circle of radius `r` about `(0, r)`, from `(0, 0)` heading `+x`:
    // the read `r - r cos a` cancels, rounding at the scale of `r`, not of
    // the point.
    for r in [5000.0, 1.0e6] {
        let curve = Curve3::Elevated(Elevated3::new(
            Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::new(0.0, r),
                    x: Vec2::new(0.0, -1.0),
                    y: Vec2::new(1.0, 0.0),
                },
                radius: r,
            }),
            ElevationLaw::level(0.0),
        ));
        let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
        let d = 0.84;
        let a = d / r;
        // Stable forms: `r sin a` and `2 r sin^2(a / 2)`.
        let point = Vec3::new(r * a.sin(), 2.0 * r * (0.5 * a).sin().powi(2), 0.0);
        let exact = reference_up(point, Vec3::new(a.cos(), a.sin(), 0.0));
        let error = law
            .error_at(d, SeamSide::Outgoing, FrameReading::Section)
            .unwrap();
        let read = station_section3_on(&curve, d, SeamSide::Outgoing).unwrap();
        covers(error, &read, exact, "a large circle");
        // The data's scale is charged: `1e-14 (r + r)`.
        assert!(error.point >= 1e-14 * r, "{error:?}");
    }
}

#[test]
fn a_vertex_after_an_ellipse_reads_the_side_named() {
    let ellipse = ellipse();
    // From where the ellipse ends: up 10 to `(40, 10)`, then left.
    let bent = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(40.0, 0.0),
            Point2::new(40.0, 10.0),
            Point2::new(30.0, 10.0),
        ],
        closed: false,
    });
    let basis = CompositeBasis::new(vec![
        StationPiece::whole(StationCurve::Two(&ellipse)).unwrap(),
        StationPiece::whole(StationCurve::Two(&bent)).unwrap(),
    ])
    .unwrap();
    let law = basis.frame_law(TOLERANCE).unwrap();
    // The joint's distance is uncertain by far more than the resolution of
    // a distance alone, so the span read about the vertex must stop at it.
    let corner = basis.pieces()[0].length() + 10.0;
    let vertex = Vec3::new(40.0, 10.0, 0.0);
    for (side, tangent) in [
        (SeamSide::Incoming, Vec3::Y),
        (SeamSide::Outgoing, -Vec3::X),
    ] {
        let error = law.error_at(corner, side, FrameReading::Section).unwrap();
        let read = basis.section_on(corner, side).unwrap();
        covers(error, &read, planar(vertex, tangent), "after an ellipse");
        assert!(
            error.point > 1e-9 && error.point <= 4.0 * TOLERANCE,
            "{error:?}"
        );
    }
}

#[test]
fn an_inexactly_placed_piece_is_not_charged() {
    let circle = Curve2::Circle(Circle2 {
        frame: frame(),
        radius: 30.0,
    });
    let rigid = Transform3::from_translation(Vec3::new(1.0, 2.0, 0.0));
    for (exact, charged) in [(true, true), (false, false)] {
        let piece = StationPiece::whole(StationCurve::Two(&circle))
            .unwrap()
            .placed(rigid, exact);
        let basis = CompositeBasis::new(vec![piece]).unwrap();
        let law = basis.frame_law(TOLERANCE).unwrap();
        let error = law.error_at(40.0, SeamSide::Outgoing, FrameReading::Section);
        assert_eq!(error.is_some(), charged, "placed exactly: {exact}");
    }
}

fn spline3() -> Curve3 {
    Curve3::BSpline(BSplineCurve {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(8.0, 2.0, 0.5),
            Point3::new(16.0, -3.0, 1.5),
            Point3::new(24.0, 2.0, 2.0),
            Point3::new(30.0, 6.0, 2.2),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

fn clothoid_plan() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
        60.0,
    ))
}

#[test]
fn what_is_not_certified_is_not_charged() {
    let ellipse = ellipse();
    let spline = spline3();
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
    let vertical = Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::intrinsic(0.0, 0.01, CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0)),
    ));
    let cases: [(&str, StationCurve<'_>); 5] = [
        ("ellipse", StationCurve::Two(&ellipse)),
        ("b-spline", StationCurve::Three(&spline)),
        (
            "banked over an intrinsic plan",
            StationCurve::Three(&banked),
        ),
        ("chain", StationCurve::Two(&chain)),
        ("an intrinsic profile", StationCurve::Three(&vertical)),
    ];
    for (name, curve) in cases {
        let law = FrameLaw::atomic(curve, None, TOLERANCE).unwrap();
        for reading in READINGS {
            assert!(
                law.error_at(10.0, SeamSide::Outgoing, reading).is_none(),
                "{name}"
            );
        }
    }
    // An offset beside a circle.
    let circle = Curve2::Circle(Circle2 {
        frame: frame(),
        radius: 30.0,
    });
    let base = StationPiece::whole(StationCurve::Two(&circle)).unwrap();
    let basis =
        CompositeBasis::new(offset_pieces(&[base], OffsetLaw::Planar { distance: 2.0 }).unwrap())
            .unwrap();
    let law = basis.frame_law(TOLERANCE).unwrap();
    assert!(law
        .error_at(10.0, SeamSide::Outgoing, FrameReading::Section)
        .is_none());
    // A distance off the basis.
    let law = FrameLaw::atomic(StationCurve::Two(&circle), None, TOLERANCE).unwrap();
    for d in [-1.0, 4.0 * PI * 30.0, Scalar::NAN] {
        assert!(law
            .error_at(d, SeamSide::Outgoing, FrameReading::Section)
            .is_none());
    }
}

/// The error reads the frame `section_on` reads: a piece whole is its
/// curve's section at the same distance and side.
#[test]
fn a_piece_reads_the_curves_own_section() {
    let a = alignment();
    let plan = Curve2::Intrinsic(a.curve.clone());
    let elevated = Curve3::Elevated(Elevated3::new(plan.clone(), profile()));
    let circle = Curve2::Circle(Circle2 {
        frame: frame(),
        radius: 30.0,
    });
    let polyline = polyline();
    let two: [&Curve2; 3] = [&plan, &circle, &polyline];
    for curve in two {
        let piece = StationPiece::whole(StationCurve::Two(curve)).unwrap();
        for u in [3.0, 40.0, 100.0] {
            for side in SIDES {
                let got = piece.section_on(u, side).unwrap();
                assert_eq!(got, station_section2_on(curve, u, side).unwrap());
            }
        }
    }
    let piece = StationPiece::whole(StationCurve::Three(&elevated)).unwrap();
    for u in [3.0, 120.0, 180.0] {
        for side in SIDES {
            let got = piece.section_on(u, side).unwrap();
            assert_eq!(got, station_section3_on(&elevated, u, side).unwrap());
        }
    }
}
