//! Certified positions of arc-length chains (#438).
//!
//! Every reference is independent of the table: compensated Simpson of the
//! closed-form heading, the cubic's binomial series, circles and lines in
//! closed form, and the equivalent piecewise `Intrinsic2` as a cross-check.

use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{
    BSplineCurve2, Chain2, ChainPiece2, CurvatureLaw, Curve2, Intrinsic2, KnotSpec, Line2,
    Sinusoid2,
};
use axiolid_evaluate::bound::{ChainPositionTable, IntrinsicPositionTable, SeamSide};
use axiolid_evaluate::chain_point;

// --- references --------------------------------------------------------------

fn rotate(v: Vec2, angle: Scalar) -> Vec2 {
    let (s, c) = angle.sin_cos();
    Vec2::new(c * v.x - s * v.y, s * v.x + c * v.y)
}

fn frame(origin: Point2, angle: Scalar) -> Frame2 {
    Frame2 {
        origin,
        x: rotate(Vec2::X, angle),
        y: rotate(Vec2::Y, angle),
    }
}

/// Neumaier-compensated sum.
#[derive(Default)]
struct Sum(Scalar, Scalar);

impl Sum {
    fn add(&mut self, x: Scalar) {
        let t = self.0 + x;
        self.1 += if self.0.abs() >= x.abs() {
            (self.0 - t) + x
        } else {
            (x - t) + self.0
        };
        self.0 = t;
    }
    fn value(&self) -> Scalar {
        self.0 + self.1
    }
}

/// `int_a^b (cos theta, sin theta)` by compensated composite Simpson with
/// panels of at most `step`.
fn simpson(theta: &dyn Fn(Scalar) -> Scalar, a: Scalar, b: Scalar, step: Scalar) -> Vec2 {
    if b <= a {
        return Vec2::ZERO;
    }
    let n = (((b - a) / step).ceil() as usize).max(1);
    let h = (b - a) / n as Scalar;
    let (mut x, mut y) = (Sum::default(), Sum::default());
    for i in 0..=n {
        let s = if i == n { b } else { a + h * i as Scalar };
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let t = theta(s);
        x.add(weight * t.cos());
        y.add(weight * t.sin());
    }
    Vec2::new(x.value() * h / 3.0, y.value() * h / 3.0)
}

fn world(frame: &Frame2, local: Vec2) -> Point2 {
    frame.origin + frame.x * local.x + frame.y * local.y
}

fn dist(a: Point2, b: Point2) -> Scalar {
    (a - b).length()
}

// --- the cubic transition (as tests/arc_length_chain.rs) ----------------------

fn cubic(a: Scalar, reach: Scalar) -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
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
        knot_spec: KnotSpec::PiecewiseBezier,
    })
}

/// Arc length of `y = a x^3` from 0 to `x` by its binomial series.
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
        let step = (cubic_length(a, x) - length) / (1.0 + 9.0 * a * a * x.powi(4)).sqrt();
        x -= step;
        if step.abs() < 1e-16 * length {
            break;
        }
    }
    x
}

const RADIUS: Scalar = 300.0;
const TRANSITION: Scalar = 60.0;
const STRAIGHT: Scalar = 50.0;
const ARC: Scalar = 40.0;

fn rate() -> Scalar {
    1.0 / (6.0 * RADIUS * TRANSITION)
}

fn layout_frame() -> Frame2 {
    frame(Point2::new(5.0, -3.0), 30f64.to_radians())
}

fn line() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Vec2::ZERO,
        direction: Vec2::X,
    })
}

/// Line -> cubic -> arc, as `tests/arc_length_chain.rs::layout`.
fn layout() -> Chain2 {
    Chain2::new(
        layout_frame(),
        vec![
            ChainPiece2::Parametric {
                curve: line(),
                start: 0.0,
                length: STRAIGHT,
            },
            ChainPiece2::Parametric {
                curve: cubic(rate(), TRANSITION),
                start: 0.0,
                length: TRANSITION,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::circular(1.0 / RADIUS),
                length: ARC,
            },
        ],
    )
}

/// The exact point, unit tangent angle and signed curvature of `layout()`
/// at `s`, in the local frame of the chain start.
fn layout_exact(s: Scalar) -> (Vec2, Scalar, Scalar) {
    let a = rate();
    // A join belongs to the piece that starts there.
    if s < STRAIGHT {
        return (Vec2::new(s, 0.0), 0.0, 0.0);
    }
    if s < STRAIGHT + TRANSITION {
        let x = cubic_abscissa(a, s - STRAIGHT);
        let slope = 3.0 * a * x * x;
        let k = 6.0 * a * x / (1.0 + slope * slope).powf(1.5);
        return (Vec2::new(STRAIGHT + x, a * x.powi(3)), slope.atan(), k);
    }
    let x = cubic_abscissa(a, TRANSITION);
    let heading = (3.0 * a * x * x).atan();
    let join = Vec2::new(STRAIGHT + x, a * x.powi(3));
    let phi = (s - STRAIGHT - TRANSITION) / RADIUS;
    let arc = Vec2::new(RADIUS * phi.sin(), RADIUS * (1.0 - phi.cos()));
    (join + rotate(arc, heading), heading + phi, 1.0 / RADIUS)
}

fn unit_angle(frame: &Frame2, angle: Scalar) -> Vec2 {
    frame.x * angle.cos() + frame.y * angle.sin()
}

fn left_normal(frame: &Frame2, angle: Scalar) -> Vec2 {
    frame.x * -angle.sin() + frame.y * angle.cos()
}

// --- tests -------------------------------------------------------------------

/// Straight 40, clothoid up, arc, clothoid down, straight 40, with the
/// closed-form heading of each piece.
const A: Scalar = 1.0 / (300.0 * 60.0);

fn alignment_pieces() -> Vec<ChainPiece2> {
    vec![
        ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::straight(),
            length: 40.0,
        },
        ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
            length: 60.0,
        },
        ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::circular(1.0 / 300.0),
            length: 80.0,
        },
        ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(1.0 / 300.0, 0.0, 60.0),
            length: 60.0,
        },
        ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::straight(),
            length: 40.0,
        },
    ]
}

fn alignment_frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    }
}

/// Heading and curvature of the alignment at `s`.
fn alignment_heading(s: Scalar) -> (Scalar, Scalar) {
    let turn_up = 0.5 * 60.0 / 300.0;
    let turn_arc = 80.0 / 300.0;
    if s <= 40.0 {
        (0.0, 0.0)
    } else if s <= 100.0 {
        (0.5 * A * (s - 40.0).powi(2), A * (s - 40.0))
    } else if s <= 180.0 {
        (turn_up + (s - 100.0) / 300.0, 1.0 / 300.0)
    } else if s <= 240.0 {
        let m = s - 180.0;
        (
            turn_up + turn_arc + m / 300.0 - 0.5 * A * m * m,
            1.0 / 300.0 - A * m,
        )
    } else {
        (turn_up + turn_arc + turn_up, 0.0)
    }
}

fn alignment_point(s: Scalar) -> Point2 {
    // Integrate piece by piece so Simpson never straddles a kink.
    let mut local = Vec2::ZERO;
    let mut from = 0.0;
    for to in [40.0, 100.0, 180.0, 240.0, 280.0] {
        if s <= from {
            break;
        }
        let upto = s.min(to);
        local += simpson(&|x| alignment_heading(x).0, from, upto, 1.0 / 2048.0);
        from = to;
    }
    world(&alignment_frame(), local)
}

#[test]
fn an_all_intrinsic_chain_reads_as_its_piecewise_curve() {
    let chain = Chain2::new(alignment_frame(), alignment_pieces());
    let table = ChainPositionTable::new(&chain, 1e-9, 280.0).unwrap();
    assert!(table.max_error() <= 2e-9, "{}", table.max_error());
    let piecewise = Intrinsic2::new(
        alignment_frame(),
        CurvatureLaw::piecewise(
            vec![40.0, 100.0, 180.0, 240.0],
            vec![
                CurvatureLaw::straight(),
                CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
                CurvatureLaw::circular(1.0 / 300.0),
                CurvatureLaw::clothoid(1.0 / 300.0, 0.0, 60.0),
                CurvatureLaw::straight(),
            ],
        ),
        280.0,
    );
    let reference = IntrinsicPositionTable::new2(&piecewise, 1e-9, 280.0).unwrap();
    for k in 0..40 {
        let s = 0.5 + 7.0 * k as Scalar;
        let (point, error) = table.point2(s).unwrap();
        assert!(error <= table.max_error());
        let exact = alignment_point(s);
        assert!(
            dist(point, exact) <= error + 1e-12,
            "s={s}: {point:?} against {exact:?}, error {error}"
        );
        let (other, other_error) = reference.point2(s).unwrap();
        assert!(dist(point, other) <= error + other_error, "s={s}");

        let (heading, kappa) = alignment_heading(s);
        let (tangent, terror) = table.tangent2_on(s, SeamSide::Outgoing).unwrap();
        let expected = unit_angle(&alignment_frame(), heading);
        assert!(
            (tangent - expected).length() <= terror + 1e-13,
            "tangent at {s}"
        );
        assert!(terror <= table.max_tangent_error());
        let (curvature, cerror) = table.curvature2_on(s, SeamSide::Outgoing).unwrap();
        let expected = left_normal(&alignment_frame(), heading) * kappa;
        assert!(
            (curvature - expected).length() <= cerror + 1e-13,
            "curvature at {s}"
        );
        assert!(cerror <= table.max_curvature_error());
    }
}

#[test]
fn a_cubic_transition_chain_matches_its_series() {
    let table = ChainPositionTable::new(&layout(), 1e-9, STRAIGHT + TRANSITION + ARC).unwrap();
    let f = layout_frame();
    for k in 0..=150 {
        let s = 150.0 * k as Scalar / 150.0;
        let (local, angle, kappa) = layout_exact(s);
        let (point, error) = table.point2(s).unwrap();
        assert!(error <= table.max_error());
        assert!(
            dist(point, world(&f, local)) <= error,
            "s={s}: off by {} against {error}",
            dist(point, world(&f, local))
        );
        let (tangent, terror) = table.tangent2_on(s, SeamSide::Outgoing).unwrap();
        assert!(
            (tangent - unit_angle(&f, angle)).length() <= terror,
            "tangent at {s}"
        );
        let (curvature, cerror) = table.curvature2_on(s, SeamSide::Outgoing).unwrap();
        assert!(
            (curvature - left_normal(&f, angle) * kappa).length() <= cerror,
            "curvature at {s}"
        );
    }
    // The enclosure is tight where it matters: nanometres, not microns.
    assert!(table.max_error() <= 1e-8, "{}", table.max_error());
}

/// A fast harmonic after a straight: `k = 0.01 sin 50 s` over 100 m, where
/// the evaluator's quadrature is off by a centimetre.
fn harmonic() -> Chain2 {
    Chain2::new(
        Frame2 {
            origin: Point2::new(20.0, 10.0),
            x: Vec2::new(0.8, 0.6),
            y: Vec2::new(-0.6, 0.8),
        },
        vec![
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::straight(),
                length: 20.0,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::Sinusoid {
                    mean: 0.0,
                    amplitude: 0.01,
                    angular_frequency: 50.0,
                    phase: 0.0,
                },
                length: 100.0,
            },
        ],
    )
}

#[test]
fn a_fast_harmonic_piece_is_read_to_its_error() {
    let chain = harmonic();
    let table = ChainPositionTable::new(&chain, 1e-6, 120.0).unwrap();
    let f = chain.start;
    for s in [116.4375, 116.5, 116.5625, 119.52] {
        let heading = |x: Scalar| 0.01 / 50.0 * (1.0 - (50.0 * x).cos());
        let local = Vec2::new(20.0, 0.0) + simpson(&heading, 0.0, s - 20.0, 1.0 / 8192.0);
        let exact = world(&f, local);
        let (point, error) = table.point2(s).unwrap();
        assert!(
            dist(point, exact) <= error,
            "s={s}: {} vs {error}",
            dist(point, exact)
        );
        // The table is the exact chain: the evaluator is the one that is
        // off, by more than the table's own error.
        if s == 116.5 {
            let read = chain_point(&chain, s).unwrap();
            assert!(dist(read, exact) > 1e-3 && dist(read, exact) > 100.0 * error);
        }
    }
}

#[test]
fn an_offset_join_reads_each_side() {
    let f = Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    };
    let chain = Chain2::new(
        f,
        vec![
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::straight(),
                length: 10.0,
            },
            ChainPiece2::Parametric {
                curve: Curve2::Line(Line2 {
                    origin: Vec2::new(0.0, 5e-10),
                    direction: Vec2::X,
                }),
                start: 0.0,
                length: 10.0,
            },
        ],
    );
    let table = ChainPositionTable::new(&chain, 1e-9, 20.0).unwrap();
    let (incoming, ei) = table.point2_on(10.0, SeamSide::Incoming).unwrap();
    let (outgoing, eo) = table.point2_on(10.0, SeamSide::Outgoing).unwrap();
    assert!(ei < 2.5e-10 && eo < 2.5e-10, "{ei} {eo}");
    assert!(dist(incoming, world(&f, Vec2::new(10.0, 0.0))) <= ei);
    assert!(dist(outgoing, world(&f, Vec2::new(10.0, 5e-10))) <= eo);
    assert!(dist(incoming, outgoing) > 4.9e-10);
    assert_eq!(table.point2(10.0).unwrap().0, outgoing);
    // Off a join the side makes no difference.
    assert_eq!(
        table.point2_on(5.0, SeamSide::Incoming).unwrap().0,
        table.point2_on(5.0, SeamSide::Outgoing).unwrap().0
    );
}

#[test]
fn errors_grow_along_the_chain() {
    // A cubic first, so its enclosure is the largest error before the lines,
    // then two lines: the second line must carry the cubic's error through
    // the first, not only the first line's own.
    let chain = Chain2::new(
        layout_frame(),
        vec![
            ChainPiece2::Parametric {
                curve: cubic(rate(), TRANSITION),
                start: 0.0,
                length: TRANSITION,
            },
            ChainPiece2::Parametric {
                curve: line(),
                start: 0.0,
                length: 30.0,
            },
            ChainPiece2::Parametric {
                curve: line(),
                start: 0.0,
                length: 30.0,
            },
        ],
    );
    let table = ChainPositionTable::new(&chain, 1e-6, 120.0).unwrap();
    let mut previous = None;
    for join in [TRANSITION, TRANSITION + 30.0] {
        let (_, end) = table.point2_on(join, SeamSide::Incoming).unwrap();
        let (_, begin) = table.point2_on(join, SeamSide::Outgoing).unwrap();
        let (_, te) = table.tangent2_on(join, SeamSide::Incoming).unwrap();
        let (_, tb) = table.tangent2_on(join, SeamSide::Outgoing).unwrap();
        assert!(
            end > 1e-9,
            "the cubic's enclosure is the error to carry: {end}"
        );
        assert!(te > 1e-12, "and its tangent's: {te}");
        assert!(begin >= end, "{begin} against {end} at {join}");
        assert!(tb >= te, "{tb} against {te} at {join}");
        // Each join carries at least what the one before it did.
        if let Some((pe, pt)) = previous {
            assert!(end >= pe && te >= pt, "{end} {te} against {pe} {pt}");
        }
        previous = Some((end, te));
    }
    let (_, far) = table.point2(120.0).unwrap();
    assert!(far >= previous.unwrap().0);
    // A line is straight whatever turns it: its curvature carries no error,
    // the cubic's does.
    let (_, ce) = table.curvature2_on(TRANSITION, SeamSide::Incoming).unwrap();
    let (_, cb) = table.curvature2_on(TRANSITION, SeamSide::Outgoing).unwrap();
    assert!(ce > 0.0 && cb == 0.0);
}

#[test]
fn what_no_table_covers_is_refused() {
    let straight = |length| ChainPiece2::Intrinsic {
        curvature: CurvatureLaw::straight(),
        length,
    };
    let f = layout_frame();
    // A family without a certified arc length.
    let sinusoid = Chain2::new(
        f,
        vec![ChainPiece2::Parametric {
            curve: Curve2::Sinusoid(Sinusoid2 {
                mean: 0.0,
                cosine: 0.0,
                sine: 0.0,
            }),
            start: 0.0,
            length: 3.0,
        }],
    );
    assert!(ChainPositionTable::new(&sinusoid, 1e-9, 3.0).is_none());
    // A piece longer than its curve (the cubic is 60 m long).
    let long = Chain2::new(
        f,
        vec![ChainPiece2::Parametric {
            curve: cubic(rate(), TRANSITION),
            start: 0.0,
            length: 61.0,
        }],
    );
    assert!(ChainPositionTable::new(&long, 1e-9, 61.0).is_none());
    // Not finite, not positive, past the end, empty.
    let chain = Chain2::new(f, vec![straight(10.0)]);
    for (error, end) in [
        (Scalar::NAN, 5.0),
        (0.0, 5.0),
        (1e-9, Scalar::INFINITY),
        (1e-9, 0.0),
        (1e-9, 10.5),
    ] {
        assert!(
            ChainPositionTable::new(&chain, error, end).is_none(),
            "{error} {end}"
        );
    }
    assert!(ChainPositionTable::new(&Chain2::new(f, Vec::new()), 1e-9, 1.0).is_none());
    // Reads outside the table.
    let table = ChainPositionTable::new(&chain, 1e-9, 6.0).unwrap();
    assert!(table.point2(6.5).is_none() && table.point2(-1.0).is_none());
    assert!(table.point2(Scalar::NAN).is_none());
    assert!(table.point2(6.0).is_some());
}
