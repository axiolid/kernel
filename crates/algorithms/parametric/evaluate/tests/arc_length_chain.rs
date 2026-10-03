//! Arc length along any curve, its inverse, and arc-length chains (#239).
//!
//! Every expectation is computed independently of the quadrature under
//! test: the cubic parabola's arc length by its binomial series, circles
//! and lines in closed form, the ellipse perimeter from its published
//! value, and chain ends by placing closed-form piece ends by hand.

use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Interval, Point2, Scalar, Vec2};
use axiolid_curve::{
    BSplineCurve2, BankConvention, Banked3, CantLaw, Chain2, ChainPiece2, Circle2, CurvatureLaw,
    Curve2, Curve3, Elevated3, ElevationLaw, Ellipse2, KnotSpec, Line2, Polyline2,
};
use axiolid_curve_evaluate_contract::{CurveEvaluator, CurveMeasure};
use axiolid_evaluate::bound::{certifies_flattening2, chord_bound2, continuity_breaks2};
use axiolid_evaluate::chain::chain_chord_bound;
use axiolid_evaluate::{
    arc_length2, chain_point, chain_tangent, derivative2, elevated_point, elevated_tangent,
    evaluate2, flatten2, parameter_at_arc_length2, ReferenceCurveEvaluator,
};

// --- references ------------------------------------------------------------

/// The cubic parabola `y = a x^3` over `x in [0, reach]` as an exact cubic
/// Bezier: control points at thirds along x, the last lifted to `a reach^3`.
/// Its parameter is `x / reach`.
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

/// Arc length of `y = a x^3` from 0 to `x` by its binomial series:
/// `sum_n C(1/2, n) (9 a^2)^n x^(4n+1) / (4n+1)`, valid for `9 a^2 x^4 < 1`.
fn cubic_length(a: Scalar, x: Scalar) -> Scalar {
    let z = 9.0 * a * a * x.powi(4);
    assert!(z < 0.5, "series reference needs 9 a^2 x^4 < 1/2, got {z}");
    let mut binomial = 1.0;
    let mut power = 1.0;
    let mut sum = 0.0;
    for n in 0..200 {
        let term = binomial * power / (4.0 * n as Scalar + 1.0);
        sum += term;
        if term.abs() < 1e-30 {
            break;
        }
        // C(1/2, n + 1) = C(1/2, n) (1/2 - n) / (n + 1).
        binomial *= (0.5 - n as Scalar) / (n as Scalar + 1.0);
        power *= z;
    }
    x * sum
}

/// Abscissa where `y = a x^3` has arc length `length`, by Newton on the
/// series: `s'(x) = sqrt(1 + 9 a^2 x^4)`.
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

fn close(a: Point2, b: Point2, tol: Scalar) -> bool {
    (a - b).length() <= tol
}

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

fn message(error: &GeomError) -> String {
    error.to_string()
}

// --- arc length and its inverse ---------------------------------------------

#[test]
fn a_cubic_parabola_trimmed_by_arc_length_ends_where_its_series_says() {
    // IFC-style CUBIC: y = x^3 / (6 R L), read to arc length L. Two cases: a
    // gentle transition (R = 300, L = 60) and a sharp one (R = L = 10) where
    // the abscissa falls well short of L.
    for (radius, length) in [(300.0, 60.0), (10.0, 10.0)] {
        let a = 1.0 / (6.0 * radius * length);
        let curve = cubic(a, length);
        let t = parameter_at_arc_length2(&curve, 0.0, length).unwrap();
        let x = cubic_abscissa(a, length);
        assert!(
            (t * length - x).abs() <= 1e-13 * length,
            "R={radius} L={length}: abscissa {} against {x}",
            t * length
        );
        assert!(x < length, "the curve is longer than its chord");

        let end = evaluate2(&curve, t).unwrap();
        assert!(
            close(end, Point2::new(x, a * x.powi(3)), 1e-13 * length),
            "R={radius}: end {end:?}"
        );
        let tangent = derivative2(&curve, t).unwrap().normalize();
        let slope = 3.0 * a * x * x;
        let expected = Vec2::new(1.0, slope) / (1.0 + slope * slope).sqrt();
        assert!(
            (tangent - expected).length() <= 1e-12,
            "R={radius}: tangent {tangent:?}"
        );

        // The quadrature itself, against the series.
        let measured = arc_length2(&curve, 0.0, 1.0).unwrap();
        let reference = cubic_length(a, length);
        assert!(
            (measured - reference).abs() <= 1e-12 * reference,
            "R={radius}: arc length {measured} against {reference}"
        );
    }
}

#[test]
fn closed_forms_hold_for_lines_circles_and_the_ellipse_perimeter() {
    // A line advances |direction| per unit parameter.
    let line = Curve2::Line(Line2 {
        origin: Vec2::new(1.0, 1.0),
        direction: Vec2::new(3.0, 4.0),
    });
    assert_eq!(parameter_at_arc_length2(&line, 0.0, 10.0).unwrap(), 2.0);
    assert_eq!(parameter_at_arc_length2(&line, 1.0, -5.0).unwrap(), 0.0);
    assert_eq!(arc_length2(&line, 2.0, 0.0).unwrap(), -10.0);

    // A circle of radius 2: half a turn is 2 pi of arc, and either way.
    let circle = Curve2::Circle(Circle2 {
        frame: frame(Point2::new(4.0, -1.0), 0.3),
        radius: 2.0,
    });
    let pi = core::f64::consts::PI;
    let t = parameter_at_arc_length2(&circle, 0.0, pi).unwrap();
    assert!((t - pi / 2.0).abs() <= 1e-12, "{t}");
    let back = parameter_at_arc_length2(&circle, 0.0, -pi).unwrap();
    assert!((back + pi / 2.0).abs() <= 1e-12, "{back}");
    // Past a full turn: a periodic curve runs on.
    let far = parameter_at_arc_length2(&circle, 0.0, 6.0 * pi).unwrap();
    assert!((far - 3.0 * pi).abs() <= 1e-11, "{far}");

    // Perimeter of the ellipse with semi-axes 2 and 1: 9.688448220547675...
    let ellipse = Curve2::Ellipse(Ellipse2 {
        frame: frame(Point2::ZERO, 0.0),
        semi_axis_x: 2.0,
        semi_axis_y: 1.0,
    });
    let perimeter = arc_length2(&ellipse, 0.0, core::f64::consts::TAU).unwrap();
    assert!(
        (perimeter - 9.688_448_220_547_675).abs() <= 1e-11,
        "{perimeter}"
    );
    let backwards = arc_length2(&ellipse, core::f64::consts::TAU, 0.0).unwrap();
    assert_eq!(backwards, -perimeter);

    // A polyline's speed jumps at its vertices: 5 + 6 + 4 = 15 of length.
    let polyline = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 4.0),
            Point2::new(3.0, 10.0),
            Point2::new(7.0, 10.0),
        ],
        closed: false,
    });
    assert!((arc_length2(&polyline, 0.0, 3.0).unwrap() - 15.0).abs() <= 1e-12);
    let t = parameter_at_arc_length2(&polyline, 0.0, 8.0).unwrap();
    assert!((t - 1.5).abs() <= 1e-12, "{t}");
    let t = parameter_at_arc_length2(&polyline, 3.0, -2.0).unwrap();
    assert!((t - 2.5).abs() <= 1e-12, "{t}");
}

#[test]
fn arc_length_refuses_what_it_cannot_honour() {
    let curve = cubic(1e-3, 10.0);
    let available = arc_length2(&curve, 0.0, 1.0).unwrap();
    let error = parameter_at_arc_length2(&curve, 0.0, available + 1.0).unwrap_err();
    assert!(message(&error).contains("exceeds"), "{error}");
    let error = parameter_at_arc_length2(&curve, 2.0, 1.0).unwrap_err();
    assert!(
        message(&error).contains("outside the curve's domain"),
        "{error}"
    );
    assert!(parameter_at_arc_length2(&curve, 0.0, Scalar::NAN).is_err());
    assert!(parameter_at_arc_length2(&curve, Scalar::INFINITY, 1.0).is_err());
    assert!(arc_length2(&curve, 0.0, 3.0).is_err());
    // Exactly the available length is the curve end, not a refusal.
    let end = parameter_at_arc_length2(&curve, 0.0, available).unwrap();
    assert!((end - 1.0).abs() <= 1e-9, "{end}");
}

// --- chains -----------------------------------------------------------------

const RADIUS: Scalar = 300.0;
const TRANSITION: Scalar = 60.0;
const STRAIGHT: Scalar = 50.0;
const ARC: Scalar = 40.0;

fn transition_rate() -> Scalar {
    1.0 / (6.0 * RADIUS * TRANSITION)
}

/// Line -> CUBIC -> arc, placed at (5, -3) heading 30 degrees.
fn layout() -> Chain2 {
    Chain2::new(
        frame(Point2::new(5.0, -3.0), 30f64.to_radians()),
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
                curve: cubic(transition_rate(), TRANSITION),
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

/// Chain point at `s`, computed by hand from closed-form piece ends.
fn layout_reference(s: Scalar) -> Point2 {
    let start = frame(Point2::new(5.0, -3.0), 30f64.to_radians());
    let world = |local: Vec2| start.origin + start.x * local.x + start.y * local.y;
    if s <= STRAIGHT {
        return world(Vec2::new(s, 0.0));
    }
    let a = transition_rate();
    if s <= STRAIGHT + TRANSITION {
        let x = cubic_abscissa(a, s - STRAIGHT);
        return world(Vec2::new(STRAIGHT + x, a * x.powi(3)));
    }
    let x = cubic_abscissa(a, TRANSITION);
    let heading = (3.0 * a * x * x).atan();
    let join = Vec2::new(STRAIGHT + x, a * x.powi(3));
    let phi = (s - STRAIGHT - TRANSITION) / RADIUS;
    let arc = Vec2::new(RADIUS * phi.sin(), RADIUS * (1.0 - phi.cos()));
    world(join + rotate(arc, heading))
}

#[test]
fn a_line_cubic_arc_chain_matches_its_hand_placed_pieces() {
    let chain = Curve2::Chain(layout());
    let total = STRAIGHT + TRANSITION + ARC;
    let domain = axiolid_evaluate::curve::domain2(&chain);
    assert_eq!(domain, Interval::new(0.0, total));
    for k in 0..=150 {
        let s = total * k as Scalar / 150.0;
        let point = evaluate2(&chain, s).unwrap();
        let expected = layout_reference(s);
        assert!(
            close(point, expected, 1e-11),
            "s={s}: {point:?} against {expected:?}"
        );
        let tangent = derivative2(&chain, s).unwrap();
        assert!(
            (tangent.length() - 1.0).abs() <= 1e-12,
            "unit tangent at {s}"
        );
    }
}

#[test]
fn a_chain_is_continuous_and_tangent_continuous_at_its_joins() {
    let chain = layout();
    for join in [STRAIGHT, STRAIGHT + TRANSITION] {
        let before = chain_point(&chain, join - 1e-7).unwrap();
        let at = chain_point(&chain, join).unwrap();
        assert!(close(before, at, 2e-7), "position jumps at {join}");
        let t_before = chain_tangent(&chain, join - 1e-7).unwrap();
        let t_at = chain_tangent(&chain, join).unwrap();
        // G1: the tangent turns by at most the curvature times the step.
        assert!(
            (t_before - t_at).length() <= 1e-9,
            "tangent jumps at {join}"
        );
    }
    // The arc starts in the direction the cubic ends in.
    let a = transition_rate();
    let x = cubic_abscissa(a, TRANSITION);
    let heading = 30f64.to_radians() + (3.0 * a * x * x).atan();
    let t = chain_tangent(&chain, STRAIGHT + TRANSITION).unwrap();
    assert!((t - rotate(Vec2::X, heading)).length() <= 1e-12, "{t:?}");
}

#[test]
fn a_chain_flattens_within_a_certified_chord_bound() {
    let chain = Curve2::Chain(layout());
    let total = STRAIGHT + TRANSITION + ARC;
    assert!(certifies_flattening2(&chain));
    assert_eq!(
        continuity_breaks2(&chain, 1),
        vec![STRAIGHT, STRAIGHT + TRANSITION]
    );
    // Across a join there is no single bound; inside a piece there is.
    assert!(chord_bound2(&chain, 40.0, 60.0).is_none());
    // The chain's own bound refuses the join too, not only through the
    // continuity breaks `chord_bound2` checks first.
    let Curve2::Chain(plain) = &chain else {
        unreachable!()
    };
    assert!(chain_chord_bound(plain, 40.0, 60.0).is_none());
    assert!(chain_chord_bound(plain, 0.0, STRAIGHT).is_some());
    let straight = chord_bound2(&chain, 0.0, STRAIGHT).unwrap();
    assert!(straight <= 1e-12, "a straight piece bows by {straight}");
    let arc = chord_bound2(&chain, STRAIGHT + TRANSITION, total).unwrap();
    let sagitta = RADIUS * (1.0 - (ARC / RADIUS / 2.0).cos());
    assert!(
        arc >= sagitta,
        "arc bound {arc} under its sagitta {sagitta}"
    );

    let tolerance = 1e-3;
    let points = flatten2(&chain, Interval::new(0.0, total), tolerance, 20).unwrap();
    // Every point of the chain lies within the tolerance of the polyline.
    for k in 0..=3000 {
        let s = total * k as Scalar / 3000.0;
        let p = layout_reference(s);
        let distance = points
            .windows(2)
            .map(|w| segment_distance(p, w[0], w[1]))
            .fold(Scalar::INFINITY, Scalar::min);
        assert!(distance <= tolerance * (1.0 + 1e-6), "s={s}: {distance}");
    }
}

#[test]
fn a_clothoid_piece_is_bounded_by_its_curvature_law() {
    // Line -> clothoid into R = 200 over 80 -> nothing else: every sub-span's
    // bound must cover the chain's measured distance from its chord.
    let (radius, length) = (200.0, 80.0);
    let chain = Curve2::Chain(Chain2::new(
        frame(Point2::new(1.0, 2.0), 0.4),
        vec![
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::circular(0.0),
                length: 10.0,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::clothoid(0.0, 1.0 / radius, length),
                length,
            },
        ],
    ));
    for (a, b) in [(10.0, 90.0), (10.0, 30.0), (50.0, 90.0), (70.0, 75.0)] {
        let bound = chord_bound2(&chain, a, b).unwrap();
        let (pa, pb) = (evaluate2(&chain, a).unwrap(), evaluate2(&chain, b).unwrap());
        let worst = (0..=400)
            .map(|k| {
                let p = evaluate2(&chain, a + (b - a) * k as Scalar / 400.0).unwrap();
                segment_distance(p, pa, pb)
            })
            .fold(0.0, Scalar::max);
        assert!(
            worst > 0.0 && bound >= worst,
            "[{a}, {b}]: bound {bound} under {worst}"
        );
        // And not wildly loose: within the h^2/8 sup|k| it claims.
        assert!(
            bound <= (b - a).powi(2) / 8.0 / radius * 1.01,
            "[{a}, {b}]: {bound}"
        );
    }
}

#[test]
fn a_piece_with_a_corner_is_not_certified() {
    // A polyline piece turns a corner: the chain is then flattened on its
    // sagitta alone, and says so.
    let chain = Chain2::new(
        frame(Point2::ZERO, 0.0),
        vec![ChainPiece2::Parametric {
            curve: Curve2::Polyline(Polyline2 {
                points: vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(3.0, 0.0),
                    Point2::new(3.0, 3.0),
                ],
                closed: false,
            }),
            start: 0.0,
            length: 5.0,
        }],
    );
    assert!(!certifies_flattening2(&Curve2::Chain(chain.clone())));
    let corner = chain_point(&chain, 3.0).unwrap();
    assert!(close(corner, Point2::new(3.0, 0.0), 1e-12), "{corner:?}");
    let up = chain_point(&chain, 5.0).unwrap();
    assert!(close(up, Point2::new(3.0, 2.0), 1e-12), "{up:?}");
}

fn segment_distance(p: Point2, a: Point2, b: Point2) -> Scalar {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

#[test]
fn an_elevated_chain_reads_heights_by_plan_distance() {
    let plan = layout();
    let curve = Elevated3::new(
        Curve2::Chain(plan.clone()),
        ElevationLaw::parabolic(12.0, 0.02, -0.01, 150.0),
    );
    let law = ElevationLaw::parabolic(12.0, 0.02, -0.01, 150.0);
    let evaluator = ReferenceCurveEvaluator::new();
    let elevated = Curve3::Elevated(curve.clone());
    for d in [
        0.0,
        25.0,
        STRAIGHT,
        80.0,
        STRAIGHT + TRANSITION,
        130.0,
        150.0,
    ] {
        let p = elevated_point(&curve, d).unwrap();
        let plan_point = layout_reference(d);
        assert!((p.x - plan_point.x).abs() <= 1e-9 && (p.y - plan_point.y).abs() <= 1e-9);
        assert_eq!(
            p.z,
            law.height_at(d).unwrap(),
            "height at plan distance {d}"
        );
        let q = evaluator
            .point_at(&elevated, CurveMeasure::Distance(d))
            .unwrap();
        assert_eq!(p, q);
        // The tangent climbs at the law's grade over the plan tangent.
        let t = elevated_tangent(&curve, d).unwrap();
        let grade = law.grade_at(d).unwrap();
        let horizontal = t.x.hypot(t.y);
        assert!((t.z / horizontal - grade).abs() <= 1e-12, "grade at {d}");
    }
}

#[test]
fn a_banked_curve_over_an_elevated_chain_reads_the_chain() {
    // #240 places its sections on an `Elevated3`; with no cant and no pivot
    // the banked centreline is the elevated chain itself.
    let total = STRAIGHT + TRANSITION + ARC;
    let base = Elevated3::new(
        Curve2::Chain(layout()),
        ElevationLaw::constant_grade(3.0, 0.01),
    );
    let banked = Curve3::Banked(Banked3::new(
        base.clone(),
        CantLaw::zero(total),
        CantLaw::zero(total),
        1.5,
        BankConvention::TangentRotation,
    ));
    for d in [0.0, 70.0, total] {
        let p = axiolid_evaluate::evaluate3(&banked, d).unwrap();
        let q = elevated_point(&base, d).unwrap();
        assert!((p - q).length() <= 1e-12, "at {d}: {p:?} against {q:?}");
    }
}

// --- refusals ---------------------------------------------------------------

fn refused(chain: Chain2, s: Scalar) -> String {
    message(&chain_point(&chain, s).unwrap_err())
}

#[test]
fn a_chain_refuses_malformed_pieces_by_name() {
    let mut chain = layout();
    chain.pieces[1] = ChainPiece2::Parametric {
        curve: cubic(transition_rate(), TRANSITION),
        start: 0.0,
        length: Scalar::NAN,
    };
    assert!(refused(chain, 10.0).contains("piece 1"));

    let mut chain = layout();
    chain.pieces[2] = ChainPiece2::Intrinsic {
        curvature: CurvatureLaw::circular(0.01),
        length: -1.0,
    };
    assert!(refused(chain, 10.0).contains("piece 2"));

    let mut chain = layout();
    chain.pieces.clear();
    assert!(refused(chain, 0.0).contains("no pieces"));

    let mut chain = layout();
    chain.start.x = Vec2::new(Scalar::INFINITY, 0.0);
    assert!(refused(chain, 0.0).contains("start frame"));
}

#[test]
fn a_piece_longer_than_its_curve_is_refused_with_the_length_available() {
    let mut chain = layout();
    chain.pieces[1] = ChainPiece2::Parametric {
        curve: cubic(transition_rate(), TRANSITION),
        start: 0.0,
        length: 2.0 * TRANSITION,
    };
    let text = refused(chain.clone(), STRAIGHT + 1.0);
    assert!(
        text.contains("piece 1") && text.contains("exceeds"),
        "{text}"
    );
    // Past the piece too: the arc is placed at its end.
    assert!(chain_point(&chain, STRAIGHT + 2.0 * TRANSITION + 1.0).is_err());
}

#[test]
fn a_piece_off_its_local_frame_is_refused_rather_than_moved() {
    let mut chain = layout();
    chain.pieces[0] = ChainPiece2::Parametric {
        curve: Curve2::Line(Line2 {
            origin: Vec2::new(0.0, 1.0),
            direction: Vec2::X,
        }),
        start: 0.0,
        length: STRAIGHT,
    };
    assert!(refused(chain, 1.0).contains("local origin"));

    let mut chain = layout();
    chain.pieces[0] = ChainPiece2::Parametric {
        curve: Curve2::Line(Line2 {
            origin: Vec2::ZERO,
            direction: Vec2::new(1.0, 0.1),
        }),
        start: 0.0,
        length: STRAIGHT,
    };
    assert!(refused(chain, 1.0).contains("tangent +x"));
}

#[test]
fn a_distance_outside_the_chain_is_refused() {
    let total = STRAIGHT + TRANSITION + ARC;
    assert!(refused(layout(), -0.5).contains("outside"));
    assert!(refused(layout(), total + 0.5).contains("outside"));
    assert!(refused(layout(), Scalar::NAN).contains("finite"));
    // A rounding step past the end is the end.
    assert!(chain_point(&layout(), total * (1.0 + 1e-15)).is_ok());
}

#[test]
fn an_elevated_b_spline_plan_is_still_refused() {
    // A B-spline read by its own parameter is not a distance (ADR 0060);
    // only a chain piece reads it by arc length.
    let curve = Elevated3::new(cubic(1e-3, 10.0), ElevationLaw::level(0.0));
    assert!(elevated_point(&curve, 1.0).is_err());
}
