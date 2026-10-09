//! Elevated and banked curves through the generic 3D curve functions
//! (#252): points, derivatives and second derivatives in plan distance
//! against finite differences, and the certified derivative and chord
//! bounds against dense sampling, over straight, arc, clothoid and chain
//! plans with line, circular-arc, parabolic and intrinsic profiles. A
//! banked curve rotating about a held rail through a Viennese bend (#279)
//! derives its point path from the cant law and is held to the same
//! checks.

use axiolid_core::{Frame2, Interval, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve2, BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, Circle2,
    CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2, KnotSpec, Line2, RailSide,
};
use axiolid_evaluate::bound::{
    certifies_flattening3, chord_bound3, continuity_breaks3, curve_derivative_bounds3,
};
use axiolid_evaluate::curve::{derivative3, domain3, evaluate3, flatten3, second_derivative3};
use axiolid_evaluate::{elevated_point, grade_corners3};

fn world() -> Frame2 {
    Frame2 {
        origin: Point2::new(10.0, -5.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    }
}

fn straight() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::new(10.0, -5.0),
        direction: Vec2::new(3.0, 4.0),
    })
}

fn arc() -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: world(),
        radius: 300.0,
    })
}

fn clothoid() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        world(),
        CurvatureLaw::clothoid(0.0, 1.0 / 250.0, 120.0),
        120.0,
    ))
}

/// Straight 40 m, clothoid 60 m, arc 40 m: an alignment's plan.
fn chain() -> Curve2 {
    Curve2::Chain(Chain2::new(
        world(),
        vec![
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::Constant { curvature: 0.0 },
                length: 40.0,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::clothoid(0.0, 1.0 / 250.0, 60.0),
                length: 60.0,
            },
            ChainPiece2::Parametric {
                curve: Curve2::Circle(Circle2 {
                    frame: Frame2 {
                        origin: Point2::new(0.0, 250.0),
                        x: Vec2::new(0.0, -1.0),
                        y: Vec2::new(1.0, 0.0),
                    },
                    radius: 250.0,
                }),
                start: 0.0,
                length: 40.0,
            },
        ],
    ))
}

fn line_profile() -> ElevationLaw {
    ElevationLaw::constant_grade(100.0, 0.02)
}

/// The consumer's case: a vertical circular sag of R = 1000.
fn sag_profile() -> ElevationLaw {
    ElevationLaw::circular_arc(100.0, -0.04, 1000.0)
}

fn parabolic_profile() -> ElevationLaw {
    ElevationLaw::parabolic(100.0, 0.03, -0.02, 120.0)
}

/// 2% for 30 m, a crest arc tangent to it, then the parabola tangent to
/// the arc's exit: G1 at both seams.
fn gradient_profile() -> ElevationLaw {
    let radius = -2000.0;
    let arc = ElevationLaw::circular_arc(100.6, 0.02, radius);
    let exit = arc.grade_at(50.0).unwrap();
    let height = arc.height_at(50.0).unwrap();
    ElevationLaw::Piecewise {
        breaks: vec![30.0, 80.0],
        laws: vec![
            ElevationLaw::constant_grade(100.0, 0.02),
            arc,
            ElevationLaw::parabolic(height, exit, exit + 0.01, 100.0),
        ],
    }
}

fn intrinsic_profile() -> ElevationLaw {
    ElevationLaw::intrinsic(
        100.0,
        -0.01,
        CurvatureLaw::clothoid(0.0, 1.0 / 1500.0, 120.0),
    )
}

fn plans() -> Vec<(&'static str, Curve2, Scalar)> {
    vec![
        ("straight", straight(), 120.0),
        ("arc", arc(), 120.0),
        ("clothoid", clothoid(), 120.0),
        ("chain", chain(), 140.0),
    ]
}

fn profiles() -> Vec<(&'static str, ElevationLaw)> {
    vec![
        ("line", line_profile()),
        ("sag", sag_profile()),
        ("parabola", parabolic_profile()),
        ("gradient", gradient_profile()),
        ("intrinsic", intrinsic_profile()),
    ]
}

fn curves() -> Vec<(String, Curve3, Scalar)> {
    let mut out = Vec::new();
    for (plan_name, plan, length) in plans() {
        for (profile_name, profile) in profiles() {
            out.push((
                format!("{plan_name}/{profile_name}"),
                Curve3::Elevated(Elevated3::new(plan.clone(), profile)),
                length,
            ));
        }
    }
    out
}

fn banked() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(clothoid(), parabolic_profile()),
        CantLaw::new(vec![
            CantPiece::sine(60.0, 0.0, 0.12),
            CantPiece::constant(60.0, 0.12),
        ]),
        CantLaw::new(vec![
            CantPiece::sine(60.0, 0.0, 0.06),
            CantPiece::cosine(60.0, 0.06, 0.02),
        ]),
        1.5,
        BankConvention::TangentRotation,
    ))
}

/// Rotation about the held right rail (#279): level 20 m, a Viennese bend
/// to 150 mm over 60 m, then held; the point path follows `D / 2`.
fn held_rail(plan: Curve2, profile: ElevationLaw, convention: BankConvention) -> Curve3 {
    let top = (0.15_f64 / 1.5).asin();
    Curve3::Banked(Banked3::new(
        Elevated3::new(plan, profile),
        CantLaw::new(vec![
            CantPiece::constant(20.0, 0.0),
            CantPiece::viennese_bend(60.0, 0.0, top),
            CantPiece::constant(40.0, 0.15),
        ]),
        CantLaw::new(vec![CantPiece::about_rail(120.0, RailSide::Right, 0.02)]),
        1.5,
        convention,
    ))
}

/// The banked curves every check runs over: rotation about a moving
/// pivot law, and about a held rail on a curved and on a straight level
/// alignment (where the pivot alone bends the point path).
fn banked_curves() -> Vec<(String, Curve3, Scalar)> {
    vec![
        ("banked".into(), banked(), 120.0),
        (
            "held rail".into(),
            held_rail(
                clothoid(),
                parabolic_profile(),
                BankConvention::VerticalRise,
            ),
            120.0,
        ),
        (
            "held rail only".into(),
            held_rail(
                straight(),
                ElevationLaw::level(0.0),
                BankConvention::TangentRotation,
            ),
            120.0,
        ),
    ]
}

/// `[0, length]` cut at the curve's `C^2` breaks.
fn smooth_spans(curve: &Curve3, length: Scalar) -> Vec<(Scalar, Scalar)> {
    let mut cuts = vec![0.0];
    cuts.extend(
        continuity_breaks3(curve, 2)
            .into_iter()
            .filter(|&b| b > 0.0 && b < length),
    );
    cuts.push(length);
    cuts.windows(2).map(|w| (w[0], w[1])).collect()
}

#[test]
fn derivatives_match_finite_differences() {
    let mut all = curves();
    all.extend(banked_curves());
    for (name, curve, length) in &all {
        for (lo, hi) in smooth_spans(curve, *length) {
            for i in 1..10 {
                let d = lo + (hi - lo) * i as Scalar / 10.0;
                let h = 1e-3;
                let p = |t| evaluate3(curve, t).unwrap();
                let fd1 = (p(d + h) - p(d - h)) / (2.0 * h);
                let d1 = derivative3(curve, d).unwrap();
                assert!(
                    (fd1 - d1).length() <= 1e-7,
                    "{name} d = {d}: c' {d1:?} vs {fd1:?}"
                );
                let q = |t| derivative3(curve, t).unwrap();
                let fd2 = (q(d + h) - q(d - h)) / (2.0 * h);
                let d2 = second_derivative3(curve, d).unwrap();
                assert!(
                    (fd2 - d2).length() <= 1e-7,
                    "{name} d = {d}: c'' {d2:?} vs {fd2:?}"
                );
            }
        }
    }
}

#[test]
fn the_point_is_the_composition_and_the_domain_is_the_plans() {
    let elevated = Elevated3::new(clothoid(), sag_profile());
    let curve = Curve3::Elevated(elevated.clone());
    for d in [0.0, 17.5, 60.0, 120.0] {
        assert_eq!(
            evaluate3(&curve, d).unwrap(),
            elevated_point(&elevated, d).unwrap()
        );
    }
    assert_eq!(domain3(&curve), Interval::new(0.0, 120.0));
    let on_line = Curve3::Elevated(Elevated3::new(straight(), line_profile()));
    assert_eq!(domain3(&on_line).start, 0.0);
    assert_eq!(domain3(&on_line).end, Scalar::INFINITY);
    let on_arc = Curve3::Elevated(Elevated3::new(arc(), line_profile()));
    assert!((domain3(&on_arc).end - core::f64::consts::TAU * 300.0).abs() < 1e-9);
    let on_chain = Curve3::Elevated(Elevated3::new(chain(), line_profile()));
    assert!((domain3(&on_chain).end - 140.0).abs() < 1e-12);
}

#[test]
fn derivative_bounds_hold_over_every_smooth_span() {
    let mut all = curves();
    all.extend(banked_curves());
    for (name, curve, length) in &all {
        for (lo, hi) in smooth_spans(curve, *length) {
            // Several sub-spans: the bound over each must cover it.
            for k in 0..4 {
                let (a, b) = (
                    lo + (hi - lo) * k as Scalar / 4.0,
                    lo + (hi - lo) * (k + 1) as Scalar / 4.0,
                );
                let bounds = curve_derivative_bounds3(curve, a, b)
                    .unwrap_or_else(|| panic!("{name}: no bound over [{a}, {b}]"));
                for i in 0..=20 {
                    // Inside the span: a seam at its end belongs to the next piece.
                    let d = a + (b - a) * (i as Scalar + 0.5) / 21.0;
                    let d1 = derivative3(curve, d).unwrap().length();
                    let d2 = second_derivative3(curve, d).unwrap();
                    assert!(d1 <= bounds.first, "{name}: |c'| {d1} > {}", bounds.first);
                    assert!(
                        d2.length() <= bounds.second,
                        "{name}: |c''| {} > {}",
                        d2.length(),
                        bounds.second
                    );
                    let h = 1e-3;
                    let (u, v) = ((d - h).max(a), (d + h).min(b));
                    let d3 = ((second_derivative3(curve, v).unwrap()
                        - second_derivative3(curve, u).unwrap())
                        / (v - u))
                        .length();
                    assert!(
                        d3 <= bounds.third * (1.0 + 1e-6) + 1e-9,
                        "{name}: |c'''| {d3} > {} at {d}",
                        bounds.third
                    );
                }
            }
        }
    }
}

/// Distance from `p` to the segment `a`-`b`.
fn to_segment(p: Point3, a: Point3, b: Point3) -> Scalar {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

#[test]
fn the_composed_chord_bound_holds_and_is_not_loose() {
    let mut all = curves();
    all.extend(banked_curves());
    // Level and straight: the pivot alone bends the point path.
    all.push((
        "pivot only".into(),
        Curve3::Banked(Banked3::new(
            Elevated3::new(straight(), ElevationLaw::level(0.0)),
            CantLaw::new(vec![CantPiece::constant(40.0, 0.1)]),
            CantLaw::new(vec![
                CantPiece::sine(20.0, 0.0, 0.5),
                CantPiece::cosine(20.0, 0.5, 0.2),
            ]),
            1.5,
            BankConvention::TangentRotation,
        )),
        40.0,
    ));
    for (name, curve, length) in &all {
        let mut cuts = vec![0.0];
        cuts.extend(
            continuity_breaks3(curve, 1)
                .into_iter()
                .filter(|&b| b > 0.0 && b < *length),
        );
        cuts.push(*length);
        for w in cuts.windows(2) {
            for steps in [1, 3, 12] {
                for k in 0..steps {
                    let a = w[0] + (w[1] - w[0]) * k as Scalar / steps as Scalar;
                    let b = w[0] + (w[1] - w[0]) * (k + 1) as Scalar / steps as Scalar;
                    let bound = chord_bound3(curve, a, b)
                        .unwrap_or_else(|| panic!("{name}: no chord bound over [{a}, {b}]"));
                    let (pa, pb) = (evaluate3(curve, a).unwrap(), evaluate3(curve, b).unwrap());
                    let mut worst: Scalar = 0.0;
                    for i in 0..=200 {
                        let d = a + (b - a) * i as Scalar / 200.0;
                        worst = worst.max(to_segment(evaluate3(curve, d).unwrap(), pa, pb));
                    }
                    assert!(
                        worst <= bound + 1e-12,
                        "{name} [{a}, {b}]: deviation {worst} > bound {bound}"
                    );
                    // The bound is h^2/8 |c''|, the deviation of a curve of
                    // that curvature: never more than a few times it on
                    // these gently varying curves. A pivot's bound is its
                    // form's supremum over the whole piece, so over part
                    // of a sine transition it is looser.
                    // Over a whole Viennese bend a held-rail pivot bends
                    // both ways, so its chord deviates less than h^2/8 |e''|.
                    let ceiling = match name.as_str() {
                        "pivot only" => 16.0,
                        "held rail only" => 8.0,
                        _ => 4.0,
                    };
                    if worst > 1e-6 {
                        assert!(
                            bound <= ceiling * worst,
                            "{name} [{a}, {b}]: bound {bound} vs deviation {worst}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn flattening_is_certified_for_closed_form_profiles() {
    for (name, curve, length) in curves() {
        let closed = !name.ends_with("intrinsic");
        assert_eq!(certifies_flattening3(&curve), closed, "{name}");
        for tol in [1e-3, 1e-4] {
            let points = flatten3(&curve, Interval::new(0.0, length), tol, 20).unwrap();
            // Dense samples of the curve, each within `tol` of the polyline.
            for i in 0..=2000 {
                let d = length * i as Scalar / 2000.0;
                let p = evaluate3(&curve, d).unwrap();
                let near = points
                    .windows(2)
                    .map(|w| to_segment(p, w[0], w[1]))
                    .fold(Scalar::INFINITY, Scalar::min);
                assert!(
                    near <= tol * (1.0 + 1e-9),
                    "{name} tol {tol}: {near} at {d}"
                );
            }
        }
    }
}

#[test]
fn breaks_name_joins_seams_and_curvature_seams() {
    let curve = Curve3::Elevated(Elevated3::new(chain(), gradient_profile()));
    let first = continuity_breaks3(&curve, 1);
    for seam in [30.0, 40.0, 80.0, 100.0] {
        assert!(
            first.iter().any(|b| (b - seam).abs() < 1e-9),
            "{seam} in {first:?}"
        );
    }
    let on_clothoid = Curve3::Elevated(Elevated3::new(
        Curve2::Intrinsic(Intrinsic2::new(
            world(),
            CurvatureLaw::piecewise(
                vec![50.0],
                vec![
                    CurvatureLaw::clothoid(0.0, 0.004, 50.0),
                    CurvatureLaw::Constant { curvature: 0.004 },
                ],
            ),
            100.0,
        )),
        line_profile(),
    ));
    assert!(continuity_breaks3(&on_clothoid, 1).is_empty());
    assert_eq!(continuity_breaks3(&on_clothoid, 2), vec![50.0]);
    // The curvature seam does not stop the chord bound, which needs C^1.
    assert!(chord_bound3(&on_clothoid, 40.0, 60.0).is_some());
    assert!(curve_derivative_bounds3(&on_clothoid, 40.0, 60.0).is_none());
    // A profile seam does stop it: the grade may jump there.
    assert!(chord_bound3(&curve, 20.0, 35.0).is_none());
}

#[test]
fn corners_are_grade_jumps_read_from_the_laws() {
    let smooth = Curve3::Elevated(Elevated3::new(straight(), gradient_profile()));
    assert!(grade_corners3(&smooth, 1e-9).unwrap().is_empty());
    let kinked = Curve3::Elevated(Elevated3::new(
        straight(),
        ElevationLaw::Piecewise {
            breaks: vec![50.0],
            laws: vec![
                ElevationLaw::constant_grade(100.0, 0.02),
                ElevationLaw::constant_grade(101.0, 0.05),
            ],
        },
    ));
    assert_eq!(grade_corners3(&kinked, 1e-9).unwrap(), vec![50.0]);
    // A pivot whose rate jumps kinks a banked curve's point path.
    let pivot_kink = Curve3::Banked(Banked3::new(
        Elevated3::new(straight(), line_profile()),
        CantLaw::new(vec![CantPiece::constant(100.0, 0.1)]),
        CantLaw::new(vec![
            CantPiece::linear(50.0, 0.0, 0.05),
            CantPiece::constant(50.0, 0.05),
        ]),
        1.5,
        BankConvention::TangentRotation,
    ));
    assert_eq!(grade_corners3(&pivot_kink, 1e-9).unwrap(), vec![50.0]);
    assert!(grade_corners3(&banked(), 1e-9).unwrap().is_empty());
    // A held rail through a Viennese bend is smooth at the cant's seams;
    // through a linear ramp its rate jumps there, where the pivot law
    // itself has no seam.
    let held = held_rail(straight(), line_profile(), BankConvention::VerticalRise);
    assert!(grade_corners3(&held, 1e-9).unwrap().is_empty());
    let ramped = Curve3::Banked(Banked3::new(
        Elevated3::new(straight(), line_profile()),
        CantLaw::new(vec![
            CantPiece::linear(50.0, 0.0, 0.1),
            CantPiece::constant(50.0, 0.1),
        ]),
        CantLaw::new(vec![CantPiece::about_rail(100.0, RailSide::Left, 0.0)]),
        1.5,
        BankConvention::TangentRotation,
    ));
    assert_eq!(grade_corners3(&ramped, 1e-9).unwrap(), vec![50.0]);
}

/// The cant's seams are the point path's where a held-rail pivot derives
/// it from the cant (#279), and only there.
#[test]
fn a_held_rail_pivot_names_the_cant_seams_it_covers() {
    let held = held_rail(
        straight(),
        ElevationLaw::level(0.0),
        BankConvention::VerticalRise,
    );
    assert_eq!(continuity_breaks3(&held, 1), vec![20.0, 80.0]);
    assert_eq!(continuity_breaks3(&held, 2), vec![20.0, 80.0]);
    assert!(chord_bound3(&held, 10.0, 30.0).is_none());
    assert!(curve_derivative_bounds3(&held, 70.0, 90.0).is_none());
    assert!(chord_bound3(&held, 20.0, 80.0).is_some());
    // Rotating about the centreline over the bend and about the rail only
    // after it: the seam at 20 m is no longer the pivot's.
    let Curve3::Banked(mut later) = held else {
        unreachable!()
    };
    later.pivot = CantLaw::new(vec![
        CantPiece::constant(50.0, 0.0),
        CantPiece::about_rail(70.0, RailSide::Right, -0.075),
    ]);
    assert_eq!(
        continuity_breaks3(&Curve3::Banked(later), 1),
        vec![50.0, 80.0]
    );
}

/// Flattening a held-rail curve against its certified chord bound keeps
/// every dense sample within the tolerance.
#[test]
fn a_held_rail_curve_flattens_within_its_tolerance() {
    for (name, curve, length) in banked_curves().into_iter().skip(1) {
        for tol in [1e-3, 1e-4] {
            let points = flatten3(&curve, Interval::new(0.0, length), tol, 20).unwrap();
            for i in 0..=2000 {
                let d = length * i as Scalar / 2000.0;
                let p = evaluate3(&curve, d).unwrap();
                let near = points
                    .windows(2)
                    .map(|w| to_segment(p, w[0], w[1]))
                    .fold(Scalar::INFINITY, Scalar::min);
                assert!(
                    near <= tol * (1.0 + 1e-9),
                    "{name} tol {tol}: {near} at {d}"
                );
            }
        }
    }
}

#[test]
fn unbounded_and_unevaluable_cases_are_refused_by_name() {
    // A B-spline plan's parameter is not a distance: no point, no domain,
    // no bound.
    let spline = Curve3::Elevated(Elevated3::new(
        Curve2::BSpline(BSplineCurve2 {
            degree: 1,
            control_points: vec![Point2::new(0.0, 0.0), Point2::new(10.0, 0.0)],
            knots: vec![0.0, 1.0],
            multiplicities: vec![2, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }),
        line_profile(),
    ));
    assert!(evaluate3(&spline, 0.5).is_err());
    assert!(second_derivative3(&spline, 0.5).is_err());
    assert_eq!(domain3(&spline), Interval::new(0.0, 0.0));
    assert!(curve_derivative_bounds3(&spline, 0.0, 1.0).is_none());
    assert!(chord_bound3(&spline, 0.0, 1.0).is_none());
    assert!(!certifies_flattening3(&spline));
    // A sag of R = 10 from -40% turns vertical before 15 m.
    let steep = Curve3::Elevated(Elevated3::new(
        straight(),
        ElevationLaw::circular_arc(0.0, -0.4, 10.0),
    ));
    assert!(second_derivative3(&steep, 20.0).is_err());
    assert!(chord_bound3(&steep, 0.0, 20.0).is_none());
    assert!(curve_derivative_bounds3(&steep, 0.0, 20.0).is_none());
    // A banked curve's pivot may not be an angle.
    let angled = Curve3::Banked(Banked3::new(
        Elevated3::new(straight(), line_profile()),
        CantLaw::new(vec![CantPiece::constant(100.0, 0.1)]),
        CantLaw::new(vec![CantPiece::viennese_bend(100.0, 0.0, 0.05)]),
        1.5,
        BankConvention::TangentRotation,
    ));
    assert!(second_derivative3(&angled, 10.0).is_err());
    assert!(curve_derivative_bounds3(&angled, 0.0, 50.0).is_none());
    let _ = Vec3::ZERO;
}
