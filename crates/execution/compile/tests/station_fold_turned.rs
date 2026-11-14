//! Folds of turned sections are decided exactly (#406, ADR 0082 amendment):
//! the rate of a section point along the turned plane's normal, `F = sigma
//! n_T + rise n_U + q'.n + w.(p x n) + x c_x + y c_y`, from the frame's
//! twist, rise and bend and the orientation's closed-form turn.
//!
//! The closed forms are a yawed section on a plan arc (`F = cos(a)(1 - l /
//! R) - x / R`), a pitched one on a vertical arc (`F = cos(t)(1 - b / R) - y
//! / R`) and one pitched against a grade in the plan frame (`F = cos(t) + G
//! sin(t)`). The oracle of the rest rebuilds the placed vertices `A(s)` from
//! the kernel's own frames, the orientation interpolated as the kernel does,
//! and differentiates them: `F = A' . N_t` by a central difference
//! (`h = 1e-5`), independently of the jets the decision rests on.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BankConvention, Banked3, CantLaw, CantPiece, Circle2, Circle3, CurvatureLaw,
    Curve2, Curve3, Elevated3, ElevationLaw, Intrinsic2, KnotSpec, Line2,
};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, SectionAtStation, SolidOperation,
    Station, StationFrame, StationOffsets, StationOrientation,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::derivative3;
use axiolid_reference::station::{station_section2, station_section3, SectionFrame};

const BUDGET: Scalar = 0.01;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(BUDGET)
        .unwrap()
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

/// What the mesh compiler makes of a run.
#[derive(Debug)]
enum Outcome {
    Certified,
    Unbounded(&'static str),
    Refused(String),
}

fn mesh(graph: &GeometryGraph, root: NodeId) -> Outcome {
    match compiler().compile_mesh_with_deviation(graph, root, &options()) {
        Ok((_, report)) => {
            let run = report
                .contributions
                .iter()
                .find(|c| c.path == DeviationPath::StationedSpine)
                .expect("the run reports its path");
            match run.bound {
                DeviationBound::Certified(_) => Outcome::Certified,
                DeviationBound::Unbounded(reason) => Outcome::Unbounded(reason),
                other => panic!("neither certified nor unbounded: {other:?}"),
            }
        }
        Err(GeomError::InvalidInput(detail)) => Outcome::Refused(detail),
        Err(other) => panic!("not a refusal as invalid input: {other:?}"),
    }
}

/// `Ok` where the exact compiler builds the run, `Err` with the detail of
/// its invalid-input refusal.
fn exact(graph: &GeometryGraph, root: NodeId) -> Result<(), String> {
    match ReferenceExactCompiler::new().compile_exact(graph, root, &options()) {
        Ok(_) => Ok(()),
        Err(GeomError::InvalidInput(detail)) => Err(detail),
        Err(other) => panic!("not built nor refused as invalid input: {other:?}"),
    }
}

fn certified(outcome: &Outcome) {
    match outcome {
        Outcome::Certified => {}
        Outcome::Unbounded(reason) => panic!("unbounded: {reason}"),
        Outcome::Refused(detail) => panic!("refused: {detail}"),
    }
}

fn refused(outcome: &Outcome) {
    match outcome {
        Outcome::Refused(detail) => assert!(detail.contains("sections fold"), "{detail}"),
        Outcome::Unbounded(reason) => panic!("unbounded: {reason}"),
        Outcome::Certified => panic!("certified"),
    }
}

fn line_x() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })
}

fn plan_circle(radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

/// A section: the profile `x` by `y` at `distance`, the offsets, an orientation.
struct Section {
    distance: Scalar,
    orientation: StationOrientation,
}

/// A run of `profile` rectangles (width `x`, height `y`) along `curve`, with
/// the same offsets throughout and the given sections.
fn run(
    curve: impl Into<GeometryNode>,
    frame: StationFrame,
    profile: (Scalar, Scalar),
    offsets: StationOffsets,
    sections: &[Section],
) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(curve.into()).unwrap();
    let rectangle = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: profile.0,
            y: profile.1,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: sections
                    .iter()
                    .map(|s| {
                        SectionAtStation::new(rectangle, Station::new(s.distance, offsets))
                            .with_orientation(s.orientation)
                    })
                    .collect(),
                frame,
            },
        ))
        .unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

/// The base orientation.
const BASE: StationOrientation = StationOrientation::new(None, None);

/// The section pitched about its lateral axis by `degrees`.
fn pitched(degrees: Scalar) -> StationOrientation {
    let t = degrees.to_radians();
    StationOrientation::new(
        Some(Vec3::new(-t.sin(), 0.0, t.cos())),
        Some(Vec3::new(t.cos(), 0.0, t.sin())),
    )
}

/// The section yawed about its up axis by `degrees`.
fn yawed(degrees: Scalar) -> StationOrientation {
    let t = degrees.to_radians();
    StationOrientation::new(None, Some(Vec3::new(t.cos(), t.sin(), 0.0)))
}

fn both(distance: Scalar, orientation: StationOrientation) -> [Section; 2] {
    [
        Section {
            distance: 0.0,
            orientation,
        },
        Section {
            distance,
            orientation,
        },
    ]
}

/// The smallest rate along the turned plane's normal of the corners `(x, y)`
/// of the section over a station interval `[0, length]` between the
/// orientations `a` and `b`, by a central difference of the placed corners.
fn sampled_min_rate(
    frame: &dyn Fn(Scalar) -> SectionFrame,
    (a, b): (StationOrientation, StationOrientation),
    offsets: StationOffsets,
    (width, height): (Scalar, Scalar),
    length: Scalar,
) -> Scalar {
    let h = 1e-5;
    let (axis_a, ref_a) = a.unit_axes().unwrap();
    let (axis_b, ref_b) = b.unit_axes().unwrap();
    let turned = |s: Scalar| {
        let u = s / length;
        frame(s)
            .oriented(Some(axis_a.lerp(axis_b, u)), Some(ref_a.lerp(ref_b, u)))
            .unwrap()
    };
    let corner = |s: Scalar, (x, y): (Scalar, Scalar)| -> Point3 {
        let turned = turned(s);
        frame(s).place(offsets.lateral, offsets.vertical, offsets.longitudinal)
            + x * turned.lateral
            + y * turned.up
    };
    let corners = [
        (-0.5 * width, -0.5 * height),
        (0.5 * width, -0.5 * height),
        (0.5 * width, 0.5 * height),
        (-0.5 * width, 0.5 * height),
    ];
    let mut worst = Scalar::INFINITY;
    for i in 0..=200 {
        let s = h + (length - 2.0 * h) * Scalar::from(i) / 200.0;
        let normal = turned(s).tangent;
        for c in corners {
            let rate = (corner(s + h, c) - corner(s - h, c)) / (2.0 * h);
            worst = worst.min(rate.dot(normal));
        }
    }
    worst
}

/// The offset `b` in `[lo, hi]` at which `minimum(b)` falls through `target`,
/// by bisection (the minimum falls as `b` grows).
fn offset_at(
    minimum: &dyn Fn(Scalar) -> Scalar,
    (mut lo, mut hi): (Scalar, Scalar),
    target: Scalar,
) -> Scalar {
    assert!(
        minimum(lo) > target && minimum(hi) < target,
        "{} {}",
        minimum(lo),
        minimum(hi)
    );
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if minimum(mid) > target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

const FRAMES: [StationFrame; 2] = [StationFrame::Section, StationFrame::Plan];

#[test]
fn a_yawed_section_folds_at_its_closed_form_radius() {
    // `F = cos(20)(1 - a / 10) - x / 10` for the corner `x = 1/2`.
    let yaw = yawed(20.0);
    for frame in FRAMES {
        let offsets = |lateral| StationOffsets::new(lateral, 0.0, 0.0);
        let case = |lateral| {
            run(
                plan_circle(10.0),
                frame,
                (1.0, 1.0),
                offsets(lateral),
                &both(10.0, yaw),
            )
        };
        let closed = |lateral: Scalar| 20.0_f64.to_radians().cos() * (1.0 - lateral / 10.0) - 0.05;
        let oracle = |lateral| {
            sampled_min_rate(
                &|s| station_section2(&plan_circle(10.0), s).unwrap(),
                (yaw, yaw),
                offsets(lateral),
                (1.0, 1.0),
                10.0,
            )
        };
        assert!((closed(9.4) - 0.0064).abs() < 1e-4);
        for lateral in [9.4, 9.48, 9.5] {
            assert!(
                (oracle(lateral) - closed(lateral)).abs() < 1e-6,
                "{lateral}"
            );
        }
        let (graph, root) = case(9.4);
        certified(&mesh(&graph, root));
        // Turned, 9.48 folds although the unturned section would not
        // (`1 - 9.98 / 10 > 0`).
        let (graph, root) = case(9.48);
        refused(&mesh(&graph, root));
        let (graph, root) = case(9.5);
        refused(&mesh(&graph, root));
        // The exact compiler refuses the fold before the turned circle.
        let detail = exact(&graph, root).expect_err("the exact compiler refuses alike");
        assert!(detail.contains("sections fold"), "{frame:?}: {detail}");
    }
}

#[test]
fn a_pitched_section_folds_on_a_vertical_arc() {
    // A circle of radius 10 in the XZ plane, centre `(0, 0, 10)`, running up
    // from the origin: the section's up points at the centre, and `F =
    // cos(30)(1 - b / 10) - y / 10` for the corner `y = 1/2`. The run stops
    // short of the vertical tangent.
    let circle = || {
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(0.0, 0.0, 10.0),
                x: -Vec3::Z,
                y: Vec3::X,
                z: -Vec3::Y,
            },
            radius: 10.0,
        })
    };
    let pitch = pitched(30.0);
    let case = |vertical| {
        run(
            circle(),
            StationFrame::Section,
            (1.0, 1.0),
            StationOffsets::new(0.0, vertical, 0.0),
            &both(5.0, pitch),
        )
    };
    let closed = |vertical: Scalar| 30.0_f64.to_radians().cos() * (1.0 - vertical / 10.0) - 0.05;
    for vertical in [9.3, 9.45, 9.5] {
        let oracle = sampled_min_rate(
            &|s| station_section3(&circle(), s).unwrap(),
            (pitch, pitch),
            StationOffsets::new(0.0, vertical, 0.0),
            (1.0, 1.0),
            5.0,
        );
        assert!((oracle - closed(vertical)).abs() < 1e-6, "{vertical}");
    }
    let (graph, root) = case(9.3);
    certified(&mesh(&graph, root));
    // Turned, 9.45 folds although the unturned section would not.
    let (graph, root) = case(9.45);
    refused(&mesh(&graph, root));
    let (graph, root) = case(9.5);
    refused(&mesh(&graph, root));
}

#[test]
fn a_section_pitched_against_the_grade_in_the_plan_frame_runs_backwards() {
    // The plan frame's up never leaves `+Z`, so the grade's climb is the
    // basis point's rise: `F = cos(t) + G sin(t)` on a grade `G = 1/2`.
    let graded = || {
        Curve3::Elevated(Elevated3::new(
            line_x(),
            ElevationLaw::Polynomial {
                coefficients: vec![0.0, 0.5],
            },
        ))
    };
    for (degrees, closed) in [(-60.0_f64, 0.067), (-70.0, -0.128)] {
        let pitch = pitched(degrees);
        let oracle = sampled_min_rate(
            &|s| station_section3(&graded(), s).unwrap().plan().unwrap(),
            (pitch, pitch),
            StationOffsets::ZERO,
            (1.0, 1.0),
            10.0,
        );
        assert!((oracle - closed).abs() < 1e-3, "{degrees}: {oracle}");
        let (graph, root) = run(
            graded(),
            StationFrame::Plan,
            (1.0, 1.0),
            StationOffsets::ZERO,
            &both(10.0, pitch),
        );
        let outcome = mesh(&graph, root);
        if closed > 0.0 {
            certified(&outcome);
        } else {
            refused(&outcome);
        }
    }
}

/// A line with a section that turns between the stations: `base` at 0 and
/// `turned` at 10, a rectangle `profile`.
fn turning(
    turned: StationOrientation,
    profile: (Scalar, Scalar),
) -> ((GeometryGraph, NodeId), Scalar) {
    let sections = [
        Section {
            distance: 0.0,
            orientation: BASE,
        },
        Section {
            distance: 10.0,
            orientation: turned,
        },
    ];
    let case = run(
        line_x(),
        StationFrame::Section,
        profile,
        StationOffsets::ZERO,
        &sections,
    );
    let oracle = sampled_min_rate(
        &|s| station_section2(&line_x(), s).unwrap(),
        (BASE, turned),
        StationOffsets::ZERO,
        profile,
        10.0,
    );
    (case, oracle)
}

#[test]
fn a_section_pitching_between_stations_folds_where_its_top_swings_back() {
    // A tall section: the corner `y = 3` of the pitching plane runs back
    // faster than the basis advances before the far station, where a middle
    // frame alone would clear it.
    let ((graph, root), oracle) = turning(pitched(80.0), (1.0, 6.0));
    assert!(oracle < -0.05, "{oracle}");
    refused(&mesh(&graph, root));
    let ((graph, root), oracle) = turning(pitched(50.0), (1.0, 6.0));
    assert!(oracle > 0.3, "{oracle}");
    certified(&mesh(&graph, root));
}

#[test]
fn a_section_swinging_past_a_right_angle_and_back_folds_at_its_station() {
    // A tiny section pitched or yawed to 91 degrees at station 4 and back to
    // the base at 10 runs backwards on the last few centimetres before
    // station 4 (`F = cos(91)` there) while the middle of the stretch before
    // it advances, and the section's own extent adds next to nothing: the
    // turn of the plane's normal over a span must widen it, or a span
    // stretching into the fold is cleared by its middle.
    for turned in [pitched(91.0), yawed(91.0)] {
        let sections = [
            Section {
                distance: 0.0,
                orientation: BASE,
            },
            Section {
                distance: 4.0,
                orientation: turned,
            },
            Section {
                distance: 10.0,
                orientation: BASE,
            },
        ];
        let (graph, root) = run(
            line_x(),
            StationFrame::Section,
            (0.02, 0.02),
            StationOffsets::ZERO,
            &sections,
        );
        let oracle = sampled_min_rate(
            &|s| station_section2(&line_x(), s).unwrap(),
            (BASE, turned),
            StationOffsets::ZERO,
            (0.02, 0.02),
            4.0,
        );
        assert!(oracle < -0.01, "{oracle}");
        refused(&mesh(&graph, root));
    }
}

#[test]
fn a_section_yawing_between_stations_folds_alike() {
    let ((graph, root), oracle) = turning(yawed(80.0), (6.0, 1.0));
    assert!(oracle < -0.05, "{oracle}");
    refused(&mesh(&graph, root));
    let ((graph, root), oracle) = turning(yawed(50.0), (6.0, 1.0));
    assert!(oracle > 0.3, "{oracle}");
    certified(&mesh(&graph, root));
}

/// A circle of radius 5 whose plane is tilted 30 degrees about `x`: its
/// reference-up frame twists.
fn tilted_circle() -> Curve3 {
    let (sin, cos) = (30.0_f64).to_radians().sin_cos();
    let (x, y) = (Vec3::X, Vec3::new(0.0, cos, sin));
    Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(0.0, 0.0, 1.0),
            x,
            y,
            z: x.cross(y),
        },
        radius: 5.0,
    })
}

/// A run whose offsets are chosen so that the oracle's smallest rate is
/// `+-margin`: the one refused and the one certified.
fn near_a_fold(
    curve: &dyn Fn() -> Curve3,
    orientation: StationOrientation,
    offsets: &dyn Fn(Scalar) -> StationOffsets,
    range: (Scalar, Scalar),
    length: Scalar,
) {
    let minimum = |b: Scalar| {
        sampled_min_rate(
            &|s| station_section3(&curve(), s).unwrap(),
            (orientation, orientation),
            offsets(b),
            (1.0, 1.0),
            length,
        )
    };
    let clear = offset_at(&minimum, range, 0.02);
    let fold = offset_at(&minimum, range, -0.02);
    let build = |b: Scalar| {
        run(
            curve(),
            StationFrame::Section,
            (1.0, 1.0),
            offsets(b),
            &both(length, orientation),
        )
    };
    let (graph, root) = build(clear);
    certified(&mesh(&graph, root));
    let (graph, root) = build(fold);
    refused(&mesh(&graph, root));
}

#[test]
fn a_turned_section_on_a_twisting_frame_is_judged_by_its_rate() {
    // The tilted circle's frame twists, so the turned plane's normal is not
    // the bend's alone: the up offset at which the oracle's smallest rate
    // is +-0.02 is cleared and refused.
    near_a_fold(
        &tilted_circle,
        yawed(30.0),
        &|b| StationOffsets::new(0.0, b, 0.0),
        (0.0, 40.0),
        6.0,
    );
}

/// The canted transition of `station_fold.rs`: a clothoid of curvature
/// `s / 1000` under a constant cant of 0.15 on a 1.5 m gauge.
fn canted() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(
            Curve2::Intrinsic(Intrinsic2::new(
                Frame2 {
                    origin: Point2::ZERO,
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                CurvatureLaw::clothoid(0.0, 0.1, 100.0),
                100.0,
            )),
            ElevationLaw::level(0.0),
        ),
        CantLaw::new(vec![CantPiece::constant(100.0, 0.15)]),
        CantLaw::zero(100.0),
        1.5,
        BankConvention::TangentRotation,
    ))
}

#[test]
fn a_pitched_section_on_a_canted_transition_is_judged_by_its_roll_rate() {
    // The roll's rate is part of the twist: a section pitched 40 degrees on
    // the rolled frame, the lateral offset chosen by the oracle.
    near_a_fold(
        &canted,
        pitched(40.0),
        &|b| StationOffsets::new(b, 0.0, 0.0),
        (0.0, 40.0),
        100.0,
    );
}

/// A quadratic space B-spline climbing at a grade of about 0.3, whose
/// curvature jumps at its interior knot.
fn climbing_quadratic() -> Curve3 {
    Curve3::BSpline(BSplineCurve {
        degree: 2,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(12.0, 0.0, 3.6),
            Point3::new(12.0, 12.0, 7.2),
            Point3::new(24.0, 14.0, 10.8),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 1, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// The arc length of `climbing_quadratic` to its knot, by Simpson's rule.
fn knot_distance() -> Scalar {
    let curve = climbing_quadratic();
    let n = 2000;
    let speed = |u: Scalar| derivative3(&curve, u).unwrap().length();
    let h = 0.5 / n as Scalar;
    (0..n)
        .map(|i| {
            let a = h * i as Scalar;
            h / 6.0 * (speed(a) + 4.0 * speed(a + 0.5 * h) + speed(a + h))
        })
        .sum()
}

#[test]
fn a_turned_section_across_a_curvature_break_is_cleared_by_the_conservative_bound() {
    // The spans across the knot have no rates, so the conservative bound
    // decides them: in the plan frame a section pitched down 50 degrees
    // advances at `rho cos(50) - |O'_z| sin(50)`, which it clears.
    let knot = knot_distance();
    let pitch = pitched(-50.0);
    let (graph, root) = run(
        climbing_quadratic(),
        StationFrame::Plan,
        (1.0, 1.0),
        StationOffsets::ZERO,
        &[
            Section {
                distance: knot - 3.0,
                orientation: pitch,
            },
            Section {
                distance: knot + 3.0,
                orientation: pitch,
            },
        ],
    );
    certified(&mesh(&graph, root));
}

/// A 1 m square run along the line `x` with a turned section at stations 0
/// and 10, whose offsets run from zero to `end`.
fn sliding(
    frame: StationFrame,
    orientation: StationOrientation,
    end: StationOffsets,
) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(line_x().into()).unwrap();
    let square = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 1.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: [(0.0, StationOffsets::ZERO), (10.0, end)]
                    .map(|(distance, offsets)| {
                        SectionAtStation::new(square, Station::new(distance, offsets))
                            .with_orientation(orientation)
                    })
                    .to_vec(),
                frame,
            },
        ))
        .unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

#[test]
fn a_turned_section_is_judged_by_the_rates_of_its_offsets() {
    // On a line the plane's normal `n` of a section yawed (pitched) by 30
    // degrees has `n_T = cos(30)` and `n_L = sin(30)` (`n_U = sin(30)`), so
    // `F = cos(30) + q_L' sin(30)` (`+ q_U' sin(30)`): -0.084 at a rate of
    // -1.9 and +0.116 at -1.5, whatever the other offset does. Dropping the
    // rate of the lateral or the vertical offset from `F` leaves a section
    // that runs backwards certified.
    let lateral = |rate: Scalar| StationOffsets::new(10.0 * rate, 0.0, 0.0);
    let vertical = |rate: Scalar| StationOffsets::new(0.0, 10.0 * rate, 0.0);
    for frame in FRAMES {
        let cases: [(&str, StationOrientation, StationOffsets, bool); 8] = [
            ("yaw, lateral -1.9", yawed(30.0), lateral(-1.9), false),
            ("yaw, lateral -1.5", yawed(30.0), lateral(-1.5), true),
            ("yaw, lateral +1.9", yawed(30.0), lateral(1.9), true),
            ("yaw, vertical -1.9", yawed(30.0), vertical(-1.9), true),
            ("pitch, vertical -1.9", pitched(30.0), vertical(-1.9), false),
            ("pitch, vertical -1.5", pitched(30.0), vertical(-1.5), true),
            ("pitch, vertical +1.9", pitched(30.0), vertical(1.9), true),
            ("pitch, lateral -1.9", pitched(30.0), lateral(-1.9), true),
        ];
        for (name, orientation, end, clear) in cases {
            let (graph, root) = sliding(frame, orientation, end);
            let outcome = mesh(&graph, root);
            if clear {
                assert!(
                    matches!(outcome, Outcome::Certified),
                    "{frame:?} {name}: {outcome:?}"
                );
            } else {
                assert!(
                    matches!(&outcome, Outcome::Refused(d) if d.contains("sections fold")),
                    "{frame:?} {name}: {outcome:?}"
                );
            }
        }
    }
}
