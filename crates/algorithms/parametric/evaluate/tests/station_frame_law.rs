//! Certified frame-law bounds for frames read by quadrature (#316): along
//! an ellipse, a plan and a 3D B-spline, offsets by distances and constant
//! offsets beside them, and pieces after joints measured by quadrature.
//!
//! The derivative bounds are checked against differences of the kernel's
//! own frames over a lattice of the span (the frames' rounding and
//! quadrature errors, about `1e-10`, sit far below the slack); the
//! measure against the arc length to the read point integrated here by a
//! fine Simpson rule; the joint error against the certified joints.

use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, CurvatureLaw, Curve2, Curve3, Ellipse2, Intrinsic2, KnotSpec, OffsetFrame,
    OffsetLaw, PathOffsets, SeamSide,
};
use axiolid_evaluate::derivative2;
use axiolid_evaluate::station::{
    offset_pieces, CompositeBasis, FrameLaw, FrameLawBounds, FrameReading, SectionFrame,
    StationCurve, StationPiece,
};

const TOLERANCE: Scalar = 1e-7;

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

/// A 3D cubic B-spline climbing gently, about 30 long.
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

/// A plan B-spline with a corner at its middle knot.
fn cornered() -> Curve2 {
    Curve2::BSpline(BSplineCurve {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 1.0),
            Point2::new(10.0, 0.0),
            Point2::new(15.0, 3.0),
            Point2::new(20.0, 3.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 2, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// A frame's axis (or point) read as a vector.
type Axis = fn(&SectionFrame) -> Vec3;

/// First and second differences of each axis and the point over a lattice
/// of `[lo, hi]`, never above the bounds.
fn assert_bounded(
    section: impl Fn(Scalar) -> SectionFrame,
    bounds: &FrameLawBounds,
    lo: Scalar,
    hi: Scalar,
) {
    let h = 2e-4 * (hi - lo);
    for i in 1..20 {
        let s = lo + (hi - lo) * Scalar::from(i) / 20.0;
        let (a, b, c) = (section(s - h), section(s), section(s + h));
        let first = |f: Axis| ((f(&c) - f(&a)) / (2.0 * h)).length();
        let second = |f: Axis| ((f(&c) - 2.0 * f(&b) + f(&a)) / (h * h)).length();
        let slack = |bound: Scalar| bound * (1.0 + 1e-3) + 1e-4;
        let axes: [(Axis, [Scalar; 2], &str); 4] = [
            (|f| f.point - Point3::ZERO, bounds.point, "point"),
            (|f| f.tangent, bounds.tangent, "tangent"),
            (|f| f.lateral, bounds.lateral, "lateral"),
            (|f| f.up, bounds.up, "up"),
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

/// A plane curve's arc length over `[t0, t1]` by Simpson's rule.
fn arc(curve: &Curve2, t0: Scalar, t1: Scalar) -> Scalar {
    let steps = 200_000;
    let h = (t1 - t0) / Scalar::from(steps);
    let speed = |t: Scalar| derivative2(curve, t).unwrap().length();
    (0..steps)
        .map(|i| {
            let t = t0 + h * Scalar::from(i);
            h / 6.0 * (speed(t) + 4.0 * speed(t + 0.5 * h) + speed(t + h))
        })
        .sum()
}

/// The ellipse's angle at its point `p`, in `[0, 2 pi)`.
fn angle_of(p: Point3) -> Scalar {
    ((p.y / 6.0).atan2((p.x - 30.0) / 10.0)).rem_euclid(core::f64::consts::TAU)
}

#[test]
fn an_ellipse_is_bounded_and_its_measure_holds_the_reading() {
    let curve = ellipse();
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    for (lo, hi) in [(2.0, 6.0), (20.0, 22.0), (40.0, 50.0)] {
        let bounds = law.bounds(lo, hi, FrameReading::Section).unwrap();
        assert_bounded(
            |s| axiolid_evaluate::station::station_section2(&curve, s).unwrap(),
            &bounds,
            lo,
            hi,
        );
        // The measure holds the true distance of each read, and stays near
        // the tolerance.
        for s in [lo, hi] {
            let read = axiolid_evaluate::station::station_section2(&curve, s).unwrap();
            let truth = arc(&curve, 0.0, angle_of(read.point));
            assert!(
                (truth - s).abs() <= bounds.measure,
                "{truth} vs {s}: {}",
                bounds.measure
            );
        }
        assert!(bounds.measure <= 10.0 * TOLERANCE, "{}", bounds.measure);
        assert!(bounds.point_error <= 1e-12 && bounds.axis_error <= 1e-8);
    }
}

/// A cubic `y = c (x - 10)^3` over `x` in `[0, 20]`, `c = 5e-4`, a
/// Bezier: at its inflection the tangent's second derivative is the
/// curvature's rate, `6 c`, above any bound from the curvature alone.
fn inflected() -> Curve2 {
    Curve2::BSpline(BSplineCurve {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, -0.5),
            Point2::new(20.0 / 3.0, 0.5),
            Point2::new(40.0 / 3.0, -0.5),
            Point2::new(20.0, 0.5),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

#[test]
fn an_inflection_bounds_the_tangent_by_its_curvature_rate() {
    let curve = inflected();
    let half = 0.5
        * axiolid_evaluate::station::station_length2(&curve)
            .unwrap()
            .unwrap();
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    let (lo, hi) = (half - 0.5, half + 0.5);
    let bounds = law.bounds(lo, hi, FrameReading::Section).unwrap();
    assert_bounded(
        |s| axiolid_evaluate::station::station_section2(&curve, s).unwrap(),
        &bounds,
        lo,
        hi,
    );
    assert!(bounds.tangent[1] >= 6.0 * 5e-4, "{:?}", bounds.tangent);
}

#[test]
fn a_reversed_ellipse_piece_measures_from_its_own_start() {
    let curve = ellipse();
    let piece = StationPiece::between(StationCurve::Two(&curve), 5.0, 25.0)
        .unwrap()
        .reversed();
    let basis = CompositeBasis::new(vec![piece]).unwrap();
    let law = basis.frame_law(TOLERANCE).unwrap();
    let bounds = law.bounds(3.0, 7.0, FrameReading::Section).unwrap();
    assert_bounded(
        |s| basis.section_on(s, SeamSide::Outgoing).unwrap(),
        &bounds,
        3.0,
        7.0,
    );
    assert!(bounds.measure <= 10.0 * TOLERANCE, "{}", bounds.measure);
}

#[test]
fn a_3d_b_spline_is_bounded_in_its_section_and_plan_frames() {
    let curve = spline3();
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let bounds = law.bounds(4.0, 5.5, reading).unwrap();
        assert!(bounds.up[0] > 0.0 || reading == FrameReading::Plan);
        assert_bounded(
            |s| {
                let f = axiolid_evaluate::station::station_section3(&curve, s).unwrap();
                match reading {
                    FrameReading::Section => f,
                    _ => f.plan().unwrap(),
                }
            },
            &bounds,
            4.0,
            5.5,
        );
    }
}

#[test]
fn a_corner_knot_reads_its_sides_exactly() {
    let curve = cornered();
    let seams = axiolid_evaluate::station::station_seams2(&curve).unwrap();
    let corner = seams.iter().find(|s| !s.smooth).unwrap().distance;
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    for (lo, hi) in [(corner - 2.0, corner), (corner, corner + 2.0)] {
        let bounds = law.bounds(lo, hi, FrameReading::Section).unwrap();
        assert!(bounds.axis_error <= 1e-8, "{}", bounds.axis_error);
    }
}

fn widening() -> OffsetLaw {
    OffsetLaw::Linear {
        start: PathOffsets::new(0.5, 0.0, 0.0),
        end: PathOffsets::new(0.75, 0.2, 0.0),
        frame: OffsetFrame::Section,
    }
}

#[test]
fn offsets_beside_an_ellipse_and_a_b_spline_are_bounded() {
    let (plan, solid) = (ellipse(), spline3());
    let cases: [(StationCurve<'_>, OffsetLaw, Scalar, Scalar); 4] = [
        (
            StationCurve::Two(&plan),
            OffsetLaw::Planar { distance: 0.5 },
            0.0,
            30.0,
        ),
        (StationCurve::Two(&plan), widening(), 2.0, 20.0),
        (StationCurve::Three(&solid), widening(), 1.0, 28.0),
        (
            StationCurve::Three(&solid),
            OffsetLaw::Directed {
                distance: 1.0,
                reference_direction: Vec3::Z,
            },
            0.0,
            28.0,
        ),
    ];
    for (k, (curve, offset_law, start, end)) in cases.into_iter().enumerate() {
        let base = StationPiece::between(curve, start, end).unwrap();
        // The second offset by distances runs along its base backwards.
        let base = if k == 1 { base.reversed() } else { base };
        let basis = CompositeBasis::new(offset_pieces(&[base], offset_law).unwrap()).unwrap();
        let law = basis.frame_law(TOLERANCE).unwrap();
        let (lo, hi) = (3.0, 4.5);
        let bounds = law.bounds(lo, hi, FrameReading::Section).unwrap();
        assert_bounded(
            |s| basis.section_on(s, SeamSide::Outgoing).unwrap(),
            &bounds,
            lo,
            hi,
        );
        assert!(bounds.measure <= 1e-6, "{}", bounds.measure);
        assert!(bounds.point_error <= 1e-9, "{}", bounds.point_error);
        assert!(bounds.axis_error <= 1e-7, "{}", bounds.axis_error);
    }
}

#[test]
fn a_piece_after_a_quadrature_carries_its_joint_and_an_intrinsic_one_is_unbounded() {
    let plan = ellipse();
    // Two spans of the ellipse end to end: the second starts where the
    // kernel sums the first's estimated length.
    let first = StationPiece::between(StationCurve::Two(&plan), 0.0, 8.0).unwrap();
    let second = StationPiece::between(StationCurve::Two(&plan), 8.0, 16.0).unwrap();
    let basis = CompositeBasis::new(vec![first, second]).unwrap();
    let joint = basis.pieces()[0].length();
    let law = basis.frame_law(TOLERANCE).unwrap();
    let bounds = law
        .bounds(joint + 1.0, joint + 4.0, FrameReading::Section)
        .unwrap();
    let certified = basis.joint_intervals(TOLERANCE).unwrap()[0];
    let error = (joint - certified.lower).max(certified.upper - joint);
    // The second piece's own reading, on its own: the joint's error comes
    // on top of it.
    let alone = CompositeBasis::new(vec![second]).unwrap();
    let own = alone
        .frame_law(TOLERANCE)
        .unwrap()
        .bounds(1.0, 4.0, FrameReading::Section)
        .unwrap()
        .measure;
    assert!(error > 0.0);
    assert!(
        bounds.measure >= own + error * (1.0 - 1e-9),
        "{} < {own} + {error}",
        bounds.measure
    );
    assert!(bounds.measure <= 10.0 * TOLERANCE, "{}", bounds.measure);

    let clothoid = Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
        60.0,
    ));
    let law = FrameLaw::atomic(StationCurve::Two(&clothoid), None, TOLERANCE).unwrap();
    assert!(law.bounds(5.0, 10.0, FrameReading::Section).is_none());
    let base = StationPiece::between(StationCurve::Two(&clothoid), 0.0, 60.0).unwrap();
    let beside =
        CompositeBasis::new(offset_pieces(&[base], OffsetLaw::Planar { distance: 2.0 }).unwrap())
            .unwrap();
    let law = beside.frame_law(TOLERANCE).unwrap();
    assert!(law.bounds(5.0, 10.0, FrameReading::Section).is_none());
}
