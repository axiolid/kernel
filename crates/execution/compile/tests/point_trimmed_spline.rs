//! A B-spline trimmed by points in the graph (#306) resolves its points
//! by the one inversion a curve path's point bound takes (#304): the trim
//! flattens to exactly the path a caller builds with `PathBound::Point`
//! bounds, and its stations read bitwise as that path's, for a 2D and a 3D
//! B-spline, with and against the trim's sense. A point off the B-spline,
//! or one it passes through twice, is refused by the #304 names.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar};
use axiolid_curve::{
    BSplineCurve, Curve2, Curve3, CurvePath, KnotSpec, PathBound, PathCurve, PathPiece, SeamSide,
};
use axiolid_mesh_compile::station::{curve_path, resolve};
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraphBuilder, GeometryNode, NodeId, OrientedCurveStation,
    Station, TrimSelector, TrimmingPreference,
};
use axiolid_reference::curve::{evaluate2, evaluate3};
use axiolid_reference::station::{CompositeBasis, PATH_POINT_AMBIGUOUS, PATH_POINT_OFF_CURVE};

const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];

fn spline<P: Clone>(points: &[P], knots: Vec<Scalar>, multiplicities: Vec<u32>) -> BSplineCurve<P> {
    BSplineCurve {
        degree: 3,
        control_points: points.to_vec(),
        knots,
        multiplicities,
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    }
}

/// A planar cubic, one span.
fn plan_spline() -> Curve2 {
    Curve2::BSpline(spline(
        &[
            Point2::new(0.0, 0.0),
            Point2::new(30.0, 0.0),
            Point2::new(60.0, 20.0),
            Point2::new(90.0, 20.0),
        ],
        vec![0.0, 1.0],
        vec![4, 4],
    ))
}

/// A space cubic with an interior knot, rising as it turns.
fn space_spline() -> Curve3 {
    Curve3::BSpline(spline(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(20.0, 5.0, 2.0),
            Point3::new(40.0, 30.0, 6.0),
            Point3::new(60.0, 35.0, 9.0),
            Point3::new(80.0, 10.0, 12.0),
        ],
        vec![0.0, 0.5, 1.0],
        vec![4, 1, 4],
    ))
}

/// A planar cubic loop whose ends meet at the origin.
fn plan_loop() -> Curve2 {
    Curve2::BSpline(spline(
        &[
            Point2::new(0.0, 0.0),
            Point2::new(40.0, 30.0),
            Point2::new(-40.0, 30.0),
            Point2::new(0.0, 0.0),
        ],
        vec![0.0, 1.0],
        vec![4, 4],
    ))
}

/// The same loop in space, tilted out of its plane.
fn space_loop() -> Curve3 {
    Curve3::BSpline(spline(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(40.0, 30.0, 10.0),
            Point3::new(-40.0, 30.0, 10.0),
            Point3::new(0.0, 0.0, 0.0),
        ],
        vec![0.0, 1.0],
        vec![4, 4],
    ))
}

fn trim(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    start: TrimSelector,
    end: TrimSelector,
    sense: bool,
    preference: TrimmingPreference,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![start],
        end: vec![end],
        sense_agreement: sense,
        preference,
    }))
    .unwrap()
}

fn at(basis: NodeId, d: Scalar, side: SeamSide) -> OrientedCurveStation {
    OrientedCurveStation::from(CurveStation::new(basis, Station::at(d))).with_seam_side(side)
}

/// The trim of `basis` between the points `start` and `end` (as
/// selectors, and as the path bounds they lift to) agrees bitwise with
/// the piece bounded by those points: its flattened path, and stations
/// along it on both sides, under both trimming preferences that read a
/// point, with and against the trim's sense. A trim against its sense
/// runs from the later point back, as a reversed bounded piece does.
fn agrees(curve: PathCurve, node: GeometryNode, points: [(TrimSelector, Point3); 2], what: &str) {
    for preference in [
        TrimmingPreference::Cartesian,
        TrimmingPreference::Unspecified,
    ] {
        for sense in [true, false] {
            let [(early, p0), (late, p1)] = points;
            let (start, end) = if sense { (early, late) } else { (late, early) };
            let mut b = GeometryGraphBuilder::new();
            let basis = b.push(node.clone()).unwrap();
            let trimmed = trim(&mut b, basis, start, end, sense, preference);
            let piece =
                PathPiece::bounded(curve.clone(), PathBound::Point(p0), PathBound::Point(p1));
            let bounded = CurvePath::from(if sense { piece } else { piece.reversed() });
            let composite = CompositeBasis::from_path(&bounded).unwrap();
            let length = composite.length();
            let distances = [0.0, 0.1 * length, 0.5 * length, 0.93 * length, length];
            let mut stations = Vec::new();
            for &d in &distances {
                for side in SIDES {
                    stations.push((d, side, b.push_value(at(trimmed, d, side)).unwrap()));
                }
            }
            let roots = stations.iter().map(|&(_, _, id)| id).collect();
            let graph = b.finish(roots).unwrap();
            let label = format!("{what}, {preference:?}, sense {sense}");
            assert_eq!(
                curve_path(&graph, trimmed).unwrap(),
                composite.path(),
                "{label}: path"
            );
            for (d, side, id) in stations {
                assert_eq!(
                    resolve(&graph, id).unwrap().section,
                    composite.section_on(d, side).unwrap(),
                    "{label}: station {side:?} at {d}"
                );
            }
        }
    }
}

#[test]
fn a_planar_b_spline_trimmed_by_points_reads_as_its_point_bounded_piece() {
    let curve = plan_spline();
    let on = |t: Scalar| evaluate2(&curve, t).unwrap();
    let selector = |t: Scalar| {
        let p = on(t);
        (TrimSelector::Point2(p), Point3::new(p.x, p.y, 0.0))
    };
    agrees(
        PathCurve::Two(curve.clone()),
        GeometryNode::Curve2(curve.clone()),
        [selector(0.2), selector(0.7)],
        "planar B-spline",
    );
}

#[test]
fn a_space_b_spline_trimmed_by_points_reads_as_its_point_bounded_piece() {
    let curve = space_spline();
    let selector = |t: Scalar| {
        let p = evaluate3(&curve, t).unwrap();
        (TrimSelector::Point3(p), p)
    };
    // Across the interior knot.
    agrees(
        PathCurve::Three(curve.clone()),
        GeometryNode::Curve3(curve.clone()),
        [selector(0.15), selector(0.8)],
        "space B-spline",
    );
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, needle: &str, what: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert_eq!(input, needle, "{what}");
        }
        other => panic!("{what}: expected a refusal naming {needle:?}, got {other:?}"),
    }
}

/// The path of a trim of `node` from `start` to `end`.
fn trimmed_path(
    node: GeometryNode,
    start: TrimSelector,
    end: TrimSelector,
    preference: TrimmingPreference,
) -> GeomResult<CurvePath> {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(node).unwrap();
    let trimmed = trim(&mut b, basis, start, end, true, preference);
    let graph = b.finish(vec![trimmed]).unwrap();
    curve_path(&graph, trimmed)
}

#[test]
fn a_point_off_a_b_spline_or_on_it_twice_is_refused_by_name() {
    let plan = plan_spline();
    let space = space_spline();
    let p = evaluate2(&plan, 0.4).unwrap();
    let q = evaluate2(&plan, 0.9).unwrap();
    let s = evaluate3(&space, 0.4).unwrap();
    let r = evaluate3(&space, 0.9).unwrap();
    for preference in [
        TrimmingPreference::Cartesian,
        TrimmingPreference::Unspecified,
    ] {
        // Off the curve by ten micrometres, beside it and lifted off it.
        refused(
            trimmed_path(
                GeometryNode::Curve2(plan.clone()),
                TrimSelector::Point2(p + axiolid_core::Vec2::new(0.0, 1e-5)),
                TrimSelector::Point2(q),
                preference,
            ),
            PATH_POINT_OFF_CURVE,
            "planar, beside",
        );
        refused(
            trimmed_path(
                GeometryNode::Curve3(space.clone()),
                TrimSelector::Point3(s),
                TrimSelector::Point3(r + axiolid_core::Vec3::new(0.0, 0.0, 1e-5)),
                preference,
            ),
            PATH_POINT_OFF_CURVE,
            "space, lifted",
        );
        // Through the loop's double point, its start and its end.
        let top = evaluate2(&plan_loop(), 0.5).unwrap();
        refused(
            trimmed_path(
                GeometryNode::Curve2(plan_loop()),
                TrimSelector::Point2(Point2::ZERO),
                TrimSelector::Point2(top),
                preference,
            ),
            PATH_POINT_AMBIGUOUS,
            "planar loop",
        );
        let top = evaluate3(&space_loop(), 0.5).unwrap();
        refused(
            trimmed_path(
                GeometryNode::Curve3(space_loop()),
                TrimSelector::Point3(top),
                TrimSelector::Point3(Point3::ZERO),
                preference,
            ),
            PATH_POINT_AMBIGUOUS,
            "space loop",
        );
    }
    // Within the joint tolerance the point is on the curve.
    assert!(trimmed_path(
        GeometryNode::Curve2(plan.clone()),
        TrimSelector::Point2(p + axiolid_core::Vec2::new(0.0, 1e-8)),
        TrimSelector::Point2(q),
        TrimmingPreference::Cartesian,
    )
    .is_ok());
}

/// A trim preferring its points falls back to a parameter selector beside
/// a point the B-spline refuses, as on every other curve; one preferring
/// its parameters never reads the point.
#[test]
fn a_refused_point_falls_back_to_a_parameter_beside_it() {
    let curve = plan_spline();
    let beside = evaluate2(&curve, 0.4).unwrap() + axiolid_core::Vec2::new(0.0, 1e-3);
    let end = evaluate2(&curve, 0.9).unwrap();
    for preference in [TrimmingPreference::Cartesian, TrimmingPreference::Parameter] {
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push(GeometryNode::Curve2(curve.clone())).unwrap();
        let trimmed = b
            .push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
                basis,
                start: vec![TrimSelector::Point2(beside), TrimSelector::Parameter(0.4)],
                end: vec![TrimSelector::Point2(end), TrimSelector::Parameter(0.9)],
                sense_agreement: true,
                preference,
            }))
            .unwrap();
        let graph = b.finish(vec![trimmed]).unwrap();
        let path = curve_path(&graph, trimmed).unwrap();
        let start = CompositeBasis::from_path(&path)
            .unwrap()
            .section_on(0.0, SeamSide::Outgoing)
            .unwrap()
            .point;
        let want = evaluate2(&curve, 0.4).unwrap();
        assert!(
            (Point2::new(start.x, start.y) - want).length() <= 1e-9,
            "{preference:?}: {start:?}"
        );
    }
}
