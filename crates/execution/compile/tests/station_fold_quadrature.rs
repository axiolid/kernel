//! Folds of frames read by quadrature are decided exactly (#406, ADR 0082
//! amendment): an ellipse, a B-spline and an offset beside them bend at the
//! curvature of the EXACT curve in its true arc length, from a certified
//! second derivative (`FrameLaw::rates`). A section reaching the sharpest
//! bend's radius is refused as invalid input, one short of it is certified,
//! and one on the outer side of a bend is certified where the conservative
//! bound could not tell.
//!
//! The closed forms are an ellipse's curvature `a b / (a^2 sin^2 t + b^2 cos^2
//! t)^(3/2)`, `1/5` at the ends of the major axis of the 20 by 10 ellipse, the
//! offset's `kappa / (1 - d kappa)` (radius 7 at `d = -2`), and a B-spline's
//! sharpest bend by differences of the kernel's own tangents.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2};
use axiolid_curve::{BSplineCurve2, Curve2, Ellipse2, KnotSpec};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, RunStation,
    SectionAtStation, SolidOperation, Station, StationFrame, StationOffsets,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::derivative2;
use axiolid_reference::station::station_section2;

const BUDGET: Scalar = 0.01;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(BUDGET)
        .unwrap()
}

/// What the mesh compiler makes of a run.
#[derive(Debug)]
enum Outcome {
    Certified,
    Unbounded(&'static str),
    Refused(String),
}

fn mesh(graph: &GeometryGraph, root: NodeId) -> Outcome {
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    match compiler.compile_mesh_with_deviation(graph, root, &options()) {
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

/// A 1 m by 1 m rectangle run along `basis` between the distances `at`,
/// `lateral` to the left.
fn run(
    basis: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId,
    lateral: Scalar,
    at: [Scalar; 2],
) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let basis = basis(&mut b);
    let profile = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 1.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let offsets = StationOffsets::new(lateral, 0.0, 0.0);
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: at
                    .iter()
                    .map(|&d| SectionAtStation::new(profile, Station::new(d, offsets)))
                    .collect(),
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

/// The ellipse of semi-axes 20 and 10 about the origin.
fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 20.0,
        semi_axis_y: 10.0,
    })
}

#[test]
fn an_ellipse_folds_where_its_sharpest_bend_meets_the_section() {
    // Its sharpest bend is `10^2 / 20 = 5` m at the ends of the major axis,
    // where the run starts: `F = 1 - (a + 1/2) / 5` there.
    let case = |lateral| run(|b| b.push_value(ellipse()).unwrap(), lateral, [0.0, 10.0]);
    let (graph, root) = case(4.6);
    refused(&mesh(&graph, root));
    let (graph, root) = case(4.3);
    certified(&mesh(&graph, root));
    // On the outer side no bend folds a section however far it stands: a
    // magnitude bound reads `1 - 20 * 0.2` and cannot tell.
    let (graph, root) = case(-20.0);
    certified(&mesh(&graph, root));
}

#[test]
fn an_offset_beside_an_ellipse_folds_at_its_own_radius() {
    // The offset by `-2` (to the right, outward) of the ellipse has the
    // radius `5 + 2 = 7` at the major axis: `F = 1 - (a + 1/2) / 7`.
    let beside = |b: &mut GeometryGraphBuilder| {
        let basis = b.push_value(ellipse()).unwrap();
        b.push(GeometryNode::CurveRelation(CurveRelation::Offset {
            basis,
            distance: -2.0,
            reference_direction: None,
        }))
        .unwrap()
    };
    let (graph, root) = run(beside, 6.6, [0.0, 10.0]);
    refused(&mesh(&graph, root));
    let (graph, root) = run(beside, 6.3, [0.0, 10.0]);
    certified(&mesh(&graph, root));
}

/// A cubic plan B-spline with an interior knot, about 25 long.
fn cubic() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(6.0, 1.0),
            Point2::new(12.0, -2.0),
            Point2::new(18.0, 3.0),
            Point2::new(24.0, 2.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// A quadratic plan B-spline whose curvature jumps at its interior knot.
fn quadratic() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(12.0, 0.0),
            Point2::new(12.0, 12.0),
            Point2::new(24.0, 14.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 1, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// The signed curvature of `curve` of largest magnitude over the distances
/// `[lo, hi]`, by differences of the kernel's own tangents.
fn sharpest(curve: &Curve2, (lo, hi): (Scalar, Scalar)) -> Scalar {
    let h = 1e-3;
    (0..=2000)
        .map(|i| lo + (hi - lo) * Scalar::from(i) / 2000.0)
        .map(|s| {
            let (before, here, after) = (
                station_section2(curve, s - h).unwrap(),
                station_section2(curve, s).unwrap(),
                station_section2(curve, s + h).unwrap(),
            );
            here.lateral
                .dot((after.tangent - before.tangent) / (2.0 * h))
        })
        .fold(
            0.0,
            |worst: Scalar, k| if k.abs() > worst.abs() { k } else { worst },
        )
}

/// The lateral offset at which a unit square reaches `scale` times the
/// radius of the sharpest bend, on its inner side.
fn reaching(kappa: Scalar, scale: Scalar) -> Scalar {
    kappa.signum() * (scale / kappa.abs() - 0.5)
}

#[test]
fn a_b_spline_folds_at_its_sampled_sharpest_bend() {
    let at = [1.0, 22.0];
    let kappa = sharpest(&cubic(), (at[0], at[1]));
    assert!(kappa.abs() > 0.05, "{kappa}");
    // A section within 5% of the sharpest bend's radius.
    let (graph, root) = run(
        |b| b.push_value(cubic()).unwrap(),
        reaching(kappa, 1.05),
        at,
    );
    refused(&mesh(&graph, root));
    let (graph, root) = run(
        |b| b.push_value(cubic()).unwrap(),
        reaching(kappa, 0.95),
        at,
    );
    certified(&mesh(&graph, root));
}

/// The arc length of `curve` from its start to the parameter `t`, by
/// Simpson's rule.
fn arc_to(curve: &Curve2, t: Scalar) -> Scalar {
    let n = 2000;
    let speed = |u: Scalar| derivative2(curve, u).unwrap().length();
    let h = t / n as Scalar;
    (0..n)
        .map(|i| {
            let a = h * i as Scalar;
            h / 6.0 * (speed(a) + 4.0 * speed(a + 0.5 * h) + speed(a + h))
        })
        .sum()
}

#[test]
fn a_degree_two_b_spline_folds_away_from_its_knot_and_is_never_refused_at_it() {
    let curve = quadratic();
    let knot = arc_to(&curve, 0.5);
    assert!(knot > 10.0, "{knot}");
    // A fold well before the knot is found.
    let before = [1.0, knot - 4.0];
    let kappa = sharpest(&curve, (before[0], before[1]));
    let (graph, root) = run(
        |b| b.push_value(quadratic()).unwrap(),
        reaching(kappa, 1.05),
        before,
    );
    refused(&mesh(&graph, root));
    let (graph, root) = run(
        |b| b.push_value(quadratic()).unwrap(),
        reaching(kappa, 0.95),
        before,
    );
    certified(&mesh(&graph, root));
    // A run across the knot, where the curvature jumps: the spans across it
    // are checked by the conservative bound, which clears a section close
    // to the curve and never refuses.
    let (graph, root) = run(
        |b| b.push_value(quadratic()).unwrap(),
        0.5,
        [knot - 3.0, knot + 3.0],
    );
    certified(&mesh(&graph, root));
}

/// The deck of openbimrs/ifc#443: the ellipse of semi-axes 10 and 6 about
/// `(30, 0)` offset by `0.5` to station 2, widening to `0.75` at station 20.
fn deck(b: &mut GeometryGraphBuilder) -> NodeId {
    let basis = b
        .push_value(Curve2::Ellipse(Ellipse2 {
            frame: Frame2 {
                origin: Point2::new(30.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            semi_axis_x: 10.0,
            semi_axis_y: 6.0,
        }))
        .unwrap();
    let lateral = |d| StationOffsets::new(d, 0.0, 0.0);
    b.push(GeometryNode::CurveRelation(
        CurveRelation::OffsetsAtStations {
            basis,
            stations: vec![
                RunStation::at_basis_start(lateral(0.5)),
                RunStation::new(Station::new(2.0, lateral(0.5))),
                RunStation::new(Station::new(20.0, lateral(0.75))),
                RunStation::at_basis_end(lateral(0.75)),
            ],
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

#[test]
fn the_deck_offset_is_decided_exactly() {
    // Its sharpest bend, at the major axis, is `36 / 10 = 3.6` m less the
    // offset to the inside, about 3.1 m: a section reaching it folds, one
    // clear of every bend is certified. The runs do not cross the mitred
    // corner at station 2.
    let (graph, root) = run(deck, 0.0, [3.0, 12.0]);
    certified(&mesh(&graph, root));
    let (graph, root) = run(deck, 1.0, [3.0, 12.0]);
    certified(&mesh(&graph, root));
    let (graph, root) = run(deck, 3.0, [0.2, 1.5]);
    refused(&mesh(&graph, root));
    let (graph, root) = run(deck, 2.0, [0.2, 1.5]);
    certified(&mesh(&graph, root));
}

/// The ellipse of semi-axes 10 and 6 about `(30, 0)`, offset by 1.5 to its
/// left (inward), as `IfcOffsetCurve2D` states it (axioval/engine).
fn offset_ellipse(b: &mut GeometryGraphBuilder) -> NodeId {
    let basis = b
        .push_value(Curve2::Ellipse(Ellipse2 {
            frame: Frame2 {
                origin: Point2::new(30.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            semi_axis_x: 10.0,
            semi_axis_y: 6.0,
        }))
        .unwrap();
    b.push(GeometryNode::CurveRelation(CurveRelation::Offset {
        basis,
        distance: 1.5,
        reference_direction: None,
    }))
    .unwrap()
}

/// A rectangle `1` m high, `first` wide at station 0 and `last` wide at
/// station 5, along the inward offset of the ellipse.
fn widening(first: Scalar, last: Scalar) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let basis = offset_ellipse(&mut b);
    let rectangle = |b: &mut GeometryGraphBuilder, x| {
        b.push_value(Profile::Rectangle(RectangleProfile {
            x,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap()
    };
    let (narrow, wide) = (rectangle(&mut b, first), rectangle(&mut b, last));
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: vec![
                    SectionAtStation::new(narrow, Station::at(0.0)),
                    SectionAtStation::new(wide, Station::at(5.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

#[test]
fn a_widening_section_along_an_inward_offset_of_an_ellipse_is_decided() {
    // The offset's sharpest bend, at the major axis, has the radius
    // `3.6 - 1.5 = 2.1`; a rectangle widening from 2 to 4 over stations 0
    // to 5 reaches `a = 2` on the inner side at station 5, where the bend
    // has eased: no fold, so the run is certified rather than unbounded.
    let (graph, root) = widening(2.0, 4.0);
    certified(&mesh(&graph, root));
}

#[test]
fn a_section_reaching_the_inward_offsets_sharpest_bend_is_refused() {
    // The same run starting 4.4 wide: its inner edge, at `a = 2.2`, is past
    // the radius 2.1 of the bend at station 0.
    let (graph, root) = widening(4.4, 4.0);
    refused(&mesh(&graph, root));
}
