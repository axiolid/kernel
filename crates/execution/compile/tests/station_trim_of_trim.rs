//! A trim of a trim reads its parameters in the inner trim's basis's
//! parameter (#424, ADR 0082 amendment), as a directrix and a polyline do.
//!
//! A circle of radius 2 trimmed `[0, pi]`, then trimmed again to `[pi / 4,
//! pi / 2]`, is the arc of the CIRCLE at those angles: a station on it
//! stands where the polyline's ends stand and where the swept directrix runs,
//! not at those distances along the inner trim.

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{Circle3, Curve3};
use axiolid_measure::volume_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::resolve;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{MeshCompiler, PolylineCompiler};
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId,
    SolidOperation, Station, TrimSelector, TrimmingPreference,
};

const RADIUS: Scalar = 2.0;
const TUBE: Scalar = 0.05;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(1e-3)
        .expect("positive chord budget")
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn circle() -> Curve3 {
    Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: RADIUS,
    })
}

fn on(t: Scalar) -> Point3 {
    Point3::new(RADIUS * t.cos(), RADIUS * t.sin(), 0.0)
}

fn trimmed(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    (start, end): (TrimSelector, TrimSelector),
    sense: bool,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![start],
        end: vec![end],
        sense_agreement: sense,
        preference: TrimmingPreference::Parameter,
    }))
    .unwrap()
}

fn p(value: Scalar) -> TrimSelector {
    TrimSelector::Parameter(value)
}

/// The graph of `circle` trimmed `inner`, then `outer`, with a station at
/// each distance of `distances` along the result.
fn stations(
    inner: (Scalar, Scalar, bool),
    outer: (Scalar, Scalar, bool),
    distances: &[Scalar],
) -> (GeometryGraph, NodeId, Vec<NodeId>) {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(circle()).unwrap();
    let inner = trimmed(&mut b, basis, (p(inner.0), p(inner.1)), inner.2);
    let outer = trimmed(&mut b, inner, (p(outer.0), p(outer.1)), outer.2);
    let nodes: Vec<NodeId> = distances
        .iter()
        .map(|&d| {
            b.push(GeometryNode::CurveStation(CurveStation::new(
                outer,
                Station::at(d),
            )))
            .unwrap()
        })
        .collect();
    let mut roots = nodes.clone();
    roots.push(outer);
    (b.finish(roots).unwrap(), outer, nodes)
}

fn close(a: Point3, b: Point3, what: &str) {
    assert!((a - b).length() <= 1e-9, "{what}: {a:?} != {b:?}");
}

/// A station on the nested trim stands on the circle at the angle of its
/// distance, and the polyline's ends are the stations at `0` and `L`.
fn agrees(
    inner: (Scalar, Scalar, bool),
    outer: (Scalar, Scalar, bool),
    from: Scalar,
    sweep: Scalar,
) {
    let length = RADIUS * sweep.abs();
    let direction = sweep.signum();
    let distances = [0.0, 0.25 * length, 0.5 * length, length];
    let (graph, curve, nodes) = stations(inner, outer, &distances);
    for (node, distance) in nodes.iter().zip(distances) {
        let station = resolve(&graph, *node).unwrap();
        close(
            station.point,
            on(from + direction * distance / RADIUS),
            "station",
        );
    }
    let polyline = compiler()
        .compile_polyline(&graph, curve, &options())
        .unwrap();
    let first = resolve(&graph, nodes[0]).unwrap().point;
    let last = resolve(&graph, nodes[3]).unwrap().point;
    close(polyline.points[0], first, "polyline start");
    close(*polyline.points.last().unwrap(), last, "polyline end");
}

#[test]
fn a_station_on_a_trim_of_a_trim_matches_the_polyline() {
    // Forwards, then reversed by the outer sense.
    agrees(
        (0.0, PI, true),
        (FRAC_PI_4, FRAC_PI_2, true),
        FRAC_PI_4,
        FRAC_PI_4,
    );
    agrees(
        (0.0, PI, true),
        (FRAC_PI_4, FRAC_PI_2, false),
        FRAC_PI_2,
        -FRAC_PI_4,
    );
}

#[test]
fn a_station_on_a_trim_of_a_reversed_trim_runs_the_inner_direction() {
    // The inner trim is the angles [-pi, 0] run downwards; the outer
    // selectors 3 pi / 2 and 7 pi / 4 are the angles -pi / 2 and -pi / 4,
    // which the inner trim reaches from -pi / 4 down to -pi / 2.
    agrees(
        (0.0, PI, false),
        (1.5 * PI, 1.75 * PI, true),
        -FRAC_PI_4,
        -FRAC_PI_4,
    );
}

/// The volume of a disk swept along the circle trimmed by each of `trims`
/// in turn.
fn volume_along(trims: &[(Scalar, Scalar, bool)]) -> Scalar {
    let mut b = GeometryGraphBuilder::new();
    let mut directrix = b.push_value(circle()).unwrap();
    for &(start, end, sense) in trims {
        directrix = trimmed(&mut b, directrix, (p(start), p(end)), sense);
    }
    let sweep = b
        .push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix,
            radius: TUBE,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        }))
        .unwrap();
    let graph = b.finish(vec![sweep]).unwrap();
    let mesh = compiler()
        .compile_mesh(&graph, sweep, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .unwrap();
    volume_properties(&mesh, Tolerance::MILLIMETRE)
        .expect("a swept disk is closed")
        .signed_volume
}

#[test]
fn the_swept_directrix_runs_the_same_arc() {
    // The tube along the arc `[pi / 4, pi / 2]` of the circle is the same
    // solid whether it is trimmed there at once or inside `[0, pi]`; were
    // the outer selectors distances along the inner trim, it would be the
    // arc `[0.785, 1.571]` metres from its start instead.
    let direct = volume_along(&[(FRAC_PI_4, FRAC_PI_2, true)]);
    let nested = volume_along(&[(0.0, PI, true), (FRAC_PI_4, FRAC_PI_2, true)]);
    assert!(
        (nested - direct).abs() / direct < 1e-9,
        "{nested} is not the arc {direct}"
    );
}

#[test]
fn a_point_or_arc_length_selector_on_a_trim_of_a_trim_is_refused() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(circle()).unwrap();
    let inner = trimmed(&mut b, basis, (p(0.0), p(PI)), true);
    let outer = trimmed(
        &mut b,
        inner,
        (TrimSelector::ArcLength(0.5), TrimSelector::ArcLength(1.5)),
        true,
    );
    let station = b
        .push(GeometryNode::CurveStation(CurveStation::new(
            outer,
            Station::at(0.1),
        )))
        .unwrap();
    let graph = b.finish(vec![station]).unwrap();
    let result: GeomResult<_> = resolve(&graph, station);
    assert!(
        matches!(
            &result,
            Err(GeomError::UnsupportedInput {
                operation: Operation::CurveEvaluation,
                input,
                ..
            }) if input.contains("a point or arc-length trim selector")
        ),
        "{result:?}"
    );
}
