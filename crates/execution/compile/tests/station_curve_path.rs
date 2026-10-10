//! A graph curve relation flattened into the neutral curve path (#290):
//! `station::curve_path` hands out the pieces a station on the relation is
//! measured along, so a consumer reading them through the curve-evaluation
//! contract gets what `station::resolve` gets. Checked bitwise against
//! resolved stations on a line and arc composite, a trim of it, and a
//! composite of segments placed at its stations, on both sides of each
//! joint.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, CurvePath, Line2, PathCurve, SeamSide};
use axiolid_mesh_compile::station::{curve_path, resolve};
use axiolid_model::{
    CurveRelation, CurveSegment, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode,
    InstanceAtStation, NodeId, OrientedCurveStation, Station, Transition, TrimSelector,
    TrimmingPreference,
};
use axiolid_reference::station::CompositeBasis;

const PI: Scalar = core::f64::consts::PI;
const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];

fn trimmed(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    start: Scalar,
    end: Scalar,
    sense: bool,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![TrimSelector::Parameter(start)],
        end: vec![TrimSelector::Parameter(end)],
        sense_agreement: sense,
        preference: TrimmingPreference::Parameter,
    }))
    .unwrap()
}

fn composite(b: &mut GeometryGraphBuilder, pieces: &[(NodeId, bool)]) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments: pieces
            .iter()
            .map(|&(curve, same_sense)| CurveSegment {
                curve,
                same_sense,
                transition: Transition::Continuous,
            })
            .collect(),
    }))
    .unwrap()
}

fn at(basis: NodeId, d: Scalar, side: SeamSide) -> OrientedCurveStation {
    OrientedCurveStation::from(CurveStation::new(basis, Station::at(d))).with_seam_side(side)
}

fn line2(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    }))
    .unwrap()
}

fn circle2(b: &mut GeometryGraphBuilder, centre: Point2) -> NodeId {
    b.push_value(Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: centre,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 5.0,
    }))
    .unwrap()
}

/// The relations under test, with the distances of their joints.
struct Bases {
    line: NodeId,
    /// A line over `[0, 10]`, then a quarter arc of radius 5 about `(15,
    /// 0)` from angle `pi`: a right-angle corner at 10.
    composite: NodeId,
    /// The composite trimmed to `[3, 12]` against its sense: the arc back
    /// to the joint at 2, then the line back to 3.
    trim: NodeId,
    /// A local line over `[0, 6]` placed at the composite's station 4 and
    /// a local quarter arc placed at its joint (outgoing): the composite
    /// again from distance 4, with its joint at 6.
    segmented: NodeId,
}

fn bases(b: &mut GeometryGraphBuilder) -> Bases {
    let b = &mut *b;
    let line = line2(b);
    let circle = circle2(b, Point2::new(15.0, 0.0));
    let straight = trimmed(b, line, 0.0, 10.0, true);
    let arc = trimmed(b, circle, PI, 1.5 * PI, true);
    let composite = composite(b, &[(straight, true), (arc, true)]);
    let trim = trimmed(b, composite, 3.0, 12.0, false);
    let local_line = line2(b);
    let local_circle = circle2(b, Point2::new(0.0, 5.0));
    let line_span = trimmed(b, local_line, 0.0, 6.0, true);
    let arc_span = trimmed(b, local_circle, -0.5 * PI, 0.0, true);
    let placed_line = b
        .push_value(InstanceAtStation::new(
            line_span,
            at(composite, 4.0, SeamSide::Outgoing),
        ))
        .unwrap();
    let placed_arc = b
        .push_value(InstanceAtStation::new(
            arc_span,
            at(composite, 10.0, SeamSide::Outgoing),
        ))
        .unwrap();
    let segmented = crate::composite(b, &[(placed_line, true), (placed_arc, true)]);
    Bases {
        line,
        composite,
        trim,
        segmented,
    }
}

/// The bases alone, as a graph.
fn graph_of_bases() -> (GeometryGraph, Bases) {
    let mut b = GeometryGraphBuilder::new();
    let bases = bases(&mut b);
    let roots = vec![bases.trim, bases.segmented, bases.line];
    (b.finish(roots).unwrap(), bases)
}

/// Stations on `basis` at `distances`, both sides, resolved through the
/// graph; and the same pieces read from the flattened path.
fn path_reads_as_resolved(which: fn(&Bases) -> NodeId, distances: &[Scalar], what: &str) {
    let mut b = GeometryGraphBuilder::new();
    let basis = which(&bases(&mut b));
    let mut stations = Vec::new();
    for &d in distances {
        for side in SIDES {
            stations.push((d, side, b.push_value(at(basis, d, side)).unwrap()));
        }
    }
    let roots: Vec<NodeId> = stations.iter().map(|&(_, _, id)| id).collect();
    let graph = b.finish(roots).unwrap();
    let path = curve_path(&graph, basis).unwrap();
    let composite = CompositeBasis::from_path(&path).unwrap();
    for (d, side, id) in stations {
        let resolved = resolve(&graph, id).unwrap();
        assert_eq!(
            resolved.section,
            composite.section_on(d, side).unwrap(),
            "{what} {side:?} at {d}"
        );
    }
}

fn around(d: Scalar) -> [Scalar; 3] {
    let hair = 0.4e-12 * d.max(1.0);
    [d - hair, d, d + hair]
}

fn tangent_at(path: &CurvePath, d: Scalar, side: SeamSide) -> Vec3 {
    CompositeBasis::from_path(path)
        .unwrap()
        .section_on(d, side)
        .unwrap()
        .tangent
}

fn close(actual: Vec3, expected: Vec3, what: &str) {
    let off = (actual - expected).length();
    assert!(off <= 1e-9, "{what}: {actual:?} != {expected:?}");
}

#[test]
fn a_flattened_relation_reads_as_its_resolved_stations() {
    let (graph, bases) = graph_of_bases();
    let mut distances = vec![0.0, 4.0, 12.0];
    distances.extend(around(10.0));
    path_reads_as_resolved(|bases| bases.composite, &distances, "composite");
    let mut distances = vec![0.0, 1.0, 5.0, 9.0];
    distances.extend(around(2.0));
    path_reads_as_resolved(|bases| bases.trim, &distances, "trim");
    let mut distances = vec![1.0, 3.0, 7.5];
    distances.extend(around(6.0));
    path_reads_as_resolved(|bases| bases.segmented, &distances, "segmented");

    // Each joint's sides, independently.
    let composite = curve_path(&graph, bases.composite).unwrap();
    close(
        tangent_at(&composite, 10.0, SeamSide::Incoming),
        Vec3::X,
        "composite in",
    );
    close(
        tangent_at(&composite, 10.0, SeamSide::Outgoing),
        -Vec3::Y,
        "composite out",
    );
    let trim = curve_path(&graph, bases.trim).unwrap();
    close(
        tangent_at(&trim, 2.0, SeamSide::Incoming),
        Vec3::Y,
        "trim in",
    );
    close(
        tangent_at(&trim, 2.0, SeamSide::Outgoing),
        -Vec3::X,
        "trim out",
    );
    let segmented = curve_path(&graph, bases.segmented).unwrap();
    close(
        tangent_at(&segmented, 6.0, SeamSide::Incoming),
        Vec3::X,
        "segments in",
    );
    close(
        tangent_at(&segmented, 6.0, SeamSide::Outgoing),
        -Vec3::Y,
        "segments out",
    );

    // The placements are carried by the pieces, exact only on the line.
    let pieces = segmented.pieces();
    assert_eq!(pieces.len(), 2);
    assert!(pieces.iter().all(|piece| piece.placement.is_some()));
    assert!(pieces[0].frame_is_exact() && !pieces[1].frame_is_exact());
    let start = CompositeBasis::from_path(&segmented)
        .unwrap()
        .section_on(0.0, SeamSide::Outgoing)
        .unwrap();
    close(
        start.point - Point3::ZERO,
        Vec3::new(4.0, 0.0, 0.0),
        "segments start at 4",
    );
}

#[test]
fn an_atomic_curve_is_one_whole_piece_and_other_relations_are_refused() {
    let (graph, bases) = graph_of_bases();
    let path = curve_path(&graph, bases.line).unwrap();
    assert_eq!(path.pieces().len(), 1);
    let piece = &path.pieces()[0];
    assert!(matches!(piece.curve, PathCurve::Two(Curve2::Line(_))));
    assert_eq!((piece.start, piece.end, piece.reversed), (0.0, 1.0, false));

    let mut b = GeometryGraphBuilder::new();
    let line = line2(&mut b);
    // An instanced curve is no station basis (an offset is, since #289).
    let offset = b
        .push(GeometryNode::Instance(axiolid_model::Instance {
            source: line,
            transform: axiolid_core::Transform3::IDENTITY,
        }))
        .unwrap();
    let graph = b.finish(vec![offset]).unwrap();
    refused(curve_path(&graph, offset), "a curve relation other than");
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, needle: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert!(input.contains(needle), "expected {needle:?} in {input:?}");
        }
        other => panic!("expected a refusal naming {needle:?}, got {other:?}"),
    }
}
