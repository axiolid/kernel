//! Stations along curve relations (#285): a composite of trimmed pieces, a
//! composite of segments placed at stations of another composite, nested
//! and trimmed composites, and the refusals. Expected frames are closed
//! forms, stations on the individual pieces at the offset distance, and
//! placements composed by hand, never the composite code under test.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_mesh_compile::station::{placement, resolve, seams, ResolvedStation};
use axiolid_model::{
    CurveRelation, CurveSegment, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode,
    InstanceAtStation, NodeId, OrientedCurveStation, SeamSide, Station, StationOffsets,
    StationOrientation, Transition, TrimSelector, TrimmingPreference,
};

const EPS: Scalar = 1e-9;
const PI: Scalar = core::f64::consts::PI;

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?} (off by {:e})",
        (actual - expected).abs().max_element()
    );
}

fn same(a: &ResolvedStation, b: &ResolvedStation, eps: Scalar, what: &str) {
    close3(a.point, b.point, eps, &format!("{what}: point"));
    close3(a.frame.x, b.frame.x, eps, &format!("{what}: x"));
    close3(a.frame.y, b.frame.y, eps, &format!("{what}: y"));
    close3(a.frame.z, b.frame.z, eps, &format!("{what}: z"));
}

fn refused_naming<T: core::fmt::Debug>(result: GeomResult<T>, needle: &str) {
    match result {
        Err(GeomError::InvalidInput(detail)) => {
            assert!(detail.contains(needle), "expected {needle:?} in {detail:?}");
        }
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert!(input.contains(needle), "expected {needle:?} in {input:?}");
        }
        other => panic!("expected a refusal naming {needle:?}, got {other:?}"),
    }
}

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

fn station(b: &mut GeometryGraphBuilder, basis: NodeId, d: Scalar, side: SeamSide) -> NodeId {
    b.push_value(CurveStation::new(basis, Station::at(d)).with_seam_side(side))
        .unwrap()
}

fn at(basis: NodeId, d: Scalar, side: SeamSide) -> OrientedCurveStation {
    OrientedCurveStation::from(CurveStation::new(basis, Station::at(d))).with_seam_side(side)
}

/// The base: a 2D line along `+x` and a circle of radius 5 about
/// `(15, 0)`, as atomic nodes.
struct Base {
    line: NodeId,
    circle: NodeId,
    /// The line over `[0, 10]`, then the circle from angle `pi` to
    /// `3 pi / 2`: a right-angle turn at distance 10.
    composite: NodeId,
}

fn base(b: &mut GeometryGraphBuilder) -> Base {
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let circle = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(15.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 5.0,
        }))
        .unwrap();
    let straight = trimmed(b, line, 0.0, 10.0, true);
    let arc = trimmed(b, circle, PI, 1.5 * PI, true);
    let composite = composite(b, &[(straight, true), (arc, true)]);
    Base {
        line,
        circle,
        composite,
    }
}

/// The closed-form section on the base's arc at angle `a`: point, tangent,
/// left lateral (up is `+Z`).
fn on_arc(a: Scalar) -> (Point3, Vec3, Vec3) {
    let (s, c) = a.sin_cos();
    let t = Vec3::new(-s, c, 0.0);
    (
        Point3::new(15.0 + 5.0 * c, 5.0 * s, 0.0),
        t,
        Vec3::new(-t.y, t.x, 0.0),
    )
}

#[test]
fn a_line_and_arc_composite_agrees_with_its_pieces_at_the_offset_distance() {
    let mut b = GeometryGraphBuilder::new();
    let base = base(&mut b);
    // Each composite station with its counterpart on the piece itself.
    let mut pairs = Vec::new();
    for d in [0.0, 3.0, 9.0] {
        pairs.push((
            station(&mut b, base.composite, d, SeamSide::Outgoing),
            station(&mut b, base.line, d, SeamSide::Outgoing),
        ));
    }
    for into in [0.5, 4.0, 2.5 * PI] {
        pairs.push((
            station(&mut b, base.composite, 10.0 + into, SeamSide::Outgoing),
            station(&mut b, base.circle, 5.0 * PI + into, SeamSide::Outgoing),
        ));
    }
    // The joint, both sides: the line's end and the arc's start.
    pairs.push((
        station(&mut b, base.composite, 10.0, SeamSide::Incoming),
        station(&mut b, base.line, 10.0, SeamSide::Outgoing),
    ));
    pairs.push((
        station(&mut b, base.composite, 10.0, SeamSide::Outgoing),
        station(&mut b, base.circle, 5.0 * PI, SeamSide::Outgoing),
    ));
    // Offsets and an orientation are read in the joint's chosen frame.
    let offset = b
        .push_value(
            OrientedCurveStation::new(
                CurveStation::new(
                    base.composite,
                    Station::new(10.0, StationOffsets::new(1.0, 2.0, 0.5)),
                ),
                StationOrientation::new(None, Some(Vec3::new(0.0, 1.0, 0.0))),
            )
            .with_seam_side(SeamSide::Incoming),
        )
        .unwrap();
    let beyond = station(
        &mut b,
        base.composite,
        10.0 + 2.5 * PI + 1e-3,
        SeamSide::Outgoing,
    );
    let roots: Vec<NodeId> = pairs
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .chain([offset, beyond])
        .collect();
    let graph = b.finish(roots).unwrap();

    for (on_composite, on_piece) in &pairs {
        same(
            &resolve(&graph, *on_composite).unwrap(),
            &resolve(&graph, *on_piece).unwrap(),
            EPS,
            "composite against piece",
        );
    }
    // Independently, closed forms: inside the arc, and the joint's sides.
    let inside = resolve(&graph, pairs[4].0).unwrap();
    let (p, t, l) = on_arc(PI + 4.0 / 5.0);
    close3(inside.point, p, EPS, "arc point");
    close3(inside.frame.x, t, EPS, "arc tangent");
    close3(inside.frame.z, -l, EPS, "arc right");
    let incoming = resolve(&graph, pairs[6].0).unwrap();
    let outgoing = resolve(&graph, pairs[7].0).unwrap();
    close3(incoming.frame.x, Vec3::X, EPS, "incoming heads +x");
    close3(outgoing.frame.x, -Vec3::Y, EPS, "outgoing heads -y");
    close3(incoming.point, Point3::new(10.0, 0.0, 0.0), EPS, "joint");
    close3(outgoing.point, Point3::new(10.0, 0.0, 0.0), EPS, "joint");

    // In the incoming frame (t = +x, l = +y): (10.5, 1, 2), turned so its
    // tangent is the lateral +y.
    let offset = resolve(&graph, offset).unwrap();
    close3(
        offset.point,
        Point3::new(10.5, 1.0, 2.0),
        EPS,
        "offset point",
    );
    close3(offset.frame.x, Vec3::Y, EPS, "turned tangent");
    refused_naming(resolve(&graph, beyond), "beyond the curve's length");

    let joints = seams(
        &graph,
        base.composite,
        &axiolid_contracts::ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
    )
    .unwrap();
    assert_eq!(joints.len(), 1);
    assert!((joints[0].distance - 10.0).abs() <= 1e-12);
    assert!(joints[0].exact && !joints[0].smooth);
}

/// Segments placed at stations of the base composite: a local line over
/// `[0, 6]` at station 4, and a local quarter arc (radius 5 about `(0,
/// 5)`, from `-pi / 2` to `0`) at the joint, outgoing. Together they
/// retrace the base from distance 4.
struct Placed {
    straight: NodeId,
    bend: NodeId,
    segmented: NodeId,
}

fn placed(b: &mut GeometryGraphBuilder, base: &Base) -> Placed {
    let local_line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let local_circle = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.0, 5.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 5.0,
        }))
        .unwrap();
    let line_span = trimmed(b, local_line, 0.0, 6.0, true);
    let arc_span = trimmed(b, local_circle, -0.5 * PI, 0.0, true);
    let straight = b
        .push_value(InstanceAtStation::new(
            line_span,
            at(base.composite, 4.0, SeamSide::Outgoing),
        ))
        .unwrap();
    let bend = b
        .push_value(InstanceAtStation::new(
            arc_span,
            at(base.composite, 10.0, SeamSide::Outgoing),
        ))
        .unwrap();
    let segmented = composite(b, &[(straight, true), (bend, true)]);
    Placed {
        straight,
        bend,
        segmented,
    }
}

#[test]
fn segments_placed_on_a_composite_compose_by_hand() {
    let mut b = GeometryGraphBuilder::new();
    let base = base(&mut b);
    let placed = placed(&mut b, &base);
    let mut along = Vec::new();
    for d in [1.0, 6.0, 7.5, 6.0 + 2.5 * PI] {
        along.push((
            station(&mut b, placed.segmented, d, SeamSide::Outgoing),
            station(&mut b, base.composite, d + 4.0, SeamSide::Outgoing),
        ));
    }
    let incoming = station(&mut b, placed.segmented, 6.0, SeamSide::Incoming);
    let on_bend = station(&mut b, placed.bend, 1.5, SeamSide::Outgoing);
    let roots: Vec<NodeId> = along
        .iter()
        .flat_map(|&(a, b)| [a, b])
        .chain([incoming, on_bend, placed.straight])
        .collect();
    let graph = b.finish(roots).unwrap();

    // The placements, by hand: at station 4 the base frame is t = +x,
    // l = +y at (4, 0, 0); at the joint, outgoing, t = -y, l = +x at
    // (10, 0, 0).
    let first = placement(&graph, placed.straight).unwrap();
    close3(
        first.transform.translation,
        Point3::new(4.0, 0.0, 0.0),
        EPS,
        "first origin",
    );
    close3(
        first.transform.transform_vector3(Vec3::X),
        Vec3::X,
        EPS,
        "first x",
    );
    assert!(first.exact, "a placement on the base's line is exact");
    let second = placement(&graph, placed.bend).unwrap();
    close3(
        second.transform.translation,
        Point3::new(10.0, 0.0, 0.0),
        EPS,
        "second origin",
    );
    close3(
        second.transform.transform_vector3(Vec3::X),
        -Vec3::Y,
        EPS,
        "second x",
    );
    close3(
        second.transform.transform_vector3(Vec3::Y),
        Vec3::X,
        EPS,
        "second y",
    );
    assert!(!second.exact, "a placement on the base's arc is not exact");

    // A station along the segmented composite inside the bend: the local
    // arc's station at 5 phi from its start, moved by the second placement.
    let phi = 1.5 / 5.0;
    let theta = -0.5 * PI + phi;
    let (s, c) = theta.sin_cos();
    let local_point = Point3::new(5.0 * c, 5.0 + 5.0 * s, 0.0);
    let local_tangent = Vec3::new(-s, c, 0.0);
    let point = second.transform.transform_point3(local_point);
    let tangent = second.transform.transform_vector3(local_tangent);
    let on_bend = resolve(&graph, on_bend).unwrap();
    close3(on_bend.point, point, EPS, "bend point by hand");
    close3(on_bend.frame.x, tangent, EPS, "bend tangent by hand");
    close3(on_bend.frame.y, Vec3::Z, EPS, "bend up");
    let on_segmented = resolve(&graph, along[2].0).unwrap();
    close3(on_segmented.point, point, EPS, "segmented point by hand");
    close3(
        on_segmented.frame.x,
        tangent,
        EPS,
        "segmented tangent by hand",
    );

    // And the segmented composite retraces the base from distance 4.
    for (on_segmented, on_base) in &along {
        same(
            &resolve(&graph, *on_segmented).unwrap(),
            &resolve(&graph, *on_base).unwrap(),
            EPS,
            "segmented against base",
        );
    }
    // Its joint at 6 reads the straight segment incoming.
    close3(
        resolve(&graph, incoming).unwrap().frame.x,
        Vec3::X,
        EPS,
        "incoming",
    );
    let joints = seams(
        &graph,
        placed.segmented,
        &axiolid_contracts::ExecutionOptions::new(axiolid_core::Tolerance::MILLIMETRE),
    )
    .unwrap();
    assert_eq!(joints.len(), 1);
    assert!((joints[0].distance - 6.0).abs() <= 1e-12);
}

#[test]
fn a_node_placed_along_a_segmented_composite_is_exact_only_on_exact_pieces() {
    let mut b = GeometryGraphBuilder::new();
    let base = base(&mut b);
    let placed = placed(&mut b, &base);
    let marker = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::Y,
        }))
        .unwrap();
    let on_line = b
        .push_value(InstanceAtStation::new(
            marker,
            at(placed.segmented, 3.0, SeamSide::Outgoing),
        ))
        .unwrap();
    let on_joint = b
        .push_value(InstanceAtStation::new(
            marker,
            at(placed.segmented, 6.0, SeamSide::Incoming),
        ))
        .unwrap();
    let on_bend = b
        .push_value(InstanceAtStation::new(
            marker,
            at(placed.segmented, 7.0, SeamSide::Outgoing),
        ))
        .unwrap();
    // A line placed on the base's arc, trimmed or whole: its placement is
    // inexact, and so is every frame along it.
    let local = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let local_span = trimmed(&mut b, local, 0.0, 2.0, true);
    let mut on_inexact = Vec::new();
    for source in [local, local_span] {
        let leaning = b
            .push_value(InstanceAtStation::new(
                source,
                at(base.composite, 12.0, SeamSide::Outgoing),
            ))
            .unwrap();
        on_inexact.push(
            b.push_value(InstanceAtStation::new(
                marker,
                at(leaning, 0.5, SeamSide::Outgoing),
            ))
            .unwrap(),
        );
    }
    let graph = b
        .finish(
            [on_line, on_joint, on_bend]
                .into_iter()
                .chain(on_inexact.iter().copied())
                .collect(),
        )
        .unwrap();
    for id in on_inexact {
        let placed = placement(&graph, id).unwrap();
        assert!(!placed.exact, "a line placed on an arc is not exact");
        let (p, t, _) = on_arc(PI + 2.0 / 5.0);
        close3(
            placed.transform.translation,
            p + 0.5 * t,
            EPS,
            "on the placed line",
        );
    }
    // On the straight segment, itself placed exactly on the base's line.
    let exact = placement(&graph, on_line).unwrap();
    assert!(exact.exact);
    close3(
        exact.transform.translation,
        Point3::new(7.0, 0.0, 0.0),
        EPS,
        "on line",
    );
    assert!(placement(&graph, on_joint).unwrap().exact);
    // On the bend: an arc, placed on the base's arc.
    let bent = placement(&graph, on_bend).unwrap();
    assert!(!bent.exact);
    let (p, t, _) = on_arc(PI + 1.0 / 5.0);
    close3(bent.transform.translation, p, EPS, "on bend");
    close3(bent.transform.transform_vector3(Vec3::X), t, EPS, "bend x");
}

#[test]
fn nested_trimmed_and_reversed_relations_measure_along_their_pieces() {
    let mut b = GeometryGraphBuilder::new();
    let base = base(&mut b);
    // The base, then 5 on along +x from its end (15, -5).
    let tail_line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::new(15.0, -5.0),
            direction: Vec2::new(2.0, 0.0),
        }))
        .unwrap();
    let tail = trimmed(&mut b, tail_line, 0.0, 2.5, true);
    let nested = composite(&mut b, &[(base.composite, true), (tail, true)]);
    let length = 10.0 + 2.5 * PI;
    let on_tail = station(&mut b, nested, length + 4.5, SeamSide::Outgoing);
    // The base between 2 and 12, both ways.
    let window = trimmed(&mut b, base.composite, 2.0, 12.0, true);
    let backwards = trimmed(&mut b, base.composite, 2.0, 12.0, false);
    let window_start = station(&mut b, window, 0.0, SeamSide::Outgoing);
    let base_at_2 = station(&mut b, base.composite, 2.0, SeamSide::Outgoing);
    let back_start = station(&mut b, backwards, 0.0, SeamSide::Outgoing);
    let base_at_12 = station(&mut b, base.composite, 12.0, SeamSide::Outgoing);
    // A whole composite reversed inside another: the arc then the line,
    // backwards; its joint at 2.5 pi heads -x on the incoming side.
    let reversed = composite(&mut b, &[(base.composite, false)]);
    let reversed_joint = station(&mut b, reversed, 2.5 * PI, SeamSide::Outgoing);
    let graph = b
        .finish(vec![
            on_tail,
            window_start,
            base_at_2,
            back_start,
            base_at_12,
            reversed_joint,
        ])
        .unwrap();
    let tail_station = resolve(&graph, on_tail).unwrap();
    close3(
        tail_station.point,
        Point3::new(19.5, -5.0, 0.0),
        EPS,
        "nested tail",
    );
    same(
        &resolve(&graph, window_start).unwrap(),
        &resolve(&graph, base_at_2).unwrap(),
        EPS,
        "trimmed composite",
    );
    let back = resolve(&graph, back_start).unwrap();
    let forward = resolve(&graph, base_at_12).unwrap();
    close3(back.point, forward.point, EPS, "reversed trim point");
    close3(back.frame.x, -forward.frame.x, EPS, "reversed trim tangent");
    close3(back.frame.z, -forward.frame.z, EPS, "reversed trim right");
    close3(back.frame.y, forward.frame.y, EPS, "reversed trim up");
    let joint = resolve(&graph, reversed_joint).unwrap();
    close3(
        joint.point,
        Point3::new(10.0, 0.0, 0.0),
        EPS,
        "reversed joint",
    );
    close3(joint.frame.x, -Vec3::X, EPS, "outgoing: the line backwards");
}

#[test]
fn a_gap_or_an_undeclared_reversed_piece_is_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let base = base(&mut b);
    let short = trimmed(&mut b, base.line, 0.0, 9.5, true);
    let arc = trimmed(&mut b, base.circle, PI, 1.5 * PI, true);
    let gapped = composite(&mut b, &[(short, true), (arc, true)]);
    // The arc from (15, 5) back to (10, 0): its end meets the line.
    let full = trimmed(&mut b, base.line, 0.0, 10.0, true);
    let back_arc = trimmed(&mut b, base.circle, 0.5 * PI, PI, true);
    let undeclared = composite(&mut b, &[(full, true), (back_arc, true)]);
    let declared = composite(&mut b, &[(full, true), (back_arc, false)]);
    let on_gapped = station(&mut b, gapped, 1.0, SeamSide::Outgoing);
    let on_undeclared = station(&mut b, undeclared, 1.0, SeamSide::Outgoing);
    let on_declared = station(&mut b, declared, 10.0, SeamSide::Outgoing);
    // An instanced curve is no station basis (an offset is, since #289).
    let offset = b
        .push(GeometryNode::Instance(axiolid_model::Instance {
            source: base.line,
            transform: axiolid_core::Transform3::IDENTITY,
        }))
        .unwrap();
    let on_offset = station(&mut b, offset, 1.0, SeamSide::Outgoing);
    let graph: GeometryGraph = b
        .finish(vec![on_gapped, on_undeclared, on_declared, on_offset])
        .unwrap();
    refused_naming(resolve(&graph, on_gapped), "a gap of 0.5");
    refused_naming(resolve(&graph, on_undeclared), "a reversed piece");
    // Declared reversed, it turns left at the joint, heading +y.
    close3(
        resolve(&graph, on_declared).unwrap().frame.x,
        Vec3::Y,
        EPS,
        "declared reversed",
    );
    refused_naming(resolve(&graph, on_offset), "other than a composite, a trim");
}
