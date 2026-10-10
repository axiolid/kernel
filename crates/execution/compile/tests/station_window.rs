//! Seam-snapping windows on oriented stations (#294): a station within its
//! window of one seam of its basis is read on that seam from its side --
//! on a composite's joint, a polyline's vertex, an offset's numerical
//! joint, and through a placement at a station -- and a window holding two
//! seams, or with a seam on its edge, is refused by name. A window of `0`
//! is today's reading.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{Circle2, Curve2, Curve3, Elevated3, ElevationLaw, Line2, Polyline2};
use axiolid_mesh_compile::station::{curve_path, placement, resolve};
use axiolid_model::{
    CurveRelation, CurveSegment, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode,
    InstanceAtStation, NodeId, OrientedCurveStation, SeamSide, Station, StationFrame,
    StationOffsets, Transition, TrimSelector, TrimmingPreference,
};
use axiolid_reference::station::{
    CompositeBasis, SEAM_WINDOW_AMBIGUOUS, SEAM_WINDOW_STRADDLED, WINDOW_CERTIFICATION,
};

const PI: Scalar = core::f64::consts::PI;
const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];
/// The model's precision: a millimetre.
const WINDOW: Scalar = 1e-3;
/// How far from the joint the station is authored: inside the window.
const NEAR: Scalar = 6e-4;

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, input: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input: found, .. }) => assert_eq!(found, input),
        other => panic!("expected a refusal naming {input:?}, got {other:?}"),
    }
}

fn trimmed(b: &mut GeometryGraphBuilder, basis: NodeId, start: Scalar, end: Scalar) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![TrimSelector::Parameter(start)],
        end: vec![TrimSelector::Parameter(end)],
        sense_agreement: true,
        preference: TrimmingPreference::Parameter,
    }))
    .unwrap()
}

fn composite(b: &mut GeometryGraphBuilder, pieces: &[NodeId]) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments: pieces
            .iter()
            .map(|&curve| CurveSegment {
                curve,
                same_sense: true,
                transition: Transition::Discontinuous,
            })
            .collect(),
    }))
    .unwrap()
}

fn line2(b: &mut GeometryGraphBuilder, origin: Point2, direction: Vec2) -> NodeId {
    b.push_value(Curve2::Line(Line2 { origin, direction }))
        .unwrap()
}

/// `(0, 0)` to `(10, 0)`, then a quarter of radius 5 about `(15, 0)`
/// turning from `+y` towards `+x`: a right-angle corner at 10.
fn corner(b: &mut GeometryGraphBuilder) -> NodeId {
    let line = line2(b, Point2::ZERO, Vec2::X);
    let circle = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(15.0, 0.0),
                x: -Vec2::X,
                y: Vec2::Y,
            },
            radius: 5.0,
        }))
        .unwrap();
    let straight = trimmed(b, line, 0.0, 10.0);
    let arc = trimmed(b, circle, 0.0, 0.5 * PI);
    composite(b, &[straight, arc])
}

fn station(basis: NodeId, d: Scalar, side: SeamSide, window: Scalar) -> OrientedCurveStation {
    CurveStation::new(basis, Station::new(d, StationOffsets::new(0.5, 0.25, 0.0)))
        .with_seam_side(side)
        .with_seam_window(window)
}

/// Resolve `stations` in one graph, built on `basis` by `build`.
fn resolved(
    build: impl Fn(&mut GeometryGraphBuilder) -> NodeId,
    stations: &[(Scalar, SeamSide, Scalar)],
) -> Vec<GeomResult<axiolid_mesh_compile::station::ResolvedStation>> {
    let mut b = GeometryGraphBuilder::new();
    let basis = build(&mut b);
    let ids: Vec<NodeId> = stations
        .iter()
        .map(|&(d, side, window)| b.push_value(station(basis, d, side, window)).unwrap())
        .collect();
    let graph = b.finish(ids.clone()).unwrap();
    ids.into_iter().map(|id| resolve(&graph, id)).collect()
}

#[test]
fn a_station_within_its_window_of_a_joint_reads_the_joint_from_its_side() {
    for side in SIDES {
        for near in [10.0 - NEAR, 10.0 + NEAR] {
            let out = resolved(
                corner,
                &[(near, side, WINDOW), (10.0, side, 0.0), (near, side, 0.0)],
            );
            let (snapped, on, unsnapped) = (
                out[0].as_ref().unwrap(),
                out[1].as_ref().unwrap(),
                out[2].as_ref().unwrap(),
            );
            // Exactly as the station placed on the joint, bitwise.
            assert_eq!(snapped, on, "{side:?} at {near}");
            // Without a window it stands where it is authored, today's
            // reading: a different point.
            assert!((unsnapped.point - on.point).length() > 0.5 * NEAR);
        }
        // Incoming reads the line, outgoing the arc's start, at the joint.
        let on = resolved(corner, &[(10.0 - NEAR, side, WINDOW)])
            .remove(0)
            .unwrap();
        let want = match side {
            SeamSide::Incoming => axiolid_core::Vec3::X,
            _ => axiolid_core::Vec3::Y,
        };
        assert!((on.frame.x - want).length() < 1e-12, "{side:?}");
    }
    // Outside the window nothing moves.
    let out = resolved(
        corner,
        &[
            (10.0 - 2.0 * WINDOW, SeamSide::Incoming, WINDOW),
            (10.0 - 2.0 * WINDOW, SeamSide::Incoming, 0.0),
        ],
    );
    assert_eq!(out[0].as_ref().unwrap(), out[1].as_ref().unwrap());
}

#[test]
fn two_joints_in_one_window_are_refused_by_name() {
    let tiny = 5e-4;
    let build = |b: &mut GeometryGraphBuilder| {
        let line = line2(b, Point2::ZERO, Vec2::X);
        let a = trimmed(b, line, 0.0, 10.0);
        let c = trimmed(b, line, 10.0, 10.0 + tiny);
        let d = trimmed(b, line, 10.0 + tiny, 20.0);
        composite(b, &[a, c, d])
    };
    let out = resolved(build, &[(10.0 + 0.5 * tiny, SeamSide::Incoming, WINDOW)]);
    refused(out.into_iter().next().unwrap(), SEAM_WINDOW_AMBIGUOUS);
    // The same station with no window reads as before.
    assert!(resolved(build, &[(10.0 + 0.5 * tiny, SeamSide::Incoming, 0.0)])[0].is_ok());
}

#[test]
fn an_atomic_polyline_snaps_onto_its_vertex() {
    let build = |b: &mut GeometryGraphBuilder| {
        b.push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::ZERO, Point2::new(3.0, 4.0), Point2::new(3.0, 10.0)],
            closed: false,
        }))
        .unwrap()
    };
    for side in SIDES {
        let out = resolved(build, &[(5.0 + NEAR, side, WINDOW), (5.0, side, 0.0)]);
        assert_eq!(out[0].as_ref().unwrap(), out[1].as_ref().unwrap());
    }
}

/// An offset by distances at stations along a road on a crest: its joint
/// at the middle station is at the offset's own plan length, a quadrature.
fn kerb(b: &mut GeometryGraphBuilder) -> NodeId {
    let road = b
        .push_value(Curve3::Elevated(Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            ElevationLaw::Polynomial {
                coefficients: vec![10.0, 0.02, -0.0002],
            },
        )))
        .unwrap();
    b.push(GeometryNode::CurveRelation(
        CurveRelation::OffsetByStations {
            basis: road,
            stations: vec![
                Station::new(0.0, StationOffsets::new(3.0, -0.2, 0.0)),
                Station::new(40.0, StationOffsets::new(3.5, 0.1, 0.0)),
                Station::new(100.0, StationOffsets::new(4.5, 0.3, 0.0)),
            ],
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

/// The kerb's joint where the kernel reads it, and its certified interval
/// at a window's certification.
fn kerb_joint() -> (Scalar, axiolid_reference::station::DistanceInterval) {
    let mut b = GeometryGraphBuilder::new();
    let id = kerb(&mut b);
    let graph = b.finish(vec![id]).unwrap();
    let path = curve_path(&graph, id).unwrap();
    let basis = CompositeBasis::from_path(&path).unwrap();
    (
        basis.pieces()[0].length(),
        basis
            .joint_intervals(WINDOW * WINDOW_CERTIFICATION)
            .unwrap()[0],
    )
}

#[test]
fn an_offsets_numerical_joint_snaps_and_its_window_edge_is_refused() {
    let (joint, certified) = kerb_joint();
    assert!(certified.contains(joint) && certified.width() > 0.0);
    for side in SIDES {
        for near in [joint - NEAR, joint + NEAR] {
            let out = resolved(kerb, &[(near, side, WINDOW), (joint, side, 0.0)]);
            assert_eq!(out[0].as_ref().unwrap(), out[1].as_ref().unwrap());
        }
    }
    // The window's edge through the joint's certified interval.
    for edge in [certified.midpoint() - WINDOW, certified.midpoint() + WINDOW] {
        let out = resolved(kerb, &[(edge, SeamSide::Incoming, WINDOW)]);
        refused(out.into_iter().next().unwrap(), SEAM_WINDOW_STRADDLED);
    }
}

/// A node placed at a windowed station is placed as at the joint.
#[test]
fn a_placement_at_a_windowed_station_is_the_placement_at_the_joint() {
    let (joint, _) = kerb_joint();
    let mut b = GeometryGraphBuilder::new();
    let basis = kerb(&mut b);
    let source = line2(&mut b, Point2::ZERO, Vec2::Y);
    let windowed = b
        .push(GeometryNode::InstanceAtStation(InstanceAtStation::new(
            source,
            station(basis, joint + NEAR, SeamSide::Incoming, WINDOW),
        )))
        .unwrap();
    let on = b
        .push(GeometryNode::InstanceAtStation(InstanceAtStation::new(
            source,
            station(basis, joint, SeamSide::Incoming, 0.0),
        )))
        .unwrap();
    let ambiguous_basis = corner(&mut b);
    let off = b
        .push(GeometryNode::InstanceAtStation(InstanceAtStation::new(
            source,
            station(
                ambiguous_basis,
                10.0 - 2.0 * WINDOW,
                SeamSide::Incoming,
                0.0,
            ),
        )))
        .unwrap();
    let graph: GeometryGraph = b.finish(vec![windowed, on, off]).unwrap();
    assert_eq!(
        placement(&graph, windowed).unwrap(),
        placement(&graph, on).unwrap()
    );
    assert!(placement(&graph, off).is_ok());
}

/// A curve placed at a windowed station, used as a basis itself, is
/// carried by the snapped placement: the window is read wherever the
/// oriented station is resolved.
#[test]
fn a_curve_placed_at_a_windowed_station_is_carried_from_the_joint() {
    let mut b = GeometryGraphBuilder::new();
    let basis = corner(&mut b);
    let source = line2(&mut b, Point2::ZERO, Vec2::X);
    let place = |b: &mut GeometryGraphBuilder, d: Scalar, window: Scalar| {
        b.push(GeometryNode::InstanceAtStation(InstanceAtStation::new(
            source,
            station(basis, d, SeamSide::Outgoing, window),
        )))
        .unwrap()
    };
    let windowed = place(&mut b, 10.0 - NEAR, WINDOW);
    let on = place(&mut b, 10.0, 0.0);
    let along = |b: &mut GeometryGraphBuilder, placed: NodeId| {
        b.push_value(CurveStation::new(placed, Station::at(0.25)))
            .unwrap()
    };
    let (a, c) = (along(&mut b, windowed), along(&mut b, on));
    let graph = b.finish(vec![a, c]).unwrap();
    assert_eq!(resolve(&graph, a).unwrap(), resolve(&graph, c).unwrap());
}
