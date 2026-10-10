//! Stations along offset curve relations (#289): a constant offset and an
//! offset by distances at stations as a station basis, measured in the
//! offset's own length, nested in composites and placements; runs of
//! sections, placements, seams and curve paths along them; and the
//! refusals. Expected values are closed forms or an independent evaluation
//! written out here, never the offset code under test.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Line2, Line3, PathCurve};
use axiolid_measure::volume_properties;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::{curve_path, placement, resolve, seams, ResolvedStation};
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveSegment, CurveStation, GeometryGraphBuilder, GeometryNode,
    InstanceAtStation, NodeId, OrientedCurveStation, SeamSide, SolidOperation, Station,
    StationFrame, StationOffsets, StationedSection, Transition, TrimSelector, TrimmingPreference,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::station::CompositeBasis;

const PI: Scalar = core::f64::consts::PI;
const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    let off = (actual - expected).length();
    assert!(
        off <= eps,
        "{what}: {actual:?} != {expected:?} (off by {off:e})"
    );
}

fn refused_naming<T: core::fmt::Debug>(result: GeomResult<T>, needle: &str) {
    match result {
        Err(GeomError::InvalidInput(detail) | GeomError::Degenerate(detail)) => {
            assert!(detail.contains(needle), "expected {needle:?} in {detail:?}");
        }
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert!(input.contains(needle), "expected {needle:?} in {input:?}");
        }
        other => panic!("expected a refusal naming {needle:?}, got {other:?}"),
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
                transition: Transition::ContinuousSameGradient,
            })
            .collect(),
    }))
    .unwrap()
}

fn offset(b: &mut GeometryGraphBuilder, basis: NodeId, distance: Scalar) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Offset {
        basis,
        distance,
        reference_direction: None,
    }))
    .unwrap()
}

fn by_stations(b: &mut GeometryGraphBuilder, basis: NodeId, stations: Vec<Station>) -> NodeId {
    b.push(GeometryNode::CurveRelation(
        CurveRelation::OffsetByStations {
            basis,
            stations,
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

fn at(basis: NodeId, d: Scalar, side: SeamSide) -> OrientedCurveStation {
    CurveStation::new(basis, Station::at(d)).with_seam_side(side)
}

fn station(b: &mut GeometryGraphBuilder, basis: NodeId, d: Scalar, side: SeamSide) -> NodeId {
    b.push_value(at(basis, d, side)).unwrap()
}

fn line2(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    }))
    .unwrap()
}

fn circle2(b: &mut GeometryGraphBuilder, centre: Point2, radius: Scalar) -> NodeId {
    b.push_value(Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: centre,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    }))
    .unwrap()
}

/// `(0, 0)` to `(10, 0)`, then a quarter of radius 5 about `(10, 5)`
/// turning left, tangent at the joint.
fn line_and_arc(b: &mut GeometryGraphBuilder) -> NodeId {
    let line = line2(b);
    let circle = circle2(b, Point2::new(10.0, 5.0), 5.0);
    let straight = trimmed(b, line, 0.0, 10.0);
    let arc = trimmed(b, circle, 1.5 * PI, 2.0 * PI);
    composite(b, &[straight, arc])
}

/// The section a resolved station presents, as `(point, tangent, lateral)`.
fn read(resolved: &ResolvedStation) -> (Point3, Vec3, Vec3) {
    (resolved.point, resolved.frame.x, -resolved.frame.z)
}

#[test]
fn a_constant_offset_of_a_line_and_arc_resolves_in_its_own_length() {
    let mut b = GeometryGraphBuilder::new();
    let base = line_and_arc(&mut b);
    let left = offset(&mut b, base, 1.0);
    let right = offset(&mut b, base, -1.0);
    let length = 10.0 + 2.0 * PI;
    let mut ids = Vec::new();
    for d in [0.0, 4.0, 10.0, 10.0 + PI, length] {
        for side in SIDES {
            ids.push((d, side, station(&mut b, left, d, side)));
        }
    }
    let beyond = station(&mut b, left, length + 1e-6, SeamSide::Outgoing);
    let on_right = station(&mut b, right, 10.0 + 3.0 * PI, SeamSide::Outgoing);
    let graph = b.finish(vec![beyond, on_right]).unwrap();
    for (d, side, id) in ids {
        let (point, tangent, lateral) = read(&resolve(&graph, id).unwrap());
        let (expected, heading) = if d <= 10.0 {
            // The line's offset, y = 1; both sides of the joint at
            // (10, 1) heading +x.
            (Point3::new(d, 1.0, 0.0), Vec3::X)
        } else {
            // The arc's offset, radius 4: the arc length is scaled by 4/5.
            let theta = 1.5 * PI + (d - 10.0) / 4.0;
            let (sin, cos) = theta.sin_cos();
            (
                Point3::new(10.0 + 4.0 * cos, 5.0 + 4.0 * sin, 0.0),
                Vec3::new(-sin, cos, 0.0),
            )
        };
        let what = format!("{d} {side:?}");
        close3(point, expected, 1e-11, &format!("{what}: point"));
        close3(tangent, heading, 1e-11, &format!("{what}: tangent"));
        close3(
            lateral,
            Vec3::new(-heading.y, heading.x, 0.0),
            1e-11,
            &format!("{what}: lateral"),
        );
    }
    refused_naming(resolve(&graph, beyond), "beyond");
    // To the right of the arc is radius 6: the base's 5 pi / 2 is 3 pi.
    let (point, _, _) = read(&resolve(&graph, on_right).unwrap());
    close3(
        point,
        Point3::new(16.0, 5.0, 0.0),
        1e-11,
        "right, at the arc's end",
    );
    // The joint is a seam, exact: both pieces are closed forms.
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let seams = seams(&graph, left, &options).unwrap();
    assert_eq!(seams.len(), 1);
    assert!((seams[0].distance - 10.0).abs() <= 1e-12 && !seams[0].smooth);
}

/// The offset of a counter-clockwise arc of radius `r` about the origin,
/// from angle 0, whose lateral offset grows linearly from 0 to `w` over
/// `span` radians: its point at `theta` and its own arc length to there,
/// by Simpson's rule.
fn spiral(r: Scalar, w: Scalar, span: Scalar, theta: Scalar) -> (Point3, Scalar) {
    let rho = |t: Scalar| r - w * t / span;
    let speed = |t: Scalar| (rho(t) * rho(t) + (w / span) * (w / span)).sqrt();
    let n = 20_000;
    let h = theta / Scalar::from(n);
    let mut sum = speed(0.0) + speed(theta);
    for k in 1..n {
        sum += speed(h * Scalar::from(k)) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    let (sin, cos) = theta.sin_cos();
    (
        Point3::new(rho(theta) * cos, rho(theta) * sin, 0.0),
        sum * h / 3.0,
    )
}

#[test]
fn an_offset_by_distances_along_a_line_and_an_arc_agrees_with_an_independent_evaluation() {
    let mut b = GeometryGraphBuilder::new();
    // Along a line, lateral 1 -> 3 and vertical 0 -> 0.5 between
    // distances 2 and 12: the chord from (2, 1, 0) to (12, 3, 0.5).
    let line = line2(&mut b);
    let widening = by_stations(
        &mut b,
        line,
        vec![
            Station::new(2.0, StationOffsets::new(1.0, 0.0, 0.0)),
            Station::new(12.0, StationOffsets::new(3.0, 0.5, 0.0)),
        ],
    );
    let chord = Vec3::new(10.0, 2.0, 0.5);
    let fractions = [0.0, 0.3, 0.75, 1.0];
    let on_line: Vec<NodeId> = fractions
        .iter()
        .map(|f| station(&mut b, widening, f * chord.length(), SeamSide::Outgoing))
        .collect();
    // Along a quarter arc of radius 10, lateral 0 -> 2 (to the centre).
    let (r, w, span) = (10.0, 2.0, 0.5 * PI);
    let circle = circle2(&mut b, Point2::ZERO, r);
    let tightening = by_stations(
        &mut b,
        circle,
        vec![
            Station::at(0.0),
            Station::new(0.5 * r * span, StationOffsets::new(0.5 * w, 0.0, 0.0)),
            Station::new(r * span, StationOffsets::new(w, 0.0, 0.0)),
        ],
    );
    let thetas = [0.2, 0.7, 1.3];
    let on_arc: Vec<(NodeId, Point3)> = thetas
        .iter()
        .map(|&theta| {
            let (point, length) = spiral(r, w, span, theta);
            (
                station(&mut b, tightening, length, SeamSide::Outgoing),
                point,
            )
        })
        .collect();
    let graph = b.finish(vec![widening, tightening]).unwrap();
    for (f, id) in fractions.iter().zip(on_line) {
        let (point, tangent, _) = read(&resolve(&graph, id).unwrap());
        close3(
            point,
            Point3::new(2.0, 1.0, 0.0) + *f * chord,
            1e-8,
            &format!("line at {f}"),
        );
        close3(tangent, chord.normalize(), 1e-8, "line tangent");
    }
    for (id, expected) in on_arc {
        let (point, _, _) = read(&resolve(&graph, id).unwrap());
        close3(point, expected, 1e-7, "arc");
    }
    // The middle station is a break of the distance law, so a seam; it
    // lies after a numerically measured piece, so its exact reading is
    // refused typed.
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    assert!(matches!(
        seams(&graph, tightening, &options),
        Err(GeomError::UnsupportedInput { input, .. }) if input.contains("exact seam distances")
    ));
    assert!(seams(&graph, widening, &options).unwrap().is_empty());
}

#[test]
fn a_3d_offset_and_an_offset_by_stations_in_plan_read_their_laws() {
    let mut b = GeometryGraphBuilder::new();
    // A 3D offset lies along V x T: +Z x +X is +Y. An untrimmed line under
    // an offset is its parameter domain, as inside a composite: 10 long.
    let line3 = b
        .push_value(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: 10.0 * Vec3::X,
        }))
        .unwrap();
    let up = b
        .push(GeometryNode::CurveRelation(CurveRelation::Offset {
            basis: line3,
            distance: 2.0,
            reference_direction: Some(Vec3::Z),
        }))
        .unwrap();
    let down = b
        .push(GeometryNode::CurveRelation(CurveRelation::Offset {
            basis: line3,
            distance: 2.0,
            reference_direction: Some(-Vec3::Z),
        }))
        .unwrap();
    // One up from a 10% grade: +Z in the plan frame, the leaning up in
    // the section frame; both a rigid lift, so the plan distance is kept.
    let grade = 0.1;
    let elevated = b
        .push_value(Curve3::Elevated(axiolid_curve::Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            axiolid_curve::ElevationLaw::constant_grade(2.0, grade),
        )))
        .unwrap();
    let lifted = |frame| {
        GeometryNode::CurveRelation(CurveRelation::OffsetByStations {
            basis: elevated,
            stations: vec![
                Station::new(0.0, StationOffsets::new(0.0, 1.0, 0.0)),
                Station::new(10.0, StationOffsets::new(0.0, 1.0, 0.0)),
            ],
            frame,
        })
    };
    let plan = b.push(lifted(StationFrame::Plan)).unwrap();
    let section = b.push(lifted(StationFrame::Section)).unwrap();
    let ids =
        [up, down, plan, section].map(|basis| station(&mut b, basis, 5.0, SeamSide::Outgoing));
    let graph = b.finish(ids.to_vec()).unwrap();
    let point = |k: usize| resolve(&graph, ids[k]).unwrap().point;
    close3(point(0), Point3::new(5.0, 2.0, 0.0), 1e-12, "along +Z x T");
    close3(point(1), Point3::new(5.0, -2.0, 0.0), 1e-12, "along -Z x T");
    let lean = (1.0 + grade * grade).sqrt();
    close3(point(2), Point3::new(5.0, 0.0, 3.5), 1e-8, "plan lift");
    close3(
        point(3),
        Point3::new(5.0 - grade / lean, 0.0, 2.5 + 1.0 / lean),
        1e-8,
        "section lift",
    );
}

#[test]
fn a_node_placed_at_a_station_of_an_offset_is_exact_only_beside_a_line() {
    let mut b = GeometryGraphBuilder::new();
    let base = line_and_arc(&mut b);
    let left = offset(&mut b, base, 1.0);
    let source = b
        .push_value(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }))
        .unwrap();
    let on_line = b
        .push_value(InstanceAtStation::new(
            source,
            at(left, 4.0, SeamSide::Outgoing),
        ))
        .unwrap();
    let on_arc = b
        .push_value(InstanceAtStation::new(
            source,
            at(left, 10.0 + PI, SeamSide::Outgoing),
        ))
        .unwrap();
    // A station along the placed line: the offset composite carried.
    let along = station(&mut b, on_line, 0.5, SeamSide::Outgoing);
    let graph = b.finish(vec![on_line, on_arc, along]).unwrap();
    let placed = placement(&graph, on_line).unwrap();
    assert!(placed.exact, "beside the line the frame is exact");
    close3(
        placed.transform.transform_point3(Point3::ZERO),
        Point3::new(4.0, 1.0, 0.0),
        1e-12,
        "origin",
    );
    close3(
        placed.transform.transform_vector3(Vec3::X),
        Vec3::X,
        1e-12,
        "x",
    );
    close3(
        placed.transform.transform_vector3(Vec3::Y),
        Vec3::Y,
        1e-12,
        "y",
    );
    let curved = placement(&graph, on_arc).unwrap();
    assert!(
        !curved.exact,
        "beside the arc the frame is read numerically"
    );
    let theta = 1.5 * PI + PI / 4.0;
    close3(
        curved.transform.transform_point3(Point3::ZERO),
        Point3::new(10.0 + 4.0 * theta.cos(), 5.0 + 4.0 * theta.sin(), 0.0),
        1e-11,
        "on the arc",
    );
    close3(
        resolve(&graph, along).unwrap().point,
        Point3::new(4.5, 1.0, 0.0),
        1e-12,
        "along the placed line",
    );
}

#[test]
fn a_curve_path_of_an_offset_reads_as_its_resolved_stations() {
    let mut b = GeometryGraphBuilder::new();
    let base = line_and_arc(&mut b);
    let left = offset(&mut b, base, 1.0);
    let ids: Vec<(Scalar, SeamSide, NodeId)> = [0.0, 10.0, 11.5]
        .into_iter()
        .flat_map(|d| SIDES.map(|side| (d, side)))
        .collect::<Vec<_>>()
        .into_iter()
        .map(|(d, side)| (d, side, station(&mut b, left, d, side)))
        .collect();
    let graph = b.finish(vec![left]).unwrap();
    let path = curve_path(&graph, left).unwrap();
    assert_eq!(path.pieces().len(), 2);
    assert!(path
        .pieces()
        .iter()
        .all(|piece| matches!(piece.curve, PathCurve::Offset(_))));
    let composite = CompositeBasis::from_path(&path).unwrap();
    for (d, side, id) in ids {
        let resolved = resolve(&graph, id).unwrap();
        assert_eq!(
            composite.section_on(d, side).unwrap(),
            resolved.section,
            "{d} {side:?}: bitwise"
        );
    }
}

#[test]
fn a_run_of_sections_along_an_offset_meshes_its_own_length() {
    // A constant rectangle at stations 0 and 10 + 2 pi along the left
    // offset of the line and arc: volume A * (10 + 2 pi), Pappus on the arc
    // about the offset's own centreline.
    let mut b = GeometryGraphBuilder::new();
    let base = line_and_arc(&mut b);
    let left = offset(&mut b, base, 1.0);
    let profile = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 0.5,
            y: 0.4,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let length = 10.0 + 2.0 * PI;
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::StationedSpine {
                directrix: left,
                sections: [0.0, 10.0, length]
                    .map(|d| StationedSection {
                        profile,
                        station: Station::at(d),
                    })
                    .to_vec(),
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, root, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .unwrap();
    let volume = volume_properties(&mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume
        .abs();
    let expected = 0.2 * length;
    assert!(
        (volume - expected).abs() <= 1e-3 * expected,
        "volume {volume} != {expected}"
    );
    for p in &mesh.positions {
        assert!(
            p.y >= 1.0 - 0.25 - 1e-9,
            "beside the offset, not the base: {p:?}"
        );
    }
}

#[test]
fn offsets_as_station_bases_refuse_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let corner = b
        .push_value(Curve2::Polyline(axiolid_curve::Polyline2 {
            points: vec![
                Point2::ZERO,
                Point2::new(10.0, 0.0),
                Point2::new(10.0, 10.0),
            ],
            closed: false,
        }))
        .unwrap();
    let across = offset(&mut b, corner, 1.0);
    let circle = circle2(&mut b, Point2::ZERO, 5.0);
    let collapsed = offset(&mut b, circle, 5.0);
    let line3 = b
        .push_value(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }))
        .unwrap();
    let planar3 = offset(&mut b, line3, 1.0);
    let line = line2(&mut b);
    let once = offset(&mut b, line, 1.0);
    let twice = offset(&mut b, once, 1.0);
    let trim = trimmed(&mut b, once, 0.0, 0.5);
    let ids = [across, collapsed, planar3, twice, trim]
        .map(|basis| station(&mut b, basis, 0.5, SeamSide::Outgoing));
    let graph = b.finish(ids.to_vec()).unwrap();
    refused_naming(resolve(&graph, ids[0]), "an offset across a corner");
    refused_naming(resolve(&graph, ids[1]), "collapses");
    refused_naming(
        resolve(&graph, ids[2]),
        "a planar (2D) offset of a 3D curve",
    );
    refused_naming(resolve(&graph, ids[3]), "an offset of an offset");
    refused_naming(resolve(&graph, ids[4]), "a trim of an offset curve");
    assert!(matches!(
        resolve(&graph, ids[4]),
        Err(GeomError::UnsupportedInput {
            operation: Operation::CurveEvaluation,
            ..
        })
    ));
}

#[test]
fn a_station_along_an_offset_round_trips_through_the_graph() {
    let mut b = GeometryGraphBuilder::new();
    let base = line_and_arc(&mut b);
    let left = offset(&mut b, base, 1.0);
    let pushed = at(left, 10.0, SeamSide::Incoming);
    let id = b.push_value(pushed).unwrap();
    let graph = b.finish(vec![id]).unwrap();
    let Some(GeometryNode::OrientedCurveStation(stored)) = graph.get(id) else {
        panic!("not an oriented station");
    };
    assert_eq!(*stored, pushed);
    assert_eq!(graph.get(id).unwrap().references(), vec![left]);
    assert_eq!(graph.get(left).unwrap().references(), vec![base]);
    close3(
        resolve(&graph, id).unwrap().point,
        Point3::new(10.0, 1.0, 0.0),
        1e-12,
        "resolved after the round trip",
    );
}
