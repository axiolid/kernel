//! Nodes placed at stations (#264): a curve or a solid in the frame of an
//! oriented curve station, its local `x`, `y`, `z` on the station's
//! tangent, left lateral and up. Frames are checked against closed forms
//! computed here, never against the evaluator under test; exactness is
//! claimed on a line basis only; stations along a placed curve work where
//! the placement keeps `+Z` and are refused by name where it tilts it.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Elevated3, ElevationLaw, Line2, Line3, Polyline2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::{placement, resolve, seams, ResolvedPlacement};
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, DeviationReport, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_model::{
    CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, InstanceAtStation, NodeId,
    OrientedCurveStation, SeamSide, SolidOperation, Station, StationFrame, StationOffsets,
    StationOrientation,
};
use axiolid_profile::{Profile, RectangleProfile};

const EPS: Scalar = 1e-9;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn compile(graph: &GeometryGraph, root: NodeId) -> GeomResult<(TriMesh, DeviationReport)> {
    compiler()
        .compile_mesh_with_deviation(graph, root, &options())
        .map(|(outcome, report)| (outcome.mesh, report))
}

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?}"
    );
}

/// The placement's columns against a frame computed by hand.
fn assert_placed(placed: &ResolvedPlacement, [t, l, u]: [Vec3; 3], point: Point3, eps: Scalar) {
    let m = placed.transform;
    close3(
        m.transform_vector3(Vec3::X),
        t,
        eps,
        "local x on the tangent",
    );
    close3(
        m.transform_vector3(Vec3::Y),
        l,
        eps,
        "local y on the lateral",
    );
    close3(m.transform_vector3(Vec3::Z), u, eps, "local z on up");
    close3(m.translation, point, eps, "local origin on the point");
    close3(placed.station.point, point, eps, "the station's point");
}

fn bounds(mesh: &TriMesh) -> (Point3, Point3) {
    mesh.positions.iter().fold(
        (
            Point3::splat(Scalar::INFINITY),
            Point3::splat(-Scalar::INFINITY),
        ),
        |(lo, hi), p| (lo.min(*p), hi.max(*p)),
    )
}

fn placement_contribution(report: &DeviationReport) -> Option<DeviationBound> {
    report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationPlacement)
        .map(|c| c.bound)
}

/// A plan line along `+y` from `(2, 1)`: tangent `+y`, lateral `-x`, up
/// `+Z`, its parameter half the distance.
fn plan_line(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::new(2.0, 1.0),
        direction: Vec2::new(0.0, 2.0),
    }))
    .unwrap()
}

/// A local 2D line along local `x`.
fn local_line(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    }))
    .unwrap()
}

/// A 1 x 2 rectangle extruded 3 along local `+z`.
fn block(b: &mut GeometryGraphBuilder) -> NodeId {
    let profile = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 1.0,
            y: 2.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 3.0,
    }))
    .unwrap()
}

fn swept(b: &mut GeometryGraphBuilder, directrix: NodeId, range: (Scalar, Scalar)) -> NodeId {
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix,
        radius: 0.1,
        inner_radius: None,
        parameter_range: Some(range),
        fillet_radius: None,
    }))
    .unwrap()
}

fn at(basis: NodeId, station: Station, frame: StationFrame) -> OrientedCurveStation {
    CurveStation {
        basis,
        station,
        frame,
    }
    .into()
}

/// Every vertex of `mesh` within `radius` (and not far inside it) of the
/// line through `point` along the unit `direction`.
fn hugs_line(mesh: &TriMesh, point: Point3, direction: Vec3, radius: Scalar) {
    for p in &mesh.positions {
        let off = *p - point;
        let distance = (off - off.dot(direction) * direction).length();
        assert!(
            distance <= radius + 1e-6 && distance >= 0.9 * radius,
            "{p:?} is {distance} from the placed line"
        );
    }
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

// --- a straight base: exact ----------------------------------------------------

#[test]
fn a_line_placed_at_a_station_on_a_straight_base_is_exact() {
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let source = local_line(&mut b);
    // Base point (2, 4, 0); offsets lateral 0.5 (along -x), vertical 0.25,
    // longitudinal 1 (along +y): (1.5, 5, 0.25).
    let station = at(
        basis,
        Station::new(3.0, StationOffsets::new(0.5, 0.25, 1.0)),
        StationFrame::Section,
    );
    let placed = b
        .push_value(InstanceAtStation::new(source, station))
        .unwrap();
    let tube = swept(&mut b, placed, (0.0, 4.0));
    let graph = b.finish(vec![placed, tube]).unwrap();

    let resolved = placement(&graph, placed).unwrap();
    assert!(resolved.exact, "a line's frame is exact");
    // Exact to the bit: every value here is representable.
    assert_eq!(
        resolved.transform,
        Transform3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z, Vec3::new(1.5, 5.0, 0.25))
    );

    // The placed line as a sweep directrix runs from the point along +y.
    let (mesh, report) = compile(&graph, tube).unwrap();
    hugs_line(&mesh, Point3::new(1.5, 5.0, 0.25), Vec3::Y, 0.1);
    let (lo, hi) = bounds(&mesh);
    assert!((lo.y - 5.0).abs() <= EPS && (hi.y - 9.0).abs() <= EPS);
    assert_eq!(placement_contribution(&report), None);
}

#[test]
fn a_solid_placed_on_a_straight_base_compiles_exactly_and_meshes_in_place() {
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let solid = block(&mut b);
    let station = at(
        basis,
        Station::new(3.0, StationOffsets::new(0.5, 0.25, 1.0)),
        StationFrame::Section,
    );
    let placed = b
        .push_value(InstanceAtStation::new(solid, station))
        .unwrap();
    let graph = b.finish(vec![placed]).unwrap();

    // Local x (width 1) along +y, local y (width 2) along -x, local z
    // (depth 3) up, around (1.5, 5, 0.25).
    let (mesh, report) = compile(&graph, placed).unwrap();
    let (lo, hi) = bounds(&mesh);
    close3(lo, Vec3::new(0.5, 4.5, 0.25), EPS, "low corner");
    close3(hi, Vec3::new(2.5, 5.5, 3.25), EPS, "high corner");
    assert_eq!(placement_contribution(&report), None);
    assert!(report.bound.is_some(), "{report:?}");

    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, placed, &options())
        .unwrap();
    let volume = exact_properties(&exact, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    assert!((volume - 6.0).abs() <= 1e-9, "{volume}");
}

// --- curved, elevated and seamed bases: an independent frame -----------------

#[test]
fn a_placement_on_an_arc_matches_an_independent_frame_and_is_not_exact() {
    // R = 30 counter-clockwise, 0.4 rad round: t = (-sin, cos, 0),
    // l = (-cos, -sin, 0), u = +Z. Offsets 2 left (towards the centre) and
    // 1 up. Axis (0, 1, 0) in (t, l, u) is the lateral, so up' = l,
    // t' = t and l' = up' x t' = l x t = -Z.
    let (radius, angle) = (30.0, 0.4_f64);
    let (s, c) = angle.sin_cos();
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }))
        .unwrap();
    let source = local_line(&mut b);
    let solid = block(&mut b);
    let station = OrientedCurveStation::new(
        CurveStation::new(
            basis,
            Station::new(radius * angle, StationOffsets::new(2.0, 1.0, 0.0)),
        ),
        StationOrientation::new(Some(Vec3::Y), None),
    );
    let placed = b
        .push_value(InstanceAtStation::new(source, station))
        .unwrap();
    let placed_solid = b
        .push_value(InstanceAtStation::new(solid, station))
        .unwrap();
    let tube = swept(&mut b, placed, (-3.0, 3.0));
    let graph = b.finish(vec![placed, placed_solid, tube]).unwrap();

    let t = Vec3::new(-s, c, 0.0);
    let l = Vec3::new(-c, -s, 0.0);
    let point = Point3::new((radius - 2.0) * c, (radius - 2.0) * s, 1.0);
    let resolved = placement(&graph, placed).unwrap();
    assert_placed(&resolved, [t, -Vec3::Z, l], point, EPS);
    assert!(!resolved.exact, "an arc's distance is read numerically");

    // The placed line, swept: along t' through the point.
    let (mesh, _) = compile(&graph, tube).unwrap();
    hugs_line(&mesh, point, t, 0.1);

    // A solid placed there meshes, with its placement named unbounded, and
    // is refused by the exact compiler by name.
    let (_, report) = compile(&graph, placed_solid).unwrap();
    assert!(
        matches!(
            placement_contribution(&report),
            Some(DeviationBound::Unbounded(_))
        ),
        "{report:?}"
    );
    assert_eq!(report.bound, None);
    refused_naming(
        ReferenceExactCompiler::new().compile_exact(&graph, placed_solid, &options()),
        "frame is not exact",
    );
}

#[test]
fn a_placement_on_an_elevated_base_leans_with_the_grade_or_stands_in_plan() {
    // 2% from height 100 along +x: at 40 the point is (40, 0, 100.8),
    // t = (1, 0, 0.02) / k, l = +y, u = (-0.02, 0, 1) / k.
    let k = (1.0_f64 + 0.02 * 0.02).sqrt();
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve3::Elevated(Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            ElevationLaw::constant_grade(100.0, 0.02),
        )))
        .unwrap();
    let source = local_line(&mut b);
    let section = b
        .push_value(InstanceAtStation::new(
            source,
            at(basis, Station::at(40.0), StationFrame::Section),
        ))
        .unwrap();
    let plan = b
        .push_value(InstanceAtStation::new(
            source,
            at(basis, Station::at(40.0), StationFrame::Plan),
        ))
        .unwrap();
    let graph = b.finish(vec![section, plan]).unwrap();
    let point = Point3::new(40.0, 0.0, 100.8);
    let leaning = placement(&graph, section).unwrap();
    assert_placed(
        &leaning,
        [
            Vec3::new(1.0, 0.0, 0.02) / k,
            Vec3::Y,
            Vec3::new(-0.02, 0.0, 1.0) / k,
        ],
        point,
        EPS,
    );
    assert!(!leaning.exact, "only a line's frame is claimed exact");
    assert_placed(
        &placement(&graph, plan).unwrap(),
        [Vec3::X, Vec3::Y, Vec3::Z],
        point,
        EPS,
    );
}

#[test]
fn a_placement_on_a_seam_reads_the_side_its_station_names() {
    // The polyline turns left at (10, 0): the incoming segment runs +x
    // (lateral +y), the outgoing one +y (lateral -x).
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::ZERO,
                Point2::new(10.0, 0.0),
                Point2::new(10.0, 10.0),
            ],
            closed: false,
        }))
        .unwrap();
    let source = local_line(&mut b);
    let on_seam = at(basis, Station::at(10.0), StationFrame::Section);
    let default = b
        .push_value(InstanceAtStation::new(source, on_seam))
        .unwrap();
    let outgoing = b
        .push_value(InstanceAtStation::new(
            source,
            on_seam.with_seam_side(SeamSide::Outgoing),
        ))
        .unwrap();
    let incoming = b
        .push_value(InstanceAtStation::new(
            source,
            on_seam.with_seam_side(SeamSide::Incoming),
        ))
        .unwrap();
    let graph = b.finish(vec![default, outgoing, incoming]).unwrap();
    let point = Point3::new(10.0, 0.0, 0.0);
    let after = [Vec3::Y, -Vec3::X, Vec3::Z];
    let before = [Vec3::X, Vec3::Y, Vec3::Z];
    assert_placed(&placement(&graph, default).unwrap(), after, point, EPS);
    assert_placed(&placement(&graph, outgoing).unwrap(), after, point, EPS);
    assert_placed(&placement(&graph, incoming).unwrap(), before, point, EPS);
}

// --- refusals --------------------------------------------------------------------

#[test]
fn out_of_range_and_degenerate_stations_are_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let polyline = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::ZERO, Point2::new(20.0, 0.0)],
            closed: false,
        }))
        .unwrap();
    let vertical = b
        .push_value(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::Z,
        }))
        .unwrap();
    let zero = b
        .push_value(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::ZERO,
        }))
        .unwrap();
    let source = local_line(&mut b);
    let solid = block(&mut b);
    let place = |b: &mut GeometryGraphBuilder, source, basis, s, frame| {
        b.push_value(InstanceAtStation::new(
            source,
            at(basis, Station::at(s), frame),
        ))
        .unwrap()
    };
    let beyond = place(&mut b, source, polyline, 25.0, StationFrame::Section);
    let upright = place(&mut b, source, vertical, 1.0, StationFrame::Section);
    let upright_plan = place(&mut b, source, vertical, 1.0, StationFrame::Plan);
    let pointless = place(&mut b, source, zero, 1.0, StationFrame::Section);
    let beyond_solid = place(&mut b, solid, polyline, 25.0, StationFrame::Section);
    let tube = swept(&mut b, beyond, (0.0, 1.0));
    let graph = b
        .finish(vec![
            beyond,
            upright,
            upright_plan,
            pointless,
            beyond_solid,
            tube,
        ])
        .unwrap();
    refused_naming(placement(&graph, beyond), "beyond the curve's length");
    refused_naming(placement(&graph, upright), "vertical");
    refused_naming(placement(&graph, upright_plan), "vertical");
    refused_naming(placement(&graph, pointless), "no direction");
    refused_naming(placement(&graph, polyline), "not an instance at a station");
    // The compilers pass the refusal on.
    refused_naming(compile(&graph, beyond_solid), "beyond the curve's length");
    refused_naming(compile(&graph, tube), "beyond the curve's length");
    refused_naming(
        ReferenceExactCompiler::new().compile_exact(&graph, beyond_solid, &options()),
        "beyond the curve's length",
    );
}

// --- stations along a placed curve -------------------------------------------

#[test]
fn a_station_along_a_curve_placed_in_plan_is_the_source_station_carried() {
    // Placed at (2, 4, 0.5) in the plan line's upright frame: the placed
    // curve runs +y. Its station at 2 with 1 to the left: source point
    // (2, 0, 0) with t = x, l = y, u = Z, carried to (2, 6, 0.5) with
    // t = +y, l = -x, u = +Z, and moved left to (1, 6, 0.5).
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let source = local_line(&mut b);
    let placed = b
        .push_value(InstanceAtStation::new(
            source,
            at(
                basis,
                Station::new(3.0, StationOffsets::new(0.0, 0.5, 0.0)),
                StationFrame::Plan,
            ),
        ))
        .unwrap();
    let station = b
        .push_value(CurveStation::new(
            placed,
            Station::new(2.0, StationOffsets::new(1.0, 0.0, 0.0)),
        ))
        .unwrap();
    // A polyline placed the same way keeps its seams.
    let corner = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::ZERO, Point2::new(3.0, 0.0), Point2::new(3.0, 4.0)],
            closed: false,
        }))
        .unwrap();
    let placed_corner = b
        .push_value(InstanceAtStation::new(
            corner,
            at(basis, Station::at(3.0), StationFrame::Section),
        ))
        .unwrap();
    let on_corner = b
        .push_value(
            CurveStation::new(placed_corner, Station::at(3.0)).with_seam_side(SeamSide::Incoming),
        )
        .unwrap();
    let graph = b.finish(vec![station, on_corner]).unwrap();
    let resolved = resolve(&graph, station).unwrap();
    close3(resolved.point, Point3::new(1.0, 6.0, 0.5), EPS, "point");
    close3(resolved.frame.x, Vec3::Y, EPS, "x tangent");
    close3(resolved.frame.y, Vec3::Z, EPS, "y up");
    close3(resolved.frame.z, Vec3::X, EPS, "z right");

    let seams = seams(&graph, placed_corner, &options()).unwrap();
    assert_eq!(seams.len(), 1);
    assert!((seams[0].distance - 3.0).abs() <= EPS);
    // On that seam, the incoming side: the source's first segment (+x),
    // carried to +y, at (2, 4 + 3, 0).
    let corner_station = resolve(&graph, on_corner).unwrap();
    close3(
        corner_station.point,
        Point3::new(2.0, 7.0, 0.0),
        EPS,
        "corner",
    );
    close3(corner_station.frame.x, Vec3::Y, EPS, "incoming tangent");
}

#[test]
fn a_station_along_a_curve_placed_in_a_tilted_frame_is_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let grade = b
        .push_value(Curve3::Elevated(Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            ElevationLaw::constant_grade(100.0, 0.02),
        )))
        .unwrap();
    let basis = plan_line(&mut b);
    let source = local_line(&mut b);
    // On a grade in its section frame, and rolled on a level line.
    let leaning = b
        .push_value(InstanceAtStation::new(
            source,
            at(grade, Station::at(10.0), StationFrame::Section),
        ))
        .unwrap();
    let rolled = b
        .push_value(InstanceAtStation::new(
            source,
            OrientedCurveStation::new(
                CurveStation::new(basis, Station::at(1.0)),
                StationOrientation::new(Some(Vec3::new(0.0, 0.1, 1.0)), None),
            ),
        ))
        .unwrap();
    let on_leaning = b
        .push_value(CurveStation::new(leaning, Station::at(1.0)))
        .unwrap();
    let on_rolled = b
        .push_value(CurveStation::new(rolled, Station::at(1.0)))
        .unwrap();
    let graph = b.finish(vec![on_leaning, on_rolled]).unwrap();
    refused_naming(resolve(&graph, on_leaning), "tilts +Z");
    refused_naming(resolve(&graph, on_rolled), "tilts +Z");
    // The placement itself is fine.
    assert!(placement(&graph, leaning).is_ok());
}
