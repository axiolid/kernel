//! Run stations at their basis's start or end, measured by the kernel
//! (#303): an offset by distances continued to its basis's ends builds
//! without any stated length, along bases whose length needs quadrature,
//! reads exactly as the run stated at the kernel's own length of the basis
//! (`curve_path` measured by `CompositeBasis`), and is refused by name
//! where its resolved end does not come after its last stated station.
//! Snapping windows apply after resolution.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2};
use axiolid_curve::{
    BSplineCurve, CurvatureLaw, Curve2, Curve3, CurvePath, Elevated3, ElevationLaw, Intrinsic2,
    KnotSpec, Line2, OffsetLaw, PathCurve, PathOffsets,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::curve_path;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId,
    RunStation, SeamSide, SolidOperation, Station, StationFrame, StationOffsets, Transition,
    TrimSelector, TrimmingPreference,
};
use axiolid_reference::station::{CompositeBasis, SEAM_WINDOW_AMBIGUOUS};

/// The offset held from the last stated station to the basis's end.
const LAST: StationOffsets = StationOffsets::new(1.25, 0.2, 0.0);
/// The offset held from the basis's start to the first stated station.
const FIRST: StationOffsets = StationOffsets::new(0.75, -0.1, 0.0);
/// Where the last stated station stands, short of every basis's end.
const SHORT: Scalar = 20.0;

/// A base a run is built along, by name, and whether the run is flattened
/// into a curve path (an offset by stations beside an offset is swept, but
/// as a station basis or a path it is an offset of an offset, which is
/// refused by name, #289).
type Base = (&'static str, fn(&mut GeometryGraphBuilder) -> NodeId, bool);

/// An offset by distances along a road on a crest (a gradient curve):
/// its own plan length is a quadrature.
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

fn offset(b: &mut GeometryGraphBuilder, basis: NodeId, distance: Scalar) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Offset {
        basis,
        distance,
        reference_direction: None,
    }))
    .unwrap()
}

fn clothoid(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
        60.0,
    )))
    .unwrap()
}

/// A constant offset beside a clothoid: no closed-form length.
fn beside_clothoid(b: &mut GeometryGraphBuilder) -> NodeId {
    let basis = clothoid(b);
    offset(b, basis, 2.0)
}

fn spline(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::BSpline(BSplineCurve {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(30.0, 0.0),
            Point2::new(60.0, 20.0),
            Point2::new(90.0, 20.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    }))
    .unwrap()
}

/// A constant offset beside a B-spline: its length is a quadrature.
fn beside_spline(b: &mut GeometryGraphBuilder) -> NodeId {
    let basis = spline(b);
    offset(b, basis, -2.0)
}

/// A gradient curve whose plan is a clothoid: plan-measured.
fn gradient(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve3::Elevated(Elevated3::new(
        Curve2::Intrinsic(Intrinsic2::new(
            Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
            60.0,
        )),
        ElevationLaw::Polynomial {
            coefficients: vec![10.0, 0.02, -0.0002],
        },
    )))
    .unwrap()
}

const BASES: [Base; 6] = [
    ("beside a gradient curve", kerb, false),
    ("a gradient curve", gradient, true),
    ("beside a clothoid", beside_clothoid, false),
    ("beside a B-spline", beside_spline, false),
    ("a clothoid", clothoid, true),
    ("a B-spline", spline, true),
];

/// A graph holding `base`, an offset by distances at `stations` along it
/// and a small disk swept along that offset: `(graph, base, run, sweep)`.
fn built(
    base: fn(&mut GeometryGraphBuilder) -> NodeId,
    stations: Vec<RunStation>,
) -> (GeometryGraph, NodeId, NodeId, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let basis = base(&mut b);
    let run = b
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetsAtStations {
                basis,
                stations,
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let swept = b
        .push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix: run,
            radius: RADIUS,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        }))
        .unwrap();
    let graph = b.finish(vec![run, swept]).unwrap();
    (graph, basis, run, swept)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> GeomResult<TriMesh> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(
        graph,
        root,
        &ExecutionOptions::new(Tolerance::MILLIMETRE),
    )
}

/// Two meshes vertex for vertex, bitwise.
fn same(a: &TriMesh, b: &TriMesh, what: &str) {
    assert_eq!(a.indices, b.indices, "{what}: indices");
    assert_eq!(a.positions, b.positions, "{what}: positions");
}

/// The swept disk's radius.
const RADIUS: Scalar = 0.05;

/// A disk swept along a run ends in a ring of vertices at `RADIUS` around
/// the run's end point, centred on it: the ring is found and its centroid
/// checked against `point`.
fn ring_around(mesh: &TriMesh, point: Point3, what: &str) {
    let ring: Vec<Point3> = mesh
        .positions
        .iter()
        .copied()
        .filter(|p| ((*p - point).length() - RADIUS).abs() <= 1e-6)
        .collect();
    assert!(
        ring.len() >= 8,
        "{what}: no ring of the sweep around {point:?}"
    );
    let centroid = ring
        .iter()
        .fold(Point3::ZERO, |sum, p| sum + (*p - Point3::ZERO))
        / ring.len() as Scalar;
    close(centroid, point, 1e-6, what);
}

fn close(actual: Point3, expected: Point3, eps: Scalar, what: &str) {
    let off = (actual - expected).length();
    assert!(
        off <= eps,
        "{what}: {actual:?} != {expected:?} (off by {off:e})"
    );
}

/// The basis's length by the kernel's own measure: its `curve_path`,
/// measured as a provider measures it.
fn kernel_length(graph: &GeometryGraph, basis: NodeId) -> Scalar {
    CompositeBasis::from_path(&curve_path(graph, basis).unwrap())
        .unwrap()
        .length()
}

fn stated(distance: Scalar, offsets: StationOffsets) -> RunStation {
    RunStation::new(Station::new(distance, offsets))
}

/// The law of an offset piece of a run's path.
fn law(path: &CurvePath, index: usize) -> OffsetLaw {
    match &path.pieces()[index].curve {
        PathCurve::Offset(offset) => offset.law,
        other => panic!("piece {index} is not an offset: {other:?}"),
    }
}

fn held(offsets: StationOffsets) -> OffsetLaw {
    let at = PathOffsets::new(offsets.lateral, offsets.vertical, offsets.longitudinal);
    OffsetLaw::Linear {
        start: at,
        end: at,
        frame: axiolid_curve::OffsetFrame::Section,
    }
}

/// The point at `distance` along `basis`, read from `side`, moved by
/// `offsets` in its section frame: what a station there resolves to.
fn point_on(
    graph: &GeometryGraph,
    basis: NodeId,
    distance: Scalar,
    side: SeamSide,
    offsets: StationOffsets,
) -> Point3 {
    CompositeBasis::from_path(&curve_path(graph, basis).unwrap())
        .unwrap()
        .section_on(distance, side)
        .unwrap()
        .place(offsets.lateral, offsets.vertical, offsets.longitudinal)
}

#[test]
fn an_offset_continued_to_the_basis_s_end_reads_as_stated_at_its_kernel_length() {
    for (name, base, flat) in BASES {
        let open = vec![
            stated(0.0, FIRST),
            stated(SHORT, LAST),
            RunStation::at_basis_end(LAST),
        ];
        let (graph, basis, run, swept) = built(base, open);
        let length = kernel_length(&graph, basis);
        assert!(length > SHORT + 10.0, "{name}: {length}");
        // The run stated at the kernel's own length is the same curve,
        // bitwise: its path and its swept mesh.
        let (twin, _, twin_run, twin_swept) = built(
            base,
            vec![
                stated(0.0, FIRST),
                stated(SHORT, LAST),
                stated(length, LAST),
            ],
        );
        let mesh = compile(&graph, swept).unwrap();
        same(&mesh, &compile(&twin, twin_swept).unwrap(), name);
        // The sweep ends at the basis's end, the last offset held.
        ring_around(
            &mesh,
            point_on(&graph, basis, length, SeamSide::Incoming, LAST),
            &format!("{name}: the sweep's end"),
        );
        if !flat {
            continue;
        }
        let path = curve_path(&graph, run).unwrap();
        assert_eq!(path, curve_path(&twin, twin_run).unwrap(), "{name}: path");
        // Its last piece holds the last offset to the basis's end, and its
        // end lies there.
        let last = path.pieces().len() - 1;
        assert_eq!(law(&path, last), held(LAST), "{name}: the last law");
        let composite = CompositeBasis::from_path(&path).unwrap();
        let end = composite
            .section_on(composite.length(), SeamSide::Incoming)
            .unwrap()
            .point;
        close(
            end,
            point_on(&graph, basis, length, SeamSide::Incoming, LAST),
            1e-6,
            &format!("{name}: the end"),
        );
    }
}

#[test]
fn an_offset_continued_from_the_basis_s_start_reads_as_stated_at_zero() {
    for (name, base, flat) in BASES {
        let open = vec![
            RunStation::at_basis_start(FIRST),
            stated(5.0, FIRST),
            stated(SHORT, LAST),
            RunStation::at_basis_end(LAST),
        ];
        let (graph, basis, run, swept) = built(base, open);
        let length = kernel_length(&graph, basis);
        let (twin, _, twin_run, twin_swept) = built(
            base,
            vec![
                stated(0.0, FIRST),
                stated(5.0, FIRST),
                stated(SHORT, LAST),
                stated(length, LAST),
            ],
        );
        let mesh = compile(&graph, swept).unwrap();
        same(&mesh, &compile(&twin, twin_swept).unwrap(), name);
        ring_around(
            &mesh,
            point_on(&graph, basis, 0.0, SeamSide::Outgoing, FIRST),
            &format!("{name}: the sweep's start"),
        );
        if !flat {
            continue;
        }
        let path = curve_path(&graph, run).unwrap();
        assert_eq!(path, curve_path(&twin, twin_run).unwrap(), "{name}: path");
        assert_eq!(law(&path, 0), held(FIRST), "{name}: the first law");
        let composite = CompositeBasis::from_path(&path).unwrap();
        close(
            composite.section_on(0.0, SeamSide::Outgoing).unwrap().point,
            point_on(&graph, basis, 0.0, SeamSide::Outgoing, FIRST),
            1e-9,
            &format!("{name}: the start"),
        );
    }
}

fn refused_out_of_order<T: core::fmt::Debug>(result: GeomResult<T>, what: &str) {
    match result {
        Err(GeomError::InvalidInput(detail)) => assert!(
            detail.contains("at or past the basis's end"),
            "{what}: {detail}"
        ),
        other => panic!("{what}: expected the order refusal, got {other:?}"),
    }
}

#[test]
fn a_run_whose_resolved_end_is_not_after_its_last_station_is_refused_by_name() {
    for (name, base, _) in BASES {
        let (graph, basis, _, _) = built(base, vec![stated(0.0, FIRST), stated(1.0, LAST)]);
        let length = kernel_length(&graph, basis);
        // Stated past the end, and exactly at it: the end comes no later.
        for last in [length + 1.0, length] {
            let (graph, _, run, swept) = built(
                base,
                vec![
                    stated(0.0, FIRST),
                    stated(last, LAST),
                    RunStation::at_basis_end(LAST),
                ],
            );
            let what = format!("{name}: last stated at {last}");
            refused_out_of_order(curve_path(&graph, run), &what);
            refused_out_of_order(compile(&graph, swept), &what);
        }
    }
}

/// Three collinear lines: 10 m, a sliver of `tiny`, the rest to `rest`.
fn sliver(b: &mut GeometryGraphBuilder, tiny: Scalar, rest: Scalar) -> NodeId {
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let mut trim = |start: Scalar, end: Scalar| {
        b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis: line,
            start: vec![TrimSelector::Parameter(start)],
            end: vec![TrimSelector::Parameter(end)],
            sense_agreement: true,
            preference: TrimmingPreference::Parameter,
        }))
        .unwrap()
    };
    let pieces = [
        trim(0.0, 10.0),
        trim(10.0, 10.0 + tiny),
        trim(10.0 + tiny, rest),
    ];
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

/// A station at the end snaps from the resolved end (#301 after #303):
/// within its window of the last joint, it stands on that joint.
#[test]
fn a_windowed_end_station_snaps_after_resolution() {
    // The basis ends 0.6 mm after its joint at 10.0005.
    fn near_end(b: &mut GeometryGraphBuilder) -> NodeId {
        sliver(b, 5e-4, 10.0011)
    }
    let window = 1e-3;
    let (graph, _, run, _) = built(
        near_end,
        vec![
            stated(2.0, FIRST),
            RunStation::at_basis_end(LAST).with_seam_window(window),
        ],
    );
    let (on, _, on_run, _) = built(near_end, vec![stated(2.0, FIRST), stated(10.0005, LAST)]);
    assert_eq!(
        curve_path(&graph, run).unwrap(),
        curve_path(&on, on_run).unwrap(),
        "snapped onto the joint"
    );
    // Without a window it is the basis's end.
    let (graph, basis, run, _) = built(
        near_end,
        vec![stated(2.0, FIRST), RunStation::at_basis_end(LAST)],
    );
    let (end, _, end_run, _) = built(
        near_end,
        vec![
            stated(2.0, FIRST),
            stated(kernel_length(&graph, basis), LAST),
        ],
    );
    assert_eq!(
        curve_path(&graph, run).unwrap(),
        curve_path(&end, end_run).unwrap(),
        "at the end"
    );
    assert_ne!(
        curve_path(&graph, run).unwrap(),
        curve_path(&on, on_run).unwrap(),
        "the window moves the end"
    );
    // A window holding two joints around the resolved end is refused as
    // any such window is.
    let (graph, _, run, _) = built(
        near_end,
        vec![
            stated(2.0, FIRST),
            RunStation::at_basis_end(LAST).with_seam_window(2e-3),
        ],
    );
    match curve_path(&graph, run) {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert_eq!(input, SEAM_WINDOW_AMBIGUOUS);
        }
        other => panic!("expected the ambiguous window refused, got {other:?}"),
    }
}

/// Along an unbounded line the basis's end is its parameter domain's,
/// `|direction|`, as inside a composite and as a path's open end (#298).
#[test]
fn along_a_line_the_basis_s_end_is_its_parameter_domain_s() {
    fn line(b: &mut GeometryGraphBuilder) -> NodeId {
        b.push_value(Curve2::Line(Line2 {
            origin: Point2::new(1.0, 2.0),
            direction: Vec2::new(3.0, 4.0),
        }))
        .unwrap()
    }
    let (graph, _, run, swept) = built(
        line,
        vec![
            RunStation::at_basis_start(FIRST),
            stated(1.0, LAST),
            RunStation::at_basis_end(LAST),
        ],
    );
    let (twin, _, twin_run, twin_swept) = built(
        line,
        vec![stated(0.0, FIRST), stated(1.0, LAST), stated(5.0, LAST)],
    );
    assert_eq!(
        curve_path(&graph, run).unwrap(),
        curve_path(&twin, twin_run).unwrap()
    );
    same(
        &compile(&graph, swept).unwrap(),
        &compile(&twin, twin_swept).unwrap(),
        "along a line",
    );
}
