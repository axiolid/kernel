//! Stations on seams through the compiler (#263): the seam side a station
//! reads, section and offset runs mitred where they cross a corner, the
//! run's ends on a seam, and seam positions read from a graph's curves.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Frame3, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, Circle3, Curve2, Curve3, KnotSpec, Line2, Line3, Polyline2, Polyline3,
};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::{resolve, seams};
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveSegment, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode,
    NodeId, OpenProfile, OrientedCurveStation, SeamSide, SectionAtStation, SolidOperation, Station,
    StationFrame, StationOffsets, StationOrientation, SurfaceRelation, Transition, TrimSelector,
    TrimmingPreference,
};
use axiolid_profile::{Profile, RectangleProfile};

const EPS: Scalar = 1e-9;

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compile(graph: &GeometryGraph, root: NodeId) -> GeomResult<TriMesh> {
    compiler().compile_mesh(graph, root, &options())
}

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?}"
    );
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .expect("a stationed spine is closed")
        .signed_volume
}

fn has_vertex(mesh: &TriMesh, point: Point3) -> bool {
    mesh.positions
        .iter()
        .any(|p| (*p - point).abs().max_element() <= EPS)
}

/// An L in the plane: 10 m along +x, a right-angled left turn at
/// (10, 0, 0), 10 m along +y.
fn ell(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
        ],
        closed: false,
    }))
    .unwrap()
}

fn rectangle(b: &mut GeometryGraphBuilder, x: Scalar, y: Scalar) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
    .unwrap()
}

fn spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<SectionAtStation>,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

// --- a station on a seam ----------------------------------------------------

#[test]
fn a_station_on_a_corner_reads_the_side_it_names() {
    let mut b = GeometryGraphBuilder::new();
    let basis = ell(&mut b);
    let station = CurveStation::new(
        basis,
        Station::new(10.0, StationOffsets::new(1.0, 0.0, 0.0)),
    );
    let plain = b.push_value(station).unwrap();
    let incoming = b
        .push(GeometryNode::OrientedCurveStation(
            station.with_seam_side(SeamSide::Incoming),
        ))
        .unwrap();
    let outgoing = b
        .push(GeometryNode::OrientedCurveStation(
            OrientedCurveStation::new(station, StationOrientation::default())
                .with_seam_side(SeamSide::Outgoing),
        ))
        .unwrap();
    let graph = b.finish(vec![plain, incoming, outgoing]).unwrap();
    let plain = resolve(&graph, plain).unwrap();
    let incoming = resolve(&graph, incoming).unwrap();
    let outgoing = resolve(&graph, outgoing).unwrap();
    // One unit to the left of the first leg is +y, of the second -x.
    close3(
        incoming.point,
        Point3::new(10.0, 1.0, 0.0),
        EPS,
        "incoming point",
    );
    close3(incoming.frame.x, Vec3::X, EPS, "incoming tangent");
    close3(
        outgoing.point,
        Point3::new(9.0, 0.0, 0.0),
        EPS,
        "outgoing point",
    );
    close3(outgoing.frame.x, Vec3::Y, EPS, "outgoing tangent");
    assert_eq!(plain, outgoing, "a plain station reads the outgoing leg");
    assert_eq!(
        OrientedCurveStation::from(station).seam,
        SeamSide::Outgoing,
        "the default side"
    );
}

// --- mitred runs --------------------------------------------------------------

/// The four corners of a 2 x 1 rectangle, lateral `x` in [-1, 1] and up
/// `y` in [-0.5, 0.5], cut by the mitre plane of the L's corner: each
/// projected along +x from the first leg's section plane.
fn mitred_rectangle() -> [Point3; 4] {
    [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]
        .map(|(x, y): (Scalar, Scalar)| Point3::new(10.0 - x, x, y))
}

#[test]
fn a_spine_across_a_corner_is_mitred_in_the_bisector_plane() {
    let mut b = GeometryGraphBuilder::new();
    let basis = ell(&mut b);
    let profile = rectangle(&mut b, 2.0, 1.0);
    let root = spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(profile, Station::at(2.0)),
            SectionAtStation::new(profile, Station::at(18.0)),
        ],
    );
    let graph = b.finish(vec![root]).unwrap();
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .unwrap();
    let mesh = outcome.mesh;
    let normal = Vec3::new(1.0, 1.0, 0.0).normalize();
    for corner in mitred_rectangle() {
        assert!(has_vertex(&mesh, corner), "mitre corner {corner:?}");
        assert!(((corner - Point3::new(10.0, 0.0, 0.0)).dot(normal)).abs() <= EPS);
    }
    // Every vertex not on a straight leg's section lies in the plane.
    for p in &mesh.positions {
        let leg = (p.x - 2.0).abs() <= EPS || (p.y - 8.0).abs() <= EPS;
        let mitre = ((*p - Point3::new(10.0, 0.0, 0.0)).dot(normal)).abs() <= EPS;
        assert!(leg || mitre, "vertex {p:?} off every section plane");
    }
    // The mitred prism holds the section along the centreline, 8 + 8 m.
    let got = volume(&mesh);
    assert!((got - 2.0 * 16.0).abs() <= 1e-9, "volume {got}");
    // The run's deviation is still reported as sampled.
    assert!(report
        .contributions
        .iter()
        .any(|c| c.path == DeviationPath::StationedSpine
            && matches!(c.bound, DeviationBound::Unbounded(_))));
}

#[test]
fn a_section_standing_on_the_corner_is_mitred_too() {
    // Three stations, the middle one on the corner, widening as it goes:
    // the corner's own section is cut in the bisector plane, offsets and
    // all, and a longitudinal offset slides it along the bisector.
    let mut b = GeometryGraphBuilder::new();
    let basis = ell(&mut b);
    let profile = rectangle(&mut b, 2.0, 1.0);
    let lift = StationOffsets::new(0.0, 0.5, 0.0);
    let root = spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(profile, Station::at(2.0)),
            SectionAtStation::new(profile, Station::new(10.0, lift)),
            SectionAtStation::new(profile, Station::at(18.0)),
        ],
    );
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    for corner in mitred_rectangle() {
        assert!(
            has_vertex(&mesh, corner + Vec3::new(0.0, 0.0, 0.5)),
            "lifted mitre corner {corner:?}"
        );
    }
}

#[test]
fn an_open_section_run_across_a_corner_is_mitred() {
    let mut b = GeometryGraphBuilder::new();
    let basis = ell(&mut b);
    let path = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(-1.0, 0.0), Point2::new(1.0, 0.0)],
            closed: false,
        }))
        .unwrap();
    let deck = b.push_value(OpenProfile::new(path)).unwrap();
    let root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::OpenSectionsAtStations {
                directrix: basis,
                sections: vec![
                    SectionAtStation::new(deck, Station::at(5.0)),
                    SectionAtStation::new(deck, Station::at(15.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    // The deck's edges meet the mitre at the inner and outer corners.
    assert!(has_vertex(&mesh, Point3::new(9.0, 1.0, 0.0)), "inner edge");
    assert!(
        has_vertex(&mesh, Point3::new(11.0, -1.0, 0.0)),
        "outer edge"
    );
    // Two trapezoids of the 2 m deck: 5 m of centreline each.
    let area: Scalar = mesh
        .indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            0.5 * (b - a).cross(c - a).length()
        })
        .sum();
    assert!((area - 20.0).abs() <= 1e-9, "area {area}");
}

#[test]
fn an_offset_run_across_a_corner_takes_the_offset_corner() {
    // One metre left of the L is (0, 1) -> (9, 1) -> (9, 10); swept with
    // a small section, nothing strays to the basis's own corner side.
    let mut b = GeometryGraphBuilder::new();
    let basis = ell(&mut b);
    let offset = b
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetByStations {
                basis,
                stations: vec![
                    Station::new(0.0, StationOffsets::new(1.0, 0.0, 0.0)),
                    Station::new(20.0, StationOffsets::new(1.0, 0.0, 0.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let profile = rectangle(&mut b, 0.2, 0.4);
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::FixedReferenceSweep {
                profile,
                directrix: offset,
                reference_direction: Vec3::Z,
                parameter_range: None,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    for p in &mesh.positions {
        assert!(
            p.x <= 9.0 + 0.3 && p.y >= 1.0 - 0.3,
            "vertex {p:?} past the offset corner"
        );
    }
}

#[test]
fn a_straight_seam_is_sampled_as_before() {
    // A collinear vertex is a seam whose tangents agree: the spine is the
    // one along the plain segment, vertex for vertex.
    let mesh_along = |points: Vec<Point3>| {
        let mut b = GeometryGraphBuilder::new();
        let basis = b
            .push_value(Curve3::Polyline(Polyline3 {
                points,
                closed: false,
            }))
            .unwrap();
        let profile = rectangle(&mut b, 2.0, 1.0);
        let root = spine(
            &mut b,
            basis,
            vec![
                SectionAtStation::new(profile, Station::at(1.0)),
                SectionAtStation::new(profile, Station::at(9.0)),
            ],
        );
        let graph = b.finish(vec![root]).unwrap();
        compile(&graph, root).unwrap()
    };
    let straight = mesh_along(vec![Point3::ZERO, Point3::new(10.0, 0.0, 0.0)]);
    let seamed = mesh_along(vec![
        Point3::ZERO,
        Point3::new(4.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
    ]);
    assert_eq!(straight.indices, seamed.indices);
    for (a, b) in straight.positions.iter().zip(&seamed.positions) {
        close3(*a, *b, 1e-12, "same vertex");
    }
}

#[test]
fn a_run_ending_on_a_corner_ends_in_the_incoming_frame() {
    let mut b = GeometryGraphBuilder::new();
    let basis = ell(&mut b);
    let profile = rectangle(&mut b, 2.0, 1.0);
    let root = spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(profile, Station::at(10.0)),
            SectionAtStation::new(profile, Station::at(20.0)),
        ],
    );
    let ending = spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(profile, Station::at(0.0)),
            SectionAtStation::new(profile, Station::at(10.0)),
        ],
    );
    let graph = b.finish(vec![root, ending]).unwrap();
    // Starting on the corner: the outgoing leg's section, square to +y.
    let starting = compile(&graph, root).unwrap();
    assert!(has_vertex(&starting, Point3::new(11.0, 0.0, -0.5)));
    assert!(has_vertex(&starting, Point3::new(9.0, 0.0, 0.5)));
    // Ending on it: the incoming leg's, square to +x; a 2 x 1 x 10 box.
    let ending = compile(&graph, ending).unwrap();
    assert!(has_vertex(&ending, Point3::new(10.0, 1.0, 0.5)));
    assert!(has_vertex(&ending, Point3::new(10.0, -1.0, -0.5)));
    assert!((volume(&ending) - 20.0).abs() <= 1e-9);
}

#[test]
fn a_run_across_a_reversal_is_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve3::Polyline(Polyline3 {
            points: vec![
                Point3::ZERO,
                Point3::new(10.0, 0.0, 0.0),
                Point3::new(2.0, 1e-9, 0.0),
            ],
            closed: false,
        }))
        .unwrap();
    let profile = rectangle(&mut b, 2.0, 1.0);
    let root = spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(profile, Station::at(5.0)),
            SectionAtStation::new(profile, Station::at(15.0)),
        ],
    );
    let graph = b.finish(vec![root]).unwrap();
    let error = compile(&graph, root).unwrap_err();
    assert!(
        matches!(&error, GeomError::Degenerate(detail) if detail.contains("turns back on itself")),
        "{error:?}"
    );
}

#[test]
fn a_section_the_mitre_would_cut_is_refused_by_name() {
    // A 2 m wide section half a metre before the corner reaches past the
    // inner corner's mitre plane (x + y = 10 there).
    let mut b = GeometryGraphBuilder::new();
    let basis = ell(&mut b);
    let profile = rectangle(&mut b, 2.0, 1.0);
    let root = spine(
        &mut b,
        basis,
        vec![
            SectionAtStation::new(profile, Station::at(2.0)),
            SectionAtStation::new(profile, Station::at(9.5)),
            SectionAtStation::new(profile, Station::at(18.0)),
        ],
    );
    let graph = b.finish(vec![root]).unwrap();
    let error = compile(&graph, root).unwrap_err();
    assert!(
        error.to_string().contains("the mitre would cut it"),
        "{error:?}"
    );
}

// --- seam positions -----------------------------------------------------------

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

#[test]
fn a_composite_of_lines_and_an_arc_reports_its_joins() {
    // 10 m along +x, a quarter turn of radius 5 to the left, 3 m along +y.
    let mut b = GeometryGraphBuilder::new();
    let first = b
        .push_value(Curve3::Line(Line3 {
            origin: Point3::ZERO,
            direction: Vec3::X,
        }))
        .unwrap();
    let circle = b
        .push_value(Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(10.0, 5.0, 0.0),
                x: -Vec3::Y,
                y: Vec3::X,
                z: Vec3::Z,
            },
            radius: 5.0,
        }))
        .unwrap();
    let last = b
        .push_value(Curve3::Line(Line3 {
            origin: Point3::new(15.0, 5.0, 0.0),
            direction: Vec3::Y,
        }))
        .unwrap();
    let quarter = 0.5 * core::f64::consts::PI;
    let segments = [
        trimmed(&mut b, first, 0.0, 10.0),
        trimmed(&mut b, circle, 0.0, quarter),
        trimmed(&mut b, last, 0.0, 3.0),
    ]
    .map(|curve| CurveSegment {
        curve,
        same_sense: true,
        transition: Transition::ContinuousSameGradient,
    })
    .to_vec();
    let composite = b
        .push(GeometryNode::CurveRelation(CurveRelation::Composite {
            segments,
        }))
        .unwrap();
    let plan = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let plan_trim = trimmed(&mut b, plan, 0.0, 1.0);
    let spline = b
        .push_value(Curve3::BSpline(BSplineCurve {
            degree: 1,
            control_points: vec![Point3::ZERO, Point3::X, Point3::new(1.0, 1.0, 0.0)],
            knots: vec![0.0, 1.0, 2.0],
            multiplicities: vec![2, 1, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }))
        .unwrap();
    let polyline = ell(&mut b);
    let graph = b
        .finish(vec![composite, plan_trim, spline, polyline])
        .unwrap();
    let joins = seams(&graph, composite, &options()).unwrap();
    let at: Vec<Scalar> = joins.iter().map(|s| s.distance).collect();
    assert_eq!(at.len(), 2);
    assert!((at[0] - 10.0).abs() <= 1e-12);
    assert!((at[1] - (10.0 + 5.0 * quarter)).abs() <= 1e-12);
    assert!(joins.iter().all(|s| s.exact && !s.smooth));
    // An atomic polyline's vertex.
    let vertex = seams(&graph, polyline, &options()).unwrap();
    assert_eq!(vertex.len(), 1);
    assert_eq!((vertex[0].distance, vertex[0].parameter), (10.0, 1.0));
    // A 2D relation is read too (#285): a trimmed line has none.
    assert!(seams(&graph, plan_trim, &options()).unwrap().is_empty());
    // Refused by name: a corner only quadrature locates.
    let error = seams(&graph, spline, &options()).unwrap_err();
    assert!(
        matches!(&error, GeomError::UnsupportedInput { input, .. }
            if input.contains("its arc length to the knot is a quadrature")),
        "{error:?}"
    );
}
