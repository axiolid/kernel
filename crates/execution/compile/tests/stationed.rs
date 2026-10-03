//! Stations through the compiler (#241): resolved points and frames, an
//! offset curve swept, station-placed spines and a sectioned surface, each
//! against a closed form; refusals by name.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Circle2, CurvatureLaw, Curve2, Curve3, Elevated3,
    ElevationLaw, Intrinsic2, Line2, Polyline2,
};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::resolve;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_mesh_compile_contract::{MeshClosure, MeshCompiler};
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId,
    OpenProfile, SolidOperation, Station, StationFrame, StationOffsets, StationedOpenSection,
    StationedSection, SurfaceRelation,
};
use axiolid_profile::{Profile, RectangleProfile};

const EPS: Scalar = 1e-9;

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?}"
    );
}

fn close(actual: Scalar, expected: Scalar, relative: Scalar, what: &str) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{what}: {actual} != {expected} (relative {:e})",
        (actual - expected).abs() / expected.abs().max(1.0)
    );
}

fn plan_line(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::new(2.0, 1.0),
        direction: Vec2::new(0.0, 2.0),
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

fn two_percent() -> Elevated3 {
    Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::constant_grade(100.0, 0.02),
    )
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .expect("a stationed spine is closed")
        .signed_volume
}

fn spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<StationedSection>,
    frame: StationFrame,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::StationedSpine {
            directrix,
            sections,
            frame,
        },
    ))
    .unwrap()
}

fn compile(graph: &GeometryGraph, root: NodeId) -> GeomResult<TriMesh> {
    compiler().compile_mesh(graph, root, &options())
}

// --- resolving a station ---------------------------------------------------

#[test]
fn a_station_on_a_plan_line_resolves_with_its_offsets() {
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    // Along +y from (2, 1): left is -x.
    let id = b
        .push_value(CurveStation::new(
            basis,
            Station::new(4.0, StationOffsets::new(1.5, 0.25, 0.5)),
        ))
        .unwrap();
    let graph = b.finish(vec![id]).unwrap();
    let resolved = resolve(&graph, id).unwrap();
    close3(
        resolved.point,
        Point3::new(2.0 - 1.5, 1.0 + 4.0 + 0.5, 0.25),
        EPS,
        "point",
    );
    close3(resolved.frame.origin, resolved.point, 0.0, "origin");
    close3(resolved.frame.x, Vec3::Y, EPS, "x tangent");
    close3(resolved.frame.y, Vec3::Z, EPS, "y up");
    close3(resolved.frame.z, Vec3::X, EPS, "z right");
}

#[test]
fn an_elevated_station_offsets_in_its_leaning_or_upright_frame() {
    let k = 1.0004_f64.sqrt();
    for (frame, up) in [
        (StationFrame::Section, Vec3::new(-0.02 / k, 0.0, 1.0 / k)),
        (StationFrame::Plan, Vec3::Z),
    ] {
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push_value(Curve3::Elevated(two_percent())).unwrap();
        let id = b
            .push_value(CurveStation {
                basis,
                station: Station::new(50.0, StationOffsets::new(-2.0, 1.0, 0.0)),
                frame,
            })
            .unwrap();
        let graph = b.finish(vec![id]).unwrap();
        let resolved = resolve(&graph, id).unwrap();
        close3(
            resolved.point,
            Point3::new(50.0, -2.0, 101.0) + up,
            EPS,
            &format!("{frame:?}"),
        );
    }
}

#[test]
fn a_banked_station_rolls_its_lateral_offset() {
    // Constant 150 mm cant on a level straight: the left rail head, b / 2
    // left, stands D / 2 above the centreline.
    let mut base = two_percent();
    base.elevation = ElevationLaw::level(0.0);
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve3::Banked(Banked3::new(
            base,
            CantLaw::new(vec![CantPiece::constant(100.0, 0.15)]),
            CantLaw::zero(100.0),
            1.5,
            BankConvention::TangentRotation,
        )))
        .unwrap();
    let id = b
        .push_value(CurveStation::new(
            basis,
            Station::new(30.0, StationOffsets::new(0.75, 0.0, 0.0)),
        ))
        .unwrap();
    let graph = b.finish(vec![id]).unwrap();
    let resolved = resolve(&graph, id).unwrap();
    let cos_psi = (1.0_f64 - 0.01).sqrt();
    close3(
        resolved.point,
        Point3::new(30.0, 0.75 * cos_psi, 0.075),
        EPS,
        "left rail head",
    );
}

#[test]
fn a_station_beyond_the_curve_or_on_a_relation_is_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let clothoid = b
        .push_value(Curve2::Intrinsic(Intrinsic2::new(
            Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
            60.0,
        )))
        .unwrap();
    let beyond = b
        .push_value(CurveStation::new(clothoid, Station::at(61.0)))
        .unwrap();
    let relation = b
        .push(GeometryNode::CurveRelation(CurveRelation::Offset {
            basis: clothoid,
            distance: 1.0,
            reference_direction: None,
        }))
        .unwrap();
    let on_relation = b
        .push_value(CurveStation::new(relation, Station::at(1.0)))
        .unwrap();
    let graph = b.finish(vec![beyond, on_relation]).unwrap();
    let error = resolve(&graph, beyond).unwrap_err();
    assert!(
        error.to_string().contains("beyond the curve's length"),
        "{error}"
    );
    assert!(matches!(
        resolve(&graph, on_relation).unwrap_err(),
        GeomError::UnsupportedInput {
            operation: Operation::CurveEvaluation,
            input: "a station along an instanced curve or a curve relation",
            ..
        }
    ));
    assert!(resolve(&graph, clothoid).is_err(), "not a station node");
}

// --- offset curve by stations -----------------------------------------------

fn swept(b: &mut GeometryGraphBuilder, directrix: NodeId) -> NodeId {
    let profile = rectangle(b, 0.2, 0.4);
    b.push(GeometryNode::SolidOperation(
        SolidOperation::FixedReferenceSweep {
            profile,
            directrix,
            reference_direction: Vec3::Z,
            parameter_range: None,
        },
    ))
    .unwrap()
}

#[test]
fn an_offset_curve_widening_linearly_is_a_straight_directrix() {
    // Lateral 1 -> 3 between distances 2 and 12 along a +y line: the
    // offset curve runs from (1, 3, 0) to (-1, 13, 0.5).
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let offset = b
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetByStations {
                basis,
                stations: vec![
                    Station::new(2.0, StationOffsets::new(1.0, 0.0, 0.0)),
                    Station::new(12.0, StationOffsets::new(3.0, 0.5, 0.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let root = swept(&mut b, offset);
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    let length = (4.0_f64 + 100.0 + 0.25).sqrt();
    close(volume(&mesh).abs(), 0.08 * length, 1e-9, "prism volume");
    let lowest = mesh
        .positions
        .iter()
        .map(|p| p.y)
        .fold(Scalar::INFINITY, Scalar::min);
    assert!(lowest > 2.8 && lowest < 3.2, "starts at the first station");
}

#[test]
fn a_constant_offset_of_an_arc_is_the_concentric_arc() {
    // 2 m left of a counter-clockwise R = 30 arc is R = 28; a quarter turn
    // of the basis is a quarter turn of the offset, swept by Pappus.
    let mut b = GeometryGraphBuilder::new();
    let radius = 30.0;
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
    let quarter = 0.5 * core::f64::consts::PI * radius;
    let offset = b
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetByStations {
                basis,
                stations: vec![
                    Station::new(0.0, StationOffsets::new(2.0, 0.0, 0.0)),
                    Station::new(quarter, StationOffsets::new(2.0, 0.0, 0.0)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let root = swept(&mut b, offset);
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    let expected = 0.08 * 0.5 * core::f64::consts::PI * 28.0;
    close(volume(&mesh).abs(), expected, 1e-4, "Pappus");
    for p in &mesh.positions {
        let r = p.x.hypot(p.y);
        assert!((r - 28.0).abs() <= 0.2 + 1e-6, "point at radius {r}");
    }
}

// --- station-placed spine ---------------------------------------------------

#[test]
fn a_constant_section_along_a_line_is_a_prism() {
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let profile = rectangle(&mut b, 1.0, 2.0);
    let at = |d: Scalar| StationedSection {
        profile,
        station: Station::new(d, StationOffsets::new(0.5, 1.0, 0.0)),
    };
    let root = spine(
        &mut b,
        basis,
        vec![at(1.0), at(4.0), at(11.0)],
        StationFrame::Section,
    );
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    close(volume(&mesh), 20.0, 1e-12, "A * L, outward");
    // A straight directrix needs no sections between the stations.
    assert_eq!(mesh.positions.len(), 3 * 4, "one ring per station");
    // Sections: lateral offset 0.5 left of +y is x = 2 - 0.5, the profile's
    // x spans +-0.5 about it; up 1 lifts y's +-1 to [0, 2].
    let (lo, hi) = bounds(&mesh);
    close3(lo, Point3::new(1.0, 2.0, 0.0), EPS, "min corner");
    close3(hi, Point3::new(2.0, 12.0, 2.0), EPS, "max corner");
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

#[test]
fn a_constant_section_along_an_arc_follows_pappus() {
    let mut b = GeometryGraphBuilder::new();
    let radius = 25.0;
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
    let profile = rectangle(&mut b, 2.0, 1.0);
    let length = 0.5 * core::f64::consts::PI * radius;
    let root = spine(
        &mut b,
        basis,
        vec![
            StationedSection {
                profile,
                station: Station::at(0.0),
            },
            StationedSection {
                profile,
                station: Station::at(length),
            },
        ],
        StationFrame::Section,
    );
    let graph = b.finish(vec![root]).unwrap();
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .unwrap();
    // The centroid rides the arc, so the exact volume is A * L.
    close(volume(&outcome.mesh), 2.0 * length, 1e-4, "Pappus");
    assert_eq!(outcome.closure, MeshClosure::Solid);
    assert!(report.bound.is_none());
    assert!(report
        .contributions
        .iter()
        .any(|c| c.path == DeviationPath::StationedSpine
            && matches!(c.bound, DeviationBound::Unbounded(reason) if reason.contains("sampled"))));
    // Every vertex keeps its section's radius band: [24, 26].
    for p in &outcome.mesh.positions {
        let r = p.x.hypot(p.y);
        assert!((24.0 - 1e-9..=26.0 + 1e-9).contains(&r), "radius {r}");
    }
}

#[test]
fn sections_are_interpolated_linearly_between_stations() {
    // 1 x 1 at distance 0 to 1 x 3 at distance 10: the area runs 1 -> 3
    // linearly, so the volume is 20; the walls are planes.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let small = rectangle(&mut b, 1.0, 1.0);
    let tall = rectangle(&mut b, 1.0, 3.0);
    let root = spine(
        &mut b,
        basis,
        vec![
            StationedSection {
                profile: small,
                station: Station::at(0.0),
            },
            StationedSection {
                profile: tall,
                station: Station::new(10.0, StationOffsets::new(0.0, 1.0, 0.0)),
            },
        ],
        StationFrame::Section,
    );
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    close(volume(&mesh), 20.0, 1e-12, "linear taper");
    assert_eq!(
        mesh.positions.len(),
        2 * 4,
        "planar walls need no refinement"
    );
    let (lo, hi) = bounds(&mesh);
    close3(lo, Point3::new(1.5, 1.0, -0.5), EPS, "start bottom");
    close3(hi, Point3::new(2.5, 11.0, 2.5), EPS, "end top");
}

#[test]
fn on_a_grade_section_frames_tilt_and_plan_frames_stand() {
    // Plan distance 10 at 2%: tilted sections sweep a right prism along
    // the 3D line, upright ones a sheared prism of plan length 10.
    for (frame, expected) in [
        (StationFrame::Section, 2.0 * 10.0 * 1.0004_f64.sqrt()),
        (StationFrame::Plan, 2.0 * 10.0),
    ] {
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push_value(Curve3::Elevated(two_percent())).unwrap();
        let profile = rectangle(&mut b, 1.0, 2.0);
        let at = |d: Scalar| StationedSection {
            profile,
            station: Station::at(d),
        };
        let root = spine(&mut b, basis, vec![at(5.0), at(15.0)], frame);
        let graph = b.finish(vec![root]).unwrap();
        let mesh = compile(&graph, root).unwrap();
        close(volume(&mesh), expected, 1e-12, &format!("{frame:?}"));
    }
}

#[test]
fn a_banked_spine_compiles_closed() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve3::Banked(Banked3::new(
            two_percent(),
            CantLaw::new(vec![CantPiece::linear(40.0, 0.0, 0.12)]),
            CantLaw::zero(40.0),
            1.5,
            BankConvention::TangentRotation,
        )))
        .unwrap();
    let profile = rectangle(&mut b, 3.0, 0.5);
    let root = spine(
        &mut b,
        basis,
        vec![
            StationedSection {
                profile,
                station: Station::at(0.0),
            },
            StationedSection {
                profile,
                station: Station::at(40.0),
            },
        ],
        StationFrame::Section,
    );
    let graph = b.finish(vec![root]).unwrap();
    let mesh = compile(&graph, root).unwrap();
    // A straight centreline: the roll twists the walls but keeps the
    // section area along the 3D length to within the twist.
    close(
        volume(&mesh),
        1.5 * 40.0 * 1.0004_f64.sqrt(),
        1e-3,
        "banked",
    );
}

#[test]
fn malformed_spines_are_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::ZERO, Point2::new(10.0, 0.0)],
            closed: false,
        }))
        .unwrap();
    let rect = rectangle(&mut b, 1.0, 1.0);
    let disk = b
        .push_value(Profile::Circle(axiolid_profile::CircleProfile {
            radius: 0.5,
            thickness: None,
        }))
        .unwrap();
    let mismatched = spine(
        &mut b,
        basis,
        vec![
            StationedSection {
                profile: rect,
                station: Station::at(0.0),
            },
            StationedSection {
                profile: disk,
                station: Station::at(10.0),
            },
        ],
        StationFrame::Section,
    );
    let beyond = spine(
        &mut b,
        basis,
        vec![
            StationedSection {
                profile: rect,
                station: Station::at(0.0),
            },
            StationedSection {
                profile: rect,
                station: Station::at(10.5),
            },
        ],
        StationFrame::Section,
    );
    let graph = b.finish(vec![mismatched, beyond]).unwrap();
    let error = compile(&graph, mismatched).unwrap_err();
    assert!(error.to_string().contains("ring structure"), "{error}");
    let error = compile(&graph, beyond).unwrap_err();
    assert!(
        error.to_string().contains("beyond the curve's length"),
        "{error}"
    );
    // The exact compiler refuses the family by name.
    let error = ReferenceExactCompiler::new()
        .compile_exact(&graph, beyond, &options())
        .unwrap_err();
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                input: "station-placed spine",
                ..
            }
        ),
        "{error:?}"
    );
}

// --- sectioned surface -------------------------------------------------------

fn open_section(b: &mut GeometryGraphBuilder, points: Vec<Point2>) -> NodeId {
    let path = b
        .push_value(Curve2::Polyline(Polyline2 {
            points,
            closed: false,
        }))
        .unwrap();
    b.push_value(OpenProfile::new(path)).unwrap()
}

fn area(mesh: &TriMesh) -> Scalar {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            0.5 * (b - a).cross(c - a).length()
        })
        .sum()
}

#[test]
fn a_sectioned_surface_joins_tagged_points_into_a_sheet() {
    // A 4 m crowned deck section (2% cross-fall each side) at distance 0,
    // widened to 6 m at distance 10: tags pin the edges and the crown.
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let narrow = open_section(
        &mut b,
        vec![
            Point2::new(-2.0, -0.04),
            Point2::new(0.0, 0.0),
            Point2::new(2.0, -0.04),
        ],
    );
    let wide = open_section(
        &mut b,
        vec![
            Point2::new(-3.0, -0.06),
            Point2::new(0.0, 0.0),
            Point2::new(3.0, -0.06),
        ],
    );
    let tags = || vec!["left".to_owned(), "crown".to_owned(), "right".to_owned()];
    let root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::SectionedSurface {
                directrix: basis,
                sections: vec![
                    StationedOpenSection {
                        profile: narrow,
                        tags: tags(),
                        station: Station::at(0.0),
                    },
                    StationedOpenSection {
                        profile: wide,
                        tags: tags(),
                        station: Station::at(10.0),
                    },
                ],
                frame: StationFrame::Plan,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!(outcome.solid_mesh().is_err(), "a sheet has no volume");
    assert!(report
        .contributions
        .iter()
        .any(|c| c.path == DeviationPath::SectionedSurface));
    // Each half is a planar trapezoid: widths 2 -> 3 (slope-corrected by
    // sqrt(1 + 0.02^2)) over 10 m.
    let expected = 2.0 * 0.5 * (2.0 + 3.0) * 10.0 * (1.0_f64 + 0.0004).sqrt();
    close(area(&outcome.mesh), expected, 1e-12, "deck area");
    let (lo, hi) = bounds(&outcome.mesh);
    close3(lo, Point3::new(-1.0, 1.0, -0.06), EPS, "min");
    close3(hi, Point3::new(5.0, 11.0, 0.0), EPS, "max");
}

#[test]
fn a_sectioned_surface_refuses_a_curved_or_mistagged_section() {
    let mut b = GeometryGraphBuilder::new();
    let basis = plan_line(&mut b);
    let three = open_section(
        &mut b,
        vec![Point2::new(-1.0, 0.0), Point2::ZERO, Point2::new(1.0, 0.0)],
    );
    let two = open_section(&mut b, vec![Point2::new(-1.0, 0.0), Point2::new(1.0, 0.0)]);
    let surface = |b: &mut GeometryGraphBuilder, first, second, tags: Vec<String>| {
        b.push(GeometryNode::SurfaceRelation(
            SurfaceRelation::SectionedSurface {
                directrix: basis,
                sections: vec![
                    StationedOpenSection {
                        profile: first,
                        tags: tags.clone(),
                        station: Station::at(0.0),
                    },
                    StationedOpenSection {
                        profile: second,
                        tags,
                        station: Station::at(5.0),
                    },
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap()
    };
    let counts = surface(&mut b, three, two, Vec::new());
    let tag_count = surface(&mut b, three, three, vec!["a".into(), "b".into()]);
    let graph = b.finish(vec![counts, tag_count]).unwrap();
    let error = compile(&graph, counts).unwrap_err();
    assert!(error.to_string().contains("same count"), "{error}");
    let error = compile(&graph, tag_count).unwrap_err();
    assert!(
        error.to_string().contains("2 tags for 3 vertices"),
        "{error}"
    );
}
