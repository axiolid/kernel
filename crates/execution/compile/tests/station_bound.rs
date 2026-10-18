//! Certified bounds and exact bodies of runs of sections at stations
//! (#247): spines along a line, an arc, a grade and a banked curve, a
//! sectioned deck, and a sectioned surface report a certified bound within
//! the budget that a dense sampling of the exact surface never exceeds;
//! the closed-form runs compile exactly with their closed-form volume and
//! area; everything else is unbounded or refused by name.

use axiolid_construct::profile::profile_rings;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Transform2, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Circle2, CurvatureLaw, Curve2, Curve3, Elevated3,
    ElevationLaw, Ellipse2, Intrinsic2, Line2, Polyline2,
};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::{exact_properties, volume_properties};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, DeviationReport, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, OpenProfile, SectionAtStation,
    SolidOperation, Station, StationFrame, StationOffsets, StationOrientation, SurfaceRelation,
};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_reference::station::{station_section2, station_section3, SectionFrame};

/// The chord budget of every run here, in metres.
const BUDGET: Scalar = 0.01;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(BUDGET)
        .unwrap()
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
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

fn plan_circle(radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

fn section(profile: NodeId, distance: Scalar, offsets: StationOffsets) -> SectionAtStation {
    SectionAtStation::new(profile, Station::new(distance, offsets))
}

fn spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<SectionAtStation>,
    frame: StationFrame,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame,
        },
    ))
    .unwrap()
}

/// What the run is, as the oracle reads it: the frame law, the stations,
/// every vertex's profile point per station, the wall edges.
struct Oracle {
    frame_at: Box<dyn Fn(Scalar) -> SectionFrame>,
    distances: Vec<Scalar>,
    locals: Vec<Vec<Point2>>,
    offsets: Vec<StationOffsets>,
    orientations: Vec<StationOrientation>,
    edges: Vec<(usize, usize)>,
}

/// The vertices and edges of a spine's sections, as the compiler matches
/// them by ring and index.
fn spine_oracle(
    graph: &GeometryGraph,
    sections: &[SectionAtStation],
    frame_at: Box<dyn Fn(Scalar) -> SectionFrame>,
) -> Oracle {
    let rings: Vec<_> = sections
        .iter()
        .map(|s| {
            let Some(GeometryNode::Profile(p)) = graph.get(s.profile) else {
                panic!("a profile")
            };
            profile_rings(p, 0.5 * BUDGET, Tolerance::MILLIMETRE).unwrap()
        })
        .collect();
    let mut edges = Vec::new();
    let mut base = 0;
    for len in std::iter::once(rings[0].outer.len()).chain(rings[0].holes.iter().map(Vec::len)) {
        edges.extend((0..len).map(|k| (base + k, base + (k + 1) % len)));
        base += len;
    }
    Oracle {
        frame_at,
        distances: sections.iter().map(|s| s.station.distance).collect(),
        locals: rings
            .iter()
            .map(|r| {
                r.outer
                    .iter()
                    .chain(r.holes.iter().flatten())
                    .copied()
                    .collect()
            })
            .collect(),
        offsets: sections.iter().map(|s| s.station.offsets).collect(),
        orientations: sections.iter().map(|s| s.orientation).collect(),
        edges,
    }
}

/// Points of the exact surface: every wall edge at `per_interval` + 1
/// distances per station interval, `per_edge` + 1 points along it.
fn exact_points(oracle: &Oracle, per_interval: usize, per_edge: usize) -> Vec<Point3> {
    let mut out = Vec::new();
    for i in 0..oracle.distances.len() - 1 {
        let (d0, d1) = (oracle.distances[i], oracle.distances[i + 1]);
        for a in 0..=per_interval {
            let u = a as Scalar / per_interval as Scalar;
            let s = d0 + (d1 - d0) * u;
            let base = (oracle.frame_at)(s);
            let (oa, ob) = (&oracle.orientations[i], &oracle.orientations[i + 1]);
            let turned = if oa.is_base() && ob.is_base() {
                base
            } else {
                let ((xa, ra), (xb, rb)) = (oa.unit_axes().unwrap(), ob.unit_axes().unwrap());
                base.oriented(Some(xa.lerp(xb, u)), Some(ra.lerp(rb, u)))
                    .unwrap()
            };
            let (p, q) = (oracle.offsets[i], oracle.offsets[i + 1]);
            let lerp = |x: Scalar, y: Scalar| x + (y - x) * u;
            let origin = base.place(
                lerp(p.lateral, q.lateral),
                lerp(p.vertical, q.vertical),
                lerp(p.longitudinal, q.longitudinal),
            );
            let vertex = |k: usize| {
                let (a, b) = (oracle.locals[i][k], oracle.locals[i + 1][k]);
                origin + lerp(a.x, b.x) * turned.lateral + lerp(a.y, b.y) * turned.up
            };
            for &(k, l) in &oracle.edges {
                let (vk, vl) = (vertex(k), vertex(l));
                for c in 0..=per_edge {
                    let v = c as Scalar / per_edge as Scalar;
                    out.push(vk + (vl - vk) * v);
                }
            }
        }
    }
    out
}

fn point_triangle(p: Point3, [a, b, c]: [Point3; 3]) -> Scalar {
    axiolid_measure::proximity::closest_point_on_triangle(p, [a, b, c])
        .map(|q| p.distance(q))
        .unwrap_or(Scalar::INFINITY)
}

/// The largest distance from `points` to the mesh's triangles.
fn measured(mesh: &TriMesh, points: &[Point3]) -> Scalar {
    let triangles: Vec<[Point3; 3]> = mesh
        .indices
        .chunks_exact(3)
        .map(|t| [0, 1, 2].map(|k| mesh.positions[t[k] as usize]))
        .collect();
    points
        .iter()
        .map(|&p| {
            triangles
                .iter()
                .map(|t| point_triangle(p, *t))
                .fold(Scalar::INFINITY, Scalar::min)
        })
        .fold(0.0, Scalar::max)
}

fn certified(report: &DeviationReport, path: DeviationPath) -> Scalar {
    let contribution = report
        .contributions
        .iter()
        .find(|c| c.path == path)
        .expect("the run reports its path");
    match contribution.bound {
        DeviationBound::Certified(value) => value,
        other => panic!("not certified: {other:?}"),
    }
}

/// Compile `root`, check the run is certified within the budget, and that
/// a dense sampling of the exact surface stays within the bound.
fn assert_certified(graph: &GeometryGraph, root: NodeId, path: DeviationPath, oracle: &Oracle) {
    let (_, report) = compiler()
        .compile_mesh_with_deviation(graph, root, &options())
        .unwrap();
    assert!(report.meets_requested(), "{report:?}");
    assert_holds(graph, root, path, oracle);
}

/// The run's certified bound, which a dense sampling of the exact surface
/// never exceeds.
fn assert_holds(
    graph: &GeometryGraph,
    root: NodeId,
    path: DeviationPath,
    oracle: &Oracle,
) -> Scalar {
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(graph, root, &options())
        .unwrap();
    let bound = certified(&report, path);
    let points = exact_points(oracle, 160, 8);
    let far = measured(&outcome.mesh, &points);
    // Rounding of the points' own coordinates aside, as every bound is.
    assert!(
        far <= bound + 1e-12 * 100.0,
        "the exact surface strays {far} from the mesh, past the bound {bound}"
    );
    bound
}

fn volume(mesh: &TriMesh) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume
}

fn close(actual: Scalar, expected: Scalar, relative: Scalar, what: &str) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{what}: {actual} != {expected}"
    );
}

// --- certified bounds -------------------------------------------------------

#[test]
fn a_taper_along_a_line_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let line = Curve2::Line(Line2 {
        origin: Point2::new(2.0, 1.0),
        direction: Vec2::new(0.6, 0.8),
    });
    let basis = b.push_value(line.clone()).unwrap();
    let (small, big) = (rectangle(&mut b, 1.0, 1.0), rectangle(&mut b, 1.0, 3.0));
    let sections = vec![
        section(small, 0.0, StationOffsets::ZERO),
        section(big, 10.0, StationOffsets::new(0.5, 0.2, 0.0)),
    ];
    let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
    let graph = b.finish(vec![root]).unwrap();
    let oracle = spine_oracle(
        &graph,
        &sections,
        Box::new(move |s| station_section2(&line, s).unwrap()),
    );
    assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
}

#[test]
fn a_constant_and_a_tapered_section_along_an_arc_are_certified() {
    for (wide, offsets) in [
        (2.0, StationOffsets::ZERO),
        (3.0, StationOffsets::new(0.4, -0.3, 0.1)),
    ] {
        let mut b = GeometryGraphBuilder::new();
        let circle = plan_circle(25.0);
        let basis = b.push_value(circle.clone()).unwrap();
        let (first, last) = (rectangle(&mut b, 2.0, 1.0), rectangle(&mut b, wide, 1.0));
        let sections = vec![
            section(first, 0.0, StationOffsets::ZERO),
            section(last, 20.0, offsets),
            section(last, 30.0, StationOffsets::ZERO),
        ];
        let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
        let graph = b.finish(vec![root]).unwrap();
        let oracle = spine_oracle(
            &graph,
            &sections,
            Box::new(move |s| station_section2(&circle, s).unwrap()),
        );
        assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
    }
}

#[test]
fn a_section_along_a_vertical_curve_is_certified() {
    for frame in [StationFrame::Section, StationFrame::Plan] {
        let mut b = GeometryGraphBuilder::new();
        let curve = Curve3::Elevated(Elevated3::new(
            plan_circle(60.0),
            ElevationLaw::parabolic(1.0, 0.04, -0.03, 40.0),
        ));
        let basis = b.push_value(curve.clone()).unwrap();
        let profile = rectangle(&mut b, 3.0, 1.5);
        let sections = vec![
            section(profile, 2.0, StationOffsets::new(0.0, 0.75, 0.0)),
            section(profile, 38.0, StationOffsets::new(0.0, 0.75, 0.0)),
        ];
        let root = spine(&mut b, basis, sections.clone(), frame);
        let graph = b.finish(vec![root]).unwrap();
        let oracle = spine_oracle(
            &graph,
            &sections,
            Box::new(move |s| {
                let f = station_section3(&curve, s).unwrap();
                match frame {
                    StationFrame::Plan => f.plan().unwrap(),
                    _ => f,
                }
            }),
        );
        assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
    }
}

#[test]
fn a_section_along_a_banked_transition_is_certified() {
    for convention in [
        BankConvention::TangentRotation,
        BankConvention::VerticalRise,
    ] {
        let mut b = GeometryGraphBuilder::new();
        let curve = Curve3::Banked(Banked3::new(
            Elevated3::new(plan_circle(150.0), ElevationLaw::constant_grade(0.0, 0.015)),
            CantLaw::new(vec![CantPiece::cosine(60.0, 0.0, 0.15)]),
            CantLaw::zero(60.0),
            1.5,
            convention,
        ));
        let basis = b.push_value(curve.clone()).unwrap();
        let profile = rectangle(&mut b, 3.0, 0.5);
        let sections = vec![
            section(profile, 0.0, StationOffsets::new(0.2, 0.25, 0.0)),
            section(profile, 60.0, StationOffsets::new(0.2, 0.25, 0.0)),
        ];
        let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
        let graph = b.finish(vec![root]).unwrap();
        let oracle = spine_oracle(
            &graph,
            &sections,
            Box::new(move |s| station_section3(&curve, s).unwrap()),
        );
        assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
    }
}

#[test]
fn a_section_turning_between_stations_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let circle = plan_circle(40.0);
    let basis = b.push_value(circle.clone()).unwrap();
    let profile = rectangle(&mut b, 2.0, 1.0);
    let sections = vec![
        section(profile, 0.0, StationOffsets::ZERO),
        section(profile, 25.0, StationOffsets::new(0.0, 0.5, 0.0)).with_orientation(
            StationOrientation::new(
                Some(Vec3::new(0.0, 0.15, 1.0)),
                Some(Vec3::new(1.0, 0.0, 0.1)),
            ),
        ),
    ];
    let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
    let graph = b.finish(vec![root]).unwrap();
    let oracle = spine_oracle(
        &graph,
        &sections,
        Box::new(move |s| station_section2(&circle, s).unwrap()),
    );
    assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
}

/// Axioval's deck (ifc-geometry, IFC4X3 sectioned solids): a polyline
/// directrix, a rectangle widening from 2 m to 4 m between two stations.
fn deck(b: &mut GeometryGraphBuilder) -> (NodeId, Vec<SectionAtStation>, Curve2) {
    let polyline = Curve2::Polyline(Polyline2 {
        points: vec![Point2::ZERO, Point2::new(4.0, 0.0), Point2::new(10.0, 0.0)],
        closed: false,
    });
    let basis = b.push_value(polyline.clone()).unwrap();
    let (narrow, wide) = (rectangle(b, 2.0, 1.0), rectangle(b, 4.0, 1.0));
    let sections = vec![
        section(narrow, 0.0, StationOffsets::ZERO),
        section(wide, 10.0, StationOffsets::ZERO),
    ];
    (
        spine(b, basis, sections.clone(), StationFrame::Section),
        sections,
        polyline,
    )
}

#[test]
fn axiovals_deck_is_certified_and_holds_thirty_cubic_metres() {
    let mut b = GeometryGraphBuilder::new();
    let (root, sections, polyline) = deck(&mut b);
    let graph = b.finish(vec![root]).unwrap();
    let oracle = spine_oracle(
        &graph,
        &sections,
        Box::new(move |s| station_section2(&polyline, s).unwrap()),
    );
    assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
    let mesh = compiler()
        .compile_mesh_with_deviation(&graph, root, &options())
        .unwrap()
        .0
        .mesh;
    close(volume(&mesh), 30.0, 1e-12, "deck volume");
}

#[test]
fn a_sectioned_surface_along_an_arc_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let circle = plan_circle(30.0);
    let basis = b.push_value(circle.clone()).unwrap();
    let open = |b: &mut GeometryGraphBuilder, points: Vec<Point2>| {
        let path = b
            .push_value(Curve2::Polyline(Polyline2 {
                points,
                closed: false,
            }))
            .unwrap();
        b.push(GeometryNode::OpenProfile(OpenProfile { path }))
            .unwrap()
    };
    let narrow = open(
        &mut b,
        vec![
            Point2::new(-2.0, -0.04),
            Point2::ZERO,
            Point2::new(2.0, -0.04),
        ],
    );
    let wide = open(
        &mut b,
        vec![
            Point2::new(-3.0, -0.06),
            Point2::ZERO,
            Point2::new(3.0, -0.06),
        ],
    );
    let sections = vec![
        section(narrow, 0.0, StationOffsets::ZERO),
        section(wide, 25.0, StationOffsets::new(0.0, 0.3, 0.0)),
    ];
    let root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::OpenSectionsAtStations {
                directrix: basis,
                sections: sections.clone(),
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let local = |points: &[Point2]| points.to_vec();
    let oracle = Oracle {
        frame_at: Box::new(move |s| station_section2(&circle, s).unwrap()),
        distances: vec![0.0, 25.0],
        locals: vec![
            local(&[
                Point2::new(-2.0, -0.04),
                Point2::ZERO,
                Point2::new(2.0, -0.04),
            ]),
            local(&[
                Point2::new(-3.0, -0.06),
                Point2::ZERO,
                Point2::new(3.0, -0.06),
            ]),
        ],
        offsets: sections.iter().map(|s| s.station.offsets).collect(),
        orientations: sections.iter().map(|s| s.orientation).collect(),
        edges: vec![(0, 1), (1, 2)],
    };
    assert_certified(&graph, root, DeviationPath::SectionedSurface, &oracle);
}

#[test]
fn a_twisted_taper_along_a_line_is_certified_by_its_walls() {
    // The second rectangle is the first turned by 0.4 rad: every vertex
    // runs on a straight line, so the walls' twist is the whole bound.
    let mut b = GeometryGraphBuilder::new();
    let line = Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::Y,
    });
    let basis = b.push_value(line.clone()).unwrap();
    let rect = RectangleProfile {
        x: 3.0,
        y: 1.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    };
    let first = b.push_value(Profile::Rectangle(rect)).unwrap();
    let turned = b
        .push_value(Profile::Derived {
            basis: Box::new(Profile::Rectangle(rect)),
            transform: Transform2::from_angle(0.4),
        })
        .unwrap();
    let sections = vec![
        section(first, 0.0, StationOffsets::ZERO),
        section(turned, 6.0, StationOffsets::ZERO),
    ];
    let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
    let graph = b.finish(vec![root]).unwrap();
    let oracle = spine_oracle(
        &graph,
        &sections,
        Box::new(move |s| station_section2(&line, s).unwrap()),
    );
    assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
}

#[test]
fn a_run_across_a_smooth_seam_is_certified() {
    // A constant grade running into a crest: the grade is continuous at
    // 20 m, its rate of change is not, so the seam splits the cells there.
    let mut b = GeometryGraphBuilder::new();
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::Piecewise {
            breaks: vec![20.0],
            laws: vec![
                ElevationLaw::constant_grade(0.0, 0.04),
                ElevationLaw::parabolic(0.8, 0.04, 0.44, 10.0),
            ],
        },
    ));
    let basis = b.push_value(curve.clone()).unwrap();
    let profile = rectangle(&mut b, 2.0, 1.0);
    let sections = vec![
        section(profile, 3.0, StationOffsets::new(0.0, 2.0, 0.0)),
        section(profile, 32.0, StationOffsets::new(0.0, 2.0, 0.0)),
    ];
    let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
    let graph = b.finish(vec![root]).unwrap();
    let oracle = spine_oracle(
        &graph,
        &sections,
        Box::new(move |s| station_section3(&curve, s).unwrap()),
    );
    assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
}

#[test]
fn sheets_turning_between_stations_along_a_line_are_certified() {
    let (sin, cos) = 0.6_f64.sin_cos();
    // A yaw about up turns a level sheet's lateral in its own plane, a
    // pitch towards the tangent turns a standing sheet's up in its own:
    // the walls stay flat, so the turn alone curves the vertex paths.
    let yaw = StationOrientation::new(None, Some(Vec3::new(cos, sin, 0.0)));
    let pitch = StationOrientation::new(Some(Vec3::new(sin, 0.0, cos)), None);
    for (points, roll) in [
        (
            vec![Point2::new(-4.0, 0.0), Point2::ZERO, Point2::new(4.0, 0.0)],
            yaw,
        ),
        (
            vec![Point2::ZERO, Point2::new(0.0, 3.0), Point2::new(0.0, 6.0)],
            pitch,
        ),
    ] {
        let mut b = GeometryGraphBuilder::new();
        let line = Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        });
        let basis = b.push_value(line.clone()).unwrap();
        let path = b
            .push_value(Curve2::Polyline(Polyline2 {
                points: points.clone(),
                closed: false,
            }))
            .unwrap();
        let open = b
            .push(GeometryNode::OpenProfile(OpenProfile { path }))
            .unwrap();
        let sections = vec![
            section(open, 0.0, StationOffsets::ZERO),
            section(open, 12.0, StationOffsets::ZERO).with_orientation(roll),
        ];
        let root = b
            .push(GeometryNode::SurfaceRelation(
                SurfaceRelation::OpenSectionsAtStations {
                    directrix: basis,
                    sections: sections.clone(),
                    frame: StationFrame::Section,
                },
            ))
            .unwrap();
        let graph = b.finish(vec![root]).unwrap();
        let oracle = Oracle {
            frame_at: Box::new(move |s| station_section2(&line, s).unwrap()),
            distances: vec![0.0, 12.0],
            locals: vec![points.clone(), points],
            offsets: sections.iter().map(|s| s.station.offsets).collect(),
            orientations: sections.iter().map(|s| s.orientation).collect(),
            edges: vec![(0, 1), (1, 2)],
        };
        assert_certified(&graph, root, DeviationPath::SectionedSurface, &oracle);
    }
}

#[test]
fn sections_moving_far_from_a_tight_circle_are_certified() {
    // Offsets drifting along the tangent of a 12 m circle, a section held
    // far ahead on it, and a section widening fast: each term of the
    // vertex paths' curvature in turn dominates.
    let cases = [
        (
            1.0,
            1.0,
            StationOffsets::new(0.0, 0.5, 0.0),
            StationOffsets::new(0.0, 0.5, 30.0),
        ),
        (
            1.0,
            1.0,
            StationOffsets::new(0.0, 0.0, 20.0),
            StationOffsets::new(0.0, 0.0, 20.0),
        ),
        (1.0, 10.0, StationOffsets::ZERO, StationOffsets::ZERO),
    ];
    for (first, last, a, z) in cases {
        let mut b = GeometryGraphBuilder::new();
        let circle = plan_circle(12.0);
        let basis = b.push_value(circle.clone()).unwrap();
        let (p, q) = (rectangle(&mut b, first, 1.0), rectangle(&mut b, last, 1.0));
        let sections = vec![section(p, 0.0, a), section(q, 15.0, z)];
        let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
        let graph = b.finish(vec![root]).unwrap();
        let oracle = spine_oracle(
            &graph,
            &sections,
            Box::new(move |s| station_section2(&circle, s).unwrap()),
        );
        assert_certified(&graph, root, DeviationPath::StationedSpine, &oracle);
    }
}

#[test]
fn a_jump_in_the_roll_is_bounded_by_the_tear_it_leaves() {
    // The cant jumps from 0 to 0.15 m at 30 m: the frame rolls by
    // asin(0.1) there, and the section's edges 1.5 m out tear by about
    // 0.15 m, which no refinement closes; the bound says so.
    let mut b = GeometryGraphBuilder::new();
    let curve = Curve3::Banked(Banked3::new(
        Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            ElevationLaw::level(0.0),
        ),
        CantLaw::new(vec![
            CantPiece::constant(30.0, 0.0),
            CantPiece::constant(30.0, 0.15),
        ]),
        CantLaw::zero(60.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    let basis = b.push_value(curve.clone()).unwrap();
    let profile = rectangle(&mut b, 3.0, 0.5);
    let sections = vec![
        section(profile, 0.0, StationOffsets::ZERO),
        section(profile, 50.0, StationOffsets::ZERO),
    ];
    let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
    let graph = b.finish(vec![root]).unwrap();
    let oracle = spine_oracle(
        &graph,
        &sections,
        Box::new(move |s| station_section3(&curve, s).unwrap()),
    );
    let bound = assert_holds(&graph, root, DeviationPath::StationedSpine, &oracle);
    assert!(bound >= 0.1, "the tear is not bounded: {bound}");
}

// --- unbounded by name ------------------------------------------------------

fn unbounded_reason(graph: &GeometryGraph, root: NodeId) -> &'static str {
    let (_, report) = compiler()
        .compile_mesh_with_deviation(graph, root, &options())
        .unwrap();
    assert!(report.bound.is_none());
    let reason = report
        .unbounded()
        .find_map(|c| match c.bound {
            DeviationBound::Unbounded(reason) if c.path == DeviationPath::StationedSpine => {
                Some(reason)
            }
            _ => None,
        })
        .expect("the run is unbounded by name");
    reason
}

#[test]
fn a_numerical_frame_and_a_curved_taper_are_unbounded_by_name() {
    let mut b = GeometryGraphBuilder::new();
    // A clothoid: its frame rests on a quadrature no bound here reaches
    // (an ellipse's is certified since #316).
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
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let profile = rectangle(&mut b, 1.0, 1.0);
    let disk = |b: &mut GeometryGraphBuilder, radius| {
        b.push_value(Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        }))
        .unwrap()
    };
    let (small, large) = (disk(&mut b, 0.5), disk(&mut b, 0.8));
    let along_clothoid = spine(
        &mut b,
        clothoid,
        vec![
            section(profile, 0.0, StationOffsets::ZERO),
            section(profile, 10.0, StationOffsets::ZERO),
        ],
        StationFrame::Section,
    );
    let curved_taper = spine(
        &mut b,
        line,
        vec![
            section(small, 0.0, StationOffsets::ZERO),
            section(large, 10.0, StationOffsets::ZERO),
        ],
        StationFrame::Section,
    );
    let one_disk = spine(
        &mut b,
        line,
        vec![
            section(small, 0.0, StationOffsets::ZERO),
            section(small, 10.0, StationOffsets::ZERO),
        ],
        StationFrame::Section,
    );
    let graph = b
        .finish(vec![along_clothoid, curved_taper, one_disk])
        .unwrap();
    let reason = unbounded_reason(&graph, along_clothoid);
    assert!(reason.contains("certified derivative bounds"), "{reason}");
    let reason = unbounded_reason(&graph, curved_taper);
    assert!(reason.contains("curved profiles"), "{reason}");
    // The exact compiler refuses that taper by name too.
    assert!(refusal(&graph, curved_taper).contains("curved section changes"));
    // One curved profile throughout is its flattening plus the run.
    let (_, report) = compiler()
        .compile_mesh_with_deviation(&graph, one_disk, &options())
        .unwrap();
    assert!(report.meets_requested(), "{report:?}");
}

// --- exact bodies -----------------------------------------------------------

fn exact(graph: &GeometryGraph, root: NodeId) -> (Scalar, Scalar) {
    let body = ReferenceExactCompiler::new()
        .compile_exact(graph, root, &options())
        .unwrap();
    let properties = exact_properties(&body, Tolerance::MILLIMETRE).unwrap();
    (properties.signed_volume, properties.area)
}

#[test]
fn axiovals_deck_compiles_exactly_with_its_closed_form() {
    let mut b = GeometryGraphBuilder::new();
    let (root, _, _) = deck(&mut b);
    let graph = b.finish(vec![root]).unwrap();
    let (volume, area) = exact(&graph, root);
    close(volume, 30.0, 1e-12, "deck volume");
    // Caps 2 and 4, deck and soffit trapezoids of 30 each, and two side
    // walls leaning out by 1 m over 10 m, 1 m high.
    close(
        area,
        2.0 + 4.0 + 60.0 + 2.0 * 101.0_f64.sqrt(),
        1e-12,
        "deck area",
    );
}

#[test]
fn a_constant_section_along_a_line_compiles_to_its_prism() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::new(1.0, 2.0),
            direction: Vec2::new(3.0, 4.0),
        }))
        .unwrap();
    let rect = rectangle(&mut b, 2.0, 0.5);
    let disk = b
        .push_value(Profile::Circle(CircleProfile {
            radius: 0.5,
            thickness: None,
        }))
        .unwrap();
    let offsets = StationOffsets::new(0.3, 1.0, 0.2);
    let rect_run = spine(
        &mut b,
        basis,
        vec![section(rect, 1.0, offsets), section(rect, 7.0, offsets)],
        StationFrame::Section,
    );
    let disk_run = spine(
        &mut b,
        basis,
        vec![section(disk, 0.0, offsets), section(disk, 8.0, offsets)],
        StationFrame::Section,
    );
    let graph = b.finish(vec![rect_run, disk_run]).unwrap();
    let (volume, area) = exact(&graph, rect_run);
    close(volume, 2.0 * 0.5 * 6.0, 1e-12, "prism volume");
    close(area, 2.0 * 1.0 + 5.0 * 6.0, 1e-12, "prism area");
    let (volume, area) = exact(&graph, disk_run);
    let pi = core::f64::consts::PI;
    close(volume, pi * 0.25 * 8.0, 1e-9, "cylinder volume");
    close(area, 2.0 * pi * 0.25 + pi * 8.0, 1e-9, "cylinder area");
    // The mesh agrees with the body.
    let mesh = compiler()
        .compile_mesh_with_deviation(&graph, rect_run, &options())
        .unwrap()
        .0
        .mesh;
    close(volume_of(&mesh), 6.0, 1e-12, "prism mesh");
}

fn volume_of(mesh: &TriMesh) -> Scalar {
    volume(mesh)
}

#[test]
fn a_constant_section_along_an_arc_compiles_to_its_revolution() {
    for (radius, lateral) in [(25.0, 0.0), (25.0, 3.0), (25.0, -2.0)] {
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push_value(plan_circle(radius)).unwrap();
        let rect = rectangle(&mut b, 2.0, 1.0);
        let offsets = StationOffsets::new(lateral, 0.5, 0.0);
        let root = spine(
            &mut b,
            basis,
            vec![section(rect, 5.0, offsets), section(rect, 35.0, offsets)],
            StationFrame::Section,
        );
        let graph = b.finish(vec![root]).unwrap();
        let (volume, area) = exact(&graph, root);
        // Pappus: the centroid runs on the circle of radius r - lateral
        // (lateral is to the left, towards the centre of a turn to the
        // left), through the angle the stations span.
        let theta = 30.0 / radius;
        let centroid = radius - lateral;
        close(volume, 2.0 * 1.0 * centroid * theta, 1e-9, "Pappus volume");
        close(
            area,
            2.0 * 2.0 + 6.0 * centroid * theta,
            1e-9,
            "Pappus area",
        );
        let mesh = compiler()
            .compile_mesh_with_deviation(&graph, root, &options())
            .unwrap()
            .0
            .mesh;
        close(volume_of(&mesh), volume, 1e-3, "mesh against body");
        // And it stands where the mesh does.
        let body = ReferenceExactCompiler::new()
            .compile_exact(&graph, root, &options())
            .unwrap();
        let at = exact_properties(&body, Tolerance::MILLIMETRE)
            .unwrap()
            .centroid;
        let meshed = volume_properties(&mesh, Tolerance::MILLIMETRE)
            .unwrap()
            .centroid;
        assert!((at - meshed).length() <= 0.05, "{at:?} != {meshed:?}");
    }
}

#[test]
fn a_mitre_on_a_curved_piece_or_with_a_changing_section_is_unbounded_by_name() {
    let mut b = GeometryGraphBuilder::new();
    // A grade break on a curved plan: the run is mitred between two
    // curved pieces.
    let crest = b
        .push_value(Curve3::Elevated(Elevated3::new(
            plan_circle(80.0),
            ElevationLaw::Piecewise {
                breaks: vec![20.0],
                laws: vec![
                    ElevationLaw::constant_grade(0.0, 0.0),
                    ElevationLaw::constant_grade(0.0, 0.08),
                ],
            },
        )))
        .unwrap();
    let bend = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![
                Point2::ZERO,
                Point2::new(10.0, 0.0),
                Point2::new(10.0, 10.0),
            ],
            closed: false,
        }))
        .unwrap();
    let (rect, wide) = (rectangle(&mut b, 2.0, 1.0), rectangle(&mut b, 3.0, 1.0));
    let curved = spine(
        &mut b,
        crest,
        vec![
            section(rect, 0.0, StationOffsets::ZERO),
            section(rect, 40.0, StationOffsets::ZERO),
        ],
        StationFrame::Section,
    );
    let tapered = spine(
        &mut b,
        bend,
        vec![
            section(rect, 2.0, StationOffsets::ZERO),
            section(wide, 18.0, StationOffsets::ZERO),
        ],
        StationFrame::Section,
    );
    let graph = b.finish(vec![curved, tapered]).unwrap();
    for root in [curved, tapered] {
        let reason = unbounded_reason(&graph, root);
        assert!(reason.contains("mitred seam"), "{reason}");
    }
}

fn refusal(graph: &GeometryGraph, root: NodeId) -> &'static str {
    match ReferenceExactCompiler::new().compile_exact(graph, root, &options()) {
        Err(GeomError::UnsupportedInput { input, .. }) => input,
        other => panic!("not refused by name: {other:?}"),
    }
}

#[test]
fn every_other_run_is_refused_exactly_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let circle = b.push_value(plan_circle(20.0)).unwrap();
    let ellipse = b
        .push_value(Curve2::Ellipse(Ellipse2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            semi_axis_x: 30.0,
            semi_axis_y: 20.0,
        }))
        .unwrap();
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let bent = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::ZERO, Point2::new(5.0, 0.0), Point2::new(9.0, 3.0)],
            closed: false,
        }))
        .unwrap();
    let rect = rectangle(&mut b, 2.0, 1.0);
    let wide = rectangle(&mut b, 3.0, 1.0);
    let twisted = b
        .push_value(Profile::Derived {
            basis: Box::new(Profile::Rectangle(RectangleProfile {
                x: 2.0,
                y: 1.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            })),
            transform: Transform2::from_angle(0.3),
        })
        .unwrap();
    let run = |b: &mut GeometryGraphBuilder, basis, first, last| {
        spine(
            b,
            basis,
            vec![
                section(first, 0.0, StationOffsets::ZERO),
                section(last, 8.0, StationOffsets::ZERO),
            ],
            StationFrame::Section,
        )
    };
    let along_ellipse = run(&mut b, ellipse, rect, rect);
    let along_bend = run(&mut b, bent, rect, rect);
    let taper_on_arc = run(&mut b, circle, rect, wide);
    let warped = run(&mut b, line, rect, twisted);
    let graph = b
        .finish(vec![along_ellipse, along_bend, taper_on_arc, warped])
        .unwrap();
    assert!(refusal(&graph, along_ellipse).contains("basis other than a line"));
    assert!(refusal(&graph, along_bend).contains("basis other than a line"));
    assert!(refusal(&graph, taper_on_arc).contains("changes along a circle"));
    assert!(refusal(&graph, warped).contains("not coplanar"));
}
