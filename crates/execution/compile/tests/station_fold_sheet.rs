//! Sheets of open sections whose section planes fold are refused (#407, ADR
//! 0082 amendment). `SurfaceRelation::SectionedSurface` and
//! `OpenSectionsAtStations` are judged as a stationed spine is (#361) at the
//! vertices of their polylines, with no flattening margin: `F = 1 - a / R`
//! on a plan arc of radius `R`, `F = 1 + c'` on a line with a longitudinal
//! offset of rate `c'`. A turned section is checked conservatively and never
//! refused: reported unbounded by name where the bound cannot tell.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Elevated3, ElevationLaw, Line2, Polyline2};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, OpenProfile, SectionAtStation,
    Station, StationFrame, StationOffsets, StationOrientation, StationedOpenSection,
    SurfaceRelation,
};

const BUDGET: Scalar = 0.01;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(BUDGET)
        .unwrap()
}

/// What the mesh compiler makes of a sheet.
#[derive(Debug)]
enum Outcome {
    Certified,
    Unbounded(&'static str),
    Refused(String),
}

fn mesh(graph: &GeometryGraph, root: NodeId) -> Outcome {
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    match compiler.compile_mesh_with_deviation(graph, root, &options()) {
        Ok((_, report)) => {
            let run = report
                .contributions
                .iter()
                .find(|c| c.path == DeviationPath::SectionedSurface)
                .expect("the sheet reports its path");
            match run.bound {
                DeviationBound::Certified(_) => Outcome::Certified,
                DeviationBound::Unbounded(reason) => Outcome::Unbounded(reason),
                other => panic!("neither certified nor unbounded: {other:?}"),
            }
        }
        Err(GeomError::InvalidInput(detail)) => Outcome::Refused(detail),
        Err(other) => panic!("not a refusal as invalid input: {other:?}"),
    }
}

/// The refusal's detail, which must name the fold.
fn refused(outcome: Outcome) -> String {
    match outcome {
        Outcome::Refused(detail) => {
            assert!(detail.contains("sections fold"), "{detail}");
            detail
        }
        other => panic!("expected a fold refusal, got {other:?}"),
    }
}

fn circle_curve(radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

fn circle(radius: Scalar) -> GeometryNode {
    GeometryNode::Curve2(circle_curve(radius))
}

fn line() -> GeometryNode {
    GeometryNode::Curve2(Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    }))
}

fn open(b: &mut GeometryGraphBuilder, points: Vec<Point2>) -> NodeId {
    let path = b
        .push_value(Curve2::Polyline(Polyline2 {
            points,
            closed: false,
        }))
        .unwrap();
    b.push(GeometryNode::OpenProfile(OpenProfile { path }))
        .unwrap()
}

fn segment(from: Scalar, to: Scalar) -> Vec<Point2> {
    vec![Point2::new(from, 0.0), Point2::new(to, 0.0)]
}

/// Which relation a sheet is built as.
#[derive(Clone, Copy, Debug)]
enum Kind {
    Sectioned,
    OpenSections,
}

const KINDS: [Kind; 2] = [Kind::Sectioned, Kind::OpenSections];
const FRAMES: [StationFrame; 2] = [StationFrame::Section, StationFrame::Plan];

/// A sheet of one polyline at stations 0 and 10, the second offset by `to`
/// and turned by `turn`.
fn sheet(
    kind: Kind,
    basis: GeometryNode,
    points: Vec<Point2>,
    frame: StationFrame,
    to: StationOffsets,
    turn: Option<StationOrientation>,
) -> Outcome {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(basis).unwrap();
    let p = open(&mut b, points);
    let relation = match kind {
        Kind::OpenSections => {
            let mut last = SectionAtStation::new(p, Station::new(10.0, to));
            if let Some(turn) = turn {
                last = last.with_orientation(turn);
            }
            SurfaceRelation::OpenSectionsAtStations {
                directrix: basis,
                sections: vec![SectionAtStation::new(p, Station::at(0.0)), last],
                frame,
            }
        }
        Kind::Sectioned => {
            assert!(turn.is_none(), "a sectioned surface has no turned section");
            SurfaceRelation::SectionedSurface {
                directrix: basis,
                sections: vec![
                    StationedOpenSection {
                        profile: p,
                        tags: vec![],
                        station: Station::at(0.0),
                    },
                    StationedOpenSection {
                        profile: p,
                        tags: vec![],
                        station: Station::new(10.0, to),
                    },
                ],
                frame,
            }
        }
    };
    let root = b.push(GeometryNode::SurfaceRelation(relation)).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    mesh(&graph, root)
}

/// A polyline from `-3` to `reach` along a plan arc of radius 10.
fn plan_arc(kind: Kind, frame: StationFrame, reach: Scalar) -> Outcome {
    sheet(
        kind,
        circle(10.0),
        segment(-3.0, reach),
        frame,
        StationOffsets::ZERO,
        None,
    )
}

fn arc_centre_cases(kind: Kind) {
    for frame in FRAMES {
        let case = format!("{kind:?} {frame:?}");
        let near = plan_arc(kind, frame, 9.9);
        assert!(matches!(near, Outcome::Certified), "{case}: {near:?}");
        for reach in [10.0, 10.5] {
            let detail = refused(plan_arc(kind, frame, reach));
            assert!(detail.contains("at distance"), "{case} {reach}: {detail}");
        }
        let clear = sheet(
            kind,
            circle(10.0),
            segment(-50.0, -40.0),
            frame,
            StationOffsets::ZERO,
            None,
        );
        assert!(matches!(clear, Outcome::Certified), "{case}: {clear:?}");
    }
}

#[test]
fn a_sectioned_surface_reaching_a_plan_arcs_centre_is_refused_and_one_clear_of_it_is_certified() {
    arc_centre_cases(Kind::Sectioned);
}

#[test]
fn an_open_section_run_reaching_a_plan_arcs_centre_is_refused_and_one_clear_of_it_is_certified() {
    arc_centre_cases(Kind::OpenSections);
}

#[test]
fn a_sheet_receding_faster_than_the_basis_advances_folds() {
    for kind in KINDS {
        for frame in FRAMES {
            let case = format!("{kind:?} {frame:?}");
            refused(sheet(
                kind,
                line(),
                segment(-3.0, 3.0),
                frame,
                StationOffsets::new(0.0, 0.0, -20.0),
                None,
            ));
            let keep = sheet(
                kind,
                line(),
                segment(-3.0, 3.0),
                frame,
                StationOffsets::new(0.0, 0.0, -9.0),
                None,
            );
            assert!(matches!(keep, Outcome::Certified), "{case}: {keep:?}");
        }
    }
}

#[test]
fn a_sheet_wholly_past_the_centre_or_crossing_a_ramps_axis_is_refused() {
    for kind in KINDS {
        for frame in FRAMES {
            refused(sheet(
                kind,
                circle(10.0),
                segment(10.5, 12.0),
                frame,
                StationOffsets::ZERO,
                None,
            ));
            let ramp = GeometryNode::Curve3(Curve3::Elevated(Elevated3::new(
                circle_curve(10.0),
                ElevationLaw::constant_grade(0.0, 0.2),
            )));
            refused(sheet(
                kind,
                ramp,
                segment(-3.0, 12.0),
                frame,
                StationOffsets::ZERO,
                None,
            ));
        }
    }
}

#[test]
fn a_reversed_tagged_sheet_is_judged_on_its_joined_vertices() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(circle(10.0)).unwrap();
    let forward = open(&mut b, segment(-3.0, 12.0));
    let backward = open(&mut b, segment(12.0, -3.0));
    let tags = |names: [&str; 2]| names.iter().map(|n| n.to_string()).collect::<Vec<_>>();
    let root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::SectionedSurface {
                directrix: basis,
                sections: vec![
                    StationedOpenSection {
                        profile: forward,
                        tags: tags(["in", "out"]),
                        station: Station::at(0.0),
                    },
                    StationedOpenSection {
                        profile: backward,
                        tags: tags(["out", "in"]),
                        station: Station::at(10.0),
                    },
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    refused(mesh(&graph, root));
}

#[test]
fn a_turned_sheet_is_checked_conservatively_and_never_refused() {
    let turn = StationOrientation::new(None, Some(Vec3::new(0.98, 0.2, 0.0)));
    for frame in FRAMES {
        let case = format!("{frame:?}");
        let clear = sheet(
            Kind::OpenSections,
            circle(10.0),
            segment(-3.0, 3.0),
            frame,
            StationOffsets::ZERO,
            Some(turn),
        );
        assert!(matches!(clear, Outcome::Certified), "{case}: {clear:?}");
        let past = sheet(
            Kind::OpenSections,
            circle(10.0),
            segment(-3.0, 11.0),
            frame,
            StationOffsets::ZERO,
            Some(turn),
        );
        match past {
            Outcome::Unbounded(reason) => {
                assert!(reason.contains("cannot be excluded"), "{case}: {reason}")
            }
            other => panic!("{case}: expected unbounded, got {other:?}"),
        }
    }
}
