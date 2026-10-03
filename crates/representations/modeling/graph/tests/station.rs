//! Station relations (#241): exact stored data, malformed stations refused
//! by name when the node is pushed, references followed.

use axiolid_core::{Point2, Scalar, Vec2};
use axiolid_curve::{Curve2, Line2, Polyline2};
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraphBuilder, GeometryNode, GraphError, NodeId,
    OpenProfile, SolidOperation, Station, StationFrame, StationOffsets, StationedOpenSection,
    StationedSection, SurfaceRelation,
};
use axiolid_profile::{Profile, RectangleProfile};

fn line(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    }))
    .unwrap()
}

fn rectangle(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x: 1.0,
        y: 2.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
    .unwrap()
}

fn open(b: &mut GeometryGraphBuilder) -> NodeId {
    let path = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(-1.0, 0.0), Point2::ZERO, Point2::new(1.0, 0.0)],
            closed: false,
        }))
        .unwrap();
    b.push_value(OpenProfile::new(path)).unwrap()
}

fn invalid(result: Result<NodeId, GraphError>, needle: &str) {
    match result {
        Err(GraphError::InvalidStation { detail }) => {
            assert!(detail.contains(needle), "expected {needle:?} in {detail:?}")
        }
        other => panic!("expected an invalid station naming {needle:?}, got {other:?}"),
    }
}

fn tags(names: &[&str]) -> Vec<String> {
    names.iter().map(|name| (*name).to_owned()).collect()
}

#[test]
fn a_curve_station_is_stored_exactly_and_references_its_basis() {
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let station = CurveStation {
        basis,
        station: Station::new(12.5, StationOffsets::new(-1.75, 0.3, 0.1)),
        frame: StationFrame::Plan,
    };
    let id = b.push_value(station).unwrap();
    let graph = b.finish(vec![id]).unwrap();
    let Some(GeometryNode::CurveStation(stored)) = graph.get(id) else {
        panic!("not a station");
    };
    assert_eq!(*stored, station);
    assert_eq!(graph.get(id).unwrap().references(), vec![basis]);
    assert_eq!(
        CurveStation::new(basis, Station::at(1.0)).frame,
        StationFrame::Section
    );
}

#[test]
fn malformed_single_stations_are_refused_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let push =
        |b: &mut GeometryGraphBuilder, station| b.push_value(CurveStation::new(basis, station));
    invalid(push(&mut b, Station::at(Scalar::NAN)), "not finite");
    invalid(push(&mut b, Station::at(Scalar::INFINITY)), "not finite");
    invalid(push(&mut b, Station::at(-0.5)), "negative");
    invalid(
        push(
            &mut b,
            Station::new(1.0, StationOffsets::new(Scalar::NAN, 0.0, 0.0)),
        ),
        "offset is not finite",
    );
    // The basis must be a curve.
    let profile = rectangle(&mut b);
    assert!(matches!(
        b.push_value(CurveStation::new(profile, Station::at(1.0))),
        Err(GraphError::InvalidReferenceType {
            expected: "curve",
            ..
        })
    ));
}

#[test]
fn an_offset_curve_needs_two_increasing_stations_and_is_three_dimensional() {
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let offset = |stations: Vec<Station>| {
        GeometryNode::CurveRelation(CurveRelation::OffsetByStations {
            basis,
            stations,
            frame: StationFrame::Section,
        })
    };
    invalid(b.push(offset(vec![Station::at(1.0)])), "at least two");
    invalid(
        b.push(offset(vec![Station::at(2.0), Station::at(2.0)])),
        "increase strictly",
    );
    invalid(
        b.push(offset(vec![Station::at(3.0), Station::at(2.0)])),
        "increase strictly",
    );
    invalid(
        b.push(offset(vec![Station::at(0.0), Station::at(Scalar::NAN)])),
        "not finite",
    );
    let curve = b
        .push(offset(vec![
            Station::new(0.0, StationOffsets::new(1.0, 0.0, 0.0)),
            Station::new(10.0, StationOffsets::new(3.0, 0.5, 0.0)),
        ]))
        .unwrap();
    // A 3D curve: a 2D slot refuses it, a 3D one takes it.
    let plane = b
        .push_value(axiolid_surface::Surface::Plane(axiolid_surface::Plane {
            frame: axiolid_core::Frame3 {
                origin: axiolid_core::Vec3::ZERO,
                x: axiolid_core::Vec3::X,
                y: axiolid_core::Vec3::Y,
                z: axiolid_core::Vec3::Z,
            },
        }))
        .unwrap();
    assert!(matches!(
        b.push(GeometryNode::CurveRelation(CurveRelation::ParameterCurve {
            basis_surface: plane,
            reference_curve: curve,
        })),
        Err(GraphError::InvalidReferenceType {
            expected: "curve2",
            ..
        })
    ));
    assert!(b
        .push(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve {
            curve_3d: curve,
            sides: axiolid_model::SurfaceSides::one(plane, basis),
            master: axiolid_model::MasterRepresentation::Curve3d,
        }))
        .is_ok());
    let disk = b
        .push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix: curve,
            radius: 0.1,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        }))
        .unwrap();
    let graph = b.finish(vec![disk]).unwrap();
    assert_eq!(graph.get(curve).unwrap().references(), vec![basis]);
}

#[test]
fn a_stationed_spine_validates_its_sections() {
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let profile = rectangle(&mut b);
    let spine = |sections: Vec<StationedSection>| {
        GeometryNode::SolidOperation(SolidOperation::StationedSpine {
            directrix: basis,
            sections,
            frame: StationFrame::Section,
        })
    };
    let at = |distance| StationedSection {
        profile,
        station: Station::at(distance),
    };
    invalid(b.push(spine(vec![at(0.0)])), "at least two");
    invalid(b.push(spine(vec![at(5.0), at(1.0)])), "increase strictly");
    invalid(b.push(spine(vec![at(-1.0), at(1.0)])), "negative");
    // A section must be an area profile.
    let open_profile = open(&mut b);
    assert!(matches!(
        b.push(spine(vec![
            at(0.0),
            StationedSection {
                profile: open_profile,
                station: Station::at(1.0),
            },
        ])),
        Err(GraphError::InvalidReferenceType {
            expected: "profile",
            ..
        })
    ));
    let id = b.push(spine(vec![at(0.0), at(4.0), at(10.0)])).unwrap();
    let graph = b.finish(vec![id]).unwrap();
    assert_eq!(
        graph.get(id).unwrap().references(),
        vec![basis, profile, profile, profile]
    );
}

#[test]
fn a_sectioned_surface_validates_its_tags() {
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let profile = open(&mut b);
    let surface = |sections: Vec<StationedOpenSection>| {
        GeometryNode::SurfaceRelation(SurfaceRelation::SectionedSurface {
            directrix: basis,
            sections,
            frame: StationFrame::Plan,
        })
    };
    let at = |distance, names: &[&str]| StationedOpenSection {
        profile,
        tags: tags(names),
        station: Station::at(distance),
    };
    invalid(
        b.push(surface(vec![
            at(0.0, &["left", "axis", "right"]),
            at(10.0, &["left", "right", "axis"]),
        ])),
        "tags are inconsistent",
    );
    invalid(
        b.push(surface(vec![at(0.0, &["a", "b"]), at(10.0, &[])])),
        "tags are inconsistent",
    );
    invalid(
        b.push(surface(vec![
            at(0.0, &["a", "a", "b"]),
            at(10.0, &["a", "a", "b"]),
        ])),
        "repeats a tag",
    );
    invalid(
        b.push(surface(vec![at(0.0, &["a"]), at(0.0, &["a"])])),
        "increase strictly",
    );
    // A section must be an open profile.
    let area = rectangle(&mut b);
    assert!(matches!(
        b.push(surface(vec![
            at(0.0, &[]),
            StationedOpenSection {
                profile: area,
                tags: Vec::new(),
                station: Station::at(5.0),
            },
        ])),
        Err(GraphError::InvalidReferenceType {
            expected: "open profile",
            ..
        })
    ));
    let id = b
        .push(surface(vec![
            at(0.0, &["left", "axis", "right"]),
            at(10.0, &["left", "axis", "right"]),
        ]))
        .unwrap();
    let graph = b.finish(vec![id]).unwrap();
    assert_eq!(
        graph.get(id).unwrap().references(),
        vec![basis, profile, profile]
    );
}
