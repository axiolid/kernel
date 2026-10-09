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

#[test]
fn an_oriented_station_and_general_sections_validate_and_reference() {
    use axiolid_core::Vec3;
    use axiolid_model::{OrientedCurveStation, SectionAtStation, StationOrientation};
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let turn = StationOrientation::new(Some(Vec3::new(0.0, 1.0, 1.0)), None);
    let id = b
        .push_value(OrientedCurveStation::new(
            CurveStation::new(basis, Station::at(2.0)),
            turn,
        ))
        .unwrap();
    invalid(
        b.push_value(OrientedCurveStation::new(
            CurveStation::new(basis, Station::at(-1.0)),
            turn,
        )),
        "negative",
    );
    invalid(
        b.push_value(OrientedCurveStation::new(
            CurveStation::new(basis, Station::at(1.0)),
            StationOrientation::new(Some(Vec3::Z), Some(Vec3::new(0.0, 0.0, -2.0))),
        )),
        "parallel",
    );
    assert!(StationOrientation::default().is_base());
    assert_eq!(
        StationOrientation::default().unit_axes(),
        Ok((Vec3::Z, Vec3::X))
    );
    let profile = rectangle(&mut b);
    let spine = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: vec![
                    SectionAtStation::new(profile, Station::at(0.0)).with_orientation(turn),
                    SectionAtStation::from(StationedSection {
                        profile,
                        station: Station::at(4.0),
                    }),
                ],
                frame: StationFrame::Plan,
            },
        ))
        .unwrap();
    // An open section may run its tags backwards, in either form.
    let path = open(&mut b);
    let reversed = |names: [&[&str]; 2]| {
        vec![
            StationedOpenSection {
                profile: path,
                tags: tags(names[0]),
                station: Station::at(0.0),
            },
            StationedOpenSection {
                profile: path,
                tags: tags(names[1]),
                station: Station::at(3.0),
            },
        ]
    };
    let forwards_back = reversed([&["l", "c", "r"], &["r", "c", "l"]]);
    b.push(GeometryNode::SurfaceRelation(
        SurfaceRelation::SectionedSurface {
            directrix: basis,
            sections: forwards_back.clone(),
            frame: StationFrame::Plan,
        },
    ))
    .unwrap();
    let general = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::OpenSectionsAtStations {
                directrix: basis,
                sections: forwards_back.into_iter().map(Into::into).collect(),
                frame: StationFrame::Plan,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![id, spine, general]).unwrap();
    assert_eq!(graph.get(id).unwrap().references(), vec![basis]);
    assert_eq!(
        graph.get(spine).unwrap().references(),
        vec![basis, profile, profile]
    );
    assert_eq!(
        graph.get(general).unwrap().references(),
        vec![basis, path, path]
    );
}

#[test]
fn a_station_names_the_side_of_a_seam_it_reads() {
    use axiolid_model::{OrientedCurveStation, SeamSide, StationOrientation};
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let plain = CurveStation::new(basis, Station::at(2.0));
    // Outgoing unless told otherwise, as every curve evaluator reads.
    assert_eq!(SeamSide::default(), SeamSide::Outgoing);
    assert_eq!(OrientedCurveStation::from(plain).seam, SeamSide::Outgoing);
    assert_eq!(
        OrientedCurveStation::new(plain, StationOrientation::default()).seam,
        SeamSide::Outgoing
    );
    let incoming = plain.with_seam_side(SeamSide::Incoming);
    assert_eq!(incoming.seam, SeamSide::Incoming);
    assert_eq!(incoming.station, plain);
    assert!(incoming.orientation.is_base());
    let id = b.push_value(incoming).unwrap();
    let graph = b.finish(vec![id]).unwrap();
    assert_eq!(
        graph.get(id),
        Some(&GeometryNode::OrientedCurveStation(incoming))
    );
    assert_eq!(graph.get(id).unwrap().references(), vec![basis]);
}

#[test]
fn an_instance_at_a_station_round_trips_and_references_source_and_basis() {
    use axiolid_core::Vec3;
    use axiolid_model::{InstanceAtStation, OrientedCurveStation, SeamSide, StationOrientation};
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let source = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::Y,
        }))
        .unwrap();
    let station = OrientedCurveStation::new(
        CurveStation {
            basis,
            station: Station::new(3.5, StationOffsets::new(0.25, -1.0, 0.5)),
            frame: StationFrame::Plan,
        },
        StationOrientation::new(Some(Vec3::new(0.0, 0.6, 0.8)), None),
    )
    .with_seam_side(SeamSide::Incoming);
    let placed = InstanceAtStation::new(source, station);
    let id = b.push_value(placed).unwrap();
    let graph = b.finish(vec![id]).unwrap();
    let Some(GeometryNode::InstanceAtStation(stored)) = graph.get(id) else {
        panic!("not an instance at a station");
    };
    // Exact stored data: the station stays symbolic, its side included.
    assert_eq!(*stored, placed);
    assert_eq!(stored.station.seam, SeamSide::Incoming);
    assert_eq!(stored.station.station.station.distance, 3.5);
    assert_eq!(graph.get(id).unwrap().references(), vec![source, basis]);
}

#[test]
fn a_curve_placed_at_a_station_is_a_3d_curve_and_a_placed_solid_a_solid() {
    use axiolid_core::Vec3;
    use axiolid_model::{CurveSegment, InstanceAtStation, Transition};
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let source = line(&mut b);
    let at = CurveStation::new(basis, Station::at(1.0)).into();
    let placed = b.push_value(InstanceAtStation::new(source, at)).unwrap();
    // A 3D slot takes it: a sweep's directrix, a composite's segment.
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix: placed,
        radius: 0.1,
        inner_radius: None,
        parameter_range: Some((0.0, 1.0)),
        fillet_radius: None,
    }))
    .unwrap();
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments: vec![CurveSegment {
            curve: placed,
            same_sense: true,
            transition: Transition::Continuous,
        }],
    }))
    .unwrap();
    // A 2D slot does not: it left its plane.
    assert!(matches!(
        b.push_value(OpenProfile::new(placed)),
        Err(GraphError::InvalidReferenceType { .. })
    ));
    let plane = b
        .push_value(axiolid_surface::Surface::Plane(axiolid_surface::Plane {
            frame: axiolid_core::Frame3 {
                origin: Vec3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
        }))
        .unwrap();
    assert!(matches!(
        b.push(GeometryNode::CurveRelation(CurveRelation::ParameterCurve {
            basis_surface: plane,
            reference_curve: placed,
        })),
        Err(GraphError::InvalidReferenceType {
            expected: "curve2",
            ..
        })
    ));
    // ... and a 3D-only slot takes it.
    b.push(GeometryNode::CurveRelation(CurveRelation::SurfaceCurve {
        curve_3d: placed,
        sides: axiolid_model::SurfaceSides::one(plane, source),
        master: axiolid_model::MasterRepresentation::Curve3d,
    }))
    .unwrap();
    // A placed profile is no longer a profile, and not a curve.
    let profile = rectangle(&mut b);
    let placed_profile = b.push_value(InstanceAtStation::new(profile, at)).unwrap();
    assert!(matches!(
        b.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: placed_profile,
            direction: Vec3::Z,
            depth: 1.0,
        })),
        Err(GraphError::InvalidReferenceType {
            expected: "profile",
            ..
        })
    ));
    assert!(matches!(
        b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix: placed_profile,
            radius: 0.1,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        })),
        Err(GraphError::InvalidReferenceType {
            expected: "curve",
            ..
        })
    ));
    // A placed solid stays a solid.
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let placed_solid = b.push_value(InstanceAtStation::new(solid, at)).unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: placed_solid,
        right: solid,
        operator: axiolid_model::BooleanOperator::Difference,
    }))
    .unwrap();
}

#[test]
fn an_instance_at_a_malformed_station_is_refused_by_name() {
    use axiolid_core::Vec3;
    use axiolid_model::{InstanceAtStation, OrientedCurveStation, StationOrientation};
    let mut b = GeometryGraphBuilder::new();
    let basis = line(&mut b);
    let source = line(&mut b);
    let push = |b: &mut GeometryGraphBuilder, station: Station, turn| {
        b.push_value(InstanceAtStation::new(
            source,
            OrientedCurveStation::new(CurveStation::new(basis, station), turn),
        ))
    };
    invalid(
        push(&mut b, Station::at(-2.0), StationOrientation::default()),
        "negative",
    );
    invalid(
        push(
            &mut b,
            Station::new(1.0, StationOffsets::new(0.0, Scalar::INFINITY, 0.0)),
            StationOrientation::default(),
        ),
        "offset is not finite",
    );
    invalid(
        push(
            &mut b,
            Station::at(1.0),
            StationOrientation::new(None, Some(Vec3::new(0.0, 0.0, 5.0))),
        ),
        "parallel",
    );
    // The basis must be a curve.
    let profile = rectangle(&mut b);
    assert!(matches!(
        b.push_value(InstanceAtStation::new(
            source,
            CurveStation::new(profile, Station::at(1.0)).into(),
        )),
        Err(GraphError::InvalidReferenceType {
            expected: "curve",
            ..
        })
    ));
}

#[test]
fn a_station_along_a_composite_of_placed_segments_round_trips() {
    // #285: a composite of a segment placed at a station of another
    // composite is a station basis: pushed, stored exactly, references
    // followed through every level.
    use axiolid_model::{
        CurveSegment, InstanceAtStation, OrientedCurveStation, SeamSide, Transition,
    };
    let mut b = GeometryGraphBuilder::new();
    let first = line(&mut b);
    let second = line(&mut b);
    let segment = |curve| CurveSegment {
        curve,
        same_sense: true,
        transition: Transition::ContinuousSameGradient,
    };
    let base = b
        .push(GeometryNode::CurveRelation(CurveRelation::Composite {
            segments: vec![segment(first), segment(second)],
        }))
        .unwrap();
    let source = line(&mut b);
    let at = OrientedCurveStation::from(CurveStation::new(base, Station::at(1.0)))
        .with_seam_side(SeamSide::Incoming);
    let placed = b.push_value(InstanceAtStation::new(source, at)).unwrap();
    let segmented = b
        .push(GeometryNode::CurveRelation(CurveRelation::Composite {
            segments: vec![segment(placed)],
        }))
        .unwrap();
    let station = CurveStation {
        basis: segmented,
        station: Station::new(0.5, StationOffsets::new(0.25, 0.0, 0.0)),
        frame: StationFrame::Plan,
    };
    let id = b.push_value(station).unwrap();
    let graph = b.finish(vec![id]).unwrap();
    assert_eq!(graph.get(id), Some(&GeometryNode::CurveStation(station)));
    assert_eq!(graph.get(id).unwrap().references(), vec![segmented]);
    assert_eq!(graph.get(segmented).unwrap().references(), vec![placed]);
    let Some(GeometryNode::InstanceAtStation(stored)) = graph.get(placed) else {
        panic!("not a placed segment");
    };
    assert_eq!(stored.station, at);
    assert_eq!(graph.get(placed).unwrap().references(), vec![source, base]);
    assert_eq!(graph.get(base).unwrap().references(), vec![first, second]);
}
