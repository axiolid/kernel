//! Seam-snapping windows and seam sides on run stations (#301): a section
//! or offset station within its window of one seam of the directrix
//! stands on that seam -- an interior one in the seam's mitre plane, the
//! run's first and last on the pieces the run lies on -- and an interior
//! station on a seam that is not mitred reads the side it names. A window
//! holding two seams, a seam on a window's edge and stations snapped out
//! of order are refused by name; a window of `0` is today's reading.

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec2};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Curve2, Curve3, Elevated3, ElevationLaw, Line2,
    Polyline2,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::{curve_path, resolve};
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveSegment, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode,
    NodeId, OpenProfile, RunStation, SeamSide, SectionAtStation, SolidOperation, Station,
    StationFrame, StationOffsets, StationedSection, SurfaceRelation, Transition, TrimSelector,
    TrimmingPreference,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::station::{
    CompositeBasis, DistanceInterval, Mitre, SEAM_WINDOW_AMBIGUOUS, SEAM_WINDOW_STRADDLED,
    WINDOW_CERTIFICATION,
};

const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];
/// The model's precision: a millimetre.
const WINDOW: Scalar = 1e-3;
/// How far from the joint a station is authored: inside the window.
const NEAR: Scalar = 6e-4;

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

fn has_vertex(mesh: &TriMesh, point: Point3) -> bool {
    mesh.positions
        .iter()
        .any(|p| (*p - point).abs().max_element() <= 1e-9)
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, input: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input: found, .. }) => assert_eq!(found, input),
        other => panic!("expected a refusal naming {input:?}, got {other:?}"),
    }
}

/// An offset by distances at stations along a road on a crest (a gradient
/// curve): its joint at the middle station is at the offset's own plan
/// length, a quadrature, and its plan and grade turn there.
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
fn kerb_joint() -> (Scalar, DistanceInterval) {
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

fn rectangle(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 1.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
    .unwrap()
}

fn deck(b: &mut GeometryGraphBuilder) -> NodeId {
    let path = b
        .push_value(Curve2::Polyline(Polyline2 {
            points: vec![Point2::new(-1.0, 0.0), Point2::new(1.0, 0.0)],
            closed: false,
        }))
        .unwrap();
    b.push_value(OpenProfile::new(path)).unwrap()
}

/// A station of a run: distance, side and window.
type At = (Scalar, SeamSide, Scalar);

/// What a run is built as.
#[derive(Clone, Copy, Debug)]
enum Kind {
    Solid,
    Sheet,
    Offset,
}

const KINDS: [Kind; 3] = [Kind::Solid, Kind::Sheet, Kind::Offset];

/// A run of `kind` along the directrix `build` makes, with `stations`,
/// compiled: a solid of rectangles, a sheet of decks, or a small solid
/// swept along an offset curve by stations (all offsets `(0.5, 0.25,
/// 0)`, so every interpolation is the same whatever the distances).
fn run(
    build: impl Fn(&mut GeometryGraphBuilder) -> NodeId,
    kind: Kind,
    stations: &[At],
) -> GeomResult<TriMesh> {
    let mut b = GeometryGraphBuilder::new();
    let directrix = build(&mut b);
    let offsets = StationOffsets::new(0.5, 0.25, 0.0);
    let root = match kind {
        Kind::Solid | Kind::Sheet => {
            let profile = match kind {
                Kind::Solid => rectangle(&mut b),
                _ => deck(&mut b),
            };
            let sections: Vec<SectionAtStation> = stations
                .iter()
                .map(|&(d, side, window)| {
                    SectionAtStation::new(profile, Station::new(d, offsets))
                        .with_seam_side(side)
                        .with_seam_window(window)
                })
                .collect();
            let node = match kind {
                Kind::Solid => GeometryNode::SolidOperation(SolidOperation::SectionsAtStations {
                    directrix,
                    sections,
                    frame: StationFrame::Section,
                }),
                _ => GeometryNode::SurfaceRelation(SurfaceRelation::OpenSectionsAtStations {
                    directrix,
                    sections,
                    frame: StationFrame::Section,
                }),
            };
            b.push(node).unwrap()
        }
        Kind::Offset => {
            let offset = b
                .push(GeometryNode::CurveRelation(
                    CurveRelation::OffsetsAtStations {
                        basis: directrix,
                        stations: stations
                            .iter()
                            .map(|&(d, side, window)| {
                                RunStation::new(Station::new(d, offsets))
                                    .with_seam_side(side)
                                    .with_seam_window(window)
                            })
                            .collect(),
                        frame: StationFrame::Section,
                    },
                ))
                .unwrap();
            b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
                directrix: offset,
                radius: 0.05,
                inner_radius: None,
                parameter_range: None,
                fillet_radius: None,
            }))
            .unwrap()
        }
    };
    let graph = b.finish(vec![root]).unwrap();
    compile(&graph, root)
}

#[test]
fn an_interior_station_within_its_window_stands_in_the_joints_mitre() {
    let (joint, _) = kerb_joint();
    for kind in KINDS {
        let on = run(
            kerb,
            kind,
            &[
                (5.0, SeamSide::Outgoing, 0.0),
                (joint, SeamSide::Outgoing, 0.0),
                (90.0, SeamSide::Outgoing, 0.0),
            ],
        )
        .unwrap();
        for side in SIDES {
            for near in [joint - NEAR, joint + NEAR] {
                let stations = |window| {
                    [
                        (5.0, SeamSide::Outgoing, 0.0),
                        (near, side, window),
                        (90.0, SeamSide::Outgoing, 0.0),
                    ]
                };
                let snapped = run(kerb, kind, &stations(WINDOW)).unwrap();
                // Exactly the run with the station on the joint, whatever
                // side it names: the mitre is both sides'.
                same(&snapped, &on, &format!("{kind:?} {side:?} at {near}"));
                // Without a window it stands where it is authored, today's
                // reading: a section that wide, 0.6 mm from the mitre, is
                // cut by it and refused; an offset point is elsewhere.
                match run(kerb, kind, &stations(0.0)) {
                    Ok(unsnapped) => {
                        assert!(matches!(kind, Kind::Offset), "{kind:?} not cut");
                        assert_ne!(unsnapped.positions, on.positions, "{kind:?} unsnapped");
                    }
                    Err(GeomError::InvalidInput(detail)) => {
                        assert!(detail.contains("the mitre would cut it"), "{detail}");
                    }
                    Err(other) => panic!("{kind:?}: {other:?}"),
                }
            }
        }
    }
}

/// The section at the joint is the mitre section the two pieces' frames
/// give there: each point placed in both and projected along its own
/// tangent onto the bisector plane.
#[test]
fn the_snapped_section_is_the_mitre_section_of_the_joint() {
    let (joint, _) = kerb_joint();
    let mut b = GeometryGraphBuilder::new();
    let basis = kerb(&mut b);
    let ids: Vec<NodeId> = SIDES
        .iter()
        .map(|&side| {
            b.push(GeometryNode::OrientedCurveStation(
                CurveStation::new(basis, Station::at(joint)).with_seam_side(side),
            ))
            .unwrap()
        })
        .collect();
    let graph = b.finish(ids.clone()).unwrap();
    let incoming = resolve(&graph, ids[0]).unwrap().section;
    let outgoing = resolve(&graph, ids[1]).unwrap().section;
    let mitre = Mitre::between(&incoming, &outgoing)
        .unwrap()
        .expect("the kerb turns at its joint");
    for side in SIDES {
        let solid = run(
            kerb,
            Kind::Solid,
            &[
                (5.0, SeamSide::Outgoing, 0.0),
                (joint + NEAR, side, WINDOW),
                (90.0, SeamSide::Outgoing, 0.0),
            ],
        )
        .unwrap();
        for (x, y) in [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)] {
            let (l, v) = (0.5 + x, 0.25 + y);
            let want = mitre.place(incoming.place(l, v, 0.0), outgoing.place(l, v, 0.0), 0.0);
            assert!(has_vertex(&solid, want), "{side:?}: mitre corner {want:?}");
            // Neither piece's plain section stands there.
            assert!(!has_vertex(&solid, incoming.place(l, v, 0.0)));
            assert!(!has_vertex(&solid, outgoing.place(l, v, 0.0)));
        }
    }
}

/// The run's first station reads the outgoing piece and its last the
/// incoming one, whatever side they name: the pieces the run lies on.
#[test]
fn the_runs_ends_snap_onto_the_joint_and_keep_the_runs_reading() {
    let (joint, _) = kerb_joint();
    for kind in KINDS {
        let starts = run(
            kerb,
            kind,
            &[
                (joint, SeamSide::Outgoing, 0.0),
                (90.0, SeamSide::Outgoing, 0.0),
            ],
        )
        .unwrap();
        let ends = run(
            kerb,
            kind,
            &[
                (5.0, SeamSide::Outgoing, 0.0),
                (joint, SeamSide::Outgoing, 0.0),
            ],
        )
        .unwrap();
        for side in SIDES {
            for near in [joint - NEAR, joint + NEAR] {
                let first = run(
                    kerb,
                    kind,
                    &[(near, side, WINDOW), (90.0, SeamSide::Outgoing, 0.0)],
                )
                .unwrap();
                same(&first, &starts, &format!("{kind:?} starts {side:?} {near}"));
                let last = run(
                    kerb,
                    kind,
                    &[(5.0, SeamSide::Outgoing, 0.0), (near, side, WINDOW)],
                )
                .unwrap();
                same(&last, &ends, &format!("{kind:?} ends {side:?} {near}"));
            }
        }
    }
}

/// A level straight whose cant jumps from 0 to 150 mm at 20: a seam whose
/// tangents agree but whose roll does not, so a run is not mitred there.
fn canted(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Curve3::Banked(Banked3::new(
        Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            ElevationLaw::constant_grade(40.0, 0.0),
        ),
        CantLaw::new(vec![
            CantPiece::constant(20.0, 0.0),
            CantPiece::constant(20.0, 0.15),
        ]),
        CantLaw::zero(40.0),
        1.5,
        BankConvention::TangentRotation,
    )))
    .unwrap()
}

/// An interior station on a seam that is not mitred stands in the frame
/// of the piece its side names, snapped or placed there; the two sides
/// differ.
#[test]
fn an_interior_station_on_a_roll_seam_reads_the_side_it_names() {
    let seam = 20.0;
    let mut b = GeometryGraphBuilder::new();
    let basis = canted(&mut b);
    let ids: Vec<NodeId> = SIDES
        .iter()
        .map(|&side| {
            b.push(GeometryNode::OrientedCurveStation(
                CurveStation::new(basis, Station::at(seam)).with_seam_side(side),
            ))
            .unwrap()
        })
        .collect();
    let graph = b.finish(ids.clone()).unwrap();
    let frames: Vec<_> = ids
        .iter()
        .map(|id| resolve(&graph, *id).unwrap().section)
        .collect();
    assert!(
        (frames[0].up - frames[1].up).length() > 0.05,
        "the roll jumps"
    );
    for kind in KINDS {
        let mut meshes = Vec::new();
        for (k, side) in SIDES.into_iter().enumerate() {
            let on = run(
                canted,
                kind,
                &[
                    (5.0, SeamSide::Outgoing, 0.0),
                    (seam, side, 0.0),
                    (35.0, SeamSide::Outgoing, 0.0),
                ],
            )
            .unwrap();
            for near in [seam - NEAR, seam + NEAR] {
                let snapped = run(
                    canted,
                    kind,
                    &[
                        (5.0, SeamSide::Outgoing, 0.0),
                        (near, side, WINDOW),
                        (35.0, SeamSide::Outgoing, 0.0),
                    ],
                )
                .unwrap();
                same(&snapped, &on, &format!("{kind:?} {side:?} at {near}"));
            }
            // The station's section stands in that side's frame: a corner
            // of the rectangle, an end of the deck, offsets `(0.5, 0.25)`.
            let local = match kind {
                Kind::Solid => Some((1.0, 0.5)),
                Kind::Sheet => Some((1.0, 0.0)),
                Kind::Offset => None,
            };
            if let Some((x, y)) = local {
                let corner = frames[k].place(0.5 + x, 0.25 + y, 0.0);
                let other = frames[1 - k].place(0.5 + x, 0.25 + y, 0.0);
                assert!(has_vertex(&on, corner), "{kind:?} {side:?}: {corner:?}");
                assert!(!has_vertex(&on, other), "{kind:?} {side:?}: {other:?}");
            }
            meshes.push(on);
        }
        assert_ne!(meshes[0].positions, meshes[1].positions, "{kind:?}");
    }
}

/// Without a window, and reading the outgoing piece, every run reads as
/// its window-free form does: the general forms default to today.
#[test]
fn a_window_free_run_reads_as_its_plain_form() {
    let (joint, _) = kerb_joint();
    let distances = [5.0, joint - 5.0, joint, 90.0];
    let general = |kind| {
        let stations: Vec<At> = distances
            .iter()
            .map(|&d| (d, SeamSide::Outgoing, 0.0))
            .collect();
        run(kerb, kind, &stations).unwrap()
    };
    let offsets = StationOffsets::new(0.5, 0.25, 0.0);
    // The plain forms: a stationed spine and an offset curve by stations.
    let mut b = GeometryGraphBuilder::new();
    let directrix = kerb(&mut b);
    let profile = rectangle(&mut b);
    let spine = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::StationedSpine {
                directrix,
                sections: distances
                    .iter()
                    .map(|&d| StationedSection {
                        profile,
                        station: Station::new(d, offsets),
                    })
                    .collect(),
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let offset = b
        .push(GeometryNode::CurveRelation(
            CurveRelation::OffsetByStations {
                basis: directrix,
                stations: distances
                    .iter()
                    .map(|&d| Station::new(d, offsets))
                    .collect(),
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let swept = b
        .push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix: offset,
            radius: 0.05,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        }))
        .unwrap();
    let graph = b.finish(vec![spine, swept]).unwrap();
    same(
        &general(Kind::Solid),
        &compile(&graph, spine).unwrap(),
        "spine",
    );
    same(
        &general(Kind::Offset),
        &compile(&graph, swept).unwrap(),
        "offset",
    );
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

/// Three collinear lines, 10 m, a sliver of `tiny`, the rest to 20 m: two
/// joints within one millimetre.
fn sliver(b: &mut GeometryGraphBuilder, tiny: Scalar) -> NodeId {
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let pieces = [
        trimmed(b, line, 0.0, 10.0),
        trimmed(b, line, 10.0, 10.0 + tiny),
        trimmed(b, line, 10.0 + tiny, 20.0),
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

#[test]
fn an_ambiguous_or_straddled_window_is_refused_by_name() {
    let tiny = 5e-4;
    let (_, certified) = kerb_joint();
    for kind in KINDS {
        let middle = 10.0 + 0.5 * tiny;
        let ambiguous = run(
            |b| sliver(b, tiny),
            kind,
            &[
                (2.0, SeamSide::Outgoing, 0.0),
                (middle, SeamSide::Incoming, WINDOW),
                (18.0, SeamSide::Outgoing, 0.0),
            ],
        );
        refused(ambiguous, SEAM_WINDOW_AMBIGUOUS);
        // The same run without a window reads as before.
        run(
            |b| sliver(b, tiny),
            kind,
            &[
                (2.0, SeamSide::Outgoing, 0.0),
                (middle, SeamSide::Incoming, 0.0),
                (18.0, SeamSide::Outgoing, 0.0),
            ],
        )
        .unwrap();
        for edge in [certified.midpoint() - WINDOW, certified.midpoint() + WINDOW] {
            let straddled = run(
                kerb,
                kind,
                &[
                    (5.0, SeamSide::Outgoing, 0.0),
                    (edge, SeamSide::Incoming, WINDOW),
                    (90.0, SeamSide::Outgoing, 0.0),
                ],
            );
            refused(straddled, SEAM_WINDOW_STRADDLED);
        }
    }
}

#[test]
fn stations_snapped_onto_one_joint_are_refused_by_name() {
    let (joint, _) = kerb_joint();
    for kind in KINDS {
        let collapsed = run(
            kerb,
            kind,
            &[
                (5.0, SeamSide::Outgoing, 0.0),
                (joint - NEAR, SeamSide::Incoming, WINDOW),
                (joint + NEAR, SeamSide::Outgoing, WINDOW),
                (90.0, SeamSide::Outgoing, 0.0),
            ],
        );
        match collapsed {
            Err(GeomError::InvalidInput(detail)) => {
                assert!(detail.contains("no longer increase strictly"), "{detail}");
            }
            other => panic!("{kind:?}: expected a refusal, got {other:?}"),
        }
        // A station snapped past an unwindowed neighbour.
        let crossed = run(
            kerb,
            kind,
            &[
                (5.0, SeamSide::Outgoing, 0.0),
                (joint - NEAR, SeamSide::Incoming, WINDOW),
                (joint - 0.5 * NEAR, SeamSide::Outgoing, 0.0),
                (90.0, SeamSide::Outgoing, 0.0),
            ],
        );
        assert!(
            matches!(crossed, Err(GeomError::InvalidInput(ref d)) if d.contains("no longer increase strictly")),
            "{kind:?}: {crossed:?}"
        );
    }
}

/// An offset curve by stations used as a station basis (#289) breaks its
/// law at the snapped distances: a station along it reads as along the
/// same offset with its station on the joint.
#[test]
fn an_offset_basis_breaks_its_law_at_the_snapped_stations() {
    let lines = |b: &mut GeometryGraphBuilder| {
        let line = b
            .push_value(Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }))
            .unwrap();
        let pieces = [trimmed(b, line, 0.0, 10.0), trimmed(b, line, 10.0, 20.0)];
        b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
            segments: pieces
                .iter()
                .map(|&curve| CurveSegment {
                    curve,
                    same_sense: true,
                    transition: Transition::Continuous,
                })
                .collect(),
        }))
        .unwrap()
    };
    let along = |middle: Scalar, window: Scalar| {
        let mut b = GeometryGraphBuilder::new();
        let basis = lines(&mut b);
        let offset = b
            .push(GeometryNode::CurveRelation(
                CurveRelation::OffsetsAtStations {
                    basis,
                    stations: vec![
                        RunStation::new(Station::new(2.0, StationOffsets::new(1.0, 0.0, 0.0))),
                        RunStation::new(Station::new(middle, StationOffsets::new(2.0, 0.0, 0.0)))
                            .with_seam_side(SeamSide::Incoming)
                            .with_seam_window(window),
                        RunStation::new(Station::new(18.0, StationOffsets::new(2.0, 0.0, 0.0))),
                    ],
                    frame: StationFrame::Section,
                },
            ))
            .unwrap();
        let ids: Vec<NodeId> = [3.0, 7.9, 8.0, 12.0]
            .iter()
            .map(|&d| {
                b.push_value(CurveStation::new(offset, Station::at(d)))
                    .unwrap()
            })
            .collect();
        let graph = b.finish(ids.clone()).unwrap();
        ids.iter()
            .map(|id| resolve(&graph, *id).unwrap().point)
            .collect::<Vec<_>>()
    };
    let on = along(10.0, 0.0);
    for near in [10.0 - NEAR, 10.0 + NEAR] {
        assert_eq!(along(near, WINDOW), on, "{near}");
        assert_ne!(along(near, 0.0), on, "{near} without a window");
    }
}
