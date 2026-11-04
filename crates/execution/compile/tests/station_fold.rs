//! Stationed runs whose section planes fold are refused (#361, ADR 0082
//! amendment): a section point on the inner side of a bend at or past its
//! centre of curvature, a vertex receding faster than the basis advances,
//! a touch within rounding. The mesh and the exact compiler refuse alike;
//! the closed forms are `F = 1 - a / R` on a plan arc of radius `R`,
//! `F = sigma (1 - b / R_v)` on a vertical arc in the section frame,
//! `F = 1 + c'` on a line with a longitudinal offset of rate `c'`. A turned
//! section and a frame read by quadrature are checked conservatively and
//! never refused: reported unbounded by name where the bound cannot tell.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame2, Point2, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    Circle2, CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Ellipse2, Intrinsic2, Line2,
};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId,
    SectionAtStation, SolidOperation, Station, StationFrame, StationOffsets, StationOrientation,
    Transition, TrimSelector, TrimmingPreference,
};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

const BUDGET: Scalar = 0.01;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(BUDGET)
        .unwrap()
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

/// What the mesh compiler makes of a run.
#[derive(Debug)]
enum Outcome {
    Certified,
    Unbounded(&'static str),
    Refused(String),
}

fn mesh(graph: &GeometryGraph, root: NodeId) -> Outcome {
    match compiler().compile_mesh_with_deviation(graph, root, &options()) {
        Ok((_, report)) => {
            let run = report
                .contributions
                .iter()
                .find(|c| c.path == DeviationPath::StationedSpine)
                .expect("the run reports its path");
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

/// `Ok` where the exact compiler builds the run, `Err` with the detail of
/// its invalid-input refusal.
fn exact(graph: &GeometryGraph, root: NodeId) -> Result<(), String> {
    match ReferenceExactCompiler::new().compile_exact(graph, root, &options()) {
        Ok(_) => Ok(()),
        Err(GeomError::InvalidInput(detail)) => Err(detail),
        Err(other) => panic!("not built nor refused as invalid input: {other:?}"),
    }
}

fn rectangle(b: &mut GeometryGraphBuilder) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x: 1.0,
        y: 1.0,
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

fn line_x() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })
}

fn at(distance: Scalar, offsets: StationOffsets) -> Station {
    Station::new(distance, offsets)
}

/// A 1 m by 1 m rectangle run along `basis` from distance 0 to `length`,
/// its offsets `from` at the start and `to` at the end.
fn run(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    frame: StationFrame,
    (from, to): (StationOffsets, StationOffsets),
    length: Scalar,
) -> NodeId {
    let profile = rectangle(b);
    run_with(b, basis, frame, profile, (from, to), length)
}

fn run_with(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    frame: StationFrame,
    profile: NodeId,
    (from, to): (StationOffsets, StationOffsets),
    length: Scalar,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix: basis,
            sections: vec![
                SectionAtStation::new(profile, at(0.0, from)),
                SectionAtStation::new(profile, at(length, to)),
            ],
            frame,
        },
    ))
    .unwrap()
}

/// One section offset sideways and up by a constant.
fn constant(lateral: Scalar, vertical: Scalar) -> (StationOffsets, StationOffsets) {
    let offsets = StationOffsets::new(lateral, vertical, 0.0);
    (offsets, offsets)
}

/// A run along a curve: its graph and root.
fn along(
    curve: impl Into<GeometryNode>,
    frame: StationFrame,
    offsets: (StationOffsets, StationOffsets),
    length: Scalar,
) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(curve.into()).unwrap();
    let root = run(&mut b, basis, frame, offsets, length);
    (b.finish(vec![root]).unwrap(), root)
}

fn refused(outcome: &Outcome) -> &str {
    match outcome {
        Outcome::Refused(detail) => {
            assert!(detail.contains("sections fold"), "{detail}");
            detail
        }
        other => panic!("not refused: {other:?}"),
    }
}

const FRAMES: [StationFrame; 2] = [StationFrame::Section, StationFrame::Plan];

#[test]
fn a_section_clear_of_a_plan_arcs_centre_is_certified_and_one_reaching_it_is_refused() {
    for frame in FRAMES {
        // A 1 m wide section centred 9 m to the left of a 10 m arc: its
        // inner edge stands 0.5 m from the centre of curvature.
        let (graph, root) = along(plan_circle(10.0), frame, constant(9.0, 0.0), 10.0);
        assert!(matches!(mesh(&graph, root), Outcome::Certified));
        assert!(exact(&graph, root).is_ok(), "{frame:?}");
        // The inner edge touches the centre (9.5), then passes it (10).
        for lateral in [9.5, 10.0] {
            let (graph, root) = along(plan_circle(10.0), frame, constant(lateral, 0.0), 10.0);
            // A witness at a distance, not an exhausted span.
            let detail = refused(&mesh(&graph, root)).to_owned();
            assert!(detail.contains("at distance"), "{detail}");
            let detail = exact(&graph, root).expect_err("the exact compiler refuses alike");
            assert!(detail.contains("sections fold"), "{detail}");
        }
        // The far side of the arc never folds.
        let (graph, root) = along(plan_circle(10.0), frame, constant(-50.0, 0.0), 10.0);
        assert!(matches!(mesh(&graph, root), Outcome::Certified));
        assert!(exact(&graph, root).is_ok());
    }
}

#[test]
fn a_reversed_arc_folds_its_right_side() {
    let mut b = GeometryGraphBuilder::new();
    let circle = b.push_value(plan_circle(10.0)).unwrap();
    let arc = b
        .push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis: circle,
            start: vec![TrimSelector::Parameter(0.0)],
            end: vec![TrimSelector::Parameter(1.5)],
            sense_agreement: true,
            preference: TrimmingPreference::Parameter,
        }))
        .unwrap();
    let basis = b
        .push(GeometryNode::CurveRelation(CurveRelation::Composite {
            segments: vec![CurveSegment {
                curve: arc,
                same_sense: false,
                transition: Transition::Continuous,
            }],
        }))
        .unwrap();
    let profile = rectangle(&mut b);
    let runs: Vec<NodeId> = [-9.5, 9.5]
        .iter()
        .map(|&lateral| {
            let (from, to) = constant(lateral, 0.0);
            run_with(
                &mut b,
                basis,
                StationFrame::Section,
                profile,
                (from, to),
                10.0,
            )
        })
        .collect();
    let graph = b.finish(runs.clone()).unwrap();
    refused(&mesh(&graph, runs[0]));
    assert!(matches!(mesh(&graph, runs[1]), Outcome::Certified));
}

/// A line along `x` under a vertical profile.
fn graded(elevation: ElevationLaw) -> Curve3 {
    Curve3::Elevated(Elevated3::new(line_x(), elevation))
}

#[test]
fn a_section_reaching_a_vertical_arcs_radius_folds_in_the_section_frame_only() {
    let sags = [
        // z = d^2 / 100: a radius of 50 m at its start.
        ElevationLaw::Polynomial {
            coefficients: vec![0.0, 0.0, 0.01],
        },
        ElevationLaw::circular_arc(0.0, 0.0, 50.0),
    ];
    for sag in sags {
        let up = |vertical| {
            along(
                graded(sag.clone()),
                StationFrame::Section,
                constant(0.0, vertical),
                10.0,
            )
        };
        let (graph, root) = up(40.0);
        assert!(matches!(mesh(&graph, root), Outcome::Certified));
        let (graph, root) = up(55.0);
        refused(&mesh(&graph, root));
        // The plan frame's up never turns toward the horizontal tangent.
        let (graph, root) = along(
            graded(sag.clone()),
            StationFrame::Plan,
            constant(0.0, 55.0),
            10.0,
        );
        assert!(matches!(mesh(&graph, root), Outcome::Certified));
    }
}

#[test]
fn a_section_below_a_crest_folds_alike() {
    let crests = [
        ElevationLaw::Polynomial {
            coefficients: vec![0.0, 0.0, -0.01],
        },
        ElevationLaw::circular_arc(0.0, 0.0, -50.0),
    ];
    for crest in crests {
        let (graph, root) = along(
            graded(crest.clone()),
            StationFrame::Section,
            constant(0.0, -55.0),
            10.0,
        );
        refused(&mesh(&graph, root));
        let (graph, root) = along(
            graded(crest),
            StationFrame::Section,
            constant(0.0, 55.0),
            10.0,
        );
        assert!(matches!(mesh(&graph, root), Outcome::Certified));
    }
}

#[test]
fn a_fold_inside_a_transition_is_found_between_the_stations() {
    // A clothoid with `A^2 = 1000` over 100 m: its curvature `s / 1000`
    // meets a lateral offset falling linearly to zero inside the run.
    let clothoid = || {
        Curve2::Intrinsic(Intrinsic2::new(
            Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            CurvatureLaw::clothoid(0.0, 0.1, 100.0),
            100.0,
        ))
    };
    let falling = |from: Scalar| {
        (
            StationOffsets::new(from, 0.0, 0.0),
            StationOffsets::new(0.0, 0.0, 0.0),
        )
    };
    let (graph, root) = along(clothoid(), StationFrame::Section, falling(30.0), 100.0);
    assert!(matches!(mesh(&graph, root), Outcome::Certified));
    let (graph, root) = along(clothoid(), StationFrame::Section, falling(50.0), 100.0);
    refused(&mesh(&graph, root));
}

#[test]
fn a_section_receding_faster_than_the_basis_advances_folds() {
    let receding = |end: Scalar| {
        (
            StationOffsets::new(0.0, 0.0, 0.0),
            StationOffsets::new(0.0, 0.0, end),
        )
    };
    let (graph, root) = along(line_x(), StationFrame::Section, receding(-20.0), 10.0);
    refused(&mesh(&graph, root));
    let detail = exact(&graph, root).expect_err("the exact compiler refuses alike");
    assert!(detail.contains("sections fold"), "{detail}");
    let (graph, root) = along(line_x(), StationFrame::Section, receding(-5.0), 10.0);
    assert!(matches!(mesh(&graph, root), Outcome::Certified));
}

#[test]
fn a_touch_at_the_end_of_a_span_is_refused_after_halving() {
    // The lateral offset grows to 9.5 at the far station: the inner edge
    // reaches the centre of curvature there and nowhere before.
    let growing = (
        StationOffsets::new(0.0, 0.0, 0.0),
        StationOffsets::new(9.5, 0.0, 0.0),
    );
    let (graph, root) = along(plan_circle(10.0), StationFrame::Section, growing, 10.0);
    let detail = refused(&mesh(&graph, root)).to_owned();
    assert!(detail.contains("between distances"), "{detail}");
}

#[test]
fn a_turned_section_is_checked_conservatively_and_never_refused() {
    let turn = StationOrientation::new(None, Some(Vec3::new(0.98, 0.2, 0.0)));
    let turned = |lateral: Scalar| {
        let mut b = GeometryGraphBuilder::new();
        let basis = b.push_value(plan_circle(10.0)).unwrap();
        let profile = rectangle(&mut b);
        let offsets = StationOffsets::new(lateral, 0.0, 0.0);
        let root = b
            .push(GeometryNode::SolidOperation(
                SolidOperation::SectionsAtStations {
                    directrix: basis,
                    sections: vec![
                        SectionAtStation::new(profile, at(0.0, offsets)),
                        SectionAtStation::new(profile, at(10.0, offsets)).with_orientation(turn),
                    ],
                    frame: StationFrame::Section,
                },
            ))
            .unwrap();
        (b.finish(vec![root]).unwrap(), root)
    };
    let (graph, root) = turned(3.0);
    assert!(matches!(mesh(&graph, root), Outcome::Certified));
    // Past the centre: the conservative bound cannot clear it, and does not
    // refuse.
    let (graph, root) = turned(11.0);
    match mesh(&graph, root) {
        Outcome::Unbounded(reason) => assert!(reason.contains("cannot be excluded"), "{reason}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_frame_read_by_quadrature_is_checked_conservatively_and_never_refused() {
    let ellipse = || {
        Curve2::Ellipse(Ellipse2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            semi_axis_x: 20.0,
            semi_axis_y: 10.0,
        })
    };
    // Its sharpest bend is `10^2 / 20 = 5` m at the ends of the major axis.
    let (graph, root) = along(ellipse(), StationFrame::Section, constant(1.5, 0.0), 10.0);
    assert!(matches!(mesh(&graph, root), Outcome::Certified));
    let (graph, root) = along(ellipse(), StationFrame::Section, constant(7.5, 0.0), 10.0);
    match mesh(&graph, root) {
        Outcome::Unbounded(reason) => assert!(reason.contains("cannot be excluded"), "{reason}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_curved_profile_is_judged_by_its_vertices_and_its_flattening_bound() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(plan_circle(10.0)).unwrap();
    let disk = b
        .push_value(Profile::Circle(CircleProfile {
            radius: 1.0,
            thickness: None,
        }))
        .unwrap();
    let clear = run_with(
        &mut b,
        basis,
        StationFrame::Section,
        disk,
        constant(8.0, 0.0),
        10.0,
    );
    let reaching = run_with(
        &mut b,
        basis,
        StationFrame::Section,
        disk,
        constant(9.2, 0.0),
        10.0,
    );
    let graph = b.finish(vec![clear, reaching]).unwrap();
    assert!(matches!(mesh(&graph, clear), Outcome::Certified));
    refused(&mesh(&graph, reaching));
    let detail = exact(&graph, reaching).expect_err("the exact compiler refuses alike");
    assert!(detail.contains("sections fold"), "{detail}");
}

fn circle3(radius: Scalar) -> Curve3 {
    Curve3::Circle(axiolid_curve::Circle3 {
        frame: axiolid_core::Frame3 {
            origin: axiolid_core::Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    })
}

#[test]
fn a_station_interval_longer_than_twice_the_radius_is_still_decided() {
    // The bend is not bounded over a span longer than twice the radius, but
    // halving the span decides it: neither a fold gets through nor a clear
    // run is lost.
    fn arc() -> GeometryNode {
        Curve3::Elevated(Elevated3::new(plan_circle(10.0), ElevationLaw::level(0.0))).into()
    }
    fn space_circle() -> GeometryNode {
        circle3(10.0).into()
    }
    type Basis = (&'static str, fn() -> GeometryNode);
    let curves: [Basis; 2] = [("space circle", space_circle), ("elevated arc", arc)];
    for (name, curve) in curves {
        for frame in FRAMES {
            let (graph, root) = along(curve(), frame, constant(10.0, 0.0), 30.0);
            refused(&mesh(&graph, root));
            let (graph, root) = along(curve(), frame, constant(9.0, 0.0), 30.0);
            assert!(
                matches!(mesh(&graph, root), Outcome::Certified),
                "{name} {frame:?}"
            );
        }
    }
}

#[test]
fn a_fold_after_a_joint_inside_a_station_interval_is_found() {
    // A leg of 10 m, then a 10 m arc, one interval of 20 m across the
    // joint: the section reaches the arc's centre just after the joint.
    let mut b = GeometryGraphBuilder::new();
    let line = b.push_value(line_x()).unwrap();
    let trim = |b: &mut GeometryGraphBuilder, basis, end| {
        b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start: vec![TrimSelector::Parameter(0.0)],
            end: vec![TrimSelector::Parameter(end)],
            sense_agreement: true,
            preference: TrimmingPreference::Parameter,
        }))
        .unwrap()
    };
    let leg = trim(&mut b, line, 10.0);
    let circle = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(10.0, 10.0),
                x: Vec2::new(0.0, -1.0),
                y: Vec2::new(1.0, 0.0),
            },
            radius: 10.0,
        }))
        .unwrap();
    let arc = trim(&mut b, circle, 1.0);
    let segment = |curve| CurveSegment {
        curve,
        same_sense: true,
        transition: Transition::Continuous,
    };
    let basis = b
        .push(GeometryNode::CurveRelation(CurveRelation::Composite {
            segments: vec![segment(leg), segment(arc)],
        }))
        .unwrap();
    let profile = rectangle(&mut b);
    let roots: Vec<NodeId> = [19.0, 19.2]
        .iter()
        .map(|&lateral| {
            let offsets = (StationOffsets::new(lateral, 0.0, 0.0), StationOffsets::ZERO);
            run_with(&mut b, basis, StationFrame::Section, profile, offsets, 20.0)
        })
        .collect();
    let graph = b.finish(roots.clone()).unwrap();
    for root in roots {
        refused(&mesh(&graph, root));
    }
}

#[test]
fn the_exact_compiler_refuses_a_fold_it_cannot_exclude_by_name() {
    // A turned section along a line whose longitudinal offset recedes: the
    // conservative bound cannot clear it, so the exact compiler does not
    // build what the mesh compiler cannot vouch for.
    let turn = StationOrientation::new(None, Some(Vec3::new(0.996, 0.087, 0.0)));
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(line_x()).unwrap();
    let profile = rectangle(&mut b);
    let section = |s: Scalar, c: Scalar| {
        SectionAtStation::new(profile, at(s, StationOffsets::new(0.0, 0.0, c)))
            .with_orientation(turn)
    };
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: vec![
                    section(0.0, 0.0),
                    section(10.0, -20.0),
                    section(20.0, -20.0),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    match ReferenceExactCompiler::new().compile_exact(&graph, root, &options()) {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert!(input.contains("cannot be excluded"), "{input}");
        }
        other => panic!("not refused by name: {other:?}"),
    }
}

fn pitched(degrees: Scalar) -> StationOrientation {
    let t = degrees.to_radians();
    StationOrientation::new(
        Some(Vec3::new(-t.sin(), 0.0, t.cos())),
        Some(Vec3::new(t.cos(), 0.0, t.sin())),
    )
}

/// 1 m rectangles at 0 and 10 along `curve`, both pitched.
fn pitched_run(
    curve: impl Into<GeometryNode>,
    degrees: Scalar,
    frame: StationFrame,
) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(curve.into()).unwrap();
    let profile = rectangle(&mut b);
    let section = |s: Scalar| {
        SectionAtStation::new(profile, at(s, StationOffsets::ZERO))
            .with_orientation(pitched(degrees))
    };
    let root = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: vec![section(0.0), section(10.0)],
                frame,
            },
        ))
        .unwrap();
    (b.finish(vec![root]).unwrap(), root)
}

#[test]
fn a_section_pitched_steeply_along_a_flat_basis_is_not_a_fold() {
    // In the plan frame `O'` is horizontal on a flat basis: no vertical part
    // of it tips a section pitched past 45 degrees over.
    for degrees in [30.0, 46.0, 60.0] {
        let (graph, root) = pitched_run(line_x(), degrees, StationFrame::Plan);
        assert!(
            matches!(mesh(&graph, root), Outcome::Certified),
            "{degrees}"
        );
        assert!(exact(&graph, root).is_ok(), "{degrees}");
        // Along a grade the plan frame leaves the climb out; still decided.
        let graded = Curve3::Elevated(Elevated3::new(
            line_x(),
            ElevationLaw::Polynomial {
                coefficients: vec![0.0, 0.05],
            },
        ));
        let (graph, root) = pitched_run(graded, degrees, StationFrame::Plan);
        assert!(
            matches!(mesh(&graph, root), Outcome::Certified),
            "graded {degrees}"
        );
    }
}

#[test]
fn a_basis_the_exact_compiler_cannot_build_is_refused_as_such_first() {
    let ellipse = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 20.0,
        semi_axis_y: 10.0,
    });
    // Past the sharpest bend, where the fold cannot be excluded.
    let (graph, root) = along(ellipse, StationFrame::Section, constant(7.5, 0.0), 10.0);
    match ReferenceExactCompiler::new().compile_exact(&graph, root, &options()) {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert!(input.contains("other than a line"), "{input}");
        }
        other => panic!("not refused as an unbuildable basis: {other:?}"),
    }
}

#[test]
fn paired_curved_profiles_are_judged_with_their_pairing_bound() {
    // A disk growing from 0.5 to 0.8 m along a 10 m arc, the section
    // centred 9.19999 m to the left: the exact boundary reaches 9.99999 m,
    // short of the centre of curvature, but the flattened vertices plus the
    // pairing bound (the sagitta of the chords, about 0.005 m) reach it, so
    // the run is refused as it is for a single curved profile. 8 m is clear.
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(plan_circle(10.0)).unwrap();
    let disk = |b: &mut GeometryGraphBuilder, radius| {
        b.push_value(Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        }))
        .unwrap()
    };
    let (small, large) = (disk(&mut b, 0.5), disk(&mut b, 0.8));
    let mut grow = |lateral: Scalar| {
        let offsets = StationOffsets::new(lateral, 0.0, 0.0);
        b.push(GeometryNode::SolidOperation(
            SolidOperation::SectionsAtStations {
                directrix: basis,
                sections: vec![
                    SectionAtStation::new(small, at(0.0, offsets)),
                    SectionAtStation::new(large, at(10.0, offsets)),
                ],
                frame: StationFrame::Section,
            },
        ))
        .unwrap()
    };
    let (clear, reaching) = (grow(8.0), grow(9.19999));
    let graph = b.finish(vec![clear, reaching]).unwrap();
    assert!(matches!(mesh(&graph, clear), Outcome::Certified));
    refused(&mesh(&graph, reaching));
    // The exact compiler builds neither: a fold is refused as such before
    // the curved section that changes.
    let detail = exact(&graph, reaching).expect_err("the exact compiler refuses");
    assert!(detail.contains("sections fold"), "{detail}");
}
