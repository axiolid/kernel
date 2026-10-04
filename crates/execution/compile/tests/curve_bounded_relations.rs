//! Curve-bounded planes whose boundaries are curve relations (#255):
//! composites of polylines, lines and arcs, and trims of a basis curve, as
//! IFC `IfcCompositeCurve` and `IfcTrimmedCurve` space boundaries are.
//! Oracles are the same boundary written as an atomic curve, and
//! closed-form areas.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame2, Frame3, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    Circle2, CurvatureLaw, Curve2, Curve3, Intrinsic2, Intrinsic3, Line2, Polyline2, Polyline3,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::{CompileOutcome, MeshClosure};
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraphBuilder, GeometryNode, NodeId, SurfaceRelation,
    Transition, TrimSelector, TrimmingPreference,
};
use axiolid_surface::{Plane, Surface};

const CHORD: Scalar = 1e-4;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(CHORD)
        .expect("positive chord budget")
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

/// Compile a curve-bounded plane over the world XY frame whose boundary
/// nodes `build` pushes.
fn compile(
    build: impl FnOnce(&mut GeometryGraphBuilder) -> Vec<NodeId>,
) -> Result<(CompileOutcome, DeviationReport), GeomError> {
    compile_on(
        Frame3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        build,
    )
}

/// [`compile`] over the plane frame `frame`.
fn compile_on(
    frame: Frame3,
    build: impl FnOnce(&mut GeometryGraphBuilder) -> Vec<NodeId>,
) -> Result<(CompileOutcome, DeviationReport), GeomError> {
    let mut b = GeometryGraphBuilder::new();
    let basis = b
        .push(GeometryNode::Surface(Surface::Plane(Plane { frame })))
        .unwrap();
    let boundaries = build(&mut b);
    let root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::CurveBounded {
                basis,
                boundaries,
                implicit_outer: false,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![root]).unwrap();
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh_with_deviation(
        &graph,
        root,
        &options(),
    )
}

/// The bound the report gives the curve-bounded plane.
fn plane_bound(report: &DeviationReport) -> DeviationBound {
    let [contribution] = report.contributions.as_slice() else {
        panic!("{report:?}");
    };
    assert_eq!(contribution.path, DeviationPath::CurveBoundedPlane);
    contribution.bound
}

/// Proven at exactly the chord budget, so the report meets it: the
/// rounding-level gap a curve's computed end leaves at a joint or where the
/// loop closes (a full turn's `sin(2 pi)`, an arc's `cos(pi / 2)`) is welded,
/// not counted.
fn assert_chord_bound(report: &DeviationReport) {
    assert_eq!(plane_bound(report), DeviationBound::Proven(CHORD));
    assert!(report.meets_requested(), "{report:?}");
}

#[test]
fn an_atomic_circle_boundary_meets_its_budget() {
    // The full turn's closing point lands `sin(2 pi) r` from its start.
    let r = 1.5;
    let (outcome, report) = compile(|b| vec![push(b, circle2((2.0, 1.5), r))]).unwrap();
    let (got, exact) = (area(&outcome.mesh), PI * r * r);
    assert!(
        got <= exact && exact - got <= TAU * r * CHORD,
        "{got} vs {exact}"
    );
    assert_chord_bound(&report);
}

#[test]
fn a_rotated_plane_still_meets_its_budget() {
    // Axes unit and orthogonal only to rounding: the frame's stretch bound
    // is a hair above one, and the boundaries are flattened to the budget
    // shrunk by it, so the reported bound stays within the request.
    let (s, c) = 0.7_f64.sin_cos();
    let x = Vec3::new(c, s, 0.0);
    let y = Vec3::new(-s * 0.6, c * 0.6, 0.8);
    let frame = Frame3 {
        origin: Point3::new(10.0, -4.0, 2.5),
        x,
        y,
        z: x.cross(y),
    };
    let (outcome, report) = compile_on(frame, |b| vec![stadium(b)]).unwrap();
    let got = area(&outcome.mesh);
    assert!(
        got <= STADIUM_AREA + 1e-9 && STADIUM_AREA - got <= STADIUM_PERIMETER * CHORD,
        "{got} vs {STADIUM_AREA}"
    );
    match plane_bound(&report) {
        DeviationBound::Proven(d) => assert!(d <= CHORD && d > 0.99 * CHORD, "{d}"),
        other => panic!("{other:?}"),
    }
    assert!(report.meets_requested(), "{report:?}");
}

fn push(b: &mut GeometryGraphBuilder, node: GeometryNode) -> NodeId {
    b.push(node).unwrap()
}

fn polyline2(points: &[(Scalar, Scalar)], closed: bool) -> GeometryNode {
    GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        closed,
    }))
}

fn polyline3(points: &[(Scalar, Scalar, Scalar)], closed: bool) -> GeometryNode {
    GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
        points: points
            .iter()
            .map(|&(x, y, z)| Point3::new(x, y, z))
            .collect(),
        closed,
    }))
}

fn circle2(centre: (Scalar, Scalar), radius: Scalar) -> GeometryNode {
    GeometryNode::Curve2(Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(centre.0, centre.1),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    }))
}

fn line2(origin: (Scalar, Scalar), direction: (Scalar, Scalar)) -> GeometryNode {
    GeometryNode::Curve2(Curve2::Line(Line2 {
        origin: Point2::new(origin.0, origin.1),
        direction: Vec2::new(direction.0, direction.1),
    }))
}

fn composite(b: &mut GeometryGraphBuilder, segments: &[(NodeId, bool)]) -> NodeId {
    push(
        b,
        GeometryNode::CurveRelation(CurveRelation::Composite {
            segments: segments
                .iter()
                .map(|&(curve, same_sense)| CurveSegment {
                    curve,
                    same_sense,
                    transition: Transition::Continuous,
                })
                .collect(),
        }),
    )
}

fn trimmed(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    start: TrimSelector,
    end: TrimSelector,
    sense_agreement: bool,
    preference: TrimmingPreference,
) -> NodeId {
    push(
        b,
        GeometryNode::CurveRelation(CurveRelation::Trimmed {
            basis,
            start: vec![start],
            end: vec![end],
            sense_agreement,
            preference,
        }),
    )
}

fn by_parameter(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    start: Scalar,
    end: Scalar,
    sense_agreement: bool,
) -> NodeId {
    trimmed(
        b,
        basis,
        TrimSelector::Parameter(start),
        TrimSelector::Parameter(end),
        sense_agreement,
        TrimmingPreference::Parameter,
    )
}

/// A non-convex L, counter-clockwise: area 4 * 3 - 2 * 2 = 8.
const L_SHAPE: [(Scalar, Scalar); 6] = [
    (0.0, 0.0),
    (4.0, 0.0),
    (4.0, 1.0),
    (2.0, 1.0),
    (2.0, 3.0),
    (0.0, 3.0),
];

/// The consumer's done-when: a one-segment composite wrapping a closed
/// polyline, against its sense, is the polyline's own surface.
#[test]
fn a_reversed_one_segment_composite_is_the_polyline_alone() {
    let (alone, alone_report) = compile(|b| vec![push(b, polyline2(&L_SHAPE, true))]).unwrap();
    assert_eq!(area(&alone.mesh), 8.0);
    for same_sense in [false, true] {
        let (wrapped, report) = compile(|b| {
            let ring = push(b, polyline2(&L_SHAPE, true));
            vec![composite(b, &[(ring, same_sense)])]
        })
        .unwrap();
        assert_eq!(wrapped.closure, MeshClosure::Surface);
        assert_eq!(area(&wrapped.mesh), 8.0, "same sense {same_sense}");
        // The same triangles, vertex for vertex, facing the same way.
        assert_eq!(wrapped.mesh, alone.mesh, "same sense {same_sense}");
        assert_eq!(plane_bound(&report), plane_bound(&alone_report));
        assert_eq!(plane_bound(&report), DeviationBound::Proven(0.0));
    }
    // A 3D polyline on the parameter plane reads the same.
    let flat: Vec<_> = L_SHAPE.iter().map(|&(x, y)| (x, y, 0.0)).collect();
    let (wrapped, _) = compile(|b| {
        let ring = push(b, polyline3(&flat, true));
        vec![composite(b, &[(ring, false)])]
    })
    .unwrap();
    assert_eq!(wrapped.mesh, alone.mesh);
}

#[test]
fn reversed_composite_holes_cut_the_same_area() {
    let outer = [(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)];
    let hole = [(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)];
    let (outcome, report) = compile(|b| {
        let o = push(b, polyline2(&outer, true));
        let h = push(b, polyline2(&hole, true));
        vec![composite(b, &[(o, false)]), composite(b, &[(h, false)])]
    })
    .unwrap();
    assert_eq!(area(&outcome.mesh), 11.0);
    for t in outcome.mesh.indices.chunks_exact(3) {
        let [a, b, c] = [0, 1, 2].map(|k| outcome.mesh.positions[t[k] as usize]);
        assert!((b - a).cross(c - a).z > 0.0);
    }
    assert_eq!(plane_bound(&report), DeviationBound::Proven(0.0));
}

/// A stadium of straight length 4 and radius 1, its four segments each
/// read a different way: a trimmed line, an arc with its circle, a
/// polyline against its sense, and an arc trimmed against its circle and
/// then run against that trim. Area `8 + pi`.
fn stadium(b: &mut GeometryGraphBuilder) -> NodeId {
    // (0, -1) -> (4, -1): a line at speed 2, parameters 0 to 2.
    let bottom_line = push(b, line2((0.0, -1.0), (2.0, 0.0)));
    let bottom = by_parameter(b, bottom_line, 0.0, 2.0, true);
    // (4, -1) -> (5, 0) -> (4, 1), counter-clockwise.
    let right_circle = push(b, circle2((4.0, 0.0), 1.0));
    let right = by_parameter(b, right_circle, -FRAC_PI_2, FRAC_PI_2, true);
    // Written (0, 1) -> (4, 1), run backwards.
    let top = push(b, polyline2(&[(0.0, 1.0), (4.0, 1.0)], false));
    // Trimmed against the circle, 3pi/2 down to pi/2: (0, -1) -> (-1, 0)
    // -> (0, 1); run backwards it closes the loop at (0, -1).
    let left_circle = push(b, circle2((0.0, 0.0), 1.0));
    let left = by_parameter(b, left_circle, 3.0 * FRAC_PI_2, FRAC_PI_2, false);
    composite(
        b,
        &[(bottom, true), (right, true), (top, false), (left, false)],
    )
}

const STADIUM_AREA: Scalar = 8.0 + PI;
const STADIUM_PERIMETER: Scalar = 8.0 + TAU;

#[test]
fn composites_of_lines_and_arcs_in_both_senses_are_certified() {
    let (outcome, report) = compile(|b| vec![stadium(b)]).unwrap();
    let got = area(&outcome.mesh);
    // Chords cut inside the arcs: smaller, by at most perimeter x budget.
    assert!(
        got <= STADIUM_AREA && STADIUM_AREA - got <= STADIUM_PERIMETER * CHORD,
        "{got} vs {STADIUM_AREA}"
    );
    // Every vertex on the exact boundary: on a straight edge or an arc.
    for p in &outcome.mesh.positions {
        let on_line = (p.y.abs() - 1.0).abs() < 1e-12 && (-1e-12..=4.0 + 1e-12).contains(&p.x);
        let on_arc = [0.0, 4.0]
            .iter()
            .any(|&cx| ((p.x - cx).hypot(p.y) - 1.0).abs() < 1e-12);
        assert!(on_line || on_arc, "{p:?}");
    }
    assert_chord_bound(&report);
    // As a hole in a 6 x 4 plate, the stadium cuts its own area.
    let (outcome, report) = compile(|b| {
        let plate = push(
            b,
            polyline2(&[(-1.5, -2.0), (5.5, -2.0), (5.5, 2.0), (-1.5, 2.0)], true),
        );
        vec![plate, stadium(b)]
    })
    .unwrap();
    let got = area(&outcome.mesh);
    let exact = 28.0 - STADIUM_AREA;
    assert!(
        got >= exact && got - exact <= STADIUM_PERIMETER * CHORD,
        "{got} vs {exact}"
    );
    assert_chord_bound(&report);
}

#[test]
fn a_joint_gap_within_tolerance_is_in_the_bound() {
    // A square whose second side starts 5e-4 above where the first ends:
    // stitched shut (the linear tolerance is 1e-3), and reported.
    let (outcome, report) = compile(|b| {
        let first = push(b, polyline2(&[(0.0, 0.0), (4.0, 0.0)], false));
        let rest = push(
            b,
            polyline2(&[(4.0, 5e-4), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)], false),
        );
        vec![composite(b, &[(first, true), (rest, true)])]
    })
    .unwrap();
    assert!((area(&outcome.mesh) - 16.0).abs() < 1e-2);
    match plane_bound(&report) {
        DeviationBound::Proven(d) => assert!((d - 5e-4).abs() < 1e-12, "{d}"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_full_turn_trimmed_by_parameter_is_its_disk() {
    let r = 1.5;
    let (outcome, report) = compile(|b| {
        let circle = push(b, circle2((2.0, 1.5), r));
        vec![by_parameter(b, circle, 0.0, TAU, true)]
    })
    .unwrap();
    let exact = PI * r * r;
    let got = area(&outcome.mesh);
    assert!(
        got <= exact && exact - got <= TAU * r * CHORD,
        "{got} vs {exact}"
    );
    assert_chord_bound(&report);
    // A closed polyline trimmed over its whole parameter range.
    let (outcome, report) = compile(|b| {
        let ring = push(
            b,
            polyline2(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)], true),
        );
        vec![by_parameter(b, ring, 0.0, 4.0, true)]
    })
    .unwrap();
    assert_eq!(area(&outcome.mesh), 12.0);
    assert_eq!(plane_bound(&report), DeviationBound::Proven(0.0));
}

#[test]
fn half_circles_trimmed_by_points_close_a_disk() {
    let r = 1.0;
    let (outcome, report) = compile(|b| {
        let circle = push(b, circle2((0.0, 0.0), r));
        let point = |x: Scalar| TrimSelector::Point2(Point2::new(x, 0.0));
        let upper = trimmed(
            b,
            circle,
            point(r),
            point(-r),
            true,
            TrimmingPreference::Cartesian,
        );
        let lower = trimmed(
            b,
            circle,
            point(-r),
            point(r),
            true,
            TrimmingPreference::Cartesian,
        );
        vec![composite(b, &[(upper, true), (lower, true)])]
    })
    .unwrap();
    let got = area(&outcome.mesh);
    assert!(got <= PI && PI - got <= TAU * r * CHORD, "{got}");
    assert_chord_bound(&report);
}

#[test]
fn arc_length_trims_resolve_on_lines_and_circles() {
    // A 3-4-5 triangle of lines at speeds 2, 5 and 1, each trimmed by arc
    // length: exact, area 6.
    let (outcome, report) = compile(|b| {
        let mut side = |origin, direction, length| {
            let line = push(b, line2(origin, direction));
            trimmed(
                b,
                line,
                TrimSelector::ArcLength(0.0),
                TrimSelector::ArcLength(length),
                true,
                TrimmingPreference::Parameter,
            )
        };
        let a = side((0.0, 0.0), (2.0, 0.0), 4.0);
        let c = side((4.0, 0.0), (-4.0, 3.0), 5.0);
        let d = side((0.0, 3.0), (0.0, -1.0), 3.0);
        vec![composite(b, &[(a, true), (c, true), (d, true)])]
    })
    .unwrap();
    assert!((area(&outcome.mesh) - 6.0).abs() < 1e-12);
    assert_eq!(plane_bound(&report), DeviationBound::Proven(0.0));
    // A circle of radius 2 trimmed by arc length over its full length 4pi,
    // and against its sense: the same disk, area 4pi.
    for sense in [true, false] {
        let (outcome, report) = compile(|b| {
            let circle = push(b, circle2((0.0, 0.0), 2.0));
            vec![trimmed(
                b,
                circle,
                TrimSelector::ArcLength(0.0),
                TrimSelector::ArcLength(2.0 * TAU),
                sense,
                TrimmingPreference::Parameter,
            )]
        })
        .unwrap();
        let (got, exact) = (area(&outcome.mesh), 4.0 * PI);
        assert!(
            got <= exact + 1e-9 && exact - got <= 2.0 * TAU * CHORD,
            "sense {sense}: {got} vs {exact}"
        );
        assert_chord_bound(&report);
    }
}

#[test]
fn an_uncertified_leaf_is_unbounded_by_name() {
    // A planar natural-equation circle: compiled, but its flattening is not
    // certified, so the report says so instead of a number.
    let r = 1.0;
    let (outcome, report) = compile(|b| {
        let circle = push(
            b,
            GeometryNode::Curve3(Curve3::Intrinsic(Intrinsic3::new(
                Frame3 {
                    origin: Point3::new(0.0, -r, 0.0),
                    x: Vec3::X,
                    y: Vec3::Y,
                    z: Vec3::Z,
                },
                CurvatureLaw::circular(1.0 / r),
                CurvatureLaw::circular(0.0),
                TAU * r,
            ))),
        );
        vec![composite(b, &[(circle, true)])]
    })
    .unwrap();
    assert!((area(&outcome.mesh) - PI).abs() < 1e-2);
    assert_eq!(
        plane_bound(&report),
        DeviationBound::Unbounded("curve-bounded plane boundary family")
    );
    assert_eq!(report.bound, None);
}

fn refusal(build: impl FnOnce(&mut GeometryGraphBuilder) -> Vec<NodeId>) -> GeomError {
    match compile(build) {
        Ok((outcome, _)) => panic!("compiled {} triangles", outcome.mesh.indices.len() / 3),
        Err(error) => error,
    }
}

#[test]
fn open_and_gapped_composites_are_refused_by_name() {
    // Two sides of a square: does not close.
    let error = refusal(|b| {
        let first = push(b, polyline2(&[(0.0, 0.0), (4.0, 0.0)], false));
        let second = push(b, polyline2(&[(4.0, 0.0), (4.0, 4.0)], false));
        vec![composite(b, &[(first, true), (second, true)])]
    });
    assert!(
        matches!(&error, GeomError::InvalidInput(m) if m.contains("open curve relation")),
        "{error:?}"
    );
    // A segment run the wrong way leaves a gap.
    let error = refusal(|b| {
        let first = push(b, polyline2(&[(0.0, 0.0), (4.0, 0.0)], false));
        let rest = push(
            b,
            polyline2(&[(4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)], false),
        );
        vec![composite(b, &[(first, true), (rest, false)])]
    });
    assert!(
        matches!(&error, GeomError::InvalidInput(m) if m.contains("gap")),
        "{error:?}"
    );
    // One half circle alone.
    let error = refusal(|b| {
        let circle = push(b, circle2((0.0, 0.0), 1.0));
        vec![by_parameter(b, circle, 0.0, PI, true)]
    });
    assert!(
        matches!(&error, GeomError::InvalidInput(m) if m.contains("open curve relation")),
        "{error:?}"
    );
}

#[test]
fn unsupported_relations_and_families_are_refused_by_name() {
    let unsupported = |error: GeomError, name: &str| match error {
        GeomError::UnsupportedInput { input, .. } => assert!(input.contains(name), "{input}"),
        other => panic!("{other:?}"),
    };
    // An offset curve, alone or inside a composite.
    let offset = |b: &mut GeometryGraphBuilder| {
        let ring = push(
            b,
            polyline2(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)], true),
        );
        push(
            b,
            GeometryNode::CurveRelation(CurveRelation::Offset {
                basis: ring,
                distance: 0.5,
                reference_direction: None,
            }),
        )
    };
    unsupported(refusal(|b| vec![offset(b)]), "offset");
    unsupported(
        refusal(|b| {
            let o = offset(b);
            vec![composite(b, &[(o, true)])]
        }),
        "offset",
    );
    // A 2D family with no 3D twin.
    unsupported(
        refusal(|b| {
            let spiral = push(
                b,
                GeometryNode::Curve2(Curve2::Intrinsic(Intrinsic2 {
                    start: Frame2 {
                        origin: Point2::new(0.0, -1.0),
                        x: Vec2::X,
                        y: Vec2::Y,
                    },
                    curvature: CurvatureLaw::circular(1.0),
                    length: TAU,
                })),
            );
            vec![composite(b, &[(spiral, true)])]
        }),
        "2D curve family",
    );
}

#[test]
fn trims_that_do_not_resolve_are_refused() {
    // A point selector off the circle is not projected onto it.
    let error = refusal(|b| {
        let circle = push(b, circle2((0.0, 0.0), 1.0));
        let point = |x: Scalar, y: Scalar| TrimSelector::Point2(Point2::new(x, y));
        let upper = trimmed(
            b,
            circle,
            point(1.0, 0.0),
            point(-1.0, 0.1),
            true,
            TrimmingPreference::Cartesian,
        );
        let lower = by_parameter(b, circle, PI, TAU, true);
        vec![composite(b, &[(upper, true), (lower, true)])]
    });
    assert!(matches!(error, GeomError::InvalidInput(_)), "{error:?}");
    // An arc length has no single curve to measure along on a relation.
    let error = refusal(|b| {
        let first = push(b, polyline2(&[(0.0, 0.0), (4.0, 0.0)], false));
        let rest = push(
            b,
            polyline2(&[(4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)], false),
        );
        let square = composite(b, &[(first, true), (rest, true)]);
        vec![trimmed(
            b,
            square,
            TrimSelector::ArcLength(0.0),
            TrimSelector::ArcLength(16.0),
            true,
            TrimmingPreference::Parameter,
        )]
    });
    assert!(
        matches!(&error, GeomError::InvalidInput(m) if m.contains("arc-length")),
        "{error:?}"
    );
    // A 3D segment off the parameter plane.
    let error = refusal(|b| {
        let ring = push(
            b,
            polyline3(
                &[
                    (0.0, 0.0, 0.0),
                    (4.0, 0.0, 0.5),
                    (4.0, 3.0, 0.0),
                    (0.0, 3.0, 0.0),
                ],
                true,
            ),
        );
        vec![composite(b, &[(ring, true)])]
    });
    assert!(
        matches!(&error, GeomError::InvalidInput(m) if m.contains("parameter plane")),
        "{error:?}"
    );
}
