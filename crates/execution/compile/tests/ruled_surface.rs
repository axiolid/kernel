//! A ruled surface compiles through `compile_mesh` as an open sheet (#402).
//!
//! `SurfaceRelation::Ruled` joins two curve sides, `S = (1 - v) A + v B`;
//! the mesh is a surface (`MeshClosure::Surface`) whose reported
//! bound is certified. The oracle here evaluates every side in closed form
//! (`cos` and `sin` for a conic, de Boor for a B-spline), never through
//! `axiolid_reference`, and the bound is checked against it: sampled points
//! of the exact surface must lie within it of the mesh, and, where the
//! construction's error is the whole of it, the densest sample comes near.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4, PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Interval, Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_curve::{BSplineCurve3, Circle3, Curve3, Ellipse3, KnotSpec, Line3, Polyline};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, DeviationReport, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_mesh_compile_contract::{CompileOutcome, MeshClosure, MeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, RuledSide, SurfaceRelation,
};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

const CHORD: Scalar = 1e-3;

fn options_with(chord: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(chord)
        .expect("positive chord budget")
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

/// A frame that is neither axis aligned nor at the origin.
fn frame() -> Frame3 {
    let z = Vec3::new(1.0, -2.0, 2.0).normalize();
    let x = Vec3::new(2.0, 1.0, 0.0).normalize();
    Frame3 {
        origin: Point3::new(5.0, -3.0, 1.5),
        x,
        y: z.cross(x),
        z,
    }
}

/// The same frame moved `height` along its axis.
fn raised(height: Scalar) -> Frame3 {
    let base = frame();
    Frame3 {
        origin: base.origin + base.z * height,
        ..base
    }
}

/// A side with its closed-form evaluation.
#[derive(Clone)]
#[allow(clippy::large_enum_variant)]
enum Spec {
    Curve {
        curve: Curve3,
        oracle: Oracle,
        interval: (Scalar, Scalar),
    },
    Apex(Point3),
}

#[derive(Clone)]
enum Oracle {
    Line(Point3, Vec3),
    Conic(Frame3, Scalar, Scalar),
    Spline {
        degree: usize,
        points: Vec<Point3>,
        weights: Option<Vec<Scalar>>,
        knots: Vec<Scalar>,
    },
}

impl Oracle {
    fn at(&self, t: Scalar) -> Point3 {
        match self {
            Self::Line(origin, direction) => *origin + *direction * t,
            Self::Conic(frame, rx, ry) => {
                frame.origin + frame.x * (rx * t.cos()) + frame.y * (ry * t.sin())
            }
            Self::Spline {
                degree,
                points,
                weights,
                knots,
            } => deboor(*degree, points, weights.as_deref(), knots, t),
        }
    }
}

/// Cox-de Boor on an expanded knot vector, in homogeneous coordinates when
/// the curve is rational.
fn deboor(
    p: usize,
    points: &[Point3],
    weights: Option<&[Scalar]>,
    kv: &[Scalar],
    t: Scalar,
) -> Point3 {
    let n = points.len();
    let mut k = p;
    while k < n - 1 && t >= kv[k + 1] {
        k += 1;
    }
    let homogeneous = |i: usize| -> [Scalar; 4] {
        let w = weights.map_or(1.0, |w| w[i]);
        [points[i].x * w, points[i].y * w, points[i].z * w, w]
    };
    let mut d: Vec<[Scalar; 4]> = (0..=p).map(|j| homogeneous(j + k - p)).collect();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = j + k - p;
            let denom = kv[i + p + 1 - r] - kv[i];
            let alpha = if denom == 0.0 {
                0.0
            } else {
                (t - kv[i]) / denom
            };
            let lower = d[j - 1];
            for (value, below) in d[j].iter_mut().zip(lower) {
                *value = below * (1.0 - alpha) + *value * alpha;
            }
        }
    }
    let [x, y, z, w] = d[p];
    Point3::new(x / w, y / w, z / w)
}

fn line(origin: Point3, direction: Vec3, interval: (Scalar, Scalar)) -> Spec {
    Spec::Curve {
        curve: Curve3::Line(Line3 { origin, direction }),
        oracle: Oracle::Line(origin, direction),
        interval,
    }
}

fn arc(frame: Frame3, radius: Scalar, interval: (Scalar, Scalar)) -> Spec {
    Spec::Curve {
        curve: Curve3::Circle(Circle3 { frame, radius }),
        oracle: Oracle::Conic(frame, radius, radius),
        interval,
    }
}

fn ellipse(frame: Frame3, rx: Scalar, ry: Scalar, interval: (Scalar, Scalar)) -> Spec {
    Spec::Curve {
        curve: Curve3::Ellipse(Ellipse3 {
            frame,
            semi_axis_x: rx,
            semi_axis_y: ry,
        }),
        oracle: Oracle::Conic(frame, rx, ry),
        interval,
    }
}

fn spline(
    degree: usize,
    points: &[Point3],
    knots: &[Scalar],
    multiplicities: &[u32],
    interval: (Scalar, Scalar),
) -> Spec {
    rational_spline(degree, points, None, knots, multiplicities, interval)
}

fn rational_spline(
    degree: usize,
    points: &[Point3],
    weights: Option<&[Scalar]>,
    knots: &[Scalar],
    multiplicities: &[u32],
    interval: (Scalar, Scalar),
) -> Spec {
    let expanded: Vec<Scalar> = knots
        .iter()
        .zip(multiplicities)
        .flat_map(|(knot, m)| std::iter::repeat_n(*knot, *m as usize))
        .collect();
    Spec::Curve {
        curve: Curve3::BSpline(BSplineCurve3 {
            degree: degree as u16,
            control_points: points.to_vec(),
            knots: knots.to_vec(),
            multiplicities: multiplicities.to_vec(),
            weights: weights.map(<[Scalar]>::to_vec),
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }),
        oracle: Oracle::Spline {
            degree,
            points: points.to_vec(),
            weights: weights.map(<[Scalar]>::to_vec),
            knots: expanded,
        },
        interval,
    }
}

impl Spec {
    fn at(&self, u: Scalar) -> Point3 {
        match self {
            Self::Curve {
                oracle, interval, ..
            } => oracle.at(interval.0 + (interval.1 - interval.0) * u),
            Self::Apex(point) => *point,
        }
    }
}

/// The exact surface.
fn surface(a: &Spec, b: &Spec, u: Scalar, v: Scalar) -> Point3 {
    a.at(u) * (1.0 - v) + b.at(u) * v
}

struct Sheet {
    graph: GeometryGraph,
    root: NodeId,
}

fn side(builder: &mut GeometryGraphBuilder, spec: &Spec) -> RuledSide {
    match spec {
        Spec::Curve {
            curve, interval, ..
        } => RuledSide::Curve {
            curve: builder.push(GeometryNode::Curve3(curve.clone())).unwrap(),
            interval: Interval {
                start: interval.0,
                end: interval.1,
            },
        },
        Spec::Apex(point) => RuledSide::Apex(*point),
    }
}

fn sheet(a: &Spec, b: &Spec) -> Sheet {
    let mut builder = GeometryGraphBuilder::new();
    let first = side(&mut builder, a);
    let second = side(&mut builder, b);
    let root = builder
        .push(GeometryNode::SurfaceRelation(SurfaceRelation::Ruled {
            first,
            second,
        }))
        .unwrap();
    Sheet {
        graph: builder.finish(vec![root]).unwrap(),
        root,
    }
}

fn compile_at(s: &Sheet, chord: Scalar) -> GeomResult<(CompileOutcome, DeviationReport)> {
    compiler().compile_mesh_with_deviation(&s.graph, s.root, &options_with(chord))
}

fn compile(s: &Sheet) -> GeomResult<(CompileOutcome, DeviationReport)> {
    compile_at(s, CHORD)
}

fn distance_to_mesh(mesh: &TriMesh, p: Point3) -> Scalar {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let triangle = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            match closest_point_on_triangle(p, triangle) {
                Ok(q) => (q - p).length(),
                // A triangle of no area is the segments of its edges.
                Err(_) => (0..3)
                    .map(|k| distance_to_segment(p, triangle[k], triangle[(k + 1) % 3]))
                    .fold(Scalar::INFINITY, Scalar::min),
            }
        })
        .fold(Scalar::INFINITY, Scalar::min)
}

fn distance_to_segment(p: Point3, a: Point3, b: Point3) -> Scalar {
    let edge = b - a;
    let length_squared = edge.dot(edge);
    let along = if length_squared > 0.0 {
        ((p - a).dot(edge) / length_squared).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (a + edge * along - p).length()
}

/// The sheet's closure is a surface, the report certified and within the
/// budget, and no sampled exact point lies farther from the mesh than the
/// bound; the largest sampled distance and the bound.
fn within_bound(
    outcome: &CompileOutcome,
    report: &DeviationReport,
    a: &Spec,
    b: &Spec,
) -> (Scalar, Scalar) {
    assert_eq!(outcome.closure, MeshClosure::Surface);
    let bound = report.bound.expect("a certified bound");
    assert!(report.meets_requested(), "{bound} > {CHORD}");
    assert!(report
        .contributions
        .iter()
        .any(|c| c.path == DeviationPath::RuledSurface
            && c.detail == "ruled"
            && matches!(c.bound, DeviationBound::Certified(_))));
    let (nu, nv) = (97, 23);
    let mut worst: Scalar = 0.0;
    for i in 0..=nu {
        for j in 0..=nv {
            let u = (i as Scalar + 0.5) / (nu as Scalar + 1.0);
            let v = (j as Scalar + 0.31) / (nv as Scalar + 1.0);
            let p = surface(a, b, u, v);
            let d = distance_to_mesh(&outcome.mesh, p);
            assert!(
                d <= bound + 1e-12 * p.abs().max_element().max(1.0),
                "({u}, {v}) at {p:?} is {d} from the mesh, beyond the bound {bound}"
            );
            worst = worst.max(d);
        }
    }
    (worst, bound)
}

/// The area-weighted sum of the mesh's triangle normals.
fn mean_normal(mesh: &TriMesh) -> Vec3 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [p, q, r] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            (q - p).cross(r - p)
        })
        .fold(Vec3::ZERO, |sum, n| sum + n)
}

/// `S_u x S_v` of the exact surface, by central differences.
fn exact_normal(a: &Spec, b: &Spec, u: Scalar, v: Scalar) -> Vec3 {
    let h = 1e-6;
    let su = surface(a, b, u + h, v) - surface(a, b, u - h, v);
    let sv = surface(a, b, u, v + h) - surface(a, b, u, v - h);
    su.cross(sv)
}

fn distinct_vertices(mesh: &TriMesh, at: Point3) -> usize {
    mesh.positions
        .iter()
        .filter(|p| (**p - at).length() < 1e-12)
        .count()
}

// --- the sheets -----------------------------------------------------------------

/// Coaxial equal arcs over an interval of length other than one: the rulings
/// are parallel and equal, there is no twist, and every vertex lies on the
/// cylinder; the chords are the whole error, the sagitta the bound.
#[test]
fn a_cylinder_patch_from_two_coaxial_arcs_meets_its_sagitta() {
    let interval = (0.3, 0.3 + FRAC_PI_2);
    let (a, b) = (arc(frame(), 2.0, interval), arc(raised(1.5), 2.0, interval));
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    let f = frame();
    for p in &outcome.mesh.positions {
        let d = *p - f.origin;
        let height = d.dot(f.z);
        let radial = (d - f.z * height).length();
        assert!((radial - 2.0).abs() < 1e-12, "{radial}");
        assert!((-1e-12..=1.5 + 1e-12).contains(&height), "{height}");
    }
    let (worst, bound) = within_bound(&outcome, &report, &a, &b);
    assert!(worst >= 0.9 * bound, "{worst} vs {bound}");
    assert!(mean_normal(&outcome.mesh).dot(exact_normal(&a, &b, 0.5, 0.5)) > 0.0);
}

/// A hyperbolic paraboloid `z = x y` from two skew lines: the rulings twist,
/// and the twist term is the bound.
#[test]
fn a_hyperbolic_paraboloid_from_two_skew_lines_meets_its_twist_bound() {
    let a = line(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        (0.0, 1.0),
    );
    let b = line(
        Point3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 1.0),
        (0.0, 1.0),
    );
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    for p in &outcome.mesh.positions {
        assert!((p.z - p.x * p.y).abs() < 1e-12, "{p:?}");
    }
    let (worst, bound) = within_bound(&outcome, &report, &a, &b);
    assert!(worst >= 0.8 * bound, "{worst} vs {bound}");
    assert!(mean_normal(&outcome.mesh).dot(exact_normal(&a, &b, 0.5, 0.5)) > 0.0);
}

#[test]
fn an_arc_and_a_line_meet_the_bound() {
    let a = arc(frame(), 2.0, (0.2, 1.4));
    let b = line(
        Point3::new(-1.0, 4.0, 3.0),
        Vec3::new(1.0, 0.3, -0.2),
        (0.0, 2.0),
    );
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
}

#[test]
fn arcs_of_unequal_radius_and_sweep_meet_the_bound() {
    let a = arc(frame(), 2.0, (0.1, 0.1 + FRAC_PI_2));
    let b = arc(raised(1.0), 3.5, (1.0, 1.0 + FRAC_PI_3));
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
}

#[test]
fn ellipses_over_different_intervals_meet_the_bound() {
    let a = ellipse(frame(), 3.0, 1.5, (0.0, 2.0));
    let b = ellipse(raised(2.0), 2.0, 1.0, (4.0, 3.0));
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
}

/// A two-span C1 cubic over knots 0, 1, 2 against a one-span cubic over 0 to
/// 1: the sides run over different knots.
#[test]
fn splines_over_different_knots_meet_the_bound() {
    let a = spline(
        3,
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.4, 0.9, 0.0),
            Point3::new(1.0, 1.2, 0.3),
            Point3::new(1.6, 0.8, 0.5),
            Point3::new(2.0, 0.2, 0.4),
            Point3::new(2.5, 0.0, 0.0),
        ],
        &[0.0, 1.0, 2.0],
        &[4, 2, 4],
        (0.0, 2.0),
    );
    let b = spline(
        3,
        &[
            Point3::new(0.0, 2.0, 1.0),
            Point3::new(1.0, 3.0, 1.5),
            Point3::new(2.0, 3.0, 0.5),
            Point3::new(3.0, 2.0, 0.0),
        ],
        &[0.0, 1.0],
        &[4, 4],
        (0.0, 1.0),
    );
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
}

/// A full-multiplicity knot at 0.3 is a kink: the grid has a column there,
/// so the mesh passes through the corner's control point, whichever way the
/// interval runs.
#[test]
fn a_kinked_spline_side_has_a_column_at_its_kink() {
    for interval in [(0.0, 1.0), (1.0, 0.0)] {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.5),
            Point3::new(2.0, 1.0, 0.5),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(3.0, 2.0, -0.5),
            Point3::new(4.0, 3.0, 0.0),
            Point3::new(5.0, 3.0, 0.5),
        ];
        let a = spline(3, &points, &[0.0, 0.3, 1.0], &[4, 3, 4], interval);
        let b = line(
            Point3::new(0.0, -1.0, 2.0),
            Vec3::new(5.0, 0.0, 0.0),
            (0.0, 1.0),
        );
        let s = sheet(&a, &b);
        let (outcome, report) = compile(&s).unwrap();
        within_bound(&outcome, &report, &a, &b);
        assert!(
            distinct_vertices(&outcome.mesh, points[3]) >= 1,
            "{interval:?}: no vertex at the kink"
        );
    }
}

/// A rational quadratic is a quarter circle exactly (weights 1, cos 45
/// degrees, 1): the bound holds over its homogeneous form.
#[test]
fn a_rational_side_meets_the_bound() {
    let f = frame();
    let r = 2.0;
    let points = [
        f.origin + f.x * r,
        f.origin + (f.x + f.y) * r,
        f.origin + f.y * r,
    ];
    let weights = [1.0, FRAC_PI_4.cos(), 1.0];
    let a = rational_spline(2, &points, Some(&weights), &[0.0, 1.0], &[3, 3], (0.0, 1.0));
    let b = line(
        Point3::new(4.0, -2.0, 4.0),
        Vec3::new(1.0, 0.5, 0.2),
        (0.0, 1.5),
    );
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
    // The oracle's arc is the circle: its middle lies on it.
    let d = a.at(0.5) - f.origin;
    assert!(((d - f.z * d.dot(f.z)).length() - r).abs() < 1e-12);
}

/// The interval is oriented: a side read from its end to its start pairs
/// its far end with the other side's start.
#[test]
fn a_reversed_interval_pairs_the_far_end() {
    // Not opposite directions: the rulings would all cross at one row.
    let a = line(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 1.0, 0.0),
        (0.0, 1.0),
    );
    let b = line(
        Point3::new(0.0, 1.0, 1.0),
        Vec3::new(1.0, 0.0, 0.0),
        (1.0, 0.0),
    );
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
    // The first column's far vertex is B at the interval's start, t = 1.
    assert_eq!(
        distinct_vertices(&outcome.mesh, Point3::new(1.0, 1.0, 1.0)),
        1
    );
    assert!(outcome.mesh.positions.contains(&Point3::new(0.0, 0.0, 0.0)));
}

/// An arc ruled to a point on its axis is a cone: one apex vertex and a
/// triangle per column.
#[test]
fn an_arc_to_an_apex_on_its_axis_is_a_cone() {
    let f = frame();
    let apex = f.origin + f.z * 3.0;
    let a = arc(f, 2.0, (0.0, PI));
    let b = Spec::Apex(apex);
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    let (worst, bound) = within_bound(&outcome, &report, &a, &b);
    assert!(worst >= 0.5 * bound, "{worst} vs {bound}");
    let mesh = &outcome.mesh;
    assert_eq!(distinct_vertices(mesh, apex), 1);
    assert_eq!(mesh.positions.len(), mesh.indices.len() / 3 + 1 + 1);
    let along: Vec<_> = mesh.positions.iter().filter(|p| **p != apex).collect();
    assert_eq!(along.len(), mesh.indices.len() / 3 + 1);
    assert!(mean_normal(mesh).dot(exact_normal(&a, &b, 0.5, 0.5)) > 0.0);

    // The apex first turns the sheet round: S_u x S_v of the swapped sides.
    let swapped = sheet(&b, &a);
    let (flipped, _) = compile(&swapped).unwrap();
    assert_eq!(distinct_vertices(&flipped.mesh, apex), 1);
    assert!(mean_normal(&flipped.mesh).dot(exact_normal(&b, &a, 0.5, 0.5)) > 0.0);
    assert!(mean_normal(&flipped.mesh).dot(mean_normal(mesh)) < 0.0);
}

/// An instance places the sheet and keeps its closure.
#[test]
fn an_instance_of_the_sheet_is_a_moved_sheet() {
    let a = arc(frame(), 2.0, (0.3, 1.5));
    let b = arc(raised(1.5), 2.0, (0.3, 1.5));
    let plain = sheet(&a, &b);
    let (base, _) = compile(&plain).unwrap();

    let mut builder = GeometryGraphBuilder::new();
    let first = side(&mut builder, &a);
    let second = side(&mut builder, &b);
    let ruled = builder
        .push(GeometryNode::SurfaceRelation(SurfaceRelation::Ruled {
            first,
            second,
        }))
        .unwrap();
    let shift = Vec3::new(10.0, 0.0, -4.0);
    let root = builder
        .push(GeometryNode::Instance(Instance {
            source: ruled,
            transform: Transform3::from_translation(shift),
        }))
        .unwrap();
    let moved = Sheet {
        graph: builder.finish(vec![root]).unwrap(),
        root,
    };
    let (outcome, report) = compile(&moved).unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!(report.meets_requested());
    assert_eq!(outcome.mesh.positions.len(), base.mesh.positions.len());
    for (moved, plain) in outcome.mesh.positions.iter().zip(&base.mesh.positions) {
        assert!((*moved - (*plain + shift)).length() < 1e-9);
    }
}

// --- refusals -------------------------------------------------------------------

fn refusal(a: &Spec, b: &Spec) -> GeomError {
    compile(&sheet(a, b)).unwrap_err()
}

#[test]
fn a_side_on_another_curve_family_is_refused_by_name() {
    let polyline = Spec::Curve {
        curve: Curve3::Polyline(Polyline {
            points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            closed: false,
        }),
        oracle: Oracle::Line(Point3::ZERO, Vec3::X),
        interval: (0.0, 2.0),
    };
    let other = line(Point3::new(0.0, 0.0, 1.0), Vec3::X, (0.0, 1.0));
    for error in [refusal(&polyline, &other), refusal(&other, &polyline)] {
        assert!(
            matches!(
                error,
                GeomError::UnsupportedInput {
                    operation: Operation::SurfaceEvaluation,
                    input: "a ruled surface side on a curve other than a line, circle, ellipse \
                            or B-spline",
                    ..
                }
            ),
            "{error:?}"
        );
    }
}

#[test]
fn a_conic_over_more_than_a_turn_is_refused() {
    let a = arc(frame(), 2.0, (0.0, TAU + 0.1));
    let b = arc(raised(1.0), 2.0, (0.0, 1.0));
    for error in [refusal(&a, &b), refusal(&b, &a)] {
        assert!(
            matches!(
                &error,
                GeomError::InvalidInput(message)
                    if message == "ruled surface side runs over more than one turn of its conic"
            ),
            "{error:?}"
        );
    }
    // A whole turn, either way round, is allowed.
    for interval in [(0.0, TAU), (TAU, 0.0)] {
        let a = arc(frame(), 2.0, interval);
        compile(&sheet(&a, &b)).unwrap();
    }
}

#[test]
fn a_spline_interval_outside_its_domain_is_refused() {
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
    ];
    let b = line(Point3::new(0.0, 0.0, 1.0), Vec3::X, (0.0, 1.0));
    for interval in [(0.0, 1.1), (-0.1, 1.0), (1.0, 2.0)] {
        let a = spline(3, &points, &[0.0, 1.0], &[4, 4], interval);
        let error = refusal(&a, &b);
        assert!(
            matches!(
                &error,
                GeomError::InvalidInput(message)
                    if message == "ruled surface side interval leaves its curve's domain"
            ),
            "{interval:?}: {error:?}"
        );
    }
}

/// Two lines over the same points under different parameterisations are
/// distinct nodes, but the sheet between them has no area.
#[test]
fn sides_that_coincide_give_no_sheet() {
    let a = line(Point3::ZERO, Vec3::X, (0.0, 1.0));
    let b = line(Point3::ZERO, Vec3::new(2.0, 0.0, 0.0), (0.0, 0.5));
    let error = refusal(&a, &b);
    assert!(
        matches!(
            &error,
            GeomError::Degenerate(message)
                if message == "the ruled surface has no area: its sides coincide"
        ),
        "{error:?}"
    );
}

#[test]
fn a_sheet_beyond_the_cell_budget_is_refused() {
    let a = arc(frame(), 1000.0, (0.0, TAU));
    let b = arc(raised(1.0), 1000.0, (0.0, TAU));
    let error = compile_at(&sheet(&a, &b), 1e-9).unwrap_err();
    assert!(
        matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "ruled surface cells"
            }
        ),
        "{error:?}"
    );
}

// --- what a surface relation is not ---------------------------------------------

/// The exact compiler has no surface relation to build: it names the family.
#[test]
fn the_exact_compiler_refuses_a_ruled_surface() {
    let a = line(Point3::ZERO, Vec3::X, (0.0, 1.0));
    let b = line(Point3::new(0.0, 1.0, 1.0), Vec3::X, (0.0, 1.0));
    let s = sheet(&a, &b);
    let error = ReferenceExactCompiler::new()
        .compile_exact(&s.graph, s.root, &options_with(CHORD))
        .unwrap_err();
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::GraphCompilation,
                input: "surface relation",
                ..
            }
        ),
        "{error:?}"
    );
}

/// A B-rep face on a surface relation has no `Surface` to be tessellated
/// over: refused as unsupported input, not as a malformed graph.
#[test]
fn a_brep_face_on_a_ruled_surface_is_refused_by_name() {
    let a = line(Point3::ZERO, Vec3::X, (0.0, 1.0));
    let b = line(Point3::new(0.0, 1.0, 1.0), Vec3::X, (0.0, 1.0));
    let mut builder = GeometryGraphBuilder::new();
    let first = side(&mut builder, &a);
    let second = side(&mut builder, &b);
    let ruled = builder
        .push(GeometryNode::SurfaceRelation(SurfaceRelation::Ruled {
            first,
            second,
        }))
        .unwrap();

    let mut brep = BRep::default();
    let corners = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(0.0, 1.0, 1.0),
    ];
    let vertices: Vec<_> = corners
        .iter()
        .map(|&position| brep.add_vertex(Vertex { position }))
        .collect();
    let edges: Vec<_> = (0..4)
        .map(|k| {
            brep.add_edge(Edge {
                start: vertices[k],
                end: vertices[(k + 1) % 4],
                curve: None,
            })
        })
        .collect();
    let wire = brep.add_loop(Loop {
        edges: edges
            .into_iter()
            .map(|edge| EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: None,
            })
            .collect(),
    });
    let face = brep.add_face(Face {
        surface: Some(ruled),
        bounds: vec![FaceBound {
            loop_id: wire,
            orientation: Orientation::Forward,
            outer: true,
        }],
        orientation: Orientation::Forward,
    });
    let shell = brep.add_shell(Shell {
        faces: vec![(face, Orientation::Forward)],
        closed: false,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    let root = builder.push(GeometryNode::BRep(brep)).unwrap();
    let graph = builder.finish(vec![root]).unwrap();

    let error = compiler()
        .compile_mesh(&graph, root, &options_with(CHORD))
        .unwrap_err();
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::Tessellation,
                input: "a B-rep face on a surface relation",
                ..
            }
        ),
        "{error:?}"
    );
}

// --- sides that meet at both ends of a cell -----------------------------------

/// A line against a quadratic that meets it at the middle and the start:
/// the first cell's rulings both vanish, yet the sides part between them.
#[test]
fn sides_that_meet_at_both_ends_of_a_cell_still_mesh_the_lens() {
    let a = line(Point3::ZERO, Vec3::new(2.0, 0.0, 0.0), (0.0, 1.0));
    let b = spline(
        2,
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 5e-4, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
            Point3::new(2.0, 0.0, 1e-3),
        ],
        &[0.0, 0.5, 1.0],
        &[3, 2, 3],
        (0.0, 1.0),
    );
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
    // The lens is meshed, not left as one triangle past the meeting point.
    assert!(outcome.mesh.indices.len() / 3 > 2);
}

/// The same line against a quadratic Bezier that meets it only at the ends:
/// not equal, so meshed rather than refused.
#[test]
fn sides_that_meet_only_at_the_ends_are_not_coincident() {
    let a = line(Point3::ZERO, Vec3::new(2.0, 0.0, 0.0), (0.0, 1.0));
    let b = spline(
        2,
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 5e-4, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        &[0.0, 1.0],
        &[3, 3],
        (0.0, 1.0),
    );
    let s = sheet(&a, &b);
    let (outcome, report) = compile(&s).unwrap();
    within_bound(&outcome, &report, &a, &b);
    assert!(outcome.mesh.indices.len() / 3 >= 2);
}

// --- the twist split and the caps -----------------------------------------------

/// Skew lines whose rulings twist by 10 per unit: the columns are split by
/// the twist estimate, about 71 by 71 cells; without the split the rows
/// alone would pass their cap.
#[test]
fn a_strongly_twisted_strip_is_split_into_columns_as_well_as_rows() {
    let a = line(Point3::ZERO, Vec3::new(1.0, 0.0, 0.0), (0.0, 1.0));
    let b = line(
        Point3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 10.0),
        (0.0, 1.0),
    );
    let (outcome, report) = compile(&sheet(&a, &b)).unwrap();
    assert!(report.meets_requested());
    let triangles = outcome.mesh.indices.len() / 3;
    assert!((9_000..=11_000).contains(&triangles), "{triangles}");
}

fn over_budget(error: GeomError) {
    assert!(
        matches!(
            error,
            GeomError::BudgetExceeded {
                resource: "ruled surface cells"
            }
        ),
        "{error:?}"
    );
}

/// Columns past the per-direction cap: two kinked pieces, each within it,
/// the second side the first moved, so there is no twist to split.
#[test]
fn columns_past_the_cap_are_refused() {
    let scale = 1e4;
    let points = |z: Scalar| {
        [
            Point3::new(0.0, 0.0, z),
            Point3::new(1.0 * scale, 1.0 * scale, z),
            Point3::new(2.0 * scale, 0.0, z),
            Point3::new(3.0 * scale, 1.0 * scale, z),
            Point3::new(4.0 * scale, 0.0, z),
        ]
    };
    let a = spline(2, &points(0.0), &[0.0, 1.0, 2.0], &[3, 2, 3], (0.0, 2.0));
    let b = spline(2, &points(1000.0), &[0.0, 1.0, 2.0], &[3, 2, 3], (0.0, 2.0));
    over_budget(compile(&sheet(&a, &b)).unwrap_err());
}

/// Columns the twist split would take past the cap: 36 arc cells, each cut
/// 167 ways.
#[test]
fn a_split_past_the_cap_is_refused() {
    let a = arc(frame(), 2.0, (0.0, FRAC_PI_2));
    let b = line(
        Point3::new(9.0, 9.0, 9.0),
        Vec3::new(0.0, 0.0, 2000.0),
        (0.0, 1.0),
    );
    over_budget(compile(&sheet(&a, &b)).unwrap_err());
}

/// Columns and rows each within their cap, their product past the total.
#[test]
fn a_sheet_past_the_total_cap_is_refused() {
    let a = line(Point3::ZERO, Vec3::new(1.0, 0.0, 0.0), (0.0, 1.0));
    let b = line(
        Point3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 2100.0),
        (0.0, 1.0),
    );
    over_budget(compile(&sheet(&a, &b)).unwrap_err());
}

// --- what a surface relation is not, as a basis or a reference ------------------

fn ruled_pair(builder: &mut GeometryGraphBuilder) -> NodeId {
    let a = line(Point3::ZERO, Vec3::X, (0.0, 1.0));
    let b = line(Point3::new(0.0, 1.0, 1.0), Vec3::X, (0.0, 1.0));
    let first = side(builder, &a);
    let second = side(builder, &b);
    builder
        .push(GeometryNode::SurfaceRelation(SurfaceRelation::Ruled {
            first,
            second,
        }))
        .unwrap()
}

/// A curve-bounded surface over a ruled one is a refusal by name, not a
/// malformed graph.
#[test]
fn a_curve_bounded_surface_over_a_surface_relation_is_refused_by_name() {
    let mut builder = GeometryGraphBuilder::new();
    let basis = ruled_pair(&mut builder);
    let boundary = builder
        .push(GeometryNode::Curve2(axiolid_curve::Curve2::Polyline(
            axiolid_curve::Polyline2 {
                points: vec![
                    axiolid_core::Point2::new(0.0, 0.0),
                    axiolid_core::Point2::new(1.0, 0.0),
                    axiolid_core::Point2::new(1.0, 1.0),
                ],
                closed: true,
            },
        )))
        .unwrap();
    let root = builder
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::CurveBounded {
                basis,
                boundaries: vec![boundary],
                implicit_outer: false,
            },
        ))
        .unwrap();
    let graph = builder.finish(vec![root]).unwrap();
    let error = compiler()
        .compile_mesh(&graph, root, &options_with(CHORD))
        .unwrap_err();
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::SurfaceEvaluation,
                input: "a curve-bounded surface over a surface relation",
                ..
            }
        ),
        "{error:?}"
    );
}

/// A sweep's reference surface that is a ruled one is refused by name.
#[test]
fn a_sweep_over_a_surface_relation_is_refused_by_name() {
    let mut builder = GeometryGraphBuilder::new();
    let reference = ruled_pair(&mut builder);
    let profile = builder
        .push(GeometryNode::Profile(axiolid_profile::Profile::Rectangle(
            axiolid_profile::RectangleProfile {
                x: 0.1,
                y: 0.1,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let directrix = builder
        .push(GeometryNode::Curve3(Curve3::Polyline(Polyline {
            points: vec![Point3::new(0.0, 0.5, 0.5), Point3::new(1.0, 0.5, 0.5)],
            closed: false,
        })))
        .unwrap();
    let root = builder
        .push(GeometryNode::SolidOperation(
            axiolid_model::SolidOperation::SurfaceCurveSweep {
                profile,
                directrix,
                reference_surface: reference,
                parameter_range: None,
            },
        ))
        .unwrap();
    let graph = builder.finish(vec![root]).unwrap();
    let error = compiler()
        .compile_mesh(&graph, root, &options_with(CHORD))
        .unwrap_err();
    assert!(
        matches!(
            error,
            GeomError::UnsupportedInput {
                operation: Operation::SurfaceEvaluation,
                input: "a surface relation as a reference surface",
                ..
            }
        ),
        "{error:?}"
    );
}

// --- no degenerate triangle, no partial coincidence ----------------------------

/// A mesh no consumer refuses: no triangle of three vertices and no area.
fn assert_no_degenerate_triangles(mesh: &TriMesh) {
    for t in mesh.indices.as_chunks::<3>().0 {
        let [p, q, r] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
        assert!(
            (q - p).cross(r - p) != Vec3::ZERO,
            "zero-area triangle {t:?}"
        );
    }
}

/// A line along x against a line along y, both from the origin: the rulings
/// vanish at the start and the sheet is a triangle fan, with no sliver.
#[test]
fn sides_that_share_an_endpoint_have_no_degenerate_triangle() {
    let a = line(Point3::ZERO, Vec3::X, (0.0, 1.0));
    let b = line(Point3::ZERO, Vec3::Y, (0.0, 1.0));
    let (outcome, report) = compile(&sheet(&a, &b)).unwrap();
    within_bound(&outcome, &report, &a, &b);
    assert_no_degenerate_triangles(&outcome.mesh);
    assert_eq!(
        axiolid_mesh::audit_mesh(&outcome.mesh, Tolerance::MILLIMETRE).degenerate_triangles,
        0
    );
    assert_eq!(distinct_vertices(&outcome.mesh, Point3::ZERO), 1);
}

#[test]
fn the_lens_and_the_cone_have_no_degenerate_triangle() {
    let a = line(Point3::ZERO, Vec3::new(2.0, 0.0, 0.0), (0.0, 1.0));
    let b = spline(
        2,
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 5e-4, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
            Point3::new(2.0, 0.0, 1e-3),
        ],
        &[0.0, 0.5, 1.0],
        &[3, 2, 3],
        (0.0, 1.0),
    );
    let (lens, _) = compile(&sheet(&a, &b)).unwrap();
    assert_no_degenerate_triangles(&lens.mesh);
    assert_eq!(
        axiolid_mesh::audit_mesh(&lens.mesh, Tolerance::MILLIMETRE).degenerate_triangles,
        0
    );

    let f = frame();
    let cone = sheet(&arc(f, 2.0, (0.0, PI)), &Spec::Apex(f.origin + f.z * 3.0));
    let (cone, _) = compile(&cone).unwrap();
    assert_no_degenerate_triangles(&cone.mesh);
    assert_eq!(
        axiolid_mesh::audit_mesh(&cone.mesh, Tolerance::MILLIMETRE).degenerate_triangles,
        0
    );
}

/// Rulings that cross at one row collapse the sheet there: refused, never
/// returned with zero-area triangles.
#[test]
fn rulings_that_all_cross_at_a_row_are_a_collapse() {
    let a = line(Point3::ZERO, Vec3::X, (0.0, 1.0));
    let b = line(Point3::new(0.0, 1.0, 1.0), Vec3::X, (1.0, 0.0));
    let error = refusal(&a, &b);
    assert!(
        matches!(
            &error,
            GeomError::Degenerate(message) if message.contains("collapses to zero area")
        ),
        "{error:?}"
    );
}

fn degenerate_over(error: GeomError, expected: &str) {
    assert!(
        matches!(&error, GeomError::Degenerate(message) if message == expected),
        "{error:?}"
    );
}

/// The two quadratics are one curve on [0, 0.5] and part after: refused as
/// coinciding over part of their span.
#[test]
fn sides_that_coincide_over_part_of_their_span_are_refused() {
    let points = |lift: Scalar| {
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.5, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.5, -0.5, 0.0),
            Point3::new(2.0, 0.0, lift),
        ]
    };
    let a = spline(2, &points(0.0), &[0.0, 0.5, 1.0], &[3, 2, 3], (0.0, 1.0));
    let b = spline(2, &points(1.0), &[0.0, 0.5, 1.0], &[3, 2, 3], (0.0, 1.0));
    degenerate_over(
        refusal(&a, &b),
        "the ruled surface's sides coincide over part of their span",
    );
}

/// One bump under two parameterisations: all of it coincides, refused as
/// such and not as over budget.
#[test]
fn one_curve_under_two_knot_vectors_is_refused_as_coincident() {
    let bump = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ];
    let a = spline(2, &bump, &[0.0, 1.0], &[3, 3], (0.0, 1.0));
    let b = spline(2, &bump, &[0.0, 2.0], &[3, 3], (0.0, 2.0));
    degenerate_over(
        refusal(&a, &b),
        "the ruled surface has no area: its sides coincide",
    );
}

/// The same, so curved that the halving reaches the column cap first.
#[test]
fn coincident_sides_too_curved_to_settle_are_over_budget() {
    let scale = 9e3;
    let bump = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(scale, scale, 0.0),
        Point3::new(2.0 * scale, 0.0, 0.0),
    ];
    let a = spline(2, &bump, &[0.0, 1.0], &[3, 3], (0.0, 1.0));
    let b = spline(2, &bump, &[0.0, 2.0], &[3, 3], (0.0, 2.0));
    over_budget(refusal(&a, &b));
}
