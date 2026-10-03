//! The deviation report is a bound: no sample of the exact surface lies
//! further from the mesh than it says (#232).
//!
//! Each case compiles a body with `compile_mesh_with_deviation`, samples
//! the exact surface densely (written out here in closed form or through
//! the reference evaluator, never read back from the compiler), measures
//! each sample's distance to the nearest mesh triangle, and asserts the
//! largest is at or below the reported bound. It also prints the ratio of
//! the bound to the measured maximum and asserts it is not absurd: a bound
//! a hundred times the truth would be sound and useless.
//!
//! Paths covered: profile flattening of ellipses, splines and sinusoids
//! (#232 (e)); curved B-rep faces on cylinders, cones, spheres, tori and
//! B-spline surfaces, with straight and curved trims (#232 (c)); disks
//! swept along B-splines and ellipses (#232 (b)); primitive cylinders and
//! cones; instances and the paths reported unbounded by name.

use std::collections::HashMap;
use std::f64::consts::{PI, TAU};

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{
    Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Transform2, Transform3, Vec2, Vec3,
};
use axiolid_curve::{
    BSplineCurve2, BSplineCurve3, Circle2, CurvatureLaw, Curve2, Curve3, Ellipse3, Intrinsic2,
    KnotSpec, Line2, Sinusoid2,
};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation};
use axiolid_primitive::Primitive;
use axiolid_profile::{Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment};
use axiolid_surface::{BSplineSurface, Cone, Cylinder, Sphere, Surface, Torus};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(budget.min(1e-3), 1e-9).unwrap())
        .with_chord_error(budget)
        .unwrap()
}

fn compile(
    build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId,
    options: &ExecutionOptions,
) -> (TriMesh, DeviationReport) {
    let mut builder = GeometryGraphBuilder::new();
    let root = build(&mut builder);
    let graph = builder.finish(vec![root]).expect("a valid graph");
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, options)
        .expect("the body compiles");
    (outcome.mesh, report)
}

/// Triangles binned on a uniform grid, each in every cell its box touches
/// once grown by `reach`.
struct TriangleGrid {
    triangles: Vec<[Point3; 3]>,
    cells: HashMap<[i64; 3], Vec<usize>>,
    cell: Scalar,
    reach: Scalar,
}

impl TriangleGrid {
    fn new(mesh: &TriMesh, reach: Scalar) -> Self {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .chunks_exact(3)
            .map(|t| [0, 1, 2].map(|i| mesh.positions[t[i] as usize]))
            .collect();
        let (lo, hi) = mesh.positions.iter().fold(
            (
                Vec3::splat(Scalar::INFINITY),
                Vec3::splat(Scalar::NEG_INFINITY),
            ),
            |(lo, hi), p| (lo.min(*p), hi.max(*p)),
        );
        let cell = ((hi - lo).max_element() / 64.0).max(4.0 * reach);
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (i, t) in triangles.iter().enumerate() {
            let a = t[0].min(t[1]).min(t[2]) - Vec3::splat(reach);
            let b = t[0].max(t[1]).max(t[2]) + Vec3::splat(reach);
            let (ka, kb) = (key(a, cell), key(b, cell));
            for x in ka[0]..=kb[0] {
                for y in ka[1]..=kb[1] {
                    for z in ka[2]..=kb[2] {
                        cells.entry([x, y, z]).or_default().push(i);
                    }
                }
            }
        }
        Self {
            triangles,
            cells,
            cell,
            reach,
        }
    }

    fn distance(&self, p: Point3) -> Scalar {
        let mut best = self.reach;
        for &i in self.cells.get(&key(p, self.cell)).into_iter().flatten() {
            if let Ok(q) = closest_point_on_triangle(p, self.triangles[i]) {
                best = best.min(p.distance(q));
            }
        }
        best
    }
}

fn key(p: Point3, cell: Scalar) -> [i64; 3] {
    [p.x, p.y, p.z].map(|v| (v / cell).floor() as i64)
}

/// The report's bound, which must exist.
fn bound_of(report: &DeviationReport) -> Scalar {
    report
        .bound
        .unwrap_or_else(|| panic!("expected a certified bound, got {report:?}"))
}

/// Assert no sample is further from the mesh than the bound, and that the
/// bound is within `ceiling` times the measured maximum. Returns the
/// measured maximum.
fn assert_bounded(
    name: &str,
    mesh: &TriMesh,
    samples: &[Point3],
    bound: Scalar,
    ceiling: Scalar,
) -> Scalar {
    let grid = TriangleGrid::new(mesh, 2.0 * bound);
    let worst = samples
        .iter()
        .map(|p| grid.distance(*p))
        .fold(0.0, Scalar::max);
    let ratio = bound / worst;
    eprintln!(
        "{name}: reported {bound:e}, measured {worst:e}, ratio {ratio:.2}, {} triangles, {} samples",
        mesh.indices.len() / 3,
        samples.len()
    );
    assert!(
        worst <= bound,
        "{name}: a surface point lies {worst:e} from the mesh, above the reported {bound:e}"
    );
    assert!(
        ratio <= ceiling,
        "{name}: bound {bound:e} is {ratio:.1}x the measured {worst:e}"
    );
    worst
}

fn span(a: Scalar, b: Scalar, n: usize) -> impl Iterator<Item = Scalar> + Clone {
    (0..n).map(move |i| a + (b - a) * ((i as Scalar) + 0.5) / n as Scalar)
}

fn frame2() -> Frame2 {
    Frame2 {
        origin: Point2::ZERO,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

fn frame3() -> Frame3 {
    Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

// --- (e) profile flattening ------------------------------------------------

fn extrude(b: &mut GeometryGraphBuilder, profile: Profile, depth: Scalar) -> NodeId {
    let profile = b.push(GeometryNode::Profile(profile)).unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth,
    }))
    .unwrap()
}

/// Wall and cap samples of a prism over a closed boundary `c(t)`,
/// `t in [0, 1]`, with an inside test for the caps.
fn prism_samples(
    boundary: impl Fn(Scalar) -> Point2,
    inside: impl Fn(Point2) -> bool,
    bounds: (Point2, Point2),
    depth: Scalar,
    n: usize,
) -> Vec<Point3> {
    let mut out = Vec::new();
    for t in span(0.0, 1.0, n) {
        let p = boundary(t);
        for z in span(0.0, depth, 3) {
            out.push(Point3::new(p.x, p.y, z));
        }
        for z in [0.0, depth] {
            out.push(Point3::new(p.x, p.y, z));
        }
    }
    let (lo, hi) = bounds;
    for x in span(lo.x, hi.x, 200) {
        for y in span(lo.y, hi.y, 200) {
            if inside(Point2::new(x, y)) {
                out.push(Point3::new(x, y, 0.0));
                out.push(Point3::new(x, y, depth));
            }
        }
    }
    out
}

#[test]
fn an_extruded_ellipse_is_within_its_reported_bound() {
    let (a, b) = (0.8, 0.15);
    for budget in [1e-3, 1e-4] {
        let (mesh, report) = compile(
            |g| {
                extrude(
                    g,
                    Profile::Ellipse(EllipseProfile {
                        semi_axis_x: a,
                        semi_axis_y: b,
                    }),
                    0.2,
                )
            },
            &options(budget),
        );
        assert_eq!(report.contributions[0].path, DeviationPath::Extrusion);
        assert!(report.meets_requested(), "{report:?}");
        let bound = bound_of(&report);
        assert_eq!(bound, budget);
        let samples = prism_samples(
            |t| Point2::new(a * (TAU * t).cos(), b * (TAU * t).sin()),
            |p| (p.x / a).powi(2) + (p.y / b).powi(2) < 1.0,
            (Point2::new(-a, -b), Point2::new(a, b)),
            0.2,
            200_000,
        );
        assert_bounded(
            &format!("ellipse extrusion at {budget:e}"),
            &mesh,
            &samples,
            bound,
            3.0,
        );
    }
}

/// A clamped cubic arch from `(0, 0)` to `(1, 0)`, closed by the x axis.
fn arch() -> (Curve2, Interval) {
    let curve = Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 0.4),
            Point2::new(0.3, 0.7),
            Point2::new(0.7, 0.1),
            Point2::new(1.0, 0.5),
            Point2::new(1.0, 0.0),
        ],
        knots: vec![0.0, 0.3, 0.7, 1.0],
        multiplicities: vec![4, 1, 1, 4],
        weights: Some(vec![1.0, 0.6, 1.4, 1.0, 0.8, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    (curve, Interval::new(0.0, 1.0))
}

fn arch_profile() -> Profile {
    let (curve, domain) = arch();
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve,
                domain,
                same_sense: true,
            },
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: Point2::new(1.0, 0.0),
                    direction: Vec2::new(-1.0, 0.0),
                }),
                domain: Interval::UNIT,
                same_sense: true,
            },
        ]),
        holes: Vec::new(),
    })
}

/// Dense samples of the arch's boundary (spline then base), and its
/// region by a crossing test against a fine polyline.
fn arch_boundary() -> (Vec<Point2>, impl Fn(Point2) -> bool) {
    let (curve, _) = arch();
    let fine: Vec<Point2> = (0..=4000)
        .map(|i| axiolid_reference::curve::evaluate2(&curve, i as Scalar / 4000.0).unwrap())
        .collect();
    let ring = fine.clone();
    let inside = move |p: Point2| {
        let mut crossings = 0;
        let n = ring.len();
        for i in 0..n {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            if (a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y) {
                crossings += 1;
            }
        }
        crossings % 2 == 1
    };
    (fine, inside)
}

#[test]
fn an_extruded_rational_spline_contour_is_within_its_reported_bound() {
    let (curve, _) = arch();
    let (_, inside) = arch_boundary();
    for budget in [1e-3, 1e-4] {
        let (mesh, report) = compile(|g| extrude(g, arch_profile(), 0.1), &options(budget));
        assert!(report.meets_requested(), "{report:?}");
        let bound = bound_of(&report);
        let samples = prism_samples(
            |t| axiolid_reference::curve::evaluate2(&curve, t).unwrap(),
            &inside,
            (Point2::new(0.0, 0.0), Point2::new(1.0, 0.6)),
            0.1,
            200_000,
        );
        assert_bounded(
            // The arch's weights run 0.6 to 1.4, and the rational chord
            // bound is 6x (1 mm) to 10x (0.1 mm) the true deviation: the flattener
            // keeps halving until the bound meets the budget, so the rings
            // come out finer than asked, and the report says the budget.
            &format!("spline contour at {budget:e}"),
            &mesh,
            &samples,
            bound,
            12.0,
        );
    }
}

#[test]
fn a_stretched_derived_profile_reports_the_stretch() {
    let (a, b) = (0.4, 0.1);
    let budget = 1e-3;
    let (mesh, report) = compile(
        |g| {
            extrude(
                g,
                Profile::Derived {
                    basis: Box::new(Profile::Ellipse(EllipseProfile {
                        semi_axis_x: a,
                        semi_axis_y: b,
                    })),
                    transform: Transform2::from_scale(Vec2::new(3.0, 1.0)),
                },
                0.1,
            )
        },
        &options(budget),
    );
    let bound = bound_of(&report);
    // The chords were cut for the basis; the placement triples them along x.
    assert!((bound - 3.0 * budget).abs() <= 1e-9, "{bound}");
    assert!(!report.meets_requested());
    let samples = prism_samples(
        |t| Point2::new(3.0 * a * (TAU * t).cos(), b * (TAU * t).sin()),
        |p| (p.x / (3.0 * a)).powi(2) + (p.y / b).powi(2) < 1.0,
        (Point2::new(-3.0 * a, -b), Point2::new(3.0 * a, b)),
        0.1,
        100_000,
    );
    assert_bounded("derived ellipse x3", &mesh, &samples, bound, 6.0);
}

#[test]
fn a_revolved_spline_profile_is_within_its_reported_bound() {
    let (curve, _) = arch();
    let budget = 1e-3;
    let (mesh, report) = compile(
        |g| {
            let profile = g.push(GeometryNode::Profile(arch_profile())).unwrap();
            g.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
                profile,
                axis_origin: Point3::new(-0.5, 0.0, 0.0),
                axis_direction: Vec3::Y,
                angle: 1.5,
            }))
            .unwrap()
        },
        &options(budget),
    );
    assert_eq!(report.contributions[0].path, DeviationPath::Revolution);
    assert!(report.meets_requested(), "{report:?}");
    let bound = bound_of(&report);
    let rotate = |p: Point2, angle: Scalar| {
        let (s, c) = angle.sin_cos();
        let x = p.x + 0.5;
        Point3::new(-0.5 + x * c, p.y, -x * s)
    };
    let mut samples = Vec::new();
    for t in span(0.0, 1.0, 2000) {
        let p = axiolid_reference::curve::evaluate2(&curve, t).unwrap();
        for angle in span(0.0, 1.5, 300) {
            samples.push(rotate(p, angle));
        }
    }
    let mesh_has_orientation = mesh.indices.len() > 3;
    assert!(mesh_has_orientation);
    assert_bounded("revolved spline profile", &mesh, &samples, bound, 3.0);
}

#[test]
fn a_profile_with_a_sinusoid_segment_is_certified() {
    // A wave top over [0, pi], closed by a line along y = 0 below it.
    let wave = Sinusoid2 {
        mean: 0.5,
        cosine: 0.1,
        sine: 0.3,
    };
    let profile = Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: Point2::new(0.0, 0.0),
                    direction: Vec2::new(PI, 0.0),
                }),
                domain: Interval::UNIT,
                same_sense: true,
            },
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: Point2::new(PI, 0.0),
                    direction: Vec2::new(0.0, wave.height(PI)),
                }),
                domain: Interval::UNIT,
                same_sense: true,
            },
            ProfileSegment {
                curve: Curve2::Sinusoid(wave),
                domain: Interval::new(0.0, PI),
                same_sense: false,
            },
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: Point2::new(0.0, wave.height(0.0)),
                    direction: Vec2::new(0.0, -wave.height(0.0)),
                }),
                domain: Interval::UNIT,
                same_sense: true,
            },
        ]),
        holes: Vec::new(),
    });
    let budget = 1e-3;
    let (mesh, report) = compile(|g| extrude(g, profile, 0.1), &options(budget));
    assert!(report.meets_requested(), "{report:?}");
    let bound = bound_of(&report);
    let mut samples = Vec::new();
    for t in span(0.0, PI, 50_000) {
        for z in span(0.0, 0.1, 3) {
            samples.push(Point3::new(t, wave.height(t), z));
        }
    }
    assert_bounded("sinusoid profile", &mesh, &samples, bound, 3.0);
}

/// An S-shaped cubic crosses its own chord at the midpoint, where the
/// sagitta test looks: measured there, one chord passes for the whole
/// curve while the arc bulges 0.09 either side of it.
#[test]
fn an_s_curve_profile_is_flattened_within_its_bound() {
    let s_curve = Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.33, -0.4),
            Point2::new(0.67, 0.4),
            Point2::new(1.0, 0.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let line = |from: Point2, to: Point2| ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    };
    let profile = Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve: s_curve.clone(),
                domain: Interval::new(0.0, 1.0),
                same_sense: true,
            },
            line(Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)),
            line(Point2::new(1.0, 1.0), Point2::new(0.0, 1.0)),
            line(Point2::new(0.0, 1.0), Point2::new(0.0, 0.0)),
        ]),
        holes: Vec::new(),
    });
    let budget = 1e-3;
    let (mesh, report) = compile(|g| extrude(g, profile, 0.1), &options(budget));
    assert!(report.meets_requested(), "{report:?}");
    let bound = bound_of(&report);
    let mut samples = Vec::new();
    for t in span(0.0, 1.0, 20_000) {
        let p = axiolid_reference::curve::evaluate2(&s_curve, t).unwrap();
        for z in span(0.0, 0.1, 3) {
            samples.push(Point3::new(p.x, p.y, z));
        }
    }
    assert_bounded("S-curve profile", &mesh, &samples, bound, 3.0);
}

#[test]
fn a_profile_with_a_clothoid_segment_is_unbounded_by_name() {
    let clothoid = Intrinsic2 {
        start: frame2(),
        curvature: CurvatureLaw::Polynomial {
            coefficients: vec![0.0, 4.0],
        },
        length: 0.5,
    };
    let end =
        axiolid_reference::curve::evaluate2(&Curve2::Intrinsic(clothoid.clone()), 0.5).unwrap();
    let profile = Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            ProfileSegment {
                curve: Curve2::Intrinsic(clothoid),
                domain: Interval::new(0.0, 0.5),
                same_sense: true,
            },
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: end,
                    direction: -end,
                }),
                domain: Interval::UNIT,
                same_sense: true,
            },
        ]),
        holes: Vec::new(),
    });
    let (_, report) = compile(|g| extrude(g, profile, 0.1), &options(1e-3));
    assert_eq!(report.bound, None);
    assert!(!report.meets_requested());
    let unbounded: Vec<_> = report.unbounded().collect();
    assert_eq!(unbounded.len(), 1);
    assert!(matches!(
        unbounded[0].bound,
        DeviationBound::Unbounded(reason) if reason.contains("clothoid")
    ));
}

// --- (c) curved B-rep faces ------------------------------------------------

/// One curved face whose rings are lists of pcurves, each a `Curve2` over
/// its own domain, consecutive end to start.
fn curved_face(
    builder: &mut GeometryGraphBuilder,
    surface: Surface,
    rings: Vec<(Vec<Curve2>, bool)>,
) -> NodeId {
    let surface_id = builder
        .push(GeometryNode::Surface(surface.clone()))
        .unwrap();
    let mut brep: BRep<NodeId> = BRep::default();
    let mut bounds = Vec::new();
    for (pcurves, outer) in rings {
        let starts: Vec<_> = pcurves
            .iter()
            .map(|c| {
                let domain = axiolid_reference::curve::domain2(c);
                let uv = axiolid_reference::curve::evaluate2(c, domain.start).unwrap();
                brep.add_vertex(Vertex {
                    position: axiolid_reference::surface::evaluate(&surface, uv.x, uv.y).unwrap(),
                })
            })
            .collect();
        let mut uses = Vec::new();
        for (i, c) in pcurves.iter().enumerate() {
            let pcurve = builder.push(GeometryNode::Curve2(c.clone())).unwrap();
            let edge = brep.add_edge(Edge {
                start: starts[i],
                end: starts[(i + 1) % starts.len()],
                curve: None,
            });
            uses.push(EdgeUse {
                edge,
                orientation: Orientation::Forward,
                pcurve: Some(pcurve),
            });
        }
        let loop_id = brep.add_loop(Loop { edges: uses });
        bounds.push(FaceBound {
            loop_id,
            orientation: Orientation::Forward,
            outer,
        });
    }
    let face = brep.add_face(Face {
        surface: Some(surface_id),
        bounds,
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
    builder.push(GeometryNode::BRep(brep)).unwrap()
}

/// The four straight pcurves of the rectangle `[u0, u1] x [v0, v1]`,
/// counter-clockwise.
fn rectangle(u: (Scalar, Scalar), v: (Scalar, Scalar)) -> Vec<Curve2> {
    let corners = [
        Point2::new(u.0, v.0),
        Point2::new(u.1, v.0),
        Point2::new(u.1, v.1),
        Point2::new(u.0, v.1),
    ];
    (0..4)
        .map(|i| {
            Curve2::Line(Line2 {
                origin: corners[i],
                direction: corners[(i + 1) % 4] - corners[i],
            })
        })
        .collect()
}

/// Samples of `surface` over the rectangle, minus an optional parameter
/// disk `(centre, radius)`.
fn patch_samples(
    surface: &Surface,
    u: (Scalar, Scalar),
    v: (Scalar, Scalar),
    hole: Option<(Point2, Scalar)>,
    n: usize,
) -> Vec<Point3> {
    let mut out = Vec::with_capacity(n * n);
    for uu in span(u.0, u.1, n) {
        for vv in span(v.0, v.1, n) {
            if let Some((centre, radius)) = hole {
                if (Point2::new(uu, vv) - centre).length() < radius {
                    continue;
                }
            }
            out.push(axiolid_reference::surface::evaluate(surface, uu, vv).unwrap());
        }
    }
    if let Some((centre, radius)) = hole {
        for t in span(0.0, TAU, 4 * n) {
            let p = centre + Point2::new(t.cos(), t.sin()) * radius;
            out.push(axiolid_reference::surface::evaluate(surface, p.x, p.y).unwrap());
        }
    }
    out
}

fn brep_case(
    name: &str,
    surface: Surface,
    u: (Scalar, Scalar),
    v: (Scalar, Scalar),
    hole: Option<(Point2, Scalar)>,
    ceiling: Scalar,
) {
    // A spline face refines on its derivative-net bounds only up to a
    // vertex cap, beyond which its bound is reported as it stands; at
    // 0.1 mm that cap is reached, so splines are checked at 1 mm.
    let budgets: &[Scalar] = if name.contains("B-spline") {
        &[1e-3]
    } else {
        &[1e-3, 1e-4]
    };
    for &budget in budgets {
        let (mesh, report) = compile(
            |g| {
                let mut rings = vec![(rectangle(u, v), true)];
                if let Some((centre, radius)) = hole {
                    // Clockwise: a hole.
                    rings.push((
                        vec![Curve2::Circle(Circle2 {
                            frame: Frame2 {
                                origin: centre,
                                x: Vec2::X,
                                y: -Vec2::Y,
                            },
                            radius,
                        })],
                        false,
                    ));
                }
                curved_face(g, surface.clone(), rings)
            },
            &options(budget),
        );
        assert_eq!(report.contributions.len(), 1, "{report:?}");
        assert_eq!(report.contributions[0].path, DeviationPath::BRepFace);
        let bound = bound_of(&report);
        let samples = patch_samples(&surface, u, v, hole, 500);
        let measured = assert_bounded(
            &format!("{name} at {budget:e}"),
            &mesh,
            &samples,
            bound,
            ceiling,
        );
        assert!(measured > 0.0);
    }
}

#[test]
fn a_cylinder_face_is_within_its_reported_bound() {
    brep_case(
        "cylinder face",
        Surface::Cylinder(Cylinder {
            frame: frame3(),
            radius: 0.3,
        }),
        (0.0, 2.0),
        (0.0, 0.5),
        None,
        3.0,
    );
}

/// A flat B-spline patch trimmed to a disk: the surface does not bend, so
/// every triangle's own bound is zero, and the exact rim's bulge past its
/// chords (the trim lens) is all the report has to carry.
#[test]
fn a_trim_lens_is_carried_by_the_report() {
    let surface = Surface::BSpline(BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![2, 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    let (centre, radius) = (Point2::new(0.5, 0.5), 0.2);
    let budget = 1e-3;
    let (mesh, report) = compile(
        |g| {
            curved_face(
                g,
                surface.clone(),
                // A round outer trim: the exact disk reaches past its
                // inscribed polygon, by the lens.
                vec![(
                    vec![Curve2::Circle(Circle2 {
                        frame: Frame2 {
                            origin: centre,
                            x: Vec2::X,
                            y: Vec2::Y,
                        },
                        radius,
                    })],
                    true,
                )],
            )
        },
        &options(budget),
    );
    let bound = bound_of(&report);
    let mut samples = Vec::new();
    for t in span(0.0, TAU, 200_000) {
        let p = centre + Point2::new(t.cos(), t.sin()) * radius;
        samples.push(axiolid_reference::surface::evaluate(&surface, p.x, p.y).unwrap());
    }
    assert_bounded("trim lens on a flat patch", &mesh, &samples, bound, 4.0);
}

#[test]
fn a_cylinder_face_with_a_round_hole_is_within_its_reported_bound() {
    brep_case(
        "cylinder face with a hole",
        Surface::Cylinder(Cylinder {
            frame: frame3(),
            radius: 0.3,
        }),
        (0.0, 2.0),
        (0.0, 0.5),
        Some((Point2::new(1.0, 0.25), 0.15)),
        6.0,
    );
}

#[test]
fn a_cone_face_is_within_its_reported_bound() {
    brep_case(
        "cone face",
        Surface::Cone(Cone {
            frame: frame3(),
            radius: 0.4,
            semi_angle: -0.4,
        }),
        (0.0, 2.5),
        (0.0, 0.6),
        None,
        4.0,
    );
}

#[test]
fn a_sphere_face_is_within_its_reported_bound() {
    brep_case(
        "sphere face",
        Surface::Sphere(Sphere {
            frame: frame3(),
            radius: 0.5,
        }),
        (0.0, 1.8),
        (-0.7, 1.0),
        None,
        6.0,
    );
}

#[test]
fn a_torus_face_is_within_its_reported_bound() {
    brep_case(
        "torus face",
        Surface::Torus(Torus {
            frame: frame3(),
            major_radius: 0.5,
            minor_radius: 0.15,
        }),
        (0.0, 1.5),
        (-2.0, 2.5),
        Some((Point2::new(0.7, 0.3), 0.3)),
        6.0,
    );
}

/// A bicubic patch with a gentle bump, trimmed by a round hole.
#[test]
fn a_bspline_face_is_within_its_reported_bound() {
    let net: Vec<Vec<Point3>> = (0..5)
        .map(|i| {
            (0..5)
                .map(|j| {
                    let (x, y) = (i as Scalar * 0.25, j as Scalar * 0.25);
                    let bump = if (1..4).contains(&i) && (1..4).contains(&j) {
                        0.12
                    } else {
                        0.0
                    };
                    Point3::new(x, y, bump + 0.05 * x * y)
                })
                .collect()
        })
        .collect();
    let surface = Surface::BSpline(BSplineSurface {
        u_degree: 3,
        v_degree: 3,
        control_points: net,
        u_knots: vec![0.0, 0.5, 1.0],
        u_multiplicities: vec![4, 1, 4],
        v_knots: vec![0.0, 0.5, 1.0],
        v_multiplicities: vec![4, 1, 4],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    brep_case(
        "B-spline face",
        surface,
        (0.0, 1.0),
        (0.0, 1.0),
        Some((Point2::new(0.4, 0.55), 0.2)),
        5.0,
    );
}

/// A wavy rational patch: sound, but the quotient's net bounds are loose.
#[test]
fn a_rational_bspline_face_is_within_its_reported_bound() {
    let net: Vec<Vec<Point3>> = (0..5)
        .map(|i| {
            (0..4)
                .map(|j| {
                    let (x, y) = (i as Scalar * 0.25, j as Scalar * 0.3);
                    Point3::new(x, y, 0.15 * ((i * 3 + j * 5) % 4) as Scalar)
                })
                .collect()
        })
        .collect();
    let surface = Surface::BSpline(BSplineSurface {
        u_degree: 3,
        v_degree: 2,
        control_points: net,
        u_knots: vec![0.0, 0.5, 1.0],
        u_multiplicities: vec![4, 1, 4],
        v_knots: vec![0.0, 0.5, 1.0],
        v_multiplicities: vec![3, 1, 3],
        weights: Some(
            (0..5)
                .map(|i| {
                    (0..4)
                        .map(|j| 0.8 + 0.15 * ((i + j) % 3) as Scalar)
                        .collect()
                })
                .collect(),
        ),
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    brep_case(
        "rational B-spline face",
        surface,
        (0.0, 1.0),
        (0.0, 1.0),
        Some((Point2::new(0.4, 0.55), 0.2)),
        12.0,
    );
}

// --- (b) disks swept along smooth directrices -------------------------------

fn spine() -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.3, 0.0, 0.1),
            Point3::new(0.5, 0.3, 0.0),
            Point3::new(0.8, 0.4, 0.2),
            Point3::new(1.0, 0.2, 0.3),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// Samples of the tube of `radius` about `curve` over `[s0, s1]`, and its
/// end disks (annuli from `inner`).
fn tube_samples(
    curve: &Curve3,
    s: (Scalar, Scalar),
    radius: Scalar,
    inner: Option<Scalar>,
    n: usize,
) -> Vec<Point3> {
    let frame = |t: Scalar| {
        let c = axiolid_reference::curve::evaluate3(curve, t).unwrap();
        let d = axiolid_reference::curve::derivative3(curve, t)
            .unwrap()
            .normalize();
        let a = if d.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
        let n = (a - d * d.dot(a)).normalize();
        (c, n, d.cross(n))
    };
    let mut out = Vec::new();
    for t in span(s.0, s.1, n) {
        let (c, nn, b) = frame(t);
        for r in std::iter::once(radius).chain(inner) {
            for theta in span(0.0, TAU, 360) {
                out.push(c + (nn * theta.cos() + b * theta.sin()) * r);
            }
        }
    }
    for t in [s.0, s.1] {
        let (c, nn, b) = frame(t);
        for r in span(inner.unwrap_or(0.0), radius, 40) {
            for theta in span(0.0, TAU, 360) {
                out.push(c + (nn * theta.cos() + b * theta.sin()) * r);
            }
        }
    }
    out
}

fn swept(
    b: &mut GeometryGraphBuilder,
    curve: Curve3,
    radius: Scalar,
    inner: Option<Scalar>,
    range: Option<(Scalar, Scalar)>,
) -> NodeId {
    let directrix = b.push(GeometryNode::Curve3(curve)).unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix,
        radius,
        inner_radius: inner,
        parameter_range: range,
        fillet_radius: None,
    }))
    .unwrap()
}

#[test]
fn a_disk_swept_along_a_bspline_is_within_its_reported_bound() {
    let curve = spine();
    for (radius, inner) in [(0.04, None), (0.06, Some(0.045))] {
        let budget = 1e-3;
        let started = std::time::Instant::now();
        let (mesh, report) = compile(
            |g| swept(g, curve.clone(), radius, inner, None),
            &options(budget),
        );
        eprintln!("compiled and certified in {:?}", started.elapsed());
        assert_eq!(report.contributions[0].path, DeviationPath::SweptDisk);
        assert_eq!(report.contributions[0].detail, "B-spline");
        assert!(matches!(
            report.contributions[0].bound,
            DeviationBound::Certified(_)
        ));
        let bound = bound_of(&report);
        let samples = tube_samples(&curve, (0.0, 1.0), radius, inner, 3000);
        assert_bounded(
            &format!("B-spline pipe r {radius} inner {inner:?}"),
            &mesh,
            &samples,
            bound,
            2.0,
        );
    }
}

#[test]
fn a_disk_swept_along_an_ellipse_arc_is_within_its_reported_bound() {
    let curve = Curve3::Ellipse(Ellipse3 {
        frame: frame3(),
        semi_axis_x: 0.6,
        semi_axis_y: 0.3,
    });
    let range = (0.3, 2.6);
    let (mesh, report) = compile(
        |g| swept(g, curve.clone(), 0.05, None, Some(range)),
        &options(1e-3),
    );
    assert_eq!(report.contributions[0].detail, "ellipse");
    let bound = bound_of(&report);
    let samples = tube_samples(&curve, range, 0.05, None, 3000);
    assert_bounded("ellipse-arc pipe", &mesh, &samples, bound, 2.0);
}

// --- primitives, instances, unbounded paths --------------------------------

/// A pipe along a filleted polyline goes through `axiolid_construct::pipe`,
/// which proves the budget; along a straight polyline its exact tube is a
/// capped cylinder, sampled here.
#[test]
fn a_pipe_along_pieces_reports_the_proven_budget() {
    let budget = 1e-3;
    let polyline = |points: Vec<Point3>| {
        Curve3::Polyline(axiolid_curve::Polyline3 {
            points,
            closed: false,
        })
    };
    let (_, report) = compile(
        |g| {
            let directrix = g
                .push(GeometryNode::Curve3(polyline(vec![
                    Point3::ZERO,
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(1.0, 1.0, 0.0),
                ])))
                .unwrap();
            g.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
                directrix,
                radius: 0.05,
                inner_radius: None,
                parameter_range: None,
                fillet_radius: Some(0.2),
            }))
            .unwrap()
        },
        &options(budget),
    );
    assert_eq!(report.contributions[0].detail, "segments and arcs");
    assert_eq!(
        report.contributions[0].bound,
        DeviationBound::Proven(budget)
    );
    assert!(report.meets_requested());

    let r = 0.05;
    let (mesh, report) = compile(
        |g| {
            swept(
                g,
                polyline(vec![Point3::ZERO, Point3::new(0.6, 0.0, 0.0)]),
                r,
                None,
                None,
            )
        },
        &options(budget),
    );
    let bound = bound_of(&report);
    assert_eq!(bound, budget);
    let mut samples = Vec::new();
    for x in span(0.0, 0.6, 600) {
        for t in span(0.0, TAU, 2000) {
            samples.push(Point3::new(x, r * t.cos(), r * t.sin()));
        }
    }
    // The pipe chords its disk to half the budget and its straight piece
    // spends none of the other half, so the proof is twice the truth here.
    assert_bounded("straight pipe", &mesh, &samples, bound, 2.5);
}

#[test]
fn primitive_cylinders_and_cones_are_within_their_reported_bound() {
    for (name, primitive, top) in [
        (
            "cylinder",
            Primitive::Cylinder {
                radius: 0.25,
                height: 0.4,
            },
            0.25,
        ),
        (
            "cone",
            Primitive::Cone {
                radius: 0.25,
                height: 0.4,
            },
            0.0,
        ),
    ] {
        let (mesh, report) = compile(
            |g| g.push(GeometryNode::Primitive(primitive)).unwrap(),
            &options(1e-3),
        );
        assert_eq!(report.contributions[0].detail, name);
        let bound = bound_of(&report);
        let mut samples = Vec::new();
        for z in span(0.0, 0.4, 200) {
            let r = 0.25 + (top - 0.25) * z / 0.4;
            for phi in span(0.0, TAU, 4000) {
                samples.push(Point3::new(r * phi.cos(), r * phi.sin(), z));
            }
        }
        for phi in span(0.0, TAU, 2000) {
            samples.push(Point3::new(0.25 * phi.cos(), 0.25 * phi.sin(), 0.0));
        }
        assert_bounded(name, &mesh, &samples, bound, 2.0);
    }
}

#[test]
fn an_instance_scales_the_bound_and_a_boolean_is_unbounded() {
    let budget = 1e-3;
    let (_, report) = compile(
        |g| {
            let body = extrude(
                g,
                Profile::Ellipse(EllipseProfile {
                    semi_axis_x: 0.3,
                    semi_axis_y: 0.1,
                }),
                0.1,
            );
            g.push(GeometryNode::Instance(Instance {
                source: body,
                transform: Transform3::from_scale(Vec3::splat(2.0)),
            }))
            .unwrap()
        },
        &options(budget),
    );
    // The local budget shrinks by the stretch, the bound grows back by it.
    let bound = bound_of(&report);
    assert!((bound - budget).abs() <= 1e-15, "{bound}");

    let (_, report) = compile(
        |g| {
            let a = g
                .push(GeometryNode::Primitive(Primitive::Block {
                    x: 1.0,
                    y: 1.0,
                    z: 1.0,
                }))
                .unwrap();
            let b = g
                .push(GeometryNode::Primitive(Primitive::Sphere { radius: 0.6 }))
                .unwrap();
            g.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: a,
                right: b,
                operator: axiolid_core::BooleanOperator::Difference,
            }))
            .unwrap()
        },
        &options(budget),
    );
    assert_eq!(report.bound, None);
    assert_eq!(report.contributions[0].path, DeviationPath::Boolean);
    // No exact result to measure against (#235): named by the refusal.
    assert_eq!(
        report.contributions[0].bound,
        DeviationBound::Unbounded("exact boolean operand that is not an extrusion")
    );
}
