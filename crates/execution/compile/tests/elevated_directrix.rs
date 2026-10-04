//! Sweeps along elevated and banked directrices (#252): a gradient curve
//! lowers to `Curve3::Elevated`, a canted one to `Curve3::Banked`.
//!
//! Swept disks and swept areas along straight, arc and clothoid plans with
//! line, circular-arc (R = 1000) and parabolic profiles compile to closed
//! two-manifold meshes; a disk's tube is certified against the exact tube,
//! and dense samples of that tube (written here from the reference
//! evaluator, never read back from the compiler) lie within the reported
//! bound at 1 mm and 0.1 mm. Volumes match `A L3D`, `L3D` the directrix's
//! 3D length (a tube's volume is its section times its centreline's length
//! while it does not overlap itself, Pappus's theorem for sections normal
//! to the curve); on a straight plan at a constant grade that is the closed
//! form `pi r^2 L sqrt(1 + g^2)`.

use std::collections::HashMap;
use std::f64::consts::{PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve2, BankConvention, Banked3, CantLaw, CantPiece, Circle2, CurvatureLaw, Curve2,
    Curve3, Elevated3, ElevationLaw, Intrinsic2, KnotSpec, Line2,
};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    exact_directrix, DeviationBound, DeviationPath, DeviationReport, ExactDirectrix,
    ReferenceMeshCompiler,
};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation, SurfaceRelation};
use axiolid_profile::{Profile, RectangleProfile};

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(budget, 1e-9).unwrap())
        .with_chord_error(budget)
        .unwrap()
}

fn compile(
    build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId,
    options: &ExecutionOptions,
) -> Result<(TriMesh, DeviationReport), GeomError> {
    let mut builder = GeometryGraphBuilder::new();
    let root = build(&mut builder);
    let graph = builder.finish(vec![root]).expect("a valid graph");
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, options)
        .map(|(outcome, report)| (outcome.mesh, report))
}

// --- the alignments ----------------------------------------------------------

fn world() -> Frame2 {
    Frame2 {
        origin: Point2::new(20.0, 10.0),
        x: Vec2::new(0.8, 0.6),
        y: Vec2::new(-0.6, 0.8),
    }
}

fn straight() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::new(20.0, 10.0),
        direction: Vec2::new(0.8, 0.6),
    })
}

fn arc() -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: world(),
        radius: 150.0,
    })
}

fn clothoid() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        world(),
        CurvatureLaw::clothoid(0.0, 1.0 / 120.0, 60.0),
        60.0,
    ))
}

fn plans() -> Vec<(&'static str, Curve2)> {
    vec![
        ("straight", straight()),
        ("arc", arc()),
        ("clothoid", clothoid()),
    ]
}

fn profiles() -> Vec<(&'static str, ElevationLaw)> {
    vec![
        ("line", ElevationLaw::constant_grade(50.0, 0.04)),
        // The consumer's pipe: a vertical circular sag of R = 1000.
        ("sag", ElevationLaw::circular_arc(50.0, -0.03, 1000.0)),
        ("parabola", ElevationLaw::parabolic(50.0, 0.05, -0.03, 60.0)),
    ]
}

fn elevated(plan: &Curve2, profile: &ElevationLaw) -> Curve3 {
    Curve3::Elevated(Elevated3::new(plan.clone(), profile.clone()))
}

fn banked() -> Curve3 {
    Curve3::Banked(Banked3::new(
        Elevated3::new(clothoid(), ElevationLaw::parabolic(50.0, 0.02, -0.01, 60.0)),
        CantLaw::new(vec![
            CantPiece::sine(30.0, 0.0, 0.12),
            CantPiece::constant(30.0, 0.12),
        ]),
        // Rotation about the low rail: the rotation point rises with the
        // cant through the transition.
        CantLaw::new(vec![
            CantPiece::sine(30.0, 0.0, 0.06),
            CantPiece::constant(30.0, 0.06),
        ]),
        1.5,
        BankConvention::TangentRotation,
    ))
}

const LENGTH: Scalar = 60.0;

// --- graph building ----------------------------------------------------------

fn disk(
    b: &mut GeometryGraphBuilder,
    curve: Curve3,
    radius: Scalar,
    range: Option<(Scalar, Scalar)>,
) -> NodeId {
    let directrix = b.push(GeometryNode::Curve3(curve)).unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix,
        radius,
        inner_radius: None,
        parameter_range: range,
        fillet_radius: None,
    }))
    .unwrap()
}

fn rectangle(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn area(b: &mut GeometryGraphBuilder, curve: Curve3, range: Option<(Scalar, Scalar)>) -> NodeId {
    let profile = b.push(GeometryNode::Profile(rectangle(0.8, 0.5))).unwrap();
    let directrix = b.push(GeometryNode::Curve3(curve)).unwrap();
    b.push(GeometryNode::SolidOperation(
        SolidOperation::FixedReferenceSweep {
            profile,
            directrix,
            reference_direction: Vec3::Z,
            parameter_range: range,
        },
    ))
    .unwrap()
}

// --- measurement -------------------------------------------------------------

/// Triangles binned on a uniform grid, each in every cell its box touches
/// once grown by `reach`.
struct TriangleGrid {
    triangles: Vec<[Point3; 3]>,
    cells: HashMap<[i64; 3], Vec<usize>>,
    cell: Scalar,
    reach: Scalar,
}

fn key(p: Point3, cell: Scalar) -> [i64; 3] {
    [p.x, p.y, p.z].map(|v| (v / cell).floor() as i64)
}

impl TriangleGrid {
    fn new(mesh: &TriMesh, reach: Scalar, cell: Scalar) -> Self {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .chunks_exact(3)
            .map(|t| [0, 1, 2].map(|i| mesh.positions[t[i] as usize]))
            .collect();
        let cell = cell.max(4.0 * reach);
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

fn span(a: Scalar, b: Scalar, n: usize) -> impl Iterator<Item = Scalar> + Clone {
    (0..n).map(move |i| a + (b - a) * ((i as Scalar) + 0.5) / n as Scalar)
}

/// Samples of the exact tube of `radius` about `curve` over `[a, b]` and
/// its two end disks, square to the curve.
fn tube_samples(curve: &Curve3, (a, b): (Scalar, Scalar), radius: Scalar, n: usize) -> Vec<Point3> {
    let frame = |t: Scalar| {
        let c = axiolid_reference::curve::evaluate3(curve, t).unwrap();
        let d = axiolid_reference::curve::derivative3(curve, t)
            .unwrap()
            .normalize();
        let r = if d.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
        let n = (r - d * d.dot(r)).normalize();
        (c, n, d.cross(n))
    };
    let mut out = Vec::new();
    for t in span(a, b, n) {
        let (c, nn, bb) = frame(t);
        for theta in span(0.0, TAU, 180) {
            out.push(c + (nn * theta.cos() + bb * theta.sin()) * radius);
        }
    }
    for t in [a, b] {
        let (c, nn, bb) = frame(t);
        for r in span(0.0, radius, 20) {
            for theta in span(0.0, TAU, 180) {
                out.push(c + (nn * theta.cos() + bb * theta.sin()) * r);
            }
        }
    }
    out
}

/// The worst distance of `samples` from the mesh; panics above `bound`.
fn assert_within(name: &str, mesh: &TriMesh, samples: &[Point3], bound: Scalar) -> Scalar {
    let grid = TriangleGrid::new(mesh, 2.0 * bound, 0.25);
    let worst = samples
        .iter()
        .map(|p| grid.distance(*p))
        .fold(0.0, Scalar::max);
    eprintln!(
        "{name}: reported {bound:e}, measured {worst:e}, ratio {:.2}, {} triangles",
        bound / worst,
        mesh.indices.len() / 3
    );
    assert!(
        worst <= bound,
        "{name}: a point of the exact tube lies {worst:e} from the mesh, above {bound:e}"
    );
    worst
}

/// The directrix's 3D length over `[a, b]`: Gauss-Legendre on `|c'|`.
fn length3(curve: &Curve3, a: Scalar, b: Scalar) -> Scalar {
    let nodes = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let panels = 2000;
    let h = (b - a) / panels as Scalar;
    (0..panels)
        .map(|k| {
            let mid = a + h * (k as Scalar + 0.5);
            nodes
                .iter()
                .map(|(x, w)| {
                    w * 0.5
                        * h
                        * axiolid_reference::curve::derivative3(curve, mid + 0.5 * h * x)
                            .unwrap()
                            .length()
                })
                .sum::<Scalar>()
        })
        .sum()
}

fn volume(mesh: &TriMesh, name: &str) -> Scalar {
    volume_properties(mesh, Tolerance::MILLIMETRE)
        .unwrap_or_else(|e| panic!("{name}: the mesh is not a closed two-manifold: {e}"))
        .signed_volume
        .abs()
}

// --- swept disks -------------------------------------------------------------

/// Sweep a disk of `radius` along `curve` over `range` at `budget`, and
/// check it: closed, the contribution named `detail`, every dense sample
/// of the exact tube within the reported bound, the bound not above 1.5
/// times the budget (the sweep's stations are placed on a second-order
/// estimate that may use the whole budget, and the certificate settles
/// once within it, so a bound a little over is honest) nor more than 2.5
/// times the measured maximum, and the volume `pi r^2 L3D` (Pappus) up to
/// the inscribed rings.
fn assert_tube(
    name: &str,
    curve: &Curve3,
    radius: Scalar,
    budget: Scalar,
    range: (Scalar, Scalar),
    detail: &str,
) -> DeviationReport {
    let started = std::time::Instant::now();
    let (mesh, report) = compile(
        |g| disk(g, curve.clone(), radius, Some(range)),
        &options(budget),
    )
    .unwrap_or_else(|e| panic!("{name}: {e}"));
    eprintln!("{name}: compiled and certified in {:?}", started.elapsed());
    let contribution = &report.contributions[0];
    assert_eq!(contribution.path, DeviationPath::SweptDisk, "{name}");
    assert_eq!(contribution.detail, detail, "{name}");
    let bound = report.bound.unwrap_or_else(|| panic!("{name}: {report:?}"));
    let samples = tube_samples(curve, range, radius, 1000);
    let worst = assert_within(name, &mesh, &samples, bound);
    // Both end sections stand square to the exact curve, not to the end
    // chords: every vertex of an end ring lies in the normal plane there.
    for end in [range.0, range.1] {
        let c = axiolid_reference::curve::evaluate3(curve, end).unwrap();
        let t = axiolid_reference::curve::derivative3(curve, end)
            .unwrap()
            .normalize();
        let ring: Vec<Point3> = mesh
            .positions
            .iter()
            .copied()
            .filter(|p| p.distance(c) <= radius * (1.0 + 1e-9))
            .collect();
        assert!(ring.len() >= 3, "{name}: no end ring at {end}");
        let lean = ring
            .iter()
            .map(|p| (*p - c).dot(t).abs())
            .fold(0.0, Scalar::max);
        assert!(lean <= 1e-9, "{name}: end ring at {end} leans {lean:e}");
    }
    assert!(bound <= 1.5 * budget, "{name}: {bound} over {budget}");
    assert!(bound <= 2.5 * worst, "{name}: {bound} against {worst}");
    let exact = PI * radius * radius * length3(curve, range.0, range.1);
    let v = volume(&mesh, name);
    assert!(
        v <= exact * (1.0 + 1e-9) && v >= exact * (1.0 - 2.0 * budget / radius),
        "{name}: volume {v} against {exact}"
    );
    report
}

/// Every plan under every profile: watertight, certified (proven on the
/// straight constant grade, a segment), within the bound when sampled
/// densely, and the tube's volume.
#[test]
fn a_disk_swept_along_an_elevated_directrix_is_certified() {
    for (plan_name, plan) in plans() {
        for (profile_name, profile) in profiles() {
            let name = format!("{plan_name}/{profile_name}");
            let straight = plan_name == "straight" && profile_name == "line";
            let detail = if straight {
                "segment"
            } else {
                "elevated curve"
            };
            let report = assert_tube(
                &name,
                &elevated(&plan, &profile),
                0.1,
                1e-3,
                (10.0, 40.0),
                detail,
            );
            let bound = &report.contributions[0].bound;
            if straight {
                assert_eq!(*bound, DeviationBound::Proven(1e-3), "{name}");
            } else {
                assert!(matches!(bound, DeviationBound::Certified(_)), "{name}");
            }
        }
    }
}

/// The consumer's case at 0.1 mm: a pipe along a clothoid plan under the
/// R = 1000 sag, and along an arc plan under the parabola.
#[test]
fn an_elevated_pipe_is_certified_at_a_tenth_of_a_millimetre() {
    for (name, plan, profile) in [
        (
            "clothoid/sag 0.1 mm",
            clothoid(),
            ElevationLaw::circular_arc(50.0, -0.03, 1000.0),
        ),
        (
            "arc/parabola 0.1 mm",
            arc(),
            ElevationLaw::parabolic(50.0, 0.05, -0.03, 60.0),
        ),
    ] {
        assert_tube(
            name,
            &elevated(&plan, &profile),
            0.05,
            1e-4,
            (20.0, 30.0),
            "elevated curve",
        );
    }
}

/// A straight plan at a constant grade: a slanted cylinder, whose volume
/// is the closed form `pi r^2 L sqrt(1 + g^2)` up to the inscribed ring,
/// and exactly the ring's area times that length.
#[test]
fn a_constant_grade_on_a_straight_plan_has_the_closed_form_volume() {
    let (radius, grade) = (0.3, 0.06);
    let curve = Curve3::Elevated(Elevated3::new(
        straight(),
        ElevationLaw::constant_grade(10.0, grade),
    ));
    let budget = 1e-4;
    let (mesh, report) = compile(
        |g| disk(g, curve.clone(), radius, Some((5.0, 45.0))),
        &options(budget),
    )
    .unwrap();
    assert_eq!(report.contributions[0].detail, "segment");
    let length = 40.0 * grade.hypot(1.0);
    let exact = PI * radius * radius * length;
    let v = volume(&mesh, "straight constant grade");
    assert!(
        (v - exact).abs() <= exact * 2.0 * budget / radius,
        "volume {v} against {exact}"
    );
    // The mesh's ring: a regular polygon, every vertex on the circle.
    let ring = mesh.positions.len() / 2;
    let polygon = 0.5 * ring as Scalar * radius * radius * (TAU / ring as Scalar).sin();
    assert!(
        (v - polygon * length).abs() <= 1e-9 * v,
        "volume {v} against the prism {}",
        polygon * length
    );
    // The exact compiler reads the same directrix as one segment.
    let mut builder = GeometryGraphBuilder::new();
    let id = builder.push(GeometryNode::Curve3(curve.clone())).unwrap();
    let graph = builder.finish(vec![id]).unwrap();
    let ExactDirectrix::Segment(a, b) =
        exact_directrix(&graph, id, Some((5.0, 45.0)), &options(budget)).unwrap()
    else {
        panic!("a segment");
    };
    assert!((a.distance(b) - length).abs() < 1e-9);
    assert!((a.z - (10.0 + 5.0 * grade)).abs() < 1e-12);
}

/// A banked directrix: the disk follows the rotation point, which rises
/// with the pivot through the transition; certified like an elevated one.
#[test]
fn a_disk_swept_along_a_banked_directrix_is_certified() {
    let curve = banked();
    assert_tube("banked", &curve, 0.1, 1e-3, (5.0, 40.0), "banked curve");
    assert_tube(
        "banked 0.1 mm",
        &curve,
        0.05,
        1e-4,
        (25.0, 35.0),
        "banked curve",
    );
    // The tube is around the rotation point, `pivot` above the profile.
    let Curve3::Banked(b) = &curve else {
        unreachable!()
    };
    let lifted = axiolid_reference::curve::evaluate3(&curve, 45.0).unwrap();
    let base = axiolid_reference::elevated_point(&b.base, 45.0).unwrap();
    assert!((lifted.z - base.z - 0.06).abs() < 1e-12);
}

// --- swept areas -------------------------------------------------------------

/// A rectangle swept with a fixed reference along every elevated
/// directrix: closed, and its volume `A L3D` (sections normal to the
/// curve, centroid on it). The straight constant grade is a segment, so
/// its deviation is proven; along a curve the frame law is not bounded,
/// which the report names.
#[test]
fn an_area_swept_along_an_elevated_directrix_is_closed() {
    let budget = 1e-3;
    for (plan_name, plan) in plans() {
        for (profile_name, profile) in profiles() {
            let name = format!("{plan_name}/{profile_name}");
            let curve = elevated(&plan, &profile);
            let range = Some((0.0, LENGTH));
            let (mesh, report) = compile(|g| area(g, curve.clone(), range), &options(budget))
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            let contribution = &report.contributions[0];
            assert_eq!(contribution.path, DeviationPath::FixedReferenceSweep);
            if plan_name == "straight" && profile_name == "line" {
                assert_eq!(contribution.bound, DeviationBound::Proven(0.0), "{name}");
            } else {
                assert!(
                    matches!(
                        contribution.bound,
                        DeviationBound::Unbounded(reason) if reason.contains("frame law")
                    ),
                    "{name}: {:?}",
                    contribution.bound
                );
            }
            let exact = 0.8 * 0.5 * length3(&curve, 0.0, LENGTH);
            let v = volume(&mesh, &name);
            assert!(
                (v - exact).abs() <= exact * 1e-4,
                "{name}: volume {v} against {exact}"
            );
        }
    }
}

/// Along a banked directrix the fixed-reference frame is the sweep's own
/// law, not the cant's roll: the area still sweeps closed.
#[test]
fn an_area_swept_along_a_banked_directrix_is_closed() {
    let curve = banked();
    let (mesh, _) = compile(|g| area(g, curve.clone(), None), &options(1e-3)).unwrap();
    let exact = 0.8 * 0.5 * length3(&curve, 0.0, LENGTH);
    let v = volume(&mesh, "banked area");
    assert!((v - exact).abs() <= exact * 1e-4, "{v} against {exact}");
}

// --- corners and refusals ----------------------------------------------------

/// A grade that jumps at a seam is a corner: the disk sweeps closed, but
/// the deviation is unbounded by name.
#[test]
fn a_grade_break_is_unbounded_by_name() {
    let kinked = Curve3::Elevated(Elevated3::new(
        clothoid(),
        ElevationLaw::Piecewise {
            breaks: vec![30.0],
            laws: vec![
                ElevationLaw::constant_grade(50.0, 0.02),
                ElevationLaw::constant_grade(50.6, 0.06),
            ],
        },
    ));
    let (mesh, report) = compile(|g| disk(g, kinked, 0.2, None), &options(1e-3)).unwrap();
    volume(&mesh, "kinked");
    assert_eq!(
        report.contributions[0].detail,
        "elevated directrix with a grade break"
    );
    assert!(report.bound.is_none());
    assert!(matches!(
        report.contributions[0].bound,
        DeviationBound::Unbounded(_)
    ));
    // A tangent-continuous seam is no corner: line, crest, line.
    let arc = ElevationLaw::circular_arc(50.4, 0.02, -1500.0);
    let (exit, height) = (arc.grade_at(30.0).unwrap(), arc.height_at(30.0).unwrap());
    let smooth = Curve3::Elevated(Elevated3::new(
        clothoid(),
        ElevationLaw::Piecewise {
            breaks: vec![20.0, 50.0],
            laws: vec![
                ElevationLaw::constant_grade(50.0, 0.02),
                arc,
                ElevationLaw::constant_grade(height, exit),
            ],
        },
    ));
    assert_tube(
        "gradient",
        &smooth,
        0.1,
        1e-3,
        (10.0, 60.0),
        "elevated curve",
    );
}

#[test]
fn an_unbounded_or_unevaluable_directrix_is_refused_by_name() {
    let on_line = elevated(&straight(), &ElevationLaw::constant_grade(0.0, 0.02));
    let error = compile(|g| disk(g, on_line.clone(), 0.1, None), &options(1e-3)).unwrap_err();
    assert!(error.to_string().contains("no end"), "{error}");
    let error = compile(|g| area(g, on_line.clone(), None), &options(1e-3)).unwrap_err();
    assert!(error.to_string().contains("no end"), "{error}");
    // Before the plan's start: outside the domain, not extrapolated.
    let error = compile(
        |g| {
            disk(
                g,
                elevated(&clothoid(), &ElevationLaw::level(0.0)),
                0.1,
                Some((-5.0, 20.0)),
            )
        },
        &options(1e-3),
    )
    .unwrap_err();
    assert!(
        error.to_string().contains("outside curve domain"),
        "{error}"
    );
    // A B-spline plan's parameter is not a distance.
    let spline = Curve3::Elevated(Elevated3::new(
        Curve2::BSpline(BSplineCurve2 {
            degree: 1,
            control_points: vec![Point2::new(0.0, 0.0), Point2::new(10.0, 0.0)],
            knots: vec![0.0, 1.0],
            multiplicities: vec![2, 2],
            weights: None,
            closed: false,
            self_intersect: None,
            knot_spec: KnotSpec::Unspecified,
        }),
        ElevationLaw::level(0.0),
    ));
    assert!(compile(|g| disk(g, spline, 0.1, Some((0.0, 1.0))), &options(1e-3)).is_err());
    // A sag that turns vertical inside the span.
    let steep = elevated(&clothoid(), &ElevationLaw::circular_arc(0.0, -0.4, 10.0));
    assert!(compile(|g| disk(g, steep, 0.1, None), &options(1e-3)).is_err());
}

/// A surface-curve sweep along a gradient curve, its reference surface the
/// curve extruded vertically (an alignment's vertical surface): closed,
/// and its volume `A L3D`; the frame law is not bounded, by name.
#[test]
fn an_area_swept_along_an_elevated_directrix_on_its_vertical_surface_is_closed() {
    let curve = elevated(
        &clothoid(),
        &ElevationLaw::circular_arc(50.0, -0.03, 1000.0),
    );
    let (mesh, report) = compile(
        |b| {
            let profile = b.push(GeometryNode::Profile(rectangle(0.8, 0.5))).unwrap();
            let directrix = b.push(GeometryNode::Curve3(curve.clone())).unwrap();
            let surface = b
                .push(GeometryNode::SurfaceRelation(
                    SurfaceRelation::LinearExtrusion {
                        swept_curve: directrix,
                        direction: Vec3::Z,
                    },
                ))
                .unwrap();
            b.push(GeometryNode::SolidOperation(
                SolidOperation::SurfaceCurveSweep {
                    profile,
                    directrix,
                    reference_surface: surface,
                    parameter_range: Some((10.0, 50.0)),
                },
            ))
            .unwrap()
        },
        &options(1e-3),
    )
    .unwrap();
    assert!(matches!(
        report.contributions[0].bound,
        DeviationBound::Unbounded(reason) if reason.contains("frame law")
    ));
    let exact = 0.8 * 0.5 * length3(&curve, 10.0, 50.0);
    let v = volume(&mesh, "surface-curve sweep");
    assert!((v - exact).abs() <= exact * 1e-4, "{v} against {exact}");
}
