//! Every point of a doubly curved surface lies within the chord budget of
//! its mesh (#231).
//!
//! A mesh whose vertices sit on the surface is not enough: the triangles
//! between them must stay within the budget too. Revolving a chorded
//! profile about an axis chords the surface twice, once across the profile
//! and once round the axis, and the two deviations add inside a triangle.
//! Each case samples the exact surface densely, measures each sample's
//! distance to the nearest mesh triangle, and asserts the largest is at or
//! below the requested budget, with no slack: the compiler's bound is a
//! proof, not an estimate, so sampling can only under-report it.
//!
//! The exact surfaces are written out here in closed form (tori, spheres,
//! the tube round an arc, the revolved rounded rectangle), never read back
//! from the compiler.

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{
    Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Transform2, Vec2, Vec3,
};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2, Line3, Polyline3};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation,
    Transition, TrimSelector, TrimmingPreference,
};
use axiolid_primitive::Primitive;
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
};

/// The two budgets the issue names: 1 mm and 0.1 mm.
fn budgets() -> [(ExecutionOptions, Scalar); 3] {
    [
        (ExecutionOptions::new(Tolerance::MILLIMETRE), 1e-3),
        (
            ExecutionOptions::new(Tolerance::new(1e-4, 1e-9).unwrap()),
            1e-4,
        ),
        // An explicit chord budget finer than the linear tolerance governs
        // the surface too, not only the profile.
        (
            ExecutionOptions::new(Tolerance::MILLIMETRE)
                .with_chord_error(1e-4)
                .unwrap(),
            1e-4,
        ),
    ]
}

fn compile(
    build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId,
    options: &ExecutionOptions,
) -> TriMesh {
    let mut builder = GeometryGraphBuilder::new();
    let root = build(&mut builder);
    let graph = builder.finish(vec![root]).expect("a valid graph");
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, root, options)
        .expect("the body compiles")
}

/// Triangles binned on a uniform grid, each in every cell its box touches
/// once grown by `reach`, so one cell lookup finds every triangle within
/// `reach` of a point in that cell.
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
        let cell = (hi - lo).max_element() / 64.0;
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

    /// Distance to the nearest triangle, or `reach` when none is closer.
    fn distance(&self, p: Point3) -> Scalar {
        let mut best = self.reach;
        for &i in self.cells.get(&key(p, self.cell)).into_iter().flatten() {
            // A zero-area triangle covers nothing its neighbours' shared
            // edges do not.
            if let Ok(q) = closest_point_on_triangle(p, self.triangles[i]) {
                best = best.min(p.distance(q));
            }
        }
        best
    }

    fn max_distance(&self, samples: &[Point3]) -> Scalar {
        samples
            .iter()
            .map(|p| self.distance(*p))
            .fold(0.0, Scalar::max)
    }
}

fn key(p: Point3, cell: Scalar) -> [i64; 3] {
    [p.x, p.y, p.z].map(|v| (v / cell).floor() as i64)
}

/// Assert the largest sampled surface-to-mesh distance is within `budget`.
fn assert_within(name: &str, mesh: &TriMesh, samples: &[Point3], budget: Scalar) -> Scalar {
    // Only distances up to twice the budget are resolved: anything beyond
    // reads as twice the budget, which fails the assertion just the same.
    let grid = TriangleGrid::new(mesh, 2.0 * budget);
    let worst = grid.max_distance(samples);
    eprintln!(
        "{name}: budget {budget:e}, max surface->mesh {worst:e}, {} triangles",
        mesh.indices.len() / 3
    );
    assert!(
        worst <= budget,
        "{name}: a surface point lies {worst:e} from the mesh, above the {budget:e} budget"
    );
    worst
}

/// `n` samples of `[a, b]`, offset off the lattice so they do not land on
/// mesh vertices, where the deviation is zero.
fn span(a: Scalar, b: Scalar, n: usize) -> impl Iterator<Item = Scalar> + Clone {
    (0..n).map(move |i| a + (b - a) * ((i as Scalar) + 0.5) / n as Scalar)
}

/// Rodrigues rotation of `p` about the line through `origin` along unit `dir`.
fn rotate(p: Point3, origin: Point3, dir: Vec3, angle: Scalar) -> Point3 {
    let v = p - origin;
    let (s, c) = angle.sin_cos();
    origin + v * c + dir.cross(v) * s + dir * (dir.dot(v) * (1.0 - c))
}

/// Samples of the surface swept by revolving profile-plane points
/// `profile` (a closed boundary sampled densely) and, for a partial turn,
/// the two caps (`cap` samples the profile's region).
fn revolved_samples(
    boundary: &[Point2],
    cap: &[Point2],
    axis_origin: Point3,
    axis: Vec3,
    angle: Scalar,
    turns: usize,
) -> Vec<Point3> {
    let lift = |p: &Point2| Point3::new(p.x, p.y, 0.0);
    let mut out = Vec::new();
    for t in span(0.0, angle, turns) {
        out.extend(
            boundary
                .iter()
                .map(|p| rotate(lift(p), axis_origin, axis, t)),
        );
    }
    if (angle.abs() - TAU).abs() > 1e-9 {
        for t in [0.0, angle] {
            out.extend(cap.iter().map(|p| rotate(lift(p), axis_origin, axis, t)));
        }
    }
    out
}

/// A circle of radius `r` about the profile origin, densely, and its disk.
fn circle_boundary(r: Scalar, n: usize) -> Vec<Point2> {
    span(0.0, TAU, n)
        .map(|a| Point2::new(r * a.cos(), r * a.sin()))
        .collect()
}

fn disk(r: Scalar, rings: usize, n: usize) -> Vec<Point2> {
    let mut out = Vec::new();
    for rho in span(0.0, r, rings) {
        out.extend(span(0.0, TAU, n).map(|a| Point2::new(rho * a.cos(), rho * a.sin())));
    }
    out
}

fn revolution(
    b: &mut GeometryGraphBuilder,
    profile: Profile,
    axis_origin: Point3,
    angle: Scalar,
) -> NodeId {
    revolution_about(b, profile, axis_origin, Vec3::Y, angle)
}

fn revolution_about(
    b: &mut GeometryGraphBuilder,
    profile: Profile,
    axis_origin: Point3,
    axis_direction: Vec3,
    angle: Scalar,
) -> NodeId {
    let profile = b.push(GeometryNode::Profile(profile)).unwrap();
    b.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
        profile,
        axis_origin,
        axis_direction,
        angle,
    }))
    .unwrap()
}

#[test]
fn a_revolution_about_an_axis_off_the_profile_plane_stays_within_the_budget() {
    // The axis leans out of the profile's plane, so the walls are
    // hyperboloidal quads, not flat trapezoids: their twist counts too.
    let r = 0.1;
    let axis_origin = Point3::new(-0.5, 0.0, 0.08);
    let axis = Vec3::new(0.0, 1.0, 0.4).normalize();
    let samples = revolved_samples(&circle_boundary(r, 360), &[], axis_origin, axis, TAU, 720);
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                revolution_about(
                    b,
                    Profile::Circle(CircleProfile {
                        radius: r,
                        thickness: None,
                    }),
                    axis_origin,
                    axis,
                    TAU,
                )
            },
            &options,
        );
        assert_within("revolution about a skew axis", &mesh, &samples, budget);
    }
}

fn rect(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The rectangle's corners in the order its rings list them.
fn corners(x: Scalar, y: Scalar) -> [Point2; 4] {
    let (hx, hy) = (x / 2.0, y / 2.0);
    [
        Point2::new(-hx, -hy),
        Point2::new(hx, -hy),
        Point2::new(hx, hy),
        Point2::new(-hx, hy),
    ]
}

#[test]
fn a_tapered_revolution_stays_within_the_budget() {
    // The section shrinks from 0.2 x 0.3 to 0.1 x 0.12 while it turns:
    // each corner runs on a spiral, not a circle, and the walls twist.
    let (start, end) = (corners(0.2, 0.3), corners(0.1, 0.12));
    let axis_origin = Point3::new(-0.3, 0.0, 0.0);
    let angle = 1.5;
    let lerp =
        |p: Point2, q: Point2, t: Scalar| Point2::new(p.x + (q.x - p.x) * t, p.y + (q.y - p.y) * t);
    let mut samples = Vec::new();
    for t in span(0.0, 1.0, 720) {
        for j in 0..4 {
            let (a, b) = (
                lerp(start[j], end[j], t),
                lerp(start[(j + 1) % 4], end[(j + 1) % 4], t),
            );
            for x in span(0.0, 1.0, 90) {
                let p = lerp(a, b, x);
                samples.push(rotate(
                    Point3::new(p.x, p.y, 0.0),
                    axis_origin,
                    Vec3::Y,
                    angle * t,
                ));
            }
        }
    }
    // The two caps: each profile's region at its own end of the turn.
    for ((x, y), turn) in [((0.2, 0.3), 0.0), ((0.1, 0.12), angle)] {
        for px in span(-0.5 * x, 0.5 * x, 30) {
            for py in span(-0.5 * y, 0.5 * y, 30) {
                samples.push(rotate(Point3::new(px, py, 0.0), axis_origin, Vec3::Y, turn));
            }
        }
    }
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                let start_profile = b.push(GeometryNode::Profile(rect(0.2, 0.3))).unwrap();
                let end_profile = b.push(GeometryNode::Profile(rect(0.1, 0.12))).unwrap();
                b.push(GeometryNode::SolidOperation(
                    SolidOperation::TaperedRevolution {
                        start_profile,
                        end_profile,
                        axis_origin,
                        axis_direction: Vec3::Y,
                        angle,
                    },
                ))
                .unwrap()
            },
            &options,
        );
        assert_within("tapered revolution", &mesh, &samples, budget);
    }
}

#[test]
fn a_section_sliding_outward_while_it_turns_stays_within_the_budget() {
    // A slab slides 0.3 m away from the axis over the turn: its walls
    // stay flat (every quad is planar), but each corner runs on a spiral
    // whose curvature its outward speed adds to. A short turn with a long
    // slide makes that speed, not the turn, set the step.
    let (x, y, slide) = (0.05, 0.2, 0.3);
    let start = corners(x, y);
    let end = start.map(|p| Point2::new(p.x + slide, p.y));
    let axis_origin = Point3::new(-0.3, 0.0, 0.0);
    let angle = 0.3;
    let lerp =
        |p: Point2, q: Point2, t: Scalar| Point2::new(p.x + (q.x - p.x) * t, p.y + (q.y - p.y) * t);
    let mut samples = Vec::new();
    for t in span(0.0, 1.0, 720) {
        for j in 0..4 {
            let (a, b) = (
                lerp(start[j], end[j], t),
                lerp(start[(j + 1) % 4], end[(j + 1) % 4], t),
            );
            for s in span(0.0, 1.0, 90) {
                let p = lerp(a, b, s);
                samples.push(rotate(
                    Point3::new(p.x, p.y, 0.0),
                    axis_origin,
                    Vec3::Y,
                    angle * t,
                ));
            }
        }
    }
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                let start_profile = b.push(GeometryNode::Profile(rect(x, y))).unwrap();
                let end_profile = b
                    .push(GeometryNode::Profile(Profile::Derived {
                        basis: Box::new(rect(x, y)),
                        transform: Transform2::from_translation(Vec2::new(slide, 0.0)),
                    }))
                    .unwrap();
                b.push(GeometryNode::SolidOperation(
                    SolidOperation::TaperedRevolution {
                        start_profile,
                        end_profile,
                        axis_origin,
                        axis_direction: Vec3::Y,
                        angle,
                    },
                ))
                .unwrap()
            },
            &options,
        );
        assert_within("sliding tapered revolution", &mesh, &samples, budget);
    }
}

#[test]
fn a_rectangle_revolved_about_a_skew_axis_stays_within_the_budget() {
    // Straight edges about an axis off their plane sweep hyperboloids:
    // the walls are twisted quads whose twist, not the sagitta, governs.
    let (x, y) = (0.4, 0.05);
    let axis_origin = Point3::new(-0.4, 0.0, 0.1);
    let axis = Vec3::new(0.0, 1.0, 0.5).normalize();
    let rim = corners(x, y);
    let mut boundary = Vec::new();
    for j in 0..4 {
        let (p, q) = (rim[j], rim[(j + 1) % 4]);
        boundary.extend(
            span(0.0, 1.0, 180).map(|s| Point2::new(p.x + (q.x - p.x) * s, p.y + (q.y - p.y) * s)),
        );
    }
    let samples = revolved_samples(&boundary, &[], axis_origin, axis, TAU, 720);
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| revolution_about(b, rect(x, y), axis_origin, axis, TAU),
            &options,
        );
        assert_within("rectangle about a skew axis", &mesh, &samples, budget);
    }
}

#[test]
fn a_budget_beyond_the_step_cap_is_refused_not_coarsened() {
    // 4096 steps round the axis cannot hold a 0.6 m equator to 1 nm; the
    // compiler must say so rather than hand back a coarser mesh.
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(1e-9)
        .unwrap();
    let mut builder = GeometryGraphBuilder::new();
    let torus = revolution(
        &mut builder,
        Profile::Circle(CircleProfile {
            radius: 0.1,
            thickness: None,
        }),
        Point3::new(-0.5, 0.0, 0.0),
        TAU,
    );
    let sphere = builder
        .push(GeometryNode::Primitive(Primitive::Sphere { radius: 1e4 }))
        .unwrap();
    let graph = builder.finish(vec![torus, sphere]).unwrap();
    let compiler = ReferenceMeshCompiler::new(BoolmeshBoolean::new());
    for (root, options) in [
        (torus, options),
        (sphere, ExecutionOptions::new(Tolerance::MILLIMETRE)),
    ] {
        let refused = compiler.compile_mesh(&graph, root, &options);
        assert!(
            matches!(refused, Err(GeomError::BudgetExceeded { .. })),
            "{refused:?}"
        );
    }
}

/// A circle profile of radius `r` revolved about the parallel axis `big`
/// away: a torus, or a torus segment for a partial turn.
fn check_revolved_torus(big: Scalar, r: Scalar, angle: Scalar) {
    let axis_origin = Point3::new(-big, 0.0, 0.0);
    let boundary = circle_boundary(r, 360);
    let cap = disk(r, 12, 90);
    let samples = revolved_samples(&boundary, &cap, axis_origin, Vec3::Y, angle, 720);
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                revolution(
                    b,
                    Profile::Circle(CircleProfile {
                        radius: r,
                        thickness: None,
                    }),
                    axis_origin,
                    angle,
                )
            },
            &options,
        );
        assert_within(
            &format!("revolved torus R {big} r {r} angle {angle:.3}"),
            &mesh,
            &samples,
            budget,
        );
    }
}

#[test]
fn a_revolved_torus_stays_within_the_budget() {
    // The issue's reproduction: R = 0.5, r = 0.1 deviated 1.44 mm at 1 mm.
    check_revolved_torus(0.5, 0.1, TAU);
}

#[test]
fn a_fat_revolved_torus_stays_within_the_budget() {
    // R close to r: the inner side has strongly negative Gaussian
    // curvature, the outer equator twice the radius of the tube centre.
    check_revolved_torus(0.12, 0.1, TAU);
}

#[test]
fn a_thin_revolved_torus_stays_within_the_budget() {
    check_revolved_torus(1.0, 0.05, TAU);
}

#[test]
fn partial_revolutions_stay_within_the_budget() {
    // Caps included: a partial turn is closed by the chorded profile.
    check_revolved_torus(0.3, 0.2, 1.0);
    check_revolved_torus(0.5, 0.1, 1.5 * PI);
}

/// A half disk of radius `r` against the axis: an arc and its diameter,
/// a general contour profile whose revolution is a sphere.
fn half_disk(r: Scalar) -> Profile {
    let arc = ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(0.0, 0.0),
                x: Vec2::new(1.0, 0.0),
                y: Vec2::new(0.0, 1.0),
            },
            radius: r,
        }),
        domain: Interval {
            start: -FRAC_PI_2,
            end: FRAC_PI_2,
        },
        same_sense: true,
    };
    let diameter = ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: Point2::new(0.0, r),
            direction: Vec2::new(0.0, -2.0 * r),
        }),
        domain: Interval::UNIT,
        same_sense: true,
    };
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![arc, diameter]),
        holes: Vec::new(),
    })
}

#[test]
fn a_revolved_arc_profile_sphere_stays_within_the_budget() {
    let r = 0.3;
    let boundary: Vec<Point2> = span(-FRAC_PI_2, FRAC_PI_2, 360)
        .map(|a| Point2::new(r * a.cos(), r * a.sin()))
        .collect();
    for angle in [TAU, 1.2] {
        let mut cap = Vec::new();
        for rho in span(0.0, r, 12) {
            cap.extend(
                span(-FRAC_PI_2, FRAC_PI_2, 45).map(|a| Point2::new(rho * a.cos(), rho * a.sin())),
            );
        }
        let samples = revolved_samples(&boundary, &cap, Point3::ZERO, Vec3::Y, angle, 720);
        for (options, budget) in budgets() {
            let mesh = compile(
                |b| revolution(b, half_disk(r), Point3::ZERO, angle),
                &options,
            );
            assert_within(
                &format!("revolved half disk r {r} angle {angle:.3}"),
                &mesh,
                &samples,
                budget,
            );
        }
    }
}

/// The boundary of an `x` by `y` rectangle with corners rounded to `f`.
fn rounded_rectangle_boundary(x: Scalar, y: Scalar, f: Scalar, n: usize) -> Vec<Point2> {
    let (hx, hy) = (x / 2.0 - f, y / 2.0 - f);
    let mut out = Vec::new();
    for (cx, cy, start) in [
        (hx, hy, 0.0),
        (-hx, hy, FRAC_PI_2),
        (-hx, -hy, PI),
        (hx, -hy, 1.5 * PI),
    ] {
        out.extend(
            span(start, start + FRAC_PI_2, n)
                .map(|a| Point2::new(cx + f * a.cos(), cy + f * a.sin())),
        );
    }
    for s in span(-1.0, 1.0, n) {
        out.push(Point2::new(x / 2.0, s * hy));
        out.push(Point2::new(-x / 2.0, s * hy));
        out.push(Point2::new(s * hx, y / 2.0));
        out.push(Point2::new(s * hx, -y / 2.0));
    }
    out
}

#[test]
fn a_revolved_rounded_rectangle_stays_within_the_budget() {
    // Fillet arcs in a general profile sweep toroidal bands.
    let (x, y, f) = (0.3, 0.2, 0.06);
    let axis_origin = Point3::new(-0.2, 0.0, 0.0);
    let boundary = rounded_rectangle_boundary(x, y, f, 90);
    let samples = revolved_samples(&boundary, &[], axis_origin, Vec3::Y, TAU, 720);
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                revolution(
                    b,
                    Profile::Rectangle(RectangleProfile {
                        x,
                        y,
                        thickness: None,
                        outer_radius: Some(f),
                        inner_radius: None,
                    }),
                    axis_origin,
                    TAU,
                )
            },
            &options,
        );
        assert_within("revolved rounded rectangle", &mesh, &samples, budget);
    }
}

fn sphere_samples(r: Scalar) -> Vec<Point3> {
    let mut out = Vec::new();
    for v in span(0.0, PI, 360) {
        for u in span(0.0, TAU, 720) {
            out.push(Point3::new(
                r * v.sin() * u.cos(),
                r * v.sin() * u.sin(),
                r * v.cos(),
            ));
        }
    }
    out
}

#[test]
fn a_primitive_sphere_stays_within_the_budget() {
    // 5 mm lands on 7 segments round at 1 mm: a stack count halved
    // downwards (3) would make the polar step outgrow its half.
    for r in [0.3, 0.05, 0.005] {
        let samples = sphere_samples(r);
        for (options, budget) in budgets() {
            let mesh = compile(
                |b| {
                    b.push(GeometryNode::Primitive(Primitive::Sphere { radius: r }))
                        .unwrap()
                },
                &options,
            );
            assert_within(&format!("primitive sphere r {r}"), &mesh, &samples, budget);
        }
    }
}

#[test]
fn a_primitive_torus_stays_within_the_budget() {
    for (big, r) in [(0.5, 0.1), (0.12, 0.1)] {
        let mut samples = Vec::new();
        for u in span(0.0, TAU, 720) {
            for v in span(0.0, TAU, 360) {
                let rho = big + r * v.cos();
                samples.push(Point3::new(rho * u.cos(), rho * u.sin(), r * v.sin()));
            }
        }
        for (options, budget) in budgets() {
            let mesh = compile(
                |b| {
                    b.push(GeometryNode::Primitive(Primitive::Torus {
                        major_radius: big,
                        minor_radius: r,
                    }))
                    .unwrap()
                },
                &options,
            );
            assert_within(
                &format!("primitive torus R {big} r {r}"),
                &mesh,
                &samples,
                budget,
            );
        }
    }
}

/// Samples of the tube of radius `r` round the arc of radius `big` in the
/// XY plane from `start` to `end`, with its two end disks.
fn tube_samples(big: Scalar, r: Scalar, start: Scalar, end: Scalar) -> Vec<Point3> {
    let centre = |a: Scalar| Point3::new(big * a.cos(), big * a.sin(), 0.0);
    let normal = |a: Scalar| Vec3::new(a.cos(), a.sin(), 0.0);
    let mut out = Vec::new();
    for a in span(start, end, 720) {
        for v in span(0.0, TAU, 360) {
            out.push(centre(a) + normal(a) * (r * v.cos()) + Vec3::Z * (r * v.sin()));
        }
    }
    for a in [start, end] {
        for rho in span(0.0, r, 12) {
            for v in span(0.0, TAU, 90) {
                out.push(centre(a) + normal(a) * (rho * v.cos()) + Vec3::Z * (rho * v.sin()));
            }
        }
    }
    out
}

fn arc(big: Scalar) -> Curve3 {
    Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: big,
    })
}

#[test]
fn a_disk_swept_along_an_arc_stays_within_the_budget() {
    // Large r / R: the outer side of the tube is r / R further from the
    // arc's centre than the directrix the chords were sized for.
    let (start, end) = (0.2, 0.2 + 1.5 * PI);
    for (big, r) in [(1.0, 0.05), (0.2, 0.1), (0.12, 0.1)] {
        let samples = tube_samples(big, r, start, end);
        for (options, budget) in budgets() {
            let mesh = compile(
                |b| {
                    let directrix = b.push(GeometryNode::Curve3(arc(big))).unwrap();
                    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
                        directrix,
                        radius: r,
                        inner_radius: None,
                        parameter_range: Some((start, end)),
                        fillet_radius: None,
                    }))
                    .unwrap()
                },
                &options,
            );
            assert_within(
                &format!("swept disk R {big} r {r}"),
                &mesh,
                &samples,
                budget,
            );
        }
    }
}

#[test]
fn a_disk_swept_against_a_trimmed_arc_stays_within_the_budget() {
    // The same arc walked backwards through a trim against the basis's
    // sense: its end tangents must turn round with it.
    let (big, r) = (0.2, 0.1);
    let (start, end) = (0.2, 0.2 + 1.5 * PI);
    let samples = tube_samples(big, r, start, end);
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                let basis = b.push(GeometryNode::Curve3(arc(big))).unwrap();
                let directrix = b
                    .push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
                        basis,
                        start: vec![TrimSelector::Parameter(end)],
                        end: vec![TrimSelector::Parameter(start)],
                        sense_agreement: false,
                        preference: TrimmingPreference::Parameter,
                    }))
                    .unwrap();
                b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
                    directrix,
                    radius: r,
                    inner_radius: None,
                    parameter_range: None,
                    fillet_radius: None,
                }))
                .unwrap()
            },
            &options,
        );
        assert_within("disk swept against a trimmed arc", &mesh, &samples, budget);
    }
}

#[test]
fn a_profile_swept_along_an_arc_stays_within_the_budget() {
    // A fixed-reference sweep of a rectangle along an arc, with the arc's
    // plane normal as reference, is a revolution of the rectangle about
    // the arc's axis: its far wall is 0.1 m outside the directrix.
    let (big, x, y) = (0.2, 0.08, 0.2);
    let (start, end) = (0.3, 0.3 + PI);
    let mut samples = Vec::new();
    // Profile x is the reference (Z), profile y the outward normal.
    let rim = corners(x, y);
    let mut boundary = Vec::new();
    for j in 0..4 {
        let (p, q) = (rim[j], rim[(j + 1) % 4]);
        boundary.extend(
            span(0.0, 1.0, 90).map(|s| Point2::new(p.x + (q.x - p.x) * s, p.y + (q.y - p.y) * s)),
        );
    }
    let place = |a: Scalar, p: Point2| {
        Point3::new(big * a.cos(), big * a.sin(), 0.0)
            + Vec3::Z * p.x
            + Vec3::new(a.cos(), a.sin(), 0.0) * p.y
    };
    for a in span(start, end, 720) {
        samples.extend(boundary.iter().map(|p| place(a, *p)));
    }
    for a in [start, end] {
        for px in span(-0.5 * x, 0.5 * x, 30) {
            for py in span(-0.5 * y, 0.5 * y, 30) {
                samples.push(place(a, Point2::new(px, py)));
            }
        }
    }
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                let profile = b
                    .push(GeometryNode::Profile(Profile::Rectangle(
                        RectangleProfile {
                            x,
                            y,
                            thickness: None,
                            outer_radius: None,
                            inner_radius: None,
                        },
                    )))
                    .unwrap();
                let directrix = b.push(GeometryNode::Curve3(arc(big))).unwrap();
                b.push(GeometryNode::SolidOperation(
                    SolidOperation::FixedReferenceSweep {
                        profile,
                        directrix,
                        reference_direction: Vec3::Z,
                        parameter_range: Some((start, end)),
                    },
                ))
                .unwrap()
            },
            &options,
        );
        assert_within("rectangle swept along an arc", &mesh, &samples, budget);
    }
}

// --- #232: pipes along chains of straight legs and bends, primitives ---

/// One smooth piece of a pipe's centreline, in closed form: a straight leg
/// or a circular bend (`from` turned by `angle` about the line through
/// `centre` along unit `axis`).
#[derive(Clone, Copy, Debug)]
enum Leg {
    Straight {
        from: Point3,
        dir: Vec3,
        length: Scalar,
    },
    Bend {
        centre: Point3,
        axis: Vec3,
        from: Point3,
        angle: Scalar,
    },
}

impl Leg {
    fn point(&self, s: Scalar) -> Point3 {
        match *self {
            Leg::Straight { from, dir, length } => from + dir * (length * s),
            Leg::Bend {
                centre,
                axis,
                from,
                angle,
            } => rotate(from, centre, axis, angle * s),
        }
    }

    fn tangent(&self, s: Scalar) -> Vec3 {
        match *self {
            Leg::Straight { dir, .. } => dir,
            Leg::Bend { centre, axis, .. } => axis.cross(self.point(s) - centre).normalize(),
        }
    }

    fn length(&self) -> Scalar {
        match *self {
            Leg::Straight { length, .. } => length,
            Leg::Bend {
                centre,
                from,
                angle,
                ..
            } => (from - centre).length() * angle,
        }
    }
}

/// A centreline drawn as a turtle walks: straight legs and bends.
struct Turtle {
    at: Point3,
    dir: Vec3,
    legs: Vec<Leg>,
}

impl Turtle {
    fn new(at: Point3, dir: Vec3) -> Self {
        Self {
            at,
            dir: dir.normalize(),
            legs: Vec::new(),
        }
    }

    fn straight(mut self, length: Scalar) -> Self {
        self.legs.push(Leg::Straight {
            from: self.at,
            dir: self.dir,
            length,
        });
        self.at += self.dir * length;
        self
    }

    /// Bend through `angle` on a centreline radius `radius`, turning
    /// toward `toward` (made perpendicular to the heading).
    fn bend(mut self, radius: Scalar, angle: Scalar, toward: Vec3) -> Self {
        let inward = (toward - self.dir * self.dir.dot(toward)).normalize();
        let centre = self.at + inward * radius;
        let axis = self.dir.cross(inward).normalize();
        let leg = Leg::Bend {
            centre,
            axis,
            from: self.at,
            angle,
        };
        self.at = leg.point(1.0);
        self.dir = leg.tangent(1.0);
        self.legs.push(leg);
        self
    }

    fn length(&self) -> Scalar {
        self.legs.iter().map(Leg::length).sum()
    }
}

/// Two unit vectors perpendicular to unit `t` and to each other.
fn normals(t: Vec3) -> (Vec3, Vec3) {
    let seed = if t.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let n1 = (seed - t * t.dot(seed)).normalize();
    (n1, t.cross(n1))
}

/// Samples of the exact tube of radius `r` round `legs`: each leg's own
/// wall, square to its own tangent, and the two end disks.
fn pipe_samples(legs: &[Leg], r: Scalar) -> Vec<Point3> {
    let mut out = Vec::new();
    for leg in legs {
        let along = match leg {
            Leg::Straight { .. } => 24,
            Leg::Bend { .. } => 720,
        };
        for s in span(0.0, 1.0, along) {
            let (c, (n1, n2)) = (leg.point(s), normals(leg.tangent(s)));
            out.extend(span(0.0, TAU, 360).map(|v| c + n1 * (r * v.cos()) + n2 * (r * v.sin())));
        }
    }
    let first = legs.first().unwrap();
    let last = legs.last().unwrap();
    for (c, t) in [
        (first.point(0.0), first.tangent(0.0)),
        (last.point(1.0), last.tangent(1.0)),
    ] {
        let (n1, n2) = normals(t);
        for rho in span(0.0, r, 12) {
            out.extend(span(0.0, TAU, 90).map(|v| c + n1 * (rho * v.cos()) + n2 * (rho * v.sin())));
        }
    }
    out
}

/// How the composite writes its legs into the graph.
#[derive(Clone, Copy, Debug)]
enum Spelling {
    /// Trims of a line and of a circle, all with the composite's sense.
    Trims,
    /// Two-point polylines for the legs, and each bend trimmed backwards
    /// against its circle and used against the composite's sense.
    Reversed,
}

/// Push `legs` as one composite curve.
fn composite(b: &mut GeometryGraphBuilder, legs: &[Leg], spelling: Spelling) -> NodeId {
    let mut segments = Vec::new();
    for leg in legs {
        let (curve, same_sense) = match (*leg, spelling) {
            (Leg::Straight { from, dir, length }, Spelling::Trims) => {
                let line = b
                    .push(GeometryNode::Curve3(Curve3::Line(Line3 {
                        origin: from,
                        direction: dir,
                    })))
                    .unwrap();
                (trim(b, line, 0.0, length, true), true)
            }
            (Leg::Straight { from, dir, length }, Spelling::Reversed) => {
                let polyline = Curve3::Polyline(Polyline3 {
                    points: vec![from, from + dir * length],
                    closed: false,
                });
                (b.push(GeometryNode::Curve3(polyline)).unwrap(), true)
            }
            (
                Leg::Bend {
                    centre,
                    axis,
                    from,
                    angle,
                },
                spelling,
            ) => {
                let radius = (from - centre).length();
                let x = (from - centre) / radius;
                let circle = b
                    .push(GeometryNode::Curve3(Curve3::Circle(Circle3 {
                        frame: Frame3 {
                            origin: centre,
                            x,
                            y: axis.cross(x),
                            z: axis,
                        },
                        radius,
                    })))
                    .unwrap();
                match spelling {
                    Spelling::Trims => (trim(b, circle, 0.0, angle, true), true),
                    Spelling::Reversed => (trim(b, circle, angle, 0.0, false), false),
                }
            }
        };
        segments.push(CurveSegment {
            curve,
            same_sense,
            transition: Transition::Continuous,
        });
    }
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments,
    }))
    .unwrap()
}

fn trim(
    b: &mut GeometryGraphBuilder,
    basis: NodeId,
    start: Scalar,
    end: Scalar,
    sense: bool,
) -> NodeId {
    b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
        basis,
        start: vec![TrimSelector::Parameter(start)],
        end: vec![TrimSelector::Parameter(end)],
        sense_agreement: sense,
        preference: TrimmingPreference::Parameter,
    }))
    .unwrap()
}

fn swept_disk(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    radius: Scalar,
    fillet_radius: Option<Scalar>,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix,
        radius,
        inner_radius: None,
        parameter_range: None,
        fillet_radius,
    }))
    .unwrap()
}

fn try_compile(
    build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId,
    options: &ExecutionOptions,
) -> Result<TriMesh, GeomError> {
    let mut builder = GeometryGraphBuilder::new();
    let root = build(&mut builder);
    let graph = builder.finish(vec![root]).expect("a valid graph");
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, root, options)
}

/// Every edge is used by exactly two triangles, once in each direction:
/// closed, two-manifold and consistently oriented.
fn assert_watertight(name: &str, mesh: &TriMesh) {
    let mut uses: HashMap<(u32, u32), (usize, i64)> = HashMap::new();
    for t in mesh.indices.chunks_exact(3) {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let entry = uses.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(
        uses.values().all(|(n, d)| *n == 2 && *d == 0),
        "{name}: the mesh is not closed and consistently oriented"
    );
}

/// The volume is positive (outward) and within `area * budget` of the
/// exact tube's, `pi r^2` times the centreline length (Pappus for each
/// bend), never above it: the mesh is inscribed.
fn assert_pipe_volume(name: &str, mesh: &TriMesh, r: Scalar, length: Scalar, budget: Scalar) {
    let got = axiolid_measure::volume_properties(mesh, Tolerance::new(1e-9, 1e-9).unwrap())
        .unwrap_or_else(|e| panic!("{name}: not a closed two-manifold: {e:?}"))
        .signed_volume;
    let want = PI * r * r * length;
    let area = TAU * r * length + 2.0 * PI * r * r;
    assert!(
        got > 0.0 && got <= want * (1.0 + 1e-12) && want - got <= area * budget,
        "{name}: volume {got} against Pappus {want} (area {area}, budget {budget:e})"
    );
}

fn check_pipe(name: &str, legs: &[Leg], length: Scalar, r: Scalar, spelling: Spelling) {
    let samples = pipe_samples(legs, r);
    for (options, budget) in budgets() {
        let mesh = compile(
            |b| {
                let path = composite(b, legs, spelling);
                swept_disk(b, path, r, None)
            },
            &options,
        );
        let name = format!("{name} r {r} ({spelling:?})");
        assert_watertight(&name, &mesh);
        assert_pipe_volume(&name, &mesh, r, length, budget);
        assert_within(&name, &mesh, &samples, budget);
    }
}

#[test]
fn pipes_along_line_arc_line_composites_stay_within_the_budget() {
    // Bend radii of one, two and three pipe diameters, a thin pipe and a
    // fat one: the outer side of the bend sets the step, not the centreline.
    for (big, r, spelling) in [
        (0.1, 0.05, Spelling::Trims),
        (0.04, 0.01, Spelling::Reversed),
        (0.06, 0.01, Spelling::Trims),
        (0.6, 0.3, Spelling::Reversed),
    ] {
        let path = Turtle::new(Point3::new(0.1, -0.2, 0.3), Vec3::X)
            .straight(0.5)
            .bend(big, FRAC_PI_2, Vec3::Y)
            .straight(0.4);
        check_pipe(
            &format!("line + bend R {big} + line"),
            &path.legs,
            path.length(),
            r,
            spelling,
        );
    }
}

#[test]
fn a_pipe_bending_out_of_its_plane_stays_within_the_budget() {
    // Three bends in three planes, one of them through more than a half
    // turn, two of them back to back; the frames carry round all of them.
    let r = 0.05;
    let path = Turtle::new(Point3::ZERO, Vec3::new(1.0, 0.2, 0.0))
        .straight(0.3)
        .bend(0.1, FRAC_PI_2, Vec3::Y)
        .straight(0.2)
        .bend(0.15, 1.2, Vec3::Z)
        .bend(0.1, 0.6 * PI, Vec3::new(1.0, -1.0, 0.5))
        .straight(0.25);
    for spelling in [Spelling::Trims, Spelling::Reversed] {
        check_pipe(
            "pipe bending out of plane",
            &path.legs,
            path.length(),
            r,
            spelling,
        );
    }
}

/// The corners of the polyline whose corners, rounded to `fillet`, give
/// `legs` (straight legs alternating with bends of radius `fillet`): each
/// corner lies where the legs either side of a bend meet when extended.
fn corners_of(legs: &[Leg], fillet: Scalar) -> Vec<Point3> {
    let mut points = vec![legs[0].point(0.0)];
    for leg in legs {
        if let Leg::Bend { angle, .. } = *leg {
            let reach = fillet * (0.5 * angle).tan();
            points.push(leg.point(0.0) + leg.tangent(0.0) * reach);
        }
    }
    points.push(legs.last().unwrap().point(1.0));
    points
}

fn filleted(b: &mut GeometryGraphBuilder, corners: &[Point3], r: Scalar, fillet: Scalar) -> NodeId {
    let polyline = b
        .push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
            points: corners.to_vec(),
            closed: false,
        })))
        .unwrap();
    swept_disk(b, polyline, r, Some(fillet))
}

#[test]
fn filleted_polylines_stay_within_the_budget() {
    // IfcSweptDiskSolidPolygonal: a polyline whose corners are rounded to
    // the fillet radius, two bends in a plane and three out of it, two of
    // them back to back where their fillets use up the whole leg between.
    let cases = [
        (
            0.05,
            0.1,
            Turtle::new(Point3::new(0.0, 0.0, 0.0), Vec3::X)
                .straight(0.4)
                .bend(0.1, FRAC_PI_2, Vec3::Y)
                .straight(0.3)
                .bend(0.1, FRAC_PI_2, -Vec3::X)
                .straight(0.4),
        ),
        (
            0.02,
            0.06,
            Turtle::new(Point3::new(0.2, 0.1, -0.1), Vec3::new(1.0, 1.0, 0.0))
                .straight(0.3)
                .bend(0.06, 1.1, Vec3::Z)
                .straight(0.2)
                .bend(0.06, FRAC_PI_2, Vec3::new(0.3, -1.0, 0.0))
                .straight(0.0)
                .bend(0.06, 0.7, Vec3::new(-1.0, 0.0, 1.0))
                .straight(0.35),
        ),
    ];
    for (r, fillet, path) in cases {
        let legs: Vec<Leg> = path
            .legs
            .iter()
            .copied()
            .filter(|leg| leg.length() > 0.0)
            .collect();
        let corners = corners_of(&path.legs, fillet);
        let samples = pipe_samples(&legs, r);
        for (options, budget) in budgets() {
            let mesh = compile(|b| filleted(b, &corners, r, fillet), &options);
            let name = format!(
                "polyline of {} corners filleted to {fillet}, r {r}",
                corners.len() - 2
            );
            assert_watertight(&name, &mesh);
            assert_pipe_volume(&name, &mesh, r, path.length(), budget);
            assert_within(&name, &mesh, &samples, budget);
        }
    }
}

/// The largest turn a joint may make and still count as tangent
/// continuous (`axiolid_construct::pipe`): its shared section, square to
/// the incoming leg, then leans off the outgoing leg's own by at most a
/// quarter of the chord budget.
fn joint_tolerance(chord: Scalar, r: Scalar) -> Scalar {
    2.0 * (chord / (8.0 * r)).min(1.0).asin()
}

/// A straight leg that meets a bend at its start turned by `kink` about
/// the bend's axis, then the bend, then a leg straight off its end.
fn kinked(kink: Scalar, r_bend: Scalar) -> Vec<Leg> {
    let bend = Turtle::new(Point3::new(0.4, 0.0, 0.0), Vec3::X)
        .bend(r_bend, FRAC_PI_2, Vec3::Y)
        .straight(0.3);
    let dir = Vec3::new(kink.cos(), -kink.sin(), 0.0);
    let mut legs = vec![Leg::Straight {
        from: Point3::new(0.4, 0.0, 0.0) - dir * 0.4,
        dir,
        length: 0.4,
    }];
    legs.extend(bend.legs);
    legs
}

#[test]
fn a_joint_within_the_stated_angular_tolerance_is_swept_within_the_budget() {
    // A joint that turns by less than the tolerance is tangent continuous:
    // every leg's own tube stays within the budget of the mesh.
    let r = 0.05;
    for (options, budget) in budgets() {
        let legs = kinked(0.9 * joint_tolerance(budget, r), 0.12);
        let length = legs.iter().map(Leg::length).sum();
        let samples = pipe_samples(&legs, r);
        let mesh = compile(
            |b| {
                let path = composite(b, &legs, Spelling::Trims);
                swept_disk(b, path, r, None)
            },
            &options,
        );
        assert_watertight("kink within tolerance", &mesh);
        assert_pipe_volume("kink within tolerance", &mesh, r, length, budget);
        assert_within("kink within tolerance", &mesh, &samples, budget);
        // Just past the tolerance it is a corner.
        let legs = kinked(1.1 * joint_tolerance(budget, r), 0.12);
        let refused = try_compile(
            |b| {
                let path = composite(b, &legs, Spelling::Trims);
                swept_disk(b, path, r, None)
            },
            &options,
        );
        assert!(is_invalid(&refused, "corner"), "{}", outcome(&refused));
    }
}

#[test]
fn a_gap_at_a_joint_is_paid_for_by_the_next_piece() {
    // A short bend sized so its two spans, alone, use nearly all of their
    // half of the 1 mm budget at the outer side ((R + r)(1 - cos) = 0.499
    // mm), and whose start lies 0.24 mm outward of the incoming leg's end,
    // within the linear tolerance the composite accepts. The shared joint
    // section sits on the leg, so the bend's own first span must make room
    // for the measured shift; otherwise its outer side misses the budget.
    let (r, big, gap): (Scalar, Scalar, Scalar) = (0.05, 0.15, 2.4e-4);
    let step = 2.0 * (1.0 - 0.499e-3 / (big + r)).acos();
    let bend = Turtle::new(Point3::new(0.4, -gap, 0.0), Vec3::X)
        .bend(big, 2.0 * step, Vec3::Y)
        .straight(0.2);
    let mut legs = vec![Leg::Straight {
        from: Point3::ZERO,
        dir: Vec3::X,
        length: 0.4,
    }];
    legs.extend(bend.legs);
    let samples = pipe_samples(&legs, r);
    let (options, budget) = budgets()[0].clone();
    let mesh = compile(
        |b| {
            let path = composite(b, &legs, Spelling::Trims);
            swept_disk(b, path, r, None)
        },
        &options,
    );
    assert_watertight("gap at a joint", &mesh);
    assert_within("gap at a joint", &mesh, &samples, budget);
}

/// A compile result for a failure message, without dumping the mesh.
fn outcome(result: &Result<TriMesh, GeomError>) -> String {
    match result {
        Ok(mesh) => format!("compiled to {} triangles", mesh.indices.len() / 3),
        Err(error) => format!("{error:?}"),
    }
}

fn is_invalid(result: &Result<TriMesh, GeomError>, words: &str) -> bool {
    matches!(result, Err(GeomError::InvalidInput(m)) if m.contains(words))
}

fn polyline(b: &mut GeometryGraphBuilder, points: &[Point3], closed: bool) -> NodeId {
    b.push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
        points: points.to_vec(),
        closed,
    })))
    .unwrap()
}

// --- #245: corners mitred at half angle ---

/// A swept disk along `directrix`, hollow when `inner` is given.
fn swept_tube(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    radius: Scalar,
    inner: Option<Scalar>,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
        directrix,
        radius,
        inner_radius: inner,
        parameter_range: None,
        fillet_radius: None,
    }))
    .unwrap()
}

/// The mitre plane's unit normal at corner `k` of `points` (between legs
/// `k - 1` and `k`), or the leg's own tangent at an open end: the
/// bisector of the two tangents.
fn mitre_normal(points: &[Point3], k: usize) -> Vec3 {
    let leg = |i: usize| (points[i + 1] - points[i]).normalize();
    let legs = points.len() - 1;
    match (k.checked_sub(1), k < legs) {
        (Some(i), true) => (leg(i) + leg(k)).normalize(),
        (Some(i), false) => leg(i),
        (None, _) => leg(0),
    }
}

/// Samples of the exact tube round the polyline `points`, mitred at half
/// angle at every corner, for each radius in `radii` (the outer wall and a
/// bore): each leg's cylinder cut by its two end planes (the cut curves,
/// the ellipses both tubes share, included), and the two end annuli.
fn mitred_samples(points: &[Point3], radii: &[Scalar]) -> Vec<Point3> {
    let mut out = Vec::new();
    for k in 0..points.len() - 1 {
        let (start, end) = (points[k], points[k + 1]);
        let u = (end - start).normalize();
        let length = (end - start).length();
        let (n0, n1) = (mitre_normal(points, k), mitre_normal(points, k + 1));
        let (a, b) = normals(u);
        for &r in radii {
            for v in span(0.0, TAU, 360) {
                let w = a * (r * v.cos()) + b * (r * v.sin());
                let s0 = -w.dot(n0) / u.dot(n0);
                let s1 = length - w.dot(n1) / u.dot(n1);
                let along = span(0.0, 1.0, 24).chain([0.0, 1.0]);
                out.extend(along.map(|t| start + w + u * (s0 + (s1 - s0) * t)));
            }
        }
    }
    let (lo, hi) = (
        radii.iter().copied().fold(Scalar::INFINITY, Scalar::min),
        radii[0],
    );
    let lo = if radii.len() > 1 { lo } else { 0.0 };
    let last = points.len() - 1;
    for (c, t) in [
        (points[0], mitre_normal(points, 0)),
        (points[last], mitre_normal(points, last)),
    ] {
        let (a, b) = normals(t);
        for rho in span(lo, hi, 12) {
            out.extend(span(0.0, TAU, 90).map(|v| c + a * (rho * v.cos()) + b * (rho * v.sin())));
        }
    }
    out
}

/// The area of the regular `n`-gon inscribed in a circle of radius `r`.
fn polygon_area(n: usize, r: Scalar) -> Scalar {
    0.5 * n as Scalar * r * r * (TAU / n as Scalar).sin()
}

/// The number of mesh vertices at distance `r` from `at`: the ring of the
/// section there.
fn ring_size(mesh: &TriMesh, at: Point3, r: Scalar) -> usize {
    mesh.positions
        .iter()
        .filter(|p| ((**p - at).length() - r).abs() <= 1e-12)
        .count()
}

/// The closed form for a mitred tube: each leg's cylinder cut by planes
/// through its two end points on the centreline holds `pi r^2 L`, since
/// a plane through the axis cuts as much off one side of the square
/// section as it adds on the other. The mesh, a prism over a regular ring
/// polygon centred on the axis cut by the same planes, holds exactly the
/// polygon's area times the centreline length; checked to rounding, with
/// one station per vertex of the polyline (the mitre rings shared). The
/// exact tube's volume is then above it by at most its area times the
/// budget.
fn assert_mitred_volume(
    name: &str,
    mesh: &TriMesh,
    points: &[Point3],
    r: Scalar,
    inner: Option<Scalar>,
    budget: Scalar,
) {
    let got = axiolid_measure::volume_properties(mesh, Tolerance::new(1e-9, 1e-9).unwrap())
        .unwrap_or_else(|e| panic!("{name}: not a closed two-manifold: {e:?}"))
        .signed_volume;
    let length: Scalar = points.windows(2).map(|p| (p[1] - p[0]).length()).sum();
    let outer = ring_size(mesh, points[0], r);
    let bore = inner.map_or(0, |ri| ring_size(mesh, points[0], ri));
    assert!(outer >= 3, "{name}: no ring of radius {r} at the start");
    assert_eq!(
        mesh.positions.len(),
        points.len() * (outer + bore),
        "{name}: one station per polyline vertex"
    );
    let polygon = polygon_area(outer, r) - inner.map_or(0.0, |ri| polygon_area(bore, ri));
    let prism = polygon * length;
    assert!(
        (got - prism).abs() <= 1e-9 * prism,
        "{name}: volume {got} against the mitred prism's {prism}"
    );
    let ri = inner.unwrap_or(0.0);
    let exact = PI * (r * r - ri * ri) * length;
    let area = TAU * (r + ri) * length + 2.0 * PI * (r * r - ri * ri);
    eprintln!("{name}: volume {got:.12e}, mitred tube {exact:.12e}, prism {prism:.12e}");
    assert!(
        got <= exact && exact - got <= area * budget,
        "{name}: volume {got} against the mitred tube's {exact} (area {area}, budget {budget:e})"
    );
}

/// The polylines of the issue: corners of 90, 30 and 150 degrees, and one
/// leaving its plane at every corner.
fn cornered_polylines() -> Vec<(&'static str, Vec<Point3>)> {
    let turn = |deg: Scalar| {
        let a = deg.to_radians();
        vec![
            Point3::new(0.1, -0.2, 0.3),
            Point3::new(0.7, -0.2, 0.3),
            Point3::new(0.7 + 0.6 * a.cos(), -0.2 + 0.6 * a.sin(), 0.3),
        ]
    };
    vec![
        ("90 degree corner", turn(90.0)),
        ("30 degree corner", turn(30.0)),
        ("150 degree corner", turn(150.0)),
        (
            "corners out of plane",
            vec![
                Point3::ZERO,
                Point3::new(0.6, 0.0, 0.0),
                Point3::new(0.6, 0.5, 0.0),
                Point3::new(0.6, 0.5, 0.5),
                Point3::new(1.1, 0.9, 0.7),
                Point3::new(0.6, 1.2, 0.2),
            ],
        ),
    ]
}

#[test]
fn polylines_with_sharp_corners_are_mitred_within_the_budget() {
    // IfcSweptDiskSolid along a polyline with no fillet radius: each corner
    // is mitred at half angle, both legs cut by the bisector plane.
    for (name, points) in cornered_polylines() {
        for (r, inner) in [(0.05, None), (0.1, Some(0.07))] {
            let radii: Vec<Scalar> = [Some(r), inner].into_iter().flatten().collect();
            let samples = mitred_samples(&points, &radii);
            for (options, budget) in budgets() {
                let mesh = compile(
                    |b| {
                        let path = polyline(b, &points, false);
                        swept_tube(b, path, r, inner)
                    },
                    &options,
                );
                let name = format!("{name}, r {r}, bore {inner:?}");
                assert_watertight(&name, &mesh);
                assert_mitred_volume(&name, &mesh, &points, r, inner, budget);
                assert_within(&name, &mesh, &samples, budget);
            }
        }
    }
}

#[test]
fn a_composite_of_two_lines_meeting_at_a_corner_is_mitred() {
    let r = 0.05;
    for degrees in [10.0, 90.0] {
        let turn = Scalar::to_radians(degrees);
        let points = [
            Point3::ZERO,
            Point3::X,
            Point3::X + Vec3::new(turn.cos(), turn.sin(), 0.0),
        ];
        let legs = [
            Leg::Straight {
                from: points[0],
                dir: Vec3::X,
                length: 1.0,
            },
            Leg::Straight {
                from: points[1],
                dir: Vec3::new(turn.cos(), turn.sin(), 0.0),
                length: 1.0,
            },
        ];
        let samples = mitred_samples(&points, &[r]);
        for spelling in [Spelling::Trims, Spelling::Reversed] {
            for (options, budget) in budgets() {
                let mesh = compile(
                    |b| {
                        let path = composite(b, &legs, spelling);
                        swept_disk(b, path, r, None)
                    },
                    &options,
                );
                let name = format!("composite corner of {degrees} degrees ({spelling:?})");
                assert_watertight(&name, &mesh);
                assert_mitred_volume(&name, &mesh, &points, r, None, budget);
                assert_within(&name, &mesh, &samples, budget);
            }
        }
    }
}

#[test]
fn a_mitre_reports_the_budget_it_is_proven_to() {
    let (options, budget) = budgets()[0].clone();
    let points = &cornered_polylines()[3].1;
    let mut builder = GeometryGraphBuilder::new();
    let path = polyline(&mut builder, points, false);
    let root = swept_tube(&mut builder, path, 0.05, Some(0.03));
    let graph = builder.finish(vec![root]).expect("a valid graph");
    let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options)
        .expect("the body compiles");
    assert!(report.meets_requested(), "{report:?}");
    assert!(report.bound.is_some_and(|b| b <= budget), "{report:?}");
    assert!(
        report
            .contributions
            .iter()
            .any(|c| c.path == DeviationPath::SweptDisk
                && matches!(c.bound, DeviationBound::Proven(b) if b <= budget)),
        "{report:?}"
    );
}

#[test]
fn a_cut_through_a_polyline_corner_keeps_the_corner() {
    // The downstream case: a pipe along an L, clipped by a plane just off
    // its centreline. The clip keeps the mitred corner: the outer wall of
    // each leg runs on to the mitre plane, `r` past the corner along each
    // leg (a section square to the bisector would stop at r / sqrt 2).
    let (options, budget) = budgets()[0].clone();
    let (r, h) = (0.05, 0.01);
    let points = [Point3::ZERO, Point3::X, Point3::new(1.0, 1.0, 0.0)];
    let mesh = compile(
        |b| {
            let path = polyline(b, &points, false);
            let pipe = swept_tube(b, path, r, None);
            let plane = b
                .push(GeometryNode::HalfSpace(axiolid_primitive::HalfSpace {
                    boundary: axiolid_core::Plane3 {
                        origin: Point3::new(0.0, 0.0, h),
                        normal: Vec3::Z,
                    },
                    agreement: true,
                }))
                .unwrap();
            b.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
                left: pipe,
                right: plane,
                operator: axiolid_core::BooleanOperator::Difference,
            }))
            .unwrap()
        },
        &options,
    );
    assert_watertight("clipped corner", &mesh);
    let (lo, hi) = mesh.positions.iter().fold(
        (
            Vec3::splat(Scalar::INFINITY),
            Vec3::splat(Scalar::NEG_INFINITY),
        ),
        |(lo, hi), p| (lo.min(*p), hi.max(*p)),
    );
    eprintln!("clipped corner: box {lo:?} .. {hi:?}");
    assert!(hi.x >= 1.0 + r - budget && hi.x <= 1.0 + r, "{hi:?}");
    assert!(lo.y <= -r + budget && lo.y >= -r, "{lo:?}");
    assert!(hi.z <= h + 1e-9 || lo.z >= h - 1e-9, "{lo:?} .. {hi:?}");
    // Each leg keeps the part of its section on one side of the plane,
    // so the volume is that segment's area times the centreline length
    // (the mitre's cut is odd across the plane of the corner).
    let cap = r * r * (h / r).acos() - h * (r * r - h * h).sqrt();
    let kept = if hi.z <= h + 1e-9 {
        PI * r * r - cap
    } else {
        cap
    };
    let want = kept * 2.0;
    let got = axiolid_measure::volume_properties(&mesh, Tolerance::new(1e-9, 1e-9).unwrap())
        .unwrap()
        .signed_volume;
    let area = TAU * r * 2.0 + 4.0 * PI * r * r;
    assert!(
        got > 0.0 && (want - got).abs() <= area * budget,
        "clipped corner: volume {got} against {want}"
    );
}

/// The refusal's message, for assertions on its words.
fn refused(points: &[Point3], r: Scalar) -> Result<TriMesh, GeomError> {
    try_compile(
        |b| {
            let path = polyline(b, points, false);
            swept_disk(b, path, r, None)
        },
        &ExecutionOptions::new(Tolerance::MILLIMETRE),
    )
}

#[test]
fn only_impossible_mitres_are_refused_by_name() {
    let r = 0.05;
    // A 150 degree corner reaches r tan(75 degrees) = 0.187 along each leg.
    let a = 150f64.to_radians();
    let sharp = |leg: Scalar| {
        [
            Point3::ZERO,
            Point3::X,
            Point3::X + Vec3::new(a.cos(), a.sin(), 0.0) * leg,
        ]
    };
    let fits = refused(&sharp(0.19), r);
    assert!(fits.is_ok(), "{}", outcome(&fits));
    let short = refused(&sharp(0.18), r);
    assert!(
        is_invalid(&short, "cut through itself"),
        "{}",
        outcome(&short)
    );
    // A U: both mitres of the middle leg cut into its inner side, 2 r in
    // all.
    let u_turn = |leg: Scalar| {
        [
            Point3::ZERO,
            Point3::X,
            Point3::new(1.0, leg, 0.0),
            Point3::new(0.0, leg, 0.0),
        ]
    };
    let fits = refused(&u_turn(0.11), r);
    assert!(fits.is_ok(), "{}", outcome(&fits));
    let short = refused(&u_turn(0.09), r);
    assert!(
        is_invalid(&short, "cut through itself"),
        "{}",
        outcome(&short)
    );
    // A reversal has no half angle to mitre at.
    let back = refused(&[Point3::ZERO, Point3::X, 0.5 * Point3::X], r);
    assert!(is_invalid(&back, "reverses"), "{}", outcome(&back));
    // A corner beside an arc: no ring lies on both the cylinder's and the
    // torus's cut, so it is refused rather than mitred.
    let bend = Turtle::new(Point3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0))
        .bend(0.2, FRAC_PI_2, -Vec3::X)
        .straight(0.5);
    let mut legs = vec![Leg::Straight {
        from: Point3::ZERO,
        dir: Vec3::X,
        length: 1.0,
    }];
    legs.extend(bend.legs);
    let beside = try_compile(
        |b| {
            let path = composite(b, &legs, Spelling::Trims);
            swept_disk(b, path, r, None)
        },
        &ExecutionOptions::new(Tolerance::MILLIMETRE),
    );
    assert!(is_invalid(&beside, "beside an arc"), "{}", outcome(&beside));
    // A closed polyline would need a mitre where it closes: not built.
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    for fillet in [None, Some(0.2)] {
        let closed = try_compile(
            |b| {
                let path = polyline(b, &sharp(1.0), true);
                swept_disk(b, path, r, fillet)
            },
            &options,
        );
        assert!(
            matches!(closed, Err(GeomError::UnsupportedInput { input, .. }) if input.contains("closed polyline")),
            "{fillet:?}: {}",
            outcome(&closed)
        );
    }
    // Collinear vertices are no corner.
    let straight = refused(&[Point3::ZERO, Point3::X, 2.0 * Point3::X], r);
    assert!(straight.is_ok(), "{}", outcome(&straight));
}

#[test]
fn a_fillet_that_does_not_fit_or_folds_the_tube_is_refused_by_name() {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE);
    let corners = [
        Point3::ZERO,
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 0.3, 0.0),
        Point3::new(2.0, 0.3, 0.0),
    ];
    let attempt =
        |r: Scalar, fillet: Scalar| try_compile(|b| filleted(b, &corners, r, fillet), &options);
    // Two right-angle corners 0.3 apart each need `fillet` of the leg
    // between them: 0.15 fits exactly, 0.2 does not.
    let fits = attempt(0.05, 0.15);
    assert!(fits.is_ok(), "{}", outcome(&fits));
    let too_large = attempt(0.05, 0.2);
    assert!(
        is_invalid(&too_large, "does not fit"),
        "{}",
        outcome(&too_large)
    );
    // A disk wider than the fillet folds the inside of the bend.
    let folded = attempt(0.12, 0.1);
    assert!(is_invalid(&folded, "fillet radius"), "{}", outcome(&folded));
    // As wide as the fillet, which the format rule permits (fillet radius
    // at least the disk radius), the bend is a horn torus whose inner wall
    // pinches to a point: refused by name, citing that rule (#245).
    let horn = attempt(0.1, 0.1);
    assert!(
        is_invalid(&horn, "equals the fillet radius") && is_invalid(&horn, "horn torus"),
        "{}",
        outcome(&horn)
    );
    // A bend in a composite as tight as the disk folds it too.
    let path = Turtle::new(Point3::ZERO, Vec3::X)
        .straight(0.3)
        .bend(0.05, FRAC_PI_2, Vec3::Y)
        .straight(0.3);
    let tight = try_compile(
        |b| {
            let directrix = composite(b, &path.legs, Spelling::Trims);
            swept_disk(b, directrix, 0.05, None)
        },
        &options,
    );
    assert!(is_invalid(&tight, "bend radius"), "{}", outcome(&tight));
}

#[test]
fn a_pipe_budget_beyond_the_step_cap_is_refused() {
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(1e-9)
        .unwrap();
    let path = Turtle::new(Point3::ZERO, Vec3::X)
        .straight(0.3)
        .bend(0.6, FRAC_PI_2, Vec3::Y)
        .straight(0.3);
    let refused = try_compile(
        |b| {
            let directrix = composite(b, &path.legs, Spelling::Trims);
            swept_disk(b, directrix, 0.3, None)
        },
        &options,
    );
    assert!(
        matches!(refused, Err(GeomError::BudgetExceeded { .. })),
        "{}",
        outcome(&refused)
    );
}

/// Samples of a cylinder (`top` = `r`) or a cone (`top` = 0) along +z
/// from z = 0 to `h`: the side and both caps.
fn frustum_samples(r: Scalar, top: Scalar, h: Scalar) -> Vec<Point3> {
    let mut out = Vec::new();
    for z in span(0.0, h, 64) {
        let rho = r + (top - r) * z / h;
        out.extend(span(0.0, TAU, 1440).map(|a| Point3::new(rho * a.cos(), rho * a.sin(), z)));
    }
    for (rho_max, z) in [(r, 0.0), (top, h)] {
        for rho in span(0.0, rho_max, 16) {
            out.extend(span(0.0, TAU, 720).map(|a| Point3::new(rho * a.cos(), rho * a.sin(), z)));
        }
    }
    out
}

#[test]
fn primitive_cylinders_and_cones_stay_within_the_budget() {
    let h = 0.5;
    for r in [0.01, 0.3, 1.0] {
        for (primitive, top, what) in [
            (
                Primitive::Cylinder {
                    radius: r,
                    height: h,
                },
                r,
                "cylinder",
            ),
            (
                Primitive::Cone {
                    radius: r,
                    height: h,
                },
                0.0,
                "cone",
            ),
        ] {
            let samples = frustum_samples(r, top, h);
            for (options, budget) in budgets() {
                let mesh = compile(
                    |b| b.push(GeometryNode::Primitive(primitive)).unwrap(),
                    &options,
                );
                assert_within(&format!("primitive {what} r {r}"), &mesh, &samples, budget);
            }
        }
    }
}

#[test]
fn a_primitive_cylinder_or_cone_beyond_the_segment_cap_is_refused() {
    // 4096 segments hold a 1 m circle to 0.3 um, not to 1 nm: refuse
    // rather than clamp to a coarser mesh.
    let options = ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(1e-9)
        .unwrap();
    for primitive in [
        Primitive::Cylinder {
            radius: 1.0,
            height: 1.0,
        },
        Primitive::Cone {
            radius: 1.0,
            height: 1.0,
        },
    ] {
        let refused = try_compile(
            |b| b.push(GeometryNode::Primitive(primitive)).unwrap(),
            &options,
        );
        assert!(
            matches!(refused, Err(GeomError::BudgetExceeded { .. })),
            "{primitive:?}: {}",
            outcome(&refused)
        );
    }
}
