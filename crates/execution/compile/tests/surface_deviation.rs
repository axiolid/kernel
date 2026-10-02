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
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Line2};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    CurveRelation, GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation, TrimSelector,
    TrimmingPreference,
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
