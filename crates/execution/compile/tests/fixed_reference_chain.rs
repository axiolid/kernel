//! A profile swept with a fixed reference along a tangent-continuous chain
//! of segments and arcs stays within the chord budget (#122).
//!
//! Each case samples the exact swept surface densely, placing the
//! profile's boundary by the fixed-reference law itself (x is the reference
//! made square to the tangent, z the tangent, y = z x x) at every sampled
//! station of every piece, and asserts the largest distance from a sample
//! to the mesh is within the requested budget. The exact surface is written
//! out here in closed form, never read back from the compiler. A swept
//! solid must also be closed and two-manifold, and a centred profile's
//! volume is its area times the centreline length. What the construction
//! refuses is checked by name, and the deviation report must call the
//! proven chains proven and an oblique arc unbounded.

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Frame3, Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{Circle3, Curve3, Polyline3};
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_model::{
    CurveRelation, CurveSegment, GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation,
    Transition, TrimSelector, TrimmingPreference,
};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

/// The budgets the sweeps are checked at: 1 mm, 0.1 mm, and an explicit
/// chord budget finer than the linear tolerance.
fn budgets() -> [(ExecutionOptions, Scalar); 3] {
    [
        (ExecutionOptions::new(Tolerance::MILLIMETRE), 1e-3),
        (
            ExecutionOptions::new(Tolerance::new(1e-4, 1e-9).unwrap()),
            1e-4,
        ),
        (
            ExecutionOptions::new(Tolerance::MILLIMETRE)
                .with_chord_error(1e-4)
                .unwrap(),
            1e-4,
        ),
    ]
}

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(budget.min(1e-3), 1e-9).unwrap())
        .with_chord_error(budget)
        .unwrap()
}

fn try_compile(
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

fn compile(
    build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId,
    options: &ExecutionOptions,
) -> (TriMesh, DeviationReport) {
    try_compile(build, options).expect("the body compiles")
}

fn refusal(build: impl FnOnce(&mut GeometryGraphBuilder) -> NodeId, needle: &str) -> GeomError {
    let error = try_compile(build, &options(1e-3)).expect_err("the sweep is refused");
    let message = error.to_string();
    assert!(message.contains(needle), "{needle}: {message}");
    error
}

// --- distance from the exact surface to the mesh ---------------------------

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
    fn new(mesh: &TriMesh, reach: Scalar) -> Self {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
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
            if let Ok(q) = closest_point_on_triangle(p, self.triangles[i]) {
                best = best.min(p.distance(q));
            }
        }
        best
    }
}

fn assert_within(name: &str, mesh: &TriMesh, samples: &[Point3], budget: Scalar) {
    let grid = TriangleGrid::new(mesh, 2.0 * budget);
    let worst = samples
        .iter()
        .map(|p| grid.distance(*p))
        .fold(0.0, Scalar::max);
    eprintln!("{name}: budget {budget:e}, max surface->mesh {worst:e}");
    assert!(
        worst <= budget,
        "{name}: a surface point lies {worst:e} from the mesh, above the {budget:e} budget"
    );
}

/// Every edge is used by exactly two triangles, once in each direction.
fn assert_watertight(name: &str, mesh: &TriMesh) {
    let mut uses: HashMap<(u32, u32), (usize, i64)> = HashMap::new();
    for t in mesh.indices.as_chunks::<3>().0 {
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

// --- the exact centreline ---------------------------------------------------

fn rotate(p: Point3, origin: Point3, dir: Vec3, angle: Scalar) -> Point3 {
    let v = p - origin;
    let (s, c) = angle.sin_cos();
    origin + v * c + dir.cross(v) * s + dir * (dir.dot(v) * (1.0 - c))
}

/// One smooth piece in closed form.
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
            Leg::Bend {
                centre,
                axis,
                angle,
                ..
            } => (axis.cross(self.point(s) - centre) * angle.signum()).normalize(),
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
            } => (from - centre).length() * angle.abs(),
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

    /// Bend on radius `radius` through the positive `angle` about `axis`.
    fn bend(mut self, axis: Vec3, radius: Scalar, angle: Scalar) -> Self {
        let centre = self.at + axis.cross(self.dir) * radius;
        let leg = Leg::Bend {
            centre,
            axis,
            from: self.at,
            angle,
        };
        self.legs.push(leg);
        self.at = leg.point(1.0);
        self.dir = leg.tangent(1.0);
        self
    }
}

fn length_of(legs: &[Leg]) -> Scalar {
    legs.iter().map(Leg::length).sum()
}

/// Push `legs` as one composite: polylines for the legs, each bend a
/// trimmed circle.
fn composite(b: &mut GeometryGraphBuilder, legs: &[Leg]) -> NodeId {
    let mut segments = Vec::new();
    for leg in legs {
        let curve = match *leg {
            Leg::Straight { from, dir, length } => b
                .push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
                    points: vec![from, from + dir * length],
                    closed: false,
                })))
                .unwrap(),
            Leg::Bend {
                centre,
                axis,
                from,
                angle,
            } => {
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
                b.push(GeometryNode::CurveRelation(CurveRelation::Trimmed {
                    basis: circle,
                    start: vec![TrimSelector::Parameter(0.0)],
                    end: vec![TrimSelector::Parameter(angle)],
                    sense_agreement: true,
                    preference: TrimmingPreference::Parameter,
                }))
                .unwrap()
            }
        };
        segments.push(CurveSegment {
            curve,
            same_sense: true,
            transition: Transition::Continuous,
        });
    }
    b.push(GeometryNode::CurveRelation(CurveRelation::Composite {
        segments,
    }))
    .unwrap()
}

// --- the profile ------------------------------------------------------------

/// A profile centred on its origin: its boundary and region sampled in
/// closed form.
struct Shape {
    profile: Profile,
    boundary: Vec<Point2>,
    region: Vec<Point2>,
    area: Scalar,
    perimeter: Scalar,
}

fn rect(x: Scalar, y: Scalar) -> Shape {
    let (hx, hy) = (x / 2.0, y / 2.0);
    let rim = [
        Point2::new(-hx, -hy),
        Point2::new(hx, -hy),
        Point2::new(hx, hy),
        Point2::new(-hx, hy),
    ];
    let mut boundary = Vec::new();
    for j in 0..4 {
        let (p, q) = (rim[j], rim[(j + 1) % 4]);
        boundary.extend((0..40).map(|i| {
            let s = (i as Scalar + 0.5) / 40.0;
            Point2::new(p.x + (q.x - p.x) * s, p.y + (q.y - p.y) * s)
        }));
    }
    let region = (0..20)
        .flat_map(|i| {
            (0..20).map(move |j| {
                Point2::new(
                    hx * (2.0 * (i as Scalar + 0.5) / 20.0 - 1.0),
                    hy * (2.0 * (j as Scalar + 0.5) / 20.0 - 1.0),
                )
            })
        })
        .collect();
    Shape {
        profile: Profile::Rectangle(RectangleProfile {
            x,
            y,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        boundary,
        region,
        area: x * y,
        perimeter: 2.0 * (x + y),
    }
}

/// An annulus: arcs for the profile and a hole in it.
fn annulus(radius: Scalar, thickness: Scalar) -> Shape {
    let inner = radius - thickness;
    let on = |r: Scalar, n: usize, offset: Scalar| {
        (0..n)
            .map(|i| {
                let a = TAU * (i as Scalar + offset) / n as Scalar;
                Point2::new(r * a.cos(), r * a.sin())
            })
            .collect::<Vec<_>>()
    };
    let mut boundary = on(radius, 360, 0.5);
    boundary.extend(on(inner, 360, 0.5));
    let region = (0..8)
        .flat_map(|k| on(inner + thickness * (k as Scalar + 0.5) / 8.0, 90, 0.5))
        .collect();
    Shape {
        profile: Profile::Circle(CircleProfile {
            radius,
            thickness: Some(thickness),
        }),
        boundary,
        region,
        area: PI * (radius * radius - inner * inner),
        perimeter: TAU * (radius + inner),
    }
}

/// Samples of the swept surface: the profile's boundary placed by the
/// fixed-reference law along every leg, and, with `caps`, its region at
/// both ends.
fn samples(legs: &[Leg], reference: Vec3, shape: &Shape, caps: bool) -> Vec<Point3> {
    let place = |leg: &Leg, s: Scalar, p: Point2| {
        let t = leg.tangent(s);
        let x = (reference - t * t.dot(reference)).normalize();
        let y = t.cross(x);
        leg.point(s) + x * p.x + y * p.y
    };
    let mut out = Vec::new();
    for leg in legs {
        let n = if matches!(leg, Leg::Bend { .. }) {
            360
        } else {
            12
        };
        for i in 0..n {
            let s = (i as Scalar + 0.5) / n as Scalar;
            out.extend(shape.boundary.iter().map(|p| place(leg, s, *p)));
        }
    }
    if caps {
        let (first, last) = (legs.first().unwrap(), legs.last().unwrap());
        out.extend(shape.region.iter().map(|p| place(first, 0.0, *p)));
        out.extend(shape.region.iter().map(|p| place(last, 1.0, *p)));
    }
    out
}

fn sweep(
    b: &mut GeometryGraphBuilder,
    profile: Profile,
    directrix: NodeId,
    reference: Vec3,
    range: Option<(Scalar, Scalar)>,
) -> NodeId {
    let profile = b.push(GeometryNode::Profile(profile)).unwrap();
    b.push(GeometryNode::SolidOperation(
        SolidOperation::FixedReferenceSweep {
            profile,
            directrix,
            reference_direction: reference,
            parameter_range: range,
        },
    ))
    .unwrap()
}

/// The volume is the profile's area times the centreline length, up to
/// the budget over the surface.
fn assert_volume(name: &str, mesh: &TriMesh, shape: &Shape, length: Scalar, budget: Scalar) {
    let got = axiolid_measure::volume_properties(mesh, Tolerance::new(1e-9, 1e-9).unwrap())
        .unwrap_or_else(|e| panic!("{name}: not a closed two-manifold: {e:?}"))
        .signed_volume;
    let want = shape.area * length;
    let surface = shape.perimeter * length + 2.0 * shape.area;
    assert!(
        got > 0.0 && (got - want).abs() <= surface * budget,
        "{name}: volume {got} against {want} (surface {surface}, budget {budget:e})"
    );
}

/// Sweep `legs` (as a composite, or a loop) at the three budgets.
fn check_chain(name: &str, legs: &[Leg], reference: Vec3, shape: &Shape, closed: bool) {
    let exact = samples(legs, reference, shape, !closed);
    for (options, budget) in budgets() {
        let (mesh, report) = compile(
            |b| {
                let path = composite(b, legs);
                sweep(b, shape.profile.clone(), path, reference, None)
            },
            &options,
        );
        assert_within(name, &mesh, &exact, budget);
        assert_watertight(name, &mesh);
        assert_volume(name, &mesh, shape, length_of(legs), budget);
        assert!(report.meets_requested(), "{name}: {report:?}");
    }
}

fn line_arc_line() -> Vec<Leg> {
    Turtle::new(Point3::ZERO, Vec3::X)
        .straight(0.6)
        .bend(Vec3::Z, 0.4, 1.5)
        .straight(0.5)
        .legs
}

#[test]
fn a_rectangle_swept_along_a_tangent_chain_stays_within_the_budget() {
    check_chain(
        "rectangle along line-arc-line",
        &line_arc_line(),
        Vec3::Z,
        &rect(0.04, 0.06),
        false,
    );
}

#[test]
fn a_profile_with_arcs_and_a_hole_stays_within_the_budget() {
    check_chain(
        "annulus along line-arc-line",
        &line_arc_line(),
        Vec3::Z,
        &annulus(0.06, 0.02),
        false,
    );
}

#[test]
fn a_curved_beam_with_an_in_plane_reference_stays_within_the_budget() {
    // One arc of 2 rad about Z; the reference is the arc's mid radial
    // direction, so the profile's +y maps to -Z along all of it.
    let legs = Turtle::new(Point3::new(2.0, 0.0, 0.0), Vec3::Y)
        .bend(Vec3::Z, 2.0, 2.0)
        .legs;
    let reference = Vec3::new(1.0f64.cos(), 1.0f64.sin(), 0.0);
    for s in [0.0, 0.3, 1.0] {
        let t = legs[0].tangent(s);
        let x = (reference - t * t.dot(reference)).normalize();
        assert!((t.cross(x).dot(Vec3::Z).abs() - 1.0).abs() < 1e-12);
    }
    check_chain(
        "curved beam, in-plane reference",
        &legs,
        reference,
        &rect(0.1, 0.08),
        false,
    );
}

#[test]
fn an_in_plane_reference_along_a_chain_is_re_seeded_at_each_piece() {
    // A quarter turn between two legs, the reference in its plane and a
    // clear 45 degrees off the tangent at both ends: the frame law restarts
    // on each piece, so the last leg's profile is not the first one's.
    let legs = Turtle::new(Point3::ZERO, Vec3::X)
        .straight(0.6)
        .bend(Vec3::Z, 0.4, FRAC_PI_2)
        .straight(0.5)
        .legs;
    check_chain(
        "in-plane reference along line-arc-line",
        &legs,
        Vec3::new(-1.0, 1.0, 0.0),
        &rect(0.04, 0.06),
        false,
    );
}

#[test]
fn an_arch_in_a_vertical_plane_stays_within_the_budget() {
    let legs = Turtle::new(Point3::new(0.5, 0.0, -0.6), Vec3::Z)
        .straight(0.6)
        .bend(-Vec3::Y, 0.5, PI)
        .straight(0.6)
        .legs;
    check_chain("arch", &legs, Vec3::Y, &rect(0.05, 0.07), false);
}

#[test]
fn a_stadium_loop_is_swept_without_caps() {
    let legs = Turtle::new(Point3::new(0.0, -0.5, 0.0), Vec3::X)
        .straight(1.0)
        .bend(Vec3::Z, 0.5, PI)
        .straight(1.0)
        .bend(Vec3::Z, 0.5, PI)
        .legs;
    check_chain("stadium", &legs, Vec3::Z, &rect(0.05, 0.05), true);
}

#[test]
fn a_composite_trimmed_by_arc_length_is_swept_over_the_trim() {
    let full = Turtle::new(Point3::ZERO, Vec3::X)
        .straight(2.0)
        .bend(Vec3::Z, 0.5, 1.5)
        .straight(1.5)
        .legs;
    let total = length_of(&full);
    // Drop a unit off the first leg and half a unit off the last.
    let trimmed: Vec<Leg> = vec![
        match full[0] {
            Leg::Straight { from, dir, .. } => Leg::Straight {
                from: from + dir,
                dir,
                length: 1.0,
            },
            bend => bend,
        },
        full[1],
        match full[2] {
            Leg::Straight { from, dir, .. } => Leg::Straight {
                from,
                dir,
                length: 1.0,
            },
            bend => bend,
        },
    ];
    let shape = rect(0.04, 0.05);
    let exact = samples(&trimmed, Vec3::Z, &shape, true);
    for (options, budget) in budgets() {
        let (mesh, _) = compile(
            |b| {
                let path = composite(b, &full);
                sweep(
                    b,
                    shape.profile.clone(),
                    path,
                    Vec3::Z,
                    Some((1.0, total - 0.5)),
                )
            },
            &options,
        );
        assert_within("trimmed composite", &mesh, &exact, budget);
        assert_watertight("trimmed composite", &mesh);
    }
}

// --- what is refused --------------------------------------------------------

fn polyline(b: &mut GeometryGraphBuilder, points: &[[Scalar; 3]], closed: bool) -> NodeId {
    b.push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
        points: points
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect(),
        closed,
    })))
    .unwrap()
}

fn circle(b: &mut GeometryGraphBuilder, radius: Scalar) -> NodeId {
    b.push(GeometryNode::Curve3(Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    })))
    .unwrap()
}

#[test]
fn a_reference_parallel_to_a_segment_is_refused() {
    refusal(
        |b| {
            let path = polyline(b, &[[0.0; 3], [1.0, 0.0, 0.0]], false);
            sweep(b, rect(0.1, 0.1).profile, path, Vec3::X, None)
        },
        "tangent",
    );
    // Against the direction of travel too: the reference is a line.
    refusal(
        |b| {
            let path = polyline(b, &[[0.0; 3], [1.0, 0.0, 0.0]], false);
            sweep(b, rect(0.1, 0.1).profile, path, -Vec3::X, None)
        },
        "tangent",
    );
}

#[test]
fn an_in_plane_arc_spanning_more_than_a_half_turn_is_refused() {
    // The mid radial direction is parallel to the tangent a quarter turn
    // either side of the middle.
    let build = |span: Scalar| {
        move |b: &mut GeometryGraphBuilder| {
            let path = circle(b, 1.0);
            let mid = 0.5 * span;
            let reference = Vec3::new(mid.cos(), mid.sin(), 0.0);
            sweep(
                b,
                rect(0.1, 0.1).profile,
                path,
                reference,
                Some((0.0, span)),
            )
        }
    };
    refusal(build(3.5), "tangent");
    let (mesh, _) = compile(build(2.5), &options(1e-3));
    assert_watertight("in-plane arc", &mesh);
}

#[test]
fn corners_are_refused_by_name() {
    refusal(
        |b| {
            let path = polyline(b, &[[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]], false);
            sweep(b, rect(0.1, 0.1).profile, path, Vec3::Z, None)
        },
        "tangent continuous",
    );
    // A bend of 0.01 rad is a corner for a profile of this reach at this
    // budget; one of 0.002 rad is tangent continuous.
    let bent = |angle: Scalar| {
        move |b: &mut GeometryGraphBuilder| {
            let path = polyline(
                b,
                &[
                    [0.0; 3],
                    [1.0, 0.0, 0.0],
                    [1.0 + angle.cos(), angle.sin(), 0.0],
                ],
                false,
            );
            sweep(b, rect(0.1, 0.1).profile, path, Vec3::Z, None)
        }
    };
    refusal(bent(0.01), "tangent continuous");
    let (mesh, _) = compile(bent(0.002), &options(1e-3));
    assert!(!mesh.indices.is_empty());
    // A closed polyline is a loop of corners.
    refusal(
        |b| {
            let path = polyline(
                b,
                &[[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                true,
            );
            sweep(b, rect(0.1, 0.1).profile, path, Vec3::Z, None)
        },
        "tangent continuous",
    );
    // A segment ending square to an arc's start.
    refusal(
        |b| {
            let legs = [
                Leg::Straight {
                    from: Point3::ZERO,
                    dir: Vec3::X,
                    length: 1.0,
                },
                Leg::Bend {
                    centre: Point3::ZERO,
                    axis: Vec3::Z,
                    from: Point3::new(1.0, 0.0, 0.0),
                    angle: 1.0,
                },
            ];
            let path = composite(b, &legs);
            sweep(b, rect(0.1, 0.1).profile, path, Vec3::Z, None)
        },
        "tangent continuous",
    );
}

#[test]
fn a_gap_between_pieces_is_refused() {
    refusal(
        |b| {
            let legs = [
                Leg::Straight {
                    from: Point3::ZERO,
                    dir: Vec3::X,
                    length: 1.0,
                },
                Leg::Straight {
                    from: Point3::new(1.0005, 0.0, 0.0),
                    dir: Vec3::X,
                    length: 1.0,
                },
            ];
            let path = composite(b, &legs);
            sweep(b, rect(0.1, 0.1).profile, path, Vec3::Z, None)
        },
        "gap",
    );
}

#[test]
fn a_profile_reaching_a_bends_axis_is_refused() {
    let radius = 0.2;
    let budget = 1e-3;
    let sweep_with = move |y: Scalar| {
        move |b: &mut GeometryGraphBuilder| {
            let path = circle(b, radius);
            sweep(b, rect(0.04, y).profile, path, Vec3::Z, Some((0.0, 1.5)))
        }
    };
    // Reaching the axis exactly, and within the half-budget slack of it.
    refusal(sweep_with(2.0 * radius), "fold through the bend's axis");
    refusal(
        sweep_with(2.0 * (radius - 0.25 * budget)),
        "fold through the bend's axis",
    );
    // Clear of both.
    let (mesh, _) = compile(sweep_with(2.0 * (radius - budget)), &options(budget));
    assert_watertight("fat profile on a tight bend", &mesh);
}

// --- the deviation report ---------------------------------------------------

#[test]
fn a_chain_is_proven_with_its_wall_term_and_one_segment_without() {
    let budget = 1e-3;
    let (_, report) = compile(
        |b| {
            let path = composite(b, &line_arc_line());
            sweep(b, rect(0.04, 0.06).profile, path, Vec3::Z, None)
        },
        &options(budget),
    );
    assert!(report.meets_requested(), "{report:?}");
    let [contribution] = report.contributions.as_slice() else {
        panic!("{report:?}")
    };
    assert_eq!(contribution.path, DeviationPath::FixedReferenceSweep);
    assert_eq!(contribution.detail, "segments and arcs");
    assert_eq!(contribution.bound, DeviationBound::Proven(0.5 * budget));

    let (_, report) = compile(
        |b| {
            let path = polyline(b, &[[0.0; 3], [1.0, 0.0, 0.0]], false);
            sweep(b, rect(0.04, 0.06).profile, path, Vec3::Z, None)
        },
        &options(budget),
    );
    assert_eq!(report.bound, Some(0.0), "{report:?}");
    assert!(report.meets_requested());
}

#[test]
fn an_arc_oblique_to_the_reference_is_meshed_and_reported_unbounded() {
    let (mesh, report) = compile(
        |b| {
            let path = circle(b, 0.5);
            sweep(
                b,
                rect(0.04, 0.04).profile,
                path,
                Vec3::new(1.0, 0.0, 1.0),
                Some((0.0, 1.5)),
            )
        },
        &options(1e-3),
    );
    assert!(!mesh.indices.is_empty());
    assert_eq!(report.bound, None);
    assert!(!report.meets_requested());
    let unbounded: Vec<_> = report.unbounded().collect();
    assert!(matches!(
        unbounded.as_slice(),
        [c] if matches!(c.bound, DeviationBound::Unbounded(reason) if reason.contains("twists"))
    ));
}

// --- a disk swept along a lone arc ------------------------------------------

#[test]
fn a_disk_wider_than_a_lone_arcs_bend_is_refused() {
    let disk = |radius: Scalar| {
        move |b: &mut GeometryGraphBuilder| {
            let directrix = circle(b, 0.1);
            b.push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
                directrix,
                radius,
                inner_radius: None,
                parameter_range: Some((0.0, 1.0)),
                fillet_radius: None,
            }))
            .unwrap()
        }
    };
    refusal(disk(0.12), "not below the bend radius");
    refusal(disk(0.1), "equals the bend radius");
    let (mesh, _) = compile(disk(0.05), &options(1e-3));
    assert_watertight("narrow disk on a lone arc", &mesh);
}
