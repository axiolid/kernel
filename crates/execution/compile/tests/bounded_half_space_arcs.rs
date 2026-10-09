//! A bounded half-space whose boundary has circular-arc segments (#277).
//!
//! The boundary is a profile node, a contour of lines and exact `Circle2`
//! arcs, and the half-space is the side of a plane within the prism of
//! that region along the plane normal, so every arc bounds it by a right
//! circular cylinder. A wall `[-L/2, L/2] x [-T/2, T/2] x [0, H]` is
//! clipped above `z = ZC` by such a half-space under a rigid placement.
//!
//! The boundary of the consumer's wall (openbimrs/ifc#398) is six segments,
//! two of them arcs over one circle of radius 1.2: here a rectangle with a
//! round bite taken out of its side, its arc split at the circle's
//! rightmost point. The removed region of the wall's plan is the strip
//! `|y| <= T/2` right of the arc (or, for the bulging boundary, left of
//! it), whose area is closed form, so the clipped volume is too:
//! `L T H - (H - ZC) area`.
//!
//! The exact compiler is checked against that volume. The mesh compiler's
//! deviation report must be certified against the exact result and hold
//! under dense sampling of the exact faces, which also checks that both
//! compilers frame the boundary alike: a footprint mirrored or turned in
//! one of them would put the cylinder wall where the other has none. A
//! polyline boundary is untouched (`exact_half_space_clip.rs`); invalid
//! profile boundaries are refused by both compilers.

use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{
    BooleanOperator, Frame2, Interval, Plane3, Point2, Point3, Scalar, Tolerance, Transform3, Vec2,
    Vec3,
};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_measure::{exact_properties, FaceDomain};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationPath, ReferenceExactCompiler, ReferenceMeshCompiler,
};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

/// Wall length, thickness and height, and the clip plane's height.
const L: f64 = 6.0;
const T: f64 = 0.3;
const H: f64 = 3.0;
const ZC: f64 = 2.0;
/// The consumer's arc radius.
const R: f64 = 1.2;

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
        .with_chord_error(budget)
        .expect("a positive budget")
}

fn line(a: Point2, b: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: b - a,
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    }
}

fn arc(centre: Point2, radius: f64, (from, to): (f64, f64)) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius,
        }),
        domain: Interval::new(from.min(to), from.max(to)),
        same_sense: from < to,
    }
}

fn contour(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

/// Six segments, two of them arcs over the circle of radius `R` about
/// `(0, cy)`: the rectangle `[0, 4] x [cy - 1, cy + 2]` less that disk,
/// the arc split where it crosses `y = cy`, run clockwise.
fn bitten(cy: f64) -> Profile {
    let c = Point2::new(0.0, cy);
    let foot = (R * R - 1.0).sqrt();
    let below = -(1.0 / R).asin();
    contour(vec![
        line(Point2::new(foot, cy - 1.0), Point2::new(4.0, cy - 1.0)),
        line(Point2::new(4.0, cy - 1.0), Point2::new(4.0, cy + 2.0)),
        line(Point2::new(4.0, cy + 2.0), Point2::new(0.0, cy + 2.0)),
        line(Point2::new(0.0, cy + 2.0), Point2::new(0.0, cy + R)),
        arc(c, R, (FRAC_PI_2, 0.0)),
        arc(c, R, (0.0, below)),
    ])
}

/// Seven segments: the rectangle `[-4, 0] x [cy - 2, cy + 2]` with the
/// right half of the same disk added, the arc split at `y = cy` and run
/// counter-clockwise.
fn bulging(cy: f64) -> Profile {
    let c = Point2::new(0.0, cy);
    contour(vec![
        line(Point2::new(-4.0, cy - 2.0), Point2::new(0.0, cy - 2.0)),
        line(Point2::new(0.0, cy - 2.0), Point2::new(0.0, cy - R)),
        arc(c, R, (-FRAC_PI_2, 0.0)),
        arc(c, R, (0.0, FRAC_PI_2)),
        line(Point2::new(0.0, cy + R), Point2::new(0.0, cy + 2.0)),
        line(Point2::new(0.0, cy + 2.0), Point2::new(-4.0, cy + 2.0)),
        line(Point2::new(-4.0, cy + 2.0), Point2::new(-4.0, cy - 2.0)),
    ])
}

/// `integral of sqrt(R^2 - u^2) du` over the strip `|y| <= T/2`, `u = y - cy`.
fn disk_strip(cy: f64) -> f64 {
    let f = |u: f64| 0.5 * u * (R * R - u * u).sqrt() + 0.5 * R * R * (u / R).asin();
    let t = T / 2.0;
    f(t - cy) - f(-t - cy)
}

/// The wall's plan area the boundary covers: right of the arc to the wall's
/// end at `x = L/2` for the bite, from `x = -L/2` to the arc for the bulge.
fn covered(profile: &str, cy: f64) -> f64 {
    match profile {
        "bitten" => L / 2.0 * T - disk_strip(cy),
        "bulging" => L / 2.0 * T + disk_strip(cy),
        _ => unreachable!(),
    }
}

fn boundary(profile: &str, cy: f64) -> Profile {
    match profile {
        "bitten" => bitten(cy),
        "bulging" => bulging(cy),
        _ => unreachable!(),
    }
}

fn expected_volume(profile: &str, cy: f64) -> f64 {
    L * T * H - (H - ZC) * covered(profile, cy)
}

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

/// How the clip is stated: the plane's normal up and its side kept, or the
/// normal down and the other side kept. The boundary frame is a quarter
/// turn about `z`, moved off the plane; the profile is authored in it, so
/// its points are the plan's turned back. Down, the in-plane `y` is
/// `normal x x`, the plan's `-y`, so `cy` is authored negated.
#[derive(Debug, Clone, Copy)]
struct Stated {
    up: bool,
}

/// The wall minus the half-space above `z = ZC` bounded by `profile` with
/// its circle about `(0, cy)` in the wall's plan, all under `placement`.
fn clipped(
    g: &mut GeometryGraphBuilder,
    boundary: Profile,
    how: Stated,
    placement: Transform3,
) -> NodeId {
    let mut push = |node| g.push(node).expect("a valid node");
    let profile = push(GeometryNode::Profile(rect(L, T)));
    let wall = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: H,
    }));
    let wall = push(GeometryNode::Instance(Instance {
        source: wall,
        transform: placement,
    }));
    let normal = if how.up { Vec3::Z } else { -Vec3::Z };
    let half_space = push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, ZC),
            normal,
        },
        agreement: how.up,
    }));
    let boundary = push(GeometryNode::Profile(boundary));
    let bounded = push(GeometryNode::SolidOperation(
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary,
            placement: Transform3::from_translation(Vec3::new(0.0, 0.0, 5.0)),
        },
    ));
    let bounded = push(GeometryNode::Instance(Instance {
        source: bounded,
        transform: placement,
    }));
    push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: wall,
        right: bounded,
        operator: BooleanOperator::Difference,
    }))
}

/// Build the clipped wall for `profile` about `(0, cy)` in the plan.
fn case(profile: &str, cy: f64, how: Stated, placement: Transform3) -> (GeometryGraph, NodeId) {
    let mut g = GeometryGraphBuilder::new();
    // Down, the profile's `y` is the plan's `-y`.
    let authored = if how.up { cy } else { 0.0 - cy };
    let root = clipped(&mut g, boundary(profile, authored), how, placement);
    (g.finish(vec![root]).expect("a valid graph"), root)
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

fn mesh_volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

#[test]
fn a_wall_clipped_by_an_arc_bounded_half_space_has_the_closed_form_volume() {
    for profile in ["bitten", "bulging"] {
        for cy in [0.0, 0.3] {
            for up in [true, false] {
                for placement in [building(), general()] {
                    let how = Stated { up };
                    let (graph, root) = case(profile, cy, how, placement);
                    let exact = ReferenceExactCompiler::new()
                        .compile_exact(&graph, root, &options(1e-3))
                        .unwrap_or_else(|e| panic!("{profile} {cy} {how:?}: {e:?}"));
                    let health = geometric_audit(&exact, Tolerance::METRE);
                    assert!(health.is_consistent(), "{:?}", health.defects());
                    let topology = axiolid_topology::audit_brep(exact.topology());
                    assert!(topology.is_closed_manifold(), "{topology:?}");
                    let expected = expected_volume(profile, cy);
                    let measured = volume(&exact);
                    assert!(
                        (measured - expected).abs() <= 1e-9 * expected,
                        "{profile} {cy} {how:?}: volume {measured}, expected {expected}"
                    );
                    // The arc is a cylinder wall of the result, not chords.
                    assert!(
                        exact
                            .surfaces()
                            .iter()
                            .any(|s| matches!(s, axiolid_surface::Surface::Cylinder(_))),
                        "{profile}: the arc must cut a cylindrical face"
                    );
                }
            }
        }
    }
}

#[test]
fn the_mesh_compiler_clips_within_its_chord_budget() {
    for profile in ["bitten", "bulging"] {
        for up in [true, false] {
            let budget = 1e-3;
            let (graph, root) = case(profile, 0.3, Stated { up }, building());
            let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
                .compile_mesh(&graph, root, &options(budget))
                .expect("the mesh compiler clips too");
            let expected = expected_volume(profile, 0.3);
            // Each chord of the arc moves the cut by at most the budget over
            // the strip's width and the clipped height.
            let slack = budget * T * (H - ZC) + 1e-9;
            let measured = mesh_volume(&mesh);
            assert!(
                (measured - expected).abs() <= slack,
                "{profile} up {up}: mesh volume {measured}, exact {expected}"
            );
        }
    }
}

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

/// A lattice of each face's parameter box where the face certainly
/// contains it, and `n` points along every pcurve span of every loop.
fn exact_samples(exact: &ExactBRep, n: usize) -> Vec<Point3> {
    let topology = exact.topology();
    let mut out = Vec::new();
    for (index, face) in topology.faces().iter().enumerate() {
        let surface = &exact.surfaces()[face.surface.expect("a surface").index()];
        let domain = FaceDomain::new(exact, topology.face_id_at(index).unwrap(), Tolerance::METRE)
            .expect("a domain")
            .expect("a monotone split");
        let (lo, hi) = domain.bounds();
        for i in 0..n {
            for j in 0..n {
                let at = Point2::new(
                    lo.x + (hi.x - lo.x) * (i as Scalar + 0.5) / n as Scalar,
                    lo.y + (hi.y - lo.y) * (j as Scalar + 0.5) / n as Scalar,
                );
                if domain.contains(at).expect("decidable") == Some(true) {
                    out.push(axiolid_reference::surface::evaluate(surface, at.x, at.y).unwrap());
                }
            }
        }
        for bound in &face.bounds {
            let wire = &topology.loops()[bound.loop_id.index()];
            for (k, edge_use) in wire.edges.iter().enumerate() {
                let curve = &exact.curves2()[edge_use.pcurve.expect("a pcurve").index()];
                let span = exact.pcurve_interval(bound.loop_id, k).expect("a span");
                for s in 0..=n {
                    let t = span.start + (span.end - span.start) * s as Scalar / n as Scalar;
                    let uv = axiolid_reference::evaluate2(curve, t).unwrap();
                    out.push(axiolid_reference::surface::evaluate(surface, uv.x, uv.y).unwrap());
                }
            }
        }
    }
    out
}

#[test]
fn the_mesh_deviation_is_certified_and_holds_under_dense_sampling() {
    for profile in ["bitten", "bulging"] {
        for up in [true, false] {
            for budget in [1e-2, 1e-3] {
                let options = options(budget);
                let (graph, root) = case(profile, 0.3, Stated { up }, general());
                let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
                    .compile_mesh_with_deviation(&graph, root, &options)
                    .expect("the clipped wall compiles");
                assert_eq!(report.contributions.len(), 1, "{report:?}");
                let contribution = report.contributions[0];
                assert_eq!(contribution.path, DeviationPath::Boolean);
                let DeviationBound::Certified(bound) = contribution.bound else {
                    panic!("{profile} up {up}: expected a certified bound, got {report:?}");
                };
                assert!(
                    bound <= 2.0 * budget,
                    "{profile} up {up}: bound {bound} for a budget of {budget}"
                );
                let exact = ReferenceExactCompiler::new()
                    .compile_exact(&graph, root, &options)
                    .expect("the exact result");
                let samples = exact_samples(&exact, 200);
                let grid = TriangleGrid::new(&outcome.mesh, 2.0 * bound);
                let worst = samples
                    .iter()
                    .map(|p| grid.distance(*p))
                    .fold(0.0, Scalar::max);
                eprintln!(
                    "{profile} up {up} budget {budget:e}: reported {bound:e}, measured \
                     {worst:e}, {} samples",
                    samples.len()
                );
                assert!(
                    worst <= bound,
                    "{profile} up {up}: an exact point lies {worst:e} from the mesh, \
                     above the reported {bound:e}"
                );
            }
        }
    }
}

/// Both compilers' answers for the wall clipped by `boundary`.
fn both(boundary: Profile) -> (Result<ExactBRep, GeomError>, Result<TriMesh, GeomError>) {
    let mut g = GeometryGraphBuilder::new();
    let root = clipped(&mut g, boundary, Stated { up: true }, building());
    let graph = g.finish(vec![root]).expect("a valid graph");
    let options = options(1e-3);
    (
        ReferenceExactCompiler::new().compile_exact(&graph, root, &options),
        ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(&graph, root, &options),
    )
}

fn refused_by_both(name: &str, boundary: Profile, part: &str) {
    let (exact, mesh) = both(boundary);
    for (path, error) in [("exact", exact.err()), ("mesh", mesh.err())] {
        let error = error.unwrap_or_else(|| panic!("{name}: the {path} compiler must refuse"));
        assert!(
            matches!(error, GeomError::InvalidInput(_)),
            "{name}: the {path} refusal must blame the input, got {error:?}"
        );
        let text = format!("{error}");
        assert!(
            text.contains(part),
            "{name}: the {path} refusal `{text}` must name `{part}`"
        );
    }
}

#[test]
fn invalid_arc_boundaries_are_refused_by_name() {
    // Open: the bite's bottom line dropped. The mesh flattener alone would
    // bridge the gap with a chord.
    let Profile::Contour(mut open) = bitten(0.0) else {
        unreachable!()
    };
    open.outer.segments.remove(0);
    refused_by_both("open", Profile::Contour(open), "does not close");

    // A gap between two segments, not at the closing joint.
    let Profile::Contour(mut gapped) = bitten(0.0) else {
        unreachable!()
    };
    gapped.outer.segments.remove(2);
    refused_by_both("gapped", Profile::Contour(gapped), "gap");

    // A zero radius.
    let c = Point2::new(0.0, 0.0);
    let point = contour(vec![
        line(c, Point2::new(2.0, 0.0)),
        line(Point2::new(2.0, 0.0), c),
        arc(c, 0.0, (0.0, 1.0)),
    ]);
    refused_by_both("zero radius", point, "radius 0");

    // Self-crossing: the half circle over `[0, 2]` crosses the slanted line
    // from its start.
    let crossing = contour(vec![
        line(Point2::new(0.0, 0.0), Point2::new(2.0, 1.0)),
        line(Point2::new(2.0, 1.0), Point2::new(2.0, 0.0)),
        arc(Point2::new(1.0, 0.0), 1.0, (0.0, FRAC_PI_2 * 2.0)),
    ]);
    refused_by_both("self-crossing", crossing, "crosses or touches itself");

    // A line crossing an arc: the upper half circle of radius `sqrt 2`
    // over `(-1, 0)` and `(1, 0)`, closed through a corner above it.
    let over = contour(vec![
        arc(
            Point2::new(0.0, -1.0),
            2.0_f64.sqrt(),
            (FRAC_PI_2 / 2.0, FRAC_PI_2 * 1.5),
        ),
        line(Point2::new(-1.0, 0.0), Point2::new(0.5, 0.5)),
        line(Point2::new(0.5, 0.5), Point2::new(1.0, 0.0)),
    ]);
    refused_by_both("line crossing an arc", over, "crosses or touches itself");

    // Two arcs touching: the upper unit half circle, and the arc over
    // `(-1, -1)` and `(1, -1)` about `(0, -0.5)`, which runs through both
    // of its ends.
    let low = Point2::new(0.0, -0.5);
    let r = 1.25_f64.sqrt();
    let from = (-0.5_f64).atan2(-1.0) + 2.0 * PI;
    let to = (-0.5_f64).atan2(1.0);
    let touching = contour(vec![
        arc(Point2::new(0.0, 0.0), 1.0, (0.0, PI)),
        line(Point2::new(-1.0, 0.0), Point2::new(-1.0, -1.0)),
        arc(low, r, (from, to)),
        line(Point2::new(1.0, -1.0), Point2::new(1.0, 0.0)),
    ]);
    refused_by_both("arcs touching", touching, "crosses or touches itself");

    // A line passing half a tolerance over the unit circle, against the
    // inside of an arc from -10 to 180 degrees, split at 85: it touches
    // within the tolerance where the exact tangency has no root.
    let lift = 1.0 + 5e-7;
    let low = (-PI / 18.0).sin_cos();
    let corner = Point2::new(low.1, low.0);
    let grazed = contour(vec![
        arc(Point2::new(0.0, 0.0), 1.0, (-PI / 18.0, PI)),
        line(Point2::new(-1.0, 0.0), Point2::new(-1.0, lift)),
        line(Point2::new(-1.0, lift), Point2::new(1.5, lift)),
        line(Point2::new(1.5, lift), Point2::new(1.5, corner.y)),
        line(Point2::new(1.5, corner.y), corner),
    ]);
    refused_by_both(
        "line grazing an arc",
        grazed,
        "edge 1 of ring 0 and edge 3 of ring 0",
    );

    // A bow tie closed by a half circle: two of its lines cross.
    let bow = contour(vec![
        line(Point2::new(0.0, 0.0), Point2::new(2.0, 2.0)),
        line(Point2::new(2.0, 2.0), Point2::new(2.0, 0.0)),
        line(Point2::new(2.0, 0.0), Point2::new(0.0, 2.0)),
        arc(Point2::new(0.0, 1.0), 1.0, (FRAC_PI_2, 3.0 * FRAC_PI_2)),
    ]);
    refused_by_both("bow tie", bow, "edge 0 of ring 0 and edge 2 of ring 0");

    // An arc folding back over the one before it on the same circle: from
    // 45 to 90 degrees, then back from 90 to 0.
    let c = Point2::new(0.0, 0.0);
    let diagonal = Point2::new(FRAC_PI_4.cos(), FRAC_PI_4.sin());
    let folded = contour(vec![
        arc(c, 1.0, (FRAC_PI_4, FRAC_PI_2)),
        arc(c, 1.0, (FRAC_PI_2, 0.0)),
        line(Point2::new(1.0, 0.0), Point2::new(1.5, 0.0)),
        line(Point2::new(1.5, 0.0), Point2::new(1.5, 1.5)),
        line(Point2::new(1.5, 1.5), diagonal),
    ]);
    refused_by_both(
        "arc folding back",
        folded,
        "edge 0 of ring 0 and edge 1 of ring 0",
    );

    // A segment of no length.
    let Profile::Contour(mut stalled) = bitten(0.0) else {
        unreachable!()
    };
    let at = Point2::new(4.0, -1.0);
    stalled.outer.segments.insert(1, line(at, at));
    refused_by_both(
        "segment of no length",
        Profile::Contour(stalled),
        "no length",
    );

    // A hole touching the outer ring.
    let Profile::Contour(mut holed) = bitten(0.0) else {
        unreachable!()
    };
    let (a, b, d) = (
        Point2::new(3.0, -0.5),
        Point2::new(4.0, 0.0),
        Point2::new(3.0, 0.5),
    );
    holed
        .holes
        .push(Contour::new(vec![line(a, b), line(b, d), line(d, a)]));
    refused_by_both("touching hole", Profile::Contour(holed), "ring 1");
}

#[test]
fn tangent_joints_and_half_turns_are_accepted() {
    // A stadium: two lines tangent to two half circles of radius `R` about
    // `(-1, 0)` and `(1, 0)`. Each half turn is lowered as two sub-arcs.
    let stadium = contour(vec![
        line(Point2::new(-1.0, -R), Point2::new(1.0, -R)),
        arc(Point2::new(1.0, 0.0), R, (-FRAC_PI_2, FRAC_PI_2)),
        line(Point2::new(1.0, R), Point2::new(-1.0, R)),
        arc(Point2::new(-1.0, 0.0), R, (FRAC_PI_2, 3.0 * FRAC_PI_2)),
    ]);
    // Over the strip, `|x| <= 1 + sqrt(R^2 - y^2)`.
    let covered = 2.0 * T + 2.0 * disk_strip(0.0);
    let expected = L * T * H - (H - ZC) * covered;
    let (exact, mesh) = both(stadium);
    let exact = exact.expect("the exact compiler clips by a stadium");
    let measured = volume(&exact);
    assert!(
        (measured - expected).abs() <= 1e-9 * expected,
        "volume {measured}, expected {expected}"
    );
    let mesh = mesh.expect("the mesh compiler clips by a stadium");
    let slack = 2.0 * 1e-3 * T * (H - ZC);
    assert!((mesh_volume(&mesh) - expected).abs() <= slack);
}

#[test]
fn a_profile_boundary_alone_is_still_an_unbounded_half_space() {
    let mut g = GeometryGraphBuilder::new();
    let mut push = |node| g.push(node).expect("a valid node");
    let half_space = push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: Point3::ZERO,
            normal: Vec3::Z,
        },
        agreement: true,
    }));
    let boundary = push(GeometryNode::Profile(bitten(0.0)));
    let bounded = push(GeometryNode::SolidOperation(
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary,
            placement: Transform3::IDENTITY,
        },
    ));
    let graph = g.finish(vec![bounded]).expect("a valid graph");
    let error = ReferenceExactCompiler::new()
        .compile_exact(&graph, bounded, &options(1e-3))
        .expect_err("refused");
    assert!(
        matches!(&error, GeomError::UnsupportedInput { input, .. } if input.contains("half-space")),
        "{error:?}"
    );
    // The mesh compiler bounds it by its slab, as a polygon's.
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, bounded, &options(1e-3))
        .expect("the slab meshes");
    assert!(mesh_volume(&mesh) > 0.0);
}
