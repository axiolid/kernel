//! Exact half-space clipping of placed solids: walls cut by roofs (#234).
//!
//! The wall is `[-L/2, L/2] x [-T/2, T/2] x [0, H]` in its own frame, under a
//! rigid placement. A roof plane `z = z0 + s x` in that frame clips it; the
//! half-space is given either in world coordinates or in the wall's frame
//! under the wall's placement. Each case checks, against closed forms, the
//! exact volume and the certified distance from the result's boundary to a
//! probe box that would cross the unclipped wall, plus a clean audit and a
//! closed two-manifold topology. The unbounded cases are also compared with
//! the mesh compiler, whose semantics (which side `agreement` keeps, how a
//! bounded half-space frames its boundary) the exact path keeps.

use std::f64::consts::{FRAC_PI_2, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Plane3, Point2, Point3, Tolerance, Transform3, Vec3};
use axiolid_curve::{Circle2, Curve2, Polyline2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::{boundary_distance, exact_properties};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{ReferenceExactCompiler, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// The linear tolerance every test here runs at.
const EPS: f64 = 1e-6;

/// Wall length, thickness and height.
const L: f64 = 6.0;
const T: f64 = 0.3;
const H: f64 = 3.0;

struct Graph {
    builder: GeometryGraphBuilder,
}

impl Graph {
    fn new() -> Self {
        Self {
            builder: GeometryGraphBuilder::new(),
        }
    }

    fn push(&mut self, node: GeometryNode) -> NodeId {
        self.builder.push(node).expect("a valid node")
    }

    fn extrusion(&mut self, profile: Profile, depth: f64) -> NodeId {
        let profile = self.push(GeometryNode::Profile(profile));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn instance(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }

    fn minus(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.boolean(left, right, BooleanOperator::Difference)
    }

    fn wall(&mut self, placement: Transform3) -> NodeId {
        let wall = self.extrusion(rect(L, T), H);
        self.instance(wall, placement)
    }

    /// The side of the plane through `origin` with `normal` (in the wall's
    /// frame) that `agreement` selects, in world coordinates when `world`,
    /// else in the wall's frame under an instance of `placement`.
    fn half_space(
        &mut self,
        (origin, normal): (Point3, Vec3),
        agreement: bool,
        placement: Transform3,
        world: bool,
    ) -> NodeId {
        if world {
            self.push(GeometryNode::HalfSpace(HalfSpace {
                boundary: Plane3 {
                    origin: placement.transform_point3(origin),
                    normal: placement.transform_vector3(normal),
                },
                agreement,
            }))
        } else {
            let local = self.push(GeometryNode::HalfSpace(HalfSpace {
                boundary: Plane3 { origin, normal },
                agreement,
            }));
            self.instance(local, placement)
        }
    }

    /// A half-space bounded by `boundary`, authored in `frame`, all in the
    /// wall's frame under an instance of `placement`.
    fn bounded(
        &mut self,
        (origin, normal): (Point3, Vec3),
        agreement: bool,
        boundary: Vec<Point2>,
        frame: Transform3,
        placement: Transform3,
    ) -> NodeId {
        let half_space = self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 { origin, normal },
            agreement,
        }));
        let boundary = self.push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
            points: boundary,
            closed: true,
        })));
        let bounded = self.push(GeometryNode::SolidOperation(
            SolidOperation::BoundedHalfSpace {
                half_space,
                boundary,
                placement: frame,
            },
        ));
        self.instance(bounded, placement)
    }

    /// A box `[min, max]` in the wall's frame, under `placement`.
    fn probe(&mut self, min: Vec3, max: Vec3, placement: Transform3) -> NodeId {
        let size = max - min;
        let block = self.extrusion(rect(size.x, size.y), size.z);
        let centre = Vec3::new((min.x + max.x) / 2.0, (min.y + max.y) / 2.0, min.z);
        self.instance(block, placement * Transform3::from_translation(centre))
    }

    /// A window `w x h` through the wall at `x = cx` from `z = sill`.
    fn window(&mut self, (w, h): (f64, f64), (cx, sill): (f64, f64), p: Transform3) -> NodeId {
        let opening = self.extrusion(rect(w, h), T + 0.2);
        let across = self.instance(
            opening,
            Transform3::from_translation(Vec3::new(cx, T / 2.0 + 0.1, sill + h / 2.0))
                * Transform3::from_rotation_x(FRAC_PI_2),
        );
        self.instance(across, p)
    }

    fn finish(self, roots: Vec<NodeId>) -> GeometryGraph {
        self.builder.finish(roots).expect("a valid graph")
    }
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

/// The roof `z = z0 + s x` in the wall's frame, its normal pointing up.
fn roof_plane(z0: f64, s: f64) -> (Point3, Vec3) {
    (Point3::new(0.0, 0.0, z0), Vec3::new(-s, 0.0, 1.0))
}

fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

fn compile(graph: &GeometryGraph, roots: &[NodeId]) -> Result<Vec<ExactBRep>, GeomError> {
    let mut out = Vec::new();
    ReferenceExactCompiler::new()
        .compile_exact_batch_into(graph, roots, &options(), &mut out)
        .map(|()| out)
}

fn audited(brep: &ExactBRep) {
    let health = geometric_audit(brep, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(brep.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
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

/// Compile `body` and each probe; check the volume, the audit, and each
/// probe's certified distance to the body. Returns the body.
fn check(graph: Graph, body: NodeId, expected: f64, probes: &[(NodeId, f64)]) -> ExactBRep {
    let mut roots = vec![body];
    roots.extend(probes.iter().map(|&(probe, _)| probe));
    let graph = graph.finish(roots.clone());
    let solids = compile(&graph, &roots).expect("exact");
    audited(&solids[0]);
    let measured = volume(&solids[0]);
    assert!(
        (measured - expected).abs() <= 1e-9 * expected,
        "volume {measured}, expected {expected}"
    );
    for (probe, &(_, distance)) in solids[1..].iter().zip(probes) {
        let bounds = boundary_distance(&solids[0], probe, 1e-9, Tolerance::METRE).expect("bounded");
        assert!(
            bounds.lower <= distance + 1e-12 && distance <= bounds.upper + 1e-12,
            "[{}, {}] must contain {distance}",
            bounds.lower,
            bounds.upper
        );
        assert!(bounds.upper - bounds.lower <= 1e-9);
    }
    solids.into_iter().next().expect("the body")
}

/// The mesh compiler's volume for `body`, which must agree with the exact
/// one: the two compilers read `agreement` and the boundary frame alike.
fn mesh_agrees(graph: Graph, body: NodeId, expected: f64) {
    let graph = graph.finish(vec![body]);
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, body, &options())
        .expect("the mesh compiler clips too");
    let measured = mesh_volume(&mesh);
    assert!(
        (measured - expected).abs() <= 1e-6 * expected,
        "mesh volume {measured}, exact {expected}"
    );
}

#[test]
fn one_roof_plane_clips_a_placed_wall() {
    let (z0, s) = (2.4, 0.1);
    for placement in [building(), general()] {
        for world in [true, false] {
            let build = |g: &mut Graph| {
                let wall = g.wall(placement);
                let roof = g.half_space(roof_plane(z0, s), true, placement, world);
                g.minus(wall, roof)
            };
            let mut g = Graph::new();
            let body = build(&mut g);
            // Above the roof over x in [-1, 1], 0.2 over the roof at x = 1,
            // where the roof is highest: the distance is that gap measured
            // across the roof.
            let probe = g.probe(
                Vec3::new(-1.0, -T, z0 + s + 0.2),
                Vec3::new(1.0, T, H + 1.0),
                placement,
            );
            let to_roof = 0.2 / (1.0 + s * s).sqrt();
            check(g, body, T * L * z0, &[(probe, to_roof)]);

            // The mesh compiler reads a half-space operand unplaced.
            if world {
                let mut g = Graph::new();
                let body = build(&mut g);
                mesh_agrees(g, body, T * L * z0);
            }
        }
    }
}

/// The gable `z0 - s |x|` under two roof planes.
fn gable(g: &mut Graph, body: NodeId, (z0, s): (f64, f64), placement: Transform3) -> NodeId {
    let left = g.half_space(roof_plane(z0, s), true, placement, true);
    let body = g.minus(body, left);
    let right = g.half_space(roof_plane(z0, -s), true, placement, false);
    g.minus(body, right)
}

#[test]
fn two_roof_planes_clip_a_gable_wall() {
    let (z0, s) = (2.4, 0.3);
    for placement in [building(), general()] {
        let mut g = Graph::new();
        let wall = g.wall(placement);
        let body = gable(&mut g, wall, (z0, s), placement);
        // Over the ridge: the ridge line is the nearest part of the wall.
        let probe = g.probe(
            Vec3::new(-0.5, -T, z0 + 0.15),
            Vec3::new(0.5, T, H + 1.0),
            placement,
        );
        check(g, body, T * (L * z0 - s * L * L / 4.0), &[(probe, 0.15)]);
    }
}

#[test]
fn a_gable_wall_with_a_window_clips_in_either_order() {
    let (z0, s) = (2.4, 0.3);
    let ((w, h), at) = ((1.0, 1.2), (0.5, 0.6));
    let expected = T * (L * z0 - s * L * L / 4.0) - w * h * T;
    let placement = general();
    for window_first in [true, false] {
        let mut g = Graph::new();
        let wall = g.wall(placement);
        let body = if window_first {
            let window = g.window((w, h), at, placement);
            let cut = g.minus(wall, window);
            gable(&mut g, cut, (z0, s), placement)
        } else {
            let clipped = gable(&mut g, wall, (z0, s), placement);
            let window = g.window((w, h), at, placement);
            g.minus(clipped, window)
        };
        let gap = 0.05;
        let probe = g.probe(
            Vec3::new(at.0 - w / 2.0 + gap, -T, at.1 + gap),
            Vec3::new(at.0 + w / 2.0 - gap, T, at.1 + h - gap),
            placement,
        );
        check(g, body, expected, &[(probe, gap)]);
    }
}

#[test]
fn an_intersection_keeps_the_selected_side() {
    // `agreement` selects the normal side; an intersection keeps it.
    let (z0, s) = (2.4, 0.1);
    let placement = general();
    for (agreement, expected) in [(false, T * L * z0), (true, L * T * H - T * L * z0)] {
        let mut g = Graph::new();
        let wall = g.wall(placement);
        let roof = g.half_space(roof_plane(z0, s), agreement, placement, true);
        let body = g.boolean(wall, roof, BooleanOperator::Intersection);
        check(g, body, expected, &[]);
        let mut g = Graph::new();
        let wall = g.wall(placement);
        let roof = g.half_space(roof_plane(z0, s), agreement, placement, true);
        let body = g.boolean(wall, roof, BooleanOperator::Intersection);
        mesh_agrees(g, body, expected);
    }
}

/// A quarter turn about `z`, moved to `at`.
fn quarter_turn(at: Vec3) -> Transform3 {
    Transform3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z, at)
}

#[test]
fn a_bounded_half_space_cuts_only_part_of_the_wall() {
    // Above z = 2 over x in [1, 4], y in [-1, 1] of the wall's frame: a step
    // down at x = 1 to the wall's end. The boundary is authored in a frame
    // turned a quarter about the plane normal and moved off the plane; the
    // frame's x runs along the wall's y. With the normal down and the other
    // side selected, the in-plane y turns over, and so does the boundary.
    let (zc, x1) = (2.0, 1.0);
    let expected = L * T * H - (L / 2.0 - x1) * T * (H - zc);
    let frame = quarter_turn(Vec3::new(2.0, 0.0, 5.0));
    let up = vec![
        Point2::new(-1.0, -2.0),
        Point2::new(1.0, -2.0),
        Point2::new(1.0, 1.0),
        Point2::new(-1.0, 1.0),
        Point2::new(-1.0, -2.0),
    ];
    // Mirrored, and reversed to stay counter-clockwise: the mesh compiler
    // sweeps the boundary as wound.
    let down = up.iter().rev().map(|p| Point2::new(p.x, -p.y)).collect();
    for placement in [building(), general()] {
        for (normal, agreement, boundary) in [(Vec3::Z, true, &up), (-Vec3::Z, false, &down)] {
            let build = |g: &mut Graph| {
                let wall = g.wall(placement);
                let plane = (Point3::new(0.0, 0.0, zc), normal);
                let tool = g.bounded(plane, agreement, boundary.clone(), frame, placement);
                g.minus(wall, tool)
            };
            let mut g = Graph::new();
            let body = build(&mut g);
            // In the notch, 0.05 from the step and 0.1 above its floor.
            let probe = g.probe(
                Vec3::new(x1 + 0.05, -T, zc + 0.1),
                Vec3::new(L, T, H + 1.0),
                placement,
            );
            check(g, body, expected, &[(probe, 0.05)]);
            let mut g = Graph::new();
            let body = build(&mut g);
            mesh_agrees(g, body, expected);
        }
    }
}

#[test]
fn a_bounded_half_space_wider_than_the_wall_clips_as_its_plane() {
    let (z0, s) = (2.4, 0.1);
    let square = vec![
        Point2::new(-10.0, -10.0),
        Point2::new(10.0, -10.0),
        Point2::new(10.0, 10.0),
        Point2::new(-10.0, 10.0),
    ];
    let frame = Transform3::from_axis_angle(Vec3::new(-s, 0.0, 1.0).normalize(), 0.4);
    let placement = general();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let tool = g.bounded(roof_plane(z0, s), true, square, frame, placement);
    let body = g.minus(wall, tool);
    check(g, body, T * L * z0, &[]);
}

/// The wall clipped by `z = H + offset` from above, compiled.
fn wall_under_a_flat_roof(offset: f64, placement: Transform3) -> ExactBRep {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let plane = (Point3::new(0.0, 0.0, H + offset), Vec3::Z);
    let roof = g.half_space(plane, true, placement, true);
    let body = g.minus(wall, roof);
    let graph = g.finish(vec![body]);
    let solid = compile(&graph, &[body]).expect("exact").remove(0);
    audited(&solid);
    solid
}

#[test]
fn a_roof_flush_with_the_wall_top_removes_nothing() {
    for placement in [Transform3::IDENTITY, building(), general()] {
        // On the top face, and a tenth of the tolerance either side of it:
        // the boolean of a roof moved by at most the tolerance, which keeps
        // the whole wall.
        for offset in [0.0, 0.1 * EPS, -0.1 * EPS] {
            let solid = wall_under_a_flat_roof(offset, placement);
            let measured = volume(&solid);
            assert!(
                (measured - L * T * H).abs() <= L * T * EPS,
                "offset {offset}: volume {measured}"
            );
            assert_eq!(solid.topology().faces().len(), 6, "offset {offset}");
        }
        // Ten tolerances below it, the layer goes.
        let solid = wall_under_a_flat_roof(-10.0 * EPS, placement);
        let expected = L * T * (H - 10.0 * EPS);
        let measured = volume(&solid);
        assert!(
            (measured - expected).abs() <= 1e-9 * expected,
            "volume {measured}, expected {expected}: the layer must go"
        );
    }
}

/// A round column of radius `r` clipped by the vertical plane `x = r -
/// inset` from outside, compiled.
fn column_by_a_tangent_plane(r: f64, inset: f64) -> Result<ExactBRep, GeomError> {
    let placement = general();
    let mut g = Graph::new();
    let column = g.extrusion(
        Profile::Circle(CircleProfile {
            radius: r,
            thickness: None,
        }),
        H,
    );
    let column = g.instance(column, placement);
    let plane = (Point3::new(r - inset, 0.0, 0.0), Vec3::X);
    let cut = g.half_space(plane, true, placement, true);
    let body = g.minus(column, cut);
    let graph = g.finish(vec![body]);
    let solid = compile(&graph, &[body])?.remove(0);
    audited(&solid);
    Ok(solid)
}

#[test]
fn a_plane_tangent_to_a_column_removes_nothing_and_ten_tolerances_in_cut_a_sliver() {
    let r = 0.3;
    let whole = PI * r * r * H;
    // Tangent, and a tenth of the tolerance outside.
    for inset in [0.0, -0.1 * EPS] {
        let measured = volume(&column_by_a_tangent_plane(r, inset).expect("exact"));
        assert!(
            (measured - whole).abs() <= 1e-9 * whole,
            "inset {inset}: volume {measured}"
        );
    }
    // A fraction of the tolerance inside reads as touching: the plane moves
    // out to the column, its chords across the caps meet the one contact
    // ruling, and the whole column is kept (#243).
    for inset in [0.1 * EPS, 0.5 * EPS, 0.9 * EPS] {
        let solid = column_by_a_tangent_plane(r, inset).expect("read as touching");
        let measured = volume(&solid);
        assert!(
            (measured - whole).abs() <= 2.0 * r * H * EPS,
            "inset {inset}: volume {measured}"
        );
    }
    // A circular segment of height 10 eps is cut off.
    let inset = 10.0 * EPS;
    let d = r - inset;
    let segment = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let expected = whole - segment * H;
    let measured = volume(&column_by_a_tangent_plane(r, inset).expect("exact"));
    assert!(
        (measured - expected).abs() <= 1e-3 * segment * H,
        "volume {measured}, expected {expected}: the sliver must go"
    );
}

#[test]
fn a_round_column_flattened_on_both_sides() {
    // Two clips `c` deep into opposite sides of the column: the envelope
    // must reach round the whole circle, not stop at its seam.
    let (r, c, h): (f64, f64, f64) = (0.4, 0.1, 1.0);
    let segment = r * r * ((r - c) / r).acos() - (r - c) * (r * r - (r - c).powi(2)).sqrt();
    let placement = building();
    let mut g = Graph::new();
    let column = g.extrusion(
        Profile::Circle(CircleProfile {
            radius: r,
            thickness: None,
        }),
        h,
    );
    let mut body = g.instance(column, placement);
    for normal in [Vec3::X, -Vec3::X] {
        let side = g.half_space((normal * (r - c), normal), true, placement, true);
        body = g.minus(body, side);
    }
    check(g, body, h * (PI * r * r - 2.0 * segment), &[]);
}

#[test]
fn a_roof_missing_the_wall_keeps_it_and_one_over_it_all_empties_it() {
    let placement = general();
    for (z, keeps) in [(H + 0.5, true), (-0.5, false)] {
        let mut g = Graph::new();
        let wall = g.wall(placement);
        let plane = (Point3::new(0.0, 0.0, z), Vec3::new(0.05, 0.0, 1.0));
        let roof = g.half_space(plane, true, placement, true);
        let body = g.minus(wall, roof);
        if keeps {
            check(g, body, L * T * H, &[]);
        } else {
            let error = refusal(g, body);
            assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
        }
    }
}

fn refusal(graph: Graph, root: NodeId) -> GeomError {
    let graph = graph.finish(vec![root]);
    compile(&graph, &[root]).expect_err("refused")
}

fn named(error: &GeomError, part: &str) {
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains(part)),
        "{error:?} must name `{part}`"
    );
}

#[test]
fn unsupported_clips_are_refused_by_name() {
    let placement = building();
    let plane = roof_plane(2.4, 0.1);

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let roof = g.half_space(plane, true, placement, true);
    let union = g.boolean(wall, roof, BooleanOperator::Union);
    named(&refusal(g, union), "union with a half-space");

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let roof = g.half_space(plane, true, placement, true);
    let swapped = g.boolean(roof, wall, BooleanOperator::Intersection);
    named(&refusal(g, swapped), "not an extrusion");

    let mut g = Graph::new();
    let roof = g.half_space(plane, true, placement, true);
    let floor = g.half_space((Point3::ZERO, Vec3::Z), false, placement, true);
    let slab = g.boolean(roof, floor, BooleanOperator::Intersection);
    named(&refusal(g, slab), "subject is a half-space");

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let roof = g.half_space(plane, true, Transform3::from_scale(Vec3::splat(2.0)), false);
    let scaled = g.minus(wall, roof);
    named(&refusal(g, scaled), "scaled");

    let mut g = Graph::new();
    let profile = g.push(GeometryNode::Profile(rect(0.5, 0.5)));
    let ring = g.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
        profile,
        axis_origin: Point3::new(-2.0, 0.0, 0.0),
        axis_direction: Vec3::Y,
        angle: 2.0 * PI,
    }));
    let roof = g.half_space((Point3::ZERO, Vec3::Z), true, Transform3::IDENTITY, true);
    let clipped = g.minus(ring, roof);
    named(&refusal(g, clipped), "not built from placed extrusions");

    let mut g = Graph::new();
    let wall = g.wall(placement);
    let half_space = g.push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: plane.0,
            normal: plane.1,
        },
        agreement: true,
    }));
    let circle = g.push(GeometryNode::Curve2(Curve2::Circle(Circle2 {
        frame: axiolid_core::Frame2 {
            origin: Point2::ZERO,
            x: axiolid_core::Vec2::X,
            y: axiolid_core::Vec2::Y,
        },
        radius: 1.0,
    })));
    let disc = g.push(GeometryNode::SolidOperation(
        SolidOperation::BoundedHalfSpace {
            half_space,
            boundary: circle,
            placement: Transform3::IDENTITY,
        },
    ));
    let clipped = g.minus(wall, disc);
    named(&refusal(g, clipped), "boundary is not a polyline");

    // On its own a half-space, bounded or not, is unbounded.
    let mut g = Graph::new();
    let roof = g.half_space(plane, true, placement, true);
    named(&refusal(g, roof), "half-space");
    let mut g = Graph::new();
    let square = vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
    ];
    let bounded = g.bounded(plane, true, square, Transform3::IDENTITY, placement);
    named(&refusal(g, bounded), "bounded half-space");
}
