//! Exact differences of placed extrusions: openings in walls and slabs
//! (#228).
//!
//! A wall is a rectangle extruded along its local `z`, `L` long in `x` and
//! `T` thick in `y`. A door or window opening is extruded across it, along
//! the wall's `-y` (perpendicular to the wall's extrusion); a shaft through
//! a slab runs along the slab's own `z` (parallel). Every operand is placed
//! by an `Instance`, and the wall's placement is shared by its openings, as
//! a building model places an opening relative to its wall.
//!
//! Each case checks, against closed forms:
//!
//! - the exact volume (`exact_properties`);
//! - the certified distance between the result's boundary and a probe box
//!   that passes through the opening without touching it. Without the cut
//!   the probe would cross the wall's boundary and the distance would be
//!   zero, so the distance pins the opening's position and size;
//! - a clean geometric audit and a closed two-manifold topology.

use std::f64::consts::{FRAC_PI_2, PI};

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Interval, Point2, Point3, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::{boundary_distance, exact_properties};
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile,
};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

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

    /// The wall, `[-L/2, L/2] x [-T/2, T/2] x [0, H]` in its own frame,
    /// under `placement`.
    fn wall(&mut self, placement: Transform3) -> NodeId {
        let wall = self.extrusion(rect(L, T), H);
        self.instance(wall, placement)
    }

    /// An opening across the wall: `profile` (its `x` along the wall, its
    /// `y` up) extruded along the wall's `-y` from `y = start` for `depth`,
    /// its profile origin at wall `x = cx`, `z = cz`, then under the wall's
    /// `placement`.
    fn opening(
        &mut self,
        profile: Profile,
        depth: f64,
        (cx, cz): (f64, f64),
        start: f64,
        placement: Transform3,
    ) -> NodeId {
        let opening = self.extrusion(profile, depth);
        let across = self.instance(
            opening,
            Transform3::from_translation(Vec3::new(cx, start, cz))
                * Transform3::from_rotation_x(FRAC_PI_2),
        );
        self.instance(across, placement)
    }

    /// A box `[x0, x1] x [y0, y1] x [z0, z1]` in the wall's frame, under the
    /// wall's `placement`.
    fn probe(&mut self, min: Vec3, max: Vec3, placement: Transform3) -> NodeId {
        let size = max - min;
        let block = self.extrusion(rect(size.x, size.y), size.z);
        let centre = Vec3::new((min.x + max.x) / 2.0, (min.y + max.y) / 2.0, min.z);
        self.instance(block, placement * Transform3::from_translation(centre))
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

fn line(a: Point2, b: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: a,
            direction: Vec2::new(b.x - a.x, b.y - a.y),
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    }
}

/// An arched opening: a `w x h` rectangle on `y = 0` under a half circle of
/// radius `w / 2` centred at `(0, h)`, as `parts` equal arc segments (one is
/// how an IFC arch commonly comes: a single semicircle).
fn arched_in(w: f64, h: f64, parts: usize) -> Profile {
    let r = w / 2.0;
    let p = |x: f64, y: f64| Point2::new(x, y);
    let step = PI / parts as f64;
    let part = |k: usize| ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: axiolid_core::Frame2 {
                origin: p(0.0, h),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: r,
        }),
        domain: Interval::new(k as f64 * step, (k + 1) as f64 * step),
        same_sense: true,
    };
    let mut segments = vec![line(p(-r, 0.0), p(r, 0.0)), line(p(r, 0.0), p(r, h))];
    segments.extend((0..parts).map(part));
    segments.push(line(p(-r, h), p(-r, 0.0)));
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

/// An arched opening whose half circle is one segment.
fn arched(w: f64, h: f64) -> Profile {
    arched_in(w, h, 1)
}

/// A wall placement as a building model has one: turned about `z` and
/// moved.
fn building() -> Transform3 {
    Transform3::from_translation(Vec3::new(12.5, -4.0, 3.2)) * Transform3::from_rotation_z(0.6)
}

/// A rigid placement with no axis left aligned.
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

/// Compile `body` and each probe; check the volume, the audit, and each
/// probe's certified distance to the body.
fn check(graph: Graph, body: NodeId, expected_volume: f64, probes: &[(NodeId, f64)]) {
    let mut roots = vec![body];
    roots.extend(probes.iter().map(|&(probe, _)| probe));
    let graph = graph.finish(roots.clone());
    let solids = compile(&graph, &roots).expect("exact");
    audited(&solids[0]);
    let measured = volume(&solids[0]);
    assert!(
        (measured - expected_volume).abs() <= 1e-9 * expected_volume,
        "volume {measured}, expected {expected_volume}"
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
}

/// A probe through a rectangular opening `[x0, x1] x [z0, z1]`, inset by
/// `gap` on every side and sticking out of the wall on both faces.
fn through_probe(
    g: &mut Graph,
    (x0, x1): (f64, f64),
    (z0, z1): (f64, f64),
    gap: f64,
    p: Transform3,
) -> NodeId {
    g.probe(
        Vec3::new(x0 + gap, -T, z0 + gap),
        Vec3::new(x1 - gap, T, z1 - gap),
        p,
    )
}

#[test]
fn a_perpendicular_window_cuts_through_a_placed_wall() {
    for placement in [building(), general()] {
        let mut g = Graph::new();
        let (w, h, cx, sill) = (1.2, 1.4, 0.8, 0.9);
        let wall = g.wall(placement);
        // Overshooting both wall faces.
        let opening = g.opening(
            rect(w, h),
            T + 0.2,
            (cx, sill + h / 2.0),
            T / 2.0 + 0.1,
            placement,
        );
        let cut = g.minus(wall, opening);
        let probe = through_probe(
            &mut g,
            (cx - w / 2.0, cx + w / 2.0),
            (sill, sill + h),
            0.05,
            placement,
        );
        check(g, cut, L * T * H - w * h * T, &[(probe, 0.05)]);
    }
}

#[test]
fn an_opening_exactly_as_deep_as_the_wall_is_thick_cuts_through() {
    // The opening's caps lie in the wall's faces: coplanar contact, which
    // two independent placements leave coplanar only up to rounding.
    for placement in [Transform3::IDENTITY, building(), general()] {
        let mut g = Graph::new();
        let (w, h, cx, sill) = (0.9, 1.1, -1.5, 1.0);
        let wall = g.wall(placement);
        let opening = g.opening(rect(w, h), T, (cx, sill + h / 2.0), T / 2.0, placement);
        let cut = g.minus(wall, opening);
        let probe = through_probe(
            &mut g,
            (cx - w / 2.0, cx + w / 2.0),
            (sill, sill + h),
            0.02,
            placement,
        );
        check(g, cut, L * T * H - w * h * T, &[(probe, 0.02)]);
    }
}

#[test]
fn a_blind_perpendicular_recess_leaves_a_floor() {
    let placement = general();
    let mut g = Graph::new();
    let (w, h, cx, sill, d) = (1.0, 0.8, 0.0, 1.2, 0.12);
    let wall = g.wall(placement);
    // From outside the wall's +y face, ending d inside it.
    let opening = g.opening(
        rect(w, h),
        0.1 + d,
        (cx, sill + h / 2.0),
        T / 2.0 + 0.1,
        placement,
    );
    let cut = g.minus(wall, opening);
    // The probe enters the recess from outside and stops 0.03 short of
    // its floor at y = T/2 - d; its sides are 0.05 from the recess walls.
    let floor = T / 2.0 - d;
    let probe = g.probe(
        Vec3::new(cx - w / 2.0 + 0.05, floor + 0.03, sill + 0.05),
        Vec3::new(cx + w / 2.0 - 0.05, T, sill + h - 0.05),
        placement,
    );
    check(g, cut, L * T * H - w * h * d, &[(probe, 0.03)]);
}

#[test]
fn an_arched_opening_cuts_through_with_a_cylindrical_soffit() {
    for placement in [building(), general()] {
        for parts in [1, 2] {
            arched_opening(placement, parts);
        }
    }
}

fn arched_opening(placement: Transform3, parts: usize) {
    let mut g = Graph::new();
    let (w, h, cx, sill) = (1.0, 1.5, 1.0, 0.3);
    let r = w / 2.0;
    let wall = g.wall(placement);
    let opening = g.opening(
        arched_in(w, h, parts),
        T + 0.2,
        (cx, sill),
        T / 2.0 + 0.1,
        placement,
    );
    let cut = g.minus(wall, opening);
    // The probe's top corners lie inside the half circle; the arc is the
    // nearest part of the opening: r - |corner - centre|.
    let (a, top) = (0.3, sill + h + 0.3);
    let probe = g.probe(
        Vec3::new(cx - a, -T, sill + 0.2),
        Vec3::new(cx + a, T, top),
        placement,
    );
    let to_arc = r - (a * a + 0.3 * 0.3_f64).sqrt();
    assert!(to_arc < 0.2);
    let area = w * h + PI * r * r / 2.0;
    check(g, cut, L * T * H - area * T, &[(probe, to_arc)]);
}

#[test]
fn a_parallel_shaft_cuts_through_a_placed_slab() {
    let placement = general();
    let mut g = Graph::new();
    let (sx, sy, st) = (4.0, 3.0, 0.25);
    let slab = g.extrusion(rect(sx, sy), st);
    let slab = g.instance(slab, placement);
    // A shaft turned in the slab's plane, flush with both slab faces.
    let (w, d, turn, at) = (0.6, 0.4, 0.3, Vec3::new(0.7, -0.5, 0.0));
    let shaft = g.extrusion(rect(w, d), st);
    let local = Transform3::from_translation(at) * Transform3::from_rotation_z(turn);
    let shaft = g.instance(shaft, placement * local);
    let cut = g.minus(slab, shaft);
    let gap = 0.05;
    let probe = g.extrusion(rect(w - 2.0 * gap, d - 2.0 * gap), 1.0);
    let probe = g.instance(
        probe,
        placement * local * Transform3::from_translation(Vec3::new(0.0, 0.0, -0.4)),
    );
    check(g, cut, (sx * sy - w * d) * st, &[(probe, gap)]);
}

#[test]
fn a_round_shaft_runs_parallel_through_a_slab() {
    let placement = building();
    let mut g = Graph::new();
    let (sx, sy, st, r) = (4.0, 3.0, 0.25, 0.3);
    let slab = g.extrusion(rect(sx, sy), st);
    let slab = g.instance(slab, placement);
    let shaft = g.extrusion(
        Profile::Circle(CircleProfile {
            radius: r,
            thickness: None,
        }),
        st + 0.5,
    );
    let local = Transform3::from_translation(Vec3::new(-1.0, 0.4, -0.25));
    let shaft = g.instance(shaft, placement * local);
    let cut = g.minus(slab, shaft);
    // A square probe in the hole: its vertical edges are r - a sqrt 2 from
    // the cylinder.
    let a = 0.15;
    let probe = g.extrusion(rect(2.0 * a, 2.0 * a), 2.0);
    let probe = g.instance(
        probe,
        placement * Transform3::from_translation(Vec3::new(-1.0, 0.4, -1.0)),
    );
    check(
        g,
        cut,
        (sx * sy - PI * r * r) * st,
        &[(probe, r - a * 2.0_f64.sqrt())],
    );
}

#[test]
fn two_and_three_openings_compose_in_one_wall() {
    let placement = general();
    let openings = [
        (rect(0.9, 2.1), (-2.0, 1.05), 0.9 * 2.1),
        (rect(1.2, 1.2), (0.2, 1.6), 1.2 * 1.2),
        (arched(0.8, 1.0), (2.0, 0.8), 0.8 + PI * 0.16 / 2.0),
    ];
    for count in [2, 3] {
        let mut g = Graph::new();
        let mut body = g.wall(placement);
        let mut removed = 0.0;
        let mut probes = Vec::new();
        for (profile, at, area) in openings.iter().take(count) {
            let opening = g.opening(profile.clone(), T + 0.2, *at, T / 2.0 + 0.1, placement);
            body = g.minus(body, opening);
            removed += area * T;
        }
        // A probe through the first opening, a door whose threshold lies in
        // the wall's base, and one through the second.
        probes.push((
            through_probe(&mut g, (-2.45, -1.55), (0.0, 2.1), 0.04, placement),
            0.04,
        ));
        probes.push((
            through_probe(&mut g, (-0.4, 0.8), (1.0, 2.2), 0.07, placement),
            0.07,
        ));
        check(g, body, L * T * H - removed, &probes);
    }
}

#[test]
fn an_opening_touching_the_wall_end_and_base_leaves_a_notch() {
    // A door at the wall's end: its side lies in the wall's end face and
    // its threshold in the wall's base.
    let placement = building();
    let mut g = Graph::new();
    let (w, h) = (1.0, 2.0);
    let cx = L / 2.0 - w / 2.0;
    let wall = g.wall(placement);
    let opening = g.opening(rect(w, h), T, (cx, h / 2.0), T / 2.0, placement);
    let cut = g.minus(wall, opening);
    // A probe in the notch: 0.1 below its lintel, 0.1 from its inner side,
    // reaching outside past the end and the base.
    let probe = g.probe(
        Vec3::new(cx - w / 2.0 + 0.1, -T, -0.5),
        Vec3::new(L / 2.0 + 0.5, T, h - 0.1),
        placement,
    );
    check(g, cut, L * T * H - w * h * T, &[(probe, 0.1)]);
}

#[test]
fn an_instanced_cut_wall_moves_with_its_placement() {
    // The difference computed in the wall's own frame, then placed whole.
    let mut g = Graph::new();
    let (w, h, cx, sill) = (1.2, 1.4, 0.8, 0.9);
    let wall = g.wall(Transform3::IDENTITY);
    let opening = g.opening(
        rect(w, h),
        T + 0.2,
        (cx, sill + h / 2.0),
        T / 2.0 + 0.1,
        Transform3::IDENTITY,
    );
    let cut = g.minus(wall, opening);
    let placed = g.instance(cut, general());
    let probe = through_probe(
        &mut g,
        (cx - w / 2.0, cx + w / 2.0),
        (sill, sill + h),
        0.05,
        general(),
    );
    check(g, placed, L * T * H - w * h * T, &[(probe, 0.05)]);
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
fn a_tool_that_is_not_an_extrusion_is_refused_by_name() {
    let mut g = Graph::new();
    let wall = g.wall(building());
    let profile = g.push(GeometryNode::Profile(Profile::Circle(CircleProfile {
        radius: 0.2,
        thickness: None,
    })));
    let ring = g.push(GeometryNode::SolidOperation(SolidOperation::Revolution {
        profile,
        axis_origin: Point3::new(-1.0, 0.0, 0.0),
        axis_direction: Vec3::Y,
        angle: 2.0 * PI,
    }));
    let ring = g.instance(ring, building());
    let cut = g.minus(wall, ring);
    named(&refusal(g, cut), "not an extrusion");
}

#[test]
fn a_tool_that_is_itself_a_boolean_is_refused_by_name() {
    let mut g = Graph::new();
    let wall = g.wall(building());
    let a = g.opening(rect(1.0, 1.0), 1.0, (0.0, 1.0), 0.5, building());
    let b = g.opening(rect(0.5, 0.5), 1.0, (0.0, 1.0), 0.5, building());
    let tool = g.minus(a, b);
    let cut = g.minus(wall, tool);
    named(&refusal(g, cut), "tool that is not a placed extrusion");
}

#[test]
fn a_scaled_opening_is_refused_by_name() {
    let mut g = Graph::new();
    let wall = g.wall(building());
    let opening = g.extrusion(rect(1.0, 1.0), 1.0);
    let scaled = g.instance(opening, Transform3::from_scale(Vec3::new(1.0, 1.0, 2.0)));
    let cut = g.minus(wall, scaled);
    named(&refusal(g, cut), "scaled");
}

#[test]
fn a_union_or_intersection_of_placed_operands_is_refused_by_name() {
    for operator in [BooleanOperator::Union, BooleanOperator::Intersection] {
        let mut g = Graph::new();
        let wall = g.wall(building());
        let opening = g.opening(rect(1.0, 1.0), 1.0, (0.0, 1.0), 0.5, building());
        let joined = g.boolean(wall, opening, operator);
        named(&refusal(g, joined), "union or intersection");
    }
}

#[test]
fn an_opening_that_removes_the_whole_wall_is_degenerate() {
    let mut g = Graph::new();
    let wall = g.wall(building());
    let opening = g.opening(
        rect(L + 1.0, H + 1.0),
        T + 1.0,
        (0.0, H / 2.0),
        T / 2.0 + 0.5,
        building(),
    );
    let cut = g.minus(wall, opening);
    let error = refusal(g, cut);
    assert!(matches!(error, GeomError::Degenerate(_)), "{error:?}");
}

#[test]
fn an_opening_that_misses_the_wall_leaves_it_whole() {
    let placement = general();
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let opening = g.opening(
        rect(1.0, 1.0),
        T + 0.2,
        (0.0, H + 2.0),
        T / 2.0 + 0.1,
        placement,
    );
    let cut = g.minus(wall, opening);
    check(g, cut, L * T * H, &[]);
}

#[test]
fn a_configuration_the_general_boolean_refuses_is_refused_by_name() {
    // An elliptical column crossed by a round bore: the general boolean
    // does not build this section, so the compiler names the refusal
    // instead of meshing.
    let mut g = Graph::new();
    let column = g.extrusion(
        Profile::Ellipse(axiolid_profile::EllipseProfile {
            semi_axis_x: 0.5,
            semi_axis_y: 0.3,
        }),
        2.0,
    );
    let column = g.instance(column, building());
    let bore = g.extrusion(
        Profile::Circle(CircleProfile {
            radius: 0.2,
            thickness: None,
        }),
        2.0,
    );
    let bore = g.instance(
        bore,
        building()
            * Transform3::from_translation(Vec3::new(-1.0, 0.0, 1.0))
            * Transform3::from_rotation_y(FRAC_PI_2),
    );
    let cut = g.minus(column, bore);
    named(&refusal(g, cut), "exact boolean");
}

/// The linear tolerance every test here runs at.
const EPS: f64 = 1e-6;

/// A wall with one window `w x h` whose far cap stops `short` before the
/// wall's far face (negative: past it), and its compiled solid.
fn window_stopping_short(short: f64, placement: Transform3) -> ExactBRep {
    let mut g = Graph::new();
    let (w, h) = (1.2, 1.4);
    let wall = g.wall(placement);
    // From 0.1 outside the near face to `short` before the far face.
    let opening = g.opening(
        rect(w, h),
        T + 0.1 - short,
        (0.8, 1.6),
        T / 2.0 + 0.1,
        placement,
    );
    let cut = g.minus(wall, opening);
    let graph = g.finish(vec![cut]);
    let solid = compile(&graph, &[cut]).expect("exact").remove(0);
    audited(&solid);
    solid
}

#[test]
fn a_gap_of_ten_tolerances_is_kept_as_a_gap() {
    // Ten tolerances is beyond every within-tolerance reading: the exact
    // predicates decide, and the window leaves a skin 1e-5 thick.
    for placement in [Transform3::IDENTITY, general()] {
        let skin = 10.0 * EPS;
        let solid = window_stopping_short(skin, placement);
        let expected = L * T * H - 1.2 * 1.4 * (T - skin);
        let measured = volume(&solid);
        assert!(
            (measured - expected).abs() <= 1e-9 * expected,
            "volume {measured}, expected {expected}: the skin must stay"
        );
        // The skin is a face of the result: more faces than the through
        // cut, which has none there.
        let through = window_stopping_short(-0.1, placement);
        assert!(solid.topology().faces().len() > through.topology().faces().len());
    }
}

#[test]
fn a_gap_within_tolerance_reads_as_flush() {
    // A tenth of the tolerance short of the far face: the boolean of the
    // window moved by at most the tolerance, a through cut. Its volume is
    // within the perturbation's (w h eps) of the through cut's, and its
    // faces are the flush cut's.
    for placement in [Transform3::IDENTITY, general()] {
        let solid = window_stopping_short(0.1 * EPS, placement);
        let flush = window_stopping_short(0.0, placement);
        let through = L * T * H - 1.2 * 1.4 * T;
        let measured = volume(&solid);
        assert!(
            (measured - through).abs() <= 1.2 * 1.4 * EPS,
            "volume {measured}, through {through}"
        );
        assert_eq!(
            solid.topology().faces().len(),
            flush.topology().faces().len()
        );
    }
}

/// A door `1 x 2` on the wall's base whose side stops `gap` before the
/// wall's end, compiled.
fn door_short_of_the_end(gap: f64) -> ExactBRep {
    let placement = general();
    let mut g = Graph::new();
    let (w, h) = (1.0, 2.0);
    let wall = g.wall(placement);
    let opening = g.opening(
        rect(w, h),
        T + 0.2,
        (L / 2.0 - w / 2.0 - gap, h / 2.0),
        T / 2.0 + 0.1,
        placement,
    );
    let cut = g.minus(wall, opening);
    let graph = g.finish(vec![cut]);
    let solid = compile(&graph, &[cut]).expect("exact").remove(0);
    audited(&solid);
    let expected = L * T * H - w * h * T;
    let measured = volume(&solid);
    assert!(
        (measured - expected).abs() <= w * h * EPS,
        "volume {measured}, expected {expected}"
    );
    solid
}

#[test]
fn a_wall_end_ten_tolerances_from_an_opening_is_kept() {
    // A door ending 1e-5 before the wall's end leaves a sliver of wall (its
    // own end face, and the base split in two); one ending a tenth of the
    // tolerance before it reads as touching the end, the notch.
    let notch = door_short_of_the_end(0.0).topology().faces().len();
    let sliver = door_short_of_the_end(10.0 * EPS).topology().faces().len();
    let touching = door_short_of_the_end(0.1 * EPS).topology().faces().len();
    assert!(sliver > notch, "{sliver} faces: the sliver must stay");
    assert_eq!(touching, notch);
}
