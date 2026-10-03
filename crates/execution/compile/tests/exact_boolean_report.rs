//! Exact placed differences at `Tolerance::ZERO`, and the report of what a
//! difference read within tolerance (#236).
//!
//! A wall `[-L/2, L/2] x [-T/2, T/2] x [0, H]` loses an opening extruded
//! across it, along the wall's `-y`. Operands placed by matrices whose
//! entries are exactly `0` and `+-1`, with dyadic offsets, meet in faces
//! that are exactly coplanar, parallel or perpendicular: the difference is
//! decided by exact predicates, succeeds at `Tolerance::ZERO`, and reports
//! no within-tolerance decision. Under a general (rounded) rotation the
//! flush opening's caps are coplanar with the wall's faces only up to
//! rounding: at a positive tolerance that is read, and reported.
//!
//! Each case checks the exact volume, the certified distance from the
//! result to a probe through the opening, a clean audit and a closed
//! two-manifold.

use std::f64::consts::FRAC_PI_2;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Mat3, Plane3, Point3, Tolerance, Transform3, Vec3};
use axiolid_measure::{boundary_distance, exact_properties};
use axiolid_mesh_compile::{
    BooleanReport, ReferenceExactCompiler, ToleranceDecisionKind, ROUNDING_FACTOR,
};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Profile, RectangleProfile};

/// Wall length, thickness and height: dyadic, so every sum below is exact.
const L: f64 = 6.0;
const T: f64 = 0.25;
const H: f64 = 3.0;
/// Opening width, height, centre along the wall, and sill.
const W: f64 = 1.25;
const OH: f64 = 1.5;
const CX: f64 = 0.75;
const SILL: f64 = 0.5;

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

    fn extrusion(&mut self, x: f64, y: f64, depth: f64) -> NodeId {
        let profile = self.push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x,
                y,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn instance(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    fn minus(&mut self, left: NodeId, right: NodeId) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator: axiolid_core::BooleanOperator::Difference,
        }))
    }

    fn wall(&mut self, placement: Transform3) -> NodeId {
        let wall = self.extrusion(L, T, H);
        self.instance(wall, placement)
    }

    /// An opening `W x OH` at wall `x = cx`, from `y = start` along `-y`
    /// for `depth`, turned across the wall by `across`.
    fn opening(
        &mut self,
        (cx, depth, start): (f64, f64, f64),
        across: Transform3,
        placement: Transform3,
    ) -> NodeId {
        let opening = self.extrusion(W, OH, depth);
        let local = self.instance(
            opening,
            Transform3::from_translation(Vec3::new(cx, start, SILL + OH / 2.0)) * across,
        );
        self.instance(local, placement)
    }

    /// A box `[min, max]` in the wall's frame, under `placement`.
    fn probe(&mut self, min: Vec3, max: Vec3, placement: Transform3) -> NodeId {
        let size = max - min;
        let block = self.extrusion(size.x, size.y, size.z);
        let centre = Vec3::new((min.x + max.x) / 2.0, (min.y + max.y) / 2.0, min.z);
        self.instance(block, placement * Transform3::from_translation(centre))
    }

    fn finish(self, roots: Vec<NodeId>) -> GeometryGraph {
        self.builder.finish(roots).expect("a valid graph")
    }
}

/// A quarter turn about `x`, every entry exactly `0` or `+-1`: the
/// opening's extrusion `+z` becomes the wall's `-y`.
fn exact_across() -> Transform3 {
    Transform3::from_mat3(Mat3::from_cols(Vec3::X, Vec3::Z, -Vec3::Y))
}

/// The same turn from sines and cosines, as a general modeller builds it:
/// `cos(pi / 2)` leaves `6e-17` in the matrix.
fn rounded_across() -> Transform3 {
    Transform3::from_rotation_x(FRAC_PI_2)
}

/// Wall placements with exact entries: none, and a quarter or half turn
/// about `z` with a dyadic offset, as a building model places walls along
/// its grid.
fn exact_placements() -> [(&'static str, Transform3); 3] {
    let quarter = Mat3::from_cols(Vec3::Y, -Vec3::X, Vec3::Z);
    let half = Mat3::from_cols(-Vec3::X, -Vec3::Y, Vec3::Z);
    [
        ("identity", Transform3::IDENTITY),
        (
            "quarter",
            Transform3::from_mat3_translation(quarter, Vec3::new(12.5, -4.0, 3.25)),
        ),
        (
            "half",
            Transform3::from_mat3_translation(half, Vec3::new(-7.75, 2.5, 0.0)),
        ),
    ]
}

/// A rigid placement with no axis left aligned.
fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

/// The three openings: through (overshooting both faces), flush (caps in
/// both faces) and blind (from outside the `+y` face to the wall's middle
/// plane), with the volume each removes.
fn openings() -> [(&'static str, (f64, f64, f64), f64); 3] {
    [
        ("through", (CX, T + 0.5, T / 2.0 + 0.25), W * OH * T),
        ("flush", (CX, T, T / 2.0), W * OH * T),
        (
            "blind",
            (CX, 0.25 + T / 2.0, T / 2.0 + 0.25),
            W * OH * T / 2.0,
        ),
    ]
}

fn compile(
    graph: &GeometryGraph,
    roots: &[NodeId],
    tolerance: Tolerance,
) -> Result<Vec<(ExactBRep, BooleanReport)>, GeomError> {
    ReferenceExactCompiler::new().compile_exact_batch_with_reports(
        graph,
        roots,
        &ExecutionOptions::new(tolerance),
    )
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

/// A wall under `placement` minus one opening, and a probe that enters the
/// opening `gap` inside its jambs, head and sill, from outside the `+y`
/// face to `stop` (wall `y`): the certified distance from the cut wall to
/// it is `gap` exactly when the opening is where it should be.
struct Case {
    graph: GeometryGraph,
    body: NodeId,
    probe: NodeId,
}

const GAP: f64 = 0.0625;

fn case(opening: (f64, f64, f64), across: Transform3, placement: Transform3, stop: f64) -> Case {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let cut = g.opening(opening, across, placement);
    let body = g.minus(wall, cut);
    let probe = g.probe(
        Vec3::new(CX - W / 2.0 + GAP, stop, SILL + GAP),
        Vec3::new(CX + W / 2.0 - GAP, T, SILL + OH - GAP),
        placement,
    );
    Case {
        graph: g.finish(vec![body, probe]),
        body,
        probe,
    }
}

/// Where a probe through the opening stops: past the far face for a
/// through or flush opening, `GAP` short of the floor of a blind one.
fn stop(name: &str) -> f64 {
    if name == "blind" {
        GAP
    } else {
        -T
    }
}

/// Compile a case and check its volume, its probe distance and its audit;
/// the report is returned.
fn checked(case: &Case, removed: f64, tolerance: Tolerance, slack: f64) -> BooleanReport {
    let mut solids =
        compile(&case.graph, &[case.body, case.probe], tolerance).expect("compiles exactly");
    let (probe, _) = solids.pop().expect("probe");
    let (body, report) = solids.pop().expect("body");
    audited(&body);
    let expected = L * T * H - removed;
    let measured = volume(&body);
    assert!(
        (measured - expected).abs() <= slack + 1e-12 * expected,
        "volume {measured}, expected {expected}"
    );
    let bounds = boundary_distance(&body, &probe, 1e-9, Tolerance::METRE).expect("bounded");
    assert!(
        bounds.lower <= GAP + 1e-12 + slack && GAP <= bounds.upper + 1e-12 + slack,
        "[{}, {}] must contain {GAP}",
        bounds.lower,
        bounds.upper
    );
    report
}

#[test]
fn exact_axis_openings_are_exact_at_zero_tolerance() {
    for (place, placement) in exact_placements() {
        for (name, opening, removed) in openings() {
            let c = case(opening, exact_across(), placement, stop(name));
            let report = checked(&c, removed, Tolerance::ZERO, 0.0);
            assert!(report.is_exact(), "{place} {name}: {report:?}");
        }
    }
}

#[test]
fn exact_axis_openings_read_nothing_within_a_positive_tolerance_either() {
    // The exact predicates decide first: a positive tolerance changes
    // nothing for operands that are exactly coplanar or perpendicular.
    for (place, placement) in exact_placements() {
        for (name, opening, removed) in openings() {
            let c = case(opening, exact_across(), placement, stop(name));
            let report = checked(&c, removed, Tolerance::METRE, 0.0);
            assert!(report.is_exact(), "{place} {name}: {report:?}");
        }
    }
}

#[test]
fn two_exact_openings_compose_exactly_at_zero_tolerance() {
    for (place, placement) in exact_placements() {
        let mut g = Graph::new();
        let mut body = g.wall(placement);
        for cx in [-1.5, 1.5] {
            let flush = g.opening((cx, T, T / 2.0), exact_across(), placement);
            body = g.minus(body, flush);
        }
        let graph = g.finish(vec![body]);
        let (solid, report) = ReferenceExactCompiler::new()
            .compile_exact_with_report(&graph, body, &ExecutionOptions::new(Tolerance::ZERO))
            .expect("exact");
        audited(&solid);
        assert!(report.is_exact(), "{place}: {report:?}");
        let expected = L * T * H - 2.0 * W * OH * T;
        assert!((volume(&solid) - expected).abs() <= 1e-12 * expected);
    }
}

#[test]
fn a_rounded_flush_opening_is_read_within_tolerance_and_reported() {
    // The quarter turn from `cos(pi / 2)`, under a general placement: the
    // flush caps miss the wall's faces by rounding. At a positive tolerance
    // they are read as coplanar, and the report says so, within it.
    for placement in [Transform3::IDENTITY, general()] {
        let (name, opening, removed) = openings()[1];
        let c = case(opening, rounded_across(), placement, stop(name));
        let eps = Tolerance::METRE.linear();
        let report = checked(&c, removed, Tolerance::METRE, W * OH * eps);
        assert!(!report.is_exact(), "{report:?}");
        assert!(report.contains(ToleranceDecisionKind::CoincidentSupports));
        assert!(report.linear() > 0.0 && report.linear() <= eps);
        assert!(report.angular() <= Tolerance::METRE.angular());
    }
}

#[test]
fn a_rounded_flush_opening_is_refused_at_zero_tolerance() {
    // Nothing may be read within a zero tolerance: the caps that miss the
    // wall's faces by rounding are refused by name, never guessed.
    let (name, opening, _) = openings()[1];
    let c = case(opening, rounded_across(), general(), stop(name));
    let error = compile(&c.graph, &[c.body], Tolerance::ZERO).expect_err("refused");
    assert!(
        matches!(&error, GeomError::UnsupportedInput { input, .. } if input.contains("exact boolean")),
        "{error:?}"
    );
}

#[test]
fn a_body_without_a_general_boolean_reports_exact() {
    let mut g = Graph::new();
    let wall = g.wall(general());
    let graph = g.finish(vec![wall]);
    let (_, report) = ReferenceExactCompiler::new()
        .compile_exact_with_report(&graph, wall, &ExecutionOptions::new(Tolerance::METRE))
        .expect("exact");
    assert!(report.is_exact());
}

#[test]
fn a_placed_cut_keeps_its_report() {
    // An instance of a cut wall carries the cut's report; so does the
    // same body compiled again from the cache.
    let mut g = Graph::new();
    let (_, opening, _) = openings()[1];
    let wall = g.wall(Transform3::IDENTITY);
    let cut = g.opening(opening, rounded_across(), general());
    let wall = g.instance(wall, general());
    let body = g.minus(wall, cut);
    let placed = g.instance(body, Transform3::from_translation(Vec3::X));
    let graph = g.finish(vec![body, placed]);
    let reports = compile(&graph, &[body, placed, body], Tolerance::METRE).expect("exact");
    assert!(!reports[0].1.is_exact());
    assert_eq!(reports[0].1, reports[1].1);
    assert_eq!(reports[0].1, reports[2].1);
}

/// The wall under `placement` minus the half-space above its top face,
/// the plane given in world coordinates: a roof flush with the wall.
fn wall_under_a_flush_roof(
    placement: Transform3,
    tolerance: Tolerance,
) -> (ExactBRep, BooleanReport) {
    let mut g = Graph::new();
    let wall = g.wall(placement);
    let roof = g.push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: placement.transform_point3(Point3::new(0.0, 0.0, H)),
            normal: placement.transform_vector3(Vec3::Z),
        },
        agreement: true,
    }));
    let body = g.minus(wall, roof);
    let graph = g.finish(vec![body]);
    let (solid, report) = ReferenceExactCompiler::new()
        .compile_exact_with_report(&graph, body, &ExecutionOptions::new(tolerance))
        .expect("exact");
    audited(&solid);
    assert!((volume(&solid) - L * T * H).abs() <= L * T * tolerance.linear() + 1e-12);
    (solid, report)
}

#[test]
fn a_clip_reports_what_it_read() {
    // Placed by exact axes, the roof plane is exactly the wall's top face:
    // nothing read, at zero tolerance too. Under a general placement it
    // is the top face only up to rounding: read as coplanar, and reported.
    for (_, placement) in exact_placements() {
        let (_, report) = wall_under_a_flush_roof(placement, Tolerance::ZERO);
        assert!(report.is_exact(), "{report:?}");
    }
    let (_, report) = wall_under_a_flush_roof(general(), Tolerance::METRE);
    assert!(
        report.contains(ToleranceDecisionKind::CoincidentSupports),
        "{report:?}"
    );
}

#[test]
fn a_body_reports_the_largest_rounding_floor_beneath_it() {
    // #244. A wall with two exact openings: the first through,
    // overshooting to `y = +-10`, so its difference's extent is `10`; the
    // second flush, cut from that result, whose extent is the wall's
    // `L / 2 = 3`. The body's floor is the larger, at zero tolerance too,
    // and an instance placed far off keeps it: a rigid placement leaves
    // distances between points unchanged.
    let mut g = Graph::new();
    let wall = g.wall(Transform3::IDENTITY);
    let through = g.opening((-1.5, 20.0, 10.0), exact_across(), Transform3::IDENTITY);
    let first = g.minus(wall, through);
    let flush = g.opening((1.5, T, T / 2.0), exact_across(), Transform3::IDENTITY);
    let body = g.minus(first, flush);
    let placed = g.instance(
        body,
        Transform3::from_translation(Vec3::new(100.0, 0.0, 0.0)),
    );
    let alone = g.minus(wall, flush);
    let graph = g.finish(vec![body, placed, alone]);
    let reports = compile(&graph, &[body, placed, alone], Tolerance::ZERO).expect("exact");
    for (solid, report) in &reports {
        audited(solid);
        assert!(report.is_exact(), "{report:?}");
    }
    let expected = L * T * H - 2.0 * W * OH * T;
    assert!((volume(&reports[0].0) - expected).abs() <= 1e-12 * expected);
    assert_eq!(reports[0].1.extent(), 10.0, "{:?}", reports[0].1);
    assert_eq!(reports[0].1.rounding_floor(), 10.0 * ROUNDING_FACTOR);
    assert_eq!(reports[1].1, reports[0].1);
    // The flush cut alone: the wall's own extent.
    assert_eq!(reports[2].1.extent(), L / 2.0, "{:?}", reports[2].1);
    assert_eq!(reports[2].1.rounding_floor(), L / 2.0 * ROUNDING_FACTOR);
}
