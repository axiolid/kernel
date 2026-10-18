//! Certified bounds of runs of sections at stations whose frames are read
//! by quadrature (#316): runs along the openbimrs/ifc#443 deck's offset by
//! distances beside an ellipse, along an atomic ellipse, a tapered run
//! along a B-spline and a run along a constant offset of an ellipse each
//! report a certified deviation within the request, and a dense sampling
//! of the EXACT solid -- its frames computed here from the curves in
//! closed form at their true arc length, by a fine quadrature independent
//! of the kernel's -- never strays past it. The bases that stay
//! uncertified are still unbounded by name, and so is the deck's cell
//! beside the mitre where its law starts to widen.

use axiolid_construct::profile::profile_rings;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{BSplineCurve2, CurvatureLaw, Curve2, Ellipse2, Intrinsic2, KnotSpec};
use axiolid_measure::volume_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    CurveRelation, GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, RunStation,
    SectionAtStation, SolidOperation, Station, StationFrame, StationOffsets,
};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_reference::{derivative2, evaluate2};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn rectangle(b: &mut GeometryGraphBuilder, x: Scalar, y: Scalar) -> NodeId {
    b.push_value(Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    }))
    .unwrap()
}

fn spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<SectionAtStation>,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

/// The deck's ellipse: semi-axes 10 and 6 about `(30, 0)`.
fn ellipse() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(30.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    })
}

/// A cubic plan B-spline with an interior knot, about 25 long.
fn spline() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(6.0, 1.0),
            Point2::new(12.0, -2.0),
            Point2::new(18.0, 3.0),
            Point2::new(24.0, 2.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

// --- the exact solid, independently -------------------------------------------

/// A plane curve in closed form: its point and velocity at a parameter.
struct Plane {
    point: Box<dyn Fn(Scalar) -> Point2>,
    velocity: Box<dyn Fn(Scalar) -> Vec2>,
}

impl Plane {
    /// A curve the evaluator reads in closed form.
    fn of(curve: Curve2) -> Self {
        let other = curve.clone();
        Self {
            point: Box::new(move |t| evaluate2(&curve, t).unwrap()),
            velocity: Box::new(move |t| derivative2(&other, t).unwrap()),
        }
    }
}

/// Gauss-Legendre nodes and weights on `[-1, 1]`, five points.
const NODES: [Scalar; 5] = [
    0.0,
    -0.538_469_310_105_683_1,
    0.538_469_310_105_683_1,
    -0.906_179_845_938_664,
    0.906_179_845_938_664,
];
const WEIGHTS: [Scalar; 5] = [
    0.568_888_888_888_888_9,
    0.478_628_670_499_366_5,
    0.478_628_670_499_366_5,
    0.236_926_885_056_189_1,
    0.236_926_885_056_189_1,
];

fn panel(speed: &dyn Fn(Scalar) -> Scalar, a: Scalar, b: Scalar) -> Scalar {
    let (half, mid) = (0.5 * (b - a), 0.5 * (a + b));
    NODES
        .iter()
        .zip(WEIGHTS)
        .map(|(x, w)| w * speed(mid + half * x))
        .sum::<Scalar>()
        * half
}

/// A curve's arc length tabulated from `t0` on a fine grid, read back at
/// any parameter or inverted at any length.
struct Lengths {
    t0: Scalar,
    step: Scalar,
    sums: Vec<Scalar>,
}

impl Lengths {
    fn new(speed: &dyn Fn(Scalar) -> Scalar, t0: Scalar, t1: Scalar, steps: usize) -> Self {
        let step = (t1 - t0) / steps as Scalar;
        let mut sums = vec![0.0];
        for i in 0..steps {
            let a = t0 + step * i as Scalar;
            let piece = panel(speed, a, a + step);
            sums.push(sums[i] + piece);
        }
        Self { t0, step, sums }
    }

    /// The arc length from `t0` to `t`.
    fn at(&self, speed: &dyn Fn(Scalar) -> Scalar, t: Scalar) -> Scalar {
        let i = (((t - self.t0) / self.step).floor().max(0.0) as usize).min(self.sums.len() - 2);
        let a = self.t0 + self.step * i as Scalar;
        self.sums[i] + panel(speed, a, t)
    }

    /// The parameter at arc length `s` from `t0`.
    fn inverse(&self, speed: &dyn Fn(Scalar) -> Scalar, s: Scalar) -> Scalar {
        let i = self
            .sums
            .partition_point(|&x| x <= s)
            .clamp(1, self.sums.len() - 1)
            - 1;
        let mut t = self.t0 + self.step * i as Scalar;
        for _ in 0..30 {
            let error = self.at(speed, t) - s;
            t -= error / speed(t);
            if error.abs() < 1e-13 {
                break;
            }
        }
        t
    }
}

/// The exact solid's wall: sections of `corners` at `stations` (with
/// their offsets), each point interpolated linearly in distance, placed
/// in the planar frame of `curve` at its true arc length from `t0`.
struct Exact {
    curve: Plane,
    lengths: Lengths,
    stations: Vec<(Scalar, StationOffsets)>,
    corners: Vec<Vec<Point2>>,
}

impl Exact {
    fn new(
        curve: Plane,
        (t0, t1): (Scalar, Scalar),
        graph: &GeometryGraph,
        sections: &[SectionAtStation],
    ) -> Self {
        let speed = |t: Scalar| (curve.velocity)(t).length();
        let lengths = Lengths::new(&speed, t0, t1, 100_000);
        Self {
            curve,
            lengths,
            stations: sections
                .iter()
                .map(|s| (s.station.distance, s.station.offsets))
                .collect(),
            corners: sections
                .iter()
                .map(|s| {
                    let Some(GeometryNode::Profile(p)) = graph.get(s.profile) else {
                        panic!("a profile")
                    };
                    profile_rings(p, 0.5e-3, Tolerance::MILLIMETRE)
                        .unwrap()
                        .outer
                })
                .collect(),
        }
    }

    /// Points of every wall edge at `per` distances between stations,
    /// `along` + 1 per edge.
    fn points(&self, per: usize, along: usize) -> Vec<Point3> {
        let speed = |t: Scalar| (self.curve.velocity)(t).length();
        let mut out = Vec::new();
        for i in 0..self.stations.len() - 1 {
            let ((d0, o0), (d1, o1)) = (self.stations[i], self.stations[i + 1]);
            for a in 0..=per {
                let u = a as Scalar / per as Scalar;
                let t = self.lengths.inverse(&speed, d0 + (d1 - d0) * u);
                let p = (self.curve.point)(t);
                let v = (self.curve.velocity)(t);
                let tangent = Vec3::new(v.x, v.y, 0.0) / v.length();
                let lateral = Vec3::new(-tangent.y, tangent.x, 0.0);
                let lerp = |x: Scalar, y: Scalar| x + (y - x) * u;
                let origin = Point3::new(p.x, p.y, 0.0)
                    + lerp(o0.lateral, o1.lateral) * lateral
                    + lerp(o0.vertical, o1.vertical) * Vec3::Z
                    + lerp(o0.longitudinal, o1.longitudinal) * tangent;
                let ring: Vec<Point3> = self.corners[i]
                    .iter()
                    .zip(&self.corners[i + 1])
                    .map(|(p, q)| origin + lerp(p.x, q.x) * lateral + lerp(p.y, q.y) * Vec3::Z)
                    .collect();
                for k in 0..ring.len() {
                    let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                    for c in 0..=along {
                        out.push(a + (b - a) * (c as Scalar / along as Scalar));
                    }
                }
            }
        }
        out
    }
}

/// The largest distance from `points` to the mesh's triangles.
fn measured(mesh: &TriMesh, points: &[Point3]) -> Scalar {
    let triangles: Vec<[Point3; 3]> = mesh
        .indices
        .chunks_exact(3)
        .map(|t| [0, 1, 2].map(|k| mesh.positions[t[k] as usize]))
        .collect();
    points
        .iter()
        .map(|&p| {
            triangles
                .iter()
                .filter_map(|&t| axiolid_measure::proximity::closest_point_on_triangle(p, t).ok())
                .map(|q| p.distance(q))
                .fold(Scalar::INFINITY, Scalar::min)
        })
        .fold(0.0, Scalar::max)
}

/// The run's own deviation in `graph`'s compilation of `root`.
fn run_bound(graph: &GeometryGraph, root: NodeId) -> (TriMesh, DeviationBound, bool) {
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(graph, root, &options())
        .unwrap();
    let bound = report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationedSpine)
        .expect("the run reports its path")
        .bound;
    (outcome.mesh, bound, report.meets_requested())
}

/// Compile `root`: its run certified within the request, a dense sampling
/// of `exact` within the bound; the mesh and the bound.
fn assert_certified(graph: &GeometryGraph, root: NodeId, exact: &Exact) -> (TriMesh, Scalar) {
    let (mesh, bound, meets) = run_bound(graph, root);
    let DeviationBound::Certified(bound) = bound else {
        panic!("not certified: {bound:?}")
    };
    assert!(meets, "{bound} is past the request");
    let far = measured(&mesh, &exact.points(400, 4));
    assert!(
        far <= bound,
        "the exact solid strays {far} from the mesh, past the bound {bound}"
    );
    (mesh, bound)
}

fn unbounded(graph: &GeometryGraph, root: NodeId) -> &'static str {
    match run_bound(graph, root).1 {
        DeviationBound::Unbounded(reason) => reason,
        other => panic!("bounded: {other:?}"),
    }
}

// --- the deck -------------------------------------------------------------------

/// The deck's directrix: the ellipse offset by `0.5` to station 2,
/// widening to `0.75` at station 20 and held to the end.
fn deck_directrix(b: &mut GeometryGraphBuilder) -> NodeId {
    let basis = b.push_value(ellipse()).unwrap();
    let lateral = |d| StationOffsets::new(d, 0.0, 0.0);
    b.push(GeometryNode::CurveRelation(
        CurveRelation::OffsetsAtStations {
            basis,
            stations: vec![
                RunStation::at_basis_start(lateral(0.5)),
                RunStation::new(Station::new(2.0, lateral(0.5))),
                RunStation::new(Station::new(20.0, lateral(0.75))),
                RunStation::at_basis_end(lateral(0.75)),
            ],
            frame: StationFrame::Section,
        },
    ))
    .unwrap()
}

/// The deck's directrix exactly: the law runs linearly in the ellipse's
/// arc length, its point beside the ellipse's along the left normal.
fn deck_curve() -> Plane {
    fn e(t: Scalar) -> Vec2 {
        Vec2::new(-10.0 * t.sin(), 6.0 * t.cos())
    }
    fn base_speed(t: Scalar) -> Scalar {
        e(t).length()
    }
    let base = std::rc::Rc::new(Lengths::new(&base_speed, 0.0, 6.3, 100_000));
    let law = |w: Scalar| {
        if w <= 2.0 {
            (0.5, 0.0)
        } else if w <= 20.0 {
            (0.5 + 0.25 * (w - 2.0) / 18.0, 0.25 / 18.0)
        } else {
            (0.75, 0.0)
        }
    };
    let lengths = base.clone();
    Plane {
        point: Box::new(move |t| {
            let (d, _) = law(lengths.at(&base_speed, t));
            let v = e(t);
            Point2::new(30.0 + 10.0 * t.cos(), 6.0 * t.sin())
                + d * Vec2::new(-v.y, v.x) / v.length()
        }),
        velocity: Box::new(move |t| {
            let (v, dv) = (e(t), Vec2::new(-10.0 * t.cos(), -6.0 * t.sin()));
            let g = v.length();
            let (d, rate) = law(base.at(&base_speed, t));
            let n = Vec2::new(-v.y, v.x);
            let dn = Vec2::new(-dv.y, dv.x);
            let turn = dn / g - n * (v.dot(dv) / (g * g * g));
            v + rate * n + d * turn
        }),
    }
}

fn deck(b: &mut GeometryGraphBuilder, at: [Scalar; 2]) -> (NodeId, Vec<SectionAtStation>) {
    let directrix = deck_directrix(b);
    let profile = rectangle(b, 2.0, 1.0);
    let sections: Vec<_> = at
        .iter()
        .map(|&d| SectionAtStation::new(profile, Station::at(d)))
        .collect();
    (spine(b, directrix, sections.clone()), sections)
}

#[test]
fn the_ifc_443_deck_is_framed_and_names_its_mitre() {
    let mut b = GeometryGraphBuilder::new();
    let (root, _) = deck(&mut b, [0.0, 5.0]);
    let graph = b.finish(vec![root]).unwrap();
    // Its frames are certified now: what remains is the cell beside the
    // mitre where the law starts to widen, at station 2, a corner of the
    // offset on a curved piece.
    let reason = unbounded(&graph, root);
    assert!(reason.contains("mitred seam"), "{reason}");
    let mesh = run_bound(&graph, root).0;
    let volume = volume_properties(&mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    assert!((volume - 10.0).abs() <= 1e-3, "{volume}");
}

#[test]
fn a_run_along_the_decks_widening_offset_is_certified() {
    for at in [[3.0, 12.0], [0.2, 1.5], [21.0, 30.0]] {
        let mut b = GeometryGraphBuilder::new();
        let (root, sections) = deck(&mut b, at);
        let graph = b.finish(vec![root]).unwrap();
        let exact = Exact::new(deck_curve(), (0.0, 4.0), &graph, &sections);
        eprintln!("run {at:?}");
        let (_, bound) = assert_certified(&graph, root, &exact);
        eprintln!("run {at:?}: {bound}");
    }
}

// --- other bases read by quadrature ---------------------------------------------

#[test]
fn a_taper_along_a_b_spline_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(spline()).unwrap();
    let (small, big) = (rectangle(&mut b, 1.0, 1.0), rectangle(&mut b, 2.0, 1.5));
    let sections = vec![
        SectionAtStation::new(small, Station::new(1.0, StationOffsets::ZERO)),
        SectionAtStation::new(big, Station::new(12.0, StationOffsets::new(0.4, -0.3, 0.1))),
        SectionAtStation::new(small, Station::new(22.0, StationOffsets::ZERO)),
    ];
    let root = spine(&mut b, basis, sections.clone());
    let graph = b.finish(vec![root]).unwrap();
    let exact = Exact::new(Plane::of(spline()), (0.0, 1.0), &graph, &sections);
    assert_certified(&graph, root, &exact);
}

#[test]
fn a_run_along_an_ellipse_and_along_its_offset_is_certified() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(ellipse()).unwrap();
    let beside = b
        .push(GeometryNode::CurveRelation(CurveRelation::Offset {
            basis,
            distance: 0.5,
            reference_direction: None,
        }))
        .unwrap();
    let profile = rectangle(&mut b, 2.0, 1.0);
    let sections = vec![
        SectionAtStation::new(profile, Station::at(4.0)),
        SectionAtStation::new(profile, Station::at(20.0)),
    ];
    let along = spine(&mut b, basis, sections.clone());
    let offset = spine(&mut b, beside, sections.clone());
    let graph = b.finish(vec![along, offset]).unwrap();
    let exact = Exact::new(Plane::of(ellipse()), (0.0, 3.0), &graph, &sections);
    assert_certified(&graph, along, &exact);
    // The constant offset: the ellipse's point moved along its normal.
    let curve = Plane {
        point: Box::new(|t| {
            let v = Vec2::new(-10.0 * t.sin(), 6.0 * t.cos());
            Point2::new(30.0 + 10.0 * t.cos(), 6.0 * t.sin())
                + 0.5 * Vec2::new(-v.y, v.x) / v.length()
        }),
        velocity: Box::new(|t| {
            let (v, dv) = (
                Vec2::new(-10.0 * t.sin(), 6.0 * t.cos()),
                Vec2::new(-10.0 * t.cos(), -6.0 * t.sin()),
            );
            let g = v.length();
            let (n, dn) = (Vec2::new(-v.y, v.x), Vec2::new(-dv.y, dv.x));
            v + 0.5 * (dn / g - n * (v.dot(dv) / (g * g * g)))
        }),
    };
    let exact = Exact::new(curve, (0.0, 3.0), &graph, &sections);
    assert_certified(&graph, offset, &exact);
}

#[test]
fn intrinsic_curves_and_their_offsets_stay_unbounded_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let clothoid = b
        .push_value(Curve2::Intrinsic(Intrinsic2::new(
            Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0),
            60.0,
        )))
        .unwrap();
    let beside = b
        .push(GeometryNode::CurveRelation(CurveRelation::Offset {
            basis: clothoid,
            distance: 2.0,
            reference_direction: None,
        }))
        .unwrap();
    let profile = rectangle(&mut b, 1.0, 1.0);
    let sections = vec![
        SectionAtStation::new(profile, Station::at(5.0)),
        SectionAtStation::new(profile, Station::at(40.0)),
    ];
    let along = spine(&mut b, clothoid, sections.clone());
    let offset = spine(&mut b, beside, sections);
    let graph = b.finish(vec![along, offset]).unwrap();
    for root in [along, offset] {
        let reason = unbounded(&graph, root);
        assert!(reason.contains("certified derivative bounds"), "{reason}");
        assert!(reason.contains("intrinsic curve"), "{reason}");
    }
}
