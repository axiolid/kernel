//! Certified bounds of runs of sections at stations whose frames are read
//! by quadrature (#316): runs along the openbimrs/ifc#443 deck's offset by
//! distances beside an ellipse, along an atomic ellipse, a tapered run
//! along a B-spline and a run along a constant offset of an ellipse each
//! report a certified deviation within the request, and a dense sampling
//! of the EXACT solid -- its frames computed here from the curves in
//! closed form at their true arc length, by a fine quadrature independent
//! of the kernel's -- never strays past it. So do the deck's cell beside
//! the mitre where its law starts to widen and a B-spline's across its
//! corner knot (#335), against the exact mitred solid of two curves in
//! closed form. The bases that stay uncertified are still unbounded by
//! name.

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

fn options_at(chord: Scalar) -> ExecutionOptions {
    options().with_chord_error(chord).unwrap()
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
        .as_chunks::<3>()
        .0
        .iter()
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
    run_bound_with(graph, root, &options())
}

fn run_bound_with(
    graph: &GeometryGraph,
    root: NodeId,
    options: &ExecutionOptions,
) -> (TriMesh, DeviationBound, bool) {
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(graph, root, options)
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

// --- the exact solid across a mitre ---------------------------------------------

/// One side of a mitre in closed form: its curve, its arc lengths from
/// the curve's parameter `0`, and where the run's knot stands on it.
struct Side {
    curve: Plane,
    lengths: Lengths,
    /// The parameter of the knot on this curve.
    knot: Scalar,
}

impl Side {
    fn new(curve: Plane, end: Scalar, knot: Scalar) -> Self {
        let speed = |t: Scalar| (curve.velocity)(t).length();
        let lengths = Lengths::new(&speed, 0.0, end, 100_000);
        Self {
            curve,
            lengths,
            knot,
        }
    }

    /// The run's distance of the knot along this curve.
    fn knot_length(&self) -> Scalar {
        let speed = |t: Scalar| (self.curve.velocity)(t).length();
        self.lengths.at(&speed, self.knot)
    }

    /// The point and the unit tangent and lateral at `s`, the run's
    /// distance, whose knot stands at `knot_at`.
    fn frame(&self, s: Scalar, knot_at: Scalar) -> (Point3, Vec3, Vec3) {
        let speed = |t: Scalar| (self.curve.velocity)(t).length();
        let t = self
            .lengths
            .inverse(&speed, self.knot_length() + (s - knot_at));
        let (p, v) = ((self.curve.point)(t), (self.curve.velocity)(t));
        let tangent = Vec3::new(v.x, v.y, 0.0) / v.length();
        (
            Point3::new(p.x, p.y, 0.0),
            tangent,
            Vec3::new(-tangent.y, tangent.x, 0.0),
        )
    }
}

/// The exact mitred solid of a constant rectangle (`2 hw` by `2 hh`) along
/// two curves meeting in a corner at the run's distance `knot_at`: each
/// side's walls kept behind the mitre plane, extruded to it along its end
/// tangent, and the faces between the two cuts.
struct Mitred {
    incoming: Side,
    outgoing: Side,
    knot_at: Scalar,
    at: [Scalar; 2],
    half: (Scalar, Scalar),
}

const RECTANGLE: [(Scalar, Scalar); 4] = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];

fn nearest(triangles: &[[Point3; 3]], p: Point3) -> Scalar {
    triangles
        .iter()
        .filter_map(|&t| axiolid_measure::proximity::closest_point_on_triangle(p, t).ok())
        .map(|q| p.distance(q))
        .fold(Scalar::INFINITY, Scalar::min)
}

impl Mitred {
    fn point(&self, side: &Side, s: Scalar, (x, y): (Scalar, Scalar)) -> Point3 {
        let (origin, _, lateral) = side.frame(s, self.knot_at);
        origin + lateral * (x * self.half.0) + Vec3::Z * (y * self.half.1)
    }

    /// The largest distance of the exact solid from `mesh`: walls with
    /// `steps` distances per side and `per_edge` points per edge, then the
    /// faces.
    fn measure(&self, mesh: &TriMesh, steps: usize, per_edge: usize) -> Scalar {
        let triangles: Vec<[Point3; 3]> = mesh
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .map(|t| [0, 1, 2].map(|k| mesh.positions[t[k] as usize]))
            .collect();
        let (origin, t_in, _) = self.incoming.frame(self.knot_at, self.knot_at);
        let (_, t_out, _) = self.outgoing.frame(self.knot_at, self.knot_at);
        let normal = (t_in + t_out).normalize();
        let height = |p: Point3| (p - origin).dot(normal);
        let along = |a: (Scalar, Scalar), b: (Scalar, Scalar), v: Scalar| {
            (a.0 + (b.0 - a.0) * v, a.1 + (b.1 - a.1) * v)
        };
        let mut worst: Scalar = 0.0;
        for e in 0..4 {
            let (a, b) = (RECTANGLE[e], RECTANGLE[(e + 1) % 4]);
            for iv in 0..=per_edge {
                let xy = along(a, b, iv as Scalar / per_edge as Scalar);
                for j in 0..=steps {
                    let u = j as Scalar / steps as Scalar;
                    let s = self.at[0] + (self.knot_at - self.at[0]) * u;
                    let p = self.point(&self.incoming, s, xy);
                    if height(p) <= 0.0 {
                        worst = worst.max(nearest(&triangles, p));
                    }
                    let s = self.knot_at + (self.at[1] - self.knot_at) * u;
                    let p = self.point(&self.outgoing, s, xy);
                    if height(p) >= 0.0 {
                        worst = worst.max(nearest(&triangles, p));
                    }
                }
                let end = self.point(&self.incoming, self.knot_at, xy);
                if height(end) < 0.0 {
                    let reach = -height(end) / t_in.dot(normal);
                    for j in 0..=40 {
                        let p = end + t_in * (reach * j as Scalar / 40.0);
                        worst = worst.max(nearest(&triangles, p));
                    }
                }
                let start = self.point(&self.outgoing, self.knot_at, xy);
                if height(start) > 0.0 {
                    let reach = height(start) / t_out.dot(normal);
                    for j in 0..=40 {
                        let p = start - t_out * (reach * j as Scalar / 40.0);
                        worst = worst.max(nearest(&triangles, p));
                    }
                }
            }
        }
        // The faces: the plane's points inside one side's cut ring only.
        let across = Vec3::new(-normal.y, normal.x, 0.0);
        let to_plane = |p: Point3| ((p - origin).dot(across), (p - origin).dot(Vec3::Z));
        let cut = |ending: bool, xy: (Scalar, Scalar)| -> Point3 {
            let (side, tangent, sign, range) = if ending {
                (&self.incoming, t_in, 1.0, (self.at[0], self.knot_at))
            } else {
                (&self.outgoing, t_out, -1.0, (self.at[1], self.knot_at))
            };
            let at_knot = self.point(side, self.knot_at, xy);
            let h = height(at_knot);
            if h * sign <= 0.0 {
                return at_knot - tangent * (h / tangent.dot(normal));
            }
            // The wall crosses the plane before the knot.
            let (mut near, mut far) = range;
            for _ in 0..100 {
                let mid = 0.5 * (near + far);
                if height(self.point(side, mid, xy)) * sign > 0.0 {
                    far = mid;
                } else {
                    near = mid;
                }
            }
            self.point(side, 0.5 * (near + far), xy)
        };
        let ring = |ending: bool| -> Vec<(Scalar, Scalar)> {
            let mut out = Vec::new();
            for e in 0..4 {
                let (a, b) = (RECTANGLE[e], RECTANGLE[(e + 1) % 4]);
                for q in 0..64 {
                    out.push(to_plane(cut(ending, along(a, b, q as Scalar / 64.0))));
                }
            }
            out
        };
        let (ring_in, ring_out) = (ring(true), ring(false));
        let inside = |ring: &[(Scalar, Scalar)], q: (Scalar, Scalar)| {
            let mut in_ring = false;
            for (i, a) in ring.iter().enumerate() {
                let b = ring[(i + 1) % ring.len()];
                if (a.1 > q.1) != (b.1 > q.1) && q.0 < (b.0 - a.0) * (q.1 - a.1) / (b.1 - a.1) + a.0
                {
                    in_ring = !in_ring;
                }
            }
            in_ring
        };
        let (mut x0, mut x1, mut y0, mut y1) = (Scalar::MAX, Scalar::MIN, Scalar::MAX, Scalar::MIN);
        for p in ring_in.iter().chain(&ring_out) {
            x0 = x0.min(p.0);
            x1 = x1.max(p.0);
            y0 = y0.min(p.1);
            y1 = y1.max(p.1);
        }
        for i in 0..=400 {
            for j in 0..=40 {
                let q = (
                    x0 + (x1 - x0) * i as Scalar / 400.0,
                    y0 + (y1 - y0) * j as Scalar / 40.0,
                );
                if inside(&ring_in, q) != inside(&ring_out, q) {
                    let p = origin + across * q.0 + Vec3::Z * q.1;
                    worst = worst.max(nearest(&triangles, p));
                }
            }
        }
        worst
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

/// The deck's ellipse: its velocity and speed at the parameter `t`.
fn deck_ellipse_velocity(t: Scalar) -> Vec2 {
    Vec2::new(-10.0 * t.sin(), 6.0 * t.cos())
}

fn deck_ellipse_speed(t: Scalar) -> Scalar {
    deck_ellipse_velocity(t).length()
}

/// The ellipse's arc length from the parameter `0`.
fn deck_base() -> std::rc::Rc<Lengths> {
    std::rc::Rc::new(Lengths::new(&deck_ellipse_speed, 0.0, 6.3, 100_000))
}

/// The offset law along the ellipse's arc length `w`: the distance beside
/// it and its rate, constant to station 2, widening to 20, held after.
fn deck_law(w: Scalar) -> (Scalar, Scalar) {
    if w <= 2.0 {
        (0.5, 0.0)
    } else if w <= 20.0 {
        (0.5 + 0.25 * (w - 2.0) / 18.0, 0.25 / 18.0)
    } else {
        (0.75, 0.0)
    }
}

/// The deck's directrix exactly: the law runs linearly in the ellipse's
/// arc length, its point beside the ellipse's along the left normal.
fn deck_curve() -> Plane {
    deck_plane(deck_law)
}

/// The ellipse offset by `law` of its arc length.
fn deck_plane(law: fn(Scalar) -> (Scalar, Scalar)) -> Plane {
    let e = deck_ellipse_velocity;
    let base = deck_base();
    let lengths = base.clone();
    Plane {
        point: Box::new(move |t| {
            let (d, _) = law(lengths.at(&deck_ellipse_speed, t));
            let v = e(t);
            Point2::new(30.0 + 10.0 * t.cos(), 6.0 * t.sin())
                + d * Vec2::new(-v.y, v.x) / v.length()
        }),
        velocity: Box::new(move |t| {
            let (v, dv) = (e(t), Vec2::new(-10.0 * t.cos(), -6.0 * t.sin()));
            let g = v.length();
            let (d, rate) = law(base.at(&deck_ellipse_speed, t));
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
fn the_ifc_443_deck_is_certified_across_its_mitre() {
    let mut b = GeometryGraphBuilder::new();
    let (root, _) = deck(&mut b, [0.0, 5.0]);
    let graph = b.finish(vec![root]).unwrap();
    let (mesh, bound, meets) = run_bound(&graph, root);
    let DeviationBound::Certified(bound) = bound else {
        panic!("not certified: {bound:?}")
    };
    // Its frames are certified, and so is the cell beside the mitre where
    // the law starts to widen, at station 2 (#335).
    assert!((bound - 4.88e-4).abs() <= 1e-6, "{bound}");
    assert!(meets, "{bound} is past the request");
    let volume = volume_properties(&mesh, Tolerance::MILLIMETRE)
        .unwrap()
        .signed_volume;
    assert!((volume - 10.0).abs() <= 1e-3, "{volume}");
    // The mitre's exact solid: the offset before the corner holds its
    // distance, the one after widens; the corner is at the offset's own
    // length where the ellipse's arc length reaches 2.
    let knot = deck_base().inverse(&deck_ellipse_speed, 2.0);
    let incoming = Side::new(deck_plane(|_| (0.5, 0.0)), 6.3, knot);
    let outgoing = Side::new(
        deck_plane(|w| (0.5 + 0.25 * (w - 2.0) / 18.0, 0.25 / 18.0)),
        6.3,
        knot,
    );
    let knot_at = incoming.knot_length();
    let far = Mitred {
        incoming,
        outgoing,
        knot_at,
        at: [0.0, 5.0],
        half: (1.0, 0.5),
    }
    .measure(&mesh, 200, 8);
    assert!(
        far <= bound,
        "the exact solid strays {far}, past the bound {bound}"
    );
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

/// A quadratic plan B-spline whose curvature jumps at its interior knot.
fn quadratic_spline() -> Curve2 {
    Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(12.0, 0.0),
            Point2::new(12.0, 12.0),
            Point2::new(24.0, 14.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 1, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// A run of a 1 by 1 rectangle between the distances `at`, `lateral` to the
/// side of `curve` or of its offset by `beside`; the run and its sections.
fn broken_run(
    curve: Curve2,
    beside: Option<Scalar>,
    lateral: Scalar,
    at: [Scalar; 2],
) -> (GeometryGraph, NodeId, Vec<SectionAtStation>) {
    let mut b = GeometryGraphBuilder::new();
    let mut basis = b.push_value(curve).unwrap();
    if let Some(distance) = beside {
        basis = b
            .push(GeometryNode::CurveRelation(CurveRelation::Offset {
                basis,
                distance,
                reference_direction: None,
            }))
            .unwrap();
    }
    let profile = rectangle(&mut b, 1.0, 1.0);
    let offsets = StationOffsets::new(0.0, lateral, 0.0);
    let sections = vec![
        SectionAtStation::new(profile, Station::new(at[0], offsets)),
        SectionAtStation::new(profile, Station::new(at[1], offsets)),
    ];
    let root = spine(&mut b, basis, sections.clone());
    (b.finish(vec![root]).unwrap(), root, sections)
}

/// Runs far to the side of a curve across a point where its frame's
/// curvature jumps stay within their certified deviation: the bound is
/// not split there, yet a dense sampling of the exact solid never strays
/// past it (checked before #406 decided folds by these frames' rates).
fn assert_across_a_break(curve: fn() -> Curve2, laterals: [Scalar; 2], exact: fn() -> Plane) {
    for lateral in laterals {
        for chord in [1.0, 0.05] {
            for at in [[1.0, 22.0], [8.0, 18.0]] {
                let (graph, root, sections) = broken_run(curve(), None, lateral, at);
                let (mesh, bound, _) = run_bound_with(&graph, root, &options_at(chord));
                let DeviationBound::Certified(bound) = bound else {
                    panic!("not certified (lateral {lateral}, chord {chord}, stations {at:?}): {bound:?}")
                };
                let exact = Exact::new(exact(), (0.0, 1.0), &graph, &sections);
                let far = measured(&mesh, &exact.points(400, 4));
                assert!(
                    far <= bound,
                    "lateral {lateral}, chord {chord}, stations {at:?}: the exact solid \
                     strays {far}, past {bound}"
                );
            }
        }
    }
}

#[test]
fn a_run_across_a_knot_of_a_quadratic_b_spline_is_certified() {
    assert_across_a_break(quadratic_spline, [5.0, -9.0], || {
        Plane::of(quadratic_spline())
    });
}

/// The cubic spline's offset by 1 along its normal, in closed form: the
/// velocity is a central difference of the point. The base's third
/// derivative enters the offset's curvature, which jumps at the base's knot.
fn cubic_offset() -> Plane {
    fn point(t: Scalar) -> Point2 {
        let base = evaluate2(&spline(), t).unwrap();
        let v = derivative2(&spline(), t).unwrap();
        base + Vec2::new(-v.y, v.x) / v.length()
    }
    Plane {
        point: Box::new(point),
        velocity: Box::new(|t| {
            let h = 1e-6;
            (point(t + h) - point(t - h)) / (2.0 * h)
        }),
    }
}

#[test]
fn a_run_across_a_knot_beside_a_cubic_offset_is_certified_or_unbounded_by_name() {
    for chord in [1.0, 0.05] {
        let (graph, root, sections) = broken_run(spline(), Some(1.0), 1.0, [1.0, 22.0]);
        let (mesh, bound, _) = run_bound_with(&graph, root, &options_at(chord));
        match bound {
            DeviationBound::Certified(bound) => {
                let exact = Exact::new(cubic_offset(), (0.0, 1.0), &graph, &sections);
                let far = measured(&mesh, &exact.points(400, 4));
                assert!(
                    far <= bound,
                    "chord {chord}: the exact solid strays {far}, past {bound}"
                );
            }
            DeviationBound::Unbounded(reason) => {
                assert!(
                    reason.contains("no certified derivative bounds"),
                    "{reason}"
                );
            }
            other => panic!("{other:?}"),
        }
    }
}

/// A quadratic Bezier arc as a plane curve in closed form.
fn bezier(p: [Point2; 3]) -> Plane {
    Plane {
        point: Box::new(move |t| {
            let (a, b, c) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t), t * t);
            Point2::new(
                a * p[0].x + b * p[1].x + c * p[2].x,
                a * p[0].y + b * p[1].y + c * p[2].y,
            )
        }),
        velocity: Box::new(move |t| 2.0 * (1.0 - t) * (p[1] - p[0]) + 2.0 * t * (p[2] - p[1])),
    }
}

#[test]
fn a_b_spline_with_a_corner_knot_is_certified_across_its_mitre() {
    // Two quadratic Beziers joined at (10, 0) with a corner: the spline
    // of degree two with an interior knot of multiplicity two.
    let spline = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 1.0),
            Point2::new(10.0, 0.0),
            Point2::new(15.0, 0.5),
            Point2::new(20.0, 2.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 2, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push_value(spline).unwrap();
    let profile = rectangle(&mut b, 2.0, 1.0);
    let at = [1.0, 19.0];
    let sections: Vec<_> = at
        .iter()
        .map(|&d| SectionAtStation::new(profile, Station::at(d)))
        .collect();
    let root = spine(&mut b, basis, sections);
    let graph = b.finish(vec![root]).unwrap();
    let incoming = Side::new(
        bezier([
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 1.0),
            Point2::new(10.0, 0.0),
        ]),
        1.0,
        1.0,
    );
    let outgoing = Side::new(
        bezier([
            Point2::new(10.0, 0.0),
            Point2::new(15.0, 0.5),
            Point2::new(20.0, 2.0),
        ]),
        1.0,
        0.0,
    );
    let knot_at = incoming.knot_length();
    let exact = Mitred {
        incoming,
        outgoing,
        knot_at,
        at,
        half: (1.0, 0.5),
    };
    for (chord, expected, meets) in [(0.01, 3.38e-3, true), (1e-3, 1.72e-3, false)] {
        let (mesh, bound, met) = run_bound_with(&graph, root, &options_at(chord));
        let DeviationBound::Certified(bound) = bound else {
            panic!("not certified at {chord}: {bound:?}")
        };
        assert!((bound - expected).abs() <= 1e-5, "{bound} at {chord}");
        assert_eq!(met, meets, "{bound} against {chord}");
        let far = exact.measure(&mesh, 200, 8);
        assert!(
            far <= bound,
            "the exact solid strays {far}, past the bound {bound}"
        );
    }
}

#[test]
fn an_offset_beside_an_intrinsic_curve_stays_unbounded_by_name() {
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
    // The clothoid itself is certified since #361 (`station_bound_clothoid`).
    let offset = spine(&mut b, beside, sections);
    let graph = b.finish(vec![offset]).unwrap();
    let reason = unbounded(&graph, offset);
    assert!(reason.contains("certified derivative bounds"), "{reason}");
    assert!(reason.contains("intrinsic curve"), "{reason}");
}
