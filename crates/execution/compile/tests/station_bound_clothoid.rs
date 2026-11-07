//! Certified bounds of runs of sections at stations along clothoid plans
//! (#361): a rectangular deck widening from 10 m to 12 m with its axis
//! drifting 1.5 m sideways, along a gradient curve over a plan of a line,
//! a clothoid, an arc, a clothoid and a line, with a parabolic sag and a
//! circular crest, in the section and the plan frame, across the plan's
//! seams, at a tight budget, and along a bare clothoid. Each run reports a
//! certified deviation within the request, and a dense sampling of the
//! EXACT solid never strays past it.
//!
//! The exact solid is independent of the kernel's quadratures: the plan's
//! points from the Fresnel integrals' power series and closed forms, its
//! headings and the profile's heights and grades in closed form, and the
//! frames built here as the section frames are defined.

use axiolid_construct::profile::profile_rings;
use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, CurvatureLaw, Curve2, Curve3,
    Elevated3, ElevationLaw, Intrinsic2,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, ReferenceMeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, SectionAtStation, SolidOperation,
    Station, StationFrame, StationOffsets,
};
use axiolid_profile::{Profile, RectangleProfile};

/// The chord budget of the runs, in metres.
const BUDGET: Scalar = 0.01;

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(budget)
        .unwrap()
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

// --- the exact curve, independently ---------------------------------------------

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(3.0, -1.0),
        x: Vec2::new(0.6, 0.8),
        y: Vec2::new(-0.8, 0.6),
    }
}

/// `int_0^s (cos(a u^2), sin(a u^2)) du` by its power series.
fn fresnel(a: Scalar, s: Scalar) -> Vec2 {
    let (mut x, mut y) = (0.0, 0.0);
    let mut factorial = 1.0;
    for n in 0..40 {
        if n > 0 {
            factorial *= n as Scalar;
        }
        let value = a.powi(n) * s.powi(2 * n + 1) / (factorial * (2 * n + 1) as Scalar);
        match n % 4 {
            0 => x += value,
            1 => y += value,
            2 => x -= value,
            _ => y -= value,
        }
    }
    Vec2::new(x, y)
}

fn rotate(angle: Scalar, v: Vec2) -> Vec2 {
    let (s, c) = angle.sin_cos();
    Vec2::new(c * v.x - s * v.y, s * v.x + c * v.y)
}

#[derive(Clone, Copy)]
enum Piece {
    Line(Scalar),
    Rising(Scalar, Scalar),
    Arc(Scalar, Scalar),
    Falling(Scalar, Scalar),
}

impl Piece {
    fn length(self) -> Scalar {
        match self {
            Piece::Line(l) | Piece::Rising(l, _) | Piece::Arc(l, _) | Piece::Falling(l, _) => l,
        }
    }

    fn law(self) -> CurvatureLaw {
        match self {
            Piece::Line(_) => CurvatureLaw::Constant { curvature: 0.0 },
            Piece::Rising(l, r) => CurvatureLaw::clothoid(0.0, 1.0 / r, l),
            Piece::Arc(_, r) => CurvatureLaw::Constant { curvature: 1.0 / r },
            Piece::Falling(l, r) => CurvatureLaw::clothoid(1.0 / r, 0.0, l),
        }
    }

    fn turn(self) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Rising(l, r) | Piece::Falling(l, r) => 0.5 * l / r,
            Piece::Arc(l, r) => l / r,
        }
    }

    fn heading(self, s: Scalar) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Rising(l, r) => 0.5 * s * s / (r * l),
            Piece::Arc(_, r) => s / r,
            Piece::Falling(l, r) => self.turn() - 0.5 * (l - s) * (l - s) / (r * l),
        }
    }

    fn displacement(self, s: Scalar) -> Vec2 {
        match self {
            Piece::Line(_) => Vec2::new(s, 0.0),
            Piece::Rising(l, r) => fresnel(0.5 / (r * l), s),
            Piece::Arc(_, r) => Vec2::new(r * (s / r).sin(), r * (1.0 - (s / r).cos())),
            Piece::Falling(l, r) => {
                let a = 0.5 / (r * l);
                let tail = fresnel(a, l) - fresnel(a, l - s);
                rotate(self.turn(), Vec2::new(tail.x, -tail.y))
            }
        }
    }
}

struct Alignment {
    curve: Intrinsic2,
    pieces: Vec<Piece>,
}

impl Alignment {
    fn new(pieces: Vec<Piece>) -> Self {
        let mut breaks = Vec::new();
        let mut length = 0.0;
        for piece in &pieces[..pieces.len() - 1] {
            length += piece.length();
            breaks.push(length);
        }
        length += pieces.last().unwrap().length();
        let law = CurvatureLaw::piecewise(breaks, pieces.iter().map(|p| p.law()).collect());
        Self {
            curve: Intrinsic2::new(frame(), law, length),
            pieces,
        }
    }

    /// The exact point and unit tangent at `s`, in the world.
    fn world(&self, s: Scalar) -> (Vec2, Vec2) {
        let (mut at, mut heading, mut begin) = (Vec2::ZERO, 0.0, 0.0);
        for piece in &self.pieces {
            let inside = (s - begin).clamp(0.0, piece.length());
            at += rotate(heading, piece.displacement(inside));
            if s <= begin + piece.length() {
                heading += piece.heading(inside);
                break;
            }
            heading += piece.turn();
            begin += piece.length();
        }
        let f = frame();
        (
            f.origin + f.x * at.x + f.y * at.y,
            f.x * heading.cos() + f.y * heading.sin(),
        )
    }
}

/// Line 40, clothoid 0 to 1/300 over 60, arc 80, clothoid back 60, line 40.
fn alignment() -> Alignment {
    Alignment::new(vec![
        Piece::Line(40.0),
        Piece::Rising(60.0, 300.0),
        Piece::Arc(80.0, 300.0),
        Piece::Falling(60.0, 300.0),
        Piece::Line(40.0),
    ])
}

/// A parabolic sag over `[0, 120]` into a circular crest of radius 4000.
fn profile() -> ElevationLaw {
    ElevationLaw::Piecewise {
        breaks: vec![120.0],
        laws: vec![
            ElevationLaw::parabolic(10.0, 0.01, 0.0148, 120.0),
            ElevationLaw::circular_arc(11.488, 0.0148, -4000.0),
        ],
    }
}

fn profile_at(d: Scalar) -> (Scalar, Scalar) {
    if d <= 120.0 {
        return (10.0 + 0.01 * d + 2e-5 * d * d, 0.01 + 4e-5 * d);
    }
    let (r, t0) = (-4000.0_f64, 0.0148_f64.atan());
    let t = (t0.sin() + (d - 120.0) / r).asin();
    (11.488 + r * (t0.cos() - t.cos()), t.tan())
}

/// An exact section frame at a plan distance.
#[derive(Clone, Copy)]
struct Reference {
    point: Vec3,
    tangent: Vec3,
    lateral: Vec3,
    up: Vec3,
}

fn reference_up(point: Vec3, tangent: Vec3) -> Reference {
    let right = tangent.cross(Vec3::Z);
    let right = right / right.length();
    Reference {
        point,
        tangent,
        lateral: -right,
        up: right.cross(tangent),
    }
}

fn plan_of(r: Reference) -> Reference {
    let t = Vec3::new(r.tangent.x, r.tangent.y, 0.0);
    let t = t / t.length();
    Reference {
        point: r.point,
        tangent: t,
        lateral: Vec3::new(-t.y, t.x, 0.0),
        up: Vec3::Z,
    }
}

fn elevated_at(a: &Alignment, d: Scalar) -> Reference {
    let (p, t) = a.world(d);
    let (z, g) = profile_at(d);
    let size = (1.0 + g * g).sqrt();
    reference_up(
        Vec3::new(p.x, p.y, z),
        Vec3::new(t.x / size, t.y / size, g / size),
    )
}

// --- the exact solid -------------------------------------------------------------

struct Oracle {
    frame_at: Box<dyn Fn(Scalar) -> Reference>,
    distances: Vec<Scalar>,
    locals: Vec<Vec<Point2>>,
    offsets: Vec<StationOffsets>,
}

fn oracle(
    graph: &GeometryGraph,
    sections: &[SectionAtStation],
    frame_at: Box<dyn Fn(Scalar) -> Reference>,
    budget: Scalar,
) -> Oracle {
    Oracle {
        frame_at,
        distances: sections.iter().map(|s| s.station.distance).collect(),
        locals: sections
            .iter()
            .map(|s| {
                let Some(GeometryNode::Profile(p)) = graph.get(s.profile) else {
                    panic!("a profile")
                };
                profile_rings(p, 0.5 * budget, Tolerance::MILLIMETRE)
                    .unwrap()
                    .outer
            })
            .collect(),
        offsets: sections.iter().map(|s| s.station.offsets).collect(),
    }
}

/// Points of the exact wall: every edge of the sections, interpolated
/// linearly in distance, at `per_interval` + 1 distances per station
/// interval and `per_edge` + 1 points along it.
fn exact_points(oracle: &Oracle, per_interval: usize, per_edge: usize) -> Vec<Point3> {
    let mut out = Vec::new();
    for i in 0..oracle.distances.len() - 1 {
        let (d0, d1) = (oracle.distances[i], oracle.distances[i + 1]);
        for a in 0..=per_interval {
            let u = a as Scalar / per_interval as Scalar;
            let f = (oracle.frame_at)(d0 + (d1 - d0) * u);
            let (p, q) = (oracle.offsets[i], oracle.offsets[i + 1]);
            let lerp = |x: Scalar, y: Scalar| x + (y - x) * u;
            let origin = f.point
                + lerp(p.lateral, q.lateral) * f.lateral
                + lerp(p.vertical, q.vertical) * f.up
                + lerp(p.longitudinal, q.longitudinal) * f.tangent;
            let ring: Vec<Point3> = oracle.locals[i]
                .iter()
                .zip(&oracle.locals[i + 1])
                .map(|(a, b)| {
                    Point3::from(origin + lerp(a.x, b.x) * f.lateral + lerp(a.y, b.y) * f.up)
                })
                .collect();
            for k in 0..ring.len() {
                let (v, w) = (ring[k], ring[(k + 1) % ring.len()]);
                for c in 0..=per_edge {
                    out.push(v + (w - v) * (c as Scalar / per_edge as Scalar));
                }
            }
        }
    }
    out
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

fn run_bound(
    graph: &GeometryGraph,
    root: NodeId,
    budget: Scalar,
) -> (TriMesh, DeviationBound, bool) {
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(graph, root, &options(budget))
        .unwrap();
    let bound = report
        .contributions
        .iter()
        .find(|c| c.path == DeviationPath::StationedSpine)
        .expect("the run reports its path")
        .bound;
    (outcome.mesh, bound, report.meets_requested())
}

/// The run is certified within the request, and a dense sampling of the
/// exact wall stays within the bound; the bound.
fn assert_certified(graph: &GeometryGraph, root: NodeId, budget: Scalar, exact: &Oracle) -> Scalar {
    let (mesh, bound, meets) = run_bound(graph, root, budget);
    let DeviationBound::Certified(bound) = bound else {
        panic!("not certified: {bound:?}")
    };
    assert!(meets, "{bound} is past the request {budget}");
    let far = measured(&mesh, &exact_points(exact, 24, 2));
    assert!(
        far <= bound + 1e-10,
        "the exact solid strays {far} from the mesh, past the bound {bound}"
    );
    bound
}

fn spine(
    b: &mut GeometryGraphBuilder,
    directrix: NodeId,
    sections: Vec<SectionAtStation>,
    frame: StationFrame,
) -> NodeId {
    b.push(GeometryNode::SolidOperation(
        SolidOperation::SectionsAtStations {
            directrix,
            sections,
            frame,
        },
    ))
    .unwrap()
}

/// A deck 10 m wide widening to 12 m, its axis drifting from 0 to 1.5 m,
/// at `stations`, alternating between the two.
fn deck_run(
    stations: &[Scalar],
    frame: StationFrame,
    budget: Scalar,
) -> (GeometryGraph, NodeId, Oracle) {
    let a = alignment();
    let mut b = GeometryGraphBuilder::new();
    let curve = Curve3::Elevated(Elevated3::new(
        Curve2::Intrinsic(a.curve.clone()),
        profile(),
    ));
    let basis = b.push_value(curve).unwrap();
    let narrow = rectangle(&mut b, 10.0, 0.4);
    let wide = rectangle(&mut b, 12.0, 0.4);
    let sections: Vec<SectionAtStation> = stations
        .iter()
        .enumerate()
        .map(|(i, &d)| {
            let (profile, drift) = if i % 2 == 0 {
                (narrow, 0.0)
            } else {
                (wide, 1.5)
            };
            SectionAtStation::new(
                profile,
                Station::new(d, StationOffsets::new(drift, 0.2, 0.0)),
            )
        })
        .collect();
    let root = spine(&mut b, basis, sections.clone(), frame);
    let graph = b.finish(vec![root]).unwrap();
    let oracle = oracle(
        &graph,
        &sections,
        match frame {
            StationFrame::Plan => Box::new(move |d| plan_of(elevated_at(&a, d))),
            _ => Box::new(move |d| elevated_at(&a, d)),
        },
        budget,
    );
    (graph, root, oracle)
}

#[test]
fn a_deck_along_the_gradient_curve_is_certified() {
    let (graph, root, exact) = deck_run(
        &[10.0, 75.0, 130.0, 205.0, 270.0],
        StationFrame::Section,
        BUDGET,
    );
    assert_certified(&graph, root, BUDGET, &exact);
}

#[test]
fn a_deck_in_the_plan_frame_is_certified() {
    let (graph, root, exact) = deck_run(
        &[10.0, 75.0, 130.0, 205.0, 270.0],
        StationFrame::Plan,
        BUDGET,
    );
    assert_certified(&graph, root, BUDGET, &exact);
}

#[test]
fn stations_straddling_the_clothoid_to_arc_seam_are_certified() {
    // The seam at 100 and the profile's at 120.
    let (graph, root, exact) = deck_run(&[90.0, 110.0, 125.0], StationFrame::Section, BUDGET);
    assert_certified(&graph, root, BUDGET, &exact);
}

#[test]
fn a_tight_budget_is_certified() {
    let tight = 1e-4;
    let (graph, root, exact) = deck_run(&[60.0, 140.0], StationFrame::Section, tight);
    let bound = assert_certified(&graph, root, tight, &exact);
    assert!(bound <= tight, "{bound}");
}

#[test]
fn a_run_along_a_bare_clothoid_is_certified() {
    let c = 1.0 / 1800.0;
    for frame_kind in [StationFrame::Section, StationFrame::Plan] {
        let mut b = GeometryGraphBuilder::new();
        let basis = b
            .push_value(Curve2::Intrinsic(Intrinsic2::new(
                frame(),
                CurvatureLaw::clothoid(0.0, 60.0 * c, 60.0),
                60.0,
            )))
            .unwrap();
        let narrow = rectangle(&mut b, 1.0, 1.0);
        let wide = rectangle(&mut b, 2.0, 1.0);
        let sections = vec![
            SectionAtStation::new(narrow, Station::new(5.0, StationOffsets::ZERO)),
            SectionAtStation::new(wide, Station::new(55.0, StationOffsets::new(0.3, 0.0, 0.0))),
        ];
        let root = spine(&mut b, basis, sections.clone(), frame_kind);
        let graph = b.finish(vec![root]).unwrap();
        let exact = oracle(
            &graph,
            &sections,
            Box::new(move |s| {
                let f = frame();
                let local = fresnel(0.5 * c, s);
                let p = f.origin + f.x * local.x + f.y * local.y;
                let t = f.x * (0.5 * c * s * s).cos() + f.y * (0.5 * c * s * s).sin();
                Reference {
                    point: Vec3::new(p.x, p.y, 0.0),
                    tangent: Vec3::new(t.x, t.y, 0.0),
                    lateral: Vec3::new(-t.y, t.x, 0.0),
                    up: Vec3::Z,
                }
            }),
            BUDGET,
        );
        assert_certified(&graph, root, BUDGET, &exact);
    }
}

/// The cant over the alignment, in closed form: a ramp on the first line, a
/// second across the first clothoid, a plateau on the arc and the way back.
const CANT_B: Scalar = 1.5;

fn cant_law() -> CantLaw {
    CantLaw::new(vec![
        CantPiece::linear(40.0, 0.0, 0.05),
        CantPiece::linear(60.0, 0.05, 0.1),
        CantPiece::constant(80.0, 0.1),
        CantPiece::linear(60.0, 0.1, 0.0),
        CantPiece::constant(40.0, 0.0),
    ])
}

fn cant_value(d: Scalar) -> Scalar {
    match d {
        d if d < 40.0 => 0.05 * d / 40.0,
        d if d < 100.0 => 0.05 + 0.05 * (d - 40.0) / 60.0,
        d if d < 180.0 => 0.1,
        d if d < 240.0 => 0.1 * (240.0 - d) / 60.0,
        _ => 0.0,
    }
}

/// The rolled section frame of the canted alignment at `d`: the horizontal
/// left normal and `u = t x n` turned by the roll the convention takes from
/// the cant.
fn canted_at(a: &Alignment, d: Scalar, convention: BankConvention) -> Reference {
    let (p, t) = a.world(d);
    let (z, g) = profile_at(d);
    let tangent = Vec3::new(t.x, t.y, g) / (1.0 + g * g).sqrt();
    let normal = Vec3::new(-t.y, t.x, 0.0);
    let square_up = tangent.cross(normal);
    let sine = match convention {
        BankConvention::TangentRotation => cant_value(d) / CANT_B,
        _ => cant_value(d) * (1.0 + g * g).sqrt() / CANT_B,
    };
    let (sin, cos) = sine.asin().sin_cos();
    Reference {
        point: Vec3::new(p.x, p.y, z),
        tangent,
        lateral: cos * normal + sin * square_up,
        up: -sin * normal + cos * square_up,
    }
}

#[test]
fn a_canted_deck_along_the_alignment_is_certified() {
    for convention in [
        BankConvention::TangentRotation,
        BankConvention::VerticalRise,
    ] {
        let a = alignment();
        let mut b = GeometryGraphBuilder::new();
        let banked = Curve3::Banked(Banked3::new(
            Elevated3::new(Curve2::Intrinsic(a.curve.clone()), profile()),
            cant_law(),
            CantLaw::zero(280.0),
            CANT_B,
            convention,
        ));
        let basis = b.push_value(banked).unwrap();
        let narrow = rectangle(&mut b, 10.0, 0.4);
        let wide = rectangle(&mut b, 12.0, 0.4);
        // Across the cant's ramps and the plan's and profile's seams.
        let sections: Vec<SectionAtStation> = [10.0, 75.0, 130.0, 205.0, 270.0]
            .iter()
            .enumerate()
            .map(|(i, &d)| {
                let (profile, drift) = if i % 2 == 0 {
                    (narrow, 0.0)
                } else {
                    (wide, 1.5)
                };
                SectionAtStation::new(
                    profile,
                    Station::new(d, StationOffsets::new(drift, 0.2, 0.0)),
                )
            })
            .collect();
        let root = spine(&mut b, basis, sections.clone(), StationFrame::Section);
        let graph = b.finish(vec![root]).unwrap();
        let exact = oracle(
            &graph,
            &sections,
            Box::new(move |d| canted_at(&a, d, convention)),
            BUDGET,
        );
        let bound = assert_certified(&graph, root, BUDGET, &exact);
        eprintln!("{convention:?}: certified {bound:e}");
    }
}

#[test]
fn a_canted_run_over_a_chain_plan_stays_unbounded_by_name() {
    let mut b = GeometryGraphBuilder::new();
    let chain = Curve2::Chain(Chain2::new(
        frame(),
        vec![ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 280.0),
            length: 280.0,
        }],
    ));
    let banked = Curve3::Banked(Banked3::new(
        Elevated3::new(chain, ElevationLaw::level(0.0)),
        CantLaw::new(vec![CantPiece::constant(280.0, 0.1)]),
        CantLaw::zero(280.0),
        CANT_B,
        BankConvention::TangentRotation,
    ));
    let basis = b.push_value(banked).unwrap();
    let profile = rectangle(&mut b, 3.0, 0.5);
    let sections = vec![
        SectionAtStation::new(profile, Station::new(10.0, StationOffsets::ZERO)),
        SectionAtStation::new(profile, Station::new(60.0, StationOffsets::ZERO)),
    ];
    let root = spine(&mut b, basis, sections, StationFrame::Section);
    let graph = b.finish(vec![root]).unwrap();
    match run_bound(&graph, root, BUDGET).1 {
        DeviationBound::Unbounded(reason) => {
            assert!(reason.contains("certified derivative bounds"), "{reason}");
            assert!(reason.contains("a banked curve over one"), "{reason}");
        }
        other => panic!("bounded: {other:?}"),
    }
}
