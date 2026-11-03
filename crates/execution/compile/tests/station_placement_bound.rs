//! Certified placements at stations on curved bases (#405): a block and a
//! cylinder placed in the frame of a station on a circle, a clothoid, an
//! alignment of lines, arcs and clothoids (inside a piece and on a seam,
//! from both sides), and a gradient curve over a line, a circle and the
//! alignment, in the section and the plan frame, with offsets and with an
//! oriented frame, at two chord budgets. Every placement reports its
//! frame's error as a certified `StationPlacement` contribution within the
//! request, the kernel's frame lies within the error of an independent
//! exact one, and a dense sampling of the EXACT placed surface never
//! strays past the report's bound from the mesh.
//!
//! The exact frames are independent of the kernel's quadratures: the
//! plan's points from the Fresnel integrals' power series and closed
//! forms, its headings and the profiles' heights and grades in closed
//! form, and the frames built here as the section frames are defined.
//! Frames whose error is not certified (an ellipse, a B-spline, an offset,
//! a chain, a banked curve over an intrinsic plan, a vertical clothoid, a
//! basis itself placed inexactly) report their placement unbounded by name.

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Frame2, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, Circle2,
    CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Ellipse2, Intrinsic2, KnotSpec, Line2,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::station::{placement, placement_accuracy, PlacementAccuracy};
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_model::{
    CurveRelation, CurveStation, GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance,
    InstanceAtStation, NodeId, OrientedCurveStation, SeamSide, SolidOperation, Station,
    StationFrame, StationOffsets, StationOrientation,
};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

fn options(chord: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(chord)
        .unwrap()
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

// --- the exact curves, independently ---------------------------------------------

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
    // `(a s^2)^n / n!`, built term by term: the powers alone overflow.
    let mut term = 1.0;
    for n in 0..120 {
        if n > 0 {
            term *= a * s * s / Scalar::from(n);
        }
        let value = s * term / Scalar::from(2 * n + 1);
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
    /// A clothoid with curvature `c s`.
    Ramp(Scalar, Scalar),
    Rising(Scalar, Scalar),
    Arc(Scalar, Scalar),
    Falling(Scalar, Scalar),
}

impl Piece {
    fn length(self) -> Scalar {
        match self {
            Piece::Line(l)
            | Piece::Ramp(l, _)
            | Piece::Rising(l, _)
            | Piece::Arc(l, _)
            | Piece::Falling(l, _) => l,
        }
    }

    fn law(self) -> CurvatureLaw {
        match self {
            Piece::Line(_) => CurvatureLaw::Constant { curvature: 0.0 },
            Piece::Ramp(_, c) => CurvatureLaw::Polynomial {
                coefficients: vec![0.0, c],
            },
            Piece::Rising(l, r) => CurvatureLaw::clothoid(0.0, 1.0 / r, l),
            Piece::Arc(_, r) => CurvatureLaw::Constant { curvature: 1.0 / r },
            Piece::Falling(l, r) => CurvatureLaw::clothoid(1.0 / r, 0.0, l),
        }
    }

    fn turn(self) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Ramp(l, c) => 0.5 * c * l * l,
            Piece::Rising(l, r) | Piece::Falling(l, r) => 0.5 * l / r,
            Piece::Arc(l, r) => l / r,
        }
    }

    /// The heading turned over the first `s`.
    fn heading(self, s: Scalar) -> Scalar {
        match self {
            Piece::Line(_) => 0.0,
            Piece::Ramp(_, c) => 0.5 * c * s * s,
            Piece::Rising(l, r) => 0.5 * s * s / (r * l),
            Piece::Arc(_, r) => s / r,
            Piece::Falling(l, r) => self.turn() - 0.5 * (l - s) * (l - s) / (r * l),
        }
    }

    fn displacement(self, s: Scalar) -> Vec2 {
        match self {
            Piece::Line(_) => Vec2::new(s, 0.0),
            Piece::Ramp(_, c) => fresnel(0.5 * c, s),
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

    fn plan_at(&self, d: Scalar) -> Reference {
        let (p, t) = self.world(d);
        planar(Vec3::new(p.x, p.y, 0.0), Vec3::new(t.x, t.y, 0.0))
    }

    fn elevated_at(&self, d: Scalar) -> Reference {
        let (p, t) = self.world(d);
        let (z, g) = profile_at(d);
        let size = (1.0 + g * g).sqrt();
        reference_up(
            Vec3::new(p.x, p.y, z),
            Vec3::new(t.x / size, t.y / size, g / size),
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

/// An exact section frame.
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

fn planar(point: Vec3, tangent: Vec3) -> Reference {
    Reference {
        point,
        tangent,
        lateral: Vec3::new(-tangent.y, tangent.x, 0.0),
        up: Vec3::Z,
    }
}

fn plan_of(r: Reference) -> Reference {
    let t = Vec3::new(r.tangent.x, r.tangent.y, 0.0);
    planar(r.point, t / t.length())
}

// --- the cases --------------------------------------------------------------------

/// A basis curve.
#[derive(Clone)]
enum Base {
    Two(Curve2),
    Three(Curve3),
}

/// A station on a basis, with the exact frame there.
struct Case {
    name: &'static str,
    base: Base,
    distance: Scalar,
    side: SeamSide,
    frame: StationFrame,
    exact: Reference,
}

fn circle2() -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(1000.0, 2000.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 300.0,
    })
}

/// The exact frame on the circle of `circle2` at the plan distance `d`:
/// its position and its plan tangent.
fn on_circle(d: Scalar) -> (Vec3, Vec3) {
    let (sin, cos) = (d / 300.0).sin_cos();
    (
        Vec3::new(1000.0 + 300.0 * cos, 2000.0 + 300.0 * sin, 0.0),
        Vec3::new(-sin, cos, 0.0),
    )
}

fn cases() -> Vec<Case> {
    let mut out = Vec::new();
    let (point, tangent) = on_circle(250.0);
    out.push(Case {
        name: "circle",
        base: Base::Two(circle2()),
        distance: 250.0,
        side: SeamSide::Outgoing,
        frame: StationFrame::Section,
        exact: planar(point, tangent),
    });

    let ramp = Alignment::new(vec![Piece::Ramp(120.0, 1.0 / 300.0)]);
    out.push(Case {
        name: "clothoid",
        base: Base::Two(Curve2::Intrinsic(ramp.curve.clone())),
        distance: 90.0,
        side: SeamSide::Outgoing,
        frame: StationFrame::Section,
        exact: ramp.plan_at(90.0),
    });

    let a = alignment();
    let plan = Curve2::Intrinsic(a.curve.clone());
    // Inside the clothoid, and on its seam to the arc from both sides.
    for (name, d, side) in [
        ("alignment, inside a clothoid", 70.0, SeamSide::Outgoing),
        ("alignment, on a seam, outgoing", 100.0, SeamSide::Outgoing),
        ("alignment, on a seam, incoming", 100.0, SeamSide::Incoming),
    ] {
        out.push(Case {
            name,
            base: Base::Two(plan.clone()),
            distance: d,
            side,
            frame: StationFrame::Section,
            exact: a.plan_at(d),
        });
    }

    // A grade of 2% from the origin over a line.
    let k = (1.0_f64 + 0.02 * 0.02).sqrt();
    out.push(Case {
        name: "gradient over a line",
        base: Base::Three(Curve3::Elevated(Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::new(10.0, 20.0),
                direction: Vec2::X,
            }),
            ElevationLaw::Polynomial {
                coefficients: vec![0.0, 0.02],
            },
        ))),
        distance: 40.0,
        side: SeamSide::Outgoing,
        frame: StationFrame::Section,
        exact: reference_up(Vec3::new(50.0, 20.0, 0.8), Vec3::new(1.0, 0.0, 0.02) / k),
    });

    // A parabolic profile over the circle: z = 5 + 0.01 d + c d^2.
    let c = (0.0148 - 0.01) / (2.0 * 200.0);
    let grade: Scalar = 0.01 + 2.0 * c * 250.0;
    let (point, tangent) = on_circle(250.0);
    let size = (1.0 + grade * grade).sqrt();
    out.push(Case {
        name: "gradient over a circle",
        base: Base::Three(Curve3::Elevated(Elevated3::new(
            circle2(),
            ElevationLaw::parabolic(5.0, 0.01, 0.0148, 200.0),
        ))),
        distance: 250.0,
        side: SeamSide::Outgoing,
        frame: StationFrame::Section,
        exact: reference_up(
            Vec3::new(point.x, point.y, 5.0 + 0.01 * 250.0 + c * 250.0 * 250.0),
            Vec3::new(tangent.x / size, tangent.y / size, grade / size),
        ),
    });

    // The gradient curve over the alignment, in both frames, inside the
    // profile's sag and crest.
    let elevated = Curve3::Elevated(Elevated3::new(plan, profile()));
    for (name, d, frame) in [
        (
            "gradient over the alignment, sag",
            75.0,
            StationFrame::Section,
        ),
        (
            "gradient over the alignment, crest",
            205.0,
            StationFrame::Section,
        ),
        (
            "gradient over the alignment, plan",
            205.0,
            StationFrame::Plan,
        ),
    ] {
        let exact = match frame {
            StationFrame::Plan => plan_of(a.elevated_at(d)),
            _ => a.elevated_at(d),
        };
        out.push(Case {
            name,
            base: Base::Three(elevated.clone()),
            distance: d,
            side: SeamSide::Outgoing,
            frame,
            exact,
        });
    }
    out
}

/// The orientations placed in: the base frame, and one turning every axis
/// (`axis` for the up, `reference` for the tangent).
const TURNED: (Vec3, Vec3) = (Vec3::new(0.3, 0.2, 1.0), Vec3::new(1.0, 0.1, 0.0));

/// The exact frame's axes (tangent, lateral, up), turned by `orientation`
/// as `SectionFrame::oriented` turns them.
fn exact_axes(base: Reference, turned: bool) -> [Vec3; 3] {
    if !turned {
        return [base.tangent, base.lateral, base.up];
    }
    let axis = TURNED.0.normalize();
    let reference = TURNED.1.normalize();
    let x = (reference - reference.dot(axis) * axis).normalize();
    let y = axis.cross(x);
    let world = |v: Vec3| v.x * base.tangent + v.y * base.lateral + v.z * base.up;
    [world(x), world(y), world(axis)]
}

const OFFSETS: StationOffsets = StationOffsets::new(2.0, 1.0, 0.5);

/// The exact station's point: the offsets read in the base frame.
fn exact_point(base: Reference) -> Vec3 {
    base.point
        + OFFSETS.lateral * base.lateral
        + OFFSETS.vertical * base.up
        + OFFSETS.longitudinal * base.tangent
}

// --- the sources --------------------------------------------------------------------

/// A 4 x 2 x 3 block, moved so its farthest vertex has `|x| + |y| + |z|`
/// 6, as local points.
fn block_corners() -> [Vec3; 8] {
    let (lo, hi) = (Vec3::new(-1.5, -0.5, -1.0), Vec3::new(2.5, 1.5, 2.0));
    [0, 1, 2, 3, 4, 5, 6, 7].map(|i| {
        Vec3::new(
            if i & 1 == 0 { lo.x } else { hi.x },
            if i & 2 == 0 { lo.y } else { hi.y },
            if i & 4 == 0 { lo.z } else { hi.z },
        )
    })
}

fn push_block(b: &mut GeometryGraphBuilder) -> NodeId {
    let profile = b
        .push_value(Profile::Rectangle(RectangleProfile {
            x: 4.0,
            y: 2.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }))
        .unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 3.0,
        }))
        .unwrap();
    // The rectangle is centred on the origin, the extrusion runs from
    // `z = 0` to `3`: moved to `[-1.5, 2.5] x [-0.5, 1.5] x [-1, 2]`.
    b.push(GeometryNode::Instance(Instance {
        source: solid,
        transform: Transform3::from_translation(Vec3::new(0.5, 0.5, -1.0)),
    }))
    .unwrap()
}

/// A cylinder of radius 1.5 and height 3, moved to `x, y` centre `(0.5,
/// 0.5)`, `z` from -1 to 2.
fn push_cylinder(b: &mut GeometryGraphBuilder) -> NodeId {
    let profile = b
        .push_value(Profile::Circle(CircleProfile {
            radius: 1.5,
            thickness: None,
        }))
        .unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 3.0,
        }))
        .unwrap();
    b.push(GeometryNode::Instance(Instance {
        source: solid,
        transform: Transform3::from_translation(Vec3::new(0.5, 0.5, -1.0)),
    }))
    .unwrap()
}

/// `n + 1` by `n + 1` points on each face of the block.
fn block_points(n: usize) -> Vec<Vec3> {
    let (lo, hi) = (Vec3::new(-1.5, -0.5, -1.0), Vec3::new(2.5, 1.5, 2.0));
    let at = |lo: Scalar, hi: Scalar, i: usize| lo + (hi - lo) * i as Scalar / n as Scalar;
    let mut out = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            for (end, other) in [(lo, hi), (hi, lo)] {
                out.push(Vec3::new(end.x, at(lo.y, hi.y, i), at(lo.z, hi.z, j)));
                out.push(Vec3::new(at(lo.x, hi.x, i), end.y, at(lo.z, hi.z, j)));
                out.push(Vec3::new(at(lo.x, hi.x, i), at(lo.y, hi.y, j), end.z));
                let _ = other;
            }
        }
    }
    out
}

/// Rings of the cylinder: `rings + 1` heights, `angles` points each.
fn cylinder_points(rings: usize, angles: usize) -> Vec<Vec3> {
    let mut out = Vec::new();
    for i in 0..=rings {
        let z = -1.0 + 3.0 * i as Scalar / rings as Scalar;
        for k in 0..angles {
            let (sin, cos) = (std::f64::consts::TAU * k as Scalar / angles as Scalar).sin_cos();
            out.push(Vec3::new(0.5 + 1.5 * cos, 0.5 + 1.5 * sin, z));
        }
    }
    out
}

// --- compiling -------------------------------------------------------------------

fn push_base(b: &mut GeometryGraphBuilder, base: &Base) -> NodeId {
    match base {
        Base::Two(curve) => b.push_value(curve.clone()).unwrap(),
        Base::Three(curve) => b.push_value(curve.clone()).unwrap(),
    }
}

fn station(
    basis: NodeId,
    distance: Scalar,
    side: SeamSide,
    frame: StationFrame,
    turned: bool,
) -> OrientedCurveStation {
    let mut curve = CurveStation::new(basis, Station::new(distance, OFFSETS));
    curve.frame = frame;
    let orientation = match turned {
        true => StationOrientation::new(Some(TURNED.0), Some(TURNED.1)),
        false => StationOrientation::new(None, None),
    };
    OrientedCurveStation::new(curve, orientation).with_seam_side(side)
}

fn contribution(report: &DeviationReport, path: DeviationPath) -> Option<DeviationBound> {
    report
        .contributions
        .iter()
        .find(|c| c.path == path)
        .map(|c| c.bound)
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

/// A source placed at the case's station, compiled at `chord`.
struct Placed {
    graph: GeometryGraph,
    node: NodeId,
}

fn placed(case: &Case, turned: bool, source: fn(&mut GeometryGraphBuilder) -> NodeId) -> Placed {
    let mut b = GeometryGraphBuilder::new();
    let basis = push_base(&mut b, &case.base);
    let source = source(&mut b);
    let at = station(basis, case.distance, case.side, case.frame, turned);
    let node = b.push_value(InstanceAtStation::new(source, at)).unwrap();
    let graph = b.finish(vec![node]).unwrap();
    Placed { graph, node }
}

/// What the placement of `case` promises, at `chord`, for one source.
fn check(
    case: &Case,
    turned: bool,
    chord: Scalar,
    source: fn(&mut GeometryGraphBuilder) -> NodeId,
    points: &[Vec3],
    block: bool,
) {
    let what = format!("{} (turned: {turned}, chord {chord})", case.name);
    let Placed { graph, node } = placed(case, turned, source);

    // The accuracy is certified, and the kernel's frame lies within it of
    // the exact one.
    let accuracy = placement_accuracy(&graph, node, 1e-3 * chord).unwrap();
    let PlacementAccuracy::Certified(bound) = accuracy else {
        panic!("{what}: {accuracy:?}")
    };
    let resolved = placement(&graph, node).unwrap();
    assert!(!resolved.exact, "{what}");
    let [t, l, u] = exact_axes(case.exact, turned);
    let point = exact_point(case.exact);
    let moved = (resolved.transform.translation - point).length();
    assert!(
        moved <= bound.point_error,
        "{what}: point {moved} > {bound:?}"
    );
    for (axis, exact) in [(Vec3::X, t), (Vec3::Y, l), (Vec3::Z, u)] {
        let off = (resolved.transform.transform_vector3(axis) - exact).length();
        assert!(off <= bound.axis_error, "{what}: axis {off} > {bound:?}");
    }

    // The report names the frame's error, certified, within the request.
    let (outcome, report) = compiler()
        .compile_mesh_with_deviation(&graph, node, &options(chord))
        .unwrap();
    let Some(DeviationBound::Certified(frame)) =
        contribution(&report, DeviationPath::StationPlacement)
    else {
        panic!("{what}: {report:?}")
    };
    assert!(report.unbounded().next().is_none(), "{what}: {report:?}");
    assert!(report.meets_requested(), "{what}: {report:?}");
    assert!(frame <= 0.01 * chord, "{what}: the frame's error {frame}");
    assert!(frame >= bound.point_error, "{what}: {frame} < {bound:?}");

    // Every other contribution grows by it, and it reaches the block's
    // farthest vertex (`|x| + |y| + |z|` 6).
    if block {
        let reach = bound.point_error + bound.axis_error * 6.0;
        assert!(frame >= reach, "{what}: {frame} < {reach}");
        let others: Vec<_> = report
            .contributions
            .iter()
            .filter(|c| {
                !matches!(
                    c.path,
                    DeviationPath::StationPlacement | DeviationPath::Placement
                )
            })
            .collect();
        assert!(!others.is_empty(), "{what}: {report:?}");
        for c in others {
            let value = c.bound.value().expect("bounded");
            assert!(value >= frame, "{what}: {c:?} is below the frame's {frame}");
        }
    }

    // The exact placed surface stays within the report's bound.
    let exact: Vec<Point3> = points
        .iter()
        .map(|q| Point3::from(point + q.x * t + q.y * l + q.z * u))
        .collect();
    let far = measured(&outcome.mesh, &exact);
    let worst = report.bound.expect("a bound");
    assert!(
        far <= worst + 1e-12,
        "{what}: the exact placement strays {far} from the mesh, past {worst}"
    );
}

#[test]
fn a_block_placed_at_a_station_is_certified_at_both_budgets() {
    let points = block_points(40);
    for case in cases() {
        for turned in [false, true] {
            for chord in [1e-3, 1e-6] {
                check(&case, turned, chord, push_block, &points, true);
            }
        }
    }
}

#[test]
fn a_cylinder_placed_at_a_station_is_certified() {
    let points = cylinder_points(40, 80);
    for case in cases() {
        for turned in [false, true] {
            check(&case, turned, 1e-3, push_cylinder, &points, false);
        }
    }
}

#[test]
fn the_block_reaches_the_lever_the_bound_assumes() {
    // The farthest vertex has `|x| + |y| + |z|` 6.
    let reach = block_corners()
        .iter()
        .map(|v| v.x.abs() + v.y.abs() + v.z.abs())
        .fold(0.0, Scalar::max);
    assert_eq!(reach, 6.0);
}

// --- unbounded by name ---------------------------------------------------------------

fn unbounded_cases() -> Vec<(&'static str, Base, Scalar)> {
    let plan = Curve2::Intrinsic(Intrinsic2::new(
        frame(),
        CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
        60.0,
    ));
    let ellipse = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(30.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    });
    let spline = Curve3::BSpline(BSplineCurve {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(8.0, 2.0, 0.5),
            Point3::new(16.0, -3.0, 1.5),
            Point3::new(24.0, 2.0, 2.0),
            Point3::new(30.0, 6.0, 2.2),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let banked = Curve3::Banked(Banked3::new(
        Elevated3::new(plan.clone(), ElevationLaw::level(0.0)),
        CantLaw::new(vec![CantPiece::constant(60.0, 0.15)]),
        CantLaw::zero(60.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    let chain = Curve2::Chain(Chain2::new(
        frame(),
        vec![ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
            length: 60.0,
        }],
    ));
    let vertical = Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::intrinsic(0.0, 0.01, CurvatureLaw::clothoid(0.0, 1.0 / 300.0, 60.0)),
    ));
    vec![
        ("ellipse", Base::Two(ellipse), 10.0),
        ("b-spline", Base::Three(spline), 10.0),
        ("banked over a clothoid", Base::Three(banked), 10.0),
        ("chain", Base::Two(chain), 10.0),
        ("vertical clothoid", Base::Three(vertical), 20.0),
    ]
}

fn assert_unbounded(graph: &GeometryGraph, node: NodeId, what: &str) {
    let accuracy = placement_accuracy(graph, node, 1e-6).unwrap();
    let PlacementAccuracy::Unbounded(reason) = accuracy else {
        panic!("{what}: {accuracy:?}")
    };
    assert!(reason.contains("no certified error"), "{what}: {reason}");
    let (_, report) = compiler()
        .compile_mesh_with_deviation(graph, node, &options(1e-3))
        .unwrap();
    assert!(
        matches!(
            contribution(&report, DeviationPath::StationPlacement),
            Some(DeviationBound::Unbounded(r)) if r == reason
        ),
        "{what}: {report:?}"
    );
    assert_eq!(report.bound, None, "{what}");
    assert!(!report.meets_requested(), "{what}");
}

#[test]
fn a_frame_without_a_certified_error_is_unbounded_by_name() {
    for (name, base, distance) in unbounded_cases() {
        let case = Case {
            name,
            base,
            distance,
            side: SeamSide::Outgoing,
            frame: StationFrame::Section,
            exact: planar(Vec3::ZERO, Vec3::X),
        };
        let Placed { graph, node } = placed(&case, false, push_block);
        assert_unbounded(&graph, node, name);
    }
}

#[test]
fn an_offset_curve_and_a_nested_placement_are_unbounded_by_name() {
    // An offset of a circle: its frame is read by quadrature.
    let mut b = GeometryGraphBuilder::new();
    let circle = b.push_value(circle2()).unwrap();
    let offset = b
        .push(GeometryNode::CurveRelation(CurveRelation::Offset {
            basis: circle,
            distance: 2.0,
            reference_direction: None,
        }))
        .unwrap();
    let block = push_block(&mut b);
    let at = CurveStation::new(offset, Station::at(20.0)).into();
    let node = b.push_value(InstanceAtStation::new(block, at)).unwrap();
    let graph = b.finish(vec![node]).unwrap();
    assert_unbounded(&graph, node, "offset");

    // A solid at a station along a line that is itself placed on an arc.
    let mut b = GeometryGraphBuilder::new();
    let arc = b
        .push_value(Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 30.0,
        }))
        .unwrap();
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let on_arc = CurveStation::new(arc, Station::at(12.0)).into();
    let placed_line = b.push_value(InstanceAtStation::new(line, on_arc)).unwrap();
    let block = push_block(&mut b);
    let at = CurveStation::new(placed_line, Station::at(3.0)).into();
    let node = b.push_value(InstanceAtStation::new(block, at)).unwrap();
    let graph = b.finish(vec![node]).unwrap();
    assert_unbounded(&graph, node, "nested");
}

#[test]
fn a_line_is_exact_and_a_bad_tolerance_or_node_is_refused() {
    let mut b = GeometryGraphBuilder::new();
    let line = b
        .push_value(Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }))
        .unwrap();
    let block = push_block(&mut b);
    let at = CurveStation::new(line, Station::at(3.0)).into();
    let node = b.push_value(InstanceAtStation::new(block, at)).unwrap();
    let graph = b.finish(vec![node, block]).unwrap();
    assert_eq!(
        placement_accuracy(&graph, node, 1e-6).unwrap(),
        PlacementAccuracy::Exact
    );
    for tolerance in [0.0, -1.0, Scalar::NAN, Scalar::INFINITY] {
        assert!(placement_accuracy(&graph, node, tolerance).is_err());
    }
    assert!(placement_accuracy(&graph, block, 1e-6).is_err());
    let (_, report) = compiler()
        .compile_mesh_with_deviation(&graph, node, &options(1e-3))
        .unwrap();
    assert_eq!(contribution(&report, DeviationPath::StationPlacement), None);
}

// --- turning an orientation, conditioned by the angle ---------------------------------

/// A double-double number: an unevaluated sum of two `f64`, about 106 bits.
#[derive(Clone, Copy)]
struct Dd(f64, f64);

fn two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    let bb = s - a;
    Dd(s, (a - (s - bb)) + (b - bb))
}

fn fast(hi: f64, lo: f64) -> Dd {
    let s = hi + lo;
    Dd(s, lo - (s - hi))
}

impl Dd {
    fn add(self, o: Dd) -> Dd {
        let s = two_sum(self.0, o.0);
        fast(s.0, s.1 + self.1 + o.1)
    }
    fn neg(self) -> Dd {
        Dd(-self.0, -self.1)
    }
    fn sub(self, o: Dd) -> Dd {
        self.add(o.neg())
    }
    fn mul(self, o: Dd) -> Dd {
        let p = self.0 * o.0;
        let e = self.0.mul_add(o.0, -p);
        fast(p, e + self.0 * o.1 + self.1 * o.0)
    }
    fn div(self, o: Dd) -> Dd {
        let q1 = self.0 / o.0;
        let r = self.sub(o.mul(Dd(q1, 0.0)));
        let q2 = r.0 / o.0;
        let r = r.sub(o.mul(Dd(q2, 0.0)));
        let q3 = r.0 / o.0;
        fast(q1, q2).add(Dd(q3, 0.0))
    }
    fn sqrt(self) -> Dd {
        let s = self.0.sqrt();
        let r = self.sub(Dd(s, 0.0).mul(Dd(s, 0.0)));
        fast(s, r.0 / (2.0 * s))
    }
}

type V = [Dd; 3];

fn dot(a: V, b: V) -> Dd {
    (0..3).fold(Dd(0.0, 0.0), |sum, i| sum.add(a[i].mul(b[i])))
}

fn cross(a: V, b: V) -> V {
    [
        a[1].mul(b[2]).sub(a[2].mul(b[1])),
        a[2].mul(b[0]).sub(a[0].mul(b[2])),
        a[0].mul(b[1]).sub(a[1].mul(b[0])),
    ]
}

fn unit(a: V) -> V {
    let n = dot(a, a).sqrt();
    a.map(|c| c.div(n))
}

fn dd(v: Vec3) -> V {
    [Dd(v.x, 0.0), Dd(v.y, 0.0), Dd(v.z, 0.0)]
}

/// The orientation's axes as `(x, y, axis)` in the base frame, by an exact
/// Gram-Schmidt of the stored inputs.
fn gram_schmidt(axis: Vec3, reference: Vec3) -> [Vec3; 3] {
    let (a, r) = (unit(dd(axis)), unit(dd(reference)));
    let along = dot(r, a);
    let w = [0, 1, 2].map(|i| r[i].sub(along.mul(a[i])));
    let x = unit(w);
    let y = cross(a, x);
    let to = |v: V| Vec3::new(v[0].0, v[1].0, v[2].0);
    [to(x), to(y), to(a)]
}

#[test]
fn a_nearly_parallel_orientation_is_charged_its_conditioning() {
    let case = cases().into_iter().next().expect("the circle");
    let base = case.exact;
    for s in [1e-3, 1e-6, 1e-8] {
        let axis = Vec3::new(0.6, 0.0, 0.8);
        let reference = Vec3::new(0.6 + 0.8 * s, 0.3 * s, 0.8 - 0.6 * s);
        let mut b = GeometryGraphBuilder::new();
        let basis = push_base(&mut b, &case.base);
        let block = push_block(&mut b);
        let curve = CurveStation::new(basis, Station::at(case.distance));
        let at =
            OrientedCurveStation::new(curve, StationOrientation::new(Some(axis), Some(reference)));
        let node = b.push_value(InstanceAtStation::new(block, at)).unwrap();
        let graph = b.finish(vec![node]).unwrap();
        let accuracy = placement_accuracy(&graph, node, 1e-9).unwrap();
        let PlacementAccuracy::Certified(bound) = accuracy else {
            panic!("s {s}: {accuracy:?}")
        };
        let [x, y, a] = gram_schmidt(axis, reference);
        let world = |v: Vec3| v.x * base.tangent + v.y * base.lateral + v.z * base.up;
        let exact = [world(x), world(y), world(a)];
        let resolved = placement(&graph, node).unwrap();
        for (k, unit) in [Vec3::X, Vec3::Y, Vec3::Z].into_iter().enumerate() {
            let off = (resolved.transform.transform_vector3(unit) - exact[k]).length();
            assert!(
                off <= bound.axis_error,
                "s {s}, axis {k}: {off} > {bound:?}"
            );
        }
    }
}
