//! A slab extruded downward from its top, cut through by openings flush
//! with both its faces, and openings flush or nearly flush with a host far
//! from the origin (#291).
//!
//! The slab is `0.3` thick, extruded downward from its top (placement axis
//! `(0, 0, -1)`, extrusion direction `(0, 0, 1)`), over an 18-vertex
//! outline about `138 x 73`; its 22 rectangular openings are each extruded
//! `0.3` downward from `z = 0.3`, so each is flush with both faces. Some
//! share a line with the outline (a notch's edge, a step, the outer
//! edges), two share a side with each other, one lies outside the outline
//! against a notch. Subtracted in the slab's frame, then placed near
//! `(643 692, 5 649 156)` turned by 2.27 degrees.
//!
//! mesh-compile 0.3.15 meshed a real slab like it; the #276 snap moved its
//! openings' vertices onto planes within rounding of the slab's faces and
//! the result was refused as a void tangent to its host's face. The snap
//! now lands only exactly, on axis-aligned faces, by a rounding residue,
//! and a snapped boolean that is refused falls back to the operands as
//! given.

use std::f64::consts::PI;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Interval, Point2, Scalar, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve2, Line2};
use axiolid_inspect::{enclosed_volume, VolumeInterval};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::deviation::SNAPPED_OPERANDS;
use axiolid_mesh_compile::{DeviationBound, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

const DEPTH: Scalar = 0.3;

/// The outline in the slab's frame, counter-clockwise.
const OUTLINE: [(Scalar, Scalar); 18] = [
    (0.0, 0.0),
    (60.0, 0.0),
    (60.0, -2.4),
    (133.4, -2.4),
    (133.4, 25.7),
    (128.05, 25.7),
    (128.05, 48.9),
    (133.4, 48.9),
    (133.4, 70.15),
    (80.3, 70.15),
    (80.3, 60.65),
    (52.7, 60.65),
    (52.7, 70.15),
    (0.0, 70.15),
    (0.0, 41.3),
    (-4.6, 41.3),
    (-4.6, 20.1),
    (0.0, 20.1),
];

/// Openings `[x0, x1] x [y0, y1]` in the slab's frame.
const OPENINGS: [(Scalar, Scalar, Scalar, Scalar); 22] = [
    (59.64, 62.47, 57.2, 60.65), // flush with the notch's edge
    (62.47, 65.3, 57.2, 60.65),  // and its neighbour, sharing a side
    (10.0, 12.35, 5.0, 7.15),
    (10.0, 12.35, 9.0, 11.15), // on the same lines as the one above
    (20.0, 21.2, 0.0, 1.8),    // flush with the outline's bottom edge
    (30.0, 31.2, 0.0, 1.8),
    (64.1, 66.6, -2.4, -0.6),    // flush with the step
    (131.2, 133.4, 10.0, 12.5),  // flush with the right edge
    (128.05, 130.0, 30.0, 33.3), // outside, against the right notch
    (126.0, 128.05, 40.0, 43.3), // inside, flush with the right notch
    (100.0, 102.9, 68.0, 70.15), // flush with the top edge
    (0.0, 1.45, 50.0, 52.75),    // flush with the left edge
    (-4.6, -2.0, 25.0, 27.5),    // flush with the left bump
    (40.0, 43.3, 30.0, 33.3),
    (43.3, 46.6, 36.0, 39.3), // on the line of the one above
    (70.0, 73.65, 20.0, 23.65),
    (70.0, 73.65, 25.0, 28.65),
    (55.0, 58.0, 55.0, 60.65), // flush with the notch's edge
    (70.0, 74.0, 57.0, 60.65), // and again
    (90.0, 92.35, 40.0, 42.35),
    (110.0, 112.35, 40.0, 42.35), // on the same lines as the one above
    (115.0, 118.0, 60.0, 63.0),
];

/// The opening outside the outline, which cuts nothing.
const OUTSIDE: usize = 8;

/// More than the slab's surface area, openings included: it bounds how far
/// moving every point of it by `d` changes its volume (`area * d`).
const SLAB_AREA: Scalar = 2.0e4;

/// How a frame turns its axis down.
#[derive(Debug, Clone, Copy)]
enum Frame {
    /// `x, -y, -z`: a half turn about `x`, exactly.
    HalfTurn,
    /// `x, y, -z`: reflected.
    Reflected,
    /// `from_rotation_x(PI)`: the half turn with its sine rounded to
    /// `1.2e-16`, so its faces are axis-aligned only within rounding.
    RoundedHalfTurn,
}

impl Frame {
    /// The frame at `origin`, its axis down.
    fn at(self, origin: Vec3) -> Transform3 {
        match self {
            Frame::HalfTurn => Transform3::from_cols(Vec3::X, Vec3::NEG_Y, Vec3::NEG_Z, origin),
            Frame::Reflected => Transform3::from_cols(Vec3::X, Vec3::Y, Vec3::NEG_Z, origin),
            Frame::RoundedHalfTurn => {
                Transform3::from_translation(origin) * Transform3::from_rotation_x(PI)
            }
        }
    }

    /// Whether the frame turns `y` over, so the profile is mirrored.
    fn mirrors(self) -> bool {
        !matches!(self, Frame::Reflected)
    }
}

/// How each opening's placement is computed.
#[derive(Debug, Clone, Copy)]
struct Openings {
    frame: Frame,
    /// At `z = (e + 0.3) - e` for a storey elevation `e`, a rounding off
    /// the slab's top.
    storey: bool,
    /// Through world coordinates and back: `site^-1 * (site * placement)`.
    round_trip: bool,
}

fn site() -> Transform3 {
    Transform3::from_translation(Vec3::new(643_692.0, 5_649_156.0, 0.0))
        * Transform3::from_rotation_z(2.27_f64.to_radians())
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

fn rectangle(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The slab less its openings, in its own frame, under `placement` when
/// one is given.
fn slab(host: Frame, openings: Openings, placement: Option<Transform3>) -> (GeometryGraph, NodeId) {
    let mut b = GeometryGraphBuilder::new();
    let mut push = |node| b.push(node).expect("a valid node");
    // The profile is counter-clockwise in its own plane: mirrored, the
    // outline is run backwards.
    let sign = if host.mirrors() { -1.0 } else { 1.0 };
    let mut points: Vec<Point2> = OUTLINE
        .iter()
        .map(|&(x, y)| Point2::new(x, sign * y))
        .collect();
    if host.mirrors() {
        points.reverse();
    }
    let segments = (0..points.len())
        .map(|i| line(points[i], points[(i + 1) % points.len()]))
        .collect();
    let profile = push(GeometryNode::Profile(Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })));
    let body = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: DEPTH,
    }));
    let mut body = push(GeometryNode::Instance(Instance {
        source: body,
        transform: host.at(Vec3::new(0.0, 0.0, DEPTH)),
    }));
    let elevations = [9.45, 3.0, 12.6, 6.15, 2.7];
    for (i, (x0, x1, y0, y1)) in OPENINGS.into_iter().enumerate() {
        let top = if openings.storey {
            let e: Scalar = elevations[i % elevations.len()];
            (e + DEPTH) - e
        } else {
            DEPTH
        };
        let profile = push(GeometryNode::Profile(rectangle(x1 - x0, y1 - y0)));
        let solid = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: DEPTH,
        }));
        let centre = Vec3::new((x0 + x1) / 2.0, (y0 + y1) / 2.0, top);
        let mut transform = openings.frame.at(centre);
        if openings.round_trip {
            transform = site().inverse() * (site() * transform);
        }
        let opening = push(GeometryNode::Instance(Instance {
            source: solid,
            transform,
        }));
        body = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: body,
            right: opening,
            operator: BooleanOperator::Difference,
        }));
    }
    if let Some(transform) = placement {
        body = push(GeometryNode::Instance(Instance {
            source: body,
            transform,
        }));
    }
    (b.finish(vec![body]).expect("a valid graph"), body)
}

/// The slab's closed-form volume.
fn net() -> Scalar {
    let n = OUTLINE.len();
    let outline: Scalar = (0..n)
        .map(|i| {
            let (a, b) = (OUTLINE[i], OUTLINE[(i + 1) % n]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<Scalar>()
        / 2.0;
    let cut: Scalar = OPENINGS
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != OUTSIDE)
        .map(|(_, &(x0, x1, y0, y1))| (x1 - x0) * (y1 - y0))
        .sum();
    (outline - cut) * DEPTH
}

/// The rounding of a coordinate placed `size` from the origin.
fn ulp(size: Scalar) -> Scalar {
    size * Scalar::EPSILON
}

fn options(tolerance: Tolerance) -> ExecutionOptions {
    ExecutionOptions::new(tolerance)
        .with_chord_error(1e-3)
        .expect("a valid budget")
}

fn compile(
    graph: &GeometryGraph,
    root: NodeId,
    tolerance: Tolerance,
) -> Result<(TriMesh, DeviationReport), GeomError> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(graph, root, &options(tolerance))
        .map(|(outcome, report)| (outcome.mesh, report))
}

/// The mesh alone: the slab's deviation report would measure it against
/// its exact result, which is slow in a debug build and not what these
/// cases are about.
fn mesh(graph: &GeometryGraph, root: NodeId) -> Result<TriMesh, GeomError> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new()).compile_mesh(
        graph,
        root,
        &options(Tolerance::MILLIMETRE),
    )
}

/// The largest snap the report names, `None` when nothing was snapped.
fn snapped(report: &DeviationReport) -> Option<Scalar> {
    report
        .contributions
        .iter()
        .filter(|c| c.detail == SNAPPED_OPERANDS)
        .map(|c| match c.bound {
            DeviationBound::Certified(moved) => moved,
            other => panic!("a snap is certified: {other:?}"),
        })
        .reduce(Scalar::max)
}

/// `expected` within `slack` of the certified interval.
fn assert_volume(volume: VolumeInterval, expected: Scalar, slack: Scalar, what: &str) {
    assert!(
        volume.lower - slack <= expected && expected <= volume.upper + slack,
        "{what}: [{}, {}] does not hold {expected} within {slack:e}",
        volume.lower,
        volume.upper
    );
}

#[test]
fn a_downward_slab_with_flush_through_openings_has_its_closed_form_volume() {
    for host in [Frame::HalfTurn, Frame::Reflected] {
        for frame in [Frame::HalfTurn, Frame::Reflected, Frame::RoundedHalfTurn] {
            for storey in [false, true] {
                let openings = Openings {
                    frame,
                    storey,
                    round_trip: false,
                };
                let what = format!("{host:?} host, {openings:?}");
                let (graph, root) = slab(host, openings, None);
                let solid = mesh(&graph, root).unwrap_or_else(|e| panic!("{what}: {e}"));
                let volume = enclosed_volume(&solid).unwrap_or_else(|e| panic!("{what}: {e:?}"));
                assert_volume(volume, net(), SLAB_AREA * 4.0 * ulp(140.0), &what);
                let (graph, root) = slab(host, openings, Some(site()));
                let placed = mesh(&graph, root).unwrap_or_else(|e| panic!("{what}, placed: {e}"));
                let volume =
                    enclosed_volume(&placed).unwrap_or_else(|e| panic!("{what}, placed: {e:?}"));
                let slack = SLAB_AREA * 4.0 * ulp(5.7e6);
                assert_volume(volume, net(), slack, &format!("{what}, placed"));
            }
        }
    }
}

#[test]
fn a_slab_under_a_rounded_half_turn_is_not_refused() {
    // The reproduction: its bottom face tilts by the rounded sine, and the
    // openings come back from world coordinates `1e-10` off. The #276 snap
    // moved their bottoms onto the tilted plane within rounding, some
    // vertices above it, some below, and the result was refused as "touches
    // itself along the edge (69.99.., 28.65.., -3.5e-15) to (73.64..,
    // 28.65.., -3.5e-15)". Nothing lands on a face axis-aligned only within
    // rounding now; the residues stay as the operands give them.
    let openings = Openings {
        frame: Frame::HalfTurn,
        storey: false,
        round_trip: true,
    };
    let (graph, root) = slab(Frame::RoundedHalfTurn, openings, None);
    let host = mesh(&graph, root).unwrap_or_else(|e| panic!("{e}"));
    let volume = enclosed_volume(&host).expect("a valid solid");
    // The round trip moved the openings by up to `1e-10`.
    assert_volume(volume, net(), SLAB_AREA * 1e-9, "rounded half turn");
    // Placed, its `1e-10` residues are below the placement's rounding:
    // not certifiable, but never refused by the snap.
    let (graph, root) = slab(Frame::RoundedHalfTurn, openings, Some(site()));
    mesh(&graph, root).unwrap_or_else(|e| panic!("placed: {e}"));
}

/// A small deterministic generator (xorshift64*).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn uniform(&mut self, lo: Scalar, hi: Scalar) -> Scalar {
        lo + (hi - lo) * (self.next() >> 11) as Scalar / (1u64 << 53) as Scalar
    }
}

/// Where the property's host stands: its coordinates round to about
/// `9.3e-10`.
const FAR: Vec3 = Vec3::new(612_345.37, 5_649_156.81, 0.0);

/// The wall of the property, `[-3, 3] x [-0.15, 0.15] x [0, 3]` about
/// [`FAR`].
const WALL: (Scalar, Scalar, Scalar) = (6.0, 0.3, 3.0);

/// An offset of an opening's face from the host's, and whether it is a
/// rounding residue: flush, a few ulps of [`FAR`], or an authored skin or
/// overlap of `1e-6` or `1e-4`.
fn offset(rng: &mut Rng) -> (Scalar, bool) {
    let u = ulp(FAR.y);
    let sign = if rng.below(2) == 0 { 1.0 } else { -1.0 };
    match rng.below(4) {
        0 => (0.0, true),
        1 => (sign * (1 + rng.below(4)) as Scalar * u, true),
        2 => (sign * 1e-6, false),
        _ => (sign * 1e-4, false),
    }
}

/// One opening of the property: `[x0, x1, y0, y1, z0, z1]` about [`FAR`],
/// and the same with every rounding residue closed.
struct Cut {
    given: [Scalar; 6],
    closed: [Scalar; 6],
}

impl Cut {
    /// What `bounds` cuts out of the wall.
    fn volume(bounds: &[Scalar; 6]) -> Scalar {
        let (lx, ly, lz) = WALL;
        let dx = bounds[1].min(lx / 2.0) - bounds[0].max(-lx / 2.0);
        let dy = bounds[3].min(ly / 2.0) - bounds[2].max(-ly / 2.0);
        let dz = bounds[5].min(lz) - bounds[4].max(0.0);
        dx * dy.max(0.0) * dz.max(0.0)
    }
}

/// Up to four openings through the wall along `y`, in separate slots along
/// `x`, each face flush with the wall's or offset from it (see
/// [`offset`]); some stand on the wall's floor face.
fn cuts(rng: &mut Rng) -> Vec<Cut> {
    let (_, ly, _) = WALL;
    let count = 1 + rng.below(4);
    (0..count)
        .map(|slot| {
            let cx = -2.25 + 1.5 * slot as Scalar;
            let w = rng.uniform(0.6, 1.2);
            let (d0, r0) = offset(rng);
            let (d1, r1) = offset(rng);
            let (z0, z0_closed, z1) = if rng.below(2) == 0 {
                (0.5, 0.5, 2.5)
            } else {
                let (dz, rz) = offset(rng);
                (dz, if rz { 0.0 } else { dz }, 2.2)
            };
            let given = [
                cx - w / 2.0,
                cx + w / 2.0,
                -ly / 2.0 + d0,
                ly / 2.0 + d1,
                z0,
                z1,
            ];
            let mut closed = given;
            if r0 {
                closed[2] = -ly / 2.0;
            }
            if r1 {
                closed[3] = ly / 2.0;
            }
            closed[4] = z0_closed;
            Cut { given, closed }
        })
        .collect()
}

/// The wall less `cuts`, each operand placed by `placement`: the boolean
/// is computed at the placement's coordinates.
fn wall(cuts: &[Cut], placement: Transform3) -> (GeometryGraph, NodeId) {
    let (lx, ly, lz) = WALL;
    let mut b = GeometryGraphBuilder::new();
    let mut push = |node| b.push(node).expect("a valid node");
    let profile = push(GeometryNode::Profile(rectangle(lx, ly)));
    let body = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: lz,
    }));
    let mut body = push(GeometryNode::Instance(Instance {
        source: body,
        transform: placement,
    }));
    for cut in cuts {
        let [x0, x1, y0, y1, z0, z1] = cut.given;
        let profile = push(GeometryNode::Profile(rectangle(x1 - x0, y1 - y0)));
        let solid = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: z1 - z0,
        }));
        let local = push(GeometryNode::Instance(Instance {
            source: solid,
            transform: Transform3::from_translation(Vec3::new(
                (x0 + x1) / 2.0,
                (y0 + y1) / 2.0,
                z0,
            )),
        }));
        let opening = push(GeometryNode::Instance(Instance {
            source: local,
            transform: placement,
        }));
        body = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: body,
            right: opening,
            operator: BooleanOperator::Difference,
        }));
    }
    (b.finish(vec![body]).expect("a valid graph"), body)
}

#[test]
fn near_flush_openings_far_from_the_origin_are_never_refused_after_the_snap() {
    let (lx, ly, lz) = WALL;
    let moved_only = Transform3::from_translation(FAR);
    let turned = moved_only * Transform3::from_rotation_z(2.27_f64.to_radians());
    // The wall's surface area, its openings' included.
    let area = 2.0 * (lx * ly + lx * lz + ly * lz) + 4.0 * 2.0 * (1.2 + 3.0) * ly;
    let mut rng = Rng(0x2910_5EED_0BAD_CAFE);
    let (mut compiled, mut snaps, mut kept) = (0, 0, 0);
    for case in 0..160 {
        let cuts = cuts(&mut rng);
        let axis_aligned = case % 4 != 3;
        let placement = if axis_aligned { moved_only } else { turned };
        let (graph, root) = wall(&cuts, placement);
        let what = format!("case {case}");
        // Cut as given: nothing moves at zero tolerance.
        let Ok((_, given_report)) = compile(&graph, root, Tolerance::ZERO) else {
            continue;
        };
        assert_eq!(snapped(&given_report), None, "{what}");
        compiled += 1;
        let (mesh, report) = compile(&graph, root, Tolerance::MILLIMETRE)
            .unwrap_or_else(|e| panic!("{what}: cut as given, refused after the snap: {e}"));
        if !axis_aligned {
            continue;
        }
        let moved = snapped(&report);
        let residues = cuts.iter().any(|c| c.given != c.closed);
        assert!(
            moved.is_none() || residues,
            "{what}: a snap without a rounding residue: {report:?}"
        );
        assert!(
            moved.is_none_or(|m| m <= 16.0 * ulp(FAR.y)),
            "{what}: {report:?}"
        );
        snaps += usize::from(moved.is_some());
        kept += cuts
            .iter()
            .filter(|c| c.closed[2] != -ly / 2.0 || c.closed[3] != ly / 2.0)
            .count();
        // Rounding residues closed, authored skins kept.
        let expected = lx * ly * lz - cuts.iter().map(|c| Cut::volume(&c.closed)).sum::<Scalar>();
        let volume = enclosed_volume(&mesh).unwrap_or_else(|e| panic!("{what}: {e:?}"));
        let slack = area * (moved.unwrap_or(0.0) + 4.0 * ulp(FAR.y));
        assert_volume(volume, expected, slack, &what);
    }
    // The property was exercised, not vacuously met.
    assert!(compiled >= 120, "{compiled} of 160 compiled as given");
    assert!(snaps >= 20 && kept >= 20, "{snaps} snapped, {kept} kept");
}
