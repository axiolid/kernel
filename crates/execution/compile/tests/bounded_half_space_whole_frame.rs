//! A bounded half-space's placement is read whole (#322): mirrored, scaled
//! and sheared frames, along either prism axis, in both compilers.
//!
//! The box `[-2, 2] x [-3, 3] x [0, 2]` loses what lies above `z = 1.5`
//! within the vertical prism of an L-shaped boundary, chiral so a mirror
//! shows, placed partly over the box's edge so that where the footprint
//! lands, not only its area, decides the volume: the area of the footprint
//! within the box's plan, times `0.5`. The footprint is the boundary placed
//! by the frame's `x` and `y` as written, which is what
//! `BoundedHalfSpaceAlongZ` means for every frame and `BoundedHalfSpace`
//! for every frame whose `z` lies on the normal's side. Along the plane
//! normal a right-handed frame whose `z` points against the normal is
//! turned over, as it always was: its boundary is laid with `y` along
//! `normal x x`. Before #322 the mesh compiler normalised the frame's axes
//! and rebuilt `y` from the normal, so a scaled frame cut its boundary at
//! unit scale and a mirrored one cut its mirror image, with no error.

use axiolid_contracts::{ExecutionOptions, GeomError};
use std::f64::consts::FRAC_PI_2;

use axiolid_core::{
    BooleanOperator, Frame2, Interval, Plane3, Point2, Point3, Scalar, Tolerance, Transform3, Vec2,
    Vec3,
};
use axiolid_curve::{Circle2, Curve2, Line2, Polyline2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{ReferenceExactCompiler, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};

/// Which variant states the clip: the prism along the plane's normal
/// (`BoundedHalfSpace`) or along the placement's `z`
/// (`BoundedHalfSpaceAlongZ`, #321).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PrismAxis {
    PlaneNormal,
    PlacementZ,
}

impl PrismAxis {
    fn operation(
        self,
        half_space: NodeId,
        boundary: NodeId,
        placement: Transform3,
    ) -> SolidOperation {
        match self {
            Self::PlaneNormal => SolidOperation::BoundedHalfSpace {
                half_space,
                boundary,
                placement,
            },
            Self::PlacementZ => SolidOperation::BoundedHalfSpaceAlongZ {
                half_space,
                boundary,
                placement,
            },
        }
    }
}

/// The L, counter-clockwise.
const L_SHAPE: [(Scalar, Scalar); 6] = [
    (0.0, 0.0),
    (2.0, 0.0),
    (2.0, 1.0),
    (1.0, 1.0),
    (1.0, 2.0),
    (0.0, 2.0),
];

/// The box's plan.
const PLAN: (Point2, Point2) = (Point2::new(-2.0, -3.0), Point2::new(2.0, 3.0));

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
        .with_chord_error(1e-3)
        .expect("a positive budget")
}

/// The box less the clip of the L under `frame` along `axis`, the plane's
/// normal `normal` (`+z` or `-z`) and the side above it kept.
fn clipped(frame: Transform3, axis: PrismAxis, normal: Vec3) -> (GeometryGraph, NodeId) {
    let mut g = GeometryGraphBuilder::new();
    let mut push = |node| g.push(node).expect("a valid node");
    let profile = push(GeometryNode::Profile(Profile::Rectangle(
        RectangleProfile {
            x: 4.0,
            y: 6.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        },
    )));
    let block = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: 2.0,
    }));
    let half_space = push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, 1.5),
            normal,
        },
        agreement: normal.z > 0.0,
    }));
    let boundary = push(GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
        points: L_SHAPE.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        closed: true,
    })));
    let tool = push(GeometryNode::SolidOperation(
        axis.operation(half_space, boundary, frame),
    ));
    let root = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: block,
        right: tool,
        operator: BooleanOperator::Difference,
    }));
    (g.finish(vec![root]).expect("a valid graph"), root)
}

/// The L placed in the plan by `origin + u x + v y`.
fn footprint(origin: Point2, x: Point2, y: Point2) -> Vec<Point2> {
    L_SHAPE
        .iter()
        .map(|&(u, v)| origin + x * u + y * v)
        .collect()
}

/// The area of `polygon` within the box's plan: the polygon clipped by
/// each side of the rectangle in turn (Sutherland-Hodgman, exact for a
/// convex window), then the shoelace.
fn area_in_plan(polygon: &[Point2]) -> Scalar {
    let (lo, hi) = PLAN;
    let sides: [(Point2, Point2); 4] = [
        (Point2::X, lo),
        (-Point2::X, hi),
        (Point2::Y, lo),
        (-Point2::Y, hi),
    ];
    let mut ring = polygon.to_vec();
    for (inward, at) in sides {
        let inside = |p: Point2| (p - at).dot(inward);
        let mut out = Vec::new();
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let (da, db) = (inside(a), inside(b));
            if da >= 0.0 {
                out.push(a);
            }
            if (da >= 0.0) != (db >= 0.0) {
                out.push(a + (b - a) * (da / (da - db)));
            }
        }
        ring = out;
    }
    0.5 * (0..ring.len())
        .map(|i| ring[i].perp_dot(ring[(i + 1) % ring.len()]))
        .sum::<Scalar>()
        .abs()
}

fn mesh_volume(mesh: &TriMesh) -> Scalar {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            a.dot(b.cross(c))
        })
        .sum::<Scalar>()
        / 6.0
}

/// Both compilers' volumes of the clipped box.
fn volumes(frame: Transform3, axis: PrismAxis, normal: Vec3) -> (Scalar, Scalar) {
    let (graph, root) = clipped(frame, axis, normal);
    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options())
        .unwrap_or_else(|error| panic!("the exact compiler builds it: {error:?}"));
    let exact = exact_properties(&exact, Tolerance::METRE)
        .expect("measurable")
        .signed_volume;
    let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh(&graph, root, &options())
        .expect("the mesh compiler clips");
    (exact, mesh_volume(&mesh))
}

/// One frame: its columns and origin, and the plan footprint each axis
/// gives its boundary.
struct Case {
    name: &'static str,
    x: Vec3,
    y: Vec3,
    z: Vec3,
    origin: Point2,
    /// The footprint's `y` along the plane normal, where it differs from
    /// the frame's own.
    plane_normal_y: Option<Point2>,
}

fn cases() -> Vec<Case> {
    let case = |name, x, y, z, origin| Case {
        name,
        x,
        y,
        z,
        origin,
        plane_normal_y: None,
    };
    vec![
        case("unit", Vec3::X, Vec3::Y, Vec3::Z, Point2::new(-1.0, 2.5)),
        // Left-handed: the L mirrored in its x axis, wholly in the plan
        // where the unmirrored one overhangs it.
        case(
            "mirrored",
            Vec3::X,
            -Vec3::Y,
            Vec3::Z,
            Point2::new(-1.0, 2.5),
        ),
        case(
            "mirrored in x",
            -Vec3::X,
            Vec3::Y,
            Vec3::Z,
            Point2::new(2.5, 2.5),
        ),
        case(
            "scaled",
            Vec3::X * 1.5,
            Vec3::Y * 0.5,
            Vec3::Z,
            Point2::new(0.0, -0.5),
        ),
        case(
            "scaled and mirrored",
            Vec3::X * 1.5,
            Vec3::Y * -0.75,
            Vec3::Z * 3.0,
            Point2::new(0.0, -2.0),
        ),
        case(
            "sheared",
            Vec3::X,
            Vec3::new(0.5, 1.0, 0.0),
            Vec3::Z,
            Point2::new(1.0, -1.0),
        ),
        // Right-handed with z against the normal: along the placement's z
        // the L is as written, mirrored in the plan; along the normal it
        // is turned over, unmirrored, as it always was.
        Case {
            plane_normal_y: Some(Point2::Y),
            ..case(
                "turned over",
                Vec3::X,
                -Vec3::Y,
                -Vec3::Z,
                Point2::new(-1.0, 2.5),
            )
        },
    ]
}

#[test]
fn a_mirrored_scaled_or_sheared_frame_places_its_boundary_in_both_compilers() {
    let flat = |v: Vec3| Point2::new(v.x, v.y);
    for case in cases() {
        let frame = Transform3::from_cols(
            case.x,
            case.y,
            case.z,
            Point3::new(case.origin.x, case.origin.y, 0.25),
        );
        for (axis, y) in [
            (PrismAxis::PlacementZ, flat(case.y)),
            (
                PrismAxis::PlaneNormal,
                case.plane_normal_y.unwrap_or(flat(case.y)),
            ),
        ] {
            let removed = 0.5 * area_in_plan(&footprint(case.origin, flat(case.x), y));
            let expected = 48.0 - removed;
            let (exact, mesh) = volumes(frame, axis, Vec3::Z);
            assert!(
                (exact - expected).abs() < 1e-9,
                "{} {axis:?}: exact {exact}, expected {expected}",
                case.name
            );
            assert!(
                (mesh - expected).abs() < 1e-9,
                "{} {axis:?}: mesh {mesh}, expected {expected}",
                case.name
            );
        }
        // Along the placement's z the plane's own orientation does not
        // matter: the same side, stated from below.
        let removed = 0.5 * area_in_plan(&footprint(case.origin, flat(case.x), flat(case.y)));
        let (exact, mesh) = volumes(frame, PrismAxis::PlacementZ, -Vec3::Z);
        for (path, volume) in [("exact", exact), ("mesh", mesh)] {
            assert!(
                (volume - (48.0 - removed)).abs() < 1e-9,
                "{} from below, {path}: {volume}",
                case.name
            );
        }
    }
}

/// The cases tell their frames apart: no two footprints remove the same
/// volume where they should differ, so a frame read wrongly shows.
#[test]
fn the_cases_tell_a_mirror_and_a_scale_apart() {
    let flat = |v: Vec3| Point2::new(v.x, v.y);
    let unit = area_in_plan(&footprint(Point2::new(-1.0, 2.5), Point2::X, Point2::Y));
    let mirrored = area_in_plan(&footprint(Point2::new(-1.0, 2.5), Point2::X, -Point2::Y));
    assert!((unit - 1.0).abs() < 1e-12 && (mirrored - 3.0).abs() < 1e-12);
    for case in cases() {
        let written = area_in_plan(&footprint(case.origin, flat(case.x), flat(case.y)));
        let unit_axes = area_in_plan(&footprint(
            case.origin,
            flat(case.x.normalize()),
            flat(Vec3::Z.cross(case.x.normalize())),
        ));
        assert!(
            case.name == "unit" || case.name == "turned over" || (written - unit_axes).abs() > 0.1,
            "{}: {written} vs {unit_axes}",
            case.name
        );
    }
}

/// Only a degenerate frame is refused, by name, in both compilers: one
/// whose axes are parallel or zero.
#[test]
fn a_degenerate_frame_is_refused_by_name_in_both_compilers() {
    for (x, y, z) in [
        (Vec3::X, Vec3::X * 2.0, Vec3::Z),
        (Vec3::X, Vec3::Y, Vec3::ZERO),
        (Vec3::ZERO, Vec3::Y, Vec3::Z),
    ] {
        let frame = Transform3::from_cols(x, y, z, Point3::ZERO);
        for axis in [PrismAxis::PlaneNormal, PrismAxis::PlacementZ] {
            let (graph, root) = clipped(frame, axis, Vec3::Z);
            let exact = ReferenceExactCompiler::new()
                .compile_exact(&graph, root, &options())
                .expect_err("the exact compiler refuses");
            let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
                .compile_mesh(&graph, root, &options())
                .expect_err("the mesh compiler refuses");
            for error in [exact, mesh] {
                assert!(
                    matches!(&error, GeomError::InvalidInput(m) if m.contains("degenerate")),
                    "{axis:?} {x} {y} {z}: {error:?}"
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// A boundary with arcs under a scaled frame
// ---------------------------------------------------------------------------

/// Wall length, thickness and height, the clip plane's height, the arc
/// radius as authored.
const WALL: (Scalar, Scalar, Scalar) = (6.0, 0.3, 3.0);
const ZC: Scalar = 2.0;
const R: Scalar = 1.2;

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

fn arc(from: Scalar, to: Scalar) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: R,
        }),
        domain: Interval::new(from.min(to), from.max(to)),
        same_sense: from < to,
    }
}

/// The rectangle `[0, 4] x [-1, 2]` less the disk of radius `R` about the
/// origin, run clockwise.
fn bitten() -> Profile {
    let foot = (R * R - 1.0).sqrt();
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            line(Point2::new(foot, -1.0), Point2::new(4.0, -1.0)),
            line(Point2::new(4.0, -1.0), Point2::new(4.0, 2.0)),
            line(Point2::new(4.0, 2.0), Point2::new(0.0, 2.0)),
            line(Point2::new(0.0, 2.0), Point2::new(0.0, R)),
            arc(FRAC_PI_2, 0.0),
            arc(0.0, -(1.0 / R).asin()),
        ]),
        holes: Vec::new(),
    })
}

/// The wall less what lies above `z = ZC` in the bite placed by `frame`.
fn bitten_wall(frame: Transform3, axis: PrismAxis, length: Scalar) -> (GeometryGraph, NodeId) {
    let (_, thickness, height) = WALL;
    let mut g = GeometryGraphBuilder::new();
    let mut push = |node| g.push(node).expect("a valid node");
    let profile = push(GeometryNode::Profile(Profile::Rectangle(
        RectangleProfile {
            x: length,
            y: thickness,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        },
    )));
    let wall = push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
        profile,
        direction: Vec3::Z,
        depth: height,
    }));
    let half_space = push(GeometryNode::HalfSpace(HalfSpace {
        boundary: Plane3 {
            origin: Point3::new(0.0, 0.0, ZC),
            normal: Vec3::Z,
        },
        agreement: true,
    }));
    let boundary = push(GeometryNode::Profile(bitten()));
    let tool = push(GeometryNode::SolidOperation(
        axis.operation(half_space, boundary, frame),
    ));
    let root = push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left: wall,
        right: tool,
        operator: BooleanOperator::Difference,
    }));
    (g.finish(vec![root]).expect("a valid graph"), root)
}

/// The wall less the bite scaled by `scale`: the plan right of the arc of
/// radius `scale R` to the wall's end, over the wall's thickness, the bite
/// symmetric about the wall's axis so a mirror keeps it.
fn bitten_volume(scale: Scalar) -> Scalar {
    let (length, thickness, height) = WALL;
    let r = scale * R;
    let f = |u: Scalar| 0.5 * u * (r * r - u * u).sqrt() + 0.5 * r * r * (u / r).asin();
    let strip = f(thickness / 2.0) - f(-thickness / 2.0);
    length * thickness * height - (height - ZC) * (length / 2.0 * thickness - strip)
}

/// A uniformly scaled frame scales the arcs, which stay arcs: exact in the
/// exact compiler (a derived profile, reflected where the frame mirrors)
/// and within the chord budget in the mesh compiler, along either axis. A
/// frame scaled unevenly makes them ellipses, which the exact compiler
/// refuses by name and the mesh compiler cuts.
#[test]
fn a_scaled_frame_scales_the_arcs_of_its_boundary() {
    let (_, thickness, height) = WALL;
    let scale = 0.8;
    for (name, y) in [("scaled", Vec3::Y), ("scaled and mirrored", -Vec3::Y)] {
        let frame = Transform3::from_cols(
            Vec3::X * scale,
            y * scale,
            Vec3::Z,
            Point3::new(0.0, 0.0, 5.0),
        );
        for axis in [PrismAxis::PlaneNormal, PrismAxis::PlacementZ] {
            let (graph, root) = bitten_wall(frame, axis, WALL.0);
            let expected = bitten_volume(scale);
            let exact = ReferenceExactCompiler::new()
                .compile_exact(&graph, root, &options())
                .unwrap_or_else(|error| panic!("{name} {axis:?}: {error:?}"));
            let measured = exact_properties(&exact, Tolerance::METRE)
                .expect("measurable")
                .signed_volume;
            assert!(
                (measured - expected).abs() <= 1e-9 * expected,
                "{name} {axis:?}: exact {measured}, expected {expected}"
            );
            let mesh = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
                .compile_mesh(&graph, root, &options())
                .expect("the mesh compiler clips");
            let volume = mesh_volume(&mesh);
            assert!(
                (volume - expected).abs() <= 2.0 * 1e-3 * thickness * height,
                "{name} {axis:?}: mesh {volume}, expected {expected}"
            );
        }
    }
    let uneven = Transform3::from_cols(Vec3::X * 0.8, Vec3::Y * 0.9, Vec3::Z, Point3::ZERO);
    for axis in [PrismAxis::PlaneNormal, PrismAxis::PlacementZ] {
        let (graph, root) = bitten_wall(uneven, axis, WALL.0);
        let error = ReferenceExactCompiler::new()
            .compile_exact(&graph, root, &options())
            .expect_err("an ellipse is refused");
        assert!(
            matches!(&error, GeomError::UnsupportedInput { .. }),
            "{axis:?}: {error:?}"
        );
        ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh(&graph, root, &options())
            .expect("the mesh compiler clips");
    }
}

/// A frame that enlarges its boundary chords its arcs finer, by its scale,
/// so the chords stay within the budget once placed: the mesh's deviation,
/// certified against the exact result, is within the budget.
#[test]
fn an_enlarging_frame_chords_its_arcs_within_the_budget() {
    let budget = 1e-3;
    let frame = Transform3::from_cols(Vec3::X * 3.0, Vec3::Y * 3.0, Vec3::Z, Point3::ZERO);
    for axis in [PrismAxis::PlaneNormal, PrismAxis::PlacementZ] {
        let (graph, root) = bitten_wall(frame, axis, 20.0);
        let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_with_deviation(&graph, root, &options())
            .expect("the clipped wall compiles");
        let bound = report.bound.expect("a certified bound");
        eprintln!("{axis:?}: bound {bound:e} for a budget of {budget:e}");
        assert!(report.meets_requested(), "{axis:?}: {report:?}");
    }
}
