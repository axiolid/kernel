//! A bounded half-space whose prism runs along its placement's `z` (#321).
//!
//! `SolidOperation::BoundedHalfSpaceAlongZ` bounds the kept side of the plane by the
//! infinite prism of the boundary along the placement's `z`, at any angle
//! to the plane: IFC's `IfcPolygonalBoundedHalfSpace`, bounded by its
//! `Position`. Both compilers build it. The cases:
//!
//! - **The gable wall** (openbimrs/ifc#515): the IFC2X3 wall of
//!   `gable_wall_polygonal_clip_ifc2x3.ifc`, its numbers as written, clipped
//!   by four sloped planes, each bounded by a rectangle in an axis-aligned
//!   `Position`, so every prism is vertical and oblique to its plane. It is
//!   49.999912042537 m³, the elevation polygon (0, 0), (9775, 0),
//!   (9775, 11809.43), (9299.31, 13536.68), (4925, 16165.03),
//!   (550.69, 13536.68), (0, 11537.10) times 350 mm, which IfcOpenShell
//!   0.8.5 also computes. Along the plane normal the four strips land beside
//!   the wall, and nothing is cut.
//! - **A sloped plane over a vertical prism** of a box, closed form, in both
//!   compilers, whatever the sign or length of the placement's `z`.
//! - **An oblique prism with arcs**: a wall under a sloped plane cut by the
//!   vertical prism of a boundary with a circular bite: the exact compiler
//!   cuts each arc's right cylinder wall in an ellipse and matches the
//!   closed form, and the mesh compiler's deviation is certified against it.
//! - **A leaning prism**, its placement sheared so `z` leans off the
//!   boundary's plane: a polygon in both compilers, closed form; arcs, whose
//!   walls are elliptical cylinders (#287), in the mesh compiler, the exact
//!   one building them or refusing by name.
//! - **Refusals**: an axis in the plane, and a degenerate placement, by name
//!   in both compilers.

// The gable's numbers are the IFC file's, digit for digit, beyond what an
// f64 holds: they read as the nearest double, as the file does.
#![allow(clippy::excessive_precision)]

use std::f64::consts::FRAC_PI_2;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{
    BooleanOperator, Frame2, Interval, Plane3, Point2, Point3, Scalar, Tolerance, Transform3, Vec2,
    Vec3,
};
use axiolid_curve::{Circle2, Curve2, Line2, Polyline2};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_measure::exact_properties;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    deviation::COINCIDENT_OPERANDS, DeviationBound, DeviationPath, ReferenceExactCompiler,
    ReferenceMeshCompiler,
};
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation,
};
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

    /// A `x` by `y` rectangle about the origin extruded `depth` along `z`.
    fn block(&mut self, (x, y): (Scalar, Scalar), depth: Scalar) -> NodeId {
        let profile = self.push(GeometryNode::Profile(rect(x, y)));
        self.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth,
        }))
    }

    fn instance(&mut self, source: NodeId, transform: Transform3) -> NodeId {
        self.push(GeometryNode::Instance(Instance { source, transform }))
    }

    /// The side of the plane through `origin` with `normal` that
    /// `agreement` keeps, bounded by `boundary` in `frame` along `axis`.
    fn bounded(
        &mut self,
        (origin, normal, agreement): (Point3, Vec3, bool),
        boundary: GeometryNode,
        frame: Transform3,
        axis: PrismAxis,
    ) -> NodeId {
        let half_space = self.push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 { origin, normal },
            agreement,
        }));
        let boundary = self.push(boundary);
        self.push(GeometryNode::SolidOperation(
            axis.operation(half_space, boundary, frame),
        ))
    }

    fn boolean(&mut self, left: NodeId, right: NodeId, operator: BooleanOperator) -> NodeId {
        self.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left,
            right,
            operator,
        }))
    }

    fn finish(self, root: NodeId) -> (GeometryGraph, NodeId) {
        (
            self.builder.finish(vec![root]).expect("a valid graph"),
            root,
        )
    }
}

fn rect(x: Scalar, y: Scalar) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// A closed polyline through `points`, as IFC writes one: the first point
/// repeated at the end.
fn polyline(points: &[(Scalar, Scalar)]) -> GeometryNode {
    let mut points: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    points.push(points[0]);
    GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
        points,
        closed: true,
    }))
}

/// A frame translated to `origin`, its axes the world's.
fn at(origin: Point3) -> Transform3 {
    Transform3::from_translation(origin)
}

fn exact_volume(graph: &GeometryGraph, root: NodeId, options: &ExecutionOptions) -> Scalar {
    let body = ReferenceExactCompiler::new()
        .compile_exact(graph, root, options)
        .unwrap_or_else(|error| panic!("the exact compiler builds it: {error:?}"));
    let topology = axiolid_topology::audit_brep(body.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    exact_properties(&body, options.tolerance())
        .expect("measurable")
        .signed_volume
}

fn mesh_volume(mesh: &TriMesh) -> Scalar {
    let base = mesh.positions[0];
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|i| mesh.positions[t[i] as usize] - base);
            a.dot(b.cross(c))
        })
        .sum::<Scalar>()
        / 6.0
}

fn mesh_compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

// ---------------------------------------------------------------------------
// The gable wall (openbimrs/ifc#515)
// ---------------------------------------------------------------------------

/// IfcOpenShell 0.8.5's volume of the gable wall, in m³.
const GABLE: Scalar = 49.999_912_042_537;

/// The uncut extrusion, 9.775 x 0.35 x 16.165033 m.
const UNCUT: Scalar = 55.304_619_372_356;

/// One `IfcPolygonalBoundedHalfSpace` of the wall: `BaseSurface` through
/// `Location` with normal `Axis`, `AgreementFlag` `.T.` (the material of
/// the half-space solid opposite the normal, so the kernel keeps the side
/// against it), `Position` at `Location` with the world's axes, and the
/// boundary `(0, 0), (wide.0, 0), (wide.1, 350), (0, 350)`, a rectangle to
/// the file's rounding.
struct Clip {
    location: Point3,
    axis: Vec3,
    wide: (Scalar, Scalar),
}

/// The wall's four clips, `#8990`, `#9002`, `#9014` and `#9026`, their
/// numbers as the file writes them.
fn gable_clips() -> [Clip; 4] {
    [
        Clip {
            location: Point3::new(0.0, -175.0, 11537.100698716456),
            axis: Vec3::new(0.96410569120432499, 0.0, -0.26551876805120728),
            wide: (550.69387809692876, 550.69387809692739),
        },
        Clip {
            location: Point3::new(550.69387809692716, -175.0, 13536.684780404561),
            axis: Vec3::new(0.51503807491005349, 0.0, -0.85716730070211278),
            wide: (4374.3061219031079, 4374.3061219031097),
        },
        Clip {
            location: Point3::new(4925.0, -175.0, 16165.033064627303),
            axis: Vec3::new(-0.51503807491005638, 0.0, -0.85716730070211111),
            wide: (4374.3061219030942, 4374.306121903096),
        },
        Clip {
            location: Point3::new(9299.3061219031333, -175.0, 13536.684780404546),
            axis: Vec3::new(-0.96410569120432588, -0.0, -0.26551876805120411),
            wide: (475.6938780968685, 475.69387809687032),
        },
    ]
}

/// The wall `#8979` minus its four clips under `axis`, placed by the
/// wall's `IfcLocalPlacement` `#8971`, in millimetres.
fn gable_wall(axis: PrismAxis) -> (GeometryGraph, NodeId) {
    let mut g = Graph::new();
    let block = g.block((9775.0000000000036, 350.00000000000057), 16165.033064627303);
    // `#8976`: centred at (4887.5, 5.7e-14), its x along -X.
    let mut body = g.instance(
        block,
        Transform3::from_cols(
            -Vec3::X,
            -Vec3::Y,
            Vec3::Z,
            Vec3::new(4887.5000000000018, 5.6843418860808015e-14, 0.0),
        ),
    );
    for clip in gable_clips() {
        let boundary = polyline(&[
            (0.0, 0.0),
            (clip.wide.0, 0.0),
            (clip.wide.1, 350.0),
            (0.0, 350.0),
        ]);
        let tool = g.bounded(
            (clip.location, clip.axis, false),
            boundary,
            at(clip.location),
            axis,
        );
        body = g.boolean(body, tool, BooleanOperator::Difference);
    }
    // `#8970`: at (-8071.07, -722.43), its x along +Y.
    let placed = g.instance(
        body,
        Transform3::from_cols(
            Vec3::Y,
            -Vec3::X,
            Vec3::Z,
            Vec3::new(-8071.0663924178298, -722.43106880235769, 0.0),
        ),
    );
    g.finish(placed)
}

fn millimetres() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
}

/// Each `Position` is axis-aligned and each plane sloped, so each prism is
/// vertical and oblique to its plane: the gable, exactly.
#[test]
fn the_gable_wall_is_clipped_to_its_gable_by_the_exact_compiler() {
    let (graph, root) = gable_wall(PrismAxis::PlacementZ);
    let volume = exact_volume(&graph, root, &millimetres()) * 1e-9;
    eprintln!("gable, exact: {volume:.12} m³");
    assert!(
        (volume - GABLE).abs() < 1e-6,
        "exact volume {volume} m³, expected the gable's {GABLE}"
    );
}

/// The mesh compiler cuts the same gable, each slab sized past the wall
/// (#325): sized by its boundary alone, the first clip's slab, 1.1 m deep
/// for a 0.55 m boundary, ended 0.13 m short of the gable's top across the
/// plane, and left a sliver. Its mesh is certified against the exact
/// result.
#[test]
fn the_gable_wall_is_clipped_to_its_gable_by_the_mesh_compiler() {
    let (graph, root) = gable_wall(PrismAxis::PlacementZ);
    let (outcome, report) = mesh_compiler()
        .compile_mesh_with_deviation(&graph, root, &millimetres())
        .expect("the gable compiles");
    let volume = mesh_volume(&outcome.mesh) * 1e-9;
    eprintln!("gable, mesh: {volume:.12} m³");
    assert!(
        (volume - GABLE).abs() < 1e-6,
        "mesh volume {volume} m³, expected the gable's {GABLE}"
    );
    assert!(report.meets_requested(), "{report:?}");
}

/// Along the plane normal the same clips are the strips IFC does not
/// mean: each lands beside the wall, and nothing is cut.
#[test]
fn along_the_plane_normal_the_gable_clips_miss_the_wall() {
    let (graph, root) = gable_wall(PrismAxis::PlaneNormal);
    let volume = exact_volume(&graph, root, &millimetres()) * 1e-9;
    assert!(
        (volume - UNCUT).abs() < 1e-6,
        "exact volume {volume} m³, expected the uncut {UNCUT}"
    );
}

// ---------------------------------------------------------------------------
// A sloped plane over a vertical prism
// ---------------------------------------------------------------------------

/// The box `[-2, 2] x [-3, 3] x [0, 2]` less the side above the plane
/// `z = 1.5 - 0.75 x` within the vertical prism of `[1, 4] x [0, 0.5]`
/// placed at `(0, 2.75)`: the box's `x in [1, 2], y in [2.75, 3]` above
/// the plane, `0.25 * integral over [1, 2] of (0.5 + 0.75 x) dx`.
fn sloped_box(operator: BooleanOperator, frame: Transform3) -> (GeometryGraph, NodeId) {
    let mut g = Graph::new();
    let block = g.block((4.0, 6.0), 2.0);
    let tool = g.bounded(
        (Point3::new(0.0, 0.0, 1.5), Vec3::new(0.6, 0.0, 0.8), true),
        polyline(&[(1.0, 0.0), (4.0, 0.0), (4.0, 0.5), (1.0, 0.5)]),
        frame,
        PrismAxis::PlacementZ,
    );
    let root = g.boolean(block, tool, operator);
    g.finish(root)
}

const BOX: Scalar = 48.0;
const SLOPED_CUT: Scalar = 0.406_25;

fn metres(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
        .with_chord_error(budget)
        .expect("a positive budget")
}

#[test]
fn a_sloped_plane_over_a_vertical_prism_cuts_the_prism_in_both_compilers() {
    let options = metres(1e-3);
    for (operator, expected) in [
        (BooleanOperator::Difference, BOX - SLOPED_CUT),
        (BooleanOperator::Intersection, SLOPED_CUT),
    ] {
        // Raised along its own z, the prism is the same: it is infinite
        // both ways. So is it with z reversed, or stretched: only its
        // direction counts, and the boundary is where x and y put it.
        for (case, z, height) in [
            ("up", Vec3::Z, 0.0),
            ("raised", Vec3::Z, 7.5),
            ("down", -Vec3::Z, -3.0),
            ("stretched", Vec3::Z * 4.0, 1.0),
        ] {
            let frame = Transform3::from_cols(Vec3::X, Vec3::Y, z, Point3::new(0.0, 2.75, height));
            let (graph, root) = sloped_box(operator, frame);
            let exact = exact_volume(&graph, root, &options);
            assert!(
                (exact - expected).abs() < 1e-9,
                "{operator:?} {case}: exact {exact}, expected {expected}"
            );
            let mesh = mesh_compiler()
                .compile_mesh(&graph, root, &options)
                .expect("the mesh compiler clips too");
            let measured = mesh_volume(&mesh);
            assert!(
                (measured - expected).abs() < 1e-9,
                "{operator:?} {case}: mesh {measured}, expected {expected}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// An oblique prism with arcs
// ---------------------------------------------------------------------------

/// Wall length, thickness and height, the clip plane's height at the
/// arc's centre and its slope along the wall, and the arc radius.
const L: Scalar = 6.0;
const T: Scalar = 0.3;
const H: Scalar = 3.0;
const ZC: Scalar = 2.0;
const SLOPE: Scalar = 0.1;
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

fn arc(centre: Point2, (from, to): (Scalar, Scalar)) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
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
/// origin, the arc split where it crosses `y = 0`, run clockwise.
fn bitten() -> Profile {
    let c = Point2::ZERO;
    let foot = (R * R - 1.0).sqrt();
    let below = -(1.0 / R).asin();
    Profile::Contour(ContourProfile {
        outer: Contour::new(vec![
            line(Point2::new(foot, -1.0), Point2::new(4.0, -1.0)),
            line(Point2::new(4.0, -1.0), Point2::new(4.0, 2.0)),
            line(Point2::new(4.0, 2.0), Point2::new(0.0, 2.0)),
            line(Point2::new(0.0, 2.0), Point2::new(0.0, R)),
            arc(c, (FRAC_PI_2, 0.0)),
            arc(c, (0.0, below)),
        ]),
        holes: Vec::new(),
    })
}

/// `integral of sqrt(R^2 - y^2) dy` over the wall's thickness `|y| <= T/2`.
fn disk_strip() -> Scalar {
    let f = |u: Scalar| 0.5 * u * (R * R - u * u).sqrt() + 0.5 * R * R * (u / R).asin();
    f(T / 2.0) - f(-T / 2.0)
}

/// The wall minus the side above the plane `z = ZC + SLOPE x` within the
/// vertical prism of the bite, all under `placement`. The prism stands
/// oblique to the plane, each arc a vertical cylinder the plane cuts in an
/// ellipse. The wall's plan the bite covers is `x` from the arc,
/// `sqrt(R^2 - y^2)`, to the wall's end `L/2`, over `|y| <= T/2`, and the
/// wall above the plane there is `H - ZC - SLOPE x` high, so the cut is
/// `(H - ZC) (L T / 2 - disk_strip) - SLOPE/2 ((L^2/4 - R^2) T + T^3/12)`.
fn sloped_bite(placement: Transform3) -> (GeometryGraph, NodeId) {
    let mut g = Graph::new();
    let wall = g.block((L, T), H);
    let wall = g.instance(wall, placement);
    let tool = g.bounded(
        (Point3::new(0.0, 0.0, ZC), Vec3::new(-SLOPE, 0.0, 1.0), true),
        GeometryNode::Profile(bitten()),
        at(Point3::new(0.0, 0.0, -4.0)),
        PrismAxis::PlacementZ,
    );
    let tool = g.instance(tool, placement);
    let root = g.boolean(wall, tool, BooleanOperator::Difference);
    g.finish(root)
}

fn sloped_bite_volume() -> Scalar {
    let cut = (H - ZC) * (L * T / 2.0 - disk_strip())
        - SLOPE / 2.0 * ((L * L / 4.0 - R * R) * T + T * T * T / 12.0);
    L * T * H - cut
}

fn general() -> Transform3 {
    Transform3::from_translation(Vec3::new(-2.0, 7.0, 1.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7)
}

#[test]
fn an_oblique_prism_with_arcs_is_exact_and_the_mesh_is_certified_against_it() {
    let expected = sloped_bite_volume();
    for placement in [Transform3::IDENTITY, general()] {
        let (graph, root) = sloped_bite(placement);
        let exact = ReferenceExactCompiler::new()
            .compile_exact(&graph, root, &metres(1e-3))
            .unwrap_or_else(|error| panic!("the exact compiler builds it: {error:?}"));
        let measured = exact_properties(&exact, Tolerance::METRE)
            .expect("measurable")
            .signed_volume;
        assert!(
            (measured - expected).abs() <= 1e-9 * expected,
            "exact volume {measured}, expected {expected}"
        );
        // Each arc is a vertical round wall, not chords.
        assert!(
            exact
                .surfaces()
                .iter()
                .any(|s| matches!(s, axiolid_surface::Surface::Cylinder(_))),
            "the arc must cut a cylindrical face"
        );
        for budget in [1e-2, 1e-3] {
            let (outcome, report) = mesh_compiler()
                .compile_mesh_with_deviation(&graph, root, &metres(budget))
                .expect("the clipped wall compiles");
            // Besides the measurement, at most how far faces read as one
            // surface within tolerance moved an operand (#324).
            let measured: Vec<_> = report
                .contributions
                .iter()
                .filter(|c| c.detail != COINCIDENT_OPERANDS)
                .collect();
            assert_eq!(measured.len(), 1, "{report:?}");
            assert!(report.meets_requested(), "{report:?}");
            let contribution = *measured[0];
            assert_eq!(contribution.path, DeviationPath::Boolean);
            let DeviationBound::Certified(bound) = contribution.bound else {
                panic!("expected a certified bound, got {report:?}");
            };
            assert!(
                bound <= 2.0 * budget,
                "bound {bound} for a budget of {budget}"
            );
            let slack = 2.0 * budget * T * H + 1e-9;
            let volume = mesh_volume(&outcome.mesh);
            assert!(
                (volume - expected).abs() <= slack,
                "mesh volume {volume}, exact {expected}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// A sheared placement: the prism leaning off the boundary's own plane
// ---------------------------------------------------------------------------

/// The lean of a sheared placement's `z` along `+x` per unit height.
const LEAN: Scalar = 0.3;

/// The sloped box clipped by its rectangle under a placement whose `z`
/// leans along `+x`: at height `z` the prism's section is the rectangle
/// moved by `LEAN z`, so over `x in [1, 2]` of the box the cut runs from
/// the plane `z = 1.5 - 0.75 x` up to `z = (x - 1) / LEAN`, capped at the
/// box's top `z = 2`, over the quarter metre of `y`.
#[test]
fn a_leaning_prism_cuts_its_own_sections_in_both_compilers() {
    // Two linear pieces of the cut height, exact by the trapezoid rule.
    let (start, knee) = (1.45 / (1.0 + 0.75 * LEAN), 1.0 + 2.0 * LEAN);
    let rising = |x: Scalar| (x - 1.0) / LEAN - 1.5 + 0.75 * x;
    let capped = |x: Scalar| 0.5 + 0.75 * x;
    let cut = 0.25
        * ((knee - start) * (rising(start) + rising(knee)) / 2.0
            + (2.0 - knee) * (capped(knee) + capped(2.0)) / 2.0);
    // The same prism with its z reversed: it is infinite both ways.
    for z in [Vec3::new(LEAN, 0.0, 1.0), Vec3::new(-LEAN, 0.0, -1.0)] {
        let frame = Transform3::from_cols(Vec3::X, Vec3::Y, z, Point3::new(0.0, 2.75, 0.0));
        let (graph, root) = sloped_box(BooleanOperator::Difference, frame);
        let options = metres(1e-3);
        let exact = exact_volume(&graph, root, &options);
        assert!(
            (exact - (BOX - cut)).abs() < 1e-9,
            "z {z}: exact {exact}, expected {}",
            BOX - cut
        );
        let mesh = mesh_compiler()
            .compile_mesh(&graph, root, &options)
            .expect("the mesh compiler clips too");
        assert!((mesh_volume(&mesh) - (BOX - cut)).abs() < 1e-9, "z {z}");
    }
}

/// The bite under a level plane, its prism leaning along the wall: at
/// height `z` the section is the bite moved by `LEAN z`, so the plan it
/// covers loses `LEAN z T` to the arc's move.
fn leaning_bite() -> (GeometryGraph, NodeId) {
    let mut g = Graph::new();
    let wall = g.block((L, T), H);
    let frame = Transform3::from_cols(Vec3::X, Vec3::Y, Vec3::new(LEAN, 0.0, 1.0), Point3::ZERO);
    let tool = g.bounded(
        (Point3::new(0.0, 0.0, ZC), Vec3::Z, true),
        GeometryNode::Profile(bitten()),
        frame,
        PrismAxis::PlacementZ,
    );
    let root = g.boolean(wall, tool, BooleanOperator::Difference);
    g.finish(root)
}

/// A leaning prism of arcs sweeps elliptical cylinders. The mesh compiler
/// cuts it within its chord budget; the exact compiler builds it or, where
/// the general boolean cannot sew the cut, refuses by name: never a
/// different solid.
#[test]
fn a_leaning_prism_with_arcs_is_cut_or_refused_by_name() {
    let cut = (H - ZC) * (L * T / 2.0 - disk_strip()) - LEAN * T * (H * H - ZC * ZC) / 2.0;
    let expected = L * T * H - cut;
    let (graph, root) = leaning_bite();
    let budget = 1e-3;
    let mesh = mesh_compiler()
        .compile_mesh(&graph, root, &metres(budget))
        .expect("the mesh compiler clips");
    let volume = mesh_volume(&mesh);
    assert!(
        (volume - expected).abs() <= 2.0 * budget * T * H,
        "mesh volume {volume}, expected {expected}"
    );
    match ReferenceExactCompiler::new().compile_exact(&graph, root, &metres(budget)) {
        Ok(exact) => {
            let measured = exact_properties(&exact, Tolerance::METRE)
                .expect("measurable")
                .signed_volume;
            assert!(
                (measured - expected).abs() <= 1e-9 * expected,
                "exact volume {measured}, expected {expected}"
            );
        }
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert!(input.contains("exact boolean"), "{input}");
        }
        Err(error) => panic!("a named refusal, not {error:?}"),
    }
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

fn refused_by_both(frame: Transform3, needle: &str) {
    let (graph, root) = sloped_box(BooleanOperator::Difference, frame);
    let options = metres(1e-3);
    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options)
        .expect_err("the exact compiler refuses");
    let mesh = mesh_compiler()
        .compile_mesh(&graph, root, &options)
        .expect_err("the mesh compiler refuses");
    for error in [exact, mesh] {
        assert!(
            matches!(&error, GeomError::InvalidInput(m) if m.contains(needle)),
            "expected a refusal naming {needle:?}, got {error:?}"
        );
    }
}

/// A `z` in the sloped plane never crosses it, nor one within the angular
/// tolerance of it; a placement whose axes do not span space has no prism.
#[test]
fn an_axis_in_the_plane_or_a_degenerate_placement_is_refused_by_name() {
    let along_slope = Vec3::new(0.8, 0.0, -0.6);
    for z in [along_slope, along_slope + Vec3::new(0.6, 0.0, 0.8) * 1e-12] {
        // The boundary stands upright, across the slope, so the placement
        // itself is sound.
        let frame = Transform3::from_cols(Vec3::Y, Vec3::new(0.6, 0.0, 0.8), z, Point3::ZERO);
        refused_by_both(frame, "parallel to the clip plane");
    }
    let flat = Transform3::from_cols(Vec3::X, Vec3::Y, Vec3::new(1.0, 1.0, 0.0), Point3::ZERO);
    refused_by_both(flat, "degenerate");
}
