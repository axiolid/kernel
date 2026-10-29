//! A boolean with a warped operand reports a bound that covers the warp
//! (#266).
//!
//! A warped authored face (#254, #261) or a warped B-rep face without a
//! surface (#257) has no unique surface; its mesh is its triangles in the
//! fit plane, within the slab its corners span of every reading of it.
//! Under a boolean the report measures the mesh against the exact boolean
//! of the operands flattened to those triangles (ADR 0080: the exact
//! boolean of operands perturbed within the reported deviation) and widens
//! the bound by the slab width. Each case compiles a warped unit box (one
//! corner lifted, or a saddle) unioned with and cut by a planar box, and
//! samples the boolean of the unflattened box densely (the bilinear patch
//! through the warped face's corners, and either diagonal triangulation
//! of it), computed here from the corners: no sample lies further from
//! the mesh than the bound.

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{BooleanOperator, Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_mesh::{PolygonFace, PolygonMesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{
    DeviationBound, DeviationContribution, DeviationPath, DeviationReport, ReferenceMeshCompiler,
};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, NodeId, SolidOperation};
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_topology::{
    BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Shell, Solid, Vertex,
};

const WARPED_AUTHORED: &str = "non-planar authored face";
const WARPED_BREP: &str = "non-planar face without a surface";
const FLATTENED: &str =
    "measured against the exact boolean of operands flattened within their slab width, plus it";
const MEASURED_EXACT: &str =
    "measured against the exact compiler's result, the exact boolean of the given operands";

/// A corner of the top face lifted 5 cm over (1, 1): slab width 2.5 cm.
const LIFTED: [Scalar; 4] = [1.0, 1.0, 1.05, 1.0];
/// A saddle, corners alternately 5 cm above and below z = 1: 10 cm.
const SADDLE: [Scalar; 4] = [1.05, 0.95, 1.05, 0.95];
const FLAT: [Scalar; 4] = [1.0; 4];

/// The planar box the warped box is combined with: x in [0.5, 1.5],
/// y in [0.25, 0.75], z in [0.5, 1.5], through the warped top face and
/// the side x = 1.
const TOOL_LO: [Scalar; 3] = [0.5, 0.25, 0.5];
const TOOL_HI: [Scalar; 3] = [1.5, 0.75, 1.5];

fn options(budget: Scalar) -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(budget)
        .unwrap()
}

/// The unit cube of quads, outward, its top corners 4..8 over (0,0),
/// (1,0), (1,1), (0,1) at `top`.
fn authored_box(top: [Scalar; 4]) -> PolygonMesh {
    let plan = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let positions = (0..8)
        .map(|i| {
            let z = if i < 4 { 0.0 } else { top[i - 4] };
            Point3::new(plan[i % 4][0], plan[i % 4][1], z)
        })
        .collect();
    let faces = QUADS
        .iter()
        .map(|q| PolygonFace {
            outer: q.to_vec(),
            holes: Vec::new(),
        })
        .collect();
    PolygonMesh { positions, faces }
}

const QUADS: [[u32; 4]; 6] = [
    [0, 3, 2, 1],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 0, 4, 7],
];

/// The same box as a B-rep of faces without surfaces (#257).
fn brep_box(top: [Scalar; 4]) -> BRep<NodeId> {
    let authored = authored_box(top);
    let mut brep = BRep::default();
    let verts: Vec<_> = authored
        .positions
        .iter()
        .map(|&position| brep.add_vertex(Vertex { position }))
        .collect();
    let mut edges = std::collections::HashMap::new();
    let mut faces = Vec::new();
    for q in QUADS {
        let mut uses = Vec::new();
        for i in 0..4 {
            let (a, b) = (q[i] as usize, q[(i + 1) % 4] as usize);
            let key = (a.min(b), a.max(b));
            let edge = *edges.entry(key).or_insert_with(|| {
                brep.add_edge(Edge {
                    start: verts[key.0],
                    end: verts[key.1],
                    curve: None,
                })
            });
            uses.push(EdgeUse {
                edge,
                orientation: if a == key.0 {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                },
                pcurve: None,
            });
        }
        let loop_id = brep.add_loop(Loop { edges: uses });
        faces.push(brep.add_face(Face {
            surface: None,
            bounds: vec![FaceBound {
                loop_id,
                orientation: Orientation::Forward,
                outer: true,
            }],
            orientation: Orientation::Forward,
        }));
    }
    let shell = brep.add_shell(Shell {
        faces: faces.iter().map(|&f| (f, Orientation::Forward)).collect(),
        closed: true,
    });
    brep.add_solid(Solid {
        outer: shell,
        voids: Vec::new(),
    });
    brep
}

/// The planar box `TOOL_LO..TOOL_HI` as a placed extrusion.
fn tool(g: &mut GeometryGraphBuilder) -> NodeId {
    let size = [0, 1, 2].map(|k| TOOL_HI[k] - TOOL_LO[k]);
    let profile = g
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: size[0],
                y: size[1],
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let solid = g
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: size[2],
        }))
        .unwrap();
    let centre = Vec3::new(
        0.5 * (TOOL_LO[0] + TOOL_HI[0]),
        0.5 * (TOOL_LO[1] + TOOL_HI[1]),
        TOOL_LO[2],
    );
    g.push(GeometryNode::Instance(Instance {
        source: solid,
        transform: Transform3::from_translation(centre),
    }))
    .unwrap()
}

fn boolean(
    g: &mut GeometryGraphBuilder,
    left: NodeId,
    right: NodeId,
    operator: BooleanOperator,
) -> NodeId {
    g.push(GeometryNode::SolidOperation(SolidOperation::Boolean {
        left,
        right,
        operator,
    }))
    .unwrap()
}

/// `box` combined with the planar tool by `operator`, compiled with its
/// report.
fn compile(
    warped: GeometryNode,
    operator: BooleanOperator,
    budget: Scalar,
) -> (TriMesh, DeviationReport) {
    let mut g = GeometryGraphBuilder::new();
    let operand = g.push(warped).unwrap();
    let tool = tool(&mut g);
    let root = boolean(&mut g, operand, tool, operator);
    let graph = g.finish(vec![root]).unwrap();
    let (outcome, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options(budget))
        .expect("the boolean compiles");
    (outcome.mesh, report)
}

fn contribution<'r>(
    report: &'r DeviationReport,
    path: DeviationPath,
    detail: &str,
) -> Option<&'r DeviationContribution> {
    report
        .contributions
        .iter()
        .find(|c| c.path == path && c.detail == detail)
}

/// The width of the slab the corners span about their fit plane (the
/// centroid, the Newell normal), worked out here from the corners.
fn slab_width(top: [Scalar; 4]) -> Scalar {
    let corners = [
        Point3::new(0.0, 0.0, top[0]),
        Point3::new(1.0, 0.0, top[1]),
        Point3::new(1.0, 1.0, top[2]),
        Point3::new(0.0, 1.0, top[3]),
    ];
    let mut normal = Vec3::ZERO;
    for (i, &a) in corners.iter().enumerate() {
        normal += a.cross(corners[(i + 1) % 4]);
    }
    let n = normal.normalize();
    let centroid = corners.iter().fold(Vec3::ZERO, |s, &c| s + c) / 4.0;
    let d: Vec<Scalar> = corners.iter().map(|&c| (c - centroid).dot(n)).collect();
    d.iter().copied().fold(Scalar::MIN, Scalar::max)
        - d.iter().copied().fold(Scalar::MAX, Scalar::min)
}

/// A reading of the warped top face: its height over `(x, y)`.
#[derive(Debug, Clone, Copy)]
enum Reading {
    Bilinear,
    /// Triangles (0, 1, 2) and (0, 2, 3).
    Diagonal02,
    /// Triangles (0, 1, 3) and (1, 2, 3).
    Diagonal13,
}

impl Reading {
    fn height(self, top: [Scalar; 4], x: Scalar, y: Scalar) -> Scalar {
        let [a, b, c, d] = top;
        match self {
            Self::Bilinear => {
                a * (1.0 - x) * (1.0 - y) + b * x * (1.0 - y) + c * x * y + d * (1.0 - x) * y
            }
            // Below the diagonal y = x: corners 0, 1, 2; above: 0, 2, 3.
            Self::Diagonal02 if y <= x => a + (b - a) * x + (c - b) * y,
            Self::Diagonal02 => a + (c - d) * x + (d - a) * y,
            // Below the diagonal x + y = 1: corners 0, 1, 3; above: 1, 2, 3.
            Self::Diagonal13 if x + y <= 1.0 => a + (b - a) * x + (d - a) * y,
            Self::Diagonal13 => b + d - c + (c - d) * x + (c - b) * y,
        }
    }
}

/// Points on the boundary of the exact boolean of the warped box, read as
/// `reading`, and the tool box, on a grid of `n` steps per face side.
fn boolean_samples(
    top: [Scalar; 4],
    reading: Reading,
    operator: BooleanOperator,
    n: usize,
) -> Vec<Point3> {
    let steps = |lo: Scalar, hi: Scalar| -> Vec<Scalar> {
        (0..=n)
            .map(|i| lo + (hi - lo) * i as Scalar / n as Scalar)
            .collect()
    };
    let in_box = |p: Point3| {
        (0..3).all(|k| {
            let v = [p.x, p.y, p.z][k];
            v > TOOL_LO[k] && v < TOOL_HI[k]
        })
    };
    let in_warped = |p: Point3| {
        p.x > 0.0
            && p.x < 1.0
            && p.y > 0.0
            && p.y < 1.0
            && p.z > 0.0
            && p.z < reading.height(top, p.x, p.y)
    };
    let on_warped_or_in = |p: Point3| {
        (0.0..=1.0).contains(&p.x)
            && (0.0..=1.0).contains(&p.y)
            && p.z >= 0.0
            && p.z <= reading.height(top, p.x, p.y)
    };
    let unit = steps(0.0, 1.0);
    // The warped box's boundary.
    let mut warped = Vec::new();
    for &u in &unit {
        for &v in &unit {
            warped.push(Point3::new(u, v, 0.0));
            warped.push(Point3::new(u, v, reading.height(top, u, v)));
            for (x, y) in [(0.0, u), (1.0, u), (u, 0.0), (u, 1.0)] {
                warped.push(Point3::new(x, y, v * reading.height(top, x, y)));
            }
        }
    }
    // The tool's boundary.
    let mut tool = Vec::new();
    let axes = [0, 1, 2].map(|k| steps(TOOL_LO[k], TOOL_HI[k]));
    for k in 0..3 {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        for side in [TOOL_LO[k], TOOL_HI[k]] {
            for &a in &axes[i] {
                for &b in &axes[j] {
                    let mut p = [0.0; 3];
                    p[k] = side;
                    p[i] = a;
                    p[j] = b;
                    tool.push(Point3::new(p[0], p[1], p[2]));
                }
            }
        }
    }
    let on_box_or_in = |p: Point3| {
        (0..3).all(|k| {
            let v = [p.x, p.y, p.z][k];
            v >= TOOL_LO[k] && v <= TOOL_HI[k]
        })
    };
    // Each operand's boundary where the other leaves it on the result's.
    let (warped, tool): (Vec<Point3>, Vec<Point3>) = match operator {
        BooleanOperator::Union => (
            warped.into_iter().filter(|&p| !in_box(p)).collect(),
            tool.into_iter().filter(|&p| !in_warped(p)).collect(),
        ),
        BooleanOperator::Difference => (
            warped.into_iter().filter(|&p| !in_box(p)).collect(),
            tool.into_iter().filter(|&p| on_warped_or_in(p)).collect(),
        ),
        _ => (
            warped.into_iter().filter(|&p| on_box_or_in(p)).collect(),
            tool.into_iter().filter(|&p| on_warped_or_in(p)).collect(),
        ),
    };
    warped.into_iter().chain(tool).collect()
}

fn distance_to_mesh(p: Point3, mesh: &TriMesh) -> Scalar {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .filter_map(|t| {
            let tri = [0, 1, 2].map(|i| mesh.positions[t[i] as usize]);
            axiolid_measure::proximity::closest_point_on_triangle(p, tri)
                .ok()
                .map(|q| (q - p).length())
        })
        .fold(Scalar::INFINITY, Scalar::min)
}

/// The worst distance from the boolean of the unflattened box, under each
/// reading of its warped face, to the mesh.
fn worst_sample(top: [Scalar; 4], operator: BooleanOperator, mesh: &TriMesh) -> Scalar {
    [Reading::Bilinear, Reading::Diagonal02, Reading::Diagonal13]
        .into_iter()
        .flat_map(|reading| boolean_samples(top, reading, operator, 60))
        .map(|p| distance_to_mesh(p, mesh))
        .fold(0.0, Scalar::max)
}

/// The issue's done-when, for one warped box, operator and operand kind.
fn assert_covers_the_warp(
    top: [Scalar; 4],
    operator: BooleanOperator,
    node: GeometryNode,
    warp_detail: (DeviationPath, &str),
) {
    let width = slab_width(top);
    for budget in [1e-3, 0.05, 0.2] {
        let (mesh, report) = compile(node.clone(), operator, budget);
        // The operand's slab width stays named, at its own value.
        let warp = contribution(&report, warp_detail.0, warp_detail.1)
            .unwrap_or_else(|| panic!("the warped operand is named: {report:?}"));
        let DeviationBound::Certified(warp) = warp.bound else {
            panic!("a warp is certified: {report:?}");
        };
        assert!(
            warp >= width && warp <= width + 1e-12,
            "warp {warp}, slab {width}"
        );
        // The boolean's own bound includes it.
        let measured = contribution(&report, DeviationPath::Boolean, FLATTENED)
            .unwrap_or_else(|| panic!("measured against the flattened boolean: {report:?}"));
        let DeviationBound::Certified(measured) = measured.bound else {
            panic!("a flattened boolean is certified: {report:?}");
        };
        assert!(measured >= width, "{measured} below the slab {width}");
        // Not loose: the boxes' meshes are their exact polyhedra, so only
        // the search's resolution adds to the slab width, and a slab that
        // fits the budget leaves the measurement the rest of it.
        if width <= budget {
            assert!(measured <= budget, "{measured} over {budget}");
        } else {
            assert!(measured <= width + budget, "{measured} for a slab {width}");
        }
        let bound = report.bound.expect("bounded");
        assert!(bound >= warp && bound >= measured, "{report:?}");
        assert_eq!(report.meets_requested(), bound <= budget, "{report:?}");
        assert_eq!(report.meets_requested(), width <= budget, "{report:?}");
        assert_eq!(report.unbounded().count(), 0, "{report:?}");

        // What it bounds: the boolean of the unflattened box, under each
        // reading of its warped face, lies within it of the mesh.
        let worst = worst_sample(top, operator, &mesh);
        assert!(
            worst <= bound,
            "a sample lies {worst} from the mesh, over {bound}"
        );
        // And the readings do differ from the mesh: the bound is not idle.
        assert!(worst >= 0.2 * width, "{worst} for a slab {width}");
    }
}

#[test]
fn a_lifted_corner_unioned_with_a_planar_box_covers_its_slab() {
    assert_covers_the_warp(
        LIFTED,
        BooleanOperator::Union,
        GeometryNode::PolygonMesh(authored_box(LIFTED)),
        (DeviationPath::AuthoredMesh, WARPED_AUTHORED),
    );
}

#[test]
fn a_lifted_corner_cut_by_a_planar_box_covers_its_slab() {
    assert_covers_the_warp(
        LIFTED,
        BooleanOperator::Difference,
        GeometryNode::PolygonMesh(authored_box(LIFTED)),
        (DeviationPath::AuthoredMesh, WARPED_AUTHORED),
    );
}

#[test]
fn a_saddle_unioned_with_a_planar_box_covers_its_slab() {
    assert_covers_the_warp(
        SADDLE,
        BooleanOperator::Union,
        GeometryNode::PolygonMesh(authored_box(SADDLE)),
        (DeviationPath::AuthoredMesh, WARPED_AUTHORED),
    );
}

#[test]
fn a_saddle_cut_by_a_planar_box_covers_its_slab() {
    assert_covers_the_warp(
        SADDLE,
        BooleanOperator::Difference,
        GeometryNode::PolygonMesh(authored_box(SADDLE)),
        (DeviationPath::AuthoredMesh, WARPED_AUTHORED),
    );
}

/// A warped B-rep face without a surface (#257) composes the same way.
#[test]
fn a_warped_brep_box_cut_and_unioned_covers_its_slab() {
    for operator in [BooleanOperator::Difference, BooleanOperator::Union] {
        for top in [LIFTED, SADDLE] {
            assert_covers_the_warp(
                top,
                operator,
                GeometryNode::BRep(brep_box(top)),
                (DeviationPath::BRepFace, WARPED_BREP),
            );
        }
    }
}

/// A warped box as the tool: a planar box cut by it carries its slab too.
#[test]
fn a_planar_box_cut_by_a_warped_box_covers_its_slab() {
    for top in [LIFTED, SADDLE] {
        let mut g = GeometryGraphBuilder::new();
        let warped = g
            .push(GeometryNode::PolygonMesh(authored_box(top)))
            .unwrap();
        let tool = tool(&mut g);
        let root = boolean(&mut g, tool, warped, BooleanOperator::Difference);
        let graph = g.finish(vec![root]).unwrap();
        let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_with_deviation(&graph, root, &options(1e-3))
            .unwrap();
        let width = slab_width(top);
        let measured = contribution(&report, DeviationPath::Boolean, FLATTENED)
            .unwrap_or_else(|| panic!("{report:?}"));
        assert!(
            matches!(measured.bound, DeviationBound::Certified(v)
                if v >= width && v <= width + 1e-3),
            "{report:?}"
        );
        assert!(contribution(&report, DeviationPath::AuthoredMesh, WARPED_AUTHORED).is_some());
        assert!(!report.meets_requested());
    }
}

/// A warped box clipped by a half-space (a roof cut) is measured the same
/// way: the clip takes the flattened box as its subject.
#[test]
fn a_clipped_warped_box_covers_its_slab() {
    for top in [LIFTED, SADDLE] {
        let mut g = GeometryGraphBuilder::new();
        let warped = g
            .push(GeometryNode::PolygonMesh(authored_box(top)))
            .unwrap();
        let roof = g
            .push(GeometryNode::HalfSpace(axiolid_primitive::HalfSpace {
                boundary: axiolid_core::Plane3 {
                    origin: Point3::new(0.6, 0.0, 0.5),
                    normal: Vec3::new(1.0, 0.0, 0.4),
                },
                agreement: true,
            }))
            .unwrap();
        let root = boolean(&mut g, warped, roof, BooleanOperator::Difference);
        let graph = g.finish(vec![root]).unwrap();
        let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
            .compile_mesh_with_deviation(&graph, root, &options(0.2))
            .unwrap();
        let width = slab_width(top);
        let measured = contribution(&report, DeviationPath::Boolean, FLATTENED)
            .unwrap_or_else(|| panic!("{report:?}"));
        assert!(
            matches!(measured.bound, DeviationBound::Certified(v) if v >= width && v <= 0.2),
            "{report:?}"
        );
        assert!(report.meets_requested(), "{report:?}");
    }
}

/// An unwarped authored box against a planar box stays exact: no slab, no
/// flattening, measured against the exact boolean of the operands as
/// given, within the budget.
#[test]
fn an_unwarped_pair_stays_exact() {
    for operator in [
        BooleanOperator::Union,
        BooleanOperator::Difference,
        BooleanOperator::Intersection,
    ] {
        for node in [
            GeometryNode::PolygonMesh(authored_box(FLAT)),
            GeometryNode::BRep(brep_box(FLAT)),
        ] {
            let (mesh, report) = compile(node, operator, 1e-3);
            assert!(report.meets_requested(), "{operator:?}: {report:?}");
            let measured = contribution(&report, DeviationPath::Boolean, MEASURED_EXACT)
                .unwrap_or_else(|| panic!("{operator:?} exact: {report:?}"));
            assert!(
                matches!(measured.bound, DeviationBound::Certified(v) if v <= 1e-3),
                "{report:?}"
            );
            assert!(
                !report.contributions.iter().any(|c| c.detail == FLATTENED
                    || c.detail == WARPED_AUTHORED
                    || c.detail == WARPED_BREP),
                "{report:?}"
            );
            let worst = boolean_samples(FLAT, Reading::Bilinear, operator, 40)
                .into_iter()
                .map(|p| distance_to_mesh(p, &mesh))
                .fold(0.0, Scalar::max);
            assert!(worst <= 1e-12, "{operator:?}: {worst}");
        }
    }
}

/// Under an instance and through a chain, the slab width travels with the
/// operand: an outer, emitted boolean measures the whole chain flattened.
#[test]
fn a_chain_of_cuts_carries_the_slab_width() {
    let mut g = GeometryGraphBuilder::new();
    let warped = g
        .push(GeometryNode::PolygonMesh(authored_box(SADDLE)))
        .unwrap();
    let placed = g
        .push(GeometryNode::Instance(Instance {
            source: warped,
            transform: Transform3::from_rotation_z(0.0),
        }))
        .unwrap();
    let first = tool(&mut g);
    let cut = boolean(&mut g, placed, first, BooleanOperator::Difference);
    let second = g
        .push(GeometryNode::Instance(Instance {
            source: first,
            transform: Transform3::from_translation(Vec3::new(-1.0, 0.0, 0.0)),
        }))
        .unwrap();
    let root = boolean(&mut g, cut, second, BooleanOperator::Difference);
    let graph = g.finish(vec![root]).unwrap();
    let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options(1e-3))
        .unwrap();
    let width = slab_width(SADDLE);
    assert!(contribution(&report, DeviationPath::AuthoredMesh, WARPED_AUTHORED).is_some());
    let measured = contribution(&report, DeviationPath::Boolean, FLATTENED)
        .unwrap_or_else(|| panic!("{report:?}"));
    assert!(
        matches!(measured.bound, DeviationBound::Certified(v) if v >= width),
        "{report:?}"
    );
    assert!(!report.meets_requested());
}

/// Where the exact boolean cannot take the flattened operand, the report
/// names why, and still names the warp: here the tool is a
/// primitive, which the exact compiler refuses.
#[test]
fn a_refused_flattened_boolean_is_unbounded_by_name() {
    let mut g = GeometryGraphBuilder::new();
    let warped = g
        .push(GeometryNode::PolygonMesh(authored_box(LIFTED)))
        .unwrap();
    let tool = g
        .push(GeometryNode::Primitive(
            axiolid_primitive::Primitive::Sphere { radius: 0.2 },
        ))
        .unwrap();
    let root = boolean(&mut g, warped, tool, BooleanOperator::Difference);
    let graph = g.finish(vec![root]).unwrap();
    let (_, report) = ReferenceMeshCompiler::new(BoolmeshBoolean::new())
        .compile_mesh_with_deviation(&graph, root, &options(1e-3))
        .unwrap();
    assert_eq!(report.bound, None, "{report:?}");
    let refused: Vec<_> = report.unbounded().collect();
    assert_eq!(refused.len(), 1, "{report:?}");
    assert!(
        matches!(refused[0].bound, DeviationBound::Unbounded(reason) if reason.contains("primitive")),
        "{report:?}"
    );
    assert!(contribution(&report, DeviationPath::AuthoredMesh, WARPED_AUTHORED).is_some());
}

/// The public exact compiler still refuses a polygon mesh under a boolean:
/// only the report flattens it.
#[test]
fn the_exact_compiler_still_refuses_a_mesh_operand() {
    use axiolid_exact_compile_contract::ExactCompiler;
    let mut g = GeometryGraphBuilder::new();
    let warped = g
        .push(GeometryNode::PolygonMesh(authored_box(FLAT)))
        .unwrap();
    let tool = tool(&mut g);
    let root = boolean(&mut g, warped, tool, BooleanOperator::Difference);
    let graph = g.finish(vec![root]).unwrap();
    let error = axiolid_mesh_compile::ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options(1e-3))
        .unwrap_err();
    assert!(
        matches!(error, GeomError::UnsupportedInput { input, .. } if input.contains("polygon mesh")),
        "{error:?}"
    );
}
