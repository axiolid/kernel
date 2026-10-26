//! Open surfaces compile through `compile_mesh` as surfaces (#311).
//!
//! A surface trimmed by a parameter rectangle and a B-rep of loose faces
//! compile to a mesh whose closure is `MeshClosure::Surface`: area, no
//! volume, never an operand of a boolean. An atomic `Surface` is unbounded and is still refused. The
//! trimmed patches' reported bound is checked against the exact surface:
//! points sampled densely over the trimmed rectangle (through
//! `axiolid_reference::surface`, the exact evaluation) must lie within it
//! of the mesh, and it within the requested budget.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use axiolid_contracts::{ExecutionOptions, GeomError, GeomResult, Operation};
use axiolid_core::{BooleanOperator, Frame3, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::KnotSpec;
use axiolid_measure::proximity::closest_point_on_triangle;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::{DeviationBound, DeviationPath, DeviationReport, ReferenceMeshCompiler};
use axiolid_mesh_compile_contract::{CompileOutcome, MeshClosure, MeshCompiler};
use axiolid_model::{
    GeometryGraph, GeometryGraphBuilder, GeometryNode, NodeId, SolidOperation, SurfaceRelation,
};
use axiolid_surface::{BSplineSurface, Cylinder, Plane, Sphere, Surface};
use axiolid_topology::{BRep, Edge, EdgeUse, Face, FaceBound, Loop, Orientation, Vertex};

const CHORD: Scalar = 1e-3;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(CHORD)
        .expect("positive chord budget")
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn frame() -> Frame3 {
    let z = Vec3::new(1.0, -2.0, 2.0).normalize();
    let x = Vec3::new(2.0, 1.0, 0.0).normalize();
    Frame3 {
        origin: Point3::new(5.0, -3.0, 1.5),
        x,
        y: z.cross(x),
        z,
    }
}

struct Trimmed {
    graph: GeometryGraph,
    root: NodeId,
}

fn trimmed(
    surface: Surface,
    u: (Scalar, Scalar),
    v: (Scalar, Scalar),
    u_sense: bool,
    v_sense: bool,
) -> Trimmed {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(GeometryNode::Surface(surface)).unwrap();
    let root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::RectangularTrimmed {
                basis,
                u,
                v,
                u_sense,
                v_sense,
            },
        ))
        .unwrap();
    Trimmed {
        graph: b.finish(vec![root]).unwrap(),
        root,
    }
}

fn compile(t: &Trimmed) -> GeomResult<(CompileOutcome, DeviationReport)> {
    compiler().compile_mesh_with_deviation(&t.graph, t.root, &options())
}

fn area(mesh: &TriMesh) -> Scalar {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            0.5 * (b - a).cross(c - a).length()
        })
        .sum()
}

fn distance_to_mesh(mesh: &TriMesh, p: Point3) -> Scalar {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let triangle = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            (closest_point_on_triangle(p, triangle).unwrap() - p).length()
        })
        .fold(Scalar::INFINITY, Scalar::min)
}

/// The report's bound is certified, within the budget, and no sampled
/// exact point lies farther from the mesh; the largest sampled distance.
fn within_bound(
    outcome: &CompileOutcome,
    report: &DeviationReport,
    surface: &Surface,
    (u0, u1): (Scalar, Scalar),
    (v0, v1): (Scalar, Scalar),
    family: &str,
) -> Scalar {
    let bound = report.bound.expect("a certified bound");
    assert!(report.meets_requested(), "{bound} > {CHORD}");
    assert!(report
        .contributions
        .iter()
        .any(|c| c.path == DeviationPath::TrimmedSurface
            && c.detail == family
            && matches!(c.bound, DeviationBound::Certified(_))));
    let n = 48;
    let mut worst: Scalar = 0.0;
    for i in 0..=n {
        for j in 0..=n {
            let u = u0 + (u1 - u0) * (i as Scalar + 0.37) / (n as Scalar + 1.0);
            let v = v0 + (v1 - v0) * (j as Scalar + 0.61) / (n as Scalar + 1.0);
            let p = axiolid_reference::surface::evaluate(surface, u, v).unwrap();
            let d = distance_to_mesh(&outcome.mesh, p);
            assert!(
                d <= bound + 1e-12 * p.abs().max_element().max(1.0),
                "({u}, {v}) at {p:?} is {d} from the mesh, beyond the bound {bound}"
            );
            worst = worst.max(d);
        }
    }
    worst
}

fn cylinder() -> Surface {
    Surface::Cylinder(Cylinder {
        frame: frame(),
        radius: 2.0,
    })
}

/// The mean of the mesh's triangle normals, area weighted.
fn mean_normal(mesh: &TriMesh) -> Vec3 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            (b - a).cross(c - a)
        })
        .fold(Vec3::ZERO, |sum, n| sum + n)
}

// --- trimmed surfaces ----------------------------------------------------------

#[test]
fn a_trimmed_cylinder_is_an_open_sheet_within_its_certified_bound() {
    let (u, v) = ((0.3, 0.3 + FRAC_PI_2), (-1.0, 2.0));
    let t = trimmed(cylinder(), u, v, true, true);
    let (outcome, report) = compile(&t).unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!(matches!(
        outcome.solid_mesh(),
        Err(GeomError::InvalidInput(message)) if message.contains("surface model")
    ));
    let worst = within_bound(&outcome, &report, &cylinder(), u, v, "cylinder");
    // A cylinder's chords are its only error, and the bound is the
    // sagitta's: the densest sample comes near it.
    assert!(
        worst >= 0.5 * report.bound.unwrap(),
        "{worst} vs {report:?}"
    );
    // Area: a quarter turn of radius 2, 3 long, less the chords' share.
    let exact = FRAC_PI_2 * 2.0 * 3.0;
    let meshed = area(&outcome.mesh);
    assert!(meshed < exact && exact - meshed < 1e-3 * exact, "{meshed}");
    // Facing the cylinder's outward normal at the patch's middle.
    let mid = axiolid_reference::surface::normal(&cylinder(), 0.3 + FRAC_PI_2 / 2.0, 0.5).unwrap();
    assert!(mean_normal(&outcome.mesh).dot(mid) > 0.0);
}

#[test]
fn a_sense_against_the_basis_turns_the_sheet_and_wraps_the_seam() {
    // Against u: from 0.3 down to -1.27, the same quarter turn, facing in.
    let along = trimmed(cylinder(), (0.3, 0.3 + FRAC_PI_2), (-1.0, 2.0), true, true);
    let against = trimmed(cylinder(), (0.3 + FRAC_PI_2, 0.3), (-1.0, 2.0), false, true);
    let (a, _) = compile(&along).unwrap();
    let (b, _) = compile(&against).unwrap();
    assert!((area(&a.mesh) - area(&b.mesh)).abs() < 1e-9);
    assert!(mean_normal(&a.mesh).dot(mean_normal(&b.mesh)) < 0.0);
    // Both against: the normal turns twice, back to the basis's.
    let both = trimmed(
        cylinder(),
        (0.3 + FRAC_PI_2, 0.3),
        (2.0, -1.0),
        false,
        false,
    );
    let (c, _) = compile(&both).unwrap();
    assert!(mean_normal(&a.mesh).dot(mean_normal(&c.mesh)) > 0.0);

    // With u from 5.5 up to 0.5 the trim runs across the seam:
    // 2 pi - 5 of the turn.
    let wrapped = trimmed(cylinder(), (5.5, 0.5), (0.0, 1.0), true, true);
    let (outcome, report) = compile(&wrapped).unwrap();
    within_bound(
        &outcome,
        &report,
        &cylinder(),
        (5.5, 0.5 + TAU),
        (0.0, 1.0),
        "cylinder",
    );
    let exact = (TAU - 5.0) * 2.0;
    assert!((area(&outcome.mesh) - exact).abs() < 1e-3 * exact);
}

#[test]
fn a_trimmed_bspline_and_a_sphere_over_its_pole_are_within_their_bounds() {
    let net: Vec<Vec<Point3>> = (0..4)
        .map(|i| {
            (0..5)
                .map(|j| {
                    let (x, y) = (i as Scalar, j as Scalar);
                    Point3::new(x, y, 0.4 * ((x * 1.3).sin() + (y * 0.9).cos()))
                })
                .collect()
        })
        .collect();
    let spline = Surface::BSpline(BSplineSurface {
        u_degree: 3,
        v_degree: 2,
        control_points: net,
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![4, 4],
        v_knots: vec![0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0],
        v_multiplicities: vec![3, 1, 1, 3],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    let (u, v) = ((0.1, 0.9), (0.2, 1.0));
    let t = trimmed(spline.clone(), u, v, true, true);
    let (outcome, report) = compile(&t).unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!(within_bound(&outcome, &report, &spline, u, v, "B-spline surface") > 0.0);

    // A knot of full multiplicity is a crease no second derivative
    // crosses: the cells over it are unbounded, by name, not guessed.
    let Surface::BSpline(mut creased) = spline else {
        unreachable!()
    };
    creased.v_knots = vec![0.0, 0.5, 1.0];
    creased.v_multiplicities = vec![3, 2, 3];
    let t = trimmed(Surface::BSpline(creased), u, v, true, true);
    let (outcome, report) = compile(&t).unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert_eq!(report.bound, None);
    assert!(report.unbounded().any(|c| matches!(
        c.bound,
        DeviationBound::Unbounded(reason) if reason.contains("knot line")
    )));

    // A cap of a sphere from 30 degrees below its pole over the top: the
    // pole row collapses, and no triangle there has two equal corners.
    let sphere = Surface::Sphere(Sphere {
        frame: frame(),
        radius: 1.5,
    });
    let (u, v) = ((0.0, PI), (FRAC_PI_2 - PI / 6.0, FRAC_PI_2));
    let t = trimmed(sphere.clone(), u, v, true, true);
    let (outcome, report) = compile(&t).unwrap();
    for tri in outcome.mesh.indices.as_chunks::<3>().0 {
        let [a, b, c] = [0, 1, 2].map(|k| outcome.mesh.positions[tri[k] as usize]);
        assert!(
            a != b && b != c && c != a,
            "degenerate triangle at the pole"
        );
    }
    within_bound(&outcome, &report, &sphere, u, v, "sphere");
}

#[test]
fn a_trim_that_cannot_be_bounded_or_read_is_refused_by_name() {
    let plane = Surface::Plane(Plane { frame: frame() });
    // Against its sense on a parameter that does not wrap.
    let t = trimmed(plane.clone(), (0.0, 1.0), (0.0, 1.0), false, true);
    assert!(matches!(
        compile(&t),
        Err(GeomError::InvalidInput(m)) if m.contains("not periodic")
    ));
    // More than a turn.
    let t = trimmed(cylinder(), (0.0, 7.0), (0.0, 1.0), true, true);
    assert!(matches!(
        compile(&t),
        Err(GeomError::InvalidInput(m)) if m.contains("longer than one period")
    ));
    // An empty range.
    let t = trimmed(plane, (0.0, 1.0), (2.0, 2.0), true, true);
    assert!(matches!(
        compile(&t),
        Err(GeomError::InvalidInput(m)) if m.contains("empty")
    ));
}

#[test]
fn an_unbounded_surface_is_still_refused() {
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(GeometryNode::Surface(cylinder())).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    assert!(matches!(
        compiler().compile_mesh(&graph, root, &options()),
        Err(GeomError::Unsupported {
            operation: Operation::SurfaceEvaluation,
            ..
        })
    ));
}

// --- open input stays out of volumes and booleans ------------------------------

#[test]
fn a_boolean_refuses_an_open_surface_and_a_collection_stays_open() {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(GeometryNode::Surface(cylinder())).unwrap();
    let patch = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::RectangularTrimmed {
                basis,
                u: (0.0, 1.0),
                v: (0.0, 1.0),
                u_sense: true,
                v_sense: true,
            },
        ))
        .unwrap();
    let solid = b
        .push(GeometryNode::Primitive(
            axiolid_primitive::Primitive::Sphere { radius: 1.0 },
        ))
        .unwrap();
    // The graph itself refuses an open surface as a boolean's operand.
    assert!(b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: solid,
            right: patch,
            operator: BooleanOperator::Difference,
        }))
        .is_err());
    let both = b
        .push(GeometryNode::Collection(vec![solid, patch]))
        .unwrap();
    let graph = b.finish(vec![both]).unwrap();
    let outcome = compiler()
        .compile_mesh_reported(&graph, both, &options())
        .unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!(outcome.solid_mesh().is_err());
}

// --- faces ---------------------------------------------------------------------

#[test]
fn loose_faces_compile_as_an_open_surface() {
    // Two unit squares that share an edge, in no shell: area 2.
    let mut brep: BRep<NodeId> = BRep::default();
    let corners: Vec<_> = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
    ]
    .iter()
    .map(|p| {
        brep.add_vertex(Vertex {
            position: Vec3::new(p[0], p[1], p[2]),
        })
    })
    .collect();
    let mut edge = |a: usize, b: usize| {
        brep.add_edge(Edge {
            start: corners[a],
            end: corners[b],
            curve: None,
        })
    };
    let edges = [
        edge(0, 1),
        edge(1, 2),
        edge(2, 3),
        edge(3, 0),
        edge(1, 4),
        edge(4, 5),
        edge(5, 2),
    ];
    let forward = |e| EdgeUse {
        edge: e,
        orientation: Orientation::Forward,
        pcurve: None,
    };
    let reversed = |e| EdgeUse {
        edge: e,
        orientation: Orientation::Reversed,
        pcurve: None,
    };
    for ring in [
        vec![
            forward(edges[0]),
            forward(edges[1]),
            forward(edges[2]),
            forward(edges[3]),
        ],
        vec![
            forward(edges[4]),
            forward(edges[5]),
            forward(edges[6]),
            reversed(edges[1]),
        ],
    ] {
        let wire = brep.add_loop(Loop { edges: ring });
        brep.add_face(Face {
            surface: None,
            bounds: vec![FaceBound {
                loop_id: wire,
                outer: true,
                orientation: Orientation::Forward,
            }],
            orientation: Orientation::Forward,
        });
    }
    let mut b = GeometryGraphBuilder::new();
    let root = b.push(GeometryNode::BRep(brep)).unwrap();
    let graph = b.finish(vec![root]).unwrap();
    let outcome = compiler()
        .compile_mesh_reported(&graph, root, &options())
        .unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!((area(&outcome.mesh) - 2.0).abs() < 1e-12);
    // The shared edge is welded: four triangles over six corners.
    assert_eq!(outcome.mesh.positions.len(), 6);
    assert_eq!(outcome.mesh.indices.len(), 12);
}
