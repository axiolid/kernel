//! Curve-bounded planes (#192): IFC `IfcCurveBoundedPlane` space-boundary
//! connection surfaces. Oracles are closed-form areas.

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::{Frame3, Point2, Point3, Tolerance, Transform3, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Polyline2, Polyline3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::{CompileOutcome, MeshClosure, MeshCompiler};
use axiolid_model::{GeometryGraphBuilder, GeometryNode, NodeId, SurfaceRelation};
use axiolid_surface::{Cylinder, Plane, Surface};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::MILLIMETRE)
        .with_chord_error(1e-4)
        .expect("positive chord budget")
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

fn area(mesh: &TriMesh) -> f64 {
    mesh.indices
        .chunks_exact(3)
        .map(|t| {
            let [a, b, c] = [0, 1, 2].map(|k| mesh.positions[t[k] as usize]);
            0.5 * (b - a).cross(c - a).length()
        })
        .sum()
}

fn polyline(points: &[(f64, f64)]) -> GeometryNode {
    GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
        closed: true,
    }))
}

fn frame() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

/// Build a curve-bounded plane over `frame` from boundary nodes, and
/// compile it (optionally placed by an instance transform).
fn compile(
    frame: Frame3,
    boundaries: Vec<GeometryNode>,
    implicit_outer: bool,
    place: Option<Transform3>,
) -> Result<CompileOutcome, GeomError> {
    compile_over(
        GeometryNode::Surface(Surface::Plane(Plane { frame })),
        boundaries,
        implicit_outer,
        place,
    )
}

fn compile_over(
    basis: GeometryNode,
    boundaries: Vec<GeometryNode>,
    implicit_outer: bool,
    place: Option<Transform3>,
) -> Result<CompileOutcome, GeomError> {
    let mut b = GeometryGraphBuilder::new();
    let basis = b.push(basis).unwrap();
    let boundaries: Vec<NodeId> = boundaries.into_iter().map(|n| b.push(n).unwrap()).collect();
    let mut root = b
        .push(GeometryNode::SurfaceRelation(
            SurfaceRelation::CurveBounded {
                basis,
                boundaries,
                implicit_outer,
            },
        ))
        .unwrap();
    if let Some(transform) = place {
        root = b
            .push(GeometryNode::Instance(axiolid_model::Instance {
                source: root,
                transform,
            }))
            .unwrap();
    }
    let graph = b.finish(vec![root]).unwrap();
    compiler().compile_mesh_reported(&graph, root, &options())
}

fn plate_with_hole() -> Vec<GeometryNode> {
    vec![
        polyline(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
        polyline(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)]),
    ]
}

#[test]
fn a_four_by_three_plane_with_a_hole_has_area_eleven() {
    let outcome = compile(frame(), plate_with_hole(), false, None).unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert_eq!(area(&outcome.mesh), 11.0);
    // Every triangle faces along the plane's normal.
    for t in outcome.mesh.indices.chunks_exact(3) {
        let [a, b, c] = [0, 1, 2].map(|k| outcome.mesh.positions[t[k] as usize]);
        assert!((b - a).cross(c - a).z > 0.0);
    }
}

#[test]
fn a_placed_plane_keeps_its_area() {
    let place = Transform3::from_translation(Vec3::new(10.0, -4.0, 2.5))
        * Transform3::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 0.7);
    let outcome = compile(frame(), plate_with_hole(), false, Some(place)).unwrap();
    assert_eq!(outcome.closure, MeshClosure::Surface);
    assert!(
        (area(&outcome.mesh) - 11.0).abs() < 1e-12,
        "{}",
        area(&outcome.mesh)
    );
}

#[test]
fn boundaries_are_in_the_planes_parameters() {
    // A plane standing up: origin (5, 0, 0), x along world y, y along
    // world z. The parameter corner (4, 3) is the point (5, 4, 3).
    let standing = Frame3 {
        origin: Point3::new(5.0, 0.0, 0.0),
        x: Vec3::Y,
        y: Vec3::Z,
        z: Vec3::X,
    };
    let outcome = compile(standing, plate_with_hole(), false, None).unwrap();
    assert!(outcome.mesh.positions.contains(&Point3::new(5.0, 4.0, 3.0)));
    assert!(outcome.mesh.positions.iter().all(|p| p.x == 5.0));
    assert_eq!(area(&outcome.mesh), 11.0);
    // 3D boundaries on the parameter plane (z = 0) read the same.
    let three_d = |pts: &[(f64, f64)]| {
        GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
            points: pts.iter().map(|&(x, y)| Point3::new(x, y, 0.0)).collect(),
            closed: true,
        }))
    };
    let outcome = compile(
        standing,
        vec![
            three_d(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            three_d(&[(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)]),
        ],
        false,
        None,
    )
    .unwrap();
    assert_eq!(area(&outcome.mesh), 11.0);
}

#[test]
fn a_curved_hole_is_chorded_within_the_budget() {
    let r = 0.5;
    let circle = GeometryNode::Curve2(Curve2::Circle(Circle2 {
        frame: axiolid_core::Frame2 {
            origin: Point2::new(2.0, 1.5),
            x: axiolid_core::Vec2::X,
            y: axiolid_core::Vec2::Y,
        },
        radius: r,
    }));
    let outcome = compile(
        frame(),
        vec![
            polyline(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            circle,
        ],
        false,
        None,
    )
    .unwrap();
    let exact = 12.0 - std::f64::consts::PI * r * r;
    let got = area(&outcome.mesh);
    // The chords cut inside the circle, so the hole is smaller and the
    // plate larger, by at most the perimeter times the chord budget.
    assert!(
        got >= exact && got - exact <= std::f64::consts::TAU * r * 1e-4,
        "{got} vs {exact}"
    );
}

#[test]
fn invalid_loops_are_refused_never_filled() {
    let bowtie = vec![polyline(&[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0)])];
    let outside = vec![
        polyline(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
        polyline(&[(5.0, 1.0), (6.0, 1.0), (6.0, 2.0), (5.0, 2.0)]),
    ];
    let degenerate = vec![polyline(&[(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)])];
    let two_points = vec![polyline(&[(0.0, 0.0), (1.0, 1.0)])];
    let open = vec![GeometryNode::Curve2(Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(4.0, 3.0),
        ],
        closed: false,
    }))];
    let off_plane = vec![GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.5),
            Point3::new(4.0, 3.0, 0.0),
        ],
        closed: true,
    }))];
    for (name, boundaries) in [
        ("bow tie", bowtie),
        ("hole outside", outside),
        ("collinear", degenerate),
        ("two points", two_points),
        ("open", open),
        ("off the plane", off_plane),
    ] {
        let result = compile(frame(), boundaries, false, None);
        assert!(
            matches!(
                result,
                Err(GeomError::InvalidInput(_) | GeomError::Degenerate(_))
            ),
            "{name}: {result:?}"
        );
    }
}

#[test]
fn unsupported_forms_name_the_missing_capability() {
    let cylinder = GeometryNode::Surface(Surface::Cylinder(Cylinder {
        frame: frame(),
        radius: 1.0,
    }));
    match compile_over(cylinder, plate_with_hole(), false, None) {
        Err(GeomError::UnsupportedInput {
            operation, input, ..
        }) => {
            assert_eq!(operation, Operation::SurfaceEvaluation);
            assert!(input.contains("non-planar"), "{input}");
        }
        other => panic!("{other:?}"),
    }
    match compile(frame(), plate_with_hole(), true, None) {
        Err(GeomError::UnsupportedInput {
            operation, input, ..
        }) => {
            assert_eq!(operation, Operation::Tessellation);
            assert!(input.contains("implicit outer"), "{input}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_clockwise_outer_loop_still_faces_along_the_normal() {
    // Exporters write loops either way round; the surface's side is the
    // plane's normal, whatever the loop's winding.
    let outcome = compile(
        frame(),
        vec![
            polyline(&[(0.0, 0.0), (0.0, 3.0), (4.0, 3.0), (4.0, 0.0)]),
            polyline(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]),
        ],
        false,
        None,
    )
    .unwrap();
    assert_eq!(area(&outcome.mesh), 11.0);
    for t in outcome.mesh.indices.chunks_exact(3) {
        let [a, b, c] = [0, 1, 2].map(|k| outcome.mesh.positions[t[k] as usize]);
        assert!((b - a).cross(c - a).z > 0.0);
    }
}
