use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::{Tolerance, Vec3};
use axiolid_exact_compile_contract::ExactCompiler;
use axiolid_mesh_compile::ReferenceExactCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, SolidOperation};
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};
use axiolid_surface::Surface;
use axiolid_topology::audit_brep;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn extrusion_graph(
    profile: Profile,
    direction: Vec3,
    depth: f64,
) -> (axiolid_model::GeometryGraph, axiolid_model::NodeId) {
    let mut builder = GeometryGraphBuilder::new();
    let profile = builder
        .push(GeometryNode::Profile(profile))
        .expect("profile node");
    let extrusion = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction,
            depth,
        }))
        .expect("extrusion node");
    let graph = builder.finish(vec![extrusion]).expect("valid graph");
    (graph, extrusion)
}

#[test]
fn exact_rectangle_extrusion_compiles_to_closed_brep() {
    let (graph, root) = extrusion_graph(
        Profile::Rectangle(RectangleProfile {
            x: 4.0,
            y: 2.0,
            thickness: Some(0.25),
            outer_radius: None,
            inner_radius: None,
        }),
        Vec3::new(0.25, 0.5, 1.0),
        3.0,
    );

    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options())
        .expect("exact rectangle result");

    assert!(audit_brep(exact.topology()).is_closed_manifold());
    assert_eq!(exact.topology().faces().len(), 10);
    assert!(exact
        .surfaces()
        .iter()
        .all(|surface| matches!(surface, Surface::Plane(_))));
}

#[test]
fn exact_circle_extrusion_leaves_compiler_as_a_cylinder() {
    let (graph, root) = extrusion_graph(
        Profile::Circle(CircleProfile {
            radius: 2.0,
            thickness: None,
        }),
        Vec3::Z,
        5.0,
    );

    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options())
        .expect("exact cylinder result");

    assert!(audit_brep(exact.topology()).is_closed_manifold());
    assert_eq!(
        exact
            .surfaces()
            .iter()
            .filter(|surface| matches!(surface, Surface::Cylinder(_)))
            .count(),
        1
    );
}

#[test]
fn exact_batch_preserves_order_and_duplicate_results() {
    let (graph, root) = extrusion_graph(
        Profile::Rectangle(RectangleProfile {
            x: 2.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Vec3::Z,
        1.0,
    );

    let exact = ReferenceExactCompiler::new()
        .compile_exact_batch(&graph, &[root, root], &options())
        .expect("shared exact batch");

    assert_eq!(exact.len(), 2);
    assert_eq!(exact[0], exact[1]);
}

#[test]
fn constructor_refusal_is_owned_by_graph_compilation() {
    let (graph, root) = extrusion_graph(
        Profile::Circle(CircleProfile {
            radius: 1.0,
            thickness: Some(0.2),
        }),
        Vec3::Z,
        1.0,
    );

    let error = ReferenceExactCompiler::new()
        .compile_exact(&graph, root, &options())
        .expect_err("annular exact extrusion is not implemented");

    assert!(matches!(
        error,
        GeomError::UnsupportedInput {
            backend,
            operation: Operation::GraphCompilation,
            input: "annular circle extrusion",
        } if backend == ReferenceExactCompiler::ID
    ));
}

#[test]
fn unsupported_exact_operation_names_its_family() {
    let mut builder = GeometryGraphBuilder::new();
    let profile = builder
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 2.0,
                y: 1.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let revolution = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Revolution {
            profile,
            axis_origin: Vec3::ZERO,
            axis_direction: Vec3::Z,
            angle: std::f64::consts::TAU,
        }))
        .unwrap();
    let graph = builder.finish(vec![revolution]).unwrap();

    // Exact revolution now exists, but this fixture revolves about the
    // profile's own local z with the axis running through it: the profile
    // crosses the axis, so the swept solid is not an annular tube. The
    // refusal must therefore name that specific geometry rather than claim
    // revolution is unimplemented.
    let error = ReferenceExactCompiler::new()
        .compile_exact(&graph, revolution, &options())
        .expect_err("this profile crosses its revolution axis");

    let GeomError::UnsupportedInput { input, .. } = error else {
        panic!("expected a typed input-family refusal, got {error:?}");
    };
    assert!(
        input.contains("axis"),
        "the refusal must name the axis geometry that blocks it, got {input:?}"
    );
}

#[test]
fn an_offset_revolution_compiles_to_an_exact_annular_tube() {
    // Same graph shape an IFC revolved-area solid produces, with the axis
    // clear of the profile so the result is an annulus.
    let mut builder = GeometryGraphBuilder::new();
    let profile = builder
        .push(GeometryNode::Profile(Profile::Rectangle(
            RectangleProfile {
                x: 2.0,
                y: 3.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            },
        )))
        .unwrap();
    let revolution = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Revolution {
            profile,
            axis_origin: axiolid_core::Point3::new(5.0, 0.0, 0.0),
            axis_direction: Vec3::Y,
            angle: std::f64::consts::TAU,
        }))
        .unwrap();
    let graph = builder.finish(vec![revolution]).unwrap();

    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, revolution, &options())
        .expect("an offset full-turn revolution is exactly constructible");

    // Analytic surfaces, not tessellation: two cylinders and two planes.
    assert_eq!(exact.topology().faces().len(), 4);
    let cylinders = exact
        .surfaces()
        .iter()
        .filter(|s| matches!(s, axiolid_surface::Surface::Cylinder(_)))
        .count();
    assert_eq!(cylinders, 2, "an annular tube has two cylindrical walls");
}

/// A general exact boolean is now reachable through the compiler (#66).
///
/// Previously every exact boolean over two extrusions was a typed refusal.
/// The capability existing in axiolid-construct is not enough -- it has to be
/// reachable through the compiler, which is the flaw #37 was opened for.
#[test]
fn an_exact_boolean_over_two_extrusions_compiles() {
    let mut builder = GeometryGraphBuilder::new();
    let make = |b: &mut GeometryGraphBuilder, x: f64, y: f64| {
        let profile = b
            .push(GeometryNode::Profile(Profile::Rectangle(
                RectangleProfile {
                    x,
                    y,
                    thickness: None,
                    outer_radius: None,
                    inner_radius: None,
                },
            )))
            .unwrap();
        b.push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 3.0,
        }))
        .unwrap()
    };
    let subject = make(&mut builder, 4.0, 4.0);
    let tool = make(&mut builder, 2.0, 2.0);
    let cut = builder
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: subject,
            right: tool,
            operator: axiolid_core::BooleanOperator::Difference,
        }))
        .unwrap();
    let graph = builder.finish(vec![cut]).unwrap();

    let exact = ReferenceExactCompiler::new()
        .compile_exact(&graph, cut, &options())
        .expect("a full-height coaxial difference is exactly constructible");

    // The opening is interior, so the cap gains a hole loop.
    let max_bounds = exact
        .topology()
        .faces()
        .iter()
        .map(|f| f.bounds.len())
        .max()
        .expect("the solid has faces");
    assert!(
        max_bounds >= 2,
        "the cut must appear as a hole loop, got {max_bounds}"
    );
}

fn rect_ring(cx: f64, cy: f64, w: f64, h: f64) -> Vec<axiolid_core::Point2> {
    let (hw, hh) = (w / 2.0, h / 2.0);
    vec![
        axiolid_core::Point2::new(cx - hw, cy - hh),
        axiolid_core::Point2::new(cx + hw, cy - hh),
        axiolid_core::Point2::new(cx + hw, cy + hh),
        axiolid_core::Point2::new(cx - hw, cy + hh),
    ]
}

fn simple_base_area(brep: &axiolid_brep::ExactBRep) -> f64 {
    let mut points: Vec<(f64, f64)> = brep
        .topology()
        .vertices()
        .iter()
        .filter(|v| v.position.z.abs() < 1e-9)
        .map(|v| (v.position.x, v.position.y))
        .collect();
    points.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9);
    (0..points.len())
        .map(|i| {
            let (x0, y0) = points[i];
            let (x1, y1) = points[(i + 1) % points.len()];
            x0 * y1 - x1 * y0
        })
        .sum::<f64>()
        .abs()
        / 2.0
}

/// Differential: the exact prism boolean agrees with the mesh oracle (#66).
///
/// This lives in the compile crate rather than beside the implementation
/// because axiolid-construct sits in the ALGORITHMS layer and boolmesh is a
/// PROVIDER: an algorithms crate depending on a provider inverts the
/// layering, and the architecture gate rejects it even as a dev-dependency.
/// The compile crate already depends on both, so the differential belongs
/// here.
///
/// The two paths share no code -- boolmesh works on triangles, the exact
/// path on the planar overlay -- so agreement is evidence, not tautology.
#[test]
fn the_exact_prism_boolean_agrees_with_the_mesh_oracle() {
    use axiolid_construct::boolean_exact::{boolean_prisms_exact, Prism};
    use axiolid_construct::extrude::extrude_profile;
    use axiolid_construct::profile::Rings;
    use axiolid_contracts::ExecutionOptions;
    use axiolid_core::Point2;
    use axiolid_core::Vec3;
    use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
    use axiolid_mesh_boolean_contract::MeshBoolean;

    let subject_ring = rect_ring(0.0, 0.0, 4.0, 4.0);
    let tool_ring = rect_ring(2.0, 0.0, 4.0, 4.0);
    let height = 3.0;

    // Exact path.
    let exact = boolean_prisms_exact(
        &Prism {
            rings: vec![subject_ring.clone()],
            bottom: 0.0,
            top: height,
        },
        &Prism {
            rings: vec![tool_ring.clone()],
            bottom: 0.0,
            top: height,
        },
        axiolid_core::BooleanOperator::Intersection,
        axiolid_core::Tolerance::METRE,
    )
    .expect("exact intersection");
    let exact_volume = simple_base_area(&exact) * height;

    // Mesh path: build both prisms as meshes and intersect with boolmesh.
    let to_mesh = |ring: &[Point2]| {
        extrude_profile(
            &Rings {
                outer: ring.to_vec(),
                holes: Vec::new(),
            },
            Vec3::Z,
            height,
            axiolid_core::Tolerance::METRE,
        )
        .expect("a rectangle extrudes")
    };
    let mesh_result = BoolmeshBoolean::new()
        .boolean(
            &to_mesh(&subject_ring),
            &to_mesh(&tool_ring),
            axiolid_core::BooleanOperator::Intersection,
            &ExecutionOptions::new(Tolerance::METRE),
        )
        .expect("the mesh oracle intersects")
        .mesh;

    // Divergence-theorem volume: triangulation-invariant.
    let mesh_volume: f64 = mesh_result
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let a = mesh_result.positions[t[0] as usize];
            let b = mesh_result.positions[t[1] as usize];
            let c = mesh_result.positions[t[2] as usize];
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0;

    assert!(
        (exact_volume - mesh_volume.abs()).abs() < 1e-6,
        "exact and mesh paths disagree: {exact_volume} vs {}",
        mesh_volume.abs()
    );
    // And both must equal the hand-computed 2x4x3.
    assert!(
        (exact_volume - 24.0).abs() < 1e-9,
        "expected 24, got {exact_volume}"
    );
}
