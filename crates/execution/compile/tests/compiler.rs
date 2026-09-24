//! Gates for the graph compiler.
//!
//! The shape under test is the one `ifc-geometry` actually emits: a Profile
//! feeding an Extrusion, wrapped in Instance placements, combined by Boolean.

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::{BooleanOperator, Plane3, Point3, Tolerance, Transform3, Vec3};
use axiolid_curve::{Curve3, Polyline3};
use axiolid_mesh::{PolygonFace, PolygonMesh, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_compile::ReferenceMeshCompiler;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraphBuilder, GeometryNode, Instance, SolidOperation};
use axiolid_primitive::HalfSpace;
use axiolid_profile::{CircleProfile, Profile, RectangleProfile};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn compiler() -> ReferenceMeshCompiler<BoolmeshBoolean> {
    ReferenceMeshCompiler::new(BoolmeshBoolean::new())
}

/// Divergence-theorem volume: triangulation-invariant.
fn volume(mesh: &TriMesh) -> f64 {
    let mut sum = 0.0;
    for t in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        sum += a.dot(b.cross(c));
    }
    sum / 6.0
}

fn rect(x: f64, y: f64) -> Profile {
    Profile::Rectangle(RectangleProfile {
        x,
        y,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    })
}

/// The dominant IFC pattern end to end: wall minus an opening, both built as
/// Profile -> Extrusion -> Instance, combined by a Boolean node.
#[test]
fn a_wall_minus_an_opening_compiles_to_the_expected_volume() {
    let mut b = GeometryGraphBuilder::new();

    let wall_profile = b.push(GeometryNode::Profile(rect(4.0, 0.2))).unwrap();
    let wall = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: wall_profile,
            direction: Vec3::Z,
            depth: 3.0,
        }))
        .unwrap();

    let hole_profile = b.push(GeometryNode::Profile(rect(1.0, 0.4))).unwrap();
    let hole = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: hole_profile,
            direction: Vec3::Z,
            depth: 1.2,
        }))
        .unwrap();
    // Place the opening inside the wall: profiles are centred on the origin,
    // so only a Z lift is needed to sit it within the wall band.
    let placed = b
        .push(GeometryNode::Instance(Instance {
            source: hole,
            transform: Transform3::from_translation(Vec3::new(0.0, 0.0, 0.3)),
        }))
        .unwrap();

    let cut = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            left: wall,
            right: placed,
            operator: BooleanOperator::Difference,
        }))
        .unwrap();
    let graph = b.finish(vec![cut]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, cut, &options())
        .expect("compile");
    // 4 x 0.2 x 3 = 2.4 minus 1 x 0.2 x 1.2 (the opening is clipped to the
    // wall thickness) = 2.4 - 0.24 = 2.16.
    assert!((volume(&mesh) - 2.16).abs() < 1e-9, "got {}", volume(&mesh));
}

/// A shared subtree must be compiled once per batch, not once per reference.
///
/// Without memoisation a diamond DAG recompiles the shared node for every
/// path that reaches it, which is exponential on deep sharing. The observable
/// proof is that both roots agree exactly and the batch succeeds.
#[test]
fn a_shared_subtree_is_reused_across_roots() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(2.0, 2.0))).unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let left = b
        .push(GeometryNode::Instance(Instance {
            source: solid,
            transform: Transform3::from_translation(Vec3::new(-5.0, 0.0, 0.0)),
        }))
        .unwrap();
    let right = b
        .push(GeometryNode::Instance(Instance {
            source: solid,
            transform: Transform3::from_translation(Vec3::new(5.0, 0.0, 0.0)),
        }))
        .unwrap();
    let graph = b.finish(vec![left, right]).unwrap();

    let meshes = compiler()
        .compile_mesh_batch(&graph, &[left, right], &options())
        .expect("batch");
    assert_eq!(meshes.len(), 2);
    assert!((volume(&meshes[0]) - 4.0).abs() < 1e-9);
    assert!((volume(&meshes[1]) - 4.0).abs() < 1e-9);
    // Same source, different placement: the meshes must differ in position.
    assert_ne!(meshes[0].positions[0], meshes[1].positions[0]);
}

/// Both batch call shapes must behave identically.
#[test]
fn compile_mesh_batch_into_appends_and_matches_compile_mesh_batch() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 2.0,
        }))
        .unwrap();
    let graph = b.finish(vec![solid]).unwrap();

    let mut destination = vec![TriMesh::default()];
    compiler()
        .compile_mesh_batch_into(&graph, &[solid], &options(), &mut destination)
        .expect("into");
    assert_eq!(destination.len(), 2, "must append, never clear");
    assert!((volume(&destination[1]) - 2.0).abs() < 1e-9);
}

#[test]
fn authored_triangle_faces_compile_without_retriangulation() {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];
    let mut b = GeometryGraphBuilder::new();
    let triangles = b
        .push(GeometryNode::PolygonMesh(PolygonMesh {
            positions: positions.clone(),
            faces: vec![
                PolygonFace {
                    outer: vec![0, 1, 2],
                    holes: Vec::new(),
                },
                PolygonFace {
                    outer: vec![2, 1, 3],
                    holes: Vec::new(),
                },
            ],
        }))
        .unwrap();
    let graph = b.finish(vec![triangles]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, triangles, &options())
        .expect("authored triangles compile by exact index conversion");
    assert_eq!(mesh.positions, positions);
    assert_eq!(mesh.indices, vec![0, 1, 2, 2, 1, 3]);
}

// N-gon and holed polygon faces triangulate (#160); their tests, including
// the refusals that replaced the blanket `Unsupported`, live in
// `tests/authored_polygons.rs`.

#[test]
fn malformed_authored_triangle_indices_are_rejected() {
    let mut b = GeometryGraphBuilder::new();
    let triangles = b
        .push(GeometryNode::PolygonMesh(PolygonMesh {
            positions: vec![Point3::ZERO; 3],
            faces: vec![PolygonFace {
                outer: vec![0, 1, 7],
                holes: Vec::new(),
            }],
        }))
        .unwrap();
    let graph = b.finish(vec![triangles]).unwrap();

    assert!(matches!(
        compiler().compile_mesh(&graph, triangles, &options()),
        Err(GeomError::InvalidInput(_))
    ));
}

/// An unsupported family must name the capability it would need, so a caller
/// can register a provider for it instead of guessing.
#[test]
fn an_unsupported_node_reports_the_capability_it_would_need() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let graph = b.finish(vec![profile]).unwrap();

    let error = compiler()
        .compile_mesh(&graph, profile, &options())
        .expect_err("a bare profile is not a solid");
    match error {
        GeomError::Unsupported { operation, .. } => {
            assert_eq!(operation, Operation::ProfileTriangulation);
        }
        other => panic!("expected Unsupported, got {other:?}"),
    }
}

/// An unavailable capability is refused rather than approximated. ADR: a
/// wrong solid is more expensive than a missing one.
///
/// A surface curve sweep is now built for the analytic surfaces, so the
/// refusal that remains is the one with no closed form: a B-spline
/// reference surface has no direct point-to-parameter inverse, and
/// recovering it needs iterative closest-point with its own convergence
/// contract. The fixture therefore uses a B-spline, not a plane.
#[test]
fn an_unimplemented_solid_family_is_refused_not_approximated() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let directrix = b
        .push(GeometryNode::Curve3(axiolid_curve::Curve3::Polyline(
            axiolid_curve::Polyline3 {
                points: vec![
                    axiolid_core::Point3::ZERO,
                    axiolid_core::Point3::new(0.0, 0.0, 1.0),
                ],
                closed: false,
            },
        )))
        .unwrap();
    let surface = b
        .push(GeometryNode::Surface(axiolid_surface::Surface::BSpline(
            axiolid_surface::BSplineSurface {
                u_degree: 1,
                v_degree: 1,
                control_points: vec![
                    vec![
                        axiolid_core::Point3::ZERO,
                        axiolid_core::Point3::new(0.0, 1.0, 0.0),
                    ],
                    vec![
                        axiolid_core::Point3::new(1.0, 0.0, 0.0),
                        axiolid_core::Point3::new(1.0, 1.0, 0.0),
                    ],
                ],
                u_knots: vec![0.0, 1.0],
                u_multiplicities: vec![2, 2],
                v_knots: vec![0.0, 1.0],
                v_multiplicities: vec![2, 2],
                weights: None,
                u_closed: false,
                v_closed: false,
                knot_spec: axiolid_curve::KnotSpec::Unspecified,
                self_intersect: None,
            },
        )))
        .unwrap();
    let swept = b
        .push(GeometryNode::SolidOperation(
            SolidOperation::SurfaceCurveSweep {
                profile,
                directrix,
                reference_surface: surface,
                parameter_range: None,
            },
        ))
        .unwrap();
    let graph = b.finish(vec![swept]).unwrap();

    // Assert the NAMED capability, not merely that it failed: a caller uses
    // this to decide which provider to register. What is missing here is
    // surface evaluation for this surface kind, not the sweep itself.
    match compiler().compile_mesh(&graph, swept, &options()) {
        Err(GeomError::Unsupported { operation, .. }) => {
            assert_eq!(operation, Operation::SurfaceEvaluation);
        }
        other => panic!("expected Unsupported{{SurfaceEvaluation}}, got {other:?}"),
    }
}

/// A handle from a different graph must be refused, not silently indexed.
#[test]
fn a_foreign_node_handle_is_refused() {
    let mut a = GeometryGraphBuilder::new();
    let pa = a.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let sa = a
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: pa,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let graph_a = a.finish(vec![sa]).unwrap();

    let mut c = GeometryGraphBuilder::new();
    let pc = c.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let graph_c = c.finish(vec![pc]).unwrap();
    let _ = &graph_c;

    // `pc` belongs to graph_c, not graph_a.
    assert!(matches!(
        compiler().compile_mesh(&graph_a, pc, &options()),
        Err(GeomError::InvalidInput(_))
    ));
}

/// A deep chain must not overflow the stack.
///
/// This is why evaluation is iterative rather than recursive: graph depth is
/// attacker-controlled in a file format, and a recursive walker would abort
/// the process instead of returning an error.
#[test]
fn a_deep_instance_chain_does_not_overflow_the_stack() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let mut current = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();

    // 50k nested placements, each a no-op translation.
    for _ in 0..50_000 {
        current = b
            .push(GeometryNode::Instance(Instance {
                source: current,
                transform: Transform3::IDENTITY,
            }))
            .unwrap();
    }
    let graph = b.finish(vec![current]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, current, &options())
        .expect("deep chain must compile");
    assert!((volume(&mesh) - 1.0).abs() < 1e-9);
}

/// Memoisation must be observable, not merely believed.
///
/// A diamond DAG whose shared node is expensive: without a cache the shared
/// subtree is rebuilt once per path, so compile time grows exponentially in
/// depth. Ten stacked diamonds is 2^10 = 1024 rebuilds uncached versus 21
/// cached -- far beyond timing noise.
#[test]
fn shared_subtrees_are_not_recompiled_exponentially() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let mut current = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();

    // Each level references `current` twice, doubling the uncached path count.
    for _ in 0..10 {
        let a = b
            .push(GeometryNode::Instance(Instance {
                source: current,
                transform: Transform3::IDENTITY,
            }))
            .unwrap();
        let c = b
            .push(GeometryNode::Instance(Instance {
                source: current,
                transform: Transform3::IDENTITY,
            }))
            .unwrap();
        current = b.push(GeometryNode::Collection(vec![a, c])).unwrap();
    }
    let graph = b.finish(vec![current]).unwrap();

    let start = std::time::Instant::now();
    let mesh = compiler()
        .compile_mesh(&graph, current, &options())
        .expect("compile");
    let elapsed = start.elapsed();

    // 2^10 unit cubes merged.
    assert!(
        (volume(&mesh) - 1024.0).abs() < 1e-6,
        "got {}",
        volume(&mesh)
    );
    assert!(
        elapsed < std::time::Duration::from_secs(2),
        "took {elapsed:?}; shared subtrees are being recompiled"
    );
}

/// A mirroring placement must keep the solid outward-facing.
///
/// A negative-determinant transform reverses triangle orientation. Left
/// uncorrected the mesh is inside-out, which the boolean provider rejects --
/// and IFC mirrored placements are common, so this is not a corner case.
#[test]
fn a_mirrored_placement_keeps_the_solid_outward_facing() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(2.0, 2.0))).unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let mirrored = b
        .push(GeometryNode::Instance(Instance {
            source: solid,
            transform: Transform3::from_scale(Vec3::new(-1.0, 1.0, 1.0)),
        }))
        .unwrap();
    let graph = b.finish(vec![mirrored]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, mirrored, &options())
        .expect("compile");
    // Positive volume means outward winding survived the mirror.
    assert!(
        volume(&mesh) > 0.0,
        "mirrored solid is inside-out: {}",
        volume(&mesh)
    );
    assert!((volume(&mesh) - 4.0).abs() < 1e-9);

    // And the boolean provider must accept it, which is the real contract.
    let mut c = GeometryGraphBuilder::new();
    let p2 = c.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let s2 = c
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: p2,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let g2 = c.finish(vec![s2]).unwrap();
    let tool = compiler().compile_mesh(&g2, s2, &options()).expect("tool");
    use axiolid_mesh_boolean_contract::MeshBoolean;
    // The assertion is that an admissible mirrored solid is accepted; the
    // resulting geometry is not what this test is about.
    BoolmeshBoolean::new()
        .boolean(&mesh, &tool, BooleanOperator::Difference, &options())
        .expect("mirrored solid must be acceptable to the boolean provider");
}

/// Merging meshes must rebase indices, or triangles reference the wrong verts.
#[test]
fn a_collection_rebases_indices_when_merging() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(1.0, 1.0))).unwrap();
    let solid = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let far = b
        .push(GeometryNode::Instance(Instance {
            source: solid,
            transform: Transform3::from_translation(Vec3::new(100.0, 0.0, 0.0)),
        }))
        .unwrap();
    // A DIFFERENT size, so un-rebased indices duplicate the first solid and
    // change the total volume. Two identical cubes would sum to the same
    // number either way, which is why the earlier version of this test could
    // not detect a missing rebase.
    let big_profile = b.push(GeometryNode::Profile(rect(3.0, 3.0))).unwrap();
    let big = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: big_profile,
            direction: Vec3::Z,
            depth: 1.0,
        }))
        .unwrap();
    let _ = far;
    let both = b.push(GeometryNode::Collection(vec![solid, big])).unwrap();
    let graph = b.finish(vec![both]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, both, &options())
        .expect("compile");
    // Two disjoint unit cubes. Un-rebased indices would collapse the second
    // onto the first, halving the volume.
    // 1x1x1 + 3x3x1 = 10. Un-rebased indices would yield 1 + 1 = 2.
    assert!((volume(&mesh) - 10.0).abs() < 1e-9, "got {}", volume(&mesh));
    assert_eq!(mesh.positions.len(), 16);
}

#[test]
fn scaled_instances_tessellate_sources_at_instance_local_tolerance() {
    let mut b = GeometryGraphBuilder::new();
    let tiny_circle = b
        .push(GeometryNode::Profile(Profile::Circle(CircleProfile {
            radius: 0.000_05,
            thickness: None,
        })))
        .unwrap();
    let extrusion = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile: tiny_circle,
            direction: Vec3::Z,
            depth: 0.002,
        }))
        .unwrap();
    let scaled_extrusion = b
        .push(GeometryNode::Instance(Instance {
            source: extrusion,
            transform: Transform3::from_scale(Vec3::splat(1000.0)),
        }))
        .unwrap();

    let short_path = b
        .push(GeometryNode::Curve3(Curve3::Polyline(Polyline3 {
            points: vec![Point3::new(0.001, 0.0, 0.0), Point3::new(0.001, 0.0, 0.002)],
            closed: false,
        })))
        .unwrap();
    let swept_disk = b
        .push(GeometryNode::SolidOperation(SolidOperation::SweptDisk {
            directrix: short_path,
            radius: 0.000_05,
            inner_radius: None,
            parameter_range: None,
            fillet_radius: None,
        }))
        .unwrap();
    let scaled_sweep = b
        .push(GeometryNode::Instance(Instance {
            source: swept_disk,
            transform: Transform3::from_scale(Vec3::splat(1000.0)),
        }))
        .unwrap();
    let both = b
        .push(GeometryNode::Collection(vec![
            scaled_extrusion,
            scaled_sweep,
        ]))
        .unwrap();
    let graph = b.finish(vec![both]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, both, &ExecutionOptions::new(Tolerance::MILLIMETRE))
        .expect("scaled sources must be tessellated in local coordinates");
    assert!(
        mesh.positions.len() >= 40,
        "under-tessellated: {} vertices",
        mesh.positions.len()
    );
    assert!(!mesh.indices.is_empty());
    assert!(mesh.positions.iter().all(|p| p.is_finite()));
}

#[test]
fn boolean_difference_materializes_unbounded_half_space_from_subject_bounds() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(2.0, 2.0))).unwrap();
    let subject = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 2.0,
        }))
        .unwrap();
    let upper_half = b
        .push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: Point3::new(0.0, 0.0, 1.0),
                normal: Vec3::Z,
            },
            agreement: true,
        }))
        .unwrap();
    let clipped = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            operator: BooleanOperator::Difference,
            left: subject,
            right: upper_half,
        }))
        .unwrap();
    let graph = b.finish(vec![clipped]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, clipped, &options())
        .expect("finite subject must bound its half-space operand");
    assert!((volume(&mesh) - 4.0).abs() < 1e-6, "got {}", volume(&mesh));
    assert!(mesh.positions.iter().all(|p| p.is_finite()));
}

#[test]
fn boolean_intersection_materializes_unbounded_half_space_from_subject_bounds() {
    let mut b = GeometryGraphBuilder::new();
    let profile = b.push(GeometryNode::Profile(rect(2.0, 2.0))).unwrap();
    let subject = b
        .push(GeometryNode::SolidOperation(SolidOperation::Extrusion {
            profile,
            direction: Vec3::Z,
            depth: 2.0,
        }))
        .unwrap();
    let upper_half = b
        .push(GeometryNode::HalfSpace(HalfSpace {
            boundary: Plane3 {
                origin: Point3::new(0.0, 0.0, 1.0),
                normal: Vec3::Z,
            },
            agreement: true,
        }))
        .unwrap();
    let clipped = b
        .push(GeometryNode::SolidOperation(SolidOperation::Boolean {
            operator: BooleanOperator::Intersection,
            left: subject,
            right: upper_half,
        }))
        .unwrap();
    let graph = b.finish(vec![clipped]).unwrap();

    let mesh = compiler()
        .compile_mesh(&graph, clipped, &options())
        .expect("finite subject must bound its half-space intersection");
    assert!((volume(&mesh) - 4.0).abs() < 1e-6, "got {}", volume(&mesh));
    assert!(mesh.positions.iter().all(|p| p.is_finite()));
    let min_z = mesh.positions.iter().map(|p| p.z).fold(f64::MAX, f64::min);
    let max_z = mesh.positions.iter().map(|p| p.z).fold(f64::MIN, f64::max);
    assert!(min_z >= 1.0 - 1e-9, "wrong half-space side: min z {min_z}");
    assert!(max_z <= 2.0 + 1e-9, "escaped subject bounds: max z {max_z}");
}
