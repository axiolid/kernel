//! Backend-neutral mesh visible-face contract tests.
//!
//! Every pre-dispatch refusal must happen without calling the provider and
//! name the mesh at fault; every way a provider can break the outcome contract
//! must be caught; only unsupported or unavailable providers permit fallback.
//! `scripts/probe_visible_faces_registry_mutants.py` breaks each registry
//! check in turn and requires one of these tests to fail.

#![cfg(feature = "mesh-hidden-lines")]

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, CancellationToken,
    ExecutionOptions, ExecutionTarget, GeomError, Operation, ScratchRequirement,
};
use axiolid_core::{Aabb2, Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_dispatch::MeshVisibleFacesRegistry;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    BehindSurface, FaceGroup, FaceSource, HiddenLineError, HiddenLineMesh, MeshId, MeshRole,
    MeshVisibleFaces, Occlusion, View, ViewProjection, VisibleFaceEvidence, VisibleFaceLimits,
    VisibleFaceOutcome, VisibleFaceRequest, VisibleRegion,
};

fn cube() -> TriMesh {
    TriMesh::new(
        vec![
            [0.0, 0.0, 0.0].into(),
            [1.0, 0.0, 0.0].into(),
            [1.0, 1.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            [0.0, 0.0, 1.0].into(),
            [1.0, 0.0, 1.0].into(),
            [1.0, 1.0, 1.0].into(),
            [0.0, 1.0, 1.0].into(),
        ],
        vec![
            0, 2, 1, 0, 3, 2, // bottom
            4, 5, 6, 4, 6, 7, // top
            0, 1, 5, 0, 5, 4, // front
            1, 2, 6, 1, 6, 5, // right
            2, 3, 7, 2, 7, 6, // back
            3, 0, 4, 3, 4, 7, // left
        ],
    )
}

fn quad() -> TriMesh {
    TriMesh::new(
        vec![Vec3::ZERO, Vec3::X, Vec3::new(1.0, 1.0, 0.0), Vec3::Y],
        vec![0, 1, 2, 0, 2, 3],
    )
}

fn frame() -> Frame3 {
    Frame3 {
        origin: Vec3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn limits() -> VisibleFaceLimits {
    VisibleFaceLimits::new(8, 1_000, 1_000, 1_000_000, 1_000, 10_000)
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn view() -> View {
    View::new(frame())
}

fn mesh<'a>(id: u64, role: MeshRole, mesh: &'a TriMesh) -> HiddenLineMesh<'a> {
    HiddenLineMesh::new(MeshId(id), role, mesh)
}

/// Whether a refusal is the expected kind.
type Accepts = fn(&GeomError) -> bool;

/// What a recorder does after it has been called.
type Behaviour =
    fn(&VisibleFaceRequest<'_>, &mut VisibleFaceOutcome) -> Result<(), HiddenLineError>;

fn well_behaved(
    _request: &VisibleFaceRequest<'_>,
    _outcome: &mut VisibleFaceOutcome,
) -> Result<(), HiddenLineError> {
    Ok(())
}

#[derive(Debug)]
struct Recorder {
    id: &'static str,
    calls: Arc<AtomicUsize>,
    scratch: ScratchRequirement,
    behaviour: Behaviour,
    perspective: bool,
    see_through: bool,
}

impl Recorder {
    fn new(id: &'static str, calls: &Arc<AtomicUsize>, behaviour: Behaviour) -> Self {
        Self {
            id,
            calls: Arc::clone(calls),
            scratch: ScratchRequirement::None,
            behaviour,
            perspective: false,
            see_through: false,
        }
    }

    fn claiming_see_through(mut self) -> Self {
        self.see_through = true;
        self
    }

    fn claiming_perspective(mut self) -> Self {
        self.perspective = true;
        self
    }
}

impl Backend for Recorder {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new(self.id), ExecutionTarget::PortableCpu)
    }
}

fn triangle_ring() -> Vec<Point2> {
    vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
    ]
}

fn region(mesh: u64, source: FaceSource, outer: Vec<Point2>, triangles: &[u32]) -> VisibleRegion {
    VisibleRegion::new(MeshId(mesh), source, outer, Vec::new(), triangles.to_vec())
}

/// One valid region on the first mesh: group 0, the top face.
fn base_outcome(request: &VisibleFaceRequest<'_>) -> VisibleFaceOutcome {
    let triangles = request.meshes.iter().map(|m| m.mesh.triangle_count()).sum();
    let mut evidence = VisibleFaceEvidence::input_mesh(request.meshes.len(), triangles);
    evidence.see_through_triangles = request
        .meshes
        .iter()
        .map(|m| m.occlusion.see_through_count())
        .sum();
    evidence.regions = 1;
    VisibleFaceOutcome::new(
        request.view,
        vec![region(
            request.meshes[0].id.0,
            FaceSource::Group(FaceGroup(0)),
            triangle_ring(),
            &[2, 3],
        )],
        evidence,
    )
}

impl MeshVisibleFaces for Recorder {
    fn supports_projection(&self, projection: ViewProjection) -> bool {
        projection == ViewProjection::Orthographic || self.perspective
    }

    fn supports_occlusion(&self, occlusion: Occlusion) -> bool {
        occlusion == Occlusion::Opaque || self.see_through
    }

    fn scratch_requirement(&self) -> ScratchRequirement {
        self.scratch
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError> {
        options.check_cancelled()?;
        self.calls.fetch_add(1, Ordering::SeqCst);
        let mut outcome = base_outcome(request);
        (self.behaviour)(request, &mut outcome)?;
        Ok(outcome)
    }
}

fn registry(calls: &Arc<AtomicUsize>, behaviour: Behaviour) -> MeshVisibleFacesRegistry {
    let mut registry = MeshVisibleFacesRegistry::new();
    registry.register(0, Recorder::new("visible-faces-recorder", calls, behaviour));
    registry
}

fn counter() -> Arc<AtomicUsize> {
    Arc::new(AtomicUsize::new(0))
}

fn calls_of(calls: &Arc<AtomicUsize>) -> usize {
    calls.load(Ordering::SeqCst)
}

#[test]
fn dispatch_returns_regions_with_mesh_provenance() {
    let calls = counter();
    let cube = cube();
    let meshes = [mesh(7, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let outcome = registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect("valid visible-face contract");

    assert_eq!(calls_of(&calls), 1);
    assert_eq!(outcome.view, view());
    assert_eq!(outcome.regions.len(), 1);
    assert_eq!(outcome.regions[0].mesh, MeshId(7));
    assert_eq!(outcome.evidence.source_meshes, 1);
    assert_eq!(outcome.evidence.source_triangles, 12);
    assert!(outcome.evidence.is_derived_from_input_mesh());
}

/// A good first mesh (id 10) and a second mesh (id 20) that may be at fault.
fn two_mesh_refusal(
    offender: &TriMesh,
    role: MeshRole,
    groups: Option<&[FaceGroup]>,
) -> HiddenLineError {
    let calls = counter();
    let good = cube();
    let mut bad = mesh(20, role, offender);
    if let Some(groups) = groups {
        bad = bad.with_groups(groups);
    }
    let meshes = [mesh(10, MeshRole::Solid, &good), bad];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let error = registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect_err("the second mesh must be refused");
    assert_eq!(calls_of(&calls), 0, "refused before dispatch");
    error
}

#[test]
fn each_per_mesh_refusal_names_the_offending_mesh() {
    let origin = Point3::ZERO;
    type Case = (&'static str, TriMesh, MeshRole, Accepts);
    let cases: Vec<Case> = vec![
        (
            "no triangles",
            TriMesh::new(vec![origin, Vec3::X, Vec3::Y], vec![]),
            MeshRole::Sheet,
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "partial triangle",
            TriMesh::new(vec![origin, Vec3::X, Vec3::Y], vec![0, 1, 2, 0]),
            MeshRole::Sheet,
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "invalid index",
            TriMesh::new(vec![origin, Vec3::X, Vec3::Y], vec![0, 1, 5]),
            MeshRole::Sheet,
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "non-finite position",
            TriMesh::new(
                vec![origin, Vec3::X, Vec3::new(f64::NAN, 1.0, 0.0)],
                vec![0, 1, 2],
            ),
            MeshRole::Sheet,
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "degenerate triangle",
            TriMesh::new(
                vec![origin, Vec3::X, Vec3::new(2.0, 0.0, 0.0)],
                vec![0, 1, 2],
            ),
            MeshRole::Sheet,
            |e| matches!(e, GeomError::Degenerate(_)),
        ),
        (
            "inconsistent winding",
            TriMesh::new(
                vec![origin, Vec3::X, Vec3::Y, Vec3::new(0.0, -1.0, 0.0)],
                vec![0, 1, 2, 0, 1, 3],
            ),
            MeshRole::Sheet,
            |e| matches!(e, GeomError::NotManifold(_)),
        ),
        ("open solid", quad(), MeshRole::Solid, |e| {
            matches!(e, GeomError::NotManifold(_))
        }),
    ];
    for (name, offender, role, accepts) in cases {
        let error = two_mesh_refusal(&offender, role, None);
        assert!(accepts(&error.error), "{name}: {error:?}");
        assert_eq!(error.mesh, Some(MeshId(20)), "{name}");
        assert!(
            error.to_string().starts_with("mesh 20: "),
            "{name}: {error}"
        );
    }
}

#[test]
fn a_group_count_other_than_the_triangle_count_names_the_mesh() {
    // A quad has two triangles.
    let sheet = quad();
    for len in [0, 1, 3] {
        let groups = vec![FaceGroup(0); len];
        let error = two_mesh_refusal(&sheet, MeshRole::Sheet, Some(&groups));
        assert!(matches!(error.error, GeomError::InvalidInput(_)), "{len}");
        assert_eq!(error.mesh, Some(MeshId(20)), "{len}");
    }
}

#[test]
fn the_right_number_of_groups_is_accepted() {
    let calls = counter();
    let (good, sheet) = (cube(), quad());
    let groups = [FaceGroup(3), FaceGroup(4)];
    let meshes = [
        mesh(10, MeshRole::Solid, &good),
        mesh(20, MeshRole::Sheet, &sheet).with_groups(&groups),
    ];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect("one group per triangle");
    assert_eq!(calls_of(&calls), 1);
}

#[test]
fn a_duplicate_mesh_id_names_the_second_occurrence() {
    let calls = counter();
    let (a, b) = (cube(), cube());
    let meshes = [
        mesh(10, MeshRole::Solid, &a),
        mesh(20, MeshRole::Solid, &b),
        mesh(20, MeshRole::Solid, &a),
    ];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let error = registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .unwrap_err();
    assert!(matches!(error.error, GeomError::InvalidInput(_)));
    assert_eq!(error.mesh, Some(MeshId(20)));
    assert_eq!(calls_of(&calls), 0);
}

/// Requests refused as a whole: the refusal names no mesh.
#[test]
fn each_whole_request_refusal_names_no_mesh() {
    type Edit = fn(&mut VisibleFaceRequest<'_>);
    let cases: Vec<(&str, Edit, Accepts)> = vec![
        (
            "empty request",
            |r| r.meshes = &[],
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "too many meshes",
            |r| r.limits.max_meshes = 1,
            |e| {
                matches!(
                    e,
                    GeomError::BudgetExceeded {
                        resource: "visible-face meshes"
                    }
                )
            },
        ),
        (
            "too many vertices",
            |r| r.limits.max_source_vertices = 15,
            |e| {
                matches!(
                    e,
                    GeomError::BudgetExceeded {
                        resource: "visible-face source vertices"
                    }
                )
            },
        ),
        (
            "too many triangles",
            |r| r.limits.max_source_triangles = 23,
            |e| {
                matches!(
                    e,
                    GeomError::BudgetExceeded {
                        resource: "visible-face source triangles"
                    }
                )
            },
        ),
        (
            "dirty frame",
            |r| r.view.frame.z = Vec3::new(0.0, 0.0, 2.0),
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "non-finite near",
            |r| r.view = view().with_near(f64::NAN),
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "non-finite far",
            |r| r.view = view().with_far(f64::NEG_INFINITY),
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "far not below near",
            |r| r.view = view().with_near(1.0).with_far(1.0),
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "inverted region",
            |r| {
                r.region = Some(Aabb2 {
                    min: Point2::new(1.0, 0.0),
                    max: Point2::new(0.0, 1.0),
                })
            },
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "non-finite region",
            |r| {
                r.region = Some(Aabb2 {
                    min: Point2::new(0.0, 0.0),
                    max: Point2::new(f64::INFINITY, 1.0),
                })
            },
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
    ];
    let (a, b) = (cube(), cube());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Solid, &b)];
    for (name, edit, accepts) in cases {
        let calls = counter();
        let mut request = VisibleFaceRequest::new(&meshes, view(), limits());
        edit(&mut request);
        let error = registry(&calls, well_behaved)
            .visible_faces(&request, &options())
            .expect_err(name);
        assert!(accepts(&error.error), "{name}: {error:?}");
        assert_eq!(error.mesh, None, "{name}");
        assert_eq!(calls_of(&calls), 0, "{name}");
    }
}

#[test]
fn boundary_values_are_accepted() {
    let calls = counter();
    let (a, b) = (cube(), cube());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Solid, &b)];
    // Exactly at every limit, and an empty-area region. The recorder's one
    // region has three ring points and no candidate pairs.
    let tight = VisibleFaceLimits::new(2, 16, 24, 0, 1, 3);
    let request = VisibleFaceRequest::new(&meshes, view(), tight).with_region(Aabb2 {
        min: Point2::new(0.5, 0.5),
        max: Point2::new(0.5, 0.5),
    });
    registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect("limits are inclusive");
    assert_eq!(calls_of(&calls), 1);
}

#[test]
fn cancellation_refuses_the_whole_request() {
    let calls = counter();
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let token = CancellationToken::new();
    token.cancel();
    let error = registry(&calls, well_behaved)
        .visible_faces(&request, &options().with_cancellation(token))
        .unwrap_err();
    assert_eq!(error, HiddenLineError::request(GeomError::Cancelled));
    assert_eq!(calls_of(&calls), 0);
}

// --- outcome validation -------------------------------------------------

fn assert_violation(request: &VisibleFaceRequest<'_>, behaviour: Behaviour, name: &str) {
    let calls = counter();
    let error = registry(&calls, behaviour)
        .visible_faces(request, &options())
        .expect_err(name);
    assert_eq!(calls_of(&calls), 1, "{name}");
    assert!(
        matches!(error.error, GeomError::BackendContractViolation { .. }),
        "{name}: {error:?}"
    );
    assert_eq!(error.mesh, None, "{name}");
}

#[test]
fn outcome_violations_are_caught() {
    let (a, b) = (cube(), cube());
    // Triangle 2 of mesh 1 is in group 1, so the base region (group 0) lists
    // a triangle of another group.
    let mut groups = vec![FaceGroup(0); 12];
    groups[2] = FaceGroup(1);
    let meshes = [
        mesh(1, MeshRole::Solid, &a).with_groups(&groups),
        mesh(2, MeshRole::Solid, &b),
    ];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let plain = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Solid, &b)];
    let plain_request = VisibleFaceRequest::new(&plain, view(), limits());

    // Against the grouped request the base region is itself a violation.
    assert_violation(&request, well_behaved, "triangle of another group");

    let cases: Vec<(&str, Behaviour)> = vec![
        ("view changed", |_, o| {
            o.view = o.view.with_near(5.0);
            Ok(())
        }),
        ("unknown mesh", |_, o| {
            o.regions[0].mesh = MeshId(99);
            Ok(())
        }),
        ("outer ring of two points", |_, o| {
            o.regions[0].outer.truncate(2);
            Ok(())
        }),
        ("equal consecutive outer points", |_, o| {
            o.regions[0].outer.insert(1, Point2::new(0.0, 0.0));
            Ok(())
        }),
        ("equal wrap-around outer points", |_, o| {
            o.regions[0].outer.push(Point2::new(0.0, 0.0));
            Ok(())
        }),
        ("equal consecutive hole points", |_, o| {
            o.regions[0].holes = vec![vec![
                Point2::new(0.2, 0.1),
                Point2::new(0.2, 0.1),
                Point2::new(0.3, 0.2),
            ]];
            Ok(())
        }),
        ("non-finite outer point", |_, o| {
            o.regions[0].outer[1] = Point2::new(f64::NAN, 0.0);
            Ok(())
        }),
        ("infinite outer point", |_, o| {
            o.regions[0].outer[1] = Point2::new(0.0, f64::INFINITY);
            Ok(())
        }),
        ("hole of two points", |_, o| {
            o.regions[0].holes = vec![vec![Point2::new(0.2, 0.1), Point2::new(0.3, 0.1)]];
            Ok(())
        }),
        ("non-finite hole point", |_, o| {
            o.regions[0].holes = vec![vec![
                Point2::new(0.2, 0.1),
                Point2::new(0.3, 0.1),
                Point2::new(f64::NAN, 0.2),
            ]];
            Ok(())
        }),
        ("triangles descending", |_, o| {
            o.regions[0].triangles = vec![3, 2];
            Ok(())
        }),
        ("triangle repeated", |_, o| {
            o.regions[0].triangles = vec![2, 2];
            Ok(())
        }),
        ("triangle out of range", |_, o| {
            o.regions[0].triangles = vec![2, 12];
            Ok(())
        }),
        ("cap without a near plane", |_, o| {
            o.regions[0].source = FaceSource::Cap;
            o.regions[0].triangles.clear();
            o.evidence.cap_regions = 1;
            Ok(())
        }),
        ("regions out of order across meshes", |_, o| {
            o.regions.insert(
                0,
                region(2, FaceSource::Group(FaceGroup(0)), triangle_ring(), &[2]),
            );
            o.evidence.regions = 2;
            Ok(())
        }),
        ("regions out of order across sources", |_, o| {
            o.regions[0].source = FaceSource::Group(FaceGroup(5));
            o.regions.push(region(
                1,
                FaceSource::Group(FaceGroup(1)),
                triangle_ring(),
                &[],
            ));
            o.evidence.regions = 2;
            Ok(())
        }),
        ("regions out of order within a source", |_, o| {
            let mut later = triangle_ring();
            later[0] = Point2::new(0.5, 0.0);
            o.regions
                .insert(0, region(1, FaceSource::Group(FaceGroup(0)), later, &[]));
            o.evidence.regions = 2;
            Ok(())
        }),
        ("one region twice", |_, o| {
            let again = o.regions[0].clone();
            o.regions.push(again);
            o.evidence.regions = 2;
            Ok(())
        }),
        ("source triangles", |_, o| {
            o.evidence.source_triangles += 1;
            Ok(())
        }),
        ("source meshes", |_, o| {
            o.evidence.source_meshes += 1;
            Ok(())
        }),
        ("region count", |_, o| {
            o.evidence.regions = 0;
            Ok(())
        }),
        ("cap count", |_, o| {
            o.evidence.cap_regions = 1;
            Ok(())
        }),
    ];
    for (name, behaviour) in cases {
        assert_violation(&plain_request, behaviour, name);
    }
}

#[test]
fn a_triangle_of_the_regions_group_is_accepted() {
    let calls = counter();
    let a = cube();
    // Triangles 2 and 3 are group 4, and the base region is group 0, so
    // retarget the region to group 4.
    let mut groups = vec![FaceGroup(0); 12];
    groups[2] = FaceGroup(4);
    groups[3] = FaceGroup(4);
    let meshes = [mesh(1, MeshRole::Solid, &a).with_groups(&groups)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    registry(&calls, |_, o| {
        o.regions[0].source = FaceSource::Group(FaceGroup(4));
        Ok(())
    })
    .visible_faces(&request, &options())
    .expect("both triangles are in group 4");
    assert_eq!(calls_of(&calls), 1);
}

#[test]
fn a_cap_of_a_solid_with_a_near_plane_is_accepted() {
    let calls = counter();
    let (a, b) = (cube(), cube());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Sheet, &b)];
    let request = VisibleFaceRequest::new(&meshes, view().with_near(0.5), limits());
    let outcome = registry(&calls, |_, o| {
        // Sorted: the group before the cap on mesh 1, then mesh 2.
        o.regions
            .push(region(1, FaceSource::Cap, triangle_ring(), &[]));
        o.regions.push(region(
            2,
            FaceSource::Group(FaceGroup(0)),
            triangle_ring(),
            &[0],
        ));
        o.evidence.regions = 3;
        o.evidence.cap_regions = 1;
        Ok(())
    })
    .visible_faces(&request, &options())
    .expect("ordered output with a cap");
    assert_eq!(outcome.regions.len(), 3);
    assert_eq!(outcome.evidence.cap_regions, 1);
}

#[test]
fn a_cap_of_a_sheet_with_a_near_plane_is_refused() {
    let b = quad();
    let meshes = [mesh(1, MeshRole::Sheet, &b)];
    let request = VisibleFaceRequest::new(&meshes, view().with_near(0.5), limits());
    assert_violation(
        &request,
        |_, o| {
            o.regions[0].source = FaceSource::Cap;
            o.regions[0].triangles.clear();
            o.evidence.cap_regions = 1;
            Ok(())
        },
        "cap of a sheet",
    );
}

#[test]
fn a_cap_that_lists_triangles_is_refused() {
    let a = cube();
    let meshes = [mesh(1, MeshRole::Solid, &a)];
    let request = VisibleFaceRequest::new(&meshes, view().with_near(0.5), limits());
    assert_violation(
        &request,
        |_, o| {
            o.regions[0].source = FaceSource::Cap;
            o.evidence.cap_regions = 1;
            Ok(())
        },
        "cap with triangles",
    );
}

#[test]
fn outcomes_over_the_output_limits_exceed_the_budget() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let calls = counter();

    // One region of three points fits exactly; a second region does not.
    let mut tight = limits();
    tight.max_regions = 1;
    tight.max_region_vertices = 3;
    let request = VisibleFaceRequest::new(&meshes, view(), tight);
    registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect("at the limits");
    let mut roomy_vertices = tight;
    roomy_vertices.max_region_vertices = 100;
    let mut no_pairs = tight;
    no_pairs.max_candidate_pairs = 0;

    let over = |behaviour: Behaviour, resource: &'static str, request: &VisibleFaceRequest<'_>| {
        let calls = counter();
        let error = registry(&calls, behaviour)
            .visible_faces(request, &options())
            .unwrap_err();
        assert_eq!(calls_of(&calls), 1, "{resource}");
        assert_eq!(error.error, GeomError::BudgetExceeded { resource });
        assert_eq!(error.mesh, None, "{resource}");
    };
    over(
        |_, o| {
            o.regions.push(region(
                1,
                FaceSource::Group(FaceGroup(1)),
                triangle_ring(),
                &[],
            ));
            o.evidence.regions = 2;
            Ok(())
        },
        "visible-face regions",
        &VisibleFaceRequest::new(&meshes, view(), roomy_vertices),
    );
    over(
        |_, o| {
            o.regions[0].holes = vec![triangle_ring()];
            Ok(())
        },
        "visible-face region vertices",
        &request,
    );
    over(
        |_, o| {
            o.evidence.candidate_pairs = 1;
            Ok(())
        },
        "visible-face candidate pairs",
        &VisibleFaceRequest::new(&meshes, view(), no_pairs),
    );
}

// --- fallback ----------------------------------------------------------

fn unsupported(
    _: &VisibleFaceRequest<'_>,
    _: &mut VisibleFaceOutcome,
) -> Result<(), HiddenLineError> {
    Err(GeomError::Unsupported {
        backend: BackendId::new("first"),
        operation: Operation::MeshVisibleFaces,
    }
    .into())
}

fn unavailable(
    _: &VisibleFaceRequest<'_>,
    _: &mut VisibleFaceOutcome,
) -> Result<(), HiddenLineError> {
    Err(GeomError::Unavailable {
        backend: BackendId::new("first"),
        reason: "device lost".into(),
    }
    .into())
}

#[test]
fn only_unsupported_or_unavailable_providers_fall_back() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());

    for retryable in [unsupported as Behaviour, unavailable] {
        let (first, second) = (counter(), counter());
        let mut registry = MeshVisibleFacesRegistry::new();
        registry.register(10, Recorder::new("first", &first, retryable));
        registry.register(0, Recorder::new("second", &second, well_behaved));
        registry
            .visible_faces(&request, &options())
            .expect("falls back to the second provider");
        assert_eq!((calls_of(&first), calls_of(&second)), (1, 1));
    }

    let refusals: [(&str, Behaviour); 4] = [
        ("invalid input", |_, _| {
            Err(HiddenLineError::mesh(
                MeshId(1),
                GeomError::InvalidInput("provider refusal".into()),
            ))
        }),
        ("degenerate", |_, _| {
            Err(HiddenLineError::mesh(
                MeshId(1),
                GeomError::Degenerate("visible region collapsed under rounding".into()),
            ))
        }),
        ("budget", |_, _| {
            Err(GeomError::BudgetExceeded {
                resource: "visible-face candidate pairs",
            }
            .into())
        }),
        ("unsupported input", |_, _| {
            Err(GeomError::UnsupportedInput {
                backend: BackendId::new("first"),
                operation: Operation::MeshVisibleFaces,
                input: "curved surface",
            }
            .into())
        }),
    ];
    for (name, refusal) in refusals {
        let (first, second) = (counter(), counter());
        let mut registry = MeshVisibleFacesRegistry::new();
        registry.register(10, Recorder::new("first", &first, refusal));
        registry.register(0, Recorder::new("second", &second, well_behaved));
        let error = registry.visible_faces(&request, &options()).unwrap_err();
        assert_eq!(calls_of(&second), 0, "{name}: must not fall back");
        assert!(
            matches!(
                error.error,
                GeomError::InvalidInput(_)
                    | GeomError::Degenerate(_)
                    | GeomError::BudgetExceeded { .. }
                    | GeomError::UnsupportedInput { .. }
            ),
            "{name}: {error:?}"
        );
    }
}

#[test]
fn a_provider_refusal_keeps_its_mesh_through_dispatch() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let calls = counter();
    let error = registry(&calls, |_, _| {
        Err(HiddenLineError::mesh(
            MeshId(1),
            GeomError::Degenerate("visible region collapsed under rounding".into()),
        ))
    })
    .visible_faces(&request, &options())
    .unwrap_err();
    assert_eq!(error.mesh, Some(MeshId(1)));
}

#[test]
fn without_a_matching_provider_the_request_is_unsupported() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let error = MeshVisibleFacesRegistry::new()
        .visible_faces(&request, &options())
        .unwrap_err();
    assert!(matches!(
        error.error,
        GeomError::Unsupported {
            operation: Operation::MeshVisibleFaces,
            ..
        }
    ));
    assert_eq!(error.mesh, None);

    // Every provider retryable: the last retryable error is returned.
    let calls = counter();
    let error = registry(&calls, unsupported)
        .visible_faces(&request, &options())
        .unwrap_err();
    assert!(matches!(error.error, GeomError::Unsupported { .. }));
    assert_eq!(calls_of(&calls), 1);
}

#[test]
fn register_conformant_rejects_a_provider_that_fails_the_suite() {
    let mut registry = MeshVisibleFacesRegistry::new();
    let calls = counter();
    let report = registry
        .register_conformant(0, Recorder::new("recorder", &calls, well_behaved))
        .expect_err("a recorder finds no faces");
    assert!(!report.is_success());
    assert_eq!(registry.providers().count(), 0);
}

// --- budgets ------------------------------------------------------------

fn budgeted(scratch: ScratchRequirement) -> (MeshVisibleFacesRegistry, Arc<AtomicUsize>) {
    let calls = counter();
    let mut registry = MeshVisibleFacesRegistry::new();
    let mut recorder = Recorder::new("budgeted", &calls, well_behaved);
    recorder.scratch = scratch;
    registry.register(0, recorder);
    (registry, calls)
}

#[test]
fn memory_is_admitted_before_the_allocating_topology_audit() {
    let (registry, calls) = budgeted(ScratchRequirement::None);
    let open = quad();
    let meshes = [mesh(1, MeshRole::Solid, &open)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let error = registry
        .visible_faces(&request, &options().with_memory_budget(0))
        .unwrap_err();
    assert_eq!(
        error.error,
        GeomError::BudgetExceeded { resource: "memory" }
    );
    assert_eq!(calls_of(&calls), 0);
}

#[test]
fn provider_scratch_counts_every_mesh_and_the_larger_of_vertices_and_triangles() {
    let (registry, calls) = budgeted(ScratchRequirement::PerElement {
        bytes_per_element: 1_000,
    });
    let (a, b) = (cube(), quad());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Sheet, &b)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    // 12 + 2 triangles, 8 + 4 vertices: 14 elements.
    registry
        .visible_faces(&request, &options().with_memory_budget(14_000))
        .expect("14 elements fit");
    assert_eq!(calls_of(&calls), 1);
    let error = registry
        .visible_faces(&request, &options().with_memory_budget(13_999))
        .unwrap_err();
    assert_eq!(
        error.error,
        GeomError::BudgetExceeded { resource: "memory" }
    );
    assert_eq!(calls_of(&calls), 1);
}

#[test]
fn a_refusal_naming_a_mesh_outside_the_request_is_a_contract_violation() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let calls = counter();
    let error = registry(&calls, |_, _| {
        Err(HiddenLineError::mesh(
            MeshId(99),
            GeomError::Degenerate("elsewhere".into()),
        ))
    })
    .visible_faces(&request, &options())
    .unwrap_err();
    assert!(
        matches!(error.error, GeomError::BackendContractViolation { .. }),
        "{error:?}"
    );
    assert_eq!(error.mesh, None);
    assert_eq!(calls_of(&calls), 1);
}

/// A planar group is geometry the registry does not recompute: it checks the
/// structure, an id that is a triangle of the mesh and no later than any
/// triangle the region lists (the id is the group's smallest triangle).
#[test]
fn a_planar_group_is_checked_for_structure_only() {
    let a = cube();
    let meshes = [mesh(1, MeshRole::Solid, &a).with_planar_groups()];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let accepted: [(&str, Behaviour); 2] = [
        ("the id is the listed triangle", |_, o| {
            o.regions[0].source = FaceSource::Group(FaceGroup(2));
            Ok(())
        }),
        ("the id is before the listed triangle", |_, o| {
            o.regions[0].source = FaceSource::Group(FaceGroup(0));
            Ok(())
        }),
    ];
    for (name, behaviour) in accepted {
        let calls = counter();
        registry(&calls, behaviour)
            .visible_faces(&request, &options())
            .expect(name);
    }
    let refused: [(&str, Behaviour); 3] = [
        ("the id is after the listed triangle", |_, o| {
            o.regions[0].source = FaceSource::Group(FaceGroup(3));
            Ok(())
        }),
        ("the id is not a triangle", |_, o| {
            o.regions[0].source = FaceSource::Group(FaceGroup(12));
            o.regions[0].triangles.clear();
            Ok(())
        }),
        ("the id is not a triangle and one is listed", |_, o| {
            o.regions[0].source = FaceSource::Group(FaceGroup(40));
            Ok(())
        }),
    ];
    for (name, behaviour) in refused {
        assert_violation(&request, behaviour, name);
    }
}

fn perspective_view() -> View {
    View::perspective(frame(), -1.0).expect("valid perspective view")
}

#[test]
fn a_perspective_view_reaches_only_providers_that_claim_it() {
    let (declining, claiming) = (counter(), counter());
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, perspective_view(), limits());

    let mut registry = MeshVisibleFacesRegistry::new();
    registry.register(1, Recorder::new("declines", &declining, well_behaved));
    registry.register(
        0,
        Recorder::new("claims", &claiming, well_behaved).claiming_perspective(),
    );
    let outcome = registry
        .visible_faces(&request, &options())
        .expect("the claiming provider answers");
    assert_eq!(outcome.view, perspective_view());
    assert_eq!((calls_of(&declining), calls_of(&claiming)), (0, 1));

    let ortho = VisibleFaceRequest::new(&meshes, view(), limits());
    registry.visible_faces(&ortho, &options()).expect("answers");
    assert_eq!((calls_of(&declining), calls_of(&claiming)), (1, 1));
}

#[test]
fn a_perspective_view_without_a_claiming_provider_is_unsupported() {
    let calls = counter();
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = VisibleFaceRequest::new(&meshes, perspective_view(), limits());
    let error = registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect_err("no provider supports the projection");
    assert_eq!(calls_of(&calls), 0);
    assert_eq!(error.mesh, None);
    assert!(matches!(
        error.error,
        GeomError::Unsupported {
            operation: Operation::MeshVisibleFaces,
            ..
        }
    ));
}

#[test]
fn a_perspective_view_needs_a_negative_near_even_when_built_by_hand() {
    let calls = counter();
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let mut registry = MeshVisibleFacesRegistry::new();
    registry.register(
        0,
        Recorder::new("claims", &calls, well_behaved).claiming_perspective(),
    );
    let good = perspective_view();
    let mut bad_near = good;
    bad_near.near = Some(0.5);
    let mut no_near = good;
    no_near.near = None;
    for view in [bad_near, no_near] {
        let request = VisibleFaceRequest::new(&meshes, view, limits());
        let error = registry
            .visible_faces(&request, &options())
            .expect_err("refused before dispatch");
        assert!(matches!(error.error, GeomError::InvalidInput(_)));
        assert_eq!(error.mesh, None);
    }
    assert_eq!(calls_of(&calls), 0);
}

#[test]
fn routing_by_projection_comes_before_the_topology_audit() {
    let calls = counter();
    let bad = TriMesh::new(
        vec![Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::new(1.0, 1.0, 0.0)],
        vec![0, 1, 2, 0, 1, 3],
    );
    let meshes = [mesh(1, MeshRole::Sheet, &bad)];
    let request = VisibleFaceRequest::new(&meshes, perspective_view(), limits());
    let error = registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect_err("unsupported");
    assert_eq!(calls_of(&calls), 0);
    assert!(
        matches!(error.error, GeomError::Unsupported { .. }),
        "{error}"
    );
    let ortho = VisibleFaceRequest::new(&meshes, view(), limits());
    let error = registry(&calls, well_behaved)
        .visible_faces(&ortho, &options())
        .expect_err("the audit refuses the mesh");
    assert_eq!(error.mesh, Some(MeshId(1)));
}

#[test]
fn perspective_conformance_registration_follows_the_claim() {
    let calls = counter();
    let mut registry = MeshVisibleFacesRegistry::new();
    let unclaimed = registry
        .register_perspective_conformant(0, Recorder::new("declines", &calls, well_behaved));
    assert!(unclaimed.is_err());
    let claimed = registry.register_perspective_conformant(
        0,
        Recorder::new("claims", &calls, well_behaved).claiming_perspective(),
    );
    assert!(claimed.expect_err("unverifiable").claimed);
    assert_eq!(registry.providers().count(), 0);
}

fn declines_perspective(
    request: &VisibleFaceRequest<'_>,
    _outcome: &mut VisibleFaceOutcome,
) -> Result<(), HiddenLineError> {
    if request.view.projection == ViewProjection::Perspective {
        return Err(GeomError::Unsupported {
            backend: BackendId::new("declines"),
            operation: Operation::MeshVisibleFaces,
        }
        .into());
    }
    Ok(())
}

#[test]
fn perspective_registration_also_runs_the_base_suite() {
    use axiolid_mesh_hidden_lines_contract::conformance::perspective::PerspectiveConformanceFailure;

    let calls = counter();
    let mut registry = MeshVisibleFacesRegistry::new();
    let report = registry
        .register_perspective_conformant(0, Recorder::new("declines", &calls, declines_perspective))
        .expect_err("the canned regions break the base contract");
    assert!(!report.claimed);
    assert!(report
        .failures
        .iter()
        .any(|f| matches!(f, PerspectiveConformanceFailure::VisibleFacesBase(_))));
    assert!(report
        .failures
        .iter()
        .all(|f| matches!(f, PerspectiveConformanceFailure::VisibleFacesBase(_))));
    assert_eq!(registry.providers().count(), 0);
}

/// A registry whose only provider claims see-through triangles.
fn see_through_registry(
    calls: &Arc<AtomicUsize>,
    behaviour: Behaviour,
) -> MeshVisibleFacesRegistry {
    let mut registry = MeshVisibleFacesRegistry::new();
    registry.register(
        0,
        Recorder::new("claims", calls, behaviour).claiming_see_through(),
    );
    registry
}

#[test]
fn a_see_through_triangle_reaches_only_providers_that_claim_it() {
    let (declining, claiming) = (counter(), counter());
    let cube = cube();
    let flags = [Occlusion::SeeThrough; 12];
    let meshes = [
        mesh(1, MeshRole::Solid, &cube),
        mesh(2, MeshRole::Solid, &cube).with_occlusion(&flags),
    ];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let mut registry = MeshVisibleFacesRegistry::new();
    registry.register(1, Recorder::new("declines", &declining, well_behaved));
    registry.register(
        0,
        Recorder::new("claims", &claiming, well_behaved).claiming_see_through(),
    );
    registry
        .visible_faces(&request, &options())
        .expect("the claiming provider answers");
    assert_eq!((calls_of(&declining), calls_of(&claiming)), (0, 1));

    // Flags that are all opaque are no see-through request.
    let opaque = [Occlusion::Opaque; 12];
    let meshes = [mesh(1, MeshRole::Solid, &cube).with_occlusion(&opaque)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    registry
        .visible_faces(&request, &options())
        .expect("answers");
    assert_eq!((calls_of(&declining), calls_of(&claiming)), (1, 1));
}

#[test]
fn a_see_through_triangle_without_a_claiming_provider_is_unsupported() {
    let calls = counter();
    let cube = cube();
    let flags = [Occlusion::SeeThrough; 12];
    let meshes = [mesh(1, MeshRole::Solid, &cube).with_occlusion(&flags)];
    let request = VisibleFaceRequest::new(&meshes, view(), limits());
    let error = registry(&calls, well_behaved)
        .visible_faces(&request, &options())
        .expect_err("no provider honours see-through triangles");
    assert_eq!(calls_of(&calls), 0);
    assert_eq!(error.mesh, None);
    assert!(matches!(
        error.error,
        GeomError::Unsupported {
            operation: Operation::MeshVisibleFaces,
            ..
        }
    ));
}

#[test]
fn an_occlusion_slice_of_the_wrong_length_names_the_mesh() {
    let (good, cube) = (cube(), cube());
    let short = [Occlusion::SeeThrough; 4];
    let opaque_long = [Occlusion::Opaque; 13];
    for flags in [&short[..], &opaque_long[..]] {
        let meshes = [
            mesh(1, MeshRole::Solid, &good),
            mesh(2, MeshRole::Solid, &cube).with_occlusion(flags),
        ];
        let request = VisibleFaceRequest::new(&meshes, view(), limits());
        let calls = counter();
        let error = see_through_registry(&calls, well_behaved)
            .visible_faces(&request, &options())
            .expect_err("a wrong occlusion slice");
        assert_eq!(error.mesh, Some(MeshId(2)));
        assert!(
            matches!(&error.error, GeomError::InvalidInput(detail) if detail.contains("occlusion")),
            "{error:?}"
        );
        assert_eq!(calls_of(&calls), 0);
    }
}

#[test]
fn a_face_group_that_mixes_occlusion_names_the_mesh() {
    let (good, cube) = (cube(), cube());
    let mut flags = [Occlusion::Opaque; 12];
    flags[4] = Occlusion::SeeThrough;
    // The top face (triangles 2 and 3) in a group of its own.
    let mut groups = [FaceGroup(0); 12];
    groups[2] = FaceGroup(1);
    groups[3] = FaceGroup(1);
    let mut top = [Occlusion::Opaque; 12];
    top[2] = Occlusion::SeeThrough;
    top[3] = Occlusion::SeeThrough;
    let mixed_single = mesh(2, MeshRole::Solid, &cube).with_occlusion(&flags);
    let mixed_groups = mesh(2, MeshRole::Solid, &cube)
        .with_groups(&groups)
        .with_occlusion(&flags);
    for mixed in [mixed_single, mixed_groups] {
        let meshes = [mesh(1, MeshRole::Solid, &good), mixed];
        let request = VisibleFaceRequest::new(&meshes, view(), limits());
        let calls = counter();
        let error = see_through_registry(&calls, well_behaved)
            .visible_faces(&request, &options())
            .expect_err("a mixed face group");
        assert_eq!(error.mesh, Some(MeshId(2)));
        assert!(
            matches!(&error.error, GeomError::InvalidInput(detail) if detail.contains("mixes")),
            "{error:?}"
        );
        assert_eq!(calls_of(&calls), 0);
    }

    // A see-through group of its own is accepted, and planar groups are the
    // provider's to split.
    let calls = counter();
    for uniform in [
        mesh(2, MeshRole::Solid, &cube)
            .with_groups(&groups)
            .with_occlusion(&top),
        mesh(2, MeshRole::Solid, &cube)
            .with_planar_groups()
            .with_occlusion(&flags),
    ] {
        let meshes = [mesh(1, MeshRole::Solid, &good), uniform];
        let request = VisibleFaceRequest::new(&meshes, view(), limits());
        see_through_registry(&calls, well_behaved)
            .visible_faces(&request, &options())
            .expect("uniform groups");
    }
    assert_eq!(calls_of(&calls), 2);
}

const G0: FaceSource = FaceSource::Group(FaceGroup(0));

fn behind(mesh: u64, source: FaceSource, layer: u32) -> BehindSurface {
    BehindSurface::new(MeshId(mesh), source, layer)
}

/// Replace the regions by one region of the see-through mesh 2 with
/// `surfaces` behind it, and count them in the evidence.
fn pane_region(outcome: &mut VisibleFaceOutcome, surfaces: Vec<BehindSurface>) {
    outcome.evidence.behind_surfaces = surfaces.len();
    outcome.regions = vec![region(2, G0, triangle_ring(), &[2, 3]).with_behind(surfaces)];
}

/// Mesh 1 opaque, mesh 2 wholly see-through, both solid cubes, under a near
/// plane when `near` is set; answered by a provider that claims see-through
/// triangles and shapes the outcome by `behaviour`.
fn see_through_faces(
    near: bool,
    max_behind: usize,
    behaviour: Behaviour,
) -> Result<VisibleFaceOutcome, HiddenLineError> {
    let cube = cube();
    let flags = [Occlusion::SeeThrough; 12];
    let meshes = [
        mesh(1, MeshRole::Solid, &cube),
        mesh(2, MeshRole::Solid, &cube).with_occlusion(&flags),
    ];
    let view = if near { view().with_near(0.5) } else { view() };
    let limits = limits().with_max_behind_surfaces(max_behind);
    let request = VisibleFaceRequest::new(&meshes, view, limits);
    let calls = counter();
    let result = see_through_registry(&calls, behaviour).visible_faces(&request, &options());
    assert_eq!(calls_of(&calls), 1);
    result
}

#[test]
fn valid_stacks_behind_a_see_through_region_are_accepted() {
    let cases: Vec<(&str, bool, Behaviour)> = vec![
        ("an opaque surface behind", false, |_, o| {
            pane_region(o, vec![behind(1, G0, 1)]);
            Ok(())
        }),
        ("a tie, then a cap", true, |_, o| {
            pane_region(o, vec![behind(1, G0, 0), behind(1, FaceSource::Cap, 0)]);
            Ok(())
        }),
        ("glass behind glass, then a cap", true, |_, o| {
            pane_region(o, vec![behind(2, G0, 1), behind(1, FaceSource::Cap, 2)]);
            Ok(())
        }),
        (
            "two regions ordered by what is behind them",
            false,
            |_, o| {
                pane_region(o, vec![behind(1, G0, 1)]);
                let second = region(2, G0, triangle_ring(), &[2, 3])
                    .with_behind(vec![behind(2, G0, 1), behind(1, G0, 2)]);
                o.regions.push(second);
                o.evidence.regions = 2;
                o.evidence.behind_surfaces = 3;
                Ok(())
            },
        ),
    ];
    for (name, near, behaviour) in cases {
        see_through_faces(near, 16, behaviour).unwrap_or_else(|error| panic!("{name}: {error:?}"));
    }
    // An empty stack needs no budget.
    see_through_faces(false, 0, |_, o| {
        pane_region(o, Vec::new());
        Ok(())
    })
    .expect("nothing behind");
}

#[test]
fn stack_violations_are_caught() {
    let cases: Vec<(&str, bool, Behaviour)> = vec![
        ("a foreign mesh behind", false, |_, o| {
            pane_region(o, vec![behind(9, G0, 1)]);
            Ok(())
        }),
        ("a cap behind without a near plane", false, |_, o| {
            pane_region(o, vec![behind(1, FaceSource::Cap, 1)]);
            Ok(())
        }),
        ("a cap before the last layer", true, |_, o| {
            pane_region(o, vec![behind(1, FaceSource::Cap, 1), behind(2, G0, 2)]);
            Ok(())
        }),
        ("an opaque mesh before the last layer", false, |_, o| {
            pane_region(o, vec![behind(1, G0, 1), behind(2, G0, 2)]);
            Ok(())
        }),
        ("a first layer past 1", false, |_, o| {
            pane_region(o, vec![behind(1, G0, 2)]);
            Ok(())
        }),
        ("a skipped layer", false, |_, o| {
            pane_region(o, vec![behind(2, G0, 1), behind(1, G0, 3)]);
            Ok(())
        }),
        ("a layer out of order", false, |_, o| {
            pane_region(o, vec![behind(2, G0, 2), behind(1, G0, 1)]);
            Ok(())
        }),
        ("meshes out of order in a layer", false, |_, o| {
            pane_region(o, vec![behind(2, G0, 1), behind(1, G0, 1)]);
            Ok(())
        }),
        ("a repeated entry", false, |_, o| {
            pane_region(o, vec![behind(1, G0, 1), behind(1, G0, 1)]);
            Ok(())
        }),
        ("its own source in layer 0", false, |_, o| {
            pane_region(o, vec![behind(2, G0, 0)]);
            Ok(())
        }),
        ("an opaque region with a stack", false, |_, o| {
            o.regions[0] = o.regions[0].clone().with_behind(vec![behind(2, G0, 1)]);
            o.evidence.behind_surfaces = 1;
            Ok(())
        }),
        ("behind surfaces miscounted", false, |_, o| {
            pane_region(o, vec![behind(1, G0, 1)]);
            o.evidence.behind_surfaces = 0;
            Ok(())
        }),
        ("see-through triangles miscounted", false, |_, o| {
            o.evidence.see_through_triangles = 0;
            Ok(())
        }),
        (
            "regions out of order by what is behind them",
            false,
            |_, o| {
                pane_region(o, vec![behind(2, G0, 1), behind(1, G0, 2)]);
                let second =
                    region(2, G0, triangle_ring(), &[2, 3]).with_behind(vec![behind(1, G0, 1)]);
                o.regions.push(second);
                o.evidence.regions = 2;
                o.evidence.behind_surfaces = 3;
                Ok(())
            },
        ),
    ];
    for (name, near, behaviour) in cases {
        let error = see_through_faces(near, 16, behaviour).expect_err(name);
        assert!(
            matches!(error.error, GeomError::BackendContractViolation { .. }),
            "{name}: {error:?}"
        );
        assert_eq!(error.mesh, None, "{name}");
    }
}

#[test]
fn behind_surfaces_over_the_limit_exceed_the_budget() {
    let two = |_: &VisibleFaceRequest<'_>, o: &mut VisibleFaceOutcome| {
        pane_region(o, vec![behind(2, G0, 1), behind(1, G0, 2)]);
        Ok(())
    };
    see_through_faces(false, 2, two).expect("at the limit");
    for max in [0, 1] {
        let error = see_through_faces(false, max, two).expect_err("over the limit");
        assert_eq!(
            error.error,
            GeomError::BudgetExceeded {
                resource: "visible-face behind surfaces"
            },
            "{max}"
        );
    }
}
