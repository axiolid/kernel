//! Backend-neutral mesh hidden-line contract tests.
//!
//! Every pre-dispatch refusal must happen without calling the provider and
//! name the mesh at fault; every way a provider can break the outcome contract
//! must be caught; only unsupported or unavailable providers permit fallback.

#![cfg(feature = "mesh-hidden-lines")]

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, CancellationToken,
    ExecutionOptions, ExecutionTarget, GeomError, Operation, ScratchRequirement,
};
use axiolid_core::{Aabb2, Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_dispatch::MeshHiddenLinesRegistry;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeKind, EdgeSelection, HiddenLineError, HiddenLineEvidence, HiddenLineLimits, HiddenLineMesh,
    HiddenLineOutcome, HiddenLineRequest, HiddenLineSegment, MeshHiddenLines, MeshId, MeshRole,
    OrthographicView, SourceEdge, Visibility,
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

fn limits() -> HiddenLineLimits {
    HiddenLineLimits::new(8, 1_000, 1_000, 1_000_000, 1_000)
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn selection() -> EdgeSelection {
    EdgeSelection::new(Some(0.5))
}

fn view() -> OrthographicView {
    OrthographicView::new(frame())
}

fn mesh<'a>(id: u64, role: MeshRole, mesh: &'a TriMesh) -> HiddenLineMesh<'a> {
    HiddenLineMesh::new(MeshId(id), role, mesh)
}

/// Whether a refusal is the expected kind.
type Accepts = fn(&GeomError) -> bool;

/// What a recorder does after it has been called.
type Behaviour = fn(&HiddenLineRequest<'_>, &mut HiddenLineOutcome) -> Result<(), HiddenLineError>;

fn well_behaved(
    _request: &HiddenLineRequest<'_>,
    _outcome: &mut HiddenLineOutcome,
) -> Result<(), HiddenLineError> {
    Ok(())
}

#[derive(Debug)]
struct Recorder {
    id: &'static str,
    calls: Arc<AtomicUsize>,
    scratch: ScratchRequirement,
    behaviour: Behaviour,
}

impl Recorder {
    fn new(id: &'static str, calls: &Arc<AtomicUsize>, behaviour: Behaviour) -> Self {
        Self {
            id,
            calls: Arc::clone(calls),
            scratch: ScratchRequirement::None,
            behaviour,
        }
    }
}

impl Backend for Recorder {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new(self.id), ExecutionTarget::PortableCpu)
    }
}

/// One valid visible silhouette on the first edge of the first mesh, shaped
/// by `behaviour`.
fn base_outcome(request: &HiddenLineRequest<'_>) -> HiddenLineOutcome {
    let triangles = request.meshes.iter().map(|m| m.mesh.triangle_count()).sum();
    let mut evidence = HiddenLineEvidence::input_mesh(request.meshes.len(), triangles);
    evidence.visible_segments = 1;
    HiddenLineOutcome::new(
        request.view,
        vec![segment(request.meshes[0].id, [0, 1])],
        evidence,
    )
}

fn segment(mesh: MeshId, vertices: [u32; 2]) -> HiddenLineSegment {
    HiddenLineSegment::new(
        SourceEdge { mesh, vertices },
        EdgeKind::Silhouette,
        Visibility::Visible,
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
    )
}

impl MeshHiddenLines for Recorder {
    fn scratch_requirement(&self) -> ScratchRequirement {
        self.scratch
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        options.check_cancelled()?;
        self.calls.fetch_add(1, Ordering::SeqCst);
        let mut outcome = base_outcome(request);
        (self.behaviour)(request, &mut outcome)?;
        Ok(outcome)
    }
}

fn registry(calls: &Arc<AtomicUsize>, behaviour: Behaviour) -> MeshHiddenLinesRegistry {
    let mut registry = MeshHiddenLinesRegistry::new();
    registry.register(0, Recorder::new("hidden-lines-recorder", calls, behaviour));
    registry
}

fn counter() -> Arc<AtomicUsize> {
    Arc::new(AtomicUsize::new(0))
}

fn calls_of(calls: &Arc<AtomicUsize>) -> usize {
    calls.load(Ordering::SeqCst)
}

#[test]
fn dispatch_returns_segments_with_mesh_provenance() {
    let calls = counter();
    let cube = cube();
    let meshes = [mesh(7, MeshRole::Solid, &cube)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let outcome = registry(&calls, well_behaved)
        .hidden_lines(&request, &options())
        .expect("valid hidden-line contract");

    assert_eq!(calls_of(&calls), 1);
    assert_eq!(outcome.view, view());
    assert_eq!(outcome.segments.len(), 1);
    assert_eq!(outcome.segments[0].source.mesh, MeshId(7));
    assert_eq!(outcome.evidence.source_meshes, 1);
    assert_eq!(outcome.evidence.source_triangles, 12);
    assert!(outcome.evidence.is_derived_from_input_mesh());
}

/// A good first mesh (id 10) and a second mesh (id 20) that may be at fault.
fn two_mesh_refusal(offender: &TriMesh, role: MeshRole) -> HiddenLineError {
    let calls = counter();
    let good = cube();
    let meshes = [mesh(10, MeshRole::Solid, &good), mesh(20, role, offender)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let error = registry(&calls, well_behaved)
        .hidden_lines(&request, &options())
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
            "solid with two indices",
            TriMesh::new(vec![origin, Vec3::X, Vec3::Y], vec![0, 1]),
            MeshRole::Solid,
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
            "three triangles on an edge",
            TriMesh::new(
                vec![origin, Vec3::X, Vec3::Y, Vec3::new(0.0, -1.0, 0.0), Vec3::Z],
                vec![0, 1, 2, 1, 0, 3, 0, 1, 4],
            ),
            MeshRole::Sheet,
            |e| matches!(e, GeomError::NotManifold(_)),
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
        let error = two_mesh_refusal(&offender, role);
        assert!(accepts(&error.error), "{name}: {error:?}");
        assert_eq!(error.mesh, Some(MeshId(20)), "{name}");
        assert!(
            error.to_string().starts_with("mesh 20: "),
            "{name}: {error}"
        );
        assert!(!error.error.to_string().contains("mesh 20"), "{name}");
    }
}

#[test]
fn an_open_sheet_is_accepted() {
    let calls = counter();
    let good = cube();
    let sheet = quad();
    let meshes = [
        mesh(10, MeshRole::Solid, &good),
        mesh(20, MeshRole::Sheet, &sheet),
    ];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    registry(&calls, well_behaved)
        .hidden_lines(&request, &options())
        .expect("a sheet may be open");
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
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let error = registry(&calls, well_behaved)
        .hidden_lines(&request, &options())
        .unwrap_err();
    assert!(matches!(error.error, GeomError::InvalidInput(_)));
    assert_eq!(error.mesh, Some(MeshId(20)));
    assert_eq!(calls_of(&calls), 0);
}

/// Requests refused as a whole: the refusal names no mesh.
#[test]
fn each_whole_request_refusal_names_no_mesh() {
    type Edit = fn(&mut HiddenLineRequest<'_>);
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
                        resource: "hidden-line meshes"
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
                        resource: "hidden-line source vertices"
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
                        resource: "hidden-line source triangles"
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
            "crease cosine above one",
            |r| r.edges = EdgeSelection::new(Some(1.5)),
            |e| matches!(e, GeomError::InvalidInput(_)),
        ),
        (
            "crease cosine NaN",
            |r| r.edges = EdgeSelection::new(Some(f64::NAN)),
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
        let mut request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
        edit(&mut request);
        let error = registry(&calls, well_behaved)
            .hidden_lines(&request, &options())
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
    // Exactly at every limit, cosine at both ends, an empty-area region.
    let tight = HiddenLineLimits::new(2, 16, 24, 10, 1);
    for cosine in [-1.0, 1.0] {
        let request =
            HiddenLineRequest::new(&meshes, view(), EdgeSelection::new(Some(cosine)), tight)
                .with_region(Aabb2 {
                    min: Point2::new(0.5, 0.5),
                    max: Point2::new(0.5, 0.5),
                });
        registry(&calls, well_behaved)
            .hidden_lines(&request, &options())
            .expect("limits are inclusive");
    }
    assert_eq!(calls_of(&calls), 2);
}

#[test]
fn cancellation_refuses_the_whole_request() {
    let calls = counter();
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let token = CancellationToken::new();
    token.cancel();
    let error = registry(&calls, well_behaved)
        .hidden_lines(&request, &options().with_cancellation(token))
        .unwrap_err();
    assert_eq!(error, HiddenLineError::request(GeomError::Cancelled));
    assert_eq!(calls_of(&calls), 0);
}

// --- outcome validation -------------------------------------------------

fn violation_of(request: &HiddenLineRequest<'_>, behaviour: Behaviour) -> HiddenLineError {
    let calls = counter();
    let error = registry(&calls, behaviour)
        .hidden_lines(request, &options())
        .expect_err("the outcome must be refused");
    assert_eq!(calls_of(&calls), 1);
    error
}

fn assert_violation(request: &HiddenLineRequest<'_>, behaviour: Behaviour, name: &str) {
    let error = violation_of(request, behaviour);
    assert!(
        matches!(error.error, GeomError::BackendContractViolation { .. }),
        "{name}: {error:?}"
    );
    assert_eq!(error.mesh, None, "{name}");
}

#[test]
fn outcome_violations_are_caught() {
    let (a, b) = (cube(), cube());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Solid, &b)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());

    let cases: Vec<(&str, Behaviour)> = vec![
        ("view changed", |_, o| {
            o.view = o.view.with_near(5.0);
            Ok(())
        }),
        ("unknown mesh", |_, o| {
            o.segments[0].source.mesh = MeshId(99);
            Ok(())
        }),
        ("descending vertices", |_, o| {
            o.segments[0].source.vertices = [1, 0];
            Ok(())
        }),
        ("repeated vertex", |_, o| {
            o.segments[0].source.vertices = [1, 1];
            Ok(())
        }),
        ("vertex out of range", |_, o| {
            o.segments[0].source.vertices = [0, 8];
            Ok(())
        }),
        ("non-finite start", |_, o| {
            o.segments[0].start = Point2::new(f64::NAN, 0.0);
            Ok(())
        }),
        ("non-finite end", |_, o| {
            o.segments[0].end = Point2::new(0.0, f64::INFINITY);
            Ok(())
        }),
        ("zero length", |_, o| {
            o.segments[0].end = o.segments[0].start;
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
        ("visible count", |_, o| {
            o.evidence.visible_segments = 0;
            Ok(())
        }),
        ("hidden count", |_, o| {
            o.evidence.hidden_segments = 1;
            Ok(())
        }),
        ("segments out of order within a mesh", |_, o| {
            o.segments.insert(0, segment(MeshId(1), [0, 2]));
            o.evidence.visible_segments = 2;
            Ok(())
        }),
        ("one edge with two kinds", |_, o| {
            let mut second = segment(MeshId(1), [0, 1]);
            second.kind = EdgeKind::Crease;
            o.segments.push(second);
            o.evidence.visible_segments = 2;
            Ok(())
        }),
        ("segments out of order across meshes", |_, o| {
            o.segments.insert(0, segment(MeshId(2), [0, 1]));
            o.evidence.visible_segments = 2;
            Ok(())
        }),
    ];
    for (name, behaviour) in cases {
        // "out of order": the base segment (mesh 1, [0,1]) comes second.
        assert_violation(&request, behaviour, name);
    }
}

#[test]
fn ordered_segments_are_accepted() {
    let calls = counter();
    let (a, b) = (cube(), cube());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Solid, &b)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let outcome = registry(&calls, |_, o| {
        // Non-decreasing: equal keys (an edge split into pieces), then a
        // larger edge, then the next mesh.
        o.segments.push(segment(MeshId(1), [0, 1]));
        o.segments.push(segment(MeshId(1), [0, 2]));
        o.segments.push(segment(MeshId(2), [0, 1]));
        o.evidence.visible_segments = 4;
        Ok(())
    })
    .hidden_lines(&request, &options())
    .expect("ordered output");
    assert_eq!(outcome.segments.len(), 4);
}

#[test]
fn unselected_kinds_are_caught() {
    let (a, b) = (cube(), cube());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Solid, &b)];
    let kinds: [(&str, EdgeSelection, EdgeKind); 3] = [
        (
            "silhouette",
            selection().with_silhouette(false),
            EdgeKind::Silhouette,
        ),
        (
            "boundary",
            selection().with_boundary(false),
            EdgeKind::Boundary,
        ),
        ("crease", EdgeSelection::new(None), EdgeKind::Crease),
    ];
    for (name, edges, kind) in kinds {
        let request = HiddenLineRequest::new(&meshes, view(), edges, limits());
        let calls = counter();
        let behaviour: Behaviour = match kind {
            EdgeKind::Silhouette => |_, _| Ok(()),
            EdgeKind::Boundary => |_, o| {
                o.segments[0].kind = EdgeKind::Boundary;
                Ok(())
            },
            _ => |_, o| {
                o.segments[0].kind = EdgeKind::Crease;
                Ok(())
            },
        };
        let error = registry(&calls, behaviour)
            .hidden_lines(&request, &options())
            .expect_err(name);
        assert!(
            matches!(error.error, GeomError::BackendContractViolation { .. }),
            "{name}: {error:?}"
        );
    }
}

#[test]
fn hidden_segments_are_refused_when_not_requested() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let wants = selection();
    let declines = selection().with_hidden(false);
    let hide: Behaviour = |_, o| {
        o.segments[0].visibility = Visibility::Hidden;
        o.evidence.visible_segments = 0;
        o.evidence.hidden_segments = 1;
        Ok(())
    };
    // Requested: the same outcome is fine.
    let calls = counter();
    registry(&calls, hide)
        .hidden_lines(
            &HiddenLineRequest::new(&meshes, view(), wants, limits()),
            &options(),
        )
        .expect("hidden pieces were requested");
    // Not requested: refused.
    assert_violation(
        &HiddenLineRequest::new(&meshes, view(), declines, limits()),
        hide,
        "hidden emitted when not requested",
    );
}

#[test]
fn too_many_segments_exceed_the_output_budget() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let mut tight = limits();
    tight.max_output_segments = 1;
    let request = HiddenLineRequest::new(&meshes, view(), selection(), tight);
    let calls = counter();
    // One segment fits, two do not.
    registry(&calls, well_behaved)
        .hidden_lines(&request, &options())
        .expect("at the limit");
    let error = violation_of(&request, |_, o| {
        o.segments.push(segment(MeshId(1), [0, 2]));
        o.evidence.visible_segments = 2;
        Ok(())
    });
    assert_eq!(
        error.error,
        GeomError::BudgetExceeded {
            resource: "hidden-line output segments"
        }
    );
    assert_eq!(error.mesh, None);
}

// --- fallback ----------------------------------------------------------

fn unsupported(
    _: &HiddenLineRequest<'_>,
    _: &mut HiddenLineOutcome,
) -> Result<(), HiddenLineError> {
    Err(GeomError::Unsupported {
        backend: BackendId::new("first"),
        operation: Operation::MeshHiddenLines,
    }
    .into())
}

fn unavailable(
    _: &HiddenLineRequest<'_>,
    _: &mut HiddenLineOutcome,
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
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());

    for retryable in [unsupported as Behaviour, unavailable] {
        let (first, second) = (counter(), counter());
        let mut registry = MeshHiddenLinesRegistry::new();
        registry.register(10, Recorder::new("first", &first, retryable));
        registry.register(0, Recorder::new("second", &second, well_behaved));
        registry
            .hidden_lines(&request, &options())
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
                GeomError::Degenerate("provider refusal".into()),
            ))
        }),
        ("budget", |_, _| {
            Err(GeomError::BudgetExceeded {
                resource: "hidden-line candidate pairs",
            }
            .into())
        }),
        ("unsupported input", |_, _| {
            Err(GeomError::UnsupportedInput {
                backend: BackendId::new("first"),
                operation: Operation::MeshHiddenLines,
                input: "near clip through a capped solid",
            }
            .into())
        }),
    ];
    for (name, refusal) in refusals {
        let (first, second) = (counter(), counter());
        let mut registry = MeshHiddenLinesRegistry::new();
        registry.register(10, Recorder::new("first", &first, refusal));
        registry.register(0, Recorder::new("second", &second, well_behaved));
        let error = registry.hidden_lines(&request, &options()).unwrap_err();
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
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let calls = counter();
    let error = registry(&calls, |_, _| {
        Err(HiddenLineError::mesh(
            MeshId(1),
            GeomError::Degenerate("triangle 3 meets the near plane".into()),
        ))
    })
    .hidden_lines(&request, &options())
    .unwrap_err();
    assert_eq!(error.mesh, Some(MeshId(1)));
}

#[test]
fn without_a_matching_provider_the_request_is_unsupported() {
    let cube = cube();
    let meshes = [mesh(1, MeshRole::Solid, &cube)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let error = MeshHiddenLinesRegistry::new()
        .hidden_lines(&request, &options())
        .unwrap_err();
    assert!(matches!(
        error.error,
        GeomError::Unsupported {
            operation: Operation::MeshHiddenLines,
            ..
        }
    ));
    assert_eq!(error.mesh, None);

    // Every provider retryable: the last retryable error is returned.
    let calls = counter();
    let error = registry(&calls, unsupported)
        .hidden_lines(&request, &options())
        .unwrap_err();
    assert!(matches!(error.error, GeomError::Unsupported { .. }));
    assert_eq!(calls_of(&calls), 1);
}

#[test]
fn register_conformant_rejects_a_provider_that_fails_the_suite() {
    let mut registry = MeshHiddenLinesRegistry::new();
    let calls = counter();
    let report = registry
        .register_conformant(0, Recorder::new("recorder", &calls, well_behaved))
        .expect_err("a recorder draws nothing right");
    assert!(!report.is_success());
    assert_eq!(registry.providers().count(), 0);
}

// --- budgets ------------------------------------------------------------

fn budgeted(scratch: ScratchRequirement) -> (MeshHiddenLinesRegistry, Arc<AtomicUsize>) {
    let calls = counter();
    let mut registry = MeshHiddenLinesRegistry::new();
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
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let error = registry
        .hidden_lines(&request, &options().with_memory_budget(0))
        .unwrap_err();
    assert_eq!(
        error.error,
        GeomError::BudgetExceeded { resource: "memory" }
    );
    assert_eq!(calls_of(&calls), 0);
}

#[test]
fn audit_scratch_peaks_at_the_largest_mesh_not_their_sum() {
    let (registry, calls) = budgeted(ScratchRequirement::None);
    let (a, b) = (cube(), cube());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Solid, &b)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let one = axiolid_mesh::audit_mesh_scratch_bytes(12).expect("small audit bound");
    registry
        .hidden_lines(&request, &options().with_memory_budget(one))
        .expect("sequential audits share one peak");
    assert_eq!(calls_of(&calls), 1);
    let error = registry
        .hidden_lines(&request, &options().with_memory_budget(one - 1))
        .unwrap_err();
    assert_eq!(
        error.error,
        GeomError::BudgetExceeded { resource: "memory" }
    );
}

#[test]
fn provider_scratch_counts_every_mesh_and_the_larger_of_vertices_and_triangles() {
    let (registry, calls) = budgeted(ScratchRequirement::PerElement {
        bytes_per_element: 1_000,
    });
    let (a, b) = (cube(), quad());
    let meshes = [mesh(1, MeshRole::Solid, &a), mesh(2, MeshRole::Sheet, &b)];
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    // 12 + 2 triangles, 8 + 4 vertices: 14 elements.
    registry
        .hidden_lines(&request, &options().with_memory_budget(14_000))
        .expect("14 elements fit");
    assert_eq!(calls_of(&calls), 1);
    let error = registry
        .hidden_lines(&request, &options().with_memory_budget(13_999))
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
    let request = HiddenLineRequest::new(&meshes, view(), selection(), limits());
    let calls = counter();
    let error = registry(&calls, |_, _| {
        Err(HiddenLineError::mesh(
            MeshId(99),
            GeomError::Degenerate("elsewhere".into()),
        ))
    })
    .hidden_lines(&request, &options())
    .unwrap_err();
    assert!(
        matches!(error.error, GeomError::BackendContractViolation { .. }),
        "{error:?}"
    );
    assert_eq!(error.mesh, None);
    assert_eq!(calls_of(&calls), 1);
}
