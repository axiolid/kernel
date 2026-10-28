//! Backend-neutral mesh plane-section contract tests.

#![cfg(feature = "mesh-section")]

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, ExecutionOptions,
    ExecutionTarget, GeomError, GeomResult, ScratchRequirement,
};
use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_dispatch::MeshPlaneSectionRegistry;
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::{
    MeshPlaneSection, SectionContour, SectionEvidence, SectionLimits, SectionOutcome,
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

fn frame() -> Frame3 {
    Frame3 {
        origin: Vec3::new(0.0, 0.0, 0.5),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn limits() -> SectionLimits {
    SectionLimits::new(1_000, 1_000, 1_000, 100)
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

#[derive(Debug)]
struct Recorder {
    calls: Arc<AtomicUsize>,
    scratch: ScratchRequirement,
}

impl Backend for Recorder {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("section-recorder"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl MeshPlaneSection for Recorder {
    fn scratch_requirement(&self) -> ScratchRequirement {
        self.scratch
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn section(
        &self,
        mesh: &TriMesh,
        frame: Frame3,
        _limits: SectionLimits,
        options: &ExecutionOptions,
    ) -> GeomResult<SectionOutcome> {
        options.check_cancelled()?;
        self.calls.fetch_add(1, Ordering::SeqCst);
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ];
        Ok(SectionOutcome::new(
            frame,
            vec![SectionContour::new(points)],
            SectionEvidence::input_mesh(mesh.triangle_count(), 4, 1),
        ))
    }
}

fn registry(calls: Arc<AtomicUsize>, scratch: ScratchRequirement) -> MeshPlaneSectionRegistry {
    let mut registry = MeshPlaneSectionRegistry::new();
    registry.register(0, Recorder { calls, scratch });
    registry
}

#[test]
fn dispatch_returns_plane_local_closed_contours_with_mesh_provenance() {
    let calls = Arc::new(AtomicUsize::new(0));
    let outcome = registry(Arc::clone(&calls), ScratchRequirement::None)
        .section(&cube(), frame(), limits(), &options())
        .expect("valid section contract");

    assert_eq!(calls.load(Ordering::SeqCst), 1);
    assert_eq!(outcome.frame, frame());
    assert_eq!(outcome.contours.len(), 1);
    assert_eq!(outcome.contours[0].points.len(), 4);
    assert!(outcome.contours[0].is_closed());
    assert_eq!(outcome.evidence.source_triangles, 12);
    assert!(outcome.evidence.is_derived_from_input_mesh());
}

#[test]
fn registry_rejects_open_or_dirty_solids_before_dispatch() {
    let calls = Arc::new(AtomicUsize::new(0));
    let registry = registry(Arc::clone(&calls), ScratchRequirement::None);
    let open_sheet = TriMesh::new(vec![Vec3::ZERO, Vec3::X, Vec3::Y], vec![0, 1, 2]);

    assert!(matches!(
        registry.section(&open_sheet, frame(), limits(), &options()),
        Err(GeomError::NotManifold(_))
    ));
    let empty = TriMesh::new(Vec::new(), Vec::new());
    assert!(matches!(
        registry.section(&empty, frame(), limits(), &options()),
        Err(GeomError::InvalidInput(_))
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

#[test]
fn registry_rejects_a_partial_triangle_before_dispatch() {
    let calls = Arc::new(AtomicUsize::new(0));
    let registry = registry(Arc::clone(&calls), ScratchRequirement::None);
    // Two indices are no triangle: the mesh must not pass as an empty solid.
    let two = TriMesh::new(vec![Vec3::ZERO, Vec3::X, Vec3::Y], vec![0, 1]);
    assert!(matches!(
        registry.section(&two, frame(), limits(), &options()),
        Err(GeomError::InvalidInput(_))
    ));
    // A cube with a stray trailing index is refused too.
    let mut cube = cube();
    cube.indices.push(0);
    assert!(matches!(
        registry.section(&cube, frame(), limits(), &options()),
        Err(GeomError::InvalidInput(_))
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

#[test]
fn registry_rejects_non_finite_signed_volume_before_dispatch() {
    let calls = Arc::new(AtomicUsize::new(0));
    let registry = registry(Arc::clone(&calls), ScratchRequirement::None);
    let magnitude = 1.0e308;
    let extreme = TriMesh::new(
        vec![
            Point3::ZERO,
            Point3::new(magnitude, 0.0, 0.0),
            Point3::new(0.0, magnitude, 0.0),
            Point3::new(0.0, 0.0, magnitude),
        ],
        vec![
            0, 2, 1, // base, outward -z
            0, 1, 3, // side, outward -y
            0, 3, 2, // side, outward -x
            1, 2, 3, // sloped face, outward +xyz
        ],
    );

    let error = registry
        .section(&extreme, frame(), limits(), &options())
        .unwrap_err();
    assert!(
        matches!(error, GeomError::InvalidInput(ref detail) if detail.contains("non-finite signed volume")),
        "extreme finite source must fail closed, got {error:?}"
    );
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

#[test]
fn memory_budget_is_admitted_before_allocating_topology_audit() {
    let calls = Arc::new(AtomicUsize::new(0));
    let registry = registry(Arc::clone(&calls), ScratchRequirement::None);
    let open_sheet = TriMesh::new(vec![Vec3::ZERO, Vec3::X, Vec3::Y], vec![0, 1, 2]);
    let constrained = options().with_memory_budget(0);

    assert!(matches!(
        registry.section(&open_sheet, frame(), limits(), &constrained),
        Err(GeomError::BudgetExceeded { resource: "memory" })
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

#[test]
fn exact_topology_audit_bound_admits_a_valid_source() {
    let calls = Arc::new(AtomicUsize::new(0));
    let registry = registry(Arc::clone(&calls), ScratchRequirement::None);
    let audit_bytes =
        axiolid_mesh::audit_mesh_scratch_bytes(cube().triangle_count()).expect("small audit bound");
    let constrained = options().with_memory_budget(audit_bytes);

    let outcome = registry
        .section(&cube(), frame(), limits(), &constrained)
        .expect("exact audit bound is sufficient");
    assert_eq!(outcome.contours.len(), 1);
    assert_eq!(calls.load(Ordering::SeqCst), 1);
}

#[derive(Debug)]
struct VertexBudgetRecorder {
    calls: Arc<AtomicUsize>,
}

impl Backend for VertexBudgetRecorder {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("vertex-budget"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl MeshPlaneSection for VertexBudgetRecorder {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::PerElement {
            bytes_per_element: 10,
        }
    }

    fn section(
        &self,
        _mesh: &TriMesh,
        _frame: Frame3,
        _limits: SectionLimits,
        _options: &ExecutionOptions,
    ) -> GeomResult<SectionOutcome> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        Err(GeomError::InvalidInput("must not dispatch".into()))
    }
}

#[test]
fn scratch_preflight_accounts_for_unused_source_positions() {
    let calls = Arc::new(AtomicUsize::new(0));
    let mut registry = MeshPlaneSectionRegistry::new();
    registry.register(
        0,
        VertexBudgetRecorder {
            calls: Arc::clone(&calls),
        },
    );
    let mut mesh = cube();
    mesh.positions
        .extend(std::iter::repeat_n(Point3::splat(9.0), 100));

    let error = registry
        .section(
            &mesh,
            frame(),
            SectionLimits::new(200, 20, 20, 4),
            &options().with_memory_budget(900),
        )
        .unwrap_err();

    assert!(matches!(
        error,
        GeomError::BudgetExceeded { resource: "memory" }
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

#[test]
fn registry_validates_the_section_frame_before_dispatch() {
    let calls = Arc::new(AtomicUsize::new(0));
    let registry = registry(Arc::clone(&calls), ScratchRequirement::None);
    let dirty = Frame3 {
        z: Vec3::new(0.0, 0.0, 2.0),
        ..frame()
    };

    assert!(matches!(
        registry.section(&cube(), dirty, limits(), &options()),
        Err(GeomError::InvalidInput(_))
    ));

    let zero_axis = Frame3 {
        z: Vec3::ZERO,
        ..frame()
    };
    let coarse = ExecutionOptions::new(Tolerance::new(1e-6, 2.0).unwrap());
    assert!(matches!(
        registry.section(&cube(), zero_axis, limits(), &coarse),
        Err(GeomError::InvalidInput(_))
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}

#[test]
fn source_and_scratch_limits_block_dispatch() {
    let calls = Arc::new(AtomicUsize::new(0));
    let registry = registry(
        Arc::clone(&calls),
        ScratchRequirement::Fixed { bytes: 4_096 },
    );
    assert!(matches!(
        registry.section(
            &cube(),
            frame(),
            limits(),
            &options().with_memory_budget(16)
        ),
        Err(GeomError::BudgetExceeded { resource: "memory" })
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);

    let tiny = SectionLimits::new(7, 12, 1_000, 100);
    assert!(matches!(
        registry.section(&cube(), frame(), tiny, &options()),
        Err(GeomError::BudgetExceeded {
            resource: "section source vertices"
        })
    ));
    assert_eq!(calls.load(Ordering::SeqCst), 0);
}
