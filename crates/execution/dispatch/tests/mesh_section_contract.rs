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

// ---- section_request (ADR 0091) -------------------------------------------

mod request {
    use super::*;
    use axiolid_mesh_section_contract::{
        ChainEnd, InPlaneRegion, InPlaneSide, PlaneSection, PlaneSectionEvidence, RegionLimits,
        SectionChain, SectionRequest, SectionRole,
    };

    fn pts(points: &[(f64, f64)]) -> Vec<Point2> {
        points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
    }

    fn request_with(limits: SectionLimits, region_limits: RegionLimits) -> SectionRequest {
        SectionRequest::new(frame(), limits, region_limits)
    }

    fn request() -> SectionRequest {
        request_with(limits(), RegionLimits::new(100, 1_000, 1_000))
    }

    /// A valid answer: a closed chain, an open chain, and a region with two holes.
    fn valid(mesh: &TriMesh, request: &SectionRequest) -> PlaneSection {
        let chains = vec![
            SectionChain::closed(pts(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])),
            SectionChain::open(
                pts(&[(0.5, 0.0), (0.5, 0.5)]),
                ChainEnd::InPlaneRegion,
                ChainEnd::InPlaneRegion,
            ),
        ];
        let regions = vec![InPlaneRegion::new(
            InPlaneSide::Below,
            pts(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]),
            vec![
                pts(&[(0.1, 0.1), (0.1, 0.2), (0.2, 0.2)]),
                pts(&[(0.5, 0.5), (0.5, 0.6), (0.6, 0.6)]),
            ],
        )];
        let boundary_edges = match request.role {
            SectionRole::Sheet => {
                axiolid_mesh::try_audit_mesh(mesh, Tolerance::METRE)
                    .expect("audit")
                    .boundary_edges
            }
            _ => 0,
        };
        let evidence = PlaneSectionEvidence::input_mesh(mesh.triangle_count())
            .with_source_counts(mesh.triangle_count().min(2), 1, 3)
            .with_role(request.role, boundary_edges)
            .with_output(&chains, &regions);
        PlaneSection::new(request.frame, chains, regions, evidence)
    }

    type Mutation = fn(&mut PlaneSection);

    #[derive(Debug)]
    struct Scripted {
        mutate: Mutation,
        calls: Arc<AtomicUsize>,
    }

    impl Backend for Scripted {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("section-scripted"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl MeshPlaneSection for Scripted {
        fn section(
            &self,
            _mesh: &TriMesh,
            _frame: Frame3,
            _limits: SectionLimits,
            _options: &ExecutionOptions,
        ) -> GeomResult<SectionOutcome> {
            Err(GeomError::Degenerate("legacy path not under test".into()))
        }

        fn request_scratch_requirement(&self) -> ScratchRequirement {
            ScratchRequirement::None
        }

        fn section_request(
            &self,
            mesh: &TriMesh,
            request: &SectionRequest,
            _options: &ExecutionOptions,
        ) -> GeomResult<PlaneSection> {
            self.calls.fetch_add(1, Ordering::SeqCst);
            let mut section = valid(mesh, request);
            (self.mutate)(&mut section);
            Ok(section)
        }
    }

    fn run_on(
        mesh: &TriMesh,
        mutate: fn(&mut PlaneSection),
        request: &SectionRequest,
    ) -> GeomResult<PlaneSection> {
        let mut registry = MeshPlaneSectionRegistry::new();
        registry.register(
            0,
            Scripted {
                mutate,
                calls: Arc::new(AtomicUsize::new(0)),
            },
        );
        registry.section_request(mesh, request, &options())
    }

    fn run(mutate: fn(&mut PlaneSection), request: &SectionRequest) -> GeomResult<PlaneSection> {
        run_on(&cube(), mutate, request)
    }

    #[test]
    fn a_valid_answer_passes_unchanged() {
        let section = run(|_| {}, &request()).expect("valid");
        assert_eq!(section, valid(&cube(), &request()));
    }

    #[test]
    fn each_broken_invariant_is_a_contract_violation() {
        let cases: Vec<(&str, Mutation)> = vec![
            ("frame", |s| s.frame.origin = Vec3::new(0.0, 0.0, 9.0)),
            ("closed chain too short", |s| s.chains[0].points.truncate(2)),
            ("open chain too short", |s| s.chains[1].points.truncate(1)),
            ("chain non-finite", |s| {
                s.chains[0].points[2] = Point2::new(f64::NAN, 1.0)
            }),
            ("chain duplicate", |s| {
                s.chains[0].points[1] = s.chains[0].points[0]
            }),
            ("closed chain wrap duplicate", |s| {
                s.chains[0].points[3] = s.chains[0].points[0]
            }),
            ("open chain duplicate", |s| {
                s.chains[1].points[1] = s.chains[1].points[0]
            }),
            ("closed chain start", |s| s.chains[0].points.rotate_left(1)),
            ("open chain direction", |s| s.chains[1].points.reverse()),
            ("chain order", |s| s.chains.reverse()),
            ("ring too short", |s| s.regions[0].holes[0].truncate(2)),
            ("ring non-finite", |s| {
                s.regions[0].outer[1] = Point2::new(f64::INFINITY, 0.0)
            }),
            ("ring duplicate", |s| {
                s.regions[0].holes[0][1] = s.regions[0].holes[0][0]
            }),
            ("ring start", |s| s.regions[0].outer.rotate_left(1)),
            ("hole order", |s| s.regions[0].holes.reverse()),
            ("region order", |s| {
                let mut second = s.regions[0].clone();
                second.outer = pts(&[(-1.0, 0.0), (0.0, -1.0), (0.0, 0.5)]);
                second.holes.clear();
                s.regions.push(second);
            }),
            ("side order", |s| {
                let mut second = s.regions[0].clone();
                second.side = InPlaneSide::Above;
                s.regions.insert(0, second);
            }),
            ("evidence source", |s| s.evidence.source_triangles += 1),
            ("evidence chains", |s| s.evidence.open_chains += 1),
            ("evidence regions", |s| s.evidence.region_vertices += 1),
            ("evidence in-plane", |s| s.evidence.in_plane_triangles = 13),
            ("evidence role", |s| s.evidence.role = SectionRole::Sheet),
            ("evidence boundary edges", |s| s.evidence.boundary_edges = 1),
            ("solid chain ends on the mesh boundary", |s| {
                s.chains[1].ends = Some([ChainEnd::MeshBoundary, ChainEnd::InPlaneRegion])
            }),
        ];
        for (name, mutate) in cases {
            let error = run(mutate, &request()).expect_err(name);
            assert!(
                matches!(error, GeomError::BackendContractViolation { .. }),
                "{name}: {error:?}"
            );
        }
    }

    #[test]
    fn budgets_are_enforced_on_the_answer() {
        let regions = RegionLimits::new(100, 1_000, 1_000);
        let chain_points = 4 + 2;
        let ring_points = 4 + 3 + 3;
        let cases = [
            (
                request_with(SectionLimits::new(1_000, 1_000, 1_000, 1), regions),
                "section contours",
            ),
            (
                request_with(
                    SectionLimits::new(1_000, 1_000, chain_points - 1, 10),
                    regions,
                ),
                "section output vertices",
            ),
            (
                request_with(limits(), RegionLimits::new(0, 1_000, 1_000)),
                "section in-plane regions",
            ),
            (
                request_with(limits(), RegionLimits::new(10, ring_points - 1, 1_000)),
                "section in-plane region vertices",
            ),
            (
                request_with(limits(), RegionLimits::new(10, 1_000, 2)),
                "section in-plane pairs",
            ),
        ];
        for (request, resource) in cases {
            match run(|_| {}, &request) {
                Err(GeomError::BudgetExceeded { resource: got }) => assert_eq!(got, resource),
                other => panic!("{resource}: {other:?}"),
            }
        }
    }

    #[test]
    fn a_solid_request_refuses_an_open_mesh_before_dispatch() {
        let calls = Arc::new(AtomicUsize::new(0));
        let mut registry = MeshPlaneSectionRegistry::new();
        registry.register(
            0,
            Scripted {
                mutate: |_| {},
                calls: Arc::clone(&calls),
            },
        );
        let open_sheet = TriMesh::new(vec![Vec3::ZERO, Vec3::X, Vec3::Y], vec![0, 1, 2]);
        assert!(matches!(
            registry.section_request(&open_sheet, &request(), &options()),
            Err(GeomError::NotManifold(_))
        ));
        assert_eq!(calls.load(Ordering::SeqCst), 0);
    }

    /// A triangle: three vertices, three boundary edges.
    fn open_triangle() -> TriMesh {
        TriMesh::new(vec![Vec3::ZERO, Vec3::X, Vec3::Y], vec![0, 1, 2])
    }

    fn sheet_request() -> SectionRequest {
        request().with_role(SectionRole::Sheet)
    }

    #[test]
    fn a_sheet_request_accepts_an_open_mesh_and_checks_its_evidence() {
        let mesh = open_triangle();
        let section = run_on(&mesh, |_| {}, &sheet_request()).expect("a sheet may be open");
        assert_eq!(section.evidence.role, SectionRole::Sheet);
        assert_eq!(section.evidence.boundary_edges, 3);
        // A chain of a sheet may end on the mesh boundary.
        run_on(
            &mesh,
            |s| s.chains[1].ends = Some([ChainEnd::MeshBoundary, ChainEnd::InPlaneRegion]),
            &sheet_request(),
        )
        .expect("a sheet chain ends on the mesh boundary");
        let cases: Vec<(&str, Mutation)> = vec![
            ("role", |s| s.evidence.role = SectionRole::Solid),
            ("boundary edges", |s| s.evidence.boundary_edges = 2),
            ("no boundary edges", |s| s.evidence.boundary_edges = 0),
        ];
        for (name, mutate) in cases {
            let error = run_on(&mesh, mutate, &sheet_request()).expect_err(name);
            assert!(
                matches!(error, GeomError::BackendContractViolation { .. }),
                "{name}: {error:?}"
            );
        }
    }

    #[test]
    fn a_sheet_request_refuses_a_defective_surface_before_dispatch() {
        let flat = |indices: Vec<u32>| {
            TriMesh::new(
                vec![
                    Vec3::ZERO,
                    Vec3::X,
                    Vec3::Y,
                    Vec3::new(1.0, 1.0, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                ],
                indices,
            )
        };
        let calls = Arc::new(AtomicUsize::new(0));
        let mut registry = MeshPlaneSectionRegistry::new();
        registry.register(
            0,
            Scripted {
                mutate: |_| {},
                calls: Arc::clone(&calls),
            },
        );
        let refusals = [
            // Both triangles run 0 -> 1: inconsistent winding.
            ("winding", flat(vec![0, 1, 2, 0, 1, 3]), false),
            // Three triangles on the edge 0-1.
            ("non-manifold", flat(vec![0, 1, 2, 1, 0, 3, 1, 0, 4]), false),
            // Collinear corners.
            (
                "degenerate",
                TriMesh::new(vec![Vec3::ZERO, Vec3::X, Vec3::X * 2.0], vec![0, 1, 2]),
                true,
            ),
        ];
        for (name, mesh, degenerate) in refusals {
            let error = registry
                .section_request(&mesh, &sheet_request(), &options())
                .expect_err(name);
            assert!(
                if degenerate {
                    matches!(error, GeomError::Degenerate(_))
                } else {
                    matches!(error, GeomError::NotManifold(_))
                },
                "{name}: {error:?}"
            );
        }
        assert_eq!(calls.load(Ordering::SeqCst), 0);
    }

    /// Implements only the legacy method, so the default refuses.
    #[derive(Debug)]
    struct Legacy;

    impl Backend for Legacy {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("section-legacy"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl MeshPlaneSection for Legacy {
        fn section(
            &self,
            _mesh: &TriMesh,
            _frame: Frame3,
            _limits: SectionLimits,
            _options: &ExecutionOptions,
        ) -> GeomResult<SectionOutcome> {
            Err(GeomError::Degenerate("legacy path not under test".into()))
        }
    }

    #[test]
    fn unsupported_falls_back_to_the_next_provider() {
        let mut registry = MeshPlaneSectionRegistry::new();
        registry.register(5, Legacy);
        assert!(matches!(
            registry.section_request(&cube(), &request(), &options()),
            Err(GeomError::Unsupported { .. })
        ));
        let calls = Arc::new(AtomicUsize::new(0));
        registry.register(
            1,
            Scripted {
                mutate: |_| {},
                calls: Arc::clone(&calls),
            },
        );
        registry
            .section_request(&cube(), &request(), &options())
            .expect("falls back past the provider without the method");
        assert_eq!(calls.load(Ordering::SeqCst), 1);
    }
}
