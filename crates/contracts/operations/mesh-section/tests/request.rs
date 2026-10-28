//! Request types, evidence counts and the region conformance suite.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionOptions, ExecutionTarget, GeomError,
    GeomResult, Operation,
};
use axiolid_core::{Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_mesh::TriMesh;
use axiolid_mesh_section_contract::conformance::{
    RegionConformanceFailure, RegionConformanceSuite,
};
use axiolid_mesh_section_contract::{
    ChainEnd, InPlaneRegion, InPlaneSide, MeshPlaneSection, PlaneSection, PlaneSectionEvidence,
    RegionLimits, SectionChain, SectionLimits, SectionOutcome, SectionRequest, SectionRole,
};

fn world() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn square() -> Vec<Point2> {
    vec![
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ]
}

#[test]
fn constructors_keep_their_fields() {
    let limits = SectionLimits::new(1, 2, 3, 4);
    let region_limits = RegionLimits::new(5, 6, 7);
    assert_eq!(
        (
            region_limits.max_regions,
            region_limits.max_region_vertices,
            region_limits.max_in_plane_pairs
        ),
        (5, 6, 7)
    );
    let request = SectionRequest::new(world(), limits, region_limits);
    assert_eq!(request.limits, limits);
    assert_eq!(request.region_limits, region_limits);
    assert_eq!(request.frame, world());
    assert_eq!(request.role, SectionRole::Solid);
    assert_eq!(
        request.with_role(SectionRole::Sheet).role,
        SectionRole::Sheet
    );

    let closed = SectionChain::closed(square());
    assert!(closed.is_closed());
    assert_eq!(closed.ends, None);
    let open = SectionChain::open(
        vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)],
        ChainEnd::InPlaneRegion,
        ChainEnd::InPlaneRegion,
    );
    assert!(!open.is_closed());
    assert_eq!(open.ends, Some([ChainEnd::InPlaneRegion; 2]));
    let ending_on_the_mesh = SectionChain::open(
        vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)],
        ChainEnd::MeshBoundary,
        ChainEnd::InPlaneRegion,
    );
    assert_eq!(
        ending_on_the_mesh.ends,
        Some([ChainEnd::MeshBoundary, ChainEnd::InPlaneRegion])
    );

    let region = InPlaneRegion::new(InPlaneSide::Below, square(), vec![]);
    assert_eq!(region.side, InPlaneSide::Below);
    assert!(InPlaneSide::Below < InPlaneSide::Above && InPlaneSide::Above < InPlaneSide::Both);
}

#[test]
fn evidence_counts_follow_the_output() {
    let hole = vec![
        Point2::new(0.2, 0.2),
        Point2::new(0.2, 0.4),
        Point2::new(0.4, 0.4),
    ];
    let chains = vec![
        SectionChain::closed(square()),
        SectionChain::open(
            vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)],
            ChainEnd::InPlaneRegion,
            ChainEnd::InPlaneRegion,
        ),
    ];
    let regions = vec![InPlaneRegion::new(InPlaneSide::Above, square(), vec![hole])];
    let evidence = PlaneSectionEvidence::input_mesh(12)
        .with_source_counts(2, 1, 3)
        .with_output(&chains, &regions);
    assert!(evidence.is_derived_from_input_mesh());
    assert_eq!(evidence.source_triangles, 12);
    assert_eq!(
        (
            evidence.in_plane_triangles,
            evidence.touching_edges,
            evidence.in_plane_pairs
        ),
        (2, 1, 3)
    );
    assert_eq!((evidence.closed_chains, evidence.open_chains), (1, 1));
    assert_eq!(evidence.chain_vertices, 6);
    assert_eq!((evidence.regions, evidence.region_vertices), (1, 7));

    let empty = PlaneSectionEvidence::input_mesh(3);
    assert_eq!(empty.with_output(&[], &[]), empty);
    assert_eq!((empty.role, empty.boundary_edges), (SectionRole::Solid, 0));
    let sheet = empty.with_role(SectionRole::Sheet, 8);
    assert_eq!((sheet.role, sheet.boundary_edges), (SectionRole::Sheet, 8));
    assert_ne!(sheet, empty);
    assert_eq!(sheet.with_output(&[], &[]), sheet);
}

/// A provider that implements only the legacy method.
#[derive(Debug)]
struct LegacyOnly;

impl Backend for LegacyOnly {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new("legacy-only"), ExecutionTarget::PortableCpu)
    }
}

impl MeshPlaneSection for LegacyOnly {
    fn section(
        &self,
        _mesh: &TriMesh,
        _frame: Frame3,
        _limits: SectionLimits,
        _options: &ExecutionOptions,
    ) -> GeomResult<SectionOutcome> {
        Err(GeomError::Degenerate("not under test".into()))
    }
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn request() -> SectionRequest {
    SectionRequest::new(
        world(),
        SectionLimits::new(10, 10, 10, 10),
        RegionLimits::new(10, 10, 10),
    )
}

#[test]
fn the_default_method_is_unsupported_and_keeps_the_legacy_scratch() {
    let mesh = TriMesh::new(vec![], vec![]);
    let error = LegacyOnly
        .section_request(&mesh, &request(), &options())
        .expect_err("default refuses");
    assert!(matches!(
        error,
        GeomError::Unsupported {
            operation: Operation::MeshPlaneSection,
            ..
        }
    ));
    assert_eq!(
        LegacyOnly.request_scratch_requirement(),
        LegacyOnly.scratch_requirement()
    );
}

#[test]
fn a_provider_without_the_method_passes_unclaimed() {
    let report = RegionConformanceSuite::run(&LegacyOnly);
    assert!(!report.claimed);
    assert!(report.is_success());
}

/// Claims the method and answers every case with an empty section.
#[derive(Debug)]
struct EmptyAnswers;

impl Backend for EmptyAnswers {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("empty-answers"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl MeshPlaneSection for EmptyAnswers {
    fn section(
        &self,
        mesh: &TriMesh,
        frame: Frame3,
        limits: SectionLimits,
        options: &ExecutionOptions,
    ) -> GeomResult<SectionOutcome> {
        LegacyOnly.section(mesh, frame, limits, options)
    }

    fn section_request(
        &self,
        mesh: &TriMesh,
        request: &SectionRequest,
        _options: &ExecutionOptions,
    ) -> GeomResult<PlaneSection> {
        Ok(PlaneSection::new(
            request.frame,
            vec![],
            vec![],
            PlaneSectionEvidence::input_mesh(mesh.triangle_count()),
        ))
    }
}

#[test]
fn a_wrong_answer_fails_with_the_case_name() {
    let report = RegionConformanceSuite::run(&EmptyAnswers);
    assert!(report.claimed);
    assert!(report.failures.iter().any(|failure| matches!(
        failure,
        RegionConformanceFailure::WrongRegions { case } if *case == "cube through its top face"
    )));
    assert!(report.failures.iter().any(|failure| matches!(
        failure,
        RegionConformanceFailure::WrongChains { case } if *case == "parapet with a post"
    )));
}
