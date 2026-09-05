#![cfg(feature = "pointcloud-reconstruction")]
//! Dispatch policy for pointcloud reconstruction.
//!
//! These check *routing*, not geometry: which provider runs, when fallback
//! happens, and — most importantly — when it must not.

use std::sync::Arc;

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, ExecutionOptions, ExecutionTarget, GeomError,
    GeomResult, ScratchRequirement,
};
use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_dispatch::PointcloudReconstructionRegistry;
use axiolid_pointcloud::PointCloud;
use axiolid_pointcloud_reconstruction_contract::{
    PointcloudReconstruction, Reconstruction, ReconstructionEvidence, ReconstructionOutcome,
    ReconstructionRefusal, ReconstructionRequest,
};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn cloud() -> PointCloud {
    let points: Vec<Point3> = (0..64)
        .map(|i| {
            let a = i as Scalar * 0.4;
            Point3::new(a.cos(), a.sin(), (i % 5) as Scalar * 0.3)
        })
        .collect();
    PointCloud::new(points).expect("finite")
}

/// A provider that always produces a marked surface, so tests can tell
/// which one ran.
#[derive(Debug, Clone, Copy)]
struct Marker(u32);

impl Backend for Marker {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new(match self.0 {
                1 => "test.marker.one",
                2 => "test.marker.two",
                _ => "test.marker.other",
            }),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl PointcloudReconstruction for Marker {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::PerElement {
            bytes_per_element: 1,
        }
    }

    fn minimum_points(&self) -> usize {
        1
    }

    fn reconstruct(
        &self,
        _cloud: &PointCloud,
        _request: &ReconstructionRequest,
        _options: &ExecutionOptions,
    ) -> GeomResult<Reconstruction> {
        let mut evidence = ReconstructionEvidence::measured();
        // The marker value rides in a counter so the test can identify the
        // provider that answered.
        evidence.used_points = self.0 as usize;
        Ok(Reconstruction::Surface(Box::new(
            ReconstructionOutcome::new(
                axiolid_mesh::TriMesh::new(Vec::new(), Vec::new()),
                evidence,
            ),
        )))
    }
}

/// A provider that cannot take the work at all.
#[derive(Debug, Clone, Copy)]
struct Unavailable;

impl Backend for Unavailable {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("test.unavailable"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl PointcloudReconstruction for Unavailable {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::PerElement {
            bytes_per_element: 1,
        }
    }

    fn minimum_points(&self) -> usize {
        1
    }

    fn reconstruct(
        &self,
        _cloud: &PointCloud,
        _request: &ReconstructionRequest,
        _options: &ExecutionOptions,
    ) -> GeomResult<Reconstruction> {
        Err(GeomError::Unavailable {
            backend: BackendId::new("test.unavailable"),
            reason: "test backend is offline".to_owned(),
        })
    }
}

/// A provider that refuses on the data rather than failing.
#[derive(Debug, Clone, Copy)]
struct Refuser;

impl Backend for Refuser {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new("test.refuser"), ExecutionTarget::PortableCpu)
    }
}

impl PointcloudReconstruction for Refuser {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::PerElement {
            bytes_per_element: 1,
        }
    }

    fn minimum_points(&self) -> usize {
        1
    }

    fn reconstruct(
        &self,
        _cloud: &PointCloud,
        _request: &ReconstructionRequest,
        _options: &ExecutionOptions,
    ) -> GeomResult<Reconstruction> {
        Ok(Reconstruction::Refused(
            ReconstructionRefusal::CannotClose {
                detail: "capture has a hole".to_owned(),
            },
        ))
    }
}

/// An empty registry refuses by name rather than panicking or returning a
/// silently empty surface.
#[test]
fn an_empty_registry_reports_unsupported() {
    let registry = PointcloudReconstructionRegistry::new();
    let error = registry
        .reconstruct(&cloud(), &ReconstructionRequest::default(), &options())
        .expect_err("an empty registry cannot reconstruct");
    assert!(matches!(error, GeomError::Unsupported { .. }));
}

/// Higher priority runs first.
#[test]
fn the_highest_priority_provider_runs() {
    let mut registry = PointcloudReconstructionRegistry::new();
    registry.register(1, Marker(1));
    registry.register(10, Marker(2));

    let result = registry
        .reconstruct(&cloud(), &ReconstructionRequest::default(), &options())
        .expect("runs");
    assert_eq!(result.outcome().expect("surface").evidence.used_points, 2);
}

/// An unavailable provider falls through to the next one.
#[test]
fn an_unavailable_provider_falls_back() {
    let mut registry = PointcloudReconstructionRegistry::new();
    registry.register(10, Unavailable);
    registry.register(1, Marker(1));

    let result = registry
        .reconstruct(&cloud(), &ReconstructionRequest::default(), &options())
        .expect("falls back");
    assert_eq!(result.outcome().expect("surface").evidence.used_points, 1);
}

/// A refusal is an answer about the data, so dispatch must return it rather
/// than shopping for a provider willing to guess. Without this, adding a
/// lenient provider would silently override an honest "this cannot be
/// closed" with fabricated surface.
#[test]
fn a_refusal_is_not_a_fallback_trigger() {
    let mut registry = PointcloudReconstructionRegistry::new();
    registry.register(10, Refuser);
    registry.register(1, Marker(1));

    let result = registry
        .reconstruct(&cloud(), &ReconstructionRequest::default(), &options())
        .expect("runs");
    assert!(
        matches!(
            result.refusal(),
            Some(ReconstructionRefusal::CannotClose { .. })
        ),
        "the refusal must be returned, not overridden by a lower-priority provider"
    );
}

/// Dropping the provider from the registry must be detectable. This is the
/// mutation the issue asks for: a registry that answers with no providers
/// registered would mean dispatch was not consulting the registry at all.
#[test]
fn a_dropped_provider_is_caught() {
    let mut registry = PointcloudReconstructionRegistry::new();
    registry.register(1, Marker(1));
    assert!(!registry.is_empty());
    assert!(registry
        .reconstruct(&cloud(), &ReconstructionRequest::default(), &options())
        .is_ok());

    let empty = PointcloudReconstructionRegistry::new();
    assert!(empty.is_empty());
    assert!(
        empty
            .reconstruct(&cloud(), &ReconstructionRequest::default(), &options())
            .is_err(),
        "an empty registry must fail; if it succeeds, dispatch is not using the registry"
    );
}

/// A provider exceeding the caller's memory budget is excluded before it
/// can allocate, and a leaner one still runs.
#[test]
fn an_over_budget_provider_is_excluded_before_dispatch() {
    #[derive(Debug, Clone, Copy)]
    struct Hungry;
    impl Backend for Hungry {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("test.hungry"), ExecutionTarget::PortableCpu)
        }
    }
    impl PointcloudReconstruction for Hungry {
        fn scratch_requirement(&self) -> ScratchRequirement {
            ScratchRequirement::Unbounded
        }
        fn minimum_points(&self) -> usize {
            1
        }
        fn reconstruct(
            &self,
            _cloud: &PointCloud,
            _request: &ReconstructionRequest,
            _options: &ExecutionOptions,
        ) -> GeomResult<Reconstruction> {
            panic!("an over-budget provider must never be dispatched");
        }
    }

    let mut registry = PointcloudReconstructionRegistry::new();
    registry.register(10, Hungry);
    registry.register(1, Marker(1));

    let bounded = ExecutionOptions::new(Tolerance::METRE).with_memory_budget(1024);
    let result = registry
        .reconstruct(&cloud(), &ReconstructionRequest::default(), &bounded)
        .expect("the lean provider runs");
    assert_eq!(result.outcome().expect("surface").evidence.used_points, 1);
}

/// Registration through the conformance door rejects a provider that
/// violates the contract, instead of admitting it and failing later.
#[test]
fn a_non_conformant_provider_is_rejected_at_registration() {
    /// Claims a huge minimum but reconstructs from anything: its advertised
    /// threshold is a lie, so pre-flight checks against it are worthless.
    #[derive(Debug, Clone, Copy)]
    struct Liar;
    impl Backend for Liar {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("test.liar"), ExecutionTarget::PortableCpu)
        }
    }
    impl PointcloudReconstruction for Liar {
        fn minimum_points(&self) -> usize {
            1000
        }
        fn reconstruct(
            &self,
            _cloud: &PointCloud,
            _request: &ReconstructionRequest,
            _options: &ExecutionOptions,
        ) -> GeomResult<Reconstruction> {
            Ok(Reconstruction::Surface(Box::new(
                ReconstructionOutcome::new(
                    axiolid_mesh::TriMesh::new(Vec::new(), Vec::new()),
                    ReconstructionEvidence::measured(),
                ),
            )))
        }
    }

    let mut registry = PointcloudReconstructionRegistry::new();
    let report = registry
        .register_conformant(1, Liar)
        .expect_err("a provider that misstates its minimum must be rejected");
    assert!(!report.is_conformant());
    assert!(
        registry.is_empty(),
        "a rejected provider must not be stored"
    );
}

/// The reference provider passes the conformance door.
#[test]
fn the_reference_provider_registers_conformantly() {
    let mut registry = PointcloudReconstructionRegistry::new();
    registry
        .register_conformant(
            1,
            axiolid_pointcloud_reconstruction_sdf::SdfReconstruction::new(),
        )
        .expect("the reference provider must be conformant");
    assert!(!registry.is_empty());
}

/// A shared trait object registers like an owned provider.
#[test]
fn a_shared_provider_registers() {
    let mut registry = PointcloudReconstructionRegistry::new();
    let provider: Arc<dyn PointcloudReconstruction> = Arc::new(Marker(2));
    registry.register_arc(5, provider);
    assert_eq!(registry.providers().count(), 1);
}
