//! Versioned, provider-neutral contracts for downstream capability discovery.
//!
//! This module describes what a built integration surface can do. Operation
//! traits remain the capability truth: a facade or FFI layer must only publish
//! a [`CapabilityDescriptor`] after constructing and exercising the concrete
//! provider behind it.

use thiserror::Error;

use crate::{BackendId, CapabilityId, Operation};

/// First downstream integration protocol shipped by Axiolid.
pub const INTEGRATION_API_VERSION: ApiVersion = ApiVersion::new(0, 4, 0);
/// Minimum supported Rust compiler for Rust integration profiles.
pub const MINIMUM_RUST_VERSION: &str = "1.99";

/// Semantic version of an integration API or ABI.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ApiVersion {
    pub major: u16,
    pub minor: u16,
    pub patch: u16,
}

impl ApiVersion {
    pub const fn new(major: u16, minor: u16, patch: u16) -> Self {
        Self {
            major,
            minor,
            patch,
        }
    }

    /// Protocol compatibility is major-version stable and backward compatible.
    pub const fn supports(self, requested: Self) -> bool {
        self.major == requested.major
            && (self.minor > requested.minor
                || (self.minor == requested.minor && self.patch >= requested.patch))
    }
}

/// Supported way an application reaches Axiolid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IntegrationProfile {
    /// A Rust application selects only the leaf packages it needs.
    RustLeaf,
    /// A Rust application uses the feature-gated `axiolid` facade.
    RustFacade,
    /// A native application uses the versioned C ABI.
    NativeC,
}

/// Portable representation families visible at an integration boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Representation {
    Scalar,
    Linear,
    Profile2d,
    AnalyticCurve,
    AnalyticSurface,
    Topology,
    ExactBrep,
    TriangleMesh,
    MeshHealth,
    Measurements,
    RayHit,
    ModelGraph,
    SampledField,
}

/// Fidelity a capability promises for its output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Exactness {
    /// Analytic identity and certified trims are preserved.
    Exact,
    /// The caller supplies an explicit tolerance and receives a bounded result.
    ToleranceBounded,
}

impl Exactness {
    const fn satisfies(self, required: Self) -> bool {
        matches!(
            (self, required),
            (Self::Exact, _) | (Self::ToleranceBounded, Self::ToleranceBounded)
        )
    }
}

/// Cross-boundary ownership model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Ownership {
    /// Normal Rust values own their allocations.
    RustValues,
    /// Native outputs are opaque handles released by matching Axiolid functions.
    OpaqueOwnedHandles,
}

/// Thread-safety promise of one integration profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ThreadSafety {
    /// Public values are `Send + Sync`; operation-local mutable state is not shared.
    SendSyncValues,
    /// A native context may move between threads but cannot be used concurrently.
    ContextSendNotSync,
}

/// Unit, coordinate, tolerance, ownership, and execution semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BoundaryContract {
    /// Coordinates are right-handed Cartesian `f64` values.
    pub right_handed_cartesian_f64: bool,
    /// Geometry is unitless; every input and tolerance uses one caller-selected unit.
    pub caller_defined_consistent_units: bool,
    /// Approximate operations require an explicit tolerance.
    pub explicit_tolerance: bool,
    pub ownership: Ownership,
    pub thread_safety: ThreadSafety,
}

impl BoundaryContract {
    pub const fn rust() -> Self {
        Self {
            right_handed_cartesian_f64: true,
            caller_defined_consistent_units: true,
            explicit_tolerance: true,
            ownership: Ownership::RustValues,
            thread_safety: ThreadSafety::SendSyncValues,
        }
    }

    pub const fn native() -> Self {
        Self {
            ownership: Ownership::OpaqueOwnedHandles,
            thread_safety: ThreadSafety::ContextSendNotSync,
            ..Self::rust()
        }
    }
}

/// One capability backed by a concrete provider in this build.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CapabilityDescriptor {
    pub id: CapabilityId,
    pub operation: Operation,
    pub provider: BackendId,
    /// Cargo feature that admitted the provider into this build.
    pub required_feature: &'static str,
    pub inputs: &'static [Representation],
    pub output: Representation,
    pub exactness: Exactness,
    pub deterministic: bool,
}

/// Minimum behavior a caller needs before it submits geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CapabilityRequirement {
    pub id: CapabilityId,
    pub output: Representation,
    pub exactness: Exactness,
    pub deterministic: bool,
}

/// Typed reason a capability handshake refused a request.
#[derive(Debug, Clone, Copy, Error, PartialEq, Eq)]
pub enum RequirementRefusal {
    #[error("integration API {requested:?} is incompatible with advertised {advertised:?}")]
    ApiVersionUnavailable {
        requested: ApiVersion,
        advertised: ApiVersion,
    },
    #[error("capability `{capability}` is not available")]
    CapabilityUnavailable { capability: CapabilityId },
    #[error("capability `{capability}` returns {advertised:?}, not required {required:?}")]
    RepresentationUnavailable {
        capability: CapabilityId,
        required: Representation,
        advertised: Representation,
    },
    #[error("capability `{capability}` advertises {advertised:?}, not required {required:?}")]
    ExactnessUnavailable {
        capability: CapabilityId,
        required: Exactness,
        advertised: Exactness,
    },
    #[error("capability `{capability}` is not deterministic in this build")]
    DeterminismUnavailable { capability: CapabilityId },
}

/// Complete handshake returned by one compiled integration surface.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntegrationDescriptor {
    pub api_version: ApiVersion,
    pub abi_version: Option<ApiVersion>,
    pub profile: IntegrationProfile,
    pub minimum_rust_version: Option<&'static str>,
    pub enabled_features: Vec<&'static str>,
    pub representations: Vec<Representation>,
    pub capabilities: Vec<CapabilityDescriptor>,
    pub boundary: BoundaryContract,
}

impl IntegrationDescriptor {
    /// Honest baseline: a profile with no advertised representation or operation.
    pub fn empty(profile: IntegrationProfile) -> Self {
        let native = matches!(profile, IntegrationProfile::NativeC);
        Self {
            api_version: INTEGRATION_API_VERSION,
            abi_version: native.then_some(INTEGRATION_API_VERSION),
            profile,
            minimum_rust_version: (!native).then_some(MINIMUM_RUST_VERSION),
            enabled_features: Vec::new(),
            representations: Vec::new(),
            capabilities: Vec::new(),
            boundary: if native {
                BoundaryContract::native()
            } else {
                BoundaryContract::rust()
            },
        }
    }

    pub fn supports_api(&self, requested: ApiVersion) -> bool {
        self.api_version.supports(requested)
    }

    pub fn require_api(&self, requested: ApiVersion) -> Result<(), RequirementRefusal> {
        if self.supports_api(requested) {
            Ok(())
        } else {
            Err(RequirementRefusal::ApiVersionUnavailable {
                requested,
                advertised: self.api_version,
            })
        }
    }

    /// Return the concrete provider descriptor or a typed refusal.
    pub fn require(
        &self,
        requirement: CapabilityRequirement,
    ) -> Result<&CapabilityDescriptor, RequirementRefusal> {
        let Some(capability) = self
            .capabilities
            .iter()
            .find(|candidate| candidate.id == requirement.id)
        else {
            return Err(RequirementRefusal::CapabilityUnavailable {
                capability: requirement.id,
            });
        };

        if capability.output != requirement.output {
            return Err(RequirementRefusal::RepresentationUnavailable {
                capability: requirement.id,
                required: requirement.output,
                advertised: capability.output,
            });
        }
        if !capability.exactness.satisfies(requirement.exactness) {
            return Err(RequirementRefusal::ExactnessUnavailable {
                capability: requirement.id,
                required: requirement.exactness,
                advertised: capability.exactness,
            });
        }
        if requirement.deterministic && !capability.deterministic {
            return Err(RequirementRefusal::DeterminismUnavailable {
                capability: requirement.id,
            });
        }
        Ok(capability)
    }
}
