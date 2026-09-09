//! Supported application boundary for the portable v0.4 workflows.
//!
//! Concrete providers stay private. Selection is explicit on
//! [`ApplicationBuilder`](crate::application::ApplicationBuilder),
//! while results and errors use only public contract and representation types.

use core::fmt;

use axiolid_brep::ExactBRep;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_contracts::{
    capability_ids, BackendId, CapabilityDescriptor, Exactness, ExecutionOptions,
    IntegrationDescriptor, Operation, Representation,
};
use axiolid_core::{BooleanOperator, Frame3, Ray3, Scalar, Tolerance, Vec3};
use axiolid_dispatch::{MeshBooleanRegistry, MeshPlaneSectionRegistry};
use axiolid_measure::{
    surface_properties, volume_properties, MeshMeasureError, SurfaceProperties, VolumeProperties,
};
use axiolid_mesh::{try_audit_mesh, MeshAuditError, MeshHealth, TriMesh};
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use axiolid_mesh_boolean_contract::{
    conformance::ConformanceReport as BooleanConformanceReport, BooleanOutcome,
};
use axiolid_mesh_section_contract::{
    conformance::ConformanceReport as SectionConformanceReport, SectionLimits, SectionOutcome,
};
use axiolid_profile::Profile;
use axiolid_ray_mesh::{nearest_hit, RayHit3, RayMeshError};
use axiolid_reference::ScalarSection;

const APPLICATION_ID: BackendId = BackendId::new("axiolid-application");
const MEASURE_ID: BackendId = BackendId::new("scalar-measure");
const RAY_MESH_ID: BackendId = BackendId::new("scalar-ray-mesh");

static MESH_INPUT: &[Representation] = &[Representation::TriangleMesh];
static MESH_PAIR_INPUT: &[Representation] =
    &[Representation::TriangleMesh, Representation::TriangleMesh];
static PROFILE_INPUT: &[Representation] = &[Representation::Profile2d];

/// Surface and signed-volume properties computed under one tolerance policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshMeasurements {
    pub surface: SurfaceProperties,
    pub volume: VolumeProperties,
}

/// Context retained by every application-level failure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallContext {
    pub operation: Operation,
    pub provider: BackendId,
    pub tolerance: Tolerance,
}

/// Typed lower-level reason an application call failed.
#[derive(Debug)]
pub enum ApplicationErrorSource {
    Geometry(axiolid_contracts::GeomError),
    MeshAudit(Box<MeshAuditError>),
    Measurement(Box<MeshMeasureError>),
    RayMesh(RayMeshError),
}

impl fmt::Display for ApplicationErrorSource {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Geometry(error) => error.fmt(formatter),
            Self::MeshAudit(error) => error.fmt(formatter),
            Self::Measurement(error) => error.fmt(formatter),
            Self::RayMesh(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for ApplicationErrorSource {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Geometry(error) => Some(error),
            Self::MeshAudit(error) => Some(error),
            Self::Measurement(error) => Some(error),
            Self::RayMesh(error) => Some(error),
        }
    }
}

/// Application failure with operation, selected provider, and tolerance intact.
#[derive(Debug)]
pub struct ApplicationError {
    pub context: CallContext,
    pub source: Box<ApplicationErrorSource>,
}

impl fmt::Display for ApplicationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{:?} via `{}` at tolerance {:?}: {}",
            self.context.operation, self.context.provider, self.context.tolerance, self.source
        )
    }
}

impl std::error::Error for ApplicationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(self.source.as_ref())
    }
}

/// Shared conformance report without provider implementation types.
#[derive(Debug, Clone)]
pub enum ProviderConformanceReport {
    MeshBoolean(Box<BooleanConformanceReport>),
    MeshSection(Box<SectionConformanceReport>),
}

impl fmt::Display for ProviderConformanceReport {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MeshBoolean(report) => report.fmt(formatter),
            Self::MeshSection(report) => report.fmt(formatter),
        }
    }
}

/// Portable-provider conformance failed before the application became usable.
#[derive(Debug, Clone)]
pub struct ApplicationBuildError {
    pub operation: Operation,
    pub provider: BackendId,
    pub report: ProviderConformanceReport,
}

impl fmt::Display for ApplicationBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "provider `{}` failed {:?} registration: {}",
            self.provider, self.operation, self.report
        )
    }
}

impl std::error::Error for ApplicationBuildError {}

/// Explicit assembly of the supported application boundary.
#[derive(Debug, Clone, Default)]
pub struct ApplicationBuilder {
    boolean: MeshBooleanRegistry,
    section: MeshPlaneSectionRegistry,
    boolean_provider: Option<BackendId>,
    section_provider: Option<BackendId>,
}

impl ApplicationBuilder {
    #[must_use]
    pub const fn new() -> Self {
        Self {
            boolean: MeshBooleanRegistry::new(),
            section: MeshPlaneSectionRegistry::new(),
            boolean_provider: None,
            section_provider: None,
        }
    }

    /// Select the pure-Rust boolmesh provider after running shared conformance.
    pub fn with_portable_boolean(mut self) -> Result<Self, ApplicationBuildError> {
        let provider = BoolmeshBoolean::new();
        self.boolean
            .register_conformant(0, provider)
            .map_err(|report| ApplicationBuildError {
                operation: Operation::MeshBoolean,
                provider: BoolmeshBoolean::ID,
                report: ProviderConformanceReport::MeshBoolean(report),
            })?;
        self.boolean_provider = Some(BoolmeshBoolean::ID);
        Ok(self)
    }

    /// Select the deterministic scalar mesh-section provider after shared conformance.
    pub fn with_portable_section(mut self) -> Result<Self, ApplicationBuildError> {
        let provider = ScalarSection::new();
        self.section
            .register_conformant(0, provider)
            .map_err(|report| ApplicationBuildError {
                operation: Operation::MeshPlaneSection,
                provider: ScalarSection::ID,
                report: ProviderConformanceReport::MeshSection(report),
            })?;
        self.section_provider = Some(ScalarSection::ID);
        Ok(self)
    }

    #[must_use]
    pub fn build(self) -> Application {
        let mut descriptor = crate::integration::descriptor();
        descriptor.representations.extend([
            Representation::MeshHealth,
            Representation::Measurements,
            Representation::RayHit,
        ]);
        descriptor.representations.sort_unstable();
        descriptor.representations.dedup();
        descriptor.capabilities.extend(direct_capabilities());
        if let Some(provider) = self.boolean_provider {
            descriptor.capabilities.push(CapabilityDescriptor {
                id: capability_ids::MESH_BOOLEAN,
                operation: Operation::MeshBoolean,
                provider,
                required_feature: "application",
                inputs: MESH_PAIR_INPUT,
                output: Representation::TriangleMesh,
                exactness: Exactness::ToleranceBounded,
                deterministic: true,
            });
        }
        if let Some(provider) = self.section_provider {
            descriptor.capabilities.push(CapabilityDescriptor {
                id: capability_ids::MESH_SECTION,
                operation: Operation::MeshPlaneSection,
                provider,
                required_feature: "application",
                inputs: MESH_INPUT,
                output: Representation::Profile2d,
                exactness: Exactness::ToleranceBounded,
                deterministic: true,
            });
        }
        descriptor
            .capabilities
            .sort_by_key(|capability| capability.id);
        Application {
            boolean: self.boolean,
            section: self.section,
            boolean_provider: self.boolean_provider,
            section_provider: self.section_provider,
            descriptor,
        }
    }
}

fn direct_capabilities() -> [CapabilityDescriptor; 4] {
    [
        CapabilityDescriptor {
            id: capability_ids::MESH_VALIDATE,
            operation: Operation::Healing,
            provider: APPLICATION_ID,
            required_feature: "application",
            inputs: MESH_INPUT,
            output: Representation::MeshHealth,
            exactness: Exactness::ToleranceBounded,
            deterministic: true,
        },
        CapabilityDescriptor {
            id: capability_ids::MESH_MEASURE,
            operation: Operation::Measurement,
            provider: MEASURE_ID,
            required_feature: "application",
            inputs: MESH_INPUT,
            output: Representation::Measurements,
            exactness: Exactness::ToleranceBounded,
            deterministic: true,
        },
        CapabilityDescriptor {
            id: capability_ids::RAY_MESH,
            operation: Operation::SpatialQuery,
            provider: RAY_MESH_ID,
            required_feature: "application",
            inputs: MESH_INPUT,
            output: Representation::RayHit,
            exactness: Exactness::ToleranceBounded,
            deterministic: true,
        },
        CapabilityDescriptor {
            id: capability_ids::EXACT_EXTRUDE,
            operation: Operation::Sweep,
            provider: axiolid_construct::BACKEND_ID,
            required_feature: "application",
            inputs: PROFILE_INPUT,
            output: Representation::ExactBrep,
            exactness: Exactness::Exact,
            deterministic: true,
        },
    ]
}

/// Coherent application boundary backed by explicitly selected providers.
#[derive(Debug, Clone)]
pub struct Application {
    boolean: MeshBooleanRegistry,
    section: MeshPlaneSectionRegistry,
    boolean_provider: Option<BackendId>,
    section_provider: Option<BackendId>,
    descriptor: IntegrationDescriptor,
}

impl Application {
    /// Build the supported pure-Rust provider set explicitly.
    pub fn portable() -> Result<Self, ApplicationBuildError> {
        Ok(ApplicationBuilder::new()
            .with_portable_boolean()?
            .with_portable_section()?
            .build())
    }

    #[must_use]
    pub fn descriptor(&self) -> &IntegrationDescriptor {
        &self.descriptor
    }

    pub fn validate_mesh(
        &self,
        mesh: &TriMesh,
        tolerance: Tolerance,
    ) -> Result<MeshHealth, ApplicationError> {
        try_audit_mesh(mesh, tolerance).map_err(|source| ApplicationError {
            context: CallContext {
                operation: Operation::Healing,
                provider: APPLICATION_ID,
                tolerance,
            },
            source: Box::new(ApplicationErrorSource::MeshAudit(Box::new(source))),
        })
    }

    pub fn measure_mesh(
        &self,
        mesh: &TriMesh,
        tolerance: Tolerance,
    ) -> Result<MeshMeasurements, ApplicationError> {
        let context = CallContext {
            operation: Operation::Measurement,
            provider: MEASURE_ID,
            tolerance,
        };
        let surface = surface_properties(mesh, tolerance).map_err(|source| ApplicationError {
            context,
            source: Box::new(ApplicationErrorSource::Measurement(Box::new(source))),
        })?;
        let volume = volume_properties(mesh, tolerance).map_err(|source| ApplicationError {
            context,
            source: Box::new(ApplicationErrorSource::Measurement(Box::new(source))),
        })?;
        Ok(MeshMeasurements { surface, volume })
    }

    pub fn boolean(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operator: BooleanOperator,
        options: &ExecutionOptions,
    ) -> Result<BooleanOutcome, ApplicationError> {
        self.boolean
            .boolean(subject, tool, operator, options)
            .map_err(|source| {
                self.geometry_error(
                    Operation::MeshBoolean,
                    self.boolean_provider.unwrap_or(APPLICATION_ID),
                    options.tolerance(),
                    source,
                )
            })
    }

    pub fn subtract_many(
        &self,
        subject: &TriMesh,
        tools: &[TriMesh],
        options: &ExecutionOptions,
    ) -> Result<BooleanOutcome, ApplicationError> {
        self.boolean
            .subtract_many(subject, tools, options)
            .map_err(|source| {
                self.geometry_error(
                    Operation::MeshBoolean,
                    self.boolean_provider.unwrap_or(APPLICATION_ID),
                    options.tolerance(),
                    source,
                )
            })
    }

    /// Union many solids in one batch.
    ///
    /// Providers may choose the reduction order; `boolmesh` reduces in a
    /// balanced tree rather than folding left, which is materially faster
    /// as the operand count grows.
    pub fn union_many(
        &self,
        solids: &[TriMesh],
        options: &ExecutionOptions,
    ) -> Result<BooleanOutcome, ApplicationError> {
        self.boolean.union_many(solids, options).map_err(|source| {
            self.geometry_error(
                Operation::MeshBoolean,
                self.boolean_provider.unwrap_or(APPLICATION_ID),
                options.tolerance(),
                source,
            )
        })
    }

    pub fn section_mesh(
        &self,
        mesh: &TriMesh,
        frame: Frame3,
        limits: SectionLimits,
        options: &ExecutionOptions,
    ) -> Result<SectionOutcome, ApplicationError> {
        self.section
            .section(mesh, frame, limits, options)
            .map_err(|source| {
                self.geometry_error(
                    Operation::MeshPlaneSection,
                    self.section_provider.unwrap_or(APPLICATION_ID),
                    options.tolerance(),
                    source,
                )
            })
    }

    pub fn nearest_mesh_hit(
        &self,
        mesh: &TriMesh,
        ray: &Ray3,
        tolerance: Tolerance,
    ) -> Result<Option<RayHit3>, ApplicationError> {
        nearest_hit(mesh, ray, tolerance).map_err(|source| ApplicationError {
            context: CallContext {
                operation: Operation::SpatialQuery,
                provider: RAY_MESH_ID,
                tolerance,
            },
            source: Box::new(ApplicationErrorSource::RayMesh(source)),
        })
    }

    pub fn extrude_profile_exact(
        &self,
        profile: &Profile,
        direction: Vec3,
        depth: Scalar,
        tolerance: Tolerance,
    ) -> Result<ExactBRep, ApplicationError> {
        extrude_profile_exact(profile, direction, depth, tolerance).map_err(|source| {
            self.geometry_error(
                Operation::Sweep,
                axiolid_construct::BACKEND_ID,
                tolerance,
                source,
            )
        })
    }

    fn geometry_error(
        &self,
        operation: Operation,
        provider: BackendId,
        tolerance: Tolerance,
        source: axiolid_contracts::GeomError,
    ) -> ApplicationError {
        ApplicationError {
            context: CallContext {
                operation,
                provider,
                tolerance,
            },
            source: Box::new(ApplicationErrorSource::Geometry(source)),
        }
    }
}
