//! Versioned C ABI for Axiolid's supported application facade.
//!
//! ABI functions never unwind. Handles are globally unique scalar tokens; ownership
//! transfers are documented per function and no Rust allocation crosses the boundary.

#![deny(unsafe_op_in_unsafe_fn)]

use std::collections::HashMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, LazyLock, Mutex, MutexGuard};

use axiolid::application::Application;
use axiolid::core::Point3;
use axiolid::mesh::TriMesh;

/// v0.4 ABI status code.
#[non_exhaustive]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidStatus {
    Ok = 0,
    NullPointer = 1,
    InvalidArgument = 2,
    InvalidHandle = 3,
    LimitExceeded = 4,
    Unsupported = 5,
    UnsupportedExact = 6,
    BufferTooSmall = 7,
    OperationFailed = 8,
    WrongResultKind = 9,
    NoError = 10,
    Panic = 255,
}

/// Semantic ABI and package version.
#[repr(C)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidVersion {
    pub abi_major: u16,
    pub abi_minor: u16,
    pub abi_patch: u16,
    pub crate_major: u16,
    pub crate_minor: u16,
    pub crate_patch: u16,
}

/// Globally unique opaque context token. Zero is always invalid.
#[repr(transparent)]
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AxiolidContextHandle(pub u64);

impl AxiolidContextHandle {
    pub const INVALID: Self = Self(0);
}

/// Stable provider bundle selection. Integer form makes unknown C values rejectable.
pub type AxiolidProviderProfile = i32;
pub const AXIOLID_PROVIDER_PORTABLE: AxiolidProviderProfile = 1;

/// Hard allocation budgets and provider selection for one context.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidContextConfig {
    pub provider_profile: AxiolidProviderProfile,
    pub max_meshes: u32,
    pub max_results: u32,
    pub max_vertices_per_mesh: u32,
    pub max_triangles_per_mesh: u32,
}

impl Default for AxiolidContextConfig {
    fn default() -> Self {
        Self {
            provider_profile: AXIOLID_PROVIDER_PORTABLE,
            max_meshes: 1_024,
            max_results: 1_024,
            max_vertices_per_mesh: 10_000_000,
            max_triangles_per_mesh: 20_000_000,
        }
    }
}

struct Context {
    application: Application,
    config: AxiolidContextConfig,
    meshes: HandleTable<TriMesh>,
    results: HandleTable<StoredResult>,
    last_error: Option<ErrorRecord>,
}

enum StoredResult {
    Mesh(TriMesh),
    Exact(Box<axiolid::brep::ExactBRep>),
    HiddenLines(axiolid::mesh_hidden_lines::HiddenLineOutcome),
    VisibleFaces(axiolid::mesh_hidden_lines::VisibleFaceOutcome),
}

struct ErrorRecord {
    status: AxiolidStatus,
    operation: AxiolidOperation,
    tolerance: AxiolidTolerance,
    provider: String,
    message: String,
    /// The caller's id of the input mesh the error names, if it names one.
    mesh: Option<u64>,
}

struct HandleTable<T> {
    values: HashMap<u64, T>,
}

impl<T> Default for HandleTable<T> {
    fn default() -> Self {
        Self {
            values: HashMap::new(),
        }
    }
}

static NEXT_HANDLE: AtomicU64 = AtomicU64::new(1);

impl<T> HandleTable<T> {
    fn insert(&mut self, value: T) -> u64 {
        let handle = loop {
            let candidate = NEXT_HANDLE.fetch_add(1, Ordering::Relaxed);
            if candidate != 0 && !self.values.contains_key(&candidate) {
                break candidate;
            }
        };
        self.values.insert(handle, value);
        handle
    }

    fn get(&self, handle: u64) -> Option<&T> {
        self.values.get(&handle)
    }

    fn get_mut(&mut self, handle: u64) -> Option<&mut T> {
        self.values.get_mut(&handle)
    }

    fn remove(&mut self, handle: u64) -> Option<T> {
        self.values.remove(&handle)
    }

    fn live_count(&self) -> usize {
        self.values.len()
    }
}

static CONTEXTS: LazyLock<Mutex<HandleTable<Arc<Mutex<Context>>>>> =
    LazyLock::new(|| Mutex::new(HandleTable::default()));

fn lock_unpoisoned<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

fn contexts() -> MutexGuard<'static, HandleTable<Arc<Mutex<Context>>>> {
    lock_unpoisoned(&CONTEXTS)
}

fn context_entry(handle: u64) -> Option<Arc<Mutex<Context>>> {
    contexts().get(handle).cloned()
}

fn boundary(operation: impl FnOnce() -> AxiolidStatus) -> AxiolidStatus {
    catch_unwind(AssertUnwindSafe(operation)).unwrap_or(AxiolidStatus::Panic)
}

/// Write the ABI/package version to caller-owned memory.
///
/// `out_version` must be null or point to writable storage for one `AxiolidVersion`.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_version(out_version: *mut AxiolidVersion) -> AxiolidStatus {
    boundary(|| {
        if out_version.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let version = AxiolidVersion {
            abi_major: 0,
            abi_minor: 4,
            abi_patch: 1,
            crate_major: env!("CARGO_PKG_VERSION_MAJOR").parse().unwrap_or(0),
            crate_minor: env!("CARGO_PKG_VERSION_MINOR").parse().unwrap_or(0),
            crate_patch: env!("CARGO_PKG_VERSION_PATCH").parse().unwrap_or(0),
        };
        // SAFETY: null was rejected; the ABI requires writable, aligned storage for this POD type.
        unsafe { out_version.write(version) };
        AxiolidStatus::Ok
    })
}

/// Create and transfer ownership of a context handle to the caller.
///
/// Both pointers must refer to readable/writable instances for the duration of this call.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_context_create(
    config: *const AxiolidContextConfig,
    out_context: *mut AxiolidContextHandle,
) -> AxiolidStatus {
    boundary(|| {
        if config.is_null() || out_context.is_null() {
            return AxiolidStatus::NullPointer;
        }
        // SAFETY: null was rejected; caller guarantees readable aligned POD storage.
        let config = unsafe { config.read() };
        if config.provider_profile != AXIOLID_PROVIDER_PORTABLE
            || config.max_meshes == 0
            || config.max_results == 0
            || config.max_vertices_per_mesh == 0
            || config.max_triangles_per_mesh == 0
        {
            return AxiolidStatus::InvalidArgument;
        }
        let Ok(application) = Application::portable() else {
            return AxiolidStatus::OperationFailed;
        };
        let handle = AxiolidContextHandle(contexts().insert(Arc::new(Mutex::new(Context {
            application,
            config,
            meshes: HandleTable::default(),
            results: HandleTable::default(),
            last_error: None,
        }))));
        // SAFETY: null was rejected; caller guarantees writable aligned POD storage.
        unsafe { out_context.write(handle) };
        AxiolidStatus::Ok
    })
}

/// Destroy a context and all child objects it owns.
///
/// A stale or repeatedly destroyed handle is rejected without dereferencing freed memory.
#[no_mangle]
pub extern "C" fn axiolid_v0_4_context_destroy(context: AxiolidContextHandle) -> AxiolidStatus {
    boundary(|| {
        if contexts().remove(context.0).is_some() {
            AxiolidStatus::Ok
        } else {
            AxiolidStatus::InvalidHandle
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn panic_is_contained() {
        assert_eq!(boundary(|| panic!("contained")), AxiolidStatus::Panic);
    }
}

/// Globally unique opaque mesh token owned by an Axiolid context.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AxiolidMeshHandle(pub u64);

impl AxiolidMeshHandle {
    pub const INVALID: Self = Self(0);
}

/// Operation identifier in a capability descriptor.
#[non_exhaustive]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidOperation {
    CurveEvaluation = 1,
    SurfaceEvaluation = 2,
    ProfileTriangulation = 3,
    Sweep = 4,
    Tessellation = 5,
    MeshBoolean = 6,
    MeshPlaneSection = 7,
    SpatialQuery = 8,
    Measurement = 9,
    Healing = 10,
    GraphCompilation = 11,
    MeshHiddenLines = 12,
    MeshVisibleFaces = 13,
    Unknown = 255,
}

/// Representation identifier in a capability descriptor.
#[non_exhaustive]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidRepresentation {
    Scalar = 1,
    Linear = 2,
    Profile2d = 3,
    AnalyticCurve = 4,
    AnalyticSurface = 5,
    Topology = 6,
    ExactBrep = 7,
    TriangleMesh = 8,
    MeshHealth = 9,
    Measurements = 10,
    RayHit = 11,
    ModelGraph = 12,
    SampledField = 13,
    Unknown = 255,
}

/// Exactness promise in a capability descriptor.
#[non_exhaustive]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidExactness {
    Exact = 1,
    ToleranceBounded = 2,
}

/// Stable, pointer-free capability description. String lengths exclude padding.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct AxiolidCapabilityDescriptor {
    pub operation: AxiolidOperation,
    pub representation: AxiolidRepresentation,
    pub exactness: AxiolidExactness,
    pub id: [u8; 64],
    pub provider: [u8; 64],
    pub required_feature: [u8; 32],
    pub id_len: u8,
    pub provider_len: u8,
    pub required_feature_len: u8,
    pub deterministic: u8,
}

impl Default for AxiolidCapabilityDescriptor {
    fn default() -> Self {
        Self {
            id: [0; 64],
            id_len: 0,
            provider: [0; 64],
            provider_len: 0,
            required_feature: [0; 32],
            required_feature_len: 0,
            operation: AxiolidOperation::Unknown,
            representation: AxiolidRepresentation::Unknown,
            exactness: AxiolidExactness::ToleranceBounded,
            deterministic: 0,
        }
    }
}

fn copy_text<const N: usize>(text: &str, destination: &mut [u8; N]) -> u8 {
    let count = text.len().min(N);
    destination[..count].copy_from_slice(&text.as_bytes()[..count]);
    count as u8
}

fn operation(value: axiolid::contracts::Operation) -> AxiolidOperation {
    use axiolid::contracts::Operation;
    match value {
        Operation::CurveEvaluation => AxiolidOperation::CurveEvaluation,
        Operation::SurfaceEvaluation => AxiolidOperation::SurfaceEvaluation,
        Operation::ProfileTriangulation => AxiolidOperation::ProfileTriangulation,
        Operation::Sweep => AxiolidOperation::Sweep,
        Operation::Tessellation => AxiolidOperation::Tessellation,
        Operation::MeshBoolean => AxiolidOperation::MeshBoolean,
        Operation::MeshPlaneSection => AxiolidOperation::MeshPlaneSection,
        Operation::SpatialQuery => AxiolidOperation::SpatialQuery,
        Operation::Measurement => AxiolidOperation::Measurement,
        Operation::Healing => AxiolidOperation::Healing,
        Operation::GraphCompilation => AxiolidOperation::GraphCompilation,
        Operation::MeshHiddenLines => AxiolidOperation::MeshHiddenLines,
        Operation::MeshVisibleFaces => AxiolidOperation::MeshVisibleFaces,
        _ => AxiolidOperation::Unknown,
    }
}

fn representation(value: axiolid::contracts::Representation) -> AxiolidRepresentation {
    use axiolid::contracts::Representation;
    match value {
        Representation::Scalar => AxiolidRepresentation::Scalar,
        Representation::Linear => AxiolidRepresentation::Linear,
        Representation::Profile2d => AxiolidRepresentation::Profile2d,
        Representation::AnalyticCurve => AxiolidRepresentation::AnalyticCurve,
        Representation::AnalyticSurface => AxiolidRepresentation::AnalyticSurface,
        Representation::Topology => AxiolidRepresentation::Topology,
        Representation::ExactBrep => AxiolidRepresentation::ExactBrep,
        Representation::TriangleMesh => AxiolidRepresentation::TriangleMesh,
        Representation::MeshHealth => AxiolidRepresentation::MeshHealth,
        Representation::Measurements => AxiolidRepresentation::Measurements,
        Representation::RayHit => AxiolidRepresentation::RayHit,
        Representation::ModelGraph => AxiolidRepresentation::ModelGraph,
        Representation::SampledField => AxiolidRepresentation::SampledField,
    }
}

fn exactness(value: axiolid::contracts::Exactness) -> AxiolidExactness {
    match value {
        axiolid::contracts::Exactness::Exact => AxiolidExactness::Exact,
        axiolid::contracts::Exactness::ToleranceBounded => AxiolidExactness::ToleranceBounded,
    }
}

fn capability(value: &axiolid::contracts::CapabilityDescriptor) -> AxiolidCapabilityDescriptor {
    let mut result = AxiolidCapabilityDescriptor::default();
    result.id_len = copy_text(value.id.as_str(), &mut result.id);
    result.provider_len = copy_text(value.provider.as_str(), &mut result.provider);
    result.required_feature_len = copy_text(value.required_feature, &mut result.required_feature);
    result.operation = operation(value.operation);
    result.representation = representation(value.output);
    result.exactness = exactness(value.exactness);
    result.deterministic = u8::from(value.deterministic);
    result
}

fn c_capabilities(
    context: &Context,
) -> impl Iterator<Item = &axiolid::contracts::CapabilityDescriptor> {
    context
        .application
        .descriptor()
        .capabilities
        .iter()
        .filter(|descriptor| {
            matches!(
                descriptor.operation,
                axiolid::contracts::Operation::Healing
                    | axiolid::contracts::Operation::Measurement
                    | axiolid::contracts::Operation::MeshBoolean
                    | axiolid::contracts::Operation::MeshHiddenLines
                    | axiolid::contracts::Operation::MeshVisibleFaces
                    | axiolid::contracts::Operation::Sweep
            )
        })
}

/// Query the number of capabilities callable through this ABI context.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_capability_count(
    context: AxiolidContextHandle,
    out_count: *mut usize,
) -> AxiolidStatus {
    boundary(|| {
        if out_count.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let count = c_capabilities(&context).count();
        // SAFETY: null was rejected; the caller contract requires one writable usize.
        unsafe { out_count.write(count) };
        AxiolidStatus::Ok
    })
}

/// Copy one capability descriptor into caller-owned storage.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_capability_get(
    context: AxiolidContextHandle,
    index: usize,
    out_descriptor: *mut AxiolidCapabilityDescriptor,
) -> AxiolidStatus {
    boundary(|| {
        if out_descriptor.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(value) = c_capabilities(&context).nth(index) else {
            return AxiolidStatus::InvalidArgument;
        };
        // SAFETY: null was rejected; the caller contract requires one writable descriptor.
        unsafe { out_descriptor.write(capability(value)) };
        AxiolidStatus::Ok
    })
}

fn tolerance(value: AxiolidTolerance) -> Result<axiolid::core::Tolerance, AxiolidStatus> {
    axiolid::core::Tolerance::new(value.linear, value.angular)
        .map_err(|_| AxiolidStatus::InvalidArgument)
}

fn record_error(
    mut context: MutexGuard<'_, Context>,
    status: AxiolidStatus,
    operation: AxiolidOperation,
    tolerance: AxiolidTolerance,
    provider: impl Into<String>,
    message: impl Into<String>,
) -> AxiolidStatus {
    context.last_error = Some(ErrorRecord {
        status,
        operation,
        tolerance,
        provider: provider.into(),
        message: message.into(),
        mesh: None,
    });
    status
}

fn record_application_error(
    context: MutexGuard<'_, Context>,
    error: axiolid::application::ApplicationError,
) -> AxiolidStatus {
    let operation = operation(error.context.operation);
    let tolerance = AxiolidTolerance {
        linear: error.context.tolerance.linear(),
        angular: error.context.tolerance.angular(),
    };
    record_error(
        context,
        AxiolidStatus::OperationFailed,
        operation,
        tolerance,
        error.context.provider.as_str(),
        error.to_string(),
    )
}

/// Structured portion of a context-owned error record.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct AxiolidErrorInfo {
    pub status: AxiolidStatus,
    pub operation: AxiolidOperation,
    pub tolerance: AxiolidTolerance,
    pub provider: [u8; 64],
    pub message_len: usize,
    pub provider_len: u8,
    /// Must be zero; reserves layout space without exposing padding bytes.
    pub reserved: [u8; 7],
}

impl Default for AxiolidErrorInfo {
    fn default() -> Self {
        Self {
            status: AxiolidStatus::NoError,
            operation: AxiolidOperation::Unknown,
            tolerance: AxiolidTolerance {
                linear: 0.0,
                angular: 0.0,
            },
            provider: [0; 64],
            message_len: 0,
            provider_len: 0,
            reserved: [0; 7],
        }
    }
}

mod diagnostics;
mod drawing;
mod faces;
mod mesh;
mod operations;

pub use diagnostics::*;
pub use drawing::*;
pub use faces::*;
pub use mesh::*;
pub use operations::*;
