//! Adapter from an API-specific GPU executor to graph compilation.
//!
//! Before anything is submitted, the adapter checks what it can check
//! without the device: device preference, f32/f64 policy, and result
//! residency (a device cannot deliver into another device's memory, and an
//! unrecognized future residency is refused rather than assumed), then that
//! every root belongs to the graph, then the executor's own
//! option-validation hook. After the executor returns, it checks one result
//! per root. Caller faults are `Unsupported` or `InvalidInput`; executor
//! output that breaks the seam's contract is `BackendContractViolation`.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, DevicePreference, ExecutionOptions, ExecutionTarget,
    GeomError, GeomResult, Operation, Precision, Residency,
};
use axiolid_mesh::TriMesh;
use axiolid_mesh_compile_contract::MeshCompiler;
use axiolid_model::{GeometryGraph, NodeId};

use crate::{GpuDeviceDescriptor, GpuGraphExecutor};

/// Graph-compiler provider backed by one initialized GPU executor.
#[derive(Debug)]
pub struct GpuCompiler<E> {
    executor: E,
}

impl<E> GpuCompiler<E> {
    /// Wrap one initialized concrete executor.
    pub const fn new(executor: E) -> Self {
        Self { executor }
    }

    /// Underlying device facts.
    pub fn device(&self) -> &GpuDeviceDescriptor
    where
        E: GpuGraphExecutor,
    {
        self.executor.device()
    }

    /// Borrow the API-specific executor for advanced operations.
    pub const fn executor(&self) -> &E {
        &self.executor
    }

    fn validate_options(&self, options: &ExecutionOptions) -> GeomResult<()>
    where
        E: GpuGraphExecutor,
    {
        let device = self.executor.device();
        let compatible_device = match options.device() {
            DevicePreference::Auto | DevicePreference::Gpu => true,
            DevicePreference::Backend(required) => required == device.id,
            DevicePreference::Cpu => false,
        };
        if !compatible_device || (options.precision() == Precision::F64 && !device.features.float64)
        {
            return Err(GeomError::Unsupported {
                backend: device.id,
                operation: Operation::GraphCompilation,
            });
        }
        // Residency is part of the plan, not an afterthought: a device without
        // unified memory cannot serve a request whose results must stay in
        // another device's memory, and saying so here beats discovering it
        // after the upload.
        let output = options.residency().output();
        let deliverable = match output {
            Residency::Host => true,
            Residency::Device(owner) | Residency::Unified(owner) => owner == device.id,
            // `Residency` is non-exhaustive; an unrecognized future location is
            // refused rather than optimistically assumed deliverable.
            _ => false,
        };
        if !deliverable {
            return Err(GeomError::Unsupported {
                backend: device.id,
                operation: Operation::GraphCompilation,
            });
        }
        self.executor.validate_options(options)
    }

    fn validate_roots(graph: &GeometryGraph, roots: &[NodeId]) -> GeomResult<()> {
        if let Some(root) = roots.iter().find(|root| graph.get(**root).is_none()) {
            return Err(GeomError::InvalidInput(format!(
                "graph compilation root {root} does not belong to the supplied graph"
            )));
        }
        Ok(())
    }

    fn validate_result_count(
        backend: BackendId,
        root_count: usize,
        result_count: usize,
    ) -> GeomResult<()> {
        if result_count != root_count {
            return Err(GeomError::BackendContractViolation {
                backend,
                detail: format!("returned {result_count} results for {root_count} roots"),
            });
        }
        Ok(())
    }
}

impl<E: GpuGraphExecutor> Backend for GpuCompiler<E> {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(self.executor.device().id, ExecutionTarget::Gpu)
    }
}

impl<E: GpuGraphExecutor> MeshCompiler for GpuCompiler<E> {
    fn compile_mesh(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<TriMesh> {
        self.validate_options(options)?;
        Self::validate_roots(graph, &[root])?;
        let results = self.executor.compile_mesh_batch(graph, &[root], options)?;
        Self::validate_result_count(self.executor.device().id, 1, results.len())?;
        results
            .into_iter()
            .next()
            .ok_or_else(|| GeomError::BackendContractViolation {
                backend: self.executor.device().id,
                detail: "returned no result for one root".to_owned(),
            })
    }

    /// Overriding the `_into` seam keeps *both* batch call shapes on the
    /// single-dispatch GPU path; overriding only `compile_mesh_batch` would leave
    /// `compile_mesh_batch_into` silently falling back to one submission per root.
    fn compile_mesh_batch_into(
        &self,
        graph: &GeometryGraph,
        roots: &[NodeId],
        options: &ExecutionOptions,
        destination: &mut Vec<TriMesh>,
    ) -> GeomResult<()> {
        self.validate_options(options)?;
        Self::validate_roots(graph, roots)?;
        let results = self.executor.compile_mesh_batch(graph, roots, options)?;
        Self::validate_result_count(self.executor.device().id, roots.len(), results.len())?;
        destination.reserve(results.len());
        destination.extend(results);
        Ok(())
    }
}
