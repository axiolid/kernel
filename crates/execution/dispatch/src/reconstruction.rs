//! Provider registration, ordering, fallback, and budget policy for
//! pointcloud reconstruction.

use std::sync::Arc;

use axiolid_contracts::{BackendId, ExecutionOptions, GeomError, GeomResult, Operation};

use crate::device::matches_device;
use axiolid_pointcloud::PointCloud;
use axiolid_pointcloud_reconstruction_contract::{
    conformance, PointcloudReconstruction, Reconstruction, ReconstructionRequest,
};

#[derive(Clone)]
struct RegisteredReconstruction {
    priority: i32,
    provider: Arc<dyn PointcloudReconstruction>,
}

impl core::fmt::Debug for RegisteredReconstruction {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("RegisteredReconstruction")
            .field("priority", &self.priority)
            .field("backend", &self.provider.descriptor().id)
            .finish()
    }
}

/// Ordered executable providers for pointcloud reconstruction.
///
/// Fallback happens only for `Unsupported` or `Unavailable`. A *refusal* is
/// not a fallback trigger: when a provider says the data cannot support the
/// request, that is an answer about the data, and asking a second provider
/// the same question would only find one willing to guess.
#[derive(Debug, Clone, Default)]
pub struct PointcloudReconstructionRegistry {
    providers: Vec<RegisteredReconstruction>,
}

impl PointcloudReconstructionRegistry {
    /// Empty registry.
    pub const fn new() -> Self {
        Self {
            providers: Vec::new(),
        }
    }

    /// Register an implementation. Higher priorities run first.
    pub fn register<B>(&mut self, priority: i32, provider: B)
    where
        B: PointcloudReconstruction + 'static,
    {
        self.register_arc(priority, Arc::new(provider));
    }

    /// Register only if the provider passes the shared conformance suite.
    ///
    /// Conformance is a *precondition* of registration rather than a test
    /// someone might remember to run: a provider that violates the contract
    /// is rejected here, with the failing report, instead of being
    /// discovered later by a caller receiving an unexplained empty mesh.
    ///
    /// # Errors
    ///
    /// Returns the report when the provider violates any obligation.
    pub fn register_conformant<B>(
        &mut self,
        priority: i32,
        provider: B,
    ) -> Result<(), Box<conformance::ConformanceReport>>
    where
        B: PointcloudReconstruction + 'static,
    {
        let report = conformance::run(&provider);
        if !report.is_conformant() {
            return Err(Box::new(report));
        }
        self.register_arc(priority, Arc::new(provider));
        Ok(())
    }

    /// Register a shared trait object.
    pub fn register_arc(&mut self, priority: i32, provider: Arc<dyn PointcloudReconstruction>) {
        self.providers
            .push(RegisteredReconstruction { priority, provider });
        self.providers
            .sort_by_key(|entry| std::cmp::Reverse(entry.priority));
    }

    /// Registered providers in dispatch order.
    pub fn providers(&self) -> impl Iterator<Item = &dyn PointcloudReconstruction> {
        self.providers.iter().map(|entry| entry.provider.as_ref())
    }

    /// Whether any provider is registered.
    pub fn is_empty(&self) -> bool {
        self.providers.is_empty()
    }

    /// Reconstruct through the first provider able to take the work.
    ///
    /// # Errors
    ///
    /// Returns `Unsupported` when no registered provider can run at all, and
    /// `BudgetExceeded` when every candidate was excluded by the caller's
    /// memory bound.
    pub fn reconstruct(
        &self,
        cloud: &PointCloud,
        request: &ReconstructionRequest,
        options: &ExecutionOptions,
    ) -> GeomResult<Reconstruction> {
        let mut last_retryable = None;
        let mut over_budget = None;

        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if !matches_device(options.device(), descriptor.id, descriptor.target) {
                continue;
            }
            // Budget is checked before dispatch, not after: a provider that
            // cannot fit the caller's memory bound must never get the chance
            // to allocate. Reconstruction is the operation most likely to
            // exhaust memory, so this matters more here than elsewhere.
            if !entry
                .provider
                .scratch_requirement()
                .fits_budget(options, cloud.len())
            {
                over_budget = Some(GeomError::BudgetExceeded { resource: "memory" });
                continue;
            }
            match entry.provider.reconstruct(cloud, request, options) {
                // A refusal is a real answer about the data. Returning it
                // immediately is what stops fallback from shopping for a
                // provider willing to fabricate a surface.
                Ok(result) => return Ok(result),
                Err(error @ (GeomError::Unsupported { .. } | GeomError::Unavailable { .. })) => {
                    last_retryable = Some(error);
                }
                Err(error) => return Err(error),
            }
        }

        Err(last_retryable
            .or(over_budget)
            .unwrap_or(GeomError::Unsupported {
                backend: BackendId::new("pointcloud-reconstruction-registry"),
                operation: Operation::PointcloudReconstruction,
            }))
    }
}
