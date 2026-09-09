//! Provider registration, ordering, fallback, and budget policy for mesh booleans.

use std::sync::Arc;

use axiolid_contracts::{BackendId, ExecutionOptions, GeomError, GeomResult, Operation};

use crate::device::matches_device;
use axiolid_core::BooleanOperator;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::{conformance, BooleanOutcome, MeshBoolean};
use axiolid_mesh_contracts::SolidRequirements;

#[derive(Debug, Clone)]
struct RegisteredBoolean {
    priority: i32,
    provider: Arc<dyn MeshBoolean>,
}

/// Ordered executable providers for one narrow operation.
///
/// Fallback happens only for `Unsupported` or `Unavailable`; numerical and data
/// failures are returned immediately rather than hidden by another algorithm.
#[derive(Debug, Clone, Default)]
pub struct MeshBooleanRegistry {
    providers: Vec<RegisteredBoolean>,
}

impl MeshBooleanRegistry {
    /// Empty registry.
    pub const fn new() -> Self {
        Self {
            providers: Vec::new(),
        }
    }

    /// Register an implementation. Higher priorities run first.
    pub fn register<B>(&mut self, priority: i32, provider: B)
    where
        B: MeshBoolean + 'static,
    {
        self.register_arc(priority, Arc::new(provider));
    }

    /// Register only if the provider passes the shared conformance suite.
    ///
    /// ADR 0017 §6 makes conformance a *precondition* of registration rather
    /// than a test someone might remember to run. A provider that violates the
    /// contract is rejected here, with the failing report, instead of being
    /// discovered later by a caller receiving wrong geometry.
    ///
    /// [`Self::register`] remains available for tests and for deliberately
    /// partial providers; this is the door production code should use.
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
        B: MeshBoolean + 'static,
    {
        let report = conformance::run(&provider);
        if !report.is_conformant() {
            return Err(Box::new(report));
        }
        self.register_arc(priority, Arc::new(provider));
        Ok(())
    }

    /// Register a shared trait object.
    pub fn register_arc(&mut self, priority: i32, provider: Arc<dyn MeshBoolean>) {
        self.providers
            .push(RegisteredBoolean { priority, provider });
        self.providers
            .sort_by_key(|entry| std::cmp::Reverse(entry.priority));
    }

    /// Registered providers in dispatch order.
    pub fn providers(&self) -> impl Iterator<Item = &dyn MeshBoolean> {
        self.providers.iter().map(|entry| entry.provider.as_ref())
    }

    fn dispatch(
        &self,
        options: &ExecutionOptions,
        elements: usize,
        execute: impl Fn(&dyn MeshBoolean) -> GeomResult<BooleanOutcome>,
    ) -> GeomResult<BooleanOutcome> {
        let mut last_retryable = None;
        let mut over_budget = None;
        for entry in &self.providers {
            let descriptor = entry.provider.descriptor();
            if !matches_device(options.device(), descriptor.id, descriptor.target) {
                continue;
            }
            // Budget is checked before dispatch, not after: a provider that
            // cannot fit the caller's memory bound must never get the chance to
            // allocate. Treated as retryable so a leaner provider can still run.
            if !entry
                .provider
                .scratch_requirement()
                .fits_budget(options, elements)
            {
                over_budget = Some(GeomError::BudgetExceeded { resource: "memory" });
                continue;
            }
            match execute(entry.provider.as_ref()) {
                Ok(outcome) => return Ok(outcome),
                Err(error @ (GeomError::Unsupported { .. } | GeomError::Unavailable { .. })) => {
                    last_retryable = Some(error);
                }
                Err(error) => return Err(error),
            }
        }
        Err(last_retryable
            .or(over_budget)
            .unwrap_or(GeomError::Unsupported {
                backend: BackendId::new("mesh-boolean-registry"),
                operation: Operation::MeshBoolean,
            }))
    }

    /// Execute according to device policy with narrow fallback semantics.
    pub fn boolean(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operation: BooleanOperator,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        // Admissibility is contract-level and checked before any provider sees
        // the operands, so dispatch cannot change which inputs are legal.
        SolidRequirements::Oriented.validate_operands(subject, &[tool])?;
        self.dispatch(options, subject.triangle_count(), |provider| {
            provider.boolean(subject, tool, operation, options)
        })
    }

    /// Subtract many tools through one provider dispatch.
    pub fn subtract_many(
        &self,
        subject: &TriMesh,
        tools: &[TriMesh],
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        let borrowed: Vec<&TriMesh> = tools.iter().collect();
        SolidRequirements::Oriented.validate_operands(subject, &borrowed)?;
        let elements =
            subject.triangle_count() + tools.iter().map(TriMesh::triangle_count).sum::<usize>();
        self.dispatch(options, elements, |provider| {
            provider.subtract_many(subject, tools, options)
        })
    }

    /// Union many solids through one provider dispatch.
    ///
    /// An empty batch is answered without consulting a provider: the union
    /// of nothing is nothing, and there is no operand to validate or size a
    /// budget against.
    pub fn union_many(
        &self,
        solids: &[TriMesh],
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        let Some((first, rest)) = solids.split_first() else {
            return Ok(BooleanOutcome::new(TriMesh::default(), Default::default()));
        };
        // Union has no privileged operand, but validation is expressed as
        // subject-plus-tools. Naming the first solid the subject validates
        // exactly the same set, in the same way the provider reports it.
        let borrowed: Vec<&TriMesh> = rest.iter().collect();
        SolidRequirements::Oriented.validate_operands(first, &borrowed)?;
        let elements = solids.iter().map(TriMesh::triangle_count).sum::<usize>();
        self.dispatch(options, elements, |provider| {
            provider.union_many(solids, options)
        })
    }
}

#[cfg(test)]
mod tests {
    use axiolid_contracts::{DevicePreference, ExecutionTarget};
    /// Outward-oriented unit cube: the minimal admissible operand.
    ///
    /// Dispatch tests need a mesh that passes contract validation, because
    /// validation now runs before any provider is consulted.
    fn admissible_cube() -> TriMesh {
        let positions = vec![
            [0.0, 0.0, 0.0].into(),
            [1.0, 0.0, 0.0].into(),
            [1.0, 1.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            [0.0, 0.0, 1.0].into(),
            [1.0, 0.0, 1.0].into(),
            [1.0, 1.0, 1.0].into(),
            [0.0, 1.0, 1.0].into(),
        ];
        let indices = vec![
            0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7,
            6, 3, 0, 4, 3, 4, 7,
        ];
        TriMesh::new(positions, indices)
    }

    use std::sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    };

    use axiolid_contracts::{Backend, BackendDescriptor};
    use axiolid_core::Tolerance;
    use axiolid_mesh_boolean_contract::BooleanEvidence;

    use super::*;

    #[derive(Debug)]
    struct EchoBoolean {
        id: BackendId,
        target: ExecutionTarget,
    }

    impl Backend for EchoBoolean {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor {
                id: self.id,
                target: self.target,
            }
        }
    }

    impl MeshBoolean for EchoBoolean {
        fn boolean(
            &self,
            subject: &TriMesh,
            _tool: &TriMesh,
            _operation: BooleanOperator,
            _options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            Ok(BooleanOutcome::new(
                subject.clone(),
                BooleanEvidence::default(),
            ))
        }
    }

    #[derive(Debug, Clone, Copy)]
    enum ProbeResult {
        Success,
        Unsupported,
        Unavailable,
        Invalid,
    }

    #[derive(Debug)]
    struct ProbeBoolean {
        id: BackendId,
        target: ExecutionTarget,
        result: ProbeResult,
        calls: Arc<AtomicUsize>,
    }

    impl Backend for ProbeBoolean {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(self.id, self.target)
        }
    }

    impl MeshBoolean for ProbeBoolean {
        fn boolean(
            &self,
            subject: &TriMesh,
            _tool: &TriMesh,
            _operation: BooleanOperator,
            _options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            match self.result {
                ProbeResult::Success => Ok(BooleanOutcome::new(
                    subject.clone(),
                    BooleanEvidence::default(),
                )),
                ProbeResult::Unsupported => Err(GeomError::Unsupported {
                    backend: self.id,
                    operation: Operation::MeshBoolean,
                }),
                ProbeResult::Unavailable => Err(GeomError::Unavailable {
                    backend: self.id,
                    reason: "probe unavailable".to_owned(),
                }),
                ProbeResult::Invalid => {
                    Err(GeomError::InvalidInput("probe rejected input".to_owned()))
                }
            }
        }
    }

    #[derive(Debug)]
    struct BatchBoolean {
        calls: Arc<AtomicUsize>,
    }

    impl Backend for BatchBoolean {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("batch"), ExecutionTarget::OptimizedCpu)
        }
    }

    impl MeshBoolean for BatchBoolean {
        fn boolean(
            &self,
            _subject: &TriMesh,
            _tool: &TriMesh,
            _operation: BooleanOperator,
            _options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            Err(GeomError::InvalidInput(
                "batch provider must use its batch override".to_owned(),
            ))
        }

        fn subtract_many(
            &self,
            subject: &TriMesh,
            _tools: &[TriMesh],
            _options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            self.calls.fetch_add(1, Ordering::Relaxed);
            Ok(BooleanOutcome::new(
                subject.clone(),
                BooleanEvidence::default(),
            ))
        }
    }

    #[test]
    fn registry_stores_executable_traits_not_capability_flags() {
        let mut registry = MeshBooleanRegistry::new();
        registry.register(
            10,
            EchoBoolean {
                id: BackendId::new("echo"),
                target: ExecutionTarget::PortableCpu,
            },
        );
        let options = ExecutionOptions::new(Tolerance::METRE);
        let mesh = admissible_cube();
        assert_eq!(
            registry
                .boolean(&mesh, &mesh, BooleanOperator::Difference, &options)
                .expect("registered provider executes"),
            BooleanOutcome::new(mesh, BooleanEvidence::default())
        );
    }

    #[test]
    fn registry_dispatches_batch_subtraction_to_the_provider_override() {
        let calls = Arc::new(AtomicUsize::new(0));
        let mut registry = MeshBooleanRegistry::new();
        registry.register(
            10,
            BatchBoolean {
                calls: calls.clone(),
            },
        );
        let mesh = admissible_cube();
        let options = ExecutionOptions::new(Tolerance::METRE);

        assert_eq!(
            registry
                .subtract_many(&mesh, &[mesh.clone(), mesh.clone()], &options)
                .expect("batch provider executes"),
            BooleanOutcome::new(mesh, BooleanEvidence::default())
        );
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn registry_falls_back_only_for_retryable_errors_and_honors_device_policy() {
        let high_calls = Arc::new(AtomicUsize::new(0));
        let unsupported_calls = Arc::new(AtomicUsize::new(0));
        let low_calls = Arc::new(AtomicUsize::new(0));
        let mut registry = MeshBooleanRegistry::new();
        registry.register(
            100,
            ProbeBoolean {
                id: BackendId::new("unavailable-gpu"),
                target: ExecutionTarget::Gpu,
                result: ProbeResult::Unavailable,
                calls: high_calls.clone(),
            },
        );
        registry.register(
            50,
            ProbeBoolean {
                id: BackendId::new("unsupported-cpu"),
                target: ExecutionTarget::OptimizedCpu,
                result: ProbeResult::Unsupported,
                calls: unsupported_calls.clone(),
            },
        );
        registry.register(
            10,
            ProbeBoolean {
                id: BackendId::new("portable-fallback"),
                target: ExecutionTarget::PortableCpu,
                result: ProbeResult::Success,
                calls: low_calls.clone(),
            },
        );
        let mesh = admissible_cube();
        let auto = ExecutionOptions::new(Tolerance::METRE);
        assert!(registry
            .boolean(&mesh, &mesh, BooleanOperator::Union, &auto)
            .is_ok());
        assert_eq!(high_calls.load(Ordering::Relaxed), 1);
        assert_eq!(unsupported_calls.load(Ordering::Relaxed), 1);
        assert_eq!(low_calls.load(Ordering::Relaxed), 1);

        let cpu = auto.clone().with_device(DevicePreference::Cpu);
        assert!(registry
            .boolean(&mesh, &mesh, BooleanOperator::Union, &cpu)
            .is_ok());
        assert_eq!(high_calls.load(Ordering::Relaxed), 1);
        assert_eq!(unsupported_calls.load(Ordering::Relaxed), 2);
        assert_eq!(low_calls.load(Ordering::Relaxed), 2);

        let invalid_calls = Arc::new(AtomicUsize::new(0));
        let skipped_calls = Arc::new(AtomicUsize::new(0));
        let mut fail_fast = MeshBooleanRegistry::new();
        fail_fast.register(
            100,
            ProbeBoolean {
                id: BackendId::new("invalid-input"),
                target: ExecutionTarget::PortableCpu,
                result: ProbeResult::Invalid,
                calls: invalid_calls.clone(),
            },
        );
        fail_fast.register(
            10,
            ProbeBoolean {
                id: BackendId::new("must-not-run"),
                target: ExecutionTarget::PortableCpu,
                result: ProbeResult::Success,
                calls: skipped_calls.clone(),
            },
        );
        assert!(matches!(
            fail_fast.boolean(&mesh, &mesh, BooleanOperator::Union, &auto),
            Err(GeomError::InvalidInput(_))
        ));
        assert_eq!(invalid_calls.load(Ordering::Relaxed), 1);
        assert_eq!(skipped_calls.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn empty_registry_returns_structured_unsupported_error() {
        let registry = MeshBooleanRegistry::new();
        let mesh = admissible_cube();
        let error = registry
            .boolean(
                &mesh,
                &mesh,
                BooleanOperator::Union,
                &ExecutionOptions::new(Tolerance::METRE),
            )
            .unwrap_err();
        assert!(matches!(
            error,
            GeomError::Unsupported {
                backend,
                operation: Operation::MeshBoolean,
            } if backend == BackendId::new("mesh-boolean-registry")
        ));
    }
}
