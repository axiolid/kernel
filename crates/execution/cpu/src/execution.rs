//! Runtime CPU execution context. Operation providers may compose this type.

use core::fmt;
use std::num::NonZeroUsize;
#[cfg(feature = "parallel")]
use std::sync::{Arc, OnceLock};

use crate::{CpuConfigError, CpuFeatures, CpuInstructionSet};

/// CPU execution context with runtime ISA detection and, with `parallel`, a
/// local pool of exactly [`Self::thread_count`] workers.
///
/// This type deliberately implements no geometry operation trait. A provider
/// only implements `MeshBoolean`, `MeshCompiler`, or another capability
/// after the algorithm exists and is verified.
#[derive(Clone)]
pub struct CpuExecution {
    instruction_set: CpuInstructionSet,
    features: CpuFeatures,
    threads: NonZeroUsize,
    /// The local pool, shared by clones. Started with the context by
    /// [`crate::CpuExecutionBuilder`]; for [`Self::portable`] and
    /// [`Self::detect`], on first use.
    #[cfg(feature = "parallel")]
    pool: Arc<OnceLock<rayon::ThreadPool>>,
}

impl fmt::Debug for CpuExecution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CpuExecution")
            .field("instruction_set", &self.instruction_set)
            .field("features", &self.features)
            .field("threads", &self.threads)
            .finish_non_exhaustive()
    }
}

impl CpuExecution {
    /// Always-available portable execution context.
    pub fn portable() -> Self {
        Self {
            instruction_set: CpuInstructionSet::Portable,
            features: CpuFeatures::detect(),
            threads: NonZeroUsize::MIN,
            #[cfg(feature = "parallel")]
            pool: Arc::default(),
        }
    }

    /// Best instruction path compiled into this binary and available at runtime.
    pub fn detect() -> Self {
        let features = CpuFeatures::detect();
        #[cfg(feature = "simd")]
        let instruction_set = features.best();
        #[cfg(not(feature = "simd"))]
        let instruction_set = CpuInstructionSet::Portable;
        Self {
            instruction_set,
            features,
            threads: NonZeroUsize::MIN,
            #[cfg(feature = "parallel")]
            pool: Arc::default(),
        }
    }

    pub(crate) fn from_configuration(
        instruction_set: CpuInstructionSet,
        features: CpuFeatures,
        threads: NonZeroUsize,
    ) -> Result<Self, CpuConfigError> {
        #[cfg(not(feature = "parallel"))]
        if threads.get() > 1 {
            return Err(CpuConfigError::ParallelFeatureDisabled);
        }

        // One worker gets a pool too: without one, `install` would run on
        // whatever pool the caller is in, and a provider sizing its scratch
        // or its work by `rayon::current_num_threads` would see that width
        // instead of the one configured (#259).
        #[cfg(feature = "parallel")]
        let pool =
            Arc::new(OnceLock::from(local_pool(threads).map_err(|error| {
                CpuConfigError::ThreadPool(error.to_string())
            })?));

        Ok(Self {
            instruction_set,
            features,
            threads,
            #[cfg(feature = "parallel")]
            pool,
        })
    }

    /// Runtime-selected instruction set.
    pub const fn instruction_set(&self) -> CpuInstructionSet {
        self.instruction_set
    }

    /// Detected machine features.
    pub const fn features(&self) -> CpuFeatures {
        self.features
    }

    /// Configured worker bound.
    pub const fn thread_count(&self) -> NonZeroUsize {
        self.threads
    }

    /// Execute one closure inside this context's local Rayon pool, which has
    /// exactly [`Self::thread_count`] workers, one included: inside it
    /// `rayon::current_num_threads()` is that count, whatever pool the
    /// caller runs in.
    ///
    /// # Panics
    ///
    /// For a context from [`Self::portable`] or [`Self::detect`], the pool
    /// is started on first use, and this panics if the operating system
    /// refuses the thread, as a first use of Rayon's global pool does. A
    /// context from [`crate::CpuExecutionBuilder`] started its pool when it
    /// was built and reported a failure there.
    #[cfg(feature = "parallel")]
    pub fn install<R: Send>(&self, operation: impl FnOnce() -> R + Send) -> R {
        self.pool
            .get_or_init(|| local_pool(self.threads).expect("start the local pool"))
            .install(operation)
    }
}

/// A pool of exactly `threads` workers.
#[cfg(feature = "parallel")]
fn local_pool(threads: NonZeroUsize) -> Result<rayon::ThreadPool, rayon::ThreadPoolBuildError> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.get())
        .thread_name(|index| format!("axiolid-cpu-{index}"))
        .build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CpuExecutionBuilder, InstructionPolicy};

    fn assert_runtime_traits<T: fmt::Debug + Clone + Send + Sync>() {}

    #[test]
    fn portable_context_is_constructible_but_claims_no_operation_trait() {
        assert_runtime_traits::<CpuExecution>();
        let backend = CpuExecutionBuilder::new()
            .instruction_policy(InstructionPolicy::Portable)
            .build()
            .expect("portable context");
        assert_eq!(backend.instruction_set(), CpuInstructionSet::Portable);
    }

    #[cfg(not(feature = "simd"))]
    #[test]
    fn auto_selection_stays_portable_without_the_simd_feature() {
        assert_eq!(
            CpuExecution::detect().instruction_set(),
            CpuInstructionSet::Portable
        );
    }

    #[cfg(feature = "simd")]
    #[test]
    fn auto_selection_uses_the_best_runtime_supported_instruction_set() {
        let backend = CpuExecution::detect();
        assert_eq!(backend.instruction_set(), backend.features().best());
        assert!(backend.features().supports(backend.instruction_set()));
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn configured_worker_bound_controls_the_local_pool() {
        let backend = CpuExecutionBuilder::new()
            .threads(NonZeroUsize::new(2).unwrap())
            .build()
            .expect("two-worker context");
        assert_eq!(backend.thread_count().get(), 2);
        assert_eq!(backend.install(rayon::current_num_threads), 2);
    }

    /// One worker is a pool of one, not the caller's pool (#259): run from
    /// inside a wider pool, falling through to the caller would show.
    #[cfg(feature = "parallel")]
    #[test]
    fn a_one_thread_context_runs_on_a_pool_of_one() {
        let caller = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .expect("caller pool");
        let built = CpuExecutionBuilder::new()
            .threads(NonZeroUsize::MIN)
            .build()
            .expect("one-worker context");
        for context in [built, CpuExecution::portable(), CpuExecution::detect()] {
            assert_eq!(context.thread_count().get(), 1);
            let width = caller.install(|| context.install(rayon::current_num_threads));
            assert_eq!(width, 1, "{context:?} ran on the caller's pool");
            // Clones share the pool rather than each starting their own.
            let worker = |c: &CpuExecution| c.install(|| std::thread::current().id());
            assert_eq!(worker(&context), worker(&context.clone()));
        }
    }

    #[cfg(not(feature = "parallel"))]
    #[test]
    fn multiple_workers_require_the_parallel_feature() {
        let error = CpuExecutionBuilder::new()
            .threads(NonZeroUsize::new(2).unwrap())
            .build()
            .unwrap_err();
        assert_eq!(error, CpuConfigError::ParallelFeatureDisabled);
    }
}
