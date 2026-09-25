//! Execution policy passed explicitly to every costly operation.

use axiolid_core::{Scalar, Tolerance};

use crate::cancel::CancellationToken;
use crate::{BackendId, GeomError, GeomResult, Precision};

/// Reproducibility requirement.
///
/// These are three genuinely different contracts, not degrees of one. A backend
/// that reduces in a deterministic order satisfies [`Self::Topological`] while
/// still failing [`Self::Bitwise`] against a differently-scheduled backend, so
/// collapsing them into one flag lets two backends both claim "deterministic"
/// and still disagree. Ordered weakest to strongest.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Determinism {
    /// No guarantee. The backend may reorder, re-associate, and reschedule.
    BestEffort,
    /// Same connectivity and same output ordering for the same inputs.
    ///
    /// Coordinates may differ within tolerance. This is what a clash result or
    /// a topology commit needs; it does not promise identical floats.
    Topological,
    /// Topological determinism plus values within the operation's stated
    /// numerical error bound.
    NumericallyBounded,
    /// Bit-for-bit identical output for the same inputs and options.
    ///
    /// The only contract that supports hashing a result or comparing artifacts
    /// across machines.
    Bitwise,
}

impl Determinism {
    /// Whether this guarantee is at least as strong as `required`.
    ///
    /// Strength is the declaration order, so a backend offering
    /// [`Self::Bitwise`] satisfies a [`Self::Topological`] request but never
    /// the reverse. Routing must use this instead of equality, or a stronger
    /// backend gets rejected for being too good.
    pub const fn satisfies(self, required: Self) -> bool {
        (self as u8) >= (required as u8)
    }
}

/// CPU scheduling preference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Parallelism {
    /// One worker.
    Serial,
    /// Backend chooses from available parallelism and workload size.
    Auto,
    /// Upper bound on worker count. Zero is rejected by the builder.
    Threads(usize),
}

/// Device selection preference. `Auto` is a policy request, not permission to
/// silently reduce precision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DevicePreference {
    /// Select from compatible registered backends.
    Auto,
    /// Require a CPU backend.
    Cpu,
    /// Require a GPU backend.
    Gpu,
    /// Require one named backend.
    Backend(BackendId),
}

/// Where an operation's data physically lives.
///
/// Device preference says *where to run*; residency says *where the bytes
/// already are*. Routing needs both: a GPU-resident batch is cheap to run on
/// the GPU and expensive to run on the CPU, and the reverse holds for a
/// host-resident one. Without this, a planner cannot see a transfer that
/// dominates the operation it is scheduling.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Residency {
    /// Host (CPU) memory.
    Host,
    /// Memory owned by the named device.
    Device(BackendId),
    /// Host-visible memory shared with the named device; no copy is needed.
    Unified(BackendId),
}

impl Residency {
    /// Whether `backend` can read this data without a host/device transfer.
    pub fn is_local_to(self, backend: BackendId) -> bool {
        match self {
            Self::Host => false,
            Self::Device(owner) | Self::Unified(owner) => owner == backend,
        }
    }

    /// Whether the host can read this data without a transfer.
    pub const fn is_host_readable(self) -> bool {
        matches!(self, Self::Host | Self::Unified(_))
    }
}

/// Where an operation's inputs live and where its outputs are wanted.
///
/// Kept as a pair because they genuinely differ: a GPU broad phase may consume
/// device-resident geometry and still have to deliver host-readable results.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DataResidency {
    input: Residency,
    output: Residency,
}

impl DataResidency {
    /// Both inputs and outputs in host memory. The portable default.
    pub const HOST: Self = Self {
        input: Residency::Host,
        output: Residency::Host,
    };

    /// Construct an explicit input/output residency plan.
    pub const fn new(input: Residency, output: Residency) -> Self {
        Self { input, output }
    }

    /// Where inputs currently live.
    pub const fn input(self) -> Residency {
        self.input
    }

    /// Where outputs are wanted.
    pub const fn output(self) -> Residency {
        self.output
    }

    /// Whether running on `backend` needs no host/device transfer either way.
    pub fn is_transfer_free_on(self, backend: BackendId) -> bool {
        self.input.is_local_to(backend) && self.output.is_local_to(backend)
    }
}

/// Scratch memory an operation needs beyond its inputs and outputs.
///
/// Declared up front so a caller can budget, pre-reserve, or refuse before any
/// work starts. `Unbounded` is deliberately representable and deliberately
/// unpleasant: an operation that cannot bound its scratch must say so rather
/// than allocating silently in a hot loop.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ScratchRequirement {
    /// Operates in place; no scratch beyond inputs and outputs.
    None,
    /// At most `bytes` of scratch, independent of input size.
    Fixed {
        /// Upper bound in bytes.
        bytes: usize,
    },
    /// At most `bytes_per_element * elements` of scratch.
    PerElement {
        /// Upper bound per input element, in bytes.
        bytes_per_element: usize,
    },
    /// Scratch cannot be bounded ahead of time.
    ///
    /// Callers must treat this as "may allocate arbitrarily"; a memory budget
    /// cannot be enforced against it.
    Unbounded,
}

impl ScratchRequirement {
    /// Upper bound for `elements` inputs, or `None` when unbounded.
    pub const fn upper_bound_bytes(self, elements: usize) -> Option<usize> {
        match self {
            Self::None => Some(0),
            Self::Fixed { bytes } => Some(bytes),
            Self::PerElement { bytes_per_element } => bytes_per_element.checked_mul(elements),
            Self::Unbounded => None,
        }
    }

    /// Whether this requirement fits `options`' memory budget for `elements`.
    ///
    /// An unbounded requirement never fits a declared budget: allowing it would
    /// make the budget advisory, which is the failure this type exists to stop.
    /// Whether this requirement fits the caller's memory budget.
    ///
    /// No longer `const`: it reads an [`ExecutionOptions`] that now owns a
    /// shared cancellation handle, so it cannot be destructured at compile
    /// time. Budget checks happen once per dispatch, not on a hot path.
    pub fn fits_budget(self, options: &ExecutionOptions, elements: usize) -> bool {
        match options.memory_budget_bytes() {
            None => true,
            Some(budget) => match self.upper_bound_bytes(elements) {
                Some(needed) => needed <= budget,
                None => false,
            },
        }
    }
}

/// Upper bound on how many outputs an operation produces per input element.
///
/// Declared so a caller can size a destination buffer *before* the operation
/// runs. That is what makes a scan-then-scatter implementation possible: with a
/// per-element bound, exclusive-prefix-summing the per-element counts gives
/// every worker a disjoint write offset, so no lock, no atomic counter, and no
/// dynamically growing vector is needed on the hot path.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OutputBound {
    /// Exactly one output per input. The output index equals the input index.
    OneToOne,
    /// At most `max` outputs per input element.
    AtMost {
        /// Inclusive upper bound per element.
        max: usize,
    },
    /// The output count cannot be bounded before running.
    ///
    /// Callers must fall back to a growable collection. Kept representable and
    /// deliberately unpleasant so an operation that could declare a bound is
    /// not tempted to shrug.
    Unbounded,
}

impl OutputBound {
    /// Worst-case total outputs for `elements` inputs, or `None` if unbounded
    /// or the product would overflow.
    pub const fn upper_bound(self, elements: usize) -> Option<usize> {
        match self {
            Self::OneToOne => Some(elements),
            Self::AtMost { max } => max.checked_mul(elements),
            Self::Unbounded => None,
        }
    }

    /// Whether a caller can preallocate an exact destination for `elements`.
    pub const fn is_preallocatable(self, elements: usize) -> bool {
        self.upper_bound(elements).is_some()
    }

    /// Exclusive prefix sum of `counts`, plus the total.
    ///
    /// This is the scan that turns "how many outputs does each element make?"
    /// into "where does each element write?". Returns `None` if any count
    /// exceeds this bound (a provider contract violation) or the total
    /// overflows. Allocation is the caller's single result buffer; the scan
    /// itself is a running total with no per-element allocation.
    pub fn write_offsets(self, counts: &[usize]) -> Option<(Vec<usize>, usize)> {
        let per_element_max = match self {
            Self::OneToOne => Some(1),
            Self::AtMost { max } => Some(max),
            Self::Unbounded => None,
        };
        let mut offsets = Vec::with_capacity(counts.len());
        let mut running = 0usize;
        for &count in counts {
            if let Some(max) = per_element_max {
                if count > max {
                    return None;
                }
            }
            offsets.push(running);
            running = running.checked_add(count)?;
        }
        Some((offsets, running))
    }
}

/// Operation policy with explicit tolerance, precision, and cancellation.
///
/// Deliberately not `Copy`: it carries a shared [`CancellationToken`], and an
/// implicitly copied cancellation handle is a footgun. Callers clone when they
/// mean to share the token and construct fresh options when they do not.
#[derive(Debug, Clone, PartialEq)]
pub struct ExecutionOptions {
    tolerance: Tolerance,
    precision: Precision,
    determinism: Determinism,
    parallelism: Parallelism,
    device: DevicePreference,
    residency: DataResidency,
    memory_budget_bytes: Option<usize>,
    cancellation: Option<CancellationToken>,
    chord_error: Option<Scalar>,
}

impl ExecutionOptions {
    /// Start from the required model-aware tolerance.
    pub const fn new(tolerance: Tolerance) -> Self {
        Self {
            tolerance,
            precision: Precision::F64,
            determinism: Determinism::NumericallyBounded,
            parallelism: Parallelism::Auto,
            device: DevicePreference::Auto,
            residency: DataResidency::HOST,
            memory_budget_bytes: None,
            cancellation: None,
            chord_error: None,
        }
    }

    /// Replace the tolerance while preserving every other execution policy.
    ///
    /// Compilers use this when evaluating geometry in a transformed local
    /// coordinate system. Cloning first keeps cancellation, budgets,
    /// determinism, residency, and provider preferences intact.
    pub fn with_tolerance(mut self, tolerance: Tolerance) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Attach a cooperative cancellation token.
    ///
    /// Absent a token, an operation runs to completion; there is no ambient
    /// cancellation source. Providers declare how finely they poll via
    /// [`crate::CancellationGranularity`], so a caller can see the real latency.
    pub fn with_cancellation(mut self, token: CancellationToken) -> Self {
        self.cancellation = Some(token);
        self
    }

    /// The attached cancellation token, if any.
    pub fn cancellation(&self) -> Option<&CancellationToken> {
        self.cancellation.as_ref()
    }

    /// `Err(GeomError::Cancelled)` if a token is attached and cancelled.
    ///
    /// The call providers make at each poll point: `options.check_cancelled()?`.
    /// Cheap (one relaxed load) and a no-op when no token is attached.
    pub fn check_cancelled(&self) -> crate::GeomResult<()> {
        match &self.cancellation {
            Some(token) => token.check(),
            None => Ok(()),
        }
    }

    /// Set required precision.
    pub fn with_precision(mut self, precision: Precision) -> Self {
        self.precision = precision;
        self
    }

    /// Set determinism requirement.
    pub fn with_determinism(mut self, value: Determinism) -> Self {
        self.determinism = value;
        self
    }

    /// Set scheduling preference. Returns `None` for zero explicit threads.
    pub fn with_parallelism(mut self, value: Parallelism) -> Option<Self> {
        if matches!(value, Parallelism::Threads(0)) {
            return None;
        }
        self.parallelism = value;
        Some(self)
    }

    /// Set device preference.
    pub fn with_device(mut self, value: DevicePreference) -> Self {
        self.device = value;
        self
    }

    /// Declare where inputs live and where outputs are wanted.
    pub fn with_residency(mut self, value: DataResidency) -> Self {
        self.residency = value;
        self
    }

    /// Bound how far flattened curves may deviate from the exact curve.
    ///
    /// Curved geometry that a provider approximates with straight chords
    /// (profile arcs, sweep directrices, curved B-rep faces) stays within this
    /// distance of the exact curve. Without it the chord budget is the linear
    /// tolerance, which is a coincidence tolerance and coarse for small radii:
    /// at `Tolerance::MILLIMETRE` a 5 mm arc gets four chords per half turn.
    /// Quantity take-off wants a tighter budget than display.
    ///
    /// Returns `None` for a non-finite or non-positive value. A provider may
    /// still refuse a budget too fine to meet within its own work limits.
    pub fn with_chord_error(mut self, value: Scalar) -> Option<Self> {
        if !(value.is_finite() && value > 0.0) {
            return None;
        }
        self.chord_error = Some(value);
        Some(self)
    }

    /// The explicit chord budget, if one was set.
    ///
    /// `None` means the provider uses its default, the linear tolerance.
    pub fn chord_error(&self) -> Option<Scalar> {
        self.chord_error
    }

    /// Bound temporary allocation.
    pub fn with_memory_budget(mut self, bytes: usize) -> Self {
        self.memory_budget_bytes = Some(bytes);
        self
    }

    /// Tolerance.
    pub fn tolerance(&self) -> Tolerance {
        self.tolerance
    }

    /// Precision.
    pub fn precision(&self) -> Precision {
        self.precision
    }

    /// Determinism requirement.
    pub fn determinism(&self) -> Determinism {
        self.determinism
    }

    /// Scheduling preference.
    pub fn parallelism(&self) -> Parallelism {
        self.parallelism
    }

    /// Device preference.
    pub fn device(&self) -> DevicePreference {
        self.device
    }

    /// Optional temporary-memory budget.
    pub fn memory_budget_bytes(&self) -> Option<usize> {
        self.memory_budget_bytes
    }

    /// Where inputs live and where outputs are wanted.
    pub fn residency(&self) -> DataResidency {
        self.residency
    }

    /// Charge `bytes` of scratch against the budget before allocating it.
    ///
    /// The budget is only real if something checks it, so this is the single
    /// enforcement point every hot path routes through. Backends must call it
    /// *before* the allocation, not after: reporting an overrun once the
    /// allocation already succeeded defeats the purpose of a budget.
    pub fn charge_scratch(&self, bytes: usize) -> GeomResult<()> {
        match self.memory_budget_bytes {
            Some(budget) if bytes > budget => Err(GeomError::BudgetExceeded { resource: "memory" }),
            _ => Ok(()),
        }
    }
}
