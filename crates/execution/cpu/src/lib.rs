#![deny(unsafe_op_in_unsafe_fn)]

//! CPU execution context shell with runtime ISA detection and optional local Rayon pool.
//!
//! Default build is portable and single-threaded. `simd` enables runtime ISA
//! selection for providers; this crate bundles no SIMD algorithm. `parallel`
//! creates a context-owned pool instead of mutating Rayon's global pool.

pub mod config;
pub mod execution;
pub mod features;
pub mod topology;

pub use config::{CpuConfigError, CpuExecutionBuilder, InstructionPolicy};
pub use execution::CpuExecution;
pub use features::{CpuFeatures, CpuInstructionSet};
pub use topology::{CacheLevel, CpuTopology};
