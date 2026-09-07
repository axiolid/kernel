#![forbid(unsafe_code)]

//! Benchmarking foundation: deterministic workloads and measured validation.
//!
//! # Why this crate exists
//!
//! Axiolid has execution seams — runtime ISA detection, an `ExecutionTarget`
//! taxonomy, a dispatch layer — but no trustworthy numbers flowing through
//! them. Optimising against absent measurement is guesswork: it cannot identify
//! which workload deserves attention, cannot quantify a regression, and cannot
//! tell an algorithmic win from noise. This crate is the measurement system
//! that has to exist first.
//!
//! # The rule every benchmark here obeys
//!
//! **A result that is not validated is not a measurement.** A kernel that
//! declines to answer, or answers wrongly, will happily look fast. Every
//! workload in this crate therefore carries a *derived* ground truth — computed
//! from the construction of the input, never from the code under test — and the
//! harness treats a mismatch as a failure rather than a win. That lesson is
//! carried directly from the cross-kernel harness in the sibling `benchmarks`
//! repository, where an unvalidated column once reported a kernel returning
//! zero volume in 0.2 ms as the fastest result in the table.
//!
//! # Determinism
//!
//! Workloads are generated from an explicit seed with a small counter-based
//! generator, never from hashing or wall-clock state. Two runs of the same
//! workload on any machine produce byte-identical inputs, which is what makes
//! instruction-count regression measurement meaningful.

pub mod dataset;
pub mod validate;
pub mod workload;

pub use validate::{Tolerance, Validated, ValidationError};
pub use workload::{Rng, Scale, Workload};
