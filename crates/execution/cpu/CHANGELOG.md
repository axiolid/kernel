# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.2] - 2026-10-13

### Fixed

- With `parallel`, a one-thread `CpuExecution` runs `install` on a local
  pool of one worker (#259). It ran the closure directly, on whatever pool
  the caller was in, so a provider reading `rayon::current_num_threads()`
  there saw the caller's width instead of the configured one, and sized
  its work or its scratch for it. A context from `CpuExecutionBuilder`
  starts the pool when it is built and reports a failure there; one from
  `CpuExecution::portable` or `CpuExecution::detect` starts it on first
  `install`. Clones share the pool.

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

