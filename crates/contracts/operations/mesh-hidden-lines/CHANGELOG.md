# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- The `MeshHiddenLines` operation contract (ADR 0090, #138): an
  orthographic view with optional near and far clipping, a 2D crop region,
  boundary, silhouette and crease edge selection, explicit limits,
  visible and hidden segments naming their source edge, evidence counts,
  `HiddenLineError` (a `GeomError` plus the mesh at fault), and a
  conformance suite with self-contained fixtures.
