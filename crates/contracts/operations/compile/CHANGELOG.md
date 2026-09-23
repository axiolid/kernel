# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `CompileOutcome` and the provided method `MeshCompiler::compile_mesh_reported` (#115): a compiled mesh with the fate of each attribute channel. The default wraps `compile_mesh` and reports `attribute_fates: None` ("not tracked", not "nothing dropped"), so existing implementations compile unchanged.
