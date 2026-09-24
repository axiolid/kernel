# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `MeshClosure` and `CompileOutcome::closure` (#161): whether a compiled mesh
  bounds a solid (`Solid`), is a surface model with area but no volume
  (`Surface`), or was not reported (`Unknown`, the default for `untracked`
  and `tracked`, so existing compilers build unchanged).
  `CompileOutcome::solid_mesh` returns the mesh only for `Solid`, so volume
  readers refuse a surface model instead of measuring a closed shell the
  source never declared a solid. `with_closure` sets it.

## [0.3.0] - 2026-09-23

### Added

- `CompileOutcome` and the provided method `MeshCompiler::compile_mesh_reported` (#115): a compiled mesh with the fate of each attribute channel. The default wraps `compile_mesh` and reports `attribute_fates: None` ("not tracked", not "nothing dropped"), so existing implementations compile unchanged.
