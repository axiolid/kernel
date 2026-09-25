# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.1] - 2026-09-25

### Added

- `ExecutionOptions::with_chord_error` and `ExecutionOptions::chord_error`
  (#165): an explicit bound on how far a provider's straight chords may sit
  from the curve they replace, separate from the linear tolerance. The
  tolerance is a coincidence test; used as a chord budget it leaves a 5 mm
  arc a few chords at `Tolerance::MILLIMETRE`, and small profiles mesh
  percent-level off. `None` (the default) keeps each provider's previous
  behaviour. A non-finite or non-positive budget is refused (`None`).
