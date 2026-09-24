# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- `arc_overlay` is exact (ADR 0070, #155). Crossings, their order along
  each edge, inside/outside/on classification, linking into rings and
  hole ownership are all exact sign decisions on the given input, via
  `axiolid-exact`; none reads the tolerance. Crossing points of two curves
  are rounded to `f64` once, in the output, and edges shorter than the
  tolerance after rounding are merged. The public API is unchanged.
- `cavalier_contours` is no longer a dependency. The arc path costs about
  90 to 150 us per boolean on typical sections instead of about 1 us
  (`benches/arc_overlay.rs`).

### Fixed

- `arc_overlay` results no longer depend on drawing units. The arc
  backend's thresholds are fixed in drawing units, so a 5 um gap survived
  a union drawn in millimetres but vanished in metres. The drawing is now
  scaled by a power of two so those thresholds sit at the caller's linear
  tolerance, capped so coordinates stay within what f64 resolves
  (ADR 0069). Superseded by the exact core above, which needs no scaling.
