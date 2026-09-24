# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- `arc_overlay` results no longer depend on drawing units. The arc
  backend's thresholds are fixed in drawing units, so a 5 um gap survived
  a union drawn in millimetres but vanished in metres. The drawing is now
  scaled by a power of two so those thresholds sit at the caller's linear
  tolerance, capped so coordinates stay within what f64 resolves
  (ADR 0069).
