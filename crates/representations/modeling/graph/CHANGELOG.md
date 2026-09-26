# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- A `Curve2::QuadraticGraph` with finite coefficients is accepted as a
  trim basis (#119).

- A `Curve2::Sinusoid` is a valid trim basis when its three coefficients
  are finite (ADR 0071).
