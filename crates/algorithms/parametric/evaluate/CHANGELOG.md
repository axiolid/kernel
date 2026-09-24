# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `Curve2::Sinusoid` evaluation: point, first and second derivative, a
  one-turn domain, and exact inversion (the parameter is the point's first
  coordinate, then its height is checked) (ADR 0071).
