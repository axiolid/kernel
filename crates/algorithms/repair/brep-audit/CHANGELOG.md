# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Changed

- An implicit pcurve (ADR 0077) is parameterised by its cells, not in
  proportion to its edge. It passes when every lifted sample projects onto
  the edge within tolerance, inside the edge's span, in the order the use
  runs, and starting and ending at the use's ends. Every other pcurve is
  still checked against the edge at proportional parameters.
