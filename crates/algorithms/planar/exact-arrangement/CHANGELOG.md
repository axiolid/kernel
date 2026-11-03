# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- The exact point store, segment arrangement and ring tracing that the plane
  section of `axiolid-reference` used (ADR 0093, #366), as a crate of their
  own, so that the visible faces of `axiolid-hidden-lines` decide on the same
  exact signs (ADR 0094, #383).
