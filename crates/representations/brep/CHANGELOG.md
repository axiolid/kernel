# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `ExactBRepBuilder::append`: copy another exact B-rep's vertices, edges,
  loops, faces and shells, with their curves, surfaces, intervals and
  names, and return the new shell handles; optionally with every face used
  reversed, which turns an outer shell into a void (#111).
