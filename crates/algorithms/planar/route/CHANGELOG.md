# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.1] - 2026-09-27

### Fixed

- A route could run along one wall, through a vertex and on across a gap
  outside the region where two walls line up (#187): visibility tested only
  proper crossings and the segment's midpoint. A segment is now cut at
  every obstacle vertex lying on it, decided exactly, and each stretch is
  either along an obstacle edge or has its midpoint in the region. Two
  rooms whose corridor is cut are `DisconnectedComponents` again.
