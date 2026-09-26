# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `section_edges` (#167, ADR 0075 stage 1): the exact intersection curves of
  two exact B-reps' faces, each trimmed to where it lies inside both faces.
  Crossings with a boundary edge are found against the adjacent face's
  surface, or across a seam against the plane through the ruling.
- `split_face` (#167): a plane or cylinder face cut along its section edges
  into regions, traced in the face's parameters with exact pcurves (lines,
  conics, rulings, circles about the axis, `Sinusoid2` for oblique cuts).
