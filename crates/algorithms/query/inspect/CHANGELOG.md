# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `intersection_volume`, `difference_volume` and `enclosed_volume` (#183):
  the volume two closed triangle meshes share, the volume of the first
  outside the second, and a mesh's own volume, each as a `VolumeInterval`
  certified to contain the true value and never negative. No boolean is
  built: every face pair whose shadows overlap contributes the integral of
  the lower of the two planes over the overlap, with sides of lines and
  the lower plane decided exactly and the arithmetic rounded outward.
  Coplanar faces and touching bodies are exact cases. Open, non-manifold,
  self-intersecting or non-finite meshes are refused with the operand
  named (`OverlapError`); either winding is accepted.
