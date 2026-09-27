# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- `triangulate`: after recovering constraint edges, unconstrained edges are
  flipped back to locally Delaunay (Lawson's flips). Recovery used to leave
  long triangles whose circumcircles held points no constraint hid, so the
  result was not constrained Delaunay as documented -- with finely sampled
  walls, triangles spanned whole rooms (#139).
