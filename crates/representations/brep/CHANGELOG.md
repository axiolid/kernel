# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.3] - 2026-10-02

### Added

- `ExactBRep::transformed`: place an exact B-rep under a rigid motion
  (#223). Every vertex, curve and surface maps onto the same family --
  planes, cylinders, elliptical cylinders, cones, spheres, tori, circles,
  ellipses, B-splines and the section curves that carry their surface --
  with the topology, structural names and intervals kept; nothing is
  tessellated or refitted. A reflection keeps every frame right-handed
  (curved surfaces then read their angle backwards, `u -> 2 pi - u`, and
  their pcurves are reflected to match) and flips every face, so a solid
  stays outward oriented. Refused with the new `TransformError`: a
  non-finite transform, a linear part not orthonormal within
  `RIGID_TOLERANCE` (a scale or shear), an elevated alignment curve, and
  under a reflection the pcurve and carrier families with no closed-form
  reflection here (quadric and torus section graphs, traced, lifted and
  intrinsic curves on curved faces).

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Added

- `ExactBRepBuilder::append`: copy another exact B-rep's vertices, edges,
  loops, faces and shells, with their curves, surfaces, intervals and
  names, and return the new shell handles; optionally with every face used
  reversed, which turns an outer shell into a void (#111).
