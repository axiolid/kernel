# 0015 — Adopt `earcut` for polygon triangulation

- **Status:** Accepted
- **Date:** 2026-08-19
- **Deciders:** Friedrich, axiolid
- **Supersedes:** —

## Context

`MeshCompiler` needs to turn an IFC profile into triangles. Profiles have
holes: `IfcArbitraryProfileDefWithVoids`, hollow rectangular and circular
sections, and every wall whose section is a ring.

A hand-rolled ear clipper was written first, driven by the certified `orient2d`
from `axiolid-reference`. It passed on simple polygons, reflex vertices, a single
hole, and correctly refused degenerate and mis-oriented rings. It **failed** the
two-hole case: after splicing a second bridge the ring contains duplicated
vertices at both the array start and interior, and the cyclic ear walk stalls
with no ear found. Three fix attempts did not hold.

The failure is not conceptual — hole bridging is well understood — it is the
kind of index-bookkeeping problem where a mature implementation has already
absorbed the edge cases.

## Decision

Adopt `earcut` (MIT OR Apache-2.0) for polygon triangulation, behind
`axiolid-mesh-compile`, which is an L3 adapter crate exactly like `axiolid-mesh-boolean-boolmesh`.

`axiolid-reference` keeps its certified `triangulate_simple` for the hole-free case.
It is **not** dead code: it is the differential oracle that audits earcut on
every hole-free polygon, so the adopted implementation is verified rather than
trusted.

This follows the rule ADR 0003 already set for the mesh boolean: adopting beats
building when a candidate passes licence, robustness, and dependency-weight
gates.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Keep debugging the hand-rolled clipper | Another session at least, to re-solve a solved problem. The remaining Goal A risk is in extrusion, transforms, and boolean dispatch — not here. |
| `earcutr` (ISC) | Also viable, 2.7M recent downloads. `earcut` was chosen for the permissive MIT/Apache-2.0 licence, compatible with the workspace, and a more recent release (2026-07). |
| `spade` (Delaunay) | Solves a different problem: constrained Delaunay is heavier than needed, and quality triangles are not required for volume-correct extrusion. |
| `lyon_tessellation` | Aimed at path/vector rendering with anti-aliasing concerns; a much larger surface than needed. |

## Consequences

**Positive**

- The two-hole case works: measured area 175 on the exact polygon that defeated
  the hand-rolled version, and hollow-section extrusion is volume-exact.
- Dependency weight is small: `earcut -> num-traits -> autocfg`, pure Rust, no
  C++, no `-sys` crate.
- MIT OR Apache-2.0 is permissive and does not extend Axiolid's MPL-2.0 obligations to the dependency.

**Negative / costs**

- A second adopted geometry dependency. Mitigated by the same seam argument:
  `axiolid-mesh-compile` is the only crate that names `earcut`, and it is not
  re-exported, so replacing it is one crate's work rather than an API break.
- earcut is f64 but not exactness-certified; it can in principle emit a wrong
  triangulation on pathological input. The differential oracle covers the
  hole-free case; the hole case is covered by area conservation only.

**Follow-ups / risks to watch**

- Extend the differential oracle to holes if `axiolid-reference` ever grows a working
  hole path, so the audit is total rather than partial.
- Corner radii on rectangle profiles are not yet approximated; they currently
  produce sharp corners. *(Amended 2026-09-28: done under kernel#111; see
  `crates/algorithms/construction/construct/tests/rounded_rectangle.rs`.)*

## Relation to existing code

- `crates/algorithms/construction/construct/src/profile.rs` — flattening and the Earcut
  call; the only solid-generation site that names the dependency.
- `crates/algorithms/construction/construct/tests/oracle.rs` — the differential gate.
- `crates/algorithms/reference/src/polygon.rs` — certified simple-polygon
  triangulation, retained as the oracle.
