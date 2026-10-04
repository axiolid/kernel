# 0083 — Profile caps use a certified ear clipper

- **Status:** Accepted
- **Date:** 2026-10-04
- **Deciders:** Friedrich, axiolid
- **Supersedes:** the `axiolid-construct` part of 0015

## Context

ADR 0015 adopted `earcut` to triangulate profiles with holes after a
hand-rolled ear clipper failed on two holes. It recorded the risk: earcut is
not exactness-certified, and the hole case was covered by area conservation
only.

That risk was real (#253). For a 4 x 4 profile with two 1 x 1 holes side by
side in one horizontal band, earcut returned 12 triangles of the right total
area where a polygon with 12 vertices and 2 holes needs 14. One triangle
edge ran along the band's bottom line straight past both holes' inner
corners: a T-junction. The plain extrusion of that profile had 8 boundary
edges and failed its volume and closure checks. earcut drops nodes where
the bridged ring runs straight on and lets only reflex nodes block an ear,
so a vertex lying on a diagonal is not seen. Area conservation cannot catch
this: a T-junction covers the right area.

## Decision

`axiolid-construct` triangulates profile rings with its own ear clipper
(`ring_triangulation`), every decision an exact `orient2d` sign, and
certifies its own output before returning it.

- **Validate.** Rings that do not bound a polygon with holes are refused
  with `InvalidInput` naming the ring: a non-finite or repeated vertex, a
  ring folding back on or crossing itself, holes that overlap or touch each
  other or the outer ring, a hole outside the outer ring or inside another
  hole. Rings may come either way round; orientation is read exactly at
  each ring's lexicographically smallest vertex.
- **Bridge.** Holes join the boundary in order of decreasing largest x,
  each from its rightmost vertex, so a mutually visible vertex exists on the
  boundary built so far. Candidates are tried nearest first; one is taken
  only when the bridge leaves both ends into the polygon and touches no
  edge of any ring or earlier bridge, a vertex on the segment included.
- **Clip.** An ear is refused if any node other than its corners lies in
  the closed triangle, so no diagonal runs through a vertex.
- **Certify.** Every triangle strictly counter-clockwise, every ring edge
  used exactly once from inside, every other edge exactly once in each
  direction. Positive triangles bounded exactly by the rings tile the
  polygon once, so this is the whole contract; a triangulation that fails
  it is refused with `Degenerate`, never returned.

The ADR 0015 lesson is kept: the bookkeeping that defeated the first
clipper is now checked by the certificate rather than trusted, and the
tests sweep random hole layouts on a shared grid, where rays, columns and
edge lines coincide constantly.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Keep earcut and repair T-junctions afterwards | Treats one symptom of an uncertified decision procedure; the next degeneracy surfaces as another open solid. |
| `axiolid-triangulate` (constrained Delaunay) | A new internal edge for `axiolid-construct`, and its incremental build rebuilds adjacency per insertion, quadratic on a finely chorded profile. Delaunay quality is not needed for caps. |
| Report the certificate failure but still use earcut | Refuses valid profiles instead of triangulating them. |

## Consequences

**Positive**

- Profiles with several holes compile to closed, consistently wound
  extrusions; polygonal input covers exactly outer minus holes.
- Invalid ring sets are refused by name instead of producing a solid.
- `axiolid-construct` no longer depends on `earcut`.

**Negative / costs**

- More code owned here. The search grid keeps a 24,576-vertex profile at
  about 0.13 s in a release build; earcut was faster on such inputs.
- Rings that touch at a single vertex, which earcut accepted, are now
  refused: their extrusion is not a two-manifold.

**Follow-ups / risks to watch**

- `axiolid-mesh-compile` still uses earcut for planar faces; ADR 0015 stands
  there.

## Relation to existing code

- `crates/algorithms/construction/construct/src/ring_triangulation.rs` and
  its `validate`, `bridge` and `clip` submodules.
- `crates/algorithms/construction/construct/src/profile.rs` —
  `triangulate` delegates here.
- `crates/algorithms/construction/construct/tests/profile_holes.rs` — the
  #253 cases and the random-layout properties.
- `crates/algorithms/construction/construct/tests/oracle.rs` — the
  hole-free differential gate against `axiolid_reference::triangulate_simple`.
