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
  each ring's lexicographically smallest vertex. Coordinates are ordered
  by value there and wherever the clipper picks or sorts by them, `-0.0`
  equal to `0.0` (#269): under `f64::total_cmp` a `-0.0` twin of the
  lowest corner's `x`, which projecting a face onto plane axes writes,
  was taken for the lowest corner, and a straight one reversed the ring.
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
  refused: their extrusion is not a two-manifold. 2D and surface callers
  accept them through `PinchPolicy::Accept` (amendment below, #262).

**Follow-ups / risks to watch**

- `axiolid-mesh-compile` still used earcut for planar faces; it moved to
  this clipper in #260 (amendment below), and ADR 0015 is superseded.

## Relation to existing code

- `crates/algorithms/construction/construct/src/ring_triangulation.rs` and
  its `validate`, `bridge` and `clip` submodules.
- `crates/algorithms/construction/construct/src/profile.rs` —
  `triangulate` delegates here.
- `crates/algorithms/construction/construct/tests/profile_holes.rs` — the
  #253 cases and the random-layout properties.
- `crates/algorithms/construction/construct/tests/oracle.rs` — the
  hole-free differential gate against `axiolid_reference::triangulate_simple`.

## Amendment 2026-10-04: pinch policy per caller (#262)

Refusing rings that touch at one vertex is right for a solid and wrong for
a region. A consumer (axioval) re-feeds overlay output as plan regions and
surfaces: two rooms meeting at a corner as one union ring, an L-shaped
room wrapped round a column corner, a corridor eroded to a point, the same
pinch reported through even-odd fill as a figure eight, a room minus door
zones touching its walls. Each bounds a valid region, and its
triangulation is a valid surface patch whose edges are the ring edges.

- **The caller states the policy.** `profile::triangulate_with(rings,
  PinchPolicy)`; `PinchPolicy` is `#[non_exhaustive]`. `Refuse` is
  `triangulate`, used by everything that must close into a two-manifold
  (extrusion and loft caps: a pinch extrudes to a wall edge shared by four
  faces). `Accept` is for 2D regions and planar surface patches, including
  the planar faces of a B-rep (a pinch inside one face leaves every edge
  of the shell shared by two faces; the shell is at most pinched at a
  vertex, which the face did not create).
- **What is accepted.** Two edges meeting at one point that is a vertex of
  at least one of them. A vertex inside another edge is inserted into it;
  coincident vertices become one, referenced by their first index in
  `outer ++ holes`. Crossing, overlapping along a stretch, a hole outside
  the outer ring or inside another hole stay refused by name.
- **How.** Each ring is split at its repeated vertices into simple loops;
  a loop's nesting depth (one probe vertex off the other loop decides,
  or, when every vertex is shared, the sector of the other loop that its
  first edge leaves into) orients it: even depths counter-clockwise, odd
  clockwise, so the region is the points inside an odd number of loops,
  and a figure eight through a vertex reads as its two lobes. Round every
  shared vertex the loops' edges, sorted by exact angle, must alternate
  leaving and arriving; each wedge from a leaving edge counter-clockwise
  to the next arriving one becomes one node. The nodes form the boundary
  cycles of the region's connected parts; a cycle bounds from outside
  when every visit of its lexicographically smallest vertex turns
  strictly left, and each hole cycle goes to the innermost outer cycle
  strictly around it. Each part is bridged and clipped as before; a
  bridge from a hole cycle that visits its rightmost vertex twice leaves
  from the visit whose sector it enters.
- **Certificate.** Unchanged in kind, over the split loops: strictly
  counter-clockwise triangles, every loop edge once from inside, every
  other edge twinned, and `n + 2h - 2c` triangles for `n` loop vertices
  (a pinch counted once per visit), `h` hole cycles and `c` parts.
- **Unchanged.** Rings touching nowhere take the original path under
  both policies and triangulate identically.

| Option | Why not |
| --- | --- |
| Accept pinches everywhere | Extrusions of pinched caps are not two-manifolds; the solid path must keep refusing by name. |
| Split the pinch into separate ring pieces only | Covers two lobes of one ring, not a hole touching the outer ring, which has to merge with it rather than split. |
| Emit two solids touching at an edge for a pinched cap | No caller asked for it; a solid caller that wants it can split its profile and extrude the parts. |

## Amendment 2026-10-04: mesh-compile's planar faces (#260)

`axiolid-mesh-compile` triangulated authored polygon faces, curve-bounded
planes, planar B-rep faces and curved faces' parameter domains with
earcut, so the T-junction above could open those meshes too, and the
planar B-rep path checked nothing at all. They now go through this clipper
(`planar::clip_projected`), over the dependency `axiolid-mesh-compile`
already had on `axiolid-construct`; no crate edge changes and the clipper
does not move.

- **Policy: `PinchPolicy::Accept` on every face.** A face is a surface
  patch. A pinch inside one face of a closed solid leaves every edge of
  the shell shared by two faces: the face's ring edges are exactly the
  edges its neighbours share, whatever its rings touch. At worst the shell
  is pinched at that vertex, which the face did not create; the mesh's
  own closure audit still reports it. Only extrusion and loft caps, which
  build the walls themselves, refuse pinches.
- **Corners.** A corner repeating its predecessor exactly (an exporter's
  closing point) is dropped first, as earcut dropped it; every other
  corner is a triangle corner, so curved faces no longer need their
  skipped trim samples put back.
- **The certificate replaces the area cross-check** the authored path ran
  on earcut's output, and closes the gap on planar B-rep faces, whose
  `Proven(0)` and slab bounds assume a cover. Rings that bound no region
  are refused by name where earcut returned a partial cover.
- **The noise-band split stays.** Export noise puts corners nanometres
  off the straight run they lie on; the clipper may cut a sliver across
  them whose diagonal a neighbouring face cuts too, using that edge four
  times. `split_invented_edges` still removes such slivers on planar faces.

## Amendment 2026-10-09: welding shared edges through inserted vertices (#265)

The #260 amendment claimed that under `Accept` "the face's ring edges are
exactly the edges its neighbours share, whatever its rings touch". That
holds for rings touching at a shared vertex, not for a vertex lying
inside another ring's edge: the clipper inserts it into that edge, so the
face runs along the edge in two pieces while the neighbouring face that
shares it keeps it whole. The mesh had a T-junction, and mesh-compile
still reported the B-rep as `Solid` (a consumer's closure audit caught
it).

- **Weld through.** `profile::ring_touches(rings)` returns every vertex
  the accepted triangulation inserts into an edge, from the same exact
  validation. Before triangulating, mesh-compile projects each planar
  B-rep or authored polygon face exactly as it will triangulate it,
  records each touch against the edge's two shared corners (topological
  vertex, position index), and inserts the recorded corners into every
  face's copy of that edge, in order along it. The touching face then
  sees the corner twice at one point, a pinch it accepts; the neighbour
  gains a corner on its edge. No position is moved or added, so the
  volume is unchanged and the shell is closed and two-manifold.
- **Check what is claimed.** A mesh labelled `Solid` (a declared B-rep
  solid, an operation's solid, a boolean result) is checked by index:
  every edge used exactly twice in opposite directions. Otherwise it is
  `MeshClosure::OpenSolid`, which `solid_mesh` and booleans refuse by
  name. This catches what the weld cannot reach, such as a curved B-rep
  face, which samples its own edges, sharing an edge with a planar face
  whose corner lies inside it.

| Option | Why not |
| --- | --- |
| Refuse vertices inside edges on B-rep faces | Valid faceted exports have them (pockets touching a face edge); refusing them drops geometry that tessellates correctly once welded. |
| Repair T-junctions after triangulating | Needs a geometric "on this edge" decision on the welded mesh, separate from the clipper's exact one; the pre-pass reuses the clipper's own. |
| Only demote to `Surface` | Loses the reason: a surface model and a broken solid are different findings for a consumer. |
