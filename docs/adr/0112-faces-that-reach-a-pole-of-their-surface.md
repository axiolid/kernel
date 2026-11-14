# 0112 — Faces that reach a pole of their surface

- **Status:** Accepted
- **Date:** 2026-11-14
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0086](./0086-curved-faces-invert-their-edges-without-pcurves.md) (the refusal of a sample at a cone apex or sphere pole), [0103](./0103-b-rep-faces-on-swept-and-ruled-surface-relations.md) (D6), [0110](./0110-edges-without-pcurves-are-inverted-onto-surface-relations.md) (D4 and D6)

## Context

Issue #463. A curved B-rep face that reaches a pole of its surface, a
point that every `u` of a line `v = v_pole` names, was refused by the
tessellator: a sphere's pole, a cone's apex, the three-quarter cone
`revolve_profile_exact` builds, and, on surface relations (ADR 0103, 0110),
a revolution whose profile reaches its axis and a ruled surface with an
apex side. Revit and IFC write such faces for domes, cones, tree canopies
and pleated luminaire shades.

ADR 0059 fixes the topology: a pole is a vertex and has no degenerate edge,
so a face's loop leaves the pole line with a gap in parameters that is no
gap on the surface. Two things were missing in the tessellator: a way to
close that gap in the face's parameter ring, and a bound the triangles
there can be certified with. Two existing bounds fail at a pole. The lens
of an edge without a pcurve (ADR 0086) scales by `|S_u|`, which vanishes
there, so its box never closes; and the triangle bound is first order at a
pole: an ideal fan triangle `(pole, X_k, X_k+1)` has a mixed term `r du /
4` where its sagitta is `r du^2 / 8`, and a flat ruled triangle into an
apex is refined into thousands although its mesh is exact.

## Decision

We close a face at a pole with a row of parameter samples welded to the
pole's vertex, read an edge's end at a pole by its `v` alone, bound the
spans and the triangles there along the meridians, and refuse by name the
cases the face's data leave undetermined.

**D1. Poles.** A support has a pole kind (`Support::pole_kind`): a sphere,
a cone and a surface of revolution turn about their poles (`Rotational`:
meridians are curves orthogonal to the parallels, `|S_v|` independent of
`u`); a ruled surface with an apex side has straight rulings into it
(`Apex`). A parameter point is a pole (`Support::is_pole`) when `|S_u|`, its
distance from the axis per radian, is within the linear tolerance, or, on a
ruled surface, when its point is within the tolerance of the apex.

**D2. Pole ends.** An edge's first or last sample at a pole is not
inverted. `SurfaceInverter::pole` reads its `v`: `+-pi/2` on a sphere, where
the cone's radius vanishes (an ulp to its side), a ruled apex's side, and a
revolution's by the certified search along the profile with the angle left
out (ADR 0110's cover, uniqueness, spread and residual). It is `None`
exactly where the inverter raises no pole refusal, and raises that refusal
for a point on the axis away from the surface. The end takes the `u` of its
neighbouring sample, so the edge's last span runs along a meridian in
parameters. An edge with one span between two poles is sampled again until
it has two. A sample at a pole that is no end of its edge still refuses:
the edge runs through the pole, and the producer splits it at the pole's
vertex (P5); one between samples turns a quarter period and refuses (P6). A
span into a pole is halved only while its middle stays off the pole: an
edge that cannot be certified there stops at its finest sampling and is
reported as it is.

**D3. Rows.** Where two uses of a ring meet at a pole, the ring gets the
samples `(u_k, v_pole)` between them, each carrying the pole's mesh index
(I1). The ring is cut into chains at its pole junctions, each localised by
continuity and, holding a pcurve use, shifted by whole periods so its first
pcurve use keeps its parameters (I2). With `K` pole junctions (I3):

- `K = 1`: the row runs from the chain's end back to its start.
- `K >= 2`, every chain anchored: the rows run as the pcurves place them.
- `K = 2`, `u` periodic, the chains a whole number of periods apart
  (including zero): one full turn, the first row `+P`, the second by
  closure.
- Otherwise: refused (P2).

A row narrower than the angular tolerance is no row; one wider than a
period is invalid (P3). A row copies the `u` values, bit for bit, of the
ring's one run of samples at one `v` spanning it (a rim), so the face stays
a lattice; otherwise it is uniform, no step above an eighth of a period or
`sqrt(8 eps / duu)` over the ring (I4).

The two-pole convention is normative: **between two poles the pcurves name
the side of the meridians a face covers.** Of a sphere's two halves between
the meridians at `u = 0` and `u = pi`, the first runs from `0` to `pi` and
the second from `pi` to `2 pi`, its far meridian at `2 pi`, not back at
`0`. Without pcurves only a full turn is determined.

**D4. Spans along the meridians (I6).** A span of an edge without a pcurve
whose end is a pole is bounded by `T = d + o + dvv dv^2 / 8`: `d` its chord
bound, `o` the pole vertex's distance from `S(u_neighbour, v_pole)`, `dvv`
the bound on `|S_vv|` over a box of every `u`, and `dv = |v_1 - v_pole| + L
/ s`, `L` the span's length and `s` a lower bound on `|S_v|` over that box,
shown positive in at most four rounds as the lens's box is. On a ruled
apex, `S_vv = 0` and `T = d + o`. The argument: take a point of the region
between the edge's preimage and the span's parameter chord (a meridian, as
the end took its neighbour's `u`), closed at the pole by the pole row. The
meridian through it is, within the region, a segment whose ends lie on the
preimage, whose surface points are edge points within `d` of the mesh
chord, or on the pole row, whose surface point is the exact pole within `o`
of the chord's end. Distance to a segment is convex, so the chord between
the two ends is within `max(d, o)` of the mesh chord, and the meridian arc
between them within `dvv dv^2 / 8` of that chord. `S_u` is orthogonal to
`S_v` on a support turning about its pole, so the edge's `v` moves at most
`|C'| / s`, and the preimage, starting at `v_1`, stays within `L / s` of
it inside the box; the box is clamped to the support's range of `v` (a
sphere's latitudes, a cone's side of its apex, a profile's domain) and to
the span's side of the pole. No quarter-turn condition is needed: the
meridians bound the region however far the preimage swings in `u` as it
reaches the pole.

The same argument bounds a span near a pole whose lens fails, where the
meridian's ends may also lie on the parameter chord, whose image is within
`e` of the mesh chord: `T = d + e + dvv dv^2 / 8`, with `dv = |v_b - v_a|
+ L / s`. It is used only where the span's ends show the edge off the
poles (their distances from the axis or apex sum to more than `L`, and
those distances change no faster than the edge moves), so no span is ever
bounded across a pole it passes. Without it, the lens failing on every span
next to a pole drives the edge sampler until a sample lands on the pole.
Where `s` cannot be shown positive or the box crosses a tangent break of
the profile, the face is unbounded by name (P9).

**D5. Fans (I7).** A ring that is a rectangle with a pole row is meshed as
a grid. Its lattice may be uneven: a row copies its rim's `u`, and a seam's
two uses share their samples' `v`, however unevenly the edge was sampled;
such a lattice is accepted only for a ring with a pole row, so no other
face changes. Each cell of a pole row is one triangle into the pole, from
the nominal triangle that does not collapse, bounded over the cell's whole
parameter quad `A, B` (the pole corners), `D, C` (the far corners): map the
unit square onto the quad bilinearly, `phi(s, t) = (1 - t) ((1 - s) A + s
B) + t ((1 - s) D + s C)`, `F = S o phi`. Interpolating `F` linearly in `s`
and then in `t` errs by at most `(sup|F_ss| + sup|F_tt|) / 8`, and since
`phi` is linear along each of `s` and `t`, `F_ss = S''[phi_s, phi_s]` with
`phi_s` a blend of `B - A` and `C - D`, so `|F_ss| <= duu d_u^2 + 2 duv d_u
d_v + dvv d_v^2` by its larger components (`F_tt` likewise from `D - A`
and `C - B`). The corners' offsets from the surface add. Both pole corners
carry the pole's vertex `P`, so the mesh's interpolant `(1 - t) P + t ((1
- s) X_D + s X_C)` lies on the fan's triangle and reaches all of it: the
quad's surface is within `Q` of the triangle both ways. For a lattice cell
this is `Q = (duu du^2 + dvv dv^2) / 8 + offsets`. A fan is split only at
the middle of its edge opposite the pole, into two fans on the quad's line
`s = 1/2`, while its term along the row is what keeps it above the budget;
its meridians are never split. A row is demoted, all of it, to the nominal
triangles when a cell's term across it plus its offsets misses the budget
or its box has no smooth bound.

**D6. The generic regime (I8).** A ring that is no grid (a teardrop, a
trim not along meridians and parallels) or a demoted row is meshed by the
clipper with the nominal triangles and bounded by the triangle bound and
the dropped-triangle rule. It is sound and first order at the pole: such a
face needs fine triangles there and may report above the budget.

**D7. Named refusals and reports.**

| Id | Error | When |
| --- | --- | --- |
| P1 | `UnsupportedInput`: "a curved B-rep face whose trim closes across a collapsed edge of a surface that is not a sphere, a cone, a surface of revolution or a ruled surface with an apex" | a trim jumps across a collapsed edge of a support without a pole kind (a B-spline's degenerate edge), #499 |
| P2 | `UnsupportedInput`: "a curved B-rep face between two poles whose edge uses carry no pcurves: the side of its meridians it covers is not determined" | D3's fall-through, #500 |
| P3 | `InvalidInput`: "a pole row of a curved B-rep face spans more than one period" | D3 |
| P4 | `UnsupportedInput`: "a hole of a curved B-rep face whose trim reaches a pole" | a row in an inner ring, #501 |
| P5 | `Degenerate`: "runs into the cone apex / sphere pole / axis / apex", unchanged | a sample at a pole that is no end of its edge |
| P6 | `Degenerate`: "turns a quarter period or more between consecutive samples", unchanged | an edge crossing a pole between samples |
| P7 | `Degenerate`: "runs into a collapsed ruling", unchanged | a sample on a collapsed ruling, now #502 |
| P8 | unchanged | a trim circling a pole with no seam edge |
| P9 | unbounded: "the surface is not certifiably regular along an edge into a pole" | D4 cannot show `s` positive or the box crosses a tangent break |

A clipper triangle joining three samples of one pole that covers parameter
area is `Degenerate`: it spans no surface.

**D8. Scope.** The collapsed ruling (ADR 0110's R3) leaves #463 for #502:
there `v` is undetermined, not `u`, and it needs a column, not a row. A
B-spline's collapsed edge (#499), a face between two poles without pcurves
(#500) and holes reaching a pole (#501) are refused by name. Known limits
without an issue of their own: a trim circling a pole with no seam edge is
refused as before, and a pole row that is no lattice is bounded only to
first order. Everything is crate-private in `axiolid-mesh-compile`: no
public API, no dependency edge, no closure change.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Require a degenerate edge with a pcurve along the pole line (OCCT's convention) | Contradicts ADR 0059, pushes the row onto every producer, and fixes neither bound. |
| Lower relations that coincide with a sphere or a cone to atomic surfaces | ADR 0103 rejected it, and atomic poles were refused too. |
| Rows only, no fan element | Sound, but first order at every pole: flat pleated faces would not meet the budget in any reasonable size. |
| Take the side between two poles from the loop's winding without pcurves | The tessellator normalises winding and has never read it; a wrong guess returns the complementary half silently (#500). |
| Root-find a pole inside an edge and split there | ADR 0059 makes the pole a vertex; producers split at it. |
| Lift every uneven lattice to a grid | Changes the mesh of every face whose rims are sampled unevenly; only rings with a pole row need it. |
| A polar chart around the pole | Rebuilds the parameter-space triangulation, its certificate and its bounds for one case. |

## Consequences

**Positive**

- Domes, cones, tree canopies and pleated shades mesh as solids or sheets
  with a certified bound, with and without pcurves, on atomic surfaces and
  on revolution and ruled relations. `revolve_profile_exact`'s cone and
  ball tessellate (#463's done-when).
- A flat ruled face into an apex is one triangle with bound 0; a cone's
  wall one cell high is one fan per rim span, bounded by the rim's sagitta.

Measured on this branch (budget and tolerance 1e-3; time is the compile's
wall clock in the optimised test build, one test thread; "measured" is the
farthest of a dense closed-form sample of the face from the mesh, read
exactly in a triangle's plane, or for the exact solids the boundary
Hausdorff lower bounds, both within the bound):

| Case | Oracle | Triangles | Time | Certified | Measured |
| --- | --- | --- | --- | --- | --- |
| cone, `revolve_profile_exact` | exact Hausdorff | 254 | 2.3 ms | 6.4e-4 | within |
| ball, `revolve_profile_exact` | exact Hausdorff | 27,904 | 79 ms | 9.8e-4 | within |
| hemisphere on a sphere, pcurves | dense sample | 3,966 | 10 ms | 9.4e-4 | 6.0e-4 |
| hemisphere on a sphere, none | dense sample | 16,382 | 77 ms | 6.7e-4 | 1.9e-4 |
| hemisphere by revolution, pcurves | dense sample | 8,318 | 95 ms | 1.003e-3 | 6.0e-4 |
| hemisphere by revolution, none | dense sample | 65,534 | 390 ms | 5.9e-4 | 4.7e-5 |
| sphere, one seam, pcurves | dense sample | 7,808 | 19 ms | 9.4e-4 | 6.0e-4 |
| sphere, one seam, none | dense sample | 16,256 | 34 ms | 1.13e-3 | 6.2e-4 |
| sphere by revolution, one seam, pcurves | dense sample | 32,000 | 251 ms | 7.6e-4 | 1.7e-4 |
| sphere by revolution, one seam, none | | | | unbounded (#498) | |
| two halves of a sphere, pcurves | dense sample | 7,808 | 23 ms | 9.4e-4 | 6.0e-4 |
| two halves by revolution, pcurves | dense sample | 32,000 | 342 ms | 7.6e-4 | 1.7e-4 |
| two halves, either, none | | | | refused (P2) | |
| cone with floor, pcurves (atomic and by revolution) | dense sample | 126 | 1.4 ms | 8.5e-4 | 6.0e-4 |
| cone with floor, none (atomic and by revolution) | dense sample | 254 | 27–35 ms | 4.5e-4 | 1.5e-4 |
| cone wall one cell high, pcurves | dense sample | 64 | 0.35 ms | 6.024e-4 | 6.013e-4 |
| cone to its apex, no pcurves (`edge_inversion.rs`) | dense sample | 64 | | 1.45e-3 | 4.8e-4 |
| three quarters of a ruled cone, pcurves | dense sample | 64 | 0.46 ms | 3.39e-4 | 3.38e-4 |
| three quarters of a ruled cone, none | dense sample | 64 | 2.3 ms | 1.02e-3 | 3.38e-4 |
| ruled apex patch, no pcurves (`relation_edge_inversion.rs`) | dense sample | 16 | | 1.47e-3 | 4.9e-4 |
| flat ruled face into an apex, either | exact in-plane | 1 | 0.3–0.6 ms | 0 | 0 |
| pleated shade, ten flat walls, either | dense sample | 18 | 0.5–1.7 ms | 5.6e-17 | 2.3e-16 (rounding) |
| shade with arc pleats, none | dense sample | 190 | 17 ms | 7.2e-4 | 2.4e-4 |
| canopy, axis to axis, one seam, pcurves | dense sample | 15,616 | 133 ms | 9.6e-4 | 5.9e-4 |
| canopy, axis to axis, one seam, none | dense sample | 32,512 | 254 ms | 9.1e-4 | 1.5e-4 |
| canopy, two half turns, pcurves | dense sample | 15,616 | 124 ms | 9.6e-4 | 5.9e-4 |
| canopy on a floor, pcurves | dense sample | 10,238 | 75 ms | 1.04e-3 | 5.4e-4 |
| canopy on a floor, none | dense sample | 32,766 | 172 ms | 6.7e-4 | 1.9e-4 |
| teardrop at a sphere's pole, none (I8) | dense sample | 2,392 | 41 ms | 1.74e-3 | 7.7e-4 |
| oblique small circle into a pole, none | dense sample | 2,284 | 18 ms | 1.73e-3 | 5.1e-4 |

The fan bound (I7) is also checked cell by cell against a dense sample of
the surface over quads next to a sphere's pole and a cone's apex, from
1.6 by 0.6 down to 0.05 by 0.4 in parameters, with corners lifted off the
surface too (`brep::pole` unit tests), and quarters as the quad halves.
The span term (I6) is checked the same way over the lens of circles tilted
0.1 to 1.2 radians from a meridian into a sphere's pole
(`brep::inversion` unit tests): the chord bound is needed there, the
meridian term is not, since at a smooth pole every curve bends at least as
much as the meridians beside it; it stays because the argument needs it.

**Negative / costs**

- On relation faces with pcurves, triangles along a seam can stay just
  above the budget (the relation's derivative bounds are looser than an
  atomic surface's, and a boundary edge is never split): the hemisphere and
  the canopy report up to 1.04 times the budget.
- Faces without pcurves report up to about 1.7 times the budget, as ADR
  0110 records for any face without pcurves, and the oblique small circle
  into a pole needs a dense edge.
- A face without pcurves whose seam crosses a tangent break of its profile
  (a semicircle as two rational quadratic arcs) is meshed but unbounded by
  name: the seam is not certified across the break (#498), so it stays
  coarse and its spans into the poles cannot be shown regular.
- The generic regime is first order at the pole and can reach the
  certified refinement depth near it.
- Unrelated to poles, found here: a partial turn from
  `revolve_profile_exact` tessellates with its base sector's triangles laid
  against the wall's rim the wrong way round (a three-quarter cylinder as a
  three-quarter cone), so it is not reported a solid; the cut faces and the
  wall close. The sector's loop runs clockwise about its plane's normal
  with the face forward, which the planar path reads by its winding and the
  curved path by the face's orientation alone.

**Follow-ups / risks to watch**

- #499: a B-spline surface's collapsed edge.
- #500: the side between two poles from the loop's winding.
- #501: holes reaching a pole.
- #502: a collapsed ruling.
- #497: edge sampling graded by local speed; a span next to a pole is still
  halved with the whole piece.
- #498: tangent breaks of the swept curve.

## Relation to existing code

`crates/execution/compile/src/brep/pole.rs` (rows, chains, refusals, fans
and their bound), `src/brep/support.rs` (`PoleKind`, `pole_kind`,
`is_pole`, `pole_distance`, `v_range`), `src/brep/inversion.rs`
(`SurfaceInverter::pole`, the pole end, `FaceBound::pole_span` and
`meridian_lens`), `src/surface_relation/invert.rs`
(`RelationInverter::pole`), `src/brep.rs` (the uneven lattice of a ring
with a pole row, the fans' refinement). Tests: the modules' unit tests,
`tests/pole_faces.rs`, `tests/edge_inversion.rs`,
`tests/relation_edge_inversion.rs`, `src/brep/exact_tests.rs`. Mutation
probe: `scripts/probe_pole_faces_mutants.py`. The ledger rows B5 and C8
record the level.
