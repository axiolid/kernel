# 0098 — A union the mesh cut cannot settle is decided by the exact compiler

- **Status:** Accepted
- **Date:** 2026-11-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0096](./0096-a-union-decides-operand-contact-at-the-linear-tolerance.md)

## Context

Issue #391, a follow-up to #356. After ADR 0096, face contact between union
operands fuses when the operands share their linear part bit for bit: the
boolean is cut in their common local frame. Operands placed by linear parts
that differ in the last bits (a turn computed as two half turns) have no such
frame. Their contact faces are exactly coplanar in exact arithmetic and a
rounding residue apart once placed, and the mesh boolean, which has no
tolerance (ADR 0047), returns two bodies (refused by D2), a pinch (refused,
#194), or one body whose faces cross with sub-rounding wedge cracks, which no
guard sees. ADR 0096, D6 left that open as `TODO(#391)`.

The exact compiler reads such faces as one surface within the tolerance and
fuses them (ADR 0080). The exact compiler is the reference ADR 0096 D8 names
for the contract, so deciding the union there needs no second contact rule.
A finite symbolic perturbation in the mesh path was rejected for that
reason: it would be a second rule to keep in agreement with ADR 0080.

A characterization of the exact compiler on the #356 column and corbel (four
faces, fifty turns, the turn composed two ways, two sites) shows it fuses
every turn at the origin as one audited shell, and refuses 132 of 200 at a
site `(6e5, 5.6e6, 200)` away, with one composition as well as with two
(`exact_placed_boolean.rs`): a `brep-boolean` limit, not one of composition.

## Decision

**D1. The referee.** A union the mesh cut cannot settle is compiled by the
exact compiler, and the exact result is tessellated. It is decided on
geometry, not on how the mesh result looks:

- *unsettled contact*: after the cut, the pair of operand meshes actually cut
  (snapped, or as given) holds triangles that nearly share a plane within the
  linear tolerance and are not exactly coplanar
  (`axiolid_inspect::near_coplanar`, D4). Exact contact and crossings at an
  angle are decided cases;
- *a contact refusal*: the cut was refused as a pinch (#194) or by D2
  (`UNMERGED_UNION_CONTACT`).

Unions only. Operands cut in a shared frame (ADR 0096, D3) are decided there:
the faces that meet are exactly coplanar, or a residue apart that the snap
closes. The frame is shared only if the linear parts of the operands'
whole instance chains agree bit for bit, not only the outermost: sources
that were themselves placed by differently composed turns leave faces
nearly sharing a plane in any common frame, so the frame yields and the
world path re-reads and referees the union.

**D2. Never worse.** The referee's result is accepted only as a closed
`MeshClosure::Solid` with no pinch and no two bodies within the tolerance.
Where the exact compiler refuses, or its result fails those guards:

- an unsettled but valid mesh result is returned as before;
- a contact refusal keeps its message and is extended with
  `; the exact referee refused: {reason}`, so it still starts with
  `UNMERGED_UNION_CONTACT` or names the pinch;
- any other refusal is returned unchanged.

**D3. The tessellator.** The B-rep tessellator reads its supports through a
resolved view (#412): an edge's 3D curve, a use's pcurve with its parameter
span, a face's surface. `from_graph` resolves handles into the graph, and
`from_exact` into an `ExactBRep`'s own catalogs, taking each use's span from
`pcurve_interval`, which may run against the pcurve's direction; a use with
no recorded interval, or a handle outside its catalog, is refused. The
tessellator and its certified `Deviation` are the same for both; a family
without a certified bound is `Unbounded` by name.

**D4. Dependency.** `axiolid-inspect` gains `near_coplanar` and
`NearCoplanar`: the first pair of triangles (ascending, first then second) of
two meshes that are within `linear` of each other, whose corners all lie
within `linear` of the other's plane, and at least one corner of the second
is not exactly on the first's plane (certified `orient3d`); degenerate
triangles are skipped; `None` for a negative or non-finite `linear`. It is a
new public item, released with `axiolid-mesh-compile`. `axiolid-heal` becomes
a dev-dependency of `axiolid-mesh-compile` (the tests check results for
self-intersections); it was already in the production closure through
`axiolid-inspect`, so no closure profile changes.

**D5. The report.** The mesh is the tessellation of the exact union. Its
deviation is the worst tessellation contribution (unbounded wins) under
`REFEREED_UNION_EXACT` when the exact report holds no reading within the
tolerance, else under `REFEREED_UNION`, with `COINCIDENT_OPERANDS`
(`Certified`, the furthest move) when the exact compiler read faces as one
surface, and the operands' `SNAPPED_OPERANDS`-class perturbations carried.
The deviation report does not measure it again: it states its own bound. The
operands' attribute channels are `Dropped(ProviderLimitation)` by name, since
no triangle of the result descends from an operand's.

**D6. Out of scope.** Mesh operands are not refereed: the exact compiler
refuses them, and flattening them is reserved for the deviation report
(ADR 0080, #266). Differences and intersections are unchanged (ADR 0096,
D7).

## Alternatives considered

| Option | Why not |
| --- | --- |
| A finite symbolic perturbation in the mesh boolean | A second contact rule to keep in agreement with ADR 0080. |
| Referee only after D2 or a pinch | Misses the one-body result with wedge cracks, which no guard sees. |
| Referee every union | Costs an exact compile where nothing is near contact. Geometry decides: only a near-coplanar pair triggers it. |
| Fall back to the referee's error where it refuses | Worse than the mesh path: a valid mesh result would turn into a refusal. |
| Flatten mesh operands into the referee | ADR 0080 reserves flattening for the report; the exact compiler's public contract refuses meshes. |

## Consequences

**Positive**

- A union of differently composed turns is one solid wherever the exact
  compiler fuses it, with the exact union's own certified bound.
- The one-body result with sub-rounding wedge cracks is replaced by a mesh
  that does not cross itself.

**Negative / costs**

- A union with a near-coplanar pair is compiled as a mesh and again exactly.
  In a chain of N differently turned unions each level re-compiles its
  growing subtree exactly, which is quadratic in N.
- The exact compiler refuses most turns far from the origin
  (`exact_placed_boolean.rs`); there the mesh answer stands, or the mesh
  refusal with the referee's reason. That is a `brep-boolean` issue.
- The two defensive guards in the referee (a closed solid, no two bodies
  within the tolerance) cannot be reached by a constructible input today:
  the exact compiler reads contact at the tolerance the guard measures.

**Follow-ups / risks to watch**

- The exact compiler's refusals at distant sites.

## Tests

`tests/exact_union_referee.rs` (the sweep on four faces: one solid with no
self-intersection at the origin, the mesh path's answer far away; a pinned
case the mesh cut cracks and the referee replaces; a shared frame is not refereed, and
one whose sources carry the differing turns is; two tolerances apart stays two bodies; a gap within the tolerance
is read as contact and says so; a sheared placement and an edge contact keep
their refusals and name the referee; a dowel certified against its closed
form; an authored mesh operand returns the mesh result), `tests/union_contact.rs`
(the flipped refusal), the `referee` unit test (channels), `brep/exact_tests.rs`
(the tessellator against the exact B-rep: prism, cylinder, quarter torus,
box and dowel, a patch with reversed spans, an uncertified family),
`axiolid-inspect`'s `tests/coplanar.rs`, and
`scripts/probe_exact_union_referee_mutants.py`.

## Relation to existing code

- `crates/execution/compile/src/referee.rs`, `src/compiler.rs`
  (`build_boolean`, `boolean_of`, `cut`), `src/brep.rs` (`tessellate_exact`),
  `src/brep/view.rs` (`from_exact`), `src/deviation.rs`.
- `crates/algorithms/query/inspect/src/coplanar.rs`.
