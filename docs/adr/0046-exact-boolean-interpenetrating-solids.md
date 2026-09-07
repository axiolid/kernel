# 0046 — Exact boolean path for interpenetrating solids

- **Status:** Accepted
- **Date:** 2026-09-06
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

`ScalarBoolean` is the oracle the boolean conformance suite is built on. It
was exact but partial: disjoint, nested and identical operands were answered,
and properly crossing surfaces were refused with `GeomError::Unsupported`.

That refusal was a property of what was implemented, not of the geometry. It
also capped what conformance could check, because the oracle and the
production provider could only be compared on cases the oracle would answer —
which excluded the interesting ones.

`boolmesh` handles interpenetration but is epsilon-tolerant, so it cannot
serve as an exactness reference for itself.

## Decision

We will resolve interpenetrating operands exactly, in three stages, and route
to them from a new `Arrangement::Interpenetrating` rather than erroring:

1. `intersection.rs` — the intersection curve, with nodes named by source
   topology (vertex index, edge/surface pair) rather than by position, and one
   canonical name per physical point.
2. `retriangulate.rs` — each cut face rebuilt against the curve, by brute-force
   constrained triangulation over exact `orient2d` signs.
3. `assemble.rs` — each resulting piece classified inside/outside by exact ray
   parity, then welded into the result.

Every decision is an exact predicate. Positions are computed only after the
crossing they represent has been proven to exist, so rounding can move a point
slightly but cannot invent or remove one.

The provider remains an ADDITION. `boolmesh` stays registered at the same
priority and remains the production path; the exact path is the oracle and the
fallback target, not a replacement.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Replace `boolmesh` with the exact path | Slower, and `boolmesh` handles shapes the exact path still refuses. Two paths with different trade-offs is the honest arrangement. |
| Keep refusing and rely on `boolmesh` alone | Leaves no exact reference for the case that matters most, so conformance cannot check the hard geometry. |
| Use a general CDT library for retriangulation | A face cut by a chain is convex-decomposable; a general CDT is a second thing to get wrong for no gain at reference scale. |
| Merge coincident points by tolerance | Welds genuinely distinct points. Identity is topological; only bit-identical coordinates from identical arithmetic are merged. |

## Consequences

**Positive**

- Interpenetrating booleans are answered exactly, verified against
  hand-computed volumes (union 1.875, intersection 0.125, difference 0.875 on
  half-offset unit cubes).
- Conformance gained a real cross-check: the exact and epsilon-tolerant
  implementations must now agree on geometry neither could be compared on
  before.
- Coplanar faces that share a plane but no area no longer block operations.
  Two walls flush in one plane and metres apart produced sixteen coplanar face
  pairs and a refusal; they now resolve.

**Negative / costs**

- The exact path is O(n·m) over face pairs with no spatial index. Acceptable
  for a reference and for the operand sizes conformance uses; not for
  production meshes.
- More code on the oracle path, which must itself stay correct for every
  conformance verdict to mean anything.

**Follow-ups / risks to watch**

- **Coplanar shared AREA is still refused.** Two solids flush over a whole
  face meet in a region whose boundary is currently derived per triangle pair,
  which yields edges interior to that region and a branching curve. Measured on
  two stacked cubes, continuing produced a `Difference` volume of 2.67 where
  the geometry says 8. Resolving it needs the overlap of the two face SETS per
  plane — a 2D polygon union — rather than of individual triangle pairs. This
  is the next piece of work on this path.
- Real IFC hits that shape constantly (stacked slabs, columns on floors), so
  the exact path is not yet a sufficient fallback for building models.
- No spatial index. Add one before this path sees production-sized input.

## Relation to existing code

- `crates/algorithms/reference/src/intersection.rs` — curve and node identity
- `crates/algorithms/reference/src/retriangulate.rs` — constrained retriangulation
- `crates/algorithms/reference/src/assemble.rs` — classification and welding
- `crates/algorithms/reference/src/coplanar.rs` — coplanar overlap clipping
- `crates/algorithms/reference/src/boolean.rs` — `Arrangement::Interpenetrating` routing
- `crates/providers/mesh/boolmesh/tests/conformance.rs` — oracle/provider agreement
