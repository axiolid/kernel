# axiolid-mesh-boolean-boolmesh plan

Owner: geometry
Depends on: axiolid-mesh-boolean-contract, axiolid-mesh-contracts, axiolid-contracts, axiolid-mesh, axiolid-core

Design notes and standing constraints. Status lives on GitHub, not here
(kernel#25).

## What this provider owns

TriMesh <-> Manifold conversion with an orientation gate on input, and
`MeshBoolean` for union/intersection/difference. Registry integration
honours budget refusal: an over-budget provider is never invoked.

`subtract_many` unions disjoint cutters before subtracting rather than
removing one cutter per boolean. The standing rule for that optimisation:
it must beat the sequential baseline recorded in ADR 0014 (n=16: 6.95 ms,
n=64: 48.68 ms). If it ever stops beating it, it does not earn its
complexity and should go.

## Gates

- Volume-conservation and winding gates; fixture issue_2019 regression.
- Fixture issue_1155 (near-degenerate halfspace). The half-space is still
  bounded in the test; moving that bounding into axiolid-model remains a
  separate concern.
- Differential test against certified `axiolid-predicates` (`orient3d`).
  Convexity, inside/outside and winding are re-decided exactly where those
  invariants hold; non-convex results stay covered by the conservation and
  structural gates instead, because "wound away from one interior point" is
  only true for convex solids.
