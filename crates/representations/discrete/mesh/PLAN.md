# axiolid-mesh implementation plan

Status: structural audit implemented; repair and richer topology remain planned.

## Standing invariants

- Crate boundary and dependency direction are executable in the layering gate.
- `TriangleMeshView` adapts foreign index storage without ownership conversion.
- `audit_mesh` reports malformed input, non-finite coordinates, tolerance-aware
  degenerate faces, boundary edges, and non-manifold edges deterministically.
  It does not mutate or reject dirty source geometry.

## Shape of the work

Add attribute channels, explicit repair plans, and richer topology diagnostics.

## Exit evidence

Targeted tests, feature-isolated compile where applicable, mutation-verified
architecture/validation gates, and benchmarks before performance claims.
