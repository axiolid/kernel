# Common base contracts plan

Status: superseded by ADR 0035. Common vocabulary remains here; operation schemas and execution policy were extracted to sibling packages.
not standing agent instruction.

## Standing invariants

- Crate boundary and dependency direction are executable in the layering gate.
- Public operation traits compile; the mesh-boolean registry stores executable
  trait objects rather than capability metadata.

## Shape of the work

Add a narrow batch trait only when a real implementation needs it; add an
operation-specific executable registry only when more than one provider exists.

## Exit evidence

Targeted tests, feature-isolated compile where applicable, mutation-verified
architecture/validation gates, and benchmarks before performance claims.
