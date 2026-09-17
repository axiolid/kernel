# axiolid-heal implementation plan

Status: diagnosis implemented; repair not started. This is planning context,
not standing agent instruction.

## Standing invariants

- Crate boundary and dependency direction are executable in the layering gate.
- Public data/contracts compile. Behavior remains scaffold unless a test names it.

## Implemented

- `diagnose` produces a `Diagnosis` from a mesh: non-manifold edges,
  inconsistent winding, boundary edges, degenerate triangles, and
  self-intersecting triangle pairs.
- `self_intersections` decides triangle-triangle crossing exactly, through
  certified `orient3d` with a BVH broad phase. `self_intersections_brute_force`
  is the exhaustive reference the accelerated path is checked against.
- `Diagnosis::blocks_boolean` answers from measured defects rather than from
  vocabulary.

## Shape of the work

Repair. Every repair must report what it changed, per defect, so a caller can
audit the difference rather than trust it. Diagnosis landed first deliberately:
a repair that cannot name what it fixed is not auditable.

Coplanar overlapping triangles are currently reported conservatively as
intersecting; deciding them needs 2D region logic this crate does not own.

## Exit evidence

Targeted tests, feature-isolated compile where applicable, mutation-verified
architecture/validation gates, and benchmarks before performance claims.

## Repair status (#74)

- Implemented, all opt-in via `RepairPlan`: `WeldVertices`,
  `DropDegenerateElements`, `UnifyOrientation`, `OrientOutward`.
- `OrientOutward` closes the gap `unify_orientation` left open: unify makes
  neighbours agree, `OrientOutward` decides which way round the agreed shell
  faces, using `axiolid-measure` for the volume convention.
- NOT implemented: stitching coincident boundary edges into shared edges,
  and self-intersection removal. Detection of the latter exists; removal
  changes topology and is deliberately deferred.
- No repair runs implicitly. Nothing in compilation or dispatch calls into
  this crate.
