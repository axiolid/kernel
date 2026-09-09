# axiolid-mesh-boolean-boolmesh instructions

Purpose: adapt the adopted `boolmesh` crate to `axiolid_contracts::MeshBoolean` (ADR 0014).
This crate owns conversion and contract enforcement; the algorithm is upstream's.

## Module ownership

convert.rs (TriMesh <-> Manifold, orientation gate); provider.rs (the trait impl,
result contract); box_detect.rs (axis-aligned box recognition); cellular.rs (the
analytic subtraction construction). Split before unrelated concerns grow together.

## Invariants

Orientation is checked on the way IN, per argument, naming which argument failed.
An inside-out mesh is structurally valid and manifold, so nothing else catches it;
`Difference` then behaves as `Union` and returns a LARGER mesh with no error. This
happened for real during the ADR 0014 evaluation.

Input faults are `InvalidInput`/`Degenerate`/`NotManifold` (caller's fault).
Result faults are `BackendContractViolation` (upstream's fault). Never blame the
caller for an upstream defect.

Scratch is `Unbounded`: `boolmesh` exposes no bound, so a caller with a hard
budget is refused rather than silently allowed past it.

Results carry no normals. `boolmesh` computes face normals; re-exporting them as
vertex normals would misrepresent the hard edges a cut creates.

## Verification

Volume conservation (`vol(a\b) + vol(a^b) == vol(a)`) is the gate, not index
comparison: it is triangulation-invariant, so it tests geometry rather than an
output buffer we do not control. Test helpers compute volume independently of the
crate's own helper, or the test would confirm the implementation with itself.

`boolmesh` must not be re-exported. It is MPL-2.0 and swappable; leaking its types
would make the adoption visible to consumers and defeat the seam.

## Batch override

`subtract_many` groups mutually disjoint cutters (AABB overlap graph, greedy
first-fit colouring) and removes each group with one boolean. Measured 9.2x at
n=64 on the IFC-dominant layout; 0.99x worst case, so it is unconditional.

Invariants, each mutation-proven in `tests/batch.rs`:

- **Only disjoint cutters may be fused.** Concatenating overlapping solids
  yields a self-intersecting mesh; subtracting it gives a wrong answer that
  still looks like a valid result. The disjointness check is load-bearing.
- **`fuse` must rebase indices.** Forgetting the offset silently duplicates the
  first mesh's triangles.
- **Every group must be subtracted**, and the single-member fast path must use
  that group's tool, not `tools[0]`.

`union_many` reduces in a BALANCED TREE rather than folding left. Same number
of booleans (n-1); the win is operand SIZE, since a fold makes step `i` union
an accumulator already holding `i` solids. Measured on a k^3 box grid
(`benches/union_many.rs`):

| n | fold | tree | speedup |
|---|---|---|---|
| 8 | 0.38 ms | 0.23 ms | 1.6x |
| 27 | 4.30 ms | 1.60 ms | 2.7x |
| 64 | 21.69 ms | 3.97 ms | 5.5x |
| 125 | 82.80 ms | 10.33 ms | 8.0x (7.5-9.5x across three runs) |

The ratio GROWS with n, which is what makes it a complexity difference rather
than a constant factor. On OVERLAPPING grids the win is smaller (1.1x at n=8,
1.9x at n=64): operands there merge into one growing solid, so the tree has
less small-operand advantage to exploit. Both numbers are reported.

Unlike `subtract_many` this has NO correctness cliff: union is associative and
commutative, so any reduction order yields the same solid, with no disjointness
precondition and no fusing. Invariants, mutation-proven in `tests/union_batch.rs`:

- **The odd trailing solid must ride to the next level.** Dropping it is
  invisible at even counts; `odd_counts_do_not_drop_the_trailing_solid` sweeps
  n=1,3,5,7,9,11. Deleting the `remainder()` push was verified to fail 4 gates.
- **Evidence must report n-1 sub-operations.** The tree does not save calls,
  and evidence claiming otherwise would misrepresent where the win comes from.
- **Reversing the operands must not change the answer** — the commutativity the
  regrouping relies on.

⚠️ A 125-box OVERLAPPING grid trips an assert inside the absorbed kernel
(`boolean45.rs`'s `pair_up`: odd edge-point count). **Pre-existing and unrelated
to reduction order** — reproduced on the sequential fold, which `union_many`
does not touch. Same class as the `hmesh.rs:67` Menger panic: an upstream assert
firing on hard geometry instead of returning a typed error. The bench caps the
overlapping sweep at n=64 to measure up to the cliff without pretending it is
absent.

Volume comparisons between the grouped and sequential paths use a RELATIVE
tolerance: the two sum a differently ordered triangle list, so the last bits
legitimately differ. Bitwise equality fails spuriously.

## Analytic box path (opt-in)

`subtract_boxes_analytic` cuts axis-aligned boxes out of an axis-aligned box in
closed form. ~25x faster than the general solver at n=64 openings.

**Opt-in, never auto-dispatched.** Unlike the batch override above (unconditional
because its worst case is 0.99x), this path changes the OUTPUT TOPOLOGY, not just
the schedule. Dispatching on shape would make triangle counts depend on whether a
wall's openings happened to be axis-aligned. The caller asks, and handles
`Ok(None)`.

Invariants, each mutation-proven in `tests/analytic_boxes.rs`:

- **Recognition is structural, never by bounding box.** Every mesh has a bounding
  box; a sphere and its enclosing cube share one. Acceptance requires an exact
  index count, all corners on the min/max lattice, and exactly 2 triangles per
  face plane.
- **The index-count check is not redundant with the plane check.** `chunks_exact(3)`
  silently drops a trailing partial triangle, so a malformed 38-index buffer
  presents a perfect box to the plane loop. Only the length check sees it.
- **The lattice check is not redundant either**, for one reason: it walks ALL
  positions, while the plane check only sees REFERENCED ones. An unused
  off-lattice vertex is invisible to the latter.
- **Refusal must stay a refusal.** Returning a wrong solid is worse than
  returning nothing; every decline case has a test.

Three independent oracles are needed, because each is blind to a different
defect:

- signed volume misses cancelling errors (an inverted face pair sums to zero);
- edge pairing misses coincident duplicate faces (each edge still balances);
- duplicate-face detection is the only one that catches an emitted interior face.

The interior-face mutant produced 96 triangles instead of 64 with IDENTICAL
volume and ZERO edge-pairing defects. Volume alone would have passed it.
