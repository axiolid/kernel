# 0014 — Adopt `boolmesh` for the mesh boolean

- **Status:** Superseded by [0047](0047-absorb-mesh-boolean.md)
- **Date:** 2026-08-19
- **Deciders:** Friedrich, axiolid
- **Supersedes:** resolves the open evaluation in [0003](0003-pure-rust-mesh-boolean.md)

## Context

ADR 0003 decided *that* the mesh boolean is a pure-Rust implementation behind
`axiolid_contracts::MeshBoolean`, and explicitly deferred *which* implementation:
"Evaluate `boolmesh` and `manifold-rust` against the fixture corpus before
writing our own; adopting beats building if one passes."

This ADR records that evaluation. It was run in a throwaway crate outside the
workspace (`/mnt/backup/build-cache/csg-eval`), so no evaluation dependency
ever entered the workspace graph.

### Candidate survey (measured 2026-08-19 from crates.io and GitHub APIs)

| crate | version | licence | last push | transitive deps |
| --- | --- | --- | --- | --- |
| `boolmesh` | 0.1.9 | MPL-2.0 | 2026-05-07 | `glam` only |
| `manifold-rust` | 0.13.1 | **Apache-2.0** | 2026-08-19 | not evaluated yet |
| `manifold3d` | 0.4.0 | Apache-2.0 OR MIT | 2026-08-08 | C++ toolchain (excluded by ADR 0003) |
| `csgrs` | 0.20.1 | MIT | 2026-07-31 | BSP-tree based, not robustness-first |

Correction to ADR 0003: it recorded `manifold-rust` as a Manifold port without
a licence note. It is **Apache-2.0**, not MPL-2.0 — permissive, and therefore
the lower-friction option if it proves equivalent.

### Measured results for `boolmesh` 0.1.9

Dependency weight (`cargo tree`): the entire transitive graph is `glam`. No
`build.rs`, no `-sys` crate, no C++ toolchain. 4,519 LOC, **zero `unsafe`**,
edition 2024 (workspace MSRV is Rust 1.88). f64 is the default precision
(`K_PRECISION = 1e-12`); `f32` is an opt-in feature.

**Fixture `issue_2019_wall_two_overlapping_openings`** — a 4x0.2x3 wall minus
three mutually overlapping openings, two of them rotated off-axis (36.87 and
-53.13 degrees) so no axis-aligned fast path applies:

```
cut A axis-aligned   -> vol 2.160000  manifold=true
cut B rot  36.87 deg -> vol 2.080000  manifold=true
cut C rot -53.13 deg -> vol 2.080000  manifold=true
conservation: a\b + a^b = 2.400000, a = 2.400000, error = 0.000e0
```

Conservation is exact to the last bit. The result was cross-checked against an
independent 4M-sample Monte-Carlo integration of the same solids: MC gives
2.0807, boolmesh gives 2.0800 — agreement within MC noise. This matters because
volume conservation alone can be satisfied by a wrong-but-self-consistent
result; the MC check is an *independent* oracle.

**Fixture `issue_1155_halfspace_flyaway`** — an `IfcHalfSpaceSolid` clip at
millimetre scale whose plane normal is `(-1, 1.99999991118124e-9, 0)`, i.e.
2e-9 off axis. The historical failure is "flyaway": the clip emits geometry far
outside the input bounds.

```
column bounds  min(-125, -125, 11940)  max(125, 125, 23880)
clipped bounds min( 124.9999995, -125, 11940)  max(125, 124.9999, 23880)
volume grew    = false
manifold       = true
```

Result bounds stay strictly inside the input. No flyaway.

**Throughput** — the IFC-dominant wall-minus-N-openings pattern, release build,
sequential subtraction on one core:

| openings | total | per opening | manifold |
| --- | --- | --- | --- |
| 1 | 0.12 ms | 0.12 ms | true |
| 4 | 0.68 ms | 0.17 ms | true |
| 16 | 6.95 ms | 0.43 ms | true |
| 64 | 48.68 ms | 0.76 ms | true |

Per-opening cost grows with accumulated result complexity, as expected for
sequential subtraction. This is the number `batch_difference` and the
`subtract_many` batch override exist to improve; it is recorded here as the
pre-optimisation baseline, not as a final figure.

## Decision

Adopt `boolmesh` as the first concrete `MeshBoolean` provider, as a **normal
cargo dependency of a dedicated adapter crate** — never vendored, never edited.

`boolmesh` is MPL-2.0, which is *file-level* copyleft: modifications to its own
files must be published under MPL-2.0, but merely depending on it does not add
obligations beyond Axiolid's own MPL-covered source. Depending preserves a clean
licence boundary; vendoring and patching must be reviewed separately.

The provider lives behind `axiolid_contracts::MeshBoolean` per ADR 0003, so this
choice remains reversible: replacing `boolmesh` with `manifold-rust` or with an
in-house implementation is a change to one crate, not an API break.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Write our own robust boolean now | ADR 0003 already judged adoption preferable if a candidate passes. One passed, on the two fixtures chosen precisely because they are hard. Building our own now would spend the project's scarcest resource re-deriving a solved problem. |
| `manifold-rust` (Apache-2.0) | Genuinely attractive: permissive licence, pushed today, targets numerical parity with Manifold v3.5.0. Not yet evaluated against the fixtures. Recorded as the primary alternative to re-test if `boolmesh` disappoints — the seam makes that swap cheap. |
| `csgrs` | BSP-tree CSG. BSP is not the robustness-first approach the kernel needs, and MIT licence does not outweigh that. |
| Vendor `boolmesh` into the workspace | Triggers MPL-2.0 file-level copyleft on our modifications, and forks us off upstream fixes. Rejected. |

## Consequences

**Positive**

- The premise of ADR 0003 — robust CSG without a C++ kernel — is now
  *demonstrated on our own fixtures*, not merely asserted.
- One added transitive dependency (`glam`), zero `unsafe`, no build script.
- The riskiest unknown in the project is retired at the cost of one evaluation,
  not one implementation.

**Negative / costs**

- MPL-2.0 constrains what we may do to those files. Depending is fine; patching
  upstream code in-tree is not, without publishing under MPL-2.0.
- `boolmesh` is at 0.1.9 with 31 stars: young, low-traffic, one maintainer. Bus
  factor is a real risk. Mitigated by the trait seam and by `manifold-rust`
  being a known, licence-compatible fallback.
- `boolmesh` uses `glam`, which our geometry crates do not. The adapter must
  convert `TriMesh` <-> `boolmesh::Manifold` at the boundary; that conversion is
  a real cost to measure, not assume.

**Follow-ups / risks to watch**

- Winding matters: `boolmesh` treats inverted winding as an inside-out solid, so
  a subtraction silently behaves like a union. The evaluation hit exactly this
  and it produced a *plausible but wrong* answer (volume grew from 2.40 to
  2.86). The adapter must assert outward orientation on input, and the test for
  it must check volume *decreases*, not merely that the call succeeded.
- Evaluate `manifold-rust` against the same fixtures before this becomes load
  bearing, so the fallback is known-good rather than assumed-good.
- The `subtract_many` batch override should beat the sequential baseline
  recorded above; if it does not, the batch API is not earning its complexity.

## Relation to existing code

- `crates/contracts/common/src/boolean.rs` — the `MeshBoolean` trait and
  `MeshBooleanRegistry` this provider registers into.
- `docs/adr/0003-pure-rust-mesh-boolean.md` — the deferred decision this
  resolves.
- `test/fixtures/ifclite-geometry/` — the fixture corpus used above.

## Batch override: measured, 2026-08-20

The sequential baseline recorded above is now reproduced in-process by
`benches/subtract_many.rs`, and compared against the grouped override on the
same machine in the same run.

| layout | n | sequential | grouped | speedup |
|---|---|---|---|---|
| disjoint openings | 16 | 5.61 ms | 1.53 ms | **3.67x** |
| disjoint openings | 64 | 66.47 ms | 7.22 ms | **9.20x** |
| nested (complete overlap graph) | 16 | 2.11 ms | 2.11 ms | 1.00x |
| nested (complete overlap graph) | 64 | 8.18 ms | 8.28 ms | 0.99x |

The override groups mutually disjoint cutters by AABB and removes each group
with a single boolean, resting on `(S \ A) \ B == S \ (A union B)` plus the
fact that a concatenation of disjoint solids *is* their union.

**The decision this settles:** the batch API earns its complexity. The win
grows with `n` (the sequential path is superlinear because each step re-cuts a
subject that keeps growing), and the worst case costs 1% -- so the override is
safe to enable unconditionally rather than behind a heuristic.

Grouping is greedy first-fit colouring over the overlap graph. Optimal
colouring is NP-hard; the partition only needs to be good, and any partition is
*correct* because disjointness is verified before fusing. Bounding boxes
over-separate (two disjoint solids may have overlapping boxes) which costs a
group, never correctness; the converse cannot happen, which is what makes the
fused union sound.

An earlier "adversarial" case in the benchmark was not adversarial: openings
that each overlap only their neighbour 2-colour into two groups, so grouping
still won 2.9x. The nested layout above is a genuine worst case -- concentric
cutters give a complete overlap graph and N groups of one.
