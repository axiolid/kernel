# 0045 — Boolean construction arithmetic is f64; predicates verify, they do not construct

- **Status:** Accepted
- **Date:** 2026-09-04
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

`axiolid-predicates` provides certified filtered `orient3d` and `orient2d`:
static filter, then exact expansion, then dyadic fallback. A reader finding
that in the tree reasonably infers that the production mesh boolean is exact.

It is not. In `crates/providers/mesh/boolmesh/Cargo.toml`, `axiolid-predicates`
is a **dev-dependency**. Predicates re-decide orientation *after the fact* in
tests; they do not drive intersection-point construction. Intersection
coordinates are computed and stored as `f64` in `TriMesh { positions: Vec<Point3>, .. }`.

This was never written down, so the gap between what the tree implies and what
the production path provides was left for a reader to discover.

## Measured behaviour

From `axiolid/benchmarks`, driving axiolid, upstream boolmesh, Manifold, CGAL
(`Exact_predicates_exact_constructions_kernel`) and OCCT through their own
APIs on byte-identical operands.

**Well-conditioned input: f64 construction costs nothing.** Relative error
against derived ground truth after N chained subtractions of 30°-rotated
cutters:

| n | axiolid | raw boolmesh | manifold | cgal | occt |
|---|---|---|---|---|---|
| 1 | 0 | 0 | 0 | 2.05e-16 | 0 |
| 8 | 3.97e-16 | 0 | 0 | 0 | 3.97e-16 |
| 64 | 5.62e-16 | 2.81e-15 | 0 | 1.12e-15 | 3.71e-15 |

Error does not grow with chain length, and axiolid is at or better than the
exact-arithmetic kernels.

**Near-degenerate input: it does cost.** Two unit cubes overlapping by `d`
along one axis, sweeping `d` toward an exactly coincident face. Exact answer
is `d`; the table is relative error. `1.00e0` means the kernel returned zero
volume — total loss, not degradation.

| overlap d | axiolid | manifold | cgal | occt |
|---|---|---|---|---|
| 1e-3 | 1.10e-13 | 1.10e-13 | 1.10e-13 | 1.10e-13 |
| 1e-6 | 9.81e-11 | 7.50e-11 | 8.43e-11 | 1.21e-10 |
| 1e-9 | 4.11e-8 | 1.80e-8 | 2.72e-8 | **1.00e0** |
| 1e-12 | 4.73e-5 | **1.00e0** | 3.34e-5 | **1.00e0** |
| 1e-15 | 8.52e-1 | **1.00e0** | 1.12e-1 | **1.00e0** |

Axiolid tracks CGAL closely to 1e-12 and never collapses; it is the second
most robust of the four on this probe. At 1e-15 CGAL is ~8x more accurate,
and that difference is exact rational constructions. OCCT collapses earliest
(1e-9) because `BRepBuilderAPI_Sewing` runs at a fixed 1e-9 tolerance —
"exact kernel" does not imply robust here.

The cost of exactness is real: rebuilding the CGAL shim with
`Exact_predicates_inexact_constructions_kernel` and changing nothing else
made it ~20% faster at n=1 and then **abort at n=4**.

## Decision

We will keep f64 construction in the production mesh boolean, and we will
stop leaving that implicit.

- Certified predicates remain **verification-only**. They re-decide
  orientation in tests and must not be read as an exactness claim about
  constructed coordinates.
- The provider **measures conditioning** and reports it in
  `BooleanEvidence`, rather than returning a badly wrong answer with no
  signal. At 1e-15 overlap the boolean previously returned 0.15x the right
  volume and said nothing.
- Conditioning is reported as **evidence, not policy**. Axiolid does not
  choose a refusal threshold on the caller's behalf; a caller decides what
  its domain can tolerate, exactly as it already does with
  `coincident_faces_encountered` and `output_components`.

Published thresholds, from the sweep above, expressed relative to operand
size:

| relative overlap | behaviour |
|---|---|
| above 1e-6 | degradation invisible; treat as trustworthy |
| 1e-6 to 1e-12 | degrades smoothly; usable with the conditioning flag observed |
| below 1e-12 | severe; results must not be trusted without independent checking |

## Alternatives considered

| Option | Why not |
| --- | --- |
| Rational/exact constructions throughout | `docs/ROADMAP.md` refuses OCCT parity and full CGAL reimplementation. Measured cost is real (the EPICK probe aborted at n=4), and the benefit appears only below 1e-12 relative overlap |
| Refuse below a fixed threshold | Picks a policy for every caller from one probe family. A structural model at millimetre scale and a mechanism at micron scale disagree about what is degenerate. Refusal remains available to the caller on the reported evidence |
| Say nothing and keep the status quo | The 1e-15 case returns a badly wrong volume silently, which contradicts the fail-closed stance in #16 |
| Snap near-coincident geometry | Silently changes the caller's input — the same class of hidden approximation the kernel exists to avoid |

## Consequences

- A reader of the tree can no longer mistake verification-only predicates for
  exact construction.
- Callers gain a conditioning signal they can act on, and Axiolid does not
  invent a domain-specific threshold.
- The gap to CGAL below 1e-12 remains open and is now documented rather than
  unmeasured. Closing it would require exact constructions, which this ADR
  declines.

## Amendment 2026-09-28: coordinate magnitude, and no RTC facility

A boolean measured at survey coordinates (1e6 to 1e7) looked badly wrong,
which was first read as a case for exact construction. Probes showed the
error was already present in the input's volume, before any boolean ran. It
came from catastrophic cancellation in a divergence-theorem sum taken about
the world origin. `volume_properties` and `surface_properties` in
`axiolid-measure` now sum about the mesh's first vertex, and
`tests/conditioning.rs` there holds thin plates at 1e7 to within 1e-12
relative. Re-basing per triangle was also measured. It was worse by up to six
orders of magnitude because it leaves one operand of `a . (b x c)` at world
magnitude, so it was not adopted. `second_moments` stays about the origin
because its value depends on that origin by contract.

Two further arguments for a relative-to-centre (RTC) coordinate facility were
measured and refuted (`axiolid-predicates` and
`axiolid-mesh-boolean-boolmesh`, `tests/survey_scale.rs`). Coincidence that
both sides author the same way survives rounding at every magnitude, so
`orient3d` still certifies `Zero` and a flush union is exact at 1e7. The
filter's escalation rate does not change with magnitude either, so RTC would
not make predicates faster. Axiolid therefore has no RTC or local-origin
facility and no offset ownership model. The input's own grid, about 1.9e-9 m
at 1e7, is a property of the caller's coordinates. Revisit this only if a
real workload fails in a way these probes do not cover.
