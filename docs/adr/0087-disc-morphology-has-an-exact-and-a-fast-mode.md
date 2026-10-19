# 0087 — Disc morphology has an exact and a fast mode

- **Status:** Accepted
- **Date:** 2026-10-19
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #330 (consumer: axioval/engine, free-floor and turning-circle rules).
`Region::erode_inner`, `erode_outer`, `dilate_inner` and `dilate_outer` in
`axiolid-overlay` build the disc morphology from the sum of a set with a
64-gon in one exact arrangement (ADR 0070), on a stated side
(`MorphologyBound`, #163). #323 merged exactly collinear corners and trimmed
each edge hull by chords proved to lie in the set, which made walls of
sub-millimetre teeth cheap for a room's erosion at up to about 1,200
corners. Profiling #330's hole and dilation cases showed the chords already
work on whichever side the operation sums, and none is refused. The cost
that remains is the tips' slivers overlapping one another. At 10,002
corners, with teeth 1.8 mm apart, each overlaps about 44 others, giving
about 905,000 pieces. On main in release, 10,002 corners take 3–11 s and
about 1.2 GB, and time and memory grow about with the square of the
corners.

The exact union of those slivers has cusps that are crossing points, not
sums of two `f64` vertices. An exact envelope would need rings with
non-`f64` vertices, which changes the exact core. A convolution offset would
put an arc of about 58° at every tip, crossing hundreds of neighbours,
which is worse. A rule engine needs a certified side and a stated
deviation, not bit-identical vertices. A CAD kernel needs results that do
not change between releases.

## Decision

The disc morphology gets two caller-selected modes, `MorphologyMode`. The
API is additive: `erode_inner_with`, `erode_outer_with`, `dilate_inner_with`
and `dilate_outer_with` take the mode, and the four existing methods keep
their behaviour as the exact mode.

- **Exact** (the default): today's algorithm, results bit for bit as
  before, pinned by hash in `tests/minkowski.rs` and by
  `cargo bench -p axiolid-benchmark --bench disc_erosion`. A CAD kernel uses
  this mode.
- **Fast** (`MorphologyMode::fast(budget)`, opt-in): before summing, the
  summed set's notches are filled by chords across its outermost corners
  (`overlay/src/fill.rs`). Each chord gets #323's exact proof, taken for the
  complement. It is used only when its *excess* is within the caller's
  budget. The excess is a proved bound on how much larger filling makes the
  sum: about the square of the spacing of the corners on the chord over
  eight times the radius, and never more than the notch's depth. The result
  keeps the inner or outer side. A circumscribed polygon (an inner erosion,
  an outer dilation) adds the excess to the deviation. An inscribed one (an
  inner dilation, an outer erosion) is drawn in by the excess. The excess
  is reported with the result: it is part of `MorphologyBound::deviation`,
  and `Region::excess` states it on its own. A rule engine uses this mode.

The budget is the caller's (metres). Fast mode refuses a budget that is not
finite and nonnegative with `OverlayError::InvalidOffsetDistance`. A chord
over the budget is split into shorter ones, down to single edges. So when no
notch fits the budget, the fast result is the exact one, with an excess of
zero.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Exact envelope of the slivers, with crossing-point vertices | Rings with non-`f64` vertices in the exact core (ADR 0070); a large change for one consumer class |
| Convolution of the boundary with the disc's arcs | An arc of about 58° at every tip crosses hundreds of neighbours: more pieces, not fewer |
| Trimming each sliver by its neighbours' hulls, exactly | Prototyped: same hashes, but only 0.57 to 0.33 s at 2,402 corners, and no gain for dilations; the kept ends of each covered run still span the sliver |
| Simplify the region by a tolerance before the morphology | No stated side unless the simplification is one-sided, and no bound tighter than the teeth's depth; the notch fill's excess is far smaller |
| Make the fast route the default | Changes results callers have pinned; a CAD kernel must keep them |

## Consequences

**Positive**

- On #330's cases (0.4 m disc, 0.5 mm teeth, release), fast mode takes
  0.036–0.045 s and 6–7 MB at 10,002 corners, against 2.9–7.8 s and about
  1.2 GB for exact mode, and grows about linearly with the corners. With a
  1 mm budget the excess is 2 µm, or 0.9 mm where a chord cuts a room
  corner.
- Exact mode is untouched, and its results stay bit for bit.

**Negative / costs**

- Two modes to document and benchmark. Fast results depend on the budget,
  and on which chords the greedy proposals and splits find, so they are not
  stable across releases the way exact results are. Only their side and
  deviation are promised.
- A notch that fits the budget is filled even where the result could have
  been exact. A room corner cut by a chord that fits the budget is an
  example (at 10,002 corners one such chord takes 0.9 mm of a 1 mm budget).

**Follow-ups / risks to watch**

- The excess is computed in `f64` with a slack of 64 ulps of the largest
  coordinate. It is a bound, not an exact certificate. The tests check the
  side and deviation pointwise against the true disc morphology.
- Splitting a refused chord is Douglas-Peucker-like: `n log n` on walls,
  `n^2` in the worst case.

## Relation to existing code

- `crates/algorithms/planar/overlay/src/fill.rs`: the notch fill, its proof
  and its excess.
- `crates/algorithms/planar/overlay/src/cover.rs`: #323's chord proposals and
  exact proof, shared.
- `crates/algorithms/planar/overlay/src/minkowski.rs`: `MorphologyMode`, the
  `*_with` methods, and how each side uses the excess.
- `tools/benchmark/benches/disc_erosion.rs`: both modes at 2,402, 5,002 and
  10,002 corners.
- `scripts/probe_fast_morphology_mutants.py`: the faults the tests must
  catch.
