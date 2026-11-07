# 0087 — Disc morphology has an exact and a fast mode

- **Status:** Accepted
- **Date:** 2026-10-19
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #330 (consumer: a downstream rule engine, free-floor and turning-circle rules).
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

## Amendment 2026-11-04: needle notches (#340)

#330's proposals offered chords whose spanned vertices lie within `t / 16`
of the line (0.0125 m for a disc of radius 0.4 m, `t` 0.2 m). A wall of
needles 0.05 m deep, under that disc, kept every notch: nothing was
filled and the sum kept its quadratic cost.

- **Proposals reach `t`.** `fill` proposes chords within `t` of their
  line. (This amendment first described two rounds, `t / 16` and then `t`
  over the edges left; the code has proposed once, at `t`, since c4314c30,
  and #428 corrects the record. On the 10,002-corner needle wall the inner
  erosion's excess is therefore 9.0e-4, not the 2.0e-6 two rounds gave;
  the time is the same. Two rounds would need their own issue.) A needle
  deeper than `t` was kept whole until #428.
- **Slab bound.** A notch whose spanned vertices all project onto the
  chord, of length `g`, has excess `e_W = sqrt(t^2 + g^2 / 4) - t`
  however deep it is. The proof is in `fill.rs`: a point within `t`
  of a point `n` of the notch lies within `sqrt(t^2 + g^2 / 4)` of one end
  of the notch's cross-section through `n` parallel to the chord, a point
  of the set. It joins the minimum of the depth and the gap bound. The
  fills of the needle walls measured here come from the gap bound (the
  tips lie on the chord), not from `e_W`, which only competes with it; for
  a notch deeper than `t` inside the slab (reached since #428), the
  minimum is of `e_W` and the depth. A notch deeper than `t` whose overhang
  beyond the chord's ends lies in a deep piece, or which has no vertex
  exactly on its chord to cut at, is refused (since #428 an overhang in a
  shallow piece is bounded by that piece's own bounds). The triage's first formula, the sagitta
  `t - sqrt(t^2 - g^2 / 4)` of a disc over the mouth, says how far a disc
  enters the notch, not how far the sum grows, which needs every centre
  within `t` of a point of the notch; it is also larger than `e_W` (by
  `g^4 / (64 t^3)` to leading order) and undefined for `g >= 2 t`.
- **Cost.** The toothed rooms and holes of #330 (2,402, 5,002 and 10,002
  corners) give the same results bit for bit in `Fast` mode, in the same
  time as before (0.046 s and 0.042 s for the hole dilation at 10,002
  corners, main and the amendment; 0.198 s with proposals at `t` alone, a
  variant not kept). The
  table of the three variants, and the needle walls' numbers, are in
  `docs/research/exact-morphology-on-needle-walls.md`.
- **`Exact` is unchanged.** Its profile on the needle shape is in that
  note; changing it stays a separate decision.

## Amendment 2026-11-07: trim the chord a deep split leaves (#345)

#345's floor (a 3 m square, teeth 0.5 mm deep on two sides, a budget of
1 mm) kept thousands of teeth under #340's proposals at `t`. The proposal
over the room's corner is deeper than the budget; it is split at the
corner, and the chord left runs from the corner, a valley of the wall,
along the wall. It is not proved, and the nearest-vertex split of a chord
not proved halves it between valleys again and again, down to single
edges. The cost was that of #330, from 13 s and 0.86 GB at 10,000 corners
to over 20 GB at 50,000.

- **Rule.** A chord not proved is tried without its start edge, without
  its end edge, then without both, before it is split as before; the first
  that is used (proved, within the disc and the budget, no ring in its
  notch) is taken. Its ends are then the wall's tips whichever the parity
  of the teeth.
- **Contract.** Unchanged: the proof, the excess and the budget check are
  those of any chord, so the side and the reported excess hold; `Exact` is
  unchanged. #330's and #340's `fast` cases give the same hashes.
- **Cost.** Linear from 50,000 to 200,000 corners: 0.39 s and 15 MB to
  1.57 s and 48 MB at r = 1.4 m. The numbers, the diagnosis and the
  hypotheses it rejected are in
  `docs/research/fast-morphology-on-toothed-floors.md`.

## Amendment 2026-11-07: needles deeper than the disc (#428)

Nothing proposed a chord over a needle deeper than `t`, and `excess`
would have bounded one by the whole wall: for a wall of needles 0.25 m deep
and 2e-4 apart, 2.5e-4 over 100 of them (2.06 m over a 4.5 m wall), where
one needle's slab bound is 2.5e-8. #340's room with a column and 124 such
needles gave 562 corners, as `Exact`, and 500 needles 2,066.

- **Reach.** Past `t`, a proposal still pops a vertex when the chord is no
  longer than `needle_width(t, budget) = 2 sqrt(2 t budget + budget^2)`,
  the widest mouth whose slab bound is within the budget, and a box
  containing what it spans projects inside it (a needle), or when the
  vertex lies exactly on the chord between its ends (a valley join).
  `chords`, the exact route, reaches none: bit for bit unchanged.
- **Bound.** A notch deeper than `t` has the least of its whole bound and
  the greatest bound of its pieces between the spanned vertices exactly on
  the chord. The pieces cover the notch, each lies against the same line
  with its spanned boundary in `S`, so `S' + B(t)` lies in
  `S + B(t + max e_i)`. Shallow notches are not cut: it would be sound but
  moves every earlier result by about `g^4 / t^3`.
- **Split.** A deep chord over the budget is split at its on-line vertex
  nearest the middle, not the farthest vertex: 7 chords and 12 corners on
  a wall of needles with one shallow tooth, against 101 and 108.
- **Contract.** Unchanged: every chord still goes through the proof, the
  budget and the check for rings in its notch, so the side and the stated
  deviation hold. The 30 earlier `fast` bench cases give the same hashes.
- **Cost.** The 124 needles give 68 and 72 corners. The timings of the
  committed `deep` bench cases are in `docs/research/fast-morphology-on-deep-needles.md`.
- **Limit.** A join needs a vertex exactly on the chord, as the exact proof
  does. Valleys collinear only up to rounding, as on a wall rotated by
  30 degrees, are not joined: each needle is filled alone but the valleys
  stay (2,500 needles: 385 to 456 corners against 68). A sound extension
  is #483.
