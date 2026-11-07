# Exact disc morphology on needle walls

Measured evidence for #340, apart from any proposal. `MorphologyMode::Fast`
fills the notches between needles since #340 (ADR 0087, amendment
2026-11-04); `MorphologyMode::Exact` is unchanged, and whether to change it
is a separate decision.

## Shape and method

`tools/benchmark/benches/disc_erosion.rs`, walls `needle`: #323's room
(4.5 m wide) whose two long walls carry teeth 0.05 m deep, under a disc of
radius 0.4 m, so `t` is 0.2 m. The depth is above `t / 16` (0.0125 m), the
reach of #330's proposals, and below `t`. A room of `n` corners has
`(n - 2) / 4` teeth on a wall, so the tips are 36 mm apart at 502 corners,
18 mm at 1,002, 9.0 mm at 2,002 and 1.8 mm at 10,002.

Every row below is one case run alone, by setting `DISC_EROSION_CASE`
(for example `needle:erode:inner:2002:0.4:exact`) on a release build of the
bench, with `DISC_EROSION_RUNS=1`, and 3 for the 10,002-corner `fast` rows;
the bench's own sweep was not used for them. The machine was shared and
busy (load average 30 to 120), so absolute times are indicative; the growth
between sizes is the finding. Peak memory is the child's `VmHWM`.

## Exact mode

Release build. Erosion, inner side; the dilation, outer side, is within 5%
of it.

| corners | tip spacing | result corners | median s | peak MB |
|--------:|------------:|---------------:|---------:|--------:|
|     502 |       36 mm |            414 |     0.42 |     125 |
|   1,002 |       18 mm |            826 |     1.28 |     343 |
|   2,002 |      9.0 mm |          1,650 |     4.85 |   1,219 |

Doubling the corners multiplies the time by 3.1 then 3.8, and the memory by
2.7 then 3.6: the quadratic growth of #323 and #330 on a wall of teeth, here
with every tooth a needle whose edge hulls overlap the neighbours'. The
result keeps about one corner per tooth, so the cost is in the sum, not in
the answer. At 10,002 corners the same growth predicts minutes and tens of
GB; it is not run.

The cause is the one #330 named: the union of the slivers at the needles'
tips has cusps that are crossing points, not vertex sums, so the exact route
cannot sum fewer terms, and chord covers (#323) are proved exact and reach
only the notches within `t / 16` of the line, which these are not. Chords
that deep are `Fast`'s to make, because it states its excess.

## Fast mode

Release build. Before #340 `Fast` filled nothing on this shape: its
proposals reached `t / 16` (0.0125 m), so it paid what `Exact` does. Now
the needles' notches are filled. The excess is that of the gap bound: the
tips lie on the chord, so it is the square of their spacing over `8 t`,
well under the 1 mm budget and orders of magnitude under the needles'
depth.

| op, side | corners | tip spacing | result corners | median s | peak MB | excess |
|---|--------:|------------:|---------------:|---------:|--------:|-------:|
| erode, inner | 502 | 36 mm | 4 | 0.008 | 6 | 8.1e-4 |
| erode, inner | 1,002 | 18 mm | 4 | 0.014 | 8 | 2.0e-4 |
| erode, inner | 2,002 | 9.0 mm | 4 | 0.028 | 11 | 5.1e-5 |
| erode, inner | 10,002 | 1.8 mm | 4 | 0.150 | 35 | 9.0e-4 (2.0e-6 as first measured, see the correction below) |
| dilate, outer | 502 | 36 mm | 68 | 0.003 | 4 | 8.1e-4 |
| dilate, outer | 1,002 | 18 mm | 68 | 0.005 | 4 | 2.0e-4 |
| dilate, outer | 2,002 | 9.0 mm | 68 | 0.009 | 4 | 5.1e-5 |
| dilate, outer | 10,002 | 1.8 mm | 68 | 0.046 | 6 | 2.0e-6 |

Needles deeper than `t` were kept by `Fast` here, because the proposals
stopped at `t`; #428 fills them (`docs/research/fast-morphology-on-deep-needles.md`).

## Why the proposals came in two rounds

Correction (#428): the code has proposed once, at `t`, since c4314c30; the
two rounds below were measured but not kept. The inner erosion of the
10,002-corner wall therefore has the excess of `t` alone, 9.0e-4, not
2.0e-6; the times are those of the rounds' measurement.

Proposing at `t` alone fills the same notches as the two rounds, with
identical results (hashes) on all 24 cases of #330's toothed rooms and
holes, but costs more on walls of teeth: every chord as deep as `t` over a
wall is tried and split before the shallow ones that fill it. The two
rounds propose at `t / 16` first, as #330 did, then at `t` over the edges
left, so a wall of teeth is filled by the first round at the cost it had.
Release build, base and variants built alike from the same bench, the
best of three runs of five (`DISC_EROSION_CASE=<case>:0.4:fast
DISC_EROSION_RUNS=5`, one case per process; for example
`hole:dilate:inner:10002`), on a busy machine:

| case | corners | main s / MB | `t` alone s / MB | two rounds s / MB | excess |
|---|--------:|------------:|-----------------:|------------------:|-------:|
| hole dilation, inner | 2,402 | 0.013 / 5 | 0.028 / 6 | 0.014 / 5 | 3.5e-5 |
| hole dilation, inner | 5,002 | 0.025 / 5 | 0.082 / 10 | 0.026 / 5 | 8.1e-6 |
| hole dilation, inner | 10,002 | 0.046 / 6 | 0.198 / 19 | 0.042 / 6 | 9.0e-4 |
| room erosion, inner | 2,402 | 0.012 / 5 | 0.022 / 8 | 0.014 / 5 | 3.5e-5 |
| room erosion, inner | 5,002 | 0.026 / 5 | 0.035 / 11 | 0.024 / 5 | 8.1e-6 |
| room erosion, inner | 10,002 | 0.051 / 6 | 0.071 / 16 | 0.049 / 6 | 9.0e-4 |

The other 18 cases of #330 (room dilations, hole erosions, the outer sides)
were within the noise of main with `t` alone (0.93 to 1.02 times) in the
sweep of the bench's own driver. The two rounds also give the smaller
excess on the needle wall at 10,002 corners (2.0e-6 against 9.0e-4 for the
inner erosion with `t` alone), because the chords they find are shorter.

## What would change Exact

Nothing here changes it. A bounded alternative would be to run the fill and
state the excess in the deviation, which is what `Fast` is; an `Exact` that
stays exact needs a different sum for overlapping slivers, which is not a
change of this size.
