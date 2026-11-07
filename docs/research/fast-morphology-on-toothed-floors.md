# Fast disc erosion of a floor with sub-budget teeth

Measured evidence for #345 (consumer: a downstream rule engine). `Region::erode_inner_with`
in `MorphologyMode::Fast` took seconds and gigabytes on a toothed floor
whose every tooth is shallower than the budget. The cause is a split policy,
not a missing rule; the fix is in `overlay/src/fill.rs` and keeps ADR 0087's
promises (side, reported excess within the budget, `Exact` unchanged).

## Shape and method

`tools/benchmark/benches/disc_erosion.rs`, walls `floor`: #345's room, a
3 m square counter-clockwise, its bottom and top sides carrying teeth
0.5 mm deep into the room, `n + 1` corners on each (the teeth `3 / n` apart,
`n = corners / 2 - 1`). The budget is 1 mm, so every tooth is within it.
An inner erosion by 1.4 m, 1.49 m or 1.6 m (the last empties the floor).
Each case is one process (`DISC_EROSION_CASE`); `Exact` does not finish on
any of them and is not in the sweep. The machine was shared and busy, so
the absolute times are indicative; the growth is the finding. "Before" is
main at 3c8bd39a (`DISC_EROSION_RUNS=1`; the cases that did not finish
within the 20 GB address-space cap of that run are marked); "after" is the
median of 5 runs. Peak memory is the child's `VmHWM`.

## Diagnosis

A copy of `fill`'s loop that counted why each chord was refused (not kept:
the counts are in the table) and `fill::diagnosis` (an `#[ignore]` test that
prints the chords, the filled ring and the excess per size and radius, run
on the fixed code the floors' erosions take 3 or 4 chords and leave 6 to 9
corners, 8 at 50,000 corners and r = 1.4, and the dilations 2 chords and 6) showed, on the complement of the floor
(an erosion sums it), at `t` half the radius and a 1 mm budget:

| corners | r (m) | chords used | filled ring (corners) | refused: unproved / over budget / ring inside |
|--------:|------:|------------:|----------------------:|-----------------------------------------------:|
| 50,000 | 1.4 | 3 | 5,891 | 5,885 / 1 / 0 |
| 50,000 | 1.49 | 1,140 | 5,138 | 5,133 / 1 / 0 |
| 100,000 | 1.4 | 4,311 | 7,467 | 7,462 / 1 / 0 |
| 100,000 | 1.49 | 4,661 | 7,889 | 7,884 / 1 / 0 |

- **A fill failure is not the cause in general; one kind of refusal is.** No
  chord had `encloses_none` false, the filled ring was never zero in area,
  and `excess` never returned `None`. A chord is proved or it is not: the
  failures are all "not proved", after one chord over budget.
- **The room-corner chord is within the budget, and does not stop the
  policy.** The chord over the room's corner at the end of the bottom side
  has an excess of one tooth step, `3 / n` (1.2e-4 at 50,000 corners,
  6.0e-5 at 100,000: the figure in the issue) and is used. It is the
  reported excess, but not the cause.
- **The complement's outer frame is not involved.** `fill` is given the
  region's own rings only; a floor has one ring, and the complement's
  frame never enters.
- **The cause.** `cover::proposals` gives the bottom side as one chord
  along the tips, then a second chord from the last bottom tip over the
  room's corner towards a point of the top side, 0.69 m deep at `t` = 0.7.
  It is over the budget and is split at its deepest vertex, the corner.
  The part left runs from the corner, a valley of the top wall (at wall
  level, outside the room), along the wall to a point `Q`. A chord from a
  valley is not outside the set: the tips between lie on its room side, so
  it is not proved. The split for a chord not proved is at the vertex
  nearest the line, which for a chord along valleys is a valley on it
  (distance 0), the middle one. Both halves are again chords between
  valleys, and again not proved, down to single edges: one used chord per
  tooth, or none. Everything between the corner and `Q` stays unfilled:
  5,900 corners at 50,000, 7,500 at 100,000 (the wall's share within `t`
  of the corner, about an eighth of the corners), and the sum of that
  many 0.5 mm teeth 1.2e-4 m apart, each sliver overlapping thousands of
  the others, has the quadratic cost of #330. The parity of the teeth
  decides whether `Q` is a valley or a tip, which is why the cost jumped
  between neighbouring sizes and radii (the issue's 1.8 s at 99k corners
  against 7 s at 200k).
- **Why it grew with #340.** #330's proposals reached `t / 16`, so the
  chord over the corner was 0.044 m deep at most and the corner's
  neighbourhood was small. With #340's second round at `t` the chord over
  the corner is as deep as the disc allows.

## The fix

A chord that is not proved is tried without its start edge, then without
its end edge, then without both, before it is split as before. The first of the three that is
proved, within the disc and the budget, and with no ring in its notch, is
used; none changes the proof, the excess or the check for a ring. The
chord's ends are then the wall's tips, whichever the parity of the teeth.
A first version trimmed only the chords a too deep notch's split left;
trimming every chord not proved gives the same hashes on all 32 `fast`
cases of the bench's #330 and #340 sweeps (before and after) and is no
slower (the 16 cases of 10,002 corners sum to 0.89 s against 1.5 s on the
busy machine), so it is the rule: simpler, and no mutant of the
restriction survives.

The result's side and deviation are those of ADR 0087: the excess of a
trimmed chord is computed as any other's, and a chord over the budget is
still refused. The result of the floor's erosion is the same region to the
bit (the hash of every row below is the same before and after, where the
before ran): only the arrangement's input changed, from thousands of
corners to six.

## Results

Release build, inner erosion, `fast(1e-3)`. Time s / peak MB. The excess is
that of the room-corner chord, `3 / n` for `n` teeth steps up to the
budget, the same before and after.

| corners | r = 1.4 before | r = 1.4 after | r = 1.49 before | r = 1.49 after | r = 1.6 before | r = 1.6 after | excess |
|--------:|---------------:|--------------:|----------------:|---------------:|---------------:|--------------:|-------:|
| 10,000 | 13.2 / 860 | 0.078 / 6 | 2.69 / 995 | 0.075 / 6 | 3.16 / 1,193 | 0.079 / 6 | 6.0e-4 |
| 20,000 | 70.8 / 4,817 | 0.157 / 8 | 8.28 / 5,639 | 0.157 / 8 | 14.4 / 6,830 | 0.157 / 8 | 3.0e-4 |
| 30,000 | 56.0 / 4,194 | 0.228 / 9 | 55.0 / 16,480 | 0.234 / 10 | over the cap | 0.236 / 10 | 2.0e-4 |
| 50,000 | over the cap | 0.386 / 15 | over the cap | 0.370 / 14 | over the cap | 0.383 / 14 | 1.2e-4 |
| 100,000 | not run | 0.769 / 25 | not run | 0.775 / 25 | not run | 0.740 / 25 | 6.0e-5 |
| 150,000 | not run | 1.098 / 40 | not run | 1.106 / 38 | not run | 1.121 / 38 | 4.0e-5 |
| 200,000 | not run | 1.570 / 48 | not run | 1.500 / 47 | not run | 1.476 / 48 | 3.0e-5 |

Growth is linear from 50,000 to 200,000 corners (4.1 times the corners, 4.1
times the time; 3.2 times the memory). A first attempt at 50,000 corners
before the fix reached 26 GB in four minutes and was stopped. What is left
is linear work in the ring's size (validation, the index, the proposals):
about 7.7 microseconds a corner. The excess stays within the budget at
every size and is reported; below about 6,000 corners the corner chord's
step exceeds the budget, is refused, and a few teeth at the corner stay.

`the_fast_mode_on_a_floor_of_sub_budget_teeth_keeps_its_side` in
`tests/minkowski.rs` checks the side against the true
morphology (the point test of the definition, not the exact route, which
does not finish) at 10,000 and 10,002 corners, for all four operations.
