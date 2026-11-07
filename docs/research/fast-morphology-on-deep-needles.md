# Fast morphology on needles deeper than the disc (#428)

`MorphologyMode::Fast` filled the notches of a wall of needles 0.05 m deep
(#340) but kept needles deeper than `t`, half the disc radius less the
margin, because no chord was proposed over them. Release build, the room of
`benches/disc_erosion.rs` with needles 0.25 m deep (`deep`), a disc of
radius 0.4 m (`t` 0.2 m), a budget of 1 mm, one run of the bench each.

## Before

`fill` proposed chords within `t` of their line (cover.rs, `proposals`), and
a needle's tip lies beyond it, so every tip stayed on the stack: 0 chords
and an unchanged ring of 203 corners for a wall of 100 needles 2e-4 wide,
0.25 and 0.5 deep. A chord over the wall would not have helped either:
`excess` bounds a notch deeper than `t` by `e_W` of its whole chord,
2.5e-4 over 100 needles of 2e-4 (2.06 m for a wall of 4.5 m), where one
needle's is 2.5e-8. Through the public API, on #340's room with a column:

| needles | corners | median s |
|---:|---:|---:|
| 124, `Fast` | 562 | 0.42 |
| 124, `Exact` | 562 | 0.37 |
| 500, `Fast` | 2,066 | 5.0 |
| 500, `Exact` | 2,066 | 4.4 |

## After

The rule is in the amendment to ADR 0087 and `src/fill.rs`: proposals reach
needles no wider than `needle_width(t, budget)`, a deep notch is bounded by
its widest piece, a deep chord over the budget is split at its on-line
vertex nearest the middle. The 124 needles give 68 corners (dilations) and
72 (erosions) in 0.005 s, with an excess of 8.2e-4, the slab bound of
needles 4.5 m / 124 apart.

| op, side | corners | result corners | median s | peak MB | excess |
|---|--------:|------:|-------:|----:|--------:|
| erode, inner | 2,002 | 4 | 0.017 | 5 | 5.1e-5 |
| erode, inner | 10,002 | 4 | 0.072 | 6 | 9.0e-4 |
| erode, outer | 2,002 | 4 | 0.017 | 5 | 5.1e-5 |
| erode, outer | 10,002 | 4 | 0.073 | 6 | 9.0e-4 |
| dilate, inner | 2,002 | 68 | 0.015 | 4 | 5.1e-5 |
| dilate, inner | 10,002 | 68 | 0.068 | 6 | 2.0e-6 |
| dilate, outer | 2,002 | 68 | 0.014 | 4 | 5.1e-5 |
| dilate, outer | 10,002 | 68 | 0.069 | 6 | 2.0e-6 |
| `Exact`, erode, inner | 2,002 | 1,650 | 4.78 | 1,456 | - |
| `Exact`, erode, outer | 2,002 | 1,650 | 4.64 | 1,455 | - |
| `Exact`, dilate, inner | 2,002 | 2,066 | 4.67 | 1,249 | - |
| `Exact`, dilate, outer | 2,002 | 2,066 | 4.61 | 1,254 | - |

The erosions' 9.0e-4 at 10,002 corners is the room-corner chord, as for
#340's needles on main (the same two values there), not the needles.

The 30 earlier `fast` cases of the bench (#330, #340, #345) give the same
hashes as main, built into separate target directories.

## Family and side

On 240 seeded needle walls (#345's family style), 43 results changed, 41
with fewer corners, none with a side violation; the total time fell from
46.2 s to 30.1 s. `tests/minkowski.rs` checks the side pointwise in bands
for the 124-needle room and a wall of mixed mouths (2e-4 and 0.12), and a
seeded family of 60 walls stays under 3,000 corners (2,851; main: 4,632).

## Limit

A join across a valley needs a vertex exactly on the chord. A wall rotated
by 30 degrees has valleys collinear only up to rounding: each needle is
filled alone, but the valleys stay. 2,500 needles gave 385 to 456 corners
in 0.07 to 0.14 s against 68 axis-aligned, with the side checks passing.
Cutting at vertices within a small distance of the chord, with the strip
between bounded by the outside gap bound plus `sqrt(delta^2 + g^2 / 4)`, is
the sound extension; it is #483.
