# Minkowski sums read from the term edges that can bound them (#344)

`Region::minkowski_sum`, `minkowski_erosion` and the disc morphology put
every term ring into one arrangement and read each face against every hull
holding it. With a shape larger than the region's features every hull
overlaps nearly every other, so the arrangement had many crossings and each
piece lay in many hulls. ADR 0108 reads the union from the term edges that
can bound it instead (`src/term_union.rs`): an edge proved, exactly, to lie
strictly inside one convex term is left out, the rest are cut once by the
exact core's partial mode and read by the same rule. Results are unchanged
bit for bit.

## Method

Release build, each case in its own process, the median of three runs, base
(main at aa73f556, `by_arrangement` route) and branch binaries run one after
the other on the same machine in separate target directories. The machine was
shared, so absolute times carry noise of tens of percent; the ratios between
the two columns of one run are the reliable part. `hash` is the FNV-1a digest
of every vertex of the result (`tests/minkowski.rs`'s `digest`), the same
before and after in every row. All baseline hashes, the `turned` rows
included, were printed by the baseline binary itself (the bench built on
aa73f556, before any library change); the message of `74932788` lists only
the `axis` ones, and `turned` rows were compared against the same run's
output file, not a number written down by hand.

## `minkowski_furniture`

Groups of a C-shaped sofa, a round table and an L-shaped desk (152 corners
each) on a 4.0 by 2.4 m grid, turned 2.3 degrees and placed at
(640 000, 5 650 000), with a 1.7 by 0.9 m rectangle, axis-aligned (`axis`) and
turned 30 degrees (`turned`); the erosion takes the free floor.

| op | shape | corners | result | before s | before MB | after s | after MB | speed-up | hash |
|---|---|--:|--:|--:|--:|--:|--:|--:|---|
| sum | axis | 152 | 60 | 0.122 | 45 | 0.004 | 5 | 30x | `f7e756025fec9c46` |
| sum | axis | 304 | 97 | 0.264 | 89 | 0.008 | 6 | 33x | `6da71d04fa14a417` |
| sum | axis | 1064 | 236 | 0.996 | 311 | 0.028 | 12 | 36x | `5e006cc15a1dd61d` |
| sum | axis | 3040 | 370 | 3.052 | 941 | 0.084 | 26 | 36x | `5ea7fb6717fa17ff` |
| sum | turned | 152 | 61 | 0.117 | 47 | 0.004 | 5 | 29x | `c94d293f0929a9d1` |
| sum | turned | 304 | 100 | 0.261 | 94 | 0.008 | 6 | 33x | `ff66819f2ff1a3c5` |
| sum | turned | 1064 | 245 | 1.012 | 338 | 0.028 | 12 | 36x | `e8e90f47fb67ce08` |
| sum | turned | 3040 | 301 | 3.312 | 1049 | 0.079 | 26 | 42x | `2afff662ff26bcf2` |
| erosion | axis | 152 | 40 | 0.127 | 47 | 0.005 | 5 | 25x | `c6eed60f1abff4f2` |
| erosion | axis | 304 | 77 | 0.284 | 95 | 0.009 | 6 | 32x | `7f3b517ed3a76422` |
| erosion | axis | 1064 | 194 | 1.124 | 322 | 0.031 | 12 | 36x | `8fc11cde44c2dfc3` |
| erosion | axis | 3040 | 281 | 3.305 | 980 | 0.090 | 26 | 37x | `202e319bf146fb31` |
| erosion | turned | 152 | 16 | 0.131 | 51 | 0.005 | 5 | 26x | `e7bf484f4c490c00` |
| erosion | turned | 304 | 29 | 0.275 | 99 | 0.009 | 6 | 31x | `5b17af88ebc0305b` |
| erosion | turned | 1064 | 106 | 1.182 | 356 | 0.034 | 12 | 35x | `b55f4a17690ebe4a` |
| erosion | turned | 3040 | 80 | 3.964 | 1093 | 0.083 | 26 | 48x | `3386a8af1613b235` |

- The branch is about 25 to 45 times faster: 25 to 48 times over the 16 cases
  of the table, 31 to 33 times at 304 corners (bar: 5). The baseline commit's
  own first run (`74932788`'s message, taken under heavier load) gives 23 to
  58 times; the table's paired rerun is the one to quote.
- The time at 3,040 corners over the time at 304 corners is 9.2 to 10.5 after
  (11.6 to 14.4 before), against the bar of 20.
- Peak memory at 3,040 corners falls from 0.94 to 1.09 GB to 26 MB.
- The structural cost test `a_furniture_cluster_is_read_from_few_pieces`
  counts the pieces cut for the 304-corner scene: 139,142 before, 1,700 after.

## `minkowski_plan` (#292)

Rooms and walls of a floor by a 0.9 m square. No case is slower; the ribs are
the case where many hulls overlap.

| finish | op | corners | before s | after s | speed-up | hash |
|---|---|--:|--:|--:|--:|---|
| pilasters | sum | 124 | 0.012 | 0.004 | 3.0x | `0c92728099aac85c` |
| pilasters | sum | 472 | 0.050 | 0.016 | 3.1x | `1c9806f37dd425c5` |
| pilasters | sum | 934 | 0.116 | 0.033 | 3.5x | `1f6c7330751b7d81` |
| pilasters | sum | 1918 | 0.263 | 0.067 | 3.9x | `32b6bc7b464be877` |
| pilasters | sum | 3370 | 0.490 | 0.135 | 3.6x | `ce18a89ab0f0f43b` |
| pilasters | erosion | 120 | 0.011 | 0.005 | 2.2x | `e07a9f362aedfb4e` |
| pilasters | erosion | 468 | 0.051 | 0.019 | 2.7x | `2d65bce653c97f5f` |
| pilasters | erosion | 930 | 0.124 | 0.036 | 3.4x | `ea8e6c5d147541ac` |
| pilasters | erosion | 1914 | 0.293 | 0.080 | 3.7x | `a60111ed31c76300` |
| pilasters | erosion | 3366 | 0.573 | 0.141 | 4.1x | `984f4e7eb0d25140` |
| ribs | sum | 124 | 0.036 | 0.007 | 5.1x | `270a6c3a753fdda6` |
| ribs | sum | 488 | 0.174 | 0.030 | 5.8x | `794aa5b0e6f1b420` |
| ribs | sum | 974 | 0.365 | 0.065 | 5.6x | `7b767a82896ccf5d` |
| ribs | sum | 2006 | 0.820 | 0.137 | 6.0x | `55edad2b29c0303e` |
| ribs | sum | 3522 | 1.453 | 0.251 | 5.8x | `31fd1e29cdcb0b5b` |
| ribs | erosion | 120 | 0.039 | 0.010 | 3.9x | `c5d131efc2c24878` |
| ribs | erosion | 484 | 0.183 | 0.034 | 5.4x | `659a9b40b95a67da` |
| ribs | erosion | 970 | 0.392 | 0.071 | 5.5x | `9cd2681b7065ab2c` |
| ribs | erosion | 2002 | 0.828 | 0.144 | 5.8x | `e9549344637c450d` |
| ribs | erosion | 3518 | 1.607 | 0.277 | 5.8x | `d13a80dbc99e8fdc` |

## `disc_erosion` (ADR 0087)

All 110 cases of `DISC_EROSION_MODES=exact,fast` with default sizes give the
same hashes. `Fast` is unchanged within noise: no case is slower than 1.05
times its former median (the 50,000-corner `floor` 0.24 s to 0.22 s). `Exact`,
the larger cases (one run each, so noisier; the 52 s of the `needle` dilation
before is an outlier of that run, its siblings take 5 s):

| walls | op | mode | corners | before s | before MB | after s | after MB |
|---|---|---|--:|--:|--:|--:|--:|
| teeth | erode inner | exact | 1202 | 0.298 | 95 | 0.108 | 25 |
| flat | erode inner | exact | 1202 | 0.016 | 9 | 0.010 | 5 |
| teeth | erode outer | exact | 1202 | 0.271 | 97 | 0.109 | 27 |
| flat | erode outer | exact | 1202 | 0.015 | 9 | 0.009 | 6 |
| teeth | erode inner | exact | 10002 | 3.608 | 1183 | 1.924 | 395 |
| teeth | erode outer | exact | 10002 | 3.973 | 1183 | 1.972 | 396 |
| teeth | dilate inner | exact | 10002 | 7.288 | 1158 | 2.293 | 457 |
| teeth | dilate outer | exact | 10002 | 3.831 | 1158 | 2.442 | 458 |
| hole | erode inner | exact | 10002 | 9.320 | 1197 | 2.225 | 462 |
| hole | erode outer | exact | 10002 | 7.411 | 1199 | 2.213 | 463 |
| hole | dilate inner | exact | 10002 | 8.282 | 1173 | 1.998 | 387 |
| hole | dilate outer | exact | 10002 | 8.365 | 1176 | 1.955 | 383 |
| needle | erode inner | exact | 2002 | 5.013 | 1516 | 2.861 | 156 |
| needle | erode outer | exact | 2002 | 4.965 | 1476 | 2.860 | 158 |
| needle | dilate inner | exact | 2002 | 52.193 | 1237 | 2.835 | 163 |
| needle | dilate outer | exact | 2002 | 5.095 | 1236 | 2.836 | 155 |
| deep | erode inner | exact | 2002 | 5.111 | 1456 | 2.444 | 161 |
| deep | erode outer | exact | 2002 | 5.127 | 1455 | 2.431 | 156 |
| deep | dilate inner | exact | 2002 | 4.648 | 1254 | 2.303 | 155 |
| deep | dilate outer | exact | 2002 | 5.307 | 1253 | 2.293 | 159 |

`Exact` on walls of teeth is still superlinear: the arcs and slivers at every
tip are not inside one convex term, so their edges are kept (ADR 0108,
negative consequences). `Fast` remains the answer for those walls.
