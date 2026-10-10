# 0070 — Exact arc booleans replace cavalier_contours

- **Status:** Accepted
- **Date:** 2026-09-24
- **Deciders:** Friedrich, axiolid
- **Relates to:** [0050](0050-arc-aware-planar-overlay.md), [0068](0068-exact-number-type.md), [0069](0069-arc-boolean-path.md); issue #155
- **Supersedes:** the "keep cavalier_contours" half of [0069](0069-arc-boolean-path.md)

## Context

ADR 0069 kept `cavalier_contours` behind `arc_overlay`, rescaled so its
fixed thresholds track the caller's tolerance, and planned an exact path once
an exact number type existed. `axiolid-exact` now has one: interval filter,
dyadic exact fallback, nested square roots (`Tower`), exact conics.

A boolean with arcs needs points where a line meets a circle and where two
circles meet. Those carry one square root each, so every question the
overlay asks (which side, which comes first along an edge, is this point on
that curve, which way does the boundary turn here) is the sign of a
polynomial in values of the form `(a + b*sqrt(d)) / w`, with dyadic `a`,
`b`, `d`, `w`. `axiolid-exact` decides such signs exactly.

## Decision

`arc_overlay` runs on an in-tree exact core (`axiolid-overlay`,
`src/exact_arc/`). `cavalier_contours` is removed from the workspace. The
public API is unchanged.

1. **Every topological decision is exact** for the given `f64` input:
   crossings, their order along each edge, whether a piece of one boundary
   is inside, outside or on the other region, linking into rings, and hole
   ownership. No decision reads the tolerance, so unit independence is
   structural (the three unit tests of ADR 0069 pass with the scaling code
   deleted).
2. **Output is rounded once.** Crossing points are irrational in general;
   output vertices are the nearest doubles (about 50 bits). A crossing that
   lies within rounding of an input vertex would leave an edge shorter than
   the tolerance after rounding; such edges are merged, and slivers below
   the tolerance are dropped. This is the only use of the tolerance on a
   result, and it is presentation, not topology.
3. **Filter before exact.** Each point caches a sound `f64` box around its
   coordinates. Predicates are decided on the boxes when possible, then by
   the `axiolid-exact` interval tier, then exactly. Point identity also tries
   plain dyadic cross-multiplication before building a tower.

### Algorithm

Split every edge at every point where it meets the other operand
(overlapping pieces contribute their end points, so they split
identically). Take an exact rational sample strictly inside each piece
(half-angle parametrization of the arc; seeded in `f64`, verified exactly).
Classify it: on the other boundary (shared, same or opposite direction), or
inside/outside by a crossing count over y-monotone parts with the
half-open rule, which needs no special case at vertices or tangencies. Keep
pieces by the operation's table, reversing those that bound the result from
the other side. Link at exact vertices, taking the leftmost turn, which
traces minimal rings. Counter-clockwise rings are outers; each clockwise
ring is a hole of the smallest outer containing it.

## Evidence

- `tests/arc_exact_oracle.rs`: 150 grid-snapped scenes (shared edges,
  identical, tangent and concentric circles, vertices on the other
  boundary) and 120 scenes on a decimal grid that binary cannot hold, each
  under all four operations. Two oracles that never call the code under
  test: area identities (`|A u B| + |A n B| = |A| + |B|`, and so on, against
  closed-form input areas) and point membership against finely tessellated
  operands, over 400,000 checks. Every output ring passes
  `validate_arc_ring`.
- `scripts/probe_arc_overlay_mutants.py`: eight deliberate faults (turn
  rank, tie-break, half-open rule, rounding cleanup, shared pieces,
  arc-side test, point-identity shortcut, difference table); each must fail
  the suite.
- All previous `arc_overlay` tests, including ADR 0069's unit tests, pass
  unchanged, as do the downstream `axiolid-construct`, `axiolid-project` and
  `axiolid-route` suites.

## Consequences

- **Cost.** Measured with `benches/arc_overlay.rs` (release, this machine):
  90 to 150 us per boolean on typical sections (a disc against a rectangle,
  two discs, a round opening in a rounded wall), about 1 ms for a disc
  passing through the corners of a wall end (crossings within rounding of
  vertices escalate to exact arithmetic). `cavalier_contours` took about
  1 us on the same scenes. With `n` arc edges against a disc: 8 edges
  0.5 ms, 32 edges 1.3 ms, 128 edges 4.2 ms, 512 edges 15 ms. Exactness is
  paid for here; the cheap float path it replaces was not correct, and the
  cost is local to the arc path (straight-edge polygons keep their integer
  backend).
- **Dependencies.** `cavalier_contours`, `static_aabb2d_index` and
  `smallvec` leave `Cargo.lock`. `axiolid-overlay` now depends on
  `axiolid-exact` and `axiolid-guarantees`.
- **Known limits.** Self-intersecting operand rings stay outside the
  contract; the linking step reports `SelfIntersection` when it meets one
  as a dead end, but does not promise to detect every such input. Two
  result pieces leaving one vertex along the same tangent (curves touching
  tangentially exactly at a result vertex) are ranked equal; rings stay
  closed, but how touching rings are grouped is not specified.
- **Broad phase (added after landing).** Each edge carries a padded
  `f64` bounding box: the chord's box grown by the sagitta
  `|bulge| * |chord| / 2`, which holds any arc, minor or major. Edge pairs
  whose boxes are apart skip the exact crossing test, shared-edge checks
  skip edges whose box cannot hold the sample, and linking looks up
  candidate pieces in a list sorted by the lower `x` bound of their start
  enclosure. Boxes only skip work; every decision is still exact. Measured
  (release, alternating runs, two rounds):

  | Scene | Before | After |
  | --- | ---: | ---: |
  | Disc crossing a rectangle | 152 us | 100 us |
  | Round opening in a rounded wall | 120-133 us | 69-76 us |
  | Disc through a wall end's corners | 1.02 ms | 0.75 ms |
  | 512 arc edges vs a disc | 16.6 ms | 8.0 ms |
  | Two 64-edge wavy rings | 53 ms | 1.9 ms |
  | Two 256-edge wavy rings | 779 ms | 6.4-7.0 ms |
  | Two 1024-edge wavy rings | (not run) | 32-33 ms |
  | 4096-edge comb vs a local disc | 51-55 ms | 17 ms |

  Growth is now close to linear in edge count; what remains per call is
  per-edge setup (exact circle coefficients, monotone splitting, piece
  samples), not pair tests.
- **Performance work left.** The tower allocates per operation; a
  small-vector coefficient store would cut most of that. Per-edge exact
  setup could be deferred to edges whose box meets the other operand.

## Addendum: straight edges on the same core (#173)

ADR 0050 left the polygon path on `i_overlay`. That backend maps the
operands' bounding box onto an `i32` grid and back, so every output
coordinate, including an input vertex the operation does not move, came
back snapped to a step of about 1.5e-8 of the extent: `[0,4]x[0,0.2]`
clipped by a box around it had area `0.800000011920929`. Consumers report
such areas as exact evidence, and threshold rules compare them with
limits.

`overlay`, `union_soup` and `Region` now run on this core: a straight ring
is an arc ring without arcs. They use the multi-ring subdivision
(`exact_arc/arrangement.rs`, built for #120) rather than the two-operand
boolean, since an operand may hold several overlapping rings with holes
and a fill rule:

1. Cut the plane by every ring of both operands at once.
2. For each piece, count each operand's winding number on either side:
   the rings containing that side, `+1` for a ring given counter-clockwise
   and `-1` for clockwise, as the integer backend counted them. The fill
   rule reads the count, the operation combines the operands.
3. Keep pieces with the result on exactly one side, facing it; link; drop
   vertices where the boundary runs exactly straight on.
4. Round once. Rational points (input vertices and segment crossings) are
   now correctly rounded (`XPoint::rounded`), in the arc path too.

The subdivision was quadratic in the ring count (every edge scanned every
ring). It now indexes edge boxes, monotone parts and ring boxes in static
box trees, and asks each edge pair once. Pairs of input segments go
through certified `f64` filters first (apart, identical, sharing one end);
points built from input vertices carry exact point boxes, so their
identity needs no arithmetic. Measured against the integer backend
(release, this machine, `union` of the same inputs):

| Scene | `i_overlay` | Exact |
| --- | ---: | ---: |
| Two 100-gons | 0.16 ms | 2.4 ms |
| Two 1000-gons | 9.2 ms | 31 ms |
| Two 5000-gons (both: O(n^2) operand validation) | 221 ms | 330 ms |
| Soup of 1600 triangles (`union_soup`) | 0.25 ms | 36 ms |
| Soup of 14,400 triangles | 3.5 ms | 0.40 s |

The route, project, construct and B-rep suites run in the same time as
before, and `benches/arc_overlay.rs` is unchanged within the noise of a
shared machine (alternating runs; the arc scenes were at or below the
previous figures in the quieter round). Evidence: `src/exact_overlay.rs` (point membership against the
integer backend over 256 random scenes, all four operations and fill
rules), `tests/exact_straight.rs` (bit-identical vertices, correctly
rounded crossings; all four tests fail on the old backend), and
`scripts/probe_exact_overlay_mutants.py`. `i_overlay` remains for
`offset.rs` only.

Note (#292): `ArcArrangement` now keeps, for each piece, the rings that
hold both of its sides (sparse, as the subdivision computes them) and its
carriers with the side each holds, not a flag per ring per side. Which
faces a predicate selects is unchanged: `regions` fills one buffer of
flags per side for its predicate, and the crate's Minkowski operations
read each side from the ascending list of rings that hold it. The flags
made memory, and the Minkowski face selection's time, grow as pieces
times rings; a sum reads one ring per boundary edge, so a 3,400-corner
floor needed more than 12 GB.
