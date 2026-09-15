# Relations over natural-equation space curves

Status: accepted
Date: 2026-09-15

Follows ADR 0061, which added `Curve3::Intrinsic` as a representable value
with an evaluator and deliberately stopped there.

## Context

ADR 0061 left four things undone: composition, trimming, offsetting, and
format lowering. They are not one feature. Each has a different answer,
and two of them are not implementable as stated.

The governing fact is that an `Intrinsic3` stores its laws against ARC
LENGTH. That is what makes some relations exact and others impossible.

## Decision

### Trimming is exact

Restricting to `[a, b]` does not approximate the shape. It re-anchors the
same law: the trimmed curve's law is `g(u) = k(a + u)`, and the family is
closed under that shift.

- a polynomial shifts by the binomial expansion of `(a + u)^i`;
- a sinusoid shifts purely in PHASE, `p -> p + w * a`;
- a piecewise law drops spent pieces, rebases the one containing `a`, and
  moves the remaining seams back.

`CurvatureLaw::shifted` implements this in closed form. The only
quadrature involved is re-anchoring the START FRAME, which needs the
position and frame at `a` -- and that is the evaluator's existing job, not
a new approximation.

### Offsetting is exact for a helix and REFUSED otherwise

This is a representational limit, not a missing feature. The normal
offset `q(s) = p(s) + d * N(s)` has speed `|1 - d*k(s)|`, so it is only
unit-speed when `k` is constant. An `Intrinsic3` is parameterised by arc
length BY CONSTRUCTION, so for a varying law the offset is not an
`Intrinsic3` at all -- storing one would require refitting a law to a
curve that has none, which is a lossy guess wearing an exact type.

Measured: on `k(s) = 0.08 + 0.01 s` at `d = 2`, the offset speed ranges
over `0.724 .. 0.844` rather than staying at 1.

For a helix the offset IS another helix, with radius `a - d`, the same
pitch, and arc length rescaled by `c2/c`. Its start frame is NOT the
base's: the offset point keeps the angular rate but changes radius, so
the tangential/axial mix of the tangent changes. Verified against the
geometric offset to 1.7e-10 at three distances including a negative one.

### Composition is a piecewise law, guarded

Joining two curves is exactly what `CurvatureLaw::Piecewise` already
models: one absolute start frame, interior anchored purely by arc length.
`join_intrinsic3` refuses unless the second curve starts where the first
ends AND their tangents agree, because a join that jumps or kinks is not
one curve.

### Graph relations come free, and that is the point

`CurveRelation::Trimmed` and `Composite` are GENERIC over curve family.
They were never the blocker: `evaluate3` refused `Intrinsic` by name, and
every generic consumer inherited that refusal. Dispatching the family in
`evaluate3`/`derivative3`/`domain3` makes the existing relation machinery,
including sweep directrix sampling, work on a torsion curve with NO change
to the relation code.

`domain3` reports `[0, length]`, not the unit interval, because the
parameter IS arc length. A consumer trimming by parameter range depends on
that being right.

### Format lowering is out of scope, and not silently

There is no format code in this repository -- no IFC, no STEP. Lowering
belongs to a consumer that owns a format, and inventing a lowering here
with no format to lower from would be speculative. What this ADR owes such
a consumer is a value that survives the graph and evaluates, which is now
the case.

## Two bugs this surfaced

**Quadrature panels straddled curvature seams.** Gauss-Legendre assumes a
smooth integrand across a panel; a piecewise law is only piecewise-smooth.
A joined curve was wrong by 3.0e-3 until panels were split at the seams of
either law. `CurvatureLaw::seams_within` reports them, including nested
ones, and both integrators now break there.

**Piecewise laws could not be evaluated over a partial span.** The
integrator refused any seam beyond the requested `s`, which made every
intermediate evaluation of a joined curve fail. Integrating `[0, s]`
legitimately stops inside whichever piece contains `s`. The stricter rule
still applies to `total_turning`, which asks about the DECLARED length and
must refuse a law whose pieces do not tile it -- a distinction that now
lives at the caller rather than in the shared integral.

## Consequences

- Trimming, joining, and helix offsetting are exact in the law.
- Offsetting a varying-curvature space curve is refused by name.
- A torsion curve works as a graph directrix, trimmed or composited,
  through machinery that does not know the family exists.
- Not done: `CurveRelation::Offset` still resolves through the generic
  point-sampling path rather than the exact helix offset, because the
  relation carries a `reference_direction` whose interaction with a Frenet
  normal is a separate decision. Named here so it is not mistaken for
  oversight.

## Verification

- 14 relation tests, 4 dispatch tests, 3 end-to-end graph tests.
- 1618 workspace tests passing, gate green.
- Trim exactness is pinned SYMBOLICALLY, comparing closed-form integrals of
  the shifted and base laws to 1e-12 -- not through the evaluator, where
  two different panel layouts limit agreement to ~1e-7. Asserting machine
  precision between two quadratures would be asserting something false.
- Helix offset checked against the geometric offset of the base curve,
  computed independently from the Darboux axis, at d = +0.5, +1.0, -0.75.
- Mutation testing: 11 mutants, 11 killed. One initially survived --
  deleting the meeting-point check in `join_intrinsic3` -- because the
  fixture also had a wrong tangent, so the kink check masked it. The
  fixture now keeps the tangent correct and varies only position.

