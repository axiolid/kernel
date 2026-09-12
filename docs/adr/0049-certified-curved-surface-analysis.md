# 0049 - Certified curved surface analysis reports coverage, not a verdict

## Status

Accepted

## Context

Certified surface/surface intersection was restricted to single-span
affine patches. That restriction looked like an unfinished shortcut. It
is not.

Transversality was certified with a single bound over the whole patch
pair: `|n1 x n2|^2 > 0`, evaluated with interval arithmetic on the
partial derivatives. On a curved patch the normal sweeps, so its
interval hull straddles the other surface normal and the bound is
**exactly zero** -- measured, on a quarter cylinder against a plane:

```
NORMAL lower=0e0   ->  refused
```

Affine patches escape only because a constant normal makes the interval
a point, so the bound is sharp. Removing the affine gate experimentally
changed nothing: the refusal moved one line earlier. No tolerance or
budget can fix a structurally zero bound.

Subdividing recovers transversality for most cells but not all:

```
n=1     1 cell    0 positive
n=2     4 cells   2 positive
n=4    16 cells  12 positive
n=8    64 cells  56 positive
n=16  256 cells 240 positive
```

The zero cells number exactly `n` -- a one-dimensional band that never
shrinks away. Inspecting the partials along that band shows why:

```
i=0 sep=784  du_y=[1.750,2.000]
i=6 sep=16   du_y=[0.250,0.500]
i=7 sep=0    du_y=[-0.000,0.250]   <- tangent
```

At the cylinder crest the u-tangent goes horizontal and the normal turns
parallel to the plane normal. That is real tangency, not a weak bound.
Subdividing it forever would be a bug, not diligence.

## Decision

Certify transversality **per cell**, and return **every** cell with what
was proven about it, rather than collapsing the pair to one verdict.

A curved pair is routinely *part* provable. Collapsing that to
"resolved" or "unresolved" is what loses geometry: the unprovable part
either disappears, or it poisons an otherwise usable result. So
`certify_surface_arcs` returns regions, each labelled:

| Kind | Meaning | Can more budget help? |
|------|---------|----------------------|
| `Empty` | Proven disjoint | n/a, already proven |
| `Transversal` | Proven regular, positive separation bound | n/a, already proven |
| `Tangential` | Geometry obstructs certification | **No** |
| `BudgetExhausted` | Policy stopped subdivision | **Yes** |

Separating the last two matters: one says "raise the depth", the other
says "raising the depth is pointless". Folding them into a single
`Unresolved` would make a caller retry work that provably cannot succeed.

`Tangential` is assigned by a probe (do any children recover a positive
bound?) and is therefore a **labelling heuristic, not a certificate**.
Both labels mean unproven; only the advice differs. The code says so.

## The honesty mechanism

The failure that matters here is invisible: a result that quietly loses
a region still looks plausible -- sensible kinds, sensible depths, no
panic, no error. Counting arcs cannot detect it.

So completeness is enforced structurally. Subdivision bisects all four
parameter axes, so a leaf at depth `d` covers exactly `16^-d` of its
root box -- a dyadic rational, never a rounded quantity. `audit_coverage`
sums `16^(MAX_AUDIT_DEPTH - d)` over the leaves in `u128` and requires
the total to equal `16^MAX_AUDIT_DEPTH` per root pair:

- too small -> `Gap`, geometry was dropped
- too large -> `Overlap`, geometry was double-counted

Integer arithmetic throughout. No tolerance, no accumulated float error,
no judgement call. `MAX_AUDIT_DEPTH = 31` because `4 * 31 <= 127` keeps
the sum exact in `u128`; a larger bound would saturate, and a saturated
total could mask a real gap. Depths beyond it are refused, not clamped.

### Why this check outlives the code it checks

`audit_coverage` constrains only *completeness*, never *content*. It
does not encode what the answer should be, so sharpening a bound,
adding a refusal kind, or changing how tangency is handled does not
invalidate it -- any such change still has to account for the domain.

That is deliberate. A test asserting "this pair yields 12 transversal
regions" becomes a maintenance burden the first time a bound improves,
and gets updated until it asserts nothing. A test asserting "the regions
tile the domain" keeps failing on exactly the bug that would otherwise
ship silently.

Mutation evidence that the mechanism works (4/4 killed):

| Mutation | Result |
|----------|--------|
| Silently drop a tangential leaf | KILLED |
| Emit a parent alongside its children | KILLED |
| Label a region transversal without positive evidence | KILLED |
| Skip the empty check, losing that leaf | KILLED |

## Consequences

**Good.** Curved pairs are analysed instead of refused wholesale. The
unprovable part is located rather than lost. Any future change to this
pipeline is checked against domain accounting by construction.

**Cost.** The region count grows as `16^depth` in the worst case, so
depth is explicit policy rather than something the algorithm chooses.
Callers must branch on `is_fully_certified()`; a caller that reads only
the transversal regions gets an incomplete picture. That is intentional
-- the type makes the incompleteness visible instead of implying a
totality that does not hold.

**Not done here.** This certifies the *structure* of a curved
intersection -- where a regular curve provably exists, and where nothing
can yet be claimed. It does not yet trace the arcs, stitch branches, or
build trimmed B-rep faces from them. Those need curved edge geometry
throughout the assembly path, which today assumes straight edges
(`add_line_edge`, `Interval::UNIT` pcurves). The affine path remains the
exact fast case, unchanged.

## References

- `crates/algorithms/parametric/nurbs/src/certified_surface_arcs.rs`
- `crates/algorithms/parametric/nurbs/tests/curved_arcs.rs`
- `docs/plans/certified-boundary-roots.md`

## Exact elementary curves (addendum)

The analysis above certifies *where* a regular intersection curve exists. It
does not produce one. Producing a curve generally means marching sampled
points and fitting a spline through them, which is an approximation carrying
an error bound.

`ExactBRep` documents itself as holding "exact 3D curve supports", and
ADR 0046 requires that every decision be an exact predicate. A fitted spline
placed in that type would make its central claim false, so fitting is not
used.

Instead, the surface pairs whose intersection has a **closed-form** answer
are derived symbolically in `exact_surface_intersection.rs`:

| pair | condition | curve | identity |
|---|---|---|---|
| plane / plane | normals not parallel | `Line3` | direction `n1 x n2` |
| cylinder / plane | normal parallel to axis | `Circle3` | radius `r` |
| cylinder / plane | oblique, not axis-parallel | `Ellipse3` | semi-axes `r`, `r / cos(theta)` |
| sphere / plane | `\|d\| < r` | `Circle3` | radius `sqrt(r^2 - d^2)` |

Every other pair returns a typed refusal. `UnsupportedPair` says only that
*this module* does not derive the case — sphere/sphere, for instance, is a
circle in principle but is not implemented, and saying so is preferable to
fitting one and calling it exact.

Degenerate configurations refuse rather than returning a degenerate curve:

- parallel planes -> `Disjoint` (coincident or never meeting)
- tangent sphere/plane -> `NotRegularCurve` (a point, not a curve)
- axis-parallel cylinder/plane -> `NotRegularCurve` (two lines, or none)

### How the algebra is checked

Asserting "an ellipse came back" would pass on a *wrong* ellipse. The tests
instead substitute sampled curve points into both operands' own defining
equations and require residuals below `1e-12`. That catches a wrong frame,
radius, or semi-axis, which a variant check cannot.

Mutation testing confirms the checks bite: stretching the ellipse by
`cos(theta)` instead of `1/cos(theta)`, using `r^2 + d^2` for the sphere
section, centring the cylinder section on the axis origin rather than the
plane, returning a zero-radius circle for a tangent plane, and dropping the
axis-parallel refusal were each introduced deliberately and each failed the
suite.

### Scope

These curves are exact and are real B-rep-ready geometry, but they cover
elementary analytic surfaces only. A spline-surface pair still yields
regions, not curves, and still refuses rather than approximating.

