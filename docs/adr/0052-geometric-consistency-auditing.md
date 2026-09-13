# Geometric consistency auditing

Status: accepted
Date: 2026-09-12

## Context

`audit_brep` is pure topology: handles and adjacency, no coordinates and no
tolerance. That makes it exact and reproducible, and it should stay that way.

It also means a whole class of defect passes silently. Every variant of
`ExactBRepError` is about PRESENCE or handle RESOLUTION -- is there a pcurve,
does the handle resolve, is the interval finite. None asks whether the stored
geometry AGREES with itself.

Measured: replacing every pcurve of a box with `Line2 { origin: (999, -999) }`
leaves `audit_brep` reporting `tessellable = true, closed = true`. With
intervals supplied, the builder accepts it.

This is not hypothetical. ADR 0050 step 5 shipped cap loops whose pcurves were
straight chords across arc edges. The solid closed, validated, audited clean,
and reported a plausible area, while the cap boundary disagreed with the wall
boundary along the same edge. It was found by reasoning, not by any check.

## Decision

Add `axiolid-brep-audit`, a separate crate providing `geometric_audit(brep,
tolerance) -> GeometricHealth`. It EVALUATES:

- every edge's start and end vertex lies on that edge's own 3D curve;
- every pcurve, lifted through its face surface, follows the 3D curve of the
  edge it trims.

It stays separate from `audit_brep` because it is necessarily
tolerance-dependent. A clean geometric audit is a statement about agreement AT
A TOLERANCE, not an absolute one, and mixing the two would weaken the exact
topological guarantee.

Both exact boolean paths are gated through it. A boolean is where pcurves are
rebuilt against surfaces they did not originally trim, so it is the operation
most able to produce a solid that is topologically perfect and geometrically
wrong.

## What the evaluation trait question turned out to be

The intended first step was "add an evaluation trait for Curve3/Curve2/
Surface". That step was unnecessary: `CurveEvaluator` and `SurfaceEvaluator`
already exist as capability contracts, and `axiolid-evaluate` already provides
2769 lines of implementation including `invert` and `project`. The audit
consumes those rather than adding a parallel evaluation path.

## Consequences

- Sampling is at five parameters per edge use: both endpoints plus three
  interior. Endpoints alone are insufficient -- a chord standing in for an arc
  agrees EXACTLY at both ends and deviates only between them.
- The worst sampled deviation is reported, not the first one over tolerance.
  An under-reported error invites someone to widen the tolerance just past it.
- The audit is cheap on analytic solids: an exact boolean returns a handful of
  faces, so this is a few evaluations per edge use, not a mesh traversal.
- `axiolid-brep-audit` cannot depend on `axiolid-construct`, even for tests:
  that is a cycle. Integration tests that need constructors live in
  `construct`; the hand-built negative test lives in `brep-audit`.

## Findings

### The audit found a real bug on its first run

Before any deliberate negative test, `geometric_audit` rejected a solid built
by the SHIPPED fillet path:

    PcurveOffCurve { loop_id: 2, use_index: 0, error: 0.06561637489641778 }
    PcurveOffCurve { loop_id: 5, use_index: 0, error: 0.0874884998618904 }

Two independent defects, both in `extrude_with_cylindrical_blends`:

1. `add_polygon_ring` gave EVERY ring edge a `Curve3::Line`, including blend
   edges. The blend face then attached an arc pcurve to that straight edge.
2. `add_cap_loop` gave every cap edge a `Line2` chord pcurve, including over
   blend arcs -- the same defect ADR 0050 hit, in a different constructor.

The error ratio confirms the diagnosis rather than assuming it:
`0.0656 / 0.0875 = 0.75 = 0.3 / 0.4`, exactly the ratio of the two fillet
radii. The deviation is the chord-versus-arc sagitta, scaling linearly with
radius as it must.

Fixed by `add_polygon_ring_with_blends`, which gives blend edges a `Circle3`,
and `add_cap_loop_with_blends`, which gives them a `Circle2` pcurve. Both
constant-radius fillet solids now audit clean.

This is the whole argument for the crate: the bug was in code that passed
every existing test, the full gate, and review.

### Endpoint-only sampling proves nothing

A chord standing in for an arc agrees EXACTLY at both endpoints. Sampling only
there reports zero error on precisely the defect this exists to catch. The
mutant that reduces `SAMPLES` to endpoints is killed by the chord test.

### Report the worst sample, not the first failure

The first draft returned on the first sample over tolerance. For a quarter arc
of radius 1 that reported `0.2187` (the quarter-point deviation) when the true
worst is `0.2929` at the midpoint -- the classical sagitta `r(1 - cos(45 deg))`,
verified independently. Under-reporting invites widening the tolerance to just
past the visible number.

### A gate that cannot be tripped is not a gate

The first version of the gate test wrapped its assertion in `if let Ok(...)`,
so disabling the gate entirely left it passing. Two mutants survived because
of it. Measuring what actually happens at a 1e-18 tolerance -- 24 defects,
worst deviation 3.14e-16 -- gave a concrete rejection to assert instead.

Mutation, 6/6 killed: audit never reports a defect; endpoint-only sampling;
first-failure instead of worst; gate disabled; gate not called on the arc
path; gate reports no measurement.
