# 0095 — Profiles are validated on their lowered arc rings

- **Status:** Accepted
- **Date:** 2026-11-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #352 (consumer: rheonda/cad#137). `ValidateProfile` in `axiolid-profile`
has no implementation. A consumer that needs to know whether a profile has a
self-contact, a gap or a badly nested hole stopped at a stopgap that chords
every arc to 1 µm and asks the straight overlay. That cannot certify a
tangency, and it is not the geometry the exact extruder builds.

Two facts shape the answer.

- **The authored data has no exact joints.** A `ProfileSegment` stores a
  `Circle2` (an `f64` frame and radius) and an `f64` angle domain, so an arc
  end is `centre + r (cos t, sin t)`. By the Lindemann-Weierstrass theorem
  `cos t` and `sin t` are transcendental for every nonzero dyadic `t`. An
  exact predicate on the literal data answers "open contour" for almost every
  arc-to-arc joint, correctly and uselessly. `axiolid-exact`'s interval
  trigonometry can prove nonzero signs only, never a zero.
- **The exact path already has an interpretation.**
  `contour_lower::contour_to_arc_ring` turns a contour into an `ArcRing`: one
  `f64` vertex per joint (the start of the segment leaving it), a bulge per
  arc, arcs of half a turn or more split into equal sub-arcs, and a joint
  closed within `max(tolerance.linear(), 8 eps * operand magnitudes)`. The
  extruder, the revolver and the half-space path build from that ring.

`axiolid-overlay` has exact arc predicates (ADR 0070) behind private
functions, but no self-contact test for one ring and no nesting check for an
arc polygon; `validate_arc_ring` names that gap in its documentation.

## Decision

Validation runs on the lowered ring, exactly, and reports what it finds.

- **One interpretation.** The validator lowers each contour with the same
  code the exact path uses (`lower_contour`, which `contour_to_arc_ring` now
  wraps) and judges that `ArcRing`. "Valid" means valid as the exact path
  will build it. Each ring edge maps back to the segment it came from, so
  findings name segments.
- **Tolerance-free contact decisions.** `axiolid_overlay::check_arc_polygon`
  decides self-contact, shared boundary and hole nesting from exact signs over
  the ring's `f64` data. No tolerance enters a contact or nesting decision.
  The tolerance enters only where it does today: joint closure, and
  `validate_arc_ring`'s per-ring degeneracies. A tangency that holds exactly
  is a contact; one that misses by `2^-40` is not.
- **The lowering is reported.** `ProfileReport::max_joint_gap` is the largest
  joint gap the lowering closed. The lowered boundary lies within
  `2 * max_joint_gap` plus rounding of the authored one (each lowered arc is
  a similarity image of the authored arc whose ends moved by at most the gap;
  a point of a sub-arc below half a turn moves by at most `|d0| + |d1 - d0|`).
  Parts of the authored boundary closer than that, without touching, are
  decided for the lowered ring.
- **Contacts are reports, not defects.** An isolated point where two rings
  meet is listed as a contact. `is_valid` accepts it (the OGC rule, and the
  straight overlay's `hole_outside` accepts a hole touching its outer ring at
  a vertex); `is_strict` refuses it. The OGC rule also needs a connected
  interior, so contacts that cut it apart (a cycle in the graph of rings and
  distinct contact points, found by union-find) are a `DisconnectedInterior`
  defect. Inside one ring, any shared point other
  than the joint two consecutive edges share is a defect: crossings, pinches,
  tangencies and fold-backs.
- **Refusals stay typed.** A contour that does not lower (family, frame,
  range, gap) yields its finding and nothing else is checked for it; nesting
  is not decided against a missing ring.
- **Placement.** The arc check is public in `axiolid-overlay`
  (`check_arc_polygon`), which owns the exact arc predicates. The validator,
  `ExactProfileValidator`, is public in `axiolid-construct`
  (`profile_validate`), which owns the lowering. The trait stays in
  `axiolid-profile`, since a representation crate may not depend on an
  algorithm. No new crate, no closure change.
- **Additive.** New items only; `#[non_exhaustive]` on the finding and
  report types.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Validate the authored `Circle2` and angle domain exactly | Joints are transcendental: it would refuse nearly every arc contour, and `FixedInterval` cannot prove a zero. |
| Chord the arcs and use the straight overlay (the stopgap) | Cannot certify a tangency, and judges a different shape than the one the exact path builds. |
| A new `axiolid-profile-check` crate | Contour, section, centre-line and derived lowering would all move out of construct: a topology change (ADR 0035, closure profiles) for a closure saving no current consumer needs. Revisit if a 2D-only consumer needs validation. |
| The validator in `axiolid-overlay` | It would add `axiolid-profile` and `axiolid-curve` to overlay's closure and duplicate `contour_lower`. |
| Contacts as defects | The OGC rule and the straight overlay both allow isolated contacts; a consumer that wants none uses `is_strict`. |
| Wire the validator into the exact builders | A scope decision of its own; a self-touching contour the builders accept is filed separately, not changed here. |

## Consequences

**Positive**

- Validation and construction cannot disagree about the geometry.
- Tangencies and near-misses are certified in either direction.
- Findings name the segments an author wrote.

**Negative / costs**

- Exact `crossings` on tangencies reaches big integers; cost is
  `O(n log n + k)` pairs through a box tree, not constant per edge.
- Closeness within the tolerance but without touching is not reported
  (#388), and polyline segments are refused like other families although
  they lower exactly (#389).

**Follow-ups / risks to watch**

- #388: report near-contacts closer than the linear tolerance.
- #389: lower polyline contour segments.

## Relation to existing code

- `crates/algorithms/planar/overlay/src/arc_check.rs`: `check_arc_polygon`,
  built on `exact_arc` (`crossings`, `same_point`, `BoxTree`) and
  `ArcArrangement`.
- `crates/algorithms/construction/construct/src/contour_lower.rs`:
  `lower_contour`, shared with `contour_to_arc_ring`.
- `crates/algorithms/construction/construct/src/profile_validate.rs`:
  `ExactProfileValidator`.
- `scripts/probe_arc_check_mutants.py`,
  `scripts/probe_profile_validate_mutants.py`: mutation probes of the two.
