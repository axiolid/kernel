# Domain interpretations stay downstream of the kernel

Status: accepted

## Context

A capability audit asked whether Axiolid should offer a `Footprint` type —
"the effect a component has when flattened to two dimensions".

The geometry already exists. `project_mesh(mesh, plane, tolerance)` folds a
triangle mesh onto a plane and unions the result, preserving holes: a mesh with
a through-hole projects to a region that still has the hole, not to an outline
or a convex hull. `Region` then supplies area, booleans, dilate/erode and
sweep. Adding `Footprint` would be naming, not implementing.

That naming is where the problem is. "Footprint" does not denote a fixed set of
points. Asked for the footprint of a building, a quantity surveyor, a fire
engineer, and a planning officer will give different answers about overhangs,
cantilevered balconies, and structure below grade — each correct for their
purpose. The geometry is identical in every case; only the inclusion rule and
the choice of reference plane differ.

A kernel function called `footprint` would have to pick one of those rules. It
would then be silently wrong for every consumer holding a different one, and
the name would imply there was no choice to make.

The kernel already draws this line elsewhere. `crates/algorithms/sampled/field/src/clearance.rs`
reports a free span and refuses to judge it: "0.9 m of free span" is geometry,
"too low" is policy. The same split applies here.

## Decision

We will keep the mechanical operation in the kernel and leave the domain
interpretation to consumers.

Axiolid provides `project_mesh` and `Region`: project, union, preserve holes,
report evidence. It does not provide a `Footprint` type, does not select the
reference plane, and does not decide which parts of a model to include.

A consumer that needs a footprint composes one — choosing the plane, applying
its own inclusion rule, and naming the result in its own vocabulary.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Add `Footprint` to the kernel with one inclusion rule | Picks a winner among incompatible domain definitions and hides the choice behind a name that implies there is none. Wrong for every consumer whose rule differs. |
| Add `Footprint` with a configuration enum for the rules | Requires enumerating every downstream discipline's convention up front, and a closed enum in a published kernel cannot be extended by the consumer who needs the next variant. |
| Rename `project_mesh` to `footprint` | Same objection, plus it narrows a general operation. The projection is equally the input to shadow analysis, clash silhouettes, and formwork outlines, none of which are footprints. |

## Consequences

**Positive**

- One operation serves footprints, shadow areas, clash silhouettes and
  formwork outlines, because it commits to none of them.
- Disagreements about inclusion rules surface in consumer code, where the
  domain expert can see and change them, instead of inside a kernel default.
- The kernel keeps its property that every output is checkable geometry rather
  than a judgement that depends on who is asking.

**Negative / costs**

- A consumer wanting a footprint writes the plane selection and inclusion rule
  themselves; the kernel does not shorten that step.
- Several consumers may write similar wrappers before a genuinely shared
  convention becomes visible. That duplication is the evidence needed to
  justify promoting one, and promoting it earlier would be guessing.

**Follow-ups / risks to watch**

- If multiple consumers converge on the same inclusion rule, the candidate for
  promotion is that rule — as a named, documented convention — not the word
  "footprint" attached to the existing projection.
- Watch for the same pressure on other domain nouns: gross floor area,
  envelope, silhouette. The answer is the same each time.

## Relation to existing code

- `crates/algorithms/planar/project/src/lib.rs` — `project_mesh` and
  `intersect_prism`; the module documents this boundary directly.
- `crates/algorithms/planar/overlay/src/region.rs` — `Region`, the 2D result
  type carrying area and boolean operations.
- `crates/algorithms/sampled/field/src/clearance.rs` — the existing precedent
  for reporting geometry and declining the verdict.
