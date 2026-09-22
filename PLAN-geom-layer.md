# Geometry capability layer — working plan

Branch `feat/geom-capability-layer`, worktree `/mnt/backup/wt/kernel-geom`.

## Goal

Close the three audited gaps against ifc-lite, and add the convex-collision
+ spatial-index layer that lets downstream apps express Solibri-style rules.

## Hard constraint: a concurrent session shares this repo

Another session is expanding `primitives2d`/`primitives3d` and landing on
`main` (HEAD `ae894d8` is theirs). Rules for this branch:

- **Never edit** `foundation/core/src/primitives2.rs`, `primitives3.rs`,
  `mat4.rs`, or anything else in `foundation/core/src/` unless unavoidable.
- New capability goes in **new crates / new files**.
- **Never `git add -A`.** Stage explicit paths only.
- Root `Cargo.toml` is the one shared file both sessions must touch
  (workspace members + version pins). Keep edits to appended lines so a
  merge is trivial.

## Workstreams

| # | Gap | Where | Status |
|---|-----|-------|--------|
| 1 | Constrained Delaunay + quality refinement | new crate `algorithms/planar/triangulate` | pending |
| 2 | Exact arithmetic in the production boolean | `providers/mesh/boolmesh/src/csg/` | pending |
| 3 | Persistent editable half-edge topology | new crate `representations/topology/dcel` | pending |
| 4 | Convex collision: SAT, GJK, EPA, MPR | new crate `algorithms/query/collide` | pending |
| 5 | Octree + k-d tree | new files in `algorithms/query/spatial` | pending |

## Validation strategy

- Every new crate carries its own tests; no capability lands untested.
- Gap 2 is the risky one: it changes a shipped numerical path. Required
  evidence before commit:
  - the existing differential corpus still passes,
  - a case that is WRONG before and RIGHT after (otherwise the change is
    unmotivated),
  - a benchmark delta, because exactness costs time and the cost must be
    stated rather than discovered later by a user.
- `scripts/gate.sh` must pass before any push.
- Mutation-test new gates: a check that cannot fail is not a check.

## Risks / rollback

- Gap 2 could regress boolean performance badly. If the cost is
  unacceptable, the fallback is a filtered cascade (fast path first, exact
  only on a straddling filter) rather than unconditional exact arithmetic.
- Drift on `main` from the concurrent session: rebase before push, re-run
  the gate on the rebased commit, never force-push shared history.

## Next concrete action

Gap 1: scaffold `algorithms/planar/triangulate` with a CDT that preserves
constraint edges, then add Ruppert/Chew refinement behind an explicit
quality target.
