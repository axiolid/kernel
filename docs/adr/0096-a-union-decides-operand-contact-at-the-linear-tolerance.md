# 0096 — A union decides operand contact at the linear tolerance

- **Status:** Accepted
- **Date:** 2026-11-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0080](./0080-exact-compiler-cuts-placed-openings-with-the-general-boolean.md)

## Context

Issues #356 and #357 (consumer: rheonda/cad). Both are unions of two placed
solids that are in contact, and both come back from the mesh compiler as two
disjoint bodies labelled `MeshClosure::Solid`.

- **#356, face contact under a general rotation.** A column and a corbel
  stand on the same plane in exact arithmetic. Placing rounds each vertex on
  its own, about `1e-16`. The residue is far inside the snap's reach, but the
  snap (ADR 0080, #291) lands only on exactly axis-aligned planes, and after
  a turn about Z only the z faces are. The mesh boolean then decides on the
  raw values: two shells, a pinch (refused, #194), or one body with ulp-thin
  wedge cracks where the faces cross.
- **#357, a gap of `0.5e-6` under a tolerance of `1e-6`.** The snap reaches
  the rounding scale only (#291: "anything above the reach is authored and
  kept, whatever the tolerance"). The mesh boolean has no tolerance (ADR
  0047) and returns two boxes. The exact compiler reads the same faces as
  `CoincidentSupports` (plane offset at most the tolerance) and fuses them, so
  the two compilers disagree about a union the caller declared in contact.
  The deviation report does not catch it: it measures from the exact surface
  to the mesh, and the inner faces of a two-shell mesh are invisible to a
  one-sided measure.

The fix belongs in `axiolid-mesh-compile`, which owns the tolerance, the
frames, the snap and the refusals. Not in the boolean provider (its contract is
tolerance-free), and not in the exact path (it already reads contact within
the tolerance; it is checked, not changed).

## Decision

**D1. Policy.** In a union, operands within the linear tolerance are in
contact: they are merged, or the union is refused by name. This changes the
#291 rule for unions only; differences and intersections keep the rounding
reach. Closing a gap in a union adds at most the tolerance of material and
never removes any, the caller declared those points equal through the
tolerance, and the exact compiler already fuses them. #291's harm was
deleting authored skins in differences and landing inexactly; neither applies.

**D2. The guard.** The result of a mesh-compiled union never holds two bodies
(index components) whose surfaces lie within the tolerance of each other.
Otherwise the union is refused with `GeomError::Degenerate`, whose message
starts with the new `UNMERGED_UNION_CONTACT`. The check runs on the result,
after the pinch check, inside `cut`, so it applies to snapped operands and to
operands as given, and #291's never-worse fallback still holds. It costs
nothing for a result of one component. It also catches a union that leaves a
cavity open to the outside within the tolerance, the tolerance generalisation
of #194's pinch. `Degenerate` carries the distance and the bodies' boxes, as
the pinch refusal does; the prefix constant keeps it matchable.

**D5. Dependency.** `axiolid-inspect` gains `surface_gap`, the distance between
surfaces with no containment test. `axiolid-mesh-compile` moves
`axiolid-inspect` from dev- to production dependencies. It was already
allowlisted; its closure gains `axiolid-inspect` and `axiolid-heal`, and
everything else they need is already in it. Every declared closure profile
that names `axiolid-mesh-compile` lists it under `forbidden_internal`, and none
forbids `axiolid-inspect` or `axiolid-heal`, so no profile changes
(`cargo xtask architecture closure check` passes). The same change fixes `min_gap` (#390), which read a bar through
a box as clearance because no edge-crosses-face test existed.

**D3, D4, D6** (booleans in the frame their operands share; the union's contact
reach; operands with different linear parts) are decided, landing in the next
change; they are completed here then.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Rounding only for unions too, treating #357's two shells as correct | The two compilers keep disagreeing, and a report certifies a topology the exact result does not have. |
| Refuse the whole band between the rounding and the tolerance | Fail-closed, but it refuses unions whose operands the caller declared equal. |
| Guard as "more components than the operand contact graph allows" | Needs point-in-mesh containment between the operands' components, and misreads a body standing in another's cavity. Distance between result components needs neither. |
| A private triangle-triangle distance in `mesh-compile` | Mesh clearance is `axiolid-inspect`'s capability; it is added there. |
| `UnsupportedInput` instead of `Degenerate` | Its message is a fixed `&'static str` and cannot carry the distance and boxes. |

## Consequences

**Positive**

- A union never returns bodies within the tolerance of each other.
- No change to the `MeshBoolean` contract or to boolmesh.

**Negative / costs**

- A union whose operand already holds two shells within the tolerance is now
  refused even when the other operand is far away. Such an input is itself
  non-manifold under the tolerance; the refusal is accepted.
- One body with sub-rounding wedge cracks, where contact faces met under
  different placements, is not seen by the guard (`TODO(#391)`).

**Follow-ups / risks to watch**

- The deviation report is one-sided and cannot see extra shells in a measured
  boolean (`TODO(#392)`).
- Open question: a difference whose tool is thinner than the tolerance cutting
  through can leave two pieces closer than the tolerance, labelled `Solid`.
  Under D1 they would be in contact, but under #291 they are authored. This
  stays as it is; extending D2 to differences and intersections reverses #291
  and is a separate decision.

## Relation to existing code

- `crates/execution/compile/src/contact.rs`, `src/compiler.rs` (`cut`).
- `crates/algorithms/query/inspect/src/clearance.rs` (`surface_gap`).
- Tests: `crates/execution/compile/tests/union_contact.rs`,
  `tests/rotated_face_contact.rs`.
