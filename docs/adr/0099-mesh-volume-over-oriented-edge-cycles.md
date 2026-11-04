# 0099 — Mesh volume over oriented edge cycles

- **Status:** Accepted
- **Date:** 2026-11-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #416 (consumer: a takeoff tool reading cached meshes). `volume_properties` requires a
closed two-manifold. Meshes cached by other tools are often closed in a
weaker sense: a T-junction where an opening was cut leaves a long edge matched
by two short ones, and a door or window of several closed shells shares edges
with its neighbours. In the consumer's model 158 of 327 measured components
are refused, although the divergence volume reproduces the stored volume in
every closed case, and the four genuinely open meshes are exactly those where
it does not.

## Decision

**D1. The audit.** `axiolid-mesh` gains `audit_mesh_cycle` and
`try_audit_mesh_cycle`, returning `MeshCycleHealth`: the unchanged
`MeshHealth` plus `unbalanced_edges` (undirected edges of usable triangles
whose net direction is not zero) and the smallest of them. A mesh is a closed
oriented cycle when it is surface-usable, has no degenerate triangle and no
unbalanced edge. The pass is the existing one (both edge sinks count, the
scratch bound is unchanged), and every closed two-manifold is a cycle.

**D2. Why that is sound.** Equal use in each direction means the triangles
form a 2-cycle, so the divergence sum does not depend on its base point and
equals the integral of the integer winding number `w`. The enclosed volume
needs `w` in {0, 1} as well, which no edge audit sees; that is the contract
`volume_properties` has always had (overlapping or duplicated shells count
twice, an inverted shell subtracts) and it is documented, not newly
introduced. It also answers the acceptance case of a box without top and
bottom: its vector area vanishes, so a translation-invariance test would
accept it, while eight rim edges are unbalanced.

**D3. The measure.** `axiolid-measure` gains `cycle_volume_properties` and
`CycleVolumeError` (`NotClosedOrientedCycle`, carrying the boxed audit, and
`ZeroSignedVolume`). It shares one private summation with `volume_properties`,
moved unchanged, so two-manifold input is bit-identical. Everything else is
refused as before. `volume_properties`, `second_moments` and
`MeshMeasureError` do not change.

**D4. The split.** T-junctions are repaired in `axiolid-heal`, which depends on
`axiolid-measure` and not the reverse, so the measure takes an already healed
mesh. `RepairAction::SplitTJunctions` and `split_t_junctions` split only
unbalanced edges, at vertices within `tolerance.linear()` of the edge line,
strictly inside it and further than the tolerance from both endpoints. No
vertex is created or moved. A triangle is declined and kept as it was when
its own third corner lies on its edge, when a piece would be degenerate or
face against the parent, or when the two ends of a split side hold different
values in a corner-indexed buffer, or when the per-vertex value of the vertex
it would reuse (a channel, or a normal) differs from the interpolation of the
side's endpoints at its parameter (compared exactly for `Linear` channels and
normals, bit-equal at both ends and the vertex for other blends); otherwise
new corners take the value of the corner they grew from. The deviation (at most the tolerance) is reported, not
refused. The action is opt-in and the recommended pipeline is weld, drop
degenerate elements, split T-junctions, then `cycle_volume_properties` with
the same tolerance.

**D5. Additive only.** `MeshHealth`, `MeshMeasureError` and `RepairReport` are
exhaustive public types and stay untouched; every change is an addition, so
each crate takes a patch release.

## Amendment

**D6. Certificate of w in {0, 1}** (#420). `axiolid-inspect` gains
`certify_unit_winding`, `WindingError` and, for the check, reuses
`VolumeInterval`. The prism decomposition of `overlap` needs only a closed
oriented 2-cycle, and over the faces as given (no inward flip, which would
change the sign of `w`) it yields both moments of the winding number:
the integral of `w` is the signed volume, the integral of `w^2` is the
unclamped pair sum of the mesh against itself. As `w` is an integer,
`w^2 - w >= 2` where `w` is not 0 or 1, so the volume of that set is at
most the winding excess `E = (int w^2 - int w) / 2 = int w(w-1)/2`. `E` is
an upper bound on that volume and equals it where `w` is in {-1..2} (a
triple-welded cube has `E` = 3 over a bad volume of 1); it is zero exactly
when `w` is 0 or 1 almost everywhere. The certificate refuses a mesh that is not a closed
oriented cycle (with the boxed audit) or not finite, takes a caller budget
`max_excess` (finite and positive, since outward rounding can never meet 0),
and returns the outward-rounded bracket of `E` (a `VolumeInterval` of `E`) when its upper
end is within the budget, `Exceeds` when its lower end is over, and
`Undecided` when the bracket straddles it. It accepts edge-touching shells
and cavities and refuses duplicated, overlapping and inverted shells. It
sits in `axiolid-inspect` because `axiolid-measure` cannot depend on it;
`overlap.rs` only extracts the shared loops, so its results are unchanged.

## Alternatives considered

| Option | Why not |
| --- | --- |
| A field on `MeshHealth` | It is exhaustive: a new field is a breaking change. |
| Let `volume_properties` accept cycles | It would change what an existing call refuses, and its callers rely on the two-manifold guarantee. |
| Split T-junctions inside the measure | A measure never edits the mesh, and heal depends on measure. |
| Insert the vertex into the other face | Creates geometry; the vertex is already in the mesh. |
| Skip the audit and test translation invariance | A box without top and bottom passes it. |

## Consequences

**Positive**

- Meshes closed up to T-junctions and edge-touching shells are measured, and
  open ones are still refused by name with the audit that says why.
- No existing result changes.

**Negative**

- `cycle_volume_properties` returns the winding-number integral, which equals
  the enclosed volume only for non-overlapping shells; callers who cannot
  exclude overlap use `certify_unit_winding` (D6).
- The ear clipping is a private copy of the one in `axiolid-mesh-compile`
  until it moves to `axiolid-mesh` (TODO(#422)).
