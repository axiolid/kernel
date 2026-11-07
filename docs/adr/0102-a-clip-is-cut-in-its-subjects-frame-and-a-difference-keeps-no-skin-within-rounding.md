# 0102 — A clip is cut in its subject's frame, and a difference keeps no skin within rounding

- **Status:** Accepted
- **Date:** 2026-11-06
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0096](./0096-a-union-decides-operand-contact-at-the-linear-tolerance.md), [0098](./0098-a-union-the-mesh-cut-cannot-settle-is-decided-by-the-exact-compiler.md)

## Context

Issue #400, found by a consumer (openbimrs/ifc#655). A wall turned by a
placement and cut by a bounded half-space clip that reaches through its whole
thickness came back as a closed, consistently wound mesh labelled `Solid`
with a lid across the slot, a skin about `1e-14` thick. Lowered in its own
frame the same cut opens the slot cleanly.

The snap ([`snap::reach`], about 16 ε·S, with the placed scale of #312)
closes a skin that thin only where the faces are axis-aligned (#291).
The clip's faces, placed by a turn, lie on planes no `f64` coordinate holds
exactly, so they land a residue off the wall's faces and the mesh boolean,
which has no tolerance (ADR 0047), keeps the lid. ADR 0096 D3 already cuts
operands that share a linear part in their common local frame, but it
leaves out a half-space clip: the clip is built from the subject, so it
never was an operand of the shared frame. ADR 0096 D7 left a related
defect "recorded, not decided".

A skin thinner than the rounding reach is residue, not authored geometry. A
skin above the reach but below the tolerance stays authored (D7, #291).

## Decision

**D1. A clip is cut in its subject's frame.** When the subject is framed (an
`Instance`, or a boolean with `Built::local`) by a similarity, the clip is
re-expressed in that frame, built against the subject's local mesh, cut there
and the result placed by the frame. This extends ADR 0096 D3 to half-space
tools. Where the linear parts are bit-equal the relative placement is the
exact translation of `relative_offset`, identity when the translations are
equal: bitwise the cut of the whole in its own frame. Otherwise it is
`frame.inverse() * chain`, rounded; the clip's faces then lie a residue off
the subject's axis-aligned local faces, and the existing snap lands them.
Unlike D3, bit-equality is not required: a clip is a constructed stand-in with
no coordinates of its own to preserve. Any error falls back to the world cut,
so it is never worse.

**D2. A difference or intersection returns no skin within rounding.**
After the pinch check, the cut reads its result with
`axiolid_inspect::thin_sheet` at the cut's rounding reach: two triangles that
face each other on one plane, closer than the reach, over more than an edge.
A result that has one is refused by name,
`SKIN_WITHIN_ROUNDING`, as `Degenerate`, naming the two triangles, their
offset, the reach and where they lie. The guard reads the rounding reach,
never the tolerance, so a skin between the two stays authored (ADR 0096 D7). A
union is not guarded here: ADR 0098 D2 returns a union the referee refuses,
and refusing it is a decision of its own (#434).

**D3. The exact referee decides the guard's refusal, once, at the root of
its chain of clips.** A cut whose result holds a skin within rounding is not
refused where it is cut. If the node is only the subject of further
differences and intersections, and no root asked for it, the cut result is
carried on with the skin pending, so the next cut works on it as it always
did; a cut that leaves no skin clears it. The decision rises to the root of
the maximal chain of clipping booleans (a union, an instance, a tool use or a
requested root ends a chain): the exact referee (ADR 0098, #423) compiles
that node once, the whole chain, bounded half-space tools included. The
referee's own extra-shell trigger of #423 leaves a half-space to the mesh
result, but a skin within rounding is decided for it too. A pinch against
the mesh a pending skin left goes to the referee at that node.

The referee's result is accepted only if what the exact compiler read as one
is within the rounding reach (`within_reach(report.linear(), reach)`, the
gate #423 uses, so an authored sub-tolerance skin elsewhere is never fused)
and its tessellation holds no thin sheet. It is then returned as
`Decided::Referee` under the operator's `REFEREED_*` detail, with
`COINCIDENT_OPERANDS` certified to at most the reach. Otherwise the boolean
is refused with `"{text}; the exact referee refused: {why}"`, the text
starting with `SKIN_WITHIN_ROUNDING`. Refereeing each cut instead costs an
exact compile of the whole subtree per cut, and the tessellation it leaves
as the next subject pinches the next opening.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Snap onto general planes | A vertex cannot land exactly on a non-axis-aligned plane in `f64`, and #291 says a move that does not land exactly is harmful. |
| `axiolid_heal::self_intersections` as the guard | It misses a parallel skin that does not cross, and costs more. |
| Referee on near contact | It reverses D7, the reason #423 narrowed its own trigger. |
| Delete the facing triangles | Substitute geometry. |
| Referee every cut of a chain | One exact compile of the subtree per cut (sixteen flush openings took 2 s against 0.6 s), and each tessellated subject pinches against the next opening. |
| Guard with the tolerance | A skin of half a millimetre is authored, and a tolerance-thick guard would refuse it. |

## Consequences

**Positive**

- The reported wall opens its slot, bitwise as the same cut lowered in its own
  frame.
- A difference or intersection never returns a `Solid` whose faces face each
  other within the reach.
- `thin_sheet` is a public query on any mesh.

**Negative / costs**

- A chain of differences whose mesh cuts keep a lid within rounding costs one
  exact compile, at the chain's root. A sweep of a turned wall and clip in world coordinates, with no
  `Instance`, goes to the referee at 19 of 40 angles and every one is
  decided; at a site 600 km out the exact compiler cannot classify a point
  near a face boundary and the skin is refused by name.
- One BVH pass over the result of every difference and intersection cut.
- A chain whose residue skin lies beside a feature above the reach but below
  the tolerance is refused by name: the round-tripped flush slab, whose
  openings the round trip also moved by `9.3e-10`. The exact compiler reads
  that feature as one at the tolerance and cannot classify below it (#457).
  mesh-compile 0.3.44 returned that slab with every opening lidded, and
  placed it self-intersecting.

**Follow-ups / risks to watch**

- Whether the guard applies to unions: #434.
- The exact compiler classifying below a feature it would read as one: #457.

## Tests

`tests/turned_clip.rs` (the reported wall placed, written in world
coordinates and at a site: closed, one body, measured, with no thin sheet, no
self-intersection and no lid; an authored half-millimetre skin stays; a tool
offset in the frame moves the slot; the cut equals the whole cut in its own
frame bitwise; a sweep of 1 to 40 degrees in world coordinates is a clean
solid or refused by `SKIN_WITHIN_ROUNDING`; a pinned angle the referee
decides, an authored skin it keeps, a pinned site where the exact compiler
refuses, and a pinned site where it reads above the reach), the
`frame::relative` unit
tests, `tests/flush_through_cut.rs` (the round-tripped slab refused by name;
the same slab placed, and with exactly turned openings, opens every opening),
`axiolid-inspect`'s `tests/sheet.rs`, and
`scripts/probe_clip_skin_mutants.py`.

## Relation to existing code

- `crates/execution/compile/src/compiler.rs` (`clip_in_subject_frame`, `cut`,
  `CutRefusal::Skin`), `src/frame.rs` (`similarity`, `relative`),
  `src/skin.rs`.
- `src/referee.rs` (`boolean`, its `rounding` parameter).
- `crates/algorithms/query/inspect/src/sheet.rs`.

## Amendment 2026-11-07: regions in the subject's frame (#326)

D1 covers a region, an intersection of half-spaces used as a clip's tool
(ADR 0080, #326): each leaf is placed relative to the subject's frame by
`frame::relative`, the stand-in is built there leaf by leaf and cut with the
subject, and a leaf with no finite form in that frame, or a region that
leaves nothing of the subject, falls back to the world cut as before. A
region of one leaf is the single clip of D1, unchanged. Without it a region
of two half-spaces on a turned wall is cut in the world frame, and its mesh
is not the whole's cut in its own frame, placed (it is, bit for bit, with
it: `tests/half_space_region.rs`, five angles from 0.3 to 2.5 rad).
