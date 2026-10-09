# 0084 — Bounded half-space boundaries may be profiles

- **Status:** Accepted
- **Date:** 2026-10-09
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #277 (consumer: openbimrs/ifc#398, axioval/engine#307).
`SolidOperation::BoundedHalfSpace` took its boundary as a closed 2D curve,
and both compilers read only a `Curve2::Polyline`. IFC allows
`IfcPolygonalBoundedHalfSpace` boundaries of line and circular-arc
segments (a composite of polylines and trimmed lines and circles, or an
indexed poly-curve with arcs), and a real wall is clipped by one of six
segments, two of them arcs over one circle of radius 1.2. The reader
refuses them rather than polygonise, since tessellation is the kernel's.

`Curve2` has no composite variant. The graph's `CurveRelation::Composite`
can chain trimmed curves, but neither compiler resolves a closed relation
into a region, and its closure, sense and trimming would need a second
lowering beside the one profiles already have. The profile `Contour` is
exactly a closed chain of bounded `Curve2` segments with a sense each,
and the exact extruder already turns its lines into plane walls and its
arcs into cylinder walls (ADR 0053), the mesh flattener chords it under a
certified bound.

## Decision

The boundary of a bounded half-space may be a `Profile` node as well as
a 2D curve. A boundary with arcs is a `Profile::Contour` of `Line2` and
`Circle2` segments.

- **Model.** The field stays a `NodeId`; graph validation accepts a 2D
  curve or a profile and refuses anything else as "curve2 or profile".
  No type, variant or field is added, so the change is additive.
- **Meaning.** The half-space is its side of the plane within the
  infinite prism of the profile's region, swept along the plane normal,
  the profile framed by the operation's placement as the polyline is. An
  arc bounds it by a right circular cylinder; the sweep is never oblique.
- **Both compilers check it alike** (`mesh-compile`'s
  `half_space_boundary`): each contour must close within tolerance, every
  arc have a positive radius, every segment a length, and no two edges
  meet within the tolerance other than neighbours at their joint. A
  boundary that fails is refused with `InvalidInput` naming the fault.
- **Exact.** The prism is `extrude_profile_exact` of the profile, mirrored
  in its `x` axis when the kept side is opposite the normal (as the
  polyline is). **Amended (#288):** the mirror was first a derived
  profile, because the general boolean refused some cuts by a reflected
  cylinder wall; with that fixed, the prism of the profile as authored is
  reflected by its placement instead.
- **Mesh.** The profile is flattened with the chord budget; a clipped
  solid's deviation is certified against the exact result, as for every
  boolean the exact compiler builds (#235).

## Alternatives considered

| Option | Why not |
| --- | --- |
| A `Curve2::Composite` variant | A second closed-chain type beside `Contour`, with its own closure, sense and lowering, and a new variant in a widely matched enum, for a shape profiles already state. |
| Resolve a `CurveRelation::Composite` of trimmed curves | Neither compiler reads a relation as a region; the trim and sense resolution would duplicate the profile lowering. Readers that build composites can build contours. |
| Polygonise arcs in the reader | Tessellation belongs to the kernel, and a polygonal boundary has no exact cylinder wall to certify against. |

## Consequences

**Positive**

- One exact boundary type for profiles and half-spaces; every profile
  family the exact extruder builds is a boundary for free.
- The clip of a wall by an arc boundary is exact, with closed-form volume,
  and its mesh is certified.

**Negative / costs**

- A boundary is now either a curve or a profile, so a consumer matching on
  the boundary node must handle both.
- The self-crossing check decides within the tolerance, not with exact
  predicates: edges closer than the tolerance are read as touching.

**Follow-ups / risks to watch**

- The general boolean's refusal of some cuts by a reflected cylinder wall
  (a downward extrusion of an arc profile) was worked around here at
  first; #288 fixed it in the boolean (a boundary cut placed on a line
  pcurve a turn of its parameter away) and removed the workaround.

## Relation to existing code

- `crates/representations/modeling/graph/src/solid_operation.rs`,
  `validation.rs` (`ExpectedReference::HalfSpaceBoundary`).
- `crates/execution/compile/src/half_space_boundary.rs`,
  `compiler.rs` (`boundary_rings`), `exact/clip.rs` (`Footprint`).
- `crates/execution/compile/tests/bounded_half_space_arcs.rs`,
  `scripts/probe_arc_bounded_half_space_mutants.py`.
