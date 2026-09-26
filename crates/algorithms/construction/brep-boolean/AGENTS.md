# axiolid-brep-boolean

General exact booleans between exact B-reps with analytic faces: the
general-fuse pipeline of ADR 0075, built in stages.

## Pipeline and module ownership

- `section.rs`: section edges -- every face pair's exact intersection curve
  (`axiolid_nurbs::exact_surface_intersection`), cut where it crosses a
  boundary edge (transversally: against the ADJACENT face's surface, or
  across a seam against the plane through the ruling), membership from
  `axiolid_measure::FaceDomain`.
- `split.rs`: one face cut along its section edges into regions, in the
  face's parameters, with exact pcurves (planes and cylinders in stage 1).
- Later slices add classification, selection and sewing, each in its own
  module.

## Rules

- Never approximate: a pair or configuration a stage cannot build exactly is
  refused by name (`BooleanError`), never meshed or fitted.
- Decisions come from exact predicates or certified membership; a point too
  close to a boundary to decide is refused, not guessed.
- Every result test checks the exact volume identity with
  `axiolid_measure::exact_properties` once solids are assembled.

## Verification

- `tests/section_edges.rs`: section edges lie on both boundaries, close into
  loops, and match closed-form curves.
- `tests/split_faces.rs`: regions cover each face exactly and match
  closed-form areas.
- Never cut a section with an exact curve/curve test against a boundary
  edge: two curves on one surface meet only up to the rounding of their
  doubles. Cut against a transverse surface.
