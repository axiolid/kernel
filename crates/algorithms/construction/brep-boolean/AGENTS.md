# axiolid-brep-boolean

General exact booleans between exact B-reps with analytic faces: the
general-fuse pipeline of ADR 0075, built in stages.

## Pipeline and module ownership

- `section.rs`: section edges -- every face pair's exact intersection curve
  (`axiolid_nurbs::exact_surface_intersection`), trimmed exactly to where it
  lies inside both faces (curve/boundary-edge crossings from
  `exact_curve_curve_intersection3`, membership from
  `axiolid_measure::FaceDomain`).
- Later stages add face splitting, classification, selection and sewing,
  each in its own module.

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
