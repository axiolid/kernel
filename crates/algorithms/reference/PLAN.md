# axiolid-reference plan

Design notes for the predicate suite. Status lives on GitHub, not here:
this file records *what the predicates must satisfy*, which does not
change when an item ships (kernel#25).

## The predicate suite
Error-free transformations (`two_sum`, `two_diff`, `two_product`) and
arbitrary-length expansion arithmetic are the shared foundation. On top
of them: `orient2d` and `orient3d` (is a point above/on/below a line or
plane), `incircle` and `insphere` (is a point inside/on/outside a
circumcircle or circumsphere).

Each is a filtered cascade: a fast floating-point path with a computed
error bound, escalating to exact expansion arithmetic only when the
bound cannot decide the sign. Static filters precompute bounds from a
coordinate magnitude limit, skipping the per-call permanent computation.

## Gates
- Differential vs an independent exact oracle (i128 rational, integer inputs
  bounded so the determinant cannot overflow).
- Measured escalation rate per degeneracy tier, asserted to stay in band.
  The degeneracy benchmark harness reports throughput AND escalation rate
  at 0%, 0.01%, 1%, 10% degenerate inputs.
- Mutation probes on every filter bound and every exact path.

## Relationship to adopted predicates
`boolmesh` carries its own predicates and is MPL-2.0, so replacing them means
forking. See ADR 0016: ours serve our own algorithms and act as an independent
audit oracle for adopted ones, rather than trying to displace them.
