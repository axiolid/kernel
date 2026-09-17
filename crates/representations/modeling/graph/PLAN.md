# axiolid-model plan

Design notes for the immutable geometry DAG.
Status lives on GitHub, not here (kernel#25).

## Standing invariants

Node handles carry a graph-owner brand. Insertion rejects foreign,
forward, and semantically invalid reference families before an immutable
graph can exist -- the brand is what makes a handle from one graph
unusable in another, so an invalid graph is unrepresentable rather than
merely undetected.

## Shape of the work

Graph visitors, budgets, provenance side tables, and complete compiler
coverage.

## Exit evidence

Targeted tests, feature-isolated compile where applicable, mutation-verified
architecture/validation gates, and benchmarks before performance claims.
