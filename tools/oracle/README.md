# axiolid-oracle

An independent check for intersection, inversion and distance results
([ADR 0037](../../docs/adr/0037-mapped-3d-verification-oracle.md)). It maps a
claimed parameter-space result back into model space through `axiolid-evaluate`
and measures the 3D deviation there: contact between two curves or surfaces over
a claimed parameter box, and a search for a point closer than a claimed global
minimum distance. It is a falsifier, not a prover. A hit disproves a claim, and
no hit proves nothing. It shares no subdivision, interval or root-isolation code
with `axiolid-nurbs`, which it checks, and its dependency allowlist keeps it
that way. It is a dev-dependency of the crates it checks and is never published.
