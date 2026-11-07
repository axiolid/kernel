# axiolid-exact-arrangement

An exact arrangement of straight segments in the plane. Points are rational (homogeneous coordinates over exact dyadic numbers), so every orientation, order and containment decision is an exact sign: segments are cut at every crossing, T-junction and collinear overlap into pieces that no other piece crosses, equal points are joined, and rings of caller-labelled half-edges are traced and nested into outer boundaries with holes. No tolerance enters. Only the points a caller finally returns are rounded, once each, or, where an edge passes through the rounding cell of another vertex, snap-rounded to the binary64 grid (`snap_round`, ADR 0100): each edge is routed through the centre of every cell it meets, decided exactly, and the result carries a displacement bound.

It is the shared source of truth for those signs: the plane section of `axiolid-reference` (in-plane faces, ADR 0093) and the visible faces of `axiolid-hidden-lines` (ADR 0094) both decide on it. It knows no mesh, no plane and no face; callers label the pieces. For an editable subdivision with persistent handles, use `axiolid-arrangement`; for polygon booleans with a tolerance, use `axiolid-overlay`.

```bash
cargo add axiolid-exact-arrangement
```

- API documentation: [docs.rs/axiolid-exact-arrangement](https://docs.rs/axiolid-exact-arrangement)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-exact-arrangement)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
