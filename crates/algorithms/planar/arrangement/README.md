# axiolid-arrangement

An editable planar subdivision: a doubly-connected edge list whose vertices, half-edges and faces keep stable handles across edits, so a caller can move a vertex or split a face without rebuilding the plane or losing track of which face is which. Orientation decisions use certified predicates, and the unbounded outer region is a real face. It is deliberately neutral: it exposes faces, boundaries, areas and adjacency, and leaves deciding that a face is a room to the caller. For one-shot polygon booleans with no retained structure, use `axiolid-overlay` instead.

```bash
cargo add axiolid-arrangement
```

- API documentation: [docs.rs/axiolid-arrangement](https://docs.rs/axiolid-arrangement)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-arrangement)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
