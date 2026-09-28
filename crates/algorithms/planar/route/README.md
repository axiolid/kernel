# axiolid-route

Exact planar shortest paths over a visibility graph, plus distance maps, farthest points and forced walks built on the same graph. Which edges exist is decided with certified `orient2d`, so the combinatorics are exact; path lengths are sums of square roots in `f64` and carry ordinary rounding. Oversized input is refused with a proven lower bound rather than truncated, and the budget is a caller parameter. It reports routes and typed unreachable reasons, never whether a route is acceptable. For grid-sampled routing over layered fields, see `axiolid-field-ops`.

```bash
cargo add axiolid-route
```

- API documentation: [docs.rs/axiolid-route](https://docs.rs/axiolid-route)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
