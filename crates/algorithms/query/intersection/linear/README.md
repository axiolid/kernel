# axiolid-linear-intersection

Certified 2D intersections for lines and segments, with a minimal dependency closure so a line-query application does not pull in curves, surfaces, meshes or B-rep (ADR 0036). Results are classifications, not optional points: crossing, endpoint contact, parallel-disjoint, coincident, collinear-disjoint and overlap are distinct variants. Topology comes from certified predicates; the tolerance only governs acceptance of the computed coordinate. Invalid input is a typed refusal naming the operand.

```bash
cargo add axiolid-linear-intersection
```

- API documentation: [docs.rs/axiolid-linear-intersection](https://docs.rs/axiolid-linear-intersection)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
