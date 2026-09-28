# axiolid-brep-audit

Geometric consistency auditing for exact B-reps. It evaluates curves and surfaces and checks that edge vertices lie on their 3D curves and that every pcurve, mapped through its face's surface, lands on the curve of the edge it trims. It complements the exact, tolerance-free topological audit in `axiolid-topology`; because it compares positions, it needs a tolerance and reports agreement to within it (ADR 0052). It diagnoses only and repairs nothing.

```bash
cargo add axiolid-brep-audit
```

- API documentation: [docs.rs/axiolid-brep-audit](https://docs.rs/axiolid-brep-audit)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
