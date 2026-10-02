# axiolid-brep

The strict exact B-rep result: an `axiolid-topology` graph bound to owned
catalogs of `Curve3` edge supports, `Curve2` pcurves and `Surface` face
supports, with every edge and pcurve interval stated explicitly and checked
before the value can exist. Faces and edges can carry persistent structural
names that survive rebuilds. `ExactBRep::transformed` places a B-rep under a
rigid motion, family for family, and refuses a scale or shear. It does not
evaluate, intersect, tessellate or traverse geometry.

```bash
cargo add axiolid-brep
```

- API documentation: [docs.rs/axiolid-brep](https://docs.rs/axiolid-brep)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-brep)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
