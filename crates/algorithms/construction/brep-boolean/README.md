# axiolid-brep-boolean

Exact union, intersection and difference of two exact B-rep solids whose
faces lie on planes, cylinders, elliptical cylinders, cones, spheres, tori
and B-spline surfaces: the general-fuse pipeline of ADR 0075 (section edges,
face splitting in each face's own parameters, certified classification,
sewing). Operands that touch rather than cross are handled. Nothing is
meshed or fitted: a configuration the pipeline cannot build exactly is
refused with a typed error. It does not tessellate its result and does not
work on meshes; mesh booleans are operation providers selected through the
execution layer.

```bash
cargo add axiolid-brep-boolean
```

- API documentation: [docs.rs/axiolid-brep-boolean](https://docs.rs/axiolid-brep-boolean)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

`axiolid-construct` keeps narrower exact booleans that predate this crate:
planar polyhedra (`polyhedron::boolean_polyhedra_exact`) and coaxial
prisms and plane cuts over one arc arrangement (`boolean_exact`). They are
independent exact paths, and this crate's tests check vertical-column
results against them.
