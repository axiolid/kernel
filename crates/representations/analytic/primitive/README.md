# axiolid-primitive

Exact parametric primitive solids and half-spaces, used as CSG leaves. The
values stay exact until someone explicitly tessellates them: constructors do
no tessellation or boolean work, and the finite margin used when a
half-space has to be clipped for meshing is an explicit parameter. It has no
mesh or kernel dependency.

```bash
cargo add axiolid-primitive
```

- API documentation: [docs.rs/axiolid-primitive](https://docs.rs/axiolid-primitive)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
