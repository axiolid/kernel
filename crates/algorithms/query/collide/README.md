# axiolid-collide

Convex collision queries by the separating axis theorem: whether two convex shapes overlap and, when they do not, how far apart they are and along which axis. It deliberately does not report penetration depth; callers who need EPA-style contact for physics want a physics engine. For clearance between triangle meshes, use `axiolid-inspect`.

```bash
cargo add axiolid-collide
```

- API documentation: [docs.rs/axiolid-collide](https://docs.rs/axiolid-collide)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
