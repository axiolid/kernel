# axiolid-inspect

Queries over triangle meshes: clearance between meshes, point containment and winding number, ray casting, line of sight, genus and per-component topology, plane detection, and intersection or difference volumes with a certified error bound. Containment reuses the exact ray-parity test of the mesh boolean, so the two cannot disagree. Every query reports a measurement or a typed refusal and leaves the verdict ("too close", "hidden") to the caller.

```bash
cargo add axiolid-inspect
```

- API documentation: [docs.rs/axiolid-inspect](https://docs.rs/axiolid-inspect)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
