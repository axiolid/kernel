# axiolid-hidden-lines

The exact scalar provider for `MeshHiddenLines` (ADR 0090): the visible and
hidden pieces of the boundary, silhouette and crease edges of triangle
meshes in an orthographic view. Visibility is decided exactly in the view's
affine frame, with the interval filter first and dyadic big integers when it
cannot decide; no tolerance enters a decision, and the output endpoints are
the exact points rounded once. It refuses what it cannot decide with a typed
error that names the mesh. Candidate pairs are found by a bounding-box scan that is quadratic in edges
times triangles; the examined pairs are bounded by `max_candidate_pairs`. A
near clipping plane that cuts a `Solid` mesh leaves a cap that occludes, decided
by exact ray parity with a symbolic perturbation of the image point; a face of a
cut solid in the plane needs no special case. A solid wholly beyond the plane is
removed.
It does not draw perspective views, exact B-rep edges or smooth silhouettes.

```bash
cargo add axiolid-hidden-lines
```

- API documentation: [docs.rs/axiolid-hidden-lines](https://docs.rs/axiolid-hidden-lines)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-hidden-lines)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
