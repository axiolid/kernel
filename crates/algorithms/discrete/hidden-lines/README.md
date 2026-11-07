# axiolid-hidden-lines

The exact scalar provider for `MeshHiddenLines` (ADR 0090): the visible and
hidden pieces of the boundary, silhouette and crease edges of triangle
meshes (and, when selected, the edges between face groups, ADR 0106) in an orthographic view. Visibility is decided exactly in the view's
affine frame, with the interval filter first and dyadic big integers when it
cannot decide; no tolerance enters a decision, and the output endpoints are
the exact points rounded once. It refuses what it cannot decide with a typed
error that names the mesh. Candidate pairs are found in an `axiolid-spatial` BVH over view-space boxes
(worst case still edges times triangles); the examined pairs are bounded by
`max_candidate_pairs`. A
near clipping plane that cuts a `Solid` mesh leaves a cap that occludes, decided
by exact ray parity with a symbolic perturbation of the image point; a face of a
cut solid in the plane needs no special case. A solid wholly beyond the plane is
removed.
It does not draw perspective views, exact B-rep edges or smooth silhouettes.

`ScalarVisibleFaces` is the provider for `MeshVisibleFaces` (ADR 0094): the
visible face regions of caller-grouped triangles in the same view, with the
same preparation, clipping, caps and crop. The visible outlines, the clip
lines, the crop and the exact intersection curves of the triangles go into
the exact planar arrangement of `axiolid-exact-arrangement`; each side of
each piece is labelled by the exact frontmost source at a symbolic side
point, and the regions are traced from those labels. Region vertices are the
exact points rounded once where no edge passes through the rounding cell of
another vertex, and snap-rounded to the binary64 grid otherwise (ADR 0100).
Sources tied in depth each receive the area. The groups are the caller's per
triangle or, as `FaceGroups::Planar`, derived exactly: triangles join across
an edge of exactly two triangles with parallel, same-facing normals, and a
group is named by its smallest triangle (ADR 0106).
A region whose area collapses under rounding is dropped and counted in the
evidence, never refused.

```bash
cargo add axiolid-hidden-lines
```

- API documentation: [docs.rs/axiolid-hidden-lines](https://docs.rs/axiolid-hidden-lines)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-hidden-lines)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
