# axiolid-mesh-compile

The scalar reference compilers for an `axiolid-model` geometry graph.
`ReferenceMeshCompiler` walks the graph and produces one `TriMesh` per root:
it resolves instances and collections, composes transforms, tessellates
profiles, sweeps, B-reps and authored meshes, and hands booleans to whichever
`MeshBoolean` provider it is given. `ReferenceExactCompiler` compiles the
families it supports to exact B-reps, places instances under rigid
transforms, unites, intersects and cuts placed extrusions, revolutions and
booleans with the general exact boolean, clips them by half-spaces
through the same boolean, and refuses the rest by name. `ReferenceMeshCompiler` also compiles curve
nodes (atomic 2D and 3D curves, their trims, composites and offsets,
open profiles, and their placements) to polylines within a certified chord
bound (`PolylineCompiler`). This crate
owns graph traversal and dispatch; the construction algorithms themselves
(profile flattening, extrusion, revolution, sweeps) belong to
`axiolid-construct`.

```bash
cargo add axiolid-mesh-compile
```

- API documentation: [docs.rs/axiolid-mesh-compile](https://docs.rs/axiolid-mesh-compile)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-mesh-compile)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

- The caller always supplies the tolerance through `ExecutionOptions`;
  there is no built-in default. Without an explicit chord budget, curves are
  flattened to the linear tolerance.
- Read volume through `CompileOutcome::solid_mesh`: a surface model can
  compile to a closed mesh without bounding a solid.
  A mesh is reported `MeshClosure::Solid` only when it is a closed,
  consistently wound two-manifold; a declared solid that tessellates open
  is `MeshClosure::OpenSolid`, and `solid_mesh` refuses it.
- Before a mesh boolean, operand coordinates a rounding residue from the
  other operand's axis-aligned faces (16 `f64::EPSILON` of the largest
  coordinate the operands reach, or their result reaches once an
  `Instance` places it, about `2e-8` at georeferenced coordinates) are set
  onto them exactly, so a tool a rounding error short of its host's face
  leaves no skin that thin, near the origin or placed far from it (#276,
  #312); a thicker skin is authored and kept whatever the tolerance
  (#291). A union reads contact at the tolerance instead: operands whose
  axis-aligned faces lie within it are set onto each other, adding at most
  the tolerance of material (#357, ADR 0096). A coordinate moves by at
  most the tolerance, none at zero, and the deviation report names the
  move (`deviation::SNAPPED_OPERANDS`). A boolean whose snapped operands are
  refused is cut as given, with no snap reported.
- A boolean of operands placed by the same linear part (a rotation and
  uniform scale), bit for bit, is cut in their shared frame and the result placed (#356,
  ADR 0096), so faces in contact under a general rotation are exactly
  coplanar to the snap, and a chain of unions stays in that frame. A
  refused shared-frame cut falls back to the world-frame cut. Operands
  placed by differing linear parts, or stretched by direction, are cut in
  the world frame
  (TODO(#391)).
- A union whose result holds two bodies within the linear tolerance of
  each other is refused with `GeomError::Degenerate`, starting with
  `UNMERGED_UNION_CONTACT`: the operands meet within the tolerance, as the
  exact compiler reads them (ADR 0080), and the tolerance-free mesh boolean
  left them apart (#356, #357; ADR 0096). Bodies further apart are
  authored and kept.
