# axiolid-mesh-compile

The scalar reference compilers for an `axiolid-model` geometry graph.
`ReferenceMeshCompiler` walks the graph and produces one `TriMesh` per root:
it resolves instances and collections, composes transforms, tessellates
profiles, sweeps, B-reps and authored meshes, and hands booleans to whichever
`MeshBoolean` provider it is given. `ReferenceExactCompiler` compiles the
families it supports to exact B-reps, places instances under rigid
transforms, cuts placed openings from placed extrusions with the general
exact boolean, clips them by half-spaces through the same boolean, and
refuses the rest by name. This crate
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
- Before a mesh boolean, operand vertices within the linear tolerance of
  the other operand's faces are moved onto them, so a tool a rounding
  error short of its host's face leaves no skin that thin (#276). The move
  is at most the tolerance, none at zero, and the deviation report names
  it (`deviation::SNAPPED_OPERANDS`).
