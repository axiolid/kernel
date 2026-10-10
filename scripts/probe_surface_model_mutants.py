"""Mutation probe for surface-model tessellation and closure (#161).

Each mutant weakens one decision: which shells a solid-less B-rep
tessellates, how closure is labelled and propagated, and the refusals built
on it. `tests/surface_models.rs` must fail on each. A build error also
counts as killed.

#265: "collection solid if any member is" re-anchored on the three-way
`combined_closure` that keeps `OpenSolid`.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/execution/compile/src/brep.rs"
CH = "crates/execution/compile/src/channels.rs"
TR = "crates/execution/compile/src/channels/transform.rs"
C = "crates/execution/compile/src/compiler.rs"
O = "crates/contracts/operations/compile/src/outcome.rs"
T = ["-p", "axiolid-mesh-compile", "--test", "surface_models"]

MUTANTS = [
    ("solid-less shells labelled solid", B,
     "(brep.shells().iter().collect(), MeshClosure::Surface)",
     "(brep.shells().iter().collect(), MeshClosure::Solid)"),
    ("only the first shell tessellated", B,
     "(brep.shells().iter().collect(), MeshClosure::Surface)",
     "(brep.shells().iter().take(1).collect(), MeshClosure::Surface)"),
    ("solid-less brep still refused", B,
     "        None if !brep.shells().is_empty() =>",
     "        None if false =>"),
    ("declared solid labelled surface", B,
     "            (vec![shell], MeshClosure::Solid)",
     "            (vec![shell], MeshClosure::Surface)"),
    ("collection solid if any member is", CH,
     "            (MeshClosure::Solid, MeshClosure::Solid) => MeshClosure::Solid,",
     "            (MeshClosure::Solid, _) | (_, MeshClosure::Solid) => MeshClosure::Solid,"),
    ("placement forgets closure", TR,
     "        closure: built.closure,",
     "        closure: axiolid_mesh_compile_contract::MeshClosure::Solid,"),
    ("surface subject accepted", C,
     '        refuse_surface_operand(subject, "subject")?;',
     "        let _ = subject;"),
    ("surface tool accepted", C,
     '                refuse_surface_operand(built, "tool")?;',
     "                let _ = built;"),
    ("authored mesh always solid", C,
     "    let closure = if axiolid_mesh::audit_mesh(&mesh, options.tolerance()).is_closed_two_manifold() {",
     "    let closure = if true {"),
    ("reported outcome drops closure", C,
     ".with_closure(built.closure)",
     ""),
    ("surface volume allowed", O,
     "            MeshClosure::Solid => Ok(&self.mesh),",
     "            MeshClosure::Solid | MeshClosure::Surface => Ok(&self.mesh),"),
    ("unknown volume allowed", O,
     "            MeshClosure::Solid => Ok(&self.mesh),",
     "            MeshClosure::Solid | MeshClosure::Unknown => Ok(&self.mesh),"),
]

probe_lib.main(MUTANTS, tests=T, timeout=900)
