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
     "                    brep.shells().iter().map(|shell| (shell, false)).collect(),\n                    MeshClosure::Surface,",
     "                    brep.shells().iter().map(|shell| (shell, false)).collect(),\n                    MeshClosure::Solid,"),
    ("only the first shell tessellated", B,
     "                brep.shells().iter().map(|shell| (shell, false)).collect(),",
     "                brep.shells().iter().take(1).map(|shell| (shell, false)).collect(),"),
    ("solid-less brep still refused", B,
     "            if brep.shells().is_empty() {\n                if loose.faces.is_empty() {",
     "            if true {\n                if true {"),
    ("declared solid labelled surface", B,
     "            (shells, MeshClosure::Solid)",
     "            (shells, MeshClosure::Surface)"),
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
     "    let closure = if authored_mesh_is_closed(&mesh) {",
     "    let closure = if true {"),
    # The deviation variant reports closure the same way; the probed one is
    # `compile_mesh_reported`, which `tests/surface_models.rs` calls.
    ("reported outcome drops closure", C,
     "        Ok(CompileOutcome::tracked(built.mesh, built.fates.into_vec()).with_closure(built.closure))",
     "        Ok(CompileOutcome::tracked(built.mesh, built.fates.into_vec()))"),
    ("surface volume allowed", O,
     "            MeshClosure::Solid => Ok(&self.mesh),",
     "            MeshClosure::Solid | MeshClosure::Surface => Ok(&self.mesh),"),
    ("unknown volume allowed", O,
     "            MeshClosure::Solid => Ok(&self.mesh),",
     "            MeshClosure::Solid | MeshClosure::Unknown => Ok(&self.mesh),"),
]

probe_lib.main(MUTANTS, tests=T, timeout=900)
