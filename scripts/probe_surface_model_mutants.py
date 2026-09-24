"""Mutation probe for surface-model tessellation and closure (#161).

Each mutant weakens one decision: which shells a solid-less B-rep
tessellates, how closure is labelled and propagated, and the refusals built
on it. `tests/surface_models.rs` must fail on each. A build error also
counts as killed.
"""
import pathlib, subprocess, sys

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
     "        .all(|part| part.closure == MeshClosure::Solid)",
     "        .any(|part| part.closure == MeshClosure::Solid)"),
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


def run(target):
    return subprocess.run(["cargo", "test", "-q", *target], cwd=ROOT,
                          capture_output=True, text=True, timeout=900).returncode


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found ({original.count(old)})"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(T)
        except subprocess.TimeoutExpired:
            code = -1
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)

print(f"{len(MUTANTS) - len(survivors)} / {len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
