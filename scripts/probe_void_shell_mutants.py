"""Mutation probe for void shells and cavity results (#120).

Each mutant weakens one decision in how the mesh compiler tessellates a
solid's void shells or how the column builder returns a cavity, and must
turn a test red.

Not listed: the column builder's refusal of a cavity in a result of
several pieces. Two prism operands cannot produce that combination, so no
input reaches it; it guards future callers.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/execution/compile/src/brep.rs"
C = "crates/algorithms/construction/construct/src/column.rs"
TESS = ["-p", "axiolid-mesh-compile", "--test", "brep_tessellation"]
CAV = ["-p", "axiolid-construct", "--test", "column_cavities"]

MUTANTS = [
    ('voids not tessellated', B, '                check_void_closed(brep, shell)?;\n                shells.push(shell);', '                check_void_closed(brep, shell)?;', TESS),
    ('open void accepted', B, '    if shell.faces.is_empty() || uses.values().any(|&count| count != 2) {', '    if shell.faces.is_empty() {', TESS),
    ('outward void kept outward', B, '    if volume > 0.0 {\n        for triangle', '    if volume > 0.0 && false {\n        for triangle', TESS),
    ('void flipped into material', B, '    if volume > 0.0 {\n        for triangle', '    if volume < 0.0 {\n        for triangle', TESS),
    ('column drops the cavity', C, '        for void in voids.drain(..) {\n            cavities.push(emit.shell(&void)?);\n        }', '        voids.clear();', CAV),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=900,
    ).returncode

survivors = []
for name, rel, old, new, target in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(target)
        except subprocess.TimeoutExpired:
            code = -1
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
