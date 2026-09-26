"""Mutation probe for exact profile breadth (#111).

Each mutant weakens one decision in revolving sections with holes and
circles, unioning composite members exactly, carrying disjoint pieces as
separate solids, or tessellating every solid, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
BREP = "crates/representations/brep/src/lib.rs"
LOWER = "crates/algorithms/construction/construct/src/profile_lower.rs"
SECTION = "crates/algorithms/construction/construct/src/section_lower.rs"
ASSEMBLE = "crates/algorithms/construction/construct/src/assemble.rs"
TESS = "crates/execution/compile/src/brep.rs"
HOLES = ["-p", "axiolid-construct", "--test", "revolve_holes"]
COMPOSITE = ["-p", "axiolid-construct", "--test", "composite_breadth"]
DERIVED = ["-p", "axiolid-construct", "--test", "derived_composite"]
MESH = ["-p", "axiolid-mesh-compile", "--test", "brep_tessellation"]

MUTANTS = [
    ('appended void kept forward', BREP, '(Orientation::Forward, true) => Orientation::Reversed,', '(Orientation::Forward, true) => Orientation::Forward,', HOLES),
    ('hollow circle loses its bore', SECTION, '            vec![ring(circle.radius - thickness)]', '            Vec::new()', HOLES),
    ('member holes ignored in union', LOWER, '.any(|(outer, holes)| inside[*outer] && !holes.iter().any(|hole| inside[*hole]))', '.any(|(outer, _)| inside[*outer])', COMPOSITE),
    ('member vertices not welded', LOWER, '    weld_ring_vertices(&mut rings, tolerance);\n', '', COMPOSITE),
    ('zero-length pieces kept', LOWER, '        if ring.vertices[index].point == ring.vertices[next].point {', '        if false {', COMPOSITE),
    ('only the first piece kept', ASSEMBLE, '    if pieces.len() == 1 {', '    if !pieces.is_empty() {', DERIVED),
    ('only the first solid tessellated', TESS, '            for solid in brep.solids() {', '            for solid in brep.solids().iter().take(1) {', MESH),
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
