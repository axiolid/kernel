"""Mutation probe for meshed sections (#193) and tangent voids (#194).

Each mutant drops the mesh path for sections or rounded rectangles, puts a
chord point back on every axis direction of a circle, or lets a pinched
boolean result through, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/construction/construct/src/profile.rs"
C = "crates/execution/compile/src/compiler.rs"
N = "crates/execution/compile/src/pinch.rs"
TESTS = ["-p", "axiolid-mesh-compile", "--test", "section_mesh", "--test", "tangent_void", "--lib"]

MUTANTS = [
    ('sections refused again', P, '        Profile::Section(section) => {\n            let contour = crate::section_lower::section_contour(section)?;\n            contour_rings(&contour, chord_error, tolerance)\n        }', '', TESTS),
    ('rounded rectangles meshed sharp', P, '        Profile::Rectangle(r) if r.outer_radius.is_some() || r.inner_radius.is_some() => {', '        Profile::Rectangle(r) if false => {', TESTS),
    ('circle chords start on the axes', P, '    let mut ring = flatten(core::f64::consts::PI / steps as Scalar)?;', '    let mut ring = flatten(0.0 * steps as Scalar)?;', TESTS),
    ('pinched results returned', C, '        if let Some(pinch) = crate::pinch::find(&outcome.mesh) {', '        if let Some(pinch) = None::<crate::pinch::Pinch> {', TESTS),
    ('an edge with three faces is fine', N, '.find(|(_, faces)| faces.len() > 2)', '.find(|(_, faces)| faces.len() > 4)', TESTS),
    ('separate fans are fine', N, '        if (1..faces.len()).any(|i| root(&mut parent, i) != first) {', '        if false && (1..faces.len()).any(|i| root(&mut parent, i) != first) {', TESTS),
    ('positions not welded', N, '            *ids.entry(key(*p)).or_insert_with(|| {', '            *ids.entry((p.x.to_bits(), p.y.to_bits(), at.len() as u64)).or_insert_with(|| {', TESTS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
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
