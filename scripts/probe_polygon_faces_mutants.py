"""Mutation probe for authored polygon triangulation (#160).

Each mutant weakens one decision in `compile_authored_polygons` or
`planar::triangulate_polygon`; `tests/authored_polygons.rs` must fail on it.
A build error also counts as killed.

First run: 9 / 11. Two survivors, neither a test gap:
- "repeated closing corner kept": earcut drops coincident consecutive
  points itself, so the compiler's own dedup was dead code. Removed; the
  repeated-corner test still pins the behaviour end to end.
- "triangle winding not normalised": earcut 0.4 always emits triangles
  counter-clockwise in the projected frame, so the normalising swap never
  fires. Kept as a guard against an earcut convention change; equivalent
  today, so not listed. "triangle winding always flipped" covers the swap.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/execution/compile/src/planar.rs"
C = "crates/execution/compile/src/compiler.rs"
T = ["-p", "axiolid-mesh-compile", "--test", "authored_polygons"]

MUTANTS = [
    ("planarity check skipped", P,
     "    if worst > linear {\n        return Err(PolygonRefusal::NotPlanar(worst));",
     "    if worst > linear && false {\n        return Err(PolygonRefusal::NotPlanar(worst));"),
    ("planarity against zero, not tolerance", P,
     "    if worst > linear {", "    if worst > 0.0 {"),
    ("triangle winding always flipped", P,
     "        if area < 0.0 {\n            triangle.swap(1, 2);",
     "        if area >= 0.0 {\n            triangle.swap(1, 2);"),
    ("area cross-check skipped", P,
     "    if expected <= slack || (covered - expected).abs() > slack {",
     "    if expected <= slack {"),
    ("holes added instead of subtracted", P,
     "        expected += if index == 0 { area } else { -area };",
     "        expected += area;"),
    ("holes not passed to earcut", P,
     "    let mut indices = earcut_projected(&flat, &hole_starts);",
     "    let mut indices = earcut_projected(&flat[..hole_starts.first().copied().unwrap_or(flat.len())], &[]);"),
    ("triangles pass through untriangulated", C,
     "            if face.outer.len() == 3 && face.holes.is_empty() {",
     "            if face.holes.is_empty() {"),
    ("hole indices not range-checked", C,
     "            .flat_map(|face| face_rings(face).flatten().copied())",
     "            .flat_map(|face| face.outer.iter().copied())"),
    ("face index in error off by one", P,
     "        _ => GeomError::InvalidInput(format!(\n            \"authored polygon face {face} cannot be triangulated: {refusal}\"\n        )),",
     "        _ => GeomError::InvalidInput(format!(\n            \"authored polygon face {} cannot be triangulated: {refusal}\", face + 1\n        )),"),
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
