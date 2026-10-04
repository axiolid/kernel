"""Mutation probe for authored polygon triangulation (#160, #254).

Each mutant weakens one decision in `compile_authored_polygons` or
`planar::triangulate_polygon` / `triangulate_authored_polygon`;
`tests/authored_polygons.rs` or the `planar` unit tests must fail on it.
A build error also counts as killed.

First run (#160): 9 / 11. Two survivors, neither a test gap:
- "repeated closing corner kept": earcut drops coincident consecutive
  points itself, so the compiler's own dedup was dead code. Removed; the
  repeated-corner test still pins the behaviour end to end.
- "triangle winding not normalised": earcut 0.4 always emits triangles
  counter-clockwise in the projected frame, so the normalising swap never
  fires. Kept as a guard against an earcut convention change; equivalent
  today, so not listed. "triangle winding always flipped" covers the swap.

#254 (warped faces triangulated, their warp reported Certified): 18 / 18,
after adding fixtures for a warped hole and for two warped faces in both
orders, written for the "outer ring only" and "worst across faces"
mutants before the first run. Not listed, as equivalent at test
resolution: dropping the 16-ulp rounding allowance (`Some(worst)`), which
no fixture can resolve.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/execution/compile/src/planar.rs"
C = "crates/execution/compile/src/compiler.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "authored_polygons"],
    ["-p", "axiolid-mesh-compile", "--lib", "planar"],
]

MUTANTS = [
    ("planarity check skipped", P,
     "    let warp = if worst > linear {",
     "    let warp = if worst > linear && false {"),
    ("planarity against zero, not tolerance", P,
     "    let warp = if worst > linear {", "    let warp = if worst > 0.0 {"),
    ("planar callers no longer refuse a warp", P,
     "    triangulate_in_fit_plane(rings, linear, true).map(|(indices, _)| indices)",
     "    triangulate_in_fit_plane(rings, linear, false).map(|(indices, _)| indices)"),
    ("authored faces refuse a warp again", P,
     "    triangulate_in_fit_plane(rings, linear, false)\n",
     "    triangulate_in_fit_plane(rings, linear, true)\n"),
    ("warp halved", P,
     "        Some(worst + 16.0 * f64::EPSILON * reach)",
     "        Some(0.5 * worst)"),
    ("rounding allowance subtracted", P,
     "        Some(worst + 16.0 * f64::EPSILON * reach)",
     "        Some(worst - 16.0 * f64::EPSILON * reach)"),
    ("warp measured on the outer ring only", P,
     "    for ring in rings {\n        for &p in *ring {\n            let offset = p - centroid;",
     "    for ring in &rings[..1] {\n        for &p in *ring {\n            let offset = p - centroid;"),
    ("warp reported as proven", C,
     "                                crate::deviation::DeviationBound::Certified(warp),",
     "                                crate::deviation::DeviationBound::Proven(warp),"),
    ("warp folded into the exact authored contribution", C,
     "                                crate::deviation::WARPED_AUTHORED_FACE,",
     "                                \"\","),
    ("warp not reported", C,
     "                            None => built,\n                            Some(warp) => {",
     "                            _ => built,\n                            Some(warp) if false => {"),
    ("worst warp across faces dropped", C,
     "                warp = Some(warp.map_or(face_warp, |w| w.max(face_warp)));",
     "                warp = Some(face_warp);"),
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
     "    let mut indices = earcut_projected(&flat, &hole_starts);\n    if indices.len() % 3 != 0 {",
     "    let mut indices = earcut_projected(&flat[..hole_starts.first().copied().unwrap_or(flat.len())], &[]);\n    if indices.len() % 3 != 0 {"),
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


def run():
    for target in TARGETS:
        code = subprocess.run(["cargo", "test", "-q", *target], cwd=ROOT,
                              capture_output=True, text=True, timeout=900).returncode
        if code != 0:
            return code
    return 0


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found ({original.count(old)})"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run()
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
