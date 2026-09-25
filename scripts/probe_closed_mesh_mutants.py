"""Mutation probe for the #170 closed-mesh fix: edge splitting and closure.

Each mutant weakens one decision in `planar::split_invented_edges`,
`planar::split_triangle`, `planar::strictly_between`, `planar::collinear_band`
or `compiler::authored_mesh_is_closed`. The compile crate's planar unit
tests plus `authored_polygons` and `surface_models` must fail on it. A build
error also counts as killed.

First run: 11 / 12. The survivor, "ear test ignores collinear ears", was
dead code: `split_triangle` clips the first corner whose removal leaves a
ring not all on one line, and a Python model of that loop over every mix of
0..4 corners per side (125 cases) gave identical triangles with and without
the per-ear collinearity test, never a degenerate ear and never the
fallback. The test was removed; the mutant with it.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/execution/compile/src/planar.rs"
C = "crates/execution/compile/src/compiler.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--lib", "planar"],
    ["-p", "axiolid-mesh-compile", "--test", "authored_polygons"],
    ["-p", "axiolid-mesh-compile", "--test", "surface_models"],
]

MUTANTS = [
    ("split pass skipped entirely", P,
     "    if triangles == 0 || triangles.saturating_mul(3).saturating_mul(n) > MAX_SPLIT_WORK {",
     "    if true || triangles == 0 || triangles.saturating_mul(3).saturating_mul(n) > MAX_SPLIT_WORK {"),
    ("authored ring edges also split", P,
     "            if authored.contains(&key) || on_edge.contains_key(&key) {",
     "            if on_edge.contains_key(&key) {"),
    ("slivers kept", P,
     "        .filter(|t| !is_sliver(t))",
     "        .filter(|_t| true)"),
    ("thin authored triangle dropped as sliver", P,
     "            !authored.contains(&(a.min(b), a.max(b)))\n                && strictly_between(",
     "            true\n                && strictly_between("),
    ("corners on an edge left unsorted", P,
     "            found.sort_by(|x, y| x.0.total_cmp(&y.0));",
     "            found.sort_by(|x, y| y.0.total_cmp(&x.0));"),
    ("run not reversed for a reversed side", P,
     "            if p > q {\n                run.reverse();",
     "            if p < q {\n                run.reverse();"),
    ("collinear band ignores tolerance", P,
     "    turn.abs() <= (1e-9 * length).max(noise * length.sqrt())",
     "    turn.abs() <= 1e-9 * length"),
    ("collinear band far too wide", P,
     "pub(crate) fn collinear_band(linear: Scalar) -> Scalar {\n    1.0e-3 * linear",
     "pub(crate) fn collinear_band(linear: Scalar) -> Scalar {\n    1.0e3 * linear"),
    ("endpoints counted as between", P,
     "    if !(along > 0.0 && along < length) {",
     "    if !(along >= 0.0 && along <= length) {"),
    ("closure from the area-thresholded audit", C,
     "    let adjacency = axiolid_mesh::EdgeAdjacency::build(mesh);\n    adjacency.edge_count() > 0 && adjacency.is_closed_two_manifold()",
     "    axiolid_mesh::audit_mesh(mesh, axiolid_core::Tolerance::MILLIMETRE).is_closed_two_manifold()"),
    ("closure ignores winding", C,
     "    adjacency.edge_count() > 0 && adjacency.is_closed_two_manifold()",
     "    adjacency.edge_count() > 0 && adjacency.boundary_edges().next().is_none() && adjacency.non_manifold_edges().next().is_none()"),
]


def killed(target):
    try:
        return subprocess.run(["cargo", "test", "-q", *target], cwd=ROOT,
                              capture_output=True, text=True, timeout=900).returncode != 0
    except subprocess.TimeoutExpired:
        return True


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found ({original.count(old)})"
    path.write_text(original.replace(old, new))
    try:
        dead = any(killed(t) for t in TARGETS)
    finally:
        path.write_text(original)
    print(f"{'killed' if dead else 'SURVIVED':8} {name}", flush=True)
    if not dead:
        survivors.append(name)

print(f"{len(MUTANTS) - len(survivors)} / {len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
