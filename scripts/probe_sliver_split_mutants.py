"""Mutation probe for the planar faces' sliver split (#278).

`planar::split_invented_edges` drops the clipper's slivers along nearly
straight runs and splits the invented edges they leave. Before #278 it
re-tested the noise band for the corners of the dropped slivers; the band
is not transitive, so a cap corner could end up in no triangle and a
closed prism meshed open. Each mutant weakens one decision of the fix:
the chained sliver run that the edge left behind takes, or the check of
the split against the clipper's certified cover. `tests/near_collinear_corner.rs`
or the `planar` unit tests must fail on it; a build error also counts as
killed.

First run (#278): 9 / 11. Two survivors: "winding not compared" (every
fixture the split folds also repeats a directed edge, now pinned by a
direct `same_cover` test) and "a directed edge used twice allowed": a
split whose triangles all wind one way around the cover's boundary
covers each point once, so it cannot repeat one; that check was
redundant and went. 10 / 10 killed after.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/execution/compile/src/planar.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--lib", "planar"],
    ["-p", "axiolid-mesh-compile", "--test", "near_collinear_corner"],
    ["-p", "axiolid-mesh-compile", "--test", "authored_polygons"],
    ["-p", "axiolid-mesh-compile", "--test", "brep_tessellation"],
    ["-p", "axiolid-mesh-compile", "--test", "shared_edge_weld"],
]

MUTANTS = [
    ("exposed edges re-test the band, as before #278", P,
     "            if slivers.contains_key(&(q, p)) {",
     "            if false && slivers.contains_key(&(q, p)) {"),
    ("exposed runs ignored when splitting", P,
     "            if let Some(run) = exposed.get(&(p, q)) {\n                return run.clone();",
     "            if let Some(run) = exposed.get(&(p, q)).filter(|_| false) {\n                return run.clone();"),
    ("sliver run not chained: the apex alone", P,
     "                    stack.push(Step::Edge(m, q));\n                    stack.push(Step::Corner(m));\n                    stack.push(Step::Edge(p, m));",
     "                    stack.push(Step::Corner(m));"),
    ("sliver run in reverse order", P,
     "                    stack.push(Step::Edge(m, q));\n                    stack.push(Step::Corner(m));\n                    stack.push(Step::Edge(p, m));",
     "                    stack.push(Step::Edge(p, m));\n                    stack.push(Step::Corner(m));\n                    stack.push(Step::Edge(m, q));"),
    ("sliver recorded against its long side reversed", P,
     "                slivers.insert((t[(i + 1) % 3], t[(i + 2) % 3]), t[i]);",
     "                slivers.insert((t[(i + 2) % 3], t[(i + 1) % 3]), t[i]);"),
    ("split kept without the cover check", P,
     "    if same_cover(flat, indices, &out) {",
     "    if true {"),
    ("cover check: winding not compared", P,
     "    wound && net(cover) == net(split)",
     "    net(cover) == net(split)"),
    ("cover check: winding assumed counter-clockwise", P,
     "        .sum::<Scalar>()\n        .signum();",
     "        .sum::<Scalar>()\n        .signum()\n        .max(1.0);"),
    ("cover check: boundary not compared", P,
     "    wound && net(cover) == net(split)",
     "    wound"),
    ("cover check: boundary compared without direction", P,
     "                *net.entry((p.min(q), p.max(q))).or_insert(0) += if p < q { 1 } else { -1 };",
     "                *net.entry((p.min(q), p.max(q))).or_insert(0) += 1;"),
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
