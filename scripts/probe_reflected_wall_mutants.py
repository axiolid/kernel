"""Mutation probe for booleans cut by a reflected curved wall (#288).

A face whose angles lie a turn away from the surface inverse's -- the
cylinder walls of a reflected (downward-extruded) prism read `2 pi - u` --
has its boundary edges cut where a section ends. The cut's pcurve
parameter must land on the cut point: a turn of a line pcurve's parameter
is another point. Each mutant undoes part of that check and must turn a
test red: brep-boolean's `tests/reflected_walls.rs` (downward arc and
circle prisms against a box, every operator, reflected spheres, cones and
tori) or mesh-compile's `bounded_half_space_arcs.rs`, whose clip by the
opposite side now reflects its arc prism instead of mirroring the profile.

Equivalent mutants, deliberately not listed: trying the turns of the
pcurve parameter in another order (at most one lands inside the span),
and comparing the landing against the cut point instead of the edge's
point there (they agree within the incidence tolerance).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/construction/brep-boolean/src/split.rs"
TESTS = [
    ["-p", "axiolid-brep-boolean", "--test", "reflected_walls"],
    ["-p", "axiolid-mesh-compile", "--test", "bounded_half_space_arcs"],
]

MUTANTS = [
    ("any turn inside the span taken (the #288 defect)", S,
     "                    if lands(c)? {",
     "                    if (plo - slack..=phi + slack).contains(&c) {"),
    ("the span not checked", S,
     "                    if c < plo - slack || c > phi + slack {\n"
     "                        return Ok(false);",
     "                    if false {\n"
     "                        return Ok(false);"),
    ("the landing demanded exactly", S,
     "                    Ok((back - on).length() <= 2.0 * report::floored(tolerance).linear())",
     "                    Ok((back - on).length() <= 0.0)"),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", "--all-features", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode

survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        code = 0
        for t in TESTS:
            try:
                code = run(t)
            except subprocess.TimeoutExpired:
                code = -1
            if code != 0:
                break
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
