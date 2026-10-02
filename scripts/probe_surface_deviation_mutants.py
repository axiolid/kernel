"""Mutation probe for the surface-to-mesh bound of doubly curved bodies (#231).

Each mutant gives a direction of a doubly curved surface the whole chord
budget again, sizes a step at the wrong radius, drops a term of the bound
(the walls' twist, the taper, the section's reach, the directrix's own
sagitta), loosens the quad lemma, tilts a swept end, or clamps a budget it
must refuse, and must turn a test red. `surface_deviation.rs` samples the
exact surfaces against the compiled meshes; the `loft` unit tests check the
quad lemma against a brute-force distance.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/execution/compile/src/compiler.rs"
D = "crates/execution/compile/src/directrix.rs"
R = "crates/algorithms/construction/construct/src/revolve.rs"
S = "crates/algorithms/construction/construct/src/sweep.rs"
L = "crates/algorithms/construction/construct/src/loft.rs"
P = "crates/algorithms/reference/src/primitive.rs"
SURFACE = ["-p", "axiolid-mesh-compile", "--test", "surface_deviation"]
QUAD = ["-p", "axiolid-construct", "--lib", "loft"]

MUTANTS = [
    ('revolution profile chorded to the whole budget', C, '`axiolid_construct::revolve`).\n                let half = 0.5 * chord_error(options);', '`axiolid_construct::revolve`).\n                let half = chord_error(options);', [SURFACE]),
    ('revolution turned to the whole budget', C, '    Tolerance::new(0.5 * chord_error(options), options.tolerance().angular())', '    Tolerance::new(chord_error(options), options.tolerance().angular())', [SURFACE]),
    ('turn sized at the tube centre', R, '    let n = bounded_steps(rings, axis_origin, dir, angle, max_r, tolerance.linear())?;', '    let n = bounded_steps(rings, axis_origin, dir, angle, 0.8 * max_r, tolerance.linear())?;', [SURFACE]),
    ('revolution wall twist ignored', R, '        if sagitta + twist <= budget {', '        if sagitta <= budget {', [SURFACE]),
    # The revolution's step cap is guarded twice (the sagitta count and the
    # growth loop both refuse), so no single mutant there can pass the
    # refusal test; the primitive's single guard is probed below.
    ('tapered taper speed ignored', S, '    let curvature = angle * angle * max_r + 2.0 * angle.abs() * travel;', '    let curvature = angle * angle * max_r;', [SURFACE]),
    ('tapered wall twist ignored', S, '        if chord + twist <= budget {', '        if chord <= budget {', [SURFACE]),
    ('swept disk chorded to the whole budget', S, '    let rings = disk_rings(radius, inner_radius, 0.5 * chord)?;', '    let rings = disk_rings(radius, inner_radius, chord)?;', [SURFACE]),
    ('section reach ignored', S, '    directrix + reach * (1.0 - (0.5 * rotation).cos())', '    directrix', [SURFACE]),
    ('directrix sagitta ignored', S, '    let directrix = budget.min(0.5 * chord * (0.25 * turn).tan());', '    let directrix = 0.0;', [SURFACE]),
    ('end sections square to the end chords', S, '                ends.map_or(path[1] - path[0], |[start, _]| start)', '                path[1] - path[0]', [SURFACE]),
    ('smooth directrix never refined', S, '        if !smooth {\n            return Ok(stations);', '        if true {\n            return Ok(stations);', [SURFACE]),
    ('reversed trim keeps its tangents', D, '    (points, ends.map(|[start, end]| [-end, -start]))', '    (points, ends.map(|[start, end]| [end, start]))', [SURFACE]),
    ('conic end tangents not reported', D, '    Some([unit(start)?, unit(end)?])', '    Some([unit(start)?, unit(end)?]).filter(|_| false)', [SURFACE]),
    ('quad flattening bound quartered', L, '            best = best.min(g.abs());', '            best = best.min(0.25 * g.abs());', [QUAD, SURFACE]),
    ('quad mid-plane bound halved', L, '            best = best.min(2.0 * e.abs());', '            best = best.min(e.abs());', [QUAD, SURFACE]),
    ('quad convexity not checked', L, '        if convex_projection([a, b, d - n * g, c], n) {', '        if true {', [QUAD, SURFACE]),
    ('quad twist bound halved', L, '    let twist = (a - b - c + d).length() / 4.0;', '    let twist = (a - b - c + d).length() / 8.0;', [QUAD, SURFACE]),
    ('sphere stacks rounded down', P, '    let stacks = n.div_ceil(2).max(2);', '    let stacks = (n / 2).max(2);', [SURFACE]),
    ('sphere at the whole budget', P, '    let n = segments_within(r, 0.5 * tol)?;', '    let n = segments_within(r, tol)?;', [SURFACE]),
    ('torus equator at the whole budget', P, '    let n = segments_within(big + r, 0.5 * tol)?;', '    let n = segments_within(big + r, tol)?;', [SURFACE]),
    ('primitive segment cap clamps', P, '    if n == 4096 && radius', '    if false && radius', [SURFACE]),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode

survivors = []
for name, rel, old, new, targets in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    code = 0
    try:
        for target in targets:
            try:
                code = run(target)
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
