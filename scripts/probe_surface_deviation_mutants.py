"""Mutation probe for the surface-to-mesh bound of curved bodies (#231, #232,
#245).

Each mutant gives a direction of a doubly curved surface the whole chord
budget again, sizes a step at the wrong radius, drops a term of the bound
(the walls' twist, the taper, the section's reach, the directrix's own
sagitta), loosens the quad lemma, tilts a swept end, or clamps a budget it
must refuse, and must turn a test red. The pipe mutants (#232) do the same
to swept disks along chains of segments and arcs: they drop a term of a
piece's bound or the joint's measured shift, stop carrying the frame across
a joint, loosen the joint tolerance, cut fillets wrongly or let them
overrun their segments, and accept folded bends, horn tori or closed
polylines. The mitre mutants (#245) chord the ring without the mitre's
stretch, skip one side's projection, tilt the mitre plane, project
obliquely, drop the reach check or one end of it, and mitre reversals and
corners beside arcs.
`surface_deviation.rs` samples the exact surfaces against the compiled
meshes; the `loft` unit tests check the quad lemma against a brute-force
distance.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/execution/compile/src/compiler.rs"
D = "crates/execution/compile/src/directrix.rs"
R = "crates/algorithms/construction/construct/src/revolve.rs"
S = "crates/algorithms/construction/construct/src/sweep.rs"
L = "crates/algorithms/construction/construct/src/loft.rs"
P = "crates/algorithms/reference/src/primitive.rs"
PIPE = "crates/algorithms/construction/construct/src/pipe.rs"
PIECES = "crates/execution/compile/src/directrix/pieces.rs"
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
    ('cylinder segments clamped', P, '"cylinder height")?;\n    let n = segments_within(r, tol)?;', '"cylinder height")?;\n    let n = segments(r, tol);', [SURFACE]),
    ('cone segments clamped', P, '"cone height")?;\n    let n = segments_within(r, tol)?;', '"cone height")?;\n    let n = segments(r, tol);', [SURFACE]),
    ('pipe disk chorded to the whole budget', PIPE, 'sweep::disk_rings(radius, inner_radius, 0.5 * chord * (0.5 * sharpest).cos())?;', 'sweep::disk_rings(radius, inner_radius, chord * (0.5 * sharpest).cos())?;', [SURFACE]),
    ('pipe bends sampled to the whole budget', PIPE, '    let share = 0.5 * chord;', '    let share = chord;', [SURFACE]),
    ('pipe section reach ignored', PIPE, 'sweep::span_bound(&f[0], &f[1], b, reach)', 'sweep::span_bound(&f[0], &f[1], b, 0.0)', [SURFACE]),
    ('pipe bends never refined', PIPE, '        if worst <= share {', '        if true {', [SURFACE]),
    ('pipe joint shift ignored', PIPE, '                let moved = if k == 0 { shift } else { 0.0 };', '                let moved = 0.0;', [SURFACE]),
    ('pipe frame not carried across a joint', PIPE, '            Some(frame) => least_rotation(frame.x, frame.x.cross(frame.y), tangent),', '            Some(_) => seed_reference(tangent),', [SURFACE]),
    ('pipe corner tolerance quadrupled', PIPE, '    2.0 * (chord / (8.0 * radius)).min(1.0).asin()', '    2.0 * (chord / (2.0 * radius)).min(1.0).asin()', [SURFACE]),
    ('fillet cut by the sine', PIPE, '        let cut = fillet * (0.5 * theta).tan();', '        let cut = fillet * (0.5 * theta).sin();', [SURFACE]),
    ('fillet overrunning its segment accepted', PIPE, '        if needed > length + slack {', '        if false {', [SURFACE]),
    ('disk wider than the fillet accepted', PIPE, '        if radius > fillet {', '        if false {', [SURFACE]),
    ('disk as wide as the fillet not named', PIPE, '        if radius == fillet {', '        if false {', [SURFACE]),
    ('bend as tight as the disk accepted', PIPE, '                if bend == radius {', '                if false {', [SURFACE]),
    ('mitre ring without the stretch', PIPE, 'sweep::disk_rings(radius, inner_radius, 0.5 * chord * (0.5 * sharpest).cos())?;', 'sweep::disk_rings(radius, inner_radius, 0.5 * chord)?;', [SURFACE]),
    ('mitre projected at the start only', PIPE, 'for (k, normal) in [(0, ends[0]), (last, ends[1])] {', 'for (k, normal) in [(0, ends[0]), (last, None)] {', [SURFACE]),
    ('mitre plane square to the incoming leg', PIPE, '        mitres.push(Some(bisector.normalize()));', '        mitres.push(Some(u));', [SURFACE]),
    ('mitre projected obliquely', PIPE, '.map(|q| *q - tangent * ((*q - frame.origin).dot(normal) / along))', '.map(|q| *q - tangent * (*q - frame.origin).dot(normal))', [SURFACE]),
    ('mitre reach not checked', PIPE, '        if reach >= length {', '        if false {', [SURFACE]),
    ('mitre reach of one end only', PIPE, '        let g = lean(ends[1]) - lean(ends[0]);', '        let g = lean(ends[1]) - u;', [SURFACE]),
    ('reversal mitred', PIPE, '        if theta >= core::f64::consts::PI - 1e-9 || bisector.length() <= 1e-9 {', '        if false {', [SURFACE]),
    ('corner beside an arc mitred', PIPE, '        if !matches!(\n            (pair[0], pair[1]),', '        if false && !matches!(\n            (pair[0], pair[1]),', [SURFACE]),
    ('reversed arc keeps its turn', PIPE, '                start: self.end_point(),\n                angle: -angle,', '                start: self.end_point(),\n                angle,', [SURFACE]),
    ('closed polyline accepted', PIECES, '            if polyline.closed {', '            if false {', [SURFACE]),
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
