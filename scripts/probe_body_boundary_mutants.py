"""Mutation probe for boundary distance and Hausdorff distance between
bodies of several exact solids (#229; `exact_bodies.rs`, its `contact.rs`,
and `exact_hausdorff/cut.rs`).

Each mutant either lets a body whose items are not shown apart, touching
without a shared patch, or in exact face contact through to the
measurement (so a refusal test goes red), refuses a layout that is sound
(edge contact, a column threaded through a frame, a column on its
footing), cuts a shared patch wrongly (not at all, or a partner's hole
with it), names the wrong item for a witness, or loses the per-item
support seed, the cut faces' bound or the two-sided combination. The
closed forms are in `construct/tests/boundary_bodies.rs`.

The sum of two touching planes is needed only where no coordinate axis
lies in the wedge of directions a shared edge is extreme in, so its
fixture is turned and tilted; a footing tilted under a block resting on
its edge is what tells "one item has no patch" from "neither has".

Not listed: a B-spline face that reaches the separating plane never read
as holding a patch of it (no fixture has a B-spline face in a contact
plane); the margin of the enclosing-box test (every nested fixture sits
well inside); and, in the cut faces' bound, `zeta` left out or the check
that a partner's rim is on the boundary skipped. Both make the bound
unsound, but only where a partner's outline moves by more than rounding
between the two bodies or lies under a third item, and every fixture
that reaches the bound is a translate (where `zeta` is rounding) with
partners whose rims are free.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/algorithms/query/measure/src/exact_bodies.rs"
C = "crates/algorithms/query/measure/src/exact_bodies/contact.rs"
K = "crates/algorithms/query/measure/src/exact_hausdorff/cut.rs"
D = "crates/algorithms/query/measure/src/exact_distance.rs"
H = "crates/algorithms/query/measure/src/exact_hausdorff.rs"
TESTS = [
    ["-p", "axiolid-construct", "--test", "boundary_bodies"],
]

MUTANTS = [
    ("touching ranges read as apart", B,
     "        if gap > 0.0 {",
     "        if gap > -slack {"),
    ("overlapping ranges read as touching", B,
     "        if gap >= -slack {",
     "        if gap >= -1.0 {"),
    ("a planar face in the plane never holds a patch", B,
     "            Surface::Plane(_) if lo >= band.0 && hi <= band.1 => {",
     "            Surface::Plane(_) if false && lo >= band.0 && hi <= band.1 => {"),
    ("a patch-free contact needs both items patch-free", B,
     "        if !may_hold_patch(a, plane)? || !may_hold_patch(b, plane)? {",
     "        if !may_hold_patch(a, plane)? && !may_hold_patch(b, plane)? {"),
    ("no plane along the sum of two touching ones", B,
     "        for j in i + 1..base {",
     "        for j in base..base {"),
    ("faces off one axis plane read as in exact contact", B,
     "                _ => off_plane = true,",
     "                _ => {}"),
    ("near contact read as an overlap", B,
     "            return Ok(Layout::NearlyShareFace { gap: plane.gap });\n        }\n    }\n    Ok(Layout::Undecided)",
     "            return Ok(Layout::Undecided);\n        }\n    }\n    Ok(Layout::Undecided)"),
    ("boundaries apart read as items apart", B,
     "    if boundaries_apart(a, b, tolerance)? && !nested(a, b)? {",
     "    if boundaries_apart(a, b, tolerance)? {"),
    ("boundaries always read as apart", B,
     "    Ok(found.bounds.lower > 0.0)",
     "    Ok(true)"),
    ("one item sticking out of the other suffices", B,
     "    Ok(!sticks_out(a, b)? || !sticks_out(b, a)?)",
     "    Ok(!sticks_out(a, b)?)"),
    ("several solids per item read as one", B,
     "    if !one_solid(a.brep) || !one_solid(b.brep) {\n        return Ok(true);",
     "    if !one_solid(a.brep) || !one_solid(b.brep) {\n        return Ok(false);"),
    ("an element named after the next item", B,
     "            .saturating_sub(1)",
     "            .saturating_sub(0)"),
    ("one support seed for the whole body", B,
     "            .map(|k| self.edges[k]..self.edges.get(k + 1).copied().unwrap_or(total))",
     "            .map(|_| 0..total)"),
    ("two-sided upper from the backward side", B,
     "    distance.bounds.upper = forward.bounds.upper.max(backward.bounds.upper);",
     "    distance.bounds.upper = backward.bounds.upper;"),
    # Face contact.
    ("a shared patch never cut", C,
     "    if !overlap {",
     "    if true {"),
    ("a shared patch left on the free region", C,
     "        .regions(|flags| inside(flags, &mine) && !covered(flags))",
     "        .regions(|flags| inside(flags, &mine))"),
    ("a face's holes ignored when cutting", C,
     "        outer.iter().any(|&k| flags[k]) && !holes.iter().any(|&k| flags[k])",
     "        outer.iter().any(|&k| flags[k])"),
    ("a cut face bounded without its partners' counterparts", K,
     "            let Some(residue) = best else {\n                return Ok(None);",
     "            let Some(residue) = best else {\n                continue;"),
    ("cut faces bounded only as they were cut", K,
     "            let bound = m + e.max(zeta) + zeta;",
     "            let bound = Scalar::INFINITY;"),
    # Witness naming.
    ("distance witnesses' elements swapped", D,
     "                best = Some((d, wa, wb, (ea.shape, eb.shape)));",
     "                best = Some((d, wa, wb, (eb.shape, ea.shape)));"),
    ("Hausdorff witnesses' elements swapped", H,
     "                        state.witness = Some((w, found.nearest, element.shape, found.on));",
     "                        state.witness = Some((w, found.nearest, found.on, element.shape));"),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
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
        for target in TESTS:
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
