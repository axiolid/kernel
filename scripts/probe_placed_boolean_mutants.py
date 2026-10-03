"""Mutation probe for exact differences of placed operands (#228).

Each mutant reintroduces a way the placed-opening path failed or could
fail -- a rounding residue read exactly, a tangent double root split into
a sliver, the compiler's dispatch or refusals loosened -- and must turn a
test red.

Equivalent mutants, deliberately not listed: dropping the periodic wrap
merge in `merge_close` (no opening here puts cuts either side of a
conic's parameter origin); swapping the parallel/perpendicular branch
order in `plane_cylinder_within_rounding` (a direction cannot be both).

Not reached by these inputs, deliberately not listed: dropping the cuts
where a section line or conic leaves the two faces' common box. They
matter only when a kept cut lies outside that box and a skipped edge's
crossing lies between it and the next kept cut; no opening here produces
one, and the pieces are then still bounded by kept cuts inside the box.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/construction/brep-boolean/src/section.rs"
F = "crates/algorithms/construction/brep-boolean/src/split.rs"
C = "crates/execution/compile/src/exact/boolean.rs"
B = "crates/algorithms/construction/brep-boolean/src/bounds.rs"
TESTS = [
    ["-p", "axiolid-brep-boolean", "--lib"],
    ["-p", "axiolid-brep-boolean", "--test", "openings"],
    ["-p", "axiolid-brep-boolean", "--test", "contact"],
    ["-p", "axiolid-mesh-compile", "--test", "exact_placed_boolean"],
]

MUTANTS = [
    ("face boxes not enlarged by the tolerance", B,
     "    let pad = tolerance.linear() + 8.0 * Scalar::EPSILON * (1.0 + size);",
     "    let pad = -tolerance.linear();"),
    ("a cylinder's box without its radius", B,
     "            (a, b, radial(c.frame.x, c.frame.y, c.radius, c.radius))",
     "            (a, b, Vec3::ZERO)"),
    ("a curve along an edge only when exactly contained", S,
     "                    && runs_along(curve, edge_curve, span, tolerance)?",
     "                    && false"),
    ("a lost tangent root not recovered from the edge", S,
     "                        if out.len() == before {",
     "                        if false {"),
    ("split tangent roots kept apart on sections", S,
     "                let mut cuts = merge_close(branch, cuts, tolerance)?;",
     "                let mut cuts = cuts;"),
    ("a near-parallel plane and cylinder read exactly", S,
     "    if extent.is_none_or(|extent| along.abs() * extent > eps) {",
     "    if along != 0.0 {"),
    ("a plane read parallel however far it turns over the faces", S,
     "    if extent.is_none_or(|extent| along.abs() * extent > eps) {",
     "    if along.abs() > eps {"),
    ("a near-perpendicular plane and cylinder read exactly", S,
     "    if cylinder.radius * axis.cross(normal).length() <= eps {",
     "    if false {"),
    ("a near-tangent plane cuts two rulings", S,
     "    let offsets: Vec<Scalar> = if (distance.abs() - r).abs() <= eps {",
     "    let offsets: Vec<Scalar> = if distance.abs() == r {"),
    ("a chain of near cuts merged past the tolerance", S,
     "                    return Err(BooleanError::NearCoincidence);",
     "                    let _ = 0;"),
    ("cuts merged up to ten tolerances apart", S,
     "            Some((sum, count, first, last)) if (point - *last).length() <= eps => {",
     "            Some((sum, count, first, last)) if (point - *last).length() <= 10.0 * eps => {"),
    ("coincident supports read at twenty tolerances", "crates/algorithms/construction/brep-boolean/src/support.rs",
     "    let eps = tolerance.linear();",
     "    let eps = 20.0 * tolerance.linear();"),
    ("an edge read along a curve from sampled points", S,
     "        _ => Ok(false),\n    }\n}",
     "        _ => Ok(true),\n    }\n}"),
    ("a boundary use split next to its own end", F,
     "        if taken.iter().any(|q: &Point3| (point - *q).length() <= eps) {",
     "        if false {"),
    ("unions of placed operands let through", C,
     "        if operator != BooleanOperator::Difference {",
     "        if false {"),
    ("a boolean tool let through", C,
     "        if tool != Placed::Extrusion {",
     "        if false {"),
    ("any operand read as an extrusion", C,
     "                _ => {\n                    return Err(unsupported(",
     "                _ if false => {\n                    return Err(unsupported("),
    ("the operator dropped for the general boolean", C,
     "        boolean(&subject, &tool, operator, self.options.tolerance())",
     "        boolean(&subject, &tool, BooleanOperator::Union, self.options.tolerance())"),
    ("operands swapped", C,
     "        boolean(&subject, &tool, operator, self.options.tolerance())",
     "        boolean(&tool, &subject, operator, self.options.tolerance())"),
    ("an emptied wall refused instead of degenerate", C,
     "        BooleanError::EmptyResult => {",
     "        BooleanError::EmptyResult if false => {"),
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
        try:
            code = max(run(t) for t in TESTS)
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
