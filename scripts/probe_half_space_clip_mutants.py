"""Mutation probe for exact half-space clipping (#234).

Each mutant breaks one decision of the clip path -- which side the
half-space is, how the finite tool is framed, placed and sized, what the
envelope decides alone, the dispatch and its refusals -- and must turn a
test red.

Equivalent mutants, deliberately not listed: dropping the shortcut for a
subject wholly inside an unbounded half-space (the general boolean then
empties it, the same `Degenerate`); bounding every boolean's envelope by
both operands (sound, only looser); keeping the in-plane `y` unturned for
the opposite side (the tool's frame is then a reflection, which
`ExactBRep::transformed` also places exactly, and the boundary is read in
the same unturned frame).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/execution/compile/src/exact/clip.rs"
D = "crates/execution/compile/src/exact/boolean.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "exact_half_space_clip"],
]

MUTANTS = [
    ("the clip not dispatched", D,
     "        if let Some(tool) = self.clip_tool(right)? {",
     "        if let Some(tool) = self.clip_tool(right)?.filter(|_| false) {"),
    ("agreement ignored", C,
     "        let side = if half_space.agreement {",
     "        let side = if true {"),
    ("the in-plane y mirrored", C,
     "        let y = normal.cross(x) * flip;",
     "        let y = x.cross(normal) * flip;"),
    ("the boundary anchored off the plane", C,
     "                let anchor = frame.origin() - normal * (frame.origin() - origin).dot(normal);",
     "                let anchor = frame.origin();"),
    ("the operator dropped", C,
     "        boolean(&subject, &solid, operator, tolerance)",
     "        boolean(&subject, &solid, BooleanOperator::Difference, tolerance)"),
    ("the tool left unplaced", C,
     "            .transformed(&(*placement * frame))",
     "            .transformed(&frame)"),
    ("the envelope not brought into the half-space's frame", C,
     "        let corners = envelope.mapped(&inverse).corners();",
     "        let corners = envelope.corners();"),
    ("a subject missing the half-space emptied", C,
     "            return Ok(Err(if keeps {\n                Decided::Empty\n            } else {\n                Decided::Subject\n            }));",
     "            return Ok(Err(if keeps {\n                Decided::Subject\n            } else {\n                Decided::Empty\n            }));"),
    ("a subject touching the plane decided without the boolean", C,
     "        if far < 0.0 {",
     "        if far < 1e-3 {"),
    ("no margin past the envelope", C,
     "        let margin = 0.25 * (envelope.hi - envelope.lo).length() + 4.0 * tolerance.linear();",
     "        let margin = 0.0;"),
    ("a circle bounded without its radius", C,
     "                let r = reach(circle.frame.x, circle.frame.y) * circle.radius;",
     "                let r = reach(circle.frame.x, circle.frame.y) * 0.0;"),
    ("vertices left out of the envelope", C,
     "        envelope.include(vertex.position);",
     "        let _ = vertex;"),
    ("unions with a half-space let through", C,
     "        if operator == BooleanOperator::Union {",
     "        if false {"),
    ("a scaled half-space not named", C,
     "                TransformError::NotRigid => {",
     "                TransformError::NotRigid if false => {"),
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
