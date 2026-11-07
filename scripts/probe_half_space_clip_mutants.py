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
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/execution/compile/src/exact/clip.rs"
D = "crates/execution/compile/src/exact/boolean.rs"
H = "crates/algorithms/construction/construct/src/half_space.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "exact_half_space_clip"],
]

MUTANTS = [
    ("the clip not dispatched", D,
     "        if let Some((subject, region)) = self.clip_of(left, right, operator)? {",
     "        if let Some((subject, region)) = self.clip_of(left, right, operator)?.filter(|_| false) {"),
    ("agreement ignored", C,
     "        let side = if half_space.agreement {",
     "        let side = if true {"),
    ("the in-plane y mirrored", C,
     "        let y = normal.cross(x) * flip;",
     "        let y = x.cross(normal) * flip;"),
    ("the boundary anchored off the plane", H,
     "        let anchor = origin - normal * (origin - plane.origin).dot(normal);",
     "        let anchor = origin;"),
    ("the operator dropped", C,
     "boolean_with_report(&subject, &solid, operator, tolerance)",
     "boolean_with_report(&subject, &solid, BooleanOperator::Difference, tolerance)"),
    ("the tool left unplaced", C,
     "            .transformed(&(*placement * turned))",
     "            .transformed(&turned)"),
    ("the envelope not brought into the half-space's frame", C,
     "        let corners = envelope.mapped(&inverse).corners();",
     "        let corners = envelope.corners();"),
    ("a subject missing the half-space emptied", C,
     "            Cover::Misses if keeps => return Ok(Body::Empty(CLIP_EMPTIED)),\n            Cover::Misses => return Ok(Body::Solid(subject)),",
     "            Cover::Misses if keeps => return Ok(Body::Solid(subject)),\n            Cover::Misses => return Ok(Body::Empty(CLIP_EMPTIED)),"),
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
    ("unions with a half-space let through", D,
     "        if self.region(left)?.is_some() || self.region(right)?.is_some() {",
     "        if (self.region(left)?.is_some() || self.region(right)?.is_some())\n            && operator != BooleanOperator::Union\n        {"),
    ("a scaled half-space not named", C,
     "        TransformError::NotRigid => {",
     "        TransformError::NotRigid if false => {"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800)
