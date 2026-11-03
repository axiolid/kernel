"""Mutation probe for the #415 cap fix: slivers flipped away, flat ones
refused.

Each mutant weakens one decision in `flip::is_sliver`, `flip::is_flat`,
`flip::unsliver`, `profile::triangulate_cap` or
`extrude::extrude_profile`. The compile crate's `cap_collinear` tests (the
issue's 23-point ring, both placements, the randomised sweep, the controls
and the refusal), the construct crate's `cap_slivers` tests (every start
vertex, the over-refusal ring), its `flip` unit tests or its `extrusion_volume` tests (a fine
chord's thin triangles must still extrude) must fail on it. A build error also counts as killed.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROFILE = "crates/algorithms/construction/construct/src/profile.rs"
EXTRUDE = "crates/algorithms/construction/construct/src/extrude.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "cap_collinear"],
    ["-p", "axiolid-construct", "--lib", "flip"],
    ["-p", "axiolid-construct", "--test", "cap_slivers"],
    ["-p", "axiolid-construct", "--test", "extrusion_volume"],
]
FLIP = "crates/algorithms/construction/construct/src/ring_triangulation/flip.rs"
RINGS = "crates/algorithms/construction/construct/src/ring_triangulation.rs"

MUTANTS = [
    ("flat band far too narrow", FLIP,
     "    double_area <= 1e-9 * longest * longest",
     "    double_area <= 1e-30 * longest * longest"),
    ("flat band far too wide", FLIP,
     "    double_area <= 1e-9 * longest * longest",
     "    double_area <= 1e-1 * longest * longest"),
    ("flat cap triangle kept", PROFILE,
     "            && crate::ring_triangulation::is_flat(a, b, c)\n        {",
     "            && crate::ring_triangulation::is_flat(a, b, c)\n            && false\n        {"),
    ("cap sliver length is zero", PROFILE,
     "    let flat = tolerance.linear();",
     "    let flat = 0.0;"),
    ("slivers are not flipped away", RINGS,
     "    flip::unsliver(points, &mut triangles, flat);",
     "    let _ = flip::unsliver;"),
    ("flip ignores convexity", FLIP,
     "if orient(at(u), at(x), at(w)) <= 0 || orient(at(x), at(v), at(w)) <= 0 {",
     "if false {"),
    ("flip may create a sliver", FLIP,
     "if is_sliver(at(u), at(x), at(w), flat) || is_sliver(at(x), at(v), at(w), flat) {",
     "if false {"),
    ("flip takes the shortest edge", FLIP,
     "if length(e) > length(best) { e } else { best }",
     "if length(e) < length(best) { e } else { best }"),
    ("flat refusal ignores the sliver test", PROFILE,
     "        if crate::ring_triangulation::is_sliver(a, b, c, flat)\n            && crate::ring_triangulation::is_flat(a, b, c)",
     "        if crate::ring_triangulation::is_flat(a, b, c)"),
    ("extrusion caps use the plain triangulation", EXTRUDE,
     "crate::profile::triangulate_cap(rings, tolerance)?",
     "crate::profile::triangulate(rings)?"),
]

probe_lib.main(MUTANTS, tests=TARGETS, timeout=900)
