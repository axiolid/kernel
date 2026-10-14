"""Mutation probe for the refusal of a repeated vertex inside the exact overlay (#309): each mutant must be caught.

Every public entry validates its rings first, so only the crate's own
tests reach these guards: the repeated-vertex check at the top of the
internal straight-edge boolean, and the refusal of a piece of no length
where its sample is taken. Without either, the tests' calls must hang
(caught by their ten-second limit) or answer with something else.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/planar/overlay/src"
TESTS = [["-p", "axiolid-overlay", "--lib", "refused_by_name"]]

MUTANTS = [
    ("no repeated-vertex check before the boolean", "exact_overlay.rs",
     "        .any(repeats)\n",
     "        .any(|ring| repeats(ring) && false)\n"),
    ("a repeat across the closing edge missed", "exact_overlay.rs",
     "        (0..n).any(|i| ring.points[i] == ring.points[(i + 1) % n])",
     "        (1..n).any(|i| ring.points[i - 1] == ring.points[i])"),
    ("a piece of no length sampled", "exact_arc.rs",
     "    if same_point(from, to) {\n        return Err(OverlayError::RepeatedVertex);\n    }\n",
     ""),
]

probe_lib.main(MUTANTS, tests=TESTS, root=SRC, timeout=300)
