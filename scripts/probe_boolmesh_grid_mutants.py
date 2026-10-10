"""Mutation probe for the boolmesh winding broad phase (#203).

Each mutant brings back a way a query on the grid's far edge is dropped,
and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/providers/mesh/boolmesh/src/csg/manifold/collider.rs"
TESTS = [
    ["-p", "axiolid-mesh-boolean-boolmesh", "--lib"],
    ["-p", "axiolid-mesh-boolean-boolmesh", "--test", "overlapping_grid"],
]

MUTANTS = [
    ("far edge from the rounded cell product", C,
     "            if px < self.min.x || py < self.min.y || px > self.max.x || py > self.max.y {",
     "            let (hx, hy) = (self.min.x + self.cell * self.dim as Real, self.min.y + self.cell * self.dim as Real);\n            if px < self.min.x || py < self.min.y || px > hx || py > hy {"),
    ("the far edge itself excluded", C,
     "            if px < self.min.x || py < self.min.y || px > self.max.x || py > self.max.y {",
     "            if px < self.min.x || py < self.min.y || px >= self.max.x || py >= self.max.y {"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800)
