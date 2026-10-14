"""Mutation probe for the result charge of boolmesh booleans (#308).

The declared scratch bound counts input triangles, so a boolean whose
operands cross many times admits on it and then outgrows it. The provider
charges the result inside the solve: `boolean45` counts the result's
half-edges before allocating it, the intersection pass stops once its
crossings alone could not fit, and every caller that holds memory around a
boolean (a composed symmetric difference, a batch step, a level of
concurrent pairs) hands it only what it does not hold.

Each mutant takes one of those back and must turn a test red:
`batch_budget` pins the charges to the byte, `scratch_bound` measures that
crossing combs, slab soups and pierced slabs stay within the budget they
admit or are refused within it, and the library tests pin the arithmetic.
Every mutant must compile: a mutant killed by a compile error proves
nothing.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
BOOLMESH = "axiolid-mesh-boolean-boolmesh"
SRC = "crates/providers/mesh/boolmesh/src/"
B = SRC + "budget.rs"
P = SRC + "provider.rs"
CSG = SRC + "csg.rs"
B45 = SRC + "csg/boolean45.rs"
K12 = SRC + "csg/boolean03/kernel12.rs"
LIB = ["-p", BOOLMESH, "--lib"]
LIB_BATCH = ["-p", BOOLMESH, "--features", "parallel-batch", "--lib"]
BATCH_BUDGET = ["-p", BOOLMESH, "--test", "batch_budget"]
SCRATCH_BOUND = ["--release", "-p", BOOLMESH, "--test", "scratch_bound"]

MUTANTS = [
    ("the result is never charged", B45,
     "    if budget.is_some_and(|budget| !budget.fits(nh)) {",
     "    if budget.is_some_and(|budget| !budget.fits(nh)) && nh == usize::MAX {",
     [BATCH_BUDGET, SCRATCH_BOUND]),
    ("a third of the slope per excess half-edge", B,
     "pub(crate) const SCRATCH_BYTES_PER_EXCESS_HALF_EDGE: usize = 448;",
     "pub(crate) const SCRATCH_BYTES_PER_EXCESS_HALF_EDGE: usize = 150;",
     [LIB, SCRATCH_BOUND]),
    ("the declaration covers twice the half-edges", B,
     "pub(crate) const RESULT_HALF_EDGES_PER_TRIANGLE: usize = 3;",
     "pub(crate) const RESULT_HALF_EDGES_PER_TRIANGLE: usize = 6;",
     [LIB, SCRATCH_BOUND]),
    ("the result's charge ignores what the declaration covers", B,
     "            .saturating_sub(self.covered)",
     "            .saturating_sub(0)",
     [LIB, BATCH_BUDGET]),
    ("the intersection pass is not capped", CSG,
     "    let cap = budget.map(ResultBudget::crossing_limit);",
     "    let cap: Option<usize> = None;",
     [SCRATCH_BOUND]),
    ("the crossing cap counts one half-edge per crossing", B,
     "        half_edges / 3\n",
     "        half_edges\n",
     [LIB, SCRATCH_BOUND]),
    ("a symmetric difference is not charged the union it holds", P,
     "        let budget = less(budget, held(&union));",
     "        let budget = less(budget, 0);",
     [BATCH_BUDGET]),
    ("a symmetric difference is not charged the intersection it holds", P,
     "        let budget = less(budget, held(&intersection));",
     "        let budget = less(budget, 0);",
     [BATCH_BUDGET]),
    ("a symmetric difference composed with the whole budget per step", P,
     "            return self.symmetric_difference(subject, tool, options, budget);",
     "            let _ = budget;\n"
     "            return axiolid_mesh_boolean_contract::symmetric_difference_via_composition(\n"
     "                self, subject, tool, options,\n"
     "            );",
     [BATCH_BUDGET]),
    ("a subtract_many step gets the whole budget", P,
     "            let budget = crate::budget::remaining(options, held);",
     "            let budget = crate::budget::remaining(options, 0);",
     [BATCH_BUDGET]),
    ("a sequential union_many step gets the whole budget", P,
     "                    let budget = crate::budget::remaining(options, held + done);",
     "                    let budget = crate::budget::remaining(options, done);",
     [BATCH_BUDGET]),
    ("concurrent pairs each get the whole budget", B,
     "    Some(own.saturating_add(surplus / booleans.max(1)))",
     "    Some(own.saturating_add(surplus))",
     [LIB_BATCH]),
    ("a parallel-batch pair gets the whole budget", P,
     "                        step(pair, share(pair)).map(|(mesh, fates)| (Cow::Owned(mesh), fates))",
     "                        let _ = share;\n"
     "                        step(pair, options.memory_budget_bytes())\n"
     "                            .map(|(mesh, fates)| (Cow::Owned(mesh), fates))",
     [["-p", BOOLMESH, "--features", "parallel-batch", "--test", "batch_budget"]]),
]

probe_lib.main(MUTANTS, timeout=3600, reject_uncompiled=True)
