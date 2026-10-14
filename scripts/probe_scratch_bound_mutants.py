"""Mutation probe for affine scratch bounds and their worker term (#226).

Each mutant brings back a way a declared scratch bound stops being an upper
bound -- the fixed term dropped, the per-worker term dropped, the worker
count ignored somewhere between the declaration and the budget check, or the
old purely per-triangle declaration -- and must turn a test red.

`scratch_bound` runs with `parallel` (the per-worker term is measured on
pools of 1 to 64 workers) and without rayon (the base term alone carries
small inputs).

The #259 mutants take back what the batch paths charge a budget step by
step (every concurrent pair of a `parallel-batch` level, the intermediates
and copies a batch holds, its bookkeeping), the budget check of a direct
boolean in a build without `parallel`, and the pool of one that a one-thread
`CpuExecution` runs on. They name their own targets: `batch_budget`, which
pins the charges, and `scratch_bound`, whose small-budget runs catch a copy
made before the check that should have refused it. Every mutant must
compile: a mutant killed by a compile error proves nothing.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
X = "crates/contracts/common/base/src/execution.rs"
D = "crates/execution/dispatch/src/boolean.rs"
P = "crates/providers/mesh/boolmesh/src/provider.rs"
B = "crates/providers/mesh/boolmesh/src/budget.rs"
C = "crates/execution/cpu/src/execution.rs"
BOOLMESH = "axiolid-mesh-boolean-boolmesh"
BATCH_BUDGET = ["-p", BOOLMESH, "--test", "batch_budget"]
BATCH_BUDGET_PARALLEL = ["-p", BOOLMESH, "--features", "parallel-batch", "--test", "batch_budget"]
SCRATCH_BOUND = ["-p", BOOLMESH, "--test", "scratch_bound"]
CPU = ["-p", "axiolid-backend-cpu", "--features", "parallel", "--lib"]
DISPATCH_BUDGET = ["-p", "axiolid-dispatch", "--all-features", "--test", "mesh_boolean_budget"]
# The measurement first, so a mutant it kills is reported as killed by it.
TESTS = [
    ["-p", BOOLMESH, "--features", "parallel", "--test", "scratch_bound"],
    ["-p", BOOLMESH, "--test", "scratch_bound"],
    ["-p", BOOLMESH, "--features", "parallel", "--test", "registry"],
    ["-p", "axiolid-dispatch", "--all-features", "--test", "mesh_boolean_budget"],
    ["-p", "axiolid-contracts", "--test", "execution_contracts"],
]

MUTANTS = [
    ("contract: base term dropped", X,
     "                base_bytes.checked_add(variable)",
     "                let _ = base_bytes;\n                Some(variable)"),
    ("contract: per-worker term dropped", X,
     "let Some(per_workers) = bytes_per_worker.checked_mul(workers) else {",
     "let Some(per_workers) = bytes_per_worker.checked_mul(0) else {"),
    ("contract: worker count ignored (one worker)", X,
     "let Some(per_workers) = bytes_per_worker.checked_mul(workers) else {",
     "let Some(per_workers) = bytes_per_worker.checked_mul(1) else {"),
    ("contract: fits_budget ignores the requested width", X,
     "options.parallelism().worker_bound()",
     "1"),
    ("dispatch: configured pool width ignored", D,
     "return scratch.fits_budget_on(options, elements, workers);",
     "return scratch.fits_budget_on(options, elements, 1);"),
    ("dispatch: boolean budgeted for the subject alone", D,
     "let elements = subject.triangle_count() + tool.triangle_count();",
     "let elements = subject.triangle_count();"),
    ("boolmesh: base term dropped", P,
     "            base_bytes: SCRATCH_BASE_BYTES,",
     "            base_bytes: 0,"),
    ("boolmesh: per-worker term dropped", P,
     "const WORKER_BYTES: usize = 16 * 1024;",
     "const WORKER_BYTES: usize = 0;"),
    ("boolmesh: running pool width ignored", B,
     "        match concurrent_bound(1, elements, running_workers()) {",
     "        match concurrent_bound(1, elements, 1) {"),
    ("boolmesh: old flat 4096 B/triangle declaration", P,
     "        ScratchRequirement::Affine {\n"
     "            base_bytes: SCRATCH_BASE_BYTES,\n"
     "            bytes_per_element: SCRATCH_BYTES_PER_TRIANGLE,\n"
     "            bytes_per_worker: SCRATCH_BYTES_PER_WORKER,\n"
     "        }",
     "        let _ = (SCRATCH_BASE_BYTES, SCRATCH_BYTES_PER_TRIANGLE, SCRATCH_BYTES_PER_WORKER);\n"
     "        ScratchRequirement::PerElement {\n"
     "            bytes_per_element: 4096,\n"
     "        }"),
    # #259: batch paths and the one-thread pool.
    ("boolmesh: a parallel-batch level is charged one boolean", P,
     "charge_step(options, held, chunks.len(), elements)?;",
     "charge_step(options, held, 1, elements)?;",
     BATCH_BUDGET_PARALLEL),
    ("boolmesh: a parallel-batch level is charged one boolean per worker", P,
     "charge_step(options, held, chunks.len(), elements)?;",
     "charge_step(options, held, chunks.len().min(crate::budget::running_workers()), elements)?;",
     BATCH_BUDGET_PARALLEL),
    ("boolmesh: a parallel-batch level is not charged", P,
     "                crate::budget::charge_step(options, held, chunks.len(), elements)?;",
     "                let _: (usize, usize) = (held, elements);",
     BATCH_BUDGET_PARALLEL),
    # Both places: since #308 the union's own budget leaves out what the
    # level holds too, so dropping it from the step's charge alone changes
    # nothing (the union is refused at its own admission instead).
    ("boolmesh: a sequential union is not charged the level it holds", P,
     [("charge_step(options, held + done, 1, triangles(pair))?;",
       "charge_step(options, done, 1, triangles(pair))?;"),
      ("let budget = crate::budget::remaining(options, held + done);",
       "let budget = crate::budget::remaining(options, done);")],
     BATCH_BUDGET),
    ("boolmesh: subtract_many is not charged its fused tool", P,
     "                _ => crate::budget::fused_bound(&members, &subject.attributes),",
     "                _ => 0,",
     BATCH_BUDGET),
    ("boolmesh: subtract_many is not charged the subject it has cut", P,
     "                Cow::Owned(mesh) => crate::budget::heap_bytes(mesh),\n"
     "            } + crate::budget::fates_bytes(&fates)",
     "                Cow::Owned(_) => 0,\n"
     "            } + crate::budget::fates_bytes(&fates)",
     BATCH_BUDGET),
    ("boolmesh: union_many copies before charging them", P,
     "            copies + solids.len() * crate::budget::fates_bytes(&seed) + bookkeeping,",
     "            solids.len() * crate::budget::fates_bytes(&seed) + bookkeeping + copies * 0,",
     SCRATCH_BOUND),
    ("boolmesh: union_many bookkeeping is not charged", P,
     "        let bookkeeping = solids.len() * UNION_BOOKKEEPING_PER_SOLID;",
     "        let bookkeeping = solids.len() * UNION_BOOKKEEPING_PER_SOLID * 0;",
     SCRATCH_BOUND),
    ("boolmesh: subtract_many bookkeeping is not charged", P,
     "        let bookkeeping = tools.len() * SUBTRACT_BOOKKEEPING_PER_TOOL;",
     "        let bookkeeping = tools.len() * SUBTRACT_BOOKKEEPING_PER_TOOL * 0;",
     SCRATCH_BOUND),
    ("boolmesh: a direct boolean is budgeted only with parallel", P,
     "        let result_budget = crate::budget::ResultBudget::admit(budget, elements)?;",
     "        #[cfg(not(feature = \"parallel\"))]\n"
     "        let result_budget = None;\n"
     "        #[cfg(feature = \"parallel\")]\n"
     "        let result_budget = crate::budget::ResultBudget::admit(budget, elements)?;",
     BATCH_BUDGET),
    ("boolmesh: no boolean still pays the per-worker term", B,
     "    if booleans == 0 {\n        return Some(0);\n    }\n",
     "",
     BATCH_BUDGET_PARALLEL),
    ("cpu: a one-thread context runs on the caller's pool", C,
     "        self.pool\n            .get_or_init(",
     "        if self.threads.get() == 1 {\n            return operation();\n        }\n"
     "        self.pool\n            .get_or_init(",
     [CPU, DISPATCH_BUDGET]),
    ("dispatch: a one-thread context is charged the options' width", D,
     "            return scratch.fits_budget_on(options, elements, workers);",
     "            if workers > 1 {\n"
     "                return scratch.fits_budget_on(options, elements, workers);\n"
     "            }",
     DISPATCH_BUDGET),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
