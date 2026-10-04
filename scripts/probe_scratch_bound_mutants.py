"""Mutation probe for affine scratch bounds and their worker term (#226).

Each mutant brings back a way a declared scratch bound stops being an upper
bound -- the fixed term dropped, the per-worker term dropped, the worker
count ignored somewhere between the declaration and the budget check, or the
old purely per-triangle declaration -- and must turn a test red.

`scratch_bound` runs with `parallel` (the per-worker term is measured on
pools of 1 to 64 workers) and without rayon (the base term alone carries
small inputs).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
X = "crates/contracts/common/base/src/execution.rs"
D = "crates/execution/dispatch/src/boolean.rs"
P = "crates/providers/mesh/boolmesh/src/provider.rs"
BOOLMESH = "axiolid-mesh-boolean-boolmesh"
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
    ("boolmesh: running pool width ignored", P,
     "                rayon::current_num_threads(),",
     "                1,"),
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
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode


def killed_by(name):
    # Stop at the first red test: a mutant only has to be killed once.
    for target in TESTS:
        try:
            if run(target) != 0:
                return " ".join(target)
        except subprocess.TimeoutExpired:
            return "timeout: " + " ".join(target)
    return None


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        killer = killed_by(name)
    finally:
        path.write_text(original)
    if killer:
        print(f"killed   {name}  [{killer}]", flush=True)
    else:
        print(f"SURVIVED {name}", flush=True)
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
