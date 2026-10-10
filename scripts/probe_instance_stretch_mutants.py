"""Mutation probe for instance stretch and the shrunk local budget (#258).

An instance's source is meshed to the world chord budget shrunk by the
transform's stretch bound, down to the largest value whose rounded product
with it fits, and the source's bound is grown back by that same bound. Each
mutant meshes the source too coarsely, grows the bound too little, or lets
the two factors disagree, and must turn a test red: a rigid or scaled
instance reporting over its budget, or a dense sample of the exact surface
further from the mesh than the reported bound.

Each mutant runs in its own process group under a timeout that kills the
whole group; the mutated file is always restored. A missing or ambiguous
anchor aborts the probe.
"""
import os, pathlib, signal, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/execution/compile/src/compiler.rs"
D = "crates/execution/compile/src/deviation.rs"
B = "crates/execution/compile/src/bounded.rs"
TESTS = ["-p", "axiolid-mesh-compile", "--all-features", "--lib", "--test", "instance_stretch"]
TIMEOUT = 1200

MUTANTS = [
    ('local budget the plain quotient, not shrunk to fit', D,
     '    while shrunk > 0.0 && shrunk * stretch > budget {',
     '    while false && shrunk * stretch > budget {', TESTS),
    ('instance source meshed to the world budget', C,
     '        chord: crate::deviation::shrunk_budget(budget.chord, stretch),',
     '        chord: budget.chord,', TESTS),
    ('instance source budget divided by the column length only', C,
     '    let stretch = crate::deviation::stretch(transform);',
     '    let stretch = max_scale;', TESTS),
    ('instance bound not grown by the stretch', D,
     '        let factor = stretch(transform);',
     '        let factor = 1.0;', TESTS),
    ('stretch margin for near-orthogonal columns dropped', D,
     '        sx.max(sy).max(sz) * (1.0 + 3.0 * eps)',
     '        sx.max(sy).max(sz)', TESTS),
    ('rigid transform read as stretch exactly one', D,
     '        sx.max(sy).max(sz) * (1.0 + 3.0 * eps)',
     '        if (sx.max(sy).max(sz) - 1.0).abs() <= eps { 1.0 } else { sx.max(sy).max(sz) * (1.0 + 3.0 * eps) }',
     TESTS),
    ('sheared columns read as the largest column', D,
     '        (sx * sx + sy * sy + sz * sz).sqrt()',
     '        sx.max(sy).max(sz)', TESTS),
    ('plane boundaries not shrunk by the frame stretch', B,
     '        .with_chord_error(crate::deviation::shrunk_budget(chord, stretch))',
     '        .with_chord_error(chord)', TESTS),
]


def run(target):
    """The test run's exit code; -1 when it timed out (its group killed)."""
    proc = subprocess.Popen(
        ["cargo", "test", "-q", *target],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    try:
        return proc.wait(timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        return -1
    except BaseException:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        raise


def main():
    for name, rel, old, _, _ in MUTANTS:
        count = (ROOT / rel).read_text().count(old)
        if count != 1:
            sys.exit(f"anchor for '{name}' found {count} times in {rel}")
    survivors = []
    for name, rel, old, new, target in MUTANTS:
        path = ROOT / rel
        original = path.read_text()
        if original.count(old) != 1:
            sys.exit(f"anchor for '{name}' not unique/found in {rel}")
        try:
            path.write_text(original.replace(old, new))
            code = run(target)
        finally:
            path.write_text(original)
        status = "killed" if code != 0 else "SURVIVED"
        suffix = " (timeout)" if code == -1 else ""
        print(f"{status:8} {name}{suffix}", flush=True)
        if code == 0:
            survivors.append(name)
    print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
    sys.exit(1 if survivors else 0)


if __name__ == "__main__":
    main()
