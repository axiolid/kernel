"""Mutation probe for homology bases with boundary or non-orientable (#144).

`homology_generators` extends tree-cotree past the closed-orientable case
by (a) allowing it for closed non-orientable components and (b) giving
boundary edges a virtual dual node per boundary loop, then completing the
basis with `b - 1` of the traced boundary loops. Each mutant must corrupt
one of those moving parts and turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
T = "crates/algorithms/query/inspect/src/topology.rs"
TESTS = ["-p", "axiolid-inspect", "--test", "topology"]

MUTANTS = [
    # A boundary edge's dual side is always its own loop's node 0, so
    # separate boundary loops collapse into one virtual node.
    (
        'boundary loops share one virtual node',
        T,
        'let node = members.len() + loop_of[&(a, b)];',
        'let node = members.len();',
        TESTS,
    ),
    # Boundary edges get no dual node at all: every leftover interior-edge
    # index is looked up as if the edge had two triangle uses.
    (
        'boundary edge dual side ignored',
        T,
        '        let (s, t) = if let [x, y] = uses[..] {\n            (slot[&x.triangle], slot[&y.triangle])\n        } else {\n            let node = members.len() + loop_of[&(a, b)];\n            (slot[&uses[0].triangle], node)\n        };',
        '        let (s, t) = (slot[&uses[0].triangle], slot[&uses[0].triangle]);',
        TESTS,
    ),
    # Boundary loop tracing never marks the starting edge used, so a loop
    # of length > 1 spins forever or double-counts its first edge; force a
    # wrong (empty) loop set instead so every boundary case misses rank.
    (
        'boundary loops never traced',
        T,
        '    loops\n}',
        '    let _ = loops;\n    Vec::new()\n}',
        TESTS,
    ),
    # The loop-membership map keys the wrong pair (unwrapped, not
    # normalised to min/max), so lookups for edges stored the other way
    # round panic or miss.
    (
        'loop_of keyed unnormalised',
        T,
        'loop_of.insert((a.min(b), a.max(b)), index);',
        'loop_of.insert((a, b), index);',
        TESTS,
    ),
    # `close_through_tree` walks only from `a`, never lifting `b` to the
    # same depth first, so loops through an ancestor come out wrong.
    (
        'lca walk skips lifting b',
        T,
        '    while depth[&y] > depth[&x] {\n        y = parent[&y];\n        up_b.push(y);\n    }',
        '',
        TESTS,
    ),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=1200,
    ).returncode


survivors = []
for name, rel, old, new, target in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(target)
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
