"""Mutation probe for the Bentley-Ottmann segment sweep (#146).

Each mutant drops a neighbour check, weakens an exact decision, breaks the
sweep order, the status tree or the degeneracy handling, or mislabels the
result, and must turn a test red.

Equivalent mutants, deliberately not listed: rotating the lower-priority
child up in `Status::remove`, or skipping the priority rotations in
`Status::insert_before`, only unbalance the treap (a speed change, not a
result change); skipping the `q > *p` guard in `Sweep::check` is equivalent
because a crossing of two neighbours never lies behind the sweep; and
accepting a touch (one zero orientation) as a crossing in `proper_crossing`
is equivalent because the touching point is an endpoint, already queued, so
the new key deduplicates into it and is not counted.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/planar/overlay/src/segment_sweep.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "segment_sweep"]

MUTANTS = [
    ('lower neighbour never checked', S, '                self.check(below, low, p);\n', '', TESTS),
    ('upper neighbour never checked', S, '                self.check(high, above, p);\n', '', TESTS),
    ('neighbours closing a gap never checked', S, '            _ => self.check(below, above, p),', '            _ => {}', TESTS),
    ('crossings scheduled behind the sweep only', S, '        if q > *p {', '        if q < *p {', TESTS),
    ('one straddle test skipped', S, '    o3 != 0 && o4 != 0 && o3 != o4', '    true', TESTS),
    ('directions ordered downwards', S, '                1 => Ordering::Less,\n                -1 => Ordering::Greater,', '                1 => Ordering::Greater,\n                -1 => Ordering::Less,', TESTS),
    ('only the first segment through a point found', S, '        while above != NIL && side(above) == 0 {', '        if above != NIL && side(above) == 0 {', TESTS),
    ('segments never end', S, '            .map(|&id| at.is_some_and(|q| segs[id as usize].hi == q))', '            .map(|_| false)', TESTS),
    ('zero-length segments enter the status', S, '            .filter(|&id| !self.segs[id as usize].degenerate)', '            .filter(|_| true)', TESTS),
    ('undecided side of a crossing read as on', S, '            exact_sign(&value)\n', '            0\n', TESTS),
    ('undecided coordinate order read as equal', S, '    exact_sign(&pn.mul(&qw).sub(&qn.mul(&pw))).cmp(&0)', '    Ordering::Equal', TESTS),
    ('input points ordered by x alone', S, '            return cmp_f64(a.x, b.x).then_with(|| cmp_f64(a.y, b.y));', '            return cmp_f64(a.x, b.x);', TESTS),
    ('ends never put in sweep order', S, '            let reversed = cmp_f64(a.x, b.x).then_with(|| cmp_f64(a.y, b.y)) == Ordering::Greater;', '            let reversed = false;', TESTS),
    ('crossings left at the rounding guess', S, '                    round_ratio(&r.x, &r.w, guess(&r.x)),', '                    guess(&r.x),', TESTS),
    ('insertion drops a subtree', S, '            let last = self.max_from(self.node(before).left);', '            let last = self.node(before).left;', TESTS),
    ('lone segments through a point reported', S, '        if starts.len() + block.len() < 2 {', '        if starts.len() + block.len() < 1 {', TESTS),
    ('overlaps not closed at an event', S, '                    piece.to = Some(p.clone());\n                    self.closed_here.push(index);', '                    self.closed_here.push(index);', TESTS),
    ('overlap pieces never merged', S, '                    .find(|&index| self.pieces[index].segments == members);', '                    .find(|_| false);', TESTS),
    ('reversed starts labelled as starts', S, '        } else if self.reversed {\n            SegmentLocation::End\n        } else {\n            SegmentLocation::Start\n        }', '        } else {\n            SegmentLocation::Start\n        }', TESTS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
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
