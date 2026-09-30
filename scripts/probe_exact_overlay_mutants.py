"""Mutation probe for the exact straight-edge overlay (#173): each mutant must be caught."""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/planar/overlay/src"

MUTANTS = [
    # Membership: fill rule, ring orientation, which side a kept piece faces.
    ("non-zero fill counts only positive windings", "exact_overlay.rs",
     "        FillRule::NonZero => winding != 0,",
     "        FillRule::NonZero => winding > 0,"),
    ("clockwise rings count positive", "exact_overlay.rs",
     "                weight: if positive { 1 } else { -1 },",
     "                weight: 1,"),
    ("kept pieces face away from the result", "exact_overlay.rs",
     "            (left != right).then_some(right)",
     "            (left != right).then_some(left)"),
    ("turning vertices dropped as straight", "exact_arc/arrangement.rs",
     "        ask(true) == Sign::Zero && ask(false) == Sign::Positive",
     "        ask(false) == Sign::Positive"),
    # Chain reduction of a soup: cancel only edges given both ways, and
    # keep only simple cycles.
    ("edges given the same way cancel too", "exact_overlay.rs",
     "            let (edge, step) = if a < b { ((a, b), 1) } else { ((b, a), -1) };",
     "            let (edge, step) = if a < b { ((a, b), 1) } else { ((b, a), 1) };"),
    ("any cycle taken as simple", "exact_overlay.rs",
     "    if cycles.iter().all(|cycle| simple(cycle)) {",
     "    if true {"),
    ("neighbours folding back taken as simple", "exact_overlay.rs",
     "                    if (du.0 == dw.0 && u.x != v.x) || (du.1 == dw.1 && u.y != v.y) {",
     "                    if false {"),
    # One touch clause alone is equivalent in `simple`: a touching vertex
    # ends two edges, and one of them presents it to another clause.
    ("touches taken as apart", "exact_overlay.rs",
     "    (o1 == Sign::Zero && within(a, b, c))\n        || (o2 == Sign::Zero && within(a, b, d))\n        || (o3 == Sign::Zero && within(c, d, a))\n        || (o4 == Sign::Zero && within(c, d, b))",
     "    false"),
    # Output rounding: crossings to the nearest double, inputs untouched.
    ("rounding never corrects the guess upwards", "exact_arc/point.rs",
     "            Sign::Positive => r = up,",
     "            Sign::Positive => break,"),
    ("vertices rounded by the quick approximation", "exact_arc/arrangement.rs",
     "        self.vertices[index].rounded()",
     "        self.vertices[index].approx()"),
    ("exact-box identity trusts one box", "exact_arc/point.rs",
     "    if exact(a) && exact(b) {",
     "    if exact(a) || exact(b) {"),
    # f64 fast paths: each must stay a certified filter.
    ("apart test takes one side for both", "exact_arc/edge.rs",
     "    if one_side(sides[0], sides[1]) || one_side(sides[2], sides[3]) {",
     "    if one_side(sides[0], sides[1]) || one_side(sides[2], sides[2]) {"),
    # Not listed: the shared-end shortcut in `segments_quick`. Wrong there,
    # it drops the far end of a collinear overlap, but that end is a vertex
    # whose next edge meets the same line and restores the split, so no
    # result can tell (an equivalent mutant; checked, it survives).
    ("reversed twin read as same way", "exact_arc/edge.rs",
     "        } else if a == d && b == c {\n            Some(false)",
     "        } else if a == d && b == c {\n            Some(true)"),
    ("part right of the sample counted backwards", "exact_arc.rs",
     "            return if up { 1 } else { -1 };",
     "            return if up { -1 } else { 1 };"),
    # Broad phase: a query that misses items skips a real answer.
    ("box tree drops a subtree", "exact_arc/boxes.rs",
     "                    stack.push(node.end);\n",
     ""),
    ("long rings wound leftwards", "exact_arc/arrangement.rs",
     "                tree.query(|b| b.may_hold((sx.0, f64::INFINITY), sy))",
     "                tree.query(|b| b.may_hold((f64::NEG_INFINITY, sx.1), sy))"),
]

def run():
    return subprocess.run(
        ["cargo", "test", "-q", "--release", "-p", "axiolid-overlay",
         "--lib", "--test", "exact_straight", "--test", "region",
         "--test", "roundtrip", "--test", "collinear", "--test", "contract",
         "--test", "arc_overlay", "--test", "arrangement"],
        cwd=ROOT, capture_output=True, text=True, timeout=900,
    ).returncode

survivors = []
for name, rel, old, new in MUTANTS:
    path = SRC / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run()
        except subprocess.TimeoutExpired:
            code = -1  # a hang is a detection
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
