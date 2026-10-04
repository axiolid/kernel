"""Mutation probe for the certified profile triangulation (#253).

Each mutant weakens one exact decision of validation, bridging, ear
clipping or the certificate, and must turn a test red.

Equivalent mutants, deliberately not listed:

- filing every node in the ear-search grid, not only those that do not
  turn strictly left: a node that turns left never blocks a valid ear
  (see `clip.rs`), so only the work changes.
- stepping back to the previous corner after a cut instead of skipping
  ahead: a different, equally valid ear order.
- reading orientation at the lexicographically largest vertex instead
  of the smallest: it is strictly convex too.
- testing the bridge's hole end locally, the edges ending at a bridge's
  ends, or the last triangle's turn: each is implied by the edge test or
  the certificate (see `bridge.rs`, `clip.rs`), so they were removed.
- dropping the certificate's per-edge use count: an interior edge used
  twice leaves an odd number of directed edges to pair, so the twin test
  fails on the same input.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/construction/construct/src/ring_triangulation"
M = D + ".rs"
V = D + "/validate.rs"
B = D + "/bridge.rs"
C = D + "/clip.rs"
TESTS = ["-p", "axiolid-construct", "--test", "profile_holes"]
LIB = ["-p", "axiolid-construct", "--lib", "ring_triangulation"]

MUTANTS = [
    # Ear clipping.
    ('ear test ignores a vertex on the triangle\'s boundary', C,
     'orient(a, p, q) >= 0 && orient(p, b, q) >= 0 && orient(b, a, q) >= 0',
     'orient(a, p, q) > 0 && orient(p, b, q) > 0 && orient(b, a, q) > 0', TESTS),
    ('a straight corner taken as an ear', C,
     '    if orient(a, p, b) <= 0 {\n        return false;',
     '    if orient(a, p, b) < 0 {\n        return false;', TESTS),
    ('bridge duplicates block their own ears', C,
     '        if q == a || q == p || q == b {',
     '        if q == p {', TESTS),
    ('straight nodes left out of the grid', C,
     '            if orient(prev, p, next) > 0 {',
     '            if orient(prev, p, next) >= 0 {', TESTS),
    ('grid row ignores y', C,
     '            let cell = grid.row(p.y) * columns + grid.column(p.x);',
     '            let cell = grid.column(p.x);', TESTS),
    # Bridging.
    ('bridge end not checked to enter the polygon at the boundary', B,
     'locally_inside(p_prev, p, p_next, m) && clear(m, p, &edges)',
     'clear(m, p, &edges)', TESTS),
    ('bridge not checked against edges', B,
     'locally_inside(p_prev, p, p_next, m) && clear(m, p, &edges)',
     'locally_inside(p_prev, p, p_next, m)', TESTS),
    ('waiting holes not checked against a bridge', B,
     '            if g == h || pending[g] {', '            if g == h {', TESTS),
    ('earlier bridges not checked against a bridge', B,
     '        let mut edges: Vec<(Point2, Point2)> = main\n            .iter()',
     '        let mut edges: Vec<(Point2, Point2)> = main\n            .iter()\n            .filter(|_| false)', TESTS),
    ('holes bridged from their leftmost vertex', B,
     '                if p.x > q.x || (p.x == q.x && p.y > q.y) {',
     '                if p.x < q.x || (p.x == q.x && p.y > q.y) {', TESTS),
    ('splice loses the two bridge nodes from the count', B,
     '        self.len += hole_len + 2;', '        self.len += hole_len;', TESTS),
    # Shared predicates.
    ('segments touching at an endpoint read as apart', M,
     '    (o1 == 0 && within(p1, p2, q1))\n        || (o2 == 0 && within(p1, p2, q2))\n        || ',
     '    ', TESTS),
    ('a point past a vertical segment read as on it', M,
     '    a.x.min(b.x) <= c.x && c.x <= a.x.max(b.x) && a.y.min(b.y) <= c.y && c.y <= a.y.max(b.y)',
     '    a.x.min(b.x) <= c.x && c.x <= a.x.max(b.x)', TESTS),
    ('a point past a horizontal segment read as on it', M,
     '    a.x.min(b.x) <= c.x && c.x <= a.x.max(b.x) && a.y.min(b.y) <= c.y && c.y <= a.y.max(b.y)',
     '    a.y.min(b.y) <= c.y && c.y <= a.y.max(b.y)', TESTS),
    ('reflex sector test reads convex', M,
     '        left_of_in || left_of_out\n', '        left_of_in && left_of_out\n', TESTS),
    # Validation.
    ('non-finite vertex accepted', V,
     '            .any(|p| !(p.x.is_finite() && p.y.is_finite()))',
     '            .any(|_| false)', TESTS),
    ('repeated vertex accepted', V,
     '            if points[ring.start + k] == points[ring.start + next] {',
     '            if false && points[ring.start + k] == points[ring.start + next] {', TESTS),
    ('fold back accepted', V,
     '            if folds_back(at(k), at(k + 1), at(k + 2)) {',
     '            if false && folds_back(at(k), at(k + 1), at(k + 2)) {', TESTS),
    ('self-intersection accepted', V,
     '            return refuse(format!("profile {} intersects itself", name(e.ring)));',
     '            return Ok(());', TESTS),
    ('touching rings accepted', V,
     '    if !segments_touch(e.a, e.b, f.a, f.b) {\n        return Ok(());',
     '    if true {\n        return Ok(());', TESTS),
    ('edge sweep stops at a shared x', V,
     '            if f.a.x.min(f.b.x) > right {', '            if f.a.x.min(f.b.x) >= right {', TESTS),
    ('edge sweep skips pairs sharing a y', V,
     '            if f.a.y.max(f.b.y) < low || f.a.y.min(f.b.y) > high {',
     '            if f.a.y.max(f.b.y) <= low || f.a.y.min(f.b.y) >= high {', TESTS),
    ('hole outside the outer ring accepted', V,
     '        if !strictly_inside(outer, probe) {', '        if false {', TESTS),
    ('nested hole accepted', V,
     '            if g != h && strictly_inside(', '            if false && strictly_inside(', TESTS),
    ('ring orientation ignored', V,
     '        if counter_clockwise != (r == 0) {', '        if r != 0 {', TESTS),
    ('ring orientation read at the first vertex', V,
     '        .unwrap_or(0);', '        .map(|_| 0)\n        .unwrap_or(0);', TESTS),
    ('point-in-ring ignores downward edges', V,
     '        if (b.y > a.y && side > 0) || (b.y < a.y && side < 0) {',
     '        if b.y > a.y && side > 0 {', TESTS),
    # Certificate (unit tests: the algorithm above never trips it).
    ('certificate skips the triangle count', M,
     '    if triangles.len() != expected {', '    if false {', LIB),
    ('certificate accepts a flat or clockwise triangle', M,
     '        if orient(corner(0), corner(1), corner(2)) <= 0 {',
     '        if orient(corner(0), corner(1), corner(2)) < -1 {', LIB),
    ('certificate skips ring edges', M,
     '            if directed.get(&edge) != Some(&1) || directed.contains_key(&(edge.1, edge.0)) {',
     '            if false {', LIB),
    ('certificate skips twins', M,
     '        if !boundary.contains_key(&(a, b)) && directed.get(&(b, a)) != Some(&1) {',
     '        if false {', LIB),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=2400,
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
