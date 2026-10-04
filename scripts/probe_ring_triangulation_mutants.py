"""Mutation probe for the certified profile triangulation (#253, #262).

Each mutant weakens one exact decision of validation, bridging, ear
clipping, the pinch path (rings touching at single points under
`PinchPolicy::Accept`) or the certificate, and must turn a test red.

#253: 31/31. #262 added 21 pinch mutants; first run 51/52: "spokes not
checked to alternate" survived because the fixture for holes crossing at
shared corners was refused earlier, as a hole inside another. The square
now starts outside the other hole, so only the spokes catch it: 52/52.

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
- skipping the bounding-box test before a pinched loop's containment
  probe: only the work changes.
- reading a cycle's outer or hole side at its lexicographically largest
  vertex: the part's wedges are convex there too, and the outside of a
  hole cycle wraps round it as well.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/construction/construct/src/ring_triangulation"
M = D + ".rs"
V = D + "/validate.rs"
B = D + "/bridge.rs"
C = D + "/clip.rs"
P = D + "/pinch.rs"
TESTS = ["-p", "axiolid-construct", "--test", "profile_holes"]
PINCHES = ["-p", "axiolid-construct", "--test", "profile_pinches"]
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
    # Touching rings under PinchPolicy::Accept (#262).
    ('a vertex inside another edge refused', V,
     '        (Some(&(_, edge, vertex)), None) => {',
     '        (Some(&(_, edge, vertex)), None) if false => {', PINCHES),
    ('edges from a shared end along one line accepted', V,
     '                if far_e == far_f || folds_back(far_e, p, far_f) {',
     '                if far_e == far_f {', PINCHES),
    ('edges lying inside each other read as a touch', V,
     '        (Some(_), Some(_)) => overlap(e, f),', '        (Some(_), Some(_)) => Ok(()),', PINCHES),
    ('pinch path never taken', P,
     '    (shared || !touches.is_empty()).then_some(canon)',
     '    (false && (shared || !touches.is_empty())).then_some(canon)', PINCHES),
    ('signed zero read as another point', P,
     '            let key = ((p.x + 0.0).to_bits(), (p.y + 0.0).to_bits());',
     '            let key = (p.x.to_bits(), p.y.to_bits());', PINCHES),
    ('touching vertex not inserted into its edge', P,
     '                walk.extend_from_slice(list);', '                let _ = &list;', PINCHES),
    ('touching vertices kept in one order whichever way the edge runs', P,
     '                    if (pa.x != pb.x && pb.x < pa.x) || (pa.x == pb.x && pb.y < pa.y) {',
     '                    if false {', PINCHES),
    ('ring not split at a repeated vertex', P,
     '            if let Some(&i) = position.get(&v) {',
     '            if let Some(&i) = position.get(&v).filter(|_| false) {', PINCHES),
    ('loops all turned counter-clockwise', P,
     '        if ring_turns_left(&ring) != (d % 2 == 0) {',
     '        if !ring_turns_left(&ring) {', PINCHES),
    ('loops kept as given', P,
     '        if ring_turns_left(&ring) != (d % 2 == 0) {', '        if false {', PINCHES),
    ('loop sharing every vertex read on the wrong side', P,
     '        locally_inside(at(before), at(v), at(after), at(w))',
     '        !locally_inside(at(before), at(v), at(after), at(w))', PINCHES),
    ('containing loop\'s winding ignored at a shared vertex', P,
     '        if !counter_clockwise[i] {', '        if false {', PINCHES),
    ('pinched hole outside the outer ring accepted', P,
     '        if !inside[j].iter().any(|&i| pieces[i].ring == 0) {', '        if false {', PINCHES),
    ('pinched hole inside another hole accepted', P,
     '            .find(|&&i| pieces[i].ring != 0 && pieces[i].ring != piece.ring)',
     '            .find(|&&i| false && pieces[i].ring != 0 && pieces[i].ring != piece.ring)', PINCHES),
    ('spokes not checked to alternate', P,
     '    if tie || !alternates {', '    if tie {', PINCHES),
    ('spoke order ignores the half-plane', P,
     '            .cmp(&half(b.0))\n',
     '            .cmp(&half(a.0))\n', PINCHES),
    ('wedge closed by the previous arriving edge', P,
     '        let (leaving, arriving) = (spokes[i], spokes[(i + 1) % count]);',
     '        let (leaving, arriving) = (spokes[i], spokes[(i + count - 1) % count]);', PINCHES),
    ('cycle read as outer at any left turn', P,
     '                .all(|&n| {', '                .any(|&n| {', PINCHES),
    ('hole cycle given to the first outer cycle around it', P,
     '            .find(|&o| candidates.iter().all(|&d| d == o || !around(o, d)));',
     '            .find(|_| true);', PINCHES),
    ('pinched hole bridged from its first visit', B,
     '            if visits.len() == 1 {', '            if true {', PINCHES),
    ('certificate counts every part as one', M,
     '    let expected = (vertex_count + 2 * holes).saturating_sub(2 * outers);',
     '    let expected = (vertex_count + 2 * holes).saturating_sub(2);', PINCHES),
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
