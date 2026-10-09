"""Mutation probe for exact ring orientation and area in the overlay (#274).

A ring's orientation and its ZeroArea check are exact signs of its area
taken from its first vertex; areas reported as values are summed from the
first vertex too. Each mutant puts back a rounded `f64` sign, or the
shoelace over the coordinates as given, at one site, and must turn the
tests red: `tests/far_from_origin.rs` holds rings translated to 1e5, 5e6
and 1e7 m and a triangle whose rounded fan cancels to zero exactly.

Equivalent mutant, deliberately not listed: the exact sign taken over the
coordinates as given instead of relative to the first vertex. Exact
arithmetic is translation-invariant, so only how often the interval filter
falls back to the exact tier changes, never a sign.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/planar/overlay/src/"
L, S, E, A = D + "lib.rs", D + "settle.rs", D + "exact_overlay.rs", D + "arc.rs"
X, O = D + "exact_arc.rs", D + "offset.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "far_from_origin", "--test", "exact_straight"]

ROUNDED_SIGN = """{
    let a = signed(&Ring { points: points.to_vec() });
    if a > 0.0 { Sign::Positive } else if a < 0.0 { Sign::Negative } else { Sign::Zero }
}"""
RAW = """    r.points
        .iter()
        .zip(r.points.iter().cycle().skip(1))
        .take(r.points.len())
        .map(|(a, b)| a.x * b.y - b.x * a.y)
        .sum::<f64>()
        * 0.5"""

MUTANTS = [
    ('orientation from the rounded fan', L,
     'pub(crate) fn orientation(points: &[Point2]) -> Sign {\n    area_sign(points, 0.0)\n}',
     'pub(crate) fn orientation(points: &[Point2]) -> Sign ' + ROUNDED_SIGN, TESTS),
    ('ZeroArea from the rounded fan', L,
     '    let bound = 2.0 * linear * linear;\n    area_sign(points, bound) != Sign::Positive && area_sign(points, -bound) != Sign::Negative',
     '    signed(&Ring { points: points.to_vec() }).abs() <= linear * linear', TESTS),
    ('ZeroArea from the raw shoelace', L,
     '    let bound = 2.0 * linear * linear;\n    area_sign(points, bound) != Sign::Positive && area_sign(points, -bound) != Sign::Negative',
     '    let r = Ring { points: points.to_vec() };\n    (' + RAW.strip() + ').abs() <= linear * linear', TESTS),
    ('area measure from the raw shoelace', L,
     """    r.points
        .windows(2)
        .skip(1)
        .map(|w| (w[0] - o).perp_dot(w[1] - o))
        .sum::<f64>()
        * 0.5""",
     '    let _ = o;\n' + RAW, TESTS),
    ('a ring on one line refused as ZeroArea', L,
     'if zero_area(&r.points, t.linear()) && !exact_overlay::on_one_line(&r.points) {',
     'if zero_area(&r.points, t.linear()) {', TESTS),
    ('canonical orientation from the rounded fan', L,
     '    if (orientation(&r.points) == Sign::Positive) != want_positive {',
     '    if (signed(&r) > 0.0) != want_positive {', TESTS),
    ('boolean winding from the rounded fan', E,
     '            let positive = orientation(&points) == Sign::Positive;',
     '            let positive = crate::signed(&Ring { points: points.clone() }) > 0.0;', TESTS),
    ('settled pieces wound by the rounded fan', S,
     '                let same = (orientation(&piece) == Sign::Positive) == sense;',
     '                let same = (area(&piece) > 0.0) == sense;', TESTS),
    ('linked rings told outer from hole by the rounded area', X,
     '        if arc_orientation(&ring) == Sign::Positive {',
     '        if arc_ring_area(&ring) > 0.0 {', TESTS),
    ('straight arc rings oriented by the rounded area', A,
     '    if ring.is_polygonal() {\n        return crate::orientation(&points(ring));\n    }',
     '', TESTS),
    ('straight arc rings ZeroArea by the rounded area', A,
     '    let zero = if ring.is_polygonal() {', '    let zero = if false {', TESTS),
    ('arc ring area from the raw shoelace', A,
     """    for index in 1..count - 1 {
        let from = ring.vertices[index].point - origin;
        let to = ring.vertices[index + 1].point - origin;
        area += from.perp_dot(to);
    }""",
     """    let _ = origin;
    for index in 0..count {
        let from = ring.vertices[index].point;
        let to = ring.vertices[(index + 1) % count].point;
        area += from.x * to.y - to.x * from.y;
    }""", TESTS),
    ('ring_area from the raw shoelace', O,
     '    crate::signed(ring).abs()',
     '    let r = ring;\n    (' + RAW.strip() + ').abs()', TESTS),
]


def run(target):
    """The test exit code; a mutant that does not compile kills nothing."""
    done = subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    )
    assert "error[E" not in done.stderr, done.stderr[-2000:]
    return done.returncode


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
