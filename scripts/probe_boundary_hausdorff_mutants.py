"""Mutation probe for the certified Hausdorff distance between exact
boundaries (#224), and its faces matched as translates trimmed in another
chart (#227; `exact_hausdorff/translate.rs`, unit-tested there, and the
prisms of `boolean_arc_prisms_exact` in the construct tests).

Each mutant makes a bound unsound -- a matched bound below the gap it
claims to hold, a face matched across a different domain or chart, a
Lipschitz or point bound without its radius -- or loses the witnesses or
the two-sided combination, and must turn a test red. The unit tests in
`exact_hausdorff.rs` sample every family's matched bound;
`construct/tests/boundary_hausdorff.rs` holds the closed-form distances.

Equivalent mutants, deliberately not listed: the point query's lower bound
ignoring the queue (a popped bound is the least queued one), the first-come
tie-break reversed (order only), and a face patch that cannot shrink
dropped instead of kept (no test fixture reaches the rounding floor).

Not killed, and so not listed: the Lipschitz bound about the centre of a
patch with no witness (one straddling its face's boundary) losing its
radius. That bound is unsound, but in every fixture tried, including a
farthest point isolated at a vertex, the sound bounds of neighbouring
patches keep the interval above the true value until it closes. Nor the
re-charted bound of a translate losing its trim residue `L delta`: the
gate keeps the residue below 1e-9 relative, and no fixture's distance
depends on the sliver of domain it covers. Nor, since #227, edges never
split for witnesses: the translates whose farthest points sit on rims
(columns moved by 0.1 mm) now get them from the seeded support point,
and in the fixtures without a seed (turned, mirrored, resized pairs)
face patches next to the edge reach the farthest point within the
accuracy as fast.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
H = "crates/algorithms/query/measure/src/exact_hausdorff.rs"
D = "crates/algorithms/query/measure/src/exact_distance.rs"
T = "crates/algorithms/query/measure/src/exact_hausdorff/translate.rs"
TESTS = [
    ["-p", "axiolid-measure", "--all-features", "--lib", "exact_hausdorff"],
    ["-p", "axiolid-construct", "--test", "boundary_hausdorff"],
]

MUTANTS = [
    ("matched bound without its spread", H,
     "        spread += wobble[k] * delta.length();",
     "        spread += 0.0 * delta.length();"),
    ("matched bound read at the origins only", H,
     "        centre += delta * phi[k];",
     "        centre += delta * 0.0;"),
    ("cylinder angle oscillation ignored", H,
     "(vec![cu, su, m.y], vec![tu, tu, hv], vec![1.0, 1.0, reach_v])",
     "(vec![cu, su, m.y], vec![0.0, 0.0, hv], vec![1.0, 1.0, reach_v])"),
    ("cone mixed term without the angle", H,
     "            let mixed = hv + m.y.abs() * tu;",
     "            let mixed = hv;"),
    ("sphere oscillation along the meridian only", H,
     "            vec![tv + tu, tv + tu, tv],\n            vec![1.0, 1.0, 1.0],",
     "            vec![tu, tu, tv],\n            vec![1.0, 1.0, 1.0],"),
    ("torus tube oscillation round the axis only", H,
     "            vec![tu, tu, tv + tu, tv + tu, tv],",
     "            vec![tu, tu, tu, tu, tv],"),
    ("spline matched by its last control point", H,
     "                reach = reach.max((*pa - *pb).length());",
     "                reach = (*pa - *pb).length();"),
    ("faces matched on any trim", H,
     "                ia.is_some() && ia == ib && ca == cb",
     "                ia.is_some()"),
    ("cones of any apex matched", H,
     "            p.radius == q.radius && p.semi_angle == q.semi_angle",
     "            true"),
    ("splines of any knots matched", H,
     "                && p.u_knots == q.u_knots\n",
     "\n"),
    ("Lipschitz bound about a witness without its radius", H,
     "                    upper = upper.min(found.upper + element.radius);",
     "                    upper = upper.min(found.upper);"),
    ("two-sided upper from one side", H,
     "    distance.upper = forward.upper.max(backward.upper);",
     "    distance.upper = forward.upper;"),
    ("point bound without the element's radius", D,
     "    let mut best = (gap - element.radius - rounding).max(0.0);",
     "    let mut best = (gap - rounding).max(0.0);"),
    ("point bound projected without its own position", D,
     "        best = best.max(lo - x - pad).max(x - hi - pad);",
     "        best = best.max(lo - pad).max(-hi - pad);"),
    ("a point prunes patches it lies within", D,
     "    critical_along(face, p - face.centre, face.radius)",
     "    critical_along(face, p - face.centre, 0.0)"),
    # Translates trimmed in another chart (#227).
    ("a turned face matched as a translate", T,
     "        .all(|(p, q)| (*p - *q).length() <= GATE)",
     "        .all(|_| true)"),
    ("a cylinder shifted round its axis", T,
     "close(p.radius, q.radius)).then_some([false, true])",
     "close(p.radius, q.radius)).then_some([true, true])"),
    ("any trim matched as a translate", T,
     "        delta = delta.max(gap);",
     "        delta = delta.max(0.0);"),
    ("a plane re-charted without the trim shift", T,
     "frame: moved(&p.frame, p.frame.x * s.x + p.frame.y * s.y),",
     "frame: moved(&p.frame, Vec3::ZERO),"),
    ("re-charted bound without |S_A - S'|", T,
     "        let bound = matched + lipschitz * self.delta + self.slack;",
     "        let bound = lipschitz * self.delta + self.slack;"),
    ("a cone's apex move not folded in", T,
     "            delta = delta.max(pole + 1e-12 * (p.radius / pa).abs());",
     "            delta = delta.max(0.0 * pole);"),
    ("a cone of another radius matched", T,
     "        if !close(p.radius, q.radius) {",
     "        if false {"),
    ("a conic edge's support taken at its ends only", T,
     "            if inside <= hi {",
     "            if false {"),
    ("support point toward t, not against it", H,
     "        let Some((p, edge)) = support_point(source.brep, -t, edges) else {",
     "        let Some((p, edge)) = support_point(source.brep, t, edges) else {"),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode

survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        code = 0
        for target in TESTS:
            try:
                code = run(target)
            except subprocess.TimeoutExpired:
                code = -1
            if code != 0:
                break
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
