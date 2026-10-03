"""Mutation probe for plan distance between exact bodies (#217) and
between bodies of several items (#237, `exact_bodies/plan.rs`).

Each mutant makes a plan bound unsound -- a lower bound above the plan
distance, an upper bound not between boundary points' projections --
certifies an overlap or a gap it has not shown, measures a body without
its placement or in space, names the wrong item, or drops the search over
level faces that shows an overlap the distance search leaves undecided,
and must turn a test red. The closed forms are in
`construct/tests/plan_distance.rs` and `construct/tests/plan_bodies.rs`.

Equivalent mutants, deliberately not listed: the sphere fallback of a
horizontal projection using the space radius (only looser), the overlap
margin scaled up tenfold (only more cautious), and vertical faces let
into the level-face search (their shadows have no area, so no clip of
theirs clears the margin).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/query/measure/src/exact_distance.rs"
B = "crates/algorithms/query/measure/src/exact_bodies.rs"
P = "crates/algorithms/query/measure/src/exact_bodies/plan.rs"
TESTS = [
    ["-p", "axiolid-construct", "--test", "plan_distance"],
    ["-p", "axiolid-construct", "--test", "plan_bodies"],
    ["-p", "axiolid-measure", "--all-features", "--lib"],
]

MUTANTS = [
    ("the plan gap measured in space", D,
     "        (flat(b.centre - a.centre), a.plan_radius + b.plan_radius)\n    } else {",
     "        (b.centre - a.centre, a.plan_radius + b.plan_radius)\n    } else {"),
    ("plan discs of no size", D,
     "        (flat(b.centre - a.centre), a.plan_radius + b.plan_radius)\n    } else {",
     "        (flat(b.centre - a.centre), 0.0)\n    } else {"),
    ("normals read whole in plan", D,
     "        if !plan {\n            return Some(n);\n        }",
     "        if true {\n            return Some(n);\n        }"),
    ("horizontal turning bounded by the larger axis alone", D,
     "        (a * flat_length(frame.x)).hypot(b * flat_length(frame.y))",
     "        (a * flat_length(frame.x)).max(b * flat_length(frame.y))"),
    ("a plan patch's vertical speed dropped", D,
     "            turning(&c.frame, c.radius.abs(), c.radius.abs()),\n            flat_length(c.frame.z),",
     "            turning(&c.frame, c.radius.abs(), c.radius.abs()),\n            0.0,"),
    ("plan witnesses measured in space", D,
     "        Metric::Plan => flat_length(p - q),",
     "        Metric::Plan => (p - q).length(),"),
    ("an overlap shown without a margin", D,
     "            length > 0.0 && cross(e0, e1, at) / length > margin",
     "            length > 0.0 && cross(e0, e1, at) / length >= 0.0"),
    ("apart when merely not overlapping", D,
     "    if found.bounds.lower > 0.0 {\n        return Ok((\n            PlanOverlap::Disjoint {",
     "    if found.bounds.lower >= 0.0 {\n        return Ok((\n            PlanOverlap::Disjoint {"),
    ("faces pruned by the offset in space", D,
     "        Metric::Plan => (\n            flat(other.centre - face.centre),",
     "        Metric::Plan => (\n            other.centre - face.centre,"),
    # The level-face search for an overlap (#237).
    ("no level-face search after an undecided distance search", D,
     "    Ok(match level_overlap(a, b, tolerance, OVERLAP_STEPS)? {",
     "    Ok(match Option::<(Point2, Pair)>::None {"),
    ("level-face pairs refined finest first", D,
     "            heap.push((Key(ea.plan_radius + eb.plan_radius), i, j));",
     "            heap.push((Key(-(ea.plan_radius + eb.plan_radius)), i, j));"),
    ("level-face pairs kept with their shadows apart", D,
     "        if lower_bound(&side_a, ea, &side_b, eb)? <= 0.0 {",
     "        if lower_bound(&side_a, ea, &side_b, eb)? <= Scalar::INFINITY {"),
    # Bodies of several items in plan (#237).
    ("a body's plan distance measured in space", P,
     "        Metric::Plan,\n        &mut |lower, upper| upper - lower <= accuracy,",
     "        Metric::Space,\n        &mut |lower, upper| upper - lower <= accuracy,"),
    ("a body's plan clearance measured in space", P,
     "        Metric::Plan,\n        &mut |lower, upper| upper < limit || lower > limit,",
     "        Metric::Space,\n        &mut |lower, upper| upper < limit || lower > limit,"),
    ("an overlap's first item named by the second body's element", P,
     "            item_a: a.item(on_a),",
     "            item_a: a.item(on_b),"),
    ("an overlap's second item named by the first body's element", P,
     "            item_b: b.item(on_b),",
     "            item_b: b.item(on_a),"),
    ("an undecided body overlap read as disjoint", P,
     "        _ => BodyPlanOverlap::Undecided,",
     "        _ => BodyPlanOverlap::Disjoint { gap: 0.0 },"),
    ("the first body measured unplaced", B,
     "        let a = Self::new(a.items, BodySide::First)?.placed(a.placement, BodySide::First)?;",
     "        let a = Self::new(a.items, BodySide::First)?;"),
    ("the second body measured unplaced", B,
     "        let b = Self::new(b.items, BodySide::Second)?.placed(b.placement, BodySide::Second)?;",
     "        let b = Self::new(b.items, BodySide::Second)?;"),
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
        try:
            code = max(run(t) for t in TESTS)
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
