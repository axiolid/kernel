"""Mutation probe for curve-relation boundaries of curve-bounded planes (#255).

Each mutant misreads a composite or a trim, closes what must be refused,
accepts a relation kind it cannot read, or reports a deviation smaller than
the mesh's, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/execution/compile/src/bounded.rs"
R = "crates/execution/compile/src/bounded/relation.rs"
D = "crates/execution/compile/src/deviation/paths.rs"
TESTS = ["-p", "axiolid-mesh-compile", "--test", "curve_bounded_relations"]

MUTANTS = [
    ('segment sense ignored', R,
     '                    ..*segment\n',
     '                    same_sense: true,\n                    ..*segment\n', TESTS),
    ('2D point selectors not lifted', R,
     '        TrimSelector::Point2(p) => TrimSelector::Point3(lift_point(*p)),\n',
     '', TESTS),
    ('circle frame lifted with its axes swapped', R,
     '        x,\n        y,\n        z: x.cross(y),',
     '        x: y,\n        y: x,\n        z: x.cross(y),', TESTS),
    ('line direction lifted at unit speed', R,
     '            direction: lift_vector(line.direction),',
     '            direction: lift_vector(line.direction).normalize(),', TESTS),
    ('offset relation read as its basis', R,
     '        Some(GeometryNode::CurveRelation(relation)) => {',
     '        Some(GeometryNode::CurveRelation(CurveRelation::Offset { basis, .. })) => {\n'
     '            return lift_node(graph, boundary, *basis, builder, copied, depth + 1)\n'
     '        }\n'
     '        Some(GeometryNode::CurveRelation(relation)) => {', TESTS),
    ('open relation boundary closed silently', B,
     '            if gap.is_nan() || gap > linear {',
     '            if false && (gap.is_nan() || gap > linear) {', TESTS),
    ('relation boundary off the parameter plane flattened', B,
     '            let points = relation::points(graph, id, options)?;',
     '            let points: Vec<Point3> = relation::points(graph, id, options)?\n'
     '                .into_iter()\n'
     '                .map(|p| Point3::new(p.x, p.y, 0.0))\n'
     '                .collect();', TESTS),
    ('outer loop reversed about its last point', B,
     '            ring[1..].reverse();',
     '            ring.reverse();', TESTS),
    ('certified leaves reported exact', R,
     '            DeviationBound::Proven(crate::compiler::chord_error(options))',
     '            DeviationBound::Proven(0.0)', TESTS),
    ('uncertified leaves reported within the budget', R,
     '            if axiolid_reference::bound::certifies_flattening3(curve) =>',
     '            if true =>', TESTS),
    ('trims reported exact', R,
     '            leaf_bound(graph, *basis, options, depth + 1)',
     '            DeviationBound::Proven(0.0)', TESTS),
    ('joint gaps left out of the bound', R,
     '                (DeviationBound::Proven(d), Some(gap)) => DeviationBound::Proven(d + gap),',
     '                (DeviationBound::Proven(d), Some(_)) => DeviationBound::Proven(d),', TESTS),
    ('relation boundaries not dispatched to their bound', D,
     '                crate::bounded::relation::bound(graph, id, options)',
     '                DeviationBound::Proven(0.0)', TESTS),
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
