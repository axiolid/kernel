"""Mutation probe for stations and station-placed sections (#241).

Each mutant breaks one step a station's point, frame or meshing rests on:
the section frame's axes (left normal, reference-up lateral and up, the
banked roll, the upright plan frame), the distance refusals, the offsets
and their linear interpolation, the frame selection, the refinement
between stations, the ring-structure check, and the model's validation of
station runs and tags. Each must turn a test red: `evaluate/tests/station.rs`
checks points and frames against closed forms and the clothoid's Fresnel
series; `compile/tests/stationed.rs` checks resolved stations, swept
offset curves, spines and a sheet against closed-form volumes, areas and
bounds; `graph/tests/station.rs` checks validation.

Not probed, because no test can tell them apart: clamping a distance a
hair past the end onto it (the evaluator reads the same point), the start
parameter of a polyline or B-spline (every one tested starts at 0), the
end cap's ring set (equal-structure rectangles triangulate alike), the
station interval search with two stations, the surface's re-check of tags
the graph already validated, and the sheet's winding (area is unsigned).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/parametric/evaluate/src/station.rs"
C = "crates/execution/compile/src/station.rs"
V = "crates/representations/modeling/graph/src/validation.rs"
EVAL = ["-p", "axiolid-evaluate", "--test", "station"]
COMPILE = ["-p", "axiolid-mesh-compile", "--test", "stationed"]
MODEL = ["-p", "axiolid-model", "--test", "station"]

MUTANTS = [
    ('planar lateral to the right', E,
     '        lateral: Vec3::new(-t.y, t.x, 0.0),',
     '        lateral: Vec3::new(t.y, -t.x, 0.0),', [EVAL]),
    ('reference-up lateral to the right', E,
     '        lateral: -right,',
     '        lateral: right,', [EVAL]),
    ('reference-up up does not lean with the grade', E,
     '        up: right.cross(tangent),',
     '        up: Vec3::Z,', [EVAL]),
    ('banked lateral unrolled the wrong way', E,
     '                lateral: section.lateral,',
     '                lateral: -section.lateral,', [EVAL]),
    ('plan frame lateral to the right', E,
     '            lateral: Vec3::new(-tangent.y, tangent.x, 0.0),',
     '            lateral: Vec3::new(tangent.y, -tangent.x, 0.0),', [EVAL]),
    ('negative distance accepted', E,
     '    if distance < 0.0 {',
     '    if false {', [EVAL]),
    ('distance beyond the length accepted', E,
     '            if distance > length + slack {',
     '            if false {', [EVAL, COMPILE]),
    ('lateral offset dropped', C,
     '        local.x + offsets.lateral,',
     '        local.x,', [COMPILE]),
    ('vertical offset dropped', C,
     '        local.y + offsets.vertical,',
     '        local.y,', [COMPILE]),
    ('longitudinal offset dropped', C,
     '        offsets.longitudinal,\n    )',
     '        0.0,\n    )', [COMPILE]),
    ('no interpolation between stations', C,
     '    a + (b - a) * u',
     '    a', [COMPILE]),
    ('plan frame ignored', C,
     '            StationFrame::Plan => section.plan(),',
     '            StationFrame::Plan => Ok(section),', [COMPILE]),
    ('midpoint test skipped', C,
     '        .all(|((p, q), m)| (*m - 0.5 * (*p + *q)).length() <= chord);',
     '        .all(|_| true);', [COMPILE]),
    ('twist test skipped', C,
     '            .all(|&(i, j)| quad_spread(lo[i], lo[j], hi[i], hi[j]) <= chord)',
     '            .all(|_| true)', [COMPILE]),
    ('planar walls refined anyway', C,
     '        twist.min(normal.dot(d - a).abs())',
     '        twist', [COMPILE]),
    ('ring structure not checked', C,
     '    if rings.iter().any(|ring| !same_structure(ring, &rings[0])) {',
     '    if false {', [COMPILE]),
    ('equal distances accepted', V,
     '        if previous.is_some_and(|before| station.distance <= before) {',
     '        if previous.is_some_and(|before| station.distance < before) {', [MODEL]),
    ('one station accepted', V,
     '    if stations.len() < 2 {',
     '    if stations.len() < 1 {', [MODEL]),
    ('non-finite offsets accepted', V,
     '    if !station.offsets.is_finite() {',
     '    if false {', [MODEL]),
    ('inconsistent tags accepted', V,
     '    if sections.iter().any(|section| section.tags != first.tags) {',
     '    if false {', [MODEL]),
    ('repeated tags accepted', V,
     '    if !first.tags.iter().all(|tag| seen.insert(tag.as_str())) {',
     '    if false {', [MODEL]),
    ('offset curve read as 2D', V,
     '            GeometryNode::CurveRelation(CurveRelation::OffsetByStations { .. }) => {\n'
     '                if dimension != CurveDimension::Three {',
     '            GeometryNode::CurveRelation(CurveRelation::OffsetByStations { .. }) => {\n'
     '                if dimension != CurveDimension::Two {', [MODEL]),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode


survivors = []
for name, rel, old, new, targets in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    code = 0
    try:
        for target in targets:
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
