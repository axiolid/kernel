"""Mutation probe for stations and station-placed sections (#241, #246).

Each mutant breaks one step a station's point, frame or meshing rests on:
the section frame's axes (left normal, reference-up lateral and up, the
banked roll, the upright plan frame), the distance refusals, the offsets
and their linear interpolation, the frame selection, the refinement
between stations, the ring-structure check, and the model's validation of
station runs and tags; and for #246 the explicit orientation (its
Gram-Schmidt, its reading in the base frame, its interpolation, the
offsets staying in the base frame) and matching by tag (ring winding,
cyclic order, outer onto outer, segment sense, derived transforms,
reversed open sections, the tag rules the graph enforces). Each must turn
a test red: `evaluate/tests/station.rs` checks points and frames against
closed forms and the clothoid's Fresnel series; `compile/tests/stationed.rs`
and `compile/tests/stationed_oriented.rs` check resolved stations, swept
offset curves, spines and sheets against closed-form volumes, areas and
bounds; `graph/tests/station.rs` checks validation.

Not probed, because no test can tell them apart: clamping a distance a
hair past the end onto it (the evaluator reads the same point), the start
parameter of a polyline or B-spline (every one tested starts at 0), the
end cap's ring set (equal-structure rectangles triangulate alike), the
station interval search with two stations, the surface's re-check of tags
the graph already validated, and the sheet's winding (area is unsigned).
For #246 likewise the compiler's re-check of mixed tagging, which the graph
refuses before any compiler sees it.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/parametric/evaluate/src/station.rs"
C = "crates/execution/compile/src/station.rs"
V = "crates/representations/modeling/graph/src/validation.rs"
G = "crates/representations/modeling/graph/src/station.rs"
EVAL = ["-p", "axiolid-evaluate", "--test", "station"]
COMPILE = ["-p", "axiolid-mesh-compile", "--test", "stationed"]
MODEL = ["-p", "axiolid-model", "--test", "station"]
ORIENTED = ["-p", "axiolid-mesh-compile", "--test", "stationed_oriented"]

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
    ('mixed tagging accepted', V,
     '        .any(|section| section.is_empty() != first.is_empty())',
     '        .any(|_| false)', [MODEL, ORIENTED]),
    ('repeated tags accepted', V,
     '        if sorted.windows(2).any(|pair| pair[0] == pair[1]) {',
     '        if false {', [MODEL]),
    ('different tag sets accepted', V,
     '        if set(section)? != expected {',
     '        if false {', [MODEL, ORIENTED]),
    ('open sections in any tag order accepted', V,
     '        if open && *section != *first && !section.iter().eq(first.iter().rev()) {',
     '        if false {', [MODEL, ORIENTED]),
    ('reversed open sections refused', V,
     '        if open && *section != *first && !section.iter().eq(first.iter().rev()) {',
     '        if open && *section != *first {', [MODEL, ORIENTED]),
    ('station orientation not validated', V,
     '            validate_orientation(&value.orientation)\n',
     '            Ok(())\n', [MODEL, ORIENTED]),
    ('section orientation not validated', V,
     '                expect_reference(nodes, section.profile, ExpectedReference::Profile)?;\n'
     '                validate_orientation(&section.orientation)?;',
     '                expect_reference(nodes, section.profile, ExpectedReference::Profile)?;',
     [MODEL, ORIENTED]),
    ('graph accepts a parallel orientation', G,
     '        if axis.cross(reference).length() <= ORIENTATION_TOLERANCE {',
     '        if false {', [MODEL, ORIENTED]),
    ('graph defaults the reference to lateral', G,
     '            self.ref_direction.unwrap_or(Vec3::X),',
     '            self.ref_direction.unwrap_or(Vec3::Y),', [MODEL, ORIENTED]),
    ('orientation read in world axes', E,
     '        let world = |v: Vec3| v.x * self.tangent + v.y * self.lateral + v.z * self.up;',
     '        let world = |v: Vec3| v;', [EVAL, ORIENTED]),
    ('reference not made perpendicular to the axis', E,
     '        let x = (reference - reference.dot(axis) * axis).normalize();',
     '        let x = reference;', [EVAL, ORIENTED]),
    ('oriented lateral reversed', E,
     '        let y = axis.cross(x);',
     '        let y = x.cross(axis);', [EVAL, ORIENTED]),
    ('oriented up ignores the axis', E,
     '            up: world(axis),',
     '            up: self.up,', [EVAL, ORIENTED]),
    ('evaluator accepts a parallel orientation', E,
     '        if axis.cross(reference).length() <= ORIENTATION_TOLERANCE {',
     '        if false {', [EVAL]),
    ('section orientation ignored', C,
     '    if a.is_base() && b.is_base() {',
     '    if true {', [ORIENTED]),
    ('orientation not interpolated', C,
     '    base.oriented(Some(axis_a.lerp(axis_b, u)), Some(ref_a.lerp(ref_b, u)))',
     '    base.oriented(Some(axis_a), Some(ref_a))', [ORIENTED]),
    ('offsets read in the turned frame', C,
     '    base.place(offsets.lateral, offsets.vertical, offsets.longitudinal)\n',
     '    turned.place(offsets.lateral, offsets.vertical, offsets.longitudinal)\n', [ORIENTED]),
    ('resolved station frame not turned', C,
     '    let mut placed = oriented(&section, &orientation, &orientation, 0.0)?.frame();',
     '    let mut placed = section.frame();', [ORIENTED]),
    ('tagged rings not rewound', C,
     '        if (area > 0.0) != (k == 0) {',
     '        if false {', [ORIENTED]),
    ('tag count not checked', C,
     '    if count != tags.len() {',
     '    if false {', [ORIENTED]),
    ('segment sense ignored', C,
     '        if !segment.same_sense {',
     '        if false {', [ORIENTED]),
    ('derived transform ignored', C,
     '                *point = transform.transform_point2(*point);',
     '                let _ = transform;', [ORIENTED]),
    ('cyclic order not checked', C,
     '        let cyclic = have.len() == len && (0..len).all(|m| have[(start + m) % len] == want[m]);',
     '        let cyclic = have.len() == len;', [ORIENTED]),
    ('outer may match a hole', C,
     '        if (r == 0) != (found == 0) || !cyclic {',
     '        if !cyclic {', [ORIENTED]),
    ('reversed open section joined forwards', C,
     '        if section.tags.iter().eq(first.iter().rev()) {\n            points.reverse();',
     '        if section.tags.iter().eq(first.iter().rev()) {', [ORIENTED]),
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
