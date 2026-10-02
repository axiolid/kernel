"""Mutation probe for rigid placement of exact B-reps and exact swept disks
(#223).

Each mutant moves a support wrongly, keeps a reflected solid inward,
misreads a pcurve on a reflected face, accepts a transform that is not
rigid, or places a swept disk off its directrix, and must turn a test red.
`brep/tests/transform.rs` checks each support; `construct/tests/placement.rs`
and `compile/tests/exact_placement.rs` audit, measure and certify distances
of placed solids against closed forms.

Equivalent mutants, deliberately not listed: dropping the swept disk's own
"reaching the axis" check (the revolution refuses the same disk by name),
and building the segment frame left-handed (the reflection path then yields
the same rotationally symmetric solid).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
T = "crates/representations/brep/src/transform.rs"
S = "crates/algorithms/construction/construct/src/swept_disk_exact.rs"
C = "crates/execution/compile/src/exact.rs"
D = "crates/execution/compile/src/directrix.rs"
TESTS = [
    ["-p", "axiolid-brep", "--test", "transform"],
    ["-p", "axiolid-construct", "--test", "placement"],
    ["-p", "axiolid-mesh-compile", "--test", "exact_placement"],
]

MUTANTS = [
    ("reflection read as a rotation", T,
     "            mirror: m.determinant() < 0.0,",
     "            mirror: false,"),
    ("reflected faces kept forward", T,
     "                (Orientation::Forward, true) => Orientation::Reversed,",
     "                (Orientation::Forward, true) => Orientation::Forward,"),
    ("reflected plane frame left-handed", T,
     "            z: if self.mirror { -z } else { z },",
     "            z,"),
    ("reflected axial frame left-handed", T,
     "            y: if self.mirror { -y } else { y },",
     "            y,"),
    ("pcurves on a reflected curved face kept", T,
     "            let reflect = axial.get(surface.index()).copied().unwrap_or(false);",
     "            let reflect = false;"),
    ("pcurves reflected without the turn offset", T,
     "    Point2::new(TAU - p.x, p.y)",
     "    Point2::new(-p.x, p.y)"),
    ("reflected conic frame left-handed", T,
     "        y: -reflect_vector(frame.y),",
     "        y: reflect_vector(frame.y),"),
    ("reflected circle pcurve keeps its interval", T,
     "                frame: reflect_conic_frame(&circle.frame),\n                ..*circle\n            }),\n            Reparam::Negated,",
     "                frame: reflect_conic_frame(&circle.frame),\n                ..*circle\n            }),\n            Reparam::Same,"),
    ("reflected wave keeps its sine", T,
     "                sine: -wave.sine,",
     "                sine: wave.sine,"),
    ("reflected wave interval negated without the turn", T,
     "            Self::Reflected => Interval::new(TAU - interval.start, TAU - interval.end),",
     "            Self::Reflected => Interval::new(-interval.start, -interval.end),"),
    ("a scale or shear accepted", T,
     "                if (columns[i].dot(columns[j]) - want).abs() > RIGID_TOLERANCE {",
     "                if false {"),
    ("a non-finite transform accepted", T,
     "        if !columns.iter().all(|c| c.is_finite()) || !transform.translation.is_finite() {",
     "        if false {"),
    ("vertices left in place", T,
     "                position: rigid.point(vertex.position),",
     "                position: vertex.position,"),
    ("a line direction translated", T,
     "                origin: rigid.point(line.origin),\n                direction: rigid.vector(line.direction),",
     "                origin: rigid.point(line.origin),\n                direction: rigid.point(line.direction),"),
    ("an unreflectable pcurve copied", T,
     "        _ => {\n            return Err(TransformError::Unsupported(\n                \"a pcurve family on a curved face under a reflection\",",
     "        _ => {\n            return Ok((curve.clone(), Reparam::Same));\n            #[allow(unreachable_code)]\n            return Err(TransformError::Unsupported(\n                \"a pcurve family on a curved face under a reflection\","),
    ("arc sweep turned the wrong way", S,
     "    let rotation = Mat3::from_cols(radial, frame.z, -tangent);",
     "    let rotation = Mat3::from_cols(radial, frame.z, tangent);"),
    ("arc sweep not moved onto the arc", S,
     "    let translation = frame.origin + radial * bend;",
     "    let translation = frame.origin;"),
    ("arc sweep started at angle zero", S,
     "    let (sin, cos) = span.start.sin_cos();",
     "    let (sin, cos) = 0.0_f64.sin_cos();"),
    ("hollow disk swept solid", S,
     "                thickness: Some(radius - inner),",
     "                thickness: None,"),
    ("segment sweep not moved to its start", S,
     "        &Transform3::from_mat3_translation(Mat3::from_cols(x, y, z), start),",
     "        &Transform3::from_mat3_translation(Mat3::from_cols(x, y, z), Vec3::ZERO),"),
    ("instance transform ignored", C,
     "                    .transformed(&instance.transform)",
     "                    .transformed(&axiolid_core::Transform3::IDENTITY)"),
    ("a scaled instance not named", C,
     "                            unsupported(\"exact instance under a scaled or sheared transform\")",
     "                            unsupported(\"instance\")"),
    ("line directrix range ignored", D,
     "                line.origin + line.direction * start,",
     "                line.origin,"),
    ("circle directrix range ignored", D,
     "                    axiolid_core::Interval::new(lo.max(0.0), hi.min(TAU))",
     "                    axiolid_core::Interval::new(0.0, TAU)"),
    ("trimmed arc read without wrapping", D,
     "                        axiolid_core::Interval::new(lo, hi),",
     "                        axiolid_core::Interval::new(a, b),"),
    ("a cornered directrix swept as its first segment", D,
     "            if polyline.points.len() == 2 && !polyline.closed && range.is_none() =>",
     "            if polyline.points.len() >= 2 && !polyline.closed && range.is_none() =>"),
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
