"""Mutation probe for banked curves: cant laws, bank conventions and section
frames (#240, ADR 0081).

Each mutant misreads a cant form, owns a seam from the wrong side,
extrapolates a law, lets a cant exceed the rail heads, confuses the two bank
conventions, drops the pivot, rolls the frame the wrong way, or wires the
curve wrongly into evaluation, and must turn a test red.
`curve/tests/banked.rs` pins law values against hand-computed ones,
`evaluate/tests/banked.rs` pins sections and frames, and
`brep/tests/transform.rs` the refusal under a rigid motion.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
L = "crates/representations/analytic/curve/src/banked.rs"
E = "crates/algorithms/parametric/evaluate/src/banked.rs"
P = "crates/algorithms/parametric/evaluate/src/provider.rs"
C = "crates/algorithms/parametric/evaluate/src/curve.rs"
B = "crates/algorithms/parametric/evaluate/src/bound.rs"
T = "crates/representations/brep/src/transform.rs"
TESTS = [
    ["-p", "axiolid-curve", "--test", "banked"],
    ["-p", "axiolid-evaluate", "--test", "banked"],
    ["-p", "axiolid-brep", "--test", "transform"],
]

MUTANTS = [
    ("bloss read as a smoothstep of the wrong weight", L,
     "coefficients: vec![start, 0.0, 3.0 * change, -2.0 * change],",
     "coefficients: vec![start, 0.0, 2.0 * change, -2.0 * change],"),
    ("helmert second half bent the wrong way", L,
     "coefficients: vec![start + 0.5 * change, change, -0.5 * change],",
     "coefficients: vec![start + 0.5 * change, change, 0.5 * change],"),
    ("cosine transition without the half", L,
     "start + change * 0.5 * (1.0 - phase.cos()),",
     "start + change * (1.0 - phase.cos()),"),
    ("cosine rate without the half", L,
     "change * 0.5 * core::f64::consts::PI * phase.sin(),",
     "change * core::f64::consts::PI * phase.sin(),"),
    ("sine transition over pi instead of two pi", L,
     "start + change * (xi - phase.sin() / core::f64::consts::TAU),",
     "start + change * (xi - phase.sin() / core::f64::consts::PI),"),
    ("sine rate sign flipped", L,
     "                    change * (1.0 - phase.cos()),",
     "                    change * (1.0 + phase.cos()),"),
    ("viennese bend coefficient misread", L,
     "(35.0 + xi * (-84.0 + xi * (70.0 - 20.0 * xi)))",
     "(35.0 + xi * (-80.0 + xi * (70.0 - 20.0 * xi)))"),
    ("viennese bend rate misread", L,
     "140.0 * xi2 * xi",
     "120.0 * xi2 * xi"),
    ("polynomial slope without the power", L,
     "accumulated * xi + c * power as Scalar",
     "accumulated * xi + c"),
    ("rate per xi, not per plan distance", L,
     "let rate = slope / length;",
     "let rate = slope;"),
    ("a seam owned by the earlier piece", L,
     "if distance < end || (index == last && distance <= end) {",
     "if distance <= end {"),
    ("the last piece extrapolated", L,
     "if distance < end || (index == last && distance <= end) {",
     "if distance < end || index == last {"),
    ("a negative distance accepted", L,
     "if !distance.is_finite() || distance < 0.0 || !self.is_well_formed() {",
     "if !distance.is_finite() || !self.is_well_formed() {"),
    ("a cant beyond the rail heads accepted", L,
     "    if cant.abs() > rail_head_distance {",
     "    if false {"),
    ("vertical rise rolled as a tangent rotation", L,
     "Ok((cant / span).asin())",
     "Ok((cant / rail_head_distance).asin())"),
    ("a cant beyond the vertical span accepted", L,
     "if span <= 0.0 || cant.abs() > span {",
     "if span <= 0.0 {"),
    ("an angle beyond a quarter turn accepted", L,
     "        if psi.abs() > core::f64::consts::FRAC_PI_2 {\n            return Err(BankError::AngleOutOfRange { angle: psi });\n        }\n        match self {",
     "        match self {"),
    ("an angle piece accepted as a pivot", L,
     "        if self.pivot.has_angle_pieces() {",
     "        if false {"),
    ("an angle piece's cant read as the angle", L,
     "Ok(self.rail_head_distance * psi.sin())",
     "Ok(psi)"),
    ("pivot not added to the point", E,
     "Ok(base + pivot * Vec3::Z)",
     "Ok(base)"),
    ("pivot rate left out of the tangent", E,
     "let derivative = Vec3::new(plan.x, plan.y, grade + pivot_rate);",
     "let derivative = Vec3::new(plan.x, plan.y, grade);"),
    ("section rolled the wrong way", E,
     "let lateral = cos * normal + sin * square_up;",
     "let lateral = cos * normal - sin * square_up;"),
    ("section up not square to the lateral", E,
     "let up = -sin * normal + cos * square_up;",
     "let up = sin * normal + cos * square_up;"),
    ("normal to the right", E,
     "let normal = Vec3::Z.cross(plan).normalize();",
     "let normal = plan.cross(Vec3::Z).normalize();"),
    ("grade ignored by the vertical rise", E,
     "let grade_cosine = 1.0 / grade.hypot(1.0);",
     "let grade_cosine = 1.0;"),
    ("frame left-handed", E,
     "            z: -self.lateral,",
     "            z: self.lateral,"),
    ("rail heads swapped", E,
     "(self.point + half, self.point - half)",
     "(self.point - half, self.point + half)"),
    ("a station off the cant law placed", E,
     "    covered(curve, d)?;\n    let (pivot, _) = curve.pivot_at(d).map_err(refused)?;\n    let base",
     "    let (pivot, _) = curve.pivot_at(d).map_err(refused)?;\n    let base"),
    ("evaluator frames a banked curve reference-up", P,
     "        if let Curve3::Banked(banked) = curve {\n            return self.banked_frame(banked, at);\n        }\n",
     ""),
    ("another reference up accepted", P,
     "        if self.up != Vec3::Z {",
     "        if false {"),
    ("banked distance convention unsupported", P,
     "Curve3::Elevated(_) | Curve3::Banked(_) => DistanceConvention::PlanDistance,",
     "Curve3::Elevated(_) => DistanceConvention::PlanDistance,"),
    ("banked domain empty", C,
     "            end: b.span(),",
     "            end: 0.0,"),
    ("banked derivative made unit", C,
     "Curve3::Banked(b) => crate::banked::banked_derivative(b, t),",
     "Curve3::Banked(b) => crate::banked::banked_tangent(b, t),"),
    ("banked flattening claimed certified", B,
     "        Curve3::Banked(_) => false,",
     "        Curve3::Banked(_) => true,"),
    ("banked curve refused under another name", T,
     "\"a banked alignment curve (centreline and cant law)\",",
     "\"an unknown curve family\","),
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
