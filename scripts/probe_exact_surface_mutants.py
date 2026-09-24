"""Mutation probe for exact decisions in exact_surface_intersection (#119)."""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PATH = ROOT / "crates/algorithms/parametric/nurbs/src/exact_surface_intersection.rs"

MUTANTS = [('sphere tangency read as a circle', '        Sign::Positive => {}\n        // Exactly tangent: a touch is a point, not a curve.\n        Sign::Zero => return Err(ExactIntersectionRefusal::NotRegularCurve),', '        Sign::Positive | Sign::Zero => {}\n        // Exactly tangent: a touch is a point, not a curve.\n        Sign::Negative => return Err(ExactIntersectionRefusal::NotRegularCurve),'), ('parallel cylinder test weakened', '    if esign(&along) == Sign::Zero {\n        // Plane parallel to the axis', '    if along.to_f64().abs() < 1e-12 {\n        // Plane parallel to the axis'), ('perpendicular cylinder test dropped', '    let perpendicular = ecross_is_zero(&exact_axis, &exact_normal);\n    // cos(theta)', '    let perpendicular = false;\n    // cos(theta)'), ('cone perpendicular test dropped', '        (Some(a), Some(n)) => ecross_is_zero(&a, &n),', '        (Some(a), Some(n)) => { let _ = (a, n); false }'), ('tangent ruling counted as two', '        Sign::Positive => false,\n        Sign::Zero => true,', '        Sign::Positive => false,\n        Sign::Zero => false,')]

def run():
    return subprocess.run(
        ["cargo", "test", "-q", "-p", "axiolid-nurbs", "--test", "exact_intersection"],
        cwd=ROOT, capture_output=True, text=True, timeout=900,
    ).returncode

survivors = []
original = PATH.read_text()
try:
    for name, old, new in MUTANTS:
        assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
        PATH.write_text(original.replace(old, new))
        try:
            code = run()
        except subprocess.TimeoutExpired:
            code = -1
        status = "killed" if code != 0 else "SURVIVED"
        print(f"{status:8} {name}", flush=True)
        if code == 0:
            survivors.append(name)
finally:
    PATH.write_text(original)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
