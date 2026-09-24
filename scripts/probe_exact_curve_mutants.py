"""Mutation probe for exact curve intersection (#119): each mutant must be caught."""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/parametric/nurbs/src/exact_curve_intersection.rs"

MUTANTS = [
    ("multiplicity stops at one",
     "    while !q.is_zero() && root.sign_of(&q) == Sign::Zero {",
     "    while false && root.sign_of(&q) == Sign::Zero {"),
    ("cone nappe condition ignored",
     "        if !allowed(&root) {\n            continue;\n        }",
     "        if false && !allowed(&root) {\n            continue;\n        }"),
    ("antipode never reported",
     "        if lost >= 1 && on_nappe {",
     "        if lost >= 2 && on_nappe {"),
    ("common roots from first equation only",
     "        .fold(polys[0].0.clone(), |acc, (p, _)| acc.gcd(p));",
     "        .fold(polys[0].0.clone(), |acc, (_p, _)| acc);"),
    ("partial overlap accepted",
     "        .any(|root| multiplicity(root, &int) % 2 == 1);",
     "        .any(|root| multiplicity(root, &int) % 2 == 7);"),
    ("theta order ignores the lower half",
     "        let band = if is_line || root.cmp_dyadic(&Dyadic::zero()) != Sign::Negative {",
     "        let band = if true || root.cmp_dyadic(&Dyadic::zero()) != Sign::Negative {"),
]

def run():
    return subprocess.run(
        ["cargo", "test", "-q", "--release", "-p", "axiolid-nurbs",
         "--test", "exact_curve_intersection"],
        cwd=ROOT, capture_output=True, text=True, timeout=900,
    ).returncode

survivors = []
original = SRC.read_text()
for name, old, new in MUTANTS:
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    SRC.write_text(original.replace(old, new))
    try:
        try:
            code = run()
        except subprocess.TimeoutExpired:
            code = -1  # a hang is a detection
    finally:
        SRC.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
