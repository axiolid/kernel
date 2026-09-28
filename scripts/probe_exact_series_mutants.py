"""Mutation probe for certified signs of series fields (#181).

Each mutant weakens the certified tier -- the fixed-point intervals'
outward rounding, the precision raised until a sign shows, the exact zero
test's common angle and imaginary parts, the harmonics' derivatives, the
Bernstein form's remainder -- or switches the tier off, and must turn a
test red.

Equivalent mutants, deliberately not listed:

- dropping the flatness gate ([`SeriesTier::flat`]): the tier then also
  answers where `f64` and subdivision already decide, which changes cell
  boundaries and costs time but no test's verdict;
- dropping the `f64` screen before the certified form: again only cost;
- dropping the Lagrange remainder of the sine and cosine series: the
  series stops at a term below one unit of a working precision 32 bits
  finer than the one returned, so the outward rounding to the returned
  precision covers it in every case a test can reach. It stays for the
  argument's sake;
- a smaller precision margin past the corner value (`+ 64` bits): the
  Bernstein coefficients are then wider, but every box the tests certify
  clears its remainder by far more.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/parametric/nurbs/src/exact_series.rs"
E = "crates/algorithms/parametric/nurbs/src/exact_field.rs"
X = "crates/algorithms/exact/src/fixed.rs"
SERIES = ["-p", "axiolid-nurbs", "--lib", "exact_series"]
TRACE = ["-p", "axiolid-nurbs", "--lib", "implicit_trace"]
FIXED = ["-p", "axiolid-exact", "--lib", "fixed"]

MUTANTS = [
    ('no series tier', E, '            Field2::Series(f) if f.is_finite() => Some(Self::Series(SeriesTier::new(f))),', '            Field2::Series(_) => None,', TRACE),
    ('products rounded inward', X, '            lo: lo >> self.bits,', '            lo: ceil_shr(lo, self.bits),', FIXED),
    ('quotients rounded inward', X, '            hi: -((-&self.hi).div_floor(&d)),', '            hi: self.hi.div_floor(&d),', FIXED),
    ('one precision only', S, '        for bits in PRECISIONS {', '        for bits in [PRECISIONS[0]] {', SERIES),
    ('sine differentiated with the wrong sign', S, '[(1, sin), (1, cos), (-1, sin), (-1, cos)][i % 4]', '[(1, sin), (-1, cos), (-1, sin), (1, cos)][i % 4]', SERIES),
    ('imaginary parts ignored', S, '.all(|(re, im)| re.sign() == Some(Sign::Zero) && im.sign() == Some(Sign::Zero)),', '.all(|(re, _)| re.sign() == Some(Sign::Zero)),', SERIES),
    ('angles not brought to one', S, '            (shift < 60).then_some(m << shift)', '            (shift < 60).then_some(m)', SERIES),
    ('Bernstein remainder left out', S, '                    remainder = (remainder + b * weight(i, j).hi).next_up();', '                    let _ = (b, remainder);', SERIES),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
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
