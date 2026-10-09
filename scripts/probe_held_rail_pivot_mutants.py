"""Mutation probe for a banked pivot rotating about a held rail (#279,
ADR 0081 amendment).

Each mutant misreads which rail is held, drops the half of the cant the
pivot follows, its rate or the `cos(psi)` of an angle piece's rate, loses
a term of the derived second derivative, misstates a Viennese bend's
derivative suprema or their critical points, halves the pivot's bounds
wrongly, forgets the cant seams the derived pivot inherits (as breaks or
as grade corners), reads a seam's rate from the wrong side, or accepts a
held-rail piece in a cant law, and must turn a test red.
`curve/tests/banked.rs` pins the derived elevation and rate against an
independently written Viennese bend, `evaluate/tests/banked.rs` the held
rail's height and the other rail's rise `b sin(psi)`, and
`evaluate/tests/elevated.rs` the derivatives against finite differences
and the bounds against dense sampling.

Not probed, because no test can tell them apart from the original on
gently canted track:

- dropping `P_1^2` from the second-derivative bound: with `psi'` of
  order `1e-3` per metre it adds a few percent at most, inside the slack
  the bound checks allow. It is kept because the bound must hold for any
  cant law, not because a test needs it;
- reading the cant's bounds across one of its seams: `banked_breaks`
  names those seams, and every caller refuses a span across a break
  before it asks for a bound.

The sweep certification along a held-rail bend
(`mesh-compile/tests/elevated_directrix.rs`) rests on the same bounds and
is checked there, not re-run per mutant.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
L = "crates/representations/analytic/curve/src/banked.rs"
BK = "crates/algorithms/parametric/evaluate/src/banked.rs"
EL = "crates/algorithms/parametric/evaluate/src/elevated.rs"
CURVE = ["-p", "axiolid-curve", "--test", "banked"]
EVB = ["-p", "axiolid-evaluate", "--test", "banked"]
EVE = ["-p", "axiolid-evaluate", "--test", "elevated"]

MUTANTS = [
    ("left rail held with the right rail's sign", L,
     "            Self::Left => -1.0,",
     "            Self::Left => 1.0,", [CURVE, EVB]),
    ("right rail held with the left rail's sign", L,
     "            Self::Right => 1.0,",
     "            Self::Right => -1.0,", [CURVE, EVB]),
    ("pivot follows the whole cant, not half", L,
     "            return Ok((elevation + half * cant, half * rate));",
     "            return Ok((elevation + 2.0 * half * cant, half * rate));", [CURVE, EVB]),
    ("held rail's elevation ignored", L,
     "            return Ok((elevation + half * cant, half * rate));",
     "            return Ok((half * cant, half * rate));", [CURVE, EVB]),
    ("derived pivot without a rate", L,
     "            return Ok((elevation + half * cant, half * rate));",
     "            return Ok((elevation + half * cant, 0.0));", [CURVE, EVB]),
    ("angle piece's cant rate without cos(psi)", L,
     "            CantValue::Angle(rate) => self.rail_head_distance * psi.cos() * rate,",
     "            CantValue::Angle(rate) => self.rail_head_distance * rate,", [CURVE, EVB]),
    ("held-rail piece accepted in a cant law", L,
     "        if self.cant.has_rail_pieces() {\n            return Err(BankError::RailInCant);",
     "        if false {\n            return Err(BankError::RailInCant);", [CURVE]),
    ("non-finite held-rail elevation well formed", L,
     "            CantForm::AboutRail { elevation, .. } => elevation.is_finite(),",
     "            CantForm::AboutRail { .. } => true,", [CURVE]),
    ("evaluator places a station on a held-rail cant law", BK,
     "    if curve.cant.has_rail_pieces() {\n        return Err(refused(BankError::RailInCant));",
     "    if false {\n        return Err(refused(BankError::RailInCant));", [EVB]),
    ("derived second derivative not halved", BK,
     "        let out = [elevation + half * cant, half * rate, half * bend];",
     "        let out = [elevation + half * cant, half * rate, bend];", [EVE]),
    ("derived second derivative drops -sin(psi) psi'^2", BK,
     "            b * (cos * accel - sin * rate * rate),",
     "            b * cos * accel,", [EVE]),
    ("viennese psi'' bent the wrong way", BK,
     "        let bend = 420.0 * xi2 * one_minus * one_minus * (1.0 - 2.0 * xi);",
     "        let bend = 420.0 * xi2 * one_minus * one_minus * (2.0 * xi - 1.0);", [EVE]),
    ("third-derivative bound drops its cross terms", BK,
     "            b * (p3 + 3.0 * p1 * p2 + p1 * p1 * p1),",
     "            b * p3,", [EVE]),
    ("derived bounds quartered", BK,
     "        return cant_bounds(curve, lo, hi).map(|bounds| bounds.map(|v| 0.5 * v));",
     "        return cant_bounds(curve, lo, hi).map(|bounds| bounds.map(|v| 0.25 * v));", [EVE]),
    ("viennese f' supremum understated", BK,
     "            1 => 140.0 * w * w * w,",
     "            1 => 100.0 * w * w * w,", [EVE]),
    ("viennese f'' supremum understated", BK,
     "            2 => 420.0 * w * w * (1.0 - 2.0 * xi),",
     "            2 => 300.0 * w * w * (1.0 - 2.0 * xi),", [EVE]),
    ("viennese f''' supremum understated", BK,
     "            _ => 840.0 * w * (1.0 - 5.0 * w),",
     "            _ => 420.0 * w * (1.0 - 5.0 * w),", [EVE]),
    ("viennese critical points ignored", BK,
     "        .filter(|&xi| xi > lo && xi < hi)",
     "        .filter(|_| false)", [EVE]),
    ("viennese f'' critical points misplaced", BK,
     "        2 => &[0.5 - 0.223_606_797_749_979, 0.5 + 0.223_606_797_749_979],",
     "        2 => &[0.5],", [EVE]),
    ("cant seams not the derived pivot's breaks", EL,
     "    out.extend(crate::banked::derived_cant_seams(curve));",
     "", [EVE]),
    ("cant seams not read for grade corners", EL,
     "        seams.extend(crate::banked::derived_cant_seams(banked));",
     "", [EVE]),
    ("every cant seam taken as the pivot's", BK,
     "        .filter(|&seam| spans.iter().any(|&(a, b)| seam >= a && seam <= b))",
     "        .filter(|_| true)", [EVE]),
    ("derived rate at a corner read after the seam", BK,
     "    let (piece, s) = piece_at(&curve.cant, d, before)?;",
     "    let (piece, s) = piece_at(&curve.cant, d, false)?;", [EVE]),
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
