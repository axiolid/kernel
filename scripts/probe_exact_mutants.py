"""Mutation probe for axiolid-exact: each mutant must fail at least one test.

Run: python3 scripts/probe_exact_mutants.py   (about a minute; not in gate.sh)

Each entry breaks one property the crate promises. A survivor means the
tests do not pin that property; add a test, never delete the mutant.
Anchors must match exactly once, so a refactor that moves code fails loudly
here instead of silently skipping a mutant.
"""
import pathlib, subprocess, sys, os

K = pathlib.Path(__file__).resolve().parents[1]
SRC = K / "crates/algorithms/exact/src"

MUTANTS = [
    ("interval: no outward widening", "interval.rs",
     "            lo: lo.next_down(),\n            hi: hi.next_up(),",
     "            lo,\n            hi,"),
    ("interval: [0,0] not recognised as zero", "interval.rs",
     "        } else if self.lo == 0.0 && self.hi == 0.0 {",
     "        } else if false {"),
    ("interval: square lower bound 0 dropped", "interval.rs",
     "        if self.lo >= 0.0 || self.hi <= 0.0 {\n            return self.mul(self);",
     "        if true {\n            return self.mul(self);"),
    ("dyadic: no normalisation", "dyadic.rs",
     "        Self { mantissa, exponent }.normalised()",
     "        Self { mantissa, exponent }"),
    ("dyadic: subnormal exponent off by one", "dyadic.rs",
     "(fraction, -1074)", "(fraction, -1073)"),
    ("certify: no exact fallback", "certify.rs",
     "    expr.sign_in::<Dyadic>().ok_or(ExactError::Undefined)",
     "    let _ = expr.sign_in::<Dyadic>();\n    Err(ExactError::Undefined)"),
    ("root: dominance sign ignored", "root.rs",
     "    let dominance = a.square().sub(&b.square().mul(c)).sign()?;\n    Some(sign_product(sa, dominance))",
     "    let _dominance = a.square().sub(&b.square().mul(c)).sign()?;\n    Some(sa)"),
    ("root: sqrt(c)=0 case treated as b", "root.rs",
     "    let sb = if radicand == Sign::Zero {\n        Sign::Zero\n    } else {\n        b.sign()?\n    };",
     "    let sb = b.sign()?;"),
    ("two roots: factor 2 dropped", "root.rs",
     "    let irrational = two.mul(p).mul(q);",
     "    let irrational = p.mul(q);\n    let _ = two;"),
    ("cmp_sign: denominator signs ignored", "root.rs",
     "        Some(sign_product(sign_product(inner, d1), d2))",
     "        Some(inner)"),
    ("crossing: denominator sign ignored", "construct.rs",
     "        Some(sign_product(scaled.sign()?, d.sign()?))",
     "        scaled.sign()"),
    ("hit orientation: branch ignored", "construct.rs",
     "        let irrational = self.hit.root_sign::<T>().mul(&k1);",
     "        let irrational = k1;"),
    ("same-circle shortcut arms swapped", "construct.rs",
     "            (Branch::Minus, Branch::Plus) => Sign::Negative,\n            (Branch::Plus, Branch::Minus) => Sign::Positive,",
     "            (Branch::Minus, Branch::Plus) => Sign::Positive,\n            (Branch::Plus, Branch::Minus) => Sign::Negative,"),
    ("tangency: zero discriminant counted as secant", "construct.rs",
     "        Sign::Zero => HitCount::Tangent(hit(Branch::Minus)),",
     "        Sign::Zero => HitCount::Secant(hit(Branch::Minus), hit(Branch::Plus)),"),
    ('interval: sqrt enclosure not widened', 'interval.rs',
     '            lo: self.lo.sqrt().next_down().max(0.0),\n            hi: self.hi.sqrt().next_up(),',
     '            lo: self.lo.sqrt(),\n            hi: self.hi.sqrt(),'),
    ('interval: sqrt of a possibly negative radicand', 'interval.rs',
     '        if self.lo.is_nan() || self.lo < 0.0 {\n            return None;\n        }',
     '        if self.lo.is_nan() {\n            return None;\n        }\n        let this = Self { lo: self.lo.max(0.0), hi: self.hi };\n        let self = &this;'),
    ('tower: product drops b*d*r', 'tower.rs',
     '        let mut out: Vec<T> = ac.iter().zip(&bdr).map(|(p, q)| p.add(q)).collect();',
     '        let _ = &bdr;\n        let mut out: Vec<T> = ac.clone();'),
    ('tower: dominance sign ignored', 'tower.rs',
     '        let dominance = self.sign_at(level - 1, &diff)?;\n        Some(sign_product(sa, dominance))',
     '        let _dominance = self.sign_at(level - 1, &diff)?;\n        Some(sa)'),
    ('tower: negative radicand given a sign', 'tower.rs',
     '        if sr == Sign::Negative {\n            return None;\n        }',
     '        let sr = if sr == Sign::Negative { Sign::Positive } else { sr };'),
    ('tower: zero radicand treated as live', 'tower.rs',
     '        let sb = if sr == Sign::Zero {\n            Sign::Zero\n        } else {\n            self.sign_at(level - 1, b)?\n        };',
     '        let sb = self.sign_at(level - 1, b)?;'),
    ('tower: depth cap removed', 'tower.rs',
     '        if level >= MAX_DEPTH {',
     '        if false {'),
    ('poly: Sturm chain normalised with a sign flip', 'poly.rs',
     'r.scaled_down().coeffs.iter().map(|c| -c).collect()',
     'r.primitive().coeffs.iter().map(|c| -c).collect()'),
    ('poly: multiple roots not reduced', 'poly.rs',
     '        let sf = self.square_free();\n',
     '        let sf = self.primitive();\n'),
    ('poly: open interval may start at a root', 'poly.rs',
     '    let lo_is_root = poly.sign_at(&lo) == Sign::Zero;',
     '    let lo_is_root = false;'),
    ('poly: refinement keeps the wrong half', 'poly.rs',
     '        if s_mid == self.poly.sign_at(&self.lo) {\n            self.lo = mid;',
     '        if s_mid != self.poly.sign_at(&self.lo) {\n            self.lo = mid;'),
    ('poly: equality by gcd skipped', 'poly.rs',
     '            if roots_in(&g.sturm(), &lo, &hi) > 0 {\n                return Sign::Zero;',
     '            if false {\n                return Sign::Zero;'),
    ('poly: cmp_dyadic interior side flipped', 'poly.rs',
     '            // No sign change between lo and x: root is above x.\n            Sign::Positive',
     '            Sign::Negative'),
]

env = {k: v for k, v in os.environ.items()
       if k not in ("CARGO_TARGET_DIR", "CARGO_HOME", "GIT_DIR", "GIT_WORK_TREE")}
env["PATH"] = str(pathlib.Path.home() / ".cargo/bin") + ":" + env["PATH"]

def tests_pass():
    try:
        r = subprocess.run(["cargo", "test", "-q", "-p", "axiolid-exact"], cwd=K, env=env,
                           capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        # A mutant that makes the tests hang (e.g. refining two equal roots
        # forever) is detected, not survived.
        return False
    return r.returncode == 0

assert tests_pass(), "baseline must pass"
survivors = []
for name, fname, old, new in MUTANTS:
    path = SRC / fname
    original = path.read_text()
    assert original.count(old) == 1, f"anchor not unique/missing: {name}"
    path.write_text(original.replace(old, new, 1))
    try:
        killed = not tests_pass()
    finally:
        path.write_text(original)
    print(f"{'killed  ' if killed else 'SURVIVED'}  {name}")
    if not killed:
        survivors.append(name)
assert tests_pass(), "tree restored and passing"
print(f"\n{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
