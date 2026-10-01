"""Mutation probe for the exact 3D Delaunay triangulation (#126).

Each mutant breaks one decision the triangulation rests on -- the
in-sphere orientation, a coefficient of the symbolic perturbation, the
coplanar hull-face circle test, the walk, duplicate detection, the
infinite cells' orientation, the cavity's re-linking, the input range
guard -- or the predicates underneath: the new `in_diametral_sphere`, the
`insphere` filter's error bound (which used to be unsound) and its narrow-grid
integer tier. Each must turn a test red (a hang past the timeout counts).

Not listed, deliberately: ranking the perturbation by ascending instead of
descending lexicographic order. Any fixed total order on the points is a
valid perturbation, so that mutant is equivalent, not a missed fault; the
insertion-order test only demands that one order is used everywhere.
Nor the in-sphere perturbation's final fallback: it is unreachable for
distinct points (two vanishing coefficients would put three collinear
points on one sphere), so changing its answer changes nothing.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/discrete/tetrahedralize/src/delaunay.rs"
S = "crates/algorithms/discrete/tetrahedralize/src/sos.rs"
P = "crates/algorithms/predicates/src/sphere.rs"
T = ["-p", "axiolid-tetrahedralize", "--test", "delaunay3"]
PT = ["-p", "axiolid-predicates", "--test", "diametral_sphere"]
PD = ["-p", "axiolid-predicates", "--test", "delaunay"]
PL = ["-p", "axiolid-predicates", "--lib"]

MUTANTS = [
    ('insphere asked with the wrong orientation', S,
     '    match decided(insphere(a, b, d, c, query)) {',
     '    match decided(insphere(a, b, c, d, query)) {', T),
    ('perturbation: the query point raised moves it inside', S,
     '            // Raising the query point\'s own lift moves it outside.\n            return false;',
     '            return true;', T),
    ('perturbation: vertex coefficient sign flipped', S,
     '        match orientation(replaced[0], replaced[1], replaced[2], replaced[3]) {\n            Sign::Positive => return true,\n            Sign::Negative => return false,',
     '        match orientation(replaced[0], replaced[1], replaced[2], replaced[3]) {\n            Sign::Positive => return false,\n            Sign::Negative => return true,', T),
    ('a point in a hull face plane never conflicts', S,
     '        _ => in_coplanar_circle(points, face, e),',
     '        _ => false,', T),
    ('coplanar circle test inverted', S,
     '    match decided(in_diametral_sphere(a, b, c, query)) {\n        Sign::Positive => return true,\n        Sign::Negative => return false,',
     '    match decided(in_diametral_sphere(a, b, c, query)) {\n        Sign::Positive => return false,\n        Sign::Negative => return true,', T),
    ('coplanar perturbation compares the wrong side', S,
     '            return query_side == side(r, s, q);',
     '            return query_side != side(r, s, q);', T),
    ('walk steps towards the point\'s far side', D,
     '                if self.orientation_with(c, i, p) == Sign::Negative {',
     '                if self.orientation_with(c, i, p) == Sign::Positive {', T),
    ('duplicates are inserted again', D,
     '                .find(|&&v| key(self.point(v)) == key(point))',
     '                .find(|&&v| v == INFINITE && key(self.point(v)) == key(point))', T),
    ('duplicates before the first tetrahedron are kept', D,
     '        if let Some(&existing) = self.seen.get(&key(p)) {',
     '        if let Some(&existing) = self.seen.get(&key(p)).filter(|_| false) {', T),
    ('collinear points accepted into the first tetrahedron', D,
     '            2 => !collinear(self.point(self.frame[0]), self.point(self.frame[1]), p),',
     '            2 => !collinear(self.point(self.frame[0]), self.point(self.frame[1]), p) || true,', T),
    ('infinite cells oriented like their finite neighbour', D,
     '            ghost.swap(others[0], others[1]);',
     '            let _ = &others;', T),
    ('hull triangles turned inside out', D,
     '                // restores it.\n                order.swap(k, 3);\n                if k != 3 {\n                    order.swap(0, 1);',
     '                // restores it.\n                order.swap(k, 3);\n                if k == 3 {\n                    order.swap(0, 1);', T),
    ('cavity faces linked one way only', D,
     '                        self.cells[b as usize].neighbors[l] = here.0;',
     '                        let _ = (b, l);', T),
    ('outside neighbour not re-linked to the new cell', D,
     '            self.cells[outside as usize].neighbors[j] = id;',
     '            let _ = (outside, j);', T),
    ('coordinate range not guarded', D,
     '        if value != 0.0 && !(MIN_COORDINATE..=MAX_COORDINATE).contains(&value.abs()) {',
     '        if false && value != 0.0 && !(MIN_COORDINATE..=MAX_COORDINATE).contains(&value.abs()) {', T),
    ('diametral sphere: exact sign inverted', P,
     '    expansion_sign(&power).flip()',
     '    expansion_sign(&power)', PT),
    ('diametral sphere: filter error bound dropped', P,
     '    let bound = power.error * (1.0 + 2f64.powi(-40)) + tiny;',
     '    let bound = tiny;', PT),
    ('diametral sphere: differences rounded in the exact path', P,
     '    let vector = |p: Point3| [diff(p.x, a.x), diff(p.y, a.y), diff(p.z, a.z)];',
     '    let vector = |p: Point3| [vec![p.x - a.x], vec![p.y - a.y], vec![p.z - a.z]];', PT),
    ('insphere filter bounded by the rounded minors', P,
     '    let permanent = (cd_plus * bz + bd_plus * cz + bc_plus * dz) * alift',
     '    let permanent = (abc.abs() * dlift + dab.abs() * clift) + (cda.abs() * blift + bcd.abs() * alift);\n    let _old = (cd_plus * bz + bd_plus * cz + bc_plus * dz) * alift', PD),
    ('insphere integer tier: sign inverted', P,
     '    Some(match total.signum() {\n        1 => Sign::Positive,\n        -1 => Sign::Negative,',
     '    Some(match total.signum() {\n        1 => Sign::Negative,\n        -1 => Sign::Positive,', PL),
    ('insphere integer tier: taken where i128 overflows', P,
     '    if highest - lowest >= 20 {',
     '    if highest - lowest >= 40 {', PL),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=300,
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
