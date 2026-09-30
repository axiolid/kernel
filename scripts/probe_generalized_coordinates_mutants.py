"""Mutation probe for Wachspress, discrete harmonic and 3D mean-value
coordinates (#143, ledger row H5).

Each mutant weakens one decision in
`crates/algorithms/query/spatial/src/barycentric.rs` -- a convexity
threshold, a formula index, a sign, a refusal -- and
`crates/algorithms/query/spatial/tests/barycentric.rs` must turn a test red.

Equivalent mutant, deliberately not listed: replacing
`ui.dot(uj.cross(uk)).signum()` with the constant `1.0` in
`mean_value_coordinates3`. Every mesh this crate's tests use (tetrahedron,
cube, octahedron) is convex with consistently outward-wound faces, so from
any interior point every face is back-facing and the sign is the same
constant for every face already; the mutant only changes behaviour for a
non-convex mesh where some faces face the query point and others do not,
which is a materially larger fixture than this probe's scope.

Also deliberately not listed: swapping `sin_theta[1] * s[2]` for
`sin_theta[2] * s[1]` in `weights[fi]`'s own denominator. Floater, Kos and
Reimers' spherical law of sines for the fan triangle `(i, j, k)` makes
`sin(theta_j) * s_k == sin(theta_k) * s_j` an identity, not a coincidence of
the fixtures tried, so the two denominators agree everywhere the formula is
defined.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/algorithms/query/spatial/src/barycentric.rs"
T = ["-p", "axiolid-spatial", "--test", "barycentric"]

MUTANTS = [
    (
        "convexity threshold direction inverted",
        B,
        "        if turn * orientation < -threshold {",
        "        if turn * orientation > -threshold {",
    ),
    (
        "convexity check skipped entirely",
        B,
        "    check_convex(v, tolerance)?;\n    if let Some(w) = vertex_or_edge_weights(v, point, tolerance) {\n        return Ok(w);\n    }\n\n    let a: Vec<Scalar> = (0..n).map(|i| area2(v[i], v[(i + 1) % n], point)).collect();",
        "    if let Some(w) = vertex_or_edge_weights(v, point, tolerance) {\n        return Ok(w);\n    }\n\n    let a: Vec<Scalar> = (0..n).map(|i| area2(v[i], v[(i + 1) % n], point)).collect();",
    ),
    (
        "wachspress corner triangle uses the wrong neighbour",
        B,
        "        let c = area2(v[prev], v[i], v[(i + 1) % n]);",
        "        let c = area2(v[i], v[(i + 1) % n], v[(i + 2) % n]);",
    ),
    (
        "wachspress denominator pairs the wrong areas",
        B,
        "        let w = c / (a[prev] * a[i]);",
        "        let w = c / (a[i] * a[i]);",
    ),
    (
        "discrete harmonic drops the unsigned-angle correction",
        B,
        "        let cot_prev = u.dot(prev) / u.perp_dot(prev).abs();\n        let cot_next = u.dot(next) / u.perp_dot(next).abs();",
        "        let cot_prev = u.dot(prev) / u.perp_dot(prev);\n        let cot_next = u.dot(next) / u.perp_dot(next);",
    ),
    (
        "discrete harmonic divides by distance instead of squared distance",
        B,
        "        let r2 = u.length_squared();",
        "        let r2 = u.length();",
    ),
    (
        "mean-value-3 vertex snap removed",
        B,
        "    if let Some(i) = (0..n)\n        .filter(|&i| d[i] <= linear)\n        .min_by(|&i, &j| d[i].total_cmp(&d[j]))\n    {\n        let mut weights = vec![0.0; n];\n        weights[i] = 1.0;\n        return Ok(weights);\n    }",
        "",
    ),
    (
        "mean-value-3 out-of-range face index not refused",
        B,
        "            if index >= n {",
        "            if index > n {",
    ),
    (
        "mean-value-3 weight update names the wrong distance",
        B,
        "        weights[fi] +=\n            (theta[0] - c[1] * theta[2] - c[2] * theta[1]) / (d[fi] * sin_theta[1] * s[2]);",
        "        weights[fi] +=\n            (theta[0] - c[1] * theta[2] - c[2] * theta[1]) / (d[fj] * sin_theta[1] * s[2]);",
    ),
    (
        "mean-value-3 coplanar-inside shortcut skipped",
        B,
        "        if core::f64::consts::PI - h < linear {",
        "        if false {",
    ),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(T)
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
