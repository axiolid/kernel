"""Mutation probe for exact extrusion against the profile normal (#275).

Each mutant drops the mirror, the orientation it restores, or the
in-plane refusal, and must turn a test red: construct's
`tests/downward_extrusion.rs` (every family, oblique directions, mirror
vertices, signed volume, refusal by name) or mesh-compile's (the exact
compiler and the certified boolean deviation).

Equivalent mutants, deliberately not listed: `offset.z > 0.0` as
`offset.z >= 0.0` (an offset within tolerance of the plane is refused
before, so `offset.z` is never zero there), and the derived-profile branch
re-entering through `extrude_profile_exact` with the forward offset (the
offset is already positive, so it takes the forward path again).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/construction/construct/src/extrude_exact.rs"
TESTS = [
    "-p", "axiolid-construct",
    "-p", "axiolid-mesh-compile",
    "--all-features",
    "--test", "downward_extrusion",
]

MUTANTS = [
    ('downward refused again (old forward-only check)', E,
     '    if offset.z.abs() <= tolerance.linear() {',
     '    if offset.z <= tolerance.linear() {', TESTS),
    ('in-plane direction not refused', E,
     '    if offset.z.abs() <= tolerance.linear() {',
     '    if offset.z.abs() < 0.0 {', TESTS),
    ('in-plane refused only above the plane', E,
     '    if offset.z.abs() <= tolerance.linear() {',
     '    if (0.0..=tolerance.linear()).contains(&offset.z) {', TESTS),
    ('in-plane refused under another name', E,
     '        return Err(unsupported("extrusion direction in the profile plane"));',
     '        return Err(unsupported("non-forward planar extrusion"));', TESTS),
    ('mirror skipped: the solid stays above the plane', E,
     '    mirror_in_profile_plane(&forward)\n}',
     '    Ok(forward)\n}', TESTS),
    ('orientation not flipped: builders fed the downward offset', E,
     '    let forward = extrude_forward(profile, Vec3::new(offset.x, offset.y, -offset.z), tolerance)?;',
     '    let forward = extrude_forward(profile, offset, tolerance)?;', TESTS),
    ('half-turn rotation instead of the reflection', E,
     'Mat3::from_diagonal(Vec3::new(1.0, 1.0, -1.0))',
     'Mat3::from_diagonal(Vec3::new(1.0, -1.0, -1.0))', TESTS),
    ('forward and mirrored branches swapped', E,
     '    if offset.z > 0.0 {\n        return extrude_forward(profile, offset, tolerance);',
     '    if offset.z < 0.0 {\n        return extrude_forward(profile, offset, tolerance);', TESTS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
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
