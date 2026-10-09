"""Mutation probe for oblique exact extrusion of arcs (#280) and the mesh
path's refusal of a direction in the profile plane (#281).

Each mutant stands an oblique arc wall upright again, gets the oblique
cylinder's frame, semi-axes or rim graphs wrong, leaves the far rim
unsheared, refuses an oblique circle again or builds it on a right
cylinder, stops refusing an oblique ellipse, or drops, narrows or
renames the mesh path's in-plane refusal. Each must turn a test red in
construct's `tests/oblique_extrusion.rs` (geometric audit,
closed-form volume, sheared vertices, walls that shear back onto the
profile's circles, the oblique cylinder's axis and semi-axes, the
ellipse refusal), `tests/downward_extrusion.rs` (the mirror of the
oblique prisms, the two paths agreeing on the profile plane) or
`tests/exact_extrusion.rs` (the oblique circle's support).

Equivalent mutants, deliberately not listed: the mesh path's `<=`
against the tolerance as `<` (the tests probe half and twice the
tolerance, as #281 asks; the exact path's boundary is pinned by its own
test), and `extrude`'s exact-zero check as `offset.z.abs() == 0.0` (the
same test on an IEEE zero).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/construction/construct/src"
A = f"{C}/extrude_arc.rs"
E = f"{C}/extrude_exact.rs"
M = f"{C}/extrude.rs"
TESTS = [
    "-p", "axiolid-construct",
    "--all-features",
    "--test", "oblique_extrusion",
    "--test", "downward_extrusion",
    "--test", "exact_extrusion",
]

MUTANTS = [
    ('oblique arc wall stands upright again (#280 repro)', A,
     '    if span.shear != Vec2::ZERO {\n        // Leaning',
     '    if false {\n        // Leaning', TESTS),
    ('far rim edges left over the plan section', A,
     '            let shift = if is_top { span.shear } else { Vec2::ZERO };',
     '            let shift = Vec2::ZERO;', TESTS),
    ('oblique section not shortened by the rise', A,
     '            semi_axis_y: radius * rise,',
     '            semi_axis_y: radius,', TESTS),
    ('rim graph leans the wrong way', A,
     '            lean: -radius * lateral,',
     '            lean: radius * lateral,', TESTS),
    ('rim graph without the slant along the axis', A,
     '            slant: length / offset.z,',
     '            slant: 1.0,', TESTS),
    ('angle origin rotated a half turn', A,
     '        let along = Vec3::new(-h.x, -h.y, 0.0);',
     '        let along = Vec3::new(h.x, h.y, 0.0);', TESTS),
    ('arc start angle read from the wrong axes', A,
     '        radial.dot(self.along).atan2(radial.dot(self.across))',
     '        radial.dot(self.across).atan2(radial.dot(self.along))', TESTS),
    ('wall axis along the normal, not the direction', A,
     '        let z = offset / length;',
     '        let z = Vec3::Z;', TESTS),
    ('oblique circle refused again', E,
     '        Some(ObliqueWall::new(Point3::ZERO, circle.radius, offset)?)',
     '        return Err(unsupported("oblique circle extrusion"));', TESTS),
    ('oblique circle wall a right cylinder', E,
     '        Some(wall) => Surface::EllipticalCylinder(wall.surface),',
     '        Some(_) => Surface::Cylinder(Cylinder {\n'
     '            frame: frame_bottom,\n'
     '            radius: circle.radius,\n'
     '        }),', TESTS),
    ('oblique ellipse no longer refused', E,
     '    if offset.x != 0.0 || offset.y != 0.0 {\n'
     '        return Err(unsupported("oblique ellipse extrusion"));',
     '    if false {\n'
     '        return Err(unsupported("oblique ellipse extrusion"));', TESTS),
    ('mesh path builds within tolerance of the plane (#281 repro)', M,
     '    if mesh_offset(direction, depth)?.z.abs() <= tolerance.linear() {',
     '    if false {', TESTS),
    ('mesh path refuses only above the plane', M,
     '    if mesh_offset(direction, depth)?.z.abs() <= tolerance.linear() {',
     '    if (0.0..=tolerance.linear()).contains(&mesh_offset(direction, depth)?.z) {', TESTS),
    ('raw extruder builds a flat solid in the plane', M,
     '    if offset.z == 0.0 {\n        return Err(in_profile_plane());',
     '    if false {\n        return Err(in_profile_plane());', TESTS),
    ('mesh refusal under another name', M,
     '        input: "extrusion direction in the profile plane",',
     '        input: "planar mesh extrusion",', TESTS),
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
