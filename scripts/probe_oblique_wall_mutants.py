"""Mutation probe for booleans cut by tilted elliptical-cylinder walls (#287).

An oblique exact extrusion sweeps each arc of its profile into a tilted
elliptical cylinder (#280). The general boolean now evaluates, inverts and
cuts such a wall: a plane meets it in an ellipse (a circle where exactly
one) or in rulings, a wall on a parallel axis in rulings, each in closed
form, and the face is split along a closed-form pcurve (a sinusoid in its
angle, or a vertical line). Each mutant breaks one of those steps and must
turn a test red: evaluate's `elliptical_cylinder.rs` (inversion, jet),
nurbs' `exact_intersection.rs` (the closed forms against both surface
equations), brep-boolean's `oblique_walls.rs` (slabs less oblique shafts and
rounded openings, rulings, parallel walls, pcurve families, closed-form
volumes) and `reflected_walls.rs` (the tests #288 left for #287), and
mesh-compile's `boolean_deviation.rs` (the mesh certified against the exact
difference, with its closed-form volume).

Equivalent mutants, deliberately not listed: the order the two rulings of
a plane are returned in, and the sign of the ellipse frame's second axis
(either orientation names the same set).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/parametric/evaluate/src/surface.rs"
N = "crates/algorithms/parametric/nurbs/src/exact_surface_intersection.rs"
R = "crates/algorithms/parametric/nurbs/src/ruled_section.rs"
S = "crates/algorithms/construction/brep-boolean/src/split.rs"
TESTS = [
    ["-p", "axiolid-evaluate", "--test", "elliptical_cylinder"],
    ["-p", "axiolid-nurbs", "--test", "exact_intersection"],
    ["-p", "axiolid-brep-boolean", "--test", "oblique_walls"],
    ["-p", "axiolid-brep-boolean", "--test", "reflected_walls"],
    ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"],
]

MUTANTS = [
    # Evaluation and inversion.
    ("inversion reads the polar angle, not the affine one", E,
     "let scaled = Vec3::new(local.x / c.semi_axis_x, local.y / c.semi_axis_y, local.z);",
     "let scaled = local;"),
    ("the jet's second u-partial bends the wrong way", E,
     "duu: direct(&c.frame, Vec3::new(-a * co, -b * s, 0.0)),",
     "duu: direct(&c.frame, Vec3::new(-a * co, b * s, 0.0)),"),
    # Plane sections in closed form.
    ("a plane and an elliptical cylinder left to the ruled section", N,
     "        (Surface::EllipticalCylinder(c), Surface::Plane(p))\n"
     "        | (Surface::Plane(p), Surface::EllipticalCylinder(c)) => elliptical_cylinder_plane(c, p),\n",
     ""),
    ("a plane along the axis not cased as rulings", N,
     "    if esign(&nz) == Sign::Zero {\n"
     "        return elliptical_cylinder_plane_rulings(",
     "    if false {\n"
     "        return elliptical_cylinder_plane_rulings("),
    ("the section centre on the wrong side", N,
     "let centre = frame.origin - z * (k.to_f64() / fz);",
     "let centre = frame.origin + z * (k.to_f64() / fz);"),
    ("the second semi-diameter not tilted into the plane", N,
     "let second = (frame.y - z * (fy / fz)) * b;",
     "let second = frame.y * b;"),
    ("the principal axes not turned", N,
     "0.5 * (2.0 * first.dot(second)).atan2(first.dot(first) - second.dot(second))",
     "0.0"),
    ("an exact circle never named", N,
     "let circle = esign(&along) == Sign::Zero && esign(&across_a.sub(&across_b)) == Sign::Zero;",
     "let circle = false;"),
    ("a tangent plane read as missing", N,
     "    let touching = match esign(&gap) {\n"
     "        Sign::Positive => false,\n"
     "        Sign::Zero => true,",
     "    let touching = match esign(&gap) {\n"
     "        Sign::Positive => false,\n"
     "        Sign::Zero => return Err(ExactIntersectionRefusal::Disjoint),"),
    ("the rulings' chord direction mirrored", N,
     "let (ws, wt) = (-my / r, mx / r);",
     "let (ws, wt) = (my / r, mx / r);"),
    # Parallel walls.
    ("parallel walls left unsectioned", R,
     "return parallel_rulings(&c, curve_carrier).map(Some);",
     "return Ok(None);"),
    ("the ruling at u = pi dropped", R,
     "        angles.push(core::f64::consts::PI);\n",
     ""),
    ("a parallel ruling lifted on the wrong semi-axis", R,
     "frame.origin + frame.x * (carrier.x_radius * co) + frame.y * (carrier.y_radius * s);",
     "frame.origin + frame.x * (carrier.x_radius * co) + frame.y * (carrier.x_radius * s);"),
    # Face splitting.
    ("an elliptical wall's sections given traced pcurves", S,
     "                Surface::Cylinder(_) | Surface::EllipticalCylinder(_),\n",
     "                Surface::Cylinder(_),\n"),
    ("the sinusoid's cosine on the wrong semi-axis", S,
     "cosine: -c.semi_axis_x * n.dot(c.frame.x) / nz,",
     "cosine: -c.semi_axis_y * n.dot(c.frame.x) / nz,"),
    ("the sinusoid's mean on the wrong side", S,
     "                mean: n.dot(origin - c.frame.origin) / nz,\n"
     "                cosine: -c.semi_axis_x",
     "                mean: -n.dot(origin - c.frame.origin) / nz,\n"
     "                cosine: -c.semi_axis_x"),
    ("a ruling's pcurve runs down the wall", S,
     "let rise = l.direction.dot(c.frame.z);",
     "let rise = -l.direction.dot(c.frame.z);"),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", "--all-features", *target],
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
        for t in TESTS:
            try:
                code = run(t)
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
