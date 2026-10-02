"""Mutation probe for hole filling, fairing and subdivision in axiolid-refine.

Run: python3 scripts/probe_hole_subdivision_mutants.py   (a few minutes; not in gate.sh)

Each entry weakens one property the crate promises (#129): a refusal, the
minimum-weight triangulation, refinement and relaxation, the fairing
system, or a subdivision or limit mask. A survivor means the tests do not
pin that property; add a test, never delete the mutant. Anchors must match
exactly once, so a refactor that moves code fails loudly here instead of
silently skipping a mutant.
"""
import pathlib, subprocess, sys, os

K = pathlib.Path(__file__).resolve().parents[1]
SRC = K / "crates/algorithms/discrete/refine/src"

MUTANTS = [
    # Hole filling: refusals.
    ("hole: density not validated", "hole.rs",
     "    if !(options.density.is_finite() && options.density > 0.0) {",
     "    if options.density.is_nan() {"),
    ("hole: boundary limit off by one", "hole.rs",
     "    if n > options.max_boundary_vertices {",
     "    if n >= options.max_boundary_vertices {"),
    ("hole: projection simplicity not checked", "hole.rs",
     "    require_simple_projection(&points)?;",
     "    let _ = (&points, require_simple_projection);"),
    ("hole: only proper crossings counted, and as either-or", "hole.rs",
     "    if o1 * o2 < 0 && o3 * o4 < 0 {",
     "    if o1 * o2 < 0 || o3 * o4 < 0 {"),
    ("hole: vertex budget off by one", "hole.rs",
     "            if new_vertices.len() >= options.max_new_vertices {",
     "            if new_vertices.len() > options.max_new_vertices {"),
    # Hole filling: triangulation.
    ("hole: diagonals that already exist allowed", "hole.rs",
     "            if !closing && mesh.find_halfedge(corners[i], corners[k]).is_some() {",
     "            if false {"),
    ("hole: degenerate triangles allowed", "hole.rs",
     "const MIN_TRIANGLE_SINE: Scalar = 1e-10;",
     "const MIN_TRIANGLE_SINE: Scalar = -1.0;"),
    ("hole: dihedral angle ignored (area only)", "hole.rs",
     "            self.cos > other.cos",
     "            self.area < other.area"),
    ("hole: closing edge judged against the wrong outside face", "hole.rs",
     "                above.push((0, Some(outside[n - 1])));",
     "                above.push((0, Some(outside[0])));"),
    ("hole: dihedral with the triangle above a diagonal ignored", "hole.rs",
     "                        cos: below.cos.min(tn.dot(normal)),",
     "                        cos: below.cos.min(if closing { tn.dot(normal) } else { 1.0 }),"),
    # Hole filling: refinement and relaxation.
    ("hole: density factor ignored", "hole.rs",
     "    let alpha = options.density;",
     "    let alpha = 0.5 * options.density;"),
    ("hole: corner scale ignored in the split test", "hole.rs",
     "                d > sigma_c && d > sigma[v.index()]",
     "                d > 0.5 * sigma_c"),
    ("hole: no relaxation", "hole.rs",
     "        relax(mesh, patch);",
     "        let _ = relax;"),
    ("hole: flips edges that are already Delaunay", "hole.rs",
     "    if angle(c) + angle(d) <= PI + FLIP_MARGIN {",
     "    if angle(c) + angle(d) > PI + FLIP_MARGIN {"),
    # Fairing.
    ("fair: harmonic diagonal sign", "fair.rs",
     "                triplets.push((a, a, -diagonal));",
     "                triplets.push((a, a, -0.5 * diagonal));"),
    ("fair: harmonic fixed neighbours subtracted", "fair.rs",
     "                            r[a] += w * p[axis];",
     "                            r[a] -= w * p[axis];"),
    ("fair: biharmonic mass dropped from the matrix", "fair.rs",
     "                        triplets.push((a, b, wa * wb / mass));",
     "                        triplets.push((a, b, wa * wb));"),
    ("fair: biharmonic right-hand side sign", "fair.rs",
     "                        r[a] -= wa * fixed_part[axis] / mass;",
     "                        r[a] += wa * fixed_part[axis] / mass;"),
    ("fair: cotangent at the wrong corner", "fair.rs",
     "    Ok(mesh.target(mesh.next(g)))",
     "    Ok(mesh.target(mesh.prev(g)))"),
    ("fair: unanchored groups accepted", "fair.rs",
     "        if !anchored {",
     "        if false {"),
    ("fair: non-convergence accepted", "fair.rs",
     "        if solution.status != axiolid_numeric::Status::Converged {",
     "        if false {"),
    ("fair: degenerate triangles accepted", "fair.rs",
     "    if cross.is_nan() || cross <= MIN_SINE * u.length() * v.length() {",
     "    if cross.is_nan() {"),
    # Subdivision masks.
    ("loop: beta from the wrong cosine weight", "subdivide.rs",
     "    let c = 3.0 / 8.0 + 0.25 * (2.0 * PI / n).cos();",
     "    let c = 3.0 / 8.0 + 0.125 * (2.0 * PI / n).cos();"),
    ("loop: edge mask is the plain midpoint", "subdivide.rs",
     "            0.375 * (a + b) + 0.125 * (c + d)",
     "            0.5 * (a + b) + 0.0 * (c + d)"),
    ("loop: boundary edges use the interior mask", "subdivide.rs",
     "        positions.push(if mesh.is_boundary_edge(e) {",
     "        positions.push(if false {"),
    ("loop: corner triangles wound backwards", "subdivide.rs",
     "            faces.push([corner(i), mid(i), mid((i + 2) % 3)]);",
     "            faces.push([corner(i), mid((i + 2) % 3), mid(i)]);"),
    ("both: boundary vertex mask", "subdivide.rs",
     "    0.75 * p + 0.125 * (b1 + b2)",
     "    0.5 * p + 0.25 * (b1 + b2)"),
    ("both: boundary neighbour taken twice", "subdivide.rs",
     "    let previous = mesh.source(mesh.prev(h));",
     "    let previous = next;"),
    ("catmull-clark: vertex mask", "subdivide.rs",
     "            (q / n + 2.0 * r / n + (n - 3.0) * p) / n",
     "            (q / n + r / n + (n - 2.0) * p) / n"),
    ("catmull-clark: edge point ignores the faces", "subdivide.rs",
     "                0.25 * (a + b + face_point[f1.index()] + face_point[f2.index()])",
     "                0.5 * (a + b) + 0.0 * (face_point[f1.index()] + face_point[f2.index()])"),
    ("limit: Loop gamma off", "subdivide.rs",
     "                    let gamma = 1.0 / (3.0 / (8.0 * loop_beta(n)) + n);",
     "                    let gamma = 1.0 / (3.0 / (8.0 * loop_beta(n)) + n + 1.0);"),
    ("limit: Catmull-Clark denominator off", "subdivide.rs",
     "                    (n * n * p + 4.0 * edges + diagonals) / (n * (n + 5.0))",
     "                    (n * n * p + 4.0 * edges + diagonals) / (n * (n + 4.0))"),
    ("limit: boundary mask", "subdivide.rs",
     "            (4.0 * p + b1 + b2) / 6.0",
     "            (6.0 * p + b1 + b2) / 8.0"),
    ("subdivide: face budget off by one", "subdivide.rs",
     "        if faces > options.max_faces {",
     "        if faces >= options.max_faces {"),
    ("subdivide: Loop accepts polygons", "subdivide.rs",
     "        if degree != 3 {",
     "        if degree < 3 {"),
]

env = {k: v for k, v in os.environ.items()
       if k not in ("CARGO_HOME", "GIT_DIR", "GIT_WORK_TREE")}
env["PATH"] = str(pathlib.Path.home() / ".cargo/bin") + ":" + env["PATH"]

def tests_pass():
    try:
        r = subprocess.run(["cargo", "test", "-q", "-p", "axiolid-refine",
                            "--test", "hole", "--test", "subdivide"],
                           cwd=K, env=env, capture_output=True, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        # A mutant that makes a test hang is detected, not survived.
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
    print(f"{'killed  ' if killed else 'SURVIVED'}  {name}", flush=True)
    if not killed:
        survivors.append(name)
assert tests_pass(), "tree restored and passing"
print(f"\n{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
