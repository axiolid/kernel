"""Mutation probe for isotropic remeshing (#149, ledger row D7).

Each mutant weakens one decision in
`crates/algorithms/discrete/refine/src/remesh.rs` or its `edits`,
`reference` submodules -- an input refusal, a length threshold, a fold or
shape guard, feature protection, patch-restricted projection, the valence
criterion, the closest-point query -- and
`crates/algorithms/discrete/refine/tests/remesh.rs` (or the closest-point
unit tests in `reference`) must turn a test red.

Deliberately not listed:

- The output budget check. Disabled, the budget test asks for ~10^12
  triangles and the probe would exhaust memory instead of failing.
- The centroid sort's tie-break in the hierarchy. `sort_by` is stable and
  the order it starts from is already ascending, so removing it changes
  nothing.
- Dropping the normal component of the relaxation step, and projecting a
  split's midpoint only in the next relaxation instead of at once. Both are
  followed, within the same iteration, by a projection onto the vertex's
  patch, which removes the difference on every surface the tests use; they
  shape the path, not the result.
- `RemeshReport::max_vertex_deviation`. Every vertex the remesher leaves is
  on the input up to rounding, so no reachable output separates the honest
  measure from a constant zero.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/discrete/refine/src/remesh.rs"
E = "crates/algorithms/discrete/refine/src/remesh/edits.rs"
P = "crates/algorithms/discrete/refine/src/remesh/reference.rs"
T = ["-p", "axiolid-refine", "--lib", "--test", "remesh"]

MUTANTS = [
    (
        "non-finite target accepted",
        R,
        "    if !target.is_finite() || target <= 0.0 {",
        "    if target <= 0.0 {",
    ),
    (
        "feature angle outside [0, pi] accepted",
        R,
        "        if !(0.0..=core::f64::consts::PI).contains(&angle) {",
        "        if false {",
    ),
    (
        "non-finite position accepted",
        R,
        "mesh.positions.iter().position(|p| !p.is_finite())",
        "mesh.positions.iter().position(|_| false)",
    ),
    (
        "degenerate input triangle accepted",
        R,
        "        if n.length() <= DEGENERATE_SINE * scale {",
        "        if false {",
    ),
    (
        "boundary edges not protected",
        R,
        "        protected[e.index()] = if mesh.is_boundary_edge(e) {",
        "        protected[e.index()] = if false {",
    ),
    (
        "sharp edges not detected",
        R,
        "(Some(n0), Some(n1)) => n0.dot(n1).clamp(-1.0, 1.0).acos() > angle,",
        "(Some(n0), Some(n1)) => n0.dot(n1).clamp(-1.0, 1.0).acos() > angle + 10.0,",
    ),
    (
        "features unprotected by default",
        R,
        "            feature_angle: Some(Self::DEFAULT_FEATURE_ANGLE),",
        "            feature_angle: None,",
    ),
    (
        "patches flood across protected edges",
        R,
        "                if protected[e.index()] {\n                    continue;",
        "                if false {\n                    continue;",
    ),
    (
        "reported band counts every short edge",
        R,
        "        if (low..=high).contains(&length) {",
        "        if length <= high {",
    ),
    (
        "projection may land on any patch",
        P,
        "            members[p as usize].push(face as u32);",
        "            for m in members.iter_mut() { let _ = p; m.push(face as u32); }",
    ),
    (
        "hierarchy prunes nodes that may hold the nearest point",
        P,
        "            if box_distance_squared(&node, point) >= limit {",
        "            if box_distance_squared(&node, point) >= limit * 0.01 + 1e-6 {",
    ),
    (
        "closest point inside a triangle swaps its barycentrics",
        P,
        "    a + ab * (vb * denom) + ac * (vc * denom)",
        "    a + ab * (vc * denom) + ac * (vb * denom)",
    ),
    (
        "long edges split only past twice the high threshold",
        E,
        "                .filter(|&e| self.length(e) > self.high)",
        "                .filter(|&e| self.length(e) > 2.0 * self.high)",
    ),
    (
        "short edges collapsed only below half the target",
        E,
        "            low: target * 4.0 / 5.0,",
        "            low: target * 0.5,",
    ),
    (
        "collapse may create long edges",
        E,
        "&& (self.mesh.position(n) - position).length() > self.high",
        "&& (self.mesh.position(n) - position).length() > Scalar::INFINITY",
    ),
    (
        "split edge's second half loses its protection",
        E,
        "        self.set_protected(self.mesh.edge(g), protected);",
        "        self.set_protected(self.mesh.edge(g), false);",
    ),
    (
        "split faces keep no patch",
        E,
        "                self.set_patch(f, p);",
        "                let _ = (f, p);",
    ),
    (
        "merged sides keep protection only if both were protected",
        E,
        "                let merged = self.is_protected(x) || self.is_protected(y);",
        "                let merged = self.is_protected(x) && self.is_protected(y);",
    ),
    (
        "a feature vertex is removed where its line bends",
        E,
        "        incoming.dot(outgoing) > 0.0 && incoming.cross(outgoing).length() <= COLLINEAR * scale",
        "        incoming.dot(outgoing) > -2.0 * scale",
    ),
    (
        "two feature vertices are joined across a patch",
        E,
        "                (false, false) => {}",
        "                (false, false) => candidates.push((a, b, pb)),",
    ),
    (
        "feature vertices treated as free",
        E,
        "        !self.mesh.is_isolated(v) && self.protected_degree(v) == 0",
        "        !self.mesh.is_isolated(v)",
    ),
    (
        "protected edges flipped",
        E,
        "        if self.is_protected(e) || self.mesh.is_boundary_edge(e) {",
        "        if self.mesh.is_boundary_edge(e) {",
    ),
    (
        "no valence flips",
        E,
        "        if after >= before {\n            return false;",
        "        if true {\n            return false;",
    ),
    (
        "boundary valence taken as six",
        E,
        "const BOUNDARY_VALENCE: i64 = 4;",
        "const BOUNDARY_VALENCE: i64 = 6;",
    ),
    (
        "collapse fold check skipped",
        E,
        "                if !replaces_without_folding(before, after) {",
        "                if false {",
    ),
    (
        "flip fold check skipped",
        E,
        "                if !replaces_without_folding(replaced, created) {",
        "                if false {",
    ),
    (
        "vertex move fold check skipped",
        E,
        "            replaces_without_folding(before, after)\n        })",
        "            let _ = (before, after);\n            true\n        })",
    ),
    (
        "edits may leave slivers",
        E,
        "    is_well_shaped(created) && normal(created).dot(normal(replaced)) > 0.0",
        "    normal(created).dot(normal(replaced)) > 0.0",
    ),
    (
        "edits may turn a normal past a right angle",
        E,
        "    is_well_shaped(created) && normal(created).dot(normal(replaced)) > 0.0",
        "    is_well_shaped(created)",
    ),
    (
        "edits may lower the smallest angle without limit",
        E,
        "    after >= before.min(QUALITY_ANGLE)",
        "    let _ = (before, after);\n    true",
    ),
    (
        "no tangential relaxation",
        E,
        "        Some(p + shift - normal * shift.dot(normal))",
        "        Some(p + 0.0 * (shift - normal * shift.dot(normal)))",
    ),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=600,
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
