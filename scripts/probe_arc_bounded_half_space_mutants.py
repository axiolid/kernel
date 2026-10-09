"""Mutation probe for bounded half-spaces whose boundary has arcs (#277).

Each mutant breaks one decision of the profile boundary -- that the graph
accepts a profile and nothing else, that both compilers read it, check it
and frame it alike, that the exact prism is mirrored rather than reflected
for the opposite side, and each refusal of the boundary check -- and must
turn a test red.

Equivalent mutants, deliberately not listed: mirroring the centred
rectangle, circle and ellipse as derived profiles (they are their own
mirror images, and the derived lowering keeps them exact); checking only
the outer ring of a composite member (the tests' composites are none).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
V = "crates/representations/modeling/graph/src/validation.rs"
M = "crates/execution/compile/src/compiler.rs"
C = "crates/execution/compile/src/exact/clip.rs"
B = "crates/execution/compile/src/half_space_boundary.rs"
TESTS = [
    ["-p", "axiolid-model", "--test", "bounded_half_space_boundary"],
    ["-p", "axiolid-mesh-compile", "--test", "bounded_half_space_arcs"],
    ["-p", "axiolid-mesh-compile", "--test", "exact_half_space_clip"],
]

MUTANTS = [
    ("the graph refuses a profile boundary", V,
     "                return Self::Curve2.accepts(reference, node, nodes)\n"
     "                    || Self::Profile.accepts(reference, node, nodes);",
     "                return Self::Curve2.accepts(reference, node, nodes);"),
    ("the graph accepts any boundary", V,
     "                return Self::Curve2.accepts(reference, node, nodes)\n"
     "                    || Self::Profile.accepts(reference, node, nodes);",
     "                return true;"),
    ("the mesh compiler refuses a profile boundary", M,
     "            GeometryNode::Profile(shape) => {",
     "            GeometryNode::Profile(shape) if false => {"),
    ("the mesh compiler skips the boundary check", M,
     "                crate::half_space_boundary::check(shape, options.tolerance())?;",
     "                let _ = crate::half_space_boundary::check;"),
    ("the mesh compiler chords the boundary coarsely", M,
     "                return profile_rings(shape, chord, options.tolerance());",
     "                return profile_rings(shape, 100.0 * chord, options.tolerance());"),
    ("the exact compiler refuses a profile boundary", C,
     "            Some(GeometryNode::Profile(profile)) => {",
     "            Some(GeometryNode::Profile(profile)) if false => {"),
    ("the exact compiler skips the boundary check", C,
     "                crate::half_space_boundary::check(profile, self.options.tolerance())?;",
     "                let _ = crate::half_space_boundary::check;"),
    ("the profile not mirrored for the opposite side", C,
     "            Some((Footprint::Profile(profile), _)) if flip > 0.0 => (profile.clone(), Vec2::ZERO),",
     "            Some((Footprint::Profile(profile), _)) if true => (profile.clone(), Vec2::ZERO),"),
    ("the mirror is the identity", C,
     "        transform: Transform2::from_scale(Vec2::new(1.0, -1.0)),",
     "        transform: Transform2::from_scale(Vec2::new(1.0, 1.0)),"),
    ("a zero radius let through to lowering", B,
     "            if !(circle.radius.is_finite() && circle.radius > 0.0) {",
     "            if false {"),
    ("closure not checked", B,
     "    let lowered = contour_to_arc_ring(contour, tolerance)?;",
     "    let lowered = contour_to_arc_ring(contour, Tolerance::new(1e9, 1e-9).unwrap())?;"),
    ("a segment of no length let through", B,
     "        if vertex.point == to {",
     "        if false {"),
    ("crossings not checked", B,
     "            uncrossed(&edges, tolerance)",
     "            Ok(())"),
    ("holes not checked", B,
     "            for (ring, contour) in core::iter::once(&contour.outer)\n"
     "                .chain(&contour.holes)",
     "            for (ring, contour) in core::iter::once(&contour.outer)\n"
     "                .chain(&contour.holes[..0])"),
    ("neighbours not excused at their joint", B,
     "                !at_joint",
     "                true"),
    ("the closing joint not a joint", B,
     "        } else if (other.index + 1) % other.count == self.index {",
     "        } else if false {"),
    ("the arc centre on the wrong side", B,
     "        let left = Vec2::new(-chord.y, chord.x) / d;",
     "        let left = Vec2::new(chord.y, -chord.x) / d;"),
    ("the whole circle taken for the arc", B,
     "        chord.perp_dot(p - self.from) * self.bulge.signum() <= 0.0",
     "        true"),
    ("a grazing line missed", B,
     "    let root = discriminant.max(0.0).sqrt();",
     "    let root = discriminant.sqrt();"),
    ("a line-line crossing missed", B,
     "        if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {",
     "        if false {"),
    ("an arc-arc meeting needs only one arc", B,
     "            if on_a(p) && on_b(p) {",
     "            if on_a(p) || on_b(p) {"),
    ("arc ends on the other arc missed", B,
     "    for p in [a.from, a.to] {\n        if on_b(p) {",
     "    for p in [a.from, a.to] {\n        if false && on_b(p) {"),
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
