"""Mutation probe for keyhole faces of closed shells (#283).

Each mutant weakens one decision that keeps a keyhole face's seam out of
its triangle edges: construct's `ring_triangulation/seal.rs` (the flip,
the Steiner point, its margin and the certificate), the seam ends it reads
(`seam.rs`), the sealed entry point (`ring_triangulation.rs`), and
mesh-compile's placement of Steiner points in space and their report
(`planar.rs`, `brep.rs`, `compiler.rs`), and the issue's related
cleanup, `pinch::find` welding `-0.0` with `0.0`. Construct's unit tests and
`profile_pinches`, or mesh-compile's `keyhole_shell` and
`planar_projected_rings`, must fail on it. Every mutant compiles.

Not listed, as equivalent at test resolution:
- dropping the test that a Steiner point lies off the seam's line: a
  candidate is the seam's midpoint moved towards a corner off that line,
  so it leaves the line only once it has rounded into the midpoint, some
  fifty halvings in, after a candidate inside has been found.
- skipping the re-check of triangles at a Steiner point once it is placed
  in space (`steiner_turns` on the clipper's cover): the extra halving
  keeps the point inside by half the seam midpoint's clearance, far beyond
  the projection's rounding, so no fixture reaches it; it guards rounding
  the margin does not cover.
- leaving a Steiner offset out of an authored operand's flattened width
  (`exact.rs`): the offset is a rounding-level slab width (6e-15 for the
  tilted rim), far below the boolean measurement's own resolution; the
  offset is still carried and named under the boolean.
- letting the noise-band split put a seam edge back, or split invented
  edges at Steiner points: the split only re-cuts triangles along corners
  within a micrometre of an invented edge, and no keyhole fixture has
  such a corner near its seam.

First run: 22 / 22 killed. The margin mutant ("the first fitting
candidate taken, no margin") is caught by the wide keyhole box as
authored polygons, whose first fitting candidate is collinear in space
with two corners and only inside by the projection's rounding.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/construction/construct/src/ring_triangulation"
SEAL = D + "/seal.rs"
SEAM = D + "/seam.rs"
RT = D + ".rs"
P = "crates/execution/compile/src/planar.rs"
B = "crates/execution/compile/src/brep.rs"
C = "crates/execution/compile/src/compiler.rs"
PI = "crates/execution/compile/src/pinch.rs"
TARGETS = [
    ["-p", "axiolid-construct", "--lib", "ring_triangulation"],
    ["-p", "axiolid-construct", "--test", "profile_pinches"],
    ["-p", "axiolid-mesh-compile", "--test", "keyhole_shell"],
    ["-p", "axiolid-mesh-compile", "--test", "planar_projected_rings"],
]

MUTANTS = [
    # Construct: taking the seam out.
    ("the sealed path seals nothing", RT,
     "    let steiner = seal::seal(&mut points, &mut triangles, &ends)?;",
     "    let steiner: Vec<Steiner> = Vec::new();\n    let _ = &ends;"),
    ("seam ends read one edge late", SEAM,
     "            [at(seam.open), at(seam.open + 1)]",
     "            [at(seam.open + 1), at(seam.open + 2)]"),
    ("every seam passed over", SEAL,
     "    for &[p, q] in seams {",
     "    for &[p, q] in &seams[..0] {"),
    ("flipped without the convexity test", SEAL,
     "            if orient(at(c), at(p), at(d)) > 0 && orient(at(d), at(q), at(c)) > 0 {",
     "            if orient(at(c), at(p), at(d)) > 0 || orient(at(d), at(q), at(c)) > 0 {"),
    ("never flipped, always a Steiner point", SEAL,
     "            if orient(at(c), at(p), at(d)) > 0 && orient(at(d), at(q), at(c)) > 0 {",
     "            if false && orient(at(c), at(p), at(d)) > 0 && orient(at(d), at(q), at(c)) > 0 {"),
    ("the Steiner point need not lie inside", SEAL,
     "            .all(|(a, b)| orient(at(a), at(b), s) > 0)",
     "            .all(|(a, b)| orient(at(a), at(b), s) > 0 || true)"),
    ("the first fitting candidate taken, no margin", SEAL,
     "                return Some(if fits(&deeper) { deeper } else { first });",
     "                return Some(if false && fits(&deeper) { deeper } else { first });"),
    ("too few halvings to reach the seam", SEAL,
     "pub(crate) const MAX_HALVINGS: u32 = 60;",
     "pub(crate) const MAX_HALVINGS: u32 = 2;"),
    ("the chain halves towards the seam's end, not its midpoint", SEAL,
     "        let m = midpoint(at(self.seam[0]), at(self.seam[1]));",
     "        let m = midpoint(at(self.seam[0]), at(self.seam[0]));"),
    # Construct: the certificate.
    ("certificate ignores the triangle count", SEAL,
     "    if after.len() != before.len() + 2 * steiner {",
     "    if after.len() != after.len() {"),
    ("certificate ignores the winding", SEAL,
     "        if orient(corner(0), corner(1), corner(2)) <= 0 {",
     "        if orient(corner(0), corner(1), corner(2)) < -1 {"),
    ("certificate ignores an edge used twice", SEAL,
     "            if *uses > 1 {",
     "            if *uses > 2 {"),
    ("certificate ignores a seam edge", SEAL,
     "                .any(|s| (s[0], s[1]) == edge || (s[1], s[0]) == edge)",
     "                .any(|s| (s[0], s[1]) == edge && (s[1], s[0]) == edge)"),
    ("certificate ignores the boundary", SEAL,
     "    if net(before) != net(after) {",
     "    if net(before) != net(before) {"),
    # mesh-compile: placing and reporting Steiner points.
    ("Steiner point placed at a corner in space", P,
     "            |a, b| (a + b) * 0.5,",
     "            |a, _b| a,"),
    ("authored faces triangulated as surface patches", P,
     "    triangulate_in_fit_plane(rings, linear, false, true)",
     "    triangulate_in_fit_plane(rings, linear, false, false)"),
    ("axis-plane test passes every face", P,
     "    if (0..3).any(|k| points().all(|&p| coordinate(p, k) == coordinate(first, k))) {",
     "    if (0..3).any(|k| points().all(|&p| coordinate(p, k) == coordinate(first, k)) || true) {"),
    ("axis-plane test passes no face", P,
     "    if (0..3).any(|k| points().all(|&p| coordinate(p, k) == coordinate(first, k))) {",
     "    if (0..3).any(|k| points().all(|&p| coordinate(p, k) == coordinate(first, k)) && false) {"),
    ("B-rep Steiner offset not reported", B,
     "    report_steiner(&rings, &cover.steiner, deviation)?;",
     "    let _ = report_steiner;"),
    ("authored Steiner offset not reported", C,
     "                        (crate::deviation::KEYHOLE_STEINER_POINT, steiner_offset),",
     "                        (crate::deviation::KEYHOLE_STEINER_POINT, steiner_offset.and(None)),"),
    ("B-rep Steiner point indexed past its position", B,
     "        local.push(mesh.positions.len() as u32);\n        mesh.positions.push(point);",
     "        mesh.positions.push(point);\n        local.push(mesh.positions.len() as u32);"),
    # mesh-compile: the flattened-exact measurement of a boolean (#266).
    ("Steiner offset not carried through a boolean", "crates/execution/compile/src/deviation.rs",
     "                    | (DeviationPath::AuthoredMesh, KEYHOLE_STEINER_POINT)\n",
     ""),
    ("B-rep with Steiner points not a polyhedron", "crates/execution/compile/src/deviation.rs",
     "                    WARPED_BREP_FACE | KEYHOLE_STEINER_POINT,",
     "                    WARPED_BREP_FACE,"),
    # Related cleanup: mesh-compile's pinch search welds by value.
    ("pinch search welds by raw bits", PI,
     "            (p.x + 0.0).to_bits(),",
     "            p.x.to_bits(),",
     ["-p", "axiolid-mesh-compile", "--lib", "pinch"]),
]

probe_lib.main(MUTANTS, tests=TARGETS, timeout=1800, reject_uncompiled=True)
