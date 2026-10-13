"""Mutation probe for authored polygon triangulation (#160, #254).

Each mutant weakens one decision in `compile_authored_polygons` or
`planar::triangulate_polygon` / `triangulate_authored_polygon`;
`tests/authored_polygons.rs` or the `planar` unit tests must fail on it.
A build error also counts as killed.

First run (#160): 9 / 11. Two survivors, neither a test gap:
- "repeated closing corner kept": earcut drops coincident consecutive
  points itself, so the compiler's own dedup was dead code. Removed; the
  repeated-corner test still pins the behaviour end to end.
- "triangle winding not normalised": earcut 0.4 always emits triangles
  counter-clockwise in the projected frame, so the normalising swap never
  fires. Kept as a guard against an earcut convention change; equivalent
  today, so not listed. "triangle winding always flipped" covers the swap.

#254 (warped faces triangulated, their warp reported Certified): 18 / 18,
after adding fixtures for a warped hole and for two warped faces in both
orders, written for the "outer ring only" and "worst across faces"
mutants before the first run. Not listed, as equivalent at test
resolution: dropping the 16-ulp rounding allowance (`Some(worst)`), which
no fixture can resolve.

#261 (the warp is the slab width, max - min signed corner distance from
the fit plane, not the largest corner distance): the #254 anchors moved
into `CornerSpread`; added "largest signed distance only", "largest
corner distance (the #254 bound)" and "lowest corner not tracked".
21 / 21 killed on the first run.

#257 (warped B-rep faces without a surface measured the same way): the
`brep_warped_faces` target and four B-rep mutants added (unmeasured,
reported proven, folded into the exact "plane" contribution, hole
corners dropped). 25 / 25 killed on the first run.

#260 (planar faces through the certified ear clipper of axiolid-construct,
rings touching at a point accepted): the area cross-check and the winding
swap went with earcut -- the clipper's certificate proves the cover and
its triangles are counter-clockwise -- and their three mutants with them.
"holes not passed to earcut" now drops the holes before the clipper, and
the `planar_certified` target and five mutants of `clip_projected` were
added (repeated and closing corners kept, triangles not mapped back past
a dropped corner, pinches refused). 26 / 26 killed on the first run.

#265 (shared edges welded): "faces refuse rings touching at a point" and
"triangles pass through untriangulated" re-anchored on the refactored
`clip_projected` and on the face's rings after the weld; the weld itself
is probed by `probe_shared_edge_weld_mutants.py`.

#266 (a boolean with a warped operand measured against the exact boolean
of its operands flattened to their triangles, the bound widened by the
slab width): sixteen mutants on the composition (the width added, the
budget shared, the detail, the width carried through instances,
booleans' subjects, authored and B-rep operands, the warped operands
named), the dispatch (mesh unions, mesh tools, clips of meshes, the
public compiler still refusing meshes) and the polyhedron (inside-out
meshes), run against the `warped_boolean` target and the polyhedron's
unit tests. Not listed, as equivalent at test resolution: rounding the
widened sum down (one ulp), and the polyhedron's own manifold test,
which the B-rep audit after assembly repeats.

#283 (keyhole faces of a shell): seven mutants re-anchored on the
authored path's `clip_shell_face` and on the warp and Steiner reports
sharing one loop; 7 / 7 still killed. The seam itself is probed by
`probe_keyhole_shell_mutants.py`. Rebased onto #266: its two slab-width
mutants re-anchored where the width now also takes the Steiner offset.
Authored and B-rep faces no longer go through `clip_projected` (now only
curved faces' parameter domains), so "triangles not mapped back past a
dropped corner" and "faces refuse rings touching at a point" survived on
it (40 / 42); they are re-anchored on the shell path the faces take
(`clip_shell_face`'s index map, construct's `triangulate_sealed`):
42 / 42 killed.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/execution/compile/src/planar.rs"
C = "crates/execution/compile/src/compiler.rs"
SEALED = "crates/algorithms/construction/construct/src/ring_triangulation.rs"
B = "crates/execution/compile/src/brep.rs"
DB = "crates/execution/compile/src/deviation/boolean.rs"
D = "crates/execution/compile/src/deviation.rs"
E = "crates/execution/compile/src/exact.rs"
EB = "crates/execution/compile/src/exact/boolean.rs"
EC = "crates/execution/compile/src/exact/clip.rs"
EP = "crates/execution/compile/src/exact/polyhedral.rs"
# #266: the warped boolean's own tests, and the polyhedron's unit tests.
W = [
    ["-p", "axiolid-mesh-compile", "--test", "warped_boolean"],
    ["-p", "axiolid-mesh-compile", "--lib", "exact::polyhedral"],
]
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "authored_polygons"],
    ["-p", "axiolid-mesh-compile", "--lib", "planar"],
    ["-p", "axiolid-mesh-compile", "--test", "brep_warped_faces"],
    ["-p", "axiolid-mesh-compile", "--test", "planar_certified"],
]

MUTANTS = [
    ("planarity check skipped", P,
     "        (self.largest > linear).then(|| self.slab_width())",
     "        (self.largest > linear && false).then(|| self.slab_width())"),
    ("planarity against zero, not tolerance", P,
     "        (self.largest > linear).then(|| self.slab_width())",
     "        (self.largest > 0.0).then(|| self.slab_width())"),
    ("planar callers no longer refuse a warp", P,
     "    triangulate_in_fit_plane(rings, linear, true, false).map(|face| face.indices)",
     "    triangulate_in_fit_plane(rings, linear, false, false).map(|face| face.indices)"),
    ("authored faces refuse a warp again", P,
     "    triangulate_in_fit_plane(rings, linear, false, true)\n",
     "    triangulate_in_fit_plane(rings, linear, true, true)\n"),
    ("warp halved", P,
     "        (self.above - self.below) + 64.0 * f64::EPSILON * self.reach",
     "        0.5 * (self.above - self.below)"),
    ("rounding allowance subtracted", P,
     "        (self.above - self.below) + 64.0 * f64::EPSILON * self.reach",
     "        (self.above - self.below) - 64.0 * f64::EPSILON * self.reach"),
    ("warp measured on the outer ring only (hole corners ignored)", P,
     "        for ring in rings {\n            for &p in *ring {\n                let offset = p - centroid;",
     "        for ring in &rings[..1] {\n            for &p in *ring {\n                let offset = p - centroid;"),
    ("slab: largest signed distance only, not max - min (#261)", P,
     "        (self.above - self.below) + 64.0 * f64::EPSILON * self.reach",
     "        self.above + 64.0 * f64::EPSILON * self.reach"),
    ("slab: largest corner distance, the #254 bound (#261)", P,
     "        (self.above - self.below) + 64.0 * f64::EPSILON * self.reach",
     "        self.largest + 64.0 * f64::EPSILON * self.reach"),
    ("slab: lowest corner not tracked (#261)", P,
     "                spread.below = spread.below.min(distance);\n",
     ""),
    ("B-rep warp unmeasured (#257)", B,
     "    if support.is_none() {\n        report_warp(ctx, &rings, deviation)?;",
     "    if support.is_none() && false {\n        report_warp(ctx, &rings, deviation)?;"),
    ("B-rep warp reported as proven (#257)", B,
     "            DeviationBound::Certified(warp),\n        );\n    }\n    Ok(())",
     "            DeviationBound::Proven(warp),\n        );\n    }\n    Ok(())"),
    ("B-rep warp folded into the exact plane contribution (#257)", B,
     "            crate::deviation::WARPED_BREP_FACE,",
     "            \"plane\","),
    ("B-rep warp measured on the outer ring only (#257)", B,
     "    let views: Vec<&[Vec3]> = points.iter().map(Vec::as_slice).collect();\n    let warp =",
     "    let views: Vec<&[Vec3]> = points.iter().take(1).map(Vec::as_slice).collect();\n    let warp ="),
    ("warp reported as proven", C,
     "                                crate::deviation::DeviationBound::Certified(bound),",
     "                                if detail == crate::deviation::WARPED_AUTHORED_FACE {\n"
     "                                    crate::deviation::DeviationBound::Proven(bound)\n"
     "                                } else {\n"
     "                                    crate::deviation::DeviationBound::Certified(bound)\n"
     "                                },"),
    ("warp folded into the exact authored contribution", C,
     "                        (crate::deviation::WARPED_AUTHORED_FACE, warp),",
     "                        (\"\", warp),"),
    ("warp not reported", C,
     "                        (crate::deviation::WARPED_AUTHORED_FACE, warp),",
     "                        (crate::deviation::WARPED_AUTHORED_FACE, warp.and(None)),"),
    ("worst warp across faces dropped", C,
     "            warp = Some(warp.map_or(face_warp, |w| w.max(face_warp)));",
     "            warp = Some(face_warp);"),
    ("holes not passed to the clipper", P,
     "        clip_shell_face(&flat, &hole_starts, linear, &corners, &projection).map_err(refusal)?;",
     "        clip_shell_face(&flat[..hole_starts.first().copied().unwrap_or(flat.len())], &[], linear, &corners, &projection).map_err(refusal)?;"),
    ("repeated corner kept (#260)", P,
     "            if ring.last() != Some(&point) {",
     "            if true {"),
    ("closing corner kept (#260)", P,
     "        while ring.len() > 1 && ring.first() == ring.last() {",
     "        while false && ring.len() > 1 && ring.first() == ring.last() {"),
    ("triangles not mapped back past a dropped corner (#260)", P,
     "            kept[v]\n",
     "            v\n"),
    ("faces refuse rings touching at a point (#260)", SEALED,
     "        let triangles = without_seams(&points, &rings, PinchPolicy::Accept)?;",
     "        let triangles = without_seams(&points, &rings, PinchPolicy::Refuse)?;"),
    ("triangles pass through untriangulated", C,
     "        if rings.len() == 1 && rings[0].len() == 3 {",
     "        if rings.len() == 1 {"),
    ("hole indices not range-checked", C,
     "        .flat_map(|face| face_rings(face).flatten().copied())",
     "        .flat_map(|face| face.outer.iter().copied())"),
    ("face index in error off by one", P,
     "        _ => GeomError::InvalidInput(format!(\n            \"authored polygon face {face} cannot be triangulated: {refusal}\"\n        )),",
     "        _ => GeomError::InvalidInput(format!(\n            \"authored polygon face {} cannot be triangulated: {refusal}\", face + 1\n        )),"),
    ("flattening not added to the boolean's bound (#266)", DB,
     "        Ok(value) => DeviationBound::Certified(widened(value, flattened)),",
     "        Ok(value) => DeviationBound::Certified(value),", W),
    ("budget not shared with the flattening (#266)", DB,
     "        let left = budget - flattened;",
     "        let left = budget;", W),
    ("flattened boolean not named (#266)", DB,
     "    let detail = if flattened > 0.0 {",
     "    let detail = if flattened < 0.0 {", W),
    ("slab width lost under an instance (#266)", E,
     "                self.set_flattened(root, self.flattened_of(instance.source));\n",
     "", W),
    ("slab width lost through a boolean's subject (#266)", E,
     "                let width = self.flattened_of(left).max(self.flattened_of(right));",
     "                let width = self.flattened_of(right);", W),
    ("authored slab width not recorded (#266)", E,
     "                let width = warp.unwrap_or(0.0).max(steiner_offset.unwrap_or(0.0));",
     "                let width = (0.0 * warp.unwrap_or(0.0)).max(steiner_offset.unwrap_or(0.0));", W),
    ("B-rep slab width not recorded (#266)", D,
     "                ) => Some(widest.max(width)),",
     "                ) => Some(widest + 0.0 * width),", W),
    ("warped authored operand not named under a boolean (#266)", D,
     "                    | (DeviationPath::AuthoredMesh, WARPED_AUTHORED_FACE)\n",
     "", W),
    ("warped B-rep operand not named under a boolean (#266)", D,
     "                    | (DeviationPath::BRepFace, WARPED_BREP_FACE)\n",
     "", W),
    ("mesh unions refused unless both operands are meshes (#266)", EB,
     "        let meshed = subject == Placed::Mesh || tool == Placed::Mesh;",
     "        let meshed = subject == Placed::Mesh && tool == Placed::Mesh;", W),
    ("mesh tool refused (#266)", EB,
     "        if !matches!(tool, Placed::Extrusion | Placed::Mesh) {",
     "        if !matches!(tool, Placed::Extrusion) {", W),
    ("meshes placed outside the report (#266)", EB,
     "                    if self.flatten =>",
     "                    if true =>", W),
    ("mesh clip subject refused (#266)", EC,
     "                if self.flatten =>",
     "                if false =>", W),
    ("public compiler flattens meshes (#266)", E,
     "            flatten: false,",
     "            flatten: true,", W),
    ("report compiles without flattening (#266)", E,
     "        compilation.flatten = true;",
     "        compilation.flatten = false;", W),
    ("inside-out mesh accepted (#266)", EP,
     "    if volume.is_nan() || volume <= 0.0 {",
     "    if volume.is_nan() {", W),
]

probe_lib.main(MUTANTS, tests=TARGETS, timeout=900)
