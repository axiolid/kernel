"""Mutation probe for welding shared edges and checking claimed solids (#265).

Each mutant weakens one decision of the weld (`weld.rs`, the pre-passes in
`brep.rs` and `compiler.rs`, `planar::edge_touches`, construct's
`profile::ring_touches`) or of the closure check (`channels::checked_solid`,
its B-rep and boolean callers, `combined_closure`, `solid_mesh`);
`tests/shared_edge_weld.rs` or construct's `profile_pinches` must fail on
it. A build error also counts as killed.

The weld mutants leave the closure check in place, so a reproducer that
stops welding is still caught: its outcome becomes `OpenSolid`, not
`Solid`. The check mutants leave the weld in place and are caught by the
fixtures that cannot be welded (a missing neighbour, a face turned round,
an open boolean result).

Not listed, as equivalent at test resolution:
- `Built::leaf` unchecked: every operation this crate builds a solid with
  closes by construction (certified caps, primitives), so no fixture can
  open one; the check guards a future regression.
- the de-duplication in `EdgeSplits::record`: `ring_touches` lists each
  vertex and edge once, and a corner recorded twice on one edge would be
  inserted twice in a row, which the clipper drops as a repeated corner.

First run: 16 / 17. The survivor, "touch corner not mapped past a dropped
corner", needed a face whose clipper input differs from its corners: the
fixture with an exporter's closing corner on the touching face was added
for it. 17 / 17 killed after.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
W = "crates/execution/compile/src/weld.rs"
B = "crates/execution/compile/src/brep.rs"
C = "crates/execution/compile/src/compiler.rs"
P = "crates/execution/compile/src/planar.rs"
CH = "crates/execution/compile/src/channels.rs"
BO = "crates/execution/compile/src/channels/boolean.rs"
O = "crates/contracts/operations/compile/src/outcome.rs"
R = "crates/algorithms/construction/construct/src/ring_triangulation.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "shared_edge_weld"],
    ["-p", "axiolid-construct", "--test", "profile_pinches"],
]

MUTANTS = [
    # The weld.
    ("B-rep faces do not take recorded corners", B,
     "    ctx.splits\n        .split(&mut rings, |corner| corner, |vertex, point| (vertex, point));\n",
     ""),
    ("authored faces do not take recorded corners", C,
     "            splits.split(&mut rings, |i| (i, mesh.positions[i as usize]), |i, _| i);\n",
     ""),
    ("B-rep pre-pass records nothing", B,
     "            splits.record(&corners, &touches);",
     "            let _ = (&corners, &touches);"),
    ("authored pre-pass records nothing", C,
     "            splits.record(&corners, &touches);",
     "            let _ = (&corners, &touches);"),
    ("B-rep pre-pass skips faces with holes", B,
     "        if rings.len() == 1 && rings[0].len() == 3 {\n            continue;\n        }",
     "        if rings.len() > 1 || rings[0].len() == 3 {\n            continue;\n        }"),
    ("authored triangle faces keep their corner order unsplit", C,
     "            if rings.len() == 1 && rings[0].len() == 3 {",
     "            if face.outer.len() == 3 && face.holes.is_empty() {"),
    ("corners inserted in reverse order along the edge", W,
     "                along.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));",
     "                along.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));"),
    ("touch edge read as its start corner twice", P,
     "                to: kept[start + (touch.edge + 1) % len],",
     "                to: kept[start + touch.edge],"),
    ("touch corner not mapped past a dropped corner", P,
     "                corner: kept[touch.vertex],",
     "                corner: touch.vertex,"),
    ("ring_touches returns nothing", R,
     "    Ok(listed)",
     "    Ok(Vec::new())"),
    ("ring_touches lists a vertex once per meeting edge", R,
     "        if !listed.contains(&entry) {",
     "        if true {"),
    # The closure check.
    ("a claimed solid is never checked", CH,
     "    if adjacency.edge_count() > 0 && adjacency.is_closed_two_manifold() {",
     "    if true {"),
    ("closure check ignores orientation", CH,
     "    if adjacency.edge_count() > 0 && adjacency.is_closed_two_manifold() {",
     "    if adjacency.edge_count() > 0 && adjacency.edges().all(|(_, uses)| uses.len() == 2) {"),
    ("B-rep keeps its declared closure", B,
     "    if declared == MeshClosure::Solid {\n        crate::channels::checked_solid(mesh)",
     "    if false {\n        crate::channels::checked_solid(mesh)"),
    ("boolean result labelled solid unchecked", BO,
     "    let closure = super::checked_solid(&mesh);",
     "    let closure = axiolid_mesh_compile_contract::MeshClosure::Solid;"),
    ("an open solid reads as a surface in a collection", CH,
     "            (MeshClosure::OpenSolid, _) | (_, MeshClosure::OpenSolid) => MeshClosure::OpenSolid,",
     "            (MeshClosure::OpenSolid, _) | (_, MeshClosure::OpenSolid) => MeshClosure::Surface,"),
    ("solid_mesh hands out an open solid", O,
     "            MeshClosure::OpenSolid => Err(GeomError::InvalidInput(",
     "            MeshClosure::OpenSolid => Ok(&self.mesh),\n            #[allow(unreachable_patterns)]\n            MeshClosure::OpenSolid => Err(GeomError::InvalidInput("),
]


def run():
    for target in TARGETS:
        code = subprocess.run(["cargo", "test", "-q", "--all-features", *target], cwd=ROOT,
                              capture_output=True, text=True, timeout=1800).returncode
        if code != 0:
            return code
    return 0


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found ({original.count(old)})"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run()
        except subprocess.TimeoutExpired:
            code = -1
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)

print(f"{len(MUTANTS) - len(survivors)} / {len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
