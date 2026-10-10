"""Mutation probe for void shells and cavity results (#120).

Each mutant weakens one decision in how the mesh compiler tessellates a
solid's void shells or how the column builder returns a cavity, and must
turn a test red.

Not listed: the column builder's refusal of a cavity in a result of
several pieces. Two prism operands cannot produce that combination, so no
input reaches it; it guards future callers.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
B = "crates/execution/compile/src/brep.rs"
C = "crates/algorithms/construction/construct/src/column.rs"
TESS = ["-p", "axiolid-mesh-compile", "--test", "brep_tessellation"]
CAV = ["-p", "axiolid-construct", "--test", "column_cavities"]

MUTANTS = [
    ('voids not tessellated', B, '                    check_void_closed(brep, shell)?;\n                    shells.push((shell, true));', '                    check_void_closed(brep, shell)?;', TESS),
    ('open void accepted', B, '    if shell.faces.is_empty() || uses.values().any(|&count| count != 2) {', '    if shell.faces.is_empty() {', TESS),
    ('outward void kept outward', B, '    if volume > 0.0 {\n        for triangle', '    if volume > 0.0 && false {\n        for triangle', TESS),
    ('void flipped into material', B, '    if volume > 0.0 {\n        for triangle', '    if volume < 0.0 {\n        for triangle', TESS),
    ('column drops the cavity', C, '        for void in voids.drain(..) {\n            cavities.push(emit.shell(&void)?);\n        }', '        voids.clear();', CAV),
]

probe_lib.main(MUTANTS, timeout=900)
