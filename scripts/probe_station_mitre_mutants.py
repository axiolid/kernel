"""Mutation probe for the bound of cells beside a mitre (#335, ADR 0082
amendment).

Each mutant breaks one step the bound rests on. The exact solid: the
longitudinal offset's shift of the plane dropped, the cell's exact section
replaced by the mesh's mitred row, the outgoing side's `toward` or `normal`
sign flipped. The terms: the in-plane kink dropped (`kink`), the general
cell's gap between the mesh row and the side's projection dropped, the
straight cell's gap dropped, a changing section bounded as a straight one,
each chord end at a mitre taken from the mesh's row. The reach past the
plane: the crossing of `reach_past` dropped. The refusals: the
transversality check disabled, the cell between two mitres bounded as one,
`refinable` returning false (in `within` and in its own definition), a cell
past the budget refined for ever whatever its floor. The refinement's rows:
the plain section at a mitre's knot given to `accept` in place of the mitred
row. Each must turn the unit tests of `station/bound/mitre.rs`,
`compile/tests/station_bound_mitre.rs` or `station_bound_quadrature.rs` red.

Not probed, because no test can tell them apart:

- the mesh row's distance from the plane (`off_plane`): the mitre row lies
  on the plane by construction, up to rounding;
- the read errors of the plane and of the extrusion: the plane's point
  error in the straight cell (`plane_error`), the tangents' (`tangent_error`,
  the `reach * tangent_error` and `extension` terms, the share `cos -
  e_t - e_n`), the other side's (`other_error`), and the read measure and
  the slide in the transversality rate. All are zero on the frames read in
  closed form and about `1e-9` on the two read by quadrature (the deck and
  the B-spline), below the tolerance their values are pinned to; pinning
  them tighter would pin the quadrature's certified error, not the bound;
- the far row's clearance of the plane (`clear`): a vertex that crosses the
  plane at the rate the other check demands also starts clear of it on every
  test cell, and a cell that does not is refused by the rate first;
- the general cell's wall (`wall`, with its `clamp` and `quad`): the cut
  dominates it on every cell the tests build, their quads being planar and
  the faces' `reach_past` reaching at least the clamp.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
MITRE = "crates/execution/compile/src/station/bound/mitre.rs"
BOUND = "crates/execution/compile/src/station/bound.rs"
STATION = "crates/execution/compile/src/station.rs"
UNIT = ["-p", "axiolid-mesh-compile", "--lib", "mitre"]
CELLS = ["-p", "axiolid-mesh-compile", "--test", "station_bound_mitre"]
QUADRATURE = ["-p", "axiolid-mesh-compile", "--test", "station_bound_quadrature"]
TIMEOUT = 1800

MUTANTS = [
    # --- the exact solid ----------------------------------------------------
    ('the longitudinal offset does not shift the plane', MITRE,
     '            q0: base + mitre.normal * long,',
     '            q0: base,', [CELLS]),
    ('the cell\'s exact section is the mesh\'s mitred row', MITRE,
     '            exact: self.exact_at(at, side)?,',
     '            exact: row.to_vec(),', [CELLS]),
    ('the outgoing side\'s toward is not reversed', MITRE,
     '            (-mitre.outgoing.tangent, -mitre.normal, mitre.outgoing.point)',
     '            (mitre.outgoing.tangent, -mitre.normal, mitre.outgoing.point)', [CELLS]),
    ('the outgoing side\'s normal is not reversed', MITRE,
     '            (-mitre.outgoing.tangent, -mitre.normal, mitre.outgoing.point)',
     '            (-mitre.outgoing.tangent, mitre.normal, mitre.outgoing.point)', [CELLS]),
    # --- the terms ------------------------------------------------------------
    ('the kink between the two cuts is dropped', MITRE,
     '        worst = worst.max(height(v0).max(height(v1)) * ratio(v0).max(ratio(v1)));',
     '        worst = worst.max(0.0 * height(v0).max(height(v1)) * ratio(v0).max(ratio(v1)));',
     [UNIT, CELLS]),
    ('the general cell drops the mesh row\'s gap', MITRE,
     '            let gap = cut.gap[k].max(cut.gap[l]);',
     '            let gap = 0.0;', [CELLS]),
    ('the straight cell drops the mesh row\'s gap', MITRE,
     '            let gap = gap[k].max(gap[l]);',
     '            let gap = 0.0;', [CELLS]),
    ('the chord end at a mitre behind is the mesh\'s row', BOUND,
     '            Some(_) => self.exact_at(lo, SeamSide::Outgoing)?,',
     '            Some(_) => at_lo.to_vec(),', [CELLS]),
    ('the chord end at a mitre ahead is the mesh\'s row', BOUND,
     '            Some(_) => self.exact_at(hi, SeamSide::Incoming)?,',
     '            Some(_) => at_hi.to_vec(),', [CELLS]),
    ('a changing section is bounded as a straight one', BOUND,
     '            if straight && !changes && inside.is_empty() {',
     '            if straight && inside.is_empty() {', [CELLS]),
    # --- the reach past the plane ----------------------------------------------
    ('the crossing of the reach is dropped', MITRE,
     '    if (en[0] - eta) * (en[1] - eta) < 0.0 {',
     '    if false {', [UNIT]),
    # --- the refusals ------------------------------------------------------------
    ('the transversality check is disabled', MITRE,
     '            if !positive(rate) || !clear {',
     '            if false {', [CELLS]),
    ('a cell between two mitres is bounded as one', BOUND,
     '                (Some(_), Some(_)) => Ok(Err(MITRE_ON_BOTH_ENDS)),',
     '                (Some(mitre), Some(_)) => self.mitre_cell(&inputs, mitre, SeamSide::Outgoing),',
     [CELLS]),
    ('no refusal is refinable', BOUND,
     '    reason == MITRE_NOT_TRANSVERSAL || reason == MITRE_ON_BOTH_ENDS',
     '    false', [CELLS]),
    ('refinement never splits a refused cell', STATION,
     '        Err(reason) => !bound::refinable(reason),',
     '        Err(_) => true,', [CELLS]),
    ('refinement ignores the floor', STATION,
     '        Ok(cell) => cell.value <= budget || cell.floor > budget,',
     '        Ok(cell) => cell.value <= budget,', [CELLS]),
    # --- the rows refinement bounds --------------------------------------------------
    ('refinement bounds the plain section at the first knot', STATION,
     '                let row_lo = if lo == a { &start_row } else { &at_lo };',
     '                let row_lo = &at_lo;', [CELLS, QUADRATURE]),
    ('refinement bounds the plain section at the last knot', STATION,
     '                let row_hi = if hi == b { &end_row } else { &at_hi };',
     '                let row_hi = &at_hi;', [CELLS, QUADRATURE]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False, reject_uncompiled=True)
