"""Mutation probe for ring validation's candidate search and the proof
regions carry (#341, `crates/algorithms/planar/overlay/src/lib.rs`,
`src/settle.rs`, `src/region.rs` and `src/minkowski.rs`).

Each mutant weakens one step -- the boxes the tree is asked with, a
candidate pair dropped or met out of order, the proof `Region::new` sets,
its tolerance, an operand's proof read or ignored, the holes still checked
under a proof, and the proof the overlay's and the disc morphology's
settled outputs carry -- and must turn a test red.

Equivalent mutants, deliberately not listed:
- comparing every pair for longer rings: the same first pair, slower;
- the settled output's proof kept when its repairs ran out of steps: no
  ring of the suite runs them out (64 plus 8 per corner), and the proof is
  then withheld only to be safe;
- a settled ring with a non-finite point proved: operands are finite,
  and a crossing rounds into its edges' boxes, so no settled point is
  not; the check only keeps the proof's statement whole;
- a candidate `j = i + 1`: it is the edge's neighbour, which no pair test
  asks.
"""
import probe_lib

L = "crates/algorithms/planar/overlay/src/lib.rs"
S = "crates/algorithms/planar/overlay/src/settle.rs"
G = "crates/algorithms/planar/overlay/src/region.rs"
M = "crates/algorithms/planar/overlay/src/minkowski.rs"
UNIT = ["-p", "axiolid-overlay", "--lib", "tests::"]
REGION = ["-p", "axiolid-overlay", "--test", "region"]
BOTH = [UNIT, REGION]

MUTANTS = [
    ('the boxes not grown by the tolerance', L,
     '            exact_arc::Bounds::around(&[a, b], linear)',
     '            exact_arc::Bounds::around(&[a, b], 0.0)', UNIT),
    ('the last candidate pair dropped', L,
     '        tree.query_into(|b| b.overlaps(own), &mut candidates, &mut stack);',
     '        tree.query_into(|b| b.overlaps(own), &mut candidates, &mut stack);\n        candidates.pop();',
     UNIT),
    ('a candidate pair dropped from the middle', L,
     '        tree.query_into(|b| b.overlaps(own), &mut candidates, &mut stack);',
     '        tree.query_into(|b| b.overlaps(own), &mut candidates, &mut stack);\n        if candidates.len() > 2 { candidates.remove(candidates.len() / 2); }',
     UNIT),
    ('the last meeting pair reported, not the first', L,
     '        if let Some(&j) = candidates.iter().find(|&&j| j > i && meets(i, j)) {',
     '        if let Some(&j) = candidates.iter().rev().find(|&&j| j > i && meets(i, j)) {',
     UNIT),
    ('Region::new sets no proof', G,
     '        region.validated = Some(validated);',
     '        let _ = validated;', BOTH),
    ('a proof holds for any tolerance', S,
     '        self.linear == tolerance.linear()',
     '        self.linear == self.linear', BOTH),
    ('operands validated whatever their proof', G,
     '            [known(self), known(other)],',
     '            [false, false],', BOTH),
    ('operands trusted without a proof', G,
     '            [known(self), known(other)],',
     '            [true, true],', BOTH),
    ('holes not checked under a proof', L,
     '        if p.holes.iter().any(|h| hole_outside(&p.outer, h, t)) {',
     '        if false {', BOTH),
    ("a boolean's result carries no proof", G,
     '        Ok(Self::from_valid(reported.result.polygons, had_input).with_validated(validated))',
     '        Ok(Self::from_valid(reported.result.polygons, had_input).with_validated(None))',
     BOTH),
    ("a disc morphology's result carries no proof", M,
     '    Ok(Region::from_valid_polygons(polygons, had_input).with_validated(validated))',
     '    Ok(Region::from_valid_polygons(polygons, had_input).with_validated(None))', BOTH),
]

probe_lib.main(MUTANTS, timeout=900)
