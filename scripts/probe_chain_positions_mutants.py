"""Mutation probe for the certified positions of arc-length chains (#438,
ADR 0082 amendment).

Each mutant breaks one step the chain table and the bounds that read it rest
on. The table: the position error not carried from one piece to the next,
the turn error not carried, the arc-length enclosure's width dropped from a
parametric read, the evaluator's parameter accepted without checking that
its enclosure holds the distance, a read at a join ignoring which side it
is on. The chain's chord bound: the evaluated ends' gap dropped, the smaller
of the two ends' gaps taken, the upper end read as the piece that starts
there, a zero-width span answered with zero, a Sinusoid piece claimed
certified for flattening. The jets over a chain plan: the grade or the
bend of the profile dropped, the tangent error reported as zero. The
elevated chord bound: the chain arm answered with the exact-chord bound. The
tube: the named refusal of a chain plan restored. Each must turn
`evaluate/tests/chain_position_table.rs`, `evaluate/tests/arc_length_chain.rs`,
`evaluate/tests/elevated_positions.rs`, the solver's unit tests
(`bound::chain_positions`) or `compile/tests/elevated_directrix.rs` red.

Not probed, because no test can tell them apart:

- the share a parametric read may charge before it is refused: the
  enclosures stay under a quarter of the share on the test curves, so a
  looser share changes nothing they can see;
- the `K delta` and `K' delta` terms of the tangent's and the curvature's
  errors, and the rounding allowances of the derivatives (`1e-14`): they
  are at most `1e-11` against references good to `1e-13` at the sizes
  tested;
- the tube's jet charge (`wall_jet_charge`) through a sampled tube: at most
  `1e-11` of a `1e-3` budget, below any sampling. The sum it feeds
  (`tube_error`) and the first-order errors it sums are held by unit tests
  (`deviation::paths::tests`: a coarse table's wall against the exact series
  wall, and against adversarially perturbed jets); that `tube_bound` calls
  `tube_error` at all is only guarded by review.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
EVAL = "crates/algorithms/parametric/evaluate"
CHAIN_POSITIONS = f"{EVAL}/src/bound/chain_positions.rs"
ELEVATED = f"{EVAL}/src/bound/elevated_positions.rs"
BOUND = f"{EVAL}/src/bound.rs"
CHAIN = f"{EVAL}/src/chain.rs"
CURVE = f"{EVAL}/src/curve.rs"
PATHS = "crates/execution/compile/src/deviation/paths.rs"

TABLE = ["-p", "axiolid-evaluate", "--all-features", "--test", "chain_position_table"]
SOLVER = ["-p", "axiolid-evaluate", "--all-features", "--lib", "chain_positions"]
CHORDS = ["-p", "axiolid-evaluate", "--all-features", "--test", "arc_length_chain"]
ELEVATED_TESTS = ["-p", "axiolid-evaluate", "--all-features", "--test", "elevated_positions"]
WALL_CHAIN = [
    "-p",
    "axiolid-mesh-compile",
    "--lib",
    "deviation::paths::tests::a_tube_wall_over_a_harmonic_chain_plan",
]
WALL_COARSE = [
    "-p",
    "axiolid-mesh-compile",
    "--lib",
    "deviation::paths::tests::a_wall_from_a_coarse_chain_table",
]
NONE_BOUND = ["-p", "axiolid-evaluate", "--all-features", "--lib", "claimed_family"]
TUBE = [
    "--release",
    "-p",
    "axiolid-mesh-compile",
    "--test",
    "elevated_directrix",
    "a_disk_along_a_cubic_chain_plan",
]

MUTANTS = [
    # The table.
    ("the position error is not carried to the next piece", CHAIN_POSITIONS,
     "            position = position\n                + local.e",
     "            position = 0.0 * position\n                + local.e", [TABLE]),
    ("the turn error is not carried to the next piece", CHAIN_POSITIONS,
     "            let end_turn = turn + local.phi + PLACE;",
     "            let end_turn = local.phi + PLACE;", [TABLE]),
    ("a parametric read drops the enclosure's width", CHAIN_POSITIONS,
     "                let e = delta + rounding;",
     "                let e = rounding;", [TABLE]),
    ("the evaluator's parameter is accepted unchecked", CHAIN_POSITIONS,
     "    let guess = guess.clamp(start, high);",
     "    let guess = guess.clamp(start, high);\n"
     "    if guess == guess {\n"
     "        return Some((guess, 0.0));\n"
     "    }", [SOLVER, TABLE]),
    ("an incoming read at a join reads the piece that starts there", CHAIN_POSITIONS,
     "        if side == SeamSide::Incoming && index > 0 && s == self.begins[index] {",
     "        if false && side == SeamSide::Incoming && index > 0 && s == self.begins[index] {",
     [TABLE, CHORDS]),
    # The chain's chord bound.
    ("the chord bound drops the evaluated ends' gap", CHAIN_POSITIONS,
     "    let out = inflate(Iv::point(law).add(Iv::point(g)).hi);",
     "    let out = inflate(Iv::point(law).add(Iv::point(0.0 * g)).hi);", [CHORDS]),
    ("the chord bound takes the smaller end gap", CHAIN_POSITIONS,
     "    let g = gap(lo, SeamSide::Outgoing)?.max(gap(hi, SeamSide::Incoming)?);",
     "    let g = gap(lo, SeamSide::Outgoing)?.min(gap(hi, SeamSide::Incoming)?);", [CHORDS]),
    ("the upper end is read as the piece that starts there", CHAIN_POSITIONS,
     "    let g = gap(lo, SeamSide::Outgoing)?.max(gap(hi, SeamSide::Incoming)?);",
     "    let g = gap(lo, SeamSide::Outgoing)?.max(gap(hi, SeamSide::Outgoing)?);", [CHORDS]),
    ("a zero-width span of a chain has no gap", BOUND,
     "            Curve2::Lifted(_) | Curve2::Intrinsic(_) | Curve2::Chain(_)",
     "            Curve2::Lifted(_) | Curve2::Intrinsic(_)", [CHORDS]),
    ("a sinusoid piece is claimed certified for flattening", CHAIN,
     "                ) && crate::bound::certifies_flattening2(curve)",
     "                ) || crate::bound::certifies_flattening2(curve)", [CHORDS]),
    ("a read's rounding grows with the square of the stretch", CHAIN_POSITIONS,
     "PLACE * (origin + self.stretch * p.length().max(1.0))",
     "PLACE * (origin + self.stretch * w.length().max(1.0))", [CHORDS]),
    ("a claimed family passes a span without a bound", CURVE,
     "        None => !claimed,",
     "        None => true,", [NONE_BOUND, CHORDS]),
    # Elevated and banked curves over a chain plan.
    ("the elevated chord bound answers a chain plan with the exact chord", ELEVATED,
     "        Curve2::Chain(chain) => (chain.start, chain.length()?),",
     "        Curve2::Chain(_) => return Some(law),", [ELEVATED_TESTS]),
    ("the chain jet drops the profile's grade", ELEVATED,
     "                    Vec3::new(tangent.x, tangent.y, grade),",
     "                    Vec3::new(tangent.x, tangent.y, 0.0),", [ELEVATED_TESTS]),
    ("the chain jet drops the profile's bend", ELEVATED,
     "                    Vec3::new(curvature.x, curvature.y, bend),",
     "                    Vec3::new(curvature.x, curvature.y, 0.0),", [ELEVATED_TESTS]),
    ("the chain jet reports no tangent error", ELEVATED,
     "                table.max_tangent_error(),",
     "                0.0,", [ELEVATED_TESTS]),
    # The tube.
    ("the tube's centre reverts to the evaluator's", PATHS,
     "                table.jet3_on(s, side)\n",
     "                table.jet3_on(s, side).map(|(_, c1, c2)| {\n"
     "                    (axiolid_reference::curve::evaluate3(self.curve, s).unwrap(), c1, c2)\n"
     "                })\n", [WALL_CHAIN]),
    ("the tube carries no jet charge", PATHS,
     "    centre + jets.max(caps)",
     "    centre + 0.0 * jets.max(caps)", [WALL_COARSE]),
    ("the tube carries no end-disk tilt", PATHS,
     "    centre + jets.max(caps)",
     "    centre + jets", [WALL_COARSE]),
    ("a chain plan is refused by name again", PATHS,
     "    has_intrinsic_profile(&base.elevation).then_some(UNCERTIFIED_DIRECTRIX_POINTS)",
     "    (matches!(base.plan.as_ref(), Curve2::Chain(_)) || has_intrinsic_profile(&base.elevation))\n"
     "        .then_some(UNCERTIFIED_DIRECTRIX_POINTS)", [TUBE]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
