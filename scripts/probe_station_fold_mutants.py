"""Mutation probe for the refusal of stationed sections that fold (#361,
ADR 0082 amendment).

Each mutant breaks one step the decision rests on. The bend's coefficients
(`FrameLaw::bend`): the lateral turning dropped, the up turning dropped, the
plan reading given the section reading's up turning, a reversed piece's
velocity not negated, the spread over the span dropped. The forward floor
(`FrameLawBounds::forward`): returned as zero, a placement's shrink ignored.
The decision (`SectionLaw::fold`): the longitudinal offset's rate dropped, the
witness refusing only far below zero, an undecided span after the last
halving accepted, a span without a bound at its width given up, only a span
with a bound decided exactly, a span at a joint read from the earlier piece,
the exact compiler building a fold it cannot exclude, the plan reading's
vertical part of the velocity taken as all of it, the conservative bound
refusing. Each must turn
`evaluate/tests/station_bend.rs`, the `frame_bound` unit tests or
`compile/tests/station_fold.rs` red.

Not probed, because no test can tell them apart:

- the profile's flattening bound in the clear test and the conservative rate
  (a margin of about the chord budget, far below any distance the tests put
  between a vertex and a fold);
- the `64 eps` rounding allowances (about `1e-14`, below the distance of any
  fold or near-fold the tests build; a touch to within them is the witness's
  case, which the touch tests exercise through the same comparison);
- the middle error's `e` terms (zero on every closed-form frame the tests
  use; an intrinsic frame's is about `1e-12`).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
EVAL = "crates/algorithms/parametric/evaluate"
BEND = f"{EVAL}/src/station/frame_bound/bend.rs"
FB = f"{EVAL}/src/station/frame_bound.rs"
COMPILE = "crates/execution/compile"
FOLD = f"{COMPILE}/src/station/bound/fold.rs"
COMPILE_STATION = f"{COMPILE}/src/station.rs"

BEND_TEST = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_bend"]
UNIT = ["-p", "axiolid-evaluate", "--all-features", "--lib", "frame_bound"]
FOLD_TEST = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_fold"]

MUTANTS = [
    # The bend's coefficients.
    ("the lateral turning is dropped", BEND,
     "            frame.lateral.dot(w) / speed,\n",
     "            0.0 * frame.lateral.dot(w) / speed,\n", [BEND_TEST, FOLD_TEST]),
    ("the up turning is dropped", BEND,
     "            if plan { 0.0 } else { frame.up.dot(w) / speed },\n",
     "            0.0,\n", [BEND_TEST, FOLD_TEST]),
    ("the plan reading takes the section reading's up turning", BEND,
     "            if plan { 0.0 } else { frame.up.dot(w) / speed },\n",
     "            frame.up.dot(w) / speed,\n", [BEND_TEST, FOLD_TEST]),
    ("a reversed piece's velocity is not negated", BEND,
     "    let (v, w) = if piece.reversed { (-v, w) } else { (v, w) };\n",
     "    let (v, w) = (v, w);\n", [BEND_TEST, FOLD_TEST]),
    ("the spread over the span is dropped", BEND,
     "            let pad = middle_error[j] + reach * lip[j];\n",
     "            let pad = middle_error[j];\n", [BEND_TEST, FOLD_TEST]),
    # The forward floor.
    ("the forward floor is zero", FB,
     "        forward: forward * shrink * (1.0 - ROUNDING),\n",
     "        forward: 0.0,\n", [UNIT, BEND_TEST]),
    ("a placement's shrink is ignored", FB,
     "        shrink = placement_shrink(rigid);\n",
     "        shrink = 1.0;\n", [UNIT]),
    # The decision.
    ("the longitudinal offset's rate is dropped", FOLD,
     "        (self.offsets[i + 1].longitudinal - self.offsets[i].longitudinal)\n"
     "            / (self.distances[i + 1] - self.distances[i])\n",
     "        0.0 * (self.offsets[i + 1].longitudinal - self.offsets[i].longitudinal)\n"
     "            / (self.distances[i + 1] - self.distances[i])\n", [FOLD_TEST]),
    ("the witness refuses only far below zero", FOLD,
     "        if floor.is_nan() || floor <= 0.0 {\n            let rate",
     "        if floor.is_nan() || floor <= -1.0 {\n            let rate", [FOLD_TEST]),
    ("an undecided span after the last halving is accepted", FOLD,
     "            if depth >= MAX_FOLD_DEPTH || !(mid > lo && mid < hi) {\n",
     "            if depth >= MAX_FOLD_DEPTH || !(mid > lo && mid < hi) {\n"
     "                continue;\n            }\n            if false {\n", [FOLD_TEST]),
    ("a span without a bound at its width gives up", FOLD,
     "                continue;\n            };\n            let ends",
     "                return Ok(Err(FOLD_NOT_EXCLUDED));\n            };\n            let ends", [FOLD_TEST]),
    ("only a span with a bound is decided exactly", FOLD,
     "        let verdict = if plain && self.fold_frames.has_bend(span.lo, span.hi) {\n",
     "        let verdict = if plain && self.fold_frames.bend(span.lo, span.hi, reading).is_some() {\n", [FOLD_TEST]),
    ("the exact compiler builds a fold it cannot exclude", COMPILE_STATION,
     "    law.fold(margin)?.map_err(unsupported)\n",
     "    let _ = (law.fold(margin)?, unsupported);\n    Ok(())\n", [FOLD_TEST]),
    ("a span at a joint takes the earlier piece", FB,
     "        (holds && mid >= start && mid <= end).then(|| Located {\n",
     "        holds.then(|| Located {\n", [FOLD_TEST]),
    ("the plan reading's vertical part is all of sup |O'|", FOLD,
     "            let speed = bounds.point[0];\n            (speed * speed - bounds.forward * bounds.forward)",
     "            let speed = bounds.point[0];\n            (speed * speed - 0.0 * bounds.forward)", [FOLD_TEST]),
    ("the conservative bound refuses", FOLD,
     "            if depth >= CONSERVATIVE_FOLD_DEPTH {\n                return Ok(Err(reason));\n",
     "            if depth >= CONSERVATIVE_FOLD_DEPTH {\n"
     "                return Err(GeomError::InvalidInput(\"sections fold\".into()));\n",
     [FOLD_TEST]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False)
