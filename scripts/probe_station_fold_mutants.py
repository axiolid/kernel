"""Mutation probe for the refusal of stationed sections that fold (#361,
ADR 0082 amendment).

Each mutant breaks one step the decision rests on. The bend's coefficients
(`FrameLaw::rates`, of which `bend` is the projection since #406): the lateral
turning dropped, the up turning dropped, the plan reading given the section
reading's up turning, a reversed piece's velocity not negated, the spread over
the span dropped. The forward floor (`FrameLawBounds::forward`): returned as
zero, a placement's shrink ignored. The decision (`SectionLaw::decide`): the
longitudinal offset's rate dropped, the witness refusing only far below zero,
an undecided span after the last halving accepted, a span without a bound at
its width given up, only a span with rates decided exactly, a span at a joint
read from the earlier piece, the exact compiler building a turned fold, the
plan reading's vertical part of the velocity taken as all of it, the
conservative bound refusing. Each must turn
`evaluate/tests/station_bend.rs`, the `frame_bound` unit tests or
`compile/tests/station_fold.rs` red. The probe of the turned and quadrature
decision itself is `probe_exact_fold_mutants.py`.

The sheets of open sections (#407) share the decision and add three steps in
`station::sectioned_surface`: a sectioned surface and an open-section run
skipping the fold check (each by its own test), and an undecided fold
reported certified (a sheet along the offset of a cubic, whose spans across
the knot the conservative bound cannot clear). Each must turn
`compile/tests/station_fold_sheet.rs` red.

Not probed, because no test can tell them apart:

- the `within` flag of a sheet (`fold.is_ok()`): refinement does not change
  when a run is unbounded, so no test can tell it from `true`;
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
RATES = f"{EVAL}/src/station/frame_bound/rates.rs"
FB = f"{EVAL}/src/station/frame_bound.rs"
COMPILE = "crates/execution/compile"
FOLD = f"{COMPILE}/src/station/bound/fold.rs"
RATE = f"{COMPILE}/src/station/bound/fold/rate.rs"
COMPILE_STATION = f"{COMPILE}/src/station.rs"

BEND_TEST = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_bend"]
UNIT = ["-p", "axiolid-evaluate", "--all-features", "--lib", "frame_bound"]
FOLD_TEST = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_fold"]
SHEET = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_fold_sheet"]
TURNED = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_fold_turned"]
BOUND_TEST = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_bound"]

SHEET_CHECK = (COMPILE_STATION,
               "    let fold = law.fold(0.0)?;\n",
               "    let fold: Result<(), &'static str> = Ok(());\n")

MUTANTS = [
    # The bend's coefficients.
    ("the lateral turning is dropped", RATES,
     "        let lateral = RateBounds::around(frame.lateral.dot(w) / speed, turn, reach, lip[1]);\n",
     "        let lateral = RateBounds::around(0.0 * frame.lateral.dot(w) / speed, turn, reach, lip[1]);\n",
     [BEND_TEST, FOLD_TEST]),
    ("the up turning is dropped", RATES,
     "            RateBounds::around(frame.up.dot(w) / speed, turn, reach, lip[2])\n",
     "            RateBounds::around(0.0, turn, reach, lip[2])\n", [BEND_TEST, FOLD_TEST]),
    ("the plan reading takes the section reading's up turning", RATES,
     "        let up = if plan {\n            RateBounds::ZERO\n",
     "        let up = if false {\n            RateBounds::ZERO\n", [BEND_TEST, FOLD_TEST]),
    ("a reversed piece's velocity is not negated", BEND,
     "    let (v, w) = if piece.reversed { (-v, w) } else { (v, w) };\n",
     "    let (v, w) = (v, w);\n", [BEND_TEST, FOLD_TEST]),
    ("the spread over the span is dropped", RATES,
     "        let pad = error + reach * lipschitz;\n",
     "        let pad = error;\n", [BEND_TEST, FOLD_TEST]),
    # The forward floor.
    ("the forward floor is zero", FB,
     "        forward: forward * shrink * (1.0 - ROUNDING),\n",
     "        forward: 0.0,\n", [UNIT, BEND_TEST]),
    ("a placement's shrink is ignored", FB,
     "        shrink = placement_shrink(rigid);\n",
     "        shrink = 1.0;\n", [UNIT]),
    # The decision.
    ("the longitudinal offset's rate is dropped", RATE,
     "        let rates = [0, 1, 2].map(|k| (components(b)[k] - components(a)[k]) / delta);\n",
     "        let rates = [0, 1, 2].map(|k| if k == 0 { 0.0 } else { (components(b)[k] - components(a)[k]) / delta });\n",
     [FOLD_TEST]),
    ("the witness refuses only far below zero", FOLD,
     "                if let Some(f) = witness.iter().find(|f| f.floor.is_nan() || f.floor <= 0.0) {\n",
     "                if let Some(f) = witness.iter().find(|f| f.floor.is_nan() || f.floor <= -1.0) {\n",
     [FOLD_TEST]),
    ("an undecided span after the last halving is accepted", FOLD,
     "            if depth >= MAX_FOLD_DEPTH || !halvable {\n",
     "            if depth >= MAX_FOLD_DEPTH || !halvable {\n"
     "                continue;\n            }\n            if false {\n", [FOLD_TEST]),
    ("a span without a bound at its width gives up", FOLD,
     "                if depth >= CONSERVATIVE_FOLD_DEPTH || !halvable {\n",
     "                if true {\n", [FOLD_TEST]),
    ("only a span with rates is decided exactly", FOLD,
     "        if self.fold_frames.has_bend(span.lo, span.hi) {\n",
     "        if self.fold_frames.rates(span.lo, span.hi, reading).is_some() {\n", [FOLD_TEST]),
    ("the exact compiler builds a turned fold", COMPILE_STATION,
     "    law.fold(margin)?.map_err(unsupported)\n",
     "    let _ = (&law, margin, unsupported);\n    Ok(())\n", [FOLD_TEST]),
    ("a span at a joint takes the earlier piece", FB,
     "        (holds && mid >= start && mid <= end).then(|| Located {\n",
     "        holds.then(|| Located {\n", [FOLD_TEST]),
    ("the plan reading's vertical part is all of sup |O'|", FOLD,
     "            let speed = bounds.point[0];\n            (speed * speed - bounds.forward * bounds.forward)",
     "            let speed = bounds.point[0];\n            (speed * speed - 0.0 * bounds.forward)", [TURNED]),
    ("the conservative bound refuses", FOLD,
     "            if depth >= CONSERVATIVE_FOLD_DEPTH {\n                return Ok(Err(reason));\n",
     "            if depth >= CONSERVATIVE_FOLD_DEPTH {\n"
     "                return Err(GeomError::InvalidInput(\"sections fold\".into()));\n",
     [BOUND_TEST]),
    # The sheets of open sections.
    ("a sectioned surface skips the fold check", *SHEET_CHECK,
     [SHEET + ["--", "a_sectioned_surface_reaching"]]),
    ("an open section run skips the fold check", *SHEET_CHECK,
     [SHEET + ["--", "an_open_section_run_reaching"]]),
    ("a sheet's undecided fold is reported certified", COMPILE_STATION,
     "    let bound = settle(&law, Ok(0.0), fold, &samples, &rows)?;\n",
     "    let bound = settle(&law, Ok(0.0), Ok(()), &samples, &rows)?;\n", [SHEET]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False)
