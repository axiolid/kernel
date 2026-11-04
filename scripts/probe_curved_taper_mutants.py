"""Mutation probe for curved profiles that change between stations (#361,
ADR 0082 amendment).

Each mutant breaks one step the certified bound of a pairing rests on. The
chord deviation reported with the rings (`paired_profile_rings`): the
`h^2 m2 / 8` divisor made 32, so the deviation falls short of the polygon's
real distance from the interpolated section; a derived profile's stretch
dropped from the second derivative's bound and the joint snap; the ring
orientation skipped, so a bore keeps the wrong sense; a mirrored profile
accepted against an unmirrored one, whose vertices then pair the opposite
side (in the pairing, and for straight-edged sections in the compiler); the
orientation taken from eight samples per segment, which misjudges a thin
band; the fold check dropping the pairing bound from its margin;
the half step of circle-only runs taken for any run holding a circle. The
compiler's use of the pairing: `profile_bound` ignoring the paired bound, so
a changing curved profile is unbounded again. Each must turn
`construct/tests/paired_profile_rings.rs` or
`compile/tests/station_bound_curved_taper.rs` red.

Not probed, because no test can tell them apart:

- the `1 + 1e-12` inflation of the deviation (far below the slack the tests
  leave between the bound and the sampled distance);
- the vertex cap and the near-duplicate refusal (their effect is a `None`
  that the tests build only through another cause: a segment count).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
PAIRED = "crates/algorithms/construction/construct/src/profile/paired.rs"
STATION = "crates/execution/compile/src/station.rs"

PAIRED_TEST = ["-p", "axiolid-construct", "--test", "paired_profile_rings"]
TAPER_TEST = [
    "-p", "axiolid-mesh-compile", "--all-features", "--test", "station_bound_curved_taper",
]

FOLD_TEST = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "station_fold"]

MUTANTS = [
    ("the chord deviation divides by 32", PAIRED,
     "chord = chord.max(step * step * m2 / 8.0);",
     "chord = chord.max(step * step * m2 / 32.0);", [PAIRED_TEST, TAPER_TEST]),
    ("a derived profile's stretch is dropped", PAIRED,
     "inner.stretch *= stretch2(transform);",
     "inner.stretch *= 1.0;", [PAIRED_TEST]),
    ("the ring orientation is skipped", PAIRED,
     "let reverse = (r == 0 && area < 0.0) || (r > 0 && area > 0.0);",
     "let reverse = false;", [PAIRED_TEST]),
    ("a mirror mismatch is accepted", PAIRED,
     "if lowered.iter().any(|l| l.reversed != lowered[0].reversed) {",
     "if false && lowered.iter().any(|l| l.reversed != lowered[0].reversed) {", [PAIRED_TEST]),
    ("the orientation comes from samples", PAIRED,
     "area += segment_area(segment)?;",
     "{\n            let pts = sample(segment, 8, 0.0)?;\n"
     "            area += 0.5 * pts.windows(2).map(|w| w[0].x * w[1].y - w[1].x * w[0].y)"
     ".sum::<Scalar>();\n        }", [PAIRED_TEST]),
    ("the straight-edged mirror mismatch is accepted", STATION,
     "if !one && !curved && profiles.iter()",
     "if false && !one && !curved && profiles.iter()", [TAPER_TEST]),
    ("the fold check drops the paired margin", STATION,
     "law.fold(margin)?.map_err(unsupported)",
     "law.fold(0.0)?.map_err(unsupported)", [FOLD_TEST]),
    ("profile_bound ignores the paired bound", STATION,
     "        return Ok(Ok(bound));\n",
     "        let _ = bound;\n", [TAPER_TEST]),
    ("the half step applies to mixed runs", PAIRED,
     "let circles = profiles.iter().all(|p| matches!(p, Profile::Circle(_)));",
     "let circles = profiles.iter().any(|p| matches!(p, Profile::Circle(_)));", [PAIRED_TEST]),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False)
