"""Mutation probe for certified lengths beside a planar 3D intrinsic curve
(#302, ADR 0082 amendment).

A `Curve3::Intrinsic` with a torsion law identically zero (read
structurally from the stored law) and a start tangent and normal with no
`z` component is a 2D intrinsic curve lifted to a level plane, so its
offsets' lengths are certified by the midpoint rule of the 2D case. Each
mutant breaks one step of that: the torsion test dropped (a twisted curve,
or one whose torsion law is zero only once evaluated, read as planar),
either half of the level-plane test dropped (a plane rolled about the
start tangent, or pitched about the start normal, read as level), the
plan's handedness dropped (a curve turning right read as turning left),
the planar curve refused again, and its offset measured in plan (a rise
not lengthening it). Each must turn
`evaluate/tests/station_certify_planar_intrinsic.rs` red.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
CERTIFY = "crates/algorithms/parametric/evaluate/src/station/certify.rs"

PLANAR = ["-p", "axiolid-evaluate", "--all-features", "--test",
          "station_certify_planar_intrinsic"]

TEST = "    if !intrinsic.is_planar() || x.z != 0.0 || y.z != 0.0 {"

MUTANTS = [
    ("a twisted intrinsic is read as planar", CERTIFY, TEST,
     "    if x.z != 0.0 || y.z != 0.0 {", PLANAR),
    ("a plane pitched about the start normal is read as level", CERTIFY, TEST,
     "    if !intrinsic.is_planar() || y.z != 0.0 {", PLANAR),
    ("a plane rolled about the start tangent is read as level", CERTIFY, TEST,
     "    if !intrinsic.is_planar() || x.z != 0.0 {", PLANAR),
    ("the plan's handedness is dropped", CERTIFY,
     "        sign: handedness(&plan)?,",
     "        sign: handedness(&plan)?.abs(),", PLANAR),
    ("a planar 3D intrinsic is refused again", CERTIFY,
     "                (level_curvature(intrinsic)?, None, true)",
     "                return Err(level_curvature(intrinsic).err().unwrap_or_else(uncertified))",
     PLANAR),
    ("its offset is measured in plan", CERTIFY,
     "                (level_curvature(intrinsic)?, None, true)",
     "                (level_curvature(intrinsic)?, None, false)", PLANAR),
]

probe_lib.main(MUTANTS, root=ROOT, timeout=TIMEOUT, baseline=False)
