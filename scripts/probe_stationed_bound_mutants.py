"""Mutation probe for the certified bound and the exact bodies of runs of
sections at stations (#247).

Each mutant breaks one step the bound or the exact path rests on: the
frame law's derivative bounds (a circle's curvature, the reference-up
frame's normalisation, the plan frame, the rolled frame, the roll's cant
term, an elevated curve's tangent and point laws), the per-cell bound
(the chord term, the wall quad, the breakpoints at seams, splitting at
seams, the mitre rules, the turning orientation, the offsets' rate), the
refinement against the bound, the profile's share, and the exact
compiler's extrusion length, revolution angle and axis, straight
polyline test, curved-taper refusal, and the loft's planarity test and
topology. Each must turn a test red: `evaluate`'s `frame_bound` unit tests
compare the bounds with second differences of the frames, and
`compile/tests/station_bound.rs` checks that every run is certified
within its budget, that a dense sampling of the exact surface never
exceeds the bound, the unbounded reasons, and the exact bodies' volumes,
areas and refusals.

Not probed, because no test can tell them apart: the arc-length inverse's
measure slack (it is `1e-12` relative, below the rounding every sampled
distance carries), the placement's stretch (`1 + 1e-14`), a turning
orientation's floor `rho` on the interpolated unit vectors (they stay far
from zero in any orientation a graph accepts), the refusal of a
plan-measured piece under a tilting placement (the evaluator refuses it
first), and the tangent's share of the reference-up frame's up axis
(`|U'| <= |R'| + |T'|`): with `T' = a R + b U` and `R' = -a T + c U` the
up turns by `sqrt(b^2 + c^2)`, which `|R'| <= k1 / rho` alone already
covers wherever a test can read the frame, so dropping `k1` does not
undercut the true derivative there.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
FB = "crates/algorithms/parametric/evaluate/src/station/frame_bound.rs"
BK = "crates/algorithms/parametric/evaluate/src/banked.rs"
SB = "crates/execution/compile/src/station/bound.rs"
SC = "crates/execution/compile/src/station.rs"
SX = "crates/execution/compile/src/station/exact.rs"
LX = "crates/algorithms/construction/construct/src/loft_exact.rs"
FRAMES = ["-p", "axiolid-evaluate", "--lib", "frame_bound"]
BOUND = ["-p", "axiolid-mesh-compile", "--test", "station_bound"]

MUTANTS = [
    # --- the frame law's derivative bounds ---------------------------------
    ('circle point without curvature', FB,
     '        point: [1.0, k],',
     '        point: [1.0, 0.0],', [FRAMES, BOUND]),
    ('circle tangent does not turn', FB,
     '        tangent: [k, k * k],',
     '        tangent: [0.0, k * k],', [FRAMES, BOUND]),
    ('reference-up lateral not normalised', FB,
     '            let [r1, r2] = normalised(k1, k2, horizontal()?);',
     '            let [r1, r2] = [k1, k2];', [FRAMES, BOUND]),
    ('rolled frame ignores the roll', FB,
     '            let axis = [p1 + a1, p2 + p1 * p1 + 2.0 * p1 * a1 + a2];',
     '            let axis = [a1, a2];', [FRAMES, BOUND]),
    ('roll without the cant rate', BK,
     '        axiolid_curve::BankConvention::TangentRotation => (d0 / b, d1 / b, d2 / b),',
     '        axiolid_curve::BankConvention::TangentRotation => (d0 / b, 0.0, 0.0),',
     [FRAMES, BOUND]),
    ('elevated tangent without its second derivative', FB,
     '            2.0 * d3 / sigma + 6.0 * d2 * d2 / (sigma * sigma),',
     '            0.0,', [FRAMES, BOUND]),
    ('elevated point without curvature', FB,
     '        point: [bounds.first, d2],',
     '        point: [bounds.first, 0.0],', [FRAMES, BOUND]),
    ('plan frame not normalised', FB,
     '            let flat = normalised(k1, k2, horizontal()?);',
     '            let flat = [k1, k2];', [FRAMES, BOUND]),
    # --- the bound per cell ------------------------------------------------
    ('chord term dropped', SB,
     '                chord[k] = chord[k].max(0.125 * span * span * second);',
     '                chord[k] = chord[k].max(0.0 * span * second);', [BOUND]),
    ('wall quad dropped', SB,
     '            worst = worst.max(along + wall + slide[k].max(slide[l]));',
     '            worst = worst.max(along + 0.0 * wall + slide[k].max(slide[l]));', [BOUND]),
    ('breakpoints dropped', SB,
     '                corner[k] = corner[k].max((exact[k] - on_chord).length());',
     '                corner[k] = corner[k].max(0.0 * (exact[k] - on_chord).length());',
     [BOUND]),
    ('cells not split at seams', SB,
     '            .filter(|&c| c > lo + seam_slack(lo) && c < hi - seam_slack(hi))',
     '            .filter(|_| false)', [BOUND]),
    ('mitre beside a curved piece bounded', SB,
     '            if mitred && !bounds.is_straight() {',
     '            if false {', [BOUND]),
    ('mitre beside a changing section bounded', SB,
     '        if mitred && changes {',
     '        if false {', [BOUND]),
    ('turning orientation adds nothing to the lateral', SB,
     '        lateral: [n1 + x1, n2 + 2.0 * n1 * x1 + x2],',
     '        lateral: [0.0, 0.0],', [BOUND]),
    ('turning orientation adds nothing to the up', SB,
     '        up: [n1, n2],',
     '        up: [0.0, 0.0],', [BOUND]),
    ('offsets moving between stations add nothing', SB,
     '                    + 2.0 * swing',
     '                    + 0.0 * swing', [BOUND]),
    ('offsets add no frame curvature', SB,
     '                    + frame(1)',
     '                    + 0.0 * frame(1)', [BOUND]),
    ('profile points moving add nothing', SB,
     '                    + 2.0 * (xr * lt[0] + yr * ut[0]);',
     '                    + 0.0 * (xr * lt[0] + yr * ut[0]);', [BOUND]),
    # --- refinement and the profile's share ----------------------------------
    ('refinement ignores the bound', SC,
     '                passes(&at_lo, &at_hi, at_mid, edges, chord) && accept(lo, hi, &at_lo, &at_hi)?',
     '                passes(&at_lo, &at_hi, at_mid, edges, chord)', [BOUND]),
    ('different curved profiles bounded', SC,
     '    Ok(if worst == 0.0 || one {',
     '    Ok(if true {', [BOUND]),
    # --- exact bodies -----------------------------------------------------
    ('extruded by the end distance', SX,
     '            let prism = extrude_profile_exact(profile(first)?, direction, s1 - s0, tolerance)',
     '            let prism = extrude_profile_exact(profile(first)?, direction, s1, tolerance)',
     [BOUND]),
    ('revolved by the end distance', SX,
     '                (s1 - s0) / circle.radius,',
     '                s1 / circle.radius,', [BOUND]),
    ('revolved about the mirrored axis', SX,
     '                Point3::new(centre.dot(base.lateral), 0.0, 0.0),',
     '                Point3::new(-centre.dot(base.lateral), 0.0, 0.0),', [BOUND]),
    ('bent polyline read as straight', SX,
     '        if (offset - direction * t).length() > linear || t < along - linear {',
     '        if t < along - linear {', [BOUND]),
    ('curved taper built exactly', SX,
     '                    ProfileDeviation::Bounded(d) if d == 0.0',
     '                    ProfileDeviation::Bounded(_)', [BOUND]),
    ('warped wall built as a plane', LX,
     '                if corners.iter().any(|p| (*p - centre).dot(z).abs() > linear) {',
     '                if false {', [BOUND]),
    ('wall joined to the wrong rail', LX,
     '                    (rails[j][r][m], Orientation::Forward, line(b, d)),',
     '                    (rails[j][r][k], Orientation::Forward, line(b, d)),', [BOUND]),
]

probe_lib.main(MUTANTS, timeout=1800)
