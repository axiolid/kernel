"""Mutation probe for certified placements at stations on curved bases
(#405, ADR 0082 amendment).

Each mutant breaks one step the bound rests on. The frame's error
(`FrameLaw::error_at`): the distance's resolution dropped from the point
and from the axes, the joint's distance dropped, the inexactly placed piece
accepted, the seam beside the station not refused, a frame read by
quadrature certified, the span not clipped to the side a station on a seam
names, and an intrinsic frame's comparison dropped. The placement
(`placement_accuracy`): the offsets' lever on the axes dropped, a turned
frame's `sqrt 3` taken as one, its turning not conditioned by the angle
between its axis and reference, the lever measured in the L2 norm, an exact
basis reported certified, and a basis placed inexactly certified. The
report: the frame's error not grown into every contribution, a certified
placement reported unbounded, and the source's reach dropped from the error. Each must turn
`evaluate/tests/station_frame_error.rs`, `accuracy::tests`,
`compile/tests/station_placement_bound.rs` or
`compile/tests/station_placement.rs` red.

Not probed, because no test can tell them apart:

- the closed-form read's `AXIS_ROUNDING` per axis and the rounding of
  placing the offsets, far below the resolution and any frame the tests
  build (the read point's rounding is probed: it scales with the curve's
  data, `1e-10` on a circle of radius 5000 about the origin of its frame);
- the room the retry leaves the source for the frame's error: a source
  meshed finer by that little reports the same bound within the budget;
- the plan reading's widening of the axis error (`2 e / rho`, as in #361).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/parametric/evaluate/src/station/frame_bound/error.rs"
A = "crates/execution/compile/src/station/accuracy.rs"
K = "crates/execution/compile/src/compiler.rs"
D = "crates/execution/compile/src/deviation.rs"
EVAL = ["-p", "axiolid-evaluate", "--test", "station_frame_error"]
UNIT = ["-p", "axiolid-mesh-compile", "--lib", "accuracy"]
DENSE = ["-p", "axiolid-mesh-compile", "--test", "station_placement_bound"]
PLACE = ["-p", "axiolid-mesh-compile", "--test", "station_placement"]

MUTANTS = [
    ('the resolution dropped from the point', E,
     '            point: (e_p + delta * s_point) * grow,',
     '            point: e_p * grow,', [EVAL]),
    ('the resolution dropped from the axes', E,
     '            axis: (e_a + delta * s_axis) * grow,',
     '            axis: e_a * grow,', [EVAL]),
    ('the joint distance dropped', E,
     '        let delta = here.measure + joint;',
     '        let delta = here.measure;', [EVAL]),
    ('an inexactly placed piece accepted', E,
     '                if !piece.placement_exact {\n                    return None;\n                }\n                let joint',
     '                if false {\n                    return None;\n                }\n                let joint', [EVAL]),
    ('an inexactly placed atomic curve accepted', E,
     '                if !piece.placement_exact {\n                    return None;\n                }\n                let at = admitted',
     '                if false {\n                    return None;\n                }\n                let at = admitted', [EVAL]),
    ('a seam beside the station not refused', E,
     '                } else if (c - u).abs() <= 2.0 * delta {',
     '                } else if false {', [EVAL]),
    ('a frame read by quadrature certified', E,
     '            Slot::Read(_) | Slot::Unbounded => return None,',
     '            Slot::Read(_) => (here.point_error, here.axis_error),\n'
     '            Slot::Unbounded => return None,', [EVAL]),
    ('the span not clipped to the side named', E,
     [('                        SeamSide::Incoming => hi = hi.min(c),',
       '                        SeamSide::Incoming => {}'),
      ('                        _ => lo = lo.max(c),',
       '                        _ => {}')], [EVAL]),
    ('the curve data dropped from the read rounding', E,
     '.max(data_scale(piece))',
     '', [EVAL]),
    ('a joint beside the station not refused', E,
     '            if ends[0] {',
     '            if false {', [EVAL]),
    ('an intrinsic frame not compared', E,
     '            Slot::Intrinsic(_) => (here.point_error, here.axis_error),',
     '            Slot::Intrinsic(_) => (0.0, 0.0),', [EVAL]),
    ('the offsets without a lever on the axes', A,
     '    let point_error = (read_point + read_axis * reach + POINT_ROUNDING * scale) * (1.0 + MARGIN);',
     '    let point_error = (read_point + POINT_ROUNDING * scale) * (1.0 + MARGIN);', [UNIT]),
    ('the turning not conditioned by the angle', A,
     '        Some(sin) => 3.0_f64.sqrt() * read_axis + TURN_ROUNDING / sin,',
     '        Some(_) => 3.0_f64.sqrt() * read_axis + TURN_ROUNDING,', [DENSE]),
    ('a turned frame bounded as the base frame', A,
     '        Some(sin) => 3.0_f64.sqrt() * read_axis + TURN_ROUNDING / sin,',
     '        Some(sin) => read_axis + TURN_ROUNDING / sin,', [UNIT]),
    ('the lever measured in the L2 norm', A,
     '    let reach = offsets.lateral.abs() + offsets.vertical.abs() + offsets.longitudinal.abs();',
     '    let reach = (offsets.lateral.powi(2) + offsets.vertical.powi(2)\n'
     '        + offsets.longitudinal.powi(2))\n'
     '    .sqrt();', [UNIT]),
    ('an exact basis reported certified', A,
     '    if parts.exact {',
     '    if false {', [DENSE]),
    ('a basis placed inexactly certified', A,
     '            placed_exactly: false,',
     '            placed_exactly: true,', [DENSE]),
    ('the frame not grown into every contribution', D,
     '                    bound: c.bound.placed(factor, rounding + frame),',
     '                    bound: c.bound.placed(factor, rounding),', [DENSE]),
    ('the source reach dropped from the frame', K,
     '        .fold(0.0, Scalar::max)\n        * (1.0 + 1e-12);',
     '        .fold(0.0, Scalar::min)\n        * (1.0 + 1e-12);', [DENSE]),
    ('a certified placement reported unbounded', K,
     '                        match frame_term(Some(accuracy), &source.mesh) {',
     '                        match frame_term(Some(accuracy), &source.mesh).and(None::<Scalar>) {',
     [DENSE, PLACE]),
]

probe_lib.main(MUTANTS, timeout=1800)
