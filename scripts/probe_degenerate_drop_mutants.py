"""Mutation probe for opt-in dropping of degenerate overlay rings (#256).

`overlay_with` and `union_soup_with` refuse a degenerate input ring by
default. With `OverlayOptions::dropping_degenerate_rings` they leave out a
ring that is too short, has a repeated vertex or has zero area under the
exact `ZeroArea` rule, report each one with operand, polygon, ring, reason
and area, and drop the holes of a dropped outer; every other fault stays a
refusal. Each mutant breaks one of those rules and must turn
`tests/degenerate_drop.rs` red.

Equivalent mutants, deliberately not listed:

- halving the area bound without the `next_up` for a subnormal result:
  only rings whose doubled area is subnormal reach it, and the bound then
  moves by one ulp of the smallest doubles.
- the doubled area summed with rounded `f64` arithmetic instead of
  `Dyadic`: on the dyadic test rings every product and sum is exact, so the
  values agree; the exact sum is kept because a general ring's rounded fan
  is not a bound.
"""
import probe_lib

D = "crates/algorithms/planar/overlay/src/"
G, L = D + "degenerate.rs", D + "lib.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "degenerate_drop"]

MUTANTS = [
    ('drop mode refuses like the default', G,
     '    if options.degenerate_rings == DegenerateRings::Refuse {\n        crate::validate(input, tolerance)?;',
     '    if true {\n        crate::validate(input, tolerance)?;', TESTS),
    ('new() drops by default', G,
     '        Self {\n            degenerate_rings: DegenerateRings::Refuse,\n        }',
     '        Self {\n            degenerate_rings: DegenerateRings::DropReported,\n        }', TESTS),
    ('Default drops', G,
     [('    #[default]\n    Refuse,', '    Refuse,'),
      ('    DropReported,\n}', '    #[default]\n    DropReported,\n}')], TESTS),
    ('frame not validated in drop mode', G,
     '    validate_frame(input, tolerance)?;\n', '', TESTS),
    ('self-intersection droppable', G,
     '        other => Err(other),', '        _ => Ok(Some(DropReason::ZeroArea)),', TESTS),
    ('repeated vertex dropped without the merged self-intersection test', G,
     '            if merged.points.len() >= 3 && self_intersects(&merged, tolerance) {',
     '            if false {', TESTS),
    ('merging keeps every vertex', G,
     '        if merged.last().is_none_or(|&q| (p - q).length() > linear) {',
     '        if merged.last().is_none_or(|_| true) {', TESTS),
    ('merging does not close the ring', G,
     '    while merged.len() > 1 && (merged[merged.len() - 1] - merged[0]).length() <= linear {',
     '    while false && (merged[merged.len() - 1] - merged[0]).length() <= linear {', TESTS),
    ('non-finite too-short ring droppable', G,
     '    if !ring.points.iter().all(|p| p.is_finite()) {\n        return Err(OverlayError::NonFinitePoint);\n    }',
     '', TESTS),
    ('ZeroArea not droppable', G,
     '        OverlayError::ZeroArea => Ok(Some(DropReason::ZeroArea)),\n', '', TESTS),
    ('too short not droppable', G,
     '        OverlayError::RingTooShort => Ok(Some(DropReason::TooShort)),\n', '', TESTS),
    ('holes of a dropped outer silently vanish', G,
     '                (None, Some(_)) => DropReason::OuterDropped,',
     '                (None, Some(_)) => continue,', TESTS),
    ('a hole of a dropped outer loses its own reason', G,
     '                (Some(reason), _) => reason,\n                (None, Some(_)) => DropReason::OuterDropped,',
     '                (_, Some(_)) => DropReason::OuterDropped,\n                (Some(reason), None) => reason,', TESTS),
    ('hole outside its kept outer not refused', G,
     '                    if hole_outside(&polygon.outer, hole, tolerance) {',
     '                    if false {', TESTS),
    ('a degenerate hole drops its whole polygon', G,
     '        if outer.is_none() {', '        if outer.is_none() && holes.len() == polygon.holes.len() {',
     TESTS),
    # A ring with a repeated vertex fed to the exact boolean does not
    # terminate (validation always precedes it), so the zero-area and
    # too-short tests, which see the dropped ring counted, run first.
    ('a dropped outer still used', G,
     '        if outer.is_none() {', '        if true {',
     [TESTS + ["--", "a_zero_area_ring_is_dropped", "a_ring_too_short_is_dropped"], TESTS]),
    ('hole index not reported', G,
     '            dropped.push(report(hole, place(RingRole::Hole(h)), reason));',
     '            dropped.push(report(hole, place(RingRole::Hole(0)), reason));', TESTS),
    ('polygon index not reported', G,
     '            polygon: index,', '            polygon: 0,', TESTS),
    ('area of a clockwise ring reported negative', G,
     '    if twice.sign() == Some(Sign::Negative) {', '    if false {', TESTS),
    ('area bound below the area', G,
     '    let bound = half_up(enclosure.hi());',
     '    let bound = half_up(enclosure.lo()) * (1.0 - 1e-15);', TESTS),
    ('clip drops reported as subject', L,
     '        clip,\n        OverlayOperand::Clip,',
     '        clip,\n        OverlayOperand::Subject,',
     TESTS),
    ('evidence counts the input rings', L,
     '        subject_rings: subject_polygons.iter().map(|p| 1 + p.holes.len()).sum(),',
     '        subject_rings: subject.polygons.iter().map(|p| 1 + p.holes.len()).sum(),', TESTS),
    ('union soup refuses in drop mode', L,
     '        if options.degenerate_rings == DegenerateRings::Refuse {\n            validate_ring(ring, tolerance)?;',
     '        if true {\n            validate_ring(ring, tolerance)?;', TESTS),
    ('union soup reports the wrong ring', L,
     '                    polygon: index,', '                    polygon: 0,', TESTS),
]

probe_lib.main(MUTANTS, reject_uncompiled=True, timeout=1800)
