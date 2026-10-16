"""Mutation probe for the merged and trimmed edge hulls of a Minkowski sum
or erosion (#323, `crates/algorithms/planar/overlay/src/cover.rs`).

Each mutant weakens one step -- merging only vertices exactly between their
neighbours, each condition proving that a chord lies in the summed set
(the spanned vertices' side, a vertex on the chord turning away from the
set, the chord's ends, the sectors at both ends, other boundaries crossing
it or touching it, its own ring past the run), trimming a hull to what the
chord's hull does not hold, the tilt of cuts parallel to the chord, the
side of the set read for each ring, and the chords being used at all --
and must turn a test red.

Equivalent mutants, deliberately not listed:
- trimming against the cover's interior instead of the closed cover:
  fewer vertices count as covered, so more of each hull is kept and the
  union is the same (only the cost grows);
- giving a chord's cover to an edge outside its run: trimming any term by
  any other term of the union leaves the union as it is, the cover being
  a subset of the sum;
- dropping the depth bound, or the deduplication of equal slivers, or the
  anchor at the origin: each changes only how many rings the subdivision
  gets, never their union;
- the guard against a chord whose ends coincide: a valid ring has no
  repeated vertex, and merging keeps the vertices distinct.
"""
import probe_lib

C = "crates/algorithms/planar/overlay/src/cover.rs"
M = "crates/algorithms/planar/overlay/src/minkowski.rs"
# The proofs and trimming are checked in the unit tests; the whole route,
# against the untrimmed hulls, there and in the Minkowski tests.
UNIT = ["-p", "axiolid-overlay", "--lib", "--test", "minkowski"]

MUTANTS = [
    ('merging drops a vertex where the boundary turns back', C,
     '        !(orient(a, b, c) == Sign::Zero && between(a, b, c))',
     '        !(orient(a, b, c) == Sign::Zero)', UNIT),
    ('merging within a tolerance', C,
     '        !(orient(a, b, c) == Sign::Zero && between(a, b, c))',
     '        !(((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)).abs() < 1e-12 && between(a, b, c))', UNIT),
    ('chords proposed over the wrong vertices', C,
     '            if orient(at(u), at(w), v) == away {',
     '            if orient(at(u), at(w), v) != away {', UNIT),
    ('a spanned vertex on the set side accepted', C,
     '            s if s == side => return false,',
     '            s if s == side => {}', UNIT),
    ('a vertex on the chord may turn towards the set', C,
     '                    || orient(p[(i + n - 1) % n], p[i], p[(i + 1) % n]) != side.flip()',
     '                    || false', UNIT),
    ("a vertex on the chord's line past its end accepted", C,
     '                if !on_segment(a, b, p[i])',
     '                if false', UNIT),
    ("the start's sector not checked", C,
     '    if !leaves_inside(a, next(chord.from), prev(chord.from), b, side)',
     '    if false', UNIT),
    ("the end's sector not checked", C,
     '        || !leaves_inside(b, prev(to), next(to), a, side.flip())',
     '        || false', UNIT),
    ('a sector read past its far edge', C,
     '        s if s == side => orient(at, prev, other) != side,',
     '        s if s == side => true,', UNIT),
    ('a spike taken as a sector', C,
     '        Sign::Zero => false,',
     '        Sign::Zero => true,', UNIT),
    ('other boundaries may cross the chord', C,
     '        !(opposite(orient(a, b, c), orient(a, b, d)) && opposite(orient(c, d, a), orient(c, d, b)))',
     '        true', UNIT),
    ('other vertices may lie on the chord', C,
     '        if (!end(e) && on_segment(a, b, c)) || (!end(f) && on_segment(a, b, d)) {',
     '        if false {', UNIT),
    ("the chord's own ring skipped past its run", C,
     '        if r == ring && (e + length - first) % length < run {',
     '        if r == ring {', UNIT),
    ('a hull the cover holds kept', C,
     '    if covered.iter().all(|&c| c) {',
     '    if false {', UNIT),
    ('trimming drops the ends of a covered run', C,
     '        .map(|i| !(covered[i] && covered[(i + n - 1) % n] && covered[(i + 1) % n]))',
     '        .map(|i| !covered[i])', UNIT),
    ('cuts along the chord not tilted', C,
     '        if (j + n - i) % n > 2 && parallel(hull[i], hull[j], a, b) {',
     '        if false {', UNIT),
    ('a trimmed ring that tidying changes kept', C,
     '        Some(same) if same == kept => Some(kept),',
     '        Some(_) => Some(kept),', UNIT),
    ("the set's side read the wrong way round", M,
     '                    || (shape.starts.binary_search(&(i - offset)).is_ok() == outer_inside);',
     '                    || (shape.starts.binary_search(&(i - offset)).is_ok() != outer_inside);', UNIT),
    ('no chord used', M,
     '        if covers && !convex.is_empty() {',
     '        if covers && convex.is_empty() {', UNIT),
]

# Removing the whole-hull exit leaves the walk over kept vertices with
# none to find: it never ends, and counts as killed by timeout.
probe_lib.main(MUTANTS, timeout=300)
