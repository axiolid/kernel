"""Mutation probe for chords whose trimmed hulls degenerate in a Minkowski
sum or erosion (#331, `crates/algorithms/planar/overlay/src/minkowski.rs`).

Along a wall that runs parallel to a side of the convex piece, with a
corner an ulp off the chord across it, trimming an edge hull leaves a
sliver of three nearly collinear vertices that the arrangement refuses
(`ZeroArea`). Such a chord is not used: its hull and all of its trimmed
hulls are dropped and the edges it spans keep their whole hulls. Each
mutant weakens that rule -- taking the sliver, measuring it against
another bound than the arrangement's, keeping the chord while one edge
falls back to its whole hull or loses the sliver, or dropping every chord
-- and must turn a test red.

Equivalent mutant, deliberately not listed:
- using a chord whose own hull has degenerated: that hull holds a
  translate of the piece, as every edge hull does, and is longer than any
  edge it spans, so it is never smaller than the edge hulls, which were
  terms before #323 and still are; the check only keeps every step of
  the route from emitting a ring of zero area.
"""
import probe_lib

M = "crates/algorithms/planar/overlay/src/minkowski.rs"
# The 0.3.15 digest of the synthetic walls is pinned in the Minkowski
# tests; that chords still trim is checked in the unit tests.
UNIT = ["-p", "axiolid-overlay", "--lib", "--test", "minkowski"]

DROP_TEST = (
    '.any(|term| matches!(term, EdgeTerm::Trimmed(kept) if degenerate(kept)))',
    '.any(|term| matches!(term, EdgeTerm::Trimmed(kept) if false && degenerate(kept)))',
)

MUTANTS = [
    ('a trimmed hull of zero area taken as a term', M, *DROP_TEST, UNIT),
    ("zero area measured against another bound than the arrangement's", M,
     '    let degenerate = |points: &[Point2]| crate::zero_area(points, tolerance.linear());',
     '    let degenerate = |points: &[Point2]| crate::zero_area(points, 0.0);', UNIT),
    ('a degenerate trim keeps its whole hull, the chord kept', M,
     '                                    tidy_points(points, tolerance)\n',
     '                                    tidy_points(points, tolerance).filter(|p| !degenerate(p))\n',
     UNIT),
    ('a degenerate trim dropped, the chord kept', M,
     [DROP_TEST,
      ('                    EdgeTerm::Trimmed(kept) => {',
       '                    EdgeTerm::Trimmed(kept) if degenerate(&kept) => {}\n'
       '                    EdgeTerm::Trimmed(kept) => {')],
     UNIT),
    ('every chord dropped', M,
     '                if hulls.iter().flatten().any(|g| degenerate(g))',
     '                if true || hulls.iter().flatten().any(|g| degenerate(g))', UNIT),
]

probe_lib.main(MUTANTS, timeout=600)
