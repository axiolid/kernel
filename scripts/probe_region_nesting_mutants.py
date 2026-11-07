"""Mutation probe for the nesting of a region's rings (#448,
`crates/algorithms/planar/overlay/src/nesting.rs` and `src/region.rs`).

Each mutant weakens one step -- the check `Region::new` and `Region::translate` run, the regions
it skips, the rings it turns or leaves out, the crossings it refuses
before the arrangement and the names it gives them, how each side of a
piece is read, and the sweep of an edge's shared stretches -- and must turn a test red.

Equivalent mutants, deliberately not listed:
- leaving out the early return for fewer than two arranged rings: the
  arrangement of one ring has every piece held by that ring alone, which
  no rule refuses;
- dropping the dedup of `held`: carriers and the rings
  holding a piece are disjoint, so no ring is listed twice.
"""
import probe_lib

N = "crates/algorithms/planar/overlay/src/nesting.rs"
G = "crates/algorithms/planar/overlay/src/region.rs"
ARR = "crates/algorithms/planar/overlay/src/exact_arc/arrangement.rs"
SWEEP = [
    ["-p", "axiolid-overlay", "--lib", "exact_arc::arrangement::tests"],
    ["-p", "axiolid-overlay", "--test", "region_nesting"],
]
TESTS = [
    ["-p", "axiolid-overlay", "--test", "region_nesting"],
    ["-p", "axiolid-overlay", "--lib", "nesting::"],
]

MUTANTS = [
    ("Region::new runs no check", G,
     "tolerance)?;\n        crate::nesting::check(&polygons)?;\n        let mut region",
     "tolerance)?;\n        let _ = crate::nesting::check(&polygons);\n        let mut region", TESTS),
    ("translate skips the check", G,
     "        crate::nesting::check(&polygons)?;\n        Ok(Self::from_valid(polygons, !self.is_empty()))",
     "        Ok(Self::from_valid(polygons, !self.is_empty()))", TESTS),
    ("a region of two rings not checked", N,
     ".sum::<usize>() < 2 {",
     ".sum::<usize>() < 3 {", TESTS),
    ("rings not turned counter-clockwise", N,
     "            if orientation(&points) == Sign::Negative {",
     "            if false && orientation(&points) == Sign::Negative {", TESTS),
    ("rings on one line arranged", N,
     "            if on_one_line(&ring.points) {",
     "            if false && on_one_line(&ring.points) {", TESTS),
    ("a hole of an outer ring on one line not refused", N,
     "        .any(|r| r.hole && outer_of[r.polygon].is_none())",
     "        .any(|r| false && r.hole && outer_of[r.polygon].is_none())", TESTS),
    ("no crossing refused before the arrangement", N,
     "            if proper_crossing(a, b, c, d) {",
     "            if false && proper_crossing(a, b, c, d) {", TESTS),
    ("every crossing named by the rings' roles", N,
     "!a.hole && !b.hole && !holed[a.polygon] && !holed[b.polygon]",
     "true", TESTS),
    ("a crossing of two holes named a hole leaving", N,
     "        OverlayError::HolesOverlap\n    } else {",
     "        OverlayError::HoleOutsideOuter\n    } else {", TESTS),
    ("a carrier's other side taken for its inside", N,
     "held.extend(edge.sources.iter().filter(|s| s.2 == left).map(|s| s.0));",
     "held.extend(edge.sources.iter().filter(|s| s.2 != left).map(|s| s.0));", TESTS),
    ("a hole's outside read for its inside", N,
     "for source in edge.sources.iter().filter(|s| s.2 == left) {",
     "for source in edge.sources.iter().filter(|s| s.2 != left) {", TESTS),
    ("a hole outside its outer ring accepted", N,
     "outside |= held.binary_search(&outer).is_err();",
     "outside |= false && held.binary_search(&outer).is_err();", TESTS),
    ("a hole inside another of its polygon accepted", N,
     "crowded |= held",
     "crowded |= false && held", TESTS),
    ("a polygon holds a side one of its holes holds", N,
     "!role.hole\n                && !held",
     "!role.hole\n                && true || false && !held", TESTS),
    ("two polygons holding one side accepted", N,
     "overlap |= polygons_holding(&held, roles) > 1;",
     "overlap |= polygons_holding(&held, roles) > 2;", TESTS),
]

MUTANTS += [
    ("sweep uses the checked translate", G,
     "let moved = Self::from_valid(self.shifted(direction)?, true);",
     "let moved = self.translate(direction)?;", TESTS),
    ("a shared stretch never leaves the sweep", ARR,
     "while next_end < by_end.len() && spans[by_end[next_end]].2 <= at {",
     "while false && next_end < by_end.len() && spans[by_end[next_end]].2 <= at {", SWEEP),
    ("a shared stretch enters the sweep a piece late", ARR,
     "while next_start < by_start.len() && spans[by_start[next_start]].1 <= at {",
     "while next_start < by_start.len() && spans[by_start[next_start]].1 + 1 <= at {", SWEEP),
]

probe_lib.main(MUTANTS, timeout=900)
