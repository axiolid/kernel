"""Mutation probe for sampling edge curves at their breaks (#342).

Each mutant weakens one decision of `crates/execution/compile/src/brep/
trim.rs`, `inversion.rs` or `edges.rs`: the fractions at which an edge
curve has a derivative break, their exactness, and the per-piece requirement
the edge pass keeps. `tests/edge_inversion.rs` or the modules' unit tests
must fail on it. A build error also counts as killed.

Not listed, as equivalent at test resolution:
- only the unfinished pieces of an edge being halved, not all of them:
  every fixture's pieces fit or fail together (a polyline's all fit at one
  span; a circle spline's four quarters are alike), so no fixture has a
  piece that finishes before its neighbour.
- the segment cap counted over all pieces of an edge, not per piece: no
  fixture reaches the cap with more than one piece.
"""
import probe_lib

I = "crates/execution/compile/src/brep/inversion.rs"
E = "crates/execution/compile/src/brep/edges.rs"
T = "crates/execution/compile/src/brep/trim.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "edge_inversion"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::trim"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::inversion"],
]

MUTANTS = [
    ("an edge curve has no breaks", T,
     "            .cuts(curve, lo, hi)\n            .into_iter()\n",
     "            .cuts(curve, lo, hi)\n            .into_iter()\n            .take(0)\n"),
    ("a break fraction is not nudged to its parameter", T,
     "            .map(|t| self.exact_fraction(t, (t - self.span.start) / length))\n",
     "            .map(|t| (t - self.span.start) / length)\n"),
    ("break fractions rounded through single precision", T,
     "            .map(|t| self.exact_fraction(t, (t - self.span.start) / length))\n",
     "            .map(|t| ((t - self.span.start) / length) as f32 as Scalar)\n"),
    ("a piece's requirement on an edge is the first face's, not the strictest", E,
     "                    *kept = (*kept).max(n);\n",
     "                    let _ = n;\n"),
    ("an edge's planned samples ignore its breaks", E,
     "        let fractions = inversion::piece_fractions(&trace.breaks(), counts);\n",
     "        let fractions = inversion::piece_fractions(&[], &counts[..1]);\n"),
    ("an edge's requirement is taken without its breaks", E,
     "                    let breaks = if uniform { Vec::new() } else { trace.breaks() };\n",
     "                    let breaks = Vec::new();\n"),
    ("the certified terms of a piece are ignored", I,
     "                fits = certified.is_none_or(|worst| worst <= budget);\n",
     "                fits = certified.is_none_or(|worst| worst <= budget) || true;\n"),
    ("pieces are joined with their shared end twice", I,
     "        let skip = usize::from(k > 0);\n",
     "        let skip = usize::from(k > usize::MAX);\n"),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
