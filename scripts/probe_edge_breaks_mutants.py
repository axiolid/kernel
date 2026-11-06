"""Mutation probe for sampling edge curves at their breaks (#342).

Each mutant weakens one decision of `crates/execution/compile/src/brep/
trim.rs`, `inversion.rs` or `edges.rs`: the fractions at which an edge
curve has a derivative break, their exactness, and the per-piece requirement
the edge pass keeps. `tests/edge_inversion.rs` or the modules' unit tests
must fail on it. A build error also counts as killed.

First run: 12 / 12 killed.

Not listed, as untested at test resolution (no fixture, not equivalent):
- the clamp of a sample's shifted parameter into its piece's native
  ends (#426): rounding is monotone, so the plain shift never leaves them
  (a sample past a knot copy maps to at least the knot, one before it to at
  most the knot; none was found across the seam), and the clamp is a
  guarantee by construction, not an observable.
- the seam mapped to the domain start by subtraction, not by identity: a
  piece starting a hair below the domain start is harmless to the bounds.
- the segment cap counted over all pieces of an edge, not per piece: no
  fixture reaches the cap with more than one piece.
- the merged per-piece requirement of an edge checked against the cap: no
  fixture has faces whose per-piece maxima sum past it.
- the Uniform arm of the requirement merge keeping the first face's count:
  no fixture shares a pcurve edge between faces asking for different counts.
- the midpoints of a halving taken at the plain mean of the fractions, not
  the doubled sampling's odd fractions: they differ by an ulp.
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
     "            .cuts(lo, hi)\n            .into_iter()\n",
     "            .cuts(lo, hi)\n            .into_iter()\n            .take(0)\n"),
    ("a break fraction is not nudged to its parameter", T,
     "            .map(|t| self.exact_fraction(t, (t - self.span.start) / length))\n",
     "            .map(|t| (t - self.span.start) / length)\n"),
    ("break fractions rounded through single precision", T,
     "            .map(|t| self.exact_fraction(t, (t - self.span.start) / length))\n",
     "            .map(|t| ((t - self.span.start) / length) as f32 as Scalar)\n"),
    ("a break's copy past the seam is shifted back inexactly (#426)", T,
     "        self.breaks.get(i).copied().filter(|b| b + period == t)\n",
     "        self.breaks.get(i).copied().filter(|b| b + period == t && false)\n"),
    ("a sample at a break's copy is shifted back inexactly (#426)", T,
     "        if let Some(exact) = self.exact_native(t) {\n            return exact;\n        }\n",
     ""),
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
    ("an edge with more pieces than the cap allows is still cut", I,
     "                if breaks.len() + 1 > MAX_CURVED_EDGE_SEGMENTS {\n",
     "                if breaks.len() + 1 > usize::MAX {\n"),
    ("every piece is halved, not only the unfinished ones", I,
     "        for &k in &unfinished {\n            pieces[k].halve(",
     "        for k in 0..counts.len() {\n            pieces[k].halve("),
    ("pieces are joined with their shared end twice", I,
     "        let skip = usize::from(k > 0);\n",
     "        let skip = usize::from(k > usize::MAX);\n"),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
