"""Mutation probe for the fast disc morphology (#330,
`crates/algorithms/planar/overlay/src/fill.rs` and the mode switch in
`minkowski.rs`).

Each mutant weakens one step -- the mode switch either way, the side of
the set whose notches are filled, each condition on a chord (its exact
proof, the depth within the disc or the slab, the budget, no ring inside its notch),
the excess (the gap term, the supports, the depth bound, the slab bound), splitting a
refused chord (at the vertex nearest the middle for a deep one, #428), the
proposals' reach over needles and its extent box (#428), the pieces of a deep
notch's bound, the trim of a chord not proved (#345),
the ring a fill would close, the inscribed polygon drawn in
by the excess and the deviation carrying it, and the budget's refusal --
and must turn a test red.

Equivalent mutants, deliberately not listed:
- the `f64` slack on the excess (64 ulps of the largest coordinate) and its
  final widening: they cover rounding of the excess's own arithmetic, far
  below what any sampled test point can resolve;
- splitting a chord too deep or over the budget, not deeper than the disc,
  at another vertex than the farthest: every split ends in single edges, so only which shorter
  chords are found changes, and with it the cost and the excess, both
  within the stated bounds (the split of a chord not proved, which
  decides whether a wall of teeth is filled at all, is listed);
- a notch deeper than the disc refused without the slab: `excess` then
  returns the slab bound or nothing either way, the guard after it being
  the same test;
- the order of the halving limits over the supports: the bound is a
  minimum over all of them.
- the box of a part in `Along::merged` with its corners on one side of the
  part only: proposals are a choice only; `excess` and the proof check
  every chord.
"""
import probe_lib

F = "crates/algorithms/planar/overlay/src/fill.rs"
C = "crates/algorithms/planar/overlay/src/cover.rs"
M = "crates/algorithms/planar/overlay/src/minkowski.rs"
UNIT = ["-p", "axiolid-overlay", "--lib", "fill::"]
ROUTE = ["-p", "axiolid-overlay", "--test", "minkowski"]
BOTH = [UNIT, ROUTE]
PROPOSALS = ["-p", "axiolid-overlay", "--lib", "cover::"]
ALL = [PROPOSALS, UNIT, ROUTE]

MUTANTS = [
    ('the exact mode fills notches', M,
     '            MorphologyMode::Exact => None,',
     '            MorphologyMode::Exact => return self.disc(radius, side, erode, MorphologyMode::fast(1e-3), tolerance),',
     ROUTE),
    ('the fast mode fills nothing', F,
     '    if !(t > 0.0 && t.is_finite() && budget >= 0.0 && budget.is_finite()) {',
     '    if true {', BOTH),
    ("an erosion fills the region's notches, not its complement's", M,
     '                        *s = s.flip();',
     '                        *s = *s;', ROUTE),
    ('an inscribed polygon not drawn in by the excess', M,
     '                    reach -= filled.excess;',
     '                    reach -= 0.0;', ROUTE),
    ("an inscribed polygon's deviation without the excess", M,
     '                    deviation = radius - reach * half.cos();',
     '                    let _ = half;', ROUTE),
    ("a circumscribed polygon's deviation without the excess", M,
     '                    deviation += filled.excess;',
     '                    deviation += 0.0;', ROUTE),
    ('a negative or non-finite budget accepted', M,
     '            if !budget.is_finite() || budget < 0.0 {',
     '            if false {', ROUTE),
    ('an unproved chord used', F,
     '            if !cover::proved(r, &boundary, chord, outside, &index) {',
     '            if false {', BOTH),
    ('a chord over the budget used', F,
     '                Some(e) if e <= budget => {',
     '                Some(e) => {', BOTH),
    ('a notch holding another ring filled', F,
     '                    if encloses_none(r, &polygon, &index) {',
     '                    if true {', BOTH),
    ('a ring inside the notch not found', F,
     '        !inside(polygon, index.rings[r][e])',
     '        true', BOTH),
    ('a notch deeper than the disc bounded by its depth and supports', F,
     '    if deepest > t {',
     '    if false {', BOTH),
    ('the slab test dropped: every notch in the slab', F,
     '    let slab = at.iter().all(|&(s, _)| s >= slack && s <= len - slack);',
     '    let slab = true;', BOTH),
    ('the slab never holds', F,
     '    let slab = at.iter().all(|&(s, _)| s >= slack && s <= len - slack);',
     '    let slab = false;', BOTH),
    ('the slab bound without the disc radius in its denominator', F,
     '        h / ((t * t + h).sqrt() + t) * (1.0 + 1e-12)',
     '        h / (t * t + h).sqrt() * (1.0 + 1e-12)', BOTH),
    ('the slab bound with the half chord as its gap', F,
     '    let g = len + slack;',
     '    let g = 0.5 * len + slack;', BOTH),
    ('the depth ceiling dropped on a notch deeper than the disc', F,
     '        return (wide.map(|w| w.min(deepest)).filter(|e| e.is_finite()), true);',
     '        return (wide.filter(|e| e.is_finite()), true);', BOTH),
    ('proposals reach only a sixteenth of the disc', F,
     '            depth: t,\n            needle: needle_width(t, budget),',
     '            depth: t / 16.0,\n            needle: needle_width(t, budget),', BOTH),
    ('the gap term halved', F,
     '        best = best.min(top + slack + g * g / (8.0 * t));',
     '        best = best.min(top + slack + g * g / (16.0 * t));', BOTH),
    ("the supports' depth left out", F,
     '        best = best.min(top + slack + g * g / (8.0 * t));',
     '        best = best.min(slack + g * g / (8.0 * t));', BOTH),
    ('supports on the line left out', F,
     '            if d <= limit {',
     '            if d < limit {', BOTH),
    ('the depth not a bound of its own', F,
     '    let mut best = deepest;',
     '    let mut best = f64::INFINITY;', BOTH),
    ('a refused chord not split', F,
     '            let k = cut.unwrap_or_else(|| split_step(points, chord, deepest));\n            pending.push(Chord {\n                from: (chord.from + k) % n,',
     '            let k = cut.unwrap_or_else(|| split_step(points, chord, deepest));\n            continue;\n            #[allow(unreachable_code)]\n            pending.push(Chord {\n                from: (chord.from + k) % n,', BOTH),
    ('a chord not proved split at its farthest vertex', F,
     '            let k = cut.unwrap_or_else(|| split_step(points, chord, deepest));',
     '            let k = cut.unwrap_or_else(|| split_step(points, chord, deepest || true));', BOTH),
    ('equally far vertices split at the first, not the middle', F,
     '            far.then(from_middle(j).cmp(&from_middle(i)))',
     '            far.then(std::cmp::Ordering::Equal)', BOTH),
    ('a ring the fill would close used', F,
     '        if filled.len() < 3 || crate::zero_area(&filled, linear) {',
     '        if false {', BOTH),
    ('a chord not proved not trimmed', F,
     '            if matches!(verdict, Verdict::Unproved) {',
     '            if false {', BOTH),
    ('the trim without its start edge only', F,
     '[(1, 0), (0, 1), (1, 1)]', '[(0, 1)]', BOTH),
    ('the trim without its end edge only', F,
     '[(1, 0), (0, 1), (1, 1)]', '[(1, 0)]', BOTH),
    ('the trim without both edges', F,
     '[(1, 0), (0, 1), (1, 1)]', '[(1, 0), (0, 1)]', BOTH),
    # #428: needles deeper than the disc.
    ('proposals reach no needle', F,
     '            needle: needle_width(t, budget),',
     '            needle: 0.0,', BOTH),
    ('the needle width halved', F,
     '    2.0 * (2.0 * t * budget + budget * budget).sqrt()',
     '    (2.0 * t * budget + budget * budget).sqrt()', BOTH),
    ('a notch deeper than the disc not cut into pieces', F,
     '    if !deep {\n        return (whole, false);',
     '    if true {\n        return (whole, deep);', BOTH),
    ('pieces bounded by the narrowest', F,
     '.map(|(m, e)| m.max(e));',
     '.map(|(m, e)| m.min(e));', BOTH),
    ('a piece without a bound left out', F,
     '            pieces = pieces\n                .zip(notch_bound(&polygon[from..=to], t).0)\n                .map(|(m, e)| m.max(e));',
     '            pieces = pieces.map(|m| {\n                notch_bound(&polygon[from..=to], t).0.map_or(m, |e| m.max(e))\n            });', BOTH),
    ('a notch bounded by its whole alone', F,
     '        (w, p) => w.or(p),',
     '        (w, p) => w,', UNIT),
    ('a notch cut at every spanned vertex', F,
     '        .filter(|&i| cover::on_segment(a, b, polygon[i]))\n        .chain([last])',
     '        .filter(|&i| i > 0)\n        .chain([last])', BOTH),
    ('a needle wall over the budget split at its farthest vertex', F,
     '                Verdict::TooBig { cut: true } => cut_step(points, chord),',
     '                Verdict::TooBig { cut: true } => None,', BOTH),
    ('a needle wall split at its first vertex on the chord', F,
     '        .min_by_key(|&i| ((2 * i).abs_diff(chord.span), i))',
     '        .min_by_key(|&i| i)', BOTH),
    ('no needle proposed past the reach', C,
     '    if reach.needle.is_nan() || reach.needle <= 0.0 {\n        return false;',
     '    if true {\n        return false;', ALL),
    ('needles not joined across a valley on the chord', C,
     '    if orient(u, w, v) == Sign::Zero && between(u, w, v) {',
     '    if false {', ALL),
    ('a needle of any width proposed', C,
     '    len <= reach.needle && spans.lo >= 0.0 && spans.hi <= len',
     '    spans.lo >= 0.0 && spans.hi <= len', ALL),
    ('a needle leaning past its mouth proposed', C,
     '    len <= reach.needle && spans.lo >= 0.0 && spans.hi <= len',
     '    len <= reach.needle', ALL),
    ('the extent without the vertex extended over', C,
     '        let (mut lo, mut hi) = (along(w), along(w));',
     '        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);', ALL),
    ("the extent without the parts' boxes", C,
     '            if part.lo > part.hi {',
     '            if true {', ALL),
]

probe_lib.main(MUTANTS, timeout=600)
