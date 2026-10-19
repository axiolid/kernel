"""Mutation probe for the fast disc morphology (#330,
`crates/algorithms/planar/overlay/src/fill.rs` and the mode switch in
`minkowski.rs`).

Each mutant weakens one step -- the mode switch either way, the side of
the set whose notches are filled, each condition on a chord (its exact
proof, the depth within the disc, the budget, no ring inside its notch),
the excess (the gap term, the supports, the depth bound), splitting a
refused chord, the ring a fill would close, the inscribed polygon drawn in
by the excess and the deviation carrying it, and the budget's refusal --
and must turn a test red.

Equivalent mutants, deliberately not listed:
- the `f64` slack on the excess (64 ulps of the largest coordinate) and its
  final widening: they cover rounding of the excess's own arithmetic, far
  below what any sampled test point can resolve;
- splitting a chord too deep or over the budget at another vertex than
  the farthest: every split ends in single edges, so only which shorter
  chords are found changes, and with it the cost and the excess, both
  within the stated bounds (the split of a chord not proved, which
  decides whether a wall of teeth is filled at all, is listed);
- the order of the halving limits over the supports: the bound is a
  minimum over all of them.
"""
import probe_lib

F = "crates/algorithms/planar/overlay/src/fill.rs"
M = "crates/algorithms/planar/overlay/src/minkowski.rs"
UNIT = ["-p", "axiolid-overlay", "--lib", "fill::"]
ROUTE = ["-p", "axiolid-overlay", "--test", "minkowski"]
BOTH = [UNIT, ROUTE]

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
     '            let taken = if cover::proved(r, &boundary, chord, outside, &index) {',
     '            let taken = if true {', BOTH),
    ('a chord over the budget used', F,
     '                    Some(e) if e <= budget => {',
     '                    Some(e) => {', BOTH),
    ('a notch holding another ring filled', F,
     '                        if encloses_none(r, &polygon, &index) {',
     '                        if true {', BOTH),
    ('a ring inside the notch not found', F,
     '        !inside(polygon, index.rings[r][e])',
     '        true', BOTH),
    ('a notch deeper than the disc taken', F,
     '    if deepest.is_nan() || deepest > t {',
     '    if deepest.is_nan() {', BOTH),
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
     '            let k = split_step(points, chord, deepest);\n            pending.push(Chord {\n                from: (chord.from + k) % n,',
     '            let k = split_step(points, chord, deepest);\n            continue;\n            #[allow(unreachable_code)]\n            pending.push(Chord {\n                from: (chord.from + k) % n,', BOTH),
    ('a chord not proved split at its farthest vertex', F,
     '            let k = split_step(points, chord, deepest);',
     '            let k = split_step(points, chord, deepest || true);', BOTH),
    ('equally far vertices split at the first, not the middle', F,
     '            far.then(from_middle(j).cmp(&from_middle(i)))',
     '            far.then(std::cmp::Ordering::Equal)', BOTH),
    ('a ring the fill would close used', F,
     '        if filled.len() < 3 || crate::zero_area(&filled, linear) {',
     '        if false {', BOTH),
]

probe_lib.main(MUTANTS, timeout=600)
