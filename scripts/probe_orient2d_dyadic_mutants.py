"""Mutation probe for the exact dyadic fallback of `orient2d` (#404).

`orient2d` evaluates the determinant, once its f64 filter declines, over the
original coordinates in a fixed-size dyadic accumulator
(`crates/algorithms/predicates/src/orient2_dyadic.rs`). Each mutant breaks one
term, one sign, the accumulator comparison or one range guard and must turn
`tests/orient2d.rs` red (the scaled collinear/left/right cases, the one-ulp
offset at 1e300, the mixed-scale and subnormal seeded comparison against an
exact dyadic reference) or, for the heal consumer, `tests/exact_intersection.rs`.

Equivalent mutants, deliberately not listed: widening or dropping the
filter's difference range (`differences`). An overflowing product makes the
filter's determinant non-finite, which `Certified::from_filter` already
declines, and an underflowing product is rounded monotonically, so the
filter can lose a sign to zero but never flip one. The range stays as the
conservative twin of `orient3d_filter`'s, so the filter's bound is only ever
asked about products in the normal range.

Runs on `scripts/probe_lib.py`: each mutant in its own process group under a
timeout that kills the group, the source always restored, and a missing or
ambiguous anchor aborts the probe. A mutant that does not compile counts as
BROKEN, not killed.
"""
import probe_lib

DYADIC = "crates/algorithms/predicates/src/orient2_dyadic.rs"
ORIENT = "crates/algorithms/predicates/src/orientation.rs"
PREDICATES = [["-p", "axiolid-predicates", "--test", "orient2d"]]
HEAL = [["-p", "axiolid-heal", "--test", "exact_intersection"]]
TIMEOUT = 1800

MUTANTS = [
    # the six terms of the expanded determinant
    ('the ax*cy term added instead of subtracted', DYADIC,
     '(true, a.x, c.y),', '(false, a.x, c.y),', PREDICATES),
    ('the cx*by term added instead of subtracted', DYADIC,
     '(true, c.x, b.y),', '(false, c.x, b.y),', PREDICATES),
    ('the ay*bx term added instead of subtracted', DYADIC,
     '(true, a.y, b.x),', '(false, a.y, b.x),', PREDICATES),
    ('the ay*cx term subtracted instead of added', DYADIC,
     '(false, a.y, c.x),', '(true, a.y, c.x),', PREDICATES),
    ('the cy*bx term dropped', DYADIC,
     '        (false, c.y, b.x),\n', '', PREDICATES),
    ('the ax*by term dropped', DYADIC,
     '        (false, a.x, b.y),\n', '', PREDICATES),
    # accumulation
    ('the sign of a negative factor ignored', DYADIC,
     'subtract ^ first_negative ^ second_negative', 'subtract ^ first_negative', PREDICATES),
    ('positive and negative sums swapped in the comparison', DYADIC,
     'match positive.compare(&negative) {', 'match negative.compare(&positive) {', PREDICATES),
    ('an equal pair of sums read as positive', DYADIC,
     '        std::cmp::Ordering::Equal => Sign::Zero,', '        std::cmp::Ordering::Equal => Sign::Positive,', PREDICATES),
    ('the high word of the product dropped', DYADIC,
     '[product as u64, (product >> 64) as u64, 0]', '[product as u64, 0, 0]', PREDICATES),
    ('the smallest product exponent raised', DYADIC,
     'const MIN_PAIR_EXPONENT: i32 = -2148;', 'const MIN_PAIR_EXPONENT: i32 = -2100;', PREDICATES),
    # routing in orient2d
    ('a failed dyadic evaluation read as certain zero', ORIENT,
     '        _ => orient2d_exact_dyadic(a, b, c).map_or(',
     '        _ => orient2d_exact_dyadic(a, b, c).or(Some(Sign::Zero)).map_or(', PREDICATES),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
