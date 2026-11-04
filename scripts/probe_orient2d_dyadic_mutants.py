"""Mutation probe for the exact dyadic fallback of `orient2d` (#404).

`orient2d` evaluates the determinant, once its f64 filter declines, over the
original coordinates in a fixed-size dyadic accumulator
(`crates/algorithms/predicates/src/orient2_dyadic.rs`). Each mutant breaks one
term, one sign, the accumulator comparison or one range guard and must turn
`tests/orient2d.rs` red (the scaled collinear/left/right cases, the one-ulp
offset at 1e300, the mixed-scale and subnormal seeded comparison against an
exact dyadic reference) or, for the heal consumer, `tests/exact_intersection.rs`.

The exact path first tries two-term expansions (`orient2d_exact_expansion`,
exact when no coordinate difference rounds), then the dyadic accumulator; the
accumulator is shared with `orient3d` and adds whole shifted words.

Equivalent mutant, deliberately not listed: widening or dropping the UPPER
bound of the filter's difference range (`differences`). An overflowing
product makes the filter's determinant non-finite, which
`Certified::from_filter` already declines. The LOWER bound is load-bearing:
an underflowing product is rounded, so without it the filter proves a wrong
sign (a Negative for a truly Positive triple near 1e-155), and its mutant is
listed.

Runs on `scripts/probe_lib.py`: each mutant in its own process group under a
timeout that kills the group, the source always restored, and a missing or
ambiguous anchor aborts the probe. A mutant that does not compile counts as
BROKEN, not killed.
"""
import probe_lib

DYADIC = "crates/algorithms/predicates/src/orient2_dyadic.rs"
ORIENT = "crates/algorithms/predicates/src/orientation.rs"
DYADIC_ACC = "crates/algorithms/predicates/src/orient3_dyadic.rs"
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
    # the filter's range
    ('the filter range lower bound dropped', ORIENT,
     '(1.0e-90..=1.0e90).contains(&v.abs())', '(0.0..=1.0e90).contains(&v.abs())', PREDICATES),
    # the expansion stage
    ('the main-products-only stage used although a tail is nonzero', ORIENT,
     '    if tails == [0.0; 4] {', '    if tails != [0.0; 4] || true {', PREDICATES),
    ('the right low product term added instead of subtracted', ORIENT,
     '[left_low, -right_low, left_high, -right_high]', '[left_low, right_low, left_high, -right_high]', PREDICATES),
    ('the right high product term added instead of subtracted', ORIENT,
     '[left_low, -right_low, left_high, -right_high]', '[left_low, -right_low, left_high, right_high]', PREDICATES),
    ('the left low product term dropped', ORIENT,
     '[left_low, -right_low, left_high, -right_high]', '[0.0, -right_low, left_high, -right_high]', PREDICATES),
    ('the sum read from its smallest component', ORIENT,
     'expansion[..len].iter().rev().find(|v| **v != 0.0)', 'expansion[..len].iter().find(|v| **v != 0.0)', PREDICATES),
    ('a left tail product (acx * tbcy) subtracted', ORIENT,
     '(acx, tbcy, false),', '(acx, tbcy, true),', PREDICATES),
    ('a right tail product (acy * tbcx) added', ORIENT,
     '(acy, tbcx, true),', '(acy, tbcx, false),', PREDICATES),
    ('a left tail product (tacx * bcy) dropped', ORIENT,
     '        (tacx, bcy, false),\n', '', PREDICATES),
    ('a right tail product (tacy * bcx) dropped', ORIENT,
     '        (tacy, bcx, true),\n', '', PREDICATES),
    ('the tail-by-tail product (tacx * tbcy) subtracted', ORIENT,
     '(tacx, tbcy, false),', '(tacx, tbcy, true),', PREDICATES),
    ('the tail-by-tail product (tacy * tbcx) dropped', ORIENT,
     '        (tacy, tbcx, true),\n', '', PREDICATES),
    ('the low part of a tail product dropped', ORIENT,
     'terms[count + 1] = sign * low;', 'terms[count + 1] = 0.0;', PREDICATES),
    ('the smallest accepted tail dropped to zero', ORIENT,
     'const MIN_TAIL: f64 = 1.0e-120;', 'const MIN_TAIL: f64 = 0.0;', PREDICATES),
    # the shared accumulator
    ('the carry out of a limb dropped', DYADIC_ACC,
     'value = u64::from(carry);', 'value = 0;', PREDICATES),
    ('the bits shifted into the next limb dropped', DYADIC_ACC,
     'if shift != 0 && !self.add_limb(', 'if false && !self.add_limb(', PREDICATES),
    ('the word index of a product ignored', DYADIC_ACC,
     'self.add_limb(base + word_index, word << shift)', 'self.add_limb(base, word << shift)', PREDICATES),
    # routing in orient2d
    ('a failed dyadic evaluation read as certain zero', ORIENT,
     '.or_else(|| orient2d_exact_dyadic(a, b, c))',
     '.or_else(|| orient2d_exact_dyadic(a, b, c).or(Some(axiolid_guarantees::Sign::Zero)))', PREDICATES),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
