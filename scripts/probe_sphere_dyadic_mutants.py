"""Mutation probe for the range gates and exact fallback of `incircle` and `insphere` (#430).

Both predicates are lifted determinants of degree 4 and 5, so their products
leave binary64's range near 1e60 and 1e-60. The filter declines on its
permanent outside [1e-250, 1e300]; the expansion stage runs only where every
difference and rounding tail is in range; every other finite input goes to the
fixed-size dyadic evaluation in
`crates/algorithms/predicates/src/sphere_dyadic.rs`. Each mutant breaks one
gate, one term, one sign, the accumulator or a routing decision and must turn
`tests/sphere_extreme_scales.rs` red (the issue's inside and outside inputs at
1e+-80 to 1e+-300, exact zeros on power-of-two scaled lattices, the seeded
comparison with an exact dyadic reference over scales 2^-1060 to 2^1010, the
non-finite inputs) or, for the ordinary range, `tests/incircle_exact.rs`.

Equivalent mutants, deliberately not listed: raising the filter's upper
permanent bound to infinity (a permanent beyond 1e300 below infinity cannot
overflow, and an infinite one makes the error bound infinite, which certifies
nothing); reading a NaN as in range in `all_in_range` (the tail check rejects
the NaN tail of the same difference); a non-finite pre-check in `insphere`
(`insphere_small_integer`, the expansion guard and the dyadic evaluation each
refuse it); a uniform extra shift of every product's offset (it scales both
sums alike).

Runs on `scripts/probe_lib.py`: each mutant in its own process group under a
timeout that kills the group, the source always restored, and a missing or
ambiguous anchor aborts the probe. A mutant that does not compile counts as
BROKEN, not killed.
"""
import probe_lib

SPHERE = "crates/algorithms/predicates/src/sphere.rs"
DYADIC = "crates/algorithms/predicates/src/sphere_dyadic.rs"
EXTREME = [["-p", "axiolid-predicates", "--test", "sphere_extreme_scales"]]
ORDINARY = [["-p", "axiolid-predicates", "--test", "incircle_exact"]]
# A mutant that breaks the permutation walk can loop for ever; the timeout kills it.
TIMEOUT = 300

MUTANTS = [
    # the filter's range
    ('the filter accepts a permanent below the underflow range', SPHERE,
     'const PERMANENT_MIN: f64 = 1.0e-250;', 'const PERMANENT_MIN: f64 = 0.0;', EXTREME),
    ('the incircle filter ignores the range', SPHERE,
     '    if !permanent_in_range(permanent) {\n        return Certified::Uncertain {\n            attempted: Precision::F64,\n        };\n    }\n    Certified::from_filter(\n        determinant,\n        INCIRCLE_ERROR_FACTOR',
     '    Certified::from_filter(\n        determinant,\n        INCIRCLE_ERROR_FACTOR', EXTREME),
    ('the insphere filter ignores the range', SPHERE,
     '    if !permanent_in_range(permanent) {\n        return Certified::Uncertain {\n            attempted: Precision::F64,\n        };\n    }\n    Certified::from_filter(\n        determinant,\n        INSPHERE_ERROR_FACTOR',
     '    Certified::from_filter(\n        determinant,\n        INSPHERE_ERROR_FACTOR', EXTREME),
    # the expansion stage's range
    ('the incircle expansion range lower bound dropped', SPHERE,
     'const INCIRCLE_MIN: f64 = 1.0e-60;', 'const INCIRCLE_MIN: f64 = 0.0;', EXTREME),
    ('the incircle expansion range upper bound dropped', SPHERE,
     'const INCIRCLE_MAX: f64 = 1.0e60;', 'const INCIRCLE_MAX: f64 = f64::INFINITY;', EXTREME),
    ('the insphere expansion range lower bound dropped', SPHERE,
     'const INSPHERE_MIN: f64 = 1.0e-45;', 'const INSPHERE_MIN: f64 = 0.0;', EXTREME),
    ('the insphere expansion range upper bound dropped', SPHERE,
     'const INSPHERE_MAX: f64 = 1.0e45;', 'const INSPHERE_MAX: f64 = f64::INFINITY;', EXTREME),
    ('a rounding tail below the range accepted', SPHERE,
     '(tail == 0.0 || tail.abs() >= min)', '(true || tail.abs() >= min)', EXTREME),
    ('a range check on the head dropped', SPHERE,
     'all_in_range(&[head], min, max) &&', 'true &&', EXTREME),
    # routing
    ('incircle always takes the expansion stage', SPHERE,
     'if differences_expandable(&pairs, INCIRCLE_MIN, INCIRCLE_MAX) {', 'if true {', EXTREME),
    ('insphere always takes the expansion stage', SPHERE,
     'if differences_expandable(pairs.as_flattened(), INSPHERE_MIN, INSPHERE_MAX) {', 'if true {', EXTREME),
    ('a failed incircle evaluation read as certain zero', SPHERE,
     '        _ => incircle_exact(a, b, c, d).map_or(\n            Certified::Uncertain {\n                attempted: Precision::Exact,\n            },\n            Certified::exact_sign,\n        ),',
     '        _ => Certified::exact_sign(incircle_exact(a, b, c, d).unwrap_or(Sign::Zero)),', EXTREME),
    ('a failed insphere evaluation read as certain zero', SPHERE,
     '        _ => insphere_exact(a, b, c, d, e).map_or(\n            Certified::Uncertain {\n                attempted: Precision::Exact,\n            },\n            Certified::exact_sign,\n        ),',
     '        _ => Certified::exact_sign(insphere_exact(a, b, c, d, e).unwrap_or(Sign::Zero)),', EXTREME),
    # the dyadic evaluation
    ('the incircle accumulator too small for the largest products', DYADIC,
     'const INCIRCLE_LIMBS: usize = 133;', 'const INCIRCLE_LIMBS: usize = 120;', EXTREME),
    ('the insphere accumulator too small for the largest products', DYADIC,
     'const INSPHERE_LIMBS: usize = 166;', 'const INSPHERE_LIMBS: usize = 150;', EXTREME),
    ('the smallest factor exponent raised', DYADIC,
     'const MIN_FACTOR_EXPONENT: i32 = -1074;', 'const MIN_FACTOR_EXPONENT: i32 = -1000;', EXTREME),
    ('the permutation sign inverted', DYADIC,
     'let mut is_negative = inversions % 2 == 1;', 'let mut is_negative = inversions % 2 == 0;', EXTREME),
    ('a zero factor not vanishing the product', DYADIC,
     '                        vanishes = true;\n                        continue;', '                        continue;', EXTREME),
    ('the sign of a negative coordinate ignored', DYADIC,
     'is_negative ^= value_negative ^ from_origin;', 'is_negative ^= from_origin;', EXTREME),
    ('the sign of the reference point ignored', DYADIC,
     'is_negative ^= value_negative ^ from_origin;', 'is_negative ^= value_negative;', EXTREME),
    ('the lift squares only the first axis', DYADIC,
     'for lift_axis in 0..D {', 'for lift_axis in 0..1 {', EXTREME),
    ('the lift squares two different axes', DYADIC,
     'factors[count + 1] = (lift_row, lift_axis);', 'factors[count + 1] = (lift_row, (lift_axis + 1) % D);', EXTREME),
    ('one choice of point or reference dropped', DYADIC,
     'for choices in 0_u32..(1 << FACTORS) - 1 {', 'for choices in 0_u32..(1 << FACTORS) - 2 {', EXTREME),
    ('the carry between significand words dropped', DYADIC,
     '        carry = wide >> 64;', '        carry = 0;', EXTREME),
    ('the positive and negative sums swapped in the comparison', DYADIC,
     'match positive.compare(&negative) {', 'match negative.compare(&positive) {', EXTREME),
    ('an equal pair of sums read as positive', DYADIC,
     '        std::cmp::Ordering::Equal => Sign::Zero,', '        std::cmp::Ordering::Equal => Sign::Positive,', EXTREME),
    ('the permutations stop after the first', DYADIC,
     '        if !next_permutation(columns) {', '        if true || !next_permutation(columns) {', EXTREME),
    ('the pivot swap of the next permutation dropped', DYADIC,
     '    items.swap(pivot, swap);', '', EXTREME),
    ('the tail of the next permutation not reversed', DYADIC,
     '    items[pivot + 1..].reverse();', '', EXTREME),
    # the ordinary range still reads the filter
    ('the incircle error bound zeroed', SPHERE,
     'INCIRCLE_ERROR_FACTOR * permanent', '0.0 * permanent', ORDINARY),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
