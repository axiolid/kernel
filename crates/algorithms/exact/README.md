# axiolid-exact

Filtered exact arithmetic for *constructions* over `f64` input (ADR 0068):
signs and orderings of values `f64` cannot hold, such as where two segments
cross, where a line meets a circle, or roots of the form
`(a + b*sqrt(c)) / d`. Each sign question is one expression evaluated first
in outward-rounded interval arithmetic and, only if that cannot decide, in
exact dyadic big-integer arithmetic (`num-bigint`). It has no division and
evaluates no square roots; approximate values are available for output
only.

```bash
cargo add axiolid-exact
```

- API documentation: [docs.rs/axiolid-exact](https://docs.rs/axiolid-exact)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-exact)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

Choosing between this crate and `axiolid-predicates`: if the question is the
sign of a polynomial in the *input* coordinates (orientation, in-circle),
use `axiolid-predicates`, which needs no big integers. If the question is
about a *constructed* point (a crossing, a line/circle hit), use this crate.
