# axiolid-exact instructions

Purpose: exact *constructions* over `f64` input (ADR 0068). Signs of values
`f64` cannot hold -- where segments cross, where a line meets a circle --
decided in two passes over one expression: outward-rounded intervals, then
dyadic big integers (`num-bigint`) only if the interval cannot decide.

Allowed internal dependencies: `axiolid-core`, `axiolid-guarantees`. The one
external is `num-bigint`, kept behind `Dyadic` so it can be swapped. Keep
`axiolid-predicates` free of it: consumers who need only certified signs of
input points must not pay for big integers.

## Module ownership

`arith.rs` the shared `Arith` trait; `interval.rs` fast tier; `dyadic.rs`
exact tier; `certify.rs` the two-tier driver and `ExactError`; `root.rs`
`(a + b*sqrt(c)) / d` signs and comparisons; `construct.rs` public
constructions (crossings, line/circle hits, ordering along a line).

## Invariants

- One expression, two tiers. Write sign questions as `SignExpr` against
  `Arith`, never as a separate f64 fast path, so the filter and the exact
  fallback cannot evaluate different polynomials.
- No division, no square roots evaluated. Clear denominators; decide
  root signs by squaring with case analysis (`root.rs`).
- An interval may say "undecided"; it must never be wrong. A decided
  interval sign equals the exact sign (property-tested over all f64
  regimes, including subnormals and overflow).
- `approx_*` methods are for output only. Decisions use sign questions.
- Structural answers beat evaluation where they are provable: the two hits
  of one circle are ordered by branch, not by computing `x - x`, which no
  interval can certify as zero.

## Gates

```bash
cargo test -p axiolid-exact
cargo bench -p axiolid-exact --bench exact   # escalation rate beside ns/call
```
